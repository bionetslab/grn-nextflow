#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------- Params ----------------
params.in_h5ad = null
params.outdir  = 'results'
params.threads = 4
params.seed    = 42


// ---------------- Workflow ----------------
workflow {
    h5ad_ch = Channel.fromPath(params.in_h5ad)

    converted_ch = convert_h5ad_to_scribe(h5ad_ch)
    scribe_ch    = run_scribe(converted_ch)
    aggregate_edges(scribe_ch)
}

// ---------------- Process: convert_h5ad_to_scribe ----------------
process convert_h5ad_to_scribe {
    tag "h5ad→SCRIBE inputs"
    publishDir "${params.outdir}/inputs", mode: 'copy'
    cpus 1
    memory '4 GB'
    stageInMode 'copy'

    input:
    path h5ad

    output:
    tuple path("expr.csv"), path("genes.txt"), path("pseudotime.csv")

    script:
    """
    python ${projectDir}/modules/grn/scribe/h5ad_to_scribe_inputs.py \
        --h5ad "${h5ad}" \
        --expr expr.csv \
        --genes genes.txt \
        --pseudotime pseudotime.csv
    """
}

// ---------------- Process: run_scribe ----------------

process run_scribe {
    tag "Run SCRIBE"
    publishDir "${params.outdir}/results", mode: 'copy'
    cpus params.threads

    input:
    tuple path(expr), path(genes), path(pseudotime)

    output:
    tuple path("edges.tsv"), path("A_averaged.tsv")

    script:
    """
    set -euo pipefail

    if command -v Rscript >/dev/null 2>&1 ; then
      echo "[INFO] Using Rscript backend"
      Rscript ${projectDir}/modules/grn/scribe/run_scribe.R \
          --expr "${expr}" \
          --genes "${genes}" \
          --pseudotime "${pseudotime}" \
          --threads ${task.cpus} \
          --seed ${params.seed} \
          --edges edges.tsv

      Rscript ${projectDir}/modules/grn/scribe/averageA.R \
          --edges edges.tsv \
          --out A_averaged.tsv

    else
      echo "[WARN] Rscript not found — using Python fallback"

      # -------- Python fallback for SCRIBE-like scoring --------
      python - <<'PY_SCRIBE'
import numpy as np, pandas as pd
from math import sqrt
import sys

expr_path = "expr.csv"
genes_path = "genes.txt"
pt_path = "pseudotime.csv"

X = pd.read_csv(expr_path, header=None).values  # cells x genes
genes = [g.strip() for g in open(genes_path)]
pt = pd.read_csv(pt_path, header=None).iloc[:,0].values.astype(float)

if X.shape[1] != len(genes):
    raise SystemExit(f"[ERROR] expr columns ({X.shape[1]}) != genes ({len(genes)})")

# z-score by gene
X = (X - X.mean(axis=0, keepdims=True))
std = X.std(axis=0, ddof=0, keepdims=True)
std[std==0] = 1.0
X = X / std

# order by pseudotime
ord_idx = np.argsort(pt)
X = X[ord_idx, :]
pt = pt[ord_idx]

n, p = X.shape
lag = 1 if n > 6 else 0

def wcorr(y, x, w):
    wy = y - np.average(y, weights=w)
    wx = x - np.average(x, weights=w)
    num = np.sum(w*wx*wy)
    den = sqrt(np.sum(w*wx*wx) * np.sum(w*wy*wy))
    return abs(num/den) if den>0 else 0.0

edges = []
if lag == 0:
    # contemporaneous, simple weighting = 1
    w = np.ones(n, dtype=float)
    for ti in range(p):
        y = X[:,ti]
        for j in range(p):
            if j==ti: continue
            sc = wcorr(y, X[:,j], w)
            if sc>0:
                edges.append((genes[j], genes[ti], sc))
else:
    Tn = n - lag
    w = np.ones(Tn, dtype=float)  # simple weights; keep fallback minimal
    for ti in range(p):
        y = X[lag:, ti]
        for j in range(p):
            if j==ti: continue
            xj = X[:-lag, j]
            sc = wcorr(y, xj, w)
            if sc>0:
                edges.append((genes[j], genes[ti], sc))

df = pd.DataFrame(edges, columns=["regulator","target","score"]).sort_values("score", ascending=False)
df.to_csv("edges.tsv", sep="\\t", index=False)
print("[OK] Wrote edges.tsv with", len(df), "edges")
PY_SCRIBE

      # -------- Python fallback for averageA (normalize to max) --------
      python - <<'PY_AVG'
import pandas as pd, numpy as np
df = pd.read_csv("edges.tsv", sep="\\t")
if df.empty:
    pd.DataFrame(columns=["regulator","target","weight"]).to_csv("A_averaged.tsv", sep="\\t", index=False)
else:
    # take max weight per (regulator,target), normalize by global max
    df['weight'] = df['score']
    df = df.groupby(['regulator','target'], as_index=False)['weight'].max()
    m = df['weight'].max()
    if m>0: df['weight'] = df['weight'] / m
    df.sort_values('weight', ascending=False).to_csv("A_averaged.tsv", sep="\\t", index=False)
print("[OK] Wrote A_averaged.tsv")
PY_AVG
    fi
    """
}

// ---------------- Process: aggregate_edges ----------------
process aggregate_edges {
    tag "Aggregate edges"
    publishDir "${params.outdir}/final", mode: 'copy'
    cpus 1
    memory '2 GB'

    input:
    tuple path(edges), path(averaged)

    output:
    path "scribe_ranked_edges.tsv"

    script:
    """
    python ${projectDir}/modules/grn/scribe/aggregate_edges.py \
        --in "${edges}" \
        --avg "${averaged}" \
        --out scribe_ranked_edges.tsv
    """
}

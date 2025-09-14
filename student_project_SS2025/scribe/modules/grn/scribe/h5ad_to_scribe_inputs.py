#!/usr/bin/env python3

import argparse
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import anndata as ad


def read_h5ad_minimal(h5ad_path):
    """Read minimal data from h5ad: expression matrix, gene names, pseudotime."""
    # Use anndata to handle variations in storage
    adata = ad.read_h5ad(h5ad_path)

    # Expression matrix
    if sparse.issparse(adata.X):
        mat = adata.X.toarray()
    else:
        mat = np.array(adata.X)

    # Gene names (robust fallback)
    if "_index" in adata.var.columns:
        genes = adata.var["_index"].astype(str).tolist()
    elif "index" in adata.var.columns:
        genes = adata.var["index"].astype(str).tolist()
    elif "gene_symbols" in adata.var.columns:
        genes = adata.var["gene_symbols"].astype(str).tolist()
    else:
        genes = adata.var.index.astype(str).tolist()

    # Pseudotime
    pt_candidates = ["pseudotime", "Pseudotime", "pt"]
    pt = None
    for col in pt_candidates:
        if col in adata.obs.columns:
            pt = adata.obs[col].to_numpy()
            break

    if pt is None:
        raise ValueError(
            f"No pseudotime column found in {pt_candidates}. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    return mat, genes, pt


def main():
    parser = argparse.ArgumentParser(description="Convert .h5ad to SCRIBE input files")
    parser.add_argument("--h5ad", required=True, help="Path to input .h5ad file")
    parser.add_argument(
        "--expr", required=True, help="Output CSV for expression matrix"
    )
    parser.add_argument("--genes", required=True, help="Output TXT for gene names")
    parser.add_argument(
        "--pseudotime", required=True, help="Output CSV for pseudotime values"
    )
    args = parser.parse_args()

    mat, genes, pt = read_h5ad_minimal(args.h5ad)

    # Save expression matrix
    pd.DataFrame(mat).to_csv(args.expr, index=False, header=False)

    # Save genes
    with open(args.genes, "w") as f:
        for g in genes:
            f.write(str(g) + "\n")

    # Save pseudotime
    pd.DataFrame(pt).to_csv(args.pseudotime, index=False, header=False)


if __name__ == "__main__":
    main()

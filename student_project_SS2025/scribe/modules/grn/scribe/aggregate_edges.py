#!/usr/bin/env python3
import argparse, pandas as pd, numpy as np, os

p = argparse.ArgumentParser()
p.add_argument("--in", dest="inp", required=True)
p.add_argument("--avg", dest="avg", required=True)
p.add_argument("--out", dest="out", required=True)
a = p.parse_args()


def read_df(path, sep="\t"):
    if os.path.exists(path) and os.path.getsize(path) > 0:
        try:
            return pd.read_csv(path, sep=sep)
        except Exception:
            pass
    return pd.DataFrame()


edges = read_df(a.inp, sep="\t")
avg = read_df(a.avg, sep="\t")

if not edges.empty and not avg.empty:
    df = edges.merge(avg, on=["regulator", "target"], how="outer")
    s = pd.to_numeric(
        df.get("score", pd.Series([0] * len(df))), errors="coerce"
    ).fillna(0.0)
    w = pd.to_numeric(
        df.get("weight", pd.Series([0] * len(df))), errors="coerce"
    ).fillna(0.0)
    if s.max() > 0:
        s = s / s.max()
    df["final_score"] = np.maximum(s.values, w.values)
    df = df[["regulator", "target", "final_score"]]
elif not avg.empty:
    df = avg.rename(columns={"weight": "final_score"})[
        ["regulator", "target", "final_score"]
    ]
else:
    df = edges.rename(columns={"score": "final_score"})[
        ["regulator", "target", "final_score"]
    ]

df = df[df["regulator"] != df["target"]].copy()
df = df.sort_values("final_score", ascending=False).drop_duplicates(
    ["regulator", "target"]
)
df["rank"] = np.arange(1, len(df) + 1)
df.to_csv(a.out, sep="\t", index=False)
print("[OK] Wrote", a.out)

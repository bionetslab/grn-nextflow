#!/usr/bin/env python3
import argparse, sys, os, csv


def has_h5py():
    try:
        import h5py  # noqa

        return True
    except Exception:
        return False


def read_h5ad_minimal(h5ad_path):
    import h5py
    import numpy as np

    with h5py.File(h5ad_path, "r") as f:
        # Try CSR
        if "X" in f and all(
            k in f["X"] for k in ["data", "indptr", "indices", "shape"]
        ):
            data = f["X"]["data"][()]
            indptr = f["X"]["indptr"][()]
            indices = f["X"]["indices"][()]
            shape = tuple(f["X"]["shape"][()])
            dense = np.zeros(shape, dtype=float)
            for i in range(shape[0]):
                start, end = indptr[i], indptr[i + 1]
                cols = indices[start:end]
                vals = data[start:end]
                dense[i, cols] = vals
            mat = dense
        elif "X" in f:
            mat = f["X"][()]
        else:
            raise RuntimeError("Cannot find X in .h5ad")

        # var names
        if "var" in f and "_index" in f["var"]:
            genes = [
                g.decode("utf-8") if isinstance(g, bytes) else str(g)
                for g in f["var"]["_index"][()]
            ]
        elif "var_names" in f:
            genes = [
                g.decode("utf-8") if isinstance(g, bytes) else str(g)
                for g in f["var_names"][()]
            ]
        else:
            genes = [f"g{i}" for i in range(mat.shape[1])]

        # pseudotime
        pt = None
        if "obs" in f:
            for key in [b"pseudotime", b"pt", b"dpt_pseudotime", b"latent_time"]:
                if key in f["obs"]:
                    col = f["obs"][key][()]
                    try:
                        pt = col.astype(float)
                    except Exception:
                        try:
                            pt = col[:, 0].astype(float)
                        except Exception:
                            pass
                    if pt is not None:
                        break
        if pt is None:
            pt = list(range(mat.shape[0]))
        return mat, genes, pt


def write_csv_matrix(path, matrix, header):
    import numpy as np

    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for i in range(matrix.shape[0]):
            writer.writerow(np.asarray(matrix[i, :]).tolist())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=False)
    ap.add_argument("--expr", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--pseudotime", required=True)
    args = ap.parse_args()

    # passthrough if CSVs already exist
    if (
        os.path.exists("expr.csv")
        and os.path.exists("genes.txt")
        and os.path.exists("pseudotime.csv")
    ):
        print("[INFO] Found existing CSV inputs; skipping .h5ad conversion.")
        sys.exit(0)

    if args.h5ad is None or not os.path.exists(args.h5ad):
        print(
            "[ERROR] No .h5ad found and no CSVs present. Provide filtered_placeholder.h5ad or expr.csv/genes.txt/pseudotime.csv.",
            file=sys.stderr,
        )
        sys.exit(2)

    if not has_h5py():
        print(
            "[ERROR] Python package 'h5py' not found. Install with: pip install h5py",
            file=sys.stderr,
        )
        print(
            "[TIP] Or provide CSVs manually: expr.csv, genes.txt, pseudotime.csv",
            file=sys.stderr,
        )
        sys.exit(3)

    mat, genes, pt = read_h5ad_minimal(args.h5ad)

    import numpy as np

    ord_idx = np.argsort(np.asarray(pt).astype(float))
    mat = np.asarray(mat)[ord_idx, :]
    pt = np.asarray(pt)[ord_idx].astype(float)

    write_csv_matrix(args.expr, mat, genes)
    with open(args.genes, "w", newline="") as gfh:
        for g in genes:
            gfh.write(f"{g}\n")
    with open(args.pseudotime, "w", newline="") as pfh:
        w = csv.writer(pfh)
        w.writerow(["pseudotime"])
        for v in pt:
            w.writerow([float(v)])
    print("[OK] Wrote expr.csv, genes.txt, pseudotime.csv")


if __name__ == "__main__":
    main()

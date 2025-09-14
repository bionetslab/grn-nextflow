import argparse
import scanpy as sc
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help="Input h5ad file")
parser.add_argument('--expr_out', default='expression.txt')
parser.add_argument('--time_out', default='pseudotime.txt')
args = parser.parse_args()

# Load data
adata = sc.read_h5ad(args.input)

# Save expression matrix (genes x cells)
pd.DataFrame(adata.X.T).to_csv(args.expr_out, sep='\t', header=False, index=False)

# Save pseudotime (column index + value)
if 'pseudotime' in adata.obs.columns:
    pseudotime = adata.obs['pseudotime'].reset_index()
else:
    raise ValueError("No 'pseudotime' found in adata.obs")
pseudotime.to_csv(args.time_out, sep='\t', header=False, index=False)

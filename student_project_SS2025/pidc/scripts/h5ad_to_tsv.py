import scanpy as sc
import pandas as pd
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

adata = sc.read_h5ad(input_file)
X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

# Transpose to get genes as rows
df = pd.DataFrame(X.T, index=adata.var_names, columns=adata.obs_names)

# Save with gene names in first column (index=True), and sample names as headers
df.to_csv(output_file, sep="\t", index=True)

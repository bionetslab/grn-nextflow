import scanpy as sc
import numpy as np

# Load the file
adata = sc.read_h5ad("input/filtered_placeholder.h5ad")

# Create pseudotime: here just a linear progression
adata.obs["pseudotime"] = np.linspace(0, 1, adata.n_obs)

# Save it back
adata.write("input/filtered_placeholder.h5ad")

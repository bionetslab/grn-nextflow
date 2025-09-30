#!/usr/bin/env python3

import anndata as ad
import pandas as pd
import os
import sys


filename = "${csv}"
prefix = "${prefix}"
df = pd.read_csv(f"{filename}", index_col=0)
df = df.T  # Transpose the DataFrame to have genes as columns and cells as rows
adata = ad.AnnData(df)
adata.write_h5ad(f"{prefix}.h5ad")
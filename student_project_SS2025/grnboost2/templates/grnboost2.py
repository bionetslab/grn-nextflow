#!/usr/bin/env python3

import os
import platform
import json
import base64

os.environ[ 'NUMBA_CACHE_DIR' ] = '/tmp/'
# os.environ['MPLCONFIGDIR'] = "/tmp/"

import scanpy as sc
import resource
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names

from threadpoolctl import threadpool_limits
threadpool_limits(int("${task.cpus}"))
sc.settings.n_jobs = int("${task.cpus}")

from distributed import Client, LocalCluster

adata = sc.read_h5ad("${sample_h5ad}")
prefix = "${prefix}"

expr_matrix = adata.to_df()

if adata.n_vars > 15000:
    adata.layers["raw"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=15000)
    expr_matrix = adata[:, adata.var.highly_variable].to_df(layer="${layer}")

if "${params.grn_tflist}" != 'all':
    tf_names = load_tf_names("${params.grn_tflist}")
else:
    tf_names = 'all'

def main(expr_matrix, tf_names, prefix):
    client = Client(LocalCluster())
    network = grnboost2(
        expression_data=expr_matrix.to_numpy(), 
        gene_names=expr_matrix.columns, 
        tf_names=tf_names, 
        verbose=False,
        client_or_address=client
    )  

    network.to_csv(f"{prefix}_grnboost2.csv", sep='\t', index=False)
    

if __name__ == '__main__':
    main(expr_matrix, tf_names, prefix)
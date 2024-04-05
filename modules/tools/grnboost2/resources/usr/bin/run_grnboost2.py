#!/usr/bin/env python
"""Usage: run_grnboost2.py -c FILE -d FILE -o FILE [-u USE_TF_LIST] [-w PROJECT_DIRECTORY]

Runs grnboost2 (https://github.com/tmoerman/arboreto/tree/master)

-h --help               show this
-c --inputfile1 FILE    specify input file
-d --inputfile2 FILE    specify input file
-o --output FILE        specify path and file name
-u --use_tf_list USE_TF_LIST Use transcription factor list?  True or False.  [default: True]
-w --project_directory PROJECT_DIRECTORY path to the nextflow directory
--verbose               print more text

"""

# import python modules
import pandas as pd
from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from docopt import docopt
import os

def run_grnboost2(file1, file2, output_folder, use_tf_list, project_dir):
    ex_matrix_cond1 = pd.read_csv(file1, sep='\t', index_col = 0).T
    ex_matrix_cond2 = pd.read_csv(file2, sep='\t', index_col = 0).T

    expr_matrix = pd.concat([ex_matrix_cond1, ex_matrix_cond2])

    if use_tf_list:
        tfnames = load_tf_names(project_dir + "/tfList_human.tsv")
    else:
        tfnames = 'all'

    network = grnboost2(expression_data=expr_matrix, gene_names = expr_matrix.columns, tf_names = tfnames, verbose = True)  

    network.to_csv(output_folder, sep='\t', index=False)

if __name__ == '__main__':
    arguments = docopt(__doc__, version='GRNBoost2 run')
    run_grnboost2(
        file1 = arguments["--inputfile1"],
        file2 = arguments["--inputfile2"],
        output_folder= arguments["--output"],
        use_tf_list = arguments['--use_tf_list'],
        project_dir = arguments['--project_directory']
    )
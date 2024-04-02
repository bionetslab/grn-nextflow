#!/usr/bin/env python
"""Usage: run_grnboost2.py -c FILE -d FILE -o FILE

Runs grnboost2 (https://github.com/tmoerman/arboreto/tree/master)

-h --help               show this
-c --inputfile1 FILE    specify input file
-d --inputfile2 FILE    specify input file
-o --output FILE        specify path and file name
--verbose               print more text

"""

# import python modules
import pandas as pd
from arboreto.algo import grnboost2
from docopt import docopt
import os

def run_grnboost2(file1, file2, output_folder):
    ex_matrix_cond1 = pd.read_csv(file1, sep='\t', index_col = 0).T
    ex_matrix_cond2 = pd.read_csv(file2, sep='\t', index_col = 0).T

    expr_matrix = pd.concat([ex_matrix_cond1, ex_matrix_cond2])

    network = grnboost2(expression_data=expr_matrix)        
    network.to_csv(output_folder, sep='\t', index=False)

if __name__ == '__main__':
    arguments = docopt(__doc__, version='GRNBoost2 run')
    run_grnboost2(
        file1 = arguments["--inputfile1"],
        file2 = arguments["--inputfile2"],
        output_folder= arguments["--output"]
    )
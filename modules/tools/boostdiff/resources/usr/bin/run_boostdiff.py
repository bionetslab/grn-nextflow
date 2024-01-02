#!/usr/bin/env python
"""Usage: run_boostdiff.py -c FILE -d FILE -o FOLDER [-n FEATURES] [-e NUM_ESTIMATORS] [-p PROCESSES] [-s SUBSAMPLES] [-u USE_TF_LIST] [-w PROJECT_DIRECTORY]

Wrapper script taking as input the case and control files.

-h --help               show this
-c --inputfile1 FILE    specify input file
-d --inputfile2 FILE    specify input file
-o --output FOLDER      specify output folder
-n --n_features FEATURES    specify number of features to use [default: 200]
-e --n_estimators NUM_ESTIMATORS    specify the number of estimators to use [default: 50]
-p --n_processes PROCESSES  number of threads [default: 1]
-s --n_subsamples SUBSAMPLES    number o. of samples in disease and control datasets used to fit a differential tree [default: 15]
-u --use_tf_list USE_TF_LIST Use transcription factor list?  True or False.  [default: True]
-w --project_directory PROJECT_DIRECTORY path to the nextflow directory
--verbose    print more text

"""

from docopt import docopt
from boostdiff.main_boostdiff import BoostDiff
import numpy as np
import pandas as pd
import sys


def run_boostdiff(file_control, file_case, output_folder, n_estimators, n_features,
                  n_subsamples, n_processes, use_tf_list, project_dir):


    keyword = "test"
    if (use_tf_list):
        tfList = pd.read_csv(project_dir + "/tfList_human.tsv", sep="\t", header = None).iloc[:, 0].values.tolist()
    else:
        tfList = "all"

    model = BoostDiff()
    model.run(file_case, file_control, output_folder, n_estimators, \
            n_features, n_subsamples=n_subsamples, keyword=keyword, n_processes=n_processes, regulators = tfList)

def cast_arguments(argument_dict):
    int_args = ['--n_estimators','--n_features', '--n_subsamples','--n_processes']
    for ia in int_args:
        try:
            argument_dict[ia] = int(argument_dict[ia])
        except ValueError:
            raise ValueError("Trying to parse argument to integer, but failed")
    return argument_dict

if __name__ == '__main__':
    arguments = docopt(__doc__, version='Boostdiff postprocessing')
    print(arguments)
    arguments = cast_arguments(argument_dict=arguments)
    print(arguments)
    run_boostdiff(file_case=arguments["--inputfile1"], 
                  file_control = arguments["--inputfile2"], 
                  output_folder = arguments["--output"],
                  n_estimators = arguments['--n_estimators'],
                  n_features = arguments['--n_features'],
                  n_subsamples= arguments['--n_subsamples'],
                  n_processes = arguments['--n_processes'],
                  use_tf_list = arguments['--use_tf_list'],
                  project_dir = arguments['--project_directory']
                  )



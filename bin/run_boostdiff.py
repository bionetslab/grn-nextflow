#!/usr/bin/env python
"""Usage: run_boostdiff.py -c FILE -d FILE -o FOLDER

Wrapper script taking as input the case and control files.

-h --help               show this
-c --inputfile1 FILE    specify input file
-d --inputfile2 FILE    specify input file
-o --output FOLDER      specify output folder
--verbose    print more text

"""

from docopt import docopt
from boostdiff.main_boostdiff import BoostDiff
import numpy as np



def run_boostdiff(file_control, file_case, output_folder):

    n_estimators = 50
    n_features = 50
    n_subsamples = 50
    keyword = "test"
    n_processes = 12

    model = BoostDiff()
    model.run(file_case, file_control, output_folder, n_estimators, \
            n_features, n_subsamples, keyword=keyword, n_processes=n_processes)


if __name__ == '__main__':
    arguments = docopt(__doc__, version='Boostdiff postprocessing')
    print(arguments)
    run_boostdiff(arguments["--inputfile1"], arguments["--inputfile2"], arguments["--output"])



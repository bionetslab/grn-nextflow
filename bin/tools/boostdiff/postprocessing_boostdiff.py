#!/usr/bin/env python
"""Usage: postprocessing_boostdiff.py -c FILE -d FILE -e FILE -o FOLDER

Wrapper script taking as input the case and control files.

-h --help               show this
-c --inputfile1 FILE    specify input file
-d --inputfile2 FILE    specify input file
-e --casecontrol FILE   specify case control log file
-o --output FOLDER      specify output folder
--verbose    print more text

"""

from docopt import docopt
import postprocessing as pp
import os.path as op
import os
import matplotlib.pyplot as plt
def postprocess_boostdiff(condition_1_path, condition_2_path, case_control_path, output_file):

    # find the true labels for case and control
    file_case_control = op.join(case_control_path)
    conditions = []
    with open(file_case_control, 'r') as handle:
        for line in handle:
            line = line.rstrip()
            line = line.split('\t')
            l = line[0].replace('out_', "")
            l = l.replace('.tsv', "")
            conditions.append(l)
    print(conditions)

    # Specify the output files containing the mean difference in prediction error after running the BoostDiff algorithm
    file_diff_dis = op.join(condition_1_path, "differences_test.txt")
    file_diff_con = op.join(condition_2_path, "differences_test.txt")

    # Specify the output file containing the raw network output after running the BoostDiff algorithm
    file_net_dis = op.join(condition_1_path, "boostdiff_network_test.txt")
    file_net_con = op.join(condition_2_path, "boostdiff_network_test.txt")

    # Filter the raw output based on the no. of top targets in the differences files
    # Then the top 50 edges for the run where the disease condition is the target condition
    # Also the top 50 edges for the run where the control condition is the target condition
    df_dis = pp.filter_network(file_net_dis, file_diff_dis, n_top_targets=10, n_top_edges=50)
    df_con = pp.filter_network(file_net_con, file_diff_con, n_top_targets=10, n_top_edges=50)

    # Example for real, large-scale datasets: filtering based on 3rd percentile with the p parameter 
    # df_filtered = pp.filter_network(file_net, file_diff, p=3, n_top_edges=100)

    # For plotting the differential network
    # Colorize by condition
    df_both = pp.colorize_by_condition(df_dis, df_con)
    # Generate and save the plot
    file_grn = op.join(output_file, "diff_grn.png")
    os.makedirs(output_file, exist_ok=True)
    pp.plot_grn(df_both, layout="graphviz_layout",show_conflicting=True, filename=file_grn, condition1=conditions[0], condition2=conditions[1])


    # Save the final differential network to file
    df_both.to_csv(op.join(output_file, 'diff_network.tsv'), sep="\t")

if __name__ == '__main__':
    arguments = docopt(__doc__, version='Boostdiff wrapper')
    print(arguments)
    # Turn interactive plotting off
    plt.ioff()
    postprocess_boostdiff(arguments["--inputfile1"], arguments["--inputfile2"], arguments['--casecontrol'], arguments["--output"])



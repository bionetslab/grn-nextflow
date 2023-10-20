#!/usr/bin/env python
"""Usage: aggregated_postprocessing_grnboost2.py -c FILE -d FILE -o FOLDER -r RUNS [-e NUM_EDGES]

Script that aggregates the result of --n_runs runs of grnboost2 and only keeps top $NUM_EDGES based on 
a union aggregation of the importance scores of the individual runs. Additionally, small linear
models are fitted on every edge to determine the regulatory (up or down regulation) of the interaction.

-h --help               show this
-o --output FOLDER      specify output folder
-c --inputfile1 FILE    specify bulk expression input file (same input as in run_grnboost2.py)
-d --inputfile2 FILE    specify bulk expression input file (same input as in run_grnboost2.py)
-r --n_runs RUNS    no. runs of tool
-e --top_n_edges NUM_EDGES    no. top edges to filter [default: 100]

--verbose    print more text

"""

from scipy.stats import linregress
from functools import reduce
from docopt import docopt

import pandas as pd
import numpy as np
import time
import os

def load_grnboost2_results(n_runs):
    """Loads the data of $n_runs grnboost2 runs.
    
    Args:
        n_runs (int): no. of total runs that were performed
    
    Retuns:
        networks (List): List that contains $n_runs inferred networks 
    """
    networks = []
    for i in range(n_runs):
        networks.append(pd.read_csv(os.path.join(f'run_{i+1}/network.tsv'), 
                                            sep='\t', names=["TF", "target", f"importance{i+1}"], header=None))
        
    return networks
    
def union_aggregation(networks):
    """Builds the aggregation of all networks by outer joining the networks and averaging over the score/importance

    Args:
        network (List): List that contains $n_runs inferred networks

    Returns:
        List: Aggregated network
    """
    result = reduce(lambda left, right: pd.merge(left, right, on=["TF", "target"], how="outer"), networks)
    result = result[1:] # merge creates first row of columns names -> dropping this row
    result = result.fillna(0)
    result["weight"] = result.iloc[:, 2:].astype(float).mean(axis=1)
    result = result.drop(result.iloc[:, 2:-1], axis=1)
    result = result.reset_index(drop=True)
    result = result.reindex(columns=["target", "TF", "weight"])
    result = result.rename(columns={"TF": "regulator"})
    
    return result

def grnboost2_filter_graphs(network, n_top_edges=200):
    """Filtering network as performed by boostdiff. Only keeping ${n_top_edges} edges to ${n_top_targets} target nodes.

    Args:
        network (Pandas df): inferred network
        n_top_edges (int, optional): Top n edges for filtering. Defaults to 200.

    Returns:
        Pandas df: filtered network
    """
    network = network.sort_values(by="weight", ascending=False)
    return network.head(n_top_edges)

def inh_or_act(file1, file2, network):
    """Fitting a small linear model onto every edge of the network to determine the regulatory effect (activator or repressor).

    Args:
        file1 (String): Path to file that contains bulk expression data for condition 1
        file2 (String): Path to file that contains bulk expression data for condition 2
        network (Pandas df): Aggregated and filtered network

    Returns:
        Pandas df: Aggregated and filtered network with additional regulatory edge information
    """
    
    
    
    ex_matrix_cond1 = pd.read_csv(file1, sep='\t')
    ex_matrix_cond2 = pd.read_csv(file2, sep='\t')

    expr_matrix = pd.concat([ex_matrix_cond1, ex_matrix_cond2.iloc[:, 1:]], axis=1)
    expr_matrix = expr_matrix.set_index('Gene')
    
    
    gene_names = list(expr_matrix.index)
    expr_matrix.index = expr_matrix.index.str.lower()

    slopes = []

    for _, row in network.iterrows():
        target = row[0].lower()
        regulator = row[1].lower()
        
        ge_regulator = np.array(list(expr_matrix.loc[regulator]))
        ge_target = np.array(list(expr_matrix.loc[target]))
        
        # outlier detection (|Z-score| > 2 -> outlier)
        target_mean = np.mean(ge_target)
        regulator_mean = np.mean(ge_regulator)
        target_sd = np.std(ge_target)
        regulator_sd = np.std(ge_regulator)
       
        target_mask = np.logical_or(ge_target > (target_mean + 2*target_sd), ge_target < (target_mean - 2*target_sd))
        regulator_mask = np.logical_or(ge_regulator > (regulator_mean + 2*regulator_sd), ge_regulator < (regulator_mean - 2*regulator_sd))
        combined_mask = np.logical_not(target_mask | regulator_mask)

        ge_regulator = ge_regulator[combined_mask]
        ge_target = ge_target[combined_mask]
        
        slope, intercept = np.polyfit(ge_regulator, ge_target, 1)
        slopes.append(slope)

    network["effect"] = slopes    
    return network

def write_output(network, output_path="", filename="aggregated_filtered_network.tsv"):
    """Writes the network to a tsv file

    Args:
        network (Pandas df): Aggregated and filtered network
        output_path (str, optional): Output path. Defaults to "".
        filename (str, optional): Output filename. Defaults to "aggregated_filtered_network.txt".
    """
    out = os.path.join(output_path, filename)
    os.makedirs(output_path, exist_ok=True)
    network.to_csv(out, index=False, sep="\t")

def main(file1, file2, output_folder, n_runs, top_n_edges=100):
    """Workflow for aggregating the runs and filtering/postprocessing the aggregated results

    Args:
        output_folder (str): Path to output folder
        n_runs (int): No. total runs
        top_n_targets (int, optional): Top n target nodes for filtering. Defaults to 20.
        top_n_edges (int, optional): Top n edges for filtering. Defaults to 100.
    """
    networks = load_grnboost2_results(n_runs=n_runs)
    
    # perform "union" aggregation
    aggregated_network = union_aggregation(networks)
    
    # filter the aggregated networks the same way as boostdiff does it
    aggregated_filtered_network = grnboost2_filter_graphs(aggregated_network, top_n_edges)
    
    # fit linear model on edges to check if edge is activator or inhibitor
    aggregated_filtered_network = inh_or_act(file1, file2, aggregated_filtered_network)
    aggregated_filtered_network["condition"] = "grnboost2"
    aggregated_filtered_network = aggregated_filtered_network.reindex(columns=["target", "regulator", "weight", "condition", "effect"])

    write_output(network=aggregated_filtered_network, output_path=output_folder, filename="aggregated_filtered_network_grnboost.txt")
    

def cast_arguments(argument_dict):
    """Casts integer arguments to int

    Args:
        argument_dict (dict): arguments dictionary

    Raises:
        ValueError: Invalid argument was passed

    Returns:
        dict: arguments dictionary with correctly casted arguments
    """
    int_args = ['--n_runs', '--top_n_edges']
    for ia in int_args:
        try:
            argument_dict[ia] = int(argument_dict[ia])
        except ValueError:
            raise ValueError("Trying to parse argument to integer, but failed")
    return argument_dict

if __name__ == '__main__':
    arguments = docopt(__doc__, version='Boostdiff postprocessing')
    arguments = cast_arguments(argument_dict=arguments)
    print(arguments)
    main(file1 = arguments["--inputfile1"],
         file2 = arguments["--inputfile2"],
        output_folder = arguments["--output"], 
        n_runs = arguments["--n_runs"],
        top_n_edges = arguments["--top_n_edges"]
    )



# N_RUNS = 9 # no. total runs of boostdiff (run 10 is somewhat strange)

# # condition 1 is equal to "Stronger in Docile"
# # condition 2 is equal to "Stronger in Armstrong"

# def load_data(n_runs=N_RUNS):
#     """
#     Loads the data of all boostDiffRuns
#     Params:

#     Returns:
#         Pandas df x4:
#             boostDiffNetworks and differencesTests for condition 1,2 for all $N_RUNS
    
#     """
#     boostDiffNetworks_condition1, boostDiffNetworks_condition2 = [], []
#     differencesTests_condition1, differencesTests_condition2 = [], []

#     for i in range(n_runs):
        
#         boostDiffNetworks_condition1.append(pd.read_csv(f'/data_slow/yb85avem/projects/NetMap/results_run{i+1}/Arm_vs_Doc_D28:Spleen/control/boostdiff_network_test.txt', 
#                                             sep='\t', names=["target", "regulator", f"weight{i+1}"], header=None))
#         boostDiffNetworks_condition2.append(pd.read_csv(f'/data_slow/yb85avem/projects/NetMap/results_run{i+1}/Arm_vs_Doc_D28:Spleen/disease/boostdiff_network_test.txt', 
#                                             sep='\t', names=["target", "regulator", f"weight{i+1}"], header=None))
#         differencesTests_condition1.append(pd.read_csv(f'/data_slow/yb85avem/projects/NetMap/results_run{i+1}/Arm_vs_Doc_D28:Spleen/control/differences_test.txt', 
#                                             sep='\t', names=["gene", f"error_diff{i+1}"], header=None))
#         differencesTests_condition2.append(pd.read_csv(f'/data_slow/yb85avem/projects/NetMap/results_run{i+1}/Arm_vs_Doc_D28:Spleen/disease/differences_test.txt', 
#                                             sep='\t', names=["gene", f"error_diff{i+1}"], header=None))

#     return boostDiffNetworks_condition1, boostDiffNetworks_condition2, differencesTests_condition1, differencesTests_condition2


# def intersection_filter(networks_condition1, networks_condition2, differencesTests_condition1, differencesTests_condition2, n_runs=N_RUNS):

#     res = []
#     for df in [networks_condition1, networks_condition2]:
#         df_tmp = reduce(lambda left, right: pd.merge(left, right, on=["target", "regulator"], how="inner"), df)
#         df_tmp["weight"] = df_tmp.iloc[:,2:].astype(float).mean(axis=1)
#         df_tmp = df_tmp.drop(df_tmp.iloc[:, 2:-1], axis=1)
#         res.append(df_tmp.reset_index(drop=True))

#     for df in [differencesTests_condition1, differencesTests_condition2]:
#         df_tmp = reduce(lambda left, right: pd.merge(left, right, on=["gene"], how="inner"), df)
#         df_tmp["error_diff"] = df_tmp.iloc[:,1:].astype(float).mean(axis=1)
#         df_tmp = df_tmp.drop(df_tmp.iloc[:, 1:-1], axis=1)
#         res.append(df_tmp.reset_index(drop=True))

#     return res



# def union_filter(networks_condition1, networks_condition2, differencesTests_condition1, differencesTests_condition2, n_runs=N_RUNS):
    
#     res = []
#     for df in [networks_condition1, networks_condition2]:
#         df_tmp = reduce(lambda left, right: pd.merge(left, right, on=["target", "regulator"], how="outer"), df)
#         df_tmp = df_tmp.fillna(0)
#         df_tmp["weight"] = df_tmp.iloc[:,2:].astype(float).mean(axis=1)
#         df_tmp = df_tmp.drop(df_tmp.iloc[:, 2:-1], axis=1)
#         res.append(df_tmp.reset_index(drop=True))


#     for df in [differencesTests_condition1, differencesTests_condition2]:
#         df_tmp = reduce(lambda left, right: pd.merge(left, right, on=["gene"], how="outer"), df)
#         df_tmp = df_tmp.fillna(0)
#         df_tmp["error_diff"] = df_tmp.iloc[:,1:].astype(float).mean(axis=1)
#         df_tmp = df_tmp.drop(df_tmp.iloc[:, 1:-1], axis=1)
#         res.append(df_tmp.reset_index(drop=True))

#     return res


# def boostDiff_filter_graphs(network_condition1, network_condition2, differencesTest_condition1, differencesTest_condition2, n_top_targets=20, n_top_edges=100):
#     differencesTest_condition1 = differencesTest_condition1.sort_values(by=["error_diff"], ascending=True)
#     differencesTest_condition2 = differencesTest_condition2.sort_values(by=["error_diff"], ascending=True)

#     top_target_genes_cond1 = list(differencesTest_condition1.head(n_top_targets)["gene"])
#     top_target_genes_cond2 = list(differencesTest_condition2.head(n_top_targets)["gene"])

#     # network_condition1 = network_condition1[network_condition1.target.isin(top_target_genes_cond1)]
#     # network_condition2 = network_condition2[network_condition2.target.isin(top_target_genes_cond2)]

#     network_condition1["condition"] = 1
#     network_condition2["condition"] = 2

#     full_network = pd.concat([network_condition1, network_condition2])
#     full_network = full_network.sort_values(by=["weight"], ascending=False)
#     # return full_network.head(n_top_edges)
#     return full_network


# def get_TOX_submodule(network_condition1, network_condition2, differencesTest_condition1, differencesTest_condition2, alpha=5 , beta=100, n_top_targets=20):

#     network_condition1["condition"] = 1
#     network_condition2["condition"] = 2

#     differencesTest_condition1 = differencesTest_condition1.sort_values(by=["error_diff"], ascending=True)
#     differencesTest_condition2 = differencesTest_condition2.sort_values(by=["error_diff"], ascending=True)


#     tox_incidentEdges = pd.concat([
#         network_condition1[(network_condition1["target"] == "TOX")
#                         | (network_condition1["regulator"] == "TOX")],
#         network_condition2[(network_condition2["target"] == "TOX")
#                         | (network_condition2["regulator"] == "TOX")]
#         ]).reset_index(drop=True)



#     top_target_genes_cond1 = list(differencesTest_condition1.head(n_top_targets)["gene"])
#     top_target_genes_cond2 = list(differencesTest_condition2.head(n_top_targets)["gene"])
#     tox_incidentEdges = tox_incidentEdges[(tox_incidentEdges.target.isin(top_target_genes_cond1)) | (tox_incidentEdges.target.isin(top_target_genes_cond2))]
    
    
#     # Idea for algorithm to find TOX submodules:
#     # 1) TOX is the seed node (distance == 0)
#     # 2) greedily take the edge with the highest weight
#     # 3) add edges of new node to search space if distance of new node is smaller than #alpha 
#     # 4) GOTO step 2) until #beta edges are added or search space is empty

#     nodes = pd.DataFrame({"nodes": ["TOX"], "distance": [0]})
#     edges = pd.DataFrame({"target": [], "regulator": [], "weight": [], "condition": []})
#     search_space = tox_incidentEdges.sort_values(by=["weight"], ascending=False)
#     while (len(edges.index) < beta) and (len(search_space.index != 0)):
#         # print(search_space)
#         new_edge = search_space.iloc[0,:]
#         edges.loc[len(edges)] = new_edge

#         target, regulator = new_edge["target"], new_edge["regulator"]

#         new_node_name = ""
#         distance_of_new_node = -1
        
#         if (target in nodes["nodes"].values) and (regulator in nodes["nodes"].values):
#             # both nodes are already contained but not the edge -> nothing to do
#             pass
#         elif target in nodes["nodes"].values: # new node is a regulator
#             distance_of_new_node = int(nodes[nodes["nodes"] == target]["distance"]) + 1
#             new_node_name = regulator
#             # new_node = {"nodes": [new_node_name], "distance": [distance_of_new_node]}
#             nodes.loc[len(nodes)] = [new_node_name, distance_of_new_node]

#         elif (regulator in nodes["nodes"].values): # new node is a target
#             distance_of_new_node = int(nodes[nodes["nodes"] == regulator]["distance"]) + 1
#             new_node_name = target
#             # new_node = {"nodes": new_node_name, "distance": distance_of_new_node}
#             nodes.loc[len(nodes)] = [new_node_name, distance_of_new_node]

#         else: # should not happen
#             raise ValueError("Error in Implementation")

#         # updating the search space:
#         if (distance_of_new_node != -1) and distance_of_new_node < alpha:
#             new_edges = pd.concat([
#                         network_condition1[(network_condition1["target"] == new_node_name)
#                         | (network_condition1["regulator"] == new_node_name)],
#                         network_condition2[(network_condition2["target"] == new_node_name)
#                         | (network_condition2["regulator"] == new_node_name)]
#             ]).reset_index(drop=True) 

#             # filtering by error_difff
#             new_edges = new_edges[(new_edges.target.isin(top_target_genes_cond1)) | (new_edges.target.isin(top_target_genes_cond2))]

#             search_space = pd.merge(search_space, new_edges, how="outer")
#             search_space = search_space.drop([0])
#             search_space = search_space.sort_values(by=["weight"], ascending=False).reset_index(drop=True)
#         else:
#             search_space = search_space.drop([0]).reset_index(drop=True)

#     return edges


# def write_output(network, output_path="", filename="network.txt"):
#     out = os.path.join(output_path, filename)
#     os.makedirs(output_path, exist_ok=True)
#     network.to_csv(out, index=False, sep="\t")


# def inh_or_act(file_cond1, file_cond2, edges):
#     df_cond1 = pd.read_csv(file_cond1, sep="\t")
#     df_cond2 = pd.read_csv(file_cond2, sep="\t")

#     slopes = []

#     for _, row in edges.iterrows():
#         target = row[0]
#         regulator = row[1]
#         condition = row[3]
        
#         ge_regulator, ge_target = None, None

#         if condition == 1.0: 
#             ge_regulator = df_cond1[df_cond1['Gene'] == regulator].iloc[:,1:].to_numpy().flatten()
#             ge_target = df_cond1[df_cond1['Gene'] == target].iloc[:,1:].to_numpy().flatten()
#         else:
#             ge_regulator = df_cond2[df_cond2['Gene'] == regulator].iloc[:,1:].to_numpy().flatten()
#             ge_target = df_cond2[df_cond2['Gene'] == target].iloc[:,1:].to_numpy().flatten()

#         slope, intercept = np.polyfit(ge_regulator, ge_target, 1)
#         slopes.append(slope)

#     edges["effect"] = slopes
#     return edges



# def main(file_condition1, file_condition2, output_folder, intersection=True, union=True):
#     """
#     """
#     # load the data
#     boostDiffNetworks_condition1, boostDiffNetworks_condition2, differencesTests_condition1, differencesTests_condition2 = load_data()

#     if intersection:
#         # perform filtering of N_RUNS
#         intersec_network_condition1, intersec_network_condition2, intersec_differencesTests_condition1, intersec_differencesTests_condition2 = intersection_filter(
#             boostDiffNetworks_condition1, boostDiffNetworks_condition2, differencesTests_condition1, differencesTests_condition2)
        
#         # perform postprocessing similar to boostDiff
#         complete_filtered_network = boostDiff_filter_graphs(
#             intersec_network_condition1, intersec_network_condition2, intersec_differencesTests_condition1, intersec_differencesTests_condition2)
#         complete_filtered_network = inh_or_act(file_condition1, file_condition2, complete_filtered_network)
#         write_output(network=complete_filtered_network, output_path=output_folder, filename="fullFilteredNet_intersec.txt")


#         # find tox submodules
#         tox_submodule = get_TOX_submodule(
#             intersec_network_condition1, intersec_network_condition2, intersec_differencesTests_condition1, intersec_differencesTests_condition2
#         )
#         tox_submodule = inh_or_act(file_condition1, file_condition2, tox_submodule)
#         write_output(network=tox_submodule, output_path=output_folder, filename="toxSubmodule_intersec.txt")
        

#     if union:
#         # perform filtering of N_RUNS
#         union_network_condition1, union_network_condition2, union_differencesTests_condition1, union_differencesTests_condition2 = union_filter(
#             boostDiffNetworks_condition1, boostDiffNetworks_condition2, differencesTests_condition1, differencesTests_condition2)

        
#         # perform postprocessing similar to boostDiff
#         complete_filtered_network = boostDiff_filter_graphs(
#             union_network_condition1, union_network_condition2, union_differencesTests_condition1, union_differencesTests_condition2)
#         complete_filtered_network = inh_or_act(file_condition1, file_condition2, complete_filtered_network)
#         write_output(network=complete_filtered_network, output_path=output_folder, filename="fullFilteredNet_union.txt")

        
#         # find tox submodules
#         tox_submodule = get_TOX_submodule(
#             union_network_condition1, union_network_condition2, union_differencesTests_condition1, union_differencesTests_condition2
#         )
#         tox_submodule = inh_or_act(file_condition1, file_condition2, tox_submodule)
        
#         write_output(network=tox_submodule, output_path=output_folder, filename="toxSubmodule_union.txt")

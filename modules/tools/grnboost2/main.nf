process RUN_TOOL {
  conda params.conda_env_path + '/grnboost2'
  maxForks 5
  publishDir params.publish_dir

  input:
  tuple val (key), path (case_file), path (control_file)
  each run

  output:
  tuple val (key), path ("${key}/grnboost2/run_${run}/network.tsv")

  script:
  """
  mkdir -p "${key}/grnboost2/run_${run}"
  mkdir -p "${key}/grnboost2/run_${run}"
  run_grnboost2.py -c $case_file -d $control_file -o ${key}/grnboost2/run_${run}/network.tsv
  """
}

process AGGREGATE_RESULTS {
  conda params.conda_env_path + '/grnboost2'
  publishDir params.publish_dir

  input:
  tuple val (key), path (gene_expr_cond1), path (gene_expr_cond2), path (network_files, stageAs: "run_?/*")
  val (runs)
  val (top_n_edges)

  output:
  tuple val (key), path ("${key}/grnboost2/aggregated_filtered_network_grnboost.txt")

  script:
  """
  aggregated_postprocessing_grnboost2.py -c $gene_expr_cond1 -d $gene_expr_cond2 -o ${key}/grnboost2/ -r $runs -e $top_n_edges
  """
}

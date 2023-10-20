process RUN_TOOL {
  maxForks 5
  publishDir params.publish_dir

  input:
  tuple val (key), path (case_file), path (control_file)
  tuple val (n_estimators), val (n_features), val (n_subsamples), val (n_processes) 
  each run
  
  output:
  tuple val (key), path ("${key}/boostdiff/run_${run}/disease"), path ("${key}/boostdiff/run_${run}/control"), path ("${key}/boostdiff/run_${run}/case_control.txt")

  script:
  """
  mkdir -p "${key}/boostdiff/run_${run}/disease"
  mkdir -p "${key}/boostdiff/run_${run}/control"
  # This is required to keep track of which file is which because boostdiff assigns disease control
  # and deletes the filename. GroupTuple may change the file order/
  touch "${key}/boostdiff/run_${run}/case_control.txt"
  echo "$case_file\t${key}/boostdiff/run_${run}/disease" >> "${key}/boostdiff/run_${run}/case_control.txt"
  echo "$control_file\t${key}/boostdiff/run_${run}/control" >> "${key}/boostdiff/run_${run}/case_control.txt"
  run_boostdiff.py -c $case_file -d $control_file -o ${key}/boostdiff/run_${run} -n $n_features -e $n_estimators -p $n_processes -s $n_subsamples
  """
}

process AGGREGATE_RESULTS {
  publishDir params.publish_dir

  input:
  tuple val (key), path (gene_expr_cond1), path (gene_expr_cond2), path (disease_files, stageAs: "run_?/*"), path (control_files, stageAs: "run_?/*"), path (case_control_files, stageAs: "/run_?/*")
  val (runs)
  tuple val (top_n_targets), val (top_n_edges)

  output:
  tuple val (key), path ("${key}/boostdiff/aggregated_filtered_network_boostdiff.txt")

  script:
  """
  aggregated_postprocessing_boostdiff.py -o ${key}/boostdiff/ -r $runs -t $top_n_targets -e $top_n_edges
  """
}

// process POSTPROCESS_BOOSTDIFF {
//   publishDir params.publish_dir

//   input:
//   val (key)
//   tuple path (case_file), path (control_file), path(case_control)
//   val(run)
  
//   output:
//   val(key)
//   path "${key}/run_${run}/diff_grn.png"
//   path "${key}/run_${run}/diff_network.tsv"

//   script:
//   """
//   postprocessing_boostdiff.py -c $case_file -d $control_file -e $case_control -o ${key}/run_${run}
//   """
// }

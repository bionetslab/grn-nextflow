process RUN_TOOL {
  conda params.conda_env_path + '/diffcoex'
  publishDir params.publish_dir 

  input:
  tuple val (key), path (case_file), path (control_file)
  val top_n_edges

  output:
  tuple val (key), path ("${key}/diffcoex/aggregated_filtered_network_diffcoex.txt")

  script:
  """
  mkdir -p "${key}/diffcoex/"
  run_diffcoex.R -c $case_file -d $control_file -o ${key}/diffcoex/ -n ${top_n_edges}
  """
}

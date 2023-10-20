process RUN_TOOL {
  publishDir params.publish_dir

  input:
  tuple val (key), path (case_file), path (control_file)
  val top_n_edges

  output:
  tuple val (key), path ("${key}/zscores/aggregated_filtered_network_zscore.txt")

  script:
  """
  mkdir -p "${key}/zscores/"
  run_zscores.R -c $case_file -d $control_file -o ${key}/zscores/ -n ${top_n_edges}
  """
}
maxcpus = Runtime.runtime.availableProcessors()

process RUN_TOOL {
  label 'chNet'
  publishDir params.publish_dir

  input:
  tuple val (key), path (case_file), path (control_file)
  val top_n_edges

  output:
  tuple val (key), path ("${key}/chNet/aggregated_filtered_network_chNet.txt")

  script:
  """
  mkdir -p "${key}/chNet/"
  run_zscores.R -c $case_file -d $control_file -o ${key}/chNet/ --n_cpus=$maxcpus
  """
}
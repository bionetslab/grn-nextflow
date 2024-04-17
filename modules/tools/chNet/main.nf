maxcpus = Runtime.runtime.availableProcessors()

process RUN_TOOL {
  label 'chnet'
  publishDir params.publish_dir

  input:
  tuple val (key), path (case_file), path (control_file)
  val lambda
  val parallel

  output:
  tuple val (key), path ("${key}/chNet/aggregated_filtered_network_chNet.txt")

  script:
  """
  mkdir -p "${key}/chNet/"
  run_zscores.R -c $case_file -d $control_file -o ${key}/chNet/ --n_cpus=$maxcpus --lambda=$lambda --run_parallel=$parallel
  """
}
process PLOT_GRN {
  conda params.conda_env_path + '/analysis'
  publishDir params.publish_dir

  input:
  tuple val (key), path (network)

  output:
  path ("${key}/network.html")

  script:
  """
  plotting_grn.R -i ${network} -o ${key} -n ${key}/network.html -p $projectDir
  """
}
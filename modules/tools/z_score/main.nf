process RUN_TOOL {
  publishDir params.publish_dir

  input:
  tuple val (key), path (case_file), path (control_file)
  each run

  output:
  tuple val (key), path ("${key}/zscores/disease"), path ("${key}/zscores/control"), path("${key}/zscores/case_control.txt")

  script:
  """
  mkdir -p "${key}/zscores/run_${run}/disease"
  mkdir -p "${key}/zscores/run_${run}/control"
  touch "${key}/zscores/run_${run}/case_control.txt"
  echo "$case_file\t${key}/zscores/run_${run}/disease" >> "${key}/zscores/run_${run}/case_control.txt"
  echo "$control_file\t${key}/zscores/run_${run}/control" >> "${key}/zscores/run_${run}/case_control.txt"
  run_zscores.py -c $case_file -d $control_file -o ${key}/zscores/run_${run}/
  """
}
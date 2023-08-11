process RUN_TOOL {
  publishDir params.publish_dir 

  input:
  tuple val (key), path (case_file), path (control_file)
  each run

  output:
  tuple val (key), path ("${key}/diffcoex/disease"), path ("${key}/diffcoex/control"), path("${key}/diffcoex/case_control.txt")

  script:
  """
  mkdir -p "${key}/diffcoex/run_${run}/disease"
  mkdir -p "${key}/diffcoex/run_${run}/control"
  touch "${key}/diffcoex/run_${run}/case_control.txt"
  echo "$case_file\t${key}/diffcoex/run_${run}/disease" >> "${key}/diffcoex/run_${run}/case_control.txt"
  echo "$control_file\t${key}/diffcoex/run_${run}/control" >> "${key}/diffcoex/run_${run}/case_control.txt"
  run_diffcoex.py -c $case_file -d $control_file -o ${key}/diffcoex/run_${run}/
  """
}

// process AGGREGATE_RESULTS {

// }
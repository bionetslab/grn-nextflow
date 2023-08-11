process RUN_TOOL {
  publishDir params.publish_dir 

  input:
  tuple val (key), path (case_file), path (control_file)
  each run

  output:
  tuple val (key), path ("${key}/ebcoexpress/disease"), path ("${key}/ebcoexpress/control"), path("${key}/ebcoexpress/case_control.txt")

  script:
  """
  mkdir -p "${key}/ebcoexpress/run_${run}/disease"
  mkdir -p "${key}/ebcoexpress/run_${run}/control"
  touch "${key}/ebcoexpress/run_${run}/case_control.txt"
  echo "$case_file\t${key}/ebcoexpress/run_${run}/disease" >> "${key}/ebcoexpress/run_${run}/case_control.txt"
  echo "$control_file\t${key}/ebcoexpress/run_${run}/control" >> "${key}/ebcoexpress/run_${run}/case_control.txt"
  run_ebcoexpress.py -c $case_file -d $control_file -o ${key}/ebcoexpress/run_${run}/
  """
}

// process AGGREGATE_RESULTS {

// }
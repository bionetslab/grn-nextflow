process CREATE_SHINY_APP {
  publishDir params.publish_dir

  input:
  val selection
  val diffgrn_tools
  val grn_tools

  output:
  path ("run_shiny.sh")

  script:
  """
  touch "run_shiny.sh"
  echo "Rscript ${projectDir}/modules/shiny_app/resources/usr/bin/shiny_app.R -r ${params.publish_dir} -p ${projectDir} -s ${selection} -f ${params.input} --dgrntools ${diffgrn_tools} --grntools ${grn_tools} " >> "run_shiny.sh"
  chmod u+x "run_shiny.sh"
  """
}

process CREATE_SHINY_APP_BCK {
  publishDir params.publish_dir

  input: 
  tuple val (key), path (case_file), path (control_file)
  tuple val (key), path (boostdiff_network)
  val selection

  output:
  path ("${key}/run_shiny.sh")
  path ("${key}/${case_file}")
  path ("${key}/${control_file}")

  script:
  """
  mkdir -p "${key}"
  touch "${key}/run_shiny.sh"
  cat ${case_file} >> "${key}/${case_file}"
  cat ${control_file} >> "${key}/${control_file}"
  echo "Rscript ${projectDir}/modules/shiny_app/resources/usr/bin/shiny_app.R -n ${boostdiff_network} -c ${case_file} -d ${control_file} -f ${params.data.seurat_object} -g ${params.data_loading_seurat.column_name} -s ${selection} -p $projectDir" >> "${key}/run_shiny.sh"
  chmod u+x ${key}/run_shiny.sh
  """
}
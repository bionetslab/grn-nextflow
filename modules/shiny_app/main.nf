process CREATE_SHINY_APP {
  publishDir params.publish_dir

  input:
  val seurat_file
  val selection
  val diffgrn_tools
  val grn_tools
  val mode
  val metacells

  output:
  path ("run_shiny.sh")

  script:
  """
  touch "run_shiny.sh"
  echo "Rscript ${projectDir}/modules/shiny_app/resources/usr/bin/shiny_app.R -r ${params.publish_dir} -p ${projectDir} -s ${selection} -f ${seurat_file} --dgrntools ${diffgrn_tools} --grntools ${grn_tools} -m ${mode}" --metacells ${metacells} >> "run_shiny.sh"
  chmod u+x "run_shiny.sh"
  """
}

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
  path ("app/run_shiny.sh")
  path ("app/style_generator.R")
  path ("app/auxiliary_functions.R")
  path ("app/shiny_app.R")
  path ("app/style.js")
  path ("data/data.rds")
  

  script:
  """
  mkdir app
  mkdir data
  touch "app/run_shiny.sh"
  echo "Rscript app/shiny_app.R -r '/' -s ${selection} -f '/data/data.rds' --dgrntools ${diffgrn_tools} --grntools ${grn_tools} -m ${mode}" --metacells ${metacells} >> "app/run_shiny.sh"
  chmod u+x "app/run_shiny.sh"
  cp ${projectDir}/modules/shiny_app/resources/usr/bin/* app
  cp ${projectDir}/environment_configs/shiny.yml app
  cp ${projectDir}/install_shiny_packages.R app
  cp ${seurat_file} data/data.rds
  """
}

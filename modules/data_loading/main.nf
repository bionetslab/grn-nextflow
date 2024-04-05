process CREATE_METACELLS_SEURAT {
  label 'big_mem'
  label 'create_metacells'

  input:
  path seurat_object
  tuple val (key), val (assay), val (selection_criteria_keys), val(name), val (selection_criteria)
  val mode
  val create_metacells

  output:
  tuple  val (key), path("${name}.tsv")

  script:
  """
  create_metacells.R -f $seurat_object -o "$name".tsv -g $selection_criteria_keys -s $selection_criteria -a $assay -m $mode --key=$key --create_metacells=$create_metacells
  """
}

// process CREATE_METACELLS_TSVFILES {
//   label 'big_mem'
//   conda params.conda_env_path + '/create_metacells'

//   input:
//   tuple val (key), path (tsv_file_cond1), path (tsv_file_cond2)

//   output:
//   tuple val (key), path("${tsv_file_cond1}")

//   script:
//   """
//   create_metacells.R -f $tsv_file -o "aggregated_${tsv_file}" -m $mode 
//   """  
// }

process CONVERT_ANNDATA_TO_SEURAT {
  publishDir params.publish_dir
  label 'big_mem'
  label 'anndata_to_seurat'

  input:
  path anndata_object

  output:
  path ("seurat_object.rds")

  script:
  """
  convert_anndata_to_seurat.R -i $anndata_object
  """
}

process CHECK_FILES {
  label 'create_metacells'
  publishDir params.publish_dir

  input:
  tuple val (key), path (files)
  
  output:
  tuple val (key), path ("${key}/out_${files[0]}"), path ("${key}/out_${files[1]}")
  
  script:
  """
  mkdir -p "${key}/"
  check_input_files.R -c ${files[0]} -d ${files[1]} -e ${key}/out_${files[0]} -f ${key}/out_${files[1]} -k ${key} -r $params.publish_dir
  """
}

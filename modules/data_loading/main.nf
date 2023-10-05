process CREATE_METACELLS_SEURAT {
  label 'big_mem'

  input:
  path seurat_object
  tuple val (key), val (assay), val (selection_criteria_keys), val(name), val (selection_criteria), val (cluster_name), val (cluster_ids)
  val mode

  output:
  tuple  val (key), path("${name}.tsv")

  script:
  """
  create_metacells.R -f $seurat_object -o "$name".tsv -g $selection_criteria_keys -s $selection_criteria -n 30 -l $cluster_name -k $cluster_ids -a $assay -m $mode
  """
}

process CREATE_METACELLS_TSVFILES {
  label 'big_mem'

  input:
  tuple val (key), path (tsv_file)
  val mode

  output:
  tuple val (key), path("aggregated_${tsv_file}")

  script:
  """
  create_metacells.R -f $tsv_file -o "aggregated_${tsv_file}" -n 30 -m $mode 
  """
  
}

process CONVERT_SCANPY_TO_SEURAT {
  label 'big_mem'

  input:
  path scanpy_object

  output:
  path ("seurat_object.rds")

  script:
  """
  convert_scanpy_to_seurat.R -f $scanpy_object
  """
}

process CHECK_FILES {
  publishDir params.publish_dir

  input:
  tuple val (key), path (files)
  
  output:
  tuple val (key), path ("${key}/out_${files[0]}"), path ("${key}/out_${files[1]}")
  // path ("${key}/out_${files[0]}")
  // path ("${key}/out_${files[1]}")

  script:
  """
  mkdir -p "${key}/"
  check_input_files.R -c ${files[0]} -d ${files[1]} -e ${key}/out_${files[0]} -f ${key}/out_${files[1]} -k ${key} -r $params.publish_dir
  """
}
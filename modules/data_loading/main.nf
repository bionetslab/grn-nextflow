process SELECT_DATA {
  label 'big_mem'

  input:
  path seurat_object
  tuple val (key), val (assay), val (selection_criteria_keys), val (selection_criteria), val(name), val (cluster_name), val (cluster_ids)


  output:
  tuple  val (key), path("${name}.tsv")

  script:
  """
  create_metacells.R -f $seurat_object -o "$name".tsv -g $selection_criteria_keys -s $selection_criteria -n 30 -l $cluster_name -k $cluster_ids -a $assay
  """
}

process CHECK_FILES {
  input:
  tuple val (key), path (files)
  
  output:
  tuple val (key), path ("out_${files[0]}"), path ("out_${files[1]}")

  script:
  """
  check_input_files.R -c ${files[0]} -d ${files[1]} -e out_${files[0]} -f out_${files[1]}
  """
}
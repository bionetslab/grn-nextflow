//params.seurat_object = "$workflow.homeDir/../shared/netmap/data/ga_an0228_10x_deepseq_filtered_smarta_merged_tissue_integrated_rep_timepoint_infection_filtered_seurat.rds"
params.seurat_object = "$projectDir/../data/ga_an0228_10x_deepseq_filtered_smarta_merged_tissue_integrated_rep_timepoint_infection_filtered_seurat.rds"

params.column_name = 'infection:tissue:subject:time'
params.cluster_name='cluster'
params.publishDir = "$projectDir/../results/"
params.condaDir = "$projectDir/boostdiff-wf/" 
params.transcription_factors = "$projectDir/../data/DatabaseExtract_v_1.01.csv"

/*
 * define the INDEX process that creates a binary index
 * given the transcriptome file
 */

process SELECT_DATA {

  label 'big_mem'

  input:
  path seurat_object
  val column_name
  tuple val (selection_criteria), val(name), val(key), val (cluster_name), val (cluster_ids)


  output:
  tuple  val (key), path("${name}.tsv")

  script:
  """
  echo $seurat_object
  echo $workflow.homeDir
  create_metacells.R -f $seurat_object -o "$name".tsv -g $column_name -s $selection_criteria -n 30 -l $cluster_name -k $cluster_ids
  """
}

process CHECK_FILES {

  input:
  tuple val (key), path (files)
  path tfdb
  
  output:
  tuple val(key), path ("out_${files[0]}"), path ("out_${files[1]}"), emit: files
  path "transcription_factors.tsv", emit: transcriptome

  script:
  """
  check_input_files.R -c ${files[0]} -d ${files[1]} -e out_${files[0]} -f out_${files[1]} -t $tfdb -u "transcription_factors.tsv"
  """
}

process RUN_BOOSTDIFF {

  conda 'boostdiff-wf.yml'

  publishDir params.publishDir

  input:
  tuple val(key), path (case_file), path (control_file)
  val n_estimators
  val n_features
  val n_subsamples
  val n_processes
  path regulators
  val regulator_column
  
  output:
  tuple val(key), path("${key}/disease"), path("${key}/control"), path("${key}/case_control.txt")
  

  script:
  """
  mkdir -p "${key}/disease"
  mkdir -p "${key}/control"
  # This is required to keep track of which file is which because boostdiff assigns disease control
  # and deletes the filename. GroupTuple may change the file order/
  touch "${key}/case_control.txt"
  echo "$case_file\t${key}/disease" >> "${key}/case_control.txt"
  echo "$control_file\t${key}/control" >> "${key}/case_control.txt"
  run_boostdiff.py -c $case_file -d $control_file -o $key -n $n_features -e $n_estimators -p $n_processes -r $regulators -q $regulator_column -s $n_subsamples
  """
}

process POSTPROCESS_BOOSTDIFF {

  publishDir params.publishDir

  input:
  tuple val(key), path (case_file), path (control_file), path(case_control)
  
  output:
  val(key)
  path "${key}/diff_grn.png"
  path "${key}/diff_network.tsv"

  script:
  """
  postprocessing_boostdiff.py -c $case_file -d $control_file -e $case_control -o $key
  """
}

process GENE_SET_ENRICHMENT {

  publishDir params.publishDir

  input:
  val(key)
  path(network_file)

  output:
  path "${key}/gse_analysis.tsv"

  script:
  """
  gene_set_enrichment.py -i $network_file -o $key
  """
}


workflow {
  // How to configure this pipeline:
 
  // GLOBAL VARIABLES
  // params.seurat_object = "$projectDir/../data/ga_an0228_10x_deepseq_filtered_smarta_merged_tissue_integrated_rep_timepoint_infection_filtered_seurat.rds"
  // params.column_name = 'infection:tissue:subject:time' : This describes the selection criteria column_names for the data points
  // params.cluster_name='cluster'
  // params.publishDir = "$projectDir/../results/"

  /* This is the selection variable:
     [Selected_covariate_configuation, 'Output_file_name', 'Key_matching_two_files_for_boostdiff', 'Output_folder', 'secondary_selection_criterion', 'categories']
      1. Selected_covariate_configuation: This refers to the params.colum_name criterion. There must be the same number of entries as the params.colum_name string (lists not allowed)
      2. Output file name for this data subset. Must be unique within the output folder
      3. Boostdiff requires two input files. To match the inputs in the pipeline the keys need to be matching and unique accross the pipeline run.
      4. Output folder for the run
      5. Secondary selection criterion allowing for multi selection. Only one column can be selected but from this column multiple categories (see 6)
      6. Categories to select from the column defined in 5.
  */

  println(params.condaDir)
  input_case_ch = Channel
    .fromList( [//['Doc:Spleen:1:d28,Doc:Spleen:2:d28,Doc:Spleen:4:d28', "Doc_Spleen_d28", 'Doc:Spleen', 'cluster', '1:2'],
                //['Doc:Spleen:1:d10,Doc:Spleen:3:d10,Doc:Spleen:5:d10', "Doc_Spleen_d10", 'Doc:Spleen', 'cluster', '1:2'],
                // ['Doc:Liver:1:d28,Doc:Liver:2:d28,Doc:Liver:4:d28', "Doc_Liver_d28", 'Doc:Liver', 'cluster', '1:2'],
                // ['Doc:Liver:1:d10,Doc:Liver:3:d10,Doc:Liver:5:d10', "Doc_Liver_d10", 'Doc:Liver', 'cluster', '1:2'],
                // ['Arm:Liver:1:d28,Arm:Liver:3:d28,Arm:Liver:5:d28', "Arm_Liver_d28", 'Arm:Liver', 'cluster', '1:2'],
                // ['Arm:Liver:2:d10,Arm:Liver:3:d10,Arm:Liver:4:d10', "Arm_Liver_d10", 'Arm:Liver', 'cluster', '1:2'],
                ['Arm:Spleen:1:d28,Arm:Spleen:3:d28,Arm:Spleen:5:d28', "Arm_Spleen_d28", 'Arm:Spleen', 'cluster', '1:2'],
                ['Arm:Spleen:2:d10,Arm:Spleen:3:d10,Arm:Spleen:4:d10', "Arm_Spleen_d10", 'Arm:Spleen', 'cluster', '1:2']
                ] )

  input_case_ch.view { "$it" }


  data_case_ch = SELECT_DATA(params.seurat_object, params.column_name , input_case_ch)
  data_case_ch.view()


  tuple_ch = data_case_ch.groupTuple()
  tuple_ch.view { "value: $it" }

  checked_ch = CHECK_FILES(tuple_ch, params.transcription_factors)
  checked_ch.files.view()
  checked_ch.transcriptome.view()
  boostdiff_ch = RUN_BOOSTDIFF(checked_ch, 50, 30, 8, checked_ch.transcriptome, 'regulators')
  boostdiff_ch.view()
  // processed_ch = POSTPROCESS_BOOSTDIFF(boostdiff_ch)
  // gse_ch = GENE_SET_ENRICHMENT(processed_ch[0], processed_ch[2])

}

  
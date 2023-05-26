params.seurat_object = "$workflow.homeDir/../shared/netmap/data/ga_an0228_10x_deepseq_filtered_smarta_merged_tissue_integrated_rep_timepoint_infection_filtered_seurat.rds"
params.column_name = 'infection:tissue:subject:time'
params.cluster_name='cluster'
params.publish_dir = "$projectDir/../results/"

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
  create_metacells.R -f $seurat_object -o "$name".tsv -g $column_name -s $selection_criteria -n 30 -l $cluster_name -k $cluster_ids
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

process RUN_BOOSTDIFF {
  publishDir params.publish_dir

  input:
  tuple val (key), path (case_file), path (control_file)
  tuple val (n_estimators), val (n_features), val (n_subsamples), val (n_processes) 
  each run
  
  output:
  tuple val (key), path ("${key}/run_${run}/disease"), path ("${key}/run_${run}/control"), path ("${key}/run_${run}/case_control.txt")

  script:
  """
  mkdir -p "${key}/run_${run}/disease"
  mkdir -p "${key}/run_${run}/control"
  # This is required to keep track of which file is which because boostdiff assigns disease control
  # and deletes the filename. GroupTuple may change the file order/
  touch "${key}/run_${run}/case_control.txt"
  echo "$case_file\t${key}/run_${run}/disease" >> "${key}/run_${run}/case_control.txt"
  echo "$control_file\t${key}/run_${run}/control" >> "${key}/run_${run}/case_control.txt"
  run_boostdiff.py -c $case_file -d $control_file -o ${key}/run_${run} -n $n_features -e $n_estimators -p $n_processes -s $n_subsamples
  """
}

process POSTPROCESS_BOOSTDIFF {
  publishDir params.publish_dir

  input:
  val (key)
  tuple path (case_file), path (control_file), path(case_control)
  val(run)
  
  output:
  val(key)
  path "${key}/run_${run}/diff_grn.png"
  path "${key}/run_${run}/diff_network.tsv"

  script:
  """
  postprocessing_boostdiff.py -c $case_file -d $control_file -e $case_control -o ${key}/run_${run}
  """
}

process AGGREGATE_POSTPROCESS_BOOSTDIFF {
  publishDir params.publish_dir

  input:
  tuple val (key), path (gene_expr_cond1), path (gene_expr_cond2), path (disease_files, stageAs: "run_?/*"), path (control_files, stageAs: "run_?/*"), path (case_control_files, stageAs: "/run_?/*")
  val (runs)
  tuple val (top_n_targets), val (top_n_edges)

  output:
  tuple val (key), path ("${key}/aggregated_filtered_network.txt"), val (top_n_targets), val (top_n_edges)

  script:
  """
  aggregated_postprocessing.py -o $key -r $runs -t $top_n_targets -e $top_n_edges
  """
}

process PLOT_GRN {
  publishDir params.publish_dir

  input:
  tuple val (key), path (network), val (top_n_targets), val (top_n_edges)

  output:
  path ("${key}/diffGRN_top${top_n_targets}targets_top${top_n_edges}edges.html")

  script:
  """
  plotting_grn.R -i ${network} -o ${key} -n "${key}/diffGRN_top${top_n_targets}targets_top${top_n_edges}edges.html"
  """
}

workflow {
  // How to configure this pipeline:
 
  /* GLOBAL VARIABLES 
  params.seurat_object = "$projectDir/../data/ga_an0228_10x_deepseq_filtered_smarta_merged_tissue_integrated_rep_timepoint_infection_filtered_seurat.rds"
    => Path to the data
  params.column_name = 'infection:tissue:subject:time' 
    => This describes the selection criteria column_names for the data points
  params.cluster_name = 'cluster'
    => This describes the cluster name
  params.publish_dir = "$projectDir/../results/"
    => Sets the path to the output results
  */ 

  /* This is the selection variable:
     [Selected_covariate_configuation, 'Output_file_name', 'Key_matching_two_files_for_boostdiff', 'Output_folder', 'secondary_selection_criterion', 'categories']
      1. Selected_covariate_configuation: This refers to the params.colum_name criterion. There must be the same number of entries as the params.colum_name string (lists not allowed)
      2. Output file name for this data subset. Must be unique within the output folder
      3. Boostdiff requires two input files. To match the inputs in the pipeline the keys (third column) need to be matching and unique accross the pipeline run.
      4. Output folder for the run
      5. Secondary selection criterion allowing for multi selection. Only one column can be selected but from this column multiple categories (see 6)
      6. Categories to select from the column defined in 5.
  */
  // selection = 
  //   [
  //     ['Doc:Spleen:1:d28,Doc:Spleen:2:d28,Doc:Spleen:4:d28', "Doc_Spleen_d28", 'Doc:Spleen', 'cluster', '1:2'],
  //     ['Doc:Spleen:1:d10,Doc:Spleen:3:d10,Doc:Spleen:5:d10', "Doc_Spleen_d10", 'Doc:Spleen', 'cluster', '1:2'],
  //     ['Doc:Liver:1:d28,Doc:Liver:2:d28,Doc:Liver:4:d28', "Doc_Liver_d28", 'Doc:Liver', 'cluster', '1:2'],
  //     ['Doc:Liver:1:d10,Doc:Liver:3:d10,Doc:Liver:5:d10', "Doc_Liver_d10", 'Doc:Liver', 'cluster', '1:2'],
  //     ['Arm:Liver:1:d28,Arm:Liver:3:d28,Arm:Liver:5:d28', "Arm_Liver_d28", 'Arm:Liver', 'cluster', '1:2'],
  //     ['Arm:Liver:2:d10,Arm:Liver:3:d10,Arm:Liver:4:d10', "Arm_Liver_d10", 'Arm:Liver', 'cluster', '1:2'],
  //     ['Arm:Spleen:1:d28,Arm:Spleen:3:d28,Arm:Spleen:5:d28', "Arm_Spleen_d28", 'Arm:Spleen', 'cluster', '1:2'],
  //     ['Arm:Spleen:2:d10,Arm:Spleen:3:d10,Arm:Spleen:4:d10', "Arm_Spleen_d10", 'Arm:Spleen', 'cluster', '1:2']
  //   ]
  
  selection = 
    [
      ['Doc:Spleen:1:d28,Doc:Spleen:2:d28,Doc:Spleen:4:d28', "Doc_Spleen_d28", 'Arm_vs_Doc_D28:Spleen', 'cluster', '1:2'],
      ['Arm:Spleen:1:d28,Arm:Spleen:3:d28,Arm:Spleen:5:d28', "Arm_Spleen_d28", 'Arm_vs_Doc_D28:Spleen', 'cluster', '1:2'],
      // ['Doc:Spleen:1:d10,Doc:Spleen:3:d10,Doc:Spleen:5:d10', "Doc_Spleen_d10", 'Arm_vs_Doc_D10:Spleen', 'cluster', '1:2'],
      // ['Arm:Spleen:2:d10,Arm:Spleen:3:d10,Arm:Spleen:4:d10', "Arm_Spleen_d10", 'Arm_vs_Doc_D10:Spleen', 'cluster', '1:2'],
      // ['Doc:Liver:1:d28,Doc:Liver:2:d28,Doc:Liver:4:d28', "Doc_Liver_d28", 'Arm_vs_Doc_D28:Liver', 'cluster', '1:2'],
      // ['Arm:Liver:1:d28,Arm:Liver:3:d28,Arm:Liver:5:d28', "Arm_Liver_d28", 'Arm_vs_Doc_D28:Liver', 'cluster', '1:2'],
      // ['Doc:Liver:1:d10,Doc:Liver:3:d10,Doc:Liver:5:d10', "Doc_Liver_d10", 'Arm_vs_Doc_D10:Liver', 'cluster', '1:2'],
      // ['Arm:Liver:2:d10,Arm:Liver:3:d10,Arm:Liver:4:d10', "Arm_Liver_d10", 'Arm_vs_Doc_D10:Liver', 'cluster', '1:2'],
    ]
  
  n_runs = 10 // number of total runs of boostdiff
  runs = (1..n_runs)
  boostdiff_run_settings = [50, 1500, 30, 8] // [n_estimators, n_features, n_subsamples, n_processes]
  boostdiff_filtering_settings = [20, 100]  // [top_n_targets, top_n_edges]


  // Workflow:
  input_case_ch = Channel.fromList(selection)  
  input_case_ch.view { "$it" }
  
  data_case_ch = SELECT_DATA(params.seurat_object, params.column_name , input_case_ch)
  data_case_ch.view()

  tuple_ch = data_case_ch.groupTuple()
  tuple_ch.view { "value: $it" }

  checked_ch = CHECK_FILES(tuple_ch)
  checked_ch.view()

  boostdiff_ch = RUN_BOOSTDIFF(checked_ch, boostdiff_run_settings, runs)
  boostdiff_grouped_ch = boostdiff_ch.groupTuple()
  boostdiff_grouped_ch.view { "value: $it"}

  aggregate_ch = AGGREGATE_POSTPROCESS_BOOSTDIFF(checked_ch.join(boostdiff_grouped_ch), n_runs, boostdiff_filtering_settings)
  aggregate_ch.view()

  plot_ch = PLOT_GRN(aggregate_ch)
  // plot_ch.view()
}
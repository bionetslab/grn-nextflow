if ( params.help ) {
  help = """boostdiff.nf: Performs --n_runs of boostdiff on filtered data of --seurat_object based on the selection criteria.
                          Check in boostdiff.nf "workflow { ... }" for more information about the configuration of the selection criteria since it is currently hardcoded there.

                          | Required arguments:
                          |   --seurat_object Location of the seurat object data file 
                          |                   [default: ${params.seurat_object}]  
                          |
                          | Optional arguments:
                          |   --publish_dir   Location to store the output files
                          |                   [default  ${params.publish_dir}]
                          |   --column_name   Names of columns in seurat object data file
                          |                   [default: ${params.column_name}]
                          |   --n_runs        No. runs of boostdiff to perform (to combat randomness of the tool)
                          |                   [default: ${params.n_runs}]
                          |   --n_estimators  No. estimators to use for each boostdiff run
                          |                   [default: ${params.n_estimators}]
                          |   --n_features    No. features to use for each boostdiff run
                          |                   [default: ${params.n_features}]
                          |   --n_subsamples  No. subsamples to use for each boostdiff run
                          |                   [default: ${params.n_subsamples}]
                          |   --n_processes   No. processes to use for each boostdiff run
                          |                   [default: ${params.n_processes}]
                          |   --top_n_targets No. top target genes to filter results
                          |                   [default: ${params.top_n_targets}]
                          |   --top_n_edges   No. top edges to filter results
                          |                   [default: ${params.top_n_edges}]
  """.stripMargin()
  println(help)
  exit(0)
}

log.info"""\
BOOSTDIFF - NF    
=======================

seurat_object : $params.seurat_object
publish_dir   : $params.publish_dir
column_name   : $params.column_name
n_runs        : $params.n_runs
n_estimators  : $params.n_estimators
n_features    : $params.n_features   
n_subsamples  : $params.n_subsamples
n_processes   : $params.n_processes
top_n_targets : $params.top_n_targets
top_n_edges   : $params.top_n_edges

Run with --help for additional information
"""

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
      ['Doc:Liver:1:d10,Doc:Liver:3:d10,Doc:Liver:5:d10', "Doc_Liver_d10", 'Arm_vs_Doc_D10:Liver', 'cluster', '1:2'],
      ['Arm:Liver:2:d10,Arm:Liver:3:d10,Arm:Liver:4:d10', "Arm_Liver_d10", 'Arm_vs_Doc_D10:Liver', 'cluster', '1:2'],
    ]
  
  runs = (1..params.n_runs)
  boostdiff_run_settings = // [n_estimators, n_features, n_subsamples, n_processes]
  [
    params.n_estimators,
    params.n_features,
    params.n_subsamples,
    params.n_processes
  ] 
  boostdiff_filtering_settings = // [top_n_targets, top_n_edges]
  [
    params.top_n_targets,
    params.top_n_edges
  ]  


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
  boostdiff_grouped_ch.view { "value: $it" }

  aggregate_ch = AGGREGATE_POSTPROCESS_BOOSTDIFF(checked_ch.join(boostdiff_grouped_ch), params.n_runs, boostdiff_filtering_settings)
  aggregate_ch.view()

  plot_ch = PLOT_GRN(aggregate_ch)
  // plot_ch.view()
}
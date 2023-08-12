// if ( params.help ) {
//   help = """
//   boostdiff.nf: Performs --n_runs of boostdiff on filtered data of --seurat_object based on the selection criteria.
//                 Check in boostdiff.nf "workflow { ... }" for more information about the configuration of the selection criteria since it is currently hardcoded there.

//                           | Arguments:
//                           |   --seurat_object Location of the seurat object data file 
//                           |                   [default: ${params.data.seurat_object}]
//                           |   --cond1_file    Location of a .tsv file to use for the pipeline (if this and cond2_file is set the CREATE_METACELLS step is skipped)
//                           |                   [default: ${params.data.cond1_file}]
//                           |   --cond2_file    Location of a .tsv file to use for the pipeline (if this and cond1_file is set the CREATE_METACELLS step is skipped)  
//                           |                   [default: ${params.data.cond2_file}]
//                           |   --publish_dir   Location to store the output files
//                           |                   [default  ${params.publish_dir}]
//                           |   --tools         Name of tools to use as comma separated list (Available: boostdiff, grnboost2, z_score, diffcoex, ebcoexpress)
//                           |                   [default: ${params.pipeline.tools}]
//                           |   --column_name   Names of columns in seurat object data file
//                           |                   [default: ${params.data_loading_seurat.column_name}]
//                           |   --n_runs        No. runs of boostdiff to perform (to combat randomness of the tool)
//                           |                   [default: ${params.pipeline.n_runs}]
//                           |   --n_estimators  No. estimators to use for each boostdiff run
//                           |                   [default: ${params.boostdiff.n_estimators}]
//                           |   --n_features    No. features to use for each boostdiff run
//                           |                   [default: ${params.boostdiff.n_features}]
//                           |   --n_subsamples  No. subsamples to use for each boostdiff run
//                           |                   [default: ${params.boostdiff.n_subsamples}]
//                           |   --n_processes   No. processes to use for each boostdiff run
//                           |                   [default: ${params.boostdiff.n_processes}]
//                           |   --top_n_targets No. top target genes to filter results
//                           |                   [default: ${params.filtering.top_n_targets}]
//                           |   --top_n_edges   No. top edges to filter results
//                           |                   [default: ${params.filtering.top_n_edges}]
//   """.stripMargin()
//   println(help)
//   exit(0)
// }

// log.info"""\
// BOOSTDIFF - NF    
// =======================

// publish_dir   : $params.publish_dir
// tools         : $params.pipeline.tools
// n_runs        : $params.pipeline.n_runs
// seurat_object : $params.data.seurat_object
// cond1_file    : $params.data.cond1_file
// cond2_file    : $params.data.cond2_file
// column_name   : $params.data_loading_seurat.column_name
// n_estimators  : $params.boostdiff.n_estimators
// n_features    : $params.boostdiff.n_features   
// n_subsamples  : $params.boostdiff.n_subsamples
// n_processes   : $params.boostdiff.n_processes
// top_n_targets : $params.filtering.top_n_targets
// top_n_edges   : $params.filtering.top_n_edges

// Run with --help for additional information
// """

include { LOAD_DATA } from './subworkflows/load_data'
include { RUN_TOOLS } from './subworkflows/run_tools'
include { ANALYSIS } from './subworkflows/analysis'
include { SHINY_APP } from './subworkflows/shiny_app'


workflow {
  
  List<String> available_GRNInference_tools = params.available_GRNInference_tools.split(",")
  List<String> available_DGRNInference_tools = params.available_DGRNInference_tools.split(",")  
  List<String> tools = params.tools.split(",")
  tools.each { element -> 
    if (!available_GRNInference_tools.contains(element) && !available_DGRNInference_tools.contains(element)) {
          throw new Exception(element + " is not available. Please choose only a selection of the following tools: " + (available_GRNInference_tools + available_DGRNInference_tools))
    }
  }

  data = LOAD_DATA(tools)
  networks = RUN_TOOLS(data, tools)
  ANALYSIS(networks)
  // // TODO: Make NF pipeline creation of shiny App work with variable no. input network files and different data input formats
  // SHINY_APP()
}
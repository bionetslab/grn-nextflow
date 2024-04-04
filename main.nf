if ( params.help ) {
  help = """Check the github https://github.com/bionetslab/grn-nextflow/tree/dgrn_nf for more information.""".stripMargin()
  println(help)
  exit(0)
}

log.info"""\
BOOSTDIFF - NF    
=======================

publish_dir             : $params.publish_dir
mode                    : $params.mode
tools                   : $params.tools
n_runs                  : $params.n_runs
create_metacells        : $params.create_metacells
boostdiff.n_estimators  : $params.boostdiff.n_estimators
boostdiff.n_features    : $params.boostdiff.n_features   
boostdiff.n_subsamples  : $params.boostdiff.n_subsamples
boostdiff.n_processes   : $params.boostdiff.n_processes
boostdiff.top_n_targets : $params.boostdiff.top_n_targets
boostdiff.top_n_edges   : $params.boostdiff.top_n_edges

Run with --help for additional information
"""

include { LOAD_DATA } from './subworkflows/load_data'
include { RUN_TOOLS } from './subworkflows/run_tools'
include { ANALYSIS } from './subworkflows/analysis'
include { SHINY_APP } from './subworkflows/shiny_app'

workflow {
  List<String> available_GRNInference_tools = ["grnboost2"]
  List<String> available_DGRNInference_tools = ["boostdiff","zscores"]  
  List<String> tools = params.tools.split(",")
  
  diffgrn_tools = []
  grn_tools = []

  tools.each { element -> 
    if(available_DGRNInference_tools.contains(element)) {
      if(params.grn_mode) {
        throw new Exception(element + " cannot be used in grn_mode! Only GRN inference tools can be used in grn_mode.")
      }
      diffgrn_tools.add(element)
    } else if(available_GRNInference_tools.contains(element)) {
      grn_tools.add(element)
    } else {
      throw new Exception(element + " is not available. Please choose only a selection of the following tools: " + (available_GRNInference_tools + available_DGRNInference_tools))
    }
  }

  data = LOAD_DATA(tools)
  networks = RUN_TOOLS(data, tools)
  networks.view()
  // ANALYSIS(networks)
  SHINY_APP(diffgrn_tools, grn_tools)

  // // // TODO: Make NF pipeline creation of shiny App work with variable no. input network files and different data input formats
  // if (params.mode == "seurat") {
  // }
}

include { CREATE_METACELLS } from './data_loading/create_metacells'

workflow LOAD_DATA {
  take:
    tools

  main:
    if (params.mode == "seurat") {
      // TODO: Add some checks that verify a correct settings of the parameters 
      // Transforms the config.yaml file into a list of the following structure [key, assay, covariate_config, covariate_config_keys, output_folder, filter_id, filter_value]
      selection = []
      covariate_configs_keys = []        
      int selec_index_start = 0
      int selec_index_end = 1
      for (comparison in params.comparison) {
        comparison.each { key, value ->

          selection.add([])
          id = value.comparison_id
          selection[selec_index_start].add(id)
          assay = value.assay
          selection[selec_index_start].add(assay)

          covariate_config = value.covariate_config
          covariate_config_asColumnSeparatedList = []
          output_files = []    
          int cov_index = 0 
        
          // Function of this: adds every combination possible of the covarite config -> duplicates list if a list is given as a value of a key and then proceeds to create every combination
          covariate_config_keys = []
          covariate_config.each { selecKey, selecVals ->
            covariate_config_keys.add([])
            covariate_config_asColumnSeparatedList.add([])
            selecVals.each { selection_id, selection_value ->
              if (selection_id == "output_file") {
                output_files.add(selection_value)
              } else {
                covariate_config_keys[cov_index].add(selection_id)
                if (selection_value instanceof List) {
                  if (covariate_config_asColumnSeparatedList[cov_index].isEmpty()) {
                    for (elem in selection_value) {
                      covariate_config_asColumnSeparatedList[cov_index].add(elem)
                    }
                  } else {
                    for (int j = 0; j <= selection_value.size()-2; j++) {
                      covariate_config_asColumnSeparatedList[cov_index].add(covariate_config_asColumnSeparatedList[cov_index][0])
                    }
                    for (int j = 0; j <= covariate_config_asColumnSeparatedList[cov_index].size()-1; j++) {
                      covariate_config_asColumnSeparatedList[cov_index][j] = covariate_config_asColumnSeparatedList[cov_index][j] + ":" + selection_value[j]
                    }
                  }
                } else {
                  if (covariate_config_asColumnSeparatedList[cov_index].isEmpty()) {
                    covariate_config_asColumnSeparatedList[cov_index].add(selection_value)
                  } else {
                    for (int j = 0; j <= covariate_config_asColumnSeparatedList[cov_index].size()-1; j++) {
                      covariate_config_asColumnSeparatedList[cov_index][j] = covariate_config_asColumnSeparatedList[cov_index][j] + ":" + selection_value
                    }
                  }
                }
              }
            }
            // Check if covariate_config keys are the same inside one comparison -> Throw error if it is not the case 
            if (covariate_config_keys.size() > 1) {
              if (!covariate_config_keys[0].equals(covariate_config_keys[1])) {
                throw new Exception("The Covariate configuration keys have to be the same inside one comparison!")
              }
            }
            cov_index = cov_index + 1
          }

          // Function of this: Duplicates list for every covariate configuration and enters the proper values for each list
          if (covariate_config_asColumnSeparatedList.size() > 1) {
            selection.add([*selection[selec_index_start]])
            selec_index_end = selec_index_end + 1
          }          
          for (int i = 0; i <= covariate_config_asColumnSeparatedList.size()-1; i++) {         
            selection[selec_index_start + i].add(covariate_config_keys[i].join(":"))
            selection[selec_index_start + i].add(covariate_config_asColumnSeparatedList[i].join(","))
            selection[selec_index_start + i].add(output_files[i])
          }

          // adding two null values for name and id of a filter to match the input cardinality if the filter value is given
          if (value.filter == null) {
            for (int i = selec_index_start; i < selec_index_end; i++) {
                selection[i].add(value.filter)
                selection[i].add(value.filter)
              }
          } else {
            filter = value.filter
            filter.each { filterKey, filterVal ->
              for (int i = selec_index_start; i < selec_index_end; i++) {
                if (filterVal instanceof List) {
                  selection[i].add(filterVal.join(":"))
                } else {
                  selection[i].add(filterVal)
                }
              }
            }
          }
        }
        selec_index_start = selec_index_end
        selec_index_end = selec_index_start + 1
      }

      input_case_ch = Channel.fromList(selection)
      data = CREATE_METACELLS(params.mode, params.input, input_case_ch)
    } else if(params.mode == "tsv"){

      // input_files = params.input.split(",")
      // if (input_files.size() == 1) {
      //   List<String> available_DGRNInference_tools = params.available_DGRNInference_tools.split(",")  
      //   assert !available_DGRNInference_tools.contains(tools) : "You have only provided one tsv file, but selected a tool for DGRN inference! Provide two tsv files or only use GRN inference."
      //   data = CHECK_INPUT_TSVFILES(params.comparison_id, params.input_file1, null)
      // }
      data = Channel.of([params.comparison_id, [params.input_file1, params.input_file2]])
      data = CREATE_METACELLS(params.mode, data, null)
    } else {
      throw new Exception('Please select one of the following modes for the provided data: "seurat", "tsv"')
    }

  emit:
    data
}
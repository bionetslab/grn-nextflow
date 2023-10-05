include { CREATE_SHINY_APP } from '../modules/shiny_app'

workflow SHINY_APP {

	take:
        diffgrn_tools
        grn_tools

    main:
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
                selection[selec_index_start + i].add(output_files[i])
                selection[selec_index_start + i].add(covariate_config_asColumnSeparatedList[i].join(","))
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
        
        for (int i = 0; i <= selection.size()-1; i++) {
			    selection[i] = selection[i].join(",")
        }
        selection = selection.join("-")

        println grn_tools
        println diffgrn_tools
        println selection
        if (grn_tools.size == 0) {
          grn_tools = ""
        } else {
          grn_tools = grn_tools.join(",")
        }

        if (diffgrn_tools.size == 0) {
          diffgrn_tools = ""
        } else {
          diffgrn_tools = diffgrn_tools.join(",")
        }
        CREATE_SHINY_APP(selection, diffgrn_tools, grn_tools)
        // selection_ids = params.data_loading_seurat.selection_ids_toUse.split(",")
        // selection = []
        // params.data_loading_seurat.selection.eachWithIndex { element, index ->    
        //     id = "s" + index
        //     if (id in selection_ids) {
        //     selection.add([])
        //     selection[selection.size-1] = params.data_loading_seurat.selection[id].values()
        //     }
        // }
        // CREATE_SHINY_APP(
//             selection.collect { it.join(",") }.join("-")
//         )

  emit:
    selection
}
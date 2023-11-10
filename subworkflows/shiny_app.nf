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
                        // duplicate list #length(list)-1 times
                        cov_before = covariate_config_asColumnSeparatedList[cov_index]
                        for (int j = 0; j <= selection_value.size()-2; j++) {
                          covariate_config_asColumnSeparatedList[cov_index] += cov_before
                        }
                        // add the entries of the selction_values list to the duplicated entries
                        for (int j = 0; j <= covariate_config_asColumnSeparatedList[cov_index].size()-1; j++) {
                          covariate_config_asColumnSeparatedList[cov_index][j] = covariate_config_asColumnSeparatedList[cov_index][j] + ":" + selection_value[j % selection_value.size()]
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
                cov_index = cov_index + 1
              }

              // Duplicates list for every covariate configuration and enters the proper values for each list
              if (covariate_config_asColumnSeparatedList.size() > 1) {
                selection.add([*selection[selec_index_start]])
                selec_index_end = selec_index_end + 1
              }          
              for (int i = 0; i <= covariate_config_asColumnSeparatedList.size()-1; i++) {         
                selection[selec_index_start + i].add(covariate_config_keys[i].join(":"))
                selection[selec_index_start + i].add(output_files[i])
                selection[selec_index_start + i].add(covariate_config_asColumnSeparatedList[i].join(","))
              }
            }
            selec_index_start = selec_index_end
            selec_index_end = selec_index_start + 1
        }
        
        for (int i = 0; i <= selection.size()-1; i++) {
			    selection[i] = selection[i].join(",")
        }
        selection = selection.join("-")
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
        if (params.mode == "seurat") {
          seurat_file = params.input
        } else if (params.mode == "anndata") {
  	  	  seurat_file = "./seurat_object.rds"
        } else if (params.mode == "tsv") {
          selection = params.comparison_id
          seurat_file = 'NA'
        }

        CREATE_SHINY_APP(seurat_file, selection, diffgrn_tools, grn_tools, params.mode)
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
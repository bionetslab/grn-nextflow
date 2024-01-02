include { SELECT_DATA } from './data_loading/create_metacells'
include { CONVERT_ANNDATA_TO_SEURAT } from '../modules/data_loading/'

workflow LOAD_DATA {
	take:
    	tools

  	main:
	// TODO: Write code for reading the configuration file in case params.grn_mode == True
    if (params.mode == "seurat") {
      // TODO: Add some checks that verify a correct settings of the parameters 
      selection = read_config()
      input_case_ch = Channel.fromList(selection)
      data = SELECT_DATA(params.mode, params.input, input_case_ch)
    } else if(params.mode == "tsv") {
      
      data_ch = Channel.of([params.comparison_id, [params.input_file1, params.input_file2]])
      data = SELECT_DATA(params.mode, data_ch, null)

    } else if(params.mode == "anndata") {
    
	  // This currently can only handle raw inputs and no filters/aggregation etc. -> not better than a tsv file
	  // TODO: Figure out how the conversion converts column names and additional information such as UMAP plots ... (need a toy dataset for this)
		selection = read_config()
		input_case_ch = Channel.fromList(selection)
	  	seurat_file = CONVERT_ANNDATA_TO_SEURAT(params.input)
      	data = SELECT_DATA(params.mode, seurat_file, input_case_ch)

    } else {
      throw new Exception('Please select one of the following modes for the provided data: "seurat", "tsv"')
    }  
	emit:
    	data
}

def read_config() {
	selection = []
	covariate_configs_keys = []        
	selec_index_start = 0
	selec_index_end = 1
	print(params.comparison)
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
						println(selection_value)
						if (selection_value instanceof List) {
							print('test')
							if (covariate_config_asColumnSeparatedList[cov_index].isEmpty()) {
								for (elem in selection_value) {
									covariate_config_asColumnSeparatedList[cov_index].add(elem)
								}
								print(covariate_config_asColumnSeparatedList)
							} else {
								// duplicate list #length(list)-1 times
								cov_before = covariate_config_asColumnSeparatedList[cov_index]
								for (int j = 0; j <= selection_value.size()-2; j++) {
									covariate_config_asColumnSeparatedList[cov_index] += cov_before
								}
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
				// print(covariate_config_asColumnSeparatedList)
				// Check if covariate_config keys are the same inside one comparison -> Throw error if it is not the case 
				// if (covariate_config_keys.size() > 1) {
				// // if (!covariate_config_keys[0].equals(covariate_config_keys[1])) {
				// // 	throw new Exception("The Covariate configuration keys have to be the same inside one comparison!")
				// // }
				// }
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
	    	// if (value.filter == null) {
	    	//   for (int i = selec_index_start; i < selec_index_end; i++) {
	    	//       selection[i].add(value.filter)
	    	//       selection[i].add(value.filter)
	    	//     }
	    	// } else {
	    	//   filter = value.filter
	    	//   filter.each { filterKey, filterVal ->
	    	//     for (int i = selec_index_start; i < selec_index_end; i++) {
	    	//       if (filterVal instanceof List) {
	    	//         selection[i].add(filterVal.join(":"))
	    	//       } else {
	    	//         selection[i].add(filterVal)
	    	//       }
      		// 	}
      		//   }
    	}
		selec_index_start = selec_index_end
    	selec_index_end = selec_index_start + 1
	}
	println selection
	return selection
}

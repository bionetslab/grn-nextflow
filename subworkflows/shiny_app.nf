include { CREATE_SHINY_APP } from '../modules/shiny_app'

workflow SHINY_APP {

    main:
        selection_ids = params.data_loading_seurat.selection_ids_toUse.split(",")
        selection = []
        params.data_loading_seurat.selection.eachWithIndex { element, index ->    
            id = "s" + index
            if (id in selection_ids) {
            selection.add([])
            selection[selection.size-1] = params.data_loading_seurat.selection[id].values()
            }
        }
        CREATE_SHINY_APP(
            selection.collect { it.join(",") }.join("-")
        )
}
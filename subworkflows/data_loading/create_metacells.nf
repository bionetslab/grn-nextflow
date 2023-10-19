include { CREATE_METACELLS_SEURAT } from '../../modules/data_loading/'
include { CREATE_METACELLS_TSVFILES } from '../../modules/data_loading/'
include { CREATE_METACELLS_ANNDATA } from '../../modules/data_loading/'
include { CHECK_FILES } from '../../modules/data_loading/'

workflow SELECT_DATA {
    take:
        mode
        input
        selection
        
    main:
        if (params.create_metacells) {
            if (mode == "tsv") {
                data_case_ch = CREATE_METACELLS_TSVFILES(input.transpose(), mode)            
            } else if (mode == "seurat") {
                data_case_ch = CREATE_METACELLS_SEURAT(input, selection, mode)
            } else if (mode = "anndata") {
                data_case_ch = CREATE_METACELLS_ANNDATA(input, mode)
            }
            data_case_ch.view()

            tuple_ch = data_case_ch.groupTuple()
            tuple_ch.view { "value: $it" }
        } else {
            tuple_ch = data
        }
        
        checked_ch = CHECK_FILES(tuple_ch)
        checked_ch.view()

    emit:
        checked_ch
}
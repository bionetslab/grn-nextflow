include { SELECT_DATA_SEURAT } from '../../modules/data_loading/'
include { SELECT_DATA_TSVFILES } from '../../modules/data_loading/'
include { CHECK_FILES } from '../../modules/data_loading/'

workflow CREATE_METACELLS {
    take:
        mode
        input
        selection
        
    main:
        if (mode == "tsv") {
            data_case_ch = SELECT_DATA_TSVFILES(input.transpose(), mode)            
        } else if (mode == "seurat") {
            data_case_ch = SELECT_DATA_SEURAT(input, selection, mode)
        }
        data_case_ch.view()

        tuple_ch = data_case_ch.groupTuple()
        tuple_ch.view { "value: $it" }

        checked_ch = CHECK_FILES(tuple_ch)
        checked_ch.view()

    emit:
        checked_ch
}
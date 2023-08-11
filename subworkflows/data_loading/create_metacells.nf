include { SELECT_DATA } from '../../modules/data_loading/'
include { CHECK_FILES } from '../../modules/data_loading/'

workflow CREATE_METACELLS {
    take: 
        seurat_object
        selection
        
    main:
        data_case_ch = SELECT_DATA(seurat_object, selection)
        data_case_ch.view()

        tuple_ch = data_case_ch.groupTuple()
        tuple_ch.view { "value: $it" }

        checked_ch = CHECK_FILES(tuple_ch)
        checked_ch.view()

    emit:
        checked_ch
}
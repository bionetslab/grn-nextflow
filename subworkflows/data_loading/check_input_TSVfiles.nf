include { CHECK_FILES } from '../../modules/data_loading/'

workflow CHECK_INPUT_TSVFILES {
    take:
        key
        file1
        file2
    
    main:
        data_ch = Channel.of([key, [file1, file2]])
        data_ch.view()
        data = CHECK_FILES(data_ch)

    emit:
        data
}

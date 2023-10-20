include { RUN_TOOL } from '../../modules/tools/diffcoex'

workflow DIFFCOEX {
    take:
        data
    
    main:
        diffcoex = RUN_TOOL(data, params.diffcoex.top_n_edges)

    emit:
        diffcoex

}
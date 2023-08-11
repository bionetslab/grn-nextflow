include { RUN_TOOL } from '../../modules/tools/diffcoex'

workflow DIFFCOEX {
    take:
        data
        runs
    
    main:
        diffcoex = RUN_TOOL(data, runs)
        diffcoex_grouped = diffcoex.groupTuple()
}
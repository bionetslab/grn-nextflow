include { RUN_TOOL } from '../../modules/tools/grnboost2'
include { AGGREGATE_RESULTS } from '../../modules/tools/grnboost2'
include { WRITE_NETWORK_TO_DIR } from '../../modules/post_processing'

workflow GRNBOOST2 {
    take:
        data
        runs
    
    main:
        grnboost2 = RUN_TOOL(data, runs)
        grnboost2_grouped = grnboost2.groupTuple()

        grnboost2_grouped.view()
        aggregated_grnboost2 = AGGREGATE_RESULTS(data.join(grnboost2_grouped), params.n_runs, params.grnboost2.top_n_edges)
        WRITE_NETWORK_TO_DIR("grnboost2", aggregated_grnboost2)

    emit:
        aggregated_grnboost2
}
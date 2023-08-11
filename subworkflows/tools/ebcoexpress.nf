include { RUN_TOOL } from '../../modules/tools/ebcoexpress'
// include { AGGREGATE_RESULTS } from '../../modules/tools/ebcoexpress'

workflow EBCOEXPRESS {
    take:
        data
        runs
    
    main:
        ebcoexpress = RUN_TOOL(data, runs)
        ebcoexpress_grouped = ebcoexpress.groupTuple()

        // aggregated_ebcoexpress = AGGREGATE_RESULTS(data.join(ebcoexpress_grouped), params.n_runs, boostdiff_filtering_settings)
}
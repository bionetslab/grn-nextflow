include { RUN_TOOL } from '../../modules/tools/boostdiff'
include { AGGREGATE_RESULTS } from '../../modules/tools/boostdiff'


workflow BOOSTDIFF {
    take:
        data
        runs
    
    main:
        boostdiff_run_settings = // [n_estimators, n_features, n_subsamples, n_processes]
        [
            params.boostdiff.n_estimators,
            params.boostdiff.n_features,
            params.boostdiff.n_subsamples,
            params.boostdiff.n_processes
        ] 
        boostdiff_filtering_settings = // [top_n_targets, top_n_edges]
        [
            params.boostdiff.top_n_targets,
            params.boostdiff.top_n_edges
        ]
        
        boostdiff = RUN_TOOL(data, boostdiff_run_settings, params.use_tf_list, runs)
        boostdiff_grouped = boostdiff.groupTuple()

        boostdiff_network = AGGREGATE_RESULTS(data.join(boostdiff_grouped), params.n_runs, boostdiff_filtering_settings)

    emit:
        boostdiff_network[1]
    
}
include { RUN_TOOL } from '../../modules/tools/z_score'

workflow Z_SCORE {
    take:
        data
    
    main:
        z_score = RUN_TOOL(data, params.z_score.top_n_edges)
    
    emit:
        z_score
}
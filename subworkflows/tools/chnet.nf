include { RUN_TOOL } from '../../modules/tools/chNet'

workflow CHNET {
    take:
        data
    
    main:
        chnet = RUN_TOOL(data, params.chnet.lambda)
    
    emit:
        chnet
}
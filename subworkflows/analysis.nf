include { PLOT_GRN } from '../modules/analysis'

workflow ANALYSIS {
    
    take:
        networks

    main:
        PLOT_GRN(networks.transpose())
        

}
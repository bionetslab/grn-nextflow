// Define parameters with default values
params.input_path = null
params.output_dir = "results"
params.grn_tflist = 'all'

include { CSV_TOH5AD } from '../toh5ad/main.nf'
include { ADATA_TORDS } from '../tords/main.nf'

process LEAP {
    conda  "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/bioconductor-dcanr_bioconductor-singlecellexperiment_r-base_r-igraph_pruned:ae4f89e3d52cade5'
    publishDir params.output_dir, mode: 'copy'
    
    input:
    tuple val(meta), path(sample_rds)

    output:
    path "*_leap.csv"                   , emit: grn_results



    script:
    prefix = meta.id
    layer = "raw"


    template 'leap.R'
}

workflow {
    // Check if required parameters are provided
    if (!params.input_path) {
        error "Please provide --input_path parameter"
    }
    
    // Create input channel from the provided path
    input_ch = Channel.fromPath(params.input_path, checkIfExists: true)
                     .map { path -> [['id': path.baseName], path] }

    toh5ad_ch = CSV_TOH5AD(input_ch)
    tords_ch = ADATA_TORDS(toh5ad_ch)
    
    // Run the process
    LEAP(tords_ch)
}

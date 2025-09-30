// Define parameters with default values
params.input_path = null
params.output_dir = "results"
params.grn_tflist = 'all'

include { CSV_TOH5AD } from '../toh5ad/main.nf'

process GENIE3 {
    conda  "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/arboreto_dask-expr_distributed_python-abi3_pruned:dcf8d03d9c2dd880'
    publishDir params.output_dir, mode: 'copy'

    input:
    tuple val(meta), path(sample_h5ad)

    output:
    path "*.csv"                   , emit: grn_results

    script:
    // Define all variables for the template
    prefix = meta.id
    layer ="raw"
    
    template 'genie3.py'
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
    
    // Run the process
    GENIE3(toh5ad_ch)
}
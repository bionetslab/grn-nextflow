process CSV_TOH5AD {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/anndata_pandas:f8bc3ee0962e32af"

    publishDir params.output_dir, mode: 'copy'

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad

    script:
    prefix = meta.id
    counts_layer = 'X'

    template 'toh5ad.py'
}

workflow {
    // Check if required parameters are provided
    if (!params.input_path) {
        error "Please provide --input_path parameter"
    }
    
    // Create input channel from the provided path
    input_ch = Channel.fromPath(params.input_path, checkIfExists: true)
                     .map { path -> [['id': path.baseName], path] }
    
    // Run the process
    CSV_TOH5AD (input_ch)
}


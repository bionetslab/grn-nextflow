process ADATA_TORDS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/anndata2ri_bioconductor-singlecellexperiment_anndata_r-seurat:693764013a2a63f5'

    publishDir params.output_dir, mode: 'copy'

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.rds"), emit: rds

    script:
    prefix = meta.id
    counts_layer = 'X'

    template 'tords.py'
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
    ADATA_TORDS(input_ch)
}


nextflow.enable.dsl=2

workflow {
  convert_h5ad_to_pidc_inputs(
    file("input/filtered_placeholder.h5ad")
  )

  run_pidc(
    convert_h5ad_to_pidc_inputs.out.expression
  )
}

process convert_h5ad_to_pidc_inputs {
  tag "PIDC: Convert Input"

  // === CONTAINER MODE === 
  container "community.wave.seqera.io/library/anndata_numpy_pandas_pip_pruned:c5c0484dba547ab3"

  // === CONDA MODE ===
  // conda "modules/grn/pidc/pidc_env.yml"

  input:
  path expression_h5ad

  output:
  path "expression.tsv", emit: expression

  script:
  def script_path = "${workflow.projectDir}/modules/grn/pidc/h5ad_to_pidc_inputs.py"
  """
  python3 $script_path \\
    --input ${expression_h5ad} \\
    --output expression.tsv
  """
}

process run_pidc {
  tag "PIDC: Run"

  // === CONTAINER MODE === 
  container "community.wave.seqera.io/library/anndata_numpy_pandas_pip_pruned:c5c0484dba547ab3"

  // === CONDA MODE ===
  // conda "modules/grn/pidc/pidc_env.yml"

  publishDir "results", mode: 'copy'

  input:
  path "expression.tsv"

  output:
  path "pidc_ranked_edges.csv"

  script:
  def script_path = "${workflow.projectDir}/modules/grn/pidc/run_pidc.py"
  """
  python3 $script_path --input expression.tsv --output pidc_ranked_edges.csv
  """
}


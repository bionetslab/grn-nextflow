nextflow.enable.dsl=2

workflow {

  convert_h5ad_to_scode_inputs(
    file("input/filtered_placeholder.h5ad")
  )

  run_scode(
    convert_h5ad_to_scode_inputs.out.expression,
    convert_h5ad_to_scode_inputs.out.pseudotime
  )
}

process convert_h5ad_to_scode_inputs {
  tag "convert_inputs"

  // === USE CONTAINER MODE === 
  container "community.wave.seqera.io/library/anndata_numpy_pandas_pip_pruned:56d43cedd30669c7"

  // === USE CONDA MODE ===
  // conda "modules/grn/scode/scode_env.yml"

  input:
  path expression_h5ad

  output:
  path "expression.txt", emit: expression
  path "pseudotime.txt", emit: pseudotime

  script:
  def script_path = "${workflow.projectDir}/modules/grn/scode/h5ad_to_scode_inputs.py"
  """
  python3 $script_path \\
    --input ${expression_h5ad} \\
    --expr_out expression.txt \\
    --time_out pseudotime.txt
  """
}

process run_scode {
  tag "SCODE_run"

  // === USE CONTAINER MODE === 
  container "community.wave.seqera.io/library/anndata_numpy_pandas_pip_pruned:56d43cedd30669c7"

  // === USE CONDA MODE ===
  // conda "modules/grn/scode/scode_env.yml"

  publishDir "results", mode: 'copy'

  input:
  path "expression.txt"
  path "pseudotime.txt"

  output:
  path "ranked_edges.csv"

  script:
  """
  ruby ${workflow.projectDir}/modules/grn/scode/run_R.rb \\
    expression.txt \\
    pseudotime.txt \\
    results \\
    100 50 10 100 5

  cp results/ranked_edges.csv ranked_edges.csv
  """
}


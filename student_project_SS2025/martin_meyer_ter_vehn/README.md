# GRN Inference Workflows with Nextflow

This project implements four algorithms for gene regulatory network (GRN) inference from single-cell transcriptomic data in standalone Nextflow workflows.  
The workflows are prepared for later integration into the nf-core [`scdownstream`](https://nf-co.re/scdownstream) pipeline.  

Implemented algorithms:
- **GENIE3**
- **GRNBoost2**
- **LEAP**
- **ppcor**

## Usage 

Run a workflow with:

```bash
cd {algorithm}
nextflow run . --input_path {path/to/expression_matrix} --output_dir {results_dir}
```

**Examples:**
```bash
# GENIE3 algorithm
cd genie3
nextflow run . --input_path "../example_data/sim_GRN_100.csv" --output_dir "./results"

# Multiple files with wildcards
cd grnboost2  
nextflow run . --input_path "../example_data/sim_GRN_100_Phyla*.csv" --output_dir "./results"
```

The workflows are run in seperate conda environments by default.
To run the workflows in docker use:

```bash
cd {algorithm}
nextflow run . --input_path {path/to/expression_matrix} --output_dir {results_dir} --conda.enabled=false --docker.enabled=true
```

or set the corresponding configuration settings in `{algorithm}/nextflow.config`

```groovy
docker.enabled = true
conda.enabled = false
```

The result will be at `{results_dir}/{expression_matrix}_{algorithm}.csv`


## Data Format Conversion Utilities

The pipeline includes two utility workflows for data format conversion that are used internally by all algorithms but can also be run standalone:

### CSV to H5AD Conversion (`toh5ad`)

Converts CSV expression matrices to AnnData H5AD format:

```bash
cd toh5ad
nextflow run . --input_path {path/to/expression_matrix.csv} --output_dir {results_dir}
```

**Example:**
```bash
cd toh5ad
nextflow run . --input_path "../example_data/sim_GRN_100.csv" --output_dir "./results"
```

### H5AD to RDS Conversion (`tords`) 

Converts H5AD files to R RDS format for R-based algorithms:

```bash
cd tords
nextflow run . --input_path {path/to/file.h5ad} --output_dir {results_dir}
```

**Example:**
```bash
cd tords
nextflow run . --input_path "../results/sim_GRN_100.h5ad" --output_dir "./results"
```

### Internal Usage

These conversion workflows are automatically executed by the GRN inference algorithms:

- **Python algorithms** (GENIE3, GRNBoost2): CSV → H5AD
- **R algorithms** (LEAP, ppcor): CSV → H5AD → RDS

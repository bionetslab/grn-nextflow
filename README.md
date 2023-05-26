# boostdiff-nextflow
A juicy, albeit not yet fully configurable nextflow workflow for Boostdiff using a Seurat object as input.

# Installation Instruction 
1) Install nextflow with the following command (can be moved to any directory you want): (requires Java version 11 or higher)   
```
curl -fsSL get.nextflow.io | bash
```
2) Install the neccesary conda environment (this will take time). This can be sped up dramatically by using mamba instead of conda:  
```
conda env create --prefix=$(pwd)/boostdiff-wf --file=./boostdiff-wf.yml
conda activate ./boostdiff-wf
```   
3) Clone the `boostdiff-nextflow` git:
```
git clone git@github.com:gihannagalindez/boostdiff_inference.git && cd boostdiff-nextflow
pip install .
cd ..
```
Now you are set to run the boostdiff nf pipeline!

# Running the boostdiff pipeline
Run the follwing command inside the `boostdiff-nextflow` directory: (Replace `$(path_to_nextflow)` with the path to the directory that contains your nextflow installation)
```
$(path_to_nextflow)/nextflow run boostdiff.nf
```

If you change something in the settings or the run is interrupted use (reuses all the previously performed computations):
```
$(path_to_nextflow)/nextflow run boostdiff.nf -resume
```

For information about how to change the pipeline settings look into `boostdiff.nf` under:
```
workflow {
  // How to configure this pipeline:
 
  /* GLOBAL VARIABLES 
  params.seurat_object = "$projectDir/../data/ga_an0228_10x_deepseq_filtered_smarta_merged_tissue_integrated_rep_timepoint_infection_filtered_seurat.rds"
    => Path to the data
  params.column_name = 'infection:tissue:subject:time' 
    => This describes the selection criteria column_names for the data points
  params.cluster_name = 'cluster'
    => This describes the cluster name
  params.publish_dir = "$projectDir/../results/"
    => Sets the path to the output results
  */ 

  /* This is the selection variable:
     [Selected_covariate_configuation, 'Output_file_name', 'Key_matching_two_files_for_boostdiff', 'Output_folder', 'secondary_selection_criterion', 'categories']
      1. Selected_covariate_configuation: This refers to the params.colum_name criterion. There must be the same number of entries as the params.colum_name string (lists not allowed)
      2. Output file name for this data subset. Must be unique within the output folder
      3. Boostdiff requires two input files. To match the inputs in the pipeline the keys (third column) need to be matching and unique accross the pipeline run.
      4. Output folder for the run
      5. Secondary selection criterion allowing for multi selection. Only one column can be selected but from this column multiple categories (see 6)
      6. Categories to select from the column defined in 5.
  */
  .
  .
  .
```
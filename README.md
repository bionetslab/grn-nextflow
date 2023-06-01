# boostdiff-nextflow
A juicy, albeit not yet fully configurable nextflow workflow for Boostdiff using a Seurat object as input.

# Installation Instruction 
1) Install nextflow with the following command (can be moved to any directory you want): (requires Java version 11 or higher)   
```
curl -fsSL get.nextflow.io | bash
```
2) Clone this repository and navigate to it:
```
git clone git@github.com:bionetslab/boostdiff-nextflow.git && cd boostdiff-nextflow
```
2) Install the neccesary conda environment (this will take time). This can be sped up dramatically by using mamba instead of conda:  
```
conda env create --prefix=$(pwd)boostdiff-wf --file=./boostdiff-wf.yml
conda activate ./boostdiff-wf
```   
3) Clone the `boostdiff_inference` tool (https://github.com/gihannagalindez/boostdiff_inference):
```
git clone git@github.com:gihannagalindez/boostdiff_inference.git && cd boostdiff_inference
pip install .
cd ..
```
Now you are set to run the boostdiff nf pipeline!

# Running the boostdiff pipeline
Run the follwing command inside the `boostdiff-nextflow` directory: (Replace `$(path_to_nextflow)` with the path to the directory that contains your nextflow installation)
```
$(path_to_nextflow)/nextflow run boostdiff.nf
```
The results are saved to `results/`. Note that only symlinks to the `work/` folder of nextflow are saved in the `results/` folder. If you want to upload your results to a git repository you need to copy the results (with following symlinks) into a folder and upload this folder. 
If you change something in the settings or the run is interrupted use (reuses all the previously performed computations):
```
$(path_to_nextflow)/nextflow run boostdiff.nf -resume
```

For information about the configuration of the pipeline run:
```
$(path_to_nextflow)/nectflow run boostdiff.nf --help
```

# Settings of the pipeline
Standard settings of this pipeline:
- Data:
  - Cluster 1, 2
  - Armstrong vs Docile, Spleen, day 28
  - Armstrong vs Docile, Liver, day 10
- No. total runs of boostdiff: 10
- Settings for individual boostdiff runs:
  - no. estimators: 50
  - no. features: 1500
  - no. subsamples: 30
  - no. processes: 8
- Settings for filtering the aggregated results:
  - Top n nodes: 20 (most differntially expressed target genes between the two conditions)
  - Top n edges: 100 (highest ranking interactions between remaining genes)

# Pipeline workflow
1) Read in data
2) Run boostdiff for a no. total runs
3) Aggregate results by creating the union of all runs and average over the scores
4) Filter aggregated results based on the settings
5) Check regulatory interaction of every edge based on small linear model that is fitted on every edge 
6) Create output .html file

# Interpreting the results
Outputs:
  - This pipeline puts out a .txt file containing the data of the inferred differntial GRN
  - This pipeline puts out a .html file containing the graph representation of the inferred differential GRN (image shows part of such a differential GRN):
    - Nodes: represent the genes (annotated with the gene name)
    - Edges: 
      - 2 colours representing condition 1,2 (pink and green, see legend)
      - 4 possible edges:
        - pink fully drawn arrow  (up regulatory interaction that is stronger in condition 1)
        - pink dashed arrow       (down regulatory interaction that is stronger in condition 1)
        - green fully drawn arrow (up regulatory interaction that is stronger in condition 2)
        - green dashed arrow      (down regulatory interaction that is stronger in condition 2)

![diff_grn](diff_grn_example.png)



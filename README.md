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

If you change something in the settings use (reuses all the previously performed computations):
```
$(path_to_nextflow)/nextflow run boostdiff.nf -resume
```
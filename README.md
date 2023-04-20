# boostdiff-nextflow
A juicy, albeit not yet fully configurable nextflow workflow for Boostdiff using a Seurat object as input.

# Installation Instruction 
1) Install nextflow with the following command (can be moved to any directory you want): (requires Java version 11 or higher)   
```
curl -fsSL get.nextflow.io | bash
```  
2) Clone the `boostdiff-nextflow` git:
```
git clone https://github.com/bionetslab/boostdiff-nextflow.git && cd boostdiff-nextflow
```
3) Install the neccesary conda environment (this will take time):  
```
conda env create --prefix=$(pwd)/boostdiff-wf --file=./boostdiff-wf.yml
``` 
4) Install boostdiff:
```
git clone https://github.com/gihannagalindez/boostdiff_inference.git  && cd boostdiff_inference
pip install .
```
Now you are set to run the boostdiff nf pipeline!

# Running the boostdiff pipeline
Run the follwing command inside the `boostdiff-nextflow` directory: (Replace `$(path_to_nextflow)` with the path to the directory that contains your nextflow installation)
```
$(path_to_nextflow)/nextflow run boostdiff.nf -with-conda $(pwd)/boostdiff-wf
```
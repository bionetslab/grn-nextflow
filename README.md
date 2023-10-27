# Inter-Net Xplorer
## Installation Instruction 
1) Install nextflow with the following command (can be moved to any directory you want): (requires Java version 11 or higher)   
```
curl -fsSL get.nextflow.io | bash
```
2) Clone the boostdiff repository and navigate to it:
```
git clone git@github.com:bionetslab/boostdiff-nextflow.git && cd boostdiff-nextflow
```
2) Install the neccesary conda environment (this will take time). This can be sped up dramatically by using mamba instead of conda:  
```
mamba env create --prefix=$(pwd)/environment --file=./environment.yml
conda activate ./environment
```   
3) Clone the `boostdiff_inference` tool (https://github.com/gihannagalindez/boostdiff_inference):
```
git clone git@github.com:gihannagalindez/boostdiff_inference.git && cd boostdiff_inference
pip install .
cd ..
```
Now you are set to run Inter-Net Xplorer!

## Interpretation of computed networks of Inter-Net Xplorer
See this pdf [Network Explanation](https://github.com/bionetslab/grn-nextflow/blob/dgrn_nf/man/network_explanation.pdf)

## Running Inter-Net Xplorer
This section describes piece by piece how to run this pipeline. Examples will be provided along the way and at the end of this section.
### 1) Running the nextflow pipeline:
To run the the nextflow pipeline use the following command and swap out the parameters to fit your needs. The next sections go over the paramters in detail 
```
${path_to_nextflow}/nextflow run main.nf --tools=${tools_to_run} --mode=${data_mode} --input=${data_input} -params-file ${config_file} --publish_dir=${output_data_path}
```

### 2) Setting the `--tools` paramter:
The `--tools` parameter needs to be set to identify the tools that are used in the pipeline. Current available tools are:
* DGRN inference tools:
  * `boostdiff` (https://github.com/gihannagalindez/boostdiff_inference)
  * `zscores` (https://doi.org/10.1371/journal.pone.0009202)
  * `diffcoex` (https://doi.org/10.1186/1471-2105-11-497)
* GRN inference tools:
  * `grnboost2` (https://academic.oup.com/bioinformatics/article/35/12/2159/5184284)

The `--tools` parameter needs to be set as comma separated list. For example, if you want to use boostdiff and grnboost2 you need to set `--tools=boostdiff,grnboost2`   

### 3) Setting the `--mode` parameter:
The `--mode` parameter needs to be set to identify the data that you are using. Currently availabe modes are `seurat`, `tsv` and `anndata`. 

### 4) Setting `--input` parameter:
The full path has to be set for all input files! ALL values in the columns that are used for selection in the configuration file must not contain ",", "-", ":".  

#### 4.1) If `--mode=seurat`:
Use the `--input` parameter to set the path to the seurat file. **The file type must be .RDS**. If you are using this mode, you need to provide a configuration file with the `-params-file` parameter that contains information about the grouping/filtering that should be done in the Seurat object for your specific needs. See `example_config.yaml` for instructions and an example on how to write a config file for your dataset.

#### 4.2) If `--mode=tsv`:
Do not set the `--input` parameter! <br />
Use the `--input_file1` parameter to set the path to the first tsv file. <br />
Use the `--input_file2` parameter to set the path to the second tsv file. <br />
<!-- If you are only using GRN inference tools, specifying one input is enough. <br /> -->
The first column of the tsv files **has to be named** `Gene` and contain all gene names. The following columns represent the samples. If you are using this mode, you need to set the `--comparison_id` parameter. This needs to be an identifiable string because the folder with results will be named after this. You **do not** need to set the `-params-file` parameter as there cannot be done any grouping/filtering on the tsv files!  

#### 4.3) If `--mode=anndata`:
Set the `--input` paramter to the path to the AnnData object. **The file type must be .h5ad**. The AnnData object will be converted to a Seurat object in the pipeline. If you are using this mode, you need to provide a configuration file with the `-params-file` parameter that contains information about the grouping/filtering that should be done in the AnnData object for your specific needs. See `example_config.yaml` for instructions and an example on how to write a config file for your dataset.

### 5) Setting the `--publish_dir` parameter:
This parameter sets the path to the results folder where the results/outputs should be written. This folder must exist!

### 6) Optional Parameters
1) `--create_metacells`: Default value: TRUE <br />
   Determines whether metacells should be created or not. If you do not use metacells, the computation runtime of the implemented tools is really long.
2) `--work`: Default value: Path to folder where you start the the nextflow pipeline. <br />
   Change this if you want the internal nextflow files to be stored somewhere else. The internal nextflow files can be quite big, so be careful if you have limited disk usage.
3) `--n_runs`: Default value: 10 <br />
   Determines how often the tools are run that rely on randomization (boostdiff, grnboost2). This is done to improve the robustness of these tools.
4) See `nextflow.config` for all tool specific and nextflow specific parameters.  

## Further Information
README is WIP: Information to come:
1) Structure of the pipeline
2) Instructions/Example on how to extend the pipeline with a tool/analyses

<!-- ## Structure of the pipeline:
The pipeline is split into 3 steps:
1) Data loading
2) Running tools
3) Analysis
Every step has concretely defined inputs and outputs.
### 1) Data loading:
Inputs:
  * tools: The tools used inside the pipeline defined via the `--tools` parameter
  
Output:
  * data: Nextflow channel with the structure [comparison_id, [file_1.tsv, file_2.tsv]]
Function:
  This step loads the data and performs some basic checks
### 2) Running Tools
Inputs:
  * data: The key and the input files loaded in the first step of the pipeline as a nextflow channel
  * tools: The tools used inside the pipeline defined via the `--tools`parameter

Outputs:
  * networks: Nextflow channel, where all networks are grouped based on the keys
### 3) Analysis
  Inputs:
    * networks: Nextflow channel, where all networks are grouped based on the keys (computed by step 2) -->

<!-- ## How to extend the pipeline
Suppose you want to extend the pipeline with a new workflow/module for any of the steps of the pipeline. This section goes over the steps needed to accomplish this. For any information on nextflow functions, we refer to the official nextflow documentation: TODO: ADD_URL
### 1) Adding a data loading workflow/module
  1) **This step is only needed if you want to integrate a new data format into the pipeline.** <br /> Go to `subworkflows/load.data.nf` -> There you will find the workflow that calls the specific data loading module based on the mode given in the `--mode` parameter. -> Add a new case for your new data mode. Write a config parser if needed for your case. Use the `-params-file` as input parameter for the config file. An exemplary config parser can be seen for `params.mode=="seurat"`. Finally, call the `CREATE_METACELLS(mode, input)` workflow. Please try to use the parameters `--input` or `--input_file1`, `--input_file2` as input parameters. Every run of the pipeline has to have a unique ID that is used as identifier and name for the run. If you are using a config file, this config file has to have such an id (named key in the pipeline). If you do not have a config file, you can use the `--comparison_id` parameter for this.
  2) Go to `subworkflows/data_loading/create_metacells.nf`. This workflow creates a nextflow channel for the key and the inputs and calls the correct process for the data_lading based on the selected `--mode`. <br /> 
     2.1 **Only needed if you are adding a new mode!**: Add another case for your newly created mode. Goto step 2.3. <br />
     2.2 Use the `--data_loading` parameter that is used to select the data loading process to add a new case for your new data loading process inside the correct mode case. Where and how to add a new process will be described in section 3. <br />
     2.3 Create a nextflow channel for the inputs (given as parameter to this workflow) and the config file if you are using one. As examples, look at already implemented cases.
     2.4 Choose a name for your new process. This name has to be in upper case and it should follow the structure `SELECT_DATA_${--mode}_${--data-loading}` The process needs to be included at the top using: 
     ```
     include { YOUR_PROCESS_NAME } from '../../modules/data_loading/' 
     ```
     The process gets the input and config, if you are using one, channel as inputs. 
  3) Go to `modules/data_loading/main.nf`. This script contains the processes for the data loading. Add your new process with the correct process name defined in step 2.4.   -->

<!-- # Settings of the pipeline
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

 -->

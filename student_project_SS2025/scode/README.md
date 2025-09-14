# 🧬 SCODE Gene Regulatory Network Inference Pipeline (Nextflow)

This project implements the **SCODE (Sparsity-promoting Contextual Orthogonal Decomposition)** algorithm for gene regulatory network (GRN) inference from single-cell RNA-seq data using a **Nextflow DSL2 pipeline**. It accepts `.h5ad` files as input and outputs a ranked list of regulator-target gene pairs.

The implementation follows the structure outlined in the [Beeline GRN Benchmarking Suite](https://github.com/Murali-group/Beeline) and combines **Python**, **R**, and **Ruby** for complete GRN reconstruction.

---

## 📁 Project Structure

```
SCODE_Nextflow_Complete/
├── main.nf                         # Main Nextflow pipeline script
├── nextflow.config                 # Configuration file (profile settings)
├── expression.txt                 # Expression matrix (auto-generated)
├── pseudotime.txt                 # Pseudotime vector (auto-generated)
├── results/                       # Output directory
├── input/
│   └── filtered_placeholder.h5ad  # Input .h5ad data file
├── modules/
│   └── grn/
│       └── scode/
│           ├── scode_env.yml             # Conda environment file
│           ├── h5ad_to_scode_inputs.py   # Converts .h5ad to SCODE input
│           ├── SCODE.R                   # R-based SCODE implementation
│           ├── averageA.R                # Averages SCODE networks
│           └── run_R.rb                  # Ruby wrapper for SCODE execution
```

---

## 🔧 Setup Instructions

### 🧬 Prerequisites

Make sure you have the following installed:

- [Nextflow](https://www.nextflow.io/) ≥ 22.x
- [Miniconda](https://docs.conda.io/en/latest/)
- Python ≥ 3.10 with:

  - `scanpy`, `anndata`, `numpy`, `pandas`

- R ≥ 4.2 with:

  - `optparse`, `R.utils`, `ggplot2`, `devtools`, `curl`

- Ruby ≥ 2.5

---

## ✅ Environment Setup (Conda Mode)

```bash
git clone https://github.com/YOUR_USERNAME/SCODE_Nextflow_Complete.git
cd SCODE_Nextflow_Complete

# Create the Conda environment
conda env create -f modules/grn/scode/scode_env.yml -p ./scode_env

# Activate it
conda activate ./scode_env
```

---

## 🚀 Run the SCODE Pipeline

### Option 1: 🧪 Using Conda (Recommended if container build fails)

```bash
NXF_CONDA_CACHEDIR=./.conda nextflow run main.nf -with-conda -resume
```

### Option 2: 📦 Using Container (Seqera Wave)

If your container is correctly built and published on [Wave](https://wave.seqera.io), run:

```bash
nextflow run main.nf -with-wave -resume
```

🛠️ **Important**: In your `main.nf`, switch between container and conda mode by commenting/uncommenting:

```groovy
// For container:
container "community.wave.seqera.io/library/your_final_image:tag"
// For conda:
conda "modules/grn/scode/scode_env.yml"
```

> ⚠️ If the container mode fails due to image not found or build failure, please verify the dependency compatibility in your Conda environment YAML. Container build errors often arise from **conflicting versions of R/Python packages** across Conda channels (e.g., `conda-forge` vs `bioconda`). It is not a fault in your Nextflow code.

---

## 📤 Output Format

The output will be a ranked gene–gene edge list:

```
Source   Target   Rank   [Score]
```

### 🧾 Example:

```
9889,100,88,0
9890,100,89,0
9891,100,90,0
```

---

## 📸 Sample Output Screenshot

After successfully running the pipeline, you will see something like this:

![SCODE Result Screenshot](results/SCODE_output_example.png)

> _After Successfully Run the project You got this screen_

---

## 🧊 Post-processing Tips

You can further explore the results by:

- Mapping index IDs back to gene names using `.h5ad`'s `var_names`
- Visualizing the GRN using network tools (e.g., NetworkX, Cytoscape)
- Plotting expression heatmaps or SCODE matrix with `matplotlib` or `ggplot2`

---

## ⚙️ nextflow\.config (Optional)

```groovy
profiles {
  conda {
    process {
      conda = true
    }
  }
  container {
    process {
      container = "community.wave.seqera.io/library/your_final_image:tag"
    }
  }
}
```

> 💡 Not using this file is okay. You can manage execution from `main.nf` directly using process-level `conda` or `container` declarations.

---

## 📚 Reference

- **Matsumoto et al. (2017)**
  _SCODE: An efficient regulatory network inference method from single-cell RNA-seq during differentiation_
  [Bioinformatics, Volume 33, Issue 15](https://doi.org/10.1093/bioinformatics/btx194)

- [Beeline GRN Benchmarking Suite](https://github.com/Murali-group/Beeline)

---

## 👨‍💻 Author

**Prosenjit Chowdhury**
M.Sc. Artificial Intelligence – FAU Erlangen-Nürnberg
Working Student @ SAP ERP PCX
Erlangen, Germany
🔗 GitHub: [@prosenjit-chowdhury](https://github.com/prosenjit-chowdhury)

---

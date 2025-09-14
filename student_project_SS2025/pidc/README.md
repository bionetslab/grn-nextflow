# 🧬 PIDC Gene Regulatory Network Inference Pipeline (Nextflow)

This project implements the **PIDC (Partial Information Decomposition and Context)** algorithm for gene regulatory network (GRN) inference from single-cell RNA-seq data using a **Nextflow DSL2 pipeline**.  
It accepts `.h5ad` input files and outputs a ranked list of regulator-target gene pairs.

The implementation follows the structure outlined in the [Beeline GRN Benchmarking Suite](https://github.com/Murali-group/Beeline) and uses Python and Julia to perform GRN inference.

---

## 📁 Project Structure

```

pidc\_nextflow/
├── main.nf                        # Main Nextflow pipeline script
├── input/
│   └── filtered\_placeholder.h5ad # Input .h5ad data file
├── results/                       # Output directory
├── modules/
│   └── grn/
│       └── pidc/
│           ├── h5ad\_to\_pidc\_inputs.py  # Converts .h5ad to PIDC input (TSV)
│           └── run\_pidc.py             # Runs PIDC algorithm on TSV

```

---

## 🔧 Setup Instructions

### 🧬 Prerequisites

Make sure you have the following installed:

- **Nextflow ≥ 22.x**
- **Docker** or **Singularity** (for container-based execution)
- **Python ≥ 3.10** with:
  - `scanpy`, `anndata`, `numpy`, `pandas`
- **Julia ≥ 1.6** with:
  - `NetworkInference.jl`, `LightGraphs.jl`

---

## ✅ Environment Setup (Optional - Conda Mode)

```bash
# Clone the repo
git clone https://github.com/YOUR_USERNAME/pidc-nextflow-grn-pipeline.git
cd pidc-nextflow-grn-pipeline

# Create the Conda environment (only if you want to run with conda)
conda env create -f modules/grn/pidc/pidc_env.yml -p ./pidc_env

# Activate it
conda activate ./pidc_env
```

---

## 🚀 Run the PIDC Pipeline

### Option 1: 🧪 Using Conda (if container build fails)

```bash
NXF_CONDA_CACHEDIR=./.conda nextflow run main.nf -with-conda -resume
```

### Option 2: 📦 Using Container (Recommended)

```bash
nextflow run main.nf -with-wave -resume
```

> ℹ️ The container used is:
> `community.wave.seqera.io/library/anndata_numpy_pandas_pip_pruned:c5c0484dba547ab3`

---

## 🛠️ Container vs Conda

In your `main.nf`, you can toggle modes:

```nextflow
// Container mode
container "community.wave.seqera.io/library/anndata_numpy_pandas_pip_pruned:c5c0484dba547ab3"

// Conda mode
// conda "modules/grn/pidc/pidc_env.yml"
```

⚠️ If the container build fails due to dependency mismatch, switch to Conda and manually install Julia + required packages.
Container errors are usually **not your fault** — they stem from version incompatibilities.

---

## 📤 Output Format

The output will be a ranked gene–gene edge list with confidence scores.

**Format:**

```
Gene1    Gene2    Score
```

**Example:**

```
DENND2D    RCN2     1.974
RCN2       DENND2D  1.974
PLEKHA1    DEXI     1.972
```

Each line represents a predicted regulatory edge.

---

## 📸 Sample Output Screenshot

📂 After successfully running the pipeline, you’ll get a file like:

![PIDC Output Screenshot](results/PIDC_Output-1-1.png)

---

## 🧊 Post-processing Tips

You can further explore the results by:

- Removing duplicate bidirectional edges (A→B and B→A)
- Visualizing the GRN using tools like **NetworkX** or **Cytoscape**
- Comparing to ground truth if available
- Overlaying expression heatmaps from `.h5ad` onto the network

---

---

## 📚 References

- **Chan et al. (2017)**
  _Gene Regulatory Network Inference from Single-Cell Data Using Multivariate Information Measures_
  [Cell Systems](https://doi.org/10.1016/j.cels.2017.08.014)

- **Beeline GRN Benchmarking Suite**
  [https://github.com/Murali-group/Beeline](https://github.com/Murali-group/Beeline)

---

## 👨‍💻 Author

**Prosenjit Chowdhury**
_M.Sc. Artificial Intelligence – FAU Erlangen-Nürnberg_
_Working Student @ SAP ERP PCX | Erlangen, Germany_
🔗 GitHub: [@prosenjit-chowdhury](https://github.com/prosenjit-chowdhury)

```

```

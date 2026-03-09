# $\ell_1$-spectral clustering simulations

This repository contains the experimental workflow, simulation scripts, and data preprocessing pipelines for the manuscript "$\ell_1$-spectral clustering algorithm: a spectral clustering method using $\ell_1$-regularization". 

The goal of this repository is to provide a transparent, end-to-end path to reproduce the figures and tables presented in the study, using the `l1spectral` R package.

---

## Repository Structure

* **`l1spectralclustering.Rproj`**: R Project file.
* **`scripts/`**: R scripts for the full analysis.
    * `01_preprocessing.R`: Implements 75% variance filtering and 0.7 correlation network construction.
    * `02_simulations.R`: Benchmarking scripts (reproduces Table 1 and Eigengap analysis).
    * `03_tcga_analysis.R`: Real-world application on TCGA-BRCA data.
* **`data/`**: Contains processed datasets.
* **`results/`**: Contains the results to generate figures and tables.
* **`DESCRIPTION`**: Lists all R dependencies and versions.

---

## Getting Started

### 1. Requirements
Ensure you have [RStudio](https://rstudio.com/) installed. This project relies on the following key packages:
- `l1spectral`
- `WGCNA`
- `igraph`

### 2. Setup
Clone the repository and open the `.Rproj` file in RStudio. Then, install dependencies using `devtools`:

```r
install.packages("devtools")
devtools::install_deps()
```

### 3. Execution

Run the scripts in the `scripts/` folder in numerical order. Each script is self-contained but assumes the previous data-cleaning steps have been completed.
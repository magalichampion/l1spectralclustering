# $\ell_1$-spectral clustering simulations

This repository contains the experimental workflow, simulation scripts, and data preprocessing pipelines for the manuscript " $\ell_1$-spectral clustering algorithm: a spectral clustering method using $\ell_1$-regularization". 
Its goal is to provide an end-to-end path to reproduce the figures and tables presented in the study.

---

## Repository Structure

* **`l1spectralclustering.Rproj`**: R Project file.
* **`scripts/`**: R scripts for the full analysis.
    * `Create_simulated_data`: Generates synthetic graphs with ground-truth community structures across varying number of nodes, densities and noise levels.
    * `Run_simulations.R`: Benchmarks the $\ell_1$-spectral clustering against state-of-the-art methods.
    * `Compute_performance`: Calculates statistical metrics including Adjusted Mutual Information (AMI).
    * `Plot_results.R`: Produces comparative visualizations to evaluate methods performance.
* **`python/`**: Python implementations.
* **`data/`**: Processed datasets.
* **`results/`**: Results to generate figures and tables.
* **`DESCRIPTION`**: All R dependencies and versions.
* **`vignettes`**: Comprehensive step-by-step tutorials.

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
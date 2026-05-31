# HNSCclassifier: An R Package to Predict Molecular Subtypes of Head and Neck Squamous Cell Carcinoma

<!-- badges: start -->
[![R-CMD-check](https://github.com/Ronlee12355/HNSCclassifier/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Ronlee12355/HNSCclassifier/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

An R package for robust molecular subtype classification of **Head and Neck Squamous Cell Carcinoma (HNSCC)** using the TCGA-derived taxonomy (Atypical, Basal, Classical, Mesenchymal). The classifier integrates pathway‑level normalisation via single‑sample gene set enrichment analysis (ssGSEA) and a pre‑trained random forest model to ensure cross‑platform portability.

## Features

- **Command‑line function**: `classifyHNSC()` for batch prediction from a gene expression matrix.
- **Shiny web interface**: `classifyHNSC_interface()` for interactive upload and exploration.
- **Flexible input**: accepts gene symbols, Ensembl, Entrez, or RefSeq identifiers (automatic conversion with `org.Hs.eg.db`).
- **Output options**: predicted subtype labels or posterior class probabilities.
- **Cross‑platform**: designed to work with RNA‑seq (TPM/FPKM) and microarray data after log2 transformation.

---

## Installation

The package depends on Bioconductor packages (`GSVA`, `org.Hs.eg.db`,`randomForest`). Install them first, then install `HNSCclassifier` from GitHub.

```r
# Install BiocManager if not already available
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c("GSVA", "org.Hs.eg.db", "AnnotationDbi", 'shiny','DT','shinythemes','randomForest'))

# Install HNSCclassifier from GitHub
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("Ronlee12355/HNSCclassifier")
```

## References
- Cancer Genome Atlas Network. Comprehensive genomic characterization of head and neck squamous cell carcinomas. Nature. 2015;517(7536):576-582. doi:10.1038/nature14129.
- Walter V, Yin X, Wilkerson MD, et al. Molecular subtypes in head and neck cancer exhibit distinct patterns of chromosomal gain and loss of canonical cancer genes. PLoS One. 2013;8(2):e56823. doi:10.1371/journal.pone.0056823.
---
title: "HNSCclassifier: Predict Molecular Subtypes of Head and Neck Squamous Cell Carcinoma"
---

# HNSCclassifier: Predict Molecular Subtypes of Head and Neck Squamous Cell Carcinoma

<!-- badges: start -->
<a href="https://www.gnu.org/licenses/gpl-3.0"><img src="https://img.shields.io/badge/License-GPLv3-blue.svg" height="20"></a>
<a href="https://www.r-project.org/"><img src="https://img.shields.io/badge/R-%E2%89%A5%204.3-brightgreen" height="20"></a>
<a href="https://github.com/Ronlee12355/HNSCclassifier/releases"><img src="https://img.shields.io/github/v/release/Ronlee12355/HNSCclassifier" height="20"></a>
<a href="https://www.bioconductor.org/"><img src="https://img.shields.io/badge/Bioconductor-GSVA%7Corg.Hs.eg.db-blue" height="20"></a>
<!-- badges: end -->

<img src="https://github.com/Ronlee12355/HNSCclassifier/blob/main/HNSCclassifier.png" height="230" align="right"/>
An R package for robust molecular subtype classification of **Head and Neck Squamous Cell Carcinoma (HNSCC)**. This classifier assigns tumour samples to one of four TCGA-defined molecular subtypes — **Atypical**, **Basal**, **Classical**, and **Mesenchymal** — using a pre-trained random forest model coupled with pathway-level normalisation via single-sample gene set enrichment analysis (ssGSEA).

> **Why pathway-level normalisation?** Gene expression data from different platforms (RNA-seq vs. microarray, different batches) can exhibit substantial technical variation at the individual gene level. By converting gene-level expression to pathway-level enrichment scores via ssGSEA, the classifier captures higher-level biological signals that are more reproducible across platforms, ensuring cross-study and cross-platform portability.
---

## Table of Contents

- [Molecular Subtypes of HNSCC](#molecular-subtypes-of-hnscc)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Features](#features)
- [Classification Workflow](#classification-workflow)
- [Built-in Datasets](#built-in-datasets)
- [Shiny Web Interface](#shiny-web-interface)
- [References](#references)

---

## Molecular Subtypes of HNSCC

Head and Neck Squamous Cell Carcinoma is a heterogeneous disease with distinct molecular subtypes that differ in their biology, clinical outcomes, and therapeutic responses. The four subtypes, originally characterised by The Cancer Genome Atlas (TCGA) and Walter et al., are:

| Subtype | Key Characteristics |
|---|---|
| **Atypical** | Enriched for HPV-positive tumours; improved prognosis; p16<sup>INK4A</sup> overexpression; CDKN2A silencing less common |
| **Basal** | Expression patterns resembling basal epithelial cells; enrichments in epidermal development and extracellular matrix organisation genes |
| **Classical** | Heavy smoking association; the most prevalent subtype; characterised by xenobiotic metabolism, KEAP1/NRF2 pathway alterations, and oxidative stress gene signatures |
| **Mesenchymal** | Epithelial–mesenchymal transition (EMT) features; invasive/migratory phenotype; TGF-β signalling activation; poorest prognosis among the four subtypes |

Subtype information can guide prognosis stratification and may inform treatment selection in the context of clinical trials and translational research.

---

## Installation

The package depends on several Bioconductor packages. Install them first, then install `HNSCclassifier` from GitHub.

### Step 1: Install Bioconductor dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GSVA", "org.Hs.eg.db", "AnnotationDbi"))
```

### Step 2: Install optional dependencies (for the Shiny app)

```r
BiocManager::install(c("shiny", "DT", "shinythemes", "shinyjs"))
```

### Step 3: Install HNSCclassifier

```r
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("Ronlee12355/HNSCclassifier")
```

### Step 4: Load the package

```r
library(HNSCclassifier)
# Startup message: 'HNSCclassifier v0.1.0'
```

---

## Quick Start

### Basic classification

```r
library(HNSCclassifier)

# Load example dataset
data(TCGA_LUSC)

# Predict molecular subtypes (output: class labels)
subtypes <- classifyHNSC(TCGA_LUSC, outputType = "class")
table(subtypes)
#   Atypical      Basal   Classical Mesenchymal
#          7         15         20          8

# Get posterior probabilities for each subtype
probs <- classifyHNSC(TCGA_LUSC, outputType = "prob")
head(probs)
```

### Using different gene identifier types

```r
# With Ensembl IDs
data(TCGA_LUSC_ENSEMBL)
res <- classifyHNSC(TCGA_LUSC_ENSEMBL, idType = "ENSEMBL")
table(res)
```

### Launch the Shiny app

```r
classifyHNSC_interface()
```

---

## Features

- **Command-line interface**: [`classifyHNSC()`](#quick-start) for batch prediction from a numeric gene expression matrix — ideal for processing large cohorts or integrating into downstream pipelines.
- **Interactive Shiny app**: [`classifyHNSC_interface()`](#shiny-web-interface) for users who prefer a graphical, point-and-click experience.
- **Multi-platform identifier support**: accepts gene symbols, Ensembl IDs, Entrez IDs, or RefSeq accession numbers. Automatic conversion is performed via `org.Hs.eg.db`.
- **Flexible output**: return predicted subtype **class labels** (`"class"`) or a **probability matrix** (`"prob"`) showing posterior probabilities for each of the four subtypes.
- **Cross-platform robustness**: designed to work with RNA-seq (TPM, FPKM, or normalised counts) and microarray expression data after log2 transformation.
- **Automatic preprocessing**: input validation, log2 transformation detection and application, and gene-level noise reduction (zero-variance gene filtering) are handled internally.
- **Curated pathway gene sets**: built with a carefully selected collection of pathway gene sets (`required.sets`) for ssGSEA, bridging any input dataset to the model's training feature space.
- **Pre-trained random forest model**: ships with a fully trained `randomForest` object (`finalModel`) — no need for users to train or tune any model.

---

## Classification Workflow

The complete classification pipeline proceeds through five well-defined steps:

```
Input Expression Matrix (genes × samples)
        │
        ▼
┌─────────────────────────────────────┐
│  1. Input Validation                 │
│     • Check for NAs, negatives,      │
│       non-numeric columns            │
│     • Filter genes with zero         │
│       variance (MAD = 0)             │
│     • Validate sample names          │
└─────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────┐
│  2. Gene ID Conversion (if needed)   │
│     ENSEMBL / ENTREZID / REFSEQ      │
│           ────────► SYMBOL           │
│     via org.Hs.eg.db; for multiple   │
│     IDs mapping to the same symbol,  │
│     retain the maximum expression    │
└─────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────┐
│  3. Log2 Transformation              │
│     Heuristic check: if data are     │
│     not already log2-transformed,    │
│     apply log2(x + 1)                │
└─────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────┐
│  4. Pathway-Level Normalisation      │
│     GSVA::gsva(ssgseaParam(...))     │
│     with curated pathway gene sets   │
│     Gene-level → Pathway-level       │
│     enrichment scores (ssGSEA)       │
│     Ensures cross-platform portability│
└─────────────────────────────────────┘
        │
        ▼
┌─────────────────────────────────────┐
│  5. Random Forest Classification     │
│     • Center/scale pathway scores    │
│     • Predict with pre-trained model │
│     • Output: class or probability   │
└─────────────────────────────────────┘
        │
        ▼
   Subtype: Atypical / Basal / Classical / Mesenchymal
```

This design mirrors the approach described in the package methodology, where the transition from gene-level to pathway-level features is the critical step for achieving platform-independent classification.

---

## Built-in Datasets

The package includes two example datasets derived from the **TCGA Lung Squamous Cell Carcinoma (LUSC)** project ([GDC portal](https://portal.gdc.cancer.gov/projects/TCGA-LUSC)) to help users get started immediately:

| Dataset | Rows | Columns | Description |
|---|---|---|---|
| `TCGA_LUSC` | 22,962 genes (SYMBOL) | 50 tumour samples | Gene expression matrix in log2(TPM + 1) scale, with official gene symbols as row names |
| `TCGA_LUSC_ENSEMBL` | 22,962 genes (ENSEMBL) | 50 tumour samples | Same expression values as `TCGA_LUSC`, but row names are converted to Ensembl gene IDs (e.g., ENSG00000141510) |

> **Note:** These datasets are provided for demonstration and testing purposes. Users should apply the classifier to their own HNSCC gene expression data for research applications. The dataset names use "LUSC" (Lung Squamous Cell Carcinoma) because the data were originally sourced from that project, but the classifier was designed and validated for HNSCC samples.

---

## Shiny Web Interface

The package provides an interactive [Shiny](https://shiny.rstudio.com/) web application for users who prefer a graphical interface:

```r
library(HNSCclassifier)
classifyHNSC_interface()
```

The app features three tabs:

### 1. Predict Tab
- **Input options**: select gene identifier type (SYMBOL, ENSEMBL, ENTREZID, REFSEQ), output type (class labels or probabilities), and file separator (comma or tab)
- **File upload**: upload a CSV/TSV file (max 50 MB) containing a gene expression matrix
- **Interactive validation**: real-time checks on the uploaded data with clear error messages (missing Gene_ID column, NAs, negative values, etc.)
- **Data preview**: a 4 × 4 preview of the uploaded matrix
- **Downloadable example**: button to download an example CSV file (`HNSCclassifier_example.csv`) to understand the required input format
- **Classification result**: results displayed in an interactive [DT](https://rstudio.github.io/DT/) table with copy, CSV, and Excel export buttons

### 2. Tutorial Tab
Step-by-step instructions guiding users through:
1. Preparing your expression matrix (format: `Gene_ID` in the first column, sample IDs in subsequent columns)
2. Uploading and configuring options
3. Running the classification and interpreting results

### 3. Contact Tab
Author contact information and a link to the [GitHub Issues](https://github.com/Ronlee12355/HNSCclassifier/issues) page for bug reports and feature requests.

---

## References

1. Cancer Genome Atlas Network. Comprehensive genomic characterization of head and neck squamous cell carcinomas. *Nature*. 2015;517(7536):576-582. doi:[10.1038/nature14129](https://doi.org/10.1038/nature14129).

2. Walter V, Yin X, Wilkerson MD, et al. Molecular subtypes in head and neck cancer exhibit distinct patterns of chromosomal gain and loss of canonical cancer genes. *PLoS One*. 2013;8(2):e56823. doi:[10.1371/journal.pone.0056823](https://doi.org/10.1371/journal.pone.0056823).

3. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics*. 2013;14:7. doi:[10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7).

4. Breiman L. Random Forests. *Machine Learning*. 2001;45(1):5-32. doi:[10.1023/A:1010933404324](https://doi.org/10.1023/A:1010933404324).

---

## Contributing

Bug reports, feature suggestions, and contributions are welcome via the [GitHub Issues](https://github.com/Ronlee12355/HNSCclassifier/issues) page.

## License

This package is licensed under the [GNU General Public License v3 (GPL ≥ 3)](https://www.gnu.org/licenses/gpl-3.0).

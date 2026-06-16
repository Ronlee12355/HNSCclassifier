#' @keywords internal
#'
#' @section Overview:
#' **HNSCclassifier** provides a ready-to-use molecular subtyping
#' framework for head and neck squamous cell carcinoma (HNSCC).  It
#' assigns tumour samples to one of four TCGA-defined subtypes —
#' **Atypical (HPV-driven)**, **Basal**, **Classical**, and
#' **Mesenchymal** — using a pipeline designed for cross-platform
#' portability.
#'
#' The core methodological innovation is **pathway-level normalisation**
#' via single-sample Gene Set Enrichment Analysis (ssGSEA).  Rather than
#' operating on individual gene expression values — which are sensitive
#' to platform and batch effects — the classifier first projects each
#' sample onto a curated space of biological pathways, then applies a
#' pre-trained random forest model in this pathway space.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{classifyHNSC}}}{Core classifier.  Accepts a gene
#'     expression matrix and returns subtype labels or posterior
#'     probabilities.}
#'   \item{\code{\link{pathway_scores}}}{Extracts the intermediate ssGSEA
#'     pathway enrichment matrix for custom downstream analyses.}
#'   \item{\code{\link{classifyHNSC_interface}}}{Launches the interactive
#'     Shiny web application.}
#' }
#'
#' @section Visualisation:
#' \describe{
#'   \item{\code{\link{plot_subtype_probabilities}}}{Stacked bar chart of
#'     posterior subtype probabilities, with low-confidence flagging.}
#'   \item{\code{\link{plot_subtype_heatmap}}}{Annotated heatmap of
#'     subtype-driver gene expression.}
#'   \item{\code{\link{plot_confusion_matrix}}}{Heatmap visualisation of
#'     a \code{caret::confusionMatrix} for assessing classification
#'     performance against reference labels.}
#' }
#'
#' @section Built-in data:
#' \describe{
#'   \item{\code{\link{TCGA_LUSC}}}{50 TCGA-LUSC tumour samples
#'     (22,962 genes, gene symbols, log2(TPM+1)).}
#'   \item{\code{\link{TCGA_LUSC_ENSEMBL}}}{Same matrix with Ensembl gene
#'     IDs for demonstrating the \code{idType} argument.}
#' }
#'
#' @section Input requirements:
#' Expression data should be **non-negative** numeric values (e.g. TPM,
#' FPKM, or normalised counts) with genes in rows and samples in
#' columns.  The pipeline automatically detects whether a log2
#' transformation is needed and applies \code{log2(x + 1)} if so.
#' Gene identifiers can be SYMBOL (default), ENSEMBL, ENTREZID, or
#' REFSEQ; non-symbol types are auto-converted via
#' \pkg{org.Hs.eg.db}.
#'
#' @section Classification pipeline:
#' \enumerate{
#'   \item Input validation (NAs, negatives, zero-variance genes)
#'   \item Gene ID conversion to gene symbols (if needed)
#'   \item Log2 transformation detection and application
#'   \item Pathway scoring via \code{\link[GSVA]{gsva}} with ssGSEA
#'   \item Random forest classification with a pre-trained model
#' }
#'
#' @references
#' Cancer Genome Atlas Network (2015). Comprehensive genomic
#' characterization of head and neck squamous cell carcinomas.
#' \emph{Nature}, 517(7536), 576–582.
#'
#' Walter V, Yin X, Wilkerson MD, et al. (2013). Molecular subtypes
#' in head and neck cancer exhibit distinct patterns of chromosomal
#' gain and loss of canonical cancer genes. \emph{PLoS ONE}, 8(2),
#' e56823.
#'
#' Hänzelmann S, Castelo R, Guinney J (2013). GSVA: gene set
#' variation analysis for microarray and RNA-seq data.
#' \emph{BMC Bioinformatics}, 14, 7.
#'
#' @seealso
#' \itemize{
#'   \item \url{https://github.com/Ronlee12355/HNSCclassifier}
#'   \item Report bugs at
#'         \url{https://github.com/Ronlee12355/HNSCclassifier/issues}
#' }
#'
#' @author
#' **Maintainer**: Jiang Li \email{ronlee12355@outlook.com}
#'
#' Institute of Cell and Gene Technology,
#' Shenzhen University of Advanced Technology
#'
"_PACKAGE"

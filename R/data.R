#' TCGA LUSC gene expression subset (log2 TPM)
#'
#' A subset of the TCGA lung squamous cell carcinoma (LUSC) gene expression
#' dataset. It contains the first 50 primary tumor
#' samples, pre-processed to **log2(TPM + 1)** values. This dataset is used
#' for code examples and the package vignette.
#'
#' @format A numeric matrix with 22,962 rows (genes) and 50 columns (samples):
#' \describe{
#'   \item{rownames}{Official gene symbols (e.g., TP53, KRAS).}
#'   \item{colnames}{TCGA sample barcodes (e.g., TCGA-XX-XXXX-01A).}
#' }
#'
#' @source \url{https://portal.gdc.cancer.gov/projects/TCGA-LUSC}
#'
#' @details
#' The original TPM values were normalized and then log2-transformed using
#' `log2(x + 1)` to stabilize variance. The subset was created by selecting
#' the first 22,962 rows and 50 columns from the full processed matrix to
#' keep the package size small. All samples correspond to primary tumor
#' tissue (sample type code '01').
#'
#' Because the data are already log2-transformed, they will pass the
#' internal log2 check in `classifyHNSC()` and will not be re-transformed.
#'
#' @usage data(TCGA_LUSC)
#'
#' @examples
#' data(TCGA_LUSC)
#' # Check dimensions
#' dim(TCGA_LUSC)
#' # First few genes and samples
#' TCGA_LUSC[1:5, 1:3]
"TCGA_LUSC"


#' TCGA LUSC gene expression subset (Ensembl IDs, log2 TPM)
#'
#' @format A numeric matrix with 22,962 rows (Ensembl IDs) and 50 columns
#'   (primary tumor samples):
#' \describe{
#'   \item{rownames}{Ensembl gene identifiers (e.g., ENSG00000141510).}
#'   \item{colnames}{TCGA sample barcodes (e.g., TCGA-XX-XXXX-01A).}
#' }
#'
#' @source Derived from `TCGA_LUSC` (see `?TCGA_LUSC`). The original data were
#'   obtained from \url{https://portal.gdc.cancer.gov/projects/TCGA-LUSC}.
#'
#' @details
#' The matrix is identical of `TCGA_LUSC`, with rownames
#' converted to Ensembl IDs using `AnnotationDbi::mapIds` with
#' `keytype = "SYMBOL"` and `column = "ENSEMBL"`. Expression values are
#' identical to the corresponding rows of `TCGA_LUSC` (log2(TPM + 1)).
#' This dataset is intended to demonstrate the `idType = "ENSEMBL"` functionality of the
#' `classifyHNSC()` function.
#'
#' This dataset is used in the vignette and examples to show how
#' `classifyHNSC()` handles non-symbol input and performs automatic ID
#' conversion.
#'
#' @usage data(TCGA_LUSC_ENSEMBL)
#'
#' @examples
#' data(TCGA_LUSC_ENSEMBL)
#' dim(TCGA_LUSC_ENSEMBL)
#' # View first few Ensembl IDs
#' head(rownames(TCGA_LUSC_ENSEMBL))
"TCGA_LUSC_ENSEMBL"

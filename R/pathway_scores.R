#' Extract Pathway-Level ssGSEA Enrichment Scores
#'
#' This function performs input validation, optional gene ID conversion,
#' log2 transformation (if needed), and single-sample gene set enrichment
#' analysis (ssGSEA) on a gene expression matrix, then returns the
#' pathway-level enrichment score matrix. It exposes the intermediate
#' output that is normally consumed internally by `classifyHNSC()`.
#'
#' Useful for downstream analyses such as pathway-level differential
#' expression, gene set enrichment analysis (GSEA), clustering, or
#' custom machine learning experiments using the same pathway
#' representation that the classifier is built upon.
#'
#' @param input_expr A numeric gene expression matrix (or data frame) with
#'   **genes in rows** and **samples in columns**. Row and column names are
#'   mandatory. Expression values must be non-negative (e.g. TPM, FPKM,
#'   or normalised counts). If the data are not on a log2 scale, they will
#'   be automatically log2-transformed via `log2(x + 1)` internally.
#' @param idType Character string specifying the gene identifier type used
#'   in the rownames of `input_expr`. One of `"SYMBOL"` (default),
#'   `"ENSEMBL"`, `"ENTREZID"`, or `"REFSEQ"`. If not `"SYMBOL"`, the
#'   function converts rownames to gene symbols via **org.Hs.eg.db**.
#'
#' @return A numeric matrix with **pathway gene sets in rows** and
#'   **samples in columns**. Each value is the ssGSEA enrichment score
#'   for that pathway–sample pair. The pathway names correspond to the
#'   internal `required.sets` gene set collection curated for HNSCC
#'   subtype classification.
#'
#' @seealso \code{\link{classifyHNSC}} for the full classification
#'   pipeline, including the random forest step applied downstream of
#'   these scores.
#'
#' @export
#'
#' @importFrom GSVA gsva ssgseaParam
#' @importFrom stats aggregate mad
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#'
#' @examples
#' \dontrun{
#' data(TCGA_LUSC)
#' pw <- pathway_scores(TCGA_LUSC)
#' head(pw[, 1:5])            # first 5 samples, all pathways
#'
#' # With Ensembl IDs
#' pw_ens <- pathway_scores(TCGA_LUSC_ENSEMBL, idType = "ENSEMBL")
#' }
pathway_scores <- function(input_expr = NULL,
                           idType = 'SYMBOL') {
  ## ====== 1. Input checking ====== ##
  idType <- match.arg(idType, c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"))

  if (check_missing_values(input_expr)) {
    stop("Input expression matrix contains missing values.")
  }
  if (is.null(rownames(input_expr)) ||
      is.null(colnames(input_expr))) {
    stop("Rownames and colnames are mandatory in gene expression profiles.")
  }
  if (sum(apply(input_expr, 2, is.numeric)) != ncol(input_expr)) {
    stop("Only numeric values in gene expression profile are accepted.")
  }
  if (any(input_expr < 0, na.rm = TRUE)) {
    stop("Gene expression profiles cannot contain any negative value(s).")
  }
  if (ncol(input_expr) < 1) {
    stop("At least one sample is required.")
  }

  ## Remove zero-variance genes
  input_expr <- input_expr[apply(input_expr, 1, function(x)
    mad(x) > 0), ]

  ## ====== 2. Gene ID conversion (if needed) ====== ##
  if (idType != "SYMBOL") {
    mapping <- suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys      = rownames(input_expr),
        keytype   = idType,
        columns   = "SYMBOL"
      )
    )
    idx <- match(rownames(input_expr), mapping[[idType]])
    symbol_vec <- mapping$SYMBOL[idx]
    keep <- !is.na(symbol_vec)

    if (sum(keep) == 0) {
      stop("No identifiers could be mapped to gene symbols.")
    }
    expr_sub <- input_expr[keep, , drop = FALSE]
    sym_sub  <- symbol_vec[keep]

    agg <- stats::aggregate(as.data.frame(expr_sub),
                            by  = list(SYMBOL = sym_sub),
                            FUN = max)
    rownames(agg) <- agg$SYMBOL
    agg$SYMBOL <- NULL
    input_expr <- as.matrix(agg)
  }

  ## ====== 3. Log2 transformation ====== ##
  if (!Log2ed(input_expr)) {
    input_expr <- log2(input_expr + 1)
  }

  ## ====== 4. ssGSEA pathway scoring ====== ##
  input_expr_gsva <- GSVA::gsva(
    GSVA::ssgseaParam(exprData = as.matrix(input_expr),
                      geneSets = required.sets),
    verbose = FALSE
  )

  return(input_expr_gsva)
}

#' Classify Head and Neck Squamous Cell Carcinoma Samples into TCGA Subtypes
#'
#' This function takes a gene expression matrix, performs quality checks,
#' optional log2 transformation, single-sample gene set enrichment (ssGSEA)
#' at the pathway level, and then predicts the TCGA-based molecular subtype
#' of each sample using a pre-trained random forest model.
#'
#' @param input_expr A numeric gene expression matrix (or data frame) with
#'   **genes in rows** and **samples in columns**. Row and column names are
#'   mandatory. The expression values must be non-negative (e.g., TPM, FPKM,
#'   or normalized counts). If the data are not on a log2 scale, they will be
#'   automatically log2-transformed via `log2(x + 1)` internally.
#' @param idType Character string specifying the gene identifier type used in
#'   the rownames of `input_expr`. Acceptable values are `"SYMBOL"` (default),
#'   `"ENSEMBL"`, `"ENTREZID"`, or `"REFSEQ"`. If any type other than
#'   `"SYMBOL"` is provided, the function will convert the rownames to gene
#'   symbols using the **org.Hs.eg.db** annotation package. This requires
#'   **org.Hs.eg.db** and **AnnotationDbi** to be installed (they are listed
#'   under `Suggests` of the package). For multiple identifiers mapping to the
#'   same gene symbol, expression values are averaged.
#' @param outputType Character indicating the type of prediction result to
#'   return. Use `"class"` (default) to obtain the assigned subtype label for
#'   each sample, or `"prob"` to get a matrix of class probabilities (one row
#'   per sample, one column per subtype).
#'
#' @return
#' \itemize{
#'   \item If `outputType = "class"`, a named character vector of predicted
#'     subtypes (e.g., "Atypical", "Basal", "Classical", "Mesenchymal").
#'   \item If `outputType = "prob"`, a numeric matrix with samples in rows,
#'     subtypes in columns, and values ranging from 0 to 1 representing the
#'     predicted probability of each subtype.
#' }
#'
#' @details
#' The function proceeds through five main steps:
#' \enumerate{
#'   \item **Input validation**: Checks for missing values (`NA`), non-numeric
#'     columns, negative expression values, missing row/column names, and
#'     correct sample number.
#'     \item **Gene ID conversion** (if \code{idType} is not \code{"SYMBOL"}):
#'     Converts row identifiers to gene symbols using \code{org.Hs.eg.db}.
#'     If multiple identifiers map to the same gene symbol, the maximum
#'     expression value is retained for each sample.
#'   \item **Log2 transformation**: If the data are not already on a log2 scale
#'     (determined by an internal heuristic), they are transformed as
#'     `log2(x + 1)`.
#'   \item **Pathway-level representation**: Single-sample gene set enrichment
#'     analysis (ssGSEA) is performed using a curated set of pathway gene sets
#'     stored internally as `required.sets`. This step harmonizes the input
#'     data with the training cohort at the biological pathway level.
#'   \item **Random forest classification**: The pathway scores are centered,
#'     scaled, and fed into a pre-trained random forest model (`finalModel`) to
#'     generate the final subtype prediction.
#' }
#' Both `required.sets` and `finalModel` are internal package data and are
#' loaded automatically. The function requires the \pkg{GSVA} package.
#'
#' @note
#' If using `idType` values other than `"SYMBOL"`, the packages
#' **org.Hs.eg.db** and **AnnotationDbi** must be installed. You can install
#' them with:
#' \preformatted{
#' BiocManager::install("org.Hs.eg.db")
#' BiocManager::install("AnnotationDbi")
#' }
#' The conversion step may remove genes that lack a unique symbol mapping;
#' a warning will indicate the number of genes retained.
#'
#' @seealso
#' \code{\link[GSVA]{gsva}} for the underlying pathway scoring engine.
#'
#' @export
#'
#' @importFrom GSVA gsva ssgseaParam
#' @importFrom stats aggregate mad predict
#' @import randomForest
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#'
#' @examples
#' \dontrun{
#' # Example with gene symbols (default)
#' data(hnsc_example_expr)
#' subtypes <- classifyHNSC(hnsc_example_expr, outputType = "class")
#'
#' # Example with Ensembl IDs
#' # Ensure org.Hs.eg.db is installed
#' subtypes_ens <- classifyHNSC(ensembl_expr_matrix,
#'                              idType = "ENSEMBL",
#'                              outputType = "prob")
#' }
classifyHNSC <- function (input_expr = NULL,
                          idType = 'SYMBOL',
                          outputType = c("class", "prob")) {
  ## ====== 1. Input checking ====== ##
  outputType <- match.arg(outputType, c("class", "prob"))
  idType <- match.arg(idType, c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"))
  if (check_missing_values(input_expr)) {
    stop("Input expression matrix contains missing values.")
  }
  if (is.null(rownames(input_expr)) ||
      is.null(colnames(input_expr))) {
    stop('Rownames and colnames are madatory in gene expression profiles.')
  }
  if (sum(apply(input_expr, 2, is.numeric)) != ncol(input_expr)) {
    stop('Only numeric values in gene expression profile is accepted.')
  }
  if (any(input_expr < 0, na.rm = T)) {
    stop('Gene expression profiles cannot contain any negative value(s).')
  }
  if (ncol(input_expr) < 1) {
    stop('At least one sample is required for classification.')
  }

  # Before ssGSEA
  input_expr <- input_expr[
    apply(input_expr, 1, function(x){mad(x)>0}),
  ]

  ## ====== 2. ID conversion (if needed) ====== ##
  if (idType != "SYMBOL") {
    mapping <- suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = rownames(input_expr),
        keytype = idType,
        columns = "SYMBOL"
      )
    )

    # Align SYMBOL to row order
    idx <- match(rownames(input_expr), mapping[[idType]])
    symbol_vec <- mapping$SYMBOL[idx]
    keep <- !is.na(symbol_vec)
    if (sum(keep) == 0) {
      stop("No identifiers could be mapped to gene symbols.")
    }
    expr_sub <- input_expr[keep, , drop = FALSE]
    sym_sub <- symbol_vec[keep]

    # Aggregate by max per symbol
    agg <- stats::aggregate(as.data.frame(expr_sub),
                            by = list(SYMBOL = sym_sub),
                            FUN = max)
    rownames(agg) <- agg$SYMBOL
    agg$SYMBOL <- NULL
    input_expr <- as.matrix(agg)
  }

  ## ====== 3. Process input data ====== ##
  if (!Log2ed(input_expr)) {
    input_expr <- log2(input_expr + 1)
  }

  ## ====== 3. GSVA computation ====== ##
  input_expr_gsva <- GSVA::gsva(GSVA::ssgseaParam(exprData = as.matrix(input_expr), geneSets = required.sets),
                                verbose = FALSE)

  ## ====== 4. Classification ====== ##
  res <- predict(finalModel,
                 newdata = scale(t(input_expr_gsva)),
                 type = outputType)

  return(res)
}

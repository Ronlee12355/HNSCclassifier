#' Extract Differentially Expressed Genes per Subtype
#'
#' For each molecular subtype, performs gene-level differential expression
#' analysis comparing samples of that subtype against all other samples
#' (one-vs-rest). Returns genes passing user-specified fold-change and
#' adjusted p-value thresholds, sorted by significance.
#'
#' @param expr A numeric gene expression matrix with genes in rows and samples
#'   in columns. Row names must be gene symbols (or identifiers matching
#'   \code{idType}). Values should be log2-transformed (e.g. log2(TPM + 1)).
#' @param subtype A named character vector or factor giving the predicted
#'   subtype for each sample, as returned by
#'   \code{classifyHNSC(outputType = "class")}.
#' @param idType Character string specifying the gene identifier type.
#'   One of \code{"SYMBOL"} (default), \code{"ENSEMBL"}, \code{"ENTREZID"},
#'   or \code{"REFSEQ"}. Non-symbol IDs are auto-converted to gene symbols
#'   via \pkg{org.Hs.eg.db}.
#' @param method Character string specifying the statistical test.
#'   One of \code{"wilcox"} (default) for Wilcoxon rank-sum test or
#'   \code{"t.test"} for Welch's t-test.
#' @param log2fc_cutoff Numeric. Minimum absolute log2 fold-change to retain
#'   a gene (default: 1, i.e. 2-fold change).
#' @param p_cutoff Numeric. Adjusted p-value (FDR) threshold (default: 0.05).
#' @param min_samples Numeric. Minimum number of samples required in each
#'   group (subtype vs rest) to run the comparison (default: 3).
#'
#' @return A named list of data frames, one per subtype. Each data frame
#'   contains genes passing the thresholds, with columns:
#'   \itemize{
#'     \item \code{gene}: gene symbol
#'     \item \code{log2FC}: log2 fold-change (positive = up in this subtype)
#'     \item \code{mean_expr}: mean expression in this subtype
#'     \item \code{mean_other}: mean expression in all other subtypes
#'     \item \code{p_value}: raw p-value from the chosen test
#'     \item \code{adj_p_value}: Benjamini-Hochberg adjusted p-value (FDR)
#'   }
#'
#' @details
#' For each subtype, the function compares in-subtype samples against all
#' remaining samples on every gene. P-values are corrected across all genes
#' within each subtype using the Benjamini-Hochberg procedure.
#'
#' \code{method = "wilcox"} uses \code{\link[stats]{wilcox.test}} which makes
#' no distributional assumptions. \code{method = "t.test"} uses Welch's
#' t-test via \code{\link[stats]{t.test}} which assumes approximate normality.
#'
#' @importFrom stats p.adjust wilcox.test t.test setNames
#' @export
#'
#' @examples
#' \dontrun{
#' data(TCGA_LUSC)
#' subtypes <- classifyHNSC(TCGA_LUSC, outputType = "class")
#' degs <- extract_top_genes(TCGA_LUSC, subtypes,
#'                            log2fc_cutoff = 1, p_cutoff = 0.05)
#' head(degs$Basal)
#' }
extract_top_genes <- function(expr,
                               subtype,
                               idType = "SYMBOL",
                               method = c("wilcox", "t.test"),
                               log2fc_cutoff = 1,
                               p_cutoff = 0.05,
                               min_samples = 3) {

  # ---- 1. Input validation ----
  idType <- match.arg(idType, c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"))
  method <- match.arg(method)
  if (!is.numeric(log2fc_cutoff) || log2fc_cutoff < 0) {
    stop("'log2fc_cutoff' must be a non-negative number.")
  }
  if (!is.numeric(p_cutoff) || p_cutoff <= 0 || p_cutoff >= 1) {
    stop("'p_cutoff' must be between 0 and 1.")
  }

  validate_id_type(expr, idType)

  # ---- 2. Coerce subtype to named character vector ----
  if (is.atomic(subtype) && !is.null(names(subtype))) {
    subtype_vec <- as.character(subtype)
    names(subtype_vec) <- names(subtype)
  } else if (is.data.frame(subtype)) {
    if (!all(c("Sample", "Subtype") %in% colnames(subtype))) {
      stop("When 'subtype' is a data frame, it must have columns 'Sample' and 'Subtype'.")
    }
    subtype_vec <- stats::setNames(as.character(subtype$Subtype),
                                   as.character(subtype$Sample))
  } else {
    stop("'subtype' must be a named character/factor vector or a data frame.")
  }

  # ---- 3. Gene ID conversion (if needed) ----
  if (idType != "SYMBOL") {
    mapping <- suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = rownames(expr),
        keytype = idType,
        columns = "SYMBOL"
      )
    )
    # Align SYMBOL to row order
    idx <- match(rownames(expr), mapping[[idType]])
    symbol_vec <- mapping$SYMBOL[idx]
    keep <- !is.na(symbol_vec)
    if (sum(keep) == 0) {
      stop("No identifiers could be mapped to gene symbols.")
    }
    expr_sub <- expr[keep, , drop = FALSE]
    sym_sub <- symbol_vec[keep]

    # Aggregate by max per symbol (consistent with classifyHNSC)
    agg <- stats::aggregate(as.data.frame(expr_sub),
                            by = list(SYMBOL = sym_sub),
                            FUN = max)
    rownames(agg) <- agg$SYMBOL
    agg$SYMBOL <- NULL
    expr <- as.matrix(agg)
  }

  # ---- 4. Align expression with subtype ----
  common <- intersect(colnames(expr), names(subtype_vec))
  if (length(common) == 0) {
    stop("No common samples between expression matrix and subtype annotations.")
  }
  expr <- expr[, common, drop = FALSE]
  subtype_vec <- subtype_vec[common]

  # ---- 5. Choose test function ----
  test_fun <- switch(method,
    wilcox = function(x, y) stats::wilcox.test(x, y)$p.value,
    t.test = function(x, y) stats::t.test(x, y, var.equal = FALSE)$p.value
  )

  # ---- 6. One-vs-rest differential analysis per subtype ----
  subtype_levels <- sort(unique(subtype_vec))

  results <- lapply(stats::setNames(nm = subtype_levels), function(lev) {
    in_group  <- names(subtype_vec)[subtype_vec == lev]
    out_group <- names(subtype_vec)[subtype_vec != lev]

    if (length(in_group) < min_samples || length(out_group) < min_samples) {
      warning("Skipping '", lev, "': fewer than ", min_samples, " samples in one group ",
              "(in = ", length(in_group), ", out = ", length(out_group), ").")
      return(NULL)
    }

    mean_in   <- rowMeans(expr[, in_group,  drop = FALSE])
    mean_out  <- rowMeans(expr[, out_group, drop = FALSE])
    log2fc    <- mean_in - mean_out

    p_vec <- sapply(seq_len(nrow(expr)), function(i) {
      x <- as.numeric(expr[i, in_group])
      y <- as.numeric(expr[i, out_group])
      tryCatch(test_fun(x, y), error = function(e) 1)
    })

    adj_p <- stats::p.adjust(p_vec, method = "BH")

    df <- data.frame(
      gene         = rownames(expr),
      log2FC       = log2fc,
      mean_expr    = mean_in,
      mean_other   = mean_out,
      p_value      = p_vec,
      adj_p_value  = adj_p,
      stringsAsFactors = FALSE
    )

    df <- df[abs(df$log2FC) >= log2fc_cutoff & df$adj_p_value <= p_cutoff, , drop = FALSE]
    df <- df[order(abs(df$log2FC), decreasing = TRUE), , drop = FALSE]
    rownames(df) <- NULL
    df
  })

  results <- Filter(Negate(is.null), results)
  return(results)
}

#' Convert Gene Identifiers in an Expression Matrix
#'
#' Maps the row names of a gene expression matrix from one gene identifier
#' system to another using the \pkg{org.Hs.eg.db} annotation database.
#' Supports flexible aggregation strategies for many-to-one mappings, and
#' optionally retains unmapped genes with \code{NA} row names.
#'
#' @param expr A numeric gene expression matrix (or data frame) with genes
#'   in rows and samples in columns. Row names must be valid identifiers
#'   of the type specified by \code{from}.
#' @param from Character string specifying the current gene identifier type.
#'   One of \code{"SYMBOL"}, \code{"ENSEMBL"}, \code{"ENTREZID"}, or
#'   \code{"REFSEQ"}. Default: \code{"ENSEMBL"}.
#' @param to Character string specifying the target gene identifier type.
#'   One of \code{"SYMBOL"}, \code{"ENSEMBL"}, \code{"ENTREZID"}, or
#'   \code{"REFSEQ"}. Default: \code{"SYMBOL"}.
#' @param agg_fun Character string specifying how to aggregate expression
#'   values when multiple source IDs map to the same target ID:
#'   \code{"max"} (default), \code{"mean"}, or \code{"median"}.
#' @param drop_unmapped Logical. If \code{TRUE} (default), genes that cannot
#'   be mapped are removed from the output. If \code{FALSE}, unmapped genes
#'   are retained with their original row names.
#' @param show_stats Logical. If \code{TRUE} (default), prints a mapping
#'   summary to the console showing the number of genes mapped, removed,
#'   and the mapping rate.
#'
#' @return A numeric matrix with the same columns as \code{expr} and row
#'   names converted to the target identifier type.
#'
#' @details
#' **Supported identifier types:**
#' \describe{
#'   \item{\code{"SYMBOL"}}{Official HGNC gene symbols (e.g. TP53, EGFR).}
#'   \item{\code{"ENSEMBL"}}{Ensembl gene IDs (e.g. ENSG00000141510).}
#'   \item{\code{"ENTREZID"}}{NCBI Gene IDs (e.g. 7157).}
#'   \item{\code{"REFSEQ"}}{RefSeq mRNA accessions (e.g. NM_000546).}
#' }
#'
#' **ENSEMBL version suffix handling:**
#' Row names with an Ensembl version suffix (e.g. \code{"ENSG00000141510.11"})
#' are automatically stripped before lookup.
#'
#' **Aggregation strategy:**
#' When multiple source IDs map to the same target identifier (common when
#' converting from transcript-level to gene-level), one of three strategies
#' is used:
#' \itemize{
#'   \item \code{"max"}: keep the row with the highest expression value per
#'     sample (default, conservative for expression).
#'   \item \code{"mean"}: average expression across all source IDs per sample.
#'   \item \code{"median"}: median expression across all source IDs per sample.
#' }
#'
#' @section Dependencies:
#' This function requires the Bioconductor packages \pkg{org.Hs.eg.db} and
#' \pkg{AnnotationDbi}. Install them with:
#' \preformatted{
#' BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi"))
#' }
#'
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#' @importFrom stats aggregate median
#' @export
#'
#' @examples
#' \dontrun{
#' # ENSEMBL to SYMBOL (most common use case)
#' data(TCGA_LUSC_ENSEMBL)
#' expr_symbol <- convert_id(TCGA_LUSC_ENSEMBL, from = "ENSEMBL", to = "SYMBOL")
#' head(rownames(expr_symbol))
#'
#' # SYMBOL to ENTREZID, using mean for aggregation
#' data(TCGA_LUSC)
#' expr_entrez <- convert_id(TCGA_LUSC,
#'                            from = "SYMBOL",
#'                            to = "ENTREZID",
#'                            agg_fun = "mean")
#'
#' # Keep unmapped genes
#' expr_keep <- convert_id(TCGA_LUSC_ENSEMBL,
#'                          from = "ENSEMBL",
#'                          to = "SYMBOL",
#'                          drop_unmapped = FALSE)
#' }
convert_id <- function(expr,
                        from          = "ENSEMBL",
                        to            = "SYMBOL",
                        agg_fun       = c("max", "mean", "median"),
                        drop_unmapped = TRUE,
                        show_stats    = TRUE) {

  # ---- 1. Input validation ----
  valid_types <- c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ")
  from <- match.arg(from, valid_types)
  to   <- match.arg(to,   valid_types)
  agg_fun <- match.arg(agg_fun)

  if (from == to) {
    message("'from' and 'to' are the same; returning input unchanged.")
    return(as.matrix(expr))
  }

  if (is.null(rownames(expr))) {
    stop("'expr' must have row names.", call. = FALSE)
  }

  # ---- 2. Validate that declared idType matches actual row names ----
  suppressMessages(validate_id_type(expr, from))

  # ---- 3. Prepare source IDs (strip ENSEMBL version suffix) ----
  src_ids <- rownames(expr)
  src_ids_clean <- sub("\\.\\d+$", "", src_ids)
  n_input <- length(src_ids)

  # ---- 4. Look up mapping via org.Hs.eg.db ----
  mapping <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys      = src_ids_clean,
      keytype   = from,
      columns   = to
    )
  )

  # ---- 5. Align mapping to expression rows ----
  idx <- match(src_ids_clean, mapping[[from]])
  target_vec <- mapping[[to]][idx]
  keep <- !is.na(target_vec)

  n_mapped <- sum(keep)
  n_dropped <- n_input - n_mapped

  if (n_mapped == 0) {
    stop("No identifiers could be mapped from '", from, "' to '", to, "'.",
         call. = FALSE)
  }

  # ---- 6. Handle unmapped genes ----
  if (drop_unmapped) {
    expr_work <- expr[keep, , drop = FALSE]
    target_work <- target_vec[keep]
  } else {
    expr_work <- expr
    target_work <- ifelse(keep, target_vec, src_ids)
  }

  # ---- 7. Aggregate duplicates ----
  expr_df <- as.data.frame(expr_work)
  expr_df$TARGET_ID <- target_work

  agg_func <- switch(agg_fun,
    max    = function(x) {
      val <- suppressWarnings(max(x, na.rm = TRUE))
      if (is.infinite(val)) NA_real_ else val
    },
    mean   = function(x) {
      val <- mean(x, na.rm = TRUE)
      if (is.nan(val)) NA_real_ else val
    },
    median = function(x) {
      val <- stats::median(x, na.rm = TRUE)
      if (is.nan(val)) NA_real_ else val
    }
  )

  groups_list <- split(expr_df, expr_df$TARGET_ID)

  result_list <- lapply(names(groups_list), function(gname) {
    sub <- groups_list[[gname]]
    sub$TARGET_ID <- NULL
    vals <- vapply(sub, function(col) agg_func(col), numeric(1))
    vals
  })

  result_mat <- do.call(rbind, result_list)
  rownames(result_mat) <- names(groups_list)
  colnames(result_mat) <- colnames(expr)

  # ---- 8. Print mapping summary ----
  if (show_stats) {
    n_multi <- sum(duplicated(target_work) | duplicated(target_work, fromLast = TRUE))
    pct <- round(n_mapped / n_input * 100, 1)

    cat(
      "\n",
      "===== Gene ID conversion summary =====\n",
      "  Input:     ", n_input,           " genes (", from,       ")\n",
      "  Mapped:    ", n_mapped,          " genes (", to,         ")  [", pct, "%]\n",
      "  Unmapped:  ", n_dropped,         " genes\n",
      "  Duplicates:", n_multi,           " source IDs with multi-mapping\n",
      "  Aggregated by:", agg_fun,        "\n",
      "=====================================\n",
      sep = ""
    )
  }

  return(result_mat)
}

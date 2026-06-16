#' Plot a subtype‑specific gene expression heatmap
#'
#' Draws a heatmap of curated subtype‑driver genes using the output of
#' \code{\link{classifyHNSC}}.  The function scales expression values per gene
#' and annotates samples by their predicted molecular subtype.
#'
#' @param expr A numeric gene expression matrix with genes in rows and samples
#'   in columns.  Row names must be gene symbols (default) or identifiers
#'   matching \code{idType}.  Column names must be sample IDs.
#' @param subtype A **named factor vector** giving the predicted subtype for
#'   each sample (e.g. \code{c(sample1 = "Basal", sample2 = "Classical")}).
#'   This is the direct output of \code{classifyHNSC(..., outputType = "class")}.
#'   A two‑column data frame with columns \code{Sample} and \code{Subtype}
#'   is also accepted for convenience, but will be coerced to a named vector.
#'   **Probability matrices are not allowed**.
#' @param idType Character string specifying the gene identifier type used in
#'   the rownames of \code{expr}.  Acceptable values are \code{"SYMBOL"}
#'   (default), \code{"ENSEMBL"}, \code{"ENTREZID"}, or \code{"REFSEQ"}.
#'   Non‑symbol identifiers are automatically converted to gene symbols.
#' @param gene_set A named list of genes to display, where names correspond to
#'   the subtype levels.  The default uses a built‑in curated set of known
#'   HNSCC subtype markers.
#' @param main Title for the heatmap.
#' @param show_rownames Logical, whether to show gene names (default \code{TRUE}).
#' @param cluster_rows Logical, whether to cluster rows (default \code{FALSE}).
#' @param cluster_cols Logical, whether to cluster columns (default \code{FALSE}).
#' @param cap Numeric value at which to cap (trim) the scaled expression values.
#'   Values above \code{cap} are set to \code{cap}, values below \code{-cap}
#'   are set to \code{-cap}.  Default is \code{NULL} (no capping).  Typical
#'   choices are \code{2} or \code{3}, but the appropriate threshold depends
#'   on the data distribution.
#' @param color A character vector of colours to use for the heatmap gradient.
#'   If \code{NULL} (default), a three‑colour gradient is built from
#'   \code{low_color}, \code{mid_color} and \code{high_color} with
#'   \code{n_color} levels.  If supplied, this vector is passed directly to
#'   \code{\link[pheatmap]{pheatmap}} and the \code{low_}/\code{mid_}/\code{high_}
#'   arguments are ignored.
#' @param low_color,mid_color,high_color Colours for the gradient when
#'   \code{color = NULL}.  \code{low_color} and \code{high_color} map to the
#'   minimum and maximum scaled values, while \code{mid_color} corresponds to
#'   zero.  Defaults are \code{"#2c7bb6"} (blue), \code{"white"}, and
#'   \code{"#d7191c"} (red), respectively.
#' @param n_color Number of colour levels to generate in the gradient
#'   (default 100).  Only used when \code{color = NULL}.
#' @param silent Logical.  If \code{FALSE} (default), the heatmap is drawn
#'   directly to the current graphics device.  If \code{TRUE}, the heatmap
#'   is built but not drawn, returning a gtable object for manual drawing
#'   with \code{grid::grid.draw()}.
#' @param ... Additional arguments passed to \code{\link[pheatmap]{pheatmap}}.
#'   Note that \code{color} is set internally and should not be passed via
#'   \code{...}.
#'
#' @return A \code{pheatmap} object (invisibly).  The heatmap is drawn
#'   as a side effect.
#'
#' @note The \pkg{pheatmap} package must be installed.
#'
#' @importFrom stats setNames
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @export
#'
#' @examples
#' \dontrun{
#' # Directly from classifyHNSC output
#' subtypes <- classifyHNSC(TCGA_LUSC, outputType = "class")
#' plot_subtype_heatmap(TCGA_LUSC, subtypes)
#'
#' # With custom capping and colours
#' plot_subtype_heatmap(TCGA_LUSC, subtypes,
#'                      cap = 2,
#'                      low_color = "green", high_color = "red")
#' }

plot_subtype_heatmap <- function(expr,
                                 subtype,
                                 idType = c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"),
                                 gene_set = NULL,
                                 main = "Subtype heatmap",
                                 show_rownames = TRUE,
                                 cluster_rows = FALSE,
                                 cluster_cols = FALSE,
                                 cap = NULL,
                                 color = NULL,
                                 low_color = "#2c7bb6",
                                 mid_color = "white",
                                 high_color = "#d7191c",
                                 n_color = 100,
                                 silent = FALSE,
                                 ...) {
  # ---------- 1. Check pheatmap availability ----------
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required. Install with: install.packages('pheatmap')")
  }

  idType <- match.arg(idType)

  # ---- Validate that the declared idType matches actual row names ----
  validate_id_type(expr, idType)


  # ---------- 2. Coerce subtype to a named character vector ----------
  if (is.atomic(subtype) && !is.null(names(subtype))) {
    # 接受命名向量（包括 factor），统一转为字符
    subtype_vec <- as.character(subtype)
    names(subtype_vec) <- names(subtype)
  } else if (is.data.frame(subtype)) {
    if (!all(c("Sample", "Subtype") %in% colnames(subtype))) {
      stop("When 'subtype' is a data frame, it must contain columns 'Sample' and 'Subtype'.")
    }
    subtype_vec <- setNames(as.character(subtype$Subtype),
                            as.character(subtype$Sample))
  } else if (is.matrix(subtype) || (is.data.frame(subtype) && ncol(subtype) > 2)) {
    stop("This function requires class labels (outputType='class'), not a probability matrix.")
  } else {
    stop("'subtype' must be a named character/factor vector or a data frame with 'Sample' and 'Subtype' columns.")
  }

  # 构建注释数据框
  anno <- data.frame(Subtype = subtype_vec,
                     row.names = names(subtype_vec))

  # ---------- 3. ID conversion (if needed) ----------
  if (idType != "SYMBOL") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop(
        "Package 'org.Hs.eg.db' is required. Install with BiocManager::install('org.Hs.eg.db')"
      )
    if (!requireNamespace("AnnotationDbi", quietly = TRUE))
      stop(
        "Package 'AnnotationDbi' is required. Install with BiocManager::install('AnnotationDbi')"
      )
    mapping <- suppressMessages(
      AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = rownames(expr),
        keytype = idType,
        columns = "SYMBOL"
      )
    )
    map_vec <- mapping$SYMBOL[match(rownames(expr), mapping[[idType]])]
    keep <- !is.na(map_vec)
    expr <- expr[keep, , drop = FALSE]
    rownames(expr) <- map_vec[keep]
    message("Converted ", sum(keep), " rows from ", idType, " to SYMBOL.")
  }

  # ---------- 4. Prepare gene set ----------
  if (is.null(gene_set)) {
    gene_set <- subtype_genes   # internal package data
  }
  all_genes <- unique(unlist(gene_set))
  found_genes <- intersect(all_genes, rownames(expr))
  if (length(found_genes) == 0) {
    stop("None of the genes in 'gene_set' were found in the expression matrix.")
  }

  # ---------- 5. Align samples and order by subtype ----------
  common_samples <- intersect(colnames(expr), names(subtype_vec))
  if (length(common_samples) == 0) {
    stop("No common samples between expression matrix and subtype annotations.")
  }
  subtype_vec <- subtype_vec[common_samples]
  sample_order <- names(sort(subtype_vec))
  expr_sub <- expr[found_genes, sample_order, drop = FALSE]
  anno <- anno[sample_order, , drop = FALSE]

  # ---------- 6. Scale and optionally cap ----------
  tmp <- t(scale(t(expr_sub)))
  if (!is.null(cap)) {
    if (!is.numeric(cap) || length(cap) != 1 || cap <= 0) {
      stop("'cap' must be a single positive number or NULL.")
    }
    tmp[tmp > cap] <- cap
    tmp[tmp < -cap] <- -cap
  }

  # ---------- 7. Build color palette ----------
  if (is.null(color)) {
    color <- colorRampPalette(c(low_color, mid_color, high_color))(n_color)
  } else {
    if (!is.character(color))
      stop("'color' must be a character vector of colours.")
  }

  # ---------- 8. Draw heatmap ----------
  pheatmap::pheatmap(
    tmp,
    annotation_col = anno,
    show_colnames   = FALSE,
    show_rownames   = show_rownames,
    cluster_rows    = cluster_rows,
    cluster_cols    = cluster_cols,
    main            = main,
    color           = color,
    silent          = silent,
    ...
  )
}

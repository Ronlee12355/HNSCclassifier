#' Summarize Molecular Subtype Classification Results
#'
#' Produces a tidy summary table from the output of
#' \code{\link{classifyHNSC}(outputType = "class")}, showing sample counts,
#' percentages, and an optional bar plot of the subtype distribution.
#'
#' @param subtype A named character vector or factor giving the predicted
#'   subtype for each sample, as returned by
#'   \code{classifyHNSC(outputType = "class")}.  Must have names (sample IDs).
#'   Alternatively, a data frame with columns \code{Sample} and \code{Subtype}.
#' @param plot Logical. If \code{TRUE} (default), a \code{ggplot2} bar chart
#'   of the subtype distribution is printed as a side effect.
#'
#' @return A data frame (invisibly, if \code{plot = TRUE}) with columns:
#'   \itemize{
#'     \item \code{Subtype}: the four TCGA molecular subtypes
#'     \item \code{N}: number of samples assigned to each subtype
#'     \item \code{Pct}: percentage of total samples (rounded to 1 decimal)
#'     \item \code{Proportion}: numeric proportion (0–1)
#'   }
#'   The rows are sorted by frequency in descending order. If a subtype
#'   has zero samples, it is still included with \code{N = 0}.
#'
#' @details
#' The function always includes all four TCGA subtypes (Atypical, Basal,
#' Classical, Mesenchymal) in the output, even if one or more have zero
#' assigned samples. This ensures consistent output across cohorts.
#'
#' @importFrom ggplot2 ggplot aes geom_col geom_text labs theme_minimal
#' @importFrom ggplot2 theme element_text element_blank scale_fill_manual
#' @importFrom ggplot2 coord_flip .data
#' @importFrom stats setNames reorder
#' @export
#'
#' @examples
#' \dontrun{
#' data(TCGA_LUSC)
#' subtypes <- classifyHNSC(TCGA_LUSC, outputType = "class")
#'
#' # Table + bar plot (default)
#' summarize_subtype(subtypes)
#'
#' # Table only, store result
#' tbl <- summarize_subtype(subtypes, plot = FALSE)
#' print(tbl)
#' }
summarize_subtype <- function(subtype,
                              plot = TRUE) {

  # ---- Check ggplot2 availability ----
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }

  # ---- 1. Coerce subtype to named character vector ----
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

  # ---- 2. Count by subtype ----
  all_subtypes <- c("Atypical", "Basal", "Classical", "Mesenchymal")
  subtype_vec <- factor(subtype_vec, levels = all_subtypes)

  counts <- table(subtype_vec)
  n_total <- length(subtype_vec)

  # ---- 3. Build summary table ----
  tbl <- data.frame(
    Subtype    = all_subtypes,
    N          = as.integer(counts),
    Pct        = round(as.numeric(counts) / n_total * 100, 1),
    Proportion = round(as.numeric(counts) / n_total, 3),
    stringsAsFactors = FALSE
  )
  tbl <- tbl[order(tbl$N, decreasing = TRUE), , drop = FALSE]
  rownames(tbl) <- NULL

  # ---- 4. Print summary to console ----
  cat("\n")
  cat("===== HNSC Molecular Subtype Classification Summary =====\n")
  cat("Total samples:", n_total, "\n\n")
  print(tbl, row.names = FALSE, right = FALSE)
  cat("\nPredominant subtype:", tbl$Subtype[1], "\n")
  cat("=========================================================\n\n")

  # ---- 5. Optional bar plot ----
  if (plot) {
    subtype_colors <- c(
      Atypical    = "#87C304",
      Basal       = "#FC7C75",
      Classical   = "#1AD3D8",
      Mesenchymal = "#D880FF"
    )

    p <- ggplot2::ggplot(
      tbl,
      ggplot2::aes(
        x = stats::reorder(.data$Subtype, .data$N),
        y = .data$N,
        fill = .data$Subtype
      )
    ) +
      ggplot2::geom_col(width = 0.65, colour = NA) +
      ggplot2::geom_text(
        ggplot2::aes(label = paste0(.data$N, " (", .data$Pct, "%)")),
        hjust = -0.1,
        size = 4
      ) +
      ggplot2::scale_fill_manual(values = subtype_colors, guide = "none") +
      ggplot2::coord_flip() +
      ggplot2::labs(
        x = NULL,
        y = "Number of Samples",
        title = "Molecular Subtype Distribution"
      ) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )

    print(p)
    return(invisible(tbl))
  }

  return(tbl)
}

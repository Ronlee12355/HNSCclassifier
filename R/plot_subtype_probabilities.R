#' Plot Posterior Subtype Probabilities for Each Sample
#'
#' Visualises the posterior probability of each TCGA molecular subtype
#' (Atypical, Basal, Classical, Mesenchymal) across samples as a stacked
#' bar chart.  Samples are ordered by predicted subtype, then by
#' confidence within each subtype group.  Low-confidence samples — where
#' the maximum probability falls below a user-specified threshold — are
#' optionally highlighted with a visual marker.
#'
#' This plot is useful for quickly assessing classification certainty
#' across a cohort and for identifying borderline cases that may warrant
#' further investigation.
#'
#' @param probs A numeric probability matrix returned by
#'   \code{\link{classifyHNSC}(outputType = "prob")}, with samples in
#'   rows and four subtype columns.
#' @param threshold Numeric value between 0 and 1. Samples whose maximum
#'   posterior probability is below this threshold are flagged as
#'   low-confidence (default: 0.5).
#' @param mark_low Logical. If \code{TRUE} (default), a red asterisk is
#'   placed above bars corresponding to low-confidence samples.
#' @param palette A named character vector of four colours for the
#'   subtypes.  Names must match the subtype column names.  If
#'   \code{NULL}, a default colour palette is used.
#'
#' @return A \code{ggplot2} object (stacked bar chart).  The plot can be
#'   further customised with standard \code{ggplot2} layers.
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_col scale_fill_manual theme
#' @importFrom ggplot2 element_text element_blank labs coord_flip
#' @importFrom ggplot2 annotate scale_y_continuous expansion .data
#' @importFrom stats reorder
#'
#' @examples
#' \dontrun{
#' data(TCGA_LUSC)
#' probs <- classifyHNSC(TCGA_LUSC, outputType = "prob")
#' plot_subtype_probabilities(probs)
#' plot_subtype_probabilities(probs, threshold = 0.5)
#' }
plot_subtype_probabilities <- function(probs,
                                       threshold = 0.5,
                                       mark_low  = TRUE,
                                       palette   = NULL) {

  # ---- validate input -------------------------------------------------------
  required_subtypes <- c("Atypical", "Basal", "Classical", "Mesenchymal")
  missing_cols <- setdiff(required_subtypes, colnames(probs))
  if (length(missing_cols) > 0) {
    stop("'probs' must contain columns: ",
         paste(required_subtypes, collapse = ", "), ".")
  }
  if (!is.numeric(probs)) {
    stop("'probs' must be a numeric matrix.")
  }
  if (any(probs < 0 | probs > 1, na.rm = TRUE)) {
    stop("All values in 'probs' must be between 0 and 1.")
  }
  if (threshold <= 0 || threshold >= 1) {
    stop("'threshold' must be between 0 and 1 (exclusive).")
  }

  # ---- prepare data ---------------------------------------------------------
  probs <- as.data.frame(probs)
  probs$Sample <- rownames(probs)
  probs$Predicted <- apply(probs[, required_subtypes], 1, which.max)
  probs$Predicted <- factor(required_subtypes[probs$Predicted],
                            levels = required_subtypes)
  probs$MaxProb <- apply(probs[, required_subtypes], 1, max)
  probs$LowConf <- probs$MaxProb < threshold

  # Order: by Predicted subtype, then descending MaxProb within each group
  probs <- probs[order(probs$Predicted, -probs$MaxProb), ]
  probs$Sample <- factor(probs$Sample, levels = probs$Sample)

  # Reshape to long format for ggplot
  long <- data.frame(
    Sample   = rep(probs$Sample, 4),
    Subtype  = factor(rep(required_subtypes, each = nrow(probs)),
                      levels = required_subtypes),
    Prob     = c(probs$Atypical, probs$Basal,
                 probs$Classical, probs$Mesenchymal),
    LowConf  = rep(probs$LowConf, 4)
  )

  # ---- default palette ------------------------------------------------------
  if (is.null(palette)) {
    palette <- c(
      Atypical    = "#87C304",
      Basal       = "#FC7C75",
      Classical   = "#1AD3D8",
      Mesenchymal = "#D880FF"
    )
  }

  # ---- build plot -----------------------------------------------------------
  p <- ggplot(long, aes(x = .data$Sample,
                        y = .data$Prob,
                        fill = .data$Subtype)) +
    geom_col(width = 1,
             colour = NA) +
    scale_fill_manual(values = palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Posterior Probability",
         fill = "Subtype") +
    theme(
      axis.text.x      = element_text(angle = 90, hjust = 1,
                                      vjust = 0.5, size = 6),
      axis.ticks.x     = element_blank(),
      panel.background = element_blank(),
      legend.position  = "right"
    )

  # ---- mark low-confidence samples ------------------------------------------
  if (mark_low) {
    low_samples <- probs$Sample[probs$LowConf]
    if (length(low_samples) > 0) {
      p <- p + annotate(
        "text",
        x = low_samples,
        y = 1.05,
        label = "*",
        colour = "red",
        size   = 3,
        fontface = "bold"
      )
    }
  }

  return(p)
}

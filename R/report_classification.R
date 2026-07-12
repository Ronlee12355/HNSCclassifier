#' Generate an HTML Classification Report
#'
#' Compiles classification results, visualizations, and optional differential
#' expression analysis into a self-contained HTML report using an internal
#' R Markdown template.
#'
#' @param expr A numeric gene expression matrix with genes in rows and
#'   samples in columns. Required for the heatmap section.
#' @param subtype A named character vector of predicted subtypes, as returned
#'   by \code{\link{classifyHNSC}(outputType = "class")}. A data frame with
#'   columns \code{Sample} and \code{Subtype} is also accepted.
#' @param probs Optional. A probability matrix as returned by
#'   \code{\link{classifyHNSC}(outputType = "prob")}. If provided, a stacked
#'   bar chart of posterior probabilities is included in the report.
#' @param degs Optional. A named list of data frames as returned by
#'   \code{\link{extract_top_genes}}. If provided, the top differentially
#'   expressed genes per subtype are included.
#' @param idType Character string specifying the gene identifier type.
#'   One of \code{"SYMBOL"} (default), \code{"ENSEMBL"}, \code{"ENTREZID"},
#'   or \code{"REFSEQ"}. Only used when \code{expr} is provided for the
#'   heatmap section.
#' @param output_file Character string for the output HTML file path.
#'   Default: \code{"HNSC_classification_report.html"} in the current
#'   working directory.
#'
#' @return The absolute path to the generated HTML file (invisibly).
#'
#' @details
#' The report is rendered from an internal R Markdown template and includes
#' the following sections:
#' \enumerate{
#'   \item **Classification Summary** — a table of subtype distribution
#'     with counts and percentages.
#'   \item **Subtype Distribution** — a bar chart visualisation of the
#'     summary table.
#'   \item **Posterior Probabilities** — a stacked bar chart of per-sample
#'     subtype probabilities (only shown if \code{probs} is provided).
#'   \item **Subtype Marker Gene Heatmap** — an annotated heatmap of
#'     curated subtype-driver genes (only shown if \code{expr} is provided).
#'   \item **Differentially Expressed Genes** — top DEGs per subtype from
#'     one-vs-rest analysis (only shown if \code{degs} is provided).
#'   \item **Sample-Level Classification** — a complete table of per-sample
#'     assignments.
#' }
#'
#' The HTML file is self-contained (all plots are embedded as base64) and
#' can be shared with collaborators or included as supplementary material.
#'
#' @note Requires the \pkg{rmarkdown} and \pkg{knitr} packages.
#'   Install with: \code{install.packages(c("rmarkdown", "knitr"))}
#'
#' @seealso
#' \code{\link{classifyHNSC}} for generating the required inputs,
#' \code{\link{extract_top_genes}} for differential expression analysis.
#'
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' \dontrun{
#' data(TCGA_LUSC)
#' subtypes <- classifyHNSC(TCGA_LUSC, outputType = "class")
#' probs <- classifyHNSC(TCGA_LUSC, outputType = "prob")
#'
#' # Minimal report (class labels only)
#' report_classification(subtype = subtypes)
#'
#' # Full report with heatmap and probabilities
#' report_classification(
#'   expr    = TCGA_LUSC,
#'   subtype = subtypes,
#'   probs   = probs
#' )
#'
#' # Include differential expression results
#' degs <- extract_top_genes(TCGA_LUSC, subtypes)
#' report_classification(
#'   expr    = TCGA_LUSC,
#'   subtype = subtypes,
#'   probs   = probs,
#'   degs    = degs
#' )
#' }
report_classification <- function(expr     = NULL,
                                   subtype,
                                   probs    = NULL,
                                   degs     = NULL,
                                   idType   = "SYMBOL",
                                   output_file = "HNSC_classification_report.html") {

  # ---- 1. Check rmarkdown availability ----
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package 'rmarkdown' is required. Install with: install.packages('rmarkdown')")
  }

  # ---- 2. Validate and coerce subtype ----
  if (is.data.frame(subtype)) {
    if (!all(c("Sample", "Subtype") %in% colnames(subtype))) {
      stop("When 'subtype' is a data frame, it must have columns 'Sample' and 'Subtype'.")
    }
    subtype_vec <- stats::setNames(as.character(subtype$Subtype),
                                   as.character(subtype$Sample))
  } else if (is.atomic(subtype) && !is.null(names(subtype))) {
    subtype_vec <- as.character(subtype)
    names(subtype_vec) <- names(subtype)
  } else {
    stop("'subtype' must be a named character/factor vector or a data frame.")
  }

  # ---- 3. Validate probs (if provided) ----
  if (!is.null(probs)) {
    required_subtypes <- c("Atypical", "Basal", "Classical", "Mesenchymal")
    missing_cols <- setdiff(required_subtypes, colnames(probs))
    if (length(missing_cols) > 0) {
      stop("'probs' must contain columns: ",
           paste(required_subtypes, collapse = ", "), ".")
    }
    if (any(probs < 0 | probs > 1, na.rm = TRUE)) {
      stop("All values in 'probs' must be between 0 and 1.")
    }
  }

  # ---- 4. Build summary table ----
  all_subtypes <- c("Atypical", "Basal", "Classical", "Mesenchymal")
  subtype_fac <- factor(subtype_vec, levels = all_subtypes)
  counts <- table(subtype_fac)
  n_total <- length(subtype_vec)

  summary_table <- data.frame(
    Subtype    = all_subtypes,
    N          = as.integer(counts),
    Pct        = round(as.numeric(counts) / n_total * 100, 1),
    Proportion = round(as.numeric(counts) / n_total, 3),
    stringsAsFactors = FALSE
  )

  # ---- 5. Build sample table ----
  sample_table <- data.frame(
    Sample  = names(subtype_vec),
    Subtype = unname(subtype_vec),
    stringsAsFactors = FALSE
  )

  # ---- 6. Locate internal template ----
  template <- system.file("report_template.Rmd", package = "HNSCclassifier")
  if (template == "") {
    stop("Internal report template not found. Re-install the package.")
  }

  # ---- 7. Render ----
  output_dir <- dirname(output_file)
  if (output_dir == ".") {
    output_dir <- getwd()
  } else {
    output_dir <- normalizePath(output_dir, mustWork = FALSE)
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  rmarkdown::render(
    input       = template,
    output_file = basename(output_file),
    output_dir  = output_dir,
    params      = list(
      n_total       = n_total,
      summary_table = summary_table,
      sample_table  = sample_table,
      probs         = probs,
      degs          = degs,
      expr          = expr,
      subtype       = subtype_vec,
      idType        = idType,
      has_probs     = !is.null(probs),
      has_heatmap   = !is.null(expr),
      has_degs      = !is.null(degs)
    ),
    envir     = new.env(parent = globalenv()),
    quiet     = TRUE
  )

  abs_path <- normalizePath(file.path(output_dir, basename(output_file)))
  message("Report generated: ", abs_path)
  invisible(abs_path)
}

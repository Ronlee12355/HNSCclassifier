#' Check for missing values in an expression matrix
#'
#' @param x A numeric matrix or data frame.
#'
#' @return `TRUE` if any missing values are present, `FALSE` otherwise.
#' @keywords internal
check_missing_values <- function(x) {
  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame.")
  }
  res <- any(is.na(x))
  return(res)
}

#' Check if an expression matrix appears to be log2-transformed
#'
#' This internal function uses simple heuristics to guess whether the input
#' matrix has already been log2-transformed.  It returns `TRUE` if the data
#' appear to be on a log2 scale, and `FALSE` otherwise.
#' @param x A numeric matrix or data frame.
#' @return Logical: `TRUE` if the matrix is already log2-transformed, `FALSE` otherwise.
#' @keywords internal
#' @importFrom stats quantile
Log2ed <- function(x) {
  qx <- as.numeric(quantile(x, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
  LogC <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    return(FALSE)
  } else{
    return(TRUE)
  }
}

#' Detect gene identifier type from row names
#'
#' Examines a sample of row names from an expression matrix and infers
#' the gene identifier system (SYMBOL, ENSEMBL, ENTREZID, or REFSEQ)
#' using regular expression patterns.
#'
#' @param rnames Character vector of row names (gene identifiers).
#'
#' @return A character string: one of \code{"SYMBOL"}, \code{"ENSEMBL"},
#'   \code{"ENTREZID"}, or \code{"REFSEQ"}.
#'
#' @keywords internal

detect_id_type <- function(rnames) {

  # Sample up to 100 row names for efficiency on large matrices.
  # Row names are typically unique enough that a moderate sample
  # is sufficient to determine the identifier type.
  n_sample <- min(length(rnames), 100L)
  ids <- as.character(rnames[seq_len(n_sample)])

  # ---- Define detection patterns ----------------------------------------------
  #
  # ENSEMBL:  Human Ensembl gene IDs are "ENSG" followed by exactly 11 digits.
  #           Versioned forms like "ENSG00000141510.11" are also common;
  #           we strip the suffix before matching.
  #
  # ENTREZID: NCBI Gene IDs are pure integers (e.g. 7157, 1956).  No letters,
  #           no punctuation, no leading zeros — just digits.
  #
  # REFSEQ:   RefSeq mRNA accessions start with NM_ / NR_ / XM_ / XR_,
  #           followed by digits.  They always contain an underscore.
  #           e.g. NM_000546, XM_011545678.
  #
  # SYMBOL:   Official HGNC symbols (e.g. TP53, EGFR, BRCA1, HLA-DRB1).
  #           This is treated as the fallback: if none of the above
  #           patterns dominate, we assume the IDs are gene symbols.
  #           Symbols are the hardest to validate with a single regex
  #           because they can contain letters, digits, and hyphens.

  ensembl_stripped <- sub("\\.\\d+$", "", ids)  # drop version suffix

  is_ensembl  <- grepl("^ENSG\\d{11}$", ensembl_stripped, ignore.case = TRUE)
  is_entrez   <- grepl("^\\d+$", ids)
  is_refseq   <- grepl("^[NX][MR]_\\d+", ids, ignore.case = TRUE)

  # ---- Compute match rates ----------------------------------------------------
  rate_ensembl  <- mean(is_ensembl)
  rate_entrez   <- mean(is_entrez)
  rate_refseq   <- mean(is_refseq)

  # ---- Decision logic ---------------------------------------------------------
  # A match rate ≥ 0.9 for a given pattern is taken as strong evidence.
  # If more than one pattern exceeds the threshold (unlikely in practice),
  # the highest rate wins.
  # If no pattern reaches 0.9, we default to SYMBOL — either the data truly
  # contain gene symbols, or the identifiers are too ambiguous to classify.

  rates <- c(ENSEMBL = rate_ensembl, ENTREZID = rate_entrez,
             REFSEQ = rate_refseq)
  best  <- which.max(rates)

  if (rates[best] >= 0.9) {
    detected <- names(rates)[best]
  } else {
    detected <- "SYMBOL"
  }

  return(detected)
}

#' Validate that the declared gene ID type matches the actual row names
#'
#' Compares the identifier type declared by the user (\code{idType})
#' against the type detected from the expression matrix row names.  If
#' a mismatch is found, the function stops with an informative error
#' message suggesting the correct \code{idType} value.  If the types
#' agree, the function returns invisibly (no output).
#'
#' @param input_expr A numeric gene expression matrix with genes in rows.
#' @param idType Character string.  The identifier type declared by the
#'   user: \code{"SYMBOL"}, \code{"ENSEMBL"}, \code{"ENTREZID"}, or
#'   \code{"REFSEQ"}.
#'
#' @return Invisibly returns \code{NULL}.  Called for its side effect
#'   (stop on mismatch, nothing on match).
#'
#' @keywords internal
#' @importFrom utils head
validate_id_type <- function(input_expr, idType) {

  # Only the four supported types are accepted at this point.
  # (The caller is responsible for match.arg() before calling this.)
  valid_types <- c("SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ")

  if (!idType %in% valid_types) {
    stop("Unsupported idType '", idType, "'. ",
         "Must be one of: ", paste(valid_types, collapse = ", "), ".")
  }

  # Detect the actual type from row names
  detected <- detect_id_type(rownames(input_expr))

  # SYMBOL detection is a fallback when no regex matches strongly.
  # Therefore, when the user says SYMBOL, we accept it without
  # complaint — we cannot prove them wrong with regex alone.
  if (idType == "SYMBOL") {
    return(invisible(NULL))
  }

  # If user declared ENSEMBL / ENTREZID / REFSEQ but the row names
  # look like something else, stop with a helpful message.
  if (detected != idType) {

    # Build an example that shows the user what we saw
    example_ids <- head(rownames(input_expr), 5L)

    msg <- paste0(
      "Gene ID mismatch detected.\n",
      "  You specified idType = \"", idType, "\", ",
      "but the row names appear to be \"", detected, "\".\n",
      "  Example row names: ", paste(example_ids, collapse = ", "), "\n\n",
      "Please either:\n",
      "  (1) set idType = \"", detected, "\" in your classifyHNSC() call, or\n",
      "  (2) convert your row names to \"", idType,
      "\" before submitting.\n\n",
      "Supported idType values: ",
      paste(valid_types, collapse = ", ")
    )
    stop(msg)
  }

  invisible(NULL)
}

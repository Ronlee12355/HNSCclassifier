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

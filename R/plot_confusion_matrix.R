#' Plot a Confusion Matrix from caret's confusionMatrix Object
#'
#' Creates an elegant heatmap visualization of a confusion matrix using ggplot2.
#' The plot displays predicted vs actual class counts with color-coded cells
#' and numeric labels for easy interpretation.
#'
#' @param cm An object of class "confusionMatrix" created by
#'   \code{\link[caret]{confusionMatrix}} from the \pkg{caret} package.
#' @param title A character string for the plot title (default: "Confusion Matrix").
#' @param xlab Label for the x-axis representing true/reference classes
#'   (default: "Reference").
#' @param ylab Label for the y-axis representing predicted classes
#'   (default: "Prediction").
#' @param low_color Colour gradient start for cells with low counts
#'   (default: "white").
#' @param high_color Colour gradient end for cells with high counts
#'   (default: "steelblue").
#' @param text_color Colour of the numeric count labels inside cells
#'   (default: "black").
#' @param text_size Size of the cell count labels in points (default: 5).
#'
#' @return A \code{ggplot} object that can be further customized or displayed.
#'
#' @note The factor levels are automatically set to match the original order
#'   from the confusion matrix table, ensuring consistent axis ordering.
#'
#' @examples
#' \dontrun{
#' # Example 1: Basic usage
#' library(caret)
#'
#' # Create sample predictions and references
#' set.seed(123)
#' actual <- factor(sample(c("A", "B", "C"), 100, replace = TRUE))
#' predicted <- factor(sample(c("A", "B", "C"), 100, replace = TRUE))
#'
#' # Generate confusion matrix
#' cm <- confusionMatrix(predicted, actual)
#'
#' # Plot with default settings
#' plot_confusion_matrix(cm)
#'
#' # Example 2: Customized appearance
#' plot_confusion_matrix(
#'   cm,
#'   title = "Model Performance: Predicted vs Actual",
#'   low_color = "#f7fbff",
#'   high_color = "#08306b",
#'   text_color = "white",
#'   text_size = 4
#' )
#' }
#'
#' @seealso
#'   \code{\link[caret]{confusionMatrix}} for creating confusion matrices
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_text scale_fill_gradient
#'   coord_fixed theme_minimal theme element_blank element_text labs .data
#' @export
plot_confusion_matrix <- function(cm,
                                  title = "Confusion Matrix",
                                  xlab = "Reference",
                                  ylab = "Prediction",
                                  low_color = "white",
                                  high_color = "steelblue",
                                  text_color = "black",
                                  text_size = 5) {

  # Validate input
  if (!inherits(cm, "confusionMatrix")) {
    stop("'cm' must be a confusionMatrix object from the caret package")
  }

  # Extract confusion matrix table and reshape to long format
  tbl <- as.data.frame(cm$table)
  names(tbl) <- c("Prediction", "Reference", "Count")

  # Convert to factors preserving original matrix order
  # Reverse Prediction levels so they display correctly on y-axis
  tbl$Prediction <- factor(tbl$Prediction, levels = rev(rownames(cm$table)))
  tbl$Reference  <- factor(tbl$Reference,  levels = colnames(cm$table))

  # Create the heatmap visualization
  ggplot(tbl, aes(x = .data$Reference, y = .data$Prediction, fill = .data$Count)) +
    geom_tile(color = "grey80", linewidth = 0.5) +
    geom_text(aes(label = .data$Count), color = text_color, size = text_size) +
    scale_fill_gradient(low = low_color, high = high_color, name = "Count") +
    coord_fixed() +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(title = title, x = xlab, y = ylab)
}

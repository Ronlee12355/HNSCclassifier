#' Launch the HNSCclassifier Shiny Web Interface
#'
#' This function starts a local Shiny application that provides an interactive
#' graphical interface for the \code{HNSCclassifier} package. Users can upload
#' a gene expression matrix, set prediction parameters, and obtain TCGA-based
#' molecular subtype classifications without writing any R code.
#'
#'
#' @return This function is called for its side effect of launching the Shiny
#'   application. It does not return a value when the app is stopped.
#'
#' @details The Shiny app is bundled inside the package installation directory
#'   (under \code{shinyApp/}). It relies on the core classification function
#'   \code{\link{classifyHNSC}} and inherits all of its input requirements.
#'   Please make sure that all suggested packages (e.g., \pkg{shiny},
#'  \pkg{DT}, etc.) are installed before launching.
#'
#' @note The first launch may take a few seconds because the app needs to load
#'   the internal pre-trained model and pathway gene sets.
#'
#' @seealso \code{\link{classifyHNSC}} for the command-line version of the
#'   classifier.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Start the Shiny app with default settings
#' classifyHNSC_interface()
#' }
classifyHNSC_interface <-  function(){
  appDir <- system.file('shinyApp', package = 'HNSCclassifier')

  if (appDir == '') {
    stop('Could not load shiny directory. Try re-install HNSCclassifier')
  }

  shiny::runApp(appDir, display.mode = 'normal')
}

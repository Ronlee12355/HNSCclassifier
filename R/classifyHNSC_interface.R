#' The Shiny interface of head and neck squamous cell carcinoma (HNSC) molecular subtypes prediction
#'
#' @description A shiny interface for molecular subtype prediction in HNSCclassifier
#'
#' @export
#' @import shiny
#'
#' @examples
#' classifyHNSC_interface()
classifyHNSC_interface <-  function(){
  appDir <- system.file('shinyApp', package = 'HNSCclassifier')

  if (appDir == '') {
    stop('Could not load shiny directory. Try re-install HNSCclassifier')
  }

  shiny::runApp(appDir, display.mode = 'normal')
}

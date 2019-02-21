#' @export
runExample <- function() {
  appDir <- system.file("shiny-examples", "mcem_ex1.R", package = "mypackage")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}
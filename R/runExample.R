#' @export
runExample <- function() {
  appDir <- system.file("Shiny-examples", "mcem_ex1.R", package = "emphasis")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `emphasis`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}
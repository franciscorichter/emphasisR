#' @export
emphasis_dd <- function() {
  appDir <- system.file("Shiny-examples", "Emphasis_dd2.R", package = "emphasis")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `emphasis`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}
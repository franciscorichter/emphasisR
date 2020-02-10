#' @export
emphasis.app <- function() {
  appDir <- system.file("Shiny-examples", "emphasis_postprocessing.R", package = "emphasis")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `emphasis`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}


#' @export
analytics <- function() {
  appDir <- system.file("Shiny-examples", "MCEM_vizualization.R", package = "emphasis")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `emphasis`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}
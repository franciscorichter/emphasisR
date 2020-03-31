#' @export
emphasis.app <- function() {
  app_dir <- system.file("Shiny-examples",
                        "emphasis_postprocessing.R",
                        package = "emphasis")
  if (app_dir == "") {
    stop("Could not find example directory. Try re-installing `emphasis`.",
          call. = FALSE)
  }

  shiny::runApp(app_dir, display.mode = "normal")
}


#' @export
analytics <- function() {
  app_dir <- system.file("Shiny-examples",
                        "MCEM_vizualization.R",
                        package = "emphasis")
  if (app_dir == "") {
    stop("Could not find example directory. Try re-installing `emphasis`.",
         call. = FALSE)
  }

  shiny::runApp(app_dir, display.mode = "normal")
}

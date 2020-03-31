#' @keywords internal
speciation_rate <- function(tm, tree, pars, model, soc, sum_lambda = FALSE) {
  speciation_r <- get(paste0("lambda.", model))
  lambda <- speciation_r(tm, tree, pars, soc = soc, sum_lambda = sum_lambda)
  return(lambda)
}

#' @keywords internal
sum_speciation_rate <- function(x, tree, pars, model, soc) {
  n <- sapply(x, n_from_time, tree = tree, soc = soc)
  speciation_r <- get(paste0("lambda.", model))
  lambda <- speciation_r(x, tree, pars, soc = soc)
  return(n * lambda)
}

# Speciations rates
#' @keywords internal
lambda.rpd1 <- function(tm, tree, pars, soc, sum_lambda = FALSE) {
  n <- sapply(tm, n_from_time, tree = tree, soc = soc)
  lambda <- max(0, pars[2] + pars[3] * n)
  if (sum_lambda) lambda <- lambda * n
  return(lambda)
}

#' @keywords internal
lambda.rpd5c <- function(tm, tree, pars, soc, sum_lambda = FALSE) {
  pd <- sapply(tm, phylodiversity, tree = tree, soc = soc) - tm
  n <- sapply(tm, n_from_time, tree = tree, soc = soc)
  lambda <- max(0, pars[2] + pars[3] * n + pars[4] * pd / n)
  if (sum_lambda) lambda <- lambda * n
  return(lambda)
}

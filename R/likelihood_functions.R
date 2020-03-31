#' @keywords internal
loglik.tree <- function(model) {
  log.lik <- get(paste0("loglik.tree.", model))
  return(log.lik)
}

# likelihood functions
#' @keywords internal
loglik.tree.rpd5 <- function(pars, tree) {
  if (pars[4] == 0) {
    log.lik <- loglik.tree.rpd1(pars, tree)
  } else {
    to <- tree$to
    to <- utils::head(to, -1)
    to[to == 2] <- 1

    mu <- max(0, pars[1])
    wt <- diff(c(0, tree$brts))

    n <- tree$n
    pt <- c(0, tree$pd[-nrow(tree)])
    pd2 <- pt + n * wt

    brts_i <- tree$brts
    brts_im1 <- c(0, brts_i[-length(brts_i)])

    lambda <- pmax(0, pars[2] + pars[3] * n + pars[4] * pd2 / n)
    rho <- pmax(lambda[-length(n)] * to + mu * (1 - to), 0)
    c1 <- pars[2] + pars[3] * n + (pars[4] / n) * (pt - brts_im1 * n)
    inte <- NULL
    for (i in seq_along(brts_im1)) {
      if ((brts_im1[i] > (-c1[i] / pars[4])) &
          (brts_i[i] < (-c1[i] / pars[4]))) {
        if (pars[4] > 0) {
          brts_im1[i] <- -c1[i] / pars[4]
        } else {
          brts_i[i] <- -c1[i] / pars[4]
        }
      }
      inte[i] <- n[i] * (mu * wt[i] + c1[i] *
                        (brts_i[i] - brts_im1[i]) +
                         pars[4] * (brts_i[i]^2 - brts_im1[i]^2) / 2)
    }
    log.lik <- -sum(inte) + sum(log(rho))
  }
  return(log.lik)
}

#' @keywords internal
loglik.tree.rpd1 <- function(pars, tree) {
  to <- tree$to
  to <- utils::head(to, -1)
  to[to == 2] <- 1

  mu <- max(0, pars[1])
  wt <- diff(c(0, tree$brts))

  n <- tree$n

  lambda <- pmax(0, pars[2] + pars[3] * n)
  rho <- pmax(lambda[-nrow(tree)] * to + mu * (1 - to), 0)

  sigma_over_tree <- n * (mu + lambda) * wt

  log.lik <- -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}

#' @keywords internal
loglik.tree.rpd5c <- function(pars, tree) {
  if (pars[4] == 0) {
    log.lik <- loglik.tree.rpd1(pars, tree)
  } else {
    to <- tree$to
    to <- utils::head(to, -1)
    to[to == 2] <- 1

    mu <- max(0, pars[1])
    wt <- diff(c(0, tree$brts))

    n <- tree$n
    pt <- c(0, tree$pd[-nrow(tree)])
    pd2 <- pt + n * wt

    brts_i <- tree$brts
    brts_im1 <- c(0, brts_i[-length(brts_i)])

    lambda <- pmax(0, pars[2] + pars[3] * n + pars[4] * (pd2 - brts_i) / n)
    rho <- pmax(lambda[-length(n)] * to + mu * (1 - to), 0)
    c1 <- pars[2] + pars[3] * n + (pars[4] / n) * (pt - brts_im1 * n)
    c2 <- pars[4] * ((n - 1) / n)
    inte <- NULL
    for (i in seq_along(brts_im1)) {
      if ((brts_im1[i] > (-c1[i] / c2[i])) & (brts_i[i] < (-c1[i] / c2[i]))) {
        if (c2[i] > 0) {
          brts_im1[i] <- -c1[i] / c2[i]
        } else {
          brts_i[i] <- -c1[i] / c2[i]
        }
      }
      inte[i] <- n[i] * (mu * wt[i] +
                        c1[i] * (brts_i[i] - brts_im1[i]) +
                        c2[i] * (brts_i[i]^2 - brts_im1[i]^2) / 2)
    }
    log.lik <- -sum(inte) + sum(log(rho))
  }
  return(log.lik)
}

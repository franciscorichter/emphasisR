# more utilities
#' @rawNamespace useDynLib(emphasis)
#' @rawNamespace import(Rcpp)
#' @keywords internal
n_from_time <- function(tm, tree, soc) {
  # return N at tm.
  to <- utils::head(tree$to, -1)
  to[to == 2] <- 1
  n <- c(soc, soc + cumsum(to) + cumsum(to - 1))
  output_n <- n[max(which(c(-1, tree$brts) < tm))]
  return(output_n)
}

#' @keywords internal
phylodiversity <- function(tm, tree, soc) {
  i1 <- tree$brts <= tm
  i2 <- tree$to == 0 & i1
  i3 <- tree$t_ext %in% tree$brts[i2]
  dt <- diff(c(0, tree$brts[i1 & !i2 & !i3], tm))
  return(sum(dt * (soc:(length(dt) + soc - 1))))
}

#' @keywords internal
emphasis_bootstrap <- function(input,
                               n_it = 100,
                               print = FALSE,
                               file = "bootstrap_temp.RData") {
  output <- NULL
  for (i in 1:n_it) {
    st <- mcE_step(
      brts = input$brts,
      pars = input$pars,
      sample_size = input$sample_size,
      model = input$model,
      no_cores = input$cores,
      parallel = FALSE,
      soc = input$soc
    )
    if (print == TRUE) {
      print(paste("iteration", i))
      print(paste("log-lik: ", log(st$fhat), sep = ""))
      print(paste("took: ", st$E_time, sep = ""))
    }
    output <- rbind(
      output,
      data.frame(
        fhat = log(st$fhat),
        eitme = st$E_time,
        ss = input$sample_size
      )
    )
    save(output, input, file = file)
  }
  return(output)
}

#' @keywords internal
data_to_table <- function(df, replicant, left, right) {
  df <- df[df$rep == replicant, ]
  df <- df[df$iteration %in% left:right, ]
  summ <- data.frame(
    lfhat = mean(df$fhat),
    sd_fhat = stats::sd(df$fhat),
    mad_fhat = stats::mad(df$fhat),
    replicant = replicant,
    par1 = stats::median(df$par1),
    par2 = stats::median(df$par2),
    par3 = stats::median(df$par3),
    par4 = stats::median(df$par4),
    E_time = sum(df$E_time) / 60,
    M_time = sum(df$M_time) / 60,
    sample_size = mean(df$sample_size)
  )
  return(summ)
}

#' @keywords internal
AIC <- function(LogLik, k) {  # nolint    capital name function
  aic <- (2 * k) - (2 * LogLik)
  return(aic)
}

#' @keywords internal
AICweights <- function(LogLik, k) { # nolint    capital name function
  ic <- AIC(LogLik, k)
  bestmodel_ic <- min(ic)
  weights <- exp(-0.5 * (ic - bestmodel_ic))
  weights <- weights / sum(weights)
  return(weights)
}

#' @keywords internal
AICw <- function(l1, l2, k1, k2) { # nolint function name
  ic <- AIC(c(l1, l2), c(k1, k2))
  bestmodel_ic <- min(ic)
  weights <- exp(-0.5 * (ic - bestmodel_ic))
  weights <- weights / sum(weights)
  return(weights[1])
}

#' @keywords internal
vectors2phylo <- function(list) {
  t <- list$wt
  e <- list$E
  s <- list$S
  ct <- sum(t)
  newick <- paste(sl[1], ";", sep = "")
  n <- 1
  identf <- data.frame(Spec = "aa", Time = 0) # Labels of species
  for (i in 1:(length(t) - 1)) {
    # speciation
    sumt <- sum(t[1:i])
    if (is.null(s)) {
      bd <- sample(1:n, 1)
      species <- as.character(identf[bd, 1])
    } else {
      species <- s[i]
    }
    if (e[i] == 1) {
      ind <- regexpr(species, newick)[1] - 1
      atm <- sumt - identf[which(identf[, 1] == species), 2]
      newick <- paste(substr(newick, 1, ind),
        "(", substr(newick, ind + 1, ind + 4),
        ",", sl[i + 1], "):",
        as.character(atm),
        substring(newick, ind + 5),
        sep = ""
      )
      identf <- rbind(identf, data.frame(
        Spec = substr(sl[i + 1], 1, 2),
        Time = sumt
      ))
      identf[identf$Spec == species, 2] <- sumt
      n <- n + 1
    }
    # extinction
    if (e[i] == 0) {
      ind <- regexpr(species, newick)[1] + 2
      atm <- sumt - identf[which(identf[, 1] == species), 2]
      identf <- identf[!identf$Spec == species, ]
      newick <- paste(substr(newick, 1, ind), as.character(atm),
                      substring(newick, ind + 2), sep = "")
      n <- n - 1
    }
  }
  newick <- compphyl(newi = newick, identf = identf, ct = ct)
  newick <- ape::read.tree(text = newick)
  return(newick)
}

#' @keywords internal
phylo2tree <- function(tree) {
  # to map newick trees into ther xxxx format
  ltt <- ape::ltt.plot.coords(tree)
  t <- diff(ltt[, 1])
  ltt <- ltt[-1, ]
  n <- ltt[, 2]
  e <- diff(n)
  e[e == -1] <- 0
  return(list(wt = t, to = e))
}

#' @keywords internal
tree2phylo <- function(tree, initspec = 1) {
  wt <- -diff(c(0, tree$brts))
  to <- tree$to
  to[to == 2] <- 1
  ct <- sum(wt)
  newick <- paste(sl[1], ";", sep = "")
  N <- 1
  identf <- data.frame(Spec = "a", Time = 0) # Labels of species
  for (i in 1:(length(wt) - 1)) {
    # speciation
    bt <- sum(wt[1:i])
    bd <- sample(1:N, 1)
    species <- as.character(identf[bd, 1])
    if (to[i] == 1) {
      ind <- regexpr(species, newick)[1] - 1
      atm <- bt - identf[which(identf[, 1] == species), 2]
      newick <- paste(substr(newick, 1, ind),
                      "(", substr(newick, ind + 1, ind + 4),
                      ",", sl[i + 1],
                      "):", as.character(atm),
                      substring(newick, ind + 5), sep = "")
      identf <- rbind(identf,
                      data.frame(Spec = substr(sl[i + 1], 1, 2), Time = bt))
      identf[identf$Spec == species, 2] <- bt
      N <- N + 1
    }
    # extinction
    if (to[i] == 0) {
      ind <- regexpr(species, newick)[1] + 2
      atm <- bt - identf[which(identf[, 1] == species), 2]
      identf <- identf[!identf$Spec == species, ]
      newick <- paste(substr(newick, 1, ind),
                      as.character(atm),
                      substring(newick, ind + 2), sep = "")
      N <- N - 1
    }
  }
  newick <- compphyl(newi = newick, identf = identf, ct = ct)
  newick <- ape::read.tree(text = newick)
  return(newick)
}

#' @keywords internal
compphyl <- function(newi, identf, ct) {
  # set to extant species to the present time
  identf[, 1] <- as.character(identf[, 1])
  identf[, 2] <- ct - identf[, 2]
  for (i in seq_along(identf[, 1])) {
    ind <- regexpr(identf[i, 1], newi)[1] + 2
    newi <- paste(substr(newi, 1, ind),
                  as.character(identf[i, 2]),
                  substring(newi, ind + 2), sep = "")
  }
  return(newi)
}

sl <- paste(letters[1], letters, ":0", sep = "")
for (i in 2:26) {
  ll <- paste(letters[i], letters, ":0", sep = "")
  sl <- c(sl, ll)
}

# time calculation
#' @keywords internal
get.time <- function(time,
                     mode = "sec") {
  dif <- proc.time() - time
  ti <- as.numeric(dif[3])
  if (mode == "min") ti <- ti / 60
  if (mode == "hou") ti <- ti / 3600
  return(ti)
}

#' @keywords internal
get.topologies <- function(M) {
  if (M == 0) {
    return(NULL)
  }
  to <- matrix(nrow = 2 * M, ncol = 1)
  to[1, 1] <- 1
  for (i in 2:(2 * M)) {
    comb <- ncol(to)
    for (j in 1:comb) {
      ns <- sum(no.na(to[, j]))
      ne <- sum(1 - no.na(to[, j]))
      if (ns < M & ns > ne) { # extinction or speciation
        to[i, j] <- 1
        to <- cbind(to, matrix(to[, j], ncol = 1))
        to[i, ncol(to)] <- 0
      }
      if (ns == M & ns > ne) { # extinction
        to[i, j] <- 0
      }
      if (ns < M & ns == ne) { # speciation
        to[i, j] <- 1
      }
    }
  }
  return(to)
}


#' @keywords internal
sim.tree <- function(pars,
                     model,
                     ct,
                     soc) {
  tree <- data.frame(brts = 0,
                     to = 1,
                     t_ext = Inf,
                     parent = 0,
                     child = 1)
  cbt <- 0
  n <- soc
  mu <- max(0, pars[1])
  ## sim waiting time,
  spec.cnt <- soc
  while ((cbt < ct) & (n > 0)) {
    tmp_tree <- rbind(tree[-1, ],
                      data.frame(brts = ct,
                                 to = 1,
                                 t_ext = Inf,
                                 parent = NA,
                                 child = NA))
    rate_max <- max(sum_speciation_rate(cbt, tmp_tree, pars, model, soc = soc),
                    sum_speciation_rate(ct, tmp_tree, pars, model, soc = soc)) +
                    mu * n
    u1 <- stats::runif(1)
    next_event_time <- cbt - log(x = u1) / rate_max

    if (next_event_time < ct) {
      u2 <- stats::runif(1)
      pt <- (sum_speciation_rate(next_event_time,
                                 tmp_tree,
                                 pars, model,
                                 soc = soc) + 
               mu * n) / rate_max
      if (u2 < pt) {
        l1 <- speciation_rate(next_event_time,
                              tmp_tree,
                              pars = pars,
                              model = model,
                              soc = soc)
        to <- sample(c(1, 0), size = 1, prob = c(l1, mu) / (l1 + mu))
        if (to == 1) {
          spec.cnt <- spec.cnt + 1
          current.spec <- tree$child[tree$to == 1 & is.infinite(tree$t_ext)]
          tree <- rbind(tree,
                        data.frame(brts = next_event_time,
                                   to = 1,
                                   t_ext = Inf,
                                   parent = sample(current.spec, 1),
                                   child = spec.cnt))
          n <- n + 1
        } else {
          n <- n - 1
          current.spec <- tree$child[tree$to == 1 & is.infinite(tree$t_ext)]
          extinction <- sample(current.spec, 1)
          tree <- rbind(tree, data.frame(brts = next_event_time,
                                         to = 0,
                                         t_ext = Inf,
                                         parent = extinction,
                                         child = NA))
          tree$t_ext[tree$child == extinction] <- next_event_time
        }
      }
    }
    cbt <- next_event_time
  }
  tree <- rbind(tree, data.frame(brts = ct,
                                 to = 1,
                                 t_ext = Inf,
                                 parent = NA,
                                 child = NA))
  return(tree)
}

#' @keywords internal
remove.extinctions <- function(tree) {
  extant_brts <- tree$brts[tree$to == 1 & is.infinite(tree$t_ext)]
  return(extant_brts)
}

#' @keywords internal
multiplot <- function(lp,
                      plotlist = NULL,
                      file, cols = 1,
                      layout = NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(lp, plotlist)

  num_plots <- length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(num_plots / cols)),
      ncol = cols, nrow = ceiling(num_plots / cols)
    )
  }

  if (num_plots == 1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(
      nrow(layout),
      ncol(layout)
    )))

    # Make each plot, in the correct location
    for (i in 1:num_plots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = grid::viewport(
        layout.pos.row = matchidx$row,
        layout.pos.col = matchidx$col
      ))
    }
  }
}

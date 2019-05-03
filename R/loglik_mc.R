sample_brts_exts <- function(S,la,mu,age,topologies)
{
  if(S > 0)
  {
    numcols <- dim(topologies)[2]
    topology <- topologies[,DDD::sample2(x = 1:numcols,size = 1)]
    age <- abs(age)
    brtsff <- sort(runif(n = 2 * S,min = -age,max = 0))
    brts <- brtsff[topology == 1]
    exts <- brtsff[topology == 0]
    rm(brtsff)
    brts_exts <- cbind(brts,exts)
  } else
  {
    brts_exts <- NA
  }
  return(brts_exts)
}

log_sampling_prob <- function(brts,brts_exts,la,mu,age)
{
  S <- compute_dim(brts_exts)
  logprob <- lgamma(S + 1) +
    dpois(x = S,lambda = 10 * mu * exp((la - mu) * age),log = TRUE) -
    2 * S * log(age) + S * log(2)
  return(logprob)
}

loglik_full_tree <- function(brts,brts_exts,la,mu,age)
{
  age <- -abs(age)
  brts <- sort(-abs(brts))
  if(any(is.na(brts_exts)))
  {
    events <- cbind(rbind(brts,1),rbind(0,0))
  } else
  {
    events <- cbind(rbind(brts,1),rbind(brts_exts[,1],1),rbind(brts_exts[,2],-1),rbind(0,0))
  }
  events <- events[,order(events[1,])]
  numspec <- sum(events[2,] == 1)
  logla <- (numspec - 1) * log(la)
  numext <- sum(events[2,] == -1)
  if(mu > 0)
  {
    logmu <- numext * log(mu)
  } else if(numext == 0)
  {
    logmu <- 0
  } else
  {
    logmu <- -Inf
  }
  loglik <- logla + logmu
  S <- 1;
  for(i in 2:dim(events)[2])
  {
    logfactor <- log(S)
    loglik <- loglik + logfactor - S * (la + mu) * (events[1,i] - events[1,i - 1])
    S <- S + events[2,i]
  }
  return(loglik)
}

log_catalan <- function(k)
{
  logc <- -log(k + 1) + lgamma(2 * k + 1) - 2 * lgamma(k + 1)
  return(logc)
}

loglik_reconstructed_tree_approx_S <- function(brts,brts_exts_list,la,mu,age)
{
  endmc <- length(brts_exts_list)
  brts <- -abs(brts)
  S <- dim(brts_exts_list[[1]])[1]
  if(is.null(S))
  {
    S <- 0
  }
  loglik_full_tree_vec <- rep(NA,endmc)
  for(mc in 1:endmc)
  {
    loglik_full_tree_vec[mc] <- loglik_full_tree(brts = brts,brts_exts = brts_exts_list[[mc]],la = la,mu = mu,age = age)
  }
  loglik_full_tree_vec_max <- max(loglik_full_tree_vec)
  if(loglik_full_tree_vec_max != -Inf)
  {
     #logfactor <-  2 * S * log(age) - S * log(2) # Works for S = 1
     logfactor <- log_catalan(S) + 2 * S * log(age) - lgamma(2 * S + 1) # Works for S = 1
     loglik <- logfactor + loglik_full_tree_vec_max + log(mean(exp(loglik_full_tree_vec - loglik_full_tree_vec_max)))
  } else
  {
     loglik <- -Inf
  }
  return(loglik)
}

expfun <- function(la,mu,t,t0 = 0)
{
  res <- exp(-(la - mu) * abs(t - t0))
  return(res)
}

Pnot0 <- function(la,mu,t)
{
  res <- (la - mu) /
    (la - mu * expfun(la,mu,t))
  return(res)
}

ut <- function(la,mu,t)
{
  res <- la * (1 - expfun(la,mu,t)) /
         (la - mu * expfun(la,mu,t))
  return(res)
}

oneminusut <- function(la,mu,t)
{
  res <- (la - mu) * expfun(la,mu,t) /
    (la - mu * expfun(la,mu,t))
  return(res)
}

logP1 <- function(la,mu,t)
{
  return( log(Pnot0(la,mu,t)) + log(oneminusut(la,mu,t)) )
}

loglik_reconstructed_tree_exact <- function(brts,la,mu,age)
{
  brts <- abs(brts)
  S <- length(brts)
  loglik <- (S - 1) * log(la) + sum(logP1(la = la,mu = mu,t = brts))
  return(loglik)
}

compute_dim <- function(brts_exts)
{
  S <- dim(brts_exts)[1]
  if(is.null(S))
  {
    S <- 0
  }
  return(S)
}

loglik_reconstructed_tree_approx <- function(brts,la,mu,age,endmc = 1000,endS = 20)
{
  loglik_approx_S <- rep(0,endS + 1)
  for(S in 0:endS)
  {
    brts_exts_list <- list()
    if(S > 0)
    {
      topologies <- get.topologies(S)
    } else
    {
      topologies <- NULL
    }
    for(mc in 1:endmc)
    {
      brts_exts_list[[mc]] <- sample_brts_exts(S,la,mu,age,topologies)
    }
    loglik_approx_S[S + 1] <- loglik_reconstructed_tree_approx_S(brts = brts,brts_exts_list = brts_exts_list,la = la,mu = mu,age = age)
  }
  loglik_approx_S_max <- max(loglik_approx_S)
  loglik_approx1 <- loglik_approx_S_max + log(sum(exp(loglik_approx_S - loglik_approx_S_max)))
  return(list(loglik_approx1 = loglik_approx1,loglik_approx_S = loglik_approx_S))
}

run_example <- function(
  brts = 10,
  la = 0.2,
  mu = 0.1,
  age = 10,
  endmc = 10000,
  endS = 9
)
{
  loglik_true1 <- loglik_reconstructed_tree_exact(brts = brts,la = la,mu = mu,age = age)
  loglik_true2 <- DDD::bd_loglik(brts = brts,pars1 = c(la = la,mu = mu),pars2 = c(tdmodel = 0,cond = 0,ph = 1,scr = 0,soc = 1),missnumspec = 0)
  
  loglik_approx_out <- loglik_reconstructed_tree_approx(
    brts = brts,
    la = la,
    mu = mu,
    age = age,
    endmc = endmc,
    endS = endS
  )
  loglik_approx1 <- loglik_approx_out$loglik_approx1
  print(paste('Approximation 1 (logmean):',loglik_approx1))
  print(paste('True 1: ',loglik_true1))
  print(paste('True 2: ',loglik_true2))
  
  pNENL <- dd_loglik_EA2009(pars = c(0.2,0.1),lx = 100,age = 10)
  pNL <- colSums(pNENL)
  result <- rbind(loglik_approx_out$loglik_approx_S,log(pNENL[1:(1 + endS),1 + length(brts)]))
  rownames(result) <- c('Approx','True')
  colnames(result) <- 0:endS
  print(result)
  return(list(loglik_approx_out = loglik_approx_out,result = result))
}  
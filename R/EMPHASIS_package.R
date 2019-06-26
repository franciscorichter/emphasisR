### EMPHASIS functions

# negative logLikelihood of a tree
nllik.tree = function(pars,tree,model="dd",initspec=1,K=FALSE){
  nl = -loglik.tree(pars,tree,model,initspec)
  return(nl)
}


loglik.tree <- function(pars,tree,model,initspec=1){
  if(model == "dd"){
    log.lik = loglik.tree.dd(pars,tree,model,initspec)
  }
  if(model == "edd"){
    log.lik = loglik.tree.edd(pars,tree,model,initspec)
  }
  if(model == "pd"){ 
    log.lik = loglik.tree.pd(pars,tree,model,initspec)
  }
  if(model == "epd"){ 
    log.lik = loglik.tree.epd(pars,tree,model,initspec)
  }
  return(log.lik)
}

loglik.tree.dd <- function(pars,tree,model,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  lambda = lambda.dd(pars,n,GLM=TRUE)
  sigma = (lambda + mu)*n
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = (sum(-sigma*wt)+sum(log(rho)))
  return(log.lik)
}

loglik.tree.edd <- function(pars,tree,model,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  lambda = lambda.edd(pars,n)
  sigma = (lambda + mu)*n
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = (sum(-sigma*wt)+sum(log(rho)))
  return(log.lik)
}


loglik.tree.pd <- function(pars,tree,model,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  n.obs = cumsum(c(1,tree$to==2))
  n.obs = head(n.obs,-1)
  pd_obs = cumsum(n.obs*wt)
  pd_miss = n.miss_vec = ext.miss_vec = NULL
  for(j in 1:nrow(tree)){
    if(is.null(n.miss_vec)){
      pd_miss[j] = 0
    }else{
      pd_miss[j] = sum(tree$brts[j] - n.miss_vec)
    }
    if(tree$to[j]==1){
      n.miss_vec = c(n.miss_vec,tree$brts[j])
      ext.miss_vec = c(ext.miss_vec,tree$t_exp[j])
    }
    if(tree$to[j]==0){
      n.miss_vec = n.miss_vec[ext.miss_vec != tree$brts[j]]
    }
  }
  pd = pd_obs + pd_miss
  lambda = pars[1] + pars[2]*pd
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  sigma_over_tree = NULL
  for(i in 1:length(wt)){
    sigma_over_tree[i] = n[i]*(pars[3]*wt[i]+integrate(lambda.pd_t,lower=0,upper=wt[i],pars=pars,pd[i],n=n[i],wt[i])$value)
  }
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}

loglik.tree.epd <- function(pars,tree,model,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  n.obs = cumsum(c(1,tree$to==2))
  n.obs = head(n.obs,-1)
  pd_obs = cumsum(n.obs*wt)
  pd_miss = n.miss_vec = ext.miss_vec = NULL
  for(j in 1:nrow(tree)){
    if(is.null(n.miss_vec)){
      pd_miss[j] = 0
    }else{
      pd_miss[j] = sum(tree$brts[j] - n.miss_vec)
    }
    if(tree$to[j]==1){
      n.miss_vec = c(n.miss_vec,tree$brts[j])
      ext.miss_vec = c(ext.miss_vec,tree$t_exp[j])
    }
    if(tree$to[j]==0){
      n.miss_vec = n.miss_vec[ext.miss_vec != tree$brts[j]]
    }
  }
  pd = pd_obs + pd_miss
  lambda = exp(pars[1] + pars[2]*pd)
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  sigma_over_tree = NULL
  for(i in 1:length(wt)){
    sigma_over_tree[i] = n[i]*(pars[3]*wt[i]+integrate(lambda.epd_t,lower=0,upper=wt[i],pars=pars,pd[i],n=n[i],wt[i])$value)
  }
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}
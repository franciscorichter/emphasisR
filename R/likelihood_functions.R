loglik.tree <- function(model){
  log.lik = get(paste0("loglik.tree.", model))
  return(log.lik)
}

# likelihood functions 

loglik.tree.dd <- function(pars,tree){
  
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  
  lambda = lambda.dd(tm = c(0,tree$brts[-nrow(tree)]),tree = tree,pars = pars)
  
  if(any((lambda[-length(lambda)] == 0) & (tree$to[-length(lambda)] != 0)) | (min(pars)<0)){
    log.lik = -Inf
  }else{
  sigma = (lambda + mu)*sapply(c(0,tree$brts[-nrow(tree)]), n_from_time,tree = tree)
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = (sum(-sigma*wt)+sum(log(rho)))
  }

  return(log.lik)
}

loglik.tree.rpd <- function(pars,tree){
  # parameters
  lambda_0 = pars[1]
  mu = pars[3]
  beta = pars[3]
  a = parrs[4]
  b = pars[5]
  ###
  n = number_of_species(tree)
  Pt = c(0,sapply(tree$brts[-length(tree$brts)], function(x) phylodiversity_t(x,tree)))
  wt = diff(c(0,tree$brts))
  # rho
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  lambda = pmax(0,pars[1]*((pd+1)^(-pars[2])))
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  # sigma 
  sigma = n*((lambda_0+mu)*wt-((lambda_*beta*(n^a))/(n*(b+1)))*((Pt+n*wt)^(b+1)-Pt^(b+1)))
  
  loglik = sum(-sigma)+sum(rho)
  return(loglik)
  
}

loglik.tree.pd <- function(pars,tree,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  Pt = c(0,sapply(tree$brts[-length(tree$brts)], function(x) phylodiversity_t(x,tree)))
  alpha1 = pars[2]*n
  brts_i = tree$brts
  brts_im1 = c(0,tree$brts[-nrow(tree)])
  alpha0 = pars[1] - pars[2]*(Pt-n*brts_im1)
  sigma_over_tree = ((alpha0+mu)*wt - (alpha1/2)*(brts_i^2-brts_im1^2))*n
  lambda = sapply(tree$brts, function(x) lambda.pd_t(x,pars=pars,tree=tree))
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}
###

loglik.tree.pd_numerical <- function(pars,tree,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  sigma_over_tree = lambda = NULL
  bt = 0
  for(i in 1:length(wt)){
    sigma_over_tree[i] = n[i]*(pars[3]*wt[i]+quad(lambda.pd_t,xa=bt,xb=bt+wt[i],pars=pars,tree=tree))
    bt = bt + wt[i]
    lambda[i] = lambda.pd_t(bt,pars,tree)
  }
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}


loglik.tree.dd_old <- function(pars,tree){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  initspec=1
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  lambda = lambda.dd(pars,n)
  sigma = (lambda + mu)*n
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = (sum(-sigma*wt)+sum(log(rho)))
  if(min(pars)<0) log.lik = -Inf
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

loglik.tree.gddx <- function(pars,tree,model,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[2])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  lambda = lambda.gddx(pars,n)
  sigma = (lambda + mu)*n
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = (sum(-sigma*wt)+sum(log(rho)))
  return(log.lik)
}

loglik.tree.epd <- function(pars,tree,initspec=1){
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

loglik.tree.gpdx <- function(pars,tree,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[2])
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
  lambda = pmax(0,pars[1]*((pd+1)^(-pars[2])))
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  sigma_over_tree = NULL
  for(i in 1:length(wt)){
    sigma_over_tree[i] = n[i]*(pars[2]*wt[i]+integrate(lambda.dpdx_t,lower=0,upper=wt[i],pars=pars,pd=pd[i],n=n[i],wt=wt[i])$value)
  }
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}

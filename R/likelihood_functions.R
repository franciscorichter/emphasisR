loglik.tree <- function(model){
  log.lik = get(paste0("loglik.tree.", model))
  return(log.lik)
}

# likelihood functions 

loglik.tree.rpd5 <- function(pars,tree){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  
  mu = max(0,pars[1])
  wt = diff(c(0,tree$brts))
  
  n = tree$n
  brts = tree$brts[-length(tree$brts)]
  Pt = c(0,tree$pd[-nrow(tree)])
  
  brts_i = tree$brts
  brts_im1 = c(0,brts)
  
  lambda = pmax(0,pars[2] + pars[4] * Pt[-1]/n[-length(n)] + pars[3]*n[-length(n)])
  rho = pmax(lambda * to + mu * (1-to),0)
  
  sigma_over_tree = n*((mu+pmax(0,pars[2]+pars[3]*n+(pars[4]/n)*(Pt-brts_im1*n)))*wt + pars[4]*(brts_i^2-brts_im1^2)/2)
  
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}

loglik.tree.rpd1 <- function(pars,tree){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  
  mu = max(0,pars[1])
  wt = diff(c(0,tree$brts))
  
  n = tree$n
  
  lambda = pmax(0, pars[2] + pars[3] * n)
  rho = pmax(lambda[-nrow(tree)] * to + mu * (1-to),0)
  
  sigma_over_tree = n*(mu+lambda)*wt 
  
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}



loglik.tree.rpd5c <- function(pars,tree){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  
  mu = max(0,pars[1])
  wt = diff(c(0,tree$brts))
  
  n = tree$n
  brts = tree$brts[-length(tree$brts)]
  Pt = c(0,tree$pd[-nrow(tree)])
  Ptmt = tree$pd[-nrow(tree)]-brts
  brts_i = tree$brts
  brts_im1 = c(0,brts)
  
  lambda = pmax(0,pars[2] + pars[4] * Ptmt/n[-length(n)] + pars[3]*n[-length(n)])
  rho = pmax(lambda * to + mu * (1-to),0)
  
  sigma_over_tree = n*((pars[1]+pars[2]+pars[3]*n+(pars[4]/n)*(Pt-brts_im1*n))*wt + pars[4]*((n-1)/n)*(brts_i^2-brts_im1^2)/2)
  
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}

############################################################





loglik.tree.rpd2 <- function(pars,tree){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  
  n = tree$n
  brts = tree$brts[-length(tree$brts)]
  Pt = c(0,tree$pd[-nrow(tree)])
  
  brts_i = tree$brts
  brts_im1 = c(0,brts)
  
  lambda = pmax(0,pars[1] - pars[2] * Pt[-1])
  rho = pmax(lambda * to + mu * (1-to),0)
  
  sigma_over_tree = n*( (pars[1]+pars[3]-(pars[2])*(Pt-brts_im1*n))*wt - pars[2]*n*(brts_i^2-brts_im1^2)/2 )
  
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}

loglik.tree.rpd3 <- function(pars,tree){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  
  n = tree$n
  brts = tree$brts[-length(tree$brts)]
  Pt = c(0,tree$pd[-nrow(tree)])

  brts_i = tree$brts
  brts_im1 = c(0,brts)

  lambda = pmax(0,pars[1] - pars[2] * Pt[-1]/n[-length(n)])
  rho = pmax(lambda * to + mu * (1-to),0)
  
  sigma_over_tree = n*((pars[1]+pars[3]-(pars[2]/n)*(Pt-brts_im1*n))*wt - pars[2]*(brts_i^2-brts_im1^2)/2)
  
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}




loglik.tree.rpd2b <- function(pars,tree){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  
  n = tree$n
  brts = tree$brts[-length(tree$brts)]
  Pt = c(0,tree$pd[-nrow(tree)])
  
  brts_i = tree$brts
  brts_im1 = c(0,brts)
  
  lambda = pmax(0,pars[1] - pars[2] * Pt[-1]/n[-length(n)])
  rho = pmax(lambda * to + mu * (1-to),0)
  
  sigma_over_tree = n*((pars[1]+pars[3]-(pars[2]/n)*(Pt-brts_im1*n))*wt - pars[2]*(brts_i^2-brts_im1^2)/2)
  
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}


loglik.tree.rpd4 <- function(pars,tree){
  # parameters
  lambda_0 = pars[1]
  gamma = pars[2]
  mu_0 = pars[3]
  N = tree$n
  Pi = tree$pd+1
  Pim1 = c(0,Pi[-nrow(tree)])+1
  wt = diff(c(0,tree$brts))
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  lambda = lambda_0 - gamma * N/Pi
  rho = pmax(lambda[-length(lambda)]*to+mu_0*(1-to),0)
  sigma = N*( (lambda_0+mu_0)*wt - gamma*(   log(Pim1+N*wt)-log(Pim1)  )  )
  if(min(c(lambda_0,mu_0,gamma))<0) log.lik = -Inf
  loglik = sum(-sigma)+sum(log(rho))
  return(loglik)
  
}



loglik.tree.rpd <- function(pars,tree){
  # parameters
  lambda_0 = pars[1]
  gamma = pars[2]
  mu_0 = pars[3]
  alpha = pars[4]
  beta = pars[5]
  
  N = tree$n
  P = c(0,tree$pd[-nrow(tree)])
  wt = diff(c(0,tree$brts))
  a = atan(alpha)
  b = atan(beta)
  if(b==-1){
    P = P+1
  }
  # rho
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  #lambda = sapply(c(0,tree$brts[-nrow(tree)]), lambda.rpd,tree=tree,pars=pars)
  lambda = pmax(0, lambda_0 - gamma * N^a * (P)^b)
  rho = pmax(lambda[-length(lambda)]*to+mu_0*(1-to),0)
  # sigma 
  sigma = N*( (lambda_0+mu_0)*wt - gamma*(N^(a-1))*((P+N*wt)^(b+1)-P^(b+1))/(b+1) )
  if(min(c(lambda_0,mu_0,gamma))<0) log.lik = -Inf
  loglik = sum(-sigma)+sum(log(rho))
  return(loglik)
  
}

########################

loglik.tree.dd2 <- function(pars,tree){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  
  lambda = lambda.dd(tm = c(0,tree$brts[-nrow(tree)]),tree = tree,pars = pars)
  sigma = (lambda + mu)*number_of_species(tree)
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = (sum(-sigma*wt)+sum(log(rho)))
  if(min(pars)<0) log.lik = -Inf
  return(log.lik)
}

loglik.tree.dd_new <- function(pars,tree){
  
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  
  lambda = lambda.dd(tm = c(0,tree$brts[-nrow(tree)]),tree = tree,pars = pars)
  
  sigma = (lambda + mu)*sapply(c(0,tree$brts[-nrow(tree)]), number_of_species,tree = tree)
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = (sum(-sigma*wt)+sum(log(rho)))
  
  return(log.lik)
}

## work in progress
loglik.tree.erpd <- function(pars,tree){
  # parameters
  lambda_0 = pars[1]
  gamma = 1
  mu_0 = pars[2]
  alpha = pars[3]
  beta = pars[4]
  ###
  n = number_of_species(tree)
  Pt = c(0,sapply(tree$brts[-length(tree$brts)], function(x) phylodiversity(x,tree)))
  wt = diff(c(0,tree$brts))
  a = atan(alpha)
  b = atan(beta)
  # rho
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  lambda = sapply(tree$brts, lambda.rpd2,tree=tree,pars=pars)
  rho = pmax(lambda[-length(lambda)]*to+mu_0*(1-to),0)
  # sigma 
  sigma = n*( (lambda_0+mu_0)*wt - gamma*(n^(a-1))*((Pt+n*wt)^(b+1)-Pt^(b+1))/(b+1) )
  if(min(c(lambda_0,mu_0,gamma))<0) log.lik = -Inf
  loglik = sum(-sigma)+sum(log(rho))
  return(loglik)
  
}

loglik.tree.rpd_old <- function(pars,tree){
  # parameters
  lambda_0 = pars[1]
  mu = pars[3]
  beta = pars[2]
  a = pars[4]
  b = pars[5]
  ###
  n = number_of_species(tree)
  Pt = c(0,sapply(tree$brts[-length(tree$brts)], function(x) phylodiversity(x,tree)))
  wt = diff(c(0,tree$brts))
  # rho
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  lambda = lambda.rpd(tree$brts,tree,pars)
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  # sigma 
  sigma = n*((lambda_0+mu)*wt-((lambda*beta*(n^a))/(n*(b+1)))*((Pt+n*wt)^(b+1)-Pt^(b+1)))
  
  loglik = sum(-sigma)+sum(rho)
  return(loglik)
  
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

loglik.tree.rpd_taylor <- function(pars,tree,prev_pars){
  # parameters
  lambda_0 = pars[1]
  gamma = pars[2]
  mu_0 = pars[3]
  alpha = pars[4]
  beta = pars[5]
  ###
  n = number_of_species(tree)
  Pt = c(0,sapply(tree$brts[-length(tree$brts)], function(x) phylodiversity(x,tree)))
  wt = diff(c(0,tree$brts))
  a = atan(alpha)
  b = atan(beta)
  # rho
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  lambda = sapply(tree$brts, lambda.rpd_taylor,tree=tree,pars=pars)
  rho = pmax(lambda[-length(lambda)]*to+mu_0*(1-to),0)
  # sigma 
  pa = atan(prev_pars[4])
  pb = atan(prev_pars[5])
  gamma0 = prev_pars[2]
  pgamma = gamma+((a-pa)*gamma0*pa)/n
  
  sigma = n*( (lambda_0+mu_0)*wt - pgamma*(n^(pa-1))*((Pt+n*wt)^(pb+1)-Pt^(pb+1))/(pb+1) + (b-pb)*gamma0*pb*(n^pa)*(1/n)*(1/(b+1))*((Pt+n*wt)^(pb+1)-Pt^(pb+1)) )
  
  if(min(c(lambda_0,mu_0,gamma))<0) log.lik = -Inf
  loglik = sum(-sigma)+sum(rho)
  return(loglik)
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

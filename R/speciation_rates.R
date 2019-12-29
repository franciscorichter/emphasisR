speciation_rate <- function(tm,tree,pars,model){
  speciation_r = get(paste0("lambda.", model))
  lambda = speciation_r(tm,tree,pars)
  return(lambda)
}

# Speciations rates 

lambda.dd <- function(tm,tree,pars){
  N = sapply(tm, n_from_time,tree=tree)
  lambdas =  pmax(0, pars[1] - pars[2]*N)
  return(lambdas)
}

lambda.dd.n <- function(pars,n){
  lambdas =  pmax(0, pars[1] - pars[2]*n)
  return(lambdas)
}

lambda.pd <- function(tm,tree,pars){
  phylodiv <- sapply(tm,phylodiversity,tree=tree)
  lambda = max(0,pars[1] - pars[2] * phylodiv)
  return(lambda)
}


lambda.rpd <- function(tm,tree,pars){
  # parameters   pars = c(l1,ga,mu,al,be)
  lambda_0 = pars[1]
  gamma = pars[2]
  #mu = pars[3]
  alpha = pars[4]
  beta = pars[5]
  ###
  a = atan(alpha)
  b = atan(beta)
  ###
  pd = sapply(tm,phylodiversity,tree=tree)
  N = sapply(tm, n_from_time,tree=tree)
  lambda = max(0, lambda_0 - gamma * N^a * (pd+1)^b)
  return(lambda)
}

lambda.rpd2 <- function(tm,tree,pars){
  phylodiv <- sapply(tm,phylodiversity,tree=tree)
  lambda = max(0,pars[1] - pars[2] * phylodiv)
  return(lambda)
}


lambda.rpd3 <- function(tm,tree,pars){
  
  lambda_0 = pars[1]
  gamma = pars[2]
  
  pd = sapply(tm,phylodiversity,tree=tree)
  N = sapply(tm, n_from_time,tree=tree)
  lambda = max(0, lambda_0 - gamma * (pd/N) )
  return(lambda)
}

lambda.rpd5 <- function(tm,tree,pars){
  
  lambda_0 = pars[1]
  gamma = pars[2]
  beta = pars[4]
  
  pd = sapply(tm,phylodiversity,tree=tree)
  N = sapply(tm, n_from_time,tree=tree)
  lambda = max(0, lambda_0 - gamma * (pd/N) - pars[4]*N)
  return(lambda)
}

lambda.rpd2b <- function(tm,tree,pars){
  
  pd = sapply(tm,phylodiversity,tree=tree)
  N = sapply(tm, n_from_time,tree=tree)
  lambda = max(0, pars[1] - pars[2] * (pd/N))
  return(lambda)
}


lambda.rpd4 <- function(tm,tree,pars){

  lambda_0 = pars[1]
  gamma = pars[2]

  pd = sapply(tm,phylodiversity,tree=tree)
  N = sapply(tm, n_from_time,tree=tree)
  lambda = max(0, lambda_0 - gamma * (N/(pd+1)) )
  return(lambda)
}

lambda.erpd <- function(tm,tree,pars){
  # work in progress
  lambda_0 = pars[1]
  #gamma = pars[2]
  #mu = pars[3]
  alpha = pars[3]
  beta = pars[4]
  ###
  a = transform_tan(alpha)
  b = transform_tan(beta)
  ###
  pd = sapply(tm,phylodiversity,tree=tree)
  N = sapply(tm, n_from_time,tree=tree)
  lambda = max(0, lambda_0 - gamma * N^a * pd^b)
  return(lambda)
}

lambda.rpd_taylor <- function(tm,tree,pars,previous_pars){
  # parameters   pars = c(l1,ga,mu,al,be)
  lambda_0 = pars[1]
  gamma = pars[2]
  #mu = pars[3]
  alpha = pars[4]
  beta = pars[5]
  
  palpha = previous_pars[4]
  pbeta = previous_pars[5]
  pgamma = previous_pars[2]
  ###
  a = transform_tan(alpha)
  b = transform_tan(beta)
  
  pa = transform_tan(alpha)
  pb = transform_tan(beta)
  ###
  pd = sapply(tm,phylodiversity,tree=tree)
  N = sapply(tm, n_from_time,tree=tree)
  lambda = max(0, lambda_0 - gamma * N^pa * pd^pb - (a-pa)*pgamma*pa*N^(pa-1)*pd^pb - (b-pb)*pgamma*pb*N^(pa)*pd^(pb-1))
  return(lambda)
}

lambda.rpd_old <- function(tm,tree,pars){
  pd = sapply(tm,phylodiversity,tree=tree)
  N = sapply(tm, n_from_time,tree=tree)
  lambda = max(0,pars[1]*(1-((N^pars[2])*(pd^pars[3]))/pars[4]))
  return(lambda)
}

n_from_time <- function(tm,tree){
  to = top = head(tree$to,-1)
  to[to==2] = 1
  initspec = 1
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  if(tm==max(tree$brts)){
    N = n[max(which(c(0,tree$brts) < tm))]
  }else{
    N = n[max(which(c(0,tree$brts) <= tm))]
  }
} 

number_of_species <- function(tree,tm=NULL){
  initspec = 1
  to = head(tree$to,-1)
  to[to==2] = 1
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  
  
  if(is.null(tm)){
    N=n
  }else{
    if(tm==max(tree$brts)){
      N = n[max(which(c(0,tree$brts) < tm))]
    }else{
      N = n[max(which(c(0,tree$brts) <= tm))]
    }
  }
  return(N)
}


######


sum_speciation_rate <- function(x,tree,pars,model){
  b = max(tree$brts)
  N = number_of_species(tree,x)
  
  speciation_r = get(paste0("lambda.", model))
  lambda = speciation_r(x,tree,pars)
  
  return(N*lambda)
}

lambda.rpd_old <- function(time_m,tree,pars){
  pd = phylodiversity_t(time_m,tree)
  n = number_of_species(tree,tm=time_m)
  lambda = max(0,pars[1]*(1-((n^pars[2])*(pd^pars[3]))/pars[4]))
  return(lambda)
}

lambda.dd_old <- function(pars,n){
  lambdas =  pmax(0, pars[1] - pars[2]*n)
  return(lambdas)
}

lambda.edd <- function(pars,n,GLM=FALSE){
  lambdas =  exp(pars[1] + pars[2]*n)
  return(lambdas)
}

lambda.gddx <- function(pars,n,GLM=FALSE){
  lambdas =  pmax(pars[1]*(n^(-pars[2])),0)
  return(lambdas)
}

lambda.dpdx_t <- function(w_time,pars,pd,n,wt){
  lambda = pmax(pars[1]*((pd - n*(wt-w_time)+1)^(-pars[3])),0)
  return(lambda)
}

lambda.epd_t <- function(w_time,pars,pd,n,wt){
  lambda = exp(pars[1]+pars[2]*(pd - n*(wt-w_time)))
  return(lambda)
}

### old ones for testing

speciation_rate_old <- function(x,tree,pars,model){
  if(stringr:::str_detect(model,"dd")){
    b = max(tree$brts)
    to = top = head(tree$to,-1)
    to[to==2] = 1
    initspec = 1
    n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
    N = n[max(which(c(0,tree$brts) <= x))]
  }
  if(model == "dd"){  # diversity-dependence model
    lambda = lambda.dd(pars,N)
  }
  if(model == "dd1.3"){
    lambda = lambda.dd.1.3(pars,N)
  }
  if(model == "edd"){
    lambda = lambda.edd(pars,N)
  }
  if(model == "edd"){
    lambda = lambda.edd(pars,N)
  }
  if(model == "pd"){
    lambda = lambda.pd_t(time_m=x,pars=pars,tree=tree)
  }
  if(model == "rpd"){
    lambda = lambda.rpd(tm=x,pars=pars,tree=tree)
  }
  return(lambda)
}


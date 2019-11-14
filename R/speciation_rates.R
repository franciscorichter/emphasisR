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
  # parameters
  lambda_0 = pars[1]
  alpha    = pars[2]
  beta     = pars[4]
  gamma    = pars[5]
  ###
  a = (2*exp(alpha))/(1+exp(alpha)) - 1
  b = (2*exp(beta))/(1+exp(beta)) - 1
  ###
  pd = sapply(tm,phylodiversity,tree=tree)
  N = sapply(tm, n_from_time,tree=tree)
  lambda = max(0, lambda_0 - gamma * N^a * pd^b)
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
    lambda = lambda.rpd(time_m=x,pars=pars,tree=tree)
  }
  return(lambda)
}


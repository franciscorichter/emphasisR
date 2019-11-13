speciation_rate <- function(x,tree,pars,model){
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

sum_speciation_rate <- function(x,tree,pars,model,initial_point=TRUE){
  b = max(tree$brts)
  to = top = head(tree$to,-1)
  to[to==2] = 1
  initspec = 1
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  N = n[max(which(c(0,tree$brts) <= x))]
  
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
  if(model == "pd"){ ## ??????????????
    lambda = lambda.pd_t(time_m=x,pars=pars,tree=tree)
  }
  return(N*lambda)
}

# Speciations rates 

lambda.rpd <- function(time_m,tree,pars){
  pd = phylodiversity_t(time_m,tree)
  n = number_of_species(tree,tm=time_m)
  lambda = max(0,pars[1]*(1-((n^pars[2])*(pd^pars[3]))/pars[4]))
  return(lambda)
}

lambda.pd_t <- function(time_m,pars,tree){
  lambda = max(0,pars[1] - pars[2] * phylodiversity_t(time_m = time_m,tree = tree))
  return(lambda)
}

lambda.dd <- function(pars,n){
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


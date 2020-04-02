speciation_rate <- function(tm,tree,pars,model,soc,sum_lambda=FALSE){
  speciation_r = get(paste0("lambda.", model))
  lambda = speciation_r(tm,tree,pars,soc=soc,sum_lambda=sum_lambda)
  return(lambda)
}


# Speciations rates 

lambda.rpd1 <- function(tm,tree,pars,soc,sum_lambda=FALSE){

  N = sapply(tm, n_from_time,tree=tree,soc=soc)
  lambda = max(0, pars[2] + pars[3]*N )
  if(sum_lambda) lambda <- lambda*N
  return(lambda)
  
}

lambda.rpd5c <- function(tm,tree,pars,soc,sum_lambda=FALSE){
  
  pd = sapply(tm,phylodiversity,tree=tree,soc=soc)-tm
  N = sapply(tm, n_from_time,tree=tree,soc=soc)
  lambda = max(0, pars[2] + pars[3]*N + pars[4] * pd/N ) 
  if(sum_lambda) lambda <- lambda*N
  return(lambda)
  
}

##########################################################
##### lineages dependent speciation rate #################
##########################################################


lambda.ldpd <- function(tm,tree,pars,soc,sum_lambda=FALSE){
  gpd = GPD2(tm,tree)
  lambdas = colSums(gpd)
  if(sum_lambda) lambdas = sum(lambdas)
  return(lambdas)
}





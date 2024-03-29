speciation_rate <- function(tm,tree,pars,model,sum_lambda=FALSE,soc=2){
  speciation_r = get(paste0("lambda.", model))
  lambda = speciation_r(tm,tree,pars,soc=soc,sum_lambda=sum_lambda)
  return(lambda)
}

extinction_rate <- function(tm=NULL,tree=NULL,pars,model=NULL,sum_rate=FALSE,soc=2){
  mu = max(0,pars[1])
  if(sum_rate){
    N = sapply(tm, n_from_time,tree=tree,soc=soc)
    mu = N*mu
  }
  return(mu)
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


nh_rate <- function(x,model,pars,tree){
  speciation_rate(tm=x,tree = tree,pars = pars,model = model,soc=tree$n[1],sum_lambda = TRUE)*(1-exp(-(max(tree$brts)-x)*pars[1]))
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





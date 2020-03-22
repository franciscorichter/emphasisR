speciation_rate <- function(tm,tree,pars,model,soc){
  speciation_r = get(paste0("lambda.", model))
  lambda = speciation_r(tm,tree,pars,soc=soc)
  return(lambda)
}

sum_speciation_rate <- function(x,tree,pars,model,soc){
  N = sapply(x, n_from_time,tree=tree,soc=soc)
  speciation_r = get(paste0("lambda.", model))
  lambda = speciation_r(x,tree,pars,soc=soc)
  return(N*lambda)
}

# Speciations rates 

lambda.rpd1 <- function(tm,tree,pars,soc){

  N = sapply(tm, n_from_time,tree=tree,soc=soc)
  lambda = max(0, pars[2] + pars[3]*N )
  return(lambda)
  
}

lambda.rpd5c <- function(tm,tree,pars,soc){
  
  pd = sapply(tm,phylodiversity,tree=tree,soc=soc)-tm
  N = sapply(tm, n_from_time,tree=tree,soc=soc)
  lambda = max(0, pars[2] + pars[3]*N + pars[4] * pd/N ) 
  return(lambda)
  
}

#############################


augment_tree <- function(
  brts,
  pars,
  model,
  soc){
#  setTimeLimit(elapse=500, trans=T)
  mu = max(0,pars[1])
  brts = cumsum(-diff(c(brts,0)))
  b = max(brts)
  cbt = 0
  
  tree = data.frame(brts = brts,
                    t_ext = rep(Inf,length(brts)),
                    to = rep(2,length(brts)))
  tree = tree[order(tree$brts),]
  
  num_missing_branches <- 0
  
  while(cbt < b){
    next_bt = min(tree$brts[tree$brts>cbt])
    
    l1 <- speciation_rate(tm = cbt,tree = tree,pars = pars,model = model,soc = soc,sum_lambda = TRUE)*(1-exp(-mu*(b-cbt)))
    l2 <- speciation_rate(tm = next_bt,tree = tree,pars = pars,model = model,soc = soc,sum_lambda = TRUE)*(1-exp(-mu*(b-next_bt)))
    
    lambda_max = max(l1, l2)
    ###
    if(lambda_max>100){
      stop("Current parameters leds to a huge speciation rate")
    }
    ####
    u1 = runif(1)
    next_speciation_time = cbt - log(x = u1)/lambda_max
    if(next_speciation_time < next_bt){  ## 
      u2 = runif(1)
      pt = speciation_rate(next_speciation_time,tree,pars,model,soc=soc,sum_lambda = TRUE)*(1-exp(-mu*(b-next_speciation_time)))/lambda_max
      if(u2<pt){
        ###
        # here a calculation of the proper limit b should be calculated, otherwise the sampler might provide trees with likelihood zero.
        
        ###
        extinction_time = next_speciation_time + truncdist::rtrunc(1,"exp",a = 0, b = (b-next_speciation_time),rate=mu)
        
        to_add_1 <- c(next_speciation_time, extinction_time, 1)
        to_add_2 <- c(extinction_time, Inf, 0)
        tree <- rbind(tree, to_add_1, to_add_2)
        tree <- tree[order(tree$brts),]
        
        
      }
      num_missing_branches <- num_missing_branches + 1
      if(num_missing_branches>10000){
        stop("Current parameters leds to a large number of species")
      }
    }
    cbt = min(next_speciation_time, next_bt)
  }
  
  tree$pd = sapply(tree$brts, phylodiversity, tree=tree,soc=soc)
  tree$n = sapply(tree$brts, n_from_time,tree=tree,soc=soc)
  
  return(list(tree=tree))
}


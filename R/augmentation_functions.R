augment_tree2 <- function(
  brts,
  pars,
  model,
  soc = 2,
  max_species=100000){  # soc parameter is in this version of emphasis but not in the next one.

  ### input set-up
  time = proc.time()
  brts = cumsum(-diff(c(brts,0)))
  ct = b = max(brts)
  cbt = 0
  
  tree = data.frame(brts = brts,
                    t_ext = rep(Inf,length(brts)),
                    to = rep(2,length(brts)))
  tree = tree[order(tree$brts),]
  
  num_missing_branches = 0
  mu = extinction_rate(tm = NULL,pars = pars,model = NULL,sum_rate = FALSE,soc = NULL,tree = NULL) #constant extinction rate case
  
  #######
  
  while(cbt < ct){
    
    next_bt = min(tree$brts[tree$brts>cbt])
    next_speciation_time = draw_miss_spe(tree,cbt,ct,pars,next_bt,model,soc)
    
    if(next_speciation_time < next_bt){  ## 
      
      extinction_time = draw_extiction(next_speciation_time,ct,mu)
      to_add_1 = c(next_speciation_time, extinction_time, 1)
      to_add_2 = c(extinction_time, Inf, 0)
      tree = rbind(tree, to_add_1, to_add_2)
      tree = tree[order(tree$brts),]
      num_missing_branches <- num_missing_branches + 1
      if(num_missing_branches>max_species){
        stop("Current parameters leds to a large number of species")
      }  
    }
    cbt = min(next_speciation_time, next_bt)
  }
  
  ####### traits augmentation 
  tree$pd = sapply(tree$brts, phylodiversity, tree=tree,soc=soc)
  tree$n = sapply(tree$brts, n_from_time,tree=tree,soc=soc)
  #######
  
  return(list(tree=tree,time=get.time(time)))
}

nh_speciation_rate <- function(time,tree,pars,model,soc=2){
  speciation_rate(time,
                  tree,
                  pars,
                  model,
                  soc=soc,
                  sum_lambda = TRUE)*(1-exp(-pars[1]*(max(tree$brts)-time)))
}

draw_miss_spe <- function(tree,cbt,ct,pars,next_bt,model,soc=2){
  key = 0
  while(key == 0 & cbt < next_bt){
  
    lambda_max = optim(cbt,
                     fn = nh_speciation_rate,
                     pars = pars,
                     tree = tree,
                     model = model,
                     soc = soc,
                     lower = cbt,
                     upper = next_bt,
                     method ="L-BFGS-B")$value
  
    u1 = runif(1)
    cbt = cbt - log(x = u1)/lambda_max
    u2 = runif(1)
    pt = nh_speciation_rate(cbt,
                       tree,
                       pars,
                       model,
                       soc=soc)/lambda_max
      if(u2<pt){
        key = 1
      }
  }
  return(cbt)
  
}

draw_extiction <- function(cbt,ct,mu){
 cbt = cbt + truncdist::rtrunc(1,"exp",a = 0, b = (ct-cbt),rate=mu)
 return(cbt)
}



##################################################################

augment_tree <- function(
  brts,
  pars,
  model,
  optim_nhpp = FALSE,
  soc){  # soc parameter is in this version of emphasis but not in the next one.
  time = proc.time()
  brts = cumsum(-diff(c(brts,0)))
  ct = max(brts)
  cbt = 0
  
  tree = data.frame(brts = brts,
                    t_ext = rep(Inf,length(brts)),
                    to = rep(2,length(brts)))
  tree = tree[order(tree$brts),]
  
  num_missing_branches <- 0
  mu = extinction_rate(tm = NULL,pars = pars,model = NULL,sum_rate = FALSE,soc = NULL,tree = NULL) #constant extinction rate case
  while(cbt < ct){
    next_bt = min(tree$brts[tree$brts>cbt])
    b = ct
    #  b = limit_bt(cbt=cbt,model=model,tree=tree,pars=pars)
    if(optim_nhpp){
      #### optim step 
      lambda_max = optim(cbt,
                         fn=nh_speciation_rate,
                         pars=pars,
                         tree=tree,
                         model=model,
                         soc = soc,
                         lower = cbt,
                         upper = next_bt,
                         method="L-BFGS-B")$value 
    }else{
      l1 <- speciation_rate(tm = cbt,tree = tree,pars = pars,model = model,soc = soc,sum_lambda = TRUE)*(1-exp(-mu*(b-cbt)))
      l2 <- speciation_rate(tm = next_bt,tree = tree,pars = pars,model = model,soc = soc,sum_lambda = TRUE)*(1-exp(-mu*(b-next_bt)))
      lambda_max = max(l1, l2)
    }
    
    ####
    u1 = runif(1)
    next_speciation_time = cbt - log(x = u1)/lambda_max
    if(next_speciation_time < next_bt){  ## 
      u2 = runif(1)
      pt = speciation_rate(next_speciation_time,tree,pars,model,soc=soc,sum_lambda = TRUE)*(1-exp(-mu*(b-next_speciation_time)))/lambda_max
      if(u2<pt){
        extinction_time = next_speciation_time + truncdist::rtrunc(1,"exp",a = 0, b = (b-next_speciation_time),rate=mu)
        
        to_add_1 <- c(next_speciation_time, extinction_time, 1)
        to_add_2 <- c(extinction_time, Inf, 0)
        tree <- rbind(tree, to_add_1, to_add_2)
        tree <- tree[order(tree$brts),]
        num_missing_branches <- num_missing_branches + 1
        if(num_missing_branches>10000){
          stop("Current parameters leds to a large number of species")
        }
        
      }

    }
    cbt = min(next_speciation_time, next_bt)
  }
  
  tree$pd = sapply(tree$brts, phylodiversity, tree=tree,soc=soc)
  tree$n = sapply(tree$brts, n_from_time,tree=tree,soc=soc)
  
  return(list(tree=tree,time=get.time(time)))
}

#########################################


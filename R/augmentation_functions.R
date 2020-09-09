augment_tree <- function(
  brts,
  pars,
  model,
  soc = 2,
  max_species = 100000,
  sampler_spe,
  sampler_ext = "emph"){  
  
  ### input set-up
  time = proc.time()
  brts = cumsum(-diff(c(brts,0)))
  ct = max(brts)
  cbt = 0
  
  tree = data.frame(brts = brts,
                    t_ext = rep(Inf,length(brts)),
                    to = rep(2,length(brts)))
  tree = tree[order(tree$brts),]
  
  num_missing_branches = 0
  
  #######
  
  while(cbt < ct){
    
    next_bt = min(tree$brts[tree$brts>cbt])
    next_speciation_time = draw_speciation(tree,
                                         cbt,
                                         ct,
                                         pars,
                                         next_bt,
                                         model,
                                         soc,
                                         sampler=sampler_spe)
    
    if(next_speciation_time < next_bt){  ## 
      
      extinction_time = draw_extiction(tree=tree,
                                       cbt=next_speciation_time,
                                       ct=ct,
                                       pars = pars,
                                       next_bt=ct,
                                       model = model,
                                       soc = soc,
                                       sampler = sampler_ext)
            
      to_add_1 = c(next_speciation_time, extinction_time, 1)
      to_add_2 = c(extinction_time, Inf, 0)
      tree = rbind(tree, to_add_1, to_add_2)
      tree = tree[order(tree$brts),]
      num_missing_branches = num_missing_branches + 1
      
      if(num_missing_branches>max_species){
        stop("Current parameters leds to a large number of species")
      }  
    }
    cbt = min(next_speciation_time, next_bt)
  }
  
  ####### traits augmentation 
  tree$pd = sapply(tree$brts, phylodiversity, tree=tree,soc=soc)
  tree$n = sapply(tree$brts, n_from_time, tree=tree, soc=soc)
  #######
  
  return(list(tree=tree,time=get.time(time)))
}


#######

draw_speciation <- function(tree,cbt,ct,pars,next_bt,model,soc=2,sampler="emph"){
  
  key = 0
  if(sampler == "emph"){ 
    nsr = nh_speciation_rate
  }
  if(sampler == "nee"){
    nsr = nh_speciation_rate_nee
  }
  if(sampler == "kendall"){
    nsr = nh_speciation_rate_kendall
  }
  while(key == 0 & cbt < next_bt){
  
    lambda_max = optim(cbt,
                     fn = nsr,
                     pars = pars,
                     tree = tree,
                     model = model,
                     soc = soc,
                     lower = cbt,
                     upper = next_bt,
                     method ="L-BFGS-B",
                     control=list(fnscale=-1))$value
  
    u1 = runif(1)
    cbt = cbt - log(x = u1)/lambda_max
    u2 = runif(1)
    pt = nsr(cbt,
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

draw_extiction <- function(tree,cbt,ct,pars,next_bt,model,soc,sampler="emph"){
  
  if(sampler=="unif"){
    cbt = runif(1,min = cbt,max = ct)
  }
  if(sampler=="emph"){
   cbt = cbt + truncdist::rtrunc(1,"exp",a = 0, b = (ct-cbt),rate=pars[1])
  }
  if(sampler=="emph2"){
    la = max(0,pars[2]+pars[3]*n_from_time(cbt,tree,soc))
    mu = pars[1]

    mu_t = mu/(1-p(tc,tp))
    cbt = cbt + truncdist::rtrunc(1,"exp",a = 0, b = (ct-cbt),rate=mu_t)
  }
  if(sampler == "kendall"){
    key = 0
    nsr = get(paste("nh_extinction_rate_",sampler))
   
  }

  
 return(cbt)
  
}

nh_speciation_rate_kendall <- function(time,tree,pars,model=0,soc=1,ti=NULL){
  
  mu <- function(pars) max(0,pars[1])
  lambda <- function(pars) max(0,pars[2]+pars[3]*n_from_time(time,tree,soc))
  rate = lambda(pars) * (1-p_survival(time = time,mu = mu,lambda = lambda,ct = max(tree$brts),pars=pars))
  
  return(rate)
  
}


nh_extinction_rate_kendall <- function(time,tree,pars,model=0,soc=1,ti){
  
  mu <- function(pars) max(0,pars[1])
  lambda <- function(pars) max(0,pars[2]+pars[3]*n_from_time(time,tree,soc))
  rate = mu(pars) / (1-p_survival(time = time,mu = mu,lambda = lambda,ct = max(tree$brts),pars=pars))
  
  return(rate)
  
}

p_survival <- function(tc, tp, model, pars, traits, ...){
  survival = 1 / ( 1 +  pracma:::quad(f = Vectorize(p_int),
                           xa = time,
                           xb = tp,
                           time = tc,
                           model = model,
                           pars = pars,
                           tree = traits$tree, 
                           soc = traits$soc, ...) )
  return(survival)
}

p_int <- function(tau,time,model,pars, tree, soc,...){
  sigma <- function(time){
    speciation_rate(tm = time,
                    tree = tree,
                    pars=pars,
                    model=model,
                    soc = soc)+
    extinction_rate(tm = time,
                    tree = tree,
                    pars = pars,
                    model = model,
                    soc = soc)
    return(sigma)
  }
  p_int_return <- function(tau){
    extinction_rate(tm = time,
                    tree = tree,
                    pars=pars,
                    model=model,
                    soc = soc)*
      exp(-pracma:::quad(f = Vectorize(sigma),
                            xa = time,
                            xb = tau))
  }
  return(p_int_return)
}


##################################################################
## Nee et all augmentation

nh_speciation_rate_nee <- function(time,tree=0,pars,model=0,soc=1){
  ct=max(tree$brts)
  la = max(0,pars[2]+pars[3]*n_from_time(time,tree,soc))
  mu = pars[1]
  
  rate = la*(mu*(1-exp(-(la-mu)*(ct-time))))/(la-mu*exp(-(la-mu)*(ct-time)))
  if(la-mu*exp(-(la-mu)*(ct-time)) < 0.0001){
    rate = 99
  }
  
  N = n_from_time(time,tree,soc)+1
  
  return(rate*N)
  
}

nh_extinction_rate_nee <- function(time,tree=0,pars,model=0,soc=0){
  
  ct=max(tree$brts)
  la = max(0,pars[2]+pars[3]*n_from_time(time,tree,soc))
  mu = pars[1]
  
  rate = (la-mu*exp(-(la-mu)*(ct-time)))/(1-exp(-(la-mu)*(ct-time)))
  if(is.na(rate)){
    rate = 0.000001
  }
  
  return(rate)
  
}

nh_speciation_rate <- function(time,tree,pars,model,soc){
  
  speciation_rate(time,
                  tree,
                  pars,
                  model,
                  soc=soc,
                  sum_lambda = TRUE)*(1-exp(-pars[1]*(max(tree$brts)-time)))
}
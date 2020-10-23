augment_tree <- function(
  brts,
  pars,
  model,
  soc,
  max_species = 100000,
  sampler_spe,
  sampler_ext = "emph2"){  
  
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

draw_speciation <- function(tree,cbt,ct,pars,next_bt,model,soc,sampler="emph"){
  
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
  if(sampler == "fullprocess"){
    spec_rates = speciation_rate(tm = next_event_time,
                                 tree = temp_tree,
                                 pars = pars,
                                 model = model,
                                 sum_lambda = FALSE)
    
    ext_rates =  extinction_rate(tm = next_event_time,
                                 tree = temp_tree,
                                 pars = pars,
                                 model = model,
                                 sum_rate = FALSE)
    
    nsr = speciation_rate() + extinction_rate(sum_rate=TRUE)
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
    if(lambda_max==0){
      cbt = Inf
    }else{
      cbt = cbt - log(x = u1)/lambda_max
    }
    
    if(cbt < next_bt){
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
  }
  
  return(cbt)
  
}

draw_extiction <- function(tree,cbt,ct,pars,model,soc,sampler="emph"){
  
  if(sampler=="unif"){
    cbt = runif(1,min = cbt,max = ct)
  }
  if(sampler=="emph"){
   cbt = cbt + truncdist::rtrunc(1,"exp",a = 0, b = (ct-cbt),rate=pars[1])
  }
  if(sampler=="emph2"){
    la = max(0,pars[2]+pars[3]*n_from_time(cbt,tree,soc))
    mu = max(0,pars[1])

    mu_t = mu/(1-p(time = cbt,ct = ct,lambda = la,mu = mu))
    cbt = cbt + truncdist::rtrunc(1,"exp",a = 0, b = (ct-cbt),rate=mu_t)
  }
  if(sampler == "kendall"){
    mu_t = nh_extinction_rate_kendall(time = cbt,
                                      tree = tree,
                                      pars = pars,
                                      model = model,
                                      soc = soc,
                                      ti = cbt)
    cbt = cbt + truncdist::rtrunc(1,"exp",a = 0, b = (ct-cbt),rate=mu_t)
    
   
  }

  
 return(cbt)
  
}

p <- function(time,ct,lambda,mu){
  if(lambda != mu){
    pv = (lambda-mu)/(lambda - mu*exp(-(lambda-mu)*(ct-time)))
  }else{
    pv = 1/(1+mu*(ct-time))
  }
  return(pv)
}

nh_speciation_rate_kendall <- function(time,tree,pars,model,soc){
  rate = speciation_rate(tm = time,tree = tree,
                         pars = pars,
                         model = model,
                         sum_lambda = FALSE,
                         soc = soc) * (1-p_survival(tc = time,
                                      model = model,
                                      pars = pars,
                                      traits = list(tree=tree,soc = soc)))
  
  N = n_from_time(time,tree,soc)+1
  
  return(rate*N)
  
}


nh_extinction_rate_kendall <- function(time,tree,pars,model=0,soc=1,ti){
  if(is.null(ti)) ti = time
  psurv = p_survival(tc = ti,
                     model = model,
                     pars = pars,
                     traits = list(tree=tree,soc = soc))
  if(psurv == 1){
    rate = 999
  }else{
    rate = extinction_rate(tm = time,tree = tree,
                         pars = pars,
                         model = model,
                         sum_rate = FALSE,
                         soc = soc) / (1-psurv)
  }
  
  return(rate)

}


p_survival <- function(tc, model, pars, traits){
  survival = 1 / ( 1 +  p_int(time = tc,
                              model = model,
                              pars = pars,
                              tree = traits$tree,
                              soc = traits$soc) )
  return(survival)
}

p_int <- function(time,model,pars, tree, soc, ...){
  integrand <- function(tau){
    extinction_rate(tm = tau,
                    tree = tree,
                    pars = pars,
                    model = model,
                    soc = soc)*
      exp(intensity.rho.numerical(tm = time,
                                  tau = tau,
                                  tree = tree,
                                  pars = pars,
                                  model = model))
  }
  val = pracma:::quad(f = Vectorize(integrand),xa = time,xb = max(tree$brts))
  
     
  return(val)
}

intensity.rho.numerical <- function(tm,tau,tree, pars, model){
  nh_rate <- function(x){
    extinction_rate(tm=x,
                    tree = tree,
                    pars = pars,
                    model = model,
                    soc=1,
                    sum_rate = FALSE)-
    speciation_rate(tm=x,
                    tree = tree,
                    pars = pars,
                    model = model,
                    soc=1,
                    sum_lambda = FALSE)
      
  }
  brts = c(tree$brts[tree$brts>tm & tree$brts<tau],tau)
  brts_i = brts
  brts_im1 = c(tm,brts_i[-length(brts_i)])
  inte = vector(mode="numeric",length = length(brts_i))
  for(i in 1:length(brts_i)){
    inte[i] = pracma:::quad(f = Vectorize(nh_rate),xa = brts_im1[i]+0.00000000001,xb = brts_i[i])
  }
  return(sum(inte))
}

##################################################################
## Nee et all augmentation

nh_speciation_rate_nee <- function(time,tree,pars,model,soc){
  ct = max(tree$brts)
  la = max(0,pars[2]+pars[3]*n_from_time(time,tree,soc))
  mu = max(0,pars[1])
  if(is.na(la)){
    print(time)
    print(soc)
    print(tree)
  }
  if(la != mu){
    rate = la*(mu*(1-exp(-(la-mu)*(ct-time))))/(la-mu*exp(-(la-mu)*(ct-time)))
  }else{
    pf = 1/(1+mu*(ct-time))
    rate = la*(1-pf)
  }
  N = n_from_time(time,tree,soc)+1
  
  return(rate*N)
  
}

nh_extinction_rate_nee <- function(time,tree=0,pars,model=0,soc=0){
  
  ct = max(tree$brts)
  la = max(0,pars[2]+pars[3]*n_from_time(time,tree,soc))
  mu = pars[1]
  if(la != mu){
    rate  = (la-mu*exp(-(la-mu)*(ct-time)))/(1-exp(-(la-mu)*(ct-time)))
  }else{
    pf = 1/(1+mu*(ct-time))
    rate = mu/(1-pf)
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
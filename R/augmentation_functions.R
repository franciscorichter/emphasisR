

mc_sample_independent_trees <- function(brts,pars,nsim=1000,model="dd",importance_sampler="emphasis",no_cores=2,pars3=NULL,maxnumspec=NULL,seed=0,method="inverse",parallel=TRUE){
  
  time=proc.time()
  if(seed>0) set.seed(seed)
  if(method == "thinning"){
    E = mc_augmentation_thinning(brts = brts,pars = pars,model = model,importance_sampler = importance_sampler,sample_size = nsim,parallel = parallel,no_cores = no_cores)
  }else{
    E = mc_augmentation_inverse(brts = brts,pars = pars,model = model,importance_sampler = importance_sampler,sample_size = nsim,parallel = parallel,no_cores = no_cores)
  }
  return(E)
}


mc_augmentation_thinning <- function(brts,pars,model,importance_sampler,sample_size,parallel=FALSE,no_cores=2){
  time = proc.time()
  if(!parallel){
    st =  lapply(1:sample_size,function(i){augment_tree_thinning(brts = brts,pars = pars,model=model)} )
  }else{
    st = mclapply(1:sample_size,function(i){augment_tree_thinning(brts = brts,pars = pars,model=model)},mc.cores = no_cores)
  }

  trees = lapply(st,function(list) list$tree)
  dim = sapply(st,function(list) nrow(list$tree))
  E_time = get.time(time)
  
  ####
  log_lik_tree <- loglik.tree(model)
  logf = sapply(trees,log_lik_tree, pars=pars)
  logg =    sapply(st,function(list) list$logg)
  log_weights = logf-logg
  w = exp(log_weights)
  ####
  
  En = list(weights=w,trees=trees,fhat=mean(w),logf=logf,logg=logg,dim=dim,E_time=E_time)
  return(En)
}

mc_augmentation_inverse <-
  function(brts,
           pars,
           model,
           importance_sampler,
           sample_size,
           parallel = TRUE,
           no_cores = 2) {
    #### parallel set-up
    time = proc.time()
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    ##
    if (importance_sampler == "emphasis") {
      trees <- foreach(i = 1:sample_size, combine = list) %dopar% {
        df = emphasis::tree_augmentation_inverse(brts = brts,
                                                 pars = pars,
                                                 model = model)
        df$pd = sapply(df$brts, function(x)
          emphasis:::phylodiversity(x, df))
        df$n = c(1,sapply(df$brts[-nrow(df)], function(x)
          emphasis:::number_of_species(tree = df, tm = x)))
        if (!is.null(df)) {
          logg.samp = emphasis:::log_sampling_prob_nh(df = df,
                                                      pars = pars,
                                                      model = model)
          tree.info = list(logg.samp = logg.samp,
                           dim = nrow(df),
                           tree = df)
        } else{
          tree.info = list(logg.samp = 0,
                           tree = 0,
                           dim = 0)
        }
        return(tree.info)
      }
    }
    if (importance_sampler == "uniform") {
      trees <- foreach(i = 1:sample_size, combine = list) %dopar% {
        df = emphasis:::sample.uniform(brts, maxnumspec = maxnumspec)
        log.samp.unif.prob = emphasis:::log.sampling.prob.uniform(df,
                                                                  maxnumspec = maxnumspec,
                                                                  initspec = 1,
                                                                  p = 0.5)
        df$t_exp = rep(Inf, nrow(df))
        missing_speciations = NULL
        for (j in 1:nrow(df)) {
          if (df$to[j] == 1) {
            missing_speciations = c(missing_speciations, df$brts[j])
          }
          if (df$to[j] == 0) {
            sample_ext_time_index = sample(length(missing_speciations), 1)
            df$t_exp[j] = missing_speciations[sample_ext_time_index]
            missing_speciations = missing_speciations[-sample_ext_time_index]
          }
        }
        ##################
        return(list(
          logg.samp = log.samp.unif.prob,
          dim = nrow(df),
          tree = df
        ))
      }
    }
    stopCluster(cl)
    tree = lapply(trees, function(list)
      list$tree)
    log_lik_tree <- loglik.tree(model)
    logf = sapply(tree, log_lik_tree, pars = pars)
    logg = sapply(trees, function(list)
      list$logg.samp)
    dim = sapply(trees, function(list)
      list$dim)
    logweights = logf - logg
    weights = exp(logweights)
    fhat = mean(weights)
    time = get.time(time)
    E = list(
      weights = weights,
      logweights = logweights,
      fhat = fhat,
      logf = logf,
      fhat.se = 1,
      logg = logg,
      trees = tree,
      dim = dim,
      E_time = time
    )
    return(E)
  }



n_at_brts <- function(tree){
  sapply(tree$brts, function(x) number_of_species(tree))
}

pd_at_brts <- function(tree){
  sapply(tree$brts, function(x) phylodiversity(x,tree))
}

## thinning method

augment_tree_thinning <- function(brts,pars,model="dd"){
  mu = pars[3]
  if(model=="rpd2"){
    mu = pars[2]
  }
  brts = cumsum(-diff(c(brts,0)))
  b = max(brts)
  missing_branches = data.frame(speciation_time=NULL,extinction_time=NULL)
  cbt=0
  
  while(cbt < b){
    
    tree = data.frame(brts = c(missing_branches$speciation_time,brts,missing_branches$extinction_time),
                      t_ext = c(missing_branches$extinction_time, rep(Inf,length(brts)+nrow(missing_branches))),
                      to = c(rep(1,nrow(missing_branches)),rep(2,length(brts)),rep(0,nrow(missing_branches))))
    tree = tree[order(tree$brts),]
    
    next_bt = min(tree$brts[tree$brts>cbt])
    lambda_max = max(sum_speciation_rate(cbt,tree,pars,model)*(1-exp(-mu*(b-cbt))),sum_speciation_rate(next_bt,tree,pars,model)*(1-exp(-mu*(b-next_bt))))
    u1 = runif(1)
    next_speciation_time = cbt-log(x = u1)/lambda_max
    if(next_speciation_time < next_bt){
      u2 = runif(1)
      pt = sum_speciation_rate(next_speciation_time,tree,pars,model)*(1-exp(-mu*(b-next_speciation_time)))/lambda_max
      if(u2<pt){
        extinction_time = next_speciation_time + truncdist::rtrunc(1,"exp",a = 0, b = (b-next_speciation_time),rate=mu)
        missing_branches = rbind(missing_branches,data.frame(speciation_time=next_speciation_time,extinction_time=extinction_time))
        # choose species to speciate, and add it to tree. 
      }
    }
    cbt = min(next_speciation_time,next_bt)
  }
  tree = data.frame(brts = c(missing_branches$speciation_time,brts,missing_branches$extinction_time),
                    t_ext = c(missing_branches$extinction_time, rep(Inf,length(brts)+nrow(missing_branches))),
                    to = c(rep(1,nrow(missing_branches)),rep(2,length(brts)),rep(0,nrow(missing_branches))))
  tree = tree[order(tree$brts),]
  
  logg = log_sampling_prob_nh(df = tree,pars = pars,model = model)
  
  tree$pd = sapply(tree$brts, function(x)
    emphasis:::phylodiversity(x, tree))
  tree$n = c(1,sapply(tree$brts[-nrow(tree)], function(x)
    emphasis:::number_of_species(tree = tree, tm = x)))

  return(list(tree=tree,logg=logg))
}


#inverse method

tree_augmentation_inverse <- function(brts,pars,model="dd"){
  brts = cumsum(-diff(c(brts,0)))
  b = max(brts)
  mu = pars[3]
  missing_branches = data.frame(speciation_time=NULL,extinction_time=NULL)
  cbt = 0 # current branching time
  while(cbt < b){
    tree = data.frame(brts = c(missing_branches$speciation_time,brts,missing_branches$extinction_time),
                      t_ext = c(missing_branches$extinction_time, rep(Inf,length(brts)+nrow(missing_branches))),
                      to = c(rep(1,nrow(missing_branches)),rep(2,length(brts)),rep(0,nrow(missing_branches))))
    tree = tree[order(tree$brts),]
    next_bt = min(tree$brts[tree$brts>cbt])
    next_speciation_time = rnhpp(time0 = cbt,time_max = next_bt,tree = tree,model = model,pars = pars)
    
    if(next_speciation_time < next_bt){ # add new species
      extinction_time = next_speciation_time + truncdist::rtrunc(1,"exp",a = 0, b = (b-next_speciation_time),rate=mu)
      missing_branches = rbind(missing_branches,data.frame(speciation_time=next_speciation_time,extinction_time=extinction_time))
      cbt = next_speciation_time
    }else{
      cbt = next_bt
    }
  }
  df = data.frame(brts = c(missing_branches$speciation_time,brts,missing_branches$extinction_time),
                  t_ext = c(missing_branches$extinction_time, rep(Inf,length(brts)+nrow(missing_branches))),
                  to = c(rep(1,nrow(missing_branches)),rep(2,length(brts)),rep(0,nrow(missing_branches))))
  df = df[order(df$brts),]
  df$t.ext = df$t_ext-df$brts
  return(df)
  
}

##
rnhpp <- function(time0,time_max,tree,model,pars){  # random non-homogenous exponential (NHPP)
  ex = rexp(1)
  rv = IntInv(ex = ex,time0 = time0,time_max = time_max,tree = tree,model = model,pars = pars)
  return(rv)
}

#DD  inverse method
rnhpp_dd <- function(r,mu,s,cbt,next_bt){  # random non-homogenous exponential (NHPP)
  ex = rexp(1)
  rv = cbt+IntInv_dd(r = r,mu = mu,s = s,u = ex)
  if(is.na(rv) | rv>next_bt){
    rv = Inf
  }
  return(rv)
}

IntInv_dd <- function(r,mu,s,u){
  t = -W(-exp(-r*mu+mu*u/s-exp(-r*mu)))/mu+u/s-exp(-r*mu)/mu
  return(t)
}


intensity <- function(x, tree, model, time0, pars){
  max_time_for_continuity = min(tree$brts[tree$brts>time0])
  if(x==time0){
    val = 0
  }else{
    nh_rate <- function(wt){
      sum_speciation_rate(x=wt,tree = tree,pars = pars,model = model)*(1-exp(-(max(tree$brts)-wt)*pars[3]))
    }
    if(x != time0) x <- x-0.00000000001
    val = pracma:::quad(f = Vectorize(nh_rate),xa = time0,xb = x)
  }
  return(val)
}

inverse = function (f, lower = 0, upper = 100,...) {
  function (y) uniroot((function (x) f(x,...) - y), lower = lower, upper = upper)[1]
}

IntInv <- function(ex,time0,time_max,tree,model,pars){
  if(sign(intensity(x = time0,tree = tree,model = model,time0 = time0,pars = pars)-ex)==sign(intensity(x = time_max, time0 = time0,tree=tree,model=model,pars = pars)-ex)){
    value = Inf
  }else{
    IntInv_inner = inverse(intensity,time0,time_max,tree=tree,model=model,time0=time0,pars=pars)
    inverse = IntInv_inner(ex)
    value = inverse$root
  }
  return(value)
}




tree_augmentation <- function(obs_brts,pars,model,importance_sampler,maxnumspec=NULL,general=TRUE){
  initspec=1
  if(obs_brts[1]==max(obs_brts)){
    wt = -diff(c(obs_brts,0))
    obs_brts = cumsum(wt)
  }
  if(importance_sampler == "emphasis"){
    if(general){
      df = tree_augmentation_inverse(brts = obs_brts,pars = pars,model = model)
    }else{
      if(model=="dd"){
        df = nh_tree_augmentation_dd(brts = obs_brts,pars = pars)
      }
      if(model == "pd"){
        df = nh_tree_augmentation_pd(brts = obs_brts,pars = pars)
      }
    }
    logg = log_sampling_prob_nh(df = df,pars = pars,model = model,initspec = initspec)
  }
  if(importance_sampler == "uniform"){
    df = uniform_tree_augmentation(brts = obs_brts,maxnumspec=maxnumspec)
    logg = log.sampling.prob.uniform(df,maxnumspec=maxnumspec,initspec=initspec,p=0.5)
  }
  logf = loglik.tree(pars=pars,tree=df,model=model,initspec = initspec)
  logweight = logf-logg
  weight = exp(logweight)
  return(list(tree=df,weight=weight,logf=logf,logg=logg))
}


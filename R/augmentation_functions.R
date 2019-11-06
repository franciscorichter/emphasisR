mc_sample_independent_trees <- function(brts,pars,nsim=1000,model="dd",importance_sampler="emphasis",no_cores=2,pars3=NULL,maxnumspec=NULL,seed=0,initspec=1,method="inverse"){
  time=proc.time()
  if(seed>0) set.seed(seed)
  
  if(method == "thinning"){
    E = mc_sample_independent_trees(brts = brts,pars = pars,nsim = nsim, model = model, importance_sampler = importance_sampler,no_cores = no_cores,pars3 = pars3,maxnumspec = maxnumspec,seed = seed,initspec = initspec)
  }else{
  
    #### parallel set-up
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    ##
    if(importance_sampler=="emphasis"){
      trees <- foreach(i = 1:nsim, combine = list) %dopar% {
        df = emphasis::tree_augmentation_inverse(observed.branching.times = brts,pars = pars,model = model)
        if(!is.null(df)){
          logg.samp = emphasis:::log_sampling_prob_nh(df = df,pars = pars,model = model,initspec = initspec)
          logf.joint = emphasis:::loglik.tree(pars=pars,tree=df,model=model,initspec = 1)
          tree.info = list(logf.joint=logf.joint,logg.samp=logg.samp,dim=nrow(df),tree=df)
        }else{
          tree.info = list(logf.joint=0,logg.samp=0,tree=0,dim=0,num.miss=0)
        }
        return(tree.info)
      }
    }
    if(importance_sampler=="uniform"){
      trees <- foreach(i = 1:nsim, combine = list) %dopar% {
        df = emphasis:::sample.uniform(brts,maxnumspec=maxnumspec)
        log.samp.unif.prob = emphasis:::log.sampling.prob.uniform(df,maxnumspec=maxnumspec,initspec=1,p=0.5)
        df$t_exp = rep(Inf,nrow(df)) 
        missing_speciations = NULL
        for(j in 1:nrow(df)){
          if(df$to[j]==1){
            missing_speciations = c(missing_speciations,df$brts[j])
          }
          if(df$to[j]==0){
            sample_ext_time_index = sample(length(missing_speciations),1)
            df$t_exp[j] = missing_speciations[sample_ext_time_index]
            missing_speciations = missing_speciations[-sample_ext_time_index]
          }
        }
      ##################
        logf.joint = emphasis:::loglik.tree(pars=pars,tree=df,model=model,initspec = 1)
        return(list(logf.joint=logf.joint,logg.samp=log.samp.unif.prob,dim=nrow(df),tree=df))
      }
    }
    stopCluster(cl)
    logf = sapply(trees,function(list) list$logf.joint)
    logg = sapply(trees,function(list) list$logg.samp)
    dim = sapply(trees,function(list) list$dim)
    logweights = logf-logg
    weights = exp(logweights)
    fhat = mean(weights)
    trees = lapply(trees, function(list) list$tree)
    time = get.time(time)
    E = list(weights=weights,logweights=logweights,fhat=fhat,logf=logf,fhat.se=1,logg=logg,trees=trees,dim=dim,E_time=time)
  }
  return(E)
}

mc_augmentation <- function(brts,pars,model,importance_sampler,sample_size){
  time = proc.time()
  st =  lapply(1:sample_size,function(i,...){tree_augmentation(...)},obs_brts = brts,pars = pars,model=model,importance_sampler=importance_sampler)
  weights = sapply(st,function(list) list$weight)
  trees = lapply(st,function(list) list$tree)
  dim = sapply(st,function(list) nrow(list$tree))
  logf =    sapply(st,function(list) list$logf)
  logg =    sapply(st,function(list) list$logg)
  E_time = get.time(time)
  En = list(weights=weights,trees=trees,fhat=mean(weights),logf=logf,logg=logg,dim=dim,E_time=E_time)
  return(En)
}

mc_augmentation_thinning <- function(brts,pars,model,importance_sampler,sample_size,parallel=FALSE,no_cores=2){
  time = proc.time()
  if(!parallel){
    st =  lapply(1:sample_size,function(i,...){augment_tree_thinning(...)},obs_brts = brts,pars = pars,model=model)
  }else{
    st = mclapply(1:sample_size,function(i,...){augment_tree_thinning(...)},obs_brts = brts,pars = pars,model=model,mc.cores = no_cores)
  }
  weights = sapply(st,function(list) list$weight)
  trees = lapply(st,function(list) list$tree)
  dim = sapply(st,function(list) nrow(list$tree))
  logf =    sapply(st,function(list) list$logf)
  logg =    sapply(st,function(list) list$logg)
  E_time = get.time(time)
  En = list(weights=weights,trees=trees,fhat=mean(weights),logf=logf,logg=logg,dim=dim,E_time=E_time)
  return(En)
}


## thinning method

augment_tree_thinning <- function(obs_brts,pars,model="dd"){
  mu = pars[3]
  b = max(obs_brts)
  missing_branches = data.frame(speciation_time=NULL,extinction_time=NULL)
  cbt=0
  
  while(cbt < b){
    
    tree = data.frame(brts = c(missing_branches$speciation_time,obs_brts,missing_branches$extinction_time),
                      t_ext = c(missing_branches$extinction_time, rep(Inf,length(obs_brts)+nrow(missing_branches))),
                      to = c(rep(1,nrow(missing_branches)),rep(2,length(obs_brts)),rep(0,nrow(missing_branches))))
    tree = tree[order(tree$brts),]
    
    lambda_max = sum_speciation_rate(cbt,tree,pars,model)*(1-exp(-mu*(b-cbt)))
    next_bt = min(tree$brts[tree$brts>cbt])
    u1 = runif(1)
    next_speciation_time = cbt-log(x = u1)/lambda_max
    if(next_speciation_time < next_bt){
      u2 = runif(1)
      pt = sum_speciation_rate(next_speciation_time,tree,pars,model)*(1-exp(-mu*(b-next_speciation_time)))/lambda_max
      if(u2<pt){
        extinction_time = next_speciation_time + truncdist::rtrunc(1,"exp",a = 0, b = (b-next_speciation_time),rate=mu)
        missing_branches = rbind(missing_branches,data.frame(speciation_time=next_speciation_time,extinction_time=extinction_time))
      }
    }
    cbt = min(next_speciation_time,next_bt)
  }
  tree = data.frame(brts = c(missing_branches$speciation_time,obs_brts,missing_branches$extinction_time),
                    t_ext = c(missing_branches$extinction_time, rep(Inf,length(obs_brts)+nrow(missing_branches))),
                    to = c(rep(1,nrow(missing_branches)),rep(2,length(obs_brts)),rep(0,nrow(missing_branches))))
  tree = tree[order(tree$brts),]
  
  logg = log_sampling_prob_nh(df = tree,pars = pars,model = model)
  logf = loglik.tree(pars=pars,tree=tree,model=model)
  logweight = logf-logg
  weight = exp(logweight)
  return(list(tree=tree,weight=weight,logf=logf,logg=logg))
}


#inverse method

tree_augmentation_inverse <- function(observed.branching.times,pars,model="dd",initspec = 1){
  
  b = max(observed.branching.times)
  mu = pars[3]
  missing_branches = data.frame(speciation_time=NULL,extinction_time=NULL)
  cbt = 0 # current branching time
  while(cbt < b){
    tree = data.frame(brts = c(missing_branches$speciation_time,observed.branching.times,missing_branches$extinction_time),
                      t_ext = c(missing_branches$extinction_time, rep(Inf,length(observed.branching.times)+nrow(missing_branches))),
                      to = c(rep(1,nrow(missing_branches)),rep(2,length(observed.branching.times)),rep(0,nrow(missing_branches))))
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
  df = data.frame(brts = c(missing_branches$speciation_time,observed.branching.times,missing_branches$extinction_time),
                  t_ext = c(missing_branches$extinction_time, rep(Inf,length(observed.branching.times)+nrow(missing_branches))),
                  to = c(rep(1,nrow(missing_branches)),rep(2,length(observed.branching.times)),rep(0,nrow(missing_branches))))
  df = df[order(df$brts),]
  df$t.ext = df$t_ext-df$brts
  return(df)
  
}

tree_augmentation <- function(obs_brts,pars,model,importance_sampler,maxnumspec=NULL,general=TRUE){
  initspec=1
  if(obs_brts[1]==max(obs_brts)){
    wt = -diff(c(obs_brts,0))
    obs_brts = cumsum(wt)
  }
  if(importance_sampler == "emphasis"){
    if(general){
      df = tree_augmentation_inverse(observed.branching.times = obs_brts,pars = pars,model = model)
    }else{
      if(model=="dd"){
        df = nh_tree_augmentation_dd(observed.branching.times = obs_brts,pars = pars)
      }
      if(model == "pd"){
        df = nh_tree_augmentation_pd(observed.branching.times = obs_brts,pars = pars)
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


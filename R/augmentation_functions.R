mc_augmentation <- function(brts,pars,model,importance_sampler,sample_size,parallel=FALSE,no_cores=2,soc){
  time = proc.time()
  if(!parallel){
    st =  lapply(1:sample_size,function(i){augment_tree(brts = brts,pars = pars,model=model,soc=soc)} )
  }else{
    st = mclapply(1:sample_size,function(i){augment_tree(brts = brts,pars = pars,model=model,soc=soc)},mc.cores = no_cores)
  }
  trees = lapply(st,function(list) list$tree)
  dim = sapply(st,function(list) nrow(list$tree))
  
  ####
  logf = sapply(trees,loglik.tree(model), pars=pars)
  logg = sapply(trees,sampling_prob, pars=pars,model=model,soc=soc)
  E_time = get.time(time)
 #logg = log_sampling_prob_nh(df = tree,pars = pars,model = model,soc=soc)
  
 # logg =    sapply(st,function(list) list$logg)
  log_weights = logf-logg
  w = exp(log_weights)
  ####
  En = list(weights=w,trees=trees,fhat=mean(w),logf=logf,logg=logg,dim=dim,E_time=E_time)
  return(En)
}

augment_tree <- function(brts,pars,model="dd",soc){
  mu = max(0,pars[1])
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
    lambda_max = max(sum_speciation_rate(cbt,tree,pars,model,soc=soc)*(1-exp(-mu*(b-cbt))),sum_speciation_rate(next_bt,tree,pars,model,soc=soc)*(1-exp(-mu*(b-next_bt))))
    u1 = runif(1)
    next_speciation_time = cbt-log(x = u1)/lambda_max
    if(next_speciation_time < next_bt){
      u2 = runif(1)
      pt = sum_speciation_rate(next_speciation_time,tree,pars,model,soc=soc)*(1-exp(-mu*(b-next_speciation_time)))/lambda_max
      if(u2<pt){
        extinction_time = next_speciation_time + truncdist::rtrunc(1,"exp",a = 0, b = (b-next_speciation_time),rate=mu)
        missing_branches = rbind(missing_branches,data.frame(speciation_time=next_speciation_time,extinction_time=extinction_time))
      }
    }
    cbt = min(next_speciation_time,next_bt)
  }
  tree = data.frame(brts = c(missing_branches$speciation_time,brts,missing_branches$extinction_time),
                    t_ext = c(missing_branches$extinction_time, rep(Inf,length(brts)+nrow(missing_branches))),
                    to = c(rep(1,nrow(missing_branches)),rep(2,length(brts)),rep(0,nrow(missing_branches))))
  tree = tree[order(tree$brts),]
  
  #logg = log_sampling_prob_nh(df = tree,pars = pars,model = model,soc=soc)
  
  tree$pd = sapply(tree$brts, function(x)
  emphasis:::phylodiversity(x, tree,soc=soc))
  
  ## check this one
  
  tree$n = sapply(c(0,tree$brts[-length(tree$brts)]), n_from_time,tree=tree,soc=soc)
  
  return(list(tree=tree))
}

##############################################





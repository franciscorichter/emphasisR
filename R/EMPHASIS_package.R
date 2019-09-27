### EMPHASIS functions

# negative logLikelihood of a tree
nllik.tree = function(pars,tree,model="dd",initspec=1){
  nl = -loglik.tree(pars,tree,model,initspec)
  return(nl)
}


loglik.tree <- function(pars,tree,model,initspec=1){
  if(model == "dd"){
    log.lik = loglik.tree.dd(pars,tree,model,initspec)
  }
  if(model == "edd"){
    log.lik = loglik.tree.edd(pars,tree,model,initspec)
  }
  if(model == "pd"){ 
    log.lik = loglik.tree.pd(pars,tree,model,initspec)
  }
  if(model == "epd"){ 
    log.lik = loglik.tree.epd(pars,tree,model,initspec)
  }
  if(model == "gddx"){ 
    if(min(pars)<0){
      log.lik = -Inf
    }else{
      log.lik = loglik.tree.gddx(pars,tree,model,initspec)
    }
  }
  if(model == "gpdx"){ 
    if(min(pars)<0){
      log.lik = -Inf
    }else{
      log.lik = loglik.tree.gpdx(pars,tree,model,initspec)
    }
  }
  return(log.lik)
}

loglik.tree.dd <- function(pars,tree,model,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  lambda = lambda.dd(pars,n)
  sigma = (lambda + mu)*n
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = (sum(-sigma*wt)+sum(log(rho)))
  if(min(pars)<0) log.lik = -Inf
  return(log.lik)
}

loglik.tree.edd <- function(pars,tree,model,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  lambda = lambda.edd(pars,n)
  sigma = (lambda + mu)*n
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = (sum(-sigma*wt)+sum(log(rho)))
  return(log.lik)
}

loglik.tree.gddx <- function(pars,tree,model,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[2])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  lambda = lambda.gddx(pars,n)
  sigma = (lambda + mu)*n
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = (sum(-sigma*wt)+sum(log(rho)))
  return(log.lik)
}


loglik.tree.pd <- function(pars,tree,model,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  sigma_over_tree = lambda = NULL
  bt = 0
  for(i in 1:length(wt)){
    sigma_over_tree[i] = n[i]*(pars[3]*wt[i]+quad(lambda.pd_t,xa=bt,xb=bt+wt[i],pars=pars,tree=tree))
    bt = bt + wt[i]
    lambda[i] = lambda.pd_t(bt,pars,tree)
  }
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}

lambda.pd_t <- function(time_m,pars,tree){
  lambda = max(0,pars[1] + pars[2] * phylodiversity_t(time_m = time_m,tree = tree))
  return(lambda)
}

loglik.tree.epd <- function(pars,tree,model,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[3])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  n.obs = cumsum(c(1,tree$to==2))
  n.obs = head(n.obs,-1)
  pd_obs = cumsum(n.obs*wt)
  pd_miss = n.miss_vec = ext.miss_vec = NULL
  for(j in 1:nrow(tree)){
    if(is.null(n.miss_vec)){
      pd_miss[j] = 0
    }else{
      pd_miss[j] = sum(tree$brts[j] - n.miss_vec)
    }
    if(tree$to[j]==1){
      n.miss_vec = c(n.miss_vec,tree$brts[j])
      ext.miss_vec = c(ext.miss_vec,tree$t_exp[j])
    }
    if(tree$to[j]==0){
      n.miss_vec = n.miss_vec[ext.miss_vec != tree$brts[j]]
    }
  }
  pd = pd_obs + pd_miss
  lambda = exp(pars[1] + pars[2]*pd)
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  sigma_over_tree = NULL
  for(i in 1:length(wt)){
    sigma_over_tree[i] = n[i]*(pars[3]*wt[i]+integrate(lambda.epd_t,lower=0,upper=wt[i],pars=pars,pd[i],n=n[i],wt[i])$value)
  }
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}



loglik.tree.gpdx <- function(pars,tree,model,initspec=1){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  mu = max(0,pars[2])
  wt = diff(c(0,tree$brts))
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  n.obs = cumsum(c(1,tree$to==2))
  n.obs = head(n.obs,-1)
  pd_obs = cumsum(n.obs*wt)
  pd_miss = n.miss_vec = ext.miss_vec = NULL
  for(j in 1:nrow(tree)){
    if(is.null(n.miss_vec)){
      pd_miss[j] = 0
    }else{
      pd_miss[j] = sum(tree$brts[j] - n.miss_vec)
    }
    if(tree$to[j]==1){
      n.miss_vec = c(n.miss_vec,tree$brts[j])
      ext.miss_vec = c(ext.miss_vec,tree$t_exp[j])
    }
    if(tree$to[j]==0){
      n.miss_vec = n.miss_vec[ext.miss_vec != tree$brts[j]]
    }
  }
  pd = pd_obs + pd_miss
  lambda = pmax(0,pars[1]*((pd+1)^(-pars[2])))
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  sigma_over_tree = NULL
  for(i in 1:length(wt)){
    sigma_over_tree[i] = n[i]*(pars[2]*wt[i]+integrate(lambda.dpdx_t,lower=0,upper=wt[i],pars=pars,pd=pd[i],n=n[i],wt=wt[i])$value)
  }
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}

speciation_rate <- function(pars,N,model){
  if(model == "dd"){  # diversity-dependence model
    lambda = lambda.dd(pars,N)
  }
  if(model == "dd1.3"){
    lambda = lambda.dd.1.3(pars,N)
  }
  if(model == "edd"){
    lambda = lambda.edd(pars,N)
  }
  return(lambda)
}

lambda.dd <- function(pars,n){
  lambdas =  pmax(0, pars[1] - pars[2]*n)
  return(lambdas)
}

lambda.edd <- function(pars,n,GLM=FALSE){
  lambdas =  exp(pars[1] + pars[2]*n)
  return(lambdas)
}

lambda.gddx <- function(pars,n,GLM=FALSE){
  lambdas =  pmax(pars[1]*(n^(-pars[2])),0)
  return(lambdas)
}

lambda.dpdx_t <- function(w_time,pars,pd,n,wt){
  lambda = pmax(pars[1]*((pd - n*(wt-w_time)+1)^(-pars[3])),0)
  return(lambda)
}

lambda.epd_t <- function(w_time,pars,pd,n,wt){
  lambda = exp(pars[1]+pars[2]*(pd - n*(wt-w_time)))
  return(lambda)
}

###############################################################
# Uniform data augmentation importance sampler#
###############################################################

lprobto <- function(to,p=0.5){
  posspec = c(0,cumsum(to==1))<(length(to)/2)
  posext = !(c(0,cumsum(to==1))==c(0,cumsum(to==0)))
  possibletotal = posspec & posext
  to_possible = to[possibletotal]
  logprob = sum(to_possible==1)*log(p)+sum(to_possible==0)*log(1-p)
  return(logprob)
}

sampletopology <- function(S,p=0.5){
  to = NULL
  if(S>0){
    for(i in 1:(2*S)){
      if(sum(to==1)==sum(to==0)){
        prob = 1
      }
      if(sum(to==1)==S){
        prob = 0
      }
      to = c(to,rbinom(n=1,size=1,prob=prob))
      prob = p
    }
  }else{
    to = NULL
  }
  return(to)
}

sample_dim_prob <- function(d,max){
  factorial(2*d)/sum(factorial(seq(from = 0,to = 2*max,by=2)))
}

sample.uniform <- function(brts,maxnumspec,single_dimension=NULL){
  if(is.null(single_dimension)){
    probs = sample_dim_prob(0:maxnumspec,maxnumspec)
    S = sample(0:maxnumspec,1)
  }else{
    S = single_dimension
  }
  mbts.events = emphasis:::sim.branchingtimes.and.events(S=S ,ct = max(brts),p=0.5)
  df = data.frame(brts=c(brts,mbts.events$brts),to=c(rep(2,length(brts)),mbts.events$to))
  df = df[order(df$brts),]
  return(df)
}

sim.branchingtimes.and.events <- function(S=S,ct,p){
  brts = sort(runif(2*S,min=0,max=ct))
  to = sampletopology(S,p = p)
  tree = list(brts=brts,to=to)
  return(tree)
}

log.factor.samp.prob <- function(to){
  top = head(to,-1)
  number.observed = c(1,1+cumsum(top==2))
  number.missing = c(0,cumsum(top==1)-cumsum(top==0))
  factor = -sum(log((2*number.observed+number.missing)[to==1]))-sum(log(number.missing[to==0]))
  return(factor)
}

log.sampling.prob.uniform <- function(df,maxnumspec,initspec=1,p=0.5){
  ct = max(df$brts)
  to = top = df$to
  tom = top[top!=2]
  to[to==2] = 1
  num.miss = 2*sum((to==0))
  loggprob <- -log((maxnumspec+1))+lgamma(num.miss+1)-num.miss*log(ct)+lprobto(tom,p = p)+log.factor.samp.prob(top) 
  return(loggprob)
}

log.samp.prob <- function(to,maxnumspec,ct,initspec=1,conf,p){
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  n = n[-length(n)]
  loggprob <- -log((maxnumspec+1))+lgamma(length(to)+1)-length(to)*log(ct)+lprobto(to,p = p)-sum(log(conf$N-conf$P))
}


#################################
## emphasis sampler
rnhpp <- function(s,mu,r){  # random non-homogenous exponential (NHPP)
  ex = rexp(1)
  rv = IntInv(r=r,mu=mu,s=s,u=ex)
  if(is.na(rv)){
    rv = Inf
  }
  return(rv)
}

IntInv <- function(r,mu,s,u){
  t = -W(-exp(-r*mu+mu*u/s-exp(-r*mu)))/mu+u/s-exp(-r*mu)/mu
  return(t)
}

time_limit_dd <- function(obs_brts,missing_speciations,missing_extinctions,max_num_spec){
 df = data.frame(bt = c(obs_brts,missing_speciations,missing_extinctions),to = c(rep(1,length(obs_brts)+length(missing_speciations)),rep(-1,length(missing_extinctions)))) 
 df = df[order(df$bt),]
 N = c(1,cumsum(df$to))
 b = tail(df$bt[N <= max_num_spec],1)
 return(b)
}


nh_tree_augmentation <- function(observed.branching.times,pars,model="dd",initspec = 1){
  
  b = max(observed.branching.times)
  mu = pars[3]
  missing_branches = data.frame(speciation_time=NULL,extinction_time=NULL)
  N = initspec # current number of species
  cbt = 0 # current branching time
  while(cbt < b){
    lambda = speciation_rate(pars,N,model)
    next_speciation_time = cbt+rnhpp(s=N*lambda,mu=mu,r=b-cbt)
    all_bt = c(observed.branching.times,missing_branches$speciation_time,missing_branches$extinction_time)
    next_bt = min(all_bt[all_bt>cbt])
    if(next_speciation_time < next_bt){ # add new species
      extinction_time = next_speciation_time + truncdist::rtrunc(1,"exp",a = 0, b = (b-next_speciation_time),rate=mu)
      missing_branches = rbind(missing_branches,data.frame(speciation_time=next_speciation_time,extinction_time=extinction_time))
      cbt = next_speciation_time
      N = N+1
    }else{
      cbt = next_bt
      if(next_bt %in% missing_branches$extinction_time){
        N = N-1
      }else{
        N = N+1
      }
    }
  }
  df = data.frame(brts = c(missing_branches$speciation_time,observed.branching.times,missing_branches$extinction_time),
                  bte = c(missing_branches$extinction_time, rep(Inf,length(observed.branching.times)+nrow(missing_branches))),
                  to = c(rep(1,nrow(missing_branches)),rep(2,length(observed.branching.times)),rep(0,nrow(missing_branches))))
  df = df[order(df$brts),]
  df$t.ext = df$bte-df$brts
  return(df)
  
}



log_sampling_prob_nh <- function(df,pars,model="dd",initspec=1){
  b = max(df$brts)
  to = top = head(df$to,-1)
  to[to==2] = 1
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  lambda = lambda.dd(pars,n)
  brts_i = df$brts
  brts_im1 = c(0,df$brts[-nrow(df)])
  # at missing speciation times
  missing_speciations = (df$to == 1)
  nb = n[missing_speciations]
  No = c(1,1+cumsum(top==2))[missing_speciations]
  Ne = c(0,cumsum(top==1)-cumsum(top==0))[missing_speciations]
  lambda_b = lambda[missing_speciations]
  text = df$t.ext[missing_speciations]
  mu = pars[3]
  logg = sum(-n * lambda * (brts_i-brts_im1-(exp(-b*mu)/mu) *  (exp(brts_i*mu)-exp(brts_im1*mu))  ))+sum(log(nb)+log(mu))-sum(mu*text)+sum(log(lambda_b))-sum(log(2*No+Ne))

  return(logg)
}


##################################
###### E step 


mc_sample_independent_trees <- function(brts,pars,nsim=1000,model="dd",importance_sampler="emphasis",no_cores=2,pars3=NULL,maxnumspec=NULL,seed=0,initspec=1,limit_on_species=NULL){
  time=proc.time()
  if(seed>0) set.seed(seed)
  if(brts[1] != max(brts)){
    stop("please add brts on descending order")
  }
  wt = -diff(c(brts,0))
  brts = cumsum(wt)  # Do I need this>>
  if(is.null(pars3)) pars3=pars
  #### parallel set-up
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  ##
  if(importance_sampler=="emphasis"){
    trees <- foreach(i = 1:nsim, combine = list) %dopar% {
      df = emphasis::nh_tree_augmentation(brts,pars = pars,model = model)
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
  return(E)
}



##############################
####### M-step 

Q_approx = function(pars,st,model="dd",initspec=1){
  get_llik <- function(tree) nllik.tree(pars=pars,tree=tree,initspec = initspec,model=model)
  l = sapply(st$trees, get_llik)
  w = st$weights/(sum(st$weights))
  Q = sum(l*w)
  return(Q)
}

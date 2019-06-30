### EMPHASIS functions

# negative logLikelihood of a tree
nllik.tree = function(pars,tree,model="dd",initspec=1,K=FALSE){
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
  lambda = lambda.dd(pars,n,GLM=TRUE)
  sigma = (lambda + mu)*n
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  log.lik = (sum(-sigma*wt)+sum(log(rho)))
 # if(min(pars)<0){nl = Inf}
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
  lambda = pmax(0,pars[1]*((pd+1)^(-pars[3])))
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  sigma_over_tree = NULL
  for(i in 1:length(wt)){
    sigma_over_tree[i] = n[i]*(pars[2]*wt[i]+integrate(lambda.dpdx_t,lower=0,upper=wt[i],pars=pars,pd=pd[i],n=n[i],wt=wt[i])$value)
  }
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}



###############################################################
# Uniform data augmentation importance sampler (Bart version) #
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
    S = sample(1:maxnumspec,1)
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
  loggprob <- -log((maxnumspec))+lgamma(num.miss+1)-num.miss*log(ct)+lprobto(tom,p = p)+log.factor.samp.prob(top) 
  return(loggprob)
}

log.samp.prob <- function(to,maxnumspec,ct,initspec=1,conf,p){
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  n = n[-length(n)]
  loggprob <- -log((maxnumspec+1))+lgamma(length(to)+1)-length(to)*log(ct)+lprobto(to,p = p)-sum(log(conf$N-conf$P))
}

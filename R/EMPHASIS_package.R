### EMPHASIS functions




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

uniform_tree_augmentation <- function(brts,maxnumspec){
  S = sample(0:maxnumspec,1)
  mbts.events = emphasis:::sim.branchingtimes.and.events(S=S ,ct = max(brts),p=0.5)
  df = data.frame(brts=c(brts,mbts.events$brts),to=c(rep(2,length(brts)),mbts.events$to))
  df = df[order(df$brts),]
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
  if(x > max_time_for_continuity){
    stop("max_time_for_continuity_reached")
  }
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


### PD models

time_limit_dd <- function(obs_brts,missing_speciations,missing_extinctions,max_num_spec){
 df = data.frame(bt = c(obs_brts,missing_speciations,missing_extinctions),to = c(rep(1,length(obs_brts)+length(missing_speciations)),rep(-1,length(missing_extinctions)))) 
 df = df[order(df$bt),]
 N = c(1,cumsum(df$to))
 b = tail(df$bt[N <= max_num_spec],1)
 return(b)
}


##############################
####### M-step 

M_step <-function(st,init_par = NULL,model="dd",exclude_proportion_trees = 0){
  
  time0 = proc.time()
  
  w = st$weights/max(st$weights)
  contributing_trees = (w > exclude_proportion_trees)
  ######################
  sub_trees = st$trees[contributing_trees]
  effective_sample_size = sum(contributing_trees)
  sub_st = list(trees = sub_trees, weights = st$weights[contributing_trees])
  loglik = get(paste0("loglik.tree.", model))
  po = subplex(par = init_par, fn = Q_approx,st = sub_st,loglik=loglik,hessian = TRUE)
  #######################

  M_time = get.time(time0)
  return(list(po=po,M_time=M_time,effective_sample_size=effective_sample_size))
}

Q_approx = function(pars,st,loglik){
  l = sapply(st$trees, loglik, pars=pars)
  w = st$weights
  Q = -sum(l*w)
  return(Q)
}


## parallel M step 

#library("optimParallel")

M_step_parallel <-function(st,init_par = NULL,model="dd",exclude_proportion_trees = 0,no_cores=2){
  
  time0 = proc.time()
  
  w = st$weights/max(st$weights)
  contributing_trees = (w > exclude_proportion_trees)
  ######################
  sub_trees = st$trees[contributing_trees]
  effective_sample_size = sum(contributing_trees)
  sub_st = list(trees = sub_trees, weights = st$weights[contributing_trees])
  loglik = get(paste0("loglik.tree.", model))
  
  cl <- makeCluster(5)     # set the number of processor cores
  setDefaultCluster(cl=cl) # set 'cl' as default cluster
  
  optimParallel(par=init_par, fn=Q_approx, st = sub_st,loglik=loglik,
                method = "L-BFGS-B", lower=c(0, 0, 0))
  
  setDefaultCluster(cl=NULL); stopCluster(cl)
  #######################
  
  M_time = get.time(time0)
  return(list(po=po,M_time=M_time,effective_sample_size=effective_sample_size))
}





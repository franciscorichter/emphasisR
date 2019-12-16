### EMPHASIS functions

mcem.tree <- function(input,max_iterations=100000,report=TRUE,file=NULL){
  pars = PARS = input$pars
  MCEM=data.frame(loglik_hat=NULL,E_time=NULL,M_time=NULL,sample_size=NULL)
  sample_size = input$sample_size
  #prev_lg = -99999
  for(i in 1:max_iterations){
    
    st = mc_sample_independent_trees(brts = input$brts,pars = pars,nsim = sample_size,model = input$model, importance_sampler = input$importance_sampler,no_cores = input$cores, method = input$method)
    lg = log(st$fhat)
    if(report){ 
      print(paste(c("parameters: ",pars)))
      print(paste("loglikelihood: ",lg))
      print(paste("sample size: ",sample_size))
    }
    M = M_step(st = st,init_par = pars,model = input$model)
    if(!is.na(M$po$value)){
      pars = M$po$par
    }
    print(M$po$value)
    PARS = rbind(pars,PARS)
    MCEM = rbind(MCEM,data.frame(loglik_hat=log(st$fhat),E_time=st$E_time,M_time=M$M_time,sample_size=sample_size))
    #if(prev_lg > lg) sample_size = sample_size*input$aceleration_rate
    #prev_lg = lg
    save(MCEM,input,PARS,file=file)
  }
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
  
  sub_st = get_contributing_trees(st)
  loglik = get(paste0("loglik.tree.", model))
  po = subplex(par = init_par, fn = Q_approx,st = sub_st, loglik=loglik, hessian = TRUE)

  M_time = get.time(time0)
  return(list(po=po,M_time=M_time))
}

get_contributing_trees <- function(st, exclude_proportion_trees=0){
  w = st$weights/sum(st$weights)
  contributing_trees = (w > exclude_proportion_trees)
  sub_trees = st$trees[contributing_trees]
  effective_sample_size = sum(contributing_trees)
  sub_st = list(trees = sub_trees, weights = st$weights[contributing_trees])
  return(sub_st)
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





### EMPHASIS functions

# negative logLikelihood of a tree
nllik.tree = function(pars,tree,model="dd",initspec=1){
  if(is.data.frame(tree)){
    tree = df2tree2(df=tree,pars=pars,model=model,initspec=initspec)
  }
  wt = tree$wt
  to = tree$to
  to[to==2] = 1
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  if(model == "cr"){
    lambda = lambda.cr(pars,n)
  }
  if(model == "dd"){
    lambda = lambda.dd(pars,n)
  }
  if(model == "dd.1.3"){
    lambda = lamda.dd.1.3(pars,n)
  }
  mu = max(0,pars[2])
  sigma = (lambda + mu)*n
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  nl = -(sum(-sigma*wt)+sum(log(rho)))#-sum(to==1)*log(2))
  if(min(pars)<0){nl = Inf}
  return(nl)
}



lambda.dd.1.3 <- function(pars,n){
  pmax(1e-99, pars[1]*(1-n/pars[3]))
}

lambda.dd <- function(pars,n){
  pmax(0, (pars[1]-(pars[1]-pars[2])*(n/pars[3])))
}

lambda.cr <- function(pars,n){
  rep(pmax(1e-99, pars[1]), length(n))
}

Q_approx = function(pars,st,model="dd",initspec=1){
  get_llik <- function(tree) nllik.tree(pars=pars,tree=tree,initspec = initspec,model=model)
  l = sapply(st$trees, get_llik)
  w = st$w/(sum(st$w))
  Q = sum(l*w)
  return(Q)
}

Q_SAEM = function(sample,previous_Q,gamma){
  
}


M_step <-function(S,init_par = NULL,model="dd"){
  po = subplex(par = init_par, fn = Q_approx,st = S,model=model,hessian = TRUE)
  return(po)
}

##to do: M-step parallel

df2tree <- function(df,pars,model="dd",initspec=1){
  dim = dim(df)[1]
  wt = diff(c(0,df$bt))
  to = df$to[-dim]
  to[to==2] = 1
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  if(model=="dd"){
    s = n*lambda.dd(pars,n)
  }
  if(model == "dd.1.3"){
    s = n*lamda.dd.1.3(pars,n)
  }
  if(model=="cr"){
    s = n*lambda.cr(pars,n)
  }
  r = df$r
  t_ext = df$t.ext
  return(list(wt=wt,to=to,n=n,s=s,r=r,pars=pars,t_ext=t_ext,protected=df$protected,treeaug=df$to))
}


###########################

###  simulation of extinct species

####

mcem_step <- function(brts,theta_0,MC_ss=10,maxnumspec=NULL,model="dd",selectBestTrees=FALSE,bestTrees=NULL,no_cores,method="emphasis",p=0.5,parallel=TRUE){
  time = proc.time()
  st = E_step(brts = brts,pars = theta_0,nsim = MC_ss,model = model,method = method,no_cores = no_cores,maxnumspec = maxnumspec,p=p,parallel=parallel)
  fhat = st$fhat
  se = st$fhat.se
  E_time = get.time(time)
  time = proc.time()
  if(selectBestTrees){
    weights = st$weights/sum(st$weights)
    max.weight = sort(weights,decreasing = TRUE)[bestTrees]
    sub_st = lapply(st, "[", weights>=max.weight)
    loglik.proportion = sum(weights[weights>=max.weight])
  }else{
    sub_st = lapply(st, "[", st$weights!=0)
    loglik.proportion = sum(st$weights[st$weights>0])
  }
  M = M_step(S = sub_st,init_par = theta_0,model = model)
  M_time = get.time(time)
  pars = M$par
  hessian_inverse = try(diag(solve(M$hessian)))
  if(!is.numeric(h1)) h1 = c(NULL,NULL,NULL)
  return(list(pars=pars,fhat=fhat,se=se,st=st,loglik.proportion=loglik.proportion,hessian_inverse=hessian_inverse,E_time=E_time,M_time=M_time))
}


####### Monte Carlo E step

mc.Estep_parallel <- function(brts,pars,pars3=NULL,nsim=1000,model="dd",method="emphasis",no_cores=NULL,maxnumspec=NULL,p=0.5,seed=0,initspec=1,single_dimension=NULL,limit_on_species=NULL){
  time=proc.time()
  if(seed>0) set.seed(seed)
  #we need brts on ascending order
  if(brts[1]==max(brts)){
    wt = -diff(c(brts,0))
    brts = cumsum(wt)
  }
  if(is.null(pars3)) pars3=pars
  #### parallel set-up
  if(is.null(no_cores)) no_cores <- detectCores()
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  if(method=="emphasis"){
    trees <- foreach(i = 1:nsim, combine = list) %dopar% {
      ## check for dendroica_norm with pars3=2,1,30
      df = emphasis::augment.tree(brts,pars = pars3,model = model,limit_on_species=limit_on_species)
      if(!is.null(df)){
        tree = emphasis::df2tree(df,pars3,model=model,initspec=1)
        logg.samp = emphasis:::log_sampling_prob_emphasis(tree = tree,pars = pars3,model = model,initspec = initspec)
        logf.joint = -emphasis:::nllik.tree(pars=pars,tree=tree,model=model,initspec = 1)
        dim = nrow(df)
        num.miss = sum(df$to==0)
        tree.info = list(logf.joint=logf.joint,logg.samp=logg.samp,tree=tree,dim=dim,num.miss=num.miss)
      }else{
        tree.info = list(logf.joint=0,logg.samp=0,tree=0,dim=0,num.miss=0)
      }
      
      return(tree.info)
    }
  }
  if(method=="emphasis2"){
    #obstree = emphasis::bt2tree(brts)
    trees <- foreach(i = 1:nsim, combine = list) %dopar% {
      df = emphasis::augment.tree2(observed.branching.times = brts,pars = pars3,model = model)
      tree = emphasis::df2tree(df,pars3,model=model,initspec=1)
      logg.samp = emphasis:::log_sampling_prob_emphasis2(tree = tree,pars = pars3,model = model,initspec = initspec)
      logf.joint = -emphasis:::nllik.tree(pars=pars,tree=tree,model=model,initspec = 1)
      dim = nrow(df)
      num.miss = sum(df$to==0)
      return(list(logf.joint=logf.joint,logg.samp=logg.samp,tree=tree,dim=dim,num.miss=num.miss))
    }
  }
    if(method=="uniform"){
      trees <- foreach(i = 1:nsim, combine = list) %dopar% {
        df = emphasis:::sample.uniform(brts,maxnumspec=maxnumspec,single_dimension = single_dimension)
        if(!is.null(single_dimension)) maxnumspec=0
        log.samp.unif.prob = emphasis:::log.sampling.prob.uniform(df,maxnumspec=maxnumspec,initspec=1,p=0.5)
        ####df2tree######
        wt = diff(c(0,df$brts))
        to = df$to
        to = head(to,-1)
        tree=list(wt=wt,to=to)
        ##################
        logf.joint = -emphasis:::nllik.tree(pars=pars,tree=tree,model=model,initspec = 1)
        return(list(logf.joint=logf.joint,logg.samp=log.samp.unif.prob,dim=nrow(df),tree=tree))
      }
    }
  stopCluster(cl)
  logf = sapply(trees,function(list) list$logf.joint)
  logg = sapply(trees,function(list) list$logg.samp)
  dim = sapply(trees,function(list) list$dim)
  diff_logs = logf-logg
  max_log = max(diff_logs) #
  fhat = mean(exp(diff_logs))
  se = sd(exp(diff_logs))/sqrt(nsim)
  weights = exp(diff_logs)#-max_log)
  logweights = diff_logs
  trees = lapply(trees, function(list) list$tree)
  time = get.time(time)
  E = list(weights=weights,logweights=logweights,fhat=fhat,fhat.se=se,logf=logf,logg=logg,trees=trees,dim=dim,Etime=time)
  pw = proper_weighting(E)
  fhat2 = log(mean(pw$W))
  return(c(E,list(fhat2=fhat2)))
}


E_step <- function(brts,pars,nsim=1000,model="dd",method="emphasis",no_cores=2,maxnumspec=NULL,p=0.5,seed=0,parallel=TRUE){
  if(parallel){
    E = mc.Estep_parallel(brts,pars,nsim=nsim,model=model,method=method,no_cores=no_cores,maxnumspec=maxnumspec,p=p,seed=seed)
  }else{
    print("non parallel Em not ready")
  #  E = mc.Estep(brts=brts,pars=pars,nsim=nsim,model=model,method=method,maxnumspec = maxnumspec,p=p,seed=seed)
  }
  return(E)
}


############################
# emphasis data augmentation importance sampler
##############################
bt2tree <- function(brts){
  if(length(brts)>1){
    list(brts=brts,df=data.frame(parent=rep("s1",length(brts)-1),child=paste("s",2:length(brts),sep="")))
  }
  if(length(brts==1)){
    list(brts=brts,df=NULL)
  }
}

rnhe <- function(lambda,mu,r,extinct=TRUE){  # random non-homogenous exponential
  ex = rexp(1)
  rv = IntInv(r=r,mu=mu,s=lambda,u=ex)
  if(is.na(rv)){
    rv = Inf
  }
  return(rv)
}

IntInv <- function(r,mu,s,u){
  t = -W(-exp(-r*mu+mu*u/s-exp(-r*mu)))/mu+u/s-exp(-r*mu)/mu
  return(t)
}


augment.tree <- function(brts,pars,model='dd',initspec=1,seed=0,limit_on_species=NULL){
  if(seed>0) set.seed(seed)
  if(brts[1]==max(brts)){
    wt = -diff(c(brts,0))
    brts = cumsum(wt)
  }else{
    wt = diff(c(0,brts))
  }
  ct = sum(wt)
  dim = length(wt)
  num.of.branches = length(brts)+1
  mu = pars[2]
  kprima = (pars[1]*pars[3])/(pars[1]-pars[2])
  bt = NULL
  bte = NULL
  to = NULL
  N = initspec
  df=1
  
  #set of protected species
  protect = 1
  
  #set of species that are going to die
  extinct = NULL
  
  #limits on extinction times
  r = NULL
  
  #number of times we need to flip a coin to know which species dies
  protected = NULL

  for(i in 1:dim){
    cwt = wt[i]
    cbt = sum(wt[0:(i-1)])
    key = 0
    if(!is.null(limit_on_species)){
      if(length(bt) >= limit_on_species){
        df = NULL
        break
      }
    }
    while(key == 0 & !is.null(df)){
      if(model == "dd"){  # diversity-dependence model
        lambda = lambda.dd(pars,N)
      }
      if(model == "dd1.3"){
        lambda = lambda.dd.1.3(pars,N)
      }
      if(model=="cr"){
        lambda = lambda.cr(pars,N)
      }
      s = N*lambda
      sbte = bte[bte>cbt]
      b = get.max.speciation.bt(sbt = brts[brts>cbt],ebt = sbte,N = N,max.simultaneous.spec = floor(kprima))-cbt
      r = c(r,b)
      t.spe = rnhe(lambda=s,mu=mu,r=b)
      t_ext = ifelse(length(sbte)>0,min(sbte),Inf) - cbt
      mint = min(t.spe,t_ext)
      if(mint < cwt & !is.null(df)){
        if(mint == t.spe){#speciation
          bt = c(bt,cbt+t.spe)
          parent_spec = sample(c(protect,extinct),1)
          if(parent_spec %in% protect){
            # sample between parent_spec and new species
            protected = c(protected,1)
          }else{
            protected = c(protected,0)
          }
          extinct = c(extinct,num.of.branches)
          num.of.branches = num.of.branches + 1
          
          # if is protected there is 1/2, save the protected and dash the dead one 
          #if is extincted, both go to extinct.
          cbt = cbt + t.spe
          text = cbt + truncdist::rtrunc(1,"exp",a = 0, b = (b-t.spe),rate=mu)
       #   print(paste("saving extinction at: ",text))
          #t_ext_limits = c(t_ext_limits,b)
          bte = c(bte,text)
          to = c(to,1)
          N = N + 1
      #    print(paste('numberof species: ', N))
          cwt = cwt - t.spe
          if(!is.null(limit_on_species)){
            if(length(bt) >= limit_on_species){
              df = NULL
              break
            }
          }
        }else{#extinction
     #     print(paste("extinction at time ",cbt+t_ext))
          # remove the "corresponding" species
          extinct = extinct[-length(extinct)]
          bt = c(bt,cbt+t_ext)
          bte = c(bte,Inf)
          #t_ext_limits = c(t_ext_limits,Inf)
          protected = c(protected,0)
          to = c(to,0)
          cwt = cwt - t_ext
          cbt = cbt + t_ext
          N = N-1
    #      print(paste('numberof species: ', N))
        }
      }else{
    #    print(paste("nothing between ",cbt," and ",cbt + cwt))
        key = 1
        protected = c(protected,0)
      }
    }
    N = N+1
   # print(paste('numberof species: ', N))
  }
  if(!is.null(df)){
  df = data.frame(bt = c(bt,brts),bte = c(bte, rep(Inf,length(wt))),to = c(to,rep(2,length(wt))))
  df = df[order(df$bt),]
  df$t.ext = df$bte-df$bt
  df = df[-nrow(df),]
  df = rbind(df,data.frame(bt=ct,bte=Inf,to=2,t.ext=Inf))
  df$r = r
  df$protected=protected
  }
  return(df)
}


########


augment.tree2 <- function(observed.branching.times,pars,model='dd',initspec=1,seed=0){
  if(seed>0) set.seed(seed)
  brts = observed.branching.times
  if(max(brts)==brts[1]){
    wt = - diff(c(brts,0))
  }else{
    wt = diff(c(0,brts))
  }
  mu = pars[2]
  kprima = (pars[1]*pars[3])/(pars[1]-pars[2])
  
  # vectors of missing speciation times, extinction times and missing event sequence
  bt = NULL
  bte = NULL
  to = NULL
  
  #number of species
  N = initspec
  
  #limits on extinction times
  r = NULL
  
  for(i in 1:length(wt)){
    cwt = wt[i]
    cbt = sum(wt[0:(i-1)])
    key = 0
    while(key == 0){
      if(model == "dd"){  # diversity-dependence model
        lambda = lambda.dd(pars,N)
      }
      if(model == "dd1.3"){
        lambda = lambda.dd.1.3(pars,N)
      }
      if(model=="cr"){
        lambda = lambda.cr(pars,N)
      }
      s = N*lambda
      sbte = bte[bte>cbt]
      b = get.max.speciation.bt(sbt = brts[brts>cbt],ebt = sbte,N = N,max.simultaneous.spec = floor(kprima))-cbt
      r = c(r,b)
      t.spe = rnhe(lambda=s,mu=mu,r=b)
      t_ext = ifelse(length(sbte)>0,min(sbte),Inf) - cbt
      mint = min(t.spe,t_ext)
      if(mint < cwt){
        if(mint == t.spe){#speciation
          bt = c(bt,cbt+t.spe)
          cbt = cbt + t.spe
          cwt = cwt - t.spe
          text = cbt + truncdist::rtrunc(1,"exp",a = 0, b = (b-t.spe),rate=mu)
          bte = c(bte,text)
          to = c(to,1)
          N = N + 1
        }else{#extinction
          bt = c(bt,cbt+t_ext)
          bte = c(bte,Inf)
          to = c(to,0)
          cwt = cwt - t_ext
          cbt = cbt + t_ext
          N = N - 1
        }
      }else{
        key = 1
      }
    }
    N = N+1
    # print(paste('numberof species: ', N))
  }
  df = data.frame(bt = c(bt,brts),bte = c(bte, rep(Inf,length(wt))),to = c(to,rep(2,length(wt))))
  df = df[order(df$bt),]
  df$t.ext = df$bte-df$bt
  df = df[-nrow(df),]
  df = rbind(df,data.frame(bt=sum(wt),bte=Inf,to=2,t.ext=Inf))
  df$r = r
  return(df)
}

log_sampling_prob_emphasis2 <- function(tree,pars,model=NULL,initspec){
  if(is.data.frame(tree)){
    tree = emphasis::df2tree(df=tree,pars,model=model,initspec=initspec)
    wt=tree$wt
    to=tree$to
    n=tree$n
    s=tree$s
    r=tree$r
    t_ext=tree$t_ext
    top_aug=head(tree$treeaug,-1)
    pars
  }else{
    list2env(setNames(tree, c("wt","to","n","s","r","pars","t_ext","treeaug")), .GlobalEnv)
    top_aug=treeaug
  }
  mu = pars[2]
  if(mu!=0){
    sa = s[is.finite(t_ext)]
    text = t_ext[is.finite(t_ext)]
    logg = sum(-s * (wt-(exp(-r*mu)/mu) *  (exp(wt*mu)-1)  ))+length(sa)*log(mu)-sum(mu*text)+sum(log(sa))+log.factor.samp.prob_emph(top_aug)#-sum(tree$protected)*log(2)#-(sum(tree$augtree==1)-1)*log(2)#-sum(tree$protected)*log(2)#-(sum(tree$augtree==2)-1)*log(2)
  }else{
    logg = 0
  }
  return(logg)
}

log.factor.samp.prob_emph <- function(to){
  top = head(to,-1)
  number.observed = c(1,1+cumsum(top==2))
  number.missing = c(0,cumsum(top==1)-cumsum(top==0))
  factor = -sum(log((2*number.observed+number.missing)[to==1]))#-sum(log(number.missing[to==0]))
  return(factor)
}


#maximun extinction branching time to create possible trees
get.max.speciation.bt <- function(sbt,ebt,N,max.simultaneous.spec){
  ebt = ebt[is.finite(ebt)]
  df = data.frame(brts=c(sbt,ebt),to=c(rep(1,length(sbt)),rep(-1,length(ebt))))
  df = df[order(df$brts),]
  number.of.species = c(N,N+cumsum(df$to))
  max.bt = which(number.of.species == max.simultaneous.spec & c(df$to,0) == 1 )
  if(length(max.bt)>0){
    max.t_ext = df$brts[min(max.bt)]
  }else{
    max.t_ext = max(sbt,ebt)
  }
  return(max.t_ext)
}


##work in progre
log_sampling_prob_emphasis <- function(tree,pars,model=NULL,initspec){
  if(is.data.frame(tree)){
    tree = emphasis::df2tree(df=tree,pars,model=model,initspec=initspec)
    wt=tree$wt
    to=tree$to
    n=tree$n
    s=tree$s
    r=tree$r
    t_ext=tree$t_ext
    #t_ext_limits=tree$t_ext_limits
    pars
  }else{
    list2env(setNames(tree, c("wt","to","n","s","r","pars","t_ext")), .GlobalEnv)
  }
  mu = pars[2]
  if(mu!=0){
    la = s/n
    la = la[is.finite(t_ext)]
    text = t_ext[is.finite(t_ext)]
    logg = sum(-s * (wt-(exp(-r*mu)/mu) *  (exp(wt*mu)-1)  ))+length(la)*log(mu)-sum(mu*text)+sum(log(la))-sum(tree$protected)*log(2)#-(sum(tree$augtree==1)-1)*log(2)#-sum(tree$protected)*log(2)#-(sum(tree$augtree==2)-1)*log(2)
  }else{
    logg = 0
  }
  return(logg)
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


sample.uniform <- function(brts,maxnumspec,single_dimension=NULL){
  if(is.null(single_dimension)){
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

####  

proper_weighting <- function(E){
  ta = table(E$dim)
  dims = as.numeric(names(ta))
  sizes = as.numeric(ta)
  weights = NULL
  for(i in 1:length(dims)){
    weights[i] = sum(E$weights[E$dim == dims[i]])
  }
  W = weights
  w_tilde = weights/length(E$trees)
  Z = weights/sizes
  #fhat2[j] = mean(W)
  return(list(W=W,w_tilde=w_tilde,Z=Z))
}
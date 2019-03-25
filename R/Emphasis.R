### EMPHASIS functions

# negative logLikelihood of a tree
nllik.tree = function(pars,tree,model="dd",initspec=1){
  if(is.data.frame(tree)){
    tree = emphasis::df2tree(df=tree,pars=pars,model=model,initspec=initspec)
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
  rho = pmax(n[-length(n)]*(lambda[-length(lambda)]*to+mu*(1-to)),0)
  nl = -(sum(-sigma*wt)+sum(log(rho)))
  if(min(pars)<0){nl = Inf}
  return(nl)
}



lambda.dd.1.3 <- function(pars,n){
  pmax(1e-99, pars[1]*(1-n/pars[3]))
}

lambda.dd <- function(pars,n){
  pmax(1e-99, (pars[1]-(pars[1]-pars[2])*(n/pars[3])))
}

lambda.cr <- function(pars,n){
  rep(pmax(1e-99, pars[1]), length(n))
}

Q_approx = function(pars,st,model="dd",initspec=1){
  get_llik <- function(tree) nllik.tree(pars=pars,tree=tree,initspec = initspec,model=model)
  l = sapply(st$trees, get_llik)
  w = st$w
  Q = sum(l*w)
  return(Q)
}


M_step <-function(S,init_par = c(0.5,0.5,100),model="dd"){
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
  r = df$bt[dim]-c(0,df$bt[-dim])
  t_ext = df$t.ext
  return(list(wt=wt,to=to,n=n,s=s,r=r,pars=pars,t_ext=t_ext,protected=df$protected))
}


###########################
rnhe <- function(lambda,mu,Ti,extinct=TRUE){  # random non-homogenous exponential
  ex = rexp(1)
  if(extinct){
    rv = IntInv(r=Ti,mu=mu,s=lambda,u=ex)
  }else{
    rv = IntInvExtant(r=Ti,mu=mu,s=lambda,u=ex)
  }
  if(is.na(rv)){
    rv = Inf
  }
  return(rv)
}

IntInv <- function(r,mu,s,u){
  t = -W(-exp(-r*mu+mu*u/s-exp(-r*mu)))/mu+u/s-exp(-r*mu)/mu
  return(t)
}

IntInvExtant <- function(r,mu,s,u){
  t = log(u*(mu/s)*exp(mu*r)+1)/mu
}

###  simulation of extinct species


####

# relative likelihood
rel.llik <- function(S1,p0,p1,model="dd"){
  m = length(S1)
  f1 = vector(mode='numeric',length = m)
  f2 = vector(mode='numeric',length = m)
  d = vector(mode='numeric',length = m)
  #S1 = S1$[S1$w>0]
  for(i in 1:m){
    s = S1[[i]]
    f1[i] = nllik.tree(pars=p1,tree=s,model=model)
    f2[i] = nllik.tree(pars=p0,tree=s,model=model)
    #d[i] = length(s$tree$wt)
    if(is.na(f1[i])) print(s)
  }
  Delta = -log(sum(f1/f2)/m)
  return(Delta)
}

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
    loglik.proportion = 1
  }
  M = M_step(S = sub_st,init_par = theta_0)
  M_time = get.time(time)
  pars = M$par
  hessian_inverse = try(diag(solve(M$hessian)))
  if(!is.numeric(h1)) h1 = c(NULL,NULL,NULL)
  return(list(pars=pars,fhat=fhat,se=se,st=st,loglik.proportion=loglik.proportion,hessian_inverse=hessian_inverse,E_time=E_time,M_time=M_time))
}




####### Monte Carlo E step
mc.Estep_parallel <- function(brts,pars,nsim=1000,model="dd",method="emphasis",no_cores=NULL,maxnumspec=NULL,p=0.5,seed=0,initspec=1){
  if(seed>0) set.seed(seed)
  if(brts[1]==max(brts)){
    wt = -diff(c(brts,0))
    brts = cumsum(wt)
  }
  #### parallel set-up
  if(is.null(no_cores)) no_cores <- detectCores()
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  if(method=="uniform"){
    trees <- foreach(i = 1:nsim, combine = list) %dopar% {
      S = emphasis:::sim.dim(i=i,nsim=nsim,maxnumspec=maxnumspec,deterministic = FALSE)
      ct <- max(brts)
      mbts.events = emphasis:::sim.branchingtimes.and.events(S=S ,ct = ct,p=p)
      conf = emphasis:::possible.configurations(miss = mbts.events,obs = brts)
      logg.samp = emphasis:::log.samp.prob(to = mbts.events$to,maxnumspec = maxnumspec,ct=ct,conf=conf,p=p)
      tree = list(wt=diff(c(0,conf$tree$brts,ct)),to=as.integer(conf$tree$event>0))
      logf.joint = -emphasis:::nllik.tree(pars=pars,tree=tree,model=model,initspec = 1)
      return(list(logf.joint=logf.joint,logg.samp=logg.samp,tree=tree))
    }
  }
  if(method=="emphasis"){
    obstree = emphasis::bt2tree(brts)
    trees <- foreach(i = 1:nsim, combine = list) %dopar% {
      df = emphasis::augment.tree(tree = obstree,pars = pars,model = model)
      tree = emphasis::df2tree(df,pars,model=model,initspec=1)
      logg.samp = emphasis:::log_sampling_prob_emphasis(tree = tree,pars = pars,model = model,initspec = initspec)
      logf.joint = -emphasis:::nllik.tree(pars=pars,tree=tree,model=model,initspec = 1)
      return(list(logf.joint=logf.joint,logg.samp=logg.samp,tree=tree))
    }
  }
  stopCluster(cl)
  logf = sapply(trees,function(list) list$logf.joint)
  logg = sapply(trees,function(list) list$logg.samp)
  diff_logs = logf-logg
  max_log = max(diff_logs) #
  fhat = mean(exp(diff_logs))
  se = sd(exp(diff_logs))/sqrt(nsim)
  weights = exp(diff_logs-max_log)
  logweights = diff_logs
  trees = lapply(trees, function(list) list$tree)
  return(list(trees=trees,weights=weights,logweights=logweights,fhat=fhat,fhat.se=se,logf=logf,logg=logg))
}


E_step <- function(brts,pars,nsim=1000,model="dd",method="emphasis",no_cores=2,maxnumspec=NULL,p=0.5,seed=0,parallel=TRUE){
  if(parallel){
    E = mc.Estep_parallel(brts,pars,nsim=nsim,model=model,method=method,no_cores=no_cores,maxnumspec=maxnumspec,p=p,seed=seed)
  }else{
    E = mc.Estep(brts=brts,pars=pars,nsim=nsim,model=model,method=method,maxnumspec = maxnumspec,p=p,seed=seed)
  }
  return(E)
}

mc.Estep <- function(brts,pars,nsim=1000,model="dd",method="emphasis",maxnumspec=NULL,p=0.5,seed=0){
  if(seed>0) set.seed(seed)
  if(brts[1]==max(brts)){
    wt = -diff(c(brts,0))
    brts = cumsum(wt)
  }
  trees=vector(mode = "list",length = nsim)
  logg.samp=vector(mode = "numeric",length = nsim)
  logf.joint=vector(mode = "numeric",length = nsim)
  if(method=="uniform"){
    for(i in 1:nsim){
      dim = emphasis:::sim.dim(i=i,nsim=nsim,maxnumspec=maxnumspec,deterministic = FALSE)
      ct <- max(brts)
      mbts.events = emphasis:::sim.branchingtimes.and.events(S=dim ,ct = ct,p=p)
      conf = emphasis:::possible.configurations(miss = mbts.events,obs = brts)
      logg.samp[i] = emphasis:::log.samp.prob(to = mbts.events$to,maxnumspec = maxnumspec,ct=ct,conf=conf,p=p)
      trees[[i]] = list(wt=diff(c(0,conf$tree$brts,ct)),to=as.integer(conf$tree$event>0))
      logf.joint[i] = -emphasis:::nllik.tree(pars=pars,tree=tree,model=model,initspec = 1)
    }
  }
  if(method=="emphasis"){

    tree[[i]] = emphasis::augment.tree(tree = emphasis::bt2tree(brts),pars = pars,model = model)
    logg.samp[i] = emphasis:::log_sampling_prob_emphasis(tree = tree)
    logf.joint[i] = -emphasis:::nllik.tree(pars=pars,tree=tree,model=model,initspec = 1)

  }
  diff_logs = logf-logg
  max_log = max(diff_logs) #
  fhat = mean(exp(diff_logs))
  se = sd(exp(diff_logs))/sqrt(nsim)
  weights = exp(diff_logs-max_log)
  logweights = diff_logs
  return(list(trees=trees,weights=weights,logweights=logweights,fhat=fhat,fhat.se=se,logf=logf,logg=logg))
}

######################
# Uniform data augmentation importance sampler
######################

log.samp.prob <- function(to,maxnumspec,ct,initspec=1,conf,p){
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  n = n[-length(n)]
  loggprob <- -log((maxnumspec+1))+lgamma(length(to)+1)-length(to)*log(ct)+lprobto(to,p = p)-sum(log(conf$N-conf$P))
}

lprobto <- function(to,p=0.5){
  posspec = c(0,cumsum(to==1))<(length(to)/2)
  posext = !(c(0,cumsum(to==1))==c(0,cumsum(to==0)))
  possibletotal = posspec & posext
  to_possible = to[possibletotal]
  logprob = sum(to_possible==1)*log(p)+sum(to_possible==0)*log(1-p)
  return(logprob)
}

sim.branchingtimes.and.events <- function(S=S,ct,p){
  brts = sort(runif(2*S,min=0,max=ct))
  to = sampletopology(S,p = p)
  tree = list(brts=brts,to=to)
  return(tree)
}

sim.dim <- function(maxnumspec,deterministic=FALSE,i=NULL,nsim=NULL){
  if(deterministic){
    S = floor((maxnumspec+1)*i/nsim)
  }else{
    S = sample(0:maxnumspec,1)
  }
  return(S)
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

############################
# emphasis data augmentation importance sampler
##############################
bt2tree <- function(brts){
  list(brts=brts,df=data.frame(parent=rep("s1",length(brts)-1),child=paste("s",2:length(brts),sep="")))
}

augment.tree <- function(tree,pars,model='dd',initspec=1,seed=0){
  if(seed>0) set.seed(seed)
  brts=tree$brts
  wt = diff(c(0,brts))
  ct = sum(wt)
  dim = length(wt)
  mu = pars[2]
  bt = NULL
  bte = NULL
  to = NULL
  N = initspec
  num.of.branches = length(brts)+1
  #set of protected species
  protect = 1
  
  #set of species that are going to die
  extinct = NULL
  
  #number of times we need to flip a coin to know which species dies
  protected = NULL

  for(i in 1:dim){
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
      s = N*lambda
      t.spe = rnhe(lambda=s,mu=mu,Ti=ct-cbt)
      sbte = bte[bte>(cbt+1e-09)]
      t_ext = ifelse(length(sbte)>0,min(sbte),Inf) - cbt
      mint = min(t.spe,t_ext)
      
      if(mint < cwt){
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
          text = truncdist::rtrunc(1,"exp",a = cbt+t.spe, b =ct,rate=mu)
          bte = c(bte,text)
          to = c(to,1)
          N = N + 1
          cwt = cwt - t.spe
          cbt = cbt + t.spe
        }
        else{#extinction
          # remove the "corresponding" species
          extinct = extinct[-length(extinct)]
          bt = c(bt,cbt+t_ext)
          bte = c(bte,Inf)
          protected = c(protected,0)
          to = c(to,0)
          cwt = cwt - t_ext
          cbt = cbt + t_ext
          N = N-1
        }
      }
      else{
        key = 1
        protected=c(protected,0)
      }
    }
    N = N+1
  }
  df = data.frame(bt = c(bt,brts),bte = c(bte, rep(Inf,length(wt))),to = c(to,rep(2,length(wt))))
  df = df[order(df$bt),]
  df$t.ext = df$bte-df$bt
  df = df[-nrow(df),]
  df = rbind(df,data.frame(bt=ct,bte=Inf,to=2,t.ext=Inf))
  df$protected=protected
  return(df)
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
    pars
  }else{
    list2env(setNames(tree, c("wt","to","n","s","r","pars","t_ext")), .GlobalEnv)
  }
  mu = pars[2]
  if(mu!=0){
    la = s/n
    la = la[is.finite(t_ext)]
    text = t_ext[is.finite(t_ext)]
    logg = sum(-s * (wt-(exp(-r*mu)/mu) *  (exp(wt*mu)-1)))+length(la)*log(mu)-sum(mu*text)+sum(log(la))-sum(tree$protected)*log(2)
  }else{
    logg = 0
  }
  return(logg)
}





possible.configurations  <- function(miss,obs){
  if (is.vector(obs)){
    tms <- obs
    to <- rep(1,length(obs))
    obs <- list(brts=tms,to=to)
  }
  mi = 1 #index for missing
  ob = 1 #index for observed
  
  # (number of) protected species
  protected = 1
  P<-NULL
  
  # (number of) current species
  currentspecies = 1
  N<-NULL
  
  # missing branching times
  brts.m = c(miss$brts,Inf)
  
  # set of sets of guardians
  guardians = list()
  
  # missing new species (starting to count at 1 above number of present species)
  n.obs = length(obs$brts)
  newspecies.m = n.obs+1
  
  # observed new species (starting to count at 1, so next one is 2)
  newspecies.o = 2
  
  # sampled tree
  tree<-list(brts=NULL,species=NULL,event=NULL)
  
  while (mi < length(brts.m) | ob < length(obs$brts)) {
    if(obs$brts[ob] < brts.m[mi]){ # observed speciation
      spec = obs$to[ob]
      #update tree
      tree$brts = c(tree$brts,obs$brts[ob])
      tree$species =c(tree$species,spec)
      tree$event = c(tree$event,newspecies.o)
      if(spec%in%protected){ # if protected, then both species become protected
        protected = c(protected,newspecies.o)
        currentspecies = c(currentspecies,newspecies.o)
      }else{ # if unprotected, then then (1) its guardian set disappears and (2) both become protected
        index = unlist(lapply(guardians,function(y,x){x%in%y},x=spec))
        if (sum(index)>0){
          index = which(index)
          n.guardians = length(guardians[[index]])
          guardians[[index]] = NULL
        }
        protected = c(protected,spec,newspecies.o)
        P = c(P,0) # it is weird, but we want to try 
        currentspecies = c(currentspecies,newspecies.o)
        N = c(N,n.guardians) # it is weird, but we want to try 
      }
      ob = ob + 1
      newspecies.o = newspecies.o + 1
    }else{ # missing event
      if(miss$to[mi]==1){ # missing speciation
        mspec = sample(c(currentspecies,currentspecies),1)
        index = which(mspec==protected)
        N = c(N,length(currentspecies))
        P=c(P,0)
        currentspecies = c(currentspecies,newspecies.m)
        #update tree
        tree$brts = c(tree$brts,miss$brts[mi])
        tree$species =c(tree$species,mspec)
        tree$event = c(tree$event,newspecies.m)
        if(sum(index)>0){ # if a protected species speciates, then it becomes unprotected and both guardians
          protected = protected[-index]
          guardians[[length(guardians)+1]] = c(mspec,newspecies.m)
        }else{ # if a unprotected species speciates, then ...
          index = unlist(lapply(guardians,function(y,x){x%in%y},x=mspec))          
          if (sum(index)>0){ #... if it is a guardian then new species becomes guardian
            index=which(index)
            guardians[[index]] = c(guardians[[index]],newspecies.m)
          } # ... if it is not a guardian then no changes to guardianship
        }
        newspecies.m = newspecies.m + 1
      } else { #missing extinction
        N = c(N,length(currentspecies))
        P = c(P,length(protected))
        available = setdiff(currentspecies,protected)
        missextinct = sample(c(available,available),1)
        if (missextinct<=n.obs){# if we selected the label of an extant species, we arbitrarily pick the label of a missing guardian
          index = which(unlist(lapply(guardians,function(y,x){x%in%y},x=missextinct)))         
          missextinct = setdiff(guardians[[index]],missextinct)[1]
        }
        #update tree
        tree$brts = c(tree$brts,miss$brts[mi])
        tree$species =c(tree$species,missextinct)
        tree$event = c(tree$event,0)
        index = unlist(lapply(guardians,function(y,x){x%in%y},x=missextinct))          
        if (sum(index)>0){ # if it is a guardian then ...
          index=which(index)
          if (length(guardians[[index]])>2){ # ... if the set is larger than 2, then take it out of guardian set
            guardians[[index]] = setdiff(guardians[[index]],missextinct)
          } else { # ... if guardian set is of size 2, remove guardian set and protect the remaining species
            protected = c(protected, setdiff(guardians[[index]],missextinct))
            guardians[[index]] = NULL
          }
        }
        currentspecies = setdiff(currentspecies,missextinct)
      }
      mi = mi + 1
    }
  }
  return(list(N=N,P=P,tree=tree))
}


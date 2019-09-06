


Q_approx = function(pars,st,model="dd",initspec=1){
  get_llik <- function(tree) nllik.tree(pars=pars,tree=tree,initspec = initspec,model=model)
  l = sapply(st$trees, get_llik)
  w = st$weights/(sum(st$weights))
  Q = sum(l*w)
  return(Q)
}

Q_SAEM = function(sample,previous_Q,gamma){
  
}

M_step <-function(st,init_par = NULL,model="dd",proportion_of_subset=1){
  time0 = proc.time()
  weights = st$weights/sum(st$weights)
  weights_sorted = sort(weights,decreasing = TRUE)
  a = which(cumsum(weights_sorted)>=proportion_of_subset)[1]
  max.weight = weights_sorted[a]
  sub_st = lapply(st, "[", weights>=max.weight)
  loglik_proportion = sum(weights[weights>=max.weight])/sum(weights)
  effective_sample_size = length(weights[weights>=max.weight])
  po = subplex(par = init_par, fn = Q_approx,st = sub_st,model=model,hessian = TRUE)
  M_time = get.time(time0)
  return(list(po=po,loglik_proportion=loglik_proportion,effective_sample_size=effective_sample_size,M_time=M_time))
}

##to do: M-step parallel

df2tree <- function(df,pars,model="dd",initspec=1){
  dim = nrow(df)
  if(is.null(df$brts)) df$brts = df$bt
  wt = diff(c(0,df$brts))
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
  return(list(wt=wt,to=to,n=n,s=s,r=r,pars=pars,t_ext=t_ext,treeaug=df$to))
}


###########################

###  simulation of extinct species

####


####### Monte Carlo E step

mc.Estep_parallel <- function(brts,pars,pars3=NULL,nsim=1000,model="dd",method="emphasis",no_cores=NULL,maxnumspec=NULL,seed=0,initspec=1,single_dimension=NULL,limit_on_species=NULL){
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
      df = emphasis::augment.tree(brts,pars = pars3,model = model,limit_on_species=limit_on_species)
      if(!is.null(df)){
        df$brts = df$bt
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
        logf.joint = -emphasis:::nllik.tree(pars=pars,tree=df,model=model,initspec = 1)
        return(list(logf.joint=logf.joint,logg.samp=log.samp.unif.prob,dim=nrow(df),tree=df))
      }
    }
  stopCluster(cl)
  logf = sapply(trees,function(list) list$logf.joint)
  logg = sapply(trees,function(list) list$logg.samp)
  dim = sapply(trees,function(list) list$dim)
  diff_logs = logf-logg
  max_log = max(diff_logs) 
  fhat = mean(exp(diff_logs))
  se = sd(exp(diff_logs))/sqrt(nsim)
  weights = exp(diff_logs)
  logweights = diff_logs
  trees = lapply(trees, function(list) list$tree)
  time = get.time(time)
  E = list(weights=weights,logweights=logweights,fhat=fhat,fhat.se=se,logf=logf,logg=logg,trees=trees,dim=dim,Etime=time)
  pw = proper_weighting(E)
  fhat2 = log(mean(pw$W))
  return(c(E,list(fhat2=fhat2)))
}


E_step <- function(brts,pars,nsim=1000,model="dd",method="emphasis",no_cores=2,maxnumspec=NULL,seed=0,parallel=TRUE){
  if(parallel){
    E = mc.Estep_parallel(brts,pars,nsim=nsim,model=model,method=method,no_cores=no_cores,maxnumspec=maxnumspec,seed=seed)
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


augment.tree_old <- function(brts,pars,model='dd',initspec=1,seed=0,limit_on_species=NULL){
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
  kprima = -pars[1]/pars[2]
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
      t.spe = rnhe(lambda=s,mu=pars[3],r=b)
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
          text = cbt + truncdist::rtrunc(1,"exp",a = 0, b = (b-t.spe),rate=pars[3])
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


###

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
  mu = pars[3]
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

####  


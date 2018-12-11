### EMPHASIS functions

# negative logLikelihood of a tree
nllik.tree = function(pars,tree,topology=T,model="dd",truncdim=F,initspec=2){
  wt = tree$wt
  to = tree$to
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
  mu = pars[2]
  sigma = (lambda + mu)*n
#  sigma[n==2] = lambda[n==2]*2
  if(topology){
    rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  }else{
    rho = pmax(n[-length(n)]*lambda[-length(lambda)]*to+mu*(1-to),0)
  }
  if(truncdim){
    sigma = sigma[-length(sigma)]
    wt = wt[-length(wt)]
  }
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

# negative logLikelihood of a set of trees
Q.approx = function(pars,st,topology,model="dd"){
  m = length(st$trees)
  l = vector(mode="numeric",length=m)
  for(i in 1:m){
    l[i] = nllik.tree(pars,tree=st$trees[[i]],topology=topology,model=model)
  }
  w = st$w
  Q = sum(l*w)
  return(Q)
}

# MLE for a set of trees
mle.st <-function(S,init_par = c(0.5,0.5,100),topology=TRUE,model="dd"){
  po = subplex(par = init_par, fn = Q.approx, topology=topology,st = S,model=model,hessian = TRUE)
  return(po)
}

lg_prob <- function(tree,topology=T){
  list2env(setNames(tree, c("wt","to","n","s","r","pars","t_ext")), .GlobalEnv)
  mu = pars[2]
  if(mu!=0){
    term1 = sum(-s*(wt+(exp(-r*mu)/mu)*(1-exp(mu*wt))))
    la = s/n
    la = la[is.finite(t_ext)]
    text = t_ext[is.finite(t_ext)]
    if(topology){
      term2 = length(la)*log(mu)-sum(mu*text)+sum(log(la))
    }else{
      term2 = length(la)*log(mu)-sum(mu*text)+sum(log(s[is.finite(t_ext)]))
    }
    logg = term1 + term2
  }else{
    logg = 0
  }
  return(logg)
}

df2tree <- function(df,pars,model="dd"){
  dim = dim(df)[1]
  wt = diff(c(0,df$bt))
  to = df$to[-dim]
  to[to==2] = 1
  n = c(2,2+cumsum(to)+cumsum(to-1))
  if(model=="dd"){
    s = n*lambda.dd(pars,n)
  }
  if(model == "dd.1.3"){
    s = n*lamda.dd.1.3(pars,n)
  }
  r = df$bt[dim]-c(0,df$bt[-dim])
  t_ext = df$t.ext
  return(list(wt=wt,to=to,n=n,s=s,r=r,pars=pars,t_ext=t_ext))
}
# Monte-Carlo sampling function / simulation of a set of complete trees
sim.sct <- function(brts,pars,m=10,print=TRUE,topology=TRUE,model="dd",truncdim=FALSE){
  no_cores <- detectCores()
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  trees <- foreach(i = 1:m, combine = list) %dopar% {
    df =  emphasis::sim.extinct2(brts = brts,pars = pars,model=model)
    tree = emphasis::df2tree(df,pars,model=model)
    lsprob = emphasis::lg_prob(tree,topology = topology)
    lw2 = we_cal(tree)
    nl = emphasis::nllik.tree(pars,tree=tree,topology = topology,model=model,truncdim = truncdim)
    lw = -nl-lsprob#-DDD:::dd_loglik(pars1 = pars, pars2 = pars2,brts = brts, missnumspec = 0)
    fms = df$bt[is.finite(df$bte)][1] #first missing speciation
    fe = df$bt[df$to==0][1] # first extinction
    return(list(nl=nl,lw=lw,tree=tree,lsprob=lsprob,fms = fms,fe=fe,lw2=lw2))
  }
  stopCluster(cl)
  lw = sapply(trees,function(list) list$lw)
  lw2 = sapply(trees,function(list) list$lw2)
  dim = sapply(trees,function(list) length(list$tree$wt))
  nl = sapply(trees,function(list) list$nl)
  lg = sapply(trees, function(list) list$lsprob)
  fms = sapply(trees, function(list) list$fms)
  fe = sapply(trees, function(list) list$fe)
  trees = lapply(trees, function(list) list$tree)
  w = exp(lw)
  eg = exp(lg)
  if(print){
    g = qplot(dim,w)
    print(g)
  }
  return(list(w=w, dim=dim, nl=nl, trees=trees,g=eg,fms=fms,fe=fe,lw2=lw2))
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
sim.extinct <- function(brts,pars,model='dd',seed=0){
  if(seed>0) set.seed(seed)
  wt = -diff(c(brts,0))
  ct = sum(wt)
  dim = length(wt)
  mu = pars[2]
  bt = NULL
  bte = NULL
  to = NULL
  N = 2
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
      sbte = bte[bte>cbt]
      t_ext = ifelse(length(sbte)>0,min(sbte),Inf) - cbt
      mint = min(t.spe,t_ext)
      if(mint < cwt){
        if(mint == t.spe){#speciation
          bt = c(bt,cbt+t.spe)
          text = truncdist::rtrunc(1,"exp",a = cbt+t.spe, b =ct,rate=mu)
          bte = c(bte,text)
          to = c(to,1)
          N = N + 1
          cwt = cwt - t.spe
          cbt = cbt + t.spe
        }
        else{#extinction
          bt = c(bt,cbt+t_ext)
          bte = c(bte,Inf)
          to = c(to,0)
          cwt = cwt - t_ext
          cbt = cbt + t_ext
          N = N-1
        }
      }
      else{
        key = 1
      }
    }
    N = N+1
  }
  df = data.frame(bt = c(bt,ct-brts),bte = c(bte, rep(Inf,length(wt))),to = c(to,rep(2,length(wt))))
  df = df[order(df$bt),]
  df$t.ext = df$bte-df$bt
  df = df[-1,]
  df = rbind(df,data.frame(bt=ct,bte=Inf,to=2,t.ext=Inf))
  return(df)
}

####

# relative likelihood
rel.llik <- function(S1,p0,p1,model="dd"){
  m = length(S1)
  f1 = vector(mode='numeric',length = m)
  f2 = vector(mode='numeric',length = m)
  d = vector(mode='numeric',length = m)
  S1 = S1$rec[S1$w>0]
  for(i in 1:m){
    s = S1[[i]]
    f1[i] = nllik.tree(pars=p1,tree=s,model=model)
    f2[i] = nllik.tree(pars=p0,tree=s,model=model)
    d[i] = length(s$tree$wt)
    if(is.na(f1[i])) print(s)
  }
  Delta = -log(sum(f1/f2)/m)
  return(Delta)
}

# Pilot study
pilot.study <- function(brts,epsilon,m1=10,printprocess=FALSE,init_par=c(1.2,0.3,60),l1=20,model="dd"){
  # pilot study suggested by Chan et. al
  pars = init_par
  M = matrix(ncol = 3,nrow = l1)
  H = matrix(ncol = 3,nrow = l1)
  DD = NULL
  LL = NULL
  for(i in 1:l1){
    S = sim.sct(brts,pars,m=m1,print = F,model =  model)
    mle =  mle.st(S = S, model=model)
#    L = obs.lik.approx(pars = pars,st = S)
    lL = log(mean(S$w))
    if(model=="dd.1.3"){pars2[2]=1.3}
    DDD <- DDD:::dd_loglik(pars1 = pars, pars2 = pars2,
                           brts = brts, missnumspec = 0)
    pars = mle$par
    H[i,] = try(diag(solve(mle$hessian))/m1)
    M[i,] = pars

    DD = c(DD,DDD)
    LL = c(LL,lL)
    save(DD,LL,M,file = "DyL.RData")
    print(paste('l:',lL,'DDDlogLik:',DDD,'pars:',pars[1],pars[2],pars[3]))
  }
  l = 10
  PM = M[1:(l1-10),]
  PH = H[1:(l1-10),]
  M = M[(l1-9):l1,]
  H = H[(l1-9):l1,]
  Q = vector(mode="numeric",length = (l1-10))
  MLE = list()
  for(i in 1:10){
    Delta = vector(mode="numeric",length = l)
    Me = matrix(ncol = 3,nrow = l)
    if(printprocess) print(paste('iteration',i))
    for(j in 1:l){
      S = sim.sct(brts,pars=M[i,],m=m1,print = F,model=model)
      mle = mle.st(S = S,model=model)
      pars = mle$par
      Me[j,] = pars
      Delta[j] = rel.llik(S1 = S,p0 = M[i,], p1 = pars,model=model)
    }
    MLE[[i]] = Me
    mD = mean(Delta)
    Q[i] = sum((Delta-mD)^2)
  }
  s2 = sum(Q)/((l-1)*(l1-10+1))
  s1 = sqrt(s2)
  m = m1*s1/epsilon
  m = floor(m) + 1
  return(list(m=m,p=M[10,],s1=s1,M=M,H=H,MLE=MLE,PM=PM,PH=PH))
}

#MCEM
mcem.tree <- function(brts,p,model="dd"){
  m = p$m
  s1 = p$s1
  sig = 100*s1/m
  tol = 2*sig*sqrt(1/5)
  D = Inf
  k = 1
  print("initializing mcem")
  pars = p$p
  PARS = pars
  H = c(NULL,NULL,NULL)
  Me = p$M
  tE = NULL
  tM = NULL
  efficiency = NULL
  prop = 1
  while(abs(D)>tol){
    if(m*prop<p$m) m = m/prop
    time = proc.time()
    S = sim.sct(brts = brts,pars=pars,m = m,model=model)
    prop = sum(S$w>0)/length(S$w)
    efficiency = c(efficiency,prop)
    tE = c(tE,get.time(time))
    time = proc.time()
    M = mle.st(S = S,model=model)
    tM = c(tM,get.time(time))
    mle = M$par
    h1 = try(diag(solve(M$hessian))/m)
    if(is.numeric(h1)) H =  rbind(H,h1)
    D = rel.llik(S1 = S,p0 = pars,p1 = mle,model=model)
    PARS = rbind(PARS,mle)
    pars = mle
    print(paste("iteration",k,"Q: ",M$value,'proportion of useful trees',prop,'sampling size',m*prop, " lambda: ", pars[1]," mu: ", pars[2], "K:", pars[3]))
    k = k+1
  }
  PARS = data.frame(it=1:(dim(Me)[1]+dim(PARS)[1]),lambda = c(Me[,1],PARS[,1]),mu=c(Me[,2],PARS[,2]),K=c(Me[,3],PARS[,3]))
  return(list(pars=pars,PARS=PARS,H=H,tE=tE,tM=tM,efficiency=efficiency))
}


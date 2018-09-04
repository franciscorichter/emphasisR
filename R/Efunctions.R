### EMPHASIS functions

#negative logLikelihood of a tree
nllik.tree = function(pars,tree){
  wt = tree$wt
  to = tree$E
  n = c(2,2+cumsum(to)+cumsum(to-1))
  lambda = (pars[1]-(pars[1]-pars[2])*(n/pars[3]))
  mu = pars[2]
  sigma = (lambda + mu)*n
  rho = pmax(lambda[-length(lambda)]*to+mu*(1-to),0)
  nl = -(sum(-sigma*wt)+sum(log(rho)))
  if(min(pars)<0){nl = Inf}
  return(nl)
}

# negative logLikelihood of a set of trees
Q.approx = function(pars, st){
  m = length(st$rec)
  l = vector(mode = 'numeric',length = m)
  w = vector(mode = 'numeric',length = m)
  for(i in 1:m){
    s = st$rec[[i]]# complete tree
    w[i] = st$w[i]# corresponding weight
    if(w[i]!=0){# if weight is non-zero, calculate likelihood
      l[i] = nllik.tree(pars,tree=s)
    }else{
      l[i] = 0
    }
  }
  L = sum(l*w)
  return(L)
}

# MLE for a set of trees
mle.st <-function(S,init_par =c(0.5,0.5,100)){
  po =subplex(par = init_par, fn = Q.approx, st=S,hessian = TRUE)
  return(po)
}

# Monte-Carlo sampling function / simulation of a set of complete trees
sim.sct <- function(brts,pars,m=10,oc=0,print=TRUE){
  no_cores <- detectCores() - oc
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  trees <- foreach(i = 1:m, combine = list) %dopar% {
    df =  emphasis::sim.extinct2(brts = brts,pars = pars)
    lw = emphasis::logweight(pars,df)
    wtT = c(diff(c(0,df$bt)))
    E = df$to[-length(df$to)]
    E[E==2] = 1
    return(list(wt=wtT,E=E,lw=lw))
  }
  stopCluster(cl)
  lw = sapply(trees,function(list) list$lw)
  dim = sapply(trees,function(list) length(list$wt))
  w = exp(lw)
  if(print){
    g = qplot(dim,w)
    print(g)
  }
  return(list(rec = trees, w=w,dim=dim))
}

###########################
extinction.processes <- function(u,inits,mu0){
  nm = length(u)
  t.ext = vector(mode='numeric',length=nm)
  if(nm > 0){
    for(i in 1:nm){
      t.ext[i] = inits[i] - log(1-u[i])/mu0  #Inverse of the intensity function for constant extinction rate
    }
  }
  return(t.ext)
}
###  simulation of extinct species
sim.extinct <- function(brts,pars,model='dd',seed=0){
  if(seed>0) set.seed(seed)
  wt = -diff(c(brts,0))
  ct = sum(wt)
  lambda0 = pars[1]
  mu0 = pars[2]
  K = pars[3]
  dim = length(wt)
  ms = NULL # missing speciations, for now we just add time. When we consider topology we do it with species as well
  me = NULL # missing extinctions (in the uniform plane)
  bt = NULL
  bte = NULL
  to = NULL
  cbt = 0
  N = 2
  for(i in 1:dim){
    cwt = wt[i]
    cbt = sum(wt[0:(i-1)])
    key = 0
    while(key == 0){
      if(model == "dd"){  # diversity-dependence model
        lambda = max(1e-99, lambda0 - (lambda0-mu0)*N/K)
        mu = mu0
        s = N*lambda
      }else{print('Model not implemented yet, try dd')}
      t.spe = rexp(1,s)
      t.ext = extinction.processes(u=me,inits=ms,mu0=mu0)
      t_ext = ifelse(length(t.ext)>0,min(t.ext),Inf)-cbt  # if is not empty gives the waiting time for next extinction
      mint = min(t.spe,t_ext)
      if(mint < cwt){
        if(mint == t.spe){#speciation
          u = runif(1)
          if(u < pexp(ct-(cbt+t.spe),mu)){
            ms = c(ms,cbt+t.spe)
            me = c(me,u)
            bt = c(bt,cbt+t.spe)
            text = extinction.processes(u=u,inits=cbt+t.spe,mu0=mu0)
            bte = c(bte,text)
            to = c(to,1)
            N = N + 1
          }
          cwt = cwt - t.spe
          cbt = cbt + t.spe
        }
        else{#extinction
          extinctone = which(t.ext == min(t.ext))
          text = t.ext[extinctone]
          bt = c(bt,text)
          bte = c(bte,Inf)
          to = c(to,0)
          ms = ms[-extinctone]
          me = me[-extinctone]
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
  df = rbind(df,data.frame(bt=ct,bte=Inf,to=0,t.ext=Inf))
  return(df)
}

#weight
logweight <- function(pars,df,ct=NULL){
  dim = dim(df)[1]
  ct = df$bt[dim]
  wtT = diff(c(0,df$bt))
  n.tree = list(wt=wtT,E=df$to[-dim])
  n.tree$E[n.tree$E==2] = 1
  E = n.tree$E
  n = c(2,2+cumsum(E)+cumsum(E-1))
  lambda = (pars[1]-(pars[1]-pars[2])*(n/pars[3]))
  mu = pars[2]
  s = lambda*n
  lsprob = da.prob2(wt=wtT,t_ext=df$t.ext,s=s,mu=pars[2],r=ct-c(0,df$bt[-length(df$bt)]),n=n)
  lrprob = -nllik.tree(pars,n.tree)
  logweight = lrprob-lsprob
  return(logweight)
}

da.prob2 <- function(wt,t_ext,s,mu,r,n){
  t1 = -s*(wt+(exp(-r*mu)/mu)*(1-exp(mu*wt)))
  la = s/n
  la = la[t_ext<999]
  text = t_ext[t_ext<999]
  t2 = length(la)*log(mu)-sum(mu*text)+sum(log(la))
  logg = sum(t1) + t2
  return(logg)
}

####

# relative likelihood
rel.llik <- function(S1,p0,p1){
  m = length(S1)
  f1 = vector(mode='numeric',length = m)
  f2 = vector(mode='numeric',length = m)
  d = vector(mode='numeric',length = m)
  S1 = S1$rec[S1$w>0]
  for(i in 1:m){
    s = S1[[i]]
    f1[i] = nllik.tree(pars=p1,tree=s)
    f2[i] = nllik.tree(pars=p0,tree=s)
    d[i] = length(s$tree$wt)
    if(is.na(f1[i])) print(s)
  }
  Delta = -log(sum(f1/f2)/m)
  return(Delta)
}

# Pilot study
pilot.study <- function(brts,epsilon,m1=10,printprocess=FALSE,init_par=c(1.2,0.3,60),l1=20){
  # pilot study suggested by Chan et. al
  pars = init_par
  M = matrix(ncol = 3,nrow = l1)
  H = matrix(ncol = 3,nrow = l1)
  DD = NULL
  LL = NULL
  for(i in 1:l1){
    S = sim.sct(brts,pars,m=m1)
    mle =  mle.st(S = S)
    L = obs.lik.approx(pars = pars,st = S)
    S1 = sim.sct(brts,realpars,m=m1)
    L1 = obs.lik.approx(pars = pars,st = S1)
    DDD <- DDD:::dd_loglik(pars1 = pars, pars2 = pars2,
                           brts = brts, missnumspec = 0)
    pars = mle$par
    H[i,] = try(diag(solve(mle$hessian))/m1)
    M[i,] = pars

    DD = c(DD,DDD)
    LL = c(LL,L)
    save(DD,LL,M,file = "DyL.RData")
    print(paste('L-log(m):',L-log(m1),'L1-log(m):',L1-log(m1),'2DDDlogLik:',2*DDD,'pars:',pars[1],pars[2],pars[3]))
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
      S = sim.sct(brts,pars=M[i,],m=m1)
      mle = mle.st(S = S)
      pars = mle$par
      Me[j,] = pars
      Delta[j] = rel.llik(S1 = S,p0 = M[i,], p1 = pars)
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
mcem.tree <- function(brts,p){
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
    S = sim.sct(brts = brts,pars=pars,m = m)
    prop = sum(S$w>0)/length(S$w)
    efficiency = c(efficiency,prop)
    tE = c(tE,get.time(time))
    time = proc.time()
    M = mle.st(S = S)
    tM = c(tM,get.time(time))
    mle = M$par
    h1 = try(diag(solve(M$hessian))/m)
    if(is.numeric(h1)) H =  rbind(H,h1)
    D = rel.llik(S1 = S,p0 = pars,p1 = mle)
    PARS = rbind(PARS,mle)
    pars = mle
    print(paste("iteration",k,"Q: ",M$value,'proportion of useful trees',prop,'sampling size',m*prop, " lambda: ", pars[1]," mu: ", pars[2], "K:", pars[3]))
    k = k+1
  }
  PARS = data.frame(it=1:(dim(Me)[1]+dim(PARS)[1]),lambda = c(Me[,1],PARS[,1]),mu=c(Me[,2],PARS[,2]),K=c(Me[,3],PARS[,3]))
  return(list(pars=pars,PARS=PARS,H=H,tE=tE,tM=tM,efficiency=efficiency))
}

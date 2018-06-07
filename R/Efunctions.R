### functions

update.tree <- function(tree,t_spe,t_ext){
  wt = tree$wt
  E = tree$E
  ct = sum(wt)
  if(t_ext > ct){
    stop('Extinction beyond present!')
  }
  if(t_ext < t_spe){
    stop('Speciation after extinction!')
  }
  # speciation
  K = length(wt)
  k = length(wt[cumsum(wt) < t_spe])
  if((k+1)<K){
    lastbit = E[(k+1):(K-1)]
  }
  else{
    lastbit = NULL
  }
  E = c(E[0:k],1,lastbit)
  # (this can be witten in one line, is convinent?)
  if(k==0){
    wt = c(t_spe, wt[1]-t_spe, wt[2:K])
  }
  if(k > 0 & k < (K-1)){
    wt = c(wt[1:k], t_spe - sum(wt[1:k]), wt[k+1]-(t_spe - sum(wt[1:k])), wt[(k+2):K])
  }
  if(k == (K-1)){
    wt = c(wt[1:(K-1)],t_spe-sum(wt[1:(K-1)]),ct-t_spe)
  }
  #extinction
  K = length(wt)
  k = length(wt[cumsum(wt) < t_ext])
  if((k+1)<K){
    lastbit = E[(k+1):(K-1)]
  }
  else{
    lastbit = NULL
  }
  E = c(E[0:k],0,lastbit)
  # (this can be witten in one line, is convinent?)
  if(k==0){
    wt = c(t_ext, wt[1]-t_ext, wt[2:K])
  }
  if(k > 0 & k < (K-1)){
    wt = c(wt[1:k], t_ext - sum(wt[1:k]), wt[k+1]-(t_ext - sum(wt[1:k])), wt[(k+2):K])
  }
  if(k == (K-1)){
    wt = c(wt[1:(K-1)],t_ext-sum(wt[1:(K-1)]),ct-t_ext)
  }
  tree = list(wt=wt,E=E)
  return(tree)
}
sampprob <- function(t,s,mu,r){  ## equation (*)
  term1 = s*(1-exp(-mu*(r-t)))
  c = exp(-(s/mu)*exp(-mu*r))
  if(c==0){
    term2 = 0
  }else{
    term2 = c^{-exp(mu*t)+1}
  }
  f = term1*exp(-s*t)*term2
  return(f)
}

# simulation of missing part [NEED TO CLEAN AND SIMPLIFY]
sim.extinct <- function(tree,pars,model='dd',seed=0, adjustment=FALSE){
  if(seed>0) set.seed(seed)
  wt = tree$wt
  ct = sum(wt)
  tree$E = rep(1,(length(wt)-1))
  lambda0=pars[1]
  mu0=pars[2]
  K=pars[3]
  limit = lambda0*K/(lambda0-mu0)
  rs = dim = length(wt)
  if (limit < (dim-1)){
    print('parameters do not make sense, observed tree implies negative rates')
  }
  ms = NULL # missing species, for now we just add time. When we consider topology we do it with species as well
  e.lims = NULL # limits on extinctions
  cbt = 0
  N = 2
  nm= 0 # number of missing species
  rprob = NULL # true probability of Missing|observed
  sprob = NULL # sampling probability of Missing|observed
  et = NULL # event type
  h = 1 # index to fill probabilities
  dif1 = vector(mode = 'numeric',length = dim)
  dif2 = vector(mode = 'numeric', length = dim)
  for(i in 1:dim){
    rs = rs-1 # this are the remaning speciations
    cwt = wt[i]
    cbt = sum(wt[0:(i-1)])
    key = 0
    last='nothing'
    while(key == 0){
      if(model == "dd"){  # diversity-dependence model
        lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
        mu = mu0
        lambda = rep(lambda,N)
      }
      if(model == 'cr'){ # constant-rate model
        lambda = rep(lambda0,N)
        mu = rep(mu0,N)
      }
      s = sum(lambda)
      if (s == 0){
        t.spe = Inf
      }
      else{
        t.spe = rexp(1,s)
      }
      if(nm > 0){ # if there are missing species simulate extinction times
        t.ext = vector(mode = 'numeric',length = nm)
        for(j in 1:nm){
          t.ext[j] = truncdist::rtrunc(1,'exp',a=0,b=(e.lims[j]-cbt),rate=mu)
        }
        extinctedone = which(t.ext == min(t.ext))
        t.ext = min(t.ext)
      }
      else{
        t.ext = Inf
      }
      if(t.ext ==Inf & t.spe == Inf){
        print('try another K!')
      }
      mint = min(t.spe,t.ext)
      if(nm > 0){
        probs = vector(mode = 'numeric',length = nm)
        for(j in 1:nm){
          probs[j] = 1-truncdist::ptrunc(mint,'exp',a=0,b=(e.lims[j]-cbt),rate=mu)
        }
      }
      else{probs = 1}
      if(mint < cwt){
        if(mint == t.spe & mint != t.ext){#speciation
          acep = runif(1)
          thre = pexp(ct-cbt,mu)
          if(acep < thre){
            ms = c(ms,cbt+t.spe)
            rprob[h] = dexp(x = t.spe, rate = (s+nm*mu))*(s/(s+nm*mu))
            sp = sampprob(t = mint,s = s,mu = mu,r = ct-cbt)#/integrate(sampprob,lower = 0, upper = ct-cbt,s=s,mu=mu,r=ct-cbt)$value
            sprob[h] = prod(probs)*sp
            et[h] = 'speciation'
            h = h + 1
            nm = nm + 1 # number of missing species
            if((N + rs-1) < limit){
              e.lims = c(e.lims,ct)
            }
            else{
              e.lims = c(e.lims,sum(wt[1:(dim - (N+rs-floor(limit)-1))])) # cual es la interpretacion de esta formula?
            }
            N = N+1
            #print(paste('at branching time',cbt+t.spe,'a missing species arises, resulting on',N,'current species, happening before the current waiting time',cwt))
            last = 'speciation'
          }
          else{last = 'nothing'}
          cwt = cwt - t.spe
          cbt = cbt + t.spe
        }
        else{#extinction
          pickone = sample(1:nm,1)
          t_spe = ms[pickone]
          t_ext = cbt + t.ext
          tree = update.tree(tree,t_spe=t_spe,t_ext=t_ext)
          rprob[h] = dexp(x = mint, rate = (s+nm*mu))*(mu/(s+nm*mu))
          et[h] = 'extinction'
          probs = probs[-extinctedone]
          sprob[h] = prod(probs)*truncdist::dtrunc(mint,'exp',a=0,b=(e.lims[extinctedone]-cbt),rate=mu)*(1-integrate(sampprob,lower = 0, upper = mint,s=s,mu=mu,r=ct-cbt)$value)#/integrate(sampprob,lower = 0, upper = ct-cbt,s=s,mu=mu,r=ct-cbt)$value)
          ms = ms[-pickone]
          cwt = cwt - t.ext
          cbt = cbt + t.ext
          N = N-1
          h = h+1
          nm = nm - 1
          e.lims = e.lims[-extinctedone]
          last = 'extinction'
        }
      }
      else{
        key = 1
        rprob[h] = pexp(q = cwt, rate = (s+nm*mu0),lower.tail = FALSE)
        et[h] = 'nothing'
        sprob[h] = prod(probs)*(1 - integrate(Vectorize(sampprob),lower = 0, upper = cwt,s=s,mu=mu,r=ct-cbt)$value)
        h = h+1
        dif1[i] = cwt/wt[i]
        dif2[i] = (mint - cwt)/wt[i]
        if(adjustment & cwt<(mint - cwt) & last=='speciation'){ #Adjusting last speciation
          ms = ms[-length(ms)]
          e.lims = e.lims[-length(e.lims)]
          nm = nm - 1
          N = N-1
        }
      }
    }
    N= N+1
  }
  tree$rprob = rprob
  tree$sprob = sprob
  tree$et = et
  tree$weight = prod(rprob)/prod(sprob)
  logweight = log(rprob)-log(sprob)
  tree$logweight = sum(logweight)
  E = tree$E
  n = c(2,2+cumsum(E)+cumsum(E-1))
  tree$n = n
  return(tree)
}
#negative logLikelihood of a tree
nllik.tree = function(pars,tree){
  b = c(pars[1],(pars[1]-pars[2])/pars[3],pars[2])
  ldt = tree$wt[length(tree$wt)]
  dt = tree$wt[1:(length(tree$wt)-1)]
  lastn = tree$n[length(tree$n)]
  n = tree$n[1:(length(tree$n)-1)]
  E = tree$E
  sigma = n*(b[1]-b[2]*n + b[3]) #n-dimentional
  rho = pmax(b[1]*E-b[2]*n*E+b[3]*(1-E),0)
  l = -sum(-sigma*dt+log(rho)-lastn*ldt)
  if(min(b)<0){l = Inf}
  return(l)
}
# negative logLikelihood of a set of trees
nllik.st = function(pars, st){
  m = length(st$rec)
  l = vector(mode = 'numeric',length = m)
  w = vector(mode = 'numeric',length = m)
  for(i in 1:m){
    s = st$rec[[i]]
    w[i] = st$w[i]
    l[i] = nllik.tree(pars,tree=s)
  }
  L = sum(l*w)
  return(L)
}
# relative likelihood
rel.llik <- function(S1,p0,p1){
  m = length(S1)
  f1 = vector(mode='numeric',length = m)
  f2 = vector(mode='numeric',length = m)
  d = vector(mode='numeric',length = m)
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
# MLE for a set of trees
mle.st <- function(S,init_par = c(0.5,0.5,100)){
  po = subplex(par = init_par, fn = nllik.st, st=S,hessian = TRUE)
  return(po)
}
# Monte-Carlo sampling / simulation of a set of complete trees
sim.sct <- function(tree,pars,m=10,printHistD=F,parallel=FALSE){
  if(parallel){
    no_cores <- detectCores()- 1
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    trees <- foreach(i = 1:m, combine = list) %dopar% {
      ct =  emphasis::sim.extinct(tree = tree,pars = pars) # complete tree
      lw = ct$logweight
      return(list(wt=ct$wt,E=ct$E,logweight=ct$logweight,lw=lw,n=ct$n))
    }
    lw = sapply(trees,function(list) list$lw)
    stopCluster(cl)
    lw = lw - max(lw)
    w = exp(lw)
    Rec = trees #esta hay que sacarla
  }
  else{
    Rec = vector(mode = 'list',length = m)
    d = vector(mode='numeric',length = m)
    lw = vector(mode='numeric',length = m)
    rp = vector(mode='numeric',length = m)
    sp = vector(mode='numeric',length = m)
    for(j in 1:m){
      rec = sim.extinct(tree = tree,pars = pars)
      Rec[[j]] = rec
      d[j] = length(rec$wt)
      lw[j] = rec$logweight
      rp[j] = prod(rec$rprob)
      sp[j] = prod(rec$sprob)
    }
    lw = lw - max(lw)
    w = exp(lw)
  }
  return(list(rec = Rec, w=w))
}
# Pilot study
pilot.study <- function(tree,epsilon,m1=10,printprocess=FALSE,init_par=c(1.2,0.3,60),l1=20,parallel = TRUE){
  # pilot study suggested by Chan et. al
  pars = init_par
  M = matrix(ncol = 3,nrow = l1)
  H = matrix(ncol = 3,nrow = l1)
  for(i in 1:l1){
    S = sim.sct(tree,pars,m=m1,printHistD = TRUE,parallel = parallel)
    mle =  mle.st(S = S)
    pars = mle$par
    H[i,] = try(diag(solve(mle$hessian))/m1)
    M[i,] = pars
    print(pars)
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
      S = sim.sct(tree,M[i,],m=m1,parallel = parallel)
      mle = mle.st(S = S)
      pars = mle$par
      Me[j,] = pars
      Delta[j] = rel.llik(S1 = S$rec,p0 = M[i,], p1 = pars)
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
mcem.tree <- function(tree,p,parallel=TRUE){
  m = p$m
  s1 = p$s1
  sig = 100*s1/m
  tol = 2*sig*sqrt(1/5)
  D = Inf
  k = 1
  print("initializing mcem")
  pars = p$p
  PARS = c(NULL,NULL,NULL)
  H = c(NULL,NULL,NULL)
  Me = p$M
  while(abs(D)>tol){
    S = sim.sct(tree,pars,m = m,parallel = parallel)
    M = mle.st(S = S)
    mle = M$par
    h1 = try(diag(solve(M$hessian))/m)
    if(is.numeric(h1)) H =  rbind(H,h1)
    D = rel.llik(S1 = S$rec,p0 = pars,p1 = mle)
    PARS = rbind(PARS,mle)
    pars = mle
    print(paste("iteration",k,"previous D: ",D, " lambda: ", pars[1]," mu: ", pars[2], "K:", pars[3]))
    k = k+1
  }
  PARS = data.frame(it=1:(dim(Me)[1]+dim(PARS)[1]),lambda = c(Me[,1],PARS[,1]),mu=c(Me[,2],PARS[,2]),K=c(Me[,3],PARS[,3]))
  return(list(pars=pars,PARS=PARS,H=H))
}
###########################

### functions

#negative logLikelihood of a tree
nllik.tree = function(pars,tree){
  b = c(pars[1],(pars[1]-pars[2])/pars[3],pars[2])
  ldt = tree$wt[length(tree$wt)]
  dt = tree$wt[1:(length(tree$wt)-1)]
  E = tree$E
  if(is.null(tree$n)){
    n = c(2,2+cumsum(E)+cumsum(E-1))
    tree$n = n
  }
  lastn = tree$n[length(tree$n)]
  n = tree$n[1:(length(tree$n)-1)]
  sigma = n*(b[1]-b[2]*n + b[3]) #n-dimentional
  lastsigma = lastn*(b[1]-b[2]*lastn + b[3])
  rho = pmax((b[1]-b[2]*n)*E+b[3]*(1-E),0)
  l = -(sum(-sigma*dt+log(rho))-lastsigma*ldt)
  if(min(b)<0){l = Inf}
  return(l)
}

#logLikelihood of a tree (2nd process)
# llik.tree = function(pars,tree){
#   b = c(pars[1],(pars[1]-pars[2])/pars[3],pars[2])
#   ldt = tree$wt[length(tree$wt)]
#   dt = tree$wt[1:(length(tree$wt)-1)]
#   E = tree$E
#   if(is.null(tree$n)){
#     tree$n = c(2,2+cumsum(E)+cumsum(E-1))
#   }
#   lastn = tree$n[length(tree$n)]
#   n = tree$n[1:(length(tree$n)-1)]
#   sigma = n*(b[1]-b[2]*n + b[3]) #n-dimentional
#   lastsigma = lastn*(b[1]-b[2]*lastn + b[3])
#   rho = pmax((b[1]-b[2]*n)*n*E+b[3]*(1-E),0)
#   l = sum(-sigma*dt+log(rho))-lastsigma*ldt
#   if(min(b)<0){l = Inf}
#   return(l)
# }

# negative logLikelihood of a set of trees
nllik.st = function(pars, st){
  m = length(st$rec)
  l = vector(mode = 'numeric',length = m)
  w = vector(mode = 'numeric',length = m)
  for(i in 1:m){
    s = st$rec[[i]]
    w[i] = st$w[i]
    if(w[i]!=0){
      l[i] = nllik.tree(pars,tree=s)
    }else{
      l[i] = 0
    }
  }
  #L = sum(l*w)/m
  w = w/sum(w)
  L = sum(l*w)
  return(L)
}

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

# MLE for a set of trees
mle.st <- function(S,init_par = c(0.5,0.5,100)){
  po = subplex(par = init_par, fn = nllik.st, st=S,hessian = TRUE)
  return(po)
}

# Monte-Carlo sampling / simulation of a set of complete trees
sim.sct <- function(brts,pars,m=10,oc=0){
    no_cores <- detectCores() - oc
    cl <- makeCluster(no_cores)
    registerDoParallel(cl)
    trees <- foreach(i = 1:m, combine = list) %dopar% {
      ct =  emphasis::sim.extinct(brts = brts,pars = pars) #complete tree
      lw = ct$logweight
      return(list(wt=ct$wt,E=ct$E,logweight=ct$logweight,lw=lw))
    }
    stopCluster(cl)
    lw = sapply(trees,function(list) list$lw)
    dim = sapply(trees,function(list) length(list$wt))
    lw = lw - max(lw)
    w = exp(lw)
  return(list(rec = trees, w=w,dim=dim))
}
# Pilot study
pilot.study <- function(brts,epsilon,m1=10,printprocess=FALSE,init_par=c(1.2,0.3,60),l1=20){
  # pilot study suggested by Chan et. al
  pars = init_par
  M = matrix(ncol = 3,nrow = l1)
  H = matrix(ncol = 3,nrow = l1)
  for(i in 1:l1){
    S = sim.sct(brts,pars,m=m1)
    mle =  mle.st(S = S)
    pars = mle$par
    H[i,] = try(diag(solve(mle$hessian))/m1)
    M[i,] = pars
    print(paste('Q:',mle$value,'pars:',pars[1],pars[2],pars[3]))
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
###########################
sampprob <- function(t,s,mu,r,N){  ## equation (8)
  term1 = (s/N)*(1-exp(-mu*(r-t)))
  c = exp(-(s/mu)*exp(-mu*r))
  if(c==0){
    term2 = 0
  }else{
    term2 = c^{-exp(mu*t)+1}
  }
  f = term1*exp(-s*t)*term2
  return(f)
}

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
###  simulation of extincted new version
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
  to = NULL
  cbt = 0
  N = 2
  sprob = NULL # sampling probability of Missing|observed
  h = 1 # index to fill probabilities
  for(i in 1:dim){
    cwt = wt[i]
    cbt = sum(wt[0:(i-1)])
    key = 0
    gosttime = 0
    #last = 'nothing'
    while(key == 0){
      if(model == "dd"){  # diversity-dependence model
        lambda = max(1e-99, lambda0 - (lambda0-mu0)*N/K)
        mu = mu0
        s = N*lambda
      }else{print('Model not implemented yet, try dd')}
      t.spe = rexp(1,s)
      t.ext = extinction.processes(u=me,inits=ms,mu0=mu0)
      #sometimes parameters does not make sense. write a warning when that happens
      t_ext = ifelse(length(t.ext)>0,min(t.ext),Inf)-cbt  # if is not empty gives the waiting time for next extinction
      mint = min(t.spe,t_ext)
      if(mint < cwt){
        if(mint == t.spe){#speciation
          u = runif(1)
          if(u < pexp(ct-(cbt+t.spe),mu)){
            ms = c(ms,cbt+t.spe)
            me = c(me,u)
            bt = c(bt,cbt+t.spe)
            to = c(to,1)
            sprob[h] = sampprob(t = t.spe+gosttime, s = s, mu = mu, r = ct-(cbt-gosttime),N=N)
            h = h + 1
            N = N + 1
          }else{gosttime = t.spe + gosttime}
          cwt = cwt - t.spe
          cbt = cbt + t.spe
        }
        else{#extinction
          extinctone = which(t.ext == min(t.ext))
          tspe = ms[extinctone]
          text = t.ext[extinctone]
          bt = c(bt,text)
          to = c(to,0)
          sprob[h] = truncdist::dtrunc(text-tspe,'exp',a=0,b=ct-tspe,rate=mu)*(1-integrate(sampprob,lower = 0, upper = t_ext+gosttime,s=s,mu=mu,r=ct-(cbt-gosttime),N=N)$value)#I dont need integrate
          ms = ms[-extinctone]
          me = me[-extinctone]
          cwt = cwt - mint
          cbt = cbt + mint
          N = N-1
          h = h+1
          gosttime = 0
        }
      }
      else{
        key = 1
        sprob[h] = (1 - integrate(Vectorize(sampprob),lower = 0, upper = cwt+gosttime,s=s,mu=mu,r=ct-cbt,N=N)$value) # test if there is difference using vectorize
        h = h+1
      }
    }
    N = N+1
  }
  df = data.frame(bt = c(bt,ct-brts),to = c(to,rep(2,length(wt))))
  df = df[order(df$bt),]
  n.tree = list(wt=c(diff(df$bt),ct-df$bt[length(df$bt)]),E=df$to[-length(df$to)])
  if(length(n.tree$E==1) != length(n.tree$E==0)) print('algo mal!!')
  n.tree$E[n.tree$E==2] = 1
  lrprob = -nllik.tree(pars,n.tree) #f
  lsprob = sum(log(sprob)) #g
  logweight = lrprob-lsprob
  n.tree$weight = exp(logweight)
  n.tree$logweight = logweight
  n.tree$f=lrprob
  n.tree$g=lsprob
  return(n.tree)
}


post.pro <-function(file,extrafile=NULL){
load(file)
#pars = DDD::dd_ML(brts = btdd, idparsopt = 1:3,soc=2,cond=0)
MLE = p$MLE
M = mcem$PARS

it = NULL
lambda = NULL
mu = NULL
K = NULL
for( i in 1:10){
  MM = MLE[[i]]
  it = c(it,rep(i+1,10))
  lambda = c(lambda,MM[,1])
  mu = c(mu,MM[,2])
  K = c(K,MM[,3])
};
col = rep('grey',length(lambda))
SDl = rep(NaN,length(lambda))
SDm = rep(NaN,length(mu))
SDk = rep(NaN,length(K))

it = c(it,1:dim(M)[1])
lambda = c(lambda,M[,2])
mu = c(mu,M[,3])
K = c(K,M[,4])
col = c(col,rep('blue',dim(M)[1]))

MCEMc = data.frame(it=it,lambda=lambda,mu=mu,K=K,col=col)

gamLambda = gam(lambda ~ s(it), data=MCEMc)
gamMu = gam(mu ~ s(it), data=MCEMc)
gamK = gam(K ~ s(it), data=MCEMc)

ex = dim(MCEMc)[1] - dim(mcem$H)[1] - 100
hessL = c(rep(NaN,ex),mcem$H[,1])
SDl = c(SDl,sqrt(hessL+gamLambda$sig2))

hessM = c(rep(NaN,ex),mcem$H[,2])
SDm = c(SDm,sqrt(hessM+gamMu$sig2))

hessK = c(rep(NaN,ex),mcem$H[,3])
SDk = c(SDk,sqrt(hessK+gamK$sig2))

MCEMc$SDl = SDl
MCEMc$SDm = SDm
MCEMc$SDk = SDk

MCEM = rbind(MCEMc,data.frame(it=(-9:0),lambda=p$PM[,1],mu=p$PM[,2],K=p$PM[,3],col=rep('blue',10),SDl=rep(NaN,10),SDm=rep(NaN,10),SDk=rep(NaN,10)))

gl=ggplot(MCEM) + geom_point(aes(it,lambda),colour=MCEM$col) + geom_errorbar(aes(x=it, y=lambda, ymin = lambda-1.96*SDl, ymax = lambda + 1.96*SDl), colour='darkgreen') + geom_hline(yintercept = pars$lambda) + ggtitle(d.name)+theme(axis.text=element_text(size=18),axis.title=element_text(size=16))+labs(x='EM iteration',y=expression(lambda[0]))
gm=ggplot(MCEM) + geom_point(aes(it,mu),colour=MCEM$col) + geom_errorbar(aes(x=it, y=mu, ymin = mu-1.96*SDm, ymax = mu + 1.96*SDm), colour='darkgreen') + geom_hline(yintercept = pars$mu)+ ggtitle(d.name)+theme(axis.text=element_text(size=18),axis.title=element_text(size=16))+labs(x='EM iteration',y=expression(mu[0]))
gk=ggplot(MCEM) + geom_point(aes(it,K),colour=MCEM$col) + geom_errorbar(aes(x=it, y=K, ymin = K-1.96*SDk, ymax = K + 1.96*SDk), colour='darkgreen') + geom_hline(yintercept = pars$K) + ggtitle(d.name)+theme(axis.text=element_text(size=18),axis.title=element_text(size=16))+labs(x='EM iteration',y=expression(K))
gLLmcem = ggplot(data=d)+geom_line(aes(x=it,y=llik))+geom_hline(yintercept = pars$loglik)
gLLmcem2 = ggplot(df, aes(lambda, mu))+ geom_contour(aes(z = llik,color=..level..),bins=100) +geom_point(data=mcem$PARS,aes(x=lambda,y=mu,size=it))+  ggtitle(d.name)
if(!is.null(extrafile)){
  load(extrafile)
  gLLmcem2 = gLLmcem2 + geom_point(data=mcem$PARS,aes(x=lambda,y=mu,size=it),col='red')
  gLLmcem = gLLmcem + geom_line(data=d,aes(x=it,y=llik),col='red')
}
return(list(gl=gl,gm=gm,gk=gk, gLLmcem=gLLmcem, gLLmcem3d=gLLmcem2))
}

get.time <- function(time,mode='sec'){
  dif = proc.time()-time
  ti = as.numeric(dif[3])
  if(mode == 'min')  ti = ti/60
  if(mode == 'hou') ti = ti/3600
  return(ti)
}


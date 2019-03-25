### remaining EMPHASIS functions


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



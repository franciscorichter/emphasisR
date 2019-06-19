# more utilities
lik.tree <- function(pars,tree,model="dd",initspec=1){
  exp(-nllik.tree(pars,tree,model=model,initspec=initspec))
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

# stochastic function

mc.llik <- function(pars,brts,m){
  S = sim.sct(brts,pars,m,print=FALSE)
  aplli = log(mean(S$w))
  return(aplli)
}

# subset of trees
sub.st <- function(S,ntrees=10){
  qu = quantile(S$w,1-ntrees/length(S$w))
  whi = which(S$w>qu)
  length(whi)
  newtrees = vector(mode="list",length=length(whi))
  i=1
  for(j in whi){
    newtrees[[i]] = S$trees[[j]]
    i = i+1
  }
  newS = list(w=S$w[whi],dim=S$dim[whi],nl=S$nl[whi],trees=newtrees,g=S$g[whi],fsm=S$fsm[whi],fe=S$fe[whi])
  return(newS)
}



phylo2tree <- function(tree){
  # to map newick trees into ther xxxx format
  ltt = ltt.plot.coords(tree)
  t = diff(ltt[,1])
  ltt = ltt[-1,]
  n = ltt[,2]
  E = diff(n)
  E[E==-1] = 0
  return(list(wt=t,to=E))
}

vectors2phylo <- function(list,initspec=1){
  wt=list$wt
  to=list$to
  to[to==2] = 1
  ct=sum(wt)
  newick = paste(sl[1],";",sep="")
  N=1
  identf = data.frame(Spec="a",Time=0) # Labels of species
  for (i in 1:(length(wt)-1)){
    # speciation
    bt = sum(wt[1:i])
    BD = sample(1:N,1)
    species = as.character(identf[BD,1])  
    if (to[i] == 1){
      ind = regexpr(species,newick)[1]-1
      atm = bt-identf[which(identf[,1]==species),2]
      newick = paste(substr(newick,1,ind),"(",substr(newick,ind+1,ind+4),",",sl[i+1],"):",as.character(atm),substring(newick,ind+5),sep="")
      identf = rbind(identf,data.frame(Spec=substr(sl[i+1],1,2),Time=bt))
      identf[identf$Spec == species,2] = bt
      N = N+1
    }
    # extinction
    if (to[i]==0){
      ind = regexpr(species,newick)[1] + 2
      atm = bt-identf[which(identf[,1]==species),2]
      identf = identf[!identf$Spec==species,]
      newick = paste(substr(newick,1,ind),as.character(atm),substring(newick,ind+2),sep="")
      N=N-1
    }
  }
  newick = compphyl(newi=newick,identf=identf,ct=ct)
  newick = read.tree(text=newick)
  return(newick)
}

compphyl <- function(newi,identf,ct){
  #set to extant species to the present time
  identf[,1] = as.character(identf[,1])
  identf[,2] = ct-identf[,2]
  for(i in 1:length(identf[,1])){
    ind = regexpr(identf[i,1],newi)[1] + 2
    newi = paste(substr(newi,1,ind),as.character(identf[i,2]),substring(newi,ind+2),sep="")
  }
  return(newi)
}

sl = paste(letters[1],letters,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],letters,":0",sep="")
  sl = c(sl,ll)
}

#post processing
post.pro <-function(file,extrafile=NULL){
  load(file)
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

# time calculation
get.time <- function(time,mode='sec'){
  dif = proc.time()-time
  ti = as.numeric(dif[3])
  if(mode == 'min')  ti = ti/60
  if(mode == 'hou') ti = ti/3600
  return(ti)
}



W <- function (z, branch = 0)
{
  stopifnot(length(branch) == 1, is.numeric(z))
  if (anyNA(z)) {
    warning("Some values of ", deparse(substitute(z)), " are NA or NaN. ",
            "Returning 'NA' for these entries.")
    non.na.z <- z[!is.na(z)]
  }
  else {
    non.na.z <- z
  }
  W.non.na.z <- rep(NA, length(non.na.z))
  if (branch == 0) {
    W.non.na.z <- lamW::lambertW0_C(non.na.z)
  }
  else if (branch == -1) {
    if (any(is.infinite(z))) {
      warning("'Inf' is not a valid argument of the non-principal branch W",
              " (branch = -1).")
    }
    W.non.na.z <- lamW::lambertWm1_C(non.na.z)
  }
  else {
    stop("Branch was ", branch, "; must be either '0' or '-1'.")
  }
  if (length(W.non.na.z) == length(z)) {
    dim(W.non.na.z) <- dim(z)
    return(W.non.na.z)
  }
  else {
    W.z <- rep(NA, length(z))
    W.z[!is.na(z)] <- W.non.na.z
    W.z[is.nan(W.z)] <- NA
    dim(W.z) <- dim(z)
    return(W.z)
  }
}


get.topologies <- function(M){
  if(M == 0)
  {
    return(NULL)
  }
  TO = matrix(nrow=2*M,ncol=1)
  TO[1,1] = 1
  for(i in 2:(2*M)){
    comb = ncol(TO)
    for(j in 1:comb){
      ns = sum(no.na(TO[,j]))
      ne = sum(1-no.na(TO[,j]))
      if(ns < M & ns > ne){ #extinction or speciation
        TO[i,j] = 1
        TO = cbind(TO,matrix(TO[,j],ncol=1))
        TO[i,ncol(TO)] = 0
      }
      if(ns == M & ns > ne){ #extinction
        TO[i,j] = 0
      }
      if(ns < M & ns == ne){ #speciation
        TO[i,j] = 1
      }
    }
  }
  return(TO)
}



norma.const <- function(lambda,mu,CT,n){
  lg = NULL
  lf = NULL
  for(i in 1:n){
    df = sim.tree.g(lambda,mu,CT)
    lg[i] = g_samp_missobs(df,c(lambda,mu))
    to = df$to[-nrow(df)]
    to[to==2] = 1
    lf[i] = -nllik.tree(pars=c(lambda,mu,Inf),tree=list(wt=diff(c(0,df$bt)),to=to))
  }
  const.approx = sum(exp(lf-lg))/n
  return(const.approx)
}


sim.tree.g <- function(lambda,mu,CT){
  Ne = rpois(1,(lambda+mu)*CT)
  Np = rpois(1,(lambda-mu)*CT)
  brts = runif(Np,min=0,max=CT)
  ms = runif(Ne,min=0,max=CT)
  me = runif(Ne,min=ms,max=CT)
  to = c(rep(1,Ne),rep(0,Ne))
  df = data.frame(bt = c(ms,me),to=to)
  df = rbind(df,data.frame(bt=c(brts,CT),to=c(rep(2,length(brts)),2)))
  df = df[order(df$bt),]
  return(df)
}

g_samp_missobs <- function(df,pars){
  last = nrow(df)
  ct = df[last,1]
  Ne = sum(df$to==0)
  Np = sum(df$to==2) - 1
  subdf_e = df[df$to==1,]
  subdf_p = df[df$to==2,]
  Ts = ct-subdf_e$bt
  to = df$to[-last]
  rto = to
  to[to==2] = 1
  n = c(2,2+cumsum(to)+cumsum(to-1))
  nm = n[c(rto,2)==1]
  log_dens_miss = dpois(Ne,lambda =(pars[1]+pars[2])*ct,log=TRUE)-Ne*log(ct)-sum(log(Ts))-sum(log(nm))+lgamma(Ne+1)
  ne = n[c(rto,0)==2]
  log_dens_obs = dpois(Np,lambda = (pars[1]-pars[2])*ct,log=TRUE)-Np*log(ct)-sum(log(ne))+lgamma(Np+1)
  log_dens = log_dens_miss+log_dens_obs
  return(log_dens)
}


we_cal <- function(tree){ # weights calculation
  list2env(setNames(tree, c("wt","to","n","s","r","pars","t_ext")), .GlobalEnv)
  mu = pars[2]
  if(mu!=0){
    if(length(t_ext) != length(to)) to = c(to,2)
    ev = (to==1 & is.infinite(t_ext))
    nobs = c(2,2+cumsum(ev))
    nobs = nobs[-length(nobs)]
    index = (!is.finite(t_ext) & to == 1)
    logg = sum(-s*(exp(-r*mu)/mu)*(exp(wt*mu)-1))+sum(log(s[index]/n[index]))-sum(mu*nobs*wt)
  }else{
    logg = 0
  }
  return(logg)
}

sim.tree <- function(pars,CT,seed=1,model="dd"){
  set.seed(seed)
  mu = pars[2]
  N = 2
  brts = 0
  to = c(1,1)
  while(max(brts)<CT & N>0){
    if(model=="dd") lambda = lambda.dd(pars,N)
    if(model=="cr") lambda = lambda.cr(pars,N)
    sigma = N*(mu+lambda)
    wt = rexp(1,rate=sigma)
    tev = rbinom(n=1,size=1,prob= (lambda/(lambda+mu)) )
    brts = c(brts,max(brts)+wt)
    to = c(to,tev)
    if(tev==1) N = N+1
    if(tev==0 & max(brts)<CT){
      N = N-1
      ext = sample((1:length(to))[to==1],1)
      to[ext] = -1
    }
  }
  to = to[-length(to)]
  if(max(brts)>CT) brts[length(brts)]=CT
  tree = list(wt=diff(brts),to=to)
  return(tree)
}



mle.tree <- function(tree,init_par= c(1,0.5,100),topology=TRUE,model="dd",truncdim=FALSE){
  if(model=="cr") init_par = c(init_par[1],init_par[2])
  po = subplex(par = init_par, fn = nllik.tree, topology=topology, truncdim=truncdim, tree = tree,model=model,hessian = TRUE)
  return(po)
}


prune.tree <- function(tree){
  brts = sum(tree$wt)-c(0,cumsum(tree$wt))
  brts = c(brts[1],brts[-length(brts)])
  nbrts = brts[tree$to==1]
  return(nbrts)
}

sim.tree.g <- function(lambda,mu,CT,seed){
  set.seed(seed)
  Ne = rpois(1,mu*CT) # number of extinct species
  Np = rpois(1,(mu+lambda)*CT) # number of present species
  brts = runif(Np,min=0,max=CT) # branching times of speciation of present species
  ms = runif(Ne,min=0,max=CT) # branching times of speciation of extinct species
  me = runif(Ne,min=ms,max=CT) # branching times of extinctions
  # return matrix with branching times and topology
  to = c(rep(1,Ne),rep(0,Ne))
  df = data.frame(bt = c(ms,me),to=to)
  df = rbind(df,data.frame(bt=c(brts,CT),to=c(rep(2,length(brts)),2)))
  df = df[order(df$bt),]
  return(df)
}

no.na <- function(vector){
  vector[!is.na(vector)]
}

pars2 = c(2.5e+02,1.0e+00, 0.0e+00, 1.0e+00, 0.0e+00, 1.0e+00)

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


#nee et. al likelihood
#' @title Pt
#' @author Giovanni Laudanno
#' @description Nee's function: pt
#' @inheritParams default_params_doc
#' @return pt
#' @export
pt  <- function (lambda, mu, t) {
  time <- t
  Lambda <- exp((mu - lambda) * time)
  out    <- (lambda == mu) * (1/(1 + lambda * time)) +
    (lambda != mu) * ((lambda - mu + (lambda == mu))/(lambda - mu *
                                                        Lambda * (lambda != mu) + (lambda == mu)))
  return(unname(out))
}

#' @title ut
#' @author Giovanni Laudanno
#' @description Nee's function: ut
#' @inheritParams default_params_doc
#' @return ut
#' @export
ut  <- function (lambda, mu, t) {
  time <- t
  Lambda <- exp((mu - lambda) * time)
  out    <- (lambda == mu) * (lambda * time/(1 + lambda * time)) +
    (lambda != mu) * ((lambda - lambda * Lambda + (lambda == mu)) /
                        (lambda - mu * Lambda * (lambda != mu) + (lambda == mu)))
  return(unname(out))
}

#' @title Pn
#' @author Giovanni Laudanno
#' @description Nee's function: pn
#' @inheritParams default_params_doc
#' @return pn
#' @export
pn <- function(lambda, mu, t, n) {
  out <- (n > 0) * sls::pt(t = t, lambda = lambda, mu = mu) *
    (1 - sls::ut(t = t, lambda = lambda, mu = mu)) *
    sls::ut(t = t, lambda = lambda, mu = mu)^(n - 1 + 2*(n == 0)) +
    (n == 0) * (1 - sls::pt(t = t, lambda = lambda, mu = mu))
  return(out)
}

Nee_likelihood <- function(lambda, mu, brts, cond)
{
  BRTS <- c(brts[1], brts)
  NN <- length(BRTS)
  out <- prod(
    pt(lambda = lambda, mu = mu, t = BRTS) * (1 - ut(lambda = lambda, mu
                                                     = mu, t = BRTS))
  ) * lambda^(NN - 1) #* factorial(NN - 1)
  out <- out/((pt(lambda = lambda, mu = mu, t = BRTS[1]))^(2 * cond))
}

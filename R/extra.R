# more utilities


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


sim.sct.alt <- function(brts,pars){
  no_cores <- detectCores()
  cl <- makeCluster(no_cores)
  registerDoParallel(cl)
  trees <- foreach(i = 1:m, combine = list) %dopar% {
    df =  emphasis::sim.extinct2(brts = brts,pars = pars,model=model)
    tree = emphasis::df2tree(df,pars,model=model)
    lsprob = emphasis::lg_prob(tree,topology = topology)
    nl = emphasis::nllik.tree(pars,tree=tree,topology = topology,model=model)
    lw = -nl-lsprob#-DDD:::dd_loglik(pars1 = pars, pars2 = pars2,brts = brts, missnumspec = 0)
    fms = df$bt[is.finite(df$bte)][1] #first missing speciation
    fe = df$bt[df$to==0][1] # first extinction
    return(list(nl=nl,lw=lw,tree=tree,lsprob=lsprob,fms = fms,fe=fe))
  }
  stopCluster(cl)
  lw = sapply(trees,function(list) list$lw)
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
  return(list(w=w, dim=dim, nl=nl, trees=trees,g=eg,fms=fms,fe=fe))
}
#post processing
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


# more utilities


# stochastic function

mc.llik <- function(pars,brts,m){
  S = sim.sct(brts,pars,m,print=FALSE)
  aplli = log(mean(S$w))
  return(aplli)
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

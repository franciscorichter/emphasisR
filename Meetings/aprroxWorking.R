brts = 1
n.gsim <- 1000
pars = c(0.2,0.1,Inf)

nee1=Nee_likelihood(lambda = pars[1], mu = pars[2], brts = brts,cond = 0)
nee1
pt(lambda = pars[1], mu = pars[2],t = 1) *(1-ut(lambda = pars[1], mu = pars[2],t = 1))
exp(DDD:::dd_loglik(pars1 = pars, pars2 = pars2,brts = brts, missnumspec = 0))
pars2[6]=1
exp(DDD:::dd_loglik(pars1 = pars, pars2 = pars2,brts = brts, missnumspec = 0))

fhat1(obs = brts,nsim = 1000000,maxnumspec = 10)


S = emphasis:::sim.sct(brts,pars,m = n.gsim,print = FALSE,initspec = 1)
mean(S$w)

S = emphasis:::sim.sct(brts,pars,m = n.gsim,print = FALSE,initspec = 1,topology = F)
mean(S$w)

ctezero = exp(-nllik.tree(pars,tree=list(wt=diff(c(0,obt)),to=rep(1,length(obt)-1)),initspec = 1,topology = F))
f.hat4(brts = brts,pars = pars,cte = 1,m = 10000, model = "cr") + ctezero


sims = lapply(re, sim.mis01,pars=pars)
listis = separate(sims)
l0=listis[[1]]
l1=listis[[2]]

mean(sapply(l0, weight01,pars=pars))*(length(l0)/(length(l0)+length(l1)))
ctezero

mean(sapply(l1, weight01,pars=pars))*(length(l1)/(length(l0)+length(l1)))
f.hat4(brts = brts,pars = pars,cte = cte-2,m = 10000, model = "cr")

###


########################

## new method 1

sim.miss1 <- function(maxnumspec=250,ct){
  S = sample(0:maxnumspec,1)
  #S=1
  brts = sort(runif(2*S,min=0,max=ct))
  to = sampletopology(S)
  tree = list(brts=brts,to=to)
  return(tree)
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

gprob1 <- function(to,maxnumspec,ct){
  gprob = factorial(length(to))*((1/ct)^(length(to)))*probto(to)
}

probto <- function(to,p=0.5){
  posspec = c(0,cumsum(to==1))<(length(to)/2)
  posext = !(c(0,cumsum(to==1))==c(0,cumsum(to==0)))
  expo = sum(posspec & posext)
  prob = p^expo
  return(prob)
}

fhat1 <- function(obs,pars,nsim=10000,maxnumspec=250,model="cr"){
  ct = max(obs)
  f.joint = vector(mode="numeric",length=nsim)
  g.samp = vector(mode="numeric",length=nsim)
  dim = vector(mode="numeric",length=nsim)
  for(i in 1:nsim){
    miss = sim.miss1(maxnumspec = maxnumspec,ct = ct)
    dim[i] = length(miss$brts)
    g.samp[i] = gprob1(to = miss$to,maxnumspec = maxnumspec,ct=ct)
    df = data.frame(bt=c(miss$brts,obs),to=c(miss$to,rep(1,length(obs))))
    df = df[order(df$bt),]
    tree = df2tree2(df)
    f.joint[i] = lik.tree(pars=pars,tree=tree,topology = F,model=model,initspec = 1)
  }
  fhat = mean(f.joint/g.samp) # pending, better on log to avoid numerical issues
  return(fhat)
}

## old method working (Method 4)

Xd.mc <- function(d,obt,nsim=10,pars,model="dd"){
  #print(paste("Getting topologies for", d," missing species"))
  top = get.topologies(d)
  sumofllik = vector(mode="numeric",length=ncol(top))
  #print(paste("Drawing uniform trees for ",ncol(top),"topologies"))
  for(j in 1:ncol(top)){
    topi = top[,j]
    topi = replicate(nsim,topi,simplify = F)
    ntree = lapply(topi,sim.unif.tree.df,obt=obt)  # simulate unif random tree given topology
    nlliks = lapply(ntree,weight.unif,pars=pars,d=d,ct=max(obt),model=model)
    sumofllik[j] = Reduce("+", nlliks)
  }
  sumofllik = sum(sumofllik)
  return(sumofllik/(ncol(top)*nsim))
}


weight.unif <- function(tree,pars,d,ct,model){
  lik.tree(pars,tree,initspec = 1,topology = F,model = model)*((catalan(d)*(ct^(2*d)))/(factorial(2*d)))
}


sim.unif.mtree <- function(top,obt){
  nmbt = length(top) #length of missing bt
  mbt = sort(runif(nmbt,min=0,max=max(obt)))
  df = data.frame(bt=c(mbt,obt),to=c(top,rep(2,length(obt))))
  df = df[order(df$bt),]
  return(df)
}

sim.unif.tree.df <- function(top,obt){
  df = df2tree2(sim.unif.mtree(top,obt))
  return(df)
}

f.hat4 <- function(brts,pars,cte=10,m=100,model="dd"){
  int.d = NULL
  for(i in 1:cte){
    int.d[i] = Xd.mc(nsim=m,obt=brts,d = i,pars=pars)
  }
  sum(int.d)
}

## Method 4.1

Xd.mc <- function(d,obt,nsim=10,pars,model="dd"){
  #print(paste("Getting topologies for", d," missing species"))
  top = get.topologies(d)
  sumofllik = vector(mode="numeric",length=ncol(top))
  #print(paste("Drawing uniform trees for ",ncol(top),"topologies"))
  for(j in 1:ncol(top)){
    topi = top[,j]
    topi = replicate(nsim,topi,simplify = F)
    ntree = lapply(topi,sim.unif.tree.df,obt=obt)  # simulate unif random tree given topology
    nlliks = lapply(ntree,weight.unif,pars=pars,d=d,ct=max(obt),model=model)
    sumofllik[j] = Reduce("+", nlliks)
  }
  sumofllik = sum(sumofllik)
  return(sumofllik/(ncol(top)*nsim))
}


weight.unif <- function(tree,pars,d,ct,model){
  lik.tree(pars,tree,initspec = 1,topology = F,model = model)*((catalan(d)*(ct^(2*d)))/(factorial(2*d)))
}


sim.unif.mtree <- function(top,obt){
  nmbt = length(top) #length of missing bt
  mbt = sort(runif(nmbt,min=0,max=max(obt)))
  df = data.frame(bt=c(mbt,obt),to=c(top,rep(2,length(obt))))
  df = df[order(df$bt),]
  return(df)
}

sim.unif.tree.df <- function(top,obt){
  df = df2tree2(sim.unif.mtree(top,obt))
  return(df)
}

f.hat4 <- function(brts,pars,cte=10,m=100,model="dd"){
  int.d = NULL
  for(i in 1:cte){
    int.d[i] = Xd.mc(nsim=m,obt=brts,d = i,pars=pars)
  }
  sum(int.d)
}

##############################

f.hat01 <- function(pars,obs,nsim=100){
  re = replicate(nsim,obs,simplify = F)
  sims = lapply(re, sim.mis01,pars=pars)
  weights = sapply(sims, weight01,pars=pars)
  return(mean(weights))
}


weight01 <- function(tree,pars,model="cr"){
  lik.tree(pars,tree,initspec = 1,topology = F,model = model)/gprob(tree,pars)
}

sim.mis01 <- function(pars,obt,initspec=1){
  ct = max(obt)
  N = initspec
  s = lambda.cr(pars,N)
  tspe = rexp(1,s)
  mu=pars[2]
  if(tspe<ct){
    text = truncdist::rtrunc(1,"exp",a = tspe, b = ct,rate=mu)
    wt = c(tspe,text-tspe,ct-text)
    to = c(1,0)
  }else{
    wt = obt
    to = rep(1,0)
  }
  tree = list(wt=wt,to=to)
  return(tree)
}

gprob <- function(tree,pars){
  lambda=pars[1]
  mu=pars[2]
  ct = sum(tree$wt)
  if(length(tree$wt)==1){
    gprob = exp(-lambda*ct)
  }else{
    gprob = lambda*exp(-tree$wt[1]*lambda)*truncdist::dtrunc(tree$wt[1]+tree$wt[2],"exp",a = tree$wt[1], b = ct,rate=mu)
  }
  return(gprob)
}

separate <- function(list){
  j = 1
  k = 1
  l0 = list()
  l1 = list()
  for(i in 1:length(list)){
    if(length(list[[i]]$wt)==1){
      l0[[j]] = list[[i]]
      j=j+1
    }else{
      l1[[k]] = list[[i]]
      k=k+1
    }
  }
  return(list(l0=l0,l1=l1))
}

########################

Xd.mc <- function(d,obt,nsim=10,pars,model="dd"){
  #print(paste("Getting topologies for", d," missing species"))
  top = get.topologies(d)
  sumofllik = vector(mode="numeric",length=ncol(top))
  #print(paste("Drawing uniform trees for ",ncol(top),"topologies"))
  for(j in 1:ncol(top)){
    topi = top[,j]
    topi = replicate(nsim,topi,simplify = F)
    ntree = lapply(topi,sim.unif.tree.df,obt=obt)  # simulate unif random tree given topology
    nlliks = lapply(ntree,weight.unif,pars=pars,d=d,ct=max(obt),model=model)
    sumofllik[j] = Reduce("+", nlliks)
  }
  sumofllik = sum(sumofllik)
  return(sumofllik/(ncol(top)*nsim))
}


weight.unif <- function(tree,pars,d,ct,model){
  lik.tree(pars,tree,initspec = 1,topology = F,model = model)*((catalan(d)*(ct^(2*d)))/(factorial(2*d)))
}


sim.unif.mtree <- function(top,obt){
  nmbt = length(top) #length of missing bt
  mbt = sort(runif(nmbt,min=0,max=max(obt)))
  df = data.frame(bt=c(mbt,obt),to=c(top,rep(2,length(obt))))
  df = df[order(df$bt),]
  return(df)
}

sim.unif.tree.df <- function(top,obt){
  df = df2tree2(sim.unif.mtree(top,obt))
  return(df)
}


f.hat4 <- function(brts,pars,cte=10,m=100,model="dd"){
  int.d = NULL
  for(i in 1:cte){
    int.d[i] = Xd.mc(nsim=m,obt=brts,d = i,pars=pars)
  }
  sum(int.d)
}
#######

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
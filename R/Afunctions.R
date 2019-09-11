# Alternative function for testing



proper_weighting <- function(E){
  ta = table(E$dim)
  dims = as.numeric(names(ta))
  sizes = as.numeric(ta)
  weights = NULL
  for(i in 1:length(dims)){
    weights[i] = sum(E$weights[E$dim == dims[i]])
  }
  W = weights
  w_tilde = weights/length(E$trees)
  Z = weights/sizes
  #fhat2[j] = mean(W)
  return(list(W=W,w_tilde=w_tilde,Z=Z))
}


sim.sct2 <- function(brts,pars,m=10,oc=0){
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


sampprob <- function(t,s,mu,r,N){  ## equation (8)
  term1 = (s)*(1-exp(-mu*(r-t)))
  c = exp(-(s/mu)*exp(-mu*r))
  if(c==0){
    term2 = 0
  }else{
    term2 = c^{-exp(mu*t)+1}
  }
  f = term1*exp(-s*t)*term2
  return(f)
}

#weight
logweight <- function(pars,df,ct=NULL){
  dim = dim(df)[1]
  ct = df$bt[dim]
  #  xi = rep(0,times=dim)
  # xi[df$bte<(ct+1)] = 1
  wtT = diff(c(0,df$bt))
  n.tree = list(wt=wtT,E=df$to[-dim])
  n.tree$E[n.tree$E==2] = 1
  E = n.tree$E
  n = c(2,2+cumsum(E)+cumsum(E-1))
  lambda = (pars[1]-(pars[1]-pars[2])*(n/pars[3]))
  mu = pars[2]
  s = lambda*n
  #sprob = da.prob(xi=xi,wt=n.tree$wt,t_ext=df$t.ext,s=s,mu=pars[2],r=ct-c(0,df$bt[-length(df$bt)]),n=n)
  lsprob = da.prob2(wt=wtT,t_ext=df$t.ext,s=s,mu=pars[2],r=ct-c(0,df$bt[-length(df$bt)]),n=n)
  lrprob = -nllik.tree(pars,n.tree)
  logweight = lrprob-lsprob
  return(logweight)
}

prob.ms <- function(wt,t_ext,s,mu,r,n){
  if(mu != 0){
    prob = (s/n)*mu*exp(-s*(wt+(exp(-mu*r)-exp(-mu*(r-wt)))/mu)-mu*t_ext)
  }else{
    prob = 0
  }
  return(prob)
}
prob.nospecies <- function(wt,s,mu,r){
  if(mu != 0){
    prob = exp(-s*(wt+(exp(-mu*r)-exp(-mu*(r-wt)))/mu))
  }else{
    prob = 1
  }
  return(prob)
}
da.prob <- function(xi,wt,t_ext,s,mu,r,n){
  g =ifelse(xi,prob.ms(wt,t_ext,s,mu,r,n),prob.nospecies(wt,s,mu,r))
  return(g)
}


## cool new way

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

sim.extinct2 <- function(brts,pars,model='dd',seed=0){
  if(seed>0) set.seed(seed)
  wt = diff(c(0,brts))
  ct = sum(wt)
  lambda0 = pars[1]
  mu0 = mu = pars[2] #mu can be removed
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
        lambda = lambda.dd(pars,N)
      }
      if(model == "dd1.3"){
        lambda = lambda.dd.1.3(pars,N)
      }
      s = N*lambda
      t.spe = rexp(1,s)
      t.ext = extinction.processes(u=me,inits=ms,mu0=mu0)
      t_ext = ifelse(length(t.ext)>0,min(t.ext),Inf)-cbt  # if is not empty gives the waiting time for next extinction
      mint = min(t.spe,t_ext)
      if(mint < cwt){
        if(mint == t.spe){#speciation
          u = runif(1)
          if(u < 2*pexp(ct-(cbt+t.spe),mu0)*(1-pexp(ct-(cbt+t.spe),mu0))){
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
  df = data.frame(bt = c(bt,brts),bte = c(bte, rep(Inf,length(wt))),to = c(to,rep(2,length(wt))))
  df = df[order(df$bt),]
  df$t.ext = df$bte-df$bt
  df = df[-nrow(df),]
  df = rbind(df,data.frame(bt=ct,bte=Inf,to=2,t.ext=Inf))
  return(df)
}

##


###  simulation of extincted old version
sim.extinct_old <- function(brts,pars,model='dd',seed=0){
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
  Li = list()
  cbt = 0
  N = 2
  sprob = NULL # sampling probability of Missing|observed
  h = 1 # index to fill probabilities
  for(i in 1:dim){
    cwt = wt[i]
    cbt = sum(wt[0:(i-1)])
    key = 0
    gosttime = 0
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
            bt = c(bt,cbt+t.spe) #SHOULD INCLUDE GOSTTIME NO?
            text = extinction.processes(u=u,inits=cbt+t.spe,mu0=mu0)
            bte = c(bte,text)
            to = c(to,1)
            tspe = cbt+t.spe
            sprob[h] = sampprob(t = t.spe+gosttime, s = s, mu = mu, r = ct-(cbt-gosttime),N=N)*truncdist::dtrunc(text-tspe,'exp',a=0,b=ct-tspe,rate=mu)*(lambda/s)
            Li[[h]] = list(xi='speciation',text=text,wt= t.spe+gosttime,gosttime=gosttime,s=s, r = ct-(cbt-gosttime),N=N)
            if(all.equal(sprob[h],prob.ms(t.spe+gosttime,text-tspe,s,mu,r =  ct-(cbt-gosttime),n = N))){
              #  print('speciacion ms')
            }else{
              stop('problem in speciation')
            }
            h = h + 1
            N = N + 1
            gosttime = 0
          }else{gosttime = t.spe + gosttime}
          cwt = cwt - t.spe
          cbt = cbt + t.spe
        }
        else{#extinction
          extinctone = which(t.ext == min(t.ext))
          tspe = ms[extinctone]
          text = t.ext[extinctone]
          bt = c(bt,text)
          bte = c(bte,Inf)
          to = c(to,0)
          sprob[h] = (1-integrate(sampprob,lower = 0, upper = t_ext+gosttime,s=s,mu=mu,r=ct-(cbt-gosttime),N=N)$value)#I dont need integrate
          if(all.equal(sprob[h],prob.nospecies(t_ext+gosttime,s=s,mu=mu,r=ct-(cbt-gosttime)))){
          }else{
            stop('problem in ext')
          }
          Li[[h]] = list(xi='extinction', wt= t_ext+gosttime,gosttime=gosttime,s=s, r = ct-(cbt-gosttime),N=N)
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
        sprob[h] = (1 - integrate(Vectorize(sampprob),lower = 0, upper = cwt+gosttime,s=s,mu=mu,r=ct-(cbt-gosttime),N=N)$value) # test if there is difference using vectorize
        if(all.equal(sprob[h],prob.nospecies(cwt+gosttime,s=s,mu=mu,r=ct-(cbt-gosttime)))){
          #  print('paso nada ms')
        }else{
          stop('problem in nada')
        }
        Li[[h]] = list(xi='nada',wt=cwt+gosttime,gosttime=gosttime,s=s, r = ct-(cbt-gosttime),N=N)
        h = h+1
      }
    }
    N = N+1
  }
  df = data.frame(bt = c(bt,ct-brts),bte = c(bte, rep(Inf,length(wt))),to = c(to,rep(2,length(wt))))
  df = df[order(df$bt),]
  df$t.ext = df$bte-df$bt
  df$xi = 0
  df$xi[df$bte<ct] = 1
  wtT = c(diff(df$bt),ct-df$bt[length(df$bt)])
  df = df[-1,]
  n.tree = list(wt=wtT,E=df$to)
  if(length(n.tree$E==1) != length(n.tree$E==0)) print('algo mal!!')
  n.tree$E[n.tree$E==2] = 1
  E = n.tree$E
  n = c(2,2+cumsum(E)+cumsum(E-1))
  lambda = (pars[1]-(pars[1]-pars[2])*(n/pars[3]))
  mu = pars[2]
  s = lambda*n
  ssprob = da.prob(xi=c(df$xi,0),wt=n.tree$wt,t_ext=c(df$t.ext,Inf),s=s,mu=mu0,r=c(ct,ct-df$bt),n=n)
  lrprob = -nllik.tree(pars,n.tree) #f
  lsprob = sum(log(sprob)) #g
  logweight = lrprob-lsprob
  lssprob = sum(log(ssprob))
  logweight2 = lrprob-lssprob
  if(logweight==Inf) logweight = -Inf
  n.tree$weight = exp(logweight)
  n.tree$weight2 = exp(logweight2)
  n.tree$logweight = logweight
  n.tree$logweight2 = logweight2
  n.tree$f=lrprob
  n.tree$g=lsprob
  return(n.tree)
}


##############


sim.missing <- function(pars,obs,initspec = 1){
  ct = max(obs)
  N = initspec
  wt = diff(c(0,obs))
  i=1
  bt = 0
  btm = NULL
  while(i < length(wt)){
    la = lambda.cr(pars,N)
    tspe = rexp(1,la)
    if(tspe < wt[i]){
      btm = c(bt,bt + tspe)
      text = truncdist::rtrunc(1,"exp",a = cbt+t.spe, b = ct,rate=mu)
      tom = c(tom,1)
    }else{
      i = i + 1
    }
    
  }
}
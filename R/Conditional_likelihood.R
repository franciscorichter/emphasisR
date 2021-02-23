# load package to fit gam
#library(mgcv)
# Simulate survival
create_grid <- function(llim,rlim,n.grid){
  p = length(llim)
  theta.range = cbind(llim,rlim)
  pars = NULL
  for (i in 1:p){
    v = rep(rep(seq(theta.range[i,1],
                    theta.range[i,2],
                    length.out = n.grid),
                each=n.grid^(p-i)),n.grid^(i-1))
    pars = cbind(pars,v)
  }
  return(pars)
}


Simulation_step <- function(grid,model,ct,printPars=T){
  #then simulate trees
  srv = n = tm = vector(mode="numeric",length=nrow(grid))
  for (i in 1:nrow(grid)){
    if(printPars)   print(i)
    time = proc.time()
    tau = try(sim_survival(diversification_model = list(pars=grid[i,],model=model),
                           ct=ct,
                           timeLimit = 2),silent = TRUE)
    tm[i] = emphasis::get.time(time)
    if(class(tau)=="try-error"){
      srv[i] = -1
      n[i] = -1
    }else{
      srv[i] = tau$srv
      n[i] = length(tau$tree$brts)
    }
  }
  ##
  df = data.frame(srv=srv,n=n,time=tm,mu=grid[,1],lambda=grid[,2],betaN=grid[,3],betaP=grid[,4])
  df$srv[df$srv==-1] = "error" 
  df$srv[df$srv==0] = "Not survived"
  df$srv[df$srv==1] = "Extant tree"
  names(df) = c("Output simulation","Number of branching times","time","mu","lambda","betaN","betaP")
  dat = df
  dat$srv = 1
  dat$srv[dat$`Output simulation`=="Not survived"] = 0
  ##
  g1 = ggplot(df)+geom_point(aes(x=mu,y=lambda,size=`Number of branching times`,colour=`Output simulation`))+theme_bw()
  g2 = ggplot(df)+geom_point(aes(x=betaN,y=betaP,size=`Number of branching times`,colour=`Output simulation`))+theme_bw()
  return(list(df=dat,g1=g1,g2=g2))
}


fit_gam_survival <- function(ct,grid,llim=c(0,0,-0.1),rlim=c(2,5,0),splines="both",model,printPars=F){
  
  sim_vals = Simulation_step(grid,model,ct)
  
  
  # gam modeling
  time_gam = proc.time()
  if(splines=="both"){
    srv.gam = gam(srv~s(mu,lambda, k=10)+s(mu,betaN, k=10)+s(lambda,betaN, k=10),family = binomial,data=dat)
    srv.gam2 = gam(srv~s(p1)+s(p3)+s(p2),family = binomial,data=dat)
  }
  if(splines=="bivariate"){
    srv.gam = gam(srv~s(p1,p2)+s(p1,p3)+s(p2,p3),family = binomial,data=dat)
    srv.gam2=NULL
  }
  if(splines=="univariate"){
    srv.gam = gam(srv~s(p1)+s(p3)+s(p2),family = binomial,data=dat)
    srv.gam2=NULL
  }
  time_gam = emphasis::get.time(time_gam)
  dens = nrow(pars)/prod(rlim-llim)
  
  return(list(gam1=srv.gam,gam2=srv.gam2,density=dens,time_sim=time_sim,time_gam=time_gam))
}

sim_survival_dd <- function (pars, ct, soc=1, max.species=100){
  cbt = 0
  N = soc
  srv = 1
  mu = max(0, pars[1])
  while ((cbt < ct) & (N >= 1) & (N<max.species)) {
    lambda_ct = max(0, pars[2] + pars[3] * N)
    rate_max = (lambda_ct + mu) * N
    u1 = runif(1)
    next_event_time = cbt - log(x = u1)/rate_max
    if (next_event_time < ct) {
      u2 = runif(1)
      to = sample(c(1, 0), size = 1, 
                  prob = c(lambda_ct,mu)/(lambda_ct + mu))
      if (to == 1) {
        N = N + 1
      }
      else {
        N = N - 1
      }
    }
    cbt = next_event_time
  }
  if (N == 0) {
    srv = 0
  }
  return(list(srv=srv))
}


sim_survival <- function (diversification_model, ct, max.species=100,timeLimit=10){
  pars = diversification_model$pars
  model = diversification_model$model
  if(model == "rpd1"){
    st = sim.tree_rpd1(pars,ct,timeLimit = timeLimit)
  }
  if(model == "rpd5"){
    st = sim.tree_rpd5(pars,ct,timeLimit = timeLimit)
  }
  if (length(st$tree) == 1) {
    srv = 0
  }else{
    srv =1
  }
  return(list(tree=st,srv=srv))
}
  

# and we can modify the loglik function to include this feature
loglik.tree.numerical <- function(pars, tree, model,cond=F,condition_polinomal=NULL){
  to = tree$to
  to = head(to,-1)
  to[to!=0] = 1
  spec_times = tree$brts[c(to,0)==1]
  #
  speciations = sapply(spec_times, speciation_rate,tree=tree,pars=pars,model=model,soc=tree$n[1])
  extinctions = rep(pars[1],sum(to==0))
  # 
  inte = intensity.numerical2(tree, pars, model)
  loglik = sum(log(speciations)) + sum(log(extinctions)) - sum(inte)
  
  if(cond){
    loglik = loglik/conditional_probability_cond1_dd(pars,condition_polinomal)
  }
  
  return(loglik)
}




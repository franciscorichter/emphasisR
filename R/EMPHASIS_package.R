### EMPHASIS functions
emphasis <- function(brts,
                     soc=2,
                     model="rpd1",
                     init_par,
                     em_tol=0.25,
                     sample_size_tol=0.005,
                     burnin_sample_size=200,
                     pilot_sample_size=c(200,600),
                     burnin_iterations = 20,
                     parallel=TRUE){

  
  msg1 = paste("Initializing emphasis...")
  msg2 = paste("Age of the tree: ",max(brts))
  msg3 = paste("Number of speciations: ",length(brts))
  msg4 = paste("Diversification model to fit:",model)
  msg5 = "######################################"
  cat(msg1,msg2,msg3,msg4,msg5,sep="\n")
  
  cat( "Performing Phase 1: burn-in",sep= "\n")
  mc = mcEM(brts = brts,
            pars = init_par,
            sample_size = burnin_sample_size,
            model = model,
            soc = soc,
            parallel = parallel,
            print_process = FALSE,
            tol = em_tol,
            burnin = burnin_iterations)
  
  M = mc$mcem
  pars = c(mean(tail(M$par1,n = nrow(M)/2)),
           mean(tail(M$par2,n = nrow(M)/2)),
           mean(tail(M$par3,n = nrow(M)/2)),
           mean(tail(M$par4,n = nrow(M)/2)))
  
  cat("\n",msg5,sep="\n")
  cat( "Phase 2: Assesing required MC sampling size \n")

  for(i in 1:length(pilot_sample_size)){
    cat(paste("\n Sampling size: ",as.character(pilot_sample_size[i]),"\n"))
    mc = mcEM(brts = brts,
              pars = pars,
              sample_size = pilot_sample_size[i],
              model = model,
              soc = soc,
              parallel = parallel,
              print_process = FALSE,
              tol = em_tol,
              burnin = 10)
    ta = tail(mc$mcem,n = nrow(M)/2)
    pars = c(mean(ta$par1),mean(ta$par2),mean(ta$par3),mean(ta$par4))
    M = rbind(M,mc$mcem)
  }
  n.r = get_required_sampling_size(M[-(1:burnin_iterations),],tol = sample_size_tol)
  sample_size = max(pilot_sample_size+2,n.r)
  n.r_old = -1
  j = 1
  while(n.r_old < n.r){
    msg6 = paste0("Required sampling size: ",n.r)
    msg7 = paste0("Phase 3: Performing metaiteration: ",j)
    cat("\n",msg5,msg7,msg6,sep="\n")
    mc = mcEM(brts = brts,
              pars = pars,
              sample_size = sample_size,
              model = model,
              soc = soc,
              parallel = parallel,
              print_process = FALSE,
              tol = em_tol,
              burnin = 2)
    M <- rbind(M,mc$mcem)
    n.r_old = n.r
    j = j+1
    n.r = get_required_sampling_size(M[-(1:burnin_iterations),],tol = sample_size_tol)
    pars = as.numeric(colMeans(mc$mcem)[1:4])
    sample_size = n.r
  }

  cat(pars)
  return(list(pars=pars,MCEM=M))
}

mcEM <- function(brts, 
                 pars, 
                 sample_size, 
                 model, 
                 soc, 
                 tol=0.01,
                 burnin=20,
                 print_process=FALSE,
                 parallel=TRUE,
                 cores=(parallel::detectCores()-1)){
  mcem = NULL
  sde = 10; i=0
  times = NULL
  while(sde > tol){
    i = i+1
    st = mcE_step(brts = brts, 
                  pars = pars, 
                  sample_size = sample_size, 
                  model = model, 
                  no_cores = cores, 
                  parallel = parallel,
                  soc = soc)
    if(max(st$weight[!is.na(st$weight)])==0){
      print(st)
      stop("Only zero likelihood trees, maybe there is underflow")
    }
    M = M_step(st = st, init_par = pars, model = model)
    if(!is.infinite(M$po$value) & !is.na(st$fhat)){ 
      pars = M$po$par
      mcem = rbind(mcem,data.frame(par1=pars[1],
                                   par2=pars[2],
                                   par3=pars[3],
                                   par4=pars[4],
                                   fhat=st$fhat,
                                   sample_size=sample_size))
    }
    if(print_process){
      print(paste("(mean of) loglikelihood estimation: ",mean(mcem$fhat)))
    }
    times = c(times,st$E_time+M$M_time)
    time_p_it = mean(times)
    if(i>burnin){
      mcem_est = mcem[floor(nrow(mcem)/2):nrow(mcem),]
      mcem_est = mcem_est[is.finite(mcem_est$fhat),]
      sde0 = sde
      sde = sd(mcem_est$fhat)/sqrt(nrow(mcem_est))
      mde = mean(mcem_est$fhat)
      msg = paste("Iteration:",i," SE of the loglikelihood: ",sde)
      cat("\r",msg) 
    }else{
      msg = paste("Remining time (burn-in): ",round(time_p_it*(burnin-i),digits = 0),"sec")
      cat("\r",msg) 
    }
  }
  return(list(mcem=mcem))
}

##############################
####### E-step 

mcE_step <- function(brts,pars,sample_size,model,no_cores=2,seed=0,parallel=TRUE,soc=2,printprocess=TRUE){
  if(seed>0) set.seed(seed)
  time = proc.time()
  if(!parallel){
    st =  lapply(1:sample_size,function(i){augment_tree(brts = brts,pars = pars,model=model,soc=soc)} )
  }else{
    st = mclapply(1:sample_size,function(i){augment_tree(brts = brts,pars = pars,model=model,soc=soc)},mc.cores = no_cores)
  }
  trees = lapply(st,function(list) list$tree)
  dim = sapply(st,function(list) nrow(list$tree))
  
  ####
  logf = sapply(trees,loglik.tree(model), pars=pars)
  logg = sapply(trees,sampling_prob, pars=pars,model=model)
  E_time = get.time(time)
  log_weights = logf-logg
  prop_const = max(log_weights)
  log_weights_norm = log_weights - prop_const
  w = exp(log_weights_norm)
  log_fhat = log(mean(w)) + prop_const
  ####
  En = list(weights=w,trees=trees,fhat=log_fhat,logf=logf,logg=logg,dim=dim,E_time=E_time)
  return(En)

}

##############################
####### M-step 

M_step <-function(st,init_par,model,reltol=0.001, gam = NULL){
  
  time0 = proc.time()
  
  sub_st = get_contributing_trees(st)
  
  po = subplex(par = init_par, fn = Q_approx,st = sub_st, gam = gam, loglik = get(paste0("loglik.tree.", model)), hessian = FALSE,control=list(reltol=reltol))

  M_time = get.time(time0)
  return(list(po=po,M_time=M_time))
}


get_contributing_trees <- function(st, exclude_proportion_trees=0){
  w = st$weights/sum(st$weights)
  contributing_trees = (w > exclude_proportion_trees)
  sub_trees = st$trees[contributing_trees]
  effective_sample_size = sum(contributing_trees)
  sub_st = list(trees = sub_trees, weights = st$weights[contributing_trees])
  return(sub_st)
}

Q_approx = function(pars,st,loglik, gam=NULL){
  
  l = sapply(st$trees, loglik, pars=pars, gam = gam)
  w = st$weights
  Q = -sum(l*w)
  return(Q)
}

####
  
  



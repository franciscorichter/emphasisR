### EMPHASIS functions
emphasis <- function(brts,soc=2,model="rpd1",init_par,tol=0.01,parallel=TRUE,name="NN",burnin_sample_size=200,pilot_sample_size=c(200,400,600)){

  input = list(brts=brts,pars = init_par,sample_size=burnin_sample_size,model=model,cores=detectCores()-2,parallel=parallel,soc=soc)
  
  msg1 = paste("Initializing emphasis...")
  msg2 = paste("Name of the clade: ",name)
  msg2 = paste("Age of the tree: ",max(input$brts))
  msg3 = paste("Number of speciations: ",length(input$brts))
  msg4 = paste("Diversification model to fit:",input$model)
  msg5 = "######################################"
  cat(msg1,msg2,msg3,msg4,msg5,sep="\n")
  
  cat( "Performing Phase 1: burn-in",sep= "\n")
  mc = mcEM(input,print_process = F,tol = 0.1,burnin = 20)
  MCEM=mc$mcem
  input$pars = c(mean(tail(mc$mcem$par1,n = 10)),mean(tail(mc$mcem$par2,n = 10)),mean(tail(mc$mcem$par3,n = 10)),mean(tail(mc$mcem$par4,n = 10)))
  
  cat("\n",msg5,sep="\n")
  cat( "Phase 2: Assesing required MC sampling size")
  MC = list()
  
  input$sample_size = pilot_sample_size1
  M = NULL
  for(i in 1:length(pilot_sample_size)){
    input$sample_size = pilot_sample_size[i]
    cat(paste("\n Sampling size: ",as.character(input$sample_size),"\n"))
    MC[[i]] = mc = mcEM(input,print_process = FALSE,burnin = 5,tol = tol)
    ta = tail(mc$mcem,n = floor(nrow(mc$mcem)/2))
    input$pars = c(mean(ta$par1),mean(ta$par2),mean(ta$par3),mean(ta$par4))
    MCEM = rbind(MCEM,mc$mcem)
    M = rbind(M,mc$mcem)
  }
  
  n.r = get_required_sampling_size(M,tol = tol*10)
  if(n.r<0) n.r = input$sample_size*2
  input$sample_size = n.r
  msg6 = paste0("Required sampling size: ",n.r)
  msg7 = "Phase 3: First estimation"
  cat("\n",msg5,msg7,msg6,sep="\n")
  mc = mcEM(input,print_process = FALSE,burnin = 10,tol = tol)
  M<-rbind(M,mc$mcem)
  n.r = get_required_sampling_size(M,tol=tol*10)
  if(n.r>input$sample_size){
    input$sample_size = n.r
    msg6 = paste0("Required sampling size: ",n.r)
    msg7 = "Last phase: Second estimation"
    cat("\n",msg5,msg7,msg6,sep="\n")
    mc = mcEM(input,print_process = FALSE,burnin = 5,tol = tol)
    M<-rbind(M,mc$mcem)
  }
  msg6 = paste0("Required sampling size: ",n.r)
  msg7 = "Last phase: Second estimation"
  cat("\n",msg5,msg7,msg6,sep="\n")
  pars = as.numeric(colMeans(mc$mcem)[1:4])
  cat(pars)
  sp=sample_size_determination(f = M$fhat,n = M$sample_size,tol = tol*10)
  return(list(pars=pars,mc=mc,MCEM=M,required_sample_size=n.r,diag1=sp,clade=name,sample_size_completition=(sp$n.r<input$sample_size)))
}

mcEM <- function(input,print_process=FALSE,tol=0.01,burnin=20){
  pars = input$pars
  mcem = NULL
  sample_size = input$sample_size
  sde = 10; i=0
  times = diffsd = NULL
  while(sde > tol){
    i = i+1
    st = mcE_step(brts = input$brts, pars = pars,sample_size=sample_size,model=input$model,no_cores=input$cores,parallel=input$parallel,soc=input$soc)
    if(max(st$weight)==0){
      print(st)
      stop("Only zero likelihood trees, maybe there is underflow")
    }
    M = M_step(st = st, init_par = pars, model = input$model)
    if(!is.infinite(M$po$value) & !is.na(log(st$fhat))){ 
      pars = M$po$par
      mcem = rbind(mcem,data.frame(par1=pars[1],par2=pars[2],par3=pars[3],par4=pars[4],fhat=st$logf,E_time=st$E_time,M_time=M$M_time,sample_size=sample_size))
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
      diffsd = c(diffsd,min(0,sde-sde0))
      param = mean(mcem_est$par1)
      mde = mean(mcem_est$fhat)
      msg = paste("Iteration:",i," SE of the loglikelihood: ",sde)
      cat("\r",msg) 
    }else{
      msg = paste("Remining time (burn-in): ",round(time_p_it*(burnin-i),digits = 0),"sec")
      cat("\r",msg) 
    }
  }
  return(list(mcem=mcem,st=st,M=M))
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
  log_weights_norm = log_weights - max(log_weights)
  w = exp(log_weights_norm)
  ####
  En = list(weights=w,trees=trees,fhat=mean(exp(log_weights)),logf=logf,logg=logg,dim=dim,E_time=E_time)
  return(En)

}

##############################
####### M-step 

M_step <-function(st,init_par,model,reltol=0.001){
  
  time0 = proc.time()
  
  sub_st = get_contributing_trees(st)
  
  po = subplex(par = init_par, fn = Q_approx,st = sub_st, loglik = get(paste0("loglik.tree.", model)), hessian = FALSE,control=list(reltol=reltol))

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

Q_approx = function(pars,st,loglik){
  
  l = sapply(st$trees, loglik, pars=pars)
  w = st$weights
  Q = -sum(l*w)
  return(Q)
}

####
  
  



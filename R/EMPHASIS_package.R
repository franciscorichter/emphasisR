### EMPHASIS functions
emphasis <- function(brts,soc=2,model="rpd1",n_up_ss=3,init_par,sample_size=1000,name="temp",parallel=TRUE){
  init_sampling_size = sample_size
  input = list(brts=brts,pars = init_par,sample_size=init_sampling_size,model=model,cores=detectCores(),parallel=parallel,n_it = 1000,soc=soc)
  
  print(paste("Optimizing the likelihood - this may take a while."))
  print(paste("Age of the tree: ",max(input$brts)))
  print(paste("Number of speciations: ",length(input$brts)))
  print(paste("Diversification model to fit:",input$model))
  print(paste("initial parameters. ","mu: ",input$pars[1],"lambda: ",input$pars[2],"beta: ",input$pars[3]))  
  
  MCEM=NULL
  for(i in 1:n_up_ss){
    print(paste("Sampling size: ",input$sample_size))
    mc = mcEM(input,file=paste(name,"_",input$model,"_",as.character(input$sample_size),".RData",sep=""),n_it=input$n_it,print_process = FALSE)
    load(file=paste(name,"_",input$model,"_",as.character(input$sample_size),".RData",sep=""))
    input$sample_size = input$sample_size*2
    ta = tail(mc$mcem,n = floor(nrow(mc$mcem)/2))
    input$pars = c(mean(ta$par1),mean(ta$par2),mean(ta$par3))
    MCEM = rbind(MCEM,mc$mcem)
  }
  return(list(mc=mc,MCEM=MCEM))
  
  
}

mcEM <- function(input,file=".RData",print_process=FALSE,n_it=NULL,tol=0.01,burnin=20){
  pars = input$pars
  mcem = NULL
  sample_size = input$sample_size
  key = 0
  sde=mde=10; i=0
  times=NULL
  while(sde > tol){
    i = i+1
    if(print_process){
      print(paste("Performing E step, iteration",i))
      print(pars)
    }
    st = mcE_step(brts = input$brts, pars = pars,sample_size=sample_size,model=input$model,no_cores=input$cores,parallel=input$parallel,soc=input$soc)
    if(print_process){
      print(paste("Performing M step, iteration",i))
      print(paste("loglikelihood random estimation: ",log(st$fhat)))
    }
    M = M_step(st = st, init_par = pars, model = input$model)
    if(!is.infinite(M$po$value) & !is.na(log(st$fhat))){ 
      pars = M$po$par
      mcem = rbind(mcem,data.frame(par1=pars[1],par2=pars[2],par3=pars[3],par4=pars[4],fhat=log(st$fhat),E_time=st$E_time,M_time=M$M_time,sample_size=sample_size))
    }
    if(print_process){
      print(paste("Q random estimation: ",log(M$po$value)))
      print(paste("(mean of) loglikelihood estimation: ",mean(mcem$fhat)))
    }
    save(input,mcem,file=file)
    if(i>burnin){
      mcem_est = mcem[floor(nrow(mcem)/2):nrow(mcem),]
      sde = sd(mcem_est$fhat)/nrow(mcem_est)
      mde = mean(mcem_est$fhat)
      msg1 = paste("Iteration:",i,"Time per iteration:",round(st$E_time+M$M_time,digits = 2))
      msg2 = paste("loglikelihood estimation:",mde,"Standard Error:",sde)
      cat("\r",msg1, msg2, sep="\n")
      #cat(msg2)
    }else{
      times = c(times,st$E_time+M$M_time)
      time_p_it = mean(times)
      msg = paste("Remining time for burn-in: ",round(time_p_it*(burnin-i),digits = 0),"sec")
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
    st =  lapply(1:sample_size,function(i){augment_tree_tj(brts = brts,pars = pars,model=model,soc=soc)} )
  }else{
    st = mclapply(1:sample_size,function(i){augment_tree_tj(brts = brts,pars = pars,model=model,soc=soc)},mc.cores = no_cores)
  }
  trees = lapply(st,function(list) list$tree)
  dim = sapply(st,function(list) nrow(list$tree))
  
  ####
  logf = sapply(trees,loglik.tree(model), pars=pars)
  logg = sapply(trees,sampling_prob, pars=pars,model=model)
  E_time = get.time(time)
  log_weights = logf-logg
  w = exp(log_weights)
  ####
  En = list(weights=w,trees=trees,fhat=mean(w),logf=logf,logg=logg,dim=dim,E_time=E_time)
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
  
  



### EMPHASIS functions


emphasis <- function(input,file=".RData",print_process=TRUE,n_it=NULL,tol=0.01){
  pars = input$pars
  mcem = NULL
  sample_size = input$sample_size
  key = 0
  for(i in 1:n_it){
    if(print_process){
      print(paste("Performing E step, iteration",i))
      print(pars)
    }
    st = mcE_step(brts = input$brts, pars = pars,sample_size=sample_size,model=input$model,no_cores=input$cores,parallel=input$parallel,soc=input$soc)
    if(print_process){
      print(paste("Performing M step, iteration",i))
      print(paste("loglikelihood random estimation: ",log(st$fhat)))
      print(paste("(mean of) loglikelihood estimation: ",mean(mcem$fhat)))
    }
    M = M_step(st = st, init_par = pars, model = input$model)
    if(!is.infinite(M$po$value) & !is.na(log(st$fhat))){ 
      pars = M$po$par
      mcem = rbind(mcem,data.frame(par1=pars[1],par2=pars[2],par3=pars[3],par4=pars[4],fhat=log(st$fhat),E_time=st$E_time,M_time=M$M_time,sample_size=sample_size))
    }
    if(print_process){
      print(paste("Q random estimation: ",M$po$value))
      print(paste("(mean of) loglikelihood estimation: ",mean(mcem$fhat)))
    }
    save(input,mcem,file=file)
    if(i>10){
      if( abs( mean(mcem$fhat)-mean(mcem$fhat[-nrow(mcem)]) ) < tol){
        key = 1 + key
        if(key==3){
          break
        }
      }else{
        key=0
      }
    }
  }
}

##############################
####### E-step 

mcE_step <- function(brts,pars,sample_size,model,no_cores=2,seed=0,parallel=TRUE,soc=2){
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
  logg = sapply(trees,sampling_prob, pars=pars,model=model,soc=soc)
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


emphasis_production <- function(input,file=".RData",print_process=TRUE,mcem=NULL,n_it=NULL){
  if(is.null(mcem)){
    pars = input$pars
  }else{
    pars = mcem[nrow(mcem),1:4]
    pars = pars[!is.na(pars)]
  }
  print("Initialization of emphasis")
  for(i in 1:n_it){
    if(print_process){
      print(paste("iteration",i))
      print(pars)
    }
    ### sample size
    if(length(input$sample_size)>1){
      if(i < 400){
        sample_size = input$sample_size[1]
      }
      if(i >= 400 & i<800){
        sample_size = input$sample_size[2]
      }
      if(i >=800){
        sample_size = input$sample_size[3]
      }
    }else{
      sample_size = input$sample_size
    }
    ####
    st = mcE_step(brts = input$brts, pars = pars,sample_size=sample_size,model=input$model,no_cores=input$cores,parallel=input$parallel,soc=input$soc)
    if(print_process){
      print(paste("loglikelihood estimation: ",log(st$fhat)))
    }
    M = M_step(st = st, init_par = pars, model = input$model)
    if(!is.infinite(M$po$value)) pars = M$po$par
    mcem = rbind(mcem,data.frame(par1=pars[1],par2=pars[2],par3=pars[3],par4=pars[4],fhat=log(st$fhat),E_time=st$E_time,M_time=M$M_time,sample_size=sample_size))
    save(input,mcem,file=file)
    remaining_time = calculate_rem_time(mcem)
  }
}

remaining_time <- function(mcem,n_it=1000){
  
}
  
  
  



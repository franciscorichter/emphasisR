### EMPHASIS functions

mcem.tree <- function(input,max_iterations=100000,report=TRUE,file=NULL){
  pars = PARS = input$pars
  MCEM=data.frame(loglik_hat=NULL,E_time=NULL,M_time=NULL,sample_size=NULL)
  for(i in 1:max_iterations){
    
    st = mc_sample_independent_trees(brts = input$brts,pars = pars,nsim = input$sample_size,model = input$model, importance_sampler = input$importance_sampler,no_cores = input$cores, method = input$method)
lg = log(st$fhat)
    if(report){ 
      print(paste(c("parameters: ",pars)))
      print(paste("loglikelihood: ",lg))
      print(paste("sample size: ",input$sample_size))
    }
    M = M_step(st = st,init_par = pars,model = input$model)
    if(!is.na(M$po$value)){
      pars = M$po$par
    }
    print(M$po$value)
    PARS = rbind(pars,PARS)
    MCEM = rbind(MCEM,data.frame(loglik_hat=log(st$fhat),E_time=st$E_time,M_time=M$M_time,sample_size=input$sample_size))
    save(MCEM,input,PARS,file=file)
  }
}


##############################
####### E-step 


mc_sample_independent_trees <- function(brts,pars,nsim=1000,model="dd",importance_sampler="emphasis",no_cores=2,pars3=NULL,maxnumspec=NULL,seed=0,method="inverse",parallel=TRUE){
  
  time=proc.time()
  if(seed>0) set.seed(seed)
  if(method == "thinning"){
    E = mc_augmentation_thinning(brts = brts,pars = pars,model = model,importance_sampler = importance_sampler,sample_size = nsim,parallel = parallel,no_cores = no_cores)
  }else{
    E = mc_augmentation_inverse(brts = brts,pars = pars,model = model,importance_sampler = importance_sampler,sample_size = nsim,parallel = parallel,no_cores = no_cores)
  }
  return(E)
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



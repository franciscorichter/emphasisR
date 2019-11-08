pars = c(1.94,0.33,0.06)
## input 
pars = c(1,0.001,0.1)
input = list(brts=brts_vangidae,pars=pars,sample_size=100,model="dd",importance_sampler="emphasis",cores=2,maxnumspec=40,method="inverse")


##  way 1 


for(i in 1:100){
  print(pars)
  time = proc.time()
  st1 = mc_sample_independent_trees(brts = input$brts,pars = pars,nsim = input$sample_size,model = input$model, importance_sampler = input$importance_sampler,no_cores = input$cores,maxnumspec = input$maxspec, method = input$method)
  get.time(time)
  print(log(st1$fhat))

  time = proc.time()
  M = M_step(st = st1,init_par = pars,model = input$model)
  print(M$po$value)
  get.time(time)
  pars = M$po$par
}




M_step <-function(st,init_par = NULL,model="dd",proportion_of_subset=1){
  time0 = proc.time()
  weights = st$weights/sum(st$weights)
  effective_sample_size = sum(weights>0)
  po = subplex(par = init_par, fn = Q_approx,st = st,model=model,hessian = FALSE)
  M_time = get.time(time0)
  return(list(po=po,M_time=M_time,effective_sample_size=effective_sample_size))
}

Q_approx = function(pars,st,model="dd",initspec=1){
  get_llik <- function(tree) nllik.tree(pars=pars,tree=tree,initspec = initspec,model=model)
  l = sapply(st$trees, get_llik)
  w = st$weights/(sum(st$weights))
  impossible_trees = which(w==0)
  if(length(impossible_trees)>0){
    l = l[-impossible_trees]
    w = w[-impossible_trees]
  }
  Q = sum(l*w)
  return(Q)
}





###  way 3

time = proc.time()
st2 = mc_augmentation_thinning(brts=input$brts,pars = pars,model = input$model,importance_sampler = input$importance_sampler,sample_size = input$sample_size,parallel = TRUE,no_cores = input$cores)
get.time(time)

time = proc.time()
M = M_step(st = st,init_par = pars,model = input$model)
get.time(time)


########################


next_speciation_time = rnhpp(time0 = cbt,time_max = next_bt,tree = tree,model = model,pars = pars)

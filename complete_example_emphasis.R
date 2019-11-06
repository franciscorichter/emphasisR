
## input 
pars = c(1,0.01,0.1)
input = list(brts=brts_vangidae,pars=pars,sample_size=1000,model="dd",importance_sampler="emphasis",cores=2,maxnumspec=40,method="inverse")


##  way 1 

time = proc.time()
st1 = mc_sample_independent_trees(brts = input$brts,pars = pars,nsim = input$sample_size,model = input$model, importance_sampler = input$importance_sampler,no_cores = input$cores,maxnumspec = input$maxspec, method = input$method)
get.time(time)

time = proc.time()
M = M_step(st = st1,init_par = pars,model = input$model)
get.time(time)


###  way 2 
  
time = proc.time()
st = mc_augmentation(brts = input$brts,pars = pars,model = input$model,importance_sampler = input$importance_sampler,sample_size = input$sample_size)
get.time(time)

time = proc.time()
M = M_step(st = st,init_par = pars,model = input$model)
get.time(time)


###  way 3

time = proc.time()
st2 = mc_augmentation_thinning(brts=input$brts,pars = pars,model = input$model,importance_sampler = input$importance_sampler,sample_size = input$sample_size,parallel = TRUE,no_cores = input$cores)
get.time(time)

time = proc.time()
M = M_step(st = st,init_par = pars,model = input$model)
get.time(time)


########################


next_speciation_time = rnhpp(time0 = cbt,time_max = next_bt,tree = tree,model = model,pars = pars)

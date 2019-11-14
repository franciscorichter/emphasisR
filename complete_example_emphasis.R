devtools::install_github("franciscorichter/emphasis")
library(emphasis)
## input 1. Simple tree, Ok initial parameters
pars = c(0.2,0.001,0.02)
input = list(brts=brts_vangidae,pars=pars,sample_size=100,model="dd",importance_sampler="emphasis",cores=2,maxnumspec=40,method="inverse")
## input 2. Big tree. 
pars = c(0.05,0.00001,0.01)
input = list(brts=brts_cetacea,pars=pars,sample_size=10,model="dd",importance_sampler="emphasis",cores=2,maxnumspec=40,method="inverse")


##  complete iterative framework 

for(i in 1:100){
  print(pars)
  time = proc.time()
  st1 = mc_sample_independent_trees(brts = input$brts,pars = pars,nsim = input$sample_size,model = input$model, importance_sampler = input$importance_sampler,no_cores = input$cores,maxnumspec = input$maxnumspec, method = input$method)
  get.time(time)
  time = proc.time()
  M = M_step(st = st1,init_par = pars,model = input$model)
  get.time(time)
  pars = M$po$par
}



#### now with rpd 




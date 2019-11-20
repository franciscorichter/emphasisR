devtools::install_github("franciscorichter/emphasis")
library(emphasis)
## input 1. Simple tree, Ok initial parameters
l0 = 3.659926
ga = 0.1479505
mu =  0.182003
al = tan(1)
be = 0
pars = c(l0,ga,mu)*0.1
input = list(brts=brts_dendroica,pars=pars,sample_size=1000,model="dd",importance_sampler="emphasis",cores=2,method="inverse")


##  complete iterative framework 
MCEM = data.frame(par1=NULL,par2=NULL,par3=NULL,loglik_hat=NULL,E_time=NULL,M_time=NULL)
for(i in 1:10000){
  print(pars)
  st = mc_sample_independent_trees(brts = input$brts,pars = pars,nsim = input$sample_size,model = input$model, importance_sampler = input$importance_sampler,no_cores = input$cores, method = input$method)
  print(paste("loglikelihood: ",log(st$fhat)))
  time = proc.time()
  M = M_step(st = st,init_par = pars,model = input$model)
  M_time = get.time(time)
  MCEM = rbind(MCEM,data.frame(par1=pars[1],par2=pars[2],par3=pars[3],loglik_hat=log(st$fhat),E_time=st$E_time,M_time=M$M_time))
  pars = M$po$par
  save(MCEM,input,file="MCEM_dendroica_v2_20Nov.RData")
}


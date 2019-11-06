

#M_step <-function(st,init_par = NULL,model="dd",proportion_of_subset=1){
#  time0 = proc.time()
#  weights = st$weights/sum(st$weights)
#  weights_sorted = sort(weights,decreasing = TRUE)
#  a = which(cumsum(weights_sorted)>=proportion_of_subset)[1]
#  max.weight = weights_sorted[a]
#  sub_st = lapply(st, "[", weights>=max.weight)
#  loglik_proportion = sum(weights[weights>=max.weight])/sum(weights)
#  effective_sample_size = length(weights[weights>=max.weight])
#  po = subplex(par = init_par, fn = Q_approx,st = sub_st,model=model,hessian = TRUE)
#  M_time = get.time(time0)
#  return(list(po=po,loglik_proportion=loglik_proportion,effective_sample_size=effective_sample_size,M_time=M_time))
#}


M_step <-function(st,init_par = NULL,model="dd",proportion_of_subset=1){
  time0 = proc.time()
  weights = st$weights/sum(st$weights)
  effective_sample_size = length(weights>0)
  po = subplex(par = init_par, fn = Q_approx,st = st,model=model,hessian = FALSE)
  M_time = get.time(time0)
  return(list(po=po,M_time=M_time,effective_sample_size=effective_sample_size))
}




###########################

###  simulation of extinct species



############################
# emphasis data augmentation importance sampler
##############################
bt2tree <- function(brts){
  if(length(brts)>1){
    list(brts=brts,df=data.frame(parent=rep("s1",length(brts)-1),child=paste("s",2:length(brts),sep="")))
  }
  if(length(brts==1)){
    list(brts=brts,df=NULL)
  }
}




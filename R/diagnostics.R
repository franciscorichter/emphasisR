 ##  Ltt expected plots 

expectedLTT2 <- function(pars, ct=15, n_it= 10, sn_it = 100, median=FALSE, dim=NULL){ # it should also include model
  # add message that checks that do not have extinct species
  BT = vector("list", n_it)
  # diversity dependence model. parameters are par=c(lambda,mu,K)

  for (i in 1:(n_it)){
    s = sim_brts(pars = pars,model = model,ct = ct)
    BT[[i]] = s$brts
  }

  n.obs <- sapply(a, length)
  #m.obs = round(mean(n.obs))
  seq.max <- seq_len(max(n.obs))
 
  mat <- t(sapply(BT, "[", i = seq.max))
  #mat = mat[,1:m.obs]
  expect = colMeans(mat,na.rm = TRUE)
  if(median){expect = colMedians(mat,na.rm = TRUE)}
  ce = cumsum(expect)
  bt=c(cumsum(expect[ce<15]),15)
  expect = list(bt=bt, Ex = 1:length(cumsum(expect[ce<15])),same_dim=same_dim)
  return(expect)
}

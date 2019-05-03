
# Method 1
weight.unif <- function(tree,pars,d,ct){
  lik.tree(pars,tree,initspec = 1,topology = F)*((catalan(d)*(ct^(2*d)))/(factorial(2*d)))
}


sim.unif.mtree <- function(top,obt){
  nmbt = length(top) #length of missing bt
  mbt = sort(runif(nmbt,min=0,max=max(obt)))
  df = data.frame(bt=c(mbt,obt),to=c(top,rep(2,length(obt))))
  df = df[order(df$bt),]
  return(df)
}

sim.unif.tree.df <- function(top,obt){
  df = df2tree2(sim.unif.mtree(top,obt))
  return(df)
}

# Method 2 

loglik_reconstructed_tree_approx_S <- function(brts,brts_exts_list,la,mu,age)
{
  endmc <- length(brts_exts_list)
  brts <- -abs(brts)
  S <- dim(brts_exts_list[[1]])[1]
  if(is.null(S))
  {
    S <- 0
  }
  loglik_full_tree_vec <- rep(NA,endmc)
  for(mc in 1:endmc)
  {
    loglik_full_tree_vec[mc] <- loglik_full_tree(brts = brts,brts_exts = brts_exts_list[[mc]],la = la,mu = mu,age = age)
  }
  loglik_full_tree_vec_max <- max(loglik_full_tree_vec)
  if(loglik_full_tree_vec_max != -Inf)
  {
    #logfactor <-  2 * S * log(age) - S * log(2) # Works for S = 1
    logfactor <- log_catalan(S) + 2 * S * log(age) - lgamma(2 * S + 1) # Works for S = 1
    loglik <- logfactor + loglik_full_tree_vec_max + log(mean(exp(loglik_full_tree_vec - loglik_full_tree_vec_max)))
  } else
  {
    loglik <- -Inf
  }
  return(loglik)
}

#comparison of two methods

S = d = 2
brts = obt = 1
la = 0.2
mu = 0.1
pars = c(la,mu,Inf)
endmc = nsim = 10000
age = ct = 1
#S1. get top
topologies = get.topologies(S)
#S2. sample trees
##M1 
topi = get.topologies(S)
topi = replicate(nsim,topi,simplify = F)
ntree = lapply(topi,sim.unif.tree.df,obt=obt) 
##M2
for(mc in 1:endmc)
{
  brts_exts_list[[mc]] <- sample_brts_exts(S,la,mu,age,topologies)
}
ntree2 = lapply(brts_exts_list,brts2tree,brts=brts)
#S3. sum lliks
##M1
liks = lapply(ntree,weight.unif,pars=pars,d=S,ct=max(obt))
log(Reduce("+", liks))

liks2 = lapply(ntree2,weight.unif,pars=pars,d=S,ct=max(obt))
log(Reduce("+", liks2))

liks2 = sapply(ntree2,weight.unif,pars=pars,d=S,ct=max(obt))
log(mean(liks))
##M2
loglik_reconstructed_tree_approx_S(brts = brts,brts_exts_list = brts_exts_list,la = la,mu = mu,age = age)

############################
#S3. 


##M1
lliks = -sapply(ntree2,nllik.tree,pars=pars,topology=F,initspec=1)
est1 = log(mean(exp(lliks)*((catalan(d)*(ct^(2*d)))/(factorial(2*d)))))

liks = sapply(ntree2,lik.tree,pars=pars,topology=F,initspec=1)
est2 = log(mean(liks)*((catalan(d)*(ct^(2*d)))/(factorial(2*d))))

##M2
for(mc in 1:endmc)
{
  loglik_full_tree_vec[mc] <- loglik_full_tree(brts = brts,brts_exts = brts_exts_list[[mc]],la = la,mu = mu,age = age)
}
all.equal(lliks,loglik_full_tree_vec)

logfactor <- log_catalan(S) + 2 * S * log(age) - lgamma(2 * S + 1) # Works for S = 1
#loglik <- logfactor + loglik_full_tree_vec_max + log(mean(exp(loglik_full_tree_vec - loglik_full_tree_vec_max)))
loglik <- logfactor  + log(mean(exp(loglik_full_tree_vec)))
all.equal(loglik,est1)






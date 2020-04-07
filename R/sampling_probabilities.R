# augmentaion (sampling) probability 

sampling_prob <- function(tree,pars,model,numerical=FALSE){

  if(numerical){
    inte = intensity.numerical(tree = tree,pars = pars,model=model)
  }else{
    intensity.temp = get(paste0("intensity.", model))
    inte = intensity.temp(tree=tree,pars=pars)
  }
  
  missing_speciations = (tree$to == 1)
  if(sum(missing_speciations)>0){
    top = head(tree$to,-1)
    N = tree$n
    soc = N[1]
    
    nb = N[missing_speciations]
    No = c(soc,soc+cumsum(top==2))[missing_speciations]
    Ne = c(0,cumsum(top==1)-cumsum(top==0))[missing_speciations]
    
    brts_miss = tree$brts[missing_speciations]
    lambda_b = sapply(brts_miss,speciation_rate,tree = tree,pars = pars,model = model,soc=soc)
    text = tree$t_ext[missing_speciations]-brts_miss
    mu = max(0,pars[1])
    logg = -sum(inte)+sum(log(nb)+log(mu)+log(lambda_b)-mu*text-log(2*No+Ne))
  }else{
    logg = -sum(inte)
  }
  return(logg)
}

intensity.rpd1 <- function(tree, pars){
  
  mu = max(0,pars[1])
  n = tree$n
  lambda = pmax(0,pars[2] + pars[3] * n)
  wt = diff(c(0,tree$brts))
  brts_i = tree$brts
  brts_im1 = c(0,tree$brts[-nrow(tree)])
  sigma_over_tree = n*(lambda)*(wt - (exp(-mu*max(tree$brts))/mu)*(exp(mu*brts_i)-exp(mu*brts_im1)) )
  
  return(sigma_over_tree)
}

intensity.numerical <- function(tree, pars, model){
  nh_rate <- function(x){
    speciation_rate(tm=x,tree = tree,pars = pars,model = model,soc=tree$n[1],sum_lambda = TRUE)*(1-exp(-(max(tree$brts)-x)*pars[1]))
  }
  brts_i = tree$brts
  brts_im1 = c(0,brts_i[-length(brts_i)])
  inte = vector(mode="numeric",length = length(brts_i))
  for(i in 1:length(brts_i)){
    inte[i] = pracma:::quad(f = Vectorize(nh_rate),xa = brts_im1[i]+0.00000000001,xb = brts_i[i])
  }
  return(inte)
}

intensity.rpd5c <- function(tree,pars){
  n = tree$n; Pt = c(0,tree$pd)
  c2 = pars[4]*((n-1)/n); c3 = exp(-pars[1]*max(tree$brts)); c4 = pars[1]
  c1 = pars[2] + pars[3]*n + pars[4]*((Pt[-nrow(tree)]-n*c(0,tree$brts[-length(tree$brts)]))/n)
  
  roots = -c1/c2

  brts_i = tree$brts
  brts_im1 = c(0,brts_i[-length(brts_i)])
  inte = vector(mode="numeric",length = length(brts_i))
  for(i in 1:nrow(tree)){
    if((roots[i] > brts_im1[i]) & (roots[i] < brts_i[i])  ){
      if((pars[2] + pars[3] * n[i] + pars[4] *(Pt[i]/n[i]))<0){
        inte[i] = (ind.rpd5(x=brts_i[i],c1[i],c2[i],c3,c4)-ind.rpd5(x=roots[i],c1[i],c2[i],c3,c4))*n[i]
      }else{
        inte[i] = (ind.rpd5(x=roots[i],c1[i],c2[i],c3,c4)-ind.rpd5(x=brts_im1[i],c1[i],c2[i],c3,c4))*n[i]
      }
    }else{
      inte[i] = (ind.rpd5(x=brts_i[i],c1[i],c2[i],c3,c4)-ind.rpd5(x=brts_im1[i],c1[i],c2[i],c3,c4))*n[i]
    }
  }
  return(inte)
}

ind.rpd5 <- function(x,c1,c2,c3,c4){
  r = (c2 * x^2)/2 + c1*x - (c3*exp(c4*x)*(c2*(c4*x-1)+c1*c4))/(c4^2) 
  return(r)
}

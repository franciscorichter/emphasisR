loglik.tree <- function(model){
  log.lik = get(paste0("loglik.tree.", model))
  return(log.lik)
}

# likelihood functions 

loglik.tree.rpd1 <- function(pars,tree){
  to = tree$to
  to = head(to,-1)
  to[to==2] = 1
  
  mu = max(0,pars[1])
  wt = diff(c(0,tree$brts))
  
  n = tree$n
  
  lambda = pmax(0, pars[2] + pars[3] * n)
  rho = pmax(lambda[-nrow(tree)] * to + mu * (1-to),0)
  
  sigma_over_tree = n*(mu+lambda)*wt 
  
  log.lik = -sum(sigma_over_tree) + sum(log(rho))
  return(log.lik)
}

loglik.tree.rpd5c <- function(pars,tree){
  if(pars[4]==0){
    log.lik = loglik.tree.rpd1(pars,tree)
  }else{
    to = tree$to
    to = head(to,-1)
    to[to==2] = 1
    
    mu = max(0,pars[1])
    wt = diff(c(0,tree$brts))
    
    n = tree$n
    Pt = c(0,tree$pd[-nrow(tree)])
    pd2 = Pt + n*wt
    
    brts_i = tree$brts
    brts_im1 = c(0,brts_i[-length(brts_i)])
    
    lambda = pmax(0,pars[2] + pars[3]*n + pars[4] * (pd2-brts_i)/n )
    rho = pmax(lambda[-length(n)] * to + mu * (1-to),0)
    c1 = pars[2]+pars[3]*n+(pars[4]/n)*(Pt-brts_im1*n)
    c2 = pars[4]*((n-1)/n)
    inte = NULL
    for(i in 1:nrow(tree)){
      if( (brts_im1[i] > (-c1[i]/c2[i])) & (brts_i[i] < (-c1[i]/c2[i])) ){
        if(c2[i]>0){
          brts_im1[i] = -c1[i]/c2[i]
        }else{
          brts_i[i] = -c1[i]/c2[i]
        }
      }
      inte[i] = n[i]*(mu*wt[i] + c1[i]*(brts_i[i]-brts_im1[i]) + c2[i]*(brts_i[i]^2-brts_im1[i]^2)/2)
    }
    log.lik = -sum(inte) + sum(log(rho))
  }
  return(log.lik)
  
}

############################################################

loglik.tree.numerical <- function(pars, tree, model,cond=F,condition_polinomal=NULL){
  to = tree$to
  to = head(to,-1)
  to[to!=0] = 1
  spec_times = tree$brts[c(to,0)==1]
  #
  speciations = sapply(spec_times, speciation_rate,tree=tree,pars=pars,model=model,soc=tree$n[1])
  extinctions = rep(pars[1],sum(to==0))
  # 
  inte = intensity.numerical2(tree, pars, model)
  loglik = sum(log(speciations)) + sum(log(extinctions)) - sum(inte)
  
  if(cond){
    # not tested yet
    log.lik = log.lik/condition_polinomal(pars,diversification_model)
  }
  
  return(loglik)
}



intensity.numerical2 <- function(tree, pars, model){
  nh_rate <- function(x){
    speciation_rate(tm=x,tree = tree,pars = pars,model = model,soc=tree$n[1],sum_lambda = TRUE)+extinction_rate(tm=x,tree = tree,pars = pars,model = model,soc=tree$n[1],sum_rate = TRUE)
  }
  brts_i = tree$brts
  brts_im1 = c(0,brts_i[-length(brts_i)])
  inte = vector(mode="numeric",length = length(brts_i))
  for(i in 1:length(brts_i)){
    inte[i] = pracma:::quad(f = Vectorize(nh_rate),xa = brts_im1[i]+0.00000000001,xb = brts_i[i])
  }
  return(inte)
}

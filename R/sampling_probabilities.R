# augmentaion (sampling) probability 

log_sampling_prob_nh <- function(df,pars,model="dd",soc,...){
  to = top = head(df$to,-1)
  to[to==2] = 1
  initspec=soc
  N = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  brts_i = df$brts
  brts_im1 = c(0,df$brts[-nrow(df)])
  missing_speciations = (df$to == 1)
  nb = N[missing_speciations]
  No = c(1,1+cumsum(top==2))[missing_speciations]
  Ne = c(0,cumsum(top==1)-cumsum(top==0))[missing_speciations]
  lambda_b = sapply(df$brts[df$to==1]-0.000000001,speciation_rate,tree = df,pars = pars,model = model,soc)
  if(length(lambda_b)==0) lambda_b = 1
  text = df$t_ext[df$to==1]-df$brts[df$to==1]
  mu = max(0,pars[1])
  inte = vector(mode = "numeric",length = length(brts_i))
  for(i in 1:length(brts_i)){
    inte[i] = intensity(x=brts_i[i],tree = df,model = model,time0 = brts_im1[i],pars = pars,soc)
  }
  logg = -sum(inte)+sum(log(nb)+log(mu))-sum(mu*text)+sum(log(lambda_b))-sum(log(2*No+Ne))
  return(logg)
}

log_sampling_prob_nh_old <- function(df,pars,model="dd",soc,...){
  if(is.null(df$t_ext)) df$t_ext  = df$bte
  b = max(df$brts)
  to = top = head(df$to,-1)
  to[to==2] = 1
  initspec=soc
  N = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  brts_i = df$brts
  brts_im1 = c(0,df$brts[-nrow(df)])
  missing_speciations = (df$to == 1)
  nb = N[missing_speciations]
  No = c(1,1+cumsum(top==2))[missing_speciations]
  Ne = c(0,cumsum(top==1)-cumsum(top==0))[missing_speciations]
  lambda_b = sapply(df$brts[df$to==1]-0.000000001,speciation_rate,tree = df,pars = pars,model = model,soc)
  if(length(lambda_b)==0) lambda_b = 1
  text = df$t_ext[df$to==1]-df$brts[df$to==1]
  mu = max(0,pars[1])
  inte = vector(mode = "numeric",length = length(brts_i))
  for(i in 1:length(brts_i)){
    inte[i] = emphasis:::intensity(x=brts_i[i],tree = df,model = model,time0 = brts_im1[i],pars = pars,soc)
  }
  logg = -sum(inte)+sum(log(nb)+log(mu))-sum(mu*text)+sum(log(lambda_b))-sum(log(2*No+Ne))
  return(logg)
}

sampling_prob <- function(df,pars,model,soc){
  to = top = head(df$to,-1)
  to[to==2] = 1
  N = c(soc,soc+cumsum(to)+cumsum(to-1))
  
  brts_i = df$brts
  brts_im1 = c(0,df$brts[-nrow(df)])
  
  missing_speciations = (df$to == 1)
  nb = N[missing_speciations]
  No = c(1,1+cumsum(top==2))[missing_speciations]
  Ne = c(0,cumsum(top==1)-cumsum(top==0))[missing_speciations]
  lambda_b = sapply(df$brts[df$to==1]-0.000000001,speciation_rate,tree = df,pars = pars,model = model,soc=soc)
  if(length(lambda_b)==0) lambda_b = 1
  text = df$t_ext[df$to==1]-df$brts[df$to==1]
  mu = max(0,pars[1])

  intensity.temp = get(paste0("intensity.", model))
  inte = intensity.temp(tree=df,pars=pars)
  
  logg = -sum(inte)+sum(log(nb)+log(mu))-sum(mu*text)+sum(log(lambda_b))-sum(log(2*No+Ne))
  return(logg)
}



intensity.rpd1 <- function(tree, pars){
  mu = max(0,pars[1])
  n = tree$n
  lambda = pmax(0,pars[2] + pars[3] * n)
  wt = diff(c(0,tree$brts))
  brts_i = df$brts
  brts_im1 = c(0,df$brts[-nrow(df)])
  sigma_over_tree = n*(lambda)*(wt - (exp(-mu*max(df$brts))/mu)*(exp(mu*brts_i)-exp(mu*brts_im1)) )
  
  return(sigma_over_tree)
}

intensity <- function(x, tree, model, time0, pars,soc){
  max_time_for_continuity = min(tree$brts[tree$brts>time0])
  if(x==time0){
    val = 0
  }else{
    nh_rate <- function(wt){
      sum_speciation_rate(x=wt,tree = tree,pars = pars,model = model,soc=soc)*(1-exp(-(max(tree$brts)-wt)*pars[1]))
    }
    if(x != time0) x <- x-0.00000000001
    val = pracma:::quad(f = Vectorize(nh_rate),xa = time0,xb = x)
  }
  return(val)
}

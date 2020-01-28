# augmentaion (sampling) probability 

log_sampling_prob_nh <- function(df,pars,model="dd",...){
  if(is.null(df$t_ext)) df$t_ext  = df$bte
  b = max(df$brts)
  to = top = head(df$to,-1)
  to[to==2] = 1
  #initspec=1
  initspec=2
  N = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  brts_i = df$brts
  brts_im1 = c(0,df$brts[-nrow(df)])
  missing_speciations = (df$to == 1)
  nb = N[missing_speciations]
  No = c(1,1+cumsum(top==2))[missing_speciations]
  Ne = c(0,cumsum(top==1)-cumsum(top==0))[missing_speciations]
  lambda_b = sapply(df$brts[df$to==1]-0.000000001,speciation_rate,tree = df,pars = pars,model = model,...)
  if(length(lambda_b)==0) lambda_b = 1
  text = df$t_ext[df$to==1]-df$brts[df$to==1]
  mu = pars[1]
  inte = vector(mode = "numeric",length = length(brts_i))
  for(i in 1:length(brts_i)){
    inte[i] = intensity(x=brts_i[i],tree = df,model = model,time0 = brts_im1[i],pars = pars,...)
  }
  logg = -sum(inte)+sum(log(nb)+log(mu))-sum(mu*text)+sum(log(lambda_b))-sum(log(2*No+Ne))
  return(logg)
}


log_sampling_prob_nh_dd <- function(df,pars,model="dd",initspec=1){
  b = max(df$brts)
  to = top = head(df$to,-1)
  to[to==2] = 1
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  lambda = lambda.dd(pars,n)
  brts_i = df$brts
  brts_im1 = c(0,df$brts[-nrow(df)])
  # at missing speciation times
  missing_speciations = (df$to == 1)
  nb = n[missing_speciations]
  No = c(1,1+cumsum(top==2))[missing_speciations]
  Ne = c(0,cumsum(top==1)-cumsum(top==0))[missing_speciations]
  lambda_b = lambda[missing_speciations]
  text = df$t.ext[missing_speciations]
  mu = pars[3]
  logg = sum(-n * lambda * (brts_i-brts_im1-(exp(-b*mu)/mu) *  (exp(brts_i*mu)-exp(brts_im1*mu))  ))+sum(log(nb)+log(mu))-sum(mu*text)+sum(log(lambda_b))-sum(log(2*No+Ne))
  return(logg)
}

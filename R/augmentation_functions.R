augment_tree <- function(brts,pars,model,soc){
  mu = max(0,pars[1])
  brts = cumsum(-diff(c(brts,0)))
  b = max(brts)
  missing_branches = data.frame(speciation_time=NULL,extinction_time=NULL)
  cbt = 0
  while(cbt < b){
    tree = data.frame(brts = c(missing_branches$speciation_time,brts,missing_branches$extinction_time),
                      t_ext = c(missing_branches$extinction_time, rep(Inf,length(brts)+nrow(missing_branches))),
                      to = c(rep(1,nrow(missing_branches)),rep(2,length(brts)),rep(0,nrow(missing_branches))))
    tree = tree[order(tree$brts),]
    next_bt = min(tree$brts[tree$brts>cbt])
    
    tree$n = sapply(tree$brts,n_from_time,tree=tree,soc=soc)
    tree$pd = sapply(tree$brts,phylodiversity,tree=tree,soc=soc)
    lambda_max = max( sum_speciation_rate(cbt,tree,pars,model,soc=soc)*(1-exp(-mu*(b-cbt))) , sum_speciation_rate(next_bt,tree,pars,model,soc=soc)*(1-exp(-mu*(b-next_bt))))
   # lambda_max = lambda_max(cbt,tree,pars,model)
    ###
    u1 = runif(1)
    next_speciation_time = cbt - log(x = u1)/lambda_max
    if(next_speciation_time < next_bt){  ## 
      u2 = runif(1)
      pt = sum_speciation_rate(next_speciation_time,tree,pars,model,soc=soc)*(1-exp(-mu*(b-next_speciation_time)))/lambda_max
      if(u2<pt){
        extinction_time = next_speciation_time + truncdist::rtrunc(1,"exp",a = 0, b = (b-next_speciation_time),rate=mu)
        missing_branches = rbind(missing_branches,data.frame(speciation_time=next_speciation_time,extinction_time=extinction_time))
        ## tree = data.frame(brts=,t_ext=,to=)
        if(nrow(missing_branches)>1000){
          stop("Current parameters leds to a large number of species")
        }
      }
    }
    cbt = min(next_speciation_time,next_bt)
  }
  tree = data.frame(brts = c(missing_branches$speciation_time,brts,missing_branches$extinction_time),
                    t_ext = c(missing_branches$extinction_time, rep(Inf,length(brts)+nrow(missing_branches))),
                    to = c(rep(1,nrow(missing_branches)),rep(2,length(brts)),rep(0,nrow(missing_branches))))
  tree = tree[order(tree$brts),]
  
  tree$pd = sapply(tree$brts, phylodiversity, tree=tree,soc=soc)
  tree$n = sapply(tree$brts, n_from_time,tree=tree,soc=soc)
  
  return(list(tree=tree))
}

##############################################

lambda_max_rpd5c <- function(cbt,tree,pars){
  brts = tree$brts
  m=min(which(tree$brts>cbt))
  brts_im1 = c(0,brts)[m]
  
  n = tree$n[m]; Pt = c(0,tree$pd)[m]
  c2 = pars[4]*((n-1)/n); c3 = exp(-pars[1]*max(tree$brts)); c4 = pars[1]
  c1 = pars[2] + pars[3]*n + pars[4]*((Pt-n*brts_im1)/n)
  
  bt1 = cbt
  bt2 = min(tree$brts[tree$brts>cbt])
  
  #d1 = exp(-c4*bt1)*(c1*c3*c4 + c2*(exp(c4*bt1) + c3*(-1 + c4*bt1)))
  #d2 = exp(-c4*bt2)*(c1*c3*c4 + c2*(exp(c4*bt2) + c3*(-1 + c4*bt2)))
  d1 = -c3 * c4 * exp(c4 * bt1) *(c1 + c2 * bt1) + c2* (-c3) * exp(c4* bt1) + c2
  d2 = -c3 * c4 * exp(c4 * bt2) *(c1 + c2 * bt2) + c2* (-c3) * exp(c4* bt2) + c2
  
  
  if(d1>0 & d2>0){
    max_lambda = sum_speciation_rate(bt2,tree,pars,model,soc=soc)*(1-exp(-mu*(b-bt2)))
  }
  if(d1<0 & d2<0){
    max_lambda =  sum_speciation_rate(cbt,tree,pars,model,soc=soc)*(1-exp(-mu*(b-cbt)))
  }
  if((d1>0 & d2<0)){
    max_lambda = (-c2 - c1*c4 + c2* lambertWp((e^(1 + (c1*c4)/c2))/c))/(c2*c4)
    if(max_lambda<bt1 | max_lambda>bt2){
      stop("maximum out of boundaries")
    }else{
      print(paste("inflection point:",max_lambda))
    }
  }
  if((d1<0 & d2>0)){
    max_lambda = max( sum_speciation_rate(cbt,tree,pars,model,soc=soc)*(1-exp(-mu*(b-cbt))) , sum_speciation_rate(next_bt,tree,pars,model,soc=soc)*(1-exp(-mu*(b-next_bt))))
  }
  return(max_lambda)
  
}

lambda_max <- function(tm,tree,pars,model){
  lambda_max = get(paste0("lambda_max_", model))
  lm = lambda_max(tm,tree,pars)
  return(lm)
}

lambda_max_rpd1 <- function(cbt,tree,pars){
  m=min(which(tree$brts>cbt))
  brts_im1 = c(0,brts_i[-length(brts_i)])[m]
  n = tree$n[m]
  max_lambda = max(0,pars[2]+pars[3]*n)
  return(max_lambda)
  
}

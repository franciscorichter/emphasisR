augment_tree <- function(brts,pars,model,soc){
  mu = max(0,pars[1])
  brts = cumsum(-diff(c(brts,0)))
  b = max(brts)
  missing_branches = data.frame(speciation_time=NULL,extinction_time=NULL)
  cbt=0
  
  while(cbt < b){
    
    tree = data.frame(brts = c(missing_branches$speciation_time,brts,missing_branches$extinction_time),
                      t_ext = c(missing_branches$extinction_time, rep(Inf,length(brts)+nrow(missing_branches))),
                      to = c(rep(1,nrow(missing_branches)),rep(2,length(brts)),rep(0,nrow(missing_branches))))
    tree = tree[order(tree$brts),]
    
    next_bt = min(tree$brts[tree$brts>cbt])
    ###  sabado 15/2/20: cambiar esta parte a funciones max_model
    lambda_max = max( sum_speciation_rate(cbt,tree,pars,model,soc=soc)*(1-exp(-mu*(b-cbt))) , sum_speciation_rate(next_bt,tree,pars,model,soc=soc)*(1-exp(-mu*(b-next_bt))))
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
          #print(tree)
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
  
  tree$pd = sapply(tree$brts, function(x) emphasis:::phylodiversity(x, tree,soc=soc))
  tree$n = sapply(c(0,tree$brts[-length(tree$brts)]), n_from_time,tree=tree,soc=soc)
  
  return(list(tree=tree))
}

##############################################





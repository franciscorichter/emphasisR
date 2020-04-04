sim_brts = function(pars,model="rpd5c",ct,summary_sims=TRUE){
  sim=0
  number_of_empty_trees=-1
  sim_tree = get(paste0("sim.tree_",model))
  while(length(sim)==1){
    sim = sim_tree(pars = pars,ct=ct,soc=2,drop_extinct=TRUE)
    number_of_empty_trees = number_of_empty_trees + 1 
  }
  brts = max(sim)-c(0,sim)
  sim = brts[-length(brts)]
  if(summary_sims){
    fwt = sim[length(sim)]
    length_s = length(sim)
    li = list(dim = length_s,fwt = fwt,brts=sim,number_of_empty_trees=number_of_empty_trees)
  }else{
    li = list(brts=sim,number_of_empty_trees=number_of_empty_trees)
  }
  return(li)
}

### simulation of trees 

sim.tree_rpd1 <- function(pars,ct,soc,drop_extinct=TRUE){
  tree = NULL
  cbt = 0 
  N = soc
  mu = max(0,pars[1])
  while((cbt < ct)  &  (N >= soc)){
    lambda_ct = pars[2] + pars[3]*N  # diversity dependence only 
    rate_max = (lambda_ct+mu)*N
    u1 = runif(1)
    next_event_time = cbt-log(x = u1)/rate_max
    
    if(next_event_time < ct){
      u2 = runif(1)
      pt = rate_max/rate_max  # constant pice wise function = 1
      if(u2<pt){
        to = sample(c(1,0),size=1,prob=c(lambda_ct,mu)/(lambda_ct+mu))
        if(to == 1){
          N = N + 1
        }else{
          N = N - 1 
          if(N>=soc){
            ext_spec = sample(which(tree$to == 1 & tree$t_ext == Inf),1)
            tree$t_ext[ext_spec] = next_event_time
          }
        }
        tree = rbind(tree,data.frame(brts=next_event_time,to=to,t_ext=Inf))
        
      }
    }
    cbt = next_event_time
  }
  if(N==(soc-1)){
    tree = 0
  }else{
    tree = rbind(tree,data.frame(brts=ct,to=1,t_ext=Inf))
    if(drop_extinct){
      brts = tree$brts[is.infinite(tree$t_ext) & tree$to==1]
      tree = brts
    }
  }
  return(tree)
}


sim.tree_rpd5c <- function(pars,ct,soc,drop_extinct=TRUE){
  tree = NULL
  cbt = 0 
  N = soc
  mu = max(0,pars[1])
  P = 0
  while((cbt < ct)  &  (N >= soc)){
    lambda_mx = max(0,pars[2] + pars[3]*N  +  ((P + N*(ct-cbt))/N -cbt)*pars[4])
    rate_max = (lambda_mx+mu)*N 
    u1 = runif(1)
    next_event_time = cbt-log(x = u1)/rate_max
    
    if(next_event_time < ct){
      u2 = runif(1)
      lambda_ct = max(0,pars[2] + pars[3]*N  +  ((P+N*(next_event_time-cbt))/N - next_event_time)*pars[4])
      pt =( (lambda_ct+mu)*N )/rate_max  
      if(u2<pt){
        to = sample(c(1,0),size=1,prob=c(lambda_ct,mu)/(lambda_ct+mu))
        P = P + N*(next_event_time-cbt)
        if(to == 1){
          N = N + 1
        }else{
          N = N - 1 
          if(N>=soc){
            ext_spec = sample(which(tree$to == 1 & tree$t_ext == Inf),1)
            tree$t_ext[ext_spec] = next_event_time
            removed_branch_length = tree$t_ext[ext_spec] - tree$brts[ext_spec]
            P = P - removed_branch_length
          }
          
        }
        tree = rbind(tree,data.frame(brts=next_event_time,to=to,t_ext=Inf))
        
      }
    }
    cbt = next_event_time
  }
  if(N==(soc-1)){
    tree = 0
  }else{
    tree = rbind(tree,data.frame(brts=ct,to=1,t_ext=Inf))
    if(drop_extinct){
      brts = tree$brts[is.infinite(tree$t_ext) & tree$to==1]
      tree = brts
    }
  }
  return(tree)
}


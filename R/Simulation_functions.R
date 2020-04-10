sim_brts_bootstrap <- function(pars,model="rpd5c",ct,bootstrap_n=1){
  
  sim_tree = get(paste0("sim.tree_",model))
  SIMS = vector(mode = "list",length = bootstrap_n)
  sim = NULL
  for(i in 1:bootstrap_n){
    sim$brts = 0
    number_of_empty_trees = -1
    while(length(sim$brts) == 1){
      sim = sim_tree(pars = pars,ct=ct,soc=2)
      number_of_empty_trees = number_of_empty_trees + 1 
    }
    SIMS[[i]] = c(sim,list(number_of_empty_trees=number_of_empty_trees))

  }
    
  return(SIMS)
}


### simulation of trees 

sim.tree_rpd1 <- function(pars,ct,soc){
  tree = NULL
  cbt = 0 
  N = soc
  mu = max(0,pars[1])
  while((cbt < ct)  &  (N >= soc)){
    lambda_ct = max(0,pars[2] + pars[3]*N)  # diversity dependence only 
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
        tree = rbind(tree,data.frame(brts=next_event_time,to=to,t_ext=Inf,lambda=lambda_ct))
        
      }
    }
    cbt = next_event_time
  }
  if(N==(soc-1)){
    tree = 0
  }else{
    tree = rbind(tree,data.frame(brts=ct,to=1,t_ext=Inf,lambda=lambda_ct))
    brts = tree$brts[is.infinite(tree$t_ext) & tree$to==1]
  }
  return(list(tree=tree,brts=brts))
}


sim.tree_rpd5c <- function(pars,ct,soc){
  tree = NULL
  cbt = 0 
  N = soc
  mu = max(0,pars[1])
  P = 0
  while((cbt < ct)  &  (N >= soc)){
    lambda_mx = max(0,pars[2] + pars[3]*N  +  ((P + N*(ct-cbt)-cbt)/N )*pars[4])
    rate_max = (lambda_mx+mu)*N 
    u1 = runif(1)
    next_event_time = cbt-log(x = u1)/rate_max
    P = P + N*(next_event_time-cbt)
    if(next_event_time < ct){
      u2 = runif(1)
      lambda_ct = max(0,pars[2] + pars[3]*N  +  ((P+N*(next_event_time-cbt) - next_event_time)/N)*pars[4])
      pt =( (lambda_ct+mu)*N )/rate_max  
      
      if(u2<pt){
        to = sample(c(1,0),size=1,prob=c(lambda_ct,mu)/(lambda_ct+mu))
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
        tree = rbind(tree,data.frame(brts=next_event_time,to=to,t_ext=Inf,lambda=lambda_ct))
        
      }
    }
    cbt = next_event_time
  }
  if(N==(soc-1)){
    tree = 0
    brts = 0
  }else{
    tree = rbind(tree,data.frame(brts=ct,to=1,t_ext=Inf,lambda=lambda_ct))
    brts = tree$brts[is.infinite(tree$t_ext) & tree$to==1]
  }
  return(list(tree=tree,brts=brts))
}


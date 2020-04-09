n_spec <- function(cutime,brts){
  1+sum(brts<cutime) 
}

##  Ltt expected plots 

simulation_analysis <- function(pars,model,ct,n_it=100,expectedLTT=TRUE,divers_rate=TRUE,gLTT = ggplot(),ggLA = ggplot()){
  S = sim_brts_bootstrap(pars = pars,model = model,ct = ct, bootstrap_n = n_it)
  MEANS_N = MEANS_l = matrix(nrow=n_it,ncol=100)
  if(model=="rpd1") color = "blue"
  if(model=="rpd5c") color = "Darkgreen"
  for (i in 1:(n_it)){
    s = S[[i]]
    df = data.frame(time=s$brts,n=2:(length(s$brts)+1))
    df$time = -(max(df$time)-c(0,df$time[-nrow(df)]))
    gLTT = gLTT+geom_line(data=df,aes(x=time,y=n),colour=color,alpha=0.1)
    MEANS_N[i,] = sapply(seq(-ct,0,length=100), n_spec,brts=df$time)
    df = data.frame(time = seq(-ct,0,length=100), mean=colMeans(MEANS_N))
    gLTT = gLTT + geom_line(data = df, aes(x=time,y=mean),colour=color)
    df = s$tree
    ggLA = ggLA + geom_line(data=df,aes(x=brts,y=lambda),colour=color,alpha=0.1)
    MEANS_l[i,] = sapply(seq(0,ct,length=100), speciation_rate,tree=df,pars=pars,model=model,soc=2)
    df = data.frame(time = seq(0,ct,length=100), mean=colMeans(MEANS_l))
    ggLA = ggLA + geom_line(data = df, aes(x=time,y=mean),colour=color)
    
  }
  return(list(gLTT=gLTT,ggLA=ggLA))
}



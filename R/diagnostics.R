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


sample_size_determination <- function(f,n,tol=0.05,median=TRUE){    
  
  if(median){
    sn = unique(n)
    sn = sort(sn,decreasing = T)[1:2]
    fs = c(median(f[n==min(sn)]),median(f[n==max(sn)]))
    ns = c(min(sn),max(sn))
    hlp<-lm(fs~I(1/ns),weights = ns)
  }else{
    hlp<-lm(f~I(1/n),weights = n)
  }
  ab<-coef(hlp)
  df1 = data.frame(n=n,f=f)
  f.r <- ab[1]-tol
  n.r <- ab[2]/(f.r-ab[1])
  
  nn <- min(n):max(n)
  ff <- ab[1]+ab[2]/nn
  df2 = data.frame(nn=nn,ff=ff)
  
  theme_set(theme_minimal())
  ggbp = ggplot()+geom_boxplot(data=df1,aes(x=n,y=f,group=n,fill=n))+
    geom_line(data=df2,aes(x=nn,y=ff),colour="Darkgreen",size=1)+
    geom_hline(yintercept = ab[1],colour="purple")+
    # theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(legend.position = "none")+
    xlab("Monte Carlo Sample Size")+
    ylab("loglikelihood Estimation")+
    #  annotate(geom="text",x=2200, y=-62, label="Required Sample Size")+
    # geom_segment(data=data.frame(x=2300, y=-63, vx=2300, vy=-76), mapping=aes(x=x, y=y, xend=vx, yend=vy), arrow=arrow(), size=1, color="blue") + 
    geom_point(data=data.frame(x=n.r, y=f.r), mapping=aes(x=x, y=y), size=1, shape=21, fill="white")
  return(list(plot=ggbp,n.r=n.r))
}
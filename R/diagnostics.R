n_spec <- function(cutime,brts,soc=2){
  soc+sum(brts<cutime) 
}

n_spec_r <- function(cutime,brts,soc=2){
  soc+sum(brts>cutime) 
}
##  Ltt expected plots 

simulation_analysis <- function(pars, model, ct, n_it=100, expectedLTT=TRUE, divers_rate=TRUE, gLTT = ggplot(), ggLA = ggplot(), brts_to_compare ){

  S = sim_brts_bootstrap(pars = pars, model = model, ct = ct, bootstrap_n = n_it)
  MEANS_N = MEANS_l = matrix(nrow=n_it,ncol=100)
  
  if(model=="rpd1") color = "blue"
  if(model=="rpd5c") color = "Darkgreen"
  
  for (i in 1:(n_it)){
    
    s = S[[i]]
    df = data.frame(time = s$brts, n = 2:(length(s$brts)+1))
    gLTT = gLTT + geom_line(data=df,aes(x=time,y=n), colour=color, alpha=0.1)
    MEANS_N[i,] = sapply(seq(0,ct,length=100), n_spec, brts=df$time)
    
    df = s$tree
    ggLA = ggLA + geom_line(data=df,aes(x=brts,y=lambda), colour=color, alpha=0.1)
    MEANS_l[i,] = sapply(seq(0,ct,length=100), speciation_rate,tree=df, pars=pars, model=model, soc=2)
    
  }
  
  dfN = data.frame(time = seq(0,ct,length=100), mean=colMeans(MEANS_N))
  gLTT = gLTT + geom_line(data = dfN, aes(x=time,y=mean),colour=color)
  
  dfR =  sapply(seq(0,ct,length=100), n_spec, brts=brts)
  ltt_stat = sum(abs(dfR-dfN$mean)*0.1)
  
  df = data.frame(time = seq(0,ct,length=100), mean=colMeans(MEANS_l))
  ggLA = ggLA + geom_line(data = df[-1,], aes(x=time,y=mean),colour=color)
  
  
  
  return(list(gLTT=gLTT,ggLA=ggLA,dfN=dfN,ltt_stat=ltt_stat))
}


ltt_stat <- function(brts,est_seq){
  pr = sapply(seq(-ct,0,length=100), n_spec,brts=-brts)
  sum(abs(est_esq$mean-pr))*0.1
}

get_required_sampling_size <- function(M,tol=.05){
  n <- M$sample_size
  f<-  M$fhat
  hlp<-MASS:::rlm(f~I(1/n),weights = n)
  ab<-coef(hlp)
  
  f.r<-ab[1]-tol
  n.r<-ceiling(ab[2]/(f.r-ab[1]))
  return(n.r)
}

sample_size_determination <- function(f,n,tol=0.05){    

  hlp <- MASS:::rlm(f~I(1/n),weights = n)

  ab <- coef(hlp)

  f.r <- ab[1]-tol
  n.r <- ceiling(ab[2]/(f.r-ab[1]))
  
  nn <- min(n):max(c(n,n.r))
  ff <- ab[1]+ab[2]/nn
  df2 = data.frame(nn=nn,ff=ff)
  df1 = data.frame(n=n,f=f)
  
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
  return(list(plot=ggbp,n.r=n.r,fhat=ab[1]))
}


# Diagnostics function which recives an "emph object", which is the output of the emphasis function. 
emphasis_diagnostics <- function(MC){

  pars = MC$pars
  clade = MC$clade
  data(branching_times,package = "emphasis")
  ct = DDD_estimations[DDD_estimations$clade==clade,]$age
  brts = MC$brts
  sim = simulation_analysis(pars=pars,model="rpd5c",ct=ct,brts=brts)
  G = sim$gLTT + geom_step(data=data.frame(time=cumsum(c(0,-diff(brts))),n=2:(length(brts)+1)),aes(x=time,y=n))
  sim_ddd = simulation_analysis(pars=MC$pars_dd,model="rpd1",ct=ct,gLTT = G,ggLA = sim$ggLA,brts=brts)
  ltt_plot = sim_ddd$gLTT  + theme_classic() + ggtitle(clade) + ylab("Number of lineages (log)") + xlab("Time") + scale_y_log10()
  sr_plot = sim_ddd$ggLA  + theme_classic() + ggtitle(clade) + ylab("Speciation Rate") + xlab("Time") #+ geom_hline(yintercept = pars[1])
  
  mcem = MC$MCEM
  
  AIC_p = 8 - 2 * mean(mcem[mcem$sample_size==max(mcem$sample_size),]$fhat)
  AIC_d = 6 - 2 * DDD_estimations[DDD_estimations$clade==clade,]$fhat
  
  test1 = mean(mcem[mcem$sample_size==max(mcem$sample_size),]$fhat) > DDD_estimations[DDD_estimations$clade==clade,]$fhat
  
  #cat("AIC for PD model:",AIC_p)
  #cat("AIC for DD model:",AIC_d)
  test2 = AIC_p < AIC_d
  
  AICw = AICw(mean(mcem[mcem$sample_size==max(mcem$sample_size),]$fhat),DDD_estimations[DDD_estimations$clade==clade,]$fhat,4,3)
  # sampling size diagnostics
  samplingSize_requirment = MC$required_sample_size < MC$diag1$n.r   
  # convergence diagnostics
  conv_gg = ggplot(MC$MCEM) + geom_line(aes(x=1:nrow(MC$MCEM),y=fhat))
  M = MC$MCEM[MC$MCEM$sample_size==max(MC$MCEM$sample_size),]
  diffs_fhat = diff(M$fhat)
  df = as.data.frame(diffs_fhat)
  #n1 = ggplot()+geom_histogram(aes(x=diffs_fhat),binwidth=0.05)
  n1 = ggplot(df, aes(x = diffs_fhat)) + 
    geom_histogram(aes(y =..density..),
                   binwidth=0.2, 
                   colour = "black", 
                   fill = "white") +
    stat_function(fun = dnorm, args = list(mean = mean(df$diffs_fhat), sd = sd(df$diffs_fhat))) + theme_bw()
  #n2 = qqnorm(diffs_fhat);qqline(diffs_fhat, col = 2)
  n2 = ggplot2:::ggplot(df, aes(sample = diffs_fhat)) +  stat_qq() + stat_qq_line() + theme_bw()
  sht = stats:::shapiro.test(diffs_fhat)$p.value > 0.05
  
  test1 = mean(mcem[mcem$sample_size==max(mcem$sample_size),]$fhat) > DDD_estimations[DDD_estimations$clade==clade,]$fhat
  test2 = AIC_p < AIC_d
  test3 = sim$ltt_stat<sim_ddd$ltt_stat
  
  
  return(list(ltt_plot=ltt_plot,sr_plot=sr_plot,test1=test1,test2=test2,AICw=AICw,normality_test=sht,hist_norm=n1,qqplot=n2,ltt_PD=sim$ltt_stat,ltt_DD=sim_ddd$ltt_stat,test3=test3))
}



### non-tested


corr_pars <- function(mcem){
  
  MP = mcem[,c("par1","par2","par3")]
  
  
  #ggpairs(data=MP, palette = "RdYlGn", name = "rho", 
  #      label = FALSE, label_color = "black")
  
  data=MP
  data$K = (data$par1-data$par2)/data$par3
  names(data) = c("mu","lambda","beta","K")
  GGally:::ggpairs(data, columns = 1:ncol(data), title = "",  
                   axisLabels = "show")
}




# more utilities

n_from_time <- function(tm,tree,soc){
  # return N at tm.
  to = head(tree$to,-1)
  to[to==2] = 1
  n = c(soc,soc+cumsum(to)+cumsum(to-1))
  N = n[max(which(c(-1,tree$brts) < tm))]
  return(N)
} 

phylodiversity <- function(tm,tree,soc){
  i1<-tree$brts<=tm 
  i2<-tree$to==0&i1
  i3<-tree$t_ext%in%tree$brts[i2]
  dt<-diff(c(0,tree$brts[i1&!i2&!i3],tm))
  return(sum(dt*(soc:(length(dt)+soc-1))))
}

emphasis_bootstrap <- function(input,n_it=100,print=FALSE,file="bootstrap_temp.RData"){
  P=NULL
  for(i in 1:n_it){
    st = mcE_step(brts = input$brts, pars = input$pars,sample_size=input$sample_size,model=input$model,no_cores=input$cores,parallel=FALSE,soc=input$soc)
    if(print==TRUE){
      print(paste("iteration",i))
      print(paste("log-lik: ",log(st$fhat),sep=""))
      print(paste("took: ",st$E_time,sep=""))
    }
    P = rbind(P,data.frame(fhat=log(st$fhat),eitme=st$E_time,ss=input$sample_size))
    save(P,input,file=file)
  }
return(P)
}

data_to_table <- function(df,replicant,left,right){
  df = df[df$rep==replicant,]
  df = df[df$iteration %in% left:right,]
  summ = data.frame(lfhat = mean(df$fhat),sd_fhat=sd(df$fhat),mad_fhat=mad(df$fhat),replicant=replicant,par1=median(df$par1),par2=median(df$par2),par3=median(df$par3),par4=median(df$par4),E_time = sum(df$E_time)/60, M_time = sum(df$M_time)/60, sample_size=mean(df$sample_size))
  return(summ)
}

AIC <- function(LogLik,k){
  aic <- (2*k)-(2*LogLik)
  return(aic)
}

AICweights <- function(LogLik,k){
  IC <- AIC(LogLik,k)
  bestmodelIC <- min(IC)
  weights <- exp(-0.5*(IC-bestmodelIC))
  weights <- weights/sum(weights)
  return(weights)
}

AICw <- function(l1,l2,k1,k2){
  IC <- AIC(c(l1,l2),c(k1,k2))
  bestmodelIC <- min(IC)
  weights <- exp(-0.5*(IC-bestmodelIC))
  weights <- weights/sum(weights)
  return(weights[1])
}
vectors2phylo <- function(list){
  t=list$wt
  n=list$n
  E=list$E
  S=list$S
  ct=sum(t)
  newick = paste(sl[1],";",sep="")
  N=1
  identf = data.frame(Spec="aa",Time=0) # Labels of species
  for (i in 1:(length(t)-1)){
    # speciation
    sumt = sum(t[1:i])
    if( is.null(S)){
      BD = sample(1:N,1)
      species = as.character(identf[BD,1])
    }else{
      species = S[i]
    }
    if (E[i] == 1){
      ind = regexpr(species,newick)[1]-1
      atm = sumt-identf[which(identf[,1]==species),2]
      newick = paste(substr(newick,1,ind),"(",substr(newick,ind+1,ind+4),",",sl[i+1],"):",as.character(atm),substring(newick,ind+5),sep="")
      identf = rbind(identf,data.frame(Spec=substr(sl[i+1],1,2),Time=sumt))
      identf[identf$Spec == species,2] = sumt
      N = N+1
    }
    # extinction
    if (E[i]==0){
      ind = regexpr(species,newick)[1] + 2
      atm = sumt-identf[which(identf[,1]==species),2]
      identf = identf[!identf$Spec==species,]
      newick = paste(substr(newick,1,ind),as.character(atm),substring(newick,ind+2),sep="")
      N=N-1
    }
  }
  newick = compphyl(newi=newick,identf=identf,ct=ct)
  newick = read.tree(text=newick)
  return(newick)
}



phylo2tree <- function(tree){
  # to map newick trees into ther xxxx format
  ltt = ltt.plot.coords(tree)
  t = diff(ltt[,1])
  ltt = ltt[-1,]
  n = ltt[,2]
  E = diff(n)
  E[E==-1] = 0
  return(list(wt=t,to=E))
}

tree2phylo <- function(tree,initspec=1){
  wt=-diff(c(0,tree$brts))
  to=tree$to
  to[to==2] = 1
  ct=sum(wt)
  newick = paste(sl[1],";",sep="")
  N=1
  identf = data.frame(Spec="a",Time=0) # Labels of species
  for (i in 1:(length(wt)-1)){
    # speciation
    bt = sum(wt[1:i])
    BD = sample(1:N,1)
    species = as.character(identf[BD,1])  
    if (to[i] == 1){
      ind = regexpr(species,newick)[1]-1
      atm = bt-identf[which(identf[,1]==species),2]
      newick = paste(substr(newick,1,ind),"(",substr(newick,ind+1,ind+4),",",sl[i+1],"):",as.character(atm),substring(newick,ind+5),sep="")
      identf = rbind(identf,data.frame(Spec=substr(sl[i+1],1,2),Time=bt))
      identf[identf$Spec == species,2] = bt
      N = N+1
    }
    # extinction
    if (to[i]==0){
      ind = regexpr(species,newick)[1] + 2
      atm = bt-identf[which(identf[,1]==species),2]
      identf = identf[!identf$Spec==species,]
      newick = paste(substr(newick,1,ind),as.character(atm),substring(newick,ind+2),sep="")
      N=N-1
    }
  }
  newick = compphyl(newi=newick,identf=identf,ct=ct)
  newick = read.tree(text=newick)
  return(newick)
}

compphyl <- function(newi,identf,ct){
  #set to extant species to the present time
  identf[,1] = as.character(identf[,1])
  identf[,2] = ct-identf[,2]
  for(i in 1:length(identf[,1])){
    ind = regexpr(identf[i,1],newi)[1] + 2
    newi = paste(substr(newi,1,ind),as.character(identf[i,2]),substring(newi,ind+2),sep="")
  }
  return(newi)
}

sl = paste(letters[1],letters,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],letters,":0",sep="")
  sl = c(sl,ll)
}

# time calculation
get.time <- function(time,mode='sec'){
  dif = proc.time()-time
  ti = as.numeric(dif[3])
  if(mode == 'min')  ti = ti/60
  if(mode == 'hou') ti = ti/3600
  return(ti)
}

get.topologies <- function(M){
  if(M == 0)
  {
    return(NULL)
  }
  TO = matrix(nrow=2*M,ncol=1)
  TO[1,1] = 1
  for(i in 2:(2*M)){
    comb = ncol(TO)
    for(j in 1:comb){
      ns = sum(no.na(TO[,j]))
      ne = sum(1-no.na(TO[,j]))
      if(ns < M & ns > ne){ #extinction or speciation
        TO[i,j] = 1
        TO = cbind(TO,matrix(TO[,j],ncol=1))
        TO[i,ncol(TO)] = 0
      }
      if(ns == M & ns > ne){ #extinction
        TO[i,j] = 0
      }
      if(ns < M & ns == ne){ #speciation
        TO[i,j] = 1
      }
    }
  }
  return(TO)
}


### simulation of trees 

sim.tree <- function(pars, model,ct,soc){
  tree = data.frame(brts=0,to=1,t_ext=Inf, parent=0, child = 1)
  cbt = 0 
  N = soc
  mu = max(0,pars[1])
  ## sim waiting time,
  spec.cnt = soc
  while((cbt < ct)  &  (N > 0)){
    tmp.tree<-rbind(tree[-1,], data.frame(brts=ct,to=1,t_ext=Inf, parent=NA, child = NA))
    rate_max = max(sum_speciation_rate(cbt,tmp.tree,pars,model,soc = soc),sum_speciation_rate(ct,tmp.tree,pars,model,soc=soc))+mu*N
    u1 = runif(1)
    next_event_time = cbt-log(x = u1)/rate_max
    
    if(next_event_time < ct){
      u2 = runif(1)
      pt = (sum_speciation_rate(next_event_time,tmp.tree,pars,model,soc=soc)+mu*N)/rate_max
      if(u2<pt){
        l1 = speciation_rate(next_event_time,tmp.tree,pars = pars, model = model,soc=soc)
        to = sample(c(1,0),size=1,prob=c(l1,mu)/(l1+mu))
        #   print(l1/(l1+mu))
        if(to == 1){
          spec.cnt = spec.cnt + 1 
          current.spec = tree$child[tree$to==1 & is.infinite(tree$t_ext)]
          tree = rbind(tree,data.frame(brts=next_event_time,to=1,t_ext=Inf,parent=sample(current.spec,1), child=spec.cnt))
          N = N + 1
        }else{
          N = N - 1 
          #extant = which(is.infinite(tree$t_ext) & (tree$to != 0) & tree$brts < next_speciation_time)
          current.spec = tree$child[tree$to==1 & is.infinite(tree$t_ext)]
          extinction = sample(current.spec,1)
          tree = rbind(tree,data.frame(brts=next_event_time,to=0,t_ext=Inf,parent=extinction, child=NA))
          tree$t_ext[tree$child==extinction] = next_event_time
        }
      }
    }
    cbt = next_event_time
  }
  tree = rbind(tree,data.frame(brts=ct,to=1,t_ext=Inf, parent=NA, child = NA))
  return(tree)
}


remove.extinctions <- function(tree){
  extant_brts = tree$brts[tree$to==1 & is.infinite(tree$t_ext)]
  return(extant_brts)
}

multiplot <- function(lp, plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(lp, plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



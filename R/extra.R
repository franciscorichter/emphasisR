
n_from_time <- function(tm,tree,soc){
  # return N at tm.
  to = head(tree$to,-1)
  to[to==2] = 1
  n = c(soc,soc+cumsum(to)+cumsum(to-1))
  if(tm>max(tree$brts)){
    N = max(n)
  }else{
    N = n[max(which(c(-1,tree$brts) < tm))]
  }
  return(N)
} 

phylodiversity <- function(tm,tree,soc){
  i1<-tree$brts<=tm 
  i2<-tree$to==0&i1
  i3<-tree$t_ext%in%tree$brts[i2]
  dt<-diff(c(0,tree$brts[i1&!i2&!i3],tm))
  return(sum(dt*(soc:(length(dt)+soc-1))))
}


AIC_llik <- function(LogLik,k){
  aic <- (2*k)-(2*LogLik)
  return(aic)
}

AICw <- function(l1,l2,k1,k2){
  IC <- AIC_llik(c(l1,l2),c(k1,k2))
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

# time calculation
get.time <- function(time,mode='sec'){
  dif = proc.time()-time
  ti = as.numeric(dif[3])
  if(mode == 'min')  ti = ti/60
  if(mode == 'hou') ti = ti/3600
  return(ti)
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


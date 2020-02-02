# more utilities

number_of_species <- function(tree,tm=NULL){
  initspec = 1
  to = head(tree$to,-1)
  to[to==2] = 1
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  
  
  if(is.null(tm)){
    N=n
  }else{
    if(tm==max(tree$brts)){
      N = n[max(which(c(0,tree$brts) < tm))]
    }else{
      N = n[max(which(c(0,tree$brts) <= tm))]
    }
  }
  return(N)
}

phylodiversity <- function(tm,tree){
  i1<-tree$brts<=tm 
  i2<-tree$to==0&i1
  i3<-tree$t_ext%in%tree$brts[i2]
  dt<-diff(c(0,tree$brts[i1&!i2&!i3],tm))
  return(sum(dt*(2:(length(dt)+1))))
}

emphasis_bootstrap <- function(input,n_it=100){
  pp=NULL
  for(i in 1:n_it){
    print(paste("iteration",i))
    st = mcE_step(brts = input$brts, pars = input$pars,sample_size=input$sample_size,model=input$model,no_cores=input$cores,parallel=input$parallel)
    pp = rbind(pp,data.frame(fhat=log(st$fhat),eitme=st$E_time,ss=input$sample_size))
  }
return(pp)
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



sim.tree <- function(pars,CT,seed=1,model="dd"){
  set.seed(seed)
  mu = pars[2]
  N = 2
  brts = 0
  to = c(1,1)
  while(max(brts)<CT & N>0){
    if(model=="dd") lambda = lambda.dd(pars,N)
    if(model=="cr") lambda = lambda.cr(pars,N)
    sigma = N*(mu+lambda)
    wt = rexp(1,rate=sigma)
    tev = rbinom(n=1,size=1,prob= (lambda/(lambda+mu)) )
    brts = c(brts,max(brts)+wt)
    to = c(to,tev)
    if(tev==1) N = N+1
    if(tev==0 & max(brts)<CT){
      N = N-1
      ext = sample((1:length(to))[to==1],1)
      to[ext] = -1
    }
  }
  to = to[-length(to)]
  if(max(brts)>CT) brts[length(brts)]=CT
  tree = list(wt=diff(brts),to=to)
  return(tree)
}

prune.tree <- function(tree){
  brts = sum(tree$wt)-c(0,cumsum(tree$wt))
  brts = c(brts[1],brts[-length(brts)])
  nbrts = brts[tree$to==1]
  return(nbrts)
}

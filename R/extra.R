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
  return(sum(dt*(1:length(dt))))
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

vectors2phylo <- function(list,initspec=1){
  wt=list$wt
  to=list$to
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

W <- function (z, branch = 0)
{
  stopifnot(length(branch) == 1, is.numeric(z))
  if (anyNA(z)) {
    warning("Some values of ", deparse(substitute(z)), " are NA or NaN. ",
            "Returning 'NA' for these entries.")
    non.na.z <- z[!is.na(z)]
  }
  else {
    non.na.z <- z
  }
  W.non.na.z <- rep(NA, length(non.na.z))
  if (branch == 0) {
    W.non.na.z <- lamW::lambertW0_C(non.na.z)
  }
  else if (branch == -1) {
    if (any(is.infinite(z))) {
      warning("'Inf' is not a valid argument of the non-principal branch W",
              " (branch = -1).")
    }
    W.non.na.z <- lamW::lambertWm1_C(non.na.z)
  }
  else {
    stop("Branch was ", branch, "; must be either '0' or '-1'.")
  }
  if (length(W.non.na.z) == length(z)) {
    dim(W.non.na.z) <- dim(z)
    return(W.non.na.z)
  }
  else {
    W.z <- rep(NA, length(z))
    W.z[!is.na(z)] <- W.non.na.z
    W.z[is.nan(W.z)] <- NA
    dim(W.z) <- dim(z)
    return(W.z)
  }
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

nh_tree_augmentation_dd <- function(observed.branching.times,pars,model="dd",initspec = 1){
  
  b = max(observed.branching.times)
  brts = sort(brts)
  mu = pars[3]
  missing_branches = data.frame(speciation_time=NULL,extinction_time=NULL)
  N = initspec # current number of species
  cbt = 0 # current branching time
  while(cbt < b){
    lambda = lambda.dd(pars,N)
    all_bt = c(observed.branching.times,missing_branches$speciation_time,missing_branches$extinction_time)
    next_bt = min(all_bt[all_bt>cbt])
    next_speciation_time = rnhpp_dd(s=N*lambda,mu=mu,r=b-cbt,cbt=cbt,next_bt=next_bt)
    if(next_speciation_time < next_bt){ # add new species
      extinction_time = next_speciation_time + truncdist::rtrunc(1,"exp",a = 0, b = (b-next_speciation_time),rate=mu)
      missing_branches = rbind(missing_branches,data.frame(speciation_time=next_speciation_time,extinction_time=extinction_time))
      cbt = next_speciation_time
      N = N+1
    }else{
      cbt = next_bt
      if(next_bt %in% missing_branches$extinction_time){
        N = N-1
      }else{
        N = N+1
      }
    }
  }
  df = data.frame(brts = c(missing_branches$speciation_time,observed.branching.times,missing_branches$extinction_time),
                  bte = c(missing_branches$extinction_time, rep(Inf,length(observed.branching.times)+nrow(missing_branches))),
                  to = c(rep(1,nrow(missing_branches)),rep(2,length(observed.branching.times)),rep(0,nrow(missing_branches))))
  df = df[order(df$brts),]
  df$t.ext = df$bte-df$brts
  return(df)
  
}


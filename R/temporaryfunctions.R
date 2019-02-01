### temporary functions 

loggprob1 <- function(to,maxnumspec,ct,initspec=1,conf){
  n = c(initspec,initspec+cumsum(to)+cumsum(to-1))
  n = n[-length(n)]
  loggprob <- -log((maxnumspec+1))+lgamma(length(to)+1)-length(to)*log(ct)+log(probto(to))-sum(log(conf$N-conf$P))
}

possible.configurations  <- function(miss,obs){
  if (is.vector(obs)){
    tms <- obs
    to <- rep(1,length(obs))
    obs <- list(brts=tms,to=to)
  }
  mi = 1 #index for missing
  ob = 1 #index for observed
  
  # (number of) protected species
  protected = 1
  P<-NULL
  
  # (number of) current species
  currentspecies = 1
  N<-NULL
  
  # missing branching times
  brts.m = c(miss$brts,Inf)
  
  # set of sets of guardians
  guardians = list()
  
  # missing new species (starting to count at 1 above number of present species)
  n.obs = length(obs$brts)
  newspecies.m = n.obs+1
  
  # observed new species (starting to count at 1, so next one is 2)
  newspecies.o = 2
  
  # sampled tree
  tree<-list(brts=NULL,species=NULL,event=NULL)
  
  while (mi < length(brts.m) | ob < length(obs$brts)) {
    if(obs$brts[ob]<brts.m[mi]){ # observed speciation
      spec = obs$to[ob]
      #update tree
      tree$brts = c(tree$brts,obs$brts[ob])
      tree$species =c(tree$species,spec)
      tree$event = c(tree$event,newspecies.o)
      if(spec%in%protected){ # if protected, then both species become protected
        protected = c(protected,newspecies.o)
        currentspecies = c(currentspecies,newspecies.o)
      }else{ # if unprotected, then then (1) its guardian set disappears and (2) both become protected
        index = unlist(lapply(guardians,function(y,x){x%in%y},x=spec))
        if (sum(index)>0){
          index = which(index)
          n.guardians = length(guardians[[index]])
          guardians[[index]] = NULL
        }
        protected = c(protected,spec,newspecies.o)
        P = c(P,0) # it is weird, but we want to try 
        currentspecies = c(currentspecies,newspecies.o)
        N = c(N,n.guardians) # it is weird, but we want to try 
      }
      ob = ob + 1
      newspecies.o = newspecies.o + 1
    }else{ # missing event
      if(miss$to[mi]==1){ # missing speciation
        mspec = sample(c(currentspecies,currentspecies),1)
        index = which(mspec==protected)
        N = c(N,length(currentspecies))
        P=c(P,0)
        currentspecies = c(currentspecies,newspecies.m)
        #update tree
        tree$brts = c(tree$brts,miss$brts[mi])
        tree$species =c(tree$species,mspec)
        tree$event = c(tree$event,newspecies.m)
        if(sum(index)>0){ # if a protected species speciates, then it becomes unprotected and both guardians
          protected = protected[-index]
          guardians[[length(guardians)+1]] = c(mspec,newspecies.m)
        }else{ # if a unprotected species speciates, then ...
          index = unlist(lapply(guardians,function(y,x){x%in%y},x=mspec))          
          if (sum(index)>0){ #... if it is a guardian then new species becomes guardian
            index=which(index)
            guardians[[index]] = c(guardians[[index]],newspecies.m)
          } # ... if it is not a guardian then no changes to guardianship
        }
        newspecies.m = newspecies.m + 1
      } else { #missing extinction
        N = c(N,length(currentspecies))
        P = c(P,length(protected))
        available = setdiff(currentspecies,protected)
        missextinct = sample(c(available,available),1)
        if (missextinct<=n.obs){# if we selected the label of an extant species, we arbitrarily pick the label of a missing guardian
          index = which(unlist(lapply(guardians,function(y,x){x%in%y},x=missextinct)))         
          missextinct = setdiff(guardians[[index]],missextinct)[1]
        }
        #update tree
        tree$brts = c(tree$brts,miss$brts[mi])
        tree$species =c(tree$species,missextinct)
        tree$event = c(tree$event,0)
        index = unlist(lapply(guardians,function(y,x){x%in%y},x=missextinct))          
        if (sum(index)>0){ # if it is a guardian then ...
          index=which(index)
          if (length(guardians[[index]])>2){ # ... if the set is larger than 2, then take it out of guardian set
            guardians[[index]] = setdiff(guardians[[index]],missextinct)
          } else { # ... if guardian set is of size 2, remove guardian set and protect the remaining species
            protected = c(protected, setdiff(guardians[[index]],missextinct))
            guardians[[index]] = NULL
          }
        }
        currentspecies = setdiff(currentspecies,missextinct)
      }
      mi = mi + 1
    }
  }
  return(list(N=N,P=P,tree=tree))
}

probto <- function(to,p=0.5){
  posspec = c(0,cumsum(to==1))<(length(to)/2)
  posext = !(c(0,cumsum(to==1))==c(0,cumsum(to==0)))
  expo = sum(posspec & posext)
  prob = p^expo
  return(prob)
}

sim.miss1 <- function(maxnumspec=250,ct){
  S = sample(0:maxnumspec,1)
  brts = sort(runif(2*S,min=0,max=ct))
  to = sampletopology(S)
  tree = list(brts=brts,to=to)
  return(tree)
}

sampletopology <- function(S,p=0.5){
  to = NULL
  if(S>0){
    for(i in 1:(2*S)){
      if(sum(to==1)==sum(to==0)){
        prob = 1
      }
      if(sum(to==1)==S){
        prob = 0
      }
      to = c(to,rbinom(n=1,size=1,prob=prob))
      prob = p
    }
  }else{
    to = NULL
  }
  return(to)
}
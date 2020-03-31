###############################################################
# Uniform data augmentation importance sampler#
###############################################################

#' @keywords internal
lprobto <- function(to,p=0.5){
  posspec = c(0, cumsum(to == 1)) < (length(to)/2)
  posext = !(c(0, cumsum(to == 1)) == c(0, cumsum(to==0)))
  possibletotal = posspec & posext
  to_possible = to[possibletotal]
  logprob = sum(to_possible==1)*log(p)+sum(to_possible==0)*log(1-p)
  return(logprob)
}

#' @keywords internal
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
      to = c(to,stats::rbinom(n=1,size=1,prob=prob))
      prob = p
    }
  }else{
    to = NULL
  }
  return(to)
}

#' @keywords internal
sample_dim_prob <- function(d,max){
  factorial(2*d)/sum(factorial(seq(from = 0,to = 2*max,by=2)))
}

#' @keywords internal
sample.uniform <- function(brts,maxnumspec,single_dimension=NULL){
  if(is.null(single_dimension)){
 #   probs = sample_dim_prob(0:maxnumspec,maxnumspec)   not used
    S = sample(0:maxnumspec,1)
  }else{
    S = single_dimension
  }
  mbts.events = sim.branchingtimes.and.events(S=S ,ct = max(brts),p=0.5)
  df = data.frame(brts = c(brts, mbts.events$brts),
                  to =   c(rep(2, length(brts)), mbts.events$to))
  df = df[order(df$brts),]
  return(df)
}

#' @keywords internal
uniform_tree_augmentation <- function(brts,maxnumspec){
  S = sample(0:maxnumspec,1)
  mbts.events = emphasis:::sim.branchingtimes.and.events(S=S ,ct = max(brts),p=0.5)
  df = data.frame(brts=c(brts,mbts.events$brts),to=c(rep(2,length(brts)),mbts.events$to))
  df = df[order(df$brts),]
  df$t_exp = rep(Inf,nrow(df)) 
  missing_speciations = NULL
  for(j in 1:nrow(df)){
    if(df$to[j]==1){
      missing_speciations = c(missing_speciations,df$brts[j])
    }
    if(df$to[j]==0){
      sample_ext_time_index = sample(length(missing_speciations),1)
      df$t_exp[j] = missing_speciations[sample_ext_time_index]
      missing_speciations = missing_speciations[-sample_ext_time_index]
    }
  }
  
  return(df)
}

#' simulate branching times and events
#' @param S number of psecies
#' @param ct unknown
#' @param p probability
#' @return tree data frame
#' @keywords export
sim.branchingtimes.and.events <- function(S = S, 
                                          ct, 
                                          p){
  brts = sort( stats::runif(2*S, min = 0, max = ct))
  to = sampletopology(S, p = p)
  tree = list(brts = brts, to = to)
  return(tree)
}

#' @keywords internal
log.factor.samp.prob <- function(to){
  top = utils::head(to, -1)
  number.observed = c(1, 1+cumsum(top == 2))
  number.missing = c(0,cumsum(top==1) - cumsum(top == 0))
  factor = -sum(log((2 * number.observed+number.missing)[to == 1])) -
          sum(log(number.missing[to == 0]))
  return(factor)
}

#' @keywords internal
log.sampling.prob.uniform <- function(df,
                                      maxnumspec,
                                      initspec = 1,
                                      p = 0.5){
  ct = max(df$brts)
  to = top = df$to
  tom = top[top!=2]
  to[to==2] = 1
  num.miss = 2*sum((to==0))
  loggprob <- -log((maxnumspec+1)) +
                lgamma(num.miss+1) - 
              num.miss * log(ct) + 
    lprobto(tom,p = p) + 
    log.factor.samp.prob(top) 
  return(loggprob)
}

#' @keywords internal
log.samp.prob <- function(to,maxnumspec,ct,initspec=1,conf,p){
  n = c(initspec, initspec + cumsum(to) + cumsum(to-1))
  n = n[-length(n)]
  loggprob <- -log((maxnumspec + 1)) + 
                lgamma(length(to) + 1) - 
                  length(to)*log(ct) + 
                lprobto(to, p = p) - 
                sum(log(conf$N - conf$P))
}

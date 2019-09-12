# log-likelihood of full tree

test.llik.full.tree <- function(pars,brts,brts_exts){
  tree = brts2tree(brts,brts_exts)
  ll1 = -nllik.tree(tree=tree,pars=pars,topology = F,model = "cr",initspec = 1)
  ll2 = loglik_full_tree(brts = brts,brts_exts = brts_exts,la = pars[1],mu = pars[2],age = brts[1])
  if(all.equal(ll1,ll2)){
    print(paste("test full.lik passed"))
  }
  else{
    print(ll1)
    print(ll2)
    warning("criteria of llik not satisfied")
  }
}


brts2tree <- function(brts,brts_exts){
  if(brts[1]<0){ 
    wt = diff(c(brts,0))
    brts = cumsum(wt)
    age = sum(brts)
  }
  if(brts_exts[1,1]<0){
    brts2 = cumsum(diff(sort(c(brts_exts[,1],-age))))
    ext = cumsum(diff(sort(c(brts_exts[,2],-age))))
  }else{
    brts2 = brts_exts[,1]
    ext = brts_exts[,2]
  }
  to = c(rep(1,length(brts)+length(brts2)-1),2,rep(0,length(ext)))
  brts = c(sort(c(brts,brts2)),ext)
  tree = data.frame(brts=brts,to=to)
  tree = tree[order(tree$brts),]
  wt = diff(c(0,tree$brts))
  tree = list(wt=wt,to=tree$to[-length(tree$to)])
  return(tree)
}

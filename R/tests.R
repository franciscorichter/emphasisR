# log-likelihood of full tree



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

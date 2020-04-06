n_spec <- function(cutime,brts){
  1+sum(brts<cutime) 
}

##  Ltt expected plots 

expectedLTT <- function(pars, ct=15, model, n_it= 100,color="blue"){ # it should also include model
  BT = vector("list", n_it)
  g = ggplot() #+ geom_line(data=data.frame(time=-input$brts,n=1:(length(input$brts))),aes(x=time,y=n) )
  NE=ne=NULL
  MEANS_N = matrix(nrow=n_it,ncol=length(input$brts))
  for (i in 1:(n_it)){
    
    s = sim_brts(pars = pars,model = model,ct = ct)
    df = data.frame(time=-s$brts,n=2:(length(s$brts)+1))
    g = g+geom_line(data=df,aes(x=time,y=n),colour=color,alpha=0.1)
    NE = c(NE,s$number_of_empty_trees)
    BT[[i]] = s$brts
    ne = c(ne,length(s$brts))
    MEANS_N[i,] = sapply(-input$brts, n_spec,brts=-s$brts)
  }
  
  expect = colMeans(MEANS_N)
  df = data.frame(time=-input$brts,n=expect)
  g = g+geom_line(data=df,aes(x=time,y=n),colour=color,alpha=0.9)
  g  + geom_line(data=data.frame(time=-input$brts,n=1:(length(input$brts))),aes(x=time,y=n) )
  return(g)
}
   
 
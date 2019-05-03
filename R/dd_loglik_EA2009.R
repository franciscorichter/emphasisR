dd_loglik_rhsEA2009 <- function(t,x,pars)
{
  la <- pars[1]
  mu <- pars[2]
  nx <- sqrt(length(x))
  dim(x) <- c(nx,nx)
  xx <- matrix(0,nx+2,nx+2)
  xx[2:(nx+1),2:(nx+1)] <- x
  nl <- t(array(0:(nx - 1), dim = c(nx,nx)))
  dx <- mu * (nl + 1) * xx[1:nx,3:(nx+2)] + la * (nl - 1) * xx[2:(nx+1),1:(nx)] - (la + mu) * nl * xx[2:(nx+1),2:(nx+1)]
  dim(dx) <- c(nx^2,1)
  return(list(dx))
}

dd_loglik_EA2009 <- function(pars,lx,age)
{
  p <- array(0,dim = c(lx,lx))
  p[1,2] <- 1
  dim(p) <- c(lx^2,1)
  methode <- 'ode45'
  rtol <- 1E-10
  atol <- 1E-10
  y <- ode(p,c(0,age),dd_loglik_rhsEA2009,parms = pars,rtol = rtol,atol = atol,method = methode)[2,2:(lx^2 + 1)]
  dim(y) <- c(lx,lx)
  return(y)
}

if(0 == 1)
{
  library(deSolve)
  lx <- 100
  age <- 10
  pars <- c(0.2,0.1)
  pNENL <- dd_loglik_EA2009(pars = pars,lx = lx,age = age)
  pNL <- colSums(pNENL)
  ENL <- (0:(lx -1)) %*% pNL
  print(sum(pNENL))
  print(ENL)
  print(exp((pars[1] - pars[2]) * age))
}

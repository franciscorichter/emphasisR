\name{mle.st}
\alias{mle.st}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximun likelihood estimation for a set of trees
%%  ~~function to do ... ~~
}
\description{
Calculation of a maximun likelihood estimation of a set of trees
}
\usage{
mle.st(st, init_par = c(2, 1, 60), impsam = F, correction = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{st}{
list with the set of trees
}
  \item{init_par}{
initial parameters on the optimization routine
}
  \item{impsam}{
set to TRUE if an importance sampling correction is needed
}
  \item{correction}{
this is the correction value, still theory pending
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
MLE for lambda,mu and K
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
F. Richter M.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
n=100
S = vector(mode='list',length = n)
for(i in 1:n){
  S[[i]] = sim.tree()$tree
}
mle.st(S)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%s\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

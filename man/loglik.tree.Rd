\name{loglik.tree}
\alias{loglik.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
log-likelihood of a phylogenetic tree
}
\description{
it calculates the log-likelihood of a phylogenetic tree under the diversity-dependance model.
}
\usage{
llik(pars, tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars}{Parameters corresponding to the diversity-dependance model
%%     ~~Describe \code{pars} here~~
}
  \item{tree}{phylogenetic tree described as the three vectors wt,E and n.
%%     ~~Describe \code{tree} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
log-likelihood value
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
##   compare log likelihoods 

#dd
tree = data.frame(brts=c(2,3,5,6),to=c(1,0,2,2),t_ext=c(3,Inf,Inf,Inf))
tree$n = sapply(tree$brts,n_from_time,tree=tree,soc=2)
pars = c(0.2,0.6,-0.005)
loglik.tree.rpd1(pars,tree)

loglik.tree.numerical(pars,tree,"rpd1")


## pd
tree$pd <- sapply(tree$brts, phylodiversity, tree=tree,soc=2)
pars = c(0.2,0.6,-0.01,0.01)
loglik.tree.rpd5c(pars,tree)
loglik.tree.numerical(pars,tree,"rpd5c")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

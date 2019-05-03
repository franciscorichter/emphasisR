\name{nllik.tree}
\alias{nllik.tree}
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
s= sim.tree()
llik(c(0.8,0.1,40),s$tree)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

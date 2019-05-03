\name{mle.tree}
\alias{mle.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximun likelihood estimation for a phylogenetic tree
}
\description{
it calculates the MLE of a phylogenetic tree under the diversity-dependance model
}
\usage{
mle.tree(tree, init_par = c(2, 1, 60))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
phylogenetic tree described as a markov process wt,E,n
}
  \item{init_par}{
initial parameters of the optimization process
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
MLE values
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
s = sim.tree()
mle.tree(s$tree)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

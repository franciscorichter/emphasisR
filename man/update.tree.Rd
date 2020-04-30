\name{mcEM}
\alias{mcEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Monte-Carlo Expectation Maximization algorithm for phylogenetic trees
%%  ~~function to do ... ~~
}
\description{MCEM routine with fixed sampling size
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mcEM(input,tol=0.01,burnin=20,print_process=FALSE){
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{phylo object to be updated}
  \item{t_spe}{speciation time of the added missing species}
  \item{t_ext}{extinction time of the added missing species}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{wt}{waiting times of updated tree}
  \item{E}{binary vector with speciation and extinction events of updated tree}
  \item{n}{number of species of updated tree}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Francisco Richter M.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

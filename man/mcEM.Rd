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
mcEM(brts, 
     pars, 
     sample_size, 
     model, 
     soc, 
     tol,
     burnin,
     print_process,
     paralell,
     cores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{brts}{branching times of the phylogenetic tree to be analyzed}
  \item{pars}{initial parameters for the MCEM routine}
  \item{sample_size}{Monte Carlo sample size}
  \item{model}{Species diversification model}
  \item{soc}{Initial number of species}
  \item{tol}{Tolerance level for the loglikelihood}
  \item{burnin}{Number of iterations to drop before performing MCEM for estimation}
  \item{print_process}{print loglikelihood estimation at every MCEM iteration}
  \item{parallel}{TRUE to perform the mcE step in paralell computing}
  \item{cores}{Number of cores to use in case of using parallel computing}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{mcem}{MCEM chain with parameter values, loglikelihood estimation and sample size at every iteration}
}
\references{
[1] Richter, F., Haegeman, B., Etienne, R. S., & Wit, E. C. (2020). Introducing a general class of species diversification models for phylogenetic trees. Statistica Neerlandica.
}
\author{
Francisco Richter
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

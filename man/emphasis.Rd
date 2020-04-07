\name{emphasis}
\alias{emphasis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MCEM routine for parameters inference under diversity-dependance model
%%  ~~function to do ... ~~
}
\description{
Performs the emphasis framework for parameter estimation of diversification models in phylogenetic trees.
}
\usage{
emphasis(brts,model="rpd1",init_par,soc=2,sample_size=200,parallel=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{brts}{
Branching times of the phylogenetic tree 
}
  \item{init_par}{
Initial parameters for emphasis routine
}
  \item{soc}{
Number of initial species
}
  \item{sample_size}{
Initial Monte Carlo sample size
}
  \item{parallel}{
Describe if par of the code your be performed using parallel computing
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
[1] Richter, F., Haegeman, B., Etienne, R. S., & Wit, E. C. (2020). Introducing a general class of species diversification models for phylogenetic trees. Statistica Neerlandica.
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

#load data
data(branching_times)

clade = "Bucconidae"
brts =  get(paste0("brts_",clade))
pars = as.numeric(DD_est[DD_est$clade==clade,5:7])

mc = emphasis(brts,model="rpd5c",init_par=c(pars,0),soc=2,sample_size=100,parallel=F)

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

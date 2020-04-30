\name{emphasis}
\alias{emphasis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expectation-Maximization in Phylogenetic Analysis with Simulation and Importance Sampling
%%  ~~function to do ... ~~
}
\description{
Performs the emphasis framework for parameter estimation of diversification models in phylogenetic trees.
}
\usage{
emphasis(brts,
         model="rpd1",
         init_par,
         soc=2,
         em_tol=0.25,
         sample_size_tol=0.005,
         burnin_sample_size=200,
         burnin_iterations = 20,
         pilot_sample_size=c(200,600))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{brts}{
    Branching times of the phylogenetic tree 
}
  \item{model}{
    Species diversification model
}
  \item{init_par}{
    Initial parameters for emphasis routine
}
  \item{soc}{
    Number of initial species
}
  \item{em_tol}{
    Tolerance levels for the loglikelihood regarding the EM convergenge
}
  \item{sample_size_tol}{
    Tolerance levels for the sample size determination 
}
\item{burnin_sample_size}{
    MC sample size for the burnIn phase
}
\item{burnin_iterations}{
    Number of initial iterations for the burnIn phase, i.e to drop.
}
\item{pilot_sample_size}{
    vector with different sample sizes for MCEM routines inside the sample determination phase
}
\item{paralell}{
    Setting if the MC sampling should be done in paarallel or not
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
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

mc = emphasis(brts,model="rpd5c",init_par=c(pars,0))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

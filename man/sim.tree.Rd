\name{sim.tree}
\alias{sim.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Phylogenetic tree simulation
%%  ~~function to do ... ~~
}
\description{Phylogenetic tree simulation, given the model and parameters.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sim_phyl(ct=15, lambda0=0.8, mu0=0.1, K=40, model="dd",printEv=FALSE,seed=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ct}{Crown time}
  \item{lambda0,mu0}{Speciation and extinction rates parameter (they depends on the model)}
  \item{K}{Carring capacity (For diversity-dependance model only)}
  \item{Model}{Diversification rates model. The available models are 'dd' for diversity dependance and 'cr' for constant rate model}
  \item{PrintEv}{Prints the events while occurring in the tree}
  \item{seed}{Seed for the simulation.}
%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{tree}{A list containing 4 elements: wt are the waiting times of the tree, E is the binary vector representing speciations or extinctions, n is the number of species at each node of the tree and br are the branching times of the nodes.}
  \item{phylo}{is the phylo object corresponding to the tree}
  \item{tree.extant}{Same as tree, but removing extinct species}
  \item{extant.phylo}{phylo of the extant-species tree}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
[1] Generalizing species diversification models: Technical report
}
\author{F. Richter M.
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
st <- sim.tree()
plot(st$phylo)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


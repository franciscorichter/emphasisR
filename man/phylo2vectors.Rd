\name{phylo2vectors}
\alias{phylo2vectors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert phylo to vectors
%%  ~~function to do ... ~~
}
\description{Convert phylo objects onto a set of vectors of markov chains that carachterize the phylogenetic tree
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
phylo2vectors(tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{phylo object
%%     ~~Describe \code{tree} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{wt}{waiting times of the tree}
  \item{E}{binary vector specifying if speciation or extinction}
  \item{S}{complementary vector specifiying which species got extincted or speciated}
  \item{n}{vector containing the number of species}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
s = sim.tree()
all.equal(s$tree$wt,phylo2vectors(s$phylo)$wt)
all.equal(s$tree$E,phylo2vectors(s$phylo)$E)
all.equal(s$tree$n,phylo2vectors(s$phylo)$n)
}

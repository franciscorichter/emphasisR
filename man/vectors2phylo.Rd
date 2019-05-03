\name{vectors2phylo}
\alias{vectors2phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Convert markov chain sequences vectors onto phylo object
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vectors2phylo(list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list}{ List containing the vectors wt,n,E,S corresponding to vectors with waiting times, number of species, binary speciations/extinctions, and speciated/extincted species.
%%     ~~Describe \code{p} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Phylogenetic tree in consistency with the given vectors.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
F. Richter. M.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
s <-sim.tree()
plot(s$phylo)
plot(vectors2phylo(s$tree))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.


\name{update.tree}
\alias{update.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Update phylogenetic tree by adding 1 missing species
%%  ~~function to do ... ~~
}
\description{Given a vectorial form of the markovian sequence that characterizes the tree, and the speciation and extinction rate of the added missings species, it will return the new vectors of the new tree.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
update.tree(tree,t_spe,t_ext)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{phylo object to be updated}
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

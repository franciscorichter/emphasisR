\name{rec.tree}
\alias{rec.tree}
%- It should be called sim.extinct
\title{Reconstruction of extant species tree
}
\description{
Given a extant-species tree, it simulates extincted species on it (It should be called sim.extinct)
}
\usage{
rec.tree(tree, pars, model = "dd", seed = 0, adjustment = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
extant species tree given by marcov process wt, n and E.
}
  \item{pars}{
parameters corresponding to diversity-dependance model
}
  \item{model}{
choose model... currently diversity-dependance only
}
  \item{seed}{
random seed for replication
}
  \item{adjustment}{
I don´t know what is that
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
it will return a full extant+extinc phylogenetic tree...
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
 rec = rec.tree(tree=s$tree.extant, pars = c(0.8,0.1,40))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

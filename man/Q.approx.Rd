\name{Q.approx}
\alias{Q.approx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
log-likelihood of a set of phylogenetic trees
}
\description{
it calculates the log-likelihood of a set of phylogenetic trees under the diversity-dependance model.
}
\usage{
Q.approx(pars, setoftrees, impsam = F, correction = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars}{
%%     ~~Describe \code{pars} here~~
}
  \item{setoftrees}{
%%     ~~Describe \code{setoftrees} here~~
}
  \item{impsam}{
%%     ~~Describe \code{impsam} here~~
}
  \item{correction}{
%%     ~~Describe \code{correction} here~~
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
%% ~put references to the literature/web site here ~
}
\author{
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (pars, setoftrees, impsam = F, correction = 0)
{
    m = length(setoftrees)
    l = vector(mode = "numeric", length = m)
    w = vector(mode = "numeric", length = m)
    D = vector(mode = "numeric", length = m)
    for (i in 1:m) {
        s = setoftrees[[i]]
        l[i] = llik(pars = pars, tree = s) + length(s$tree$wt) *
            correction
        w[i] = s$weight
        D[i] = length(s$wt)
    }
    if (impsam) {
        weight = exp(w - correction * D)/sum(exp(w - correction *
            D))
        l = l * weight
    }
    L = sum(l)
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

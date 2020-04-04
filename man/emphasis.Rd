\name{emphasis}
\alias{emphasis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MCEM routine for parameters inference under diversity-dependance model
%%  ~~function to do ... ~~
}
\description{
We performs...
}
\usage{
emphasis(brts,model="rpd1",init_par,soc=2,sample_size=200,parallel=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
%%     ~~Describe \code{tree} here~~
}
  \item{init_par}{
%%     ~~Describe \code{init_par} here~~
}
  \item{cutTime}{
%%     ~~Describe \code{cutTime} here~~
}
  \item{printpar}{
%%     ~~Describe \code{printpar} here~~
}
  \item{impsam}{
%%     ~~Describe \code{impsam} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
}
  \item{ips}{
%%     ~~Describe \code{ips} here~~
}
  \item{n_trees}{
%%     ~~Describe \code{n_trees} here~~
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
[1] Flexible modelling framework for species assemblages
in XXX
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (tree, init_par, cutTime = 10, printpar = TRUE, impsam = FALSE,
    parallel = F, ips = 20, n_trees = 50)
{
    print("MCEM computations initiated. This might take a while.")
    TIME = proc.time()
    n_pars = length(init_par)
    pars = init_par
    Pars = matrix(nrow = (1000), ncol = n_pars)
    i = 1
    qt = 0
    q = vector(mode = "numeric", length = 1000)
    Pars[i, ] = pars
    while (get.time(TIME, mode = "min") < cutTime) {
        time = proc.time()
        if (printpar) {
            print(paste("iteration #", i, ":"))
            print(pars)
            print(paste("previous iteration took", qt, "secs"))
        }
        p = proc.time()
        trees <- sim.srt(tree = tree, pars = pars, n_trees = n_trees,
            parallel = parallel)
        pars = mle.st(trees, impsam = impsam)
        Pars[(i + 1), ] = pars
        qt = get.time(time)
        q[i] = qt
        i = i + 1
        if (i\%\%ips == 0) {
            n_trees = n_trees + 100
        }
    }
    Pars = Pars[1:i, ]
    q = q[1:i]
    return(list(pars = Pars, times = q))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

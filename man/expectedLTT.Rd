\name{expectedLTT}
\alias{expectedLTT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Maximun likelihood estimation for a set of trees
%%  ~~function to do ... ~~
}
\description{
Bootstrapping to simulate trees and obtain expected LTT plot given a model and its parameter values. 
}
\usage{
expectedLTT(pars, ct=15, model, n_it= 100,color="blue",g=ggplot(),all_trees=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pars}{
Parameters of the diversification rates. 
}
  \item{ct}{
Age of the clade
}
  \item{model}{
Species Miversification Model
}
  \item{n_it}{
Number of simulated trees
}
  \item{color}{
Color of the lines of the ltt plot
}
  \item{g}{
Plot to use to add new simulations
}
  \item{all_trees}{
True if want to add all trees ltt plots. False if only want to add the expected ltt plot 
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{g}{
A plot with the expected values. ggplot class.
}
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
#clade = "Caviidae"
#brts =  get(paste0("brts_",clade))
#pars = as.numeric(DD_est[DD_est$clade==clade,5:7])
#mc = emphasis(brts,model="rpd5c",init_par=c(pars,0),soc=2,sample_size=200,parallel=T)

g = ggplot() + geom_line(data=data.frame(time=-input$brts,n=1:(length(input$brts))),aes(x=time,y=n) )

g1 = expectedLTT(mc$pars,ct=input$brts[1],model = "rpd5c",n_it = 100,color = "Darkgreen",g=g,all_trees = F)
g2 = expectedLTT(as.numeric(DD_est[2,5:7]),ct=input$brts[1],model = "rpd1",n_it = 100,color = "blue",g=g1,all_trees = F)
g2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%s\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

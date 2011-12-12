\name{SpawMLResidMoran}
\alias{SpawMLResidMoran}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Testing the (residual) spatial dependency
}
\description{
  Test of spatial dependency in contextual level residuals from (spatially weighted)
  multilevel models. The spatial dependency is operationalized by Moran's I coefficient 
  for spatial auto-correlation. 
}
\usage{
SpawMLResidMoran(srawe.obj, distance.matrix, bandwidths,
  kernel = NULL, confidence.intervals = c(0.95)) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{srawe.obj}{
    \code{ResampledSpawMLOutput}, or a vector of residuals from a
    \code{SpawMLOutput} or matrix of residuals from a
    \code{ResampleExactSpawMLOutput}
}
\item{distance.matrix}{
  square matrix of dimension \eqn{n\times  n}{n by n}, where \eqn{n} is the number of
  contextual units.
}
\item{bandwidths}{
\code{vector} of bandwidths for which to compute Moran's I.
}
\item{kernel}{
      function applied to the distance matrix. By default \code{NULL}, in which
    case the kernel function \deqn{w_{ij}=f(d,h) =
      \left(\frac{1}{2}\right)^{d_{ij}^2/h^2}}{w_ij = f(d, h) =
  (1/2)^((d_ij/h)^2)} is used, where \eqn{w_{ij}, d_{ij}, h}{w_ij, d_ij,
  h} are elements of the weight matrix \eqn{\mathbf{W}}{W}, of the
    distance matrix \eqn{\mathbf{D}}{W} and the bandwidth \eqn{h}.
    User-supplied kernel functions have to take 2 arguments and return a
    matrix of the same dimension as the first argument.
}
\item{confidence.intervals}{
  \code{vector} of conficence intervals. Defaults  to \code{c(.95)}
  which corresponds to 95 \%.
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
Till Junge, Sandra Penic, Guy Elcheroth
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Residual spatial autocorrelation

## for residuals from ResampleAggregateSpawML
# load individual level data, remove collective guilt assignment for the
# data frame, and remove NA's 
data(traces_ind) 

traces_ind <- traces_ind[,-7]
traces_ind <- na.exclude(traces_ind)

# load contextual indicator for aggregation
data(traces_event) 

# load precise contextual indicator
data(homog_census) 

# load distance matrix and create weights
data(d_geo) 
geow_100 <- WeightMatrix(d_geo, bandwidth=100)

# perform ResampleAggregateSpawML
\donttest{acc_w_homog_100 <-
ResampleAggregateSpawML(individual.level.data=traces_ind,
  contextual.names=c("w_all", "Homog_00"), context.id="area",
  formula=cg_acc ~ victim_d + comb_d + male + age_1990 + high_school +
  higher_edu + (1|area), contextual.data=traces_event,
  precise.data=homog_census, contextual.weight.matrices=list(geow_100,
  geow_100), individual.weight.names="weight",
  aggregation.functions="weighted.mean", nb.resamples=10)  }


# perform SpawMLResidMoran for bandwidths 25, 50, 100, 200 

\donttest{MI_acc <- SpawMLResidMoran(srawe.obj=acc_w_homog_100,
  distance.matrix=d_geo, bandwidths=c(25,50,100,200)) }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

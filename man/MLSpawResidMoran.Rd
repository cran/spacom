\name{MLSpawResidMoran}
\alias{MLSpawResidMoran}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Testing the (residual) spatial dependency
}
\description{
  This step-4 function uses the contextual residuals of (spatially
  weighted) multilevel models to estimate the spatial dependency left
  unexplained by the model. The spatial dependency is operationalised by
  Moran's I coefficient for spatial auto-correlation, which can be
  recalculated for a range of different bandwidth values (which, in a
  similar way to spatial weights used for the creation of the context
  data, allow to parametrise the scale at which spatial dependency is
  being estimated). The part of explained spatial dependency can be
  obtained by comparison with the spatial dependency of the
  intercept-only, or with any other reference model, at the same
  bandwidth value.
}
\usage{
MLSpawResidMoran(ml.spaw.obj,
                 distance.matrix,
                 bandwidths,
                 kernel = NULL,
                 confidence.intervals = c(0.95)) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ml.spaw.obj}{
    \code{ResampleMLSpawOutput}, 
    \code{MLSpawExactOutput},
    or matrix of residuals}
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
\value{
  A \code{matrix} containing Moran's I's
}
\references{
Elcheroth, G., Penic, S., Fasel, R., Giudici, F., Glaeser, S., Joye, D.,
Le Goff, J.-M., Morselli, D., & Spini, D. (2012). Spatially weighted
context data: a new approach for modelling the impact of collective
experiences. \emph{LIVES Working Papers, 19}.
}
\author{
Till Junge, Sandra Penic, Guy Elcheroth
}
\examples{
# Residual spatial autocorrelation
## It is step-4 function
## for residuals from ResampleAggregateSpawML

## Data preparation

## load individual level data, remove collective guilt assignment for the
## data frame, and remove NA's
data(traces_ind)
traces_ind <- traces_ind[,-7]
traces_ind <- na.exclude(traces_ind)

## load contextual indicator for aggregation
data(traces_event)
\dontshow{traces_event=traces_event[seq(1, nrow(traces_event), by=100),]}

## load precise contextual indicator
data(homog_census)

## load distance matrix
data(d_geo)

## Step 1: Create spatial weights
geow.100 <- WeightMatrix(d_geo, bandwidth=100)


## Step 2: Compute spatially weighted aggregated contextual indicator
wv.agg.100 <- SpawAggregate(
              contextual.data = traces_event,
              context.id="area.name",
              contextual.names = "w_all",
              contextual.weight.matrices=geow.100,
              nb.resamples=5,
              aggregation.functions="weighted.mean",
              design.weight.names="weight",
              sample.seed=1)

## rename weighted variable name so it reflects the used weighting
## matrix
names(wv.agg.100) <- "w_all.100"

## Step 3: Perform ResampleMLSpawAggregate
acc_w_100 <-
ResampleMLSpawAggregate(
  individual.level.data=traces_ind,
  context.id="area.name",
  formula=cg_acc ~ victim_d + comb_d + male + age_1990 + high_school +
  higher_edu + (1|area.name) + w_all.100,
  aggregates=wv.agg.100,
  precise.data=NULL)

## Step 4: Perform MLSpawResidMoran for bandwidths 25, 50, 100, 200 

\dontrun{
MI_acc <- MLSpawResidMoran(ml.spaw.obj=acc_w_100,
  distance.matrix=d_geo,
  bandwidths=c(25,50,100,200)) }
\dontshow{MI_acc <- MLSpawResidMoran(ml.spaw.obj=acc_w_100,
  distance.matrix=d_geo, bandwidths=c(25,50)) }
  
## The results can be used for plotting spatial variogram
## See plot() and par() for details 
\dontrun{
plot(MI_acc[,1], xaxt="n", xlab="Bandwidth values", ylab="Moran's I", type="b")
axis(side=1, at=1:nrow(MI_acc), labels=rownames(MI_acc)) }
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "spatial" }% __ONLY ONE__ keyword per line

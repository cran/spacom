\name{spacom-package}
\alias{spacom-package}
\alias{spacom}
\docType{package}
\title{
Spatially weighted context data for multilevel modelling
%~~ package title ~~
}
\description{
The package provides tools to construct and exploit spatially weighted context data. Spatial weights are derived by a Kernel function from a user-defined matrix of distances between contextual units. Spatial weights can then be applied either to precise contextual measures or to aggregate estimates based on micro-level survey data, to compute spatially weighted context data. Available aggregation functions include indicators of central tendency, dispersion, or inter-group variability, and take into account survey design weights. The package further allows combining the resulting spatially weighted context data with individual-level predictor and outcome variables, for the purposes of multilevel modelling. An ad hoc stratified bootstrap resampling procedure generates robust point estimates for multilevel regression coefficients and model fit indicators, and computes confidence intervals adjusted for measurement dependency and measurement error of aggregate estimates. As an additional feature, residual and explained spatial dependency can be estimated for the tested models.

}
\details{
    \tabular{ll}{
    Package: \tab spacom\cr
    Type: \tab Package\cr
    Version: \tab 1.0-0\cr
    Date: \tab 2012-12-21\cr
    Authors: \tab Till Junge\cr
    \tab Sandra Penic\cr
    \tab Mathieu Cossuta\cr
    \tab Guy Elcheroth\cr
    Contributors: \tab Davide Morselli\cr
    \tab Stephanie Glaeser\cr
    Maintainer: \tab Till Junge <till.junge@gmail.com>\cr
    License: \tab GPL (>= 2)\cr
    LazyLoad: \tab yes\cr
    Depends: \tab methods, R (>= 2.14-0), spdep, foreach, lme4, boot, classInt,
    iterators, rgeos\cr
  }

  The package is structured by seven core functions, organised along four
logical steps. In each of the first three steps three types of input data
required for multilevel modelling with spatially weighted context data are
successively imported and transformed, to generate output objects required in
subsequent steps:

Step 1: Creation of spatial weighting matrices. User-specified
Kernel functions and bandwidth values (function \code{\link{WeightMatrix}}) are
applied to one or several user-defined distance matrices.

Step 2: Construction of spatially weighted context data. The spatial weights
generated in step 1 are combined with (a) imported precise contextual measures
(function \code{\link{SpawExact}}) and/or (b) imported micro-level survey data
to compute contextual estimates through a user-specified aggregation function
(function \code{\link{SpawAggregate}}).

Step 3: Multilevel modelling with spatially weighted context data. Imported
individual-level predictor and outcome variables are combined with the precise
(functions \code{\link{MLSpawExact}} and  \code{\link{ResampleMLSpawExact}}) or
aggregate (function \code{\link{ResampleMLSpawAggregate}}) context data
generated in step 2.

Finally, contextual residuals from multilevel models can be re-used in an
additional, fourth step to compute estimates of explained spatial dependency
(function \code{\link{MLSpawResidMoran}}).

The package furthermore provides a fast-track procedure, which condenses steps 1
to 3 and directly generates point estimates for regression coefficients across a
range of different bandwidth values, in order to explore the scale-dependency of
contextual effects prior to more intensive analyses (function
\code{\link{ExploreMLSpawExact}})

}

\author{
Till Junge, 
  Sandra Penic,
  Mathieu Cossuta,
  Guy Elcheroth


  Maintainer: Till Junge <till.junge@gmail.com>
  Contributors:   Stephanie Glaeser, Davide Morselli
}
\references{
Elcheroth, G., Penic, S., Fasel, R., Giudici, F., Glaeser, S., Joye, D., Le Goff, J.-M., Morselli, D., & Spini, D. (2012). Spatially weighted context data: a new approach for modelling the impact of collective experiences.\emph{LIVES Working Papers, 19}.
}
\keyword{ "package" }
\keyword{"spatial"}
  \seealso{

%~~ Optional links to other man pages, e.g. ~~
\code{\link[geospacom:geospacom-package]{geospacom}}
}
\examples{

## Testing the model described in Elcheroth et al (2012), page 25, Equation 9.4
## Collective guilt assignment (outcome) is predicted by a series of individual level predictors,
## and two contextual predictors: risk of war victimization weighted by geographical proximity and
## common territorial identification

## DATA PREPARATION
## Import data
## load individual level data
data(traces_ind)

## exclude variable "cg_acc" and missing values
traces_ind <- na.exclude(traces_ind[,-6])

## load contextual indicator for aggregation
data(traces_event)
\dontshow{traces_event=traces_event[seq(1, nrow(traces_event), by=100),]}

## load distance matrices
data(d_geo)
data(d_ident)

## EXPLORATORY ANALYSIS

## Examining the impact of the risk of war victimization weighted by geographical
## proximity weights with bandwidths of 25, 50 and 100 km

explo.wv <- ExploreMLSpawAggregate(individual.level.data=traces_ind,
 								   contextual.name="w_all",
 								   contextual.data=traces_event,
 								   context.id="area.name",
 								   formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+
 								   higher_edu+(1|area.name),
 								   distance.matrix=d_geo,
 								   multilevel.bandwidths=c(25, 50, 100),
							       design.weight.names = "weight",
 								   aggregation.function = "weighted.mean",
								   kernel = NULL)

## STEP 1: Create spatial weights

## geographical proximity, bandwidth=50
geow.50 <- WeightMatrix(d_geo, bandwidth=50)

## common territorial identification, bandwidth=2
w.id <- WeightMatrix(d_ident,bandwidth=2)

## STEP 2: Create spatially weighted aggregated indicator of risk of war victimization

wv.agg <- SpawAggregate(contextual.data=traces_event,
                        context.id="area.name",
                        contextual.names=c("w_all", "w_all"),
                        contextual.weight.matrices=list(geow.50, w.id),
                        aggregation.functions="weighted.mean",
                        design.weight.names="weight",
                        nb.resamples=2)
## If wanted, the names of aggregated indicators in wv.agg can be changed
names(wv.agg) <- c("wv.50", "wv.id")

## STEP 3: Perform multilevel analysis (Eq 9.4, Elcheroth et al (2012), p. 25)

mla.results <- ResampleMLSpawAggregate(individual.level.data=traces_ind,
									   context.id="area.name",
									   formula=cg_ass~victim_d+comb_d+male+age_1990+high_school+higher_edu+
									   (1|area.name)+wv.50+wv.id,
									   aggregates=wv.agg,
									   precise.data=NULL)

## STEP 4: Test the residual spatial dependency of the previous model
mla.results.moran <- MLSpawResidMoran(ml.spaw.obj=mla.results,
							          distance.matrix=d_geo,
							          bandwidths=c(25,50,100,200))

## Plotting spatial variogram with the values in mla.results.moran
## See plot() and par() for details
\dontrun{
plot(mla.results.moran[,1], xaxt="n", xlab="Bandwidth values", ylab="Moran's I", type="b")
axis(side=1, at=1:nrow(mla.results.moran), labels=rownames(mla.results.moran)) }

}

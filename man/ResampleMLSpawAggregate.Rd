\name{ResampleMLSpawAggregate}
\alias{ResampleMLSpawAggregate}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Spatially weighted multilevel analysis with stratified resampling for
  aggregated context data
}
\description{

  This step-3 function performs multilevel analyses with spatially weighted
  context data based on aggregate survey estimates. The context data created in
  step 2 are combined with additional individual outcome and (optional)
  predictor variables, to test a user-defined model. An ad hoc stratified
  bootstrap resampling procedure generates robust point estimates for regression
  coefficients and model fit indicators, and computes confidence intervals
  adjusted for measurement dependency and measurement error of the aggregate
  estimates. For each tested model, contextual residuals can be stored for later
  re-use.

}
\usage{
ResampleMLSpawAggregate(individual.level.data,
 context.id,
 formula,
 aggregates,
 precise.data = NULL,
 confidence.intervals = c(0.95),
 individual.sample.seed = NULL,
 ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{individual.level.data}{
   A \code{data.frame} containing the individual level predictor 
   and outcome variables. The data must be numeric (\code{logical},
    \code{integer}, or plain \code{numeric}). There must not be any 
    missing values, \code{NaN}s, \code{NULL}s or \code{NA}s. 
  }
  \item{context.id}{
    The name of the context ID variable. This variable allows matching
  contextual units from different data sets (\code{precise.data},
  \code{individual.level.data}). Values can by of
  any type but need to be defined consitently across data sets.
  }
  \item{formula}{
   Formula description of the model.The formula is handed down to \code{lmer} 
   function of the package \code{lme4}. Refer to the documentation of \code{lmer}
   for details.
  }  
  \item{aggregates}{
    A \code{\link{SpawAggregateOutput-class}} object generated using
    \code{\link{SpawAggregate}}, see examples.
  }
  \item{precise.data}{
    A \code{data.frame} containing the precise contextual data.There must 
    not be any missing values, \code{NaN}s, \code{NULL}s or \code{NA}s. 
    By default set to \code{NULL}.
  }
  \item{confidence.intervals}{
    \code{vector} of confidence intervals. Defaults  to \code{c(.95)}
    which corresponds to 95 \%.
  }
  \item{individual.sample.seed}{
    Seed used to generate the random sampling of the individual data
    Is one of three things
    \enumerate{
      \item \code{NULL}, in which case whatever the current random seed is
      is used
      \item an \code{integer}, which will be used to set the random
      seed. This allows reproducible random samples
      \item a saved \code{.Random.seed}  which allows reproducible random
      samples as well. The reason why both 2 and 3 are present is because
      \code{.Random.seed} can be saved a posteriori.
    }
    Defaults to \code{NULL}. 
  }
  \item{...}{
    All additional named arguments are handed through to the function
    \code{lme:lmer} (for example, \code{REML=TRUE}).
  }
}
\value{
  An object of class \code{SpawAggregateOutput-class}.
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
## Spatially weighted multilevel analysis, with resampled both individual
## level indicators and contextual predictors for aggregation. It may
## include non-resampled precise contextual indicator. 
## It is step-3 function
## Predicting collective guilt acceptance. 

## Data preparation
## load individual level data, remove collective guilt assignment from the
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


## Step 2: Create spatially weighted aggregated and precise indicators

wv.agg.100 <- SpawAggregate(
              contextual.data = traces_event,
              context.id="area.name",
              contextual.names = "w_all",
              contextual.weight.matrices=geow.100,
              nb.resamples=5,
              aggregation.functions="weighted.mean",
              design.weight.names="weight",
              sample.seed=1)
            
homog.100 <- SpawExact(precise.data=homog_census,
                       context.id="area.name",
                       contextual.names="Homog_00",
                       contextual.weight.matrices=geow.100)

## rename weighted variable names so they reflect the used weighting
## matrix
names(wv.agg.100) <- "w_all.100"
names(homog.100)[2] <- "Homog.100"


##  Step 3: Perform ResampleMLSpawAggregate
acc_w_homog_100 <-
ResampleMLSpawAggregate(
  individual.level.data=traces_ind,
  context.id="area.name",
  formula=cg_acc ~ victim_d + comb_d + male + age_1990 + high_school +
  higher_edu + (1|area.name) + w_all.100 + Homog.100,
  aggregates=wv.agg.100,
  precise.data=homog.100)  

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
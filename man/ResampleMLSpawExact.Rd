\name{ResampleMLSpawExact}
\alias{ResampleMLSpawExact}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Spatially weighted multilevel analysis with stratified resampling for
  exact context data
}
\description{
  This step-3 function performs multilevel analyses with spatially
  weighted context data based on precise macro-level measures. The
  context data created in step 2 are combined with additional individual
  outcome and (optional) predictor variables, to test a user-defined
  model. An ad hoc stratified bootstrap resampling procedure generates
  robust point estimates for regression coefficients and model fit
  indicators, and computes confidence intervals adjusted for measurement
  dependency. For each tested model, contextual residuals can be stored
  for later re-use.
}
\usage{
ResampleMLSpawExact(individual.level.data,
 context.id,
 formula,
 precise.data,
 confidence.intervals = c(0.95),
 nb.resamples = 1000,
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
    This variable allows matching
  contextual units from different data sets (\code{precise.data},
  \code{individual.level.data}). Values can by of
  any type but need to be defined consitently across data sets.

  }
  \item{formula}{
   Formula description of the model.The formula is handed down to \code{lmer} 
   function of the package \code{lme4}. Refer to the documentation of \code{lmer}
   for details.
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
  \item{nb.resamples}{
    number of resamples to be evaluated. By default set to 1000.
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
  An object of class
  \code{SpawAggregateOutput-class}.
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
# Spatially weighted multilevel analysis, with resampled individual
# level indicators and precise contextual indicator. 

## It is step-2 function

## Data preparation

##  load individual level data, remove collective guilt assignment from the
## data frame, and remove NA's 
data(traces_ind)
traces_ind <- traces_ind[,-7]
traces_ind <- na.exclude(traces_ind)

## load precise contextual indicator
data(homog_census)

## load distance matrix 
data(d_geo)

## Step 1: Create spatial weights
geow.100 <- WeightMatrix(d_geo, bandwidth=100)

## Step 2: Create spatially weighted precise contextual indicator
homog.100 <- SpawExact(precise.data=homog_census,
                       context.id="area.name",
                       contextual.names="Homog_00",
                       contextual.weight.matrices=geow.100)
## rename weighted variable names so they reflect the used weighting
## matrix
names(homog.100)[2] <- "Homog.100"                       

## Step 3: Perform ResampleMLSpawExact
acc_homog100 <-
  ResampleMLSpawExact(
    individual.level.data=traces_ind,
    context.id="area.name",
    formula=cg_acc ~ victim_d + comb_d + male + age_1990 + high_school +
    higher_edu + Homog.100 + (1|area.name), precise.data=homog.100,
    nb.resamples=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}% __ONLY ONE__ keyword per line

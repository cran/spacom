\name{MLSpawExact}
\alias{MLSpawExact}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Standard multilevel analysis applied to spatially weighted context data
}
\description{
  Integrates spatially weighted context data as higher-level predictors in multilevel analysis and produces conventional statistical estimates.
}
\usage{
MLSpawExact(individual.level.data,
 context.id,
 formula,
 precise.data=NULL,
 verbose=TRUE,
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

  \item{precise.data}{
       A \code{data.frame} containing the precise contextual data.There must
    not be any missing values, \code{NaN}s, \code{NULL}s or \code{NA}s.
    By default set to \code{NULL}.
  }

  \item{verbose}{
  if \code{TRUE}, progress information and an estimate of the remaining
  computing time is shown at the screen. Defaults to \code{TRUE}
  }

  \item{...}{
    All additional named arguments are handed through to the function
    \code{lme:lmer} (for example, \code{REML=TRUE}).
  }

}
\value{
  An object of class \code{MLSpawExactOutput-class}
}
\references{
  Elcheroth, G., Penic, S., Fasel, R., Giudici, F., Glaeser, S., Joye, D.,
  Le Goff, J.-M., Morselli, D., & Spini, D. (2012). Spatially weighted
  context data: a new approach for modelling the impact of collective
  experiences. \emph{LIVES Working Papers, 19}.
}
\author{
  Till Junge, Sandra Penic, Mathieu Cossuta, Guy Elcheroth
}
\examples{
## Spatially weighted multilevel analysis, with standard estimates of
## standard errors.
## It is step-3 function
## Model with two contextual predictors, treated as precise,
## predicting collective guilt acceptance.

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

## Step 2: Create spatially weighted contextual indicators

## from micro-level survey data (treated as precise)
aggregate <- SpawAggregate(contextual.data=traces_event,
                           context.id="area.name",
                           contextual.names='w_all',
                           contextual.weight.matrices=geow.100,
                           aggregation.functions="weighted.mean",
                           design.weight.names="weight",
                           nb.resamples=0)

## from precise indicator
weighted.homog.census <- SpawExact(precise.data=homog_census,
                                   context.id="area.name",
                                   contextual.names="Homog_00",
                                   contextual.weight.matrices=geow.100)

## merge aggregated and weighted data
context.data <- merge(aggregate, weighted.homog.census, by="area.name")


## Step 3: Perform MLSpawExact with two spatially weighted indicators
acc_w_homog_100 <-
  MLSpawExact(individual.level.data=traces_ind,
              context.id="area.name",
              formula=cg_acc ~ victim_d + comb_d + male + age_1990 +
              high_school + higher_edu + (1|area.name) + w_all.1 + Homog_00.1,
              precise.data=context.data)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}% __ONLY ONE__ keyword per line

\name{ResampleMLSpawOutput-class}
\Rdversion{1.1}
\docType{class}
\alias{ResampleMLSpawOutput-class}
\alias{print,ResampleMLSpawOutput-method}
\alias{show,ResampleMLSpawOutput-method}

\title{Class \code{"ResampleMLSpawOutput"}}
\description{
  Output object of the function \code{ResampleMLSpaw}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ResampleMLSpawOutput", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{individual.sample.seed}:}{
      Object of class \code{"integer"} }
    \item{\code{fixed}:}{Object of class \code{"data.frame"}}
    \item{\code{random.var}:}{Object of class \code{"data.frame"}}
    \item{\code{model.fit}:}{Object of class \code{"data.frame"}}
    \item{\code{ranefs}:}{Object of class \code{"matrix"}}
    \item{\code{nb.resamples}:}{Object of class \code{"integer"}}
    \item{\code{betas}:}{Object of class \code{"data.frame"}}
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "ResampleMLSpawOutput")}}
    \item{show}{\code{signature(object = "ResampleMLSpawOutput")}}
  }
}
\author{
  Till Junge
}
\examples{
## Perform ResampleMLSpawExact

## Data preparation (see ResampleMLSpawExact)

data(traces_ind)
traces_ind <- traces_ind[,-7]
traces_ind <- na.exclude(traces_ind)
data(homog_census)
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

## acc_homog100 is an object of class ResampleMLSpawOutput
class(acc_homog100)

## to assess standardized fixed effects coefficients
acc_homog100@betas

## to assess non-standardized fixed effects coefficients
acc_homog100@fixed

## to assess only median of non-standardized fixed effects coefficients
\dontrun{
acc_homog100@fixed["50\%"] }

## to assess random effects
acc_homog100@random.var
}
\keyword{classes}
\keyword{spatial}
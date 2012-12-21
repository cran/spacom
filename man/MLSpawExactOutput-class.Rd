\name{MLSpawExactOutput-class}
\Rdversion{1.1}
\docType{class}
\alias{MLSpawExactOutput-class}
\alias{AIC,MLSpawExactOutput-method}
\alias{BIC,MLSpawExactOutput-method}
\alias{fixef,MLSpawExactOutput-method}
\alias{print,MLSpawExactOutput-method}
\alias{ranef,MLSpawExactOutput-method}
\alias{show,MLSpawExactOutput-method}
\alias{summary,MLSpawExactOutput-method}
\alias{VarCorr,MLSpawExactOutput-method}

\title{Class \code{"MLSpawExactOutput"}}
\description{
  Output object of the function \code{MLSpawExact}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MLSpawExactOutput", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{lme}:}{Object of class \code{"mer"} }
    \item{\code{beta}:}{Object of class \code{"numeric"} }
  }
}
\section{Methods}{
  \describe{
    \item{AIC}{\code{signature(object = "MLSpawExactOutput")}:
      \code{\link[lme4:mer-class]{mer}} of \pkg{lme4} }
    \item{BIC}{\code{signature(object = "MLSpawExactOutput")}: 
      \code{\link[lme4:mer-class]{mer}} of \pkg{lme4}  }
    \item{fixef}{\code{signature(object = "MLSpawExactOutput")}: 
      \code{\link[lme4:mer-class]{mer}} of \pkg{lme4}  }
    \item{print}{\code{signature(x = "MLSpawExactOutput")}}
    \item{ranef}{\code{signature(object = "MLSpawExactOutput")}: 
      \code{\link[lme4:mer-class]{mer}} of \pkg{lme4}  }
    \item{show}{\code{signature(object = "MLSpawExactOutput")}}
    \item{summary}{\code{signature(object = "MLSpawExactOutput")}}
    \item{VarCorr}{\code{signature(x = "MLSpawExactOutput")}: 
      \code{\link[lme4:mer-class]{mer}} of \pkg{lme4}  }
	 }
}
\author{
  Till Junge
}
\examples{
## Perform MLSpawExact

## Import and prepare data (see details in MLSpawExact)

data(traces_ind)
traces_ind <- traces_ind[,-7]
traces_ind <- na.exclude(traces_ind)
data(traces_event)
\dontshow{traces_event=traces_event[seq(1, nrow(traces_event), by=100),]}
data(homog_census)
data(d_geo)

## Step 1: Create spatial weights
geow.100 <- WeightMatrix(d_geo, bandwidth=100)

## Step 2: Create spatially weighted contextual indicators
aggregate <- SpawAggregate(contextual.data=traces_event,
                           context.id="area.name",
                           contextual.names='w_all',
                           contextual.weight.matrices=geow.100,
                           aggregation.functions="weighted.mean",
                           design.weight.names="weight",
                           nb.resamples=0)

## Step 3: Perform MLSpawExact with the spatially weighted indicator
acc_w100 <- MLSpawExact(individual.level.data=traces_ind,
              			context.id="area.name",
              			formula=cg_acc ~ victim_d + comb_d + male + age_1990 +
              			high_school + higher_edu + (1|area.name) + w_all.1,
              			precise.data=aggregate)
              			
## acc_w100 is an object of class MLSpawExactOutput              			
class(acc_w100)

## to assess standardized coefficients
acc_w100@beta

## to assess non-standardized fixed coefficients
acc_w100@lme@fixef

## to assess the deviance of the fitted model
acc_w100@lme@deviance

## for details about slots within @lme see mer-class{lme4}
}
\keyword{classes}
\keyword{spatial}
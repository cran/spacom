\name{ExploreMLSpawExact}
\alias{ExploreMLSpawExact}
\alias{ExploreMLSpawAggregate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Exploratory multilevel analysis with spatially weighted context data
}
\description{
  This exploratory function provides a fast-track procedure, which
  condenses steps 1 to 3. On the basis of precise contextual measures
  and/or micro-level survey data for estimating contextual indicators, a
  distance matrix between contextual units, and individual outcomes and
  (optional) predictor variables, it directly generate point estimates
  for regression coefficients of a user-defined multilevel model, for a
  range of different bandwidth values. This function can be used to
  first explore the scale of relevant contextual effects, before
  parametrising the model more precisely, calculating (computationally
  expensive) confidence intervals and estimating spatial dependency for
  specific models, on a step-by-step basis.
}
\usage{
ExploreMLSpawExact(individual.level.data,
 contextual.name,
 context.id,
 formula,
 distance.matrix,
 multilevel.bandwidths,
 precise.data,
 kernel = NULL)

ExploreMLSpawAggregate(individual.level.data,
 contextual.name,
 contextual.data,
 context.id, formula,
 distance.matrix,
 multilevel.bandwidths,
 design.weight.names = NULL,
 aggregation.function = "weighted.mean",
 kernel = NULL,
 additional.args = NULL)
}
\arguments{
  \item{individual.level.data}{
A \code{data.frame} containing the individual level predictor and outcome variables.
   The data must be numeric (\code{logical}, \code{integer}, or plain
    \code{numeric}). There must not be any missing values, \code{NaN}s,
    \code{NULL}s or \code{NA}s.

}
  \item{contextual.name}{
  A name of contextual variable to be weighted.
}
  \item{contextual.data}{
 A \code{data.frame} consisting of micro-level data to generate contextual
    indicators by aggregation and containing a column named \code{context.id}
    with the context ID variable. This is an individual level dataframe.
     The data must be numeric (\code{logical}, \code{integer}, or plain
    \code{numeric}). There must not be any missing values, \code{NaN}s,
    \code{NULL}s or \code{NA}s. By default set to \code{NULL}.
}
  \item{context.id}{
    The name of the context ID variable. This variable allows matching
  contextual units from different data sets (\code{contextual.data},
  \code{individual.level.data}, \code{distance.matrix}). Values can by of
  any type but need to be defined consitently across data sets.

}
  \item{formula}{
Formula description of the model.The formula is handed down to \code{lmer} function of
   the package \code{lme4}. Refer to the documentation of \code{lmer} for details.
}
  \item{distance.matrix}{
 A square matrix of dimension \eqn{n\times  n}{n by n}, where \eqn{n} is the number of
    contextual units.
}
  \item{multilevel.bandwidths}{
A \code{vector} of scalar numeric values specifying the bandwidth \eqn{h}
}
  \item{precise.data}{
A \code{data.frame} containing the precise contextual data.There must
    not be any missing values, \code{NaN}s, \code{NULL}s or \code{NA}s.
    By default set to \code{NULL}.
}
\item{design.weight.names}{
  A name of optional design weight at the individual level
  used for aggregation (for example, for a weighted mean). By default set to \code{NULL}.

}
\item{aggregation.function}{
  A name of aggregation function. Function takes either
  \enumerate{
    \item 1 argument in which case the corresponding individual design
    weight is \code{NULL},
    \item 2 arguments in which case the second argument is taken
    from the corresponding individual design weight.
    Defaults to \code{"mean"}.
  }
}
\item{kernel}{
  A function applied to the distance matrix. By default \code{NULL}, in which
    case the kernel function \deqn{w_{ij}=f(d,h) =
      \left(\frac{1}{2}\right)^{d_{ij}^2/h^2}}{w_ij = f(d, h) =
  (1/2)^((d_ij/h)^2)} is used, where \eqn{w_{ij}, d_{ij}, h}{w_ij, d_ij,
  h} are elements of the weight matrix \eqn{\mathbf{W}}{W}, of the
    distance matrix \eqn{\mathbf{D}}{D} and the bandwidth \eqn{h}.
    User-supplied kernel functions have to take 2 arguments and return a
    matrix of the same dimension as the first argument.
  }
  \item{additional.args}{
   For aggregation functions which take additional arguments (that is in
  addition to the data to aggregate and design weights), they can be
  specified here.
  \code{additional.args} defaults to \code{NULL}
  }

}
\details{
\code{ExploreMLSpawExact} performs exploratory multilevel analysis with a precise
spatially weighted contextual indicator. \code{ExploreMLSpawAggregate} performs analysis
with an aggregated spatially weighted contextual indicator. Both functions provide only
conventional statistical estimates and accept only one contextual predictor.
}
\value{
  A \code{list} of
  \code{\link{MLSpawExactOutput-class}} objects
}
\references{
  Elcheroth, G., Penic, S., Fasel, R., Giudici, F., Glaeser, S., Joye,
  D., Le Goff, J.-M., Morselli, D., & Spini, D. (2012). Spatially
  weighted context data: a new approach for modelling the impact of
  collective experiences. \emph{LIVES Working Papers, 19}. 
}
\author{
Till Junge, Sandra Penic, Guy Elcheroth
}
\examples{
#### ExploreMLSpaw
# Exploratory spatially weighted multilevel with standard estimates of
# standard errors. Accepts only one contextual predictor. Predicting
# collective guilt acceptance. Precise contextual predictor (ethnic homogeneity)
# weighted with geographical proximity weights, h=50,100,200.

# load individual level data, remove collective guilt assignment from the
# data frame, and remove NA's
data(traces_ind)

traces_ind <- traces_ind[,-7]
traces_ind <- na.exclude(traces_ind)

# load precise contextual indicator
data(homog_census)

# load distance matrix
data(d_geo)

# perform ExploreMLSpawExact
acc_homog <- ExploreMLSpawExact(
  individual.level.data=traces_ind,
  contextual.name="Homog_00",
  context.id="area.name",
  formula=cg_acc ~ victim_d + comb_d + male + age_1990 + high_school +
  higher_edu + (1|area.name),
  precise.data=homog_census,
  distance.matrix=d_geo,
  multilevel.bandwidths=c(50, 100))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }

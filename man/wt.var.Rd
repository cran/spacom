\name{aggregation functions}
\alias{wt.var}
\alias{wt.Atkinson}
\alias{wt.RS}
\alias{wt.Theil}
\alias{wt.gini}
\alias{wt.gini.categ}
\alias{wt.mean}
\alias{wt.sd}
\alias{wt.gini.group}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Aggregation functions provided by spacom
}
\description{
  Spacom provides a set of aggregation functions which can be used with
  the functions \code{\link{SpawAggregate}} and
  \code{\link{ExploreMLSpawAggregate}}. See \bold{Details} for
  descriptions of each function. Note that you are not restricted
  to these.

  User supplied aggregation functions can be used if they have the form
  \deqn{x_w = f(x, w, ...)} where \eqn{x} is the data to be aggregated,
  \eqn{w} the weights to be applied. Any number of additional arguments
  may also be used (see argument \code{additional.args} in
  \code{\link{SpawAggregate}} for details on additional arguments. The
  return value \eqn{x_w} must be a scalar value (i.e. not a vector).
}
\usage{
wt.sd(data, weights=rep(1, length(data)))
wt.var(data, weights=rep(1, length(data)))
wt.gini(data, weights=rep(1, length(data)))
wt.gini.categ(data, weights=rep(1, length(data)))
wt.Theil(data, weights=rep(1, length(data)))
wt.RS(data, weights=rep(1, length(data)))
wt.Atkinson(data, weights=rep(1, length(data)))
wt.gini.group(data, weights=rep(1, length(data)), groups)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A \code{vector} of data to be aggregated
}
  \item{weights}{
    A \code{vector} of weights. These weights are computed by combining
    the spatial weights and design weights for each context.
}
\item{groups}{
  name of the column used to group data (only \code{wt.group.gini})
}
}

\details{
  \enumerate{
    \item{\code{wt.var(data, weights)}}{
      Computes the weighted variance according to
      \deqn{{\rm Var}(x) = \frac{\sum_{i=1}^n (w_i\cdot(x_i - \mu)^2)}
        {\sum_{i=1}^n w_i}}{Var(x) = sum(w*(data-mu))/sum(w)}
      where the weighted mean is \eqn{\mu}{mu}
      \deqn{\mu = \frac{\sum_{i=1}^n w_i\cdot x_i}{\sum_{i=1}^n w_i}}{mu = weighted.mean(data, weights)}
    }
    \item{\code{wt.sd(data, weights)}}{
      Computes the weighted standard deviation according to
      \deqn{\sigma = \sqrt{{\rm Var}(x)}}{sd = sqrt(wt.var(data, weights))}
    }
    \item{\code{wt.gini}}{
      The weighted Gini coefficient is given by
      \deqn{G = \frac{1}{2\mu}\sum_{i,j}w_iw_j\left|x_i-x_j\right|}{}

      where \eqn{\mu}{mu} is
      \deqn{\mu = \sum_{i}w_ix_i}{mu = sum(w*x)}
    }
    \item{\code{wt.gini.categ}}{
      Spacom also allows for the computation of inequality indicators for
      categorical variables.
      In this case, one has no \eqn{x_i} value but only the weighted
      frequencies  \eqn{f_i} for each category. For that case, the Gini
      becomes
      \deqn{
        G = 1-\sum_if_i^2
      }{
        G = 1-sum(f^2)
      }
      The weighted frequencies are computed from \code{data} and \code{weights}
    }
    \item{\code{wt.gini.group}}{
      get from guy
    }
  }
}

\value{
  A scalar value of type \code{numeric}. For user-supplied aggregation
  functions, this can potentially be of a differenet type.
}
\references{
Elcheroth, G., Penic, S., Fasel, R., Giudici, F., Glaeser, S., Joye, D.,
Le Goff, J.-M., Morselli, D., & Spini, D. (2012). Spatially weighted
context data: a new approach for modelling the impact of collective
experiences. \emph{LIVES Working Papers, 19}.
}
\author{
Mathieu Cossuta, Till Junge, Sandra Penic, Guy Elcheroth
}
\examples{
## Creating indicators of spatially weighted inequalities 

## load individual level data and remove missings
data(traces_ind)
traces_ind <- na.exclude(traces_ind)

## create variable "simulated income" in the dataset
## the variable is created only in demonstrative purposes, 
## and it is not part of the TRACES dataset

traces_ind$sim_inc <- rnorm(nrow(traces_ind), 1000, 100)

## Step 1: Load distance matrix and create weights
data(d_geo)
geow.50 <- WeightMatrix(d_geo, bandwidth=50)

## Step 2: Spatially weighted gini for simulated income
si.gini <- SpawAggregate(contextual.data=traces_ind,
                         context.id="area.name",
                         contextual.names="sim_inc",
                         contextual.weight.matrices=geow.50,
                         aggregation.functions="wt.gini",
                         design.weight.names=NULL,
                         nb.resamples=5)

## Step 2: Spatially weighted gini for groups (Spatially weighted inequalities
## in simulated income for men and women)  

si.gini.gr <- SpawAggregate(contextual.data=traces_ind,
                            context.id="area.name",
                            contextual.names="sim_inc",
                            contextual.weight.matrices=geow.50,
                            aggregation.functions="wt.gini.group",
                            additional.args="male",
                            design.weight.names=NULL,
                            nb.resamples=5)   
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

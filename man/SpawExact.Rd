\name{SpawExact}
\alias{SpawExact}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create and describe spatially weighted precise contextual indicators 
}
\description{

This step-2 function should be used to create spatially weighted context data
when contextual indicators are directly available as precise macro-level
measures (e.g., population values based on census or registry data, or precise
geographic measures). The generated context data can be adjusted by user-defined
population weights.

}
\usage{
SpawExact(precise.data,
          context.id,
          contextual.names,
          contextual.weight.matrices,
          population.weight.names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{precise.data}{
    A \code{data.frame}  consisting of precise contextual indicators and
    a context ID variable
  }
  \item{context.id}{
      The name of the context ID variable. This variable allows matching
  contextual units from different data sets (\code{precise.data},
  \code{contextual.weight.matrix}). Values can by of
  any type but need to be defined consitently across data sets.

  }
  \item{contextual.names}{
    A \code{list} of column names of contextual variables to be weighted
  }
  \item{contextual.weight.matrices}{
    A \code{list} of square weight
    matrices of dimension \eqn{n\times n}{n by n}, where \eqn{n} is the number of
    contextual units to be applied to the variable specified in
    \code{contextual.names}, respectively.
  }
  \item{population.weight.names}{
    A \code{list} of column names to be used as population weights. By
    default \code{NULL} which means that no population weighting is
    performed. If only one name is provided, the corresponding column is
    used as a weight for all contextual names
  }
}
\value{
  A \code{data.frame} containing the spatially weighted \code{contextual.names}. 
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
## Creating spatially weighted precise contextual indicator - spatially weighted 
## ethnic homogeneity across TRACES areas.
## This is step-2 function

## Step 1: Load distance matrix and create weights

data(d_geo)
geow.50 <- WeightMatrix(d_geo, bandwidth=50)

## Step 2: Load dataframe with precise contextual indicator and 
## create spatially weighted indicator

data(homog_census)

homog.g50 <- SpawExact(precise.data=homog_census,
                       context.id="area.name",
                       contextual.names="Homog_00",
                       contextual.weight.matrices=geow.50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

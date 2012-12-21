\name{WeightMatrix}
\alias{WeightMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create spatial weighting matrices}

\description{

  This step-1 function creates a matrix of spatial weights on the basis of a
  user-defined distance matrix, Kernel function, and bandwidth value. The
  distance matrix needs to specify a value for each of the possible \eqn{n\times
  n}{n by n} binomials that correspond to n contextual units. It can be either
  symmetric or asymmetric. In principle, its diagonal, corresponding to the
  distance of each unit with itself, should be composed of zero values. A Kernel
  function proposed by default generates spatial weights that tend toward 1 for
  distances substantially lower than the bandwidth value, toward 0 for distances
  substantially higher than the bandwidth value and toward 0.5 for distances
  approaching the bandwidth value.

}

\usage{
WeightMatrix(distance.matrix, bandwidth, kernel = NULL,
  moran = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance.matrix}{
    %%     ~~Describe \code{distance.matrix} here~~
    square matrix of dimension \eqn{n\times  n}{n by n}, where \eqn{n} is the number of
    contextual units.
}
  \item{bandwidth}{
    %%     ~~Describe \code{bandwidth} here~~
    scalar numeric value specifying the bandwidth \eqn{h}
}
  \item{kernel}{
    %%     ~~Describe \code{kernel} here~~
    function applied to the distance matrix. By default \code{NULL}, in which
    case the kernel function \deqn{w_{ij}=f(d,h) =
      \left(\frac{1}{2}\right)^{d_{ij}^2/h^2}}{w_ij = f(d, h) =
  (1/2)^((d_ij/h)^2)} is used, where \eqn{w_{ij}, d_{ij}, h}{w_ij, d_ij,
  h} are elements of the weight matrix \eqn{\mathbf{W}}{W}, of the
    distance matrix \eqn{\mathbf{D}}{D} and the bandwidth \eqn{h}.
    User-supplied kernel functions have to take 2 arguments and return a
    matrix of the same dimension as the first argument.
}
  \item{moran}{
    %%     ~~Describe \code{moran} here~~
    a logical value specifying whether the proximity weights matrix
    should have zeros in the diagonal. By default set to \code{FALSE}.
}
}
\value{
  A weights matrix of the same dimension as \code{distance.matrix}.
}
\references{
Elcheroth, G., Penic, S., Fasel, R., Giudici, F., Glaeser, S., Joye, D.,
Le Goff, J.-M., Morselli, D., & Spini, D. (2012). Spatially weighted
context data: a new approach for modelling the impact of collective
experiences. \emph{LIVES Working Papers, 19}.
}
\author{
  %%  ~~who you are~~
  Till Junge, Sandra Penic, Guy Elcheroth
}
\examples{
# creating geographical proximity weight, with bandwidth h=50

data(d_geo) 
geow_50 <- WeightMatrix(distance.matrix=d_geo, bandwidth=50)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spacom }
\keyword{ spatial }% __ONLY ONE__ keyword per line

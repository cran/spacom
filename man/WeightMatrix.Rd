\name{WeightMatrix}
\alias{WeightMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create weight matrices}
\description{
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Creates a weight matrix corresponding to a given distance matrix, bandwidth and
  kernel function for spatially weighted multilevel analysis.
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
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A weights matrix of the same dimension as \code{distance.matrix}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  %%  ~~who you are~~
  Till Junge, Sandra Penic, Guy Elcheroth
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# creating geographical proximity weight, with bandwidth h=50

data(d_geo) 
\donttest{geow_50 <- WeightMatrix(distance.matrix=d_geo, bandwidth=50)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spacom }
\keyword{ SpawML }% __ONLY ONE__ keyword per line

\name{ExploreSpawML}
\alias{ExploreSpawML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exploratory multilevel analysis with spatially weighted context data
}
\description{
Performs exploratory multilevel analysis with spatially weighted context data for a series of different bandwidth values.
}
\usage{
ExploreSpawML(individual.level.data, 
 contextual.name, 
 contextual.data, 
 context.id, formula, 
 distance.matrix, 
 multilevel.bandwidths, 
 precise.data = NULL, 
 individual.weight.names = NULL, 
 aggregation.function = "mean", 
 kernel = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{individual.level.data}{
A \code{data.frame} containing the individual level data. The data must
    be numeric (\code{logical}, \code{integer}, or plain
    \code{numeric}). There must not be any missing values, \code{NaN}s,
    \code{NULL}s or \code{NA}s. 
}
  \item{contextual.name}{
  A name of contextual variable to be weighted.
}
  \item{contextual.data}{
 A \code{data.frame} containing the micro-level data to generate contextual 
    indicators by aggregation. May be \code{NULL}, in which case the individual data is
    used. The data must be numeric (\code{logical}, \code{integer}, or plain
    \code{numeric}). There must not be any missing values, \code{NaN}s,
    \code{NULL}s or \code{NA}s. By default set to \code{NULL}. 
}
  \item{context.id}{
The name of the context ID variable.
}
  \item{formula}{
Formula description of the model.
}
  \item{distance.matrix}{
 A square matrix of dimension \eqn{n\times  n}{n by n}, where \eqn{n} is the number of
    contextual units.
}
  \item{multilevel.bandwidths}{
A \code{vector} of scalar numeric values specifying the bandwidth \eqn{h}
}
  \item{precise.data}{
   A \code{data.frame} containing the precise contextual data. May be
    \code{NULL}, but may not contain any missing values, like \code{NaN}s,
    \code{NULL}s or \code{NA}s. By default set to \code{NULL}, except for
    \code{ResampleExactSpawML} which requires precise data to be specified.
}
\item{individual.weight.names}{
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
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
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
#### ExploreSpawML
# Exploratory spatially weighted multilevel with standard estimates of
# standard errors. Accepts only one contextual predictor. Predicting
# collective guilt acceptance. Contextual predictor weighted with
# geographical proximity weights, h=50,100,200. 

# load individual level data, remove collective guilt assignment for the
# data frame, and remove NA's 
data(traces_ind) 

traces_ind <- traces_ind[,-7]
traces_ind <- na.exclude(traces_ind)

# load precise contextual indicator
data(homog_census) 

# load distance matrix
data(d_geo) 

# perform ExploreSpawML
\donttest{acc_homog <- ExploreSpawML(individual.level.data=traces_ind,
  contextual.name="Homog_00", contextual.data=NULL, context.id="area",
  formula=cg_acc ~ victim_d + comb_d + male + age_1990 + high_school +
  higher_edu + (1|area), precise.data=homog_census, distance.matrix=d_geo,
  multilevel.bandwidths=c(50,100,200)) }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

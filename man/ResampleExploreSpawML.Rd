\name{ResampleExploreSpawML}
\alias{ResampleExploreSpawML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exploratory spatially weighted multilevel analysis with stratified resampling
}
\description{
Performs exploratory multilevel analysis with spatially weighted context data for 
a series of different bandwidth values. Provides robust point estimates and 
standard errors based on stratified bootstrap resampling.
}
\usage{
ResampleExploreSpawML(individual.level.data, contextual.name,
  contextual.data, context.id, nb.resamples, formula, distance.matrix,
  multilevel.bandwidths, individual.weight.name = NULL,
  aggregation.function = "mean", confidence.intervals = c(0.95),
  individual.sample.seed = NULL, contextual.sample.seed = NULL, kernel =
  NULL) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{individual.level.data}{
A \code{data.frame} containing the individual level data. The data must
    be numeric (\code{logical}, \code{integer}, or plain
    \code{numeric}). There may not be any missing values, \code{NaN}s,
    \code{NULL}s or \code{NA}s. 
}
  \item{contextual.name}{
A name of contextual variable to be weighted.
}
  \item{contextual.data}{
A \code{data.frame} containing the micro-level data to generate contextual 
    indicators by aggregation. May be \code{NULL}, in which case the individual data is
    used. The data must be numeric (\code{logical}, \code{integer}, or plain
    \code{numeric}). There may not be any missing values, \code{NaN}s,
    \code{NULL}s or \code{NA}s. By default set to \code{NULL}. 
}
  \item{context.id}{
The name of the context ID variable.
}
  \item{nb.resamples}{
A number of resamples to be evaluated. By default set to 1000.
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
  \item{individual.weight.name}{
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
  \item{confidence.intervals}{
    A \code{vector} of confidence intervals. Defaults  to \code{c(.95)}
    which corresponds to 95 \%.
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
  \item{contextual.sample.seed}{
    Seed used to generate the random sampling of the contextual data
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
    Defaults to \code{NULL}

}
  \item{kernel}{
 A function applied to the distance matrix. By default \code{NULL}, in which
    case the kernel function \deqn{w_{ij}=f(d,h) =
  \left(\frac{1}{2}\right)^{d_{ij}^2/h^2}}{w_ij = f(d, h) = 
  (1/2)^((d_ij/h)^2)} is used, where \eqn{w_{ij}, d_{ij}, h}{w_ij, d_ij,
  h}are elements of the weight matrix \eqn{\mathbf{W}}{W}, of the
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
\examples{#### ResampleExploreSpawML
# Exploratory spatially weighted multilevel analysis, with resampled
# both individual level indicators and contextual predictors for
# aggregation. Accepts only one contextual predictor. Predicting
# collective guilt acceptance. Contextual predictor weighted with
# geographical proximity weights, h=50,100,200. 

# load individual level data, remove collective guilt assignment for the
# data frame, and remove NA's 
data(traces_ind) 

traces_ind <- traces_ind[,-7]
traces_ind <- na.exclude(traces_ind)

# load contextual indicator for aggregation
data(traces_event) 

# load distance matrix 
data(d_geo) 

\donttest{acc_homog <- ResampleExploreSpawML(individual.level.data=traces_ind,
  contextual.name="w_all", contextual.data=traces_event,
  context.id="area", nb.resamples=5,formula=cg_acc ~ victim_d + comb_d +
  male + age_1990 + high_school + higher_edu + (1|area),
  distance.matrix=d_geo, multilevel.bandwidths=c(50,100,200)) }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

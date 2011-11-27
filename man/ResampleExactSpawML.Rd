\name{ResampleExactSpawML}
\alias{ResampleExactSpawML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Spatially weighted multilevel analysis with stratified resampling for exact context data
}
\description{
  Performs multilevel analysis with spatially weighted context data based on precise contextual measures. Provides robust point estimates and adjusted standard errors, obtained by stratified bootstrap resampling.
}
\usage{
ResampleExactSpawML(individual.level.data,
 contextual.names,
 context.id,
 formula,
 precise.data,
 contextual.weight.matrices = NULL,
 confidence.intervals = c(0.95),
 nb.resamples = 1000,
 individual.sample.seed = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{individual.level.data}{
    A \code{data.frame} containing the individual level data. The data must
    be numeric (\code{logical}, \code{integer}, or plain
    \code{numeric}). There may not be any missing values, \code{NaN}s,
    \code{NULL}s or \code{NA}s. 
  }
  \item{contextual.names}{
    A \code{list} of names of contextual variables to be weighted.
  }
  \item{context.id}{
    The name of the context ID variable.
  }
  \item{formula}{
    Formula description of the model.
  }
  \item{precise.data}{
    A \code{data.frame} containing the precise contextual data. May be
    \code{NULL}, but may not contain any missing values, like \code{NaN}s,
    \code{NULL}s or \code{NA}s. By default set to \code{NULL}, except for
    \code{ResampleExactSpawML} which requires precise data to be specified.
  }
  \item{contextual.weight.matrices}{
    A \code{list} of weights to be applied to each variable specified in
    \code{contextual.names}. A weight may be a weights matrix as, for
    instance computed by \code{\link{WeightMatrix}}, or \code{NULL}, in
    which case the corresponding contextual variable is not weighted. If
    only one weight is defined (instead of a list) it is applied to all
    contextual variables. Defaults to \code{NULL}, which means that none
    of contextual indicators are weighted.     
  }
  \item{confidence.intervals}{
    \code{vector} of confidence intervals. Defaults  to \code{c(.95)}
    which corresponds to 95 \%.
  }
  \item{nb.resamples}{
    number of resamples to be evaluated. By default set to 1000
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
# Spatially weighted multilevel analysis, with resampled individual
# level indicators and precise contextual indicator. Predicting collective
#guilt acceptance.


# load individual level data, remove collective guilt assignment for the
#data frame, and remove NA's 
data(traces_ind)

traces_ind <- traces_ind[,-7]
traces_ind <- na.exclude(traces_ind)


# load precise contextual indicator
data(homog_census)

# load distance matrix and create weights
data(d_geo)
geow_100 <- WeightMatrix(d_geo, bandwidth=100)

# perform ResampleExactSpawML
\donttest{acc_homog00_100 <- ResampleExactSpawML(individual.level.data=traces_ind,
  contextual.names="Homog_00", context.id="area", formula=cg_acc ~
  victim_d + comb_d + male + age_1990 + high_school + higher_edu +
  (1|area), precise.data=homog_census,
  contextual.weight.matrices=geow_100, nb.resamples=10) }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

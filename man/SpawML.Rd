\name{SpawML}
\alias{SpawML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Standard multilevel analysis applied to spatially weighted context data
}
\description{
  Integrates spatially weighted context data as higher-level predictors in multilevel analysis and produces conventional statistical estimates.
}
\usage{
SpawML(individual.level.data,
 contextual.names,
 context.id,
 formula,
 contextual.data = NULL,
 precise.data = NULL,
 contextual.weight.matrices = NULL,
 individual.weight.names = NULL,
 aggregation.functions = "mean")
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
  
  \item{contextual.data}{
    A \code{data.frame} containing the micro-level data to generate contextual 
    indicators by aggregation. May be \code{NULL}, in which case the individual data is
    used. The data must be numeric (\code{logical}, \code{integer}, or plain
    \code{numeric}). There may not be any missing values, \code{NaN}s,
    \code{NULL}s or \code{NA}s. By default set to \code{NULL}. 
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
  \item{individual.weight.names}{
    A \code{list} of optional design weights at the individual level
    used for aggregation (for example, for a weighted mean). List must
    have same length as contextual.names. May contain \code{NULL}s for
    variables which should not be weighted at the individual level. If
    only one individual weight is defined (instead of a list) it is
    applied to all contextual variables. By default set to \code{NULL}.
  }
  \item{aggregation.functions}{
    A \code{list} of aggregation functions. Functions take either
    \enumerate{
      \item 1 argument in which case the corresponding individual design
      weight is \code{NULL}, 
      \item 2 arguments in which case the second argument is taken
      from the corresponding individual design weight.
      Defaults to \code{"mean"}.
    }
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
# Spatially weighted multilevel analysis, with standard estimates of
# standard errors. Model with two contextual predictors, treated as
# precise, predicting collective guilt acceptance. 
# load individual level data, remove collective guilt assignment for the
# data frame, and remove NA's 
data(traces_ind)

traces_ind <- traces_ind[,-7]
traces_ind <- na.exclude(traces_ind)

# load contextual indicator for aggregation
data(traces_event)

# load precise contextual indicator
data(homog_census)

# load distance matrix and create weights
data(d_geo)
geow_100 <- WeightMatrix(d_geo, bandwidth=100)

# perform SpawML
\donttest{acc_w_homog_100 <- SpawML(individual.level.data=traces_ind,
  contextual.names=c("w_all", "Homog_00"), context.id="area",
  formula=cg_acc ~ victim_d + comb_d + male + age_1990 + high_school +
  higher_edu + (1|area), contextual.data=traces_event,
  precise.data=homog_census, contextual.weight.matrices=geow_100,
  individual.weight.names="weight", aggregation.functions="weighted.mean")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

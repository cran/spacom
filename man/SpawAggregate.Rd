\name{SpawAggregate}
\alias{SpawAggregate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create and describe spatially weighted survey estimates for contextual
  indicators
}
\description{

This step-2 function should be used to create spatially weighted context data
when contextual indicators are aggregate estimates derived from available
micro-level survey data. Aggregation functions can be user-defined or chosen
among predefined options, which include indicators of central tendency (e.g.,
spatially weighted mean), dispersion (e.g., spatially weighted standard
deviation) and inter-group variability (e.g., spatially weighted Gini
coefficient for group inequality). Point estimates can be adjusted by
user-defined design weights. Descriptive confidence intervals are computed by
way of an ad hoc bootstrap resampling procedure.

}
\usage{
SpawAggregate(contextual.data,
 context.id,
 contextual.names,
 contextual.weight.matrices=NULL,
 nb.resamples = 1000,
 aggregation.functions = "weighted.mean",
 confidence.intervals = 0.95,
 design.weight.names = NULL,
 sample.seed = NULL,
 additional.args = NULL,
 verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contextual.data}{
    A \code{data.frame} consisting of micro-level data to generate contextual 
    indicators by aggregation and containing a column named \code{context.id} 
    with the context ID variable. This is an individual level dataframe. 
    The data must be numeric (\code{logical}, \code{integer}, or plain
    \code{numeric}). There must not be any missing values, \code{NaN}s,
    \code{NULL}s or \code{NA}s. 
  }
  \item{context.id}{
      The name of the context ID variable. This variable allows matching
  contextual units from different data sets (\code{contextual.data},
  \code{contexutal.weight.matrix}). Values can by of
  any type but need to be defined consitently across data sets.
  }
  \item{contextual.names}{
    A \code{list} of names of contextual variables to be aggregated.
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
  \item{nb.resamples}{
    A number of resamples to be evaluated. By default set to 1000.
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
  \item{confidence.intervals}{
    A \code{vector} of confidence intervals. Defaults  to \code{0.95}
    which corresponds to 95 \%.
  }
  \item{design.weight.names}{
    a \code{list} of optional design weights at the individual level
    used for aggregation (for example, for a weighted mean). List must
    have same length as contextual.names. May contain \code{NULL}s for
    variables which should not be weighted at the individual level. If
    only one design.weight is defined (instead of a list) it is
    applied to all contextual variables. By default set to \code{NULL}.
  }
  \item{sample.seed}{
    Is one of three things
    \enumerate{
      \item \code{NULL}, in which case whatever the current random seed is
      is used
      \item an \code{integer}, which will be used to set the random
      seed. This allows reproducible random samples
      \item a saved \code{.Random.seed}  which allows reproducible random
      samples as well. The reason why both 2. and 3. are present is because
      \code{.Random.seed} can be saved a posteriori.
    }
    Defaults to \code{NULL}. 
  }
  \item{additional.args}{
   For aggregation functions which take additional arguments (that is in
  addition to the data to aggregate and design weights), they can be
  specified here.
  \code{additional.args} defaults to \code{NULL}
  }

\item{verbose}{
  if \code{TRUE}, progress information and an estimate of the remaining
  computing time is shown at the screen. Defaults to \code{TRUE}
  }
}
\details{
  \code{SpawAggregate} can be used for two similar, yet distinct,
  purposes:
  \enumerate{
    \item Aggregate contextual data and weight it spatially. For this,
    make sure to set \code{nb.resamples=0}
    \item To resample thecontextual data repeatedly (bootstrap) and
    generate a series of aggregated and spatially weighted
    \code{data.frames}.
    For this, make sure to use a large number of resamples.
  }
  The function \code{SpawAggregate} works exactly the same way in either
  case, but the output changes. See return value for this
  
}
\value{
  An object of class \code{SpawAggregateOutput-class} if
  \code{nb.resamples} is different than 0, a \code{data.frame} in the
  contrary case.
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  %%  \item{comp1 }{Description of 'comp1'}
  %%  \item{comp2 }{Description of 'comp2'}
  %% ...
}
%\references{
%  %% ~put references to the literature/web site here ~
%}
\author{
  Till Junge, Sandra Penic, Mathieu Cossuta, Guy Elcheroth
}
%\note{
  %%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
  %% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## Creating spatially weighted aggregated contextual indicator - spatially weighted 
## risk of war victimization across TRACES areas.
## This is step-2 function

## Step 1: Load distance matrix and create weights
data(d_geo)
geow.50 <- WeightMatrix(d_geo, bandwidth=50)

## Step 2: Load dataframe with micro-level data for aggregation and create 
## spatially weighted contextual indicator (risk of war victimization weighted by geow.50)

data(traces_event)
\dontshow{traces_event=traces_event[seq(1, nrow(traces_event), by=100),]}

wv.g50 <- SpawAggregate(contextual.data=traces_event,
                           context.id="area.name",
                           contextual.names="w_all",
                           contextual.weight.matrices=geow.50,
                           aggregation.functions="weighted.mean",
                           design.weight.names="weight",
                           nb.resamples=5)

## To create non-weighted and spatially weighted (by geow.50) risk of war victimization 

wv.nw.g50 <- SpawAggregate(contextual.data=traces_event,
                           context.id="area.name",
                           contextual.names=c("w_all","w_all"),
                           contextual.weight.matrices=list(NULL,geow.50),
                           aggregation.functions="weighted.mean",
                           design.weight.names="weight",
                           nb.resamples=5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ spatial }

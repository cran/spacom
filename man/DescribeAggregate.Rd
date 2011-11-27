\name{DescribeAggregate}
\alias{DescribeAggregate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Descriptives for aggregated contextual indicators
}
\description{
  Computes descriptives of spatially weighted aggregated contextual indicators.
}
\usage{
DescribeAggregate(contextual.data,
 context.id,
 contextual.names,
 contextual.weight.matrices,
 nb.resamples = 1000,
 aggregation.functions = "mean",
 confidence.intervals = 0.95,
 individual.weight.names = NULL,
 sample.seed = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contextual.data}{
    A \code{data.frame} consisting of micro-level data to generate contextual 
    indicators by aggregation and containing a column named \code{context.id} 
    with the context ID variable. This is an individual level dataframe. 
}
  \item{context.id}{
    A name of the context ID variable.
}
  \item{contextual.names}{
A \code{list} of names of contextual variables to be weighted.
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
    A \code{vector} of confidence intervals. Defaults  to \code{c(.95)}
  which corresponds to 95 \%.
}
  \item{individual.weight.names}{
    a \code{list} of optional design weights at the individual level
  used for aggregation (for example, for a weighted mean). List must
  have same length as contextual.names. May contain \code{NULL}s for
  variables which should not be weighted at the individual level. If
  only one individual weight is defined (instead of a list) it is
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
  samples as well. The reason why both b) and c) are present is because
  \code{.Random.seed} can be saved a posteriori.
}
Defaults to \code{NULL}. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An object of class \code{DescribeAggregateOutput-class}
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
# creating spatially weighted (with geographical proximity weights, h=50
# and h=200) contextual indicator of risk of war victimization through
# stratified resampling 

# load distance matrix and create weights
data(d_geo)
geow_50 <- WeightMatrix(d_geo, bandwidth=50)
geow_200 <- WeightMatrix(d_geo, bandwidth=200)

# load contextual indicator for aggregation
data(traces_event)

# perform DescribeAggregate
\donttest{wv_g50_200 <- DescribeAggregate(
   contextual.data=traces_event,
   context.id="area",
   contextual.names=c("w_all", "w_all"),
   contextual.weight.matrices=list(geow_50, geow_200),
   aggregation.functions="weighted.mean",
   individual.weight.names="weight",
   nb.resamples=5)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

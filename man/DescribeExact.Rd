\name{DescribeExact}
\alias{DescribeExact}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Descriptives for exact context data
}
\description{
  Computes descriptives of spatially weighted exact context data

}
\usage{
DescribeExact(contextual.data, context.id, contextual.names,
 contextual.weight.matrices)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contextual.data}{
    %%     ~~Describe \code{contextual.data} here~~
    A \code{data.frame}  consisting of precise contextual indicators and
    a context ID variable
}
  \item{context.id}{
    %%     ~~Describe \code{context.id} here~~
    A \code{character} object containing the column name of the context
    ID variable
}
  \item{contextual.names}{
    %%     ~~Describe \code{contextual.names} here~~
    A \code{list} of column names of contextual variables to be weighted
}
  \item{contextual.weight.matrices}{
    %%     ~~Describe \code{contextual.weight.matrices} here~~
    A \code{list} of square weight
    matrices of dimension \eqn{n\times n}{n by n}, where \eqn{n} is the number of
    contextual units to be applied to the variable specified in
    \code{contextual.names}, respectively.
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
A \code{data.frame} containing the spatially weighted \code{contextual.names}. 
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
# creating precise contextual indicator Homog_00 weighted by geow_50
# (geographical proximity weights, h=50) 

# load distance matrix and create weights
data(d_geo)
geow_50 <- WeightMatrix(d_geo, bandwidth=50)

# load precise contextual indicator
data(homog_census)

# create geographically weighted precise contextual indicator
\donttest{homog_g50 <- DescribeExact(contextual.data=homog_census,
                           context.id="area",
                           contextual.names="Homog_00",
                           contextual.weight.matrices=geow_50)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

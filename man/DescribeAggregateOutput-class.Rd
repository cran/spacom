\name{DescribeAggregateOutput-class}
\Rdversion{1.1}
\docType{class}
\alias{DescribeAggregateOutput-class}
\alias{show,DescribeAggregateOutput-method}

\title{Class DescribeAggregateOutput}
\description{
  %%  ~~ A concise (1-5 lines) description of what the class is. ~~
  Output object of the function \code{DescribeAggregate}. The user
  should never have to create an object of this class.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DescribeAggregateOutput", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{seed}:}{Object of class \code{"integer"}}

    contains the random seed which was used in \code{DescribeAggregate}.
    Use this seed if you want to reproduce the same resamples.
    
    \item{\code{aggregated.samples}:}{Object of class \code{"list"}}
    
    matrix of dimension \eqn{n\times m}{n by m} where \eqn{n} is the number of
    contextual units and \eqn{m} is the number of resamples containing
    the aggregated context data for each resample.
    
    \item{\code{frames}:}{Object of class \code{"list"}  }

    \code{list} of \code{data.frame}s. One \code{data.frame} per
    contextual name, each containing the mean estimate, standard
    deviation, median and confidence intervals of the corresponding
    contextual variable by contextual unit.
    
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "DescribeAggregateOutput")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
  Till Junge, Sandra Penic, Guy Elcheroth
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{DescribeAggregate}}
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("DescribeAggregateOutput")
}
\keyword{classes}

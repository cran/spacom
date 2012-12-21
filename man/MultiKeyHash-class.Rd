\name{MultiKeyHash-class}
\Rdversion{1.1}
\docType{class}
\alias{MultiKeyHash-class}
\alias{[<-,MultiKeyHash,ANY,ANY,ANY-method}
\alias{[,MultiKeyHash,ANY,ANY,ANY-method}
\alias{[[<-,MultiKeyHash,ANY,ANY-method}
\alias{[[,MultiKeyHash,ANY,ANY-method}
\alias{length,MultiKeyHash-method}

\title{Class \code{"MultiKeyHash"}}
\description{
  for internal use only
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("MultiKeyHash", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{indices}:}{Object of class \code{"list"}}
    \item{\code{data}:}{Object of class \code{"list"}}
  }
}
\section{Methods}{
  \describe{
    \item{[<-}{\code{signature(x = "MultiKeyHash", i = "ANY", j = "ANY", value = "ANY")}: ... }
    \item{[}{\code{signature(x = "MultiKeyHash", i = "ANY", j = "ANY", drop = "ANY")}: ... }
    \item{[[<-}{\code{signature(x = "MultiKeyHash", i = "ANY", j = "ANY")}: ... }
    \item{[[}{\code{signature(x = "MultiKeyHash", i = "ANY", j = "ANY")}: ... }
    \item{length}{\code{signature(x = "MultiKeyHash")}: ... }
	 }
}
\author{
  Till Junge
}
\keyword{classes}
\keyword{internal}

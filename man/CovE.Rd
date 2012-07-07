\name{CovE}

\alias{CovE}
\alias{ICovE}
\alias{CovE.Scanlda}
\alias{ICovE.Scanlda}
\alias{CovE.RFcanlda}
\alias{ICovE.RFcanlda}
\alias{CovE.Scllda}
\alias{ICovE.Scllda}
\alias{CovE.RFcllda}
\alias{ICovE.RFcllda}

\title{Generic methods for extracting covariance and inverse covariance matrices from objects storing the results of a Linear Discriminant Analysis}

\description{

\sQuote{CovE} Extracts an object with an appropriate representation of a whithin groups covariance matrix from a \sQuote{Scanlda}, \sQuote{Scllda}, \sQuote{RFcanlda} or a \sQuote{RFcllda} object.

\sQuote{ICovE} Extracts an object with an appropriate representation of a whithin groups inverse covariance matrix from a \sQuote{Scanlda}, \sQuote{Scllda}, \sQuote{RFcanlda} or a \sQuote{RFcllda} object.

}

\usage{
\method{CovE}{Scanlda}(object)
\method{ICovE}{Scanlda}(object)
\method{CovE}{RFcanlda}(object)
\method{ICovE}{RFcanlda}(object)
\method{CovE}{Scllda}(object)
\method{ICovE}{Scllda}(object)
\method{CovE}{RFcllda}(object)
\method{ICovE}{RFcllda}(object)
}

\arguments{

  \item{object}{An object of class \sQuote{Scanlda}, \sQuote{Scllda}, \sQuote{RFcanlda} or \sQuote{RFcllda}.}

}

\value{An object with an appropriate representation of the matrix extracted.}

\seealso{\code{\link{Slda}}, \code{\link{RFlda}}, \code{\link{canldaRes}}, \code{\link{clldaRes}} }



\name{clldaRes}

\alias{clldaRes}
\alias{predict.clldaRes}
\alias{print.clldaRes}
\alias{coef.clldaRes}
\alias{CovE.clldaRes}
\alias{ICovE.clldaRes}

\title{Class object used for storing the results of a  high-dimensional linear discriminant analysis routine 
	(with \sQuote{ldafun} argument set to \dQuote{classification}).}

\description{\sQuote{predict.clldaRes} Classifies multivariate observations in conjunction with a \sQuote{clldaRes} object.

             \sQuote{print.clldaRes} is the S3 print method for \sQuote{clldaRes} objects.

             \sQuote{coef.clldaRes} is the S3 coef method for \sQuote{clldaRes} objects. }

\usage{
\method{predict}{clldaRes}(object, newdata, prior=object$prior, grpcodes=NULL,\dots)
}

\arguments{

  \item{object}{An object of class \sQuote{clldaRes}.}

  \item{newdata}{Matrix of cases to be classified.}

  \item{prior}{The prior probabilities used.} 

  \item{grpcodes}{Factor with the class codes. Set to 0:k-1 (k being the number of different classes) by default.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
     A list with components

\item{class}{The MAP classification (a factor)}

\item{x}{The classification scores of the test cases.}
}

\seealso{\code{\link{Dlda}}, \code{\link{Mlda}}, \code{\link{Slda}}, \code{\link{RFlda}}, \code{\link{print}}, \code{\link{coef}} }



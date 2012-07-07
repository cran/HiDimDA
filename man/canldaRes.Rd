\name{canldaRes}

\alias{canldaRes}
\alias{predict.canldaRes}
\alias{print.canldaRes}
\alias{coef.canldaRes}
\alias{CovE.canldaRes}
\alias{ICovE.canldaRes}

\title{Class object used for storing the results of a canonical high-dimensional linear discriminant analysis.}

\description{\sQuote{predict.canldaRes} classifies multivariate observations in conjunction with a \sQuote{canldaRes} object.

             \sQuote{print.canldaRes} is the S3 print method for \sQuote{canldaRes} objects.

             \sQuote{coef.canldaRes} is the S3 coef method for \sQuote{canldaRes} objects.
}

\usage{
\method{predict}{canldaRes}(object, newdata, prior=object$prior, grpcodes=NULL, 
nbvrs=ncol(object$scaling), \dots)
}

\arguments{

  \item{object}{An object of class \sQuote{canldaRes}.}

  \item{newdata}{Matrix of cases to be classified.}

  \item{prior}{The prior probabilities used.} 

  \item{grpcodes}{Factor with the class codes. Set to 0:k-1 (k being the number of different classes) by default.}

  \item{nbvrs}{Number of canonical discriminant variables used for prediction.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
     A list with components

\item{class}{The MAP classification (a factor)}

\item{ZsqDistances}{A matrix with the squared Euclidean distance, in the discriminant space, of each new observation to the group centroids.}

\item{prior}{The prior probabilities used.} 

\item{ZsqDprioradj}{The adjustments to squared Euclidean distance, induced by the chosen (or estimated) prior probabilities.}

\item{Z}{A matrix with the values of the canonical discriminant variates.}

\item{Zmeans}{A matrix with the group centroids in the canonical space.}
}

\seealso{ \code{\link{Dlda}}, \code{\link{Mlda}}, \code{\link{Slda}}, \code{\link{RFlda}},  \code{\link{print}}, \code{\link{coef}} }



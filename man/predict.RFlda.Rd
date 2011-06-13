\name{predict.RFlda}

\alias{predict.RFlda}

\title{Classify Multivariate Observations by Linear Discrimination based on Factor models.}

\description{Classifies multivariate observations in conjunction with \sQuote{RFlda}.}

\usage{
\method{predict}{RFlda}(object, newdata, prior=object$prior, grpcodes=NULL, \dots)
}

\arguments{

  \item{object}{An object of class \sQuote{RFlda}.}

  \item{newdata}{Matrix of cases to be classified.}

  \item{prior}{The prior probabilities used.} 

  \item{grpcodes}{Factor with the class codes. Set to 0:k-1 (k being the number of different classes) by default.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{
     A list with components

\item{class}{The MAP classification (a factor)}

\item{posterior}{The posterior probabilities for the classes}

\item{x}{The scores of test cases}
}

\seealso{\code{\link{RFlda}}}



\name{SigFq}

\alias{SigFq}
\alias{is.SigFq}
\alias{as.matrix.SigFq}
\alias{print.SigFq}

\title{SigFq objects: covariance matrices associated with a q-factor model}

\description{Creates a \sQuote{SigFq} object}

\usage{
SigFq(D, B, p, q)
}

\arguments{
  \item{D}{A p-dimensional vector with the specific variances of the assumed factor model.}

  \item{B}{A p by q matrix with the loadings of assumed factor model.}

  \item{p}{The full dimensionality (i.e., the number of original variables) of the process being modelled.}

  \item{q}{The number of factors in the assumed factor model.} 
}

\value{ An object of class \sQuote{SigFq} for which the generic methods \sQuote{as.matrix} (converting to a traditional numeric matrix), \sQuote{is}, as well as specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations, are available.}

\seealso{\code{\link{ForbSigap}}, \code{\link{solve.SigFq}}, \code{\link{LeftMult.SigFq}}, \code{\link{RightMult.SigFq}}, \code{\link{SigFqInv}}}


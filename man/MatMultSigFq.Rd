\name{MatMultSigFq}

\alias{LeftMult}
\alias{RightMult}
\alias{LeftMult.SigFq}
\alias{RightMult.SigFq}
\alias{LeftMult.SigFqInv}
\alias{RightMult.SigFqInv}

\title{MatMultSigFq: Specialized matrix multiplication of \sQuote{SigFq} and \sQuote{SigFqInv} objects}
\description{
     \sQuote{LeftMult} multiplies, on the left, a vector or matrix of compatible dimensions, by a \sQuote{SigFq} or \sQuote{SigFqInv} object.

     \sQuote{RightMult} multiplies, on the right, a vector or matrix of compatible dimensions, by a \sQuote{SigFq} or \sQuote{SigFqInv} object.
}

\usage{
\method{LeftMult}{SigFq}(x, a)
\method{RightMult}{SigFq}(x, a)
\method{LeftMult}{SigFqInv}(x, a)
\method{RightMult}{SigFqInv}(x, a)
}

\arguments{
  \item{x}{An object of class \sQuote{SigFq} or class \sQuote{SigFqInv}, representing the covariance or precision matrix associated with a q-factor model.} 

  \item{a}{A vector, or matrix, by which \sQuote{x} is to be multiplied.}
}

\value{A vector or a (traditional numeric) matrix with the result of the matrix product.}

\seealso{\code{\link{FrobSigAp}}, \code{\link{SigFq}}, \code{\link{SigFqInv}}, \code{\link{solve.SigFq}}, \code{\link{solve.SigFqInv}}}


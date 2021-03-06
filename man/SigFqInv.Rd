\name{SigFqInv}

\alias{SigFqInv}
\alias{is.SigFqInv}
\alias{as.matrix.SigFqInv}
\alias{print.SigFqInv}
\alias{+.SigFqInv}
\alias{-.SigFqInv}
\alias{*.SigFqInv}
\alias{/.SigFqInv}

\title{SigFqInv objects: precision (inverse of covariance) matrices associated with a q-factor model}

\description{
Creates a \sQuote{SigFqInv} object.
}

\usage{
SigFqInv(D, B, p, q, optres=NULL)
}

\arguments{
  \item{D}{A p-dimensional vector with the specific variances of the assumed factor model.}

  \item{B}{A p by q matrix with the loadings of assumed factor model.}

  \item{p}{The full dimensionality (i.e., the number of original variables) of the process being modelled.}

  \item{q}{The number of factors in the assumed factor model.} 

  \item{optres}{The full optimization results provided by the error minimization algorithm used to find SigFq parameters.}
}

\value{ An object of class \sQuote{SigFqInv} for which the generic method \sQuote{as.matrix} (converting to a traditional numeric matrix), as well as specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations, are available.}

\seealso{\code{\link{FrobSigAp}}, \code{\link{solve.SigFqInv}}, \code{\link{LeftMult.SigFqInv}}, \code{\link{RightMult.SigFqInv}}, \code{\link{SigFq}}}

 

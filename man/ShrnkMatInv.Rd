\name{ShrnkMatInv}

\alias{ShrnkMatInv}
\alias{is.ShrnkMatInv}
\alias{as.matrix.ShrnkMatInv}
\alias{print.ShrnkMatInv}
\alias{+.ShrnkMatInv}
\alias{-.ShrnkMatInv}
\alias{*.ShrnkMatInv}
\alias{/.ShrnkMatInv}

\title{ShrnkMatInv objects: precision (inverse of covariance) matrices associated with shrunken estimates of a covariance}

\description{
Creates a \sQuote{ShrnkMatInv} object.
}

\usage{
ShrnkMatInv(U, D, p, q, Intst, Trgt="Idntty")
}

\arguments{
  \item{U}{A p by q matrix with the orthonomal eigenvectors of the original (unshrunken) covariance estimate.}

  \item{D}{A p-dimensional vector with the eigenvalues of the original (unshrunken) covariance estimate.}

  \item{p}{The dimension of the covariance matrix.}

  \item{q}{The rank of the original (unshrunken) covariance estimate.} 

  \item{Intst}{The target intensity used by the shunk estimator.}

  \item{Trgt}{The target used by the shunk estimator. If set to the string \dQuote{Idntty} (default)
a p-dimensional matrix identity target will be assumed. Otherwise a matrix-type object representing a symmetric matrix target.}
}

\value{ An object of class \sQuote{ShrnkMatInv} for which the generic method \sQuote{as.matrix} (converting to a traditional numeric matrix), as well as specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations, are available.}

\seealso{

\code{\link{ShrnkMat}}, \code{\link{solve.ShrnkMatInv}}, \code{\link{LeftMult.ShrnkMatInv}}, 

\code{\link{RightMult.ShrnkMatInv}}

}

 

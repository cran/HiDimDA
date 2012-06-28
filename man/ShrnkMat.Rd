\name{ShrnkMat}

\alias{ShrnkMat}
\alias{is.ShrnkMat}
\alias{as.matrix.ShrnkMat}
\alias{print.ShrnkMat}
\alias{+.ShrnkMat}
\alias{-.ShrnkMat}
\alias{*.ShrnkMat}
\alias{/.ShrnkMat}

\title{ShrnkMat objects: shrunken matrix estimates of a covariance}

\description{
Creates a \sQuote{ShrnkMat} object.
}

\usage{
ShrnkMat(U, D, p, q, Intst, Trgt="Idntty")
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

\value{ An object of class \sQuote{ShrnkMat} for which the generic method \sQuote{as.matrix} (converting to a traditional numeric matrix), as well as specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations, are available.}

\seealso{ \code{\link{solve.ShrnkMat}}, \code{\link{LeftMult.ShrnkMat}}, \code{\link{RightMult.ShrnkMat}}, \code{\link{ShrnkMatInv}} }

 

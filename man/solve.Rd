\name{solve}

\alias{solve.DMat}
\alias{solve.ShrnkMat}
\alias{solve.ShrnkMatInv}
\alias{solve.SigFq}
\alias{solve.SigFqInv}

\title{Solve methods for \sQuote{DMat}, \sQuote{ShrnkMat}, \sQuote{ShrnkMatInv}, \sQuote{SigFq} and \sQuote{SigFqInv} objects.}

\description{
     These methods solve the equation \code{a \%*\% x = b} for \sQuote{x},
     where \sQuote{a} is a specialized square symmetric matrix represented by a 
     \sQuote{DMat}, \sQuote{ShrnkMat}, \sQuote{ShrnkMatInv}, \sQuote{SigFq} 
      or \sQuote{SigFqInv} object, and \sQuote{b} can be either a vector or a matrix.
}

\usage{
\method{solve}{DMat}(a, b = NULL, \dots)
\method{solve}{ShrnkMat}(a, b = NULL, \dots)
\method{solve}{ShrnkMatInv}(a, b = NULL, \dots)
\method{solve}{SigFq}(a, b = NULL, \dots)
\method{solve}{SigFqInv}(a, b = NULL, \dots)
}

\arguments{
  \item{a}{An object of type \sQuote{DMat}, \sQuote{ShrnkMat}, \sQuote{ShrnkMatInv}, \sQuote{SigFq} or \sQuote{SigFqInv} representing a specialized symmetric square matrix.}
  \item{b}{A numeric vector or matrix giving the right-hand side(s) of the linear system. If missing, \sQuote{b} is taken to be an identity matrix and solve will return an object representing the inverse of the matrix associated with \sQuote{a}.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{The result returned depends on the values of the arguments. When \sQuote{b} is not NULL, both functions return a numeric vector or matrix with the solution of the system. When \sQuote{b} is NULL, \sQuote{solve.DMat} returns a \sQuote{DMat} object, \sQuote{solve.ShrnkMat} a \sQuote{ShrnkMatInv} object, \sQuote{solve.ShrnkMatInv} a \sQuote{ShrnkMat} object, \sQuote{solve.SigFq} a \sQuote{SigFqInv} object and \sQuote{solve.SigFqInv} returns a \sQuote{SigFq} object.}  

\seealso{\code{\link{DMat}}, \code{\link{ShrnkMat}}, \code{\link{ShrnkMatInv}}, \code{\link{SigFq}}, \code{\link{SigFqInv}}, \code{\link{FrobSigAp}}, \code{\link{solve}}}


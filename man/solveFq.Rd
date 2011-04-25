\name{solveFq}

\alias{solve.SigFq}
\alias{solve.SigFqInv}

\title{Solve methods for \sQuote{SigFq} and \sQuote{SigFqInv} objects}

\description{
     These methods solve the equation \code{a \%*\% x = b} for \sQuote{x},
     where \sQuote{a} is a covariance (or precision, i.e., the inverse of a covariance) 
     matrix represented by a \sQuote{SigFq} (or \sQuote{SigFqInv}) object, and \sQuote{b} 
     can be either a vector or a matrix.
}

\usage{
solve.SigFq(a, b = NULL, \dots)
solve.SigFqInv(a, b = NULL, \dots)
}

\arguments{
  \item{a}{An object of type \sQuote{SigFq} or \sQuote{SigFqInv} representing the covariance, or precision, matrix associated with a q-factor model.}
  \item{b}{A numeric vector or matrix giving the right-hand side(s) of the linear system. If missing, \sQuote{b} is taken to be an identity matrix 
and solve will return an object representing the inverse of the matrix associated with \sQuote{a}.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{The result returned depends on the values of the arguments. When \sQuote{b} is not NULL, both functions return a numeric vector or matrix with
the solution of the system. When \sQuote{b} is NULL, \sQuote{solve.SigFq} returns a \sQuote{SigFqInv} object and \sQuote{solve.SigFqInv} returns a \sQuote{SigFq} object.}  

\seealso{\code{\link{SigFq}}, \code{\link{SigFqInv}}, \code{\link{ForbSigap}}, \code{\link{solve}}}


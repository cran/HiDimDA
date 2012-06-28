\name{MatMult}

\alias{LeftMult}
\alias{RightMult}
\alias{LeftMult.matrix}
\alias{RightMult.matrix}
\alias{LeftMult.DMat}
\alias{RightMult.DMat}
\alias{LeftMult.ShrnkMat}
\alias{RightMult.ShrnkMat}
\alias{LeftMult.ShrnkMatInv}
\alias{RightMult.ShrnkMatInv}
\alias{LeftMult.SigFq}
\alias{RightMult.SigFq}
\alias{LeftMult.SigFqInv}
\alias{RightMult.SigFqInv}

\title{MatMult: Specialized matrix multiplication of \sQuote{DMat}, \sQuote{ShrnkMat}, \sQuote{ShrnkMatInv}, \sQuote{SigFq} and \sQuote{SigFqInv} objects.}
\description{
     \sQuote{LeftMult} multiplies, on the left, a vector or matrix of compatible dimensions, by a  \sQuote{DMat}, \sQuote{ShrnkMat}, \sQuote{ShrnkMatInv}, \sQuote{SigFq} or \sQuote{SigFqInv} object.

     \sQuote{RightMult} multiplies, on the right, a vector or matrix of compatible dimensions, by a  \sQuote{DMat}, \sQuote{ShrnkMat}, \sQuote{ShrnkMatInv}, \sQuote{SigFq} or \sQuote{SigFqInv} object.
}

\usage{
\method{LeftMult}{DMat}(x, a)
\method{RightMult}{DMat}(x, a)
\method{LeftMult}{ShrnkMat}(x, a)
\method{RightMult}{ShrnkMat}(x, a)
\method{LeftMult}{ShrnkMatInv}(x, a)
\method{RightMult}{ShrnkMatInv}(x, a)
\method{LeftMult}{SigFq}(x, a)
\method{RightMult}{SigFq}(x, a)
\method{LeftMult}{SigFqInv}(x, a)
\method{RightMult}{SigFqInv}(x, a)
}

\arguments{
  \item{x}{An object of class \sQuote{DMat}, \sQuote{ShrnkMat}, \sQuote{ShrnkMatInv}, \sQuote{SigFq} or class \sQuote{SigFqInv}, with a compact representation of a specialized square symmetric matrix.} 

  \item{a}{A vector, or matrix, by which \sQuote{x} is to be multiplied.}
}

\value{A vector or a (traditional numeric) matrix with the result of the matrix product.}

\seealso{\code{\link{DMat}}, \code{\link{FrobSigAp}}, \code{\link{ShrnkMat}}, \code{\link{ShrnkMatInv}}, \code{\link{SigFq}}, \code{\link{SigFqInv}}, \code{\link{solve.DMat}}, \code{\link{solve.ShrnkMat}}, \code{\link{solve.ShrnkMatInv}}, \code{\link{solve.SigFq}}, \code{\link{solve.SigFqInv}}}


\name{DMat}

\alias{DMat}
\alias{is.DMat}
\alias{as.matrix.DMat}
\alias{print.DMat}
\alias{+.DMat}
\alias{-.DMat}
\alias{*.DMat}
\alias{/.DMat}

\title{DMat objects: diagonal matrices}

\description{
Creates a \sQuote{DMat} object.
}

\usage{
DMat(D)
}

\arguments{
  \item{D}{A vector with the diagonal elements of the matrix.}
}

\value{ An object of class \sQuote{DMat} for which the generic method \sQuote{as.matrix} (converting to a traditional numeric matrix), as well as specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations, are available.}

\seealso{\code{\link{solve.DMat}}, \code{\link{LeftMult.DMat}}, \code{\link{RightMult.DMat}}}


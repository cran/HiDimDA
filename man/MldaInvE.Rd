\name{MldaInvE}
\alias{MldaInvE}

\title{Maximum uncertainty Linear Discriminant Analysis inverse matrix estimator.}

\description{
Builds a well-conditioned estimator for the inverse of a symmetric positive definite matrix, with a 
bad-conditioned or singular estimate, based on the \dQuote{Maximum Uncertainty Linear Discriminant Analysis} 
approach of Thomaz, Kitani and Gillies (2006).
}

\usage{  
MldaInvE(M, check=TRUE, onlyMinv=TRUE, numtol=sqrt(Machine$double.eps))  
}

\arguments{
  \item{M}{Singular or bad-conditioned estimate of the matrix for which a well-conditioned inverse estimate is sought.}

  \item{check}{Boolean flag indicating if the symmetry of M and the sign of its eigenvalues should be check upfront.}

  \item{onlyMinv}{Boolean flag indicating if only an estimate of the matrix inverse is sought, or if a well-conditioned 
approximation to the matrix that M estimates should be returned as well.}

  \item{numtol}{Numerical tolerance.}
}

\value{If onlyMinv is set to true a matrix with the inverse estimate sought. Otherwise a list with components ME and MInvE,
with a well-conditioned approximation to the matrix that M estimates and its inverse.}

\references{Thomaz, Kitani and Gillies (2006) \dQuote{A maximum uncertainty LDA-based approach for limited sample size problems - with application to face recognition}, \emph{Journal of the Brazilian Computer Society}, 12 (2), 7-18}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{Mlda}}}

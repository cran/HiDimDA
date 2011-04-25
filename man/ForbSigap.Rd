\name{ForbSigap}

\title{Approximation of Covariance Matrices from q-factor models}

\alias{ForbSigap}
\alias{ForbSigap1}

\description{
 \sQuote{ForbSigap} finds the parameters of a q-factor model whose covariance is closest to the matrix argument \sQuote{Sigma}, according to the Forbenius norm.

 \sQuote{ForbSigap1} finds the parameters of a q-factor model whose covariance is closest to the matrix square \code{SigmaSr \%*\% t(SigmaSr)} of the argument \sQuote{SigmaSr}, according to the Forbenius norm. 
}

\details{The minimization of the error Forbenius norm is performed by the \sQuote{nlminb} PORT optimization routine. The actual computation of
the errors, and their analytical gradients and hessians, are coded in C in order to speed up the algorithm.\cr The minimization procedure takes the loadings (B) of the factor model as arguments (on optimality, the specific variances (D) are always given by \code{D = diag(Sigma) - diag(B \%*\% t(B))},
and uses the first q eigenvectors of \sQuote{Sigma} as starting points of the optimization.\cr
For small values of q (1 or 2), this procedure seems to converge  relatively fast (usually less than one minute when \eqn{q=1} and \eqn{ncol(Sigma) \leq 1000}) to the global minimum of the approximation error. 
For larger values of q, the computational time can be much higher and multiple random starting points (use the argument \sQuote{nstarts}) may be required 
in order to escape local optima.
}

\usage{
ForbSigap(Sigma, q, nstarts = 1, k0 = NULL, penF = NULL)
ForbSigap1(SigmaSr, SigmaRank, q, nstarts = 1, k0=NULL, penF=NULL)
}

\arguments{
  \item{Sigma}{Square, symmetric and positive-definite matrix to be approximated}
  \item{q}{Number of factors in the assumed factor model.}
  \item{nstarts}{Number of different randomly generated starting points used in the optimization.}
  \item{k0}{Lower bound for the elements of the specific variances. When NULL (default), k0 is set to 0.01 times the minimum diagonal element of \sQuote{Sigma}.}
  \item{penF}{Penalty factor, used to forbid specific variances below the k0 bound. When NULL (default), is set to 100 time the maximum diagonal element of \sQuote{Sigma}.}
  \item{SigmaSr}{Matrix square root of the covariance to be approximated.}
  \item{SigmaRank}{Rank of the covariance matrix to be approximated.}
}

\value{ An object of class \sQuote{SigFq} representing the covariance assumed by the closest q-factor model.

\sQuote{SigFq} objects have specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations.}

\references{Duarte Silva, A.P. (Forthcoming) Two Group Classification with High-Dimensional Correlated Data: A Factor Model Approach. 
\emph{Computational Statistics and Data Analysis.}}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{RFlda}}, \code{\link{SigFq}}, \code{\link{SigFqInv}}, \code{\link{solve.SigFq}}, \code{\link{solve.SigFqInv}}}


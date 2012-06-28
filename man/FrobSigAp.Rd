\name{FrobSigAp}

\title{Approximation of Covariance Matrices from q-factor models}

\alias{FrobSigAp}
\alias{FrobSigAp1}

\description{
 \sQuote{FrobSigAp} finds the parameters of a q-factor model whose covariance is closest to the matrix argument \sQuote{Sigma}, according to the Frobenius norm.

 \sQuote{FrobSigAp1} finds the parameters of a q-factor model whose covariance is closest to the matrix square \code{SigmaSr \%*\% t(SigmaSr)} of the argument \sQuote{SigmaSr}, according to the Frobenius norm. 
}

\details{The minimization of the error Frobenius norm is performed by the \sQuote{nlminb} PORT optimization routine. The actual computations of the errors, and their analytical gradients and hessians, are coded in C in order to speed up the algorithm.\cr The minimization procedure takes the loadings (B) of the factor model as arguments, computes the optimal specific variances by the analytical formula \code{D = diag(Sigma) - diag(B \%*\% t(B))}, and uses the first q eigenvectors of \sQuote{Sigma} as starting points of the optimization.\cr
For small values of q (1 or 2), this procedure seems to quickly converge to the global minimum of the approximation error. 
For larger values of q, the computational time can be much higher and multiple random starting points (that can be specified by the argument \sQuote{nstarts}) may be required in order to escape local optima.
}

\usage{

FrobSigAp(Sigma, q, nstarts = 1, k0 = NULL, penF = NULL, 
atol = 1e-20, rtol = sqrt(.Machine$double.eps))

FrobSigAp1(SigmaSr, SigmaRank, q, nstarts = 1, k0=NULL, penF=NULL, 
atol = 1e-20, rtol = 100 * sqrt(.Machine$double.eps))
}

\arguments{

  \item{Sigma}{Square, symmetric and positive-definite matrix to be approximated}

  \item{q}{Number of factors in the assumed factor model.}

  \item{nstarts}{Number of different randomly generated starting points used in the optimization.}

  \item{k0}{Lower bound for the elements of the specific variances. When NULL (default), k0 is set to 0.01 times the minimum diagonal element of \sQuote{Sigma}.}

  \item{penF}{Penalty factor, used to forbid specific variances below the k0 bound. When set to NULL (default), a penalty equal to 100 times the maximum diagonal element of \sQuote{Sigma} is used.}

  \item{atol}{The absolute convergence tolerance of the local optimizer.}

  \item{rtol}{The relative convergence tolerance of the local optimizer. The local optimizer stops if it is unable to reduce the approximation error (err) by a factor of \sQuote{reltol *(abs(err) + reltol)} at a step.}

  \item{SigmaSr}{Matrix square root of the covariance to be approximated.}

  \item{SigmaRank}{Rank of the covariance matrix to be approximated.}
}

\value{ An object of class \sQuote{SigFq} representing the covariance assumed by the closest q-factor model.
\sQuote{SigFq} objects have specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations.}

\references{Pedro Duarte Silva, A. (2006) \dQuote{Two Group Classification with High-Dimensional Correlated Data: A Factor Model Approach}, 
\emph{Computational Statistics and Data Analysis}, 55 (1), 2975-2990.}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{RFlda}}, \code{\link{SigFq}}, \code{\link{SigFqInv}}, \code{\link{solve.SigFq}}, \code{\link{solve.SigFqInv}}}


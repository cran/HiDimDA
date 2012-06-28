\name{ShrnkSigE}
\alias{ShrnkSigE}

\title{Shrunken Covariance Estimate.}

\description{
Builds a well-conditioned shrunken estimate of a covariance matrix based on Fisher and Sun's (2011) estimates
and generalizations of Ledoit and Wolf's (2004) optimal optimal shrunken covariance matrices.
}

\usage{  
ShrnkSigE(df, p , SigmaRank, Sigma=NULL, SigmaSr=NULL, check=TRUE, Trgt,
minp=20, numtol=sqrt(.Machine$double.eps), \dots)  
}

\arguments{

  \item{df}{Degrees of freedom of the original (unshrunken) covariance estimate.}

  \item{p}{Dimension of the covariance matrix.}

  \item{SigmaRank}{Rank of the original (unshrunken) covariance estimate.}

  \item{Sigma}{Original (unshrunken) covariance estimate.}

  \item{SigmaSr}{Matrix square-root of the original (unshrunken) covariance estimate, i.e. a matrix, SigmaSr, such that SigmaSr^T SigmaSr  = Sigma, where Sigma is the original unshrunken covariance estimate.}

  \item{Trgt}{A string code with the target type used by the shrunken estimator. The alternatives are \dQuote{CnstDiag} for a Ledoit-Wolf constant diagonal target, \dQuote{Idntty} for a p-dimensional identity, and \dQuote{VarDiag} for a diagonal target of empirical variances.} 

  \item{check}{Boolean flag indicating if the symmetry and the sign of the Sigma eigenvalues should be check upfront.}

  \item{minp}{Minimum number of variables required for the estimation of the target intensity to be considered reliable. If the dimension of Sigma is below pmin, no shrunken estimate is computed and the original sample covariance is returned.}

  \item{numtol}{Numerical tolerance. All computed eigenvalues below numtol are considered equal to zero, and the rank of original shrunken estimate is adjusted acordingly.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{ShrnkSigE can take as input an original unshrunken estimate of the covariance matrix or, in alternative, one
matrix square-root, SigmaSr (e.g. the original, centered and scaled, data matrix), such that 
\eqn{SigmaSr^T SigmaSr = Sigma}. 
In problems with more variables than observations it is preferable to use a matrix square-root for reasons of memory and 
computational efficiency.}

\value{An object of class \sQuote{ShrnkMat} with a compact representation of the shrunken covariance estimator.}

\references{

Ledoit, O. and Wolf, M. (2004) \dQuote{A well-conditioned estimator for large-dimensional covariance matrices.}, \emph{Journal of Multivariate Analysis}, 88 (2), 365-411. 

Fisher, T.J. and Sun, X. (2011) \dQuote{Improved Stein-type shrinkage estimators for the high-dimensional multivariate normal covariance matrix}, \emph{Computational Statistics and Data Analysis}, 55 (1), 1909-1918. 

}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{ShrnkMat}}, \code{\link{Slda}}}

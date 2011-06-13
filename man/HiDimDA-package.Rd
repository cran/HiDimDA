\name{HiDimDA-package}

\alias{HiDimDA-package}
\alias{HiDimDA}

\docType{package}

\title{High Dimensional Discriminant Analysis}

\description{Performs Linear Discriminant Analysis in High Dimensional problems based on covariance estimators derived from low dimensional factor models. Includes routines for classifier training, prediction, cross-validation and variable selection.}

\details{
\tabular{ll}{
Package: \tab HiDimDA\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2011-06-13\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
}
}

\author{Antonio Pedro Duarte Silva  <psilva@porto.ucp.pt>

Maintainer: Antonio Pedro Duarte Silva  <psilva@porto.ucp.pt>
}

\seealso{\code{\link{RFlda}}, \code{\link{predict.RFlda}}, \code{\link{SinghDS}}}

\references{Pedro Duarte Silva, A. (2011) \dQuote{Two Group Classification with High-Dimensional Correlated Data: A Factor Model Approach}, 
\emph{Computational Statistics and Data Analysis}, doi:10.1016/j.csda.2011.05.002.}

\keyword{HiDimDA}

\examples{

#train classifier on Singh's Prostate Cancer Data set, 
# selecting genes by the Extended HC scheme 

ldarule <- RFlda(SinghDS[,-1],SinghDS[,1])     

# get in-sample classification results

predict(ldarule,SinghDS[,-1],grpcodes=levels(SinghDS[,1]))$class         	       

# compare classifications with true assignments

cat("Original classes:\n")
print(SinghDS[,1])             		 

# show set of selected genes

cat("Genes kept in discrimination rule:\n")
print(colnames(SinghDS)[ldarule$vkpt])             		 
cat("Number of selected genes =",ldarule$nvkpt,"\n")
}


\name{Slda}

\alias{Slda}
\alias{Slda.default}
\alias{Slda.data.frame}
\alias{is.Slda}

\title{Shrunken Linear Discriminant Analysis.}

\description{
 \sQuote{Slda} finds the coefficients of a linear discriminant rule based on Fisher and Sun's (2011) estimate and generalizations of Ledoit and Wolf's (2004) optimal shrunken covariance matrix.
}

\usage{

\method{Slda}{default}(data, grouping, prior = "proportions", StddzData=TRUE, VSelfunct = SelectV,
Trgt=c("CnstDiag","Idntty","VarDiag"), minp=20, ldafun=c("canonical","classification"), \dots)

\method{Slda}{data.frame}(data, \dots)
}

\arguments{
  \item{data}{Matrix or data frame of observations.}

  \item{grouping}{Factor specifying the class for each observation.}

  \item{prior}{The prior probabilities of class membership.  If unspecified, the class proportions for the training set are used.  If
present, the probabilities should be specified in the order of the factor levels.}

  \item{StddzData}{A boolean flag indicating whether the data should be standardized first (default) or used in their original scales.}

  \item{VSelfunct}{Variable selection function. Either the string \dQuote{none} (no selection is to be performed) or a function that takes \sQuote{data} and \sQuote{grouping} as its first two arguments and returns a list with two components: (i) \sQuote{nvkpt} - the number of variables to be used in the Discriminant rule; and (ii) \sQuote{vkptInd} - the indices of the variables to be used in the Discriminant rule. The default is the \sQuote{SelectV} function that, by default, selects variables by the Expanded HC scheme described in Duarte Silva (2011).} 

  \item{Trgt}{A string code with the target type used by the shrunken estimator. The alternatives are \dQuote{CnstDiag} for a Ledoit-Wolf constant diagonal target, \dQuote{Idntty} for a p-dimensional identity, and \dQuote{VarDiag} for a diagonal target of empirical variances.} 

  \item{minp}{Minimum number of variables required for the estimation of the target intensity to be considered reliable. If the dimension of Sigma is below pmin, no shrunken estimate is computed and the original sample covariance is employed.}

  \item{ldafun}{Type of discriminant linear functions computed. The alternatives are \dQuote{canonical} for maximum-discrimination canonical
functions and  \dQuote{classification} for direct-classification functions.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{

If algument \sQuote{ldafun} is set to \dQuote{canonical} an object of class \sQuote{Scanlda}, which extends class \sQuote{canldaRes}, with the following components:

\item{prior}{The prior probabilities used.} 
\item{means}{The class means.} 
\item{scaling}{A matrix which transforms observations to discriminant functions, normalized so that the within groups covariance matrix is spherical.}
\item{svd}{The singular values, which give the ratio of the between- and within-group standard deviations on the linear discriminant variables.  Their squares are the canonical F-statistics.}
\item{vkpt}{A vector with the indices of the variables kept in the discriminant rule if the number of variables kept is less than \sQuote{ncol(data)}. NULL otherwise.} 
\item{nvkpt}{The number of variables kept in the discriminant rule if this number is less than\sQuote{ncol(data)}. NULL otherwise.} 
\item{SSig}{An object of class \sQuote{ShrnkMat} with a compact representation of the within groups covariance matrix.
\sQuote{ShrnkMat} objects have specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations.}
\item{SSigInv}{An object of class \sQuote{ShrnkMatInv} with a compact representation of the within groups precision (inverse covariance) matrix. \sQuote{ShrnkMatInv} objects have specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations.}
\item{N}{The number of observations used.} 
\item{call}{The (matched) function call.} 

If algument \sQuote{ldafun} is set to \dQuote{classification} an object of class \sQuote{Scllda}, which extends class \sQuote{clldaRes}, with the following components:

\item{prior}{The prior probabilities used.} 
\item{means}{The class means.} 
\item{coef}{A matrix with the coefficients of the k-1 classification functions.}
\item{cnst}{A vector with the thresholds (2nd members of linear classification rules) used in classification rules that assume equal priors.}
\item{vkpt}{A vector with the indices of the variables kept in the discriminant rule if the number of variables kept is less than \sQuote{ncol(data)}. NULL otherwise.} 
\item{nvkpt}{The number of variables kept in the discriminant rule if this number is less than\sQuote{ncol(data)}. NULL, otherwise.} 
\item{SSig}{An object of class \sQuote{ShrnkMat} with a compact representation of the within groups covariance matrix.
\sQuote{ShrnkMat} objects have specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations.}
\item{SSigInv}{An object of class \sQuote{ShrnkMatInv} with a compact representation of the within groups precision (inverse covariance) matrix. \sQuote{ShrnkMatInv} objects have specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations.}
\item{N}{The number of observations used.} 
\item{call}{The (matched) function call.} 

}


\references{

Ledoit, O. and Wolf, M. (2004) \dQuote{A well-conditioned estimator for large-dimensional covariance matrices.}, \emph{Journal of Multivariate Analysis}, 88 (2), 365-411. 

Fisher, T.J. and Sun, X. (2011) \dQuote{Improved Stein-type shrinkage estimators for the high-dimensional multivariate normal covariance matrix}, \emph{Computational Statistics and Data Analysis}, 55 (1), 1909-1918. 

Pedro Duarte Silva, A. (2011) \dQuote{Two Group Classification with High-Dimensional Correlated Data: A Factor Model Approach}, 
\emph{Computational Statistics and Data Analysis}, 55 (1), 2975-2990.
}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{SelectV}}, \code{\link{ShrnkSigE}}, \code{\link{ShrnkMat}}, \code{\link{ShrnkMatInv}}, \code{\link{predict.canldaRes}}, \code{\link{predict.clldaRes}}, \code{\link{AlonDS}}}

\examples{

# train classifier on Alon's Colon Cancer Data set after a logarithmic transformation 
# (selecting genes by the Expanded HC scheme). 

ldarule <- Slda(log10(AlonDS[,-1]),AlonDS$grouping)     

# show classification rule

print(ldarule)

# get in-sample classification results

predict(ldarule,log10(AlonDS[,-1]),grpcodes=levels(AlonDS$grouping))$class           	       

# compare classifications with true assignments

cat("Original classes:\n")
print(AlonDS[,1])             		 


# Estimate error rates by four-fold cross-validation.
# (Note: In cross-validation analysis it is recommended to set the argument 
# 'ldafun' to "classification", in order to speed up computations by avoiding 
# unecessary eigen-decompositions) 

CrosValRes <- DACrossVal(log10(AlonDS[,-1]),AlonDS$grouping,TrainAlg=Slda,ldafun="classification",
kfold=4,CVrep=1)
summary(CrosValRes[,,"Clerr"])
         	       
}

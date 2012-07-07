\name{RFlda}

\alias{RFlda}
\alias{RFlda.default}
\alias{RFlda.data.frame}
\alias{is.RFlda}

\title{ High-Dimensional Factor-based Linear Discriminant Analysis.}

\description{
 \sQuote{RFlda} finds the coefficients of a linear discriminant rule based on a correlation (or covariance)
 matrix estimator that tries to approximate the true correlation (covariance) by the closest (according to 
 a Frobenius norm) correlation (covariance) compatible with a q-factor model.
}

\usage{

\method{RFlda}{default}(data, grouping, q = 1, prior = "proportions", 
CorrAp = TRUE, maxq=5, VSelfunct = SelectV,
ldafun=c("canonical","classification"), nstarts = 1, 
CVqtrials=1:3, CVqfolds=3, CVqrep=1, CVqStrt=TRUE, \dots)

\method{RFlda}{data.frame}(data, \dots)
}

\arguments{
  \item{data}{Matrix or data frame of observations.}

  \item{grouping}{Factor specifying the class for each observation.}

  \item{q}{Number of factors assumed by the model. This argument can be set to a fixed number between 1 and the argument of \sQuote{maxq}, or to the string
\dQuote{CVq}. In the latter case the number of factors is chosen amongst the values of the argument \sQuote{CVqtrials}, by minimizing a cross-validated estimate of the error rate.}

  \item{prior}{The prior probabilities of class membership.  If unspecified, the class proportions for the training set are used.  If
present, the probabilities should be specified in the order of the factor levels.}

  \item{CorrAp}{A boolean flag indicating whether the approximation error of the correlation (default), or of the covariance matrix , should be minimized.}

  \item{maxq}{Upper limit on the values allowed for argument \sQuote{q}.}

  \item{VSelfunct}{Variable selection function. Either the string \dQuote{none} (no selection is to be performed) or a function that takes \sQuote{data} and \sQuote{grouping} as its first two arguments and returns a list with two components: (i) \sQuote{nvkpt} - the number of variables to be used in the Discriminant rule; and (ii) \sQuote{vkptInd} - the indices of the variables to be used in the Discriminant rule. The default is the \sQuote{SelectV} function that, by default, selects variables by the Expanded HC scheme described in Duarte Silva (2011).} 

  \item{ldafun}{Type of discriminant linear functions computed. The alternatives are \dQuote{canonical} for maximum-discrimination canonical
linear functions and  \dQuote{classification} for direct-classification linear functions.}

  \item{nstarts}{Number of different randomly generated starting points used in the minimization of the Frobenius norm of the correlation (or covariance) matrix approximation.}

  \item{CVqtrials}{Vector of values to be tried for the number of factors assumed by the model, when argument \sQuote{q} is set to \dQuote{CVq}.}

  \item{CVqfolds}{Number of training sample folds to be created in each replication of the cross-validation procedure for choosing the number of factors, when argument \sQuote{q} is set to \dQuote{CVq}.}

  \item{CVqrep}{Number of replications to be performed in the cross-validation procedure for choosing the number of factors, when argument \sQuote{q} is set to \dQuote{CVq}.}

  \item{CVqStrt}{Boolean flag indicating if, in the cross-validation procedure for choosing the number of factors when argument \sQuote{q} is set to \dQuote{CVq}, the folds should be stratified according to the original class proportions (default), or randomly generated from the whole training sample ignoring class membership.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{

If algument \sQuote{ldafun} is set to \dQuote{canonical} an object of class \sQuote{RFcanlda}, which extends class \sQuote{canldaRes}, with the following components:

\item{prior}{The prior probabilities used.} 
\item{means}{The class means.} 
\item{scaling}{A matrix which transforms observations to discriminant functions, normalized so that the within groups covariance matrix is spherical.}
\item{svd}{The singular values, which give the ratio of the between- and within-group standard deviations on the linear discriminant variables.  Their squares are the canonical F-statistics.}
\item{vkpt}{A vector with the indices of the variables kept in the discriminant rule if the number of variables kept is less than \sQuote{ncol(data)}. NULL otherwise.} 
\item{nvkpt}{The number of variables kept in the discriminant rule if this number is less than\sQuote{ncol(data)}. NULL otherwise.} 
\item{q}{The number of o factors used in the factor model chosen.} 
\item{SigFq}{An object of class \sQuote{SigFq} with the q-factor model approximation to the within groups covariance matrix. \sQuote{SigFq} objects have specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations.}
\item{SigFqInv}{An object of class \sQuote{SigFqInv} with the q-factor model approximation to the within groups precision (inverse covariance) matrix. \sQuote{SigFqInv} objects have specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations.}
\item{N}{The number of observations used.} 
\item{call}{The (matched) function call.} 

If algument \sQuote{ldafun} is set to \dQuote{classification} an object of class \sQuote{RFcllda}, which extends class \sQuote{clldaRes}, with the following components:

\item{prior}{The prior probabilities used.} 
\item{means}{The class means.} 
\item{coef}{A matrix with the coefficients of the k-1 classification functions.}
\item{cnst}{A vector with the thresholds (2nd members of linear classification rules) used in classification rules that assume equal priors.}
\item{vkpt}{A vector with the indices of the variables kept in the discriminant rule if the number of variables kept is less than \sQuote{ncol(data)}. NULL otherwise.} 
\item{nvkpt}{The number of variables kept in the discriminant rule if this number is less than\sQuote{ncol(data)}. NULL, otherwise.} 
\item{q}{The number of o factors used in the factor model chosen.} 
\item{SigFq}{An object of class \sQuote{SigFq} with the q-factor model approximation to the within groups covariance matrix. \sQuote{SigFq} objects have specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations.}
\item{SigFqInv}{An object of class \sQuote{SigFqInv} with the q-factor model approximation to the within groups precision (inverse covariance) matrix. \sQuote{SigFqInv} objects have specialized methods for matrix inversion, multiplication, and element-wise arithmetic operations.}
\item{N}{The number of observations used.} 
\item{call}{The (matched) function call.} 

}

\references{Pedro Duarte Silva, A. (2006) \dQuote{Two Group Classification with High-Dimensional Correlated Data: A Factor Model Approach}, 
\emph{Computational Statistics and Data Analysis}, 55 (1), 2975-2990.}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{FrobSigAp}}, \code{\link{SelectV}}, \code{\link{SigFq}}, \code{\link{SigFqInv}}, \code{\link{predict.canldaRes}}, \code{\link{predict.clldaRes}},
\code{\link{AlonDS}}}

\examples{

#train classifier with 10 genes (after a logarithmic transformation) 
# on Alon's Colon Cancer Data set. 

log10genes <- log10(AlonDS[,-1])

ldarule1 <- RFlda(log10genes,AlonDS$grouping,Selmethod="fixedp",maxp=10)     

# get in-sample classification results

predict(ldarule1,log10genes,grpcodes=levels(AlonDS$grouping))$class           	       

# compare classifications with true assignments

cat("Original classes:\n")
print(AlonDS$grouping)             		 

# Estimate error rates by four-fold cross-validation.
# Note: In cross-validation analysis it is recommended to set 
# the argument 'ldafun' to "classification", in order to speed up 
# computations by avoiding unecessary eigen-decompositions 

CrosValRes1 <- DACrossVal(log10genes,AlonDS$grouping,TrainAlg=RFlda,
Selmethod="fixedp",ldafun="classification",maxp=10,kfold=4,CVrep=1)
summary(CrosValRes1[,,"Clerr"])

# Find the best factor model amongst the choices q=1 or 2

ldarule2 <- RFlda(log10genes,AlonDS$grouping,q="CVq",CVqtrials=1:2,
Selmethod="fixedp",ldafun="classification",maxp=10)
cat("Best error rate estimate found with q =",ldarule2$q,"\n")

# Perform the analysis finding the number of selected genes 
# by the Expanded HC scheme 

ldarule3 <- RFlda(log10genes,AlonDS$grouping,q=ldarule2$q)     
cat("Number of selected genes =",ldarule3$nvkpt,"\n")

# get classification results

predict(ldarule3,log10genes,grpcodes=levels(AlonDS$grouping))$class           	       
}

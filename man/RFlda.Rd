\name{RFlda}

\alias{RFlda}
\alias{RFlda.default}
\alias{RFlda.data.frame}
\alias{print.RFlda}
\alias{coef.RFlda}

\title{ High-Dimensional Factor-based Linear Discriminant Analysis.}

\description{
 The function \sQuote{RFlda} finds the coefficients of a linear discriminant rule based on covariance matrix estimator
 that tries to approximate the true covariance by the closest (according to a Forbenius norm) covariance 
 compatible with a q-factor model.
}

\usage{

\method{RFlda}{default}(x, grouping, q = 1, prior = "proportions", CorrAp = TRUE, 
maxq=5, VSelfunct = SelectV, nstarts = 1, 
CVqtrials=1:3, CVqfolds=3, CVqrep=1, CVqStrt=TRUE, \dots)


\method{RFlda}{data.frame}(x, \dots)

}

\arguments{
  \item{x}{Matrix or data frame of observations.}
  \item{grouping}{Factor specifying the class for each observation.}
  \item{q}{Number of factors assumed by the model. This argument can be set to a fixed number between 1 and the argument of \sQuote{maxq}, or to the string
\dQuote{CVq}. In the latter case the number of factors is chosen amongst the values of the argument \sQuote{CVqtrials}, by minimizing a cross-validated estimate of the error rate.}
  \item{prior}{The prior probabilities of class membership.  If unspecified, the class proportions for the training set are used.  If
present, the probabilities should be specified in the order of the factor levels.}
  \item{CorrAp}{A boolean flag indicating whether the approximation error of the correlation (default), or of the covariance matrix , should be minimized.}
  \item{maxq}{Upper limit on the values allowed for argument \sQuote{q}.}
  \item{VSelfunct}{Variable selection function. Either the string \dQuote{none} (no selection is to be performed) or a function that takes \sQuote{x} and \sQuote{grouping} as its first two arguments and returns a list with two components: (i) \sQuote{nvkpt} - the number of variables to be used in the Discriminant rule; and (ii) \sQuote{vkptInd} - the indices of the variables to be used in the Discriminant rule. The default is the \sQuote{SelectV} function that, by default, selects variables by the Extended HC scheme described in Duarte Silva (forthcoming).} 
  \item{nstarts}{Number of different randomly generated starting points used in the minimization of the Forbenius norm of the correlation (or covariance) matrix approximation.}
  \item{CVqtrials}{Vector of values to be tried for the number of factors assumed by the model when argument \sQuote{q} is set to \dQuote{CVq}.}
  \item{CVqfolds}{Number of training sample folds to be created in each replication of the cross-validation procedure for choosing the number of factors when argument \sQuote{q} is set to \dQuote{CVq}.}
  \item{CVqrep}{Number of replications to be performed in the cross-validation procedure for choosing the number of factors when argument \sQuote{q} is set to \dQuote{CVq}.}
  \item{CVqStrt}{Boolean flag indicating if, in the cross-validation procedure for choosing the number of factors when argument \sQuote{q} is set to \dQuote{CVq}, the folds should be stratified according to the original class proportions (default), or randomly generated from the whole training sample ignoring class membership.}
  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{An object of class \sQuote{RFlda} with the following components:

\item{prior}{The prior probabilities used.} 
\item{means}{The class means.} 
\item{coef}{A matrix with the coefficients of the k-1 discriminant rules.}
\item{cnst}{A vector with the thresholds (2nd members of linear discriminants) used in classification rules that assume equal priors.}
\item{vkpt}{A vector with the indices of the variables kept in the discriminant rule if total number of variables is less than ncol(x). NULL, otherwise.} 
\item{nvkpt}{The number of variables kept in the discriminant rule if total number of variables is less than ncol(x). NULL, otherwise.} 
\item{q}{The number of o factors used in the factor model chosen.} 
\item{SigFq}{An object of class \sQuote{SigFq} with the q-factor model approximation to the within groups covariance matrix.}
\item{SigFqInv}{An object of class \sQuote{SigFqInv} with the q-factor model approximation to the within groups precision (inverse covariance) matrix.}
\item{N}{The number of observations used.} 
\item{call}{The (matched) function call.} 
}

\references{Duarte Silva, A.P. (Forthcoming) Two Group Classification with High-Dimensional Correlated Data: A Factor Model Approach. 
\emph{Computational Statistics and Data Analysis.}}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{ForbSigap}}, \code{\link{SelectV}}, \code{\link{SigFq}}, \code{\link{SigFqInv}}, \code{\link{predict.RFlda}},
\code{\link{AlonDS}}}

\examples{

#train classifier with 10 genes on Alon's Colon Cancer Data set. 
ldarule1 <- RFlda(AlonDS[,-1],AlonDS[,1],Selmethod="fixedp",maxp=10)     

# get in-sample classification results
predict(ldarule1,AlonDS[,-1],grpcodes=levels(AlonDS[,1]))$class           	       

# compare classifications with true assignments
cat("Original classes:\n")
print(AlonDS[,1])             		 

# Estimate error rates by five-fold cross-validation replicated twice.
CrosValRes1 <- DACrossVal(AlonDS[,-1],AlonDS[,1],TrainAlg=RFlda,
Selmethod="fixedp",maxp=10,kfold=5,CVrep=2)
summary(CrosValRes1[,,"Clerr"])

#Find the best factor model amongst the choices q=1 or 2
ldarule2 <- RFlda(AlonDS[,-1],AlonDS[,1],q="CVq",CVqtrials=1:2,
Selmethod="fixedp",maxp=10)
cat("Best error rate estimate found with q =",ldarule2$q,"\n")

# Find the number of selected genes by the Extended HC scheme 
ldarule3 <- RFlda(AlonDS[,-1],AlonDS[,1],q=ldarule2$q)     
cat("Number of selected genes =",ldarule3$nvkpt,"\n")

# get classification results
predict(ldarule3,AlonDS[,-1],grpcodes=levels(AlonDS[,1]))$class           	       
}

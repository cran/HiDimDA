\name{DACrossVal}

\title{Cross Validation for Discriminant Analysis Classification Rules}

\alias{DACrossVal}

\description{
\sQuote{DACrossVal} evaluates the performance of a Discriminant Analysis training sample algorithm by kfold 
Cross-Validation.  
}

\usage{DACrossVal(x, grouping, TrainAlg, EvalAlg=EvalClrule, 
Strfolds=TRUE, kfold=10, CVrep=20, prior="proportions", \dots)}

\arguments{
  \item{x}{Matrix of observations.}
  \item{grouping}{Factor specifying the class for each observation.}
  \item{TrainAlg}{A function with the training algorithm. It should return a list that can be used as input to the argument of \sQuote{EValAlg}.}
  \item{EvalAlg}{A function with the evaluation algorithm. By default set to \sQuote{EvalClrule} which returns a list with components \dQuote{err} (estimates of error rates by class) and \dQuote{Nk} (number of out-sample observations by class). This default can be used for all \sQuote{TrainAlg} 
arguments that return an object with a predict method returning a list with a \sQuote{class} component (a factor) containing the classification results.}
  \item{Strfolds}{Boolean flag indicating if the folds should be stratified according to the original class
proportions (default), or randomly generated from the whole training sample, ignoring class membership.}
  \item{kfold}{Number of training sample folds to be created in each replication.}
  \item{CVrep}{Number of replications to be performed.}
  \item{prior}{The prior probabilities of class membership.  If unspecified, the class proportions for the training set are used.  If present, the probabilities should be specified in the order of the factor levels.}
  \item{\dots}{Further arguments to be passed to \sQuote{TrainAlg} and \sQuote{EvalAlg}.} 
}

\value{A three dimensional array with the number of tested observations, and estimated classification errors for each combination of fold and replication tried. The array dimensions are defined as follows.\cr
The first dimension runs through the different fold-replication combinations, the second dimension represents the classes, and the third dimension has two named levels representing respectively the number of observations tested (\dQuote{Nk}), and the estimated classification errors (\dQuote{Clerr}).}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{lda}}, \code{\link{RFlda}}, \code{\link{predict.lda}}, \code{\link{predict.RFlda}}, \code{\link{AlonDS}}}

\examples{

# Evaluate the performance of traditional (Fisher's) linear discriminant
# analysis on the iris data set, by ten-fold cross-validation replicated 
# three times.

library(MASS)
CrosValRes1 <- DACrossVal(iris[1:4],iris$Species,TrainAlg=lda,CVrep=3)
summary(CrosValRes1[,,"Clerr"])

# Evaluate the performance on Alon's Colon Cancer Data set, of a 
# one-factor linear discriminant rule after variable selection 
# by the Extended HC scheme, by five-fold cross-validation.

CrosValRes2 <- DACrossVal(AlonDS[,-1],AlonDS[,1],TrainAlg=RFlda,
kfold=5,CVrep=1)
summary(CrosValRes2[,,"Clerr"])

}

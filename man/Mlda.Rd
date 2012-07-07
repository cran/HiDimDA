\name{Mlda}

\alias{Mlda}
\alias{Mlda.default}
\alias{Mlda.data.frame}
\alias{is.Mlda}

\title{Maximum uncertainty Linear Discriminant Analysis.}

\description{
 \sQuote{Mlda} finds the coefficients of a linear discriminant rule based on the \dQuote{Maximum uncertainty Linear Discriminant Analysis} approach of Thomaz, Kitani and Gillies (2006).
}


\usage{

\method{Mlda}{default}(data, grouping, prior = "proportions", StddzData=TRUE, 
VSelfunct = SelectV, ldafun=c("canonical","classification"), 
PCAstep=FALSE, \dots)

\method{Mlda}{data.frame}(data, \dots)
}

\arguments{
  \item{data}{Matrix or data frame of observations.}

  \item{grouping}{Factor specifying the class for each observation.}

  \item{prior}{The prior probabilities of class membership.  If unspecified, the class proportions for the training set are used.  If
present, the probabilities should be specified in the order of the factor levels.}

  \item{StddzData}{A boolean flag indicating whether the data should be standardized first (default) or used in their original scales.}

  \item{VSelfunct}{Variable selection function. Either the string \dQuote{none} (no selection is to be performed) or a function that takes \sQuote{data} and \sQuote{grouping} as its first two arguments and returns a list with two components: (i) \sQuote{nvkpt} - the number of variables to be used in the Discriminant rule; and (ii) \sQuote{vkptInd} - the indices of the variables to be used in the Discriminant rule. The default is the \sQuote{SelectV} function that, by default, selects variables by the Expanded HC scheme described in Duarte Silva (2011).} 

  \item{ldafun}{Type of discriminant linear functions computed. The alternatives are \dQuote{canonical} for maximum-discrimination canonical
linear functions and  \dQuote{classification} for direct-classification linear functions.}

  \item{PCAstep}{A flag indicating if data should be first projected into the space spanned by its first nrow(data)-1 Principal Components in problems where nrow(data)-1 is less than the number of selected variables. In applications with a very large number of useful variables seting PCAstep to TRUE avoids many potential memory problems and tends to substantially increase the size of the data sets that can be analyzed by Mlda.}

  \item{\dots}{Further arguments passed to or from other methods.}
}

\value{

If algument \sQuote{ldafun} is set to \dQuote{canonical} an object of class \sQuote{canldaRes} with the following components:

\item{prior}{The prior probabilities used.} 
\item{means}{The class means.} 
\item{scaling}{A matrix which transforms observations to discriminant functions, normalized so that the within groups covariance matrix is spherical.}
\item{svd}{The singular values, which give the ratio of the between- and within-group standard deviations on the linear discriminant variables.  Their squares are the canonical F-statistics.}
\item{vkpt}{A vector with the indices of the variables kept in the discriminant rule if the number of variables kept is less than \sQuote{ncol(data)}. NULL otherwise.} 
\item{nvkpt}{The number of variables kept in the discriminant rule if this number is less than\sQuote{ncol(data)}. NULL otherwise.} 
\item{N}{The number of observations used.} 
\item{call}{The (matched) function call.} 

If algument \sQuote{ldafun} is set to \dQuote{classification} an object of class \sQuote{clldaRes} with the following components:

\item{prior}{The prior probabilities used.} 
\item{means}{The class means.} 
\item{coef}{A matrix with the coefficients of the k-1 classification functions.}
\item{cnst}{A vector with the thresholds (2nd members of linear classification rules) used in classification rules that assume equal priors.}
\item{vkpt}{A vector with the indices of the variables kept in the discriminant rule if the number of variables kept is less than \sQuote{ncol(data)}. NULL otherwise.} 
\item{nvkpt}{The number of variables kept in the discriminant rule if this number is less than\sQuote{ncol(data)}. NULL, otherwise.} 
\item{N}{The number of observations used.} 
\item{call}{The (matched) function call.} 

}

\references{

Pedro Duarte Silva, A. (2011) \dQuote{Two Group Classification with High-Dimensional Correlated Data: A Factor Model Approach}, 
\emph{Computational Statistics and Data Analysis}, 55 (1), 2975-2990.

Thomaz, Kitani and Gillies (2006) \dQuote{A maximum uncertainty LDA-based approach for limited sample size problems - with application to face recognition}, \emph{Journal of the Brazilian Computer Society}, 12 (2), 7-18.

}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{SelectV}}, \code{\link{MldaInvE}}, \code{\link{predict.canldaRes}}, \code{\link{predict.clldaRes}}, \code{\link{AlonDS}}}

\examples{

# train classifier on Alon's Colon Cancer Data Set 
# (after a logarithmic transformation). 

log10genes <- log10(AlonDS[,-1])

ldarule <- Mlda(log10genes,AlonDS$grouping)     

# show classification rule

print(ldarule)

# get in-sample classification results

predict(ldarule,log10genes,grpcodes=levels(AlonDS$grouping))$class           	       

# compare classifications with true assignments

cat("Original classes:\n")
print(AlonDS$grouping)             		 

# Estimate error rates by four-fold cross-validation.
# Note: In cross-validation analysis it is recommended to set 
# the argument 'ldafun' to "classification", in order to speed up 
# computations by avoiding unecessary eigen-decompositions 

CrosValRes <- DACrossVal(log10genes,AlonDS$grouping,TrainAlg=Mlda,
ldafun="classification",kfold=4,CVrep=1)
summary(CrosValRes[,,"Clerr"])
        	       
}

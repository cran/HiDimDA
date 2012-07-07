\name{HiDimDA-package}

\alias{HiDimDA-package}
\alias{HiDimDA}

\docType{package}

\title{High Dimensional Discriminant Analysis}

\description{Performs Linear Discriminant Analysis in High Dimensional problems based on reliable covariance estimators for problems with
 (many) more variables than observations. Includes routines for classifier training, prediction, cross-validation and variable selection.}

\details{
\tabular{ll}{
Package: \tab HiDimDA\cr
Type: \tab Package\cr
Version: \tab 0.2-1\cr
Date: \tab 2012-07-07\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
LazyData: \tab yes\cr
}

HiDimDA is a package for High-Dimensional Discriminant Analysis aimed at problems with many variables, possibly much more
than the number of available observations. Its core consists of the four Linear Discriminant Analyis routines:

\tabular{ll}{
 Dlda: \tab Diagonal Linear Discriminant Analysis\cr
Slda: \tab Shrunken Linear Discriminant Analysis\cr
Mlda: \tab Maximum-uncertainty Linear Discriminant Analysis\cr
RFlda: \tab Factor-model  Linear Discriminant Analysis\cr
}

and the variable selection routine:

\tabular{ll}{
SelectV: \tab High-Dimensional variable selection for supervised classification\cr 
}

that selects variables to be used in a Discriminant classification rule by
ranking them according to two-sample t-scores (problems with two-groups),
or ANOVA F-scores (problems wih more that two groups), and discarding those
with scores below a threshold defined by the Higher Criticism (HC) approach 
of Donoho and Jin (2008), the Expanded Higher Criticism scheme 
proposed by Duarte Silva (2011), False Discovery Rate (Fdr) control as suggested by 
Benjamini and Hochberg (1995), the FAIR approach of Fan and Fan (2008), or simply by 
fixing the number of retained variables to some pre-defined constant.


All four discriminant routines, \sQuote{Dlda}, \sQuote{Slda}, \sQuote{Mlda} and \sQuote{RFlda}, compute Linear
Discriminant Functions, by default after a preliminary variable selection step, based on alternative estimators of 
a within-groups covariance matrix that leads to reliable allocation rules in problems where the number of selected
variables is close to, or larger than, the number of available observations.  


Consider a Discriminant Analysis problem with \eqn{k} groups, \eqn{p} selected variables, a training sample consisting
of \eqn{N = \sum_{g=1}^{k}n_g} observations with group and overall means, 
\eqn{\bar{X}_g} and \eqn{\bar{X}_.}, and a between-groups scatter (scaled by degrees of freedom) 
matrix, \eqn{S_B = \frac{1}{N-k} \sum_{g=1}^{k} n_g (\bar{X}_g -\bar{X}_.)(\bar{X}_g -\bar{X}_.)^T }

Following the two main classical approaches to Linear Discrimant Analysis, the Discriminant Functions returned by HiDimDA discriminant
routines are either based on the canonical linear discriminants given by the normalized eigenvectors

\deqn{LD_j = Egvct_j (S_B \hat{\Sigma}_W^{-1})}  
\deqn{j = 1,...,r=min(p,k-1)}
\deqn{[LD_1, ..., LD_r]^T \hat{\Sigma}_W [LD_1, ..., LD_r] = I_r } 

or the classification functions 

\deqn{CF_g = (\bar{X}_g - \bar{X}_1) \hat{\Sigma}_W^{-1}} 
\deqn{g = 2,...,k} 

where \eqn{\hat{\Sigma}_W^{-1}} is an estimate of  the inverse within-groups covariance. 

It is well known that these two approaches are equivalent, in the sense that classification rules that assign new observations to
the group with the closest (according to the Euclidean distance) centroid in the space of the canonical variates, 
\eqn{Z = [LD_1 ... LD_r]^T X }, give the same results as the rule that assigns a new observation to group 1 if all classification scores, 
\eqn{Clscr_g = CF_g^T X - CF_g^T  \frac{(\bar{X}_1 + \bar{X}_g)}{2} }, are negative, and to the group with the highest classification 
score otherwise.

The discriminant routines of HiDimDA compute canonical linear discriminant functions by default, and classification functions when
the argument \sQuote{ldafun} is set to \dQuote{classification}. However, unlike traditional linear discriminant analysis where
\eqn{\Sigma_W^{-1}} is estimated by the inverse of the sample covariance, 
which is not well-defined when \eqn{p \geq N-k} and is unreliable if \eqn{p} is close to \eqn{N-k}, the routines of HiDimDA use 
four alternative well-conditioned estimators of \eqn{\Sigma_W^{-1}} that lead to reliable classification rules if \eqn{p} is larger than, 
or close to, \eqn{N-k}.      

In particular, \sQuote{Dlda} estimates \eqn{\Sigma_W^{-1}} by the diagonal matrix of inverse sample variances, \sQuote{Slda} by
the inverse of an optimally shrunken Ledoit and Wolf's (2004) covariance estimate with the targets and optimal
target intensity estimators proposed by Fisher and Sun (2011), \sQuote{Mlda} uses a regularized inverse
covariance that deemphasizes the importance given to the last eigenvectors of the sample covariance (see Thomaz, Kitani
and Gillies (2006) for details), and \sQuote{RFlda} uses a factor model estimate of the true inverse correlation (or covariance)
matrix based on the approach of Duarte Silva (2011).

The HiDimDA package also includes predict methods for all discriminant routines implemented, a routine (\sQuote{DACrossVal}) for asssessing
the quality of the classification results by kfold cross-validation, and utilities for storing, extracting and efficiently handling  specialized high-dimensional covariance and inverse covariance matrix estimates. 

}

\author{Antonio Pedro Duarte Silva  <psilva@porto.ucp.pt>

Maintainer: Antonio Pedro Duarte Silva  <psilva@porto.ucp.pt>
}

\seealso{\code{\link{Dlda}}, \code{\link{Mlda}}, \code{\link{Slda}},\code{\link{RFlda}}, \code{\link{predict.canldaRes}}, \code{\link{predict.clldaRes}}, \code{\link{AlonDS}}}

\references{

Benjamini, Y. and Hochberg, Y. (1995) \dQuote{Controling the false discovery rate: A practical and powerful
approach to multiple testing}, \emph{Journal of the Royal Statistical Society} B, 57, 289-300.

Donoho, D. and Jin, J. (2008) \dQuote{Higher criticism thresholding: Optimal feature selection when useful
features are rare and weak}, In: \emph{Proceedings National Academy of Sciences}, USA 105, 14790-14795.

Fan, J. and Fan, Y. (2008) \dQuote{High-dimensional classification using features annealed independence rules},
\emph{Annals of Statistics}, 36 (6), 2605-2637.

Fisher, T.J. and Sun, X. (2011) \dQuote{Improved Stein-type shrinkage estimators for the high-dimensional multivariate normal covariance matrix}, \emph{Computational Statistics and Data Analysis}, 55 (1), 1909-1918. 

Ledoit, O. and Wolf, M. (2004) \dQuote{A well-conditioned estimator for large-dimensional covariance matrices.}, \emph{Journal of Multivariate Analysis}, 88 (2), 365-411. 

Pedro Duarte Silva, A. (2011) \dQuote{Two Group Classification with High-Dimensional Correlated Data: A Factor Model Approach}, 
\emph{Computational Statistics and Data Analysis}, 55 (1), 2975-2990.

Thomaz, C.E. Kitani, E.C. and Gillies, D.F. (2006) \dQuote{A maximum uncertainty LDA-based approach for limited sample size problems - with application to face recognition}, \emph{Journal of the Brazilian Computer Society}, 12 (2), 7-18

}

\keyword{HiDimDA}

\examples{

# train the four main classifiers with their default setings 
# on Alon's colon data set (after a logarithmic transformation), 
# selecting genes by the Expanded HC scheme 

# Pre-process and select the genes to be used in the classifiers

log10genes <- log10(AlonDS[,-1]) 
SelectionRes <- SelectV(log10genes,AlonDS$grouping)
genesused <- log10genes[SelectionRes$vkpt]

# Train classifiers

DiaglldaRule <- Dlda(genesused,AlonDS$grouping)     
FactldaRule <- RFlda(genesused,AlonDS$grouping)     
MaxUldaRule <- Mlda(genesused,AlonDS$grouping)     
ShrkldaRule <- Slda(genesused,AlonDS$grouping)     

# Get in-sample classification results

predict(DiaglldaRule,genesused,grpcodes=levels(AlonDS$grouping))$class         	       
predict(FactldaRule,genesused,grpcodes=levels(AlonDS$grouping))$class         	       
predict(MaxUldaRule,genesused,grpcodes=levels(AlonDS$grouping))$class         	       
predict(ShrkldaRule,genesused,grpcodes=levels(AlonDS$grouping))$class         	       

# Compare classifications with true assignments

cat("Original classes:\n")
print(AlonDS$grouping)             		 

# Show set of selected genes

cat("Genes kept in discrimination rule:\n")
print(colnames(genesused))             		 
cat("Number of selected genes =",SelectionRes$nvkpt,"\n")
}


\name{SelectV}

\alias{SelectV}

\title{Variable Selection for High-Dimensional Supervised Classification.}

\description{Selects variables to be used in a Discriminant Analysis classification rule.}

\details{The function \sQuote{SelectV} selects variables to be used in a Discriminant classification rule by the 
 Higher Criticism (HC) approach of Donoho and Jin (2004, 2008), the Expanded Higher Criticism scheme 
 proposed by Duarte Silva (2011), False Discovery Rate (Fdr) control as suggested by 
 Benjamini and Hochberg (1995), the FAIR (Features Annealed Independence Rules) approach of Fan and Fan (2008),
 or simply by fixing the number of selected variables to some pre-defined constant.

 The Fdr method is, by default, based on simple p-values  derived from t-scores (problems with two groups) or ANOVA F-scores 
 (problems with more than two groups). 
 When the argument \sQuote{NullDist} is set to \dQuote{Theoretical} these values are also used in the HC method. Otherwise, 
 the HC p-values are derived from an empirical Null of z-scores estimated by the maximum likelihood 
 approach of Efron (2004). 

The variable rankings are based on absolute-value t-scores or ANOVA F-scores.
}

\usage{SelectV(data, grouping, 
Selmethod=c("ExpHC","HC","Fdr","Fair","fixedp"),
NullDist=c("locfdr","Theoretical"), uselocfdr=c("onlyHC","always"), 
minlocfdrp=200, comvar=TRUE, Fdralpha=0.5, 
ExpHCalpha=0.5, HCalpha0=0.1, maxp=ncol(data), tol=1E-12, \dots)}

\arguments{

  \item{data}{Matrix or data frame of observations.}

  \item{grouping}{Factor specifying the class for each observation.}

  \item{Selmethod}{The method used to choose the number of variables selected. Current alternatives are:

	\sQuote{ExpHC} (default) for the Expanded Higher Criticism scheme of Duarte Silva (2011)

	\sQuote{HC} for the Higher Criticism (HC) approach of Donoho and Jin (2004, 2008)

	\sQuote{Fdr} for the False Discovery Rate control approach of Benjamini and Hochberg (1995)

	\sQuote{Fair} for the FAIR (Features Annealed Independence Rules) approach of Fan and Fan (2008). This option is only available for two-group classification problems. 

	\sQuote{fixedp} for a constant chosen by the user.
   }

  \item{NullDist}{The Null distribution used to compute pvalues from t-scores or F-scores. Current alternatives are \dQuote{Theoretical} for the corresponding theoretical distributions, and \dQuote{locfdr} for an empirical Null of z-scores estimated by the maximum likelihood approach of Efron (2004).}

  \item{uselocfdr}{Flag indicating the statistics for which the Null empirical distribution estimated by the locfdr approach should be used. Current alternatives are \dQuote{onlyHC} (default) and \dQuote{always}.}

  \item{minlocfdrp}{Minimum number of variables required to estimate empirical Null distributions by the locfdr method. When the number of variables is below \sQuote{minlocfdrp}, theoretical Nulls are always employed.} 

  \item{comvar}{Boolean flag indicating if a common group variance is to be assumed (default) in the computation of the t-scores used for problems with two groups.}

  \item{Fdralpha}{Control level for variable selection based on False Discovery Rate Control (see Benjamini and Hochberg (1995)).}

  \item{ExpHCalpha}{Control level for the first step of the Extended Higher Criticism scheme (see Duarte Silva (2011)).}

  \item{HCalpha0}{Proportion of pvalues used to compute the HC statistic (see Donoho and Jin (2004, 2008)).}

  \item{maxp}{Maximum number of variables to be used in the discriminant rule.} 

  \item{tol}{Numerical precision for distinguishing pvalues from 0 and 1. Computed pvalues below \sQuote{tol} are set to \sQuote{tol}, and those above 1-\sQuote{tol} are set to 1-\sQuote{tol}.}

  \item{\dots}{Arguments passed from other methods.}
}

\value{A list with two components: 

\item{nvkpt}{the number of variables to be used in the Discriminant rule} 
\item{vkptInd}{ the indices of the variables to be used in the Discriminant rule} 
}

\references{

Benjamini, Y. and Hochberg, Y. (1995) \dQuote{Controling the false discovery rate: A practical and powerful
approach to multiple testing}, \emph{Journal of the Royal Statistical Society} B, 57, 289-300.

Donoho, D. and Jin, J. (2004) \dQuote{Higher criticism for detecting sparse heterogeneous mixtures},
\emph{Annals of Statistics} 32, 962-964.

Donoho, D. and Jin, J. (2008) \dQuote{Higher criticism thresholding: Optimal feature selection when useful
features are rare and weak}, In: \emph{Proceedings National Academy of Sciences}, USA 105, 14790-14795.

Efron, B. (2004) \dQuote{Large-scale simultaneous hypothesis testing: the choice of a null hypothesis},
\emph{Journal of the American Statistical Association} 99, 96-104.

Fan, J. and Fan, Y. (2008) \dQuote{High-dimensional classification using features annealed independence rules},
\emph{Annals of Statistics}, 36 (6), 2605-2637.

Pedro Duarte Silva, A. (2011) \dQuote{Two Group Classification with High-Dimensional Correlated Data: A Factor Model Approach}, 
\emph{Computational Statistics and Data Analysis}, 55 (1), 2975-2990.

}

\author{A. Pedro Duarte Silva}

\seealso{\code{\link{Dlda}}, \code{\link{Mlda}}, \code{\link{Slda}}, \code{\link{RFlda}}, \code{\link{AlonDS}}}

\examples{

\dontrun{

# Compare the number of variables selected by the four methods 
# currently available  on Alon's Colon Cancer Data set 
# after a logarithmic transformation. 

log10genes <- log10(AlonDS[,-1])

Res <- array(dim=4)
names(Res) <- c("ExpHC","HC","Fdr","Fair")
Res[1] <- SelectV(log10genes,AlonDS[,1])$nvkpt
Res[2] <- SelectV(log10genes,AlonDS[,1],Selmethod="HC")$nvkpt
Res[3] <- SelectV(log10genes,AlonDS[,1],Selmethod="Fdr")$nvkpt
Res[4] <- SelectV(log10genes,AlonDS[,1],Selmethod="Fair")$nvkpt

print(Res)
}

}




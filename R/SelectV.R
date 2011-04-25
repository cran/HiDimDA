### SelectV.R  (2011-04-25)
###    
###
### Copyright 2011 A. Pedro Duarte Silva
###
### This file is part of the `HiDimDA' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

SelectV <- function(x,grouping,Selmethod=c("ExtHC","HC","Fdr","fixedp"),
NullDist=c("locfdr","Theoretical"),uselocfdr=c("onlyHC","always"),minlocfdrp=200,comvar=TRUE,alpha=0.1,alpha0=0.1,maxp=ncol(x),tol=1E-12,...)
{
   Selmethod <- match.arg(Selmethod)
   NullDist <-match.arg(NullDist)
   uselocfdr <- match.arg(uselocfdr)
   if (NullDist != "locfdr" && uselocfdr == "always") stop("Error: uselocfdr argument can only be used when NullDist is set to locfdr")

   p <- ncol(x)
   if (p < minlocfdrp) NullDist <- "Theoretical"
   nk <- table(grouping)
   k <- nrow(nk)
   nk <- as.vector(nk)
   n <- sum(nk)
   if (k==2) {
	tscr <- tscores(x,grouping,n,nk,comvar=comvar)
	scores <- abs(tscr$st)
   }
   else {
	fscr <- fscores(x,grouping,n,nk,k)
	scores <- fscr$st
   }
   if (Selmethod=="fixedp")  {
      	sortedscr <- sort(scores,decreasing=TRUE,index.return=TRUE)
	return(list(nvkpt=maxp,vkptInd=sort(sortedscr$ix[1:maxp])))
   }
   if (k==2) pvalues <- 1 - pt(scores,tscr$df)
   else pvalues <- 1 - pf(scores,k-1,fscr$df)
   pvalues[pvalues<tol] <- tol   	# ensure that no pvalue is numerically identical to zero.
   pvalues[pvalues>1-tol] <- 1 - tol   	# ensure that no pvalue is numerically identical to one.
   if (NullDist=="locfdr" && uselocfdr == "always") pvalues <- locfdrpval(pvalues)
   if (Selmethod=="ExtHC" || Selmethod=="Fdr")  
   {
       	sortedpv <- sort(pvalues,index.return=TRUE)
	if (Selmethod=="ExtHC") usefullpv <- sortedpv$x[sortedpv$x<alpha*1:p/(p*sum(1/1:p))]
	else usefullpv <- sortedpv$x[sortedpv$x<alpha*1:p/p]
	if (length(usefullpv)==0) Fdrnvar <- 1
	else {
		maxpv <- max(usefullpv)
		Fdrnvar <- min(maxp,which(sortedpv$x==maxpv))
	}
   }
   if (NullDist=="locfdr" && uselocfdr == "onlyHC") pvalues <- locfdrpval(pvalues)
   {
 	if (Selmethod == "ExtHC" || Selmethod == "HC" || Selmethod == "HCplus") {
		if (Selmethod== "ExtHC") minvar <- Fdrnvar
		else minvar <- 1
 		HCres <- HC(p,pvalues,minvkpt=minvar,HCplus=FALSE,alpha0=min(alpha0,maxp/p))
		return(list(nvkpt=HCres$nkptvar,vkptInd=HCres$varkept))
	}
	else {
		if (Selmethod == "Fdr") nkptvar = Fdrnvar
		return(list(nvkpt=nkptvar,vkptInd=sortedpv$ix[1:nkptvar]))
  	}
   } 
}

locfdrpval <- function(pvalues)
{
	zscores <- qnorm(pvalues)	#  Note:  This is diferent from defining zscores directly from tscores
	empnull <- mylocfdr(zscores,plot=0,silently=TRUE)
	if (class(empnull)=="error1") empnull <- mylocfdr(zscores,plot=0,nulltype=2,silently=TRUE)
	if (class(empnull)=="error3") empnull <- mylocfdr(zscores,plot=0,nulltype=1,silently=TRUE)
	if (class(empnull)!="error2")  { 
		zscores <- (zscores-empnull$fp0[3,1])/empnull$fp0[3,2]
		pvalues <- pnorm(zscores)
	}
	return(pvalues)
}

HC <- function(p,pvalues,HCplus=FALSE,minvkpt=1,alpha0=0.1)
{
#    Computes ANOVA Donoho and Jin's Higher Criticism threshold

#    Arguments:

#        p       -  the original number of variables 
#        pvalues -  a set of p pvalues  
#        HCplus  -  a boolean flag indicating if the HCplus version (always keep variables with pvalues below 1/p)
#                   of the HC criterion should be used 
#        minvkpt -  a minimum number of variables to be kept in the analysis
#        alpha0  -  the maximum percentage of variables to be kept in the analysis (by default 10%)

# Value:  a list with the three components:

#        threshold  -  the value of the Higher Criticism threshold (measured on the pvalue scale)
#        varkept    -  a vector with the indices of the variables to be kept in the analysis 
#        nkptvar     - the number of variables to be kept in the analysis 

        sortedpv <- sort(pvalues,index.return=TRUE)
        if (HCplus) p0 <- max(minvkpt,length(sortedpv$x[sortedpv$x<=1/p])+1)
        else p0 <- minvkpt
        p1 <- floor(alpha0*p)
	if (p0 >= p1) nkptvar <- p0
        else  {
		unifq <- (p0:p1)/p
		HC <- p * (unifq-sortedpv$x[p0:p1]) / sqrt( (p0:p1)*(1-unifq) )
		if (max(HC)>0.) nkptvar <- which.max(HC)+p0-1
		else nkptvar <- p0
	}
  	XPind <- sort(sortedpv$ix[1:nkptvar])

        list(threshold=sortedpv$x[nkptvar],varkept=XPind,nkptvar=nkptvar) # return(list(threshold=sortedpv$x[nkptvar],varkept=XPind,nkptvar=nkptvar))
}

tscores <- function(x,grouping,n,nk,comvar)
{
  # Computes two-group t-scores

  Xbark <- apply(x,2,grpmeans,grp=grouping)
  vark <- apply(x,2,grpvar,grp=grouping)
  if (comvar==TRUE) {
	df <- n-2
	denom <- sqrt( (1/nk[1]+1/nk[2]) * ((nk[1]-1)*vark[1,]+(nk[2]-1)*vark[2,]) / df )
  }
  else {
	tmp1 <- vark[1,]/nk[1]
	tmp2 <- vark[2,]/nk[2]
	tmps <- tmp1 + tmp2
	df <- round( tmps^2/ ( tmp1^2/(nk[1]-1)+tmp2^2/(nk[2]-1) ) )
  	denom <- sqrt(tmps)
   }
   list(st=(Xbark[1,]-Xbark[2,])/denom,df=df)  #  return( list(st=(Xbark[1,]-Xbark[2,])/denom,df=df) )
}

fscores <- function(x,grouping,n,nk,k)
{
  # Computes ANOVA f-scores

  df <- n - k
  vark <- apply(x,2,grpvar,grp=grouping)
  W <- apply((nk-1)*vark,2,sum)
  B <- (n-1)*apply(x,2,var) - W
  list(st=(B/(k-1))/(W/df),df=df)   # return(list(st=(B/(k-1))/(W/df),df=df))
}



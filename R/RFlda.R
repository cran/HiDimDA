### RFlda.R  (2011-06-13)
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

RFlda <- function(data,...) {
  if (is.null(class(data))) class(data) <- data.class(data)
  UseMethod("RFlda") 
}   

RFlda.default <- function(data,grouping,q=1,prior="proportions",CorrAp=TRUE,maxq=5,VSelfunct=SelectV,nstarts=1,
				CVqtrials=1:3,CVqfolds=3,CVqrep=1,CVqStrt=TRUE,...)
{
  if (!is.matrix(data)) stop("'data' is not a matrix")
  if (!is.factor(grouping)) stop("'grouping' is not a factor")
  n <- nrow(data)
  if ( n != length(grouping)) stop("nrow(data) and length(grouping) are different")
  maxcol <- 2000   # Maximum dimensionality allowed for square matrices. Matrix square roots are ued when this limit is surpassed.
  if (is.numeric(q)) {
	if (q > maxq) stop("The number of factors exceeds its upper limit of",maxq," . This limit can be increased by changing the value of the 'maxq' argument, but the resulting classification rule may take too long to compute.")
	return(RFqlda(data,grouping,q,prior,CorrAp,VSelfunct,maxcol,nstarts,...))
  }
  else  {

  	TAl <- function(clrule,...) clrule

	if (q!="CVq") stop("Argument q must be numeric or equal to the string 'CVq'")
	if (max(CVqtrials) > maxq) stop("The maximum number of factors tp be tested exceeds the upper limit of",maxq," . This limit can be increased by changing the value of the 'maxq' argument, but the resulting classification rule may take too long to compute.")
  	nbqtrials <- length(CVqtrials)
  	Results <- vector("list",nbqtrials)
  	CVResults <- array(dim=nbqtrials)
  	totrep <- CVqrep*CVqfolds
  	gcds <- levels(grouping)
  	for (q in CVqtrials)  Results[[q]] <- RFqlda(data,grouping,q,prior,CorrAp,VSelfunct,maxcol,nstarts,...)
  	for (q in CVqtrials)  {   
		tmp <- DACrossVal(data,grouping,TrainAlg=TAl,kfold=CVqfolds,CVrep=CVqrep,Strfolds=CVqStrt,grpcodes=gcds,clrule=Results[[q]],...)
		errates <- apply(tmp[,,1]*tmp[,,2],2,sum)/n
		CVResults[q] <- mean(errates)
    	}
    	bestq <- which.min(CVResults)	
    	Results[[bestq]]  # return(Results[[bestq]])
    }
}

RFlda.data.frame <- function(data,...)
{
   res <- RFlda.default(as.matrix(data),...)
   res$call <- match.call()
   res
}

RFqlda <- function(data,grouping,q,prior,CorrAp,VSelfunct,maxcol,nstarts,...)
{

#  Compute the matrices Xw and Xdelta of within and  between group deviations

  if  ( !is.matrix(data) || !is.factor(grouping) || !is.numeric(q) ) stop("Arguments of wrong type")
  grp <- factor(grouping,levels<-sort(unique(grouping))) 
  nk <- table(grp)
  n <- sum(nk)
  p <- ncol(data)
  k <- nrow(nk)
  if (n != length(grouping) || n != nrow(data) ) stop("Argument dimensions do not match")

  if (is.character(VSelfunct)) {
	if (VSelfunct=="none") {
		m <- p
		vkpt <- 1:p
	}
	else stop("Invalid value for VSelfunct argument\n")
  }
  else  {
	if (!is.function(VSelfunct)) stop("Invalid value for VSelfunct argument\n")
        SelV <- VSelfunct(data,grouping,...)
	if (!is.list(SelV)) stop("Invalid value for VSelfunct argument\n")
	m <- SelV$nvkpt
	vkpt <- SelV$vkptInd
	data <- matrix(data[,vkpt],n,m)
  }  
  if (CorrAp)  {
  	sclres <- scalebygrps(data,grouping,k=k,nk=nk,n=n,p=m)
  	data <- sclres$Xscld
  }
  if (q > m) {
	warning("Number of assumed factors reduced from",q,"to the number of selected variables, that was equal to",m,"\n")
	q <- m
  }

  u  <- apply(data,2,mean)
  uk <- apply(data,2,grpmeans,grp=grp)
  Xw <- matrix(0.,n,m)
  Xdelta <- matrix(0.,k-1,m)
  Xavg <- matrix(0.,k-1,m)
  for (i in 2:k) {
	Xdelta[i-1,] <- uk[i,] - uk[1,]
	Xavg[i-1,] <- (uk[i,]+uk[1,])/2
  }
  for (i in 1:k)  Xw[grp==levels(grp)[i],] <- data[grp==levels(grp)[i],] - matrix(rep(uk[i,],nk[i]),nk[i],m,byrow=TRUE)

#  Compute the matrices Sig and SigFq of total and Factor-q approximation within variances and covariances

   if (m <= maxcol)  {
	Sig <- t(Xw) %*% Xw / (n - k)
	SigFq <- FrobSigAp(Sig,q,nstarts)
   }
   else  {
	SigFq <- FrobSigAp1(Xw/sqrt(n-k),min(n,m),q,nstarts)
  }

#  Compute the coefficients of the linear discriminant rule and return the results

   SigkptInv <- solve(SigFq)
   Coef <- t(LeftMult(SigkptInv,Xdelta)) 
   colnames(Coef) <- paste("LD",1:(k-1),sep="")
   if (prior[1] == "proportions") prior <- as.numeric(nk/n)
   cnst <- array(dim=k-1)
   for (grp in 2:k) cnst[grp-1] <-  matrix(Xavg[grp-1,],1,m) %*% Coef[,grp-1]
   if (CorrAp)
     result <- list(prior=prior,means=uk,coef=Coef/sclres$stdev,cnst=cnst,vkpt=vkpt,nvkpt=m,q=q,SigFq=SigFq,SigFqInv=SigkptInv,N=n,call=match.call())
   else result <- list(prior=prior,means=uk,coef=Coef,cnst=cnst,vkpt=vkpt,nvkpt=m,q=q,SigFq=SigFq,SigFqInv=SigkptInv,N=n,call=match.call())
   class(result) <- "RFlda"
   result    # return(result)
}

is.RFlda <- function(x)  inherits(x,"RFlda")

predict.RFlda <- function(object,newdata,prior=object$prior,grpcodes=NULL,...)
{
   assigntomax <- function(scr) return(grpcodes[which.max(scr)])

   if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
   n <- nrow(newdata)
   k <- ncol(object$coef) + 1
   if (is.null(grpcodes)) grpcodes <- 0:(k-1)
   cnst <- object$cnst + log(prior[1]/prior[-1])   # Adjust lda's treshold by the chosen (or estimated) prior probabilities
   if (k>2) {
 	if (is.null(object$vkpt)) scores <- cbind(rep(0.,n),newdata%*%object$coef - matrix(rep(cnst,n),n,k-1,byrow=TRUE))
 	else scores <- cbind(rep(0.,n),newdata[,object$vkpt]%*%object$coef - matrix(rep(cnst,n),n,k-1,byrow=TRUE))
  	res <- apply(scores,1,assigntomax)
  }
  else {
	if (is.null(object$vkpt)) scores <- newdata %*% matrix(object$coef,object$nvkpt,1)
	else scores <- newdata[,object$vkpt] %*% matrix(object$coef,object$nvkpt,1)
      	res <- rep(grpcodes[1],n)
	res[scores>drop(cnst)] <- grpcodes[2]
	scores <- cbind(rep(0.,n),scores)
   }
   escores <- exp(scores)
   posterior <- escores/apply(escores,1,sum)
   posterior[which(is.infinite(escores))] <- 1.
   colnames(posterior) <- grpcodes
   list(class=factor(res),posterior=posterior,x=scores[,-1])
}

print.RFlda <- function(x,...)
{
	cat("Call:\n") ; print(x$call)  
	cat("\nPrior probabilities of groups:\n") ; print(x$prior)
	cat("\nGroup means:\n") ; print(x$means)
	cat("\nCoefficients of linear discriminants:\n") ; print(x$coef)
	if (!(is.null(x$vkpt))) {
		cat("\nVariables kept in discriminant rule:\n") ; print(x$vkpt)
		cat("Number of variables kept:",x$nvkpt,"\n")
	}
}

coef.RFlda <- function(object,...)
{
	cat("\nCoefficients of linear discriminants:\n") 
	print(object$coef)
}



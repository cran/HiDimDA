### SigFqInv.R  (2011-06-13)
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

SigFqInv <- function(D,B,p,q,optres=NULL) solve(SigFq(D,B,p,q,optres))

is.SigFqInv <- function(x)  inherits(x,"SigFqInv")

as.matrix.SigFqInv <- function(x,...)  diag(x$D)-x$B%*%t(x$B)

print.SigFqInv <- function(x,...)
{
	xi <- solve(x)
	cat(paste("Precision matrix for a ",x$q,"-factor model\n",sep=""))
	cat("\nLoadings:\n") ; print(xi$B)
	cat("\nSpecific Variances:\n",xi$D,"\n")
}

"+.SigFqInv" <- function(x,a)
{	
   if (x$p==1)  {
   	if (is.SigFq(a)) return( matrix(x$D-x$B^2+a$D+a$B^2,1,1) )
   	if (is.SigFqInv(a)) return( matrix(x$D-x$B^2+a$D-a$B^2,1,1) )
   	return( matrix(x$D-x$B^2+a,1,1) )
   }
   if (is.SigFq(a)) return( diag(x$D)-x$B%*%t(x$B)+diag(a$D)+a$B%*%t(a$B) )
   if (is.SigFqInv(a)) return( diag(x$D)-x$B%*%t(x$B)+diag(a$D)-a$B%*%t(a$B) )
   diag(x$D)-x$B%*%t(x$B)+a   # return( diag(x$D)-x$B%*%t(x$B)+a )
}

"-.SigFqInv" <- function(x,a)
{	
   if (x$p==1)  {
   	if (is.SigFq(a)) return( matrix(x$D-x$B^2-a$D-a$B^2,1,1) )
   	if (is.SigFqInv(a)) return( matrix(x$D-x$B^2-a$D+a$B^2,1,1) )
   	return( matrix(x$D-x$B^2-a,1,1) )
   }
   if (is.SigFq(a)) return( diag(x$D)-x$B%*%t(x$B)-diag(a$D)-a$B%*%t(a$B) )
   if (is.SigFqInv(a)) return( diag(x$D)-x$B%*%t(x$B)-diag(a$D)+a$B%*%t(a$B) )
   diag(x$D)-x$B%*%t(x$B)-a # return( diag(x$D)-x$B%*%t(x$B)-a )
}

"*.SigFqInv" <- function(x,a)
{
   if (x$p==1)  {
   	if (is.SigFq(a)) return( matrix((x$D-x$B^2)*(a$D+a$B^2),1,1) )
   	if (is.SigFqInv(a)) return( matrix((x$D-x$B^2)*(a$D-a$B^2),1,1) )
   	if (length(a)>1) return( matrix((x$D-x$B^2)*a,1,1) )
   }
   if (is.SigFq(a)) return( (diag(x$D)-x$B%*%t(x$B))*(diag(a$D)+a$B%*%t(a$B)) )
   if (is.SigFqInv(a)) return( (diag(x$D)-x$B%*%t(x$B))*(diag(a$D)-a$B%*%t(a$B)) )
   if (length(a)>1) return( (diag(x$D)-x$B%*%t(x$B))*a )
   result <- list(p=x$p,q=x$q,B=sqrt(a)*x$B,D=a*x$D,optres=NULL)
   class(result) <- "SigFqInv"
   result  # return(result) 
}

"/.SigFqInv" <- function(x,a)
{
   if (x$p==1)  {
   	if (is.SigFq(a)) return( matrix((x$D-x$B^2)/(a$D+a$B^2),1,1) )
   	if (is.SigFqInv(a)) return( matrix((x$D-x$B^2)/(a$D-a$B^2),1,1) )
   	if (length(a)>1) return( matrix((x$D-x$B^2)/a,1,1) )
   }
   if (is.SigFq(a)) return( (diag(x$D)-x$B%*%t(x$B))/(diag(a$D)+a$B%*%t(a$B)) )
   if (is.SigFqInv(a)) return( (diag(x$D)-x$B%*%t(x$B))/(diag(a$D)-a$B%*%t(a$B)) )
   if (length(a)>1) return( (diag(x$D)-x$B%*%t(x$B))/a )
   result <- list(p=x$p,q=x$q,B=x$B/sqrt(a),D=x$D/a,optres=NULL)
   class(result) <- "SigFqInv"
   result  # return(result) 
}

LeftMult.SigFqInv <- function(x,a)
{
   if (is.SigFq(a) || is.SigFqInv(a)) {
	if (x$p>1) tmp1 <- diag(x$D*a$D)
	else tmp1 <- x$D*a$D
	tmp2 <- a$D * x$B 
	tmp3 <- x$D * a$B 
	tmp4 <- t(a$B) %*% x$B
	if (is.SigFq(a)) result <- tmp1 - tmp2%*%t(x$B) + a$B%*%t(tmp3) - a$B%*%tmp4%*%t(x$B) 
	if (is.SigFqInv(a)) result <- tmp1 - tmp2%*%t(x$B) - a$B%*%t(tmp3) + a$B%*%tmp4%*%t(x$B)
	return(result)
   }
   tmp <- a %*% x$B
   result <- t(t(a)*x$D) - tmp%*%t(x$B)
   if (!is.matrix(a)) dim(result) <- c(1,x$p)
   result  # return(result)
}

RightMult.SigFqInv <- function(x,a)
{
   if (is.SigFq(a) || is.SigFqInv(a)) {
	if (x$p>1) tmp1 <- diag(x$D*a$D)
	else tmp1 <- x$D*a$D
	tmp2 <- x$D * a$B 
	tmp3 <- a$D * x$B 
	tmp4 <- t(x$B) %*% a$B
	if (is.SigFq(a)) result <- tmp1 + tmp2%*%t(a$B) - x$B%*%t(tmp3) - x$B%*%tmp4%*%t(a$B)  
	if (is.SigFqInv(a)) result <- tmp1 - tmp2%*%t(a$B) - x$B%*%t(tmp3) + x$B%*%tmp4%*%t(a$B)  
	return(result)
   }
   tmp <- t(x$B) %*% a
   result <- a*x$D - x$B%*%tmp
   if (!is.matrix(a)) dim(result) <- c(x$p,1)
   result  #  return(result)
}


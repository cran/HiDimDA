### miscelanea.R  (2011-04-20)
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

grpmeans <- function(x,grp) tapply(x,INDEX=grp,mean)
grpvar <- function(x,grp) tapply(x,INDEX=grp,var)
l2vnorm <- function(v)  sum(v^2)

scalebygrps <- function(x,grouping,k=2,nk=NULL,n=NULL,p=NULL) 
{
	if (is.null(nk))  nk <- as.vector(table(grouping))
	if (is.null(n))  n <- sum(nk)
	if (is.null(p)) p <- ncol(x)
	vark <-  apply(x,2,grpvar,grp=grouping)
	globals <- sqrt(apply(matrix(rep(nk-1,p),k,p,byrow=FALSE)*vark,2,sum)/(n-k))
	list(Xscld=x/matrix(globals,n,p,byrow=TRUE),stdev=globals)
}


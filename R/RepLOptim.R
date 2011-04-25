### RepLOptim.R  (2011-04-20)
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

RepLOptim <- function(parmean,parsd,fr,nrep,niter,gr=NULL,inphess=NULL,method="nlminb",lower=NULL,upper=NULL,rethess=FALSE,
allrep=NULL,neval=NULL,objbnd=Inf,parmstder=FALSE,tol=1E-8,...) 

#   RepLOptim -- Repeated local optimization 

#   Tries to minimize a function calling local optimizers several times from different random starting points
#   generated from multivariate normal distributions of independent variates. The standard deviations of the 
#   generating distributions are kept fixed, but their means are updated as better candidates for the global 
#   minimun are discovered

#   Arguments:

#   parmean --   vector of means for the parameter distribution generating starting points for the
#                local optimizer. Also used as stating point of the first call to the optimizer
#   parsd   --   vector of standard deviations for the parameter distribution generating starting points 
#                for the local optimizer. Also used as stating point of the first call to the optimizer
#   fr      --   the function to be minimized. If method is neither "nlminb" or "L-BFGS-B", fr should
#                accept a lbound and an ubound arguments for the parameters bounds, and should enforce
#                these bounds before calling the local optimization routine
#   nrep    --   number of replications (different calls to the local optimizer leading to valid solutions) 
#                to be performed
#   niter   --   maximum number of iterations performed in each call to the local optimizer.
#   gr      --   A function to return the gradient for the "nlminb", "BFGS", '"CG"'and '"L-BFGS-B"' methods.  
#                If it is 'NULL', a finite-difference approximation will be used. For the '"SANN"' method 
#                it specifies a function to generate a new candidate point.  If it is 'NULL' a default Gaussian
#                Markov kernel is used.
#   inphes  --   A function that takes the same arguments as 'objective' and evaluates the hessian of 'objective' at its 
#                first argument. Must return a square matrix of order 'length(start)' with the different hessian elements
#                in its lower triangle (the upper triangle is ignored). At present time only available for the nlminb method. 
#   method  --   local optimizer to be employed. Current alternatives are:
#                "nlminb" for the nlminb port routine, "nlm" for the nlm function and
#                "Nelder-Mead", "L-BFGS-B",  "CG", "L-BFGS-B" and "SANN" for the corresponding methods of
#                the optim function 
#   lower   --   vector of parameter lower bounds. Set to -Inf (no bounds) by default
#   upper   --   vector of parameter upper bounds. Set to Inf (no bounds) by default
#   rethess  --  Bolean flag indicating wether a numerically evaluated hessian matrix at the optimum 
#                should be computed and returned. Not available for the nlminb method.
#   allrep  --   maximum number of replications (including those leading to non-valid solutions)
#                performed. By default equals ten times nrep. Ignored when objbnd is set to Inf 
#   neval    --  maximum number of function evaluations (nlminb method only) performed in each call to 
#                the nlminb optimizer
#   objbnd   --  Upper bound for the objective. Only solutions leading to objective values below objbnd 
#                are considered as valid
#   parmstder --  Bolean flag indicating wether parameter assymptotic standard errors based on the 
#                inverse hessian approximation to the Fisher information matrix should be computed and returned.
#                Only available if hessian is set to true and if a local miminum with a positive-definite hessian 
#                was indeed found. For dificult problems this requirement may fail if nrep and niter (and maybe 
#                neval) are not large enough.
#                Note that the consistency of these estimates can only be assured when fr is based on a true 
#                likelihood (as oposed to a quasi-likelihood) function
#   tol      --  Numerical tolerance used to ensure that the hessian is non-singular. If the last eigenvalue 
#                of the hessian is positive but below tol, the hessian is considered to be semi-definite and the
#                the parameter assymptotic standard errors are not computed.
#   ...      --  Further arguments to be passed to fr, gr or to the local optimization routine selected

#   Value:  A list with the following components
#
#         par            --   the best result found for the parameter vector
#         val            --   the best value (minimum) found for the function fr
#         vallist        --   a vector with the best values found for each starting point
#         iterations     --   number the iterations performed by the local optimizer in the call
#                             that generated the best result
#         counts         --   number of times the function fr was evaluated in the call that 
#                             generated the result returned
#         convergence    --   code with the convergence status returned by the local optimizer
#         message        --   message generated by the local optimizer
#         hessian        --   Numerically evaluated hessian of fr at the result returned 
#                             Only returned when the parameter hessian is set to true 
#         hessegval      --   Eigenvalues of the hessian matrix. Used to confirm if a local minimum was indeed found.
#                             Only returned when the parameter hessian is set to true 
#         stderrors      --   Assymptotic standard deviations of the parameters based on Fisher Information 
#                             matrix. Only returned when the parse parameter is set to true and the hessian is indeed
#                             positive definite.

{
   npar <- length(parmean)
    values <- array(dim=nrep)
    if (is.null(lower)) lower <- rep(-Inf,npar)
    if (is.null(upper)) upper <- rep(Inf,npar)
    if (is.finite(objbnd))  { if (is.null(allrep)) allrep <- 10*nrep }
    else allrep <- nrep

    bestres <- NULL	
    bestval <- Inf
    bestpar <- parmean
    initpar <- parmean
    cnt <- 0
    for (i in 1:nrep)  {
       if (cnt > allrep) break
       value <- Inf
       while (value >= objbnd && cnt < allrep)
      {
       	 if (method == "nlminb")
	    tmpres <- nlminb(start=initpar,fr,gradient=gr,hessian=inphess,lower=lower,upper=upper,control=list(iter.max=niter,eval.max=neval),...)
       	 else if (method == "nlm") 
	    tmpres <- nlm(fr,p=initpar,lbound=lower,ubound=upper,iterlim=niter,...)
       	 else if (method == "L-BFGS-B")
	    tmpres <- optim(initpar,fr,gr=gr,method=method,lower=lower,upper=upper,control=list(maxit=niter),hessian=rethess,...)
       	 else  tmpres <-
            optim(initpar,fr,gr=gr,method=method,control=list(maxit=niter),lbound=lower,ubound=upper,hessian=rethess,...)
       	 if (method == "nlminb") value <- tmpres$objective
       	 else if (method == "nlm") value <- tmpres$minimum
       	 else value<- tmpres$value
         if (is.na(value)) value <- objbnd
         cnt <- cnt+1
	 u <- runif(n=npar)     # generate npar uniform random numbers
       	 initpar <- qnorm(u,mean=bestpar,sd=parsd) #  generate new parameters from a normal distribution
	 lbndind<- initpar < lower   #  identify indices of parameters that fell below their lower bounds
	 ubndind <- initpar > upper   #  identify indices of parameters that fell above their upper bounds
	 initpar[lbndind] <- lower[lbndind] + u[lbndind] * (bestpar[lbndind]-lower[lbndind]) # and correct those parameters
	 initpar[ubndind] <- upper[ubndind] - u[ubndind] * (upper[ubndind]-bestpar[ubndind]) 
       } 
      values[i] <- value
       if (!is.na(value)) if (value < bestval)  {
          bestval <- value
          if (method != "nlm") bestpar <- tmpres$par
          else bestpar <- tmpres$estimate
          bestres <- tmpres
       }
    }
    if (method == "nlminb")  {
        iterations <- bestres$iterations
        if (method != "nlm") counts <- bestres$evaluations
        else counts <- NULL
        hess <- NULL
        egval <- NULL
        parstd <- NULL
    } 
    else  {
    	iterations <- NULL
    	counts <- bestres$counts
        if (rethess==TRUE) {
           hess <- bestres$hessian
           egval <- eigen(hess,sym=TRUE,only.values=TRUE)$values
       	   if (parmstder==TRUE)
              if (egval[npar] < tol) parstd <- "Not computed because the hessian is not positive definite"
              else parstd <- sqrt(diag(solve(hess)))
           else parstd <- NULL
        }
        else {
	    hess <- NULL
            egval <- NULL
            parstd <- NULL
       }
    } 
    if (!is.null(bestres))
    	return(list(par=bestpar,val=bestval,iterations=iterations,vallist=values,counts=counts,convergence=bestres$convergence,message=bestres$message,hessian=hess,hessegval=egval,stderrors=parstd))
    else
    	return(list(par=NULL,val=Inf,iterations=NULL,vallist=NULL,counts=NULL,convergence=NULL,message="RepLOptim was unable to find any valid solution",hessian=NULL,hessegval=NULL,stderrors=NULL))
}

rghtsngv <- function(x, nv=min(nrow(x),ncol(x)), maxsqmatdim=10000)
{ 
	if (!is.matrix(x)) x <- as.matrix(x)
 	n <- nrow(x)
	p <- ncol(x)
	minnp <- min(n,p)	
	nv <- min(nv,minnp)
	d <- numeric(minnp)
	vt <- numeric(minnp*p)
	u <- numeric(n*p)
	retcode <- 9999

	iwork <- integer(8*minnp)
	jobz <- 'S'
	work<- numeric(1)
       	lwork <- max( 5*minnp+2*max(n,p),
		.Fortran("dgesdd",
			as.character(jobz),as.integer(n),as.integer(p),as.double(x),as.integer(n),
			as.double(d),as.double(u),as.integer(n),as.double(vt),as.integer(minnp),
			as.double(work),as.integer(-1),as.integer(iwork),as.integer(retcode))[[11]]
		)
	work <- numeric(lwork)	
       	Fortout <- .Fortran("dgesdd",
			as.character(jobz),as.integer(n),as.integer(p),as.double(x),as.integer(n),
			as.double(d),as.double(u),as.integer(n),as.double(vt),as.integer(minnp),
			as.double(work),as.integer(lwork),as.integer(iwork),as.integer(retcode)
		)
	if (Fortout[[14]]==0) return(list(d=drop(Fortout[[6]]),u=NULL,v=t(matrix(Fortout[[9]],minnp,p))[,1:nv]))
	else {
		dgesddrtc <- Fortout[[14]]  
      		jobu <- 'N'
      		jobvt <- 'S'
       		Fortout <- .Fortran("dgesvd",
			as.character(jobu),as.character(jobvt),as.integer(n),as.integer(p),
			as.double(x),as.integer(n),as.double(d),as.double(u),as.integer(1),
			as.double(vt),as.integer(minnp),as.double(work),as.integer(lwork),
			as.integer(retcode)
		)
		if (Fortout[[14]]==0) return(list(d=drop(Fortout[[7]]),u=NULL,v=t(matrix(Fortout[[10]],minnp,p))[,1:nv]))
	   	else  {
			dgesvdrtc <- Fortout[[14]]  
			if (p>maxsqmatdim) 
			   stop(paste("Computation of right singular vectors failed for a",n,"by",p," matrix,\n  with return codes",dgesddrtc,"and",dgesvdrtc,"for the Lapack dgesdd and dgesvd routines.\n"))
			else {
			   SpecDec <- eigen(t(x)%*%x,symmetric=TRUE)
			   return(list(d=ifelse(SpecDec$values>0.,sqrt(abs(SpecDec$values)),0.)[1:minnp],u=NULL,v=SpecDec$vectors[,1:nv]))
			}
	   	}
	}
	
}


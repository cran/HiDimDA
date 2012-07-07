**********************************************************************
       SUBROUTINE lworkgesdd(x,d,vt,n,p,minnp,u,work,iwork,info)
**********************************************************************
  
       integer n,p,minnp,info,lwork,iwork(*)
       double precision x(n,p),d(minnp),vt(minnp,p)
       double precision u(n,minnp),work(1)
       character*1 jobz

       external dgesdd

       jobz = 'S'
       call dgesdd(jobz,n,p,x,n,d,u,n,vt,minnp,
     +  work,-1,iwork,info)

	return
	end   

**********************************************************************
       SUBROUTINE rsgvdgesdd(x,d,vt,n,p,minnp,u,work,lwork,iwork,info)
**********************************************************************
  
       integer n,p,minnp,info,lwork,iwork(*)
       double precision x(n,p),d(minnp),vt(minnp,p)
       double precision u(n,minnp),work(lwork)
       character*1 jobz

       external dgesdd

       jobz = 'S'
       call dgesdd(jobz,n,p,x,n,d,u,n,vt,minnp,
     +  work,lwork,iwork,info)

	return
	end   

**********************************************************************
       SUBROUTINE rsgvdgesvd(x,d,vt,n,p,minnp,work,lwork,info)
**********************************************************************
  
       integer n,p,minnp,lwork,info
       double precision x(n,p),d(minnp),vt(minnp,p),u,work(lwork)
       character*1 jobu,jobvt

       external dgesvd

      jobu = 'N'
      jobvt = 'S'
      call dgesvd(jobu,jobvt,n,p,x,n,d,u,1,vt,minnp,
     +  work,lwork,info)

	return
	end   

**********************************************************************
      SUBROUTINE rsgvdsyev(x,d,vt,x2,d2,n,p,minnp,work,lwork,info)
**********************************************************************
  
       integer i,j,k,n,p,info,lwork
       double precision x(n,p),d(minnp),vt(minnp,p)
       double precision x2(p,p),d2(p),d2i,work(lwork)
       character*1 jobz,uplo

       external dsyev

       do i=1,p
	 do j=1,i
	    x2(i,j) = x(1,i)*x(1,j)
	    if (n.ge.1) then	
	    	do k=2,n  
		   x2(i,j) = x2(i,j) + x(k,i)*x(k,j)
		end do
	    end if
	 end do
       end do

       jobz = 'V'
       uplo = 'L'
       call dsyev(jobz,uplo,p,x2,p,d2,work,lwork,info)
	
       if (info.eq.0) then
	   do i=1,minnp
		d2i = d2(p-i+1)
		if (d2i.le.0.) then 
		   d(i) = 0.
		else
		   d(i) = sqrt(d2i)
		end if	
		do j=1,p
		  vt(i,j) = x2(j,p-i+1)
		 end do
	    end do
         end if	

       return
       end   

**********************************************************************


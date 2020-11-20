!-----------------------------------------------------------------------
! Numerical Recipes subroutines for non-linear least squares fit.
!
! Contents:
! MRQMIN
! MRQCOF
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Levenberg-Marquardt method.
!
! Attempting to reduce the value chi^2 of a fit between a set of data
! points x(ndata), y(ndata) with individual standard deviations
! sig(ndata), and a non-linear function dependent on ma coefficients
! a(ma).
!
! The input array ia(ma) indicates by non-zero entries those components
! of a(ma) that should be fitted for, and by zero entries those that
! should be held fixed at their input values. The program returns
! current best-fit values for the parameters a(ma), and chi^2 = chisq.
!
! The arrays covar(nca,nca), alpha(nca,nca) with physical dimension nca
! (>= the number of fitted parameters) are used as working space during
! most iterations.
!
! Supply a subroutine FUNCS(x,a,yfit,dyda,ma) that evaluates the fitting
! function yfit, and its derivatives dyda with respect to the fitting
! parameters a at x.
!
! On the first call provide an initial guess for the parameters a, and
! set alamda < 0 for initialization (which then sets alamda=.001). If a
! step succeeds chisq becomes smaller and alamda decreases by a factor
! of 10. If a step fails alamda grows by a factor of 10. You must call
! this routine repeatedly until convergence is achieved. Then, make one
! final call with alamda=0, so that covar(ma,ma) returns the covariance
! matrix, and alpha the curvature matrix. (Parameters held fixed will
! return zero covariances.)
!
! Uses:
! COVSRT
! GAUSSJ
! MRQCOF
!-----------------------------------------------------------------------
	subroutine MRQMIN(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,
     &	chisq,FUNCS,alamda)

	implicit none
	integer ma,nca,ndata,ia(ma),MMAX
	real*8 alamda,chisq,a(ma)
	real*8 alpha(nca,nca),covar(nca,nca)
	real*8 sig(ndata),x(ndata),y(ndata)
	external FUNCS

! Set the largest number of fit parameters.
	parameter(MMAX=20)

! Internally used variables.
	integer j,k,l,mfit
	real*8 ochisq,atry(MMAX),beta(MMAX),da(MMAX)
	save ochisq,atry,beta,da,mfit

	if(alamda.lt.0d0)then
	   mfit=0
	   do j=1,ma
	      if(ia(j).ne.0) mfit=mfit+1
	   enddo
	   alamda=0.001d0
	   call MRQCOF(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,FUNCS)
	   ochisq=chisq
	   do j=1,ma
	      atry(j)=a(j)
	   enddo
	endif

	do j=1,mfit
	   do k=1,mfit
	      covar(j,k)=alpha(j,k)
	   enddo
	   covar(j,j)=alpha(j,j)*(1d0+alamda)
	   da(j)=beta(j)
	enddo

	call GAUSSJ(covar,mfit,nca,da,1,1)

	if(alamda.eq.0d0)then
	   call COVSRT(covar,nca,ma,ia,mfit)
	   call COVSRT(alpha,nca,ma,ia,mfit)
	   return
	endif

	j=0
	do l=1,ma
	   if(ia(l).ne.0)then
	      j=j+1
	      atry(l)=a(l)+da(j)
	   endif
	enddo

	call MRQCOF(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,FUNCS)

	if(chisq.lt.ochisq)then
	   alamda=0.1d0*alamda
	   ochisq=chisq
	   do j=1,mfit
	      do k=1,mfit
		 alpha(j,k)=covar(j,k)
	      enddo
	      beta(j)=da(j)
	   enddo
	   do l=1,ma
	      a(l)=atry(l)
	   enddo
	else
	   alamda=10d0*alamda
	   chisq=ochisq
	endif

	return
	end
!-----------------------------------------------------------------------
! Subroutine used by MRQMIN to evaluate the linearized fitting matrix
! alpha, and vector beta as in (15.5.8), and calculate chi^2.
!-----------------------------------------------------------------------
	subroutine MRQCOF(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,
     &	FUNCS)

	implicit none
	integer ma,nalp,ndata,ia(ma),MMAX
	real*8 chisq,a(ma),alpha(nalp,nalp),beta(ma)
	real*8 sig(ndata),x(ndata),y(ndata)
	external FUNCS

! Set the largest number of fit parameters.
	parameter(MMAX=20)

! Internally used variables.
	integer mfit,i,j,k,l,m
	real*8 dy,sig2i,wt,ymod,dyda(MMAX)

	mfit=0
	do j=1,ma
	   if(ia(j).ne.0) mfit=mfit+1
	enddo
	do j=1,mfit
	   do k=1,j
	      alpha(j,k)=0d0
	   enddo
	   beta(j)=0d0
	enddo

	chisq=0d0
	do i=1,ndata
	   call FUNCS(x(i),a,ymod,dyda,ma)
	   sig2i=1d0/(sig(i)*sig(i))
	   dy=y(i)-ymod
	   j=0
	   do l=1,ma
	      if(ia(l).ne.0)then
		 j=j+1
		 wt=dyda(l)*sig2i
		 k=0
		 do m=1,l
		    if(ia(m).ne.0)then
		       k=k+1
		       alpha(j,k)=alpha(j,k)+wt*dyda(m)
		    endif
		 enddo
		 beta(j)=beta(j)+dy*wt
	      endif
	   enddo
	   chisq=chisq+dy*dy*sig2i
	enddo

	do j=2,mfit
	   do k=1,j-1
	      alpha(k,j)=alpha(j,k)
	   enddo
	enddo

	return
	end

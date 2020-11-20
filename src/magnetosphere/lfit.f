!-----------------------------------------------------------------------
! Numerical Recipes subroutine for general linear least squares fit.
!
! Contents:
! LFIT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Given a set of data points x(1:ndat), y(1:ndat) with individual
! standard deviations sig(1:ndat), use chi^2 minimization to fit for
! some or all of the coefficients a(1:ma) of a function that depends
! linearly on a, y = a_i * afunc_i(x). The input array ia(1:ma)
! indicates by nonzero entries those components of a that should be
! fitted for, and by zero entries those components that should be held
! fixed at their input values. The program returns values for a(1:ma),
! chi^2 = chisq, and the covariance matrix covar(1:ma,1:ma). (Parameters
! held fixed will return zero covariances.) npc is the physical
! dimension of covar(npc,npc) in the calling routine. The user supplies
! a subroutine FUNCS(x,afunc,ma) that returns the ma basis functions
! evaluated at x in the array afunc.
!
! Changes:
! - Converted to double precision.
!
! Uses:
! COVSRT
! GAUSSJ
!-----------------------------------------------------------------------
	subroutine LFIT(x,y,sig,ndat,a,ia,ma,covar,npc,chisq,FUNCS)

	implicit none
	integer ma,ia(ma),npc,ndat,MMAX
	double precision chisq,a(ma),covar(npc,npc)
	double precision sig(ndat),x(ndat),y(ndat)
	external FUNCS

! Set the maximum number of coefficients ma.
	parameter(MMAX=50)

	integer i,j,k,l,m,mfit
	double precision sig2i,sum,wt,ym,afunc(MMAX),beta(MMAX)

	mfit=0
	do j=1,ma
	   if(ia(j).ne.0) mfit=mfit+1
	enddo
	if(mfit.eq.0) stop "LFIT: no parameters to be fitted!"

! Initialize the (symmetric) matrix.
	do j=1,mfit
	   do k=1,mfit
	      covar(j,k)=0d0
	   enddo
	   beta(j)=0d0
	enddo

! Loop over data to accumulate coefficients of the normal equations.
	do i=1,ndat
	   call FUNCS(x(i),afunc,ma)
	   ym=y(i)
! Subtract off dependences on known pieces of the fitting function.
	   if(mfit.lt.ma)then
	      do j=1,ma
		 if(ia(j).eq.0) ym=ym-a(j)*afunc(j)
	      enddo
	   endif
	   sig2i=1d0/sig(i)**2
	   j=0
	   do l=1,ma
	      if(ia(l).ne.0)then
		 j=j+1
		 wt=afunc(l)*sig2i
		 k=0
		 do m=1,l
		    if(ia(m).ne.0)then
		       k=k+1
		       covar(j,k)=covar(j,k)+wt*afunc(m)
		    endif
		 enddo
		 beta(j)=beta(j)+ym*wt
	      endif
	   enddo
	enddo

! Fill in above the diagonal from symmetry.
	do j=2,mfit
	   do k=1,j-1
	      covar(k,j)=covar(j,k)
	   enddo
	enddo

! Matrix solution.
	call GAUSSJ(covar,mfit,npc,beta,1,1)

	j=0
	do l=1,ma
! Partition solution to appropriate coefficients a.
	   if(ia(l).ne.0)then
	      j=j+1
	      a(l)=beta(j)
	   endif
	enddo

! Evaluate χ 2 of the fit.
	chisq=0d0
	do i=1,ndat
	   call FUNCS(x(i),afunc,ma)
	   sum=0d0
	   do j=1,ma
	      sum=sum+a(j)*afunc(j)
	   enddo
	   chisq=chisq+((y(i)-sum)/sig(i))**2
	enddo

! Sort covariance matrix to true order of fitting coefficients.
	call COVSRT(covar,npc,ma,ia,mfit)

	return
	end

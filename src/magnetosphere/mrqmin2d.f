!-----------------------------------------------------------------------
! Numerical Recipes subroutines for non-linear least squares fit.
!
! Modified versions of MRQMIN and MRQCOF for functions of two variables
! (cf. MRQMIN and MRQCOF).
!
! Contents:
! MRQMIN2D
! MRQCOF2D 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Modified version of MRQMIN for functions of two variables.
!
! Changes with respect to the original version of MRQMIN:
! - The subroutine has been updated for functions with up to two
! variables, y(x1,x2), instead of just one, y(x).
! - A new flag has been added, which controls whether FUNCS returns the
! function or its derivative.
!
! Uses:
! COVSRT
! GAUSSJ
! MRQCOF2D
!-----------------------------------------------------------------------
	subroutine MRQMIN2D(ii,x1,x2,y,sig,ndata,a,ia,ma,covar,alpha,
     &	nca,chisq,FUNCS,alamda)

	implicit none
	integer ma,nca,ndata,ii(ndata),ia(ma),MMAX
	double precision alamda,chisq,a(ma),alpha(nca,nca)
	double precision covar(nca,nca),sig(ndata)
	double precision x1(ndata),x2(ndata),y(ndata)
	external FUNCS

! Set the largest number of fit parameters.
	parameter(MMAX=50)

! Internally used variables.
	integer j,k,l,mfit
	double precision ochisq,atry(MMAX),beta(MMAX),da(MMAX)
	save ochisq,atry,beta,da,mfit

	if(alamda.lt.0d0)then
	   mfit=0
	   do j=1,ma
	      if(ia(j).ne.0) mfit=mfit+1
	   enddo
	   alamda=0.001d0
	   call MRQCOF2D(ii,x1,x2,y,sig,ndata,a,ia,ma,alpha,beta,nca,
     &	   chisq,FUNCS)
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

	call MRQCOF2D(ii,x1,x2,y,sig,ndata,atry,ia,ma,covar,da,nca,
     &	chisq,FUNCS)

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
! Modified version of MRQCOF for functions of two variables.
!
! Changes:
! - The subroutine has been updated for functions with up to two
! variables, y(x1,x2), instead of just one, y(x).
!
! Uses:
! FUNCS
!-----------------------------------------------------------------------
	subroutine MRQCOF2D(ii,x1,x2,y,sig,ndata,a,ia,ma,alpha,beta,
     &	nalp,chisq,FUNCS)

	implicit none
	integer ma,nalp,ndata,ii(ndata),ia(ma),MMAX
	double precision chisq,a(ma),alpha(nalp,nalp),beta(ma)
	double precision sig(ndata),x1(ndata),x2(ndata),y(ndata)
	double precision dydx1,dydx2
	external FUNCS

! Set the largest number of fit parameters.
	parameter(MMAX=50)

	integer mfit,i,j,k,l,m
	double precision dy,sig2i,wt,ymod,dyda(MMAX)

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
	   call FUNCS(ii(i),x1(i),x2(i),a,ymod,dyda,dydx1,dydx2,ma)
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

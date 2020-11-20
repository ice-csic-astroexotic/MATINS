!-----------------------------------------------------------------------
! Numerical Recipes subroutine.
!
! Contents:
! COVSRT
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Expand in storage the covariance matrix covar, so as to take into
! account parameters that are being held fixed. (For the latter, return
! zero covariances.)
!-----------------------------------------------------------------------
	subroutine COVSRT(covar,npc,ma,ia,mfit)

	implicit none
	integer ma,mfit,npc,ia(ma)
	double precision covar(npc,npc)

	integer i,j,k
	double precision swap

	do i=mfit+1,ma
	   do j=1,i
	      covar(i,j)=0d0
	      covar(j,i)=0d0
	   enddo
	enddo

	k=mfit
	do j=ma,1,-1
	   if(ia(j).ne.0)then
	      do i=1,ma
		 swap=covar(i,k)
		 covar(i,k)=covar(i,j)
		 covar(i,j)=swap
	      enddo
	      do i=1,ma
		 swap=covar(k,i)
		 covar(k,i)=covar(j,i)
		 covar(j,i)=swap
	      enddo
	      k=k-1
	   endif
	enddo

	return
	end

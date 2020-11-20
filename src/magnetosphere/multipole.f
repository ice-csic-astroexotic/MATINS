!-----------------------------------------------------------------------
! Multipole expansion at some radius.
!-----------------------------------------------------------------------
	subroutine MULTIPOLE(Pm,cth,dpl,lmax,nz,al)

	implicit none

! Input variables.
	integer lmax,nz
	real*8 Pm(nz),cth(nz),dpl(lmax,nz)

! Output variables.
	real*8 al(lmax)

! Internally used variables.
	integer l,n
	real*8 factor

	do l=1,lmax
	   al(l)=0d0
	   factor=dble(2*l+1)/dble(2*l*(l+1))
	   do n=2,nz
! Trapezoidal integration. (From mu=1 to mu=-1.)
	      al(l)=al(l)+factor
     &		 *(Pm(n-1)*dpl(l,n-1)+Pm(n)*dpl(l,n))
     &		 *(cth(n-1)-cth(n))/2d0
	   enddo
	enddo

	return
	end

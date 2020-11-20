!-----------------------------------------------------------------------
! Calculate a velocity field proportional to the Lorentz force in a
! transition layer between the surface and the magnetosphere.
! The surface is at rb(2*lmax) and the magnetosphere starts at rb(np).
! (See dim2.h for definitions of the dimensions.)
!
! Called by BEVOL.
!-----------------------------------------------------------------------
	subroutine SURFACE_TRANSITION
  use grid, only : br, bth, bphi, jr, jth, jphi

	implicit none
	include '../decl/dim2.h'

! Input: needs br, bth, bphi, jr, jth, jphi.
! Output: returns vr_tran, vth_tran, vphi_tran.

! Internally used variables.
	integer i,j,k,l
	real*8 jrave,jthave

! Lorentz force.
	real*8, dimension (0:nang+1,0:np+2) :: fr,fth,fphi
      real*8, dimension (0:nang+1,0:np+2) :: vr_amb,vth_amb,vphi_amb

! Amplification factor.
	real*8 amp
	parameter(amp=0.1d0)

!-----------------------------------------------------------------------
! Lorentz force density.
!-----------------------------------------------------------------------
	fr=0d0
	fth=0d0
	fphi=0d0

	do i=1,nang
	do j=1,np
	      fr(i,j)=jr(i,j)*bphi(i,j)-jphi(i,j)*bth(i,j)
	      fth(i,j)=jphi(i,j)*br(i,j)-jr(i,j)*bphi(i,j)
	      fphi(i,j)=jr(i,j)*bth(i,j)-jth(i,j)*br(i,j)
	enddo
	enddo

!-----------------------------------------------------------------------
! Velocity.
!-----------------------------------------------------------------------
! Amplification factor.
	vr_tran=amp*fr
	vth_tran=amp*fth
	vphi_tran=amp*fphi

!-----------------------------------------------------------------------
! Manually set values at the boundaries.
!-----------------------------------------------------------------------
! Along the north axis.
	vth_tran(1,:)=0d0
	vphi_tran(1,:)=0d0

! Along the south axis.
	vth_tran(nang,:)=0d0
	vphi_tran(nang,:)=0d0

	return
	end

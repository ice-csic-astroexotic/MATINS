!-----------------------------------------------------------------------
! Declarations.
!-----------------------------------------------------------------------
! Constants added here for legacy purposes.
      use grid, only: NANG, NP
! Magnetic field.
	real*8, dimension (0:nang+1,0:np+2) :: br,bth,bphi,aphi,bm
	common /B_field/ br,bth,bphi,aphi,bm

!-----------------------------------------------------------------------
! Current and electric field.
!-----------------------------------------------------------------------
	real*8, dimension (0:nang+1,0:np+2) :: jr,jth,jphi
	real*8, dimension (0:nang+1,0:np+2) :: er,eth,ephi
	common /EJ_fields/ jr,jth,jphi,er,eth,ephi
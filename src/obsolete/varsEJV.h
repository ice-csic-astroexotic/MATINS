!-----------------------------------------------------------------------
! Current and electric field.
!-----------------------------------------------------------------------
	real*8, dimension (0:nang+1,0:np+2) :: jr,jth,jphi
	real*8, dimension (0:nang+1,0:np+2) :: er,eth,ephi
	common /EJ_fields/ jr,jth,jphi,er,eth,ephi



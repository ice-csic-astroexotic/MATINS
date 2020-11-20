!-----------------------------------------------------------------------
! Declarations for radial profiles of physical variables.
! USED BY:
! input.f where they are read and defined
! ambipolar.f, emissivity.f, cvf.f, difdrv.f
!-----------------------------------------------------------------------
	real*8, dimension (np+2) :: rho,xh,ye,yn,yp,aa,zz
	common /ns/ rho,xh,ye,yn,yp,aa,zz

	real*8, dimension (np+2) :: tcn,tcp,tccru 
	common /tcrit/ tccru,tcn,tcp  

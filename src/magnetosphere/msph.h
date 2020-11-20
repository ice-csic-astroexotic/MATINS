!-----------------------------------------------------------------------
! Magnetospheric grid definitions.
!-----------------------------------------------------------------------
! Radial grid in the magnetosphere.
! (nrad is defined in dim2.h.)
	real*8, dimension (nrad) :: r_msph
	common /magnetosphere_grid/ r_msph

! Geometric factor for an uneven grid.
	real*8 geo
	parameter (geo=1.05d0)

! Radial grid boundaries.
	real*8 rin,rout
	parameter (rin=1d0,rout=10d0)

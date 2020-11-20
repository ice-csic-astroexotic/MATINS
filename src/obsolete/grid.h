!-----------------------------------------------------------------------
! Declarations related to the grid. (Defined in BGRID.)
!-----------------------------------------------------------------------
! Radial and angular grid.
	real*8, dimension (0:np+2) :: rb
	real*8, dimension (0:nang+1) :: theta,sth,cth
	common /staggered_grid/ rb,theta,sth,cth

! Length, area and volume elements.
	real*8, dimension (0:np+2) :: lr,lth
	real*8, dimension (0:nang+1,0:np+2) :: lphi
	real*8, dimension (0:np+2) :: areaphi
	real*8, dimension (0:nang+1,0:np+2) :: arear,areath,vol
	common /grid_elements/ lr,lth,lphi,arear,areath,areaphi,vol

! Relativistic factors.
	real*8, dimension (0:np+2) :: belam,benu
	common /relativistic_factors/ belam,benu

! Relativistic corrections for the vacuum boundary conditions.
	real*8, dimension (0:nleg) :: frel
	common /relativistic_corrections/ frel

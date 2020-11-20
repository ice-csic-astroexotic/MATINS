!-----------------------------------------------------------------------
! Declarations.
!
! Thermopower.
! Units of c*Q*gradT, with Q expressed in units of k_b/e,
! velocity in km/Myr, electric field in 1e12 G*km/Myr, gradT in 1e8 K/km
!-----------------------------------------------------------------------
	real*8, dimension (0:nang+1,0:np+2) :: qpar,qperp,qhall
	real*8, dimension (0:nang+1,0:np+2) :: gradtr,gradtth
	real*8, dimension (0:nang+1,0:np+2) :: vrte,vthte,vphite
	real*8 qfac
	parameter(qfac=2.72d3)

	common /thermopower/ qpar,qperp,qhall,gradtr,gradtth
	common /thermovel/ vrte,vthte,vphite

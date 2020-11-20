!-----------------------------------------------------------------------
! Obsolete declarations.
!
! Used by:
! input.f (Old version.)
! bfluxtube.f
!-----------------------------------------------------------------------
	integer ia,icore,imic
	common /vel_index/ ia,icore,imic

	real*8, dimension (0:nang+1,0:np+2) :: vrexp,vthexp,vphiexp
	real*8, dimension (ld) :: cs2
	real*8, dimension (0:np+1) :: mic,rhop,rp
	common /core_evolution/ cs2,mic,rhop,rp,vrexp,vthexp,vphiexp

	real*8 mdot,taccr
	real*8, dimension (0:nang+1,0:np+2) :: vradv,vthadv,vphiadv
	common /accretion/ mdot,taccr,vradv,vthadv,vphiadv

! Parameters.
	real*8, parameter :: msun=1.9891d33
	real*8, parameter :: facvel=1d-15*1d6*msun

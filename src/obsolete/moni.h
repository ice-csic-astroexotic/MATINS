!-----------------------------------------------------------------------
! Declarations.
!-----------------------------------------------------------------------
	real*8 btot,btortot,bout
	real*8 heli,jtot,etot,divbtot
	real*8 dbtot,dbout
	real*8 joule,joutot
	real*8 poynting,poytot
	real*8 joucor,joucortot
	common /moni/ btot,btortot,bout,heli,
     &  jtot,etot,divbtot,dbtot,dbout,joule,joutot,
     &  poynting,poytot,joucor,joucortot

! Parameters.
	real*8 clight,erg40,ke,kj
	parameter(clight=2.99792458d10)
! Conversion factor to 1e40 erg, erg40 = 1d24*1d15*1d-40 (in BMONITORS).
	parameter(erg40=1d-1)
	parameter(ke=erg40*(1d12/(clight*1d6*yrs))**2/4d0*pi)
	parameter(kj=erg40/(4d0*pi*1d6*yrs))

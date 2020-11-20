!-----------------------------------------------------------------------
! Contents:
! BWRITE
!-----------------------------------------------------------------------
subroutine global_quantities(iterb,dtbold,tbyear)

  use constants, only : UNIT_B, CLIGHT, T_YEAR, PI, KJ_ERG40_KM3, KE_ERG40_KM3
  use grid, only : belam, benu
  use grid, only: qcj, qcjc, j2
  use legpol, only : nleg, blout
  use grid, only : aphi, br, bth, bphi, er, eth, ephi
  use grid, only: kmax, lmax, np, jcore, rb, arear, areath, vol
  implicit none

! Subroutine arguments -------------------------------------------------------

      integer iterb
      real*8 dtbold,tbyear,btotold,boutold
      real*8, dimension(kmax,lmax) :: poy, divb
      real*8, parameter :: erg40=1.d-1
      real*8, parameter ::  = erg40*(UNIT_B/(CLIGHT*1d6*T_YEAR))**2/(4d0*PI)
	real*8 btot,btortot,bout,heli
	real*8 jtot,etot,divbtot,dbtot,dbout,joule
	real*8 poynting,joule_shock
	real*8 poyd,poyu,poyl,poyr

! Local variables ------------------------------------------------------------
	real*8 facleg
! Internally used variables.
	integer i,j,k,l
!-----------------------------------------------------------------------
! Notes:
!
! Bout = 2pi*int(rns->infty,mu[-1,1]) (Br^2+Bth^2)/8pi*e^nu
!      = rns^3 sum_l b_l^2 (l+1)/2(2l+1)*e^nu
! (see Viganò's notes)
! bout and btot in units  [10^40 erg]:
! factors of 	1d15 for integrating vol (km^3 -> cm^3)
!		1d24 for B^2 in (10^12 G)^2
!		1d-40 (erg -> 1d40 erg)
! to be used for check energy conservation
! (Pons et al. 2009, eq. 17)
! Electric field in statvolts/cm (=G) E = KE_ERG40_KM3 e^-\nu Enum
! Electric energy in 10^40 erg: erg40*KE_ERG40_KM3^2 e^-2\nu
!
! KE_ERG40_KM3 = (1e12G*km/(Myr*c))^2/4pi
! tauh [Km^2/Myr*(1d12 G)]
!-----------------------------------------------------------------------
      btot=0d0
      btortot=0d0
      bout=0d0
      heli=0d0
      jtot=0d0
      etot=0d0
      divbtot=0d0

! Relativistic corrections to vacuum energy for M=1.4 Msun, R=11.6 km.
! PARTIAL CORRECTION, TO BE CHECKED!
cc	corr=1d0
cc	corr(1)=1.513d0
cc	corr(2)=2.0689d0
cc	corr(3)=2.763d0

	do l=1,nleg
	   facleg=(dble(l)+1d0)/(2d0*(2d0*dble(l)+1d0))
	   bout=bout+erg40*facleg*blout(l)**2*rb(np)**3*benu(np)
cc	   boutrel=boutrel+corr(i)*erg40*facleg*bl(i)**2*rb(np)**3*benu(np)
	enddo


! Global quantities.
      do k=2,kmax
      do l=2,lmax
	      i=2*k-2
	      j=2*l-1

	      btot=btot+erg40*vol(i,j)*(br(i,j)**2+bth(i,j)**2+bphi(i,j)**2)
	      btortot=btortot+erg40*vol(i,j)*bphi(i,j)**2

	      divb(k,l)=1d0/(vol(i,j))
     &	      *(br(i,j+1)*arear(i,j+1)-br(i,j-1)*arear(i,j-1)
     &	      +bth(i+1,j)*areath(i+1,j)-bth(i-1,j)*areath(i-1,j))

	      divbtot=divbtot+dabs(divb(k,l))

	      jtot=jtot+erg40*vol(i,j)*j2(k,l)/(16d0*pi**2*benu(j))
	      etot=etot+KE_ERG40_KM3*vol(i,j)/benu(j)**2*(er(i,j)**2+eth(i,j)**2+ephi(i,j)**2)/(2d0)
	      heli=heli+erg40*vol(i,j)*aphi(i,j)*bphi(i,j)*1d5
      enddo
      enddo
  
!-----------------------------------------------------------------------
! JOULE HEATING:
! Jnum = e^nu 4pi J / c
! Qcj = -4*pi*c^2*eta*J^2*e^(2nu)
!     = -4*pi*eta*(Jnum/4*pi)^2 = -KJ_ERG40_KM3 eta J^2
! local joule rate (qcj(k,l)) is defined negative
! total joulerate (joule) is defined positive
! dissipated joule energy (joutot) is joulerate integrated in t
!
! with KJ_ERG40_KM3 = 1e24*1d15/(4*pi*1e6*year[s]) / [10^40] = 2.522e(-16)
!
! because
! qcj in units [1e40 erg/km^3 s = 1e25 erg/cm^3/s]
! eta in units [km^2/Myr]
! J in units   [1e12 G/km]
!
! Total energy rate loss by heating: qcj*vol because:
! joule        [1e40 erg/s]
! vol          [km^3 = 1e15 cm^3]
!
! Total energy lost by heating: joule*dtb*3.1536d13 because:
! joutot       [1e40 erg]
! dtbold       [yr = 3.16d7 s]
!-----------------------------------------------------------------------
      joule=0d0
      joule_shock=0d0
      do k=2,kmax
        do l=2,lmax
          i=2*k-2
          j=2*l-1
          joule=joule-vol(i,j)*qcj(k,l)
          joule_shock=joule_shock-vol(i,j)*qcjc(k,l)
        enddo
      enddo

!-----------------------------------------------------------------------
! POYNTING FLUX: S=int(Surface) [c(ExB)exp(2nu)/4pi]/vol
! local poy (poy(k,l)) is defined as -div(S)
! total poy (poynting) is defined as +div(S)
! 
! cE           [1e12 G*km/Myr = 3.1645569 G km/s ]
! B            [1e12 G]
! Surface/Vol  [1/km]
! S            [1e40 erg/km^3/s]
! one factor exp(nu) is already contained in E
!
! We choose the wind-selected B components at each inteface
!-----------------------------------------------------------------------

	 do k=2,kmax
		do l=2,lmax

		   i=2*k-2
		   j=2*l-1
 
		  poyr=-bphi(i+1,j)*er(i+1,j)*benu(j)
     &		 +0.5d0*(br(i+1,j+1)*ephi(i+1,j+1)*benu(j+1)
     &		 +br(i+1,j-1)*ephi(i+1,j-1)*benu(j-1))
 
		  poyl=-bphi(i-1,j)*er(i-1,j)*benu(j)
     &		 +0.5d0*(br(i-1,j+1)*ephi(i-1,j+1)*benu(j+1)
     &		 +br(i-1,j-1)*ephi(i-1,j-1)*benu(j-1))
 
		  poyu=benu(j+1)*(bphi(i,j+1)*eth(i,j+1)
     &		 -0.5d0*(bth(i+1,j+1)*ephi(i+1,j+1)
     &		 +bth(i-1,j+1)*ephi(i-1,j+1)))
 
		  poyd=benu(j-1)*(bphi(i,j-1)*eth(i,j-1)
     &		 -0.5d0*( bth(i+1,j-1)*ephi(i+1,j-1)
     &		 +bth(i-1,j-1)*ephi(i-1,j-1)))
 
	      poy(k,l)=-KJ_ERG40_KM3*(poyu*arear(i,j+1)-poyd*arear(i,j-1)
     &		 +(poyr*areath(i+1,j)-poyl*areath(i-1,j)))/vol(i,j)
 
		enddo
      enddo
   
 ! Calculation of the Poynting flux, considering only the boundaries.
      poynting=0d0
      do k=2,kmax
        i=2*k-2
        j=2*lmax
        poynting = poynting + KJ_ERG40_KM3*bphi(i,j)*eth(i,j)*benu(j)
     &    -0.5d0*benu(j)*(bth(i+1,j)*ephi(i+1,j)
     &    +bth(i-1,j)*ephi(i-1,j))*arear(i,j)
        j=jcore
        poynting = poynting - KJ_ERG40_KM3*bphi(i,j)*eth(i,j)*benu(j)
     &    -0.5d0*benu(j)*(bth(i+1,j)*ephi(i+1,j)
     &    +bth(i-1,j)*ephi(i-1,j))*arear(i,j)
      enddo
      dbout=(bout-boutold)/(dtbold*T_YEAR)
      dbtot=(btot-btotold)/(dtbold*T_YEAR)

      boutold=bout
      btotold=btot
!-----------------------------------------------------------------------
! Output.
!-----------------------------------------------------------------------

! General quantities and energy conservation.
      !  There is something wrong with the units in this part, needs to be fixed
	if(tbyear.eq.0d0)then
	   open(21,file='outb/bcons.dat')
	   open(22,file='outb/btot.dat')
	   write(21,*)
	else
	   open(21,file='outb/bcons.dat',access='append')
	   open(22,file='outb/btot.dat',access='append')
	   write(21,120)tbyear,dbtot,dbout,joule,dbtot+joule+poynting,
     &	   poynting,dtbold,joule_shock
	endif
	write(22,130)iterb,tbyear,btot,btortot,
     &	jtot,etot,heli,divbtot,bout
	close(21)
	close(22)
! Vacuum b.c. weights.
       if(sum(dabs(blout)).ne.0d0)then
       if(tbyear.eq.0d0)then
        open(20,file='outb/bl_out.yg')
       else
         open(20,file='outb/bl_out.yg',access='append')
       endif
	   write(20,*)'"Time =',tbyear
	   write(20,*)'"Label = bl(out)'
	   do l=1,nleg
	      write(20,'(i6,e12.4)')l,blout(l)
	   enddo
	   write(20,*)
	   close(20)
	endif

! Format statements.
110	format(a36,12e13.4)
120	format(8e13.4)
130	format(i9,11e15.6)

	return
	end



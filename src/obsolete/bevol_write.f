!-----------------------------------------------------------------------
! Write out field data and perform additional operations.
!
! Included by bevol.f
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! MONITOR AND WRITE: OLD B, J, E; 
! Factors contained in moni.h
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! JOULE HEATING:
! Jnum = e^nu 4pi J / c
! Qcj = -4*pi*c^2*eta*J^2*e^(2nu)
!     = -4*pi*eta*(Jnum/4*pi)^2 = -kj eta J^2
! local joule rate (qcj(k,l)) is defined negative
! total joulerate (joule) is defined positive
! dissipated joule energy (joutot) is joulerate integrated in t
!
! with kj = 1e24*1d15/(4*pi*1e6*year[s]) / [10^40] = 2.522e(-16)
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
! dtb          [Myr = 3.16d13 s]
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! POYNTING FLUX: S=int(Surface) [c(ExB)exp(2nu)/4pi]/vol
! local poy (poy(k,l)) is defined as -div(S)
! total poy (poynting) is defined as +div(S)
! accumulated poynting flux (poytot) is poynting flux integrated in t
! factor kj (see moni.h)
! cE           [1e12 G*km/Myr = 3.1645569 G km/s ]
! B            [1e12 G]
! Surface/Vol  [1/km]
! S            [1e40 erg/km^3/s]
! one factor exp(nu) is already contained in E
!
! We choose the wind-selected B components at each inteface
!-----------------------------------------------------------------------

! Initialize.
	if(tbyear.eq.0d0)then
	   poy=0d0
	   qcj=0d0
	   qcjc=0d0
	   qcjcr=0d0
	   qcjcth=0d0
	   enbold=0d0
	   joutot=0d0
	   joucortot=0d0
	   poytot=0d0
	endif

	if(iwrite.eq.1.or.iwrite.eq.2)then
	   denbdt=0d0
	   do k=2,kmax
	      do l=lc+1,lmax
		 i=2*k-2
		 j=2*l-1

		 enb(k,l)=bmed(k,l)**2*benu(j)/(8d0*pi)*erg40
		 enbtor(k,l)=bphi(i,j)**2*benu(j)/(8d0*pi)*erg40

		 if(tbyear.ne.0d0)
     &		 denbdt(k,l)=(enb(k,l)-enbold(k,l))/(dtbold*1.d6*yrs)
	      enddo
	   enddo
	endif

	joule=0d0
	joucor=0d0
	do k=2,kmax
	   do l=lc+1,lmax
	      i=2*k-2
	      j=2*l-1

	      joule=joule-vol(i,j)*qcj(k,l)
	      joucor=joucor-vol(i,j)*qcjc(k,l)
	   enddo
	enddo

	joutot=joutot+joule*dtb*1d6*yrs
	joucortot=joucortot+joucor*dtb*1d6*yrs
	poytot=poytot+poynting*dtb*1d6*yrs

	if(iwrite.eq.2)then
	   call BMONITORS
	   boutold=bout
	   btotold=btot
	else if(iwrite.eq.1)then
! Unused variables (denbdt, poy, qcj, qcjcr and qcjcth) removed.
	   call BWRITE(iterb,dtbold,tbyear,btotold,boutold)
	endif

	enbold=enb
	qcjcr=0d0
	qcjcth=0d0
	poynting=0d0

	do k=2,kmax
	   do l=lc+1,lmax
	      i=2*k-2
	      j=2*l-1

	      j2(k,l)=(jr(i,j)**2+jth(i,j)**2+jphi(i,j)**2)
	      qcj(k,l)=-ires*etab(i,j)*j2(k,l)*kj

!-----------------------------------------------------------------------
! Shock correction.
!-----------------------------------------------------------------------
	      jump=0.5d0*(bphi(i+1,j)-bphi(i-1,j))
              f2=fh(j+1)/(rb(j+1)**2*benu(j+1))
              f1=fh(j-1)/(rb(j-1)**2*benu(j-1))
              lamth(i,j)=-rb(j)**2*benu(j)**2*(f2-f1)/lr(j)
	      if(i.ne.2.and.i.ne.nang.and.jump*lamth(i,j).lt.0d0)
     &	      qcjcth(k,l)=kj*2d0/3d0*lamth(i,j)*jump**3*benu(j)/lth(j)

	      jump=0.5d0*(bphi(i,j+1)-bphi(i,j-1))
              lamr(i,j)=-2d0*cth(i)/(rb(j)*sth(i))
	      if(jump*lamr(i,j)*fh(j).lt.0d0)
     &	      qcjcr(k,l)=kj*2d0/3d0*lamr(i,j)*fh(j)*jump**3
     &	      *benu(j)/lr(j)

	      qcjc(k,l)=qcjcr(k,l)+qcjcth(k,l)

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
	      if(iwrite.eq.2)then
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

		 poy(k,l)=-kj*(poyu*arear(i,j+1)-poyd*arear(i,j-1)
     &		 +(poyr*areath(i+1,j)-poyl*areath(i-1,j)))/vol(i,j)
	      endif

! Calculation of the Poynting flux, considering only the boundaries.
	      if(l.eq.lmax)then
		 poyu=bphi(i,j+1)*eth(i,j+1)*benu(j+1)
     &		 -0.5d0*benu(j+1)*(bth(i+1,j+1)*ephi(i+1,j+1)
     &		 +bth(i-1,j+1)*ephi(i-1,j+1))

		 poynting = poynting + kj*poyu*arear(i,j+1)
	      else if(l.eq.lc+1)then
		 poyd=bphi(i,j-1)*eth(i,j-1)*benu(j-1)
     &		 -0.5d0*benu(j-1)*(bth(i+1,j-1)*ephi(i+1,j-1)
     &		 +bth(i-1,j-1)*ephi(i-1,j-1))

		 poynting = poynting - kj*poyd*arear(i,j-1)
	      endif

	   enddo
	enddo

	qcj(:,lmax)=0d0
	qcjc(:,lmax)=0d0
	qcjcr(:,lmax)=0d0
	qcjcth(:,lmax)=0d0

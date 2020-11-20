!-----------------------------------------------------------------------
! Contents:
! TWIST_DATA
! TWIST_MAX
! TWIST
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculate the twist for a set of magnetic field lines.
!-----------------------------------------------------------------------
	subroutine TWIST_DATA(P,Pc,a,ma,fit_type,rd,td,nx,nz,tbyear,Po,
     &	Ro)

	implicit none

! Input variables.
	integer ma,fit_type,nx,nz
	real*8 a(ma),P(nx,nz),Pc,rd(nx),td(nz)
	real*8 tbyear,Po,Ro

! Internally used variables.
	integer i,j,ns
	integer npi,npivac
	parameter(npi=20,npivac=10)
	real*8 Pi,rs(nx*nz),ts(nx*nz),fs(nx*nz)

! Quantities for the voltage and the time derivative of the twist.
	integer na
	parameter(na=3)
	real*8, dimension(na) :: a1,a2
	real*8 dft,dlnTdt,dV,P1,P2,phi2,Po1,Po2,T1,T2,tb1,tb2,V,Vo
	real*8, dimension(npi) :: dVdP,dphidt,phi,phi1
	real*8, dimension(npi) :: Psurf,rmax,Tsurf,theta
	save a1,phi1,Po1,tb1

! Call counter.
	integer, save :: ncall=0

!-----------------------------------------------------------------------
! Save previous call data needed for the calculation of the voltage.
!-----------------------------------------------------------------------
! Check dimensions.
! (The dimension of a1 has to be specified explicitly in order to save.)
	if(ma.ne.na)then
	   write(*,*)"TWIST_DATA: Wrong dimension of a1(na)!"
	endif
! Discard the first call, which comes through BINIT, before BEVOL.
! (Because the time stamp for the first two calls is tbyear=0.)
	if(ncall.eq.1)then
	   a1=a
	   Po1=Po
	   tb1=tbyear
	else if(ncall.eq.2)then
	   a2=a
	   Po2=Po
	   tb2=tbyear
	   if(tb2.eq.tb1)then
	      write(*,*)"TWIST_DATA: Error in the voltage calculation!"
	      stop
	   endif
	endif
! The amplitude of V includes a factor of 2*pi*c and Po.
! Note:
! B is in 10^12 G, R is in km, and tb is in years.
! 2*pi*1 year/s = 1.98*10^8.
! 1 statV = (c/10^8 cm/s) V ~ 300 V.
! The units of P are, [P] = 10^22 [G cm^2].
	Vo=Po/1.98d0

!-----------------------------------------------------------------------
! Output files.
!-----------------------------------------------------------------------
	if(ncall.eq.0)then
	   open(2,file="outb/twist.dat")
	   open(3,file="outb/twist_3d.dat")
	   open(4,file="outb/voltage.dat")

	   write(2,'(a)')"# Twist."
	   write(2,'(a1,a16,5a17,2a6)')"#","Pi","phi (twist)",
     &	   "rs(1)","ts(1)","rs(ns)","ts(ns)","ns","index"

	   write(4,'(a)')"# Voltage and twist."
	   write(4,'(a1,a16,6a17,a6)')"#","P (normalized)",
     &	   "Voltage (MV)","twist (rad)","dphidt","T (Po/Ro)",
     &	   "rmax (Ro)","theta (rad)","i"
	else
	   open(2,file="outb/twist.dat",access="append")
cc	   open(3,file="outb/twist_3d.dat",access="append")
	   open(3,file="outb/twist_3d.dat")
	   open(4,file="outb/voltage.dat",access="append")
	endif

	write(2,'(a6,e17.8)')'"time=',tbyear
	write(3,'(a1,e17.8)')"# tbyear:",tbyear

!-----------------------------------------------------------------------
! Current-free lines (no twist, T=0).
!-----------------------------------------------------------------------
	do i=1,npivac
	   Pi=0d0+Pc*dble(i)/dble(npivac)

! If ns=1 then the subroutine twist has failed.
	   call TWIST(P,Pi,a,ma,fit_type,rd,td,nx,nz,rs,ts,fs,ns)

! Write twist data,
	   if(ns.ne.1)then
	      write(2,'(2e17.8,f17.8,e17.8,f17.8,e17.8,2i6)')Pi,
     &	      fs(ns),rs(1),ts(1),rs(ns),ts(ns),ns,i

! Three-dimensional field line in spherical coordinates.
	      do j=1,ns
		 write(3,*)rs(j),ts(j),fs(j)
	      enddo
	      write(3,*)
	   endif
	enddo

!-----------------------------------------------------------------------
! Force-free lines (with twist).
!-----------------------------------------------------------------------
	do i=1,npi
	   Pi=Pc+(1d0-Pc)*dble(i-1)/dble(npi-1)

! If ns=1 then the subroutine twist has failed.
	   call TWIST(P,Pi,a,ma,fit_type,rd,td,nx,nz,rs,ts,fs,ns)

! Current twist at each footpoint P(i)=Pi.
	   phi2=fs(ns)

! Derivative of the voltage wrt P.
	   if(ncall.eq.2)then
	      P2=Pi
	      P1=Pi*Po2/Po1

! Maximum radial distance of the field line.
	      rmax(i)=MAXVAL(rs(1:ns))

! The toroidal function at the two times for the relevant values of P.
	      call FTOR(P1,a1,T1,dft,ma,fit_type)
	      call FTOR(P2,a2,T2,dft,ma,fit_type)

! The amplitude of T is Po/Ro, but Ro is constant in time and drops out.
	      if(T2.eq.0d0)then
		 dlnTdt=0d0
	      else
		 dlnTdt=(Po2*T2-Po1*T1)/(tb2-tb1)/(Po2*T2)
	      endif

! Derivative of the voltage wrt P.
	      dVdP(i)=dlnTdt*fs(ns)

! Twist and its time derivative.
! Warning: Note that this is calculated here at a fixed P, but more
! generally it should be calculated at a fixed angle (theta).
! When P as a function of theta changes significantly, this will have to
! be taken into account.
	      phi(i)=phi2
	      dphidt(i)=(phi2-phi1(i))/(tb2-tb1)

! Poloidal and toroidal functions and angle.
	      Psurf(i)=Pi
	      Tsurf(i)=T2
	      theta(i)=ts(1)
	   endif

! Update the twist for the previous time step.
	   phi1(i)=phi2

! Write twist data.
! Include the last point for P=1 where the twist is 0.
	   write(2,'(2e17.8,f17.8,e17.8,f17.8,e17.8,2i6)')Pi,
     &	      fs(ns),rs(1),ts(1),rs(ns),ts(ns),ns,i

! Analytic check for a dipole. (Pass on P for a dipole.)
! The toroidal field is assumed to be of the form T=s(P-Pc).
! The twist then would be: phi=2*T*sqrt(1-P)/P^2.
! Very good agreement for sufficiently high resolution.
cc	   write(*,*)"TWIST_DATA: Analytic check for the twist:",
cc     &	      fs(ns),2d0*a(1)*(Pi-Pc)*sqrt(1d0-Pi)/Pi**2

! Three-dimensional field line in spherical coordinates.
	   do j=1,ns
	      write(3,*)rs(j),ts(j),fs(j)
	   enddo
	   write(3,*)
	enddo
	write(2,*)

! Integrate over P for the voltage.
	if(ncall.eq.2)then
	   V=0d0
	   write(4,'(a6,e17.8)')'"time=',tbyear
	   write(4,'(7e17.8,i6)')Psurf(npi),V,phi(npi),dphidt(npi),
     &	   Tsurf(npi),1d0,theta(npi),npi
	   do i=npi-1,1,-1
	      dV=-Vo*(dVdP(i+1)+dVdP(i))*(Psurf(i+1)-Psurf(i))/2d0
	      V=V+dV
	      write(4,'(7e17.8,i6)')Psurf(i),V,phi(i),dphidt(i),
     &	      Tsurf(i),rmax(i),theta(i),i
	   enddo
	   write(4,*)

! Update the previous time step data.
	   a1=a2
	   Po1=Po2
	   tb1=tb2
	endif

	close(2)
	close(3)
	close(4)

! Update ncall.
	if(ncall.eq.1) ncall=2
	if(ncall.eq.0) ncall=1

	return
	end
!-----------------------------------------------------------------------
! Find the maximum twist.
!-----------------------------------------------------------------------
	subroutine TWIST_MAX(P,Pc,a,ma,fit_type,rd,td,nx,nz,fmax,imax)

	implicit none

! Input variables.
	integer ma,fit_type,nx,nz
	real*8 a(ma),P(nx,nz),Pc
	real*8 rd(nx),td(nz)

! Output variables.
	integer imax
	real*8 fmax

! Coordinates of the field line.
	integer ns
	real*8 rs(nx*nz),ts(nx*nz),fs(nx*nz)

! Internally used variables.
	integer icount,ilimit
	parameter(ilimit=20)
	real*8 f1,f2,f3,fmid
	real*8 P1,P2,P3,Pf,Pmax,Pmid,Prange

!-----------------------------------------------------------------------
! This routine will fail if there are multiple toroidal domains.
!-----------------------------------------------------------------------

! Return analytic result for T=0.
cc	if(s.eq.0d0)then
cc	   fmax=0d0
cc	   Pmax=0d0
cc	   imax=1
cc
cc	   return
cc	endif

!-----------------------------------------------------------------------
! Initial values.
!-----------------------------------------------------------------------
! Define the maximum (final) value of P (on the equator).
! Should be 1 (P and T given through FORCE_FREE are normalized).
cc	Pf=1d0
	Pf=maxval(P(1,:))

! Initial values.
	P1=Pc
	P2=(Pc+Pf)/2d0
	P3=Pf

! Calculate the twist for P2.
	call TWIST(P,P2,a,ma,fit_type,rd,td,nx,nz,rs,ts,fs,ns)

! If ns=1 then the subroutine twist has failed.
	if(ns.eq.1) goto 10

! Initial values of the twist.
	f1=0d0
	f2=fs(ns)
	f3=0d0

!-----------------------------------------------------------------------
! Refine through bisection.
!-----------------------------------------------------------------------
	icount=0
	Prange=1d0
	do while(icount.lt.ilimit.and.Prange.gt.0.01d0)

! Bisection from the left.
	   Pmid=(P1+P2)/2d0
	   call TWIST(P,Pmid,a,ma,fit_type,rd,td,nx,nz,rs,ts,fs,ns)
! If ns=1 then the subroutine twist has failed.
	   if(ns.eq.1) goto 10
	   fmid=fs(ns)
	   if(fmid.lt.f2)then
	      P1=Pmid
	      f1=fmid
	   else if(fmid.gt.f2)then
	      P3=P2
	      f3=f2
	      P2=Pmid
	      f2=fmid
	   endif

! Bisection from the right.
	   Pmid=(P2+P3)/2d0
	   call TWIST(P,Pmid,a,ma,fit_type,rd,td,nx,nz,rs,ts,fs,ns)
! If ns=1 then the subroutine twist has failed.
	   if(ns.eq.1) goto 10
	   fmid=fs(ns)
	   if(fmid.lt.f2)then
	      P3=Pmid
	      f3=fmid
	   else if(fmid.gt.f2)then
	      P1=P2
	      f1=f2
	      P2=Pmid
	      f2=fmid
	   endif

! Update loop control variables.
	   icount=icount+1
	   Prange=(P3-P1)/(Pf-Pc)
	enddo

!-----------------------------------------------------------------------
! Results.
!-----------------------------------------------------------------------
! Check for maximum number of iterations.
	if(icount.eq.ilimit)then
	   write(*,*)"TWIST_MAX: Maximum iterations reached!"
	   stop
	endif

! Normal termination of subroutine.
	fmax=f2
	Pmax=P2
	imax=icount
	return

! Abnormal termination of subroutine.
! Set imax equal to -1 when the subroutine has failed.
10	continue
	fmax=-1d0
	Pmax=-1d0
	imax=-1
	return

	end
!-----------------------------------------------------------------------
! Numerical calculation of the field line and twist for a given P.
!-----------------------------------------------------------------------
	subroutine TWIST(P,Pi,a,ma,fit_type,rd,td,nx,nz,rs,ts,fs,ns)

	implicit none

! Input variables.
	integer fit_type,ma,nx,nz
	real*8 a(ma),P(nx,nz),Pi,rd(nx),td(nz)

! Output variables.
	integer ns
	real*8 rs(nx*nz),ts(nx*nz),fs(nx*nz)

! Internally used variables.
	integer l,m,n
	integer m0,m1,m2,m3,n0,n1,n2,n3,npts1,npts2
	integer icell(4),is(nx*nz)
	real*8 dPs(nx*nz),dPdr(nx*nz),dPdt(nx*nz),fs_alt(nx*nz)
	real*8 dr,dt,ds,dr_check,dt_check,r1,t1,left,right,down,up
	real*8 factor,f11,f12,f21,f22,pol,tor
	real*8 phi,phi1,phi2,phi_trapezoid
	real*8 ft,dft

! Numerical error.
	real*8 err
	parameter(err=1d-15)

! Initialize the count of points on the field line.
	ns=0

! Initialize arrays.
	dPs=0d0
	dPdr=0d0
	dPdt=0d0
	rs=0d0
	ts=0d0
	fs=0d0
	fs_alt=0d0
	is=0

! Find the first point on the surface.
! It cannot be at n=nz.
	m=1
	do n=1,nz-1
! Find the range where the first point is. 
	   factor=(Pi-P(m,n))*(Pi-P(m,n+1))
	   if(factor.le.0d0)then

! If the two points have the same value (as Pi), then stop.
	      if(P(m,n).eq.P(m,n+1))then
cc		 write(*,*)"TWIST: Can't locate footprint!"
		 goto 30
cc		 stop "TWIST: Can't locate footprint!"
	      endif

! Index of the field line.
	      ns=ns+1
! Radius and angle of the field line.
	      factor=(Pi-P(m,n))/(P(m,n+1)-P(m,n))
	      rs(ns)=rd(m)
	      ts(ns)=td(n)+(td(n+1)-td(n))*factor
! Define corners of the search cell.
	      m1=m
	      m2=m+1
	      n1=n
	      n2=n+1
! Corner of the lower cell. (No lower cell at the surface.)
	      m0=m1
! Corner of the left cell.
	      n0=n1-1
	      if(n0.eq.0) n0=1
! Corner of the right cell.
	      n3=n2+1
	      if(n3.eq.nz+1) n3=nz
! Define the sides to be searched in the next round.
! Sides: 1 refers to bottom, 2 to left, 3 to top, and 4 to right.
! Values: 0 for not searching, 1 for searching.
! The side from where the cell is entered should not be searched.
	      icell(1)=0
	      icell(2)=1
	      icell(3)=1
	      icell(4)=1
! Direction of crossing (radial=1 or angular=2).
	      is(ns)=2
! Angular derivative in the interval.
	      dPs(ns)=(P(m1,n2)-P(m1,n1))/(td(n2)-td(n1))
! Angular and radial derivatives at each point on the line.
! Forward difference at the surface for the radial derivative.
! The derivatives need to be interpolated.
	      left=(P(m2,n1)-P(m0,n1))/(rd(m2)-rd(m0))
	      right=(P(m2,n2)-P(m0,n2))/(rd(m2)-rd(m0))
	      dPdr(ns)=left+(ts(ns)-td(n1))*(right-left)/(td(n2)-td(n1))
! Constant slope for the angular derivative.
cc	      dPdt(ns)=dPs(ns)
! Centered difference for the angular derivative.
! The derivatives need to be interpolated.
	      left=(P(m1,n2)-P(m1,n0))/(td(n2)-td(n0))
	      right=(P(m1,n3)-P(m1,n1))/(td(n3)-td(n1))
	      dPdt(ns)=left+(ts(ns)-td(n1))*(right-left)/(td(n2)-td(n1))

	      goto 10
	   endif
	enddo

! No point found.
cc	write(*,*)"TWIST: No point on surface found!"
	goto 30
cc	stop "TWIST: No point on surface found!"

! First point successfully found.
10	continue

! Trace field line.
	do l=1,nx*nz-1
	   f11=sign(1d0,Pi-P(m1,n1))
	   f12=sign(1d0,Pi-P(m1,n2))
	   f21=sign(1d0,Pi-P(m2,n1))
	   f22=sign(1d0,Pi-P(m2,n2))

! Check the bottom side.
	   if(icell(1).ne.0.and.f11*f12.lt.0d0)then
! The point is between (m1,n1) and (m1,n2).
! Index of the field line.
	      ns=ns+1
! Radius and angle of the field line.
	      factor=(Pi-P(m1,n1))/(P(m1,n2)-P(m1,n1))
	      rs(ns)=rd(m1)
	      ts(ns)=td(n1)+(td(n2)-td(n1))*factor
! Define corners of the search cell.
	      m1=m1-1
	      m2=m2-1
	      n1=n1
	      n2=n2
! At stellar surface m1=m2=1.
	      if(m1.eq.0) m1=1
! Corner of the upper cell.
	      m3=m2+1
! Corner of the left cell.
	      n0=n1-1
	      if(n0.eq.0) n0=1
! Corner of the right cell.
	      n3=n2+1
	      if(n3.eq.nz+1) n3=nz
! Define the sides to be searched in the next round.
	      icell(1)=1
	      icell(2)=1
	      icell(3)=0
	      icell(4)=1
! Direction of crossing (radial=1 or angular=2).
	      is(ns)=2
! Angular derivative in the interval.
	      dPs(ns)=(P(m1,n2)-P(m1,n1))/(td(n2)-td(n1))
! Angular and radial derivatives at each point on the line.
! Centered difference for the radial derivative.
! Backward difference at the surface (for m1=m2=1).
! The derivatives need to be interpolated.
	      left=(P(m3,n1)-P(m1,n1))/(rd(m3)-rd(m1))
	      right=(P(m3,n2)-P(m1,n2))/(rd(m3)-rd(m1))
	      dPdr(ns)=left+(ts(ns)-td(n1))*(right-left)/(td(n2)-td(n1))
! Constant slope for the angular derivative.
cc	      dPdt(ns)=dPs(ns)
! Centered difference for the angular derivative.
! The derivatives need to be interpolated.
	      left=(P(m2,n2)-P(m2,n0))/(td(n2)-td(n0))
	      right=(P(m2,n3)-P(m2,n1))/(td(n3)-td(n1))
	      dPdt(ns)=left+(ts(ns)-td(n1))*(right-left)/(td(n2)-td(n1))

! Check the left side.
! THIS SECTION NEEDS TO BE CHECKED ANALYTICALLY! (A dipole won't do!)
	   else if(icell(2).ne.0.and.f11*f21.lt.0d0)then
! The point is between (m1,n1) and (m2,n1).
! Index of the field line.
	      ns=ns+1
! Radius and angle of the field line.
	      factor=(Pi-P(m1,n1))/(P(m2,n1)-P(m1,n1))
	      rs(ns)=rd(m1)+(rd(m2)-rd(m1))*factor
	      ts(ns)=td(n1)
! Define corners of the search cell.
	      m1=m1
	      m2=m2
	      n1=n1-1
	      n2=n2-1
! Near the surface m0=m1=1.
	      m0=m1-1
	      if(m0.eq.0) m0=1
	      m3=m2+1
	      if(m3.eq.nx+1) m3=nx
! Corner of the right cell.
	      n3=n2+1
	      if(n3.eq.nz+1) n3=nz
! Define the sides to be searched in the next round.
	      icell(1)=1
	      icell(2)=1
	      icell(3)=1
	      icell(4)=0
! Direction of crossing (radial=1 or angular=2).
	      is(ns)=1
! Radial derivative in the interval.
	      dPs(ns)=(P(m2,n2)-P(m1,n2))/(rd(m2)-rd(m1))
! Angular and radial derivatives at each point on the line.
! Constant slope for the radial derivative.
cc	      dPdr(ns)=dPs(ns)
! Centered difference for the radial derivative.
! The derivatives need to be interpolated.
	      down=(P(m2,n2)-P(m0,n2))/(rd(m2)-rd(m0))
	      up=(P(m3,n2)-P(m1,n2))/(rd(m3)-rd(m1))
	      dPdr(ns)=down+(rs(ns)-rd(m1))*(up-down)/(rd(m2)-rd(m1))
! Centered difference for the angular derivative.
! The derivatives need to be interpolated.
	      down=(P(m1,n3)-P(m1,n1))/(td(n3)-td(n1))
	      up=(P(m2,n3)-P(m2,n1))/(td(n3)-td(n1))
	      dPdt(ns)=down+(rs(ns)-rd(m1))*(up-down)/(rd(m2)-rd(m1))

! Check the top side.
	   else if(icell(3).ne.0.and.f21*f22.lt.0d0)then
! The point is between (m2,n1) and (m2,n2).
! Index of the field line.
	      ns=ns+1
! Radius and angle of the field line.
	      factor=(Pi-P(m2,n1))/(P(m2,n2)-P(m2,n1))
	      rs(ns)=rd(m2)
	      ts(ns)=td(n1)+(td(n2)-td(n1))*factor
! Define corners of the search cell.
	      m1=m1+1
	      m2=m2+1
	      n1=n1
	      n2=n2
! Corner of the lower cell.
	      m0=m1-1
! Corner of the left cell.
	      n0=n1-1
	      if(n0.eq.0) n0=1
! Corner of the right cell.
	      n3=n2+1
	      if(n3.eq.nz+1) n3=nz
! Define the sides to be searched in the next round.
	      icell(1)=0
	      icell(2)=1
	      icell(3)=1
	      icell(4)=1
! Direction of crossing (radial=1 or angular=2).
	      is(ns)=2
! Angular derivative in the interval.
	      dPs(ns)=(P(m2,n2)-P(m2,n1))/(td(n2)-td(n1))
! Angular and radial derivatives at each point on the line.
! Centered difference for the radial derivative.
! The derivatives need to be interpolated.
	      left=(P(m2,n1)-P(m0,n1))/(rd(m2)-rd(m0))
	      right=(P(m2,n2)-P(m0,n2))/(rd(m2)-rd(m0))
	      dPdr(ns)=left+(ts(ns)-td(n1))*(right-left)/(td(n2)-td(n1))
! Constant slope for the angular derivative.
cc	      dPdt(ns)=dPs(ns)
! Centered difference for the angular derivative.
! The derivatives need to be interpolated.
	      left=(P(m1,n2)-P(m1,n0))/(td(n2)-td(n0))
	      right=(P(m1,n3)-P(m1,n1))/(td(n3)-td(n1))
	      dPdt(ns)=left+(ts(ns)-td(n1))*(right-left)/(td(n2)-td(n1))

! Check the right side.
	   else if(icell(4).ne.0.and.f22*f12.lt.0d0)then
! The point is between (m2,n2) and (m1,n2).
! Index of the field line.
	      ns=ns+1
! Radius and angle of the field line.
	      factor=(Pi-P(m1,n2))/(P(m2,n2)-P(m1,n2))
	      rs(ns)=rd(m1)+(rd(m2)-rd(m1))*factor
	      ts(ns)=td(n2)
! Define corners of the search cell.
	      m1=m1
	      m2=m2
	      n1=n1+1
	      n2=n2+1
! Near the surface m0=m1=1.
	      m0=m1-1
	      if(m0.eq.0) m0=1
	      m3=m2+1
	      if(m3.eq.nx+1) m3=nx
! Corner of the left cell.
	      n0=n1-1
	      if(n0.eq.0) n0=1
! Define the sides to be searched in the next round.
	      icell(1)=1
	      icell(2)=0
	      icell(3)=1
	      icell(4)=1
! Direction of crossing (radial=1 or angular=2).
	      is(ns)=1
! Radial derivative in the interval.
	      dPs(ns)=(P(m2,n1)-P(m1,n1))/(rd(m2)-rd(m1))
! Angular and radial derivatives at each point on the line.
! Constant slope for the radial derivative.
cc	      dPdr(ns)=dPs(ns)
! Centered difference for the radial derivative.
! The derivatives need to be interpolated.
	      down=(P(m2,n1)-P(m0,n1))/(rd(m2)-rd(m0))
	      up=(P(m3,n1)-P(m1,n1))/(rd(m3)-rd(m1))
	      dPdr(ns)=down+(rs(ns)-rd(m1))*(up-down)/(rd(m2)-rd(m1))
! Centered difference for the angular derivative.
! The derivatives need to be interpolated.
	      down=(P(m1,n2)-P(m1,n0))/(td(n2)-td(n0))
	      up=(P(m2,n2)-P(m2,n0))/(td(n2)-td(n0))
	      dPdt(ns)=down+(rs(ns)-rd(m1))*(up-down)/(rd(m2)-rd(m1))

	   else
cc	      write(*,*)"TWIST: Field line lost!"
	      goto 30
cc	      stop "TWIST: Field line lost!"
	   endif

! Check for surface and grid boundaries.
! Do until you come back the stellar surface (m1=m2=1).
	   if(m1.eq.m2)then
	      goto 20
	   endif

! Error if the field line traced reaches the outer boundary.
	   if(m2.eq.nx+1)then
cc	      write(*,*)"TWIST: Field line reaches outer boundary!"
	      goto 30
cc	      stop "TWIST: Field line reaches outer boundary!"
	   endif

! Error if either of the axes is reached.
	   if(n1.eq.0.or.n2.eq.nz+1)then
cc	      write(*,*)"TWIST: An axis was reached!"
	      goto 30
cc	      stop "TWIST: An axis was reached!"
	   endif

	enddo

! Error if the stellar surface is not reached.
! Something must be wrong with P.
cc	write(*,*)"TWIST: Field line does not come back to the surface!"
	goto 30
cc	stop "TWIST: Field line does not come back to the surface!"

! Successful tracing of field line.
20	continue

! Calculate twist.
	phi=0d0
	phi1=0d0
	phi2=0d0
	phi_trapezoid=0d0
	npts1=0
	npts2=0

! Initial points in order to calculate the displacements.
	r1=rs(1)
	t1=ts(1)
	dr_check=rd(2)-rd(1)
	dt_check=td(2)-td(1)
! Toroidal field.
	call FTOR(Pi,a,ft,dft,ma,fit_type)
	tor=ft
! Skip the first point! (The zero of phi.)
	do l=2,ns
! Integration with respect to the poloidal line element.
! Step integration. (Most accurate.)
	   dr=rs(l)-rs(l-1)
	   dt=ts(l)-ts(l-1)
	   ds=sqrt(dr**2+(rs(l)*dt)**2)

	   pol=sqrt(dPdr(l)**2+(dPdt(l)/rs(l))**2)
	   phi=phi+tor*ds/pol/rs(l)/sin(ts(l))
	   fs(l)=phi
	enddo

! Normal termination of subroutine.
	return

! Abnormal termination of subroutine.
! Set ns equal to 1 when the subroutine has failed.
30	continue
	rs=0d0
	ts=0d0
	fs=0d0
	ns=1
	return

	end

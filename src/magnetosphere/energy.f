!-----------------------------------------------------------------------
! Contents:
! ENERGY
! ENERGY_SURFACE
! ENERGY_VOLUME
! MAGNETIC_FIELD
!-----------------------------------------------------------------------
	subroutine ENERGY(P,T,rd,td,nx,nz,E,Ealt,Etor,H)

	implicit none

! Input variables.
	integer nx,nz
	real*8 P(nx,nz),T(nx,nz),rd(nx),td(nz)

! Output variables.
	real*8 E,Ealt,Etor,H

! Magnetic field components.
	real*8 Br(nx,nz),Bt(nx,nz),Bf(nx,nz)

! Internally used variables.
	integer m,n,ms
	real*8 Brs(nz),Bts(nz),Bfs(nz)
	real*8 Epol,Etot,Evac,Evol,Hvol
	real*8 vol

! Magnetic field components.
	call MAGNETIC_FIELD(P,T,rd,td,nx,nz,Br,Bt,Bf)

! Forward and backward differences are first order.
! (Doubling the radial grid halves the error.)
! Centered difference is second order.

! Volume integrations for energy and helicity, from rin to rout.
! Helicity is zero for vacuum, so Hvol is the total helicity.
	call ENERGY_VOLUME(P,rd,td,Br,Bt,Bf,nx,nz,Epol,Etor,Evol,Hvol)

! Surface integral for the energy at rout. (For the external vacuum.)
	ms=nx
	Brs=Br(ms,:)
	Bts=Bt(ms,:)
	Bfs=Bf(ms,:)
	call ENERGY_SURFACE(rd(ms),td,Brs,Bts,Bfs,nz,Evac)

! Total energy.
	Etot=Evac+Evol

! Surface integral for the energy at rin.
! If there are no surface currents, this should give the total energy.
	ms=1
	Brs=Br(ms,:)
	Bts=Bt(ms,:)
	Bfs=Bf(ms,:)
	call ENERGY_SURFACE(rd(ms),td,Brs,Bts,Bfs,nz,Ealt)

! Return output.
	E=Etot
	H=Hvol

	return
	end
!-----------------------------------------------------------------------
! Surface integration for energy.
!-----------------------------------------------------------------------
	subroutine ENERGY_SURFACE(r,td,Brs,Bts,Bfs,nz,Esurf)

	implicit none

! Input variables.
	integer nz
	real*8 r,td(nz),Brs(nz),Bts(nz),Bfs(nz)

! Output variables.
	real*8 Esurf

! Internally used variables.
	integer n
	real*8 Int1,Int2

! Trapezoidal integration.
	Esurf=0d0
	Int1=(Brs(1)**2-Bts(1)**2-Bfs(1)**2)*sin(td(1))
	do n=2,nz
	   Int2=(Brs(n)**2-Bts(n)**2-Bfs(n)**2)*sin(td(n))
	   Esurf=Esurf+(Int1+Int2)*(td(n)-td(n-1))/2d0
	   Int1=Int2
	enddo
	Esurf=Esurf*r**3/4d0

	return
	end
!-----------------------------------------------------------------------
! Volume integration for energy and helicity.
!-----------------------------------------------------------------------
	subroutine ENERGY_VOLUME(P,rd,td,Br,Bt,Bf,nx,nz,Epol,Etor,Evol,
     &	Hvol)

	implicit none

! Input variables.
	integer nx,nz
	real*8 P(nx,nz),rd(nx),td(nz)
	real*8 Br(nx,nz),Bt(nx,nz),Bf(nx,nz)

! Output variables.
	real*8 Epol,Etor,Evol,Hvol

! Internally used variables.
	integer m,n
	real*8 dr,dt,pi
	real*8 Epol1,Epol2,Etor1,Etor2,Ipol1,Ipol2,Itor1,Itor2
	real*8 Hvol1,Hvol2,Ihel1,Ihel2

! Pi.
	pi=2d0*asin(1d0)

! Volume integrations for the energy and helicity.
! Trapezoidal integration.
! Integrations over phi are implicit.
	Epol=0d0
	Etor=0d0
	Evol=0d0
	Hvol=0d0
	do m=1,nx
! Integration over theta.
	   Epol2=0d0
	   Etor2=0d0
	   Hvol2=0d0
	   Ipol1=(Br(m,1)**2+Bt(m,1)**2)*sin(td(1))
	   Itor1=Bf(m,1)**2*sin(td(1))
	   Ihel1=P(m,1)*Bf(m,1)
	   do n=2,nz
	      Ipol2=(Br(m,n)**2+Bt(m,n)**2)*sin(td(n))
	      Itor2=Bf(m,n)**2*sin(td(n))
	      Ihel2=P(m,n)*Bf(m,n)
	      Epol2=Epol2+(Ipol1+Ipol2)*(td(n)-td(n-1))/2d0
	      Etor2=Etor2+(Itor1+Itor2)*(td(n)-td(n-1))/2d0
	      Hvol2=Hvol2+(Ihel1+Ihel2)*(td(n)-td(n-1))/2d0
	      Ipol1=Ipol2
	      Itor1=Itor2
	      Ihel1=Ihel2
	   enddo
	   Epol2=Epol2*rd(m)**2/4d0
	   Etor2=Etor2*rd(m)**2/4d0
	   Hvol2=Hvol2*rd(m)*4d0*pi
! Integration over radius.
	   if(m.gt.1)then
	      Epol=Epol+(Epol1+Epol2)*(rd(m)-rd(m-1))/2d0
	      Etor=Etor+(Etor1+Etor2)*(rd(m)-rd(m-1))/2d0
	      Hvol=Hvol+(Hvol1+Hvol2)*(rd(m)-rd(m-1))/2d0
	   endif
	   Epol1=Epol2
	   Etor1=Etor2
	   Hvol1=Hvol2
	enddo
	Evol=Epol+Etor

	return
	end
!-----------------------------------------------------------------------
! Magnetic field components.
!-----------------------------------------------------------------------
	subroutine MAGNETIC_FIELD(P,T,rd,td,nx,nz,Br,Bt,Bf)

	implicit none

! Input variables.
	integer nx,nz
	real*8 P(nx,nz),T(nx,nz),rd(nx),td(nz)

! Output variables.
	real*8 Br(nx,nz),Bt(nx,nz),Bf(nx,nz)

! Internally used variables.
	integer m,n
	real*8 dPdr(nx,nz),dPdm(nx,nz)
	real*8 a,b,c,P1,P2,P3,r1,r2,r3

! Auxiliary definitions.
! Radial derivative of the poloidal function.
	dPdr=0d0
	do m=2,nx-1
	   dPdr(m,:)=(P(m+1,:)-P(m-1,:))/(rd(m+1)-rd(m-1))
	enddo
! Forward and backward difference.
	dPdr(1,:)=(P(2,:)-P(1,:))/(rd(2)-rd(1))
	dPdr(nx,:)=(P(nx,:)-P(nx-1,:))/(rd(nx)-rd(nx-1))
! Angular derivative of the poloidal function (with respect to mu).
	dPdm=0d0
	do n=2,nz-1
	   dPdm(:,n)=(P(:,n+1)-P(:,n-1))/(cos(td(n+1))-cos(td(n-1)))
	enddo
! Forward and backward difference.
	dPdm(:,1)=(P(:,2)-P(:,1))/(cos(td(2))-cos(td(1)))
	dPdm(:,nz)=(P(:,nz)-P(:,nz-1))/(cos(td(nz))-cos(td(nz-1)))

! Extrapolation of the radial derivative at grid boundaries.
! Using P(1,:), P(2,:) and P(3,:) construct the quadratic polynomial
! f(x) = a*x^2 + b*x + c, and calculate the derivative f'(x) = 2a*x + b.
! Derivative at the inner boundary (m=1).
	r1=rd(1)
	r2=rd(2)
	r3=rd(3)
	do n=1,nz
	   P1=P(1,n)
	   P2=P(2,n)
	   P3=P(3,n)
! Coefficients.
	   a=(P1-P2)/(r1-r2)/(r2-r3)-(P1-P3)/(r1-r3)/(r2-r3)
	   b=(P1-P2)/(r1-r2)-a*(r1+r2)
	   c=P1-a*r1**2-b*r1
! Derivative.
	   dPdr(1,n)=2d0*a*r1+b
	enddo
! Derivative at the outer boundary (m=nx).
	r1=rd(nx-2)
	r2=rd(nx-1)
	r3=rd(nx)
	do n=1,nz
	   P1=P(nx-2,n)
	   P2=P(nx-1,n)
	   P3=P(nx,n)
! Coefficients.
	   a=(P1-P2)/(r1-r2)/(r2-r3)-(P1-P3)/(r1-r3)/(r2-r3)
	   b=(P1-P2)/(r1-r2)-a*(r1+r2)
	   c=P1-a*r1**2-b*r1
! Derivative.
	   dPdr(nx,n)=2d0*a*r3+b
	enddo

! Magnetic field components.
	Br=0d0
	Bt=0d0
	Bf=0d0
	do m=1,nx
! Radial component.
	   do n=1,nz
	      Br(m,n)=-dPdm(m,n)/rd(m)**2
	   enddo
! Theta and phi components must be zero along the axis.
	   do n=2,nz-1
	      Bt(m,n)=-dPdr(m,n)/rd(m)/sin(td(n))
	      Bf(m,n)=T(m,n)/rd(m)/sin(td(n))
	   enddo
	enddo

	return
	end

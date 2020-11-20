!-----------------------------------------------------------------------
! Test for energy, helicity and twist in the magnetosphere.
!-----------------------------------------------------------------------
	program TEST_MSPH

	implicit none

	integer m,n,nx,nz
	parameter(nx=100,nz=100)
	real*8 dr,dt,mu,Pc,Pmax,pi,rin,rout,s,sigma
	parameter(rin=1d0,rout=10d0)
	real*8, dimension(nx) :: rd
	real*8, dimension(nz) :: cth,sth,td
	real*8, dimension(nx,nz) :: Br,Bt,Bf,P,Pd,T

! Multipole expansion.
	integer ell,ellmax,l
	parameter(ellmax=3)
	real*8 dpl(ellmax,nz),DPLGNDR

! Output returned by subroutines being tested.
	integer imax
	real*8 E,Ealt,Epol,Etor,fmax,H

! Variables needed by FTOR.
	integer fit_type,fn_type,ma
	parameter(fit_type=2,ma=3)
	real*8 a(ma),f,dfdx,x
	common /function_type/ fn_type

!-----------------------------------------------------------------------
! Definitions.
!-----------------------------------------------------------------------
! Needed by FTOR_NONL (called by FTOR).
	fn_type=2

! Definition of pi.
	pi=2d0*asin(1d0)

! Radial grid.
	rd=0d0
	dr=(rout-rin)/dble(nx-1)
	do m=1,nx
	   rd(m)=rin+dble(m-1)*dr
	enddo

! Angular grid.
	td=0d0
	dt=pi/dble(nz-1)
	do n=1,nz
	   td(n)=dble(n-1)*dt
	   cth(n)=cos(td(n))
	   sth(n)=sin(td(n))
	enddo

! Legendre polynomials.
	do l=1,ellmax
	   do n=2,nz-1
	      mu=cth(n)
	      dpl(l,n)=DPLGNDR(l,0,mu)
	   enddo
! Derivatives at endpoints. (The axis is not implemented in DPLGNDR.)
	   dpl(l,1)=dble(l*(l+1))/2d0
	   dpl(l,nz)=(-1)**(l+1)*dpl(l,1)
	enddo

! Multipole index for the initial guess.
	ell=1

! Initial guess.
	Pd=0d0
	do n=1,nz
	   mu=cth(n)
	   do m=1,nx
! Single multipole.
	      Pd(m,n)=(1d0-mu**2)*dpl(ell,n)/rd(m)**ell
! Sum of several multipoles.
cc	      Pd(m,n)=0.5d0*(1d0-mu**2)*dpl(ell,n)/rd(m)**ell
cc     &	      +0.5d0*(1d0-mu**2)*dpl(ell+1,n)/rd(m)**(ell+1)
cc     &	      +0.5d0*(1d0-mu**2)*dpl(ell+2,n)/rd(m)**(ell+2)
	   enddo
	enddo

!-----------------------------------------------------------------------
! Checks.
!-----------------------------------------------------------------------
! Screen output.
	write(*,*)
	write(*,'(a10,2i6)')"Grid size:",nx,nz

!-----------------------------------------------------------------------
! Volume a spherical shell.
!-----------------------------------------------------------------------
	Br=1d0
	Bt=0d0
	Bf=0d0
	P=0d0
	call ENERGY_VOLUME(P,rd,td,Br,Bt,Bf,nx,nz,Epol,Etor,E,H)
	write(*,*)
	write(*,'(a)')
     &	"ENERGY_VOLUME: Testing for volume of a spherical shell."
	write(*,'(a24,e15.6)')"Relative error:         ",
     &	abs(6d0*E/(rd(nx)**3-rd(1)**3)-1d0)

!-----------------------------------------------------------------------
! Vacuum energy.
!-----------------------------------------------------------------------
	T=0d0
	call ENERGY(Pd,T,rd,td,nx,nz,E,Ealt,Etor,H)
	write(*,*)
	write(*,'(a)')"ENERGY: Testing for vacuum energy."
	write(*,'(a24,2e15.6)')"Numerical (E,Ealt):     ",E,Ealt
	write(*,'(a24,e15.6)')"Analytical (dipole):    ",1d0/3d0
cc	write(*,'(a24,e15.6)')"Analytical (quadrupole):",6d0/5d0

!-----------------------------------------------------------------------
! Helicity (for a fictitious mathematical construct).
!-----------------------------------------------------------------------
	s=1d0
	Pc=0.5d0
	a(1)=s
	a(2)=Pc
	a(3)=0d0
! Construction of the toroidal function.
! See FTOR and FTOR_NONL for the definition of the toroidal function.
	T=0d0
	do m=1,nx
	   do n=1,nz
	      x=Pd(m,n)
	      call FTOR(x,a,f,dfdx,ma,fit_type)
	      T(m,n)=f
	   enddo
	enddo
! Test.
	call ENERGY(Pd,T,rd,td,nx,nz,E,Ealt,Etor,H)
	write(*,*)
	write(*,'(a)')"ENERGY: Testing for helicity."
	write(*,'(a24,2e15.6)')"Numerical (H):          ",H
	H=(16d0*pi*s*Pc/3d0)*((1d0+2d0*Pc)*sqrt(1-Pc)/Pc
     &	   -3d0*log(1d0+sqrt(1d0-Pc))+3d0*log(Pc)/2d0)
	write(*,'(a24,2e15.6)')"Analytical:             ",H

!-----------------------------------------------------------------------
! Maximum twist (for a fictitious mathematical construct).
!-----------------------------------------------------------------------
	s=1d0
	Pc=0.5d0
	sigma=1.5d0
! Parameters needed by FTOR and FTOR_NONL.
	fn_type=1
	a(1)=s
	a(2)=Pc
	a(3)=sigma
! Test.
	call TWIST_MAX(Pd,Pc,a,ma,fit_type,rd,td,nx,nz,fmax,imax)
	write(*,*)
	write(*,'(a)')"TWIST_MAX: Testing for maximum twist."
	write(*,'(a24,2e15.6)')"Numerical (fmax):       ",fmax
! For sigma=1.
cc	Pmax=(3d0*Pc+2d0)/2d0-sqrt((3d0*Pc+2d0)**2/4d0-4d0*Pc)
cc	fmax=2d0*s*(Pmax-Pc)*sqrt(1d0-Pmax)/Pmax**2
! For general sigma.
	Pmax=(3d0*Pc-2d0*sigma+4d0
     &	   -sqrt((3d0*Pc-2d0*sigma+4d0)**2-16d0*Pc*(3d0-2d0*sigma)))/
     &	   (2d0*(3d0-2d0*sigma))
! Special case for sigma=3/2.
	if(sigma.eq.1.5d0) Pmax=4d0*Pc/(3d0*Pc+1d0)
	fmax=2d0*s*(Pmax-Pc)**sigma*sqrt(1d0-Pmax)/Pmax**2
	write(*,'(a24,2e15.6)')"Analytical:             ",fmax
	write(*,*)

	stop
	end
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
	include "../nr_routines/legendre.f"
	include "../magnetosphere/energy.f"
	include "../magnetosphere/toroidal_fn.f"
	include "../magnetosphere/twist.f"

!-----------------------------------------------------------------------
! Force-free magnetosphere.
! Interface subroutine linking BOUNDARY to FORCE_FREE.
!
! Contents:
! INTERFACE
!-----------------------------------------------------------------------
	subroutine INTERFACE(aphi_surf,bphi_surf,r_surf,t_surf,tbyear,
     &	iwrite,al,bpdip)

	implicit none
	include '../decl/dim2.h'
	include './msph.h'

! Input & output variables.
	integer iwrite
	real*8, dimension(3,nang) :: aphi_surf,bphi_surf
	real*8, dimension(3) :: r_surf
	real*8, dimension (nang):: t_surf
	real*8 tbyear

! Internally used variables.
	integer m,n
	real*8 mu,sinmed
	real*8 Pmax,Tmax,Tmed,Tmed1,Tmed2
	real*8 P,Pc
	real*8 r1,r2,r3
! Normalization constants.
	real*8 Po,Ro
! Limiting value of T (below that it is set equal to zero).
	real*8 Tlim
	real*8 bpdip, al(0:nleg)
	parameter(Tlim=1d-2)

! Grid.
! nrad is defined in dim2.h.
	integer, save :: nx,nz
	real*8, dimension(nrad), save :: rd
	real*8, dimension(nang), save :: td

! Poloidal and toroidal functions.
	real*8, dimension(nrad,nang) :: Pin,Pout,Tin,Tout
	real*8, dimension(nang) :: Psurf,Tsurf

! Multipole expansion.
	integer ell,ellmax,l
	parameter(ellmax=3)
	real*8 dpl(ellmax,nang),DPLGNDR

! Save initial guess for subsequent calls.
	integer, save :: ncall=0
	real*8, dimension(nrad,nang), save :: Pd=0d0

! Linear and non-linear fit for the function T(P) at the surface.
! ma is the number of parameters.
! npc should be equal to or larger than ma (see LFIT or MRQMIN).
	integer ma,npc,n1,n2,nfit
	parameter(ma=3,npc=ma)
	integer ia(ma)
	real*8 alamda,chisq
	real*8 a(ma),f(ma),sig(nang)
	real*8 alpha(npc,npc),covar(npc,npc)
	real*8 ft,dft
	external FTOR_LIN,FTOR_NONL,FTOR_NONL_SIGMA1
! Save the parameters of the toroidal function for subsequent calls.
	real*8, dimension(ma), save :: aini=1d0
! Type of the fitting function (linear=1, non-linear=2).
	integer fit_type
	parameter(fit_type=2)
! Number of iterations for the non-linear fit.
	integer i,niter
	parameter(niter=40)
! Parameters passed on to the fitting function routines.
	integer fn_type
	common /function_type/ fn_type
	common /p_critical/ Pc

! Magnetic field in the magnetosphere.
! Output passed on to BWRITE through a common block.
	real*8, dimension (nang,nrad) :: aphi_msph,bphi_msph
	real*8 Tfit(nang)
	common /magnetosphere/ aphi_msph,bphi_msph,Psurf,Tfit

! Derived quantities.
! Dipole content, energy, helicity and twist.
! Additional quantitites passed from FORCE_FREE through a common block.
	integer imax,kiter
	real*8 a1,E,Ealt,Etor,H,fmax,corr
	common /force_free_output/ a1,E,Ealt,Etor,H,fmax,corr,imax,kiter


! Derivative of the energy.
	real*8, save :: Eold=0d0,tbold=0d0
	real*8 dEdt

! Geometric factor for an uneven grid.
	real*8 summ

! Parameter for the function type (see FTOR_LIN and FTOR_NONL).
	fn_type=2

! Definition of the three radii.
! r1=rb(np), r2=rb(np+1) and r3=rb(np+2).
! These now coincide with the first three points of r_msph.
	r1=r_surf(1)
	r2=r_surf(2)
	r3=r_surf(3)

!-----------------------------------------------------------------------
! Section only performed on the first call.
!-----------------------------------------------------------------------
! Initial guess for the first call.
! ncall is updated in a subsequent if statement below.
	if(ncall.eq.0)then

! Grid construction.
! Dimensions.
	   nx=nrad
	   nz=nang

! Spherical coordinates r and t (theta).
! Radial grid.
! Construct a dimensionless radial grid (rescaled by r1).
! In the rescaled dimensionless units, we have r/r1=rin=1.
	   rd=r_msph/r1

! Angular grid.
! td is the same as the array z defined in MDGEOM.
	   do n=1,nz
	      td(n)=t_surf(n)
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
cc		 Pd(m,n)=0.5d0*(1d0-mu**2)*dpl(ell,n)/rd(m)**ell
cc     &		 +0.5d0*(1d0-mu**2)*dpl(ell+1,n)/rd(m)**(ell+1)
cc     &		 +0.5d0*(1d0-mu**2)*dpl(ell+2,n)/rd(m)**(ell+2)
	      enddo
	   enddo
	endif

!-----------------------------------------------------------------------
! Initial guess, renormalization and boundary conditions.
!-----------------------------------------------------------------------
! Provide an initial guess for the magnetosphere.
! Pd is calculated on the first call, and updated on subsequent calls.
	Pin=Pd

! Specify P and T at the surface.
! These are used as boundary conditions.
	Psurf=0d0
	Tsurf=0d0
	do n=1,nz
	   Psurf(n)=r1*sth(n)*aphi_surf(1,n)
	   Tsurf(n)=r1*sth(n)*bphi_surf(1,n)
	enddo

! Renormalize the stream functions.
! P is given in units of Po, R in units of Ro, and T in units of Po/Ro.
! B is then given in units of Po/Ro^2.
	Po=maxval(Psurf)
	Ro=r1
	Psurf=Psurf/Po
	Tsurf=Tsurf*Ro/Po

! Maximum values of the renormalized stream functions.
	Pmax=1d0
	Tmax=maxval(abs(Tsurf))

! Skip best fits when Tsurf is zero.
	if(Tmax.eq.0d0)then
	   a=0d0
	   Pc=10d0
	   goto 30
	endif

! Boundary conditions are passed on as the values at the surface.
! This should be commented out if carrying out an analytic test.
	Pin(1,:)=Psurf(:)

!-----------------------------------------------------------------------
! Crop the toroidal function.
!-----------------------------------------------------------------------
! Find the range of indices where the toroidal function is non-zero.
! Multiple domains are not accounted for.

! Lower index, n1.
	n1=1
	do n=2,nz
	   if(abs(Tsurf(n))/Tmax.gt.Tlim)then
	      n1=n-1
	      goto 10
	   endif
	enddo
10	continue

! Upper index, n2.
	n2=1
	do n=nz-1,1,-1
	   if(abs(Tsurf(n))/Tmax.gt.Tlim)then
	      n2=n+1
	      goto 20
	   endif
	enddo
20	continue

! Explicitly set T equal to zero outside the interval.
	Tsurf(1:n1)=0d0
	Tsurf(n2:nz)=0d0

! Critical value of the renormalized poloidal function.
cc	Pc=max(Psurf(n1),0.1d0)
cc	Pc=Psurf(n1)
	Pc=min(Psurf(n1),Psurf(n2))

!-----------------------------------------------------------------------
! Best fit for the toroidal function.
!-----------------------------------------------------------------------
	sig=1d0
	covar=0d0

! Linear fit.
	if(fit_type.eq.1)then

! The parameters that should be kept fixed are indicated by a zero.
	   ia(1)=1
	   ia(2)=1
	   ia(3)=1

! Specify the values of the parameters that are fixed.
	   a(1)=1d0
	   a(2)=1d0
	   a(3)=1d0

	   nfit=n2-n1+1
	   call LFIT(Psurf(n1:n2),Tsurf(n1:n2),sig,nfit,
     &	   a,ia,ma,covar,npc,chisq,FTOR_LIN)

! Non-linear fit.
	else if(fit_type.eq.2)then
	   alamda=-1d0

! The coefficients that should be kept fixed are indicated by a zero.
	   ia(1)=1
	   ia(2)=1
	   ia(3)=1

! Initial guess.
	   if(ncall.eq.0)then
! On first call set a(2)=Pc.
	      a(1)=1d0
	      a(2)=Pc
	      a(3)=1d0
	   else
! On subsequent calls use the previous result as an initial guess.
	      a=aini
	   endif

	   nfit=n2-n1-1
	   do i=1,niter
! Fit for the points excluding P=Pc.
cc	      call MRQMIN(Psurf(n1+1:n2-1),Tsurf(n1+1:n2-1),sig,nfit,
cc     &	      a,ia,ma,covar,alpha,npc,chisq,FTOR_NONL,alamda)
! Fit for all points.
	      call MRQMIN(Psurf,Tsurf,sig,nz,
     &	      a,ia,ma,covar,alpha,npc,chisq,FTOR_NONL,alamda)
	   enddo

! Special case for the non-linear function of fn_type=1 (cf. FTOR_NONL).
! If sigma is less than 1, fix sigma=1 and do a linear fit for s.
cc	   if(a(3).lt.1d0)then
cc	      ia(3)=0
cc	      a(3)=1d0

cc	      call LFIT(Psurf(n1+1:n2-1),Tsurf(n1+1:n2-1),sig,nfit,
cc     &	      a,ia,ma,covar,npc,chisq,FTOR_NONL_SIGMA1)
cc	   endif

! Update the value of Pc for subsequent use.
	   Pc=a(2)

! Update initial guess.
	   aini=a

! Error.
	else
	   write(*,'(a)')"INTERFACE: Unrecognized value of fit_type."
	   stop
	endif

! Skip to here when Tsurf is zero.
30	continue

!-----------------------------------------------------------------------
! Solution of the Grad-Shafranov equation for the magnetosphere.
!-----------------------------------------------------------------------
! General force-free in spherical coordinates.
	call FORCE_FREE(Pin,Pout,Tout,rd,td,geo,nx,nz,Pc,a,ma,fit_type,nleg,al)

!-----------------------------------------------------------------------
! Derived quantities.
!-----------------------------------------------------------------------
! Dipole content, energy, helicity and maximum twist are calculated by
! FORCE-FREE and passed on through a common block.

! Twist for a set of magnetic field lines.
cc	if(iwrite.eq.1)then
cc	   call TWIST_DATA(Pout,Pc,a,ma,fit_type,rd,td,nx,nz,tbyear,Po,
cc     &	   Ro)
cc	endif

! Dipole test for the twist (cf. TWIST_DATA).
! Pd should be the dipole field. (True only on the first call.)
cc	a(1)=1d0
cc	a(2)=Pc
cc	a(3)=0d0
cc	call TWIST_DATA(Pd,Pc,a,ma,fit_type,rd,td,nx,nz,tbyear,Po,Ro)
cc	stop

!-----------------------------------------------------------------------
! Update, output and renormalization.
!-----------------------------------------------------------------------
! Update the initial guess for the next call.
	Pd=Pout

! Update the values of aphi and bphi returned to BOUNDARY.
! Normalize the output.
	aphi_surf(2:3,:)=0d0
	bphi_surf(2:3,:)=0d0
	do m=2,3
	   do n=2,nz-1
	      aphi_surf(m,n)=Po*Pout(m,n)/(r_msph(m)*sin(td(n)))
	      bphi_surf(m,n)=(Po/Ro)*Tout(m,n)/(r_msph(m)*sin(td(n)))
	   enddo
	enddo

! aphi and bphi for the output written by BWRITE.
! The angular and radial indices are reversed.
	aphi_msph=0d0
	do m=1,nx
	   do n=2,nz-1
	      aphi_msph(n,m)=Po*Pout(m,n)/(r_msph(m)*sin(td(n)))
	      bphi_msph(n,m)=(Po/Ro)*Tout(m,n)/(r_msph(m)*sin(td(n)))
	   enddo
	enddo

! Tfit for the output written by BWRITE.
	Tfit=0d0
	do n=1,nz
	   P=Psurf(n)
	   if(P.gt.Pc)then
	      call FTOR(P,a,ft,dft,ma,fit_type)
	      Tfit(n)=ft
	   endif
	enddo

! Renormalize for the output.
	Psurf=Po*Psurf
	Tfit=(Po/Ro)*Tfit

!-----------------------------------------------------------------------
! INTERFACE output.
!-----------------------------------------------------------------------
! Renormalize the output calculated above.
! As noted above, B is given in units of Po/Ro^2.
! E (B^2 R^3) is given in units of Po^2/Ro.
! H (A B R^3) is given in units of Po^2.
! (B is in 10^12 G, and R is in km.)
	E=E*Po**2/Ro
	Ealt=Ealt*Po**2/Ro
	Etor=Etor*Po**2/Ro
	H=H*Po**2

! Derivative of the energy.
	if(tbyear.eq.tbold)then
	   dEdt=0d0
	else
	   dEdt=(E-Eold)/(tbyear-tbold)
	endif
	Eold=E
	tbold=tbyear

! Relative dipole content at the surface.
	a1=a1*rout/rin

! Amplitude of the dipolar component at the pole (passed on to MAIN).
	bpdip=2d0*a1*Po/Ro**2

! Create a new file on the first call.
	if(ncall.eq.0)then
	   open(2,file="outb/interface.d")

	   write(2,'(a)')"# -------------------------------------------"
	   write(2,'(a)')"# Output written by INTERFACE."
	   write(2,'(a)')"#"
	   write(2,'(a)')"# Parameters:"
	   write(2,'(a22,2i10)')"# Grid (nrad,kmax):   ",nrad,kmax
	   write(2,'(a22,2f10.2)')"# rin, rout:          ",rin,rout
	   write(2,'(a22,i10)')"# ma:                 ",ma
	   write(2,'(a22,i10)')"# fit_type:           ",fit_type
	   write(2,'(a22,i10)')"# fn_type:            ",fn_type
	   write(2,'(a22,f10.2)')"# Tlim:                ",Tlim
	   write(2,'(a)')"# -------------------------------------------"

	   write(2,100)"#","I","II","III","IV","V","VI","VII","VIII",
     &	   "IX","X","XI","XII","XIII","XIV","XV","XVI"
	   write(2,100)"#","tbyear","a1","E","Ealt","Etor","H","fmax",
     &	   "corr","imax","kiter","a_1","a_2","a_3","chisq","Po","dEdt"

! Update ncall.
	   ncall=1

! Append existing file on subsequent calls.
	else
	   open(2,file="outb/interface.d",access="append")
	endif

	write(2,200)tbyear,a1,E,Ealt,Etor,H,fmax,corr,imax,kiter,a,
     &	chisq,Po,dEdt
	close(2)

100	format(a1,a13,7a14,2a6,6a14)
200	format(e14.7,7e14.6,2i6,6e14.6)

	return
	end

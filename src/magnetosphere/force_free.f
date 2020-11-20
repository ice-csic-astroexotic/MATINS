!-----------------------------------------------------------------------
! Solution of the Grad-Shafranov equation for force-free magnetospheres.
! Imported from gs_v3.2.
!-----------------------------------------------------------------------
	subroutine FORCE_FREE(Pin,Pout,Tout,rd,td,geo,nx,nz,Pc,at,ma,
     &	fit_type,nleg,al)

	implicit none

! Input variables.
	integer nx,nz,nleg
	real*8 Pin(nx,nz),rd(nx),td(nz)
! Geometric factor for an uneven grid.
	real*8 geo

! Output variables.
	real*8 Pout(nx,nz),Tout(nx,nz)

! Dimension indices.
	integer m,n

! Radial and angular arrays, and related coordinates.
	real*8 cth(nz),sth(nz)
	real*8 dr,dt,mu

! Poloidal function.
	real*8 Pold(nx,nz),Pnew(nx,nz)

! Number of multipoles and related variables.
	integer l,lmax
	parameter(lmax=10)
	real*8 al(0:nleg),dpl(lmax,nz)

! External functions.
	real*8 DPLGNDR

! Parameters of the magnetosphere model.
	real*8 Pc
! Number of iterations.
	integer k,kiter,niter
	parameter(niter=40)
! Maximum value for the acceptable correction.
	real*8 corr,corrmax
	parameter(corrmax=1d-6)

! Toroidal function.
	integer fit_type,ma
	real*8 at(ma),dft,ft,x

! Dipole content, energy, helicity and twist.
	integer imax
	real*8 a1,E,Ealt,Etor,H,fmax

! Additional quantitites passed on to INTERFACE through a common block.
	common /force_free_output/ a1,E,Ealt,Etor,H,fmax,corr,imax,kiter

! Relativistic case.
	real*8 zstar
	parameter(zstar=0d0)

! Source (in the Grad-Shafranov equation).
	real*8 factor
	parameter(factor=0d0)

!-----------------------------------------------------------------------
! Uniform grid spacing is assumed!
	dr=rd(2)-rd(1)
	dt=td(2)-td(1)

! Trigonometric functions.
	do n=1,nz
	   cth(n)=cos(td(n))
	   sth(n)=sin(td(n))
	enddo

! Legendre polynomials.
	do l=1,lmax
	   do n=2,nz-1
	      mu=cth(n)
	      dpl(l,n)=DPLGNDR(l,0,mu)
	   enddo
! Derivatives at endpoints. (The axis is not implemented in DPLGNDR.)
	   dpl(l,1)=dble(l*(l+1))/2d0
	   dpl(l,nz)=(-1)**(l+1)*dpl(l,1)
	enddo

! Starting guess.
	Pold=Pin

!-----------------------------------------------------------------------
! Start of iterations.
	do k=1,niter

! Set up the operator and source, define the boundary conditions, and
! solve the linear system.
	   call OPERATOR(Pold,Pnew,rd,td,cth,sth,dr,dt,geo,nx,nz,         
     &	      factor,zstar,corr,dpl,lmax,Pc,at,ma,fit_type)

! Check for convergence at each iteration.
	   if(corr.le.corrmax)then
	      goto 10
	   endif

! Update guess (including some relaxation factor).
	   Pold=Pnew
cc	   Pold=Pold+(Pnew-Pold)*0.3d0

! End of iterations.
	enddo

! Stop after failure of convergence.
	write(*,*)"FORCE-FREE: Increase the number of iterations!"
	stop

!-----------------------------------------------------------------------
! Continue after convergence.
10	continue

! Return after successful completion.
! Poloidal function.
	Pout=Pnew

! Toroidal function.
	Tout=0d0
	do m=1,nx
	   do n=1,nz
	      x=Pout(m,n)
	      if(x.gt.Pc)then
		 call FTOR(x,at,ft,dft,ma,fit_type)
		 Tout(m,n)=ft
	      endif
	   enddo
	enddo

! Check if the toroidal field has reached the outer radius.
	do n=1,nz
	   if(Tout(nx,n).ne.0d0)then
	      write(*,*)"FORCE-FREE: Outer radius reached!"
	      stop
	   endif
	enddo

! Number of iterations.
	kiter=k

!-----------------------------------------------------------------------
! Warning: twist and energy not adapted to uneven grids!
!-----------------------------------------------------------------------
! Maximum twist.
! imax is returned as -1 when the subroutine fails.
! TWIST_MAX will give an error if the toroidal field is zero everywhere.
	fmax=0d0
	imax=-1d0
cc	call TWIST_MAX(Pout,Pc,at,ma,fit_type,rd,td,nx,nz,fmax,imax)

! Energy and helicity.
	E=0d0
	Ealt=0d0
	Etor=0d0
	H=0d0
cc	call ENERGY(Pout,Tout,rd,td,nx,nz,E,Ealt,Etor,H)

! Multipole expansion and dipole content.
	al=0d0
	a1=0d0
	call MULTIPOLE(Pout(nx,:),cth,dpl,nleg,nz,al)
	a1=al(1)

	return
	end

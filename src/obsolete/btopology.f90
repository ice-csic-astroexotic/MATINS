!> @brief This is an auxiliary routine to obtain 
!>        the initial magnetic field configuration
!>
!> @param[in] phi     Poloidal field-related potential
!> @param[in] psi     Toroidal field-related potential
!!
!! Uses:
!! funa
!! getmu
!!-----------------------------------------------------------------------
subroutine btopology(phi, psi)

  use input_params, only : bgeom, bpolin, btorin
  use grid, only: np, sth, cth, rb, aphi, bphi, jcore
  use legpol, only: nleg
  use constants, only: PI
  implicit none

  ! Input variables.
  real*8, dimension(0:np+2,nleg), intent(out) :: phi, psi

  ! Internally used variables.
  integer j, n
  real*8 rcore

  ! Variables used for model 1.
  real*8 mmu,xr,xhat

  ! External function.
  real*8 funa, funda

  ! Variables used for model 2.
  real*8 f2,f4,f6,P,Pc,sigma,x,norm

  real*8, dimension(nleg) :: initial_multipoles_phi, initial_multipoles_psi

  ! -----------------------------------------------------------------------
  ! Initialize.
  aphi = 0d0
  bphi = 0d0
  phi = 0d0
  psi = 0d0

  initial_multipoles_phi = 0d0
  initial_multipoles_psi = 0d0
! Radii of the core and star (as defined in InpUT).
  rcore=rb(jcore)
  rb(np)=rb(np)
!-----------------------------------------------------------------------
! Choice 1: Purely crustal field.
!
! Force-free poloidal field. (cf. Aguilera et al. 2008)
! The toroidal field is confined, but is not force-free.
!-----------------------------------------------------------------------
  if (bgeom == 1) then
! Poloidal field (force-free solution).
    call getmu(rb(np),rcore,mmu)
    xr=mmu*rb(np)
    do j=jcore+1,np+2
      xhat=mmu*rb(j)
      n=1               !! Only dipole initially
      phi(i,n) = rb(i)*funa(xhat,xr)
      n=2              !! Only quadrupole initially
      psi(i,n) = - (rb(i)-rb(jcore))*(rb(i)-rb(np))**2    
    enddo

    n=1
    fac = n*(n+1)*dsqrt((2.d0*n+1.d0)/(4.d0*PI))
    norm = bpolin*rb(np)/funa(xrb(np),xrb(np))/fac   !normalization constant
    phi = norm*phi
    psi = maxval(psi)*psi
  
    call potentials_to_b(phi,psi)


!-----------------------------------------------------------------------
! Choice 2: Dipolar plus toroidal field for a non-barotropic star.
! (cf. Akgun et al. 2013)
!-----------------------------------------------------------------------
  else if (bgeom == 0) then

    ! The field becomes current-free at the surface rb(np)
    ! Polynomial profile.
    f2=(35d0/8d0)/rb(np)**3
    f4=(-21d0/4d0)/rb(np)**5
    f6=(15d0/8d0)/rb(np)**7
    ! Define the vector potential aphi.
    ! Polynomial for non-barotropic star up to rb(2*lvac).
    do j=0,np
      x=rb(j)
      aphi(:,j)=(f2*x+f4*x**3+f6*x**5)*sth(:)
    enddo
    ! Vacuum poloidal field beyond rb(np) .
    aphi(:,np+1)=sth(:)/rb(np+1)**2
    aphi(:,np+2)=sth(:)/rb(np+2)**2

      ! Boundary condition at the center.
  ! (jcenter should be the same as in BEVOL.)
  ! TODO: understand this
    do j=0,jcore-1
      aphi(:,j)=aphi(3,jcore)/(sth(3)*rb(jcore))*rb(j)*sth(:)
    enddo

! Define the toroidal field bphi.
! Toroidal field confined within a critical field line.
!	   Pc=0.5d0/xc
!	   sigma=1d0
!	   do i=0,nang+1
!	      do j=0,np+2
!		 P=aphi(i,j)*rb(j)*sth(i)
!		 if(P.gt.Pc)then
!		    bphi(i,j)=(P-Pc)**sigma/(rb(j)*sth(i))
!		 endif
!	      enddo
!	   enddo

! Arbitrary initial toroidal field.
    rcore=rb(jcore)
    rb(np)=rb(np)
    do j=jcore,np
! Dipole.
      bphi(:,j)=(rb(np)-rb(j))**2*(rb(j)-rcore)**2*sth(:)
! Quadrupole.
!      bphi(:,j)=-(rb(np)-rb(j))**2*(rb(j)-rcore)**2*sth(:)*cth(:)/rb(j)
    enddo
  else

    write(*,'(a)')"BINIT_OPTIONS: Invalid value of bgeom (set 1 or 0)!"
    stop

  endif

  end subroutine btopology

!-----------------------------------------------------------------------
! Calculates the Bessel function A(x), linear combination of Bessel functions
! as defined in equation (8) of Aguilera et al. (2008).
! where x=mmu*r (variable) and xr=mmu*rb(np)
!-----------------------------------------------------------------------
    function funa(x,xr)
      implicit none
      real*8 x,xr,j1,n1,funa
! Spherical Bessel functions j1 and n1 (fist and second kind).
      j1=dsin(x)/x**2-dcos(x)/x
      n1=-dcos(x)/x**2-dsin(x)/x
      funa=x*(j1+tan(xr)*n1)
      return
    end function funa

    real*8 function funda(x,xR)
    implicit none
    real*8  x,xR,a,b,dj1,dn1
    a=1.d0
    b=dsin(xR)/dcos(xR)
    dj1=dsin(x)-dsin(x)/x**2+dcos(x)/x
    dn1=-dcos(x)+dcos(x)/x**2+dsin(x)/x
    funda = a*dj1 + b*dn1 
    return
    end function funda

!-----------------------------------------------------------------------
!      Uses the Newton-Raphson method to solve the equation 
!      sin(mu*(R_core-R_ns)) - mu*R_core*cos(mu*(R_core-R_ns)) = 0
!      to obtain mu
!     (equation 10 of Aguilera et al. 2008) 
!-----------------------------------------------------------------------
    subroutine getmu(rsurface,rcore,mu)
      implicit none
      real*8 rsurface,rcore,mu
      real*8 dr,dfun,fun
        
      dr = rcore - rsurface
      ! Initial guess, this works for a mass between 1.10 and 1.76.
      mu=2.5d0
      fun = 1.d0
      do while (dabs(fun) > 1d-5)
        fun=dsin(mu*dr)-mu*rcore*dcos(mu*dr)
        dfun=dr*dcos(mu*dr)-rcore*dcos(mu*dr)+mu*rcore*dr*dsin(mu*dr)
        mu=mu-fun/dfun
      enddo
    end subroutine getmu

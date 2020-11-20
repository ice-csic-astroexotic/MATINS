!-------------------------------------------------------------------------------
! Magneto Thermal 2D
!-------------------------------------------------------------------------------
! Module: Initialization of magnetic field
!
!> @author
!> Daniele Viganò
!-------------------------------------------------------------------------------

module initial_magnetic

  ! Module imports -------------------------------------------------------------
  use input_params, only: bpolin, magnetic_advance_method
  use grid, only: kmax, lmax, nang, np, cth, jmin, jevol
  use grid, only: aphi, br, bth, bphi, bm, bmed, benu
  use grid, only: jr, jth, jphi, er, eth, ephi, theta
  use legpol, only: getbl, nleg, blout
  use magnetic_evolution, only: curl_phi, curl_pol, compute_Epol, compute_Etor
  use magnetic_evolution, only: magnetic_bc

  implicit none

  real*8, dimension(:,:), allocatable, save :: phi, psi, phiold, psiold
  real*8, dimension(:), allocatable, save :: bcg

  contains

    !---------------------------------------------------------------------------
    !> Magnetic field initialization.
    !> @brief It defines the initial topology and strength
    !
    !>  Code owners:
    !>  Daniele Viganò
    !---------------------------------------------------------------------------
    subroutine binit

      implicit none

      ! Internally used variables.
      integer i, k, l
    
      ! Auxiliary fields for inermediate calculations.
      real*8, dimension (0:nang+1,0:np+2) :: xr,xth,xphi
    
      allocate(phi(0:np + 2,nleg))
      allocate(psi(0:np + 2,nleg))
      allocate(phiold(0:np + 2,nleg))
      allocate(psiold(0:np + 2,nleg))
      allocate(bcg(0:3*nleg))

      ! Glebsch-Gordon-related coefficients used in spectral methods only
      if (magnetic_advance_method == 'SPEC') then
        bcg(0) = 1.d0
        do i=1,3*nleg
          bcg(i) = (1.d0-0.5d0/dble(i))*bcg(i-1)
        enddo
      endif

      !-----------------------------------------------------------------------
      ! Initial magnetic field. (Returns aphi and bphi at all grid points.)
      !-----------------------------------------------------------------------
      call btopology
      !-----------------------------------------------------------------------
      ! Calculate Br and Btheta from A_phi
      call curl_phi(aphi,jmin,br,bth)

      ! Apply boundary conditions and
      call magnetic_bc(aphi,br,bth,bphi)

      ! Calculation of new modulus and its average at the center of the thermal cell
      bm=dsqrt(br**2+bth**2+bphi**2)
      do k=2,kmax
        do l=1,lmax
          bmed(k,l) = 0.25d0*bm(2*k-2,2*l-1) + &
    &       0.125*(bm(2*k-1,2*l-1) + bm(2*k-3,2*l-1) + bm(2*k-2,2*l) + bm(2*k-2,2*l-2)) + &
    &       0.0625*(bm(2*k-1,2*l-2) + bm(2*k-1,2*l) + bm(2*k-3,2*l-2) + bm(2*k-3,2*l-2)) 
        enddo
      enddo

      ! Calculate initial electric current
      xr   = 0d0
      xth  = 0d0
      xphi = 0d0
      do i=0,nang+1
        xr(i,jmin-1:) = br(i,jmin-1:)*benu
        xth(i,jmin-1:) = bth(i,jmin-1:)*benu
        xphi(i,jmin-1:) = bphi(i,jmin-1:)*benu
      enddo
      call curl_pol(xr,xth,jmin,jphi)
      call curl_phi(xphi,jmin,jr,jth)

      ! Calculate initial electric field
      call compute_Epol(br,bth,bphi,jr,jth,jphi,er,eth)
      call compute_Etor(br,bth,jr,jth,jphi,ephi)
            
      ! Screen output.
      write(*,'(a,1pe12.2)')"BINIT: Poloidal field at the polar surface [G]:", -2d0*blout(1)*1.e12
      write(*,'(a,1pe12.2)')"BINIT: Toroidal field maximum [G]:             ",maxval(dabs(bphi))*1.e12
      write(*,*)
    
    end subroutine binit


    !!-----------------------------------------------------------------------
    !> @brief This defines the topology of the initial magnetic field configuration
    !>
    !! Uses:
    !! funa
    !! getmu
    !!-----------------------------------------------------------------------
    subroutine btopology
    
      use input_params, only : bgeom, bpolin, btorin
      use input_params, only : n_initial_multipoles_phi, initial_multipoles_phi
      use input_params, only : n_initial_multipoles_psi, initial_multipoles_psi
      use grid, only: np, sth, cth, rb, aphi, bphi, jcore
      use legpol, only: nleg, frel
      use constants, only: PI

      implicit none
        
      ! Internally used variables.
      integer i, j, n

      ! Variables used for model 1.
      real*8 mu
      real*8, dimension(nleg) :: multipoles_phi, multipoles_psi
    
      ! Variables used for model 2.
      real*8 f2,f4,f6,P,Pc,sigma
         
    
      ! -----------------------------------------------------------------------
      ! Initialize.
      aphi = 0d0
      bphi = 0d0
      phi = 0d0
      psi = 0d0

      multipoles_phi = 0d0
      multipoles_psi = 0d0

      !-----------------------------------------------------------------------
      ! Choice 1: Purely crustal field.
      ! Modify the initial multipolar weights
      ! Radial function for phi: matching the potential solution outside
      ! The toroidal field is confined, but is not force-free.
      !-----------------------------------------------------------------------

      if (bgeom == 1) then

        ! jevol is the first index to be evolved in the magnetic evolution
        jevol = jcore + 1
        ! jmin is the first radial index to be calculated by curl operators for currents
        jmin = jcore - 1

        ! Define the initial multipolar weights of phi and psi
        multipoles_phi(:n_initial_multipoles_phi) = initial_multipoles_phi(:)
        multipoles_psi(:n_initial_multipoles_psi) = initial_multipoles_psi(:)

        call getmu(rb(np),rb(jcore),mu)
        do j=jcore+1,np
          phi(j,:) = multipoles_phi*funa(mu*rb(j),mu*rb(np))
          psi(j,:) = -multipoles_psi*(rb(np)-rb(j))**2*(rb(j)-rb(jcore))**2
        enddo

        do n=1,nleg
          phi(np+1,n) = phi(np,n) - (rb(np+1)-rb(np))*n*frel(n)*phi(np,n)/rb(np)
          phi(np+2,n) = phi(np+1,n) - (rb(np+2)-rb(np+1))*n*frel(n)*phi(np+1,n)/rb(np+1)
          psi(np+1:np+2,n) = 0d0
        end do        

        phi = bpolin*rb(np)**2/(2d0*phi(np,1))*phi

        call potentials_to_b

        if (maxval(abs(bphi)) /= 0) then
          psi = btorin*psi/maxval(abs(bphi))
          bphi = btorin*bphi/maxval(abs(bphi))
        end if

      !-----------------------------------------------------------------------
      ! Choice 2: Core-extended field
      ! Dipolar plus toroidal field for a non-barotropic star.
      ! (cf. Akgun et al. 2013)
      !-----------------------------------------------------------------------
      else if (bgeom == 2) then

        jevol = 3
        jmin = 1
  
        ! The field becomes current-free at the surface rb(np)
        ! Polynomial profile.
        f2=(35d0/8d0)/rb(np)**3
        f4=(-21d0/4d0)/rb(np)**5
        f6=(15d0/8d0)/rb(np)**7
        ! Define the vector potential aphi.
        ! Polynomial for non-barotropic star up to rb(2*lvac).
        do j=0,np
          aphi(:,j)=(f2*rb(j)+f4*rb(j)**3+f6*rb(j)**5)*sth(:)
        enddo
        ! Vacuum poloidal field beyond rb(np) .
        aphi(:,np+1)=sth(:)/rb(np+1)**2
        aphi(:,np+2)=sth(:)/rb(np+2)**2
    
        ! Normalization valid for the dipole
        aphi = aphi*0.5d0*bpolin*rb(np)/aphi((nang+1)/2,np)
    
        ! Define the toroidal field bphi.
        ! Toroidal field confined within a critical field line.
        Pc=aphi((nang+1)/2,np)*rb(np)
        sigma=1d0
        do i=0,nang+1
          do j=0,np+2
            P = aphi(i,j)*rb(j)*sth(i)
            if (P > Pc) then
              bphi(i,j)=(P-Pc)**sigma/(rb(j)*sth(i))
            endif
          end do
        end do
        bphi = bphi*btorin/maxval(dabs(bphi))
    
      else
    
        write(*,'(a)')"BINIT_OPTIONS: Invalid value of bgeom: choose 1 or 2 (in/input.dat)"
        stop
    
      endif

    end subroutine btopology
    
    !-----------------------------------------------------------------------
    ! Calculates the Bessel function A(x), linear combination of Bessel functions
    ! as defined in equation (8) of Aguilera et al. (2008).
    ! where x=mu*r (variable) and xr=mu*rb(np)
    ! It allows a smooth matching (function and derivative)
    ! with a non-relativistic dipole only
    !-----------------------------------------------------------------------
    real*8 function funa(x,xr)
      implicit none
      real*8 x,xr,j1,n1
      
      ! Spherical Bessel functions j1 and n1 (fist and second kind).
      j1=dsin(x)/x**2-dcos(x)/x
      n1=-dcos(x)/x**2-dsin(x)/x
      funa=x*(j1+tan(xr)*n1)
      return
    end function funa
        
    !-----------------------------------------------------------------------
    !> @brief Definition of the radial function of poloidal field
    !         for crust-confined field, matching with potential solution
    !
    ! It uses the Newton-Raphson method to solve the equation 
    ! sin(mu*(R_core-R_ns)) - mu*R_core*cos(mu*(R_core-R_ns)) = 0
    ! to obtain mu (eq.10 of Aguilera et al. 2008, A&A) 
    !-----------------------------------------------------------------------
    subroutine getmu(rsurface,rcore,mu)
      implicit none
      real*8 rsurface,rcore,mu
      real*8 dr,dfun,fun
        
      dr = rcore - rsurface
      ! Initial guess, this works for a mass between 1.10 and 1.76 for sly4 EoS.
      mu=2.5d0
      fun = 1.d0
      do while (dabs(fun) > 1d-5)
        fun=dsin(mu*dr)-mu*rcore*dcos(mu*dr)
        dfun=dr*dcos(mu*dr)-rcore*dcos(mu*dr)+mu*rcore*dr*dsin(mu*dr)
        mu=mu-fun/dfun
      enddo
    end subroutine getmu


    !-----------------------------------------------------------------------
    !> @brief Conversion from potentials phi and psi to Aphi and B components
    !
    !-----------------------------------------------------------------------
    subroutine potentials_to_b

      use grid, only: rb
      use legpol, only: nleg, dpln, blout
      use constants, only: PI

      implicit none
    
      integer :: j,n

      aphi = 0d0
      br = 0d0
      bth = 0d0
      bphi = 0d0

      do n = 1,nleg
        do j = 1,np+2
          aphi(1:nang,j) = aphi(1:nang,j) - dpln(1:nang,n)*phi(j,n)/rb(j)
          bphi(1:nang,j) = bphi(1:nang,j) - dpln(1:nang,n)*psi(j,n)/rb(j)
        end do
        blout(n) = - n*phi(np,n)/rb(np)**2
      end do
      aphi(0,:) = - aphi(2,:)
      aphi(nang+1,:) = - aphi(nang-1,:)

      ! Curl(Aphi) to obtain the poloidal magnetic field components
      call curl_phi(aphi,jmin,br,bth)

    end subroutine potentials_to_b
  
end module initial_magnetic

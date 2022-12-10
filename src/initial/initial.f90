!-------------------------------------------------------------------------------
! initial conditions
!@brief In this subroutine, we are defining the initial conditions for the magnetic 
! the thermal parts. 
! 
! To define the initial magnetic topology, we use select case
!>          It calls one of the following subroutines:
!>              bessel 
!>              test_pa
!>              hotspot
!>              pureTQ
!>              DP+TQ
!>              ScalarFunc
!>              Whister
! and they return the initial topology of the magnetic field in cubed-sphere coordinates
!
!
!> @author
!>  Clara Dehman
!>  Stefano Ascenzi
!>  Daniele Viganò
!-------------------------------------------------------------------------------


subroutine initial_condition()

  use input_params
  use constants, only: PI
  use grid, only: nr, nang, nrt, nangt
  use grid, only: r, theta, phi, tem0, temp, temp_surf, T_core
  use grid, only: br, beta, bxi, b2, bm, lmax
  use grid, only: jr, jeta, jxi, j2
  use grid, only: er, exi, eeta
  use grid, only: ievol, enu, etab
  use grid, only: f_spherical_to_cs
  use grid, only: curl_fnvol, fghost
  use grid, only: dot_prod
  use initial_magnetic, only: binit, binit_axi, DP_TQ, whistler, pureTQ, bessel, test_pa, hotspot
  use magnetic_evolution, only: compute_E, compute_joule

  implicit none

  real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: xr, xxi, xeta, bth, bphi
  

  integer p, j, k
  
  ! Temperature definition
  temp = T_init     ! Uniform temperature
  tem0 = T_init
  T_core = T_init   ! Temperature of the core

  !********* non uniform temperature (for debugging) *********
  !do j = 1, nangt
  !  do k = 1,nangt
  !    do p = 1,6
  !      temp(:,j,k,p) = T_init*(1.d0 + dcos(theta(2*j, 2*k, p)**2))
  !    enddo
  !  enddo
  !enddo
  !T_core = 1.5*T_init
  ! ********** END ******************

  call envelope_model()

  !   Choose the initial magnetic topology.
  if (bpolmax == 0. .and. btormax == 0.) then
    br = 0d0
    bxi = 0d0
    beta = 0d0
    b2 = 0d0
    bm = 0d0
    jr = 0d0
    jxi = 0d0
    jeta = 0d0
    er = 0d0
    exi = 0d0
    eeta = 0d0    
    j2 = 0d0

  else

    select case (init_mag_top)
      case("Bessel")
        call bessel(bpolmax,btormax)
      case ("PerAzo") 
        call test_pa
      case ("Hotspot") 
        call hotspot
      case ("PureTQ")
        call pureTQ(btormax)
      case ("DP+TQ")
        call  DP_TQ(bpolmax,btormax)
      case ("Scalar")
        call  binit(bpolmax,btormax)
      case ("Scalaxi")
        call  binit_axi(bpolmax,btormax)
      case ("Whistler")
        call  whistler(bpolmax,btormax)
      case default
        write(*,*) "<error>", &
                 & "[B_initial]", &
                 & "Invalid initial magnetic topology: ", &
                 & init_mag_top
      stop
    end select
  
   ! Calculation of electric currents
    xr   = 0d0
    xxi = 0d0
    xeta  = 0d0
    do p = 1, 6
      do j = 0, nang+1
        do k = 0, nang+1
          xr(ievol-1:,j,k,p) = br(ievol-1:,j,k,p)*enu(ievol-1:)
          xxi(ievol-1:, j,k,p) = bxi(ievol-1:,j,k,p)*enu(ievol-1:)
          xeta(ievol-1:,j,k,p) = beta(ievol-1:,j,k,p)*enu(ievol-1:)
        enddo
      enddo
    enddo
  
    ! Calculation of electrical currents
    call curl_fnvol(xr,xxi,xeta,jr,jxi,jeta,ievol)
    call fghost(jr,jxi,jeta)
  
    ! Calculation of electric fields
    call compute_E(br,bxi,beta,jr,jxi,jeta,er,exi,eeta)
    
    ! Calculation of magnetic field intensity
    call dot_prod(br,br,bxi,bxi,beta,beta,b2)
    bm = sqrt(b2)
    ! Calculate J**2 used for the Joule dissipation
    call dot_prod(jr,jr,jxi,jxi,jeta,jeta,j2)
    call compute_joule

  endif

  
end subroutine initial_condition

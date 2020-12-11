!-------------------------------------------------------------------------------
! Magneto Thermal 3D
!-------------------------------------------------------------------------------
! Module: Grid
!
!> @author
! Clara Dehman 
! Jose A. Pons
! Daniele Vigano
!
!> Grid.
!> @brief This module is responsible for handling the grid and structure.
!
!-------------------------------------------------------------------------------
module grid

! Clara Dehman: 
! Ghost cells:
! In the radial direction, we need at least one ghost cell in order to calculate the magnetic field at the surface
! In xi and eta physical directions, our sphere is complete so in reality we don't need the ghost cell but taking into account that we have six patches, 
! then to compute the magnetic field at the surface, we will use a cell from another patch, which is in different coordinate system, 
! thus it is not straight forward because we will need to change the coordinate system
! another alternative could be by introducing additional cells which goes smoothly with the same patch, meaning in the same coordinate system and so we don't need 
! to change the coordinate system when computing the magnetic field at the surface
! However, I don't believe we will win much because we will need to work out the transformation of the coordinate system from one patch to another so then we can direclty use it


! Clara Dehman
! The fourth dimension:
! We have introduce a fourth dimension ipatch, since we are dividing our sphere that represents the NS
! into 6 patches identical in size. Thus ipatch must identify which patch we are dealing with


! Clara Dehman
! I have considered eta and xi physical direction to describe each patch because the infinitesimal displacements are 
! described in term of eta and xi and not X/Y physical direction, then for the rest of the calculation, we can easily  
! find X and Y since X=tan(xi) and Y=tan(eta). Moreover, the transformation matrix between the spherical coordinates 
! system and the eta/xi coordinate system is given in term of eta/xi and not X/Y. 

  ! Modules --------------------------------------------------------------------
  use input_params, only: enable_hall_effect
  use input_params, only: superfluid_n_crust, superfluid_n_core, superfluid_p_core

  ! Module constants -----------------------------------------------------------
  ! None

  ! Module variables -----------------------------------------------------------

  ! Radial and angular dimensions.
  ! Default values are -1 to issue an error if not initialized.
  integer, save :: jmax = -1  
  integer, save :: kmax = -1
  integer, save :: lmax = -1
  ! Cubed sphere grid dimensions in r, eta, and xi physical directions --------- 
  integer, save :: np = -1
  integer, save :: neta = -1 
  integer, save :: nxi = -1 
  !-----------------------------------------------------------------------------
  integer, save :: nleg = -1
  integer, save :: mleg = -1
  integer, save :: ng = -1

  ! ipatch is an index used to indicate in which patches we are working in the cubed sphere grid 
  ! it goes from 1 to 6, which corresponds to the number of patches in our grid
  integer, save :: ipatch

  ! Radial and angular grids. 
  real*8, dimension(:), allocatable, save :: rb ! Radial grid.
  real*8, dimension(:), allocatable, save :: eta ! eta direction.
  real*8, dimension(:), allocatable, save :: xi  ! xi direction.
  !-----------------------------------------------------------------------------
  real*8, dimension(:), allocatable, save :: X  ! tan(xi).
  real*8, dimension(:), allocatable, save :: Y  ! tan(eta).
  !-----------------------------------------------------------------------------
  real*8, dimension(:, :),  allocatable, save :: delta ! delta=1+x^2+y^2
  real*8, dimension(:),  allocatable, save :: C !C=(1+x^2)^{1/2}
  real*8, dimension(:),  allocatable, save :: D !D=(1+y^2)^{1/2}

  real*8, dimension(:), allocatable, save :: ceta ! cos(eta).
  real*8, dimension(:), allocatable, save :: cxi ! cos(xi).
  !-----------------------------------------------------------------------------  
  ! Ghost cells in xi and eta directions (used for interpolation)
  real*8, dimension(:), allocatable, save :: xi_ghost  
  real*8, dimension(:), allocatable, save :: eta_ghost 
  ! Length, area, and volume elements. 
  real*8, dimension(:), allocatable, save :: lr ! radial length element of the cell.
  real*8, dimension(:, :, :), allocatable, save :: lxi ! length element of the cell in xi physical direction.
  real*8, dimension(:, :, :), allocatable, save :: leta ! length element of the cell in eta physical direction.
  real*8, dimension(:, :, :), allocatable, save :: arear ! Area of the radial interface
  real*8, dimension(:, :, :), allocatable, save :: areaxi 
  real*8, dimension(:, :, :), allocatable, save :: areaeta 
  real*8, dimension(:, :, :), allocatable, save :: vol ! Volume of each cell
  ! Covariant components 
  real*8, dimension(:), allocatable, save :: l_r 
  real*8, dimension(:, :, :), allocatable, save :: l_xi
  real*8, dimension(:, :, :), allocatable, save :: l_eta
  real*8, dimension(:, :, :), allocatable, save :: area_r 
  real*8, dimension(:, :, :), allocatable, save :: area_xi 
  real*8, dimension(:, :, :), allocatable, save :: area_eta

  ! Relativistic factors 
  real*8, dimension(:), allocatable, save :: belam  ! e^lambda, length correction
  real*8, dimension(:), allocatable, save :: benu   ! e^nu, lapse function
  ! metric tensor of the cubed-sphere formalism
  real*8, dimension(:, :, :, :, :), allocatable, save :: g 


  ! Star's structure quantities
  ! Radial profiles of:
  !   density (rho)
  !   nuclei fraction (xh)
  !   electron fraction (ye)
  !   neutron fraction (yn)
  !   proton fraction (yp)
  !   atomic mass number (aa)
  !   atomic number (zz)
  real*8, dimension (:), allocatable, save :: rho, xh, ye, yn, yp, aa, zz
  ! Flags to swith on/off the hall and ambipolar terms.
  integer, dimension (:), allocatable, save :: ia_amb, ia_hall
  real*8, save :: g14   ! Surface gravity [10^14 gr*cm/s^2].
  real*8, save :: moment_inertia   ! Moment of Inertia [g cm^2].
  real*8, save :: spindown_prefactor   ! Spin-down torque factor (used in pevol).
  real*8, dimension(:), allocatable, save :: mass, press   ! Mass and pressure.

  ! Microphysics
  real*8, dimension(:, :), allocatable, save :: c_v      ! Specific heat [10^40 erg/(km^3*10^8 K)]
  real*8, dimension(:, :, :, :), allocatable, save :: etab     ! Magnetic diffusivity [km^2/Myr]                 
  real*8, dimension(:), allocatable, save :: fh          ! Hall coefficient c/4*pi*n_e [km^2/(Myr*10^12 G)]
  real*8, dimension(:, :), allocatable, save :: omegatau ! Magnetization parameter [adimensional]
  real*8, dimension(:, :), allocatable, save :: dfc2, dfc3
  real*8, dimension(:, :, :, :), allocatable, save :: anis2, anis3 ! Anisotropy factors
  real*8, dimension(:, :), allocatable, save :: q_neutrino ! Neutrino emissivity [10^40 erg/(km^3*s)], no redshift included
  real*8, dimension(:, :), allocatable, save :: rmc      ! Removal coefficient (implicit part of neutrino emissivity) in units [10^40 erg/(km^3*s*10**8 K)]
  real*8, dimension(:, :), allocatable, save :: gammac   ! Coulomb parameter [adimensional]

  ! Superfluid and superconducting gaps: critical temperatures in units of 10**8 K.
  real*8, dimension(:), allocatable, save :: tcn, tcp, tccru

  ! Temperatures
  ! Redshifted temperature T*e^nu, which is the evolved variable
  real*8, dimension(:, :), allocatable, save :: tem      ! [K]
  ! Physical temperature T, entering in microphysical calculations
  real*8, dimension(:, :), allocatable, save :: tem0     ! [K]
  ! Surface boundary condition for heat diffusion
  real*8, dimension(:), allocatable, save :: tss      ! Surface temperature [K]
  real*8, dimension(:), allocatable, save :: sfluxb   ! Surface thermal flux [10^40 erg/km^2/s]
  real*8, dimension(:), allocatable, save :: cfluxb   ! Implicit component of the surface thermal flux (only for Gudmundsson)

  ! Magnetic evolution. 
  real*8, dimension(:, :, :, :), allocatable, save :: aphi        ! Toroidal vector potential [10^12 G*km] 
  real*8, dimension(:, :, :, :), allocatable, save :: br,bxi,beta! Magnetic field components [10^12 G]  
  real*8, dimension(:, :, :, :), allocatable, save :: bm          ! Modulus of B [10^12 G]
  real*8, dimension(:, :), allocatable, save :: bmed        ! Modulus of B averaged on the thermal cell [10^12 G]


  ! The numerical currents are defined as Jnum = curl [e^nu B] = e^nu*4*pi/c times physical currents:
  ! The numerical electric field is defined as Enum = c*e^nu times physical Electric field
  real*8, dimension(:, :, :, :), allocatable, save :: jr,jxi,jeta ! Electric currents components [10^12 G/km]
  real*8, dimension(:, :, :, :), allocatable, save :: er,exi,eeta  ! Electric field components
  real*8, dimension(:, :, :, :), allocatable, save :: j2          ! Square of the currents J^2 [10^24 (G/km)^2]



! To be modified according  to the output choice -------------------------------
  real*8, dimension(:, :, :, :), allocatable, save :: bth,bphi   ! Magnetic field components [10^12 G] 
  real*8, dimension(:, :, :, :), allocatable, save :: jth,jphi   ! Electric currents components [10^12 G/km] 
  real*8, dimension(:, :, :, :), allocatable, save :: eth,ephi   ! Electric field components
! ------------------------------------------------------------------------------

  ! Radial indexes marking the relevant radial points to start the calculations with
  ! jcore is the last point of the core (set in grid.f90)
  integer, save :: jcore
  ! jevol is the first point to be evolved (set in binit.f90)
  ! jmin minimum value to calculate in the curl operators (set in binit.f90)
  ! jmin = jevol - 2, and jevol depends on the geometry of the initial magnetic field
  integer, save :: jevol
  integer, save :: jmin


  ! Radial profile of Courant magnetic timestep estimation
  real*8, dimension(:), allocatable, save :: dtb_courant_profile  ! [years]

  ! Parameters of the Burgers equation (used in magnetic_evolution module)
  real*8, dimension(:, :), allocatable, save :: lamr  ! Radial part
  real*8, dimension(:), allocatable, save :: lamth    ! Meridional part

  ! Joule dissipation, they include the e^2nu factor, in units [10^40 erg/(km^3*s)]
  real*8, dimension(:, :), allocatable, save :: q_joule          ! Instantaneous dissipation
  real*8, dimension(:, :), allocatable, save :: q_joule_shock    ! Shock correction to the dissipation
  real*8, dimension(:, :), allocatable, save :: q_joule_average  ! Joule rate averaged over dt

  ! Crustal magnetic stress quantities
  real*8, dimension(:, :), allocatable, save :: shearModulus     ! Shear modulus
  real*8, dimension(:, :), allocatable, save :: shearMaximum     ! Maximum shear allowed
  ! Components of the equilibrium value from which stresses are calculated
  real*8, dimension(:, :), allocatable, save :: freq, ftheq, fphieq


  contains

    !---------------------------------------------------------------------------
    !> Allocate grid arrays.
    !> @brief Allocates all dynamically-sized arrays that are needed for the
    !>        grid. The subroutine also performs precondition checks to ensure
    !>        that the angular and radial dimensions have appropriate
    !>        (non-negative) values.
    !---------------------------------------------------------------------------
    subroutine allocate_grid()

      ! Modules ----------------------------------------------------------------
      ! None.

      implicit none

      ! Subroutine arguments ---------------------------------------------------
      ! None.

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
      ! None.

      ! Preconditions ----------------------------------------------------------

      if (jmax <= 0) then
        print *, "<Error>[GRID] Invalid jmax value = ", jmax
        call abort()
      end if

      if (kmax <= 0) then
        print *, "<Error>[GRID] Invalid kmax value = ", kmax
        call abort()
      end if

      if (lmax <= 0) then
        print *, "<Error>[GRID] Invalid LMAX value = ", lmax
        call abort()
      end if

      ! if (nang <= 0) then
      !  print *, "<Error>[GRID] Invalid nang value = ", nang
      !  call abort()
      ! end if

      if (np <= 0) then
        print *, "<Error>[GRID] Invalid np value = ", np
        call abort()
      end if

      if (neta <= 1) then
        print *, "<Error>[GRID] Invalid np value = ", neta
        call abort()
      end if

      if (nxi <= 1) then
        print *, "<Error>[GRID] Invalid np value = ", nxi
        call abort()
      end if
      ! ------------------------------------------------------------------------

      allocate(rb(0:np+2))
      allocate(xi(0:nxi+1))
      allocate(eta(0:neta+1))
      allocate(X(0:nxi+1))
      allocate(Y(0:neta+1))
      allocate(delta(0:nxi+1, 0:neta+1))  
      allocate(C(0:nxi+1))
      allocate(D(0:neta+1))

      allocate(xi_ghost(0:nxi+1))
      allocate(eta_ghost(0:neta+1))

      allocate(g(0:np+2, 0:nxi+1, 0:neta+1, 1:3, 1:3))

      allocate(lr(0:np+2))
      allocate(lxi(0:np+2, 0:nxi+1, 0:neta+1))
      allocate(leta(0:np+2, 0:nxi+1, 0:neta+1))
      allocate(areaxi(0:np+2, 0:nxi+1, 0:neta+1))
      allocate(arear(0:np+2, 0:nxi+1, 0:neta+1))
      allocate(areaeta(0:np+2, 0:nxi+1, 0:neta+1))
      allocate(vol(0:np+2, 0:nxi+1, 0:neta+1))
      
      ! Covariant components 
      allocate(l_r(0:np+2))
      allocate(l_xi(0:np+2, 0:nxi+1, 0:neta+1))
      allocate(l_eta(0:np+2, 0:nxi+1, 0:neta+1))
      allocate(area_xi(0:np+2, 0:nxi+1, 0:neta+1))
      allocate(area_r(0:np+2, 0:nxi+1, 0:neta+1))
      allocate(area_eta(0:np+2, 0:nxi+1, 0:neta+1))

      allocate(belam(0:np+2))
      allocate(benu(0:np+2))

     !--------------------------------------------------------------------------

      allocate(rho(np + 2))
      allocate(xh(np + 2))
      allocate(ye(np + 2))
      allocate(yn(np + 2))
      allocate(yp(np + 2))
      allocate(aa(np + 2))
      allocate(zz(np + 2))

      allocate(fh(0:np + 2))
      allocate(ia_amb(0:np + 2))
      allocate(ia_hall(0:np + 2))
      allocate(dtb_courant_profile(0:np + 2))

      allocate(tcn(np + 2))
      allocate(tcp(np + 2))
      allocate(tccru(np + 2))
      allocate(mass(np + 2))
      allocate(press(np + 2))

     !--------------------------------------------------------------------------

      allocate(br(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(bxi(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(beta(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(bm(0:np+2, 0:nxi+1, 0:neta+1, 1:6))

      allocate(jr(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(jxi(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(jeta(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(j2(0:np+2, 0:nxi+1, 0:neta+1, 1:6))

      allocate(er(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(exi(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(eeta(0:np+2, 0:nxi+1, 0:neta+1, 1:6))

    ! to be modified according to the output choice ---------------------------- 
      allocate(aphi(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(bth(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(bphi(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(jth(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(jphi(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(eth(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(ephi(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
    !---------------------------------------------------------------------------
      allocate(freq(kmax,jcore/2:lmax))
      allocate(ftheq(kmax,jcore/2:lmax))
      allocate(fphieq(kmax,jcore/2:lmax))
      allocate(gammac(kmax,jcore/2:lmax))  

      allocate(etab(0:np+2, 0:nxi+1, 0:neta+1, 1:6))
      allocate(dfc2(kmax, lmax))
      allocate(dfc3(kmax, lmax))
      allocate(anis2(2, 2, kmax, lmax))
      allocate(anis3(2, 2, kmax, lmax))

      allocate(tem(kmax, lmax))
      allocate(tem0(kmax, lmax))
      allocate(q_neutrino(kmax, lmax))
      allocate(q_joule(kmax, lmax))
      allocate(q_joule_shock(kmax, lmax))
      allocate(q_joule_average(kmax, lmax))
      allocate(rmc(kmax, lmax))
      allocate(c_v(kmax, lmax))
      allocate(bmed(kmax, lmax))
      allocate(cfluxb(kmax))
      allocate(sfluxb(kmax))
      allocate(tss(kmax))
      allocate(shearModulus(kmax,jcore/2:lmax))
      allocate(shearMaximum(kmax,jcore/2:lmax))
      allocate(omegatau(kmax, lmax))

    end subroutine allocate_grid

    !---------------------------------------------------------------------------
    !> Grid size setter.
    !> @brief Setter subroutine for the grid angular (xi and eta) and radial dimensions.
    !>        Additionally, the subroutine calculates the magnetic grid
    !>        dimensions (nxi, neta and np).
    !> @param[in] angular_dimension Number of angular cells in xi and eta physical directions.
    !> @param[in] radial_dimension Number of radial cells.
    !---------------------------------------------------------------------------
    subroutine set_grid_size(xi_dimension,eta_dimension,radial_dimension)

      ! Modules ----------------------------------------------------------------
      ! None.

      ! Subroutine arguments ---------------------------------------------------
      integer, intent(in) :: xi_dimension
      integer, intent(in) :: eta_dimension
      integer, intent(in) :: radial_dimension

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
      ! None.

      ! ------------------------------------------------------------------------

      jmax = xi_dimension
      kmax = eta_dimension
      lmax = radial_dimension
      ! Clara Dehman
      ! magnetic grid: In order to conserve the shape of the thermal grid cell at the border of each patch, 
      ! we must consider an even number of cells in the magnetic grid in the xi and eta angular directions
      nxi = 2 * jmax + 1
      neta = 2 * kmax + 1
      np = 2 * lmax
      nleg = nxi/2
      mleg = neta/2
      ng = nleg

    end subroutine

    !---------------------------------------------------------------------------
    !> Grid builder.
    !> Infinitesimal geometrical elements.
    !> @brief Construction the whole grid by initializing the arrays in the 
    !>  cubed-sphere formalism.
    !> We are defining the geometrical elements in one patch one since in the 
    !> other patches we will have the same expressions for the infinitesimal 
    !> geometrical elements
    !> 
    !> code owner
    !> Clara Dehman
    !---------------------------------------------------------------------------
    subroutine build_grid()

      ! Modules ----------------------------------------------------------------
      use constants, only : PI

      implicit none

      ! Subroutine arguments ---------------------------------------------------
      ! None.

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
       integer :: i, j, k
       real*8 :: dxi ! Angular step in xi direction, constant
       real*8 :: deta ! Angular step in eta direction, constant
      ! ------------------------------------------------------------------------

      ! xi/eta grid, no ghost cells --------------------------------------------
      ! The angular variables xi and eta are spanning the range [-Pi/4 : Pi/4] 
       dxi = PI / (2*dble(nxi-1))
       deta = PI / (2*dble(neta-1))

       do j = 0, nxi+1
        xi(j) = dble(j) * dxi
       end do
       do k = 0, neta+1
        eta(k) = dble(k) * deta
       end do
    
       X(:) = dtan(xi(:))
       Y(:) = dtan(eta(:))
       C(:) = (1 + X(:)**2)**(1/2)
       D(:) = (1 + Y(:)**2)**(1/2)

       do j = 0, nxi+1
       do k = 0, neta+1
        delta(j, k) = 1 + X(j)**2 + Y(k)**2
       end do
       end do

      ! Define cos in the interval [1 : neta/nxi] for eta and xi
        cxi(:) = dcos(xi(:))
        ceta(:) = dcos(eta(:))
       ! Contravariant components
        lr = 0d0
        lxi = 0d0
        leta = 0d0
        arear = 0d0
        areaxi = 0d0
        areaeta = 0d0
       ! Covariant components 
        l_r = 0d0
        l_xi = 0d0
        l_eta = 0d0
        area_r = 0d0
        area_xi = 0d0
        area_eta = 0d0
       
        vol = 0d0

      ! lr, areaxi, areaeta, vol
      ! Not defined at rb(0) and rb(np+2).
       lr(1:np+1) = (rb(2:np+2) - rb(0:np)) * belam(1:np+1)
       do j = 0, nxi+1
       do k = 0, neta+1
          lxi(:, j, k) = rb(:)*D(k)*dxi/(delta(j,k)*cxi(j)**2)
          leta(:, j, k) = rb(:)*C(j)*deta/(delta(j,k)*ceta(k)**2)  
       end do
       end do

       ! Covariant components of the length elements
       l_r(1:np+1) = lr(1:np+1)
       do j = 0, nxi+1
       do k = 0, neta+1
          l_xi(:, j, k) = lxi(:,j,k) - X(j)*Y(k)/(C(j)*D(K))* leta(:, j, k)
          l_eta(:, j, k) = leta(:,j,k) - X(j)*Y(k)/(C(j)*D(K))* lxi(:, j, k)
       end do
       end do

      ! Only the covariant surface element are defined 
       do j = 0, nxi+1
       do k = 0, neta+1  
         area_r(:, j, k) = rb(:)**2*deta*dxi/(delta(j,k)**(3/2)  &
          &             * cxi(j)**2*ceta(k)**2)      
         area_xi(1:np+1, j, k) = deta*rb(:)*(rb(2:np+2) - rb(0:np))  &
          &             * belam(1:np+1)/(D(k)*delta(j,k)**(1/2)*ceta(k)**2)
         area_eta(1:np+1, j, k) = dxi*rb(:)*(rb(2:np+2) - rb(0:np)) &
          &             * belam(1:np+1)/(C(j)*delta(j,k)**(1/2)*cxi(j)**2)
         vol(1:np+1, j, k) = deta*dxi*rb(:)**2*(rb(2:np+2) - rb(0:np)) &
          &             *belam(1:np+1)/(delta(j,k)**(3/2)*cxi(j)**2*ceta(k)**2)
       end do
       end do
 
    end subroutine build_grid

    ! --------------------------------------------------------------------------
    !> Main structure calculation routine.
    !> @brief Main interface to calculate the structure. This subroutine:
    !>          1) Defines the radial grid (irregular, finer in the crust)
    !>          2) Calls the TOV solver to calculate the star structure
    !>          3) Calls the EOS to obtain needed variables (composition)
    !>          4) Calculates and stores the supercon/superfluid gaps (uses fits,
    !              there are different models implemented)
    !---------------------------------------------------------------------------
    subroutine build_structure()

      ! Modules ----------------------------------------------------------------
      use constants, only : RHO_TO_N, PI, UNIT_R

      implicit none

      ! Subroutine arguments ---------------------------------------------------
      ! None.

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
      ! Baryon density
      real*8 n_b(np+2)
      ! variable radial grid size
      real*8 dr
      ! supefluid gaps and critical temperatures
      real*8 gapT0, Tc
      ! Auxiliary variables
      integer j, l, index

      ! ------------------------------------------------------------------------

      ! Build the radial grid (for a spherically symmetric NS model).
      ! This has to be manually adapted if resolution or EOS is changed !
      ! To avoid many grid point where resolution is not needed, the radial grid
      ! is coarse in the center (0.25) and fine in the outer crust (0.01).
      ! A smooth transition with a sort of Fermi-Dirac function 1/(1+exp(a*(x-x0)))
      ! enforces that the grid size does not jump abruptly. The numbers are tunes 
      ! for ths particular EOS and star model used (x0=9.41, close to the crust/core 
      ! interface, and a=6, for a typical lmax=99). A quick fix simply rescales the
      ! grid size everywhere when lmax is modified.

      rb(0) = 0.d0
      do l = 1, np+2
        dr = 0.01d0 + (0.25d0 - 0.01d0) / &
          &  (1.d0 + dexp(6.0 * (rb(l-1) - 9.41)))
        !  quick fix to make it run with variable resolutions
        dr = dr*100.d0/(lmax+1)
        rb(l) = rb(l-1) + dr
      end do ! l

      call init_eos_tab()

      call solve_tov_structure(np+2, rb(1:np+2), mass, press, g14)

      ! Call EOS to obtain the composition through the star.
      index = 0
      moment_inertia = 0d0
      do j = 1, np+2

        call geteost(press(j), rho(j), n_b(j), zz(j), aa(j), yn(j), xh(j))

        moment_inertia = moment_inertia + rho(j)*4d0*PI*rb(j)**4*(rb(j)-rb(j-1))*UNIT_R**5
        if (xh(j) == 0d0) then
          ! Charge neutrality in the core.
          yp(j) = 1.d0-yn(j)
          ye(j) = yp(j)
        else
          ! Zero proton fraction in the crust.
          yp(j) = 0d0
          ye(j) = xh(j) * zz(j) / aa(j)

          if (index == 0) then
            ! Label jcore as the last point of the core
            jcore = j - 1
            index = 1
          end if

        end if

      end do ! j

      ! Spin-down torque factor.
      ! PPdot = spindown_prefactor Bp12^2, assuming vacuum orthogonal rotator,
      ! where k = 2 pi^2 R^6/(3 I c^3) (see section 2.2 of Daniele's thesis).
      ! The stellar radius (rns) is rb(np) 
      spindown_prefactor = 2.44d-40 * (rb(np) / 10d0)**6 / (moment_inertia / 1d45) * 1d24

      ! Define the crust and core gaps.
      tccru = 0.d0
      tcn = 0.d0
      tcp = 0.d0

      ! For neutron 1s0 crust.
      if (superfluid_n_crust /= 0) then
        do l = 1, np+2
          if(xh(l) > 0.d0)then
            call gapmodel(superfluid_n_crust, 1, rho(l), yn(l), gapT0, Tc)
            tccru(l)=Tc*1.d-8
          end if
        end do ! l
      end if

      ! For neutron 3p2 core.
      if (superfluid_n_core /= 0) then
        do l = 1, np+2
          if(xh(l) == 0.d0) then
            call gapmodel(superfluid_n_core, 2, rho(l), yn(l), gapT0, Tc)
            tcn(l) = Tc*1.d-8
          end if
        end do ! l
      end if

      ! For proton 1s0 core.
      if (superfluid_p_core /= 0) then
        do l = 1, np+2
          if(xh(l) == 0.d0)then
            call gapmodel(superfluid_p_core, 1, rho(l), yp(l), gapT0, Tc) ! Check index 1 or 3 (was 1 in 2013)
            tcp(l)=Tc*1.d-8
          end if
        end do ! l
      end if

      ! Construction of the Hall factors fh and ia_hall.
      !
      ! Hall factor = c*B0 / 4*pi*e*n_e
      !             = c*1.d12 / (4*pi*1.6d-19*1 statC*(rho*Ye/m_b)
      !             = 3d10*1d12*1.67d-27 / (4*pi*1.6d-19*3d10*rho*Ye)
      !             = 1d8*1.673 / (1.602*4*pi)
      ! Multiplied by conversion factor k=1.d-10*3.16d13
      ! cm^2*(B0/1e12 G)/Myr = k Km^2*(B0/1e12 G)/Myr
      !
      ! Below the conversion factors are obtained with nuclear physics
      ! constants (fine structure, hc, and so on).
      !
      ! The final factor is
      ! fh [Km^2/Myr*(1d12 G)] ~ 2.6d10/(Ye*rho[g/cm3])
      fh=0d0
      ia_amb=0
      ia_hall=0

      do l = 1, np+2

        fh(l) = 1.564d-5 / ((rho(l) / RHO_TO_N) * ye(l))

        if(xh(l) /= 0d0) then
          ! Hall effect in the crust.
          ia_hall(l) = 1
        else
          ! Ambipolar diffusion in the core.
          ia_amb(l) = 1
        end if

      end do ! l

      ! Manually set the first two points (to avoid spurious derivatives).
      fh(0) = fh(1)
      ia_hall(0) = ia_hall(1)
      ia_amb(0) = ia_amb(1)

      ! Apply flag to turn on/off Hall effect.
      if (enable_hall_effect .eqv. .false.) then
        ia_hall = 0
      end if

      ! Output structure grid and gaps to files, write screen output information.
      call output_structure()

    end subroutine build_structure

    !!--------------------------------------------------------------------------
    !> @brief Subroutine metric_tensor
    !!
    !! In this subroutine we are defining our non-orthogonal metric tensor
    !! in the versor basis
    !! Code owners:
    !!    Clara Dehman
    !!--------------------------------------------------------------------------
    subroutine metric_tensor()
    
    implicit none
  
  ! Local variables ------------------------------------------------------------
    integer :: i,j,k

  ! initialize to zero
    g = 0d0  
     
  ! diagonal terms
     g(:, :, :, 1, 1) = 1d0
     g(:, :, :, 2, 2) = 1d0
     g(:, :, :, 3, 3) = 1d0

  ! off-diagonal terms
    do j = 1, nxi
    do k = 1, neta
     g(:, j, k, 2, 3) = -X(j)*Y(k)/(C(j)*D(k)) 
    end do
    end do 
    g(:, :, :, 3, 2) =  g(:, :, :, 2, 3) 

    g(:, :, :, 1, 2) = 0d0
    g(:, :, :, 1, 3) = 0d0
    g(:, :, :, 2, 1) = 0d0
    g(:, :, :, 3, 1) = 0d0

   end subroutine metric_tensor

    !!--------------------------------------------------------------------------
    !> @brief Subroutine dot_prod
    !!
    !! In this subroutine we express the dot product in term of the non-orthogonal
    !! metric tensor
    !! We are working in versor basis
    !! Code owners:
    !!    Clara Dehman
    !!--------------------------------------------------------------------------
    subroutine dot_prod(mr,nr,meta,neta,mxi,nxi,dprod)

    implicit none
    integer :: i, j
    real*8, dimension (0:np+2, 0:nxi+1, 0:neta+1), intent(in) :: mr, meta, mxi
    real*8, dimension (0:np+2, 0:nxi+1, 0:neta+1), intent(in) :: nr, neta, nxi
    real*8, dimension (0:np+2, 0:nxi+1, 0:neta+1), intent(out) :: dprod 

  ! initialize to zero
    dprod = 0d0     
  ! Writing the dot product in term of the metric g
     do i = 0, np+2 
     do j = 0, nxi+1
      dprod(i, j, :) =  g(i, j, :, 1, 1)*mr(i, j, :)*nr(i, j, :) & 
          & + g(i, j, :, 2, 2)*mxi(i, j, :)*nxi(i, j, :) &
          & + g(i, j, :, 3, 3)*meta(i, j, :)*neta(i, j, :)  &
          & + g(i, j, :, 2, 3)*(mxi(i, j, :)*neta(i, j, :) + meta(i, j, :)*nxi(i, j, :))
     end do 
     end do 
    end subroutine dot_prod

    !!--------------------------------------------------------------------------
    !> @brief Subroutine crossprod_cont
    !!
    !! In this subroutine we define the contravariant components of the cross product
    !! in term of the contravariant components of the vectors
    !! using the versor basis for our non-orthogonal metric
    !! Code owners:
    !!    Clara Dehman
    !!--------------------------------------------------------------------------
    !! Note: Xprodr refers to the contravariant component of the cross product 
    !! in r physical direction
    !!--------------------------------------------------------------------------
    subroutine crossprod_cont(mr,nr,meta,neta,mxi,nxi,Xprodxi, Xprodeta, Xprodr)

    implicit none
    integer :: i, j
    real*8, dimension (0:np+2, 0:nxi+1, 0:neta+1), intent(in) :: mr, meta, mxi
    real*8, dimension (0:np+2, 0:nxi+1, 0:neta+1), intent(in) :: nr, neta, nxi
    real*8, dimension (0:np+2, 0:nxi+1, 0:neta+1), intent(out) :: Xprodxi, Xprodeta, Xprodr
    
  ! initialize to zero
    Xprodxi = 0.d0
    Xprodeta = 0.d0
    Xprodr = 0.d0

    do i = 0, np+2
    do j = 0, nxi+1
      Xprodr(i, j, :) = delta(j, :)^(1/2)/(C(j)*D(:)) &
          &   * (mxi(i, j, :)*neta(i, j, :) - meta(i, j, :)*nxi(i, j, :))
      Xprodxi(i, j, :) = (C(j)*D(:)*(meta(i, j, :)*nr(i, j, :) - mr(i, j, :)*neta(i, j, :)) &
          &   + X(j)*Y(:)*(mr(i, j, :)*nxi(i, j, :) - mxi(i, j, :)*nr(i, j, :)))/(delta(j, :)^(1/2))      
      Xprodeta(i, j, :) = (X(j)*Y(:)*(meta(i, j, :)*nr(i, j, :) - mr(i, j, :)*neta(i, j, :)) &
          &   + C(j)*D(:)*(mr(i, j, :)*nxi(i, j, :) - mxi(i, j, :)*nr(i, j, :)))/(delta(j, :)^(1/2))
    end do 
    end do 

    end subroutine crossprod_cont


    !!--------------------------------------------------------------------------
    !> @brief Subroutine crossprod_cov
    !!
    !! In this subroutine we define the covariant components of the cross product
    !! in term of the contravariant components of the vectors
    !! using the versor basis for our non-orthogonal metric
    !! Code owners:
    !!    Clara Dehman
    !!--------------------------------------------------------------------------
    !! Note: Xprod_r refers to the covariant component of the cross product 
    !! in r physical direction
    !!--------------------------------------------------------------------------
    subroutine crossprod_cov(mr,nr,meta,neta,mxi,nxi,Xprod_xi, Xprod_eta, Xprod_r)

    implicit none
    integer :: i, j
    real*8, dimension (0:np+2, 0:nxi+1, 0:neta+1), intent(in) :: mr, meta, mxi
    real*8, dimension (0:np+2, 0:nxi+1, 0:neta+1), intent(in) :: nr, neta, nxi
    real*8, dimension (0:np+2, 0:nxi+1, 0:neta+1), intent(out) :: Xprod_xi, Xprod_eta, Xprod_r
    
  ! initialize to zero
    Xprod_xi = 0.d0
    Xprod_eta = 0.d0
    Xprod_r = 0.d0

    do i = 0, np+2
    do j = 0, nxi+1
      Xprod_r(i, j, :) = delta(j, :)^(1/2)/(C(j)*D(:)) &
          &   * (mxi(i, j, :)*neta(i, j, :) - meta(i, j, :)*nxi(i, j, :))
      Xprod_xi(i, j, :) = delta(j, :)^(1/2)/(C(j)*D(:)) &
          &   * (meta(i, j, :)*nr(i, j, :) - mr(i, j, :)*neta(i, j, :))
      Xprod_eta(i, j, :) = delta(j, :)^(1/2)/(C(j)*D(:)) &
          &   * (mr(i, j, :)*nxi(i, j, :) - mxi(i, j, :)*nr(i, j, :))
    end do 
    end do 

    end subroutine crossprod_cov

    ! --------------------------------------------------------------------------
    !> The subroutine xighost_position return the value of j1 and j2 
    !> such as xi(j1) and xi(j2) are surrounding xi_ghost value.
    !> xi_ghost are the ghost cells along the xi direction and for a given
    !> eta and r value 
    !> Clara Dehman 
    !---------------------------------------------------------------------------
    subroutine xighost_position(j1,j2)
      
    integer j 
    integer, intent(out) :: j1,j2

    do j = 0, nxi
    if(xi_ghost > xi(j) .and. xi_ghost < xi(j+1))
    j = j1
    j+1 = j2
    end if
    end do
    end subroutine xighost_position

    ! --------------------------------------------------------------------------
    !> The subroutine etaghost_position return the value of k1 and k2 
    !> such as eta(k1) and eta(k2) are surrounding eta_ghost value.
    !> eta_ghost are the ghost cells along the eta direction and for a given 
    !> xi and r value 
    !> Clara Dehman 
    !---------------------------------------------------------------------------
    subroutine etaghost_position(k1,k2)

    integer k 
    integer, intent(out) :: k1,k2

    do k = 0, neta
    if(eta_ghost > eta(k) .and. eta_ghost < eta(k+1))
    k = k1
    k+1 = k2
    end if
    end do
    end subroutine etaghost_position


    ! --------------------------------------------------------------------------
    !> Output struture information to files and screen.
    !> @brief Outputs structure grid and gaps information to files and also
    !>        useful information (surface density, core radius...) to screen.
    !---------------------------------------------------------------------------
    subroutine output_structure()

      ! Modules ----------------------------------------------------------------
      use utils, only: get_free_unit
      use constants, only: T8_TO_MEV, GRAV_CONSTANT, MSUN, CLIGHT, UNIT_R

      implicit none

      ! Subroutine arguments ---------------------------------------------------
      ! None.

      ! Local constants --------------------------------------------------------
      character(len=15), parameter :: GRID_FORMAT = "(11(1pe14.7))"
      character(len=14), parameter :: GAPS_FORMAT = "(11(1x,e15.7))"

      ! Local variables --------------------------------------------------------
      integer :: j, l
      integer :: unit
      ! Dimensionless parameter GM/c^2R
      real*8 compactness

      ! ------------------------------------------------------------------------

      unit = get_free_unit()

      open(unit, file="out/grid.d")
      do j = 1, np+2
        write(unit, GRID_FORMAT) &
          & rb(j), rho(j), ye(j), aa(j), zz(j), &
          & xh(j), yn(j), benu(j), belam(j), mass(j), press(j)
      end do ! j
      close(unit)

      open(unit, file="out/gaps.d")
      do l = 1, np
        write(unit, GAPS_FORMAT) &
          & rho(l), T8_TO_MEV*tccru(l), T8_TO_MEV*tcn(l), T8_TO_MEV*tcp(l)
      end do
      close(unit)
      
      compactness = GRAV_CONSTANT*MSUN/(CLIGHT**2*UNIT_R)*mass(np)/rb(np)

      ! Screen output.
      write(*, "(a,1pe11.3)") "<info>[STRUCTURE] Surface density [g/cm^3]:", rho(np)
      write(*, "(a,f8.3)") "<info>[STRUCTURE] Core radius [km]:", rb(jcore)
      write(*, "(a,f8.3)") "<info>[STRUCTURE] Star radius [km]:", rb(np)
      write(*, "(a,f8.3)") "<info>[STRUCTURE] Mass [Msun]:", mass(np)
      write(*, "(a,f8.3)") "<info>[STRUCTURE] Compactness GM/Rc^2 :", compactness
      write(*, "(a,f8.3)") "<info>[STRUCTURE] Surface gravity [10^14 cm/s^2] :", g14
      write(*, "(a,f8.3)") "<info>[STRUCTURE] Surface redshift (1+z) :", 1.d0/benu(np)
      write(*, "(a,1pe11.3)") "<info>[STRUCTURE] Moment of Inertia [g*cm^2]:", moment_inertia
      write(*, "(a)") "<info>[STRUCTURE] Various quantities initialized."
      write(*, *)

    end subroutine output_structure

end module grid




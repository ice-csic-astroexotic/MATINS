!-------------------------------------------------------------------------------
! Module: Grid
!
!> @brief This module is responsible for handling the grid and structure.
!
!> @author
! Clara Dehman 
! Daniele Viganò
! Stefano Ascenzi
!-------------------------------------------------------------------------------
module grid

  use OMP_LIB 
  use input_params
  use constants, only : PI

  implicit none

  integer, save :: nr, nang     ! Number of points in the full grid (magnetic field)
  integer, save :: nrt, nangt  ! Number of points in the thermal grid
  integer, save :: lmax
  real*8, parameter :: rmax = 11.d0 ! maximum radius
  real*8, parameter :: rmin = 10.d0  ! minimum radius
  real*8, dimension(:,:,:,:,:), allocatable, save :: g        ! Metric
  real*8, dimension(:,:,:,:), allocatable, save :: jac_eq     ! Jacobian equatorial patches (2x2 angular directions)
  real*8, dimension(:,:,:,:), allocatable, save :: jacinv_eq  ! Inverse Jacobian equatorial patches
  real*8, dimension(:,:,:,:), allocatable, save :: jac_n      ! Jacobian north patch
  real*8, dimension(:,:,:,:), allocatable, save :: jacinv_n   ! Inverse Jacobian north patch
  real*8, dimension(:,:,:,:), allocatable, save :: jac_s      ! Jacobian south patch
  real*8, dimension(:,:,:,:), allocatable, save :: jacinv_s   ! Inverse Jacobian south patch
  real*8, dimension(:,:,:,:), allocatable, save :: jac_pp     ! Jacobian: patch to patch - used for the ghost pointss
  real*8, dimension(:,:,:,:), allocatable, save :: jac_avg    ! Jacobian: patch to patch - used for the average of the field components
  real*8, dimension(:,:,:,:), allocatable, save :: xc, yc, zc ! Cartesian coordinates in each patch
  real*8, dimension(:,:,:), allocatable, save :: theta, phi   ! Spherical coordinates in each patch
  !-----------------------------------------------------------------------------
  real*8, dimension(:), allocatable, save :: r   ! radial direction
  real*8, dimension(:), allocatable, save :: eta ! eta direction
  real*8, dimension(:), allocatable, save :: xi  ! xi direction
  real*8, dimension(:), allocatable, save :: X  ! tan(xi)
  real*8, dimension(:), allocatable, save :: Y  ! tan(eta)
  real*8, dimension(:, :), allocatable, save :: delta ! delta = 1+X^2+Y^2
  real*8, dimension(:), allocatable, save :: C ! C=(1+X^2)^{1/2}
  real*8, dimension(:), allocatable, save :: D ! D=(1+Y^2)^{1/2}
  real*8, dimension(:), allocatable, save :: pm ! interpolated xi or eta
  real*8, dimension(:), allocatable, save :: pmt ! interpolated xi or eta for the thermal grid
  real*8, dimension(:), allocatable, save :: elambda ! relativistic factor of the radial component
  real*8, dimension(:), allocatable, save :: enu, enu_core ! redshift factor
  real*8, dimension(:), allocatable, save :: rtot, elambda_tot, enu_tot
  ! Length, area, and volume elements. 
  real*8, dimension(:, :, :), allocatable, save :: vol ! Volume of each cell
  real*8, dimension(:), allocatable, save :: vol_shell   ! Volume of shells in the core
  real*8, dimension(:), allocatable, save :: lr ! radial length element of the cell (covariant and contravariant are the same)
  ! contravariant components
  ! DV: Are they used? If not, let's not save them
  real*8, dimension(:, :, :), allocatable, save :: lxi ! length element of the cell in xi physical direction.
  real*8, dimension(:, :, :), allocatable, save :: leta ! length element of the cell in eta physical direction.
!  real*8, dimension(:, :, :), allocatable, save :: arear ! Area of the radial interface
!  real*8, dimension(:, :, :), allocatable, save :: areaxi 
!  real*8, dimension(:, :, :), allocatable, save :: areaeta 
  ! Covariant components 
  real*8, dimension(:, :, :), allocatable, save :: l_xi
  real*8, dimension(:, :, :), allocatable, save :: l_eta
  real*8, dimension(:, :, :), allocatable, save :: area_r 
  real*8, dimension(:, :, :), allocatable, save :: area_xi 
  real*8, dimension(:, :, :), allocatable, save :: area_eta

  ! Values of theta at the longitudes 0, 90, 180, 270 used in the 1D output
  real*8, dimension(:), allocatable, save :: theta_meridian, theta_meridian_2PI
  ! Values of phi at the equator used in the 1D output
  real*8, dimension(:), allocatable, save :: phi_equator
  ! Factors related to the grid edges
  real*8, dimension(:,:), allocatable, save :: wint   ! Weight used for integrals accounting for edges
  real*8, dimension(:,:,:), allocatable, save :: pos_pm      ! Set of the projected position at each side of each edge, parallel coordinate
  real*8, dimension(:,:), allocatable, save :: pos_qm        ! Projected position at each side of each edge, pseudo-perpendicular coordinate
  real*8, dimension(:), allocatable, save :: edge_w, edge_wt    ! Linear interpolation weights at the edges for magnetic and thermal grids


  ! Star's structure quantities
  ! Radial profiles of:
  !   energy density (rho)
  !   baryon density (nb)
  !   nuclei fraction (xh)
  !   electron fraction (ye)
  !   free neutrons fraction (yn)
  !   free protons fraction (yp)
  !   atomic mass number (aa)
  !   atomic number (zz)
  !   free neutrons density (nn)
  !   free protons density (npr)
  !   electron density (ne)
  !   Fermi momentum for e, p, n, mu (kFe, kFp, kFn, kFmu)
  !   effective masses of e, n, p (in units of their rest masses)
  real*8, dimension (:), allocatable, save :: rho, nb, xh, ye, yn, yp, ymu, aa, zz
  real*8, dimension (:), allocatable, save :: nn, npr, ne, nmu, kFn, kFp, kFe, kFmu
  real*8, dimension (:), allocatable, save :: effme, effmn, effmp
  real*8, dimension (:), allocatable, save :: Zimp
  real*8, save :: g14   ! Surface gravity [10^14 gr*cm/s^2].
  real*8, save :: moment_inertia   ! Moment of Inertia [g cm^2].
  real*8, dimension(:), allocatable, save :: mass, press   ! Enclosed mass and pressure.

  ! Cooling quantities
  real*8, save :: T_core      ! Redshifted temperature of the core
  real*8, save :: cv_core_tot, qnu_core_tot, qnu_core_tot_der, cv_core_tot_der
  real*8, dimension(:, :, :, :), allocatable, save :: temp, tem0 ! Temperature (redshifted and physical)
  real*8, dimension(:,:,:), allocatable, save :: temp_surf !Surface Temperature
  real*8, dimension(:,:,:), allocatable, save :: bb_flux, sfluxb, cfluxb !BB flux and flux derivative
  real*8, dimension(:, :, :, :), allocatable, save :: q_neutrino, q_neutrino_der    ! Neutrino emissivity and its T-derivative (used in the implicit part of neutrino emissivity) in units [10^40 erg/(km^3*s*10**8 K)]
  real*8, dimension(:, :, :, :), allocatable, save :: omegatau_arr, cv, kappa_perp_arr
  real*8, dimension(:), allocatable :: cv_core, q_neutrino_core
  real*8, dimension(:,:,:,:), allocatable :: flux_r_out, flux_xi_xip, flux_eta_etap
  ! Superfluid and superconducting gaps: critical temperatures in units of 10**8 K.
  real*8, dimension(:), allocatable, save :: tcn, tcp, tccru
  real*8, dimension(:), allocatable :: gapn_crust, gapn_core, gapp_core

  ! Magnetic quantities
  real*8, dimension (:), allocatable, save :: fh
  real*8, dimension (:,:,:,:), allocatable, save :: etab
  real*8, dimension(:, :, :, :), allocatable, save :: br,bxi,beta ! Contravariant magnetic field components [10^12 G] 
  real*8, dimension(:, :, :, :), allocatable, save :: er,exi,eeta ! Contravariant electric field components    
  real*8, dimension(:, :, :, :), allocatable, save :: jr,jxi,jeta ! Contravariant current components [10^12 G]  
  real*8, dimension(:, :, :,:), allocatable, save :: bm ! Modulus of B 
  real*8, dimension(:, :, :,:), allocatable, save :: b2 ! Square of B
  ! Associated legendre polynomial, spherical harmonic and related quantities
  real*8, dimension(:,:,:,:,:), allocatable, save :: y_lm, dyth_lm, dyphi_lm!, dy2phi_lm
  real*8, dimension(:), allocatable, save :: frel
  real*8, dimension(:,:), allocatable, save :: blm
  real*8, dimension(:,:), allocatable, save :: espec_vol, espec_pol, espec_tor ! Energy spectrum 
  real*8, dimension(:,:,:,:), allocatable, save ::  phi_scalar, psi_scalar ! Phi and Psi scalar functions
  real*8, save :: bpdip

  ! Joule dissipation, they include the e^2nu factor, in units [10^40 erg/(km^3*s)]
  real*8, dimension(:, :, :, :), allocatable, save :: q_joule  ! Instantaneous dissipation
  real*8, dimension(:, :, :, :), allocatable, save :: j2       ! Square of the currents J^2 [10^24 (G/km)^2]

  real*8, save :: en_joule_star_tot, poynting_star_tot, poynting_star_tot_surface, poynting_star_tot_interior 

  ! ievol is the first value to calculate curl operator and to evolve 
  integer, save :: ievol, ncore
!  real*8, save :: alpha

  contains
  
    subroutine allocate_grid()

      nrt = radial_dimension
      nangt = angular_dimension
      ncore = radial_dimension_core

      write(6,'(a,i3,a3,i2,a2)'),"Allocating grid with resolution (thermal grid): 6 x ",nrt," x ",nangt,"^2"

      ! Dimensions of the grid
      nr = 2*nrt     ! i = 2*it - 1 for the thermal cell centers
      nang = 2*nangt + 1 ! j = 2*jt, k = 2*kt for the thermal cell centers
      lmax = 30       ! Maximum degree of spherical harmonics

      allocate(r(0:nr+1))
      ! Cubed sphere (CS) coordinates and related quantities (Ronchi 1996)
      allocate(xi(0:nang+1))
      allocate(eta(0:nang+1))
      allocate(X(0:nang+1))
      allocate(Y(0:nang+1))
      allocate(delta(0:nang+1, 0:nang+1))  
      allocate(C(0:nang+1))
      allocate(D(0:nang+1))
      allocate(pm(0:nang+1))
      allocate(pmt(0:nang+1))
      ! Metric and Jacobians spherical to CS
      allocate(g(0:nr+1,0:nang+1,0:nang+1,3,3))
      allocate(jac_eq(0:nang+1,0:nang+1,2,2))
      allocate(jacinv_eq(0:nang+1,0:nang+1,2,2))
      allocate(jac_n(0:nang+1,0:nang+1,2,2))
      allocate(jacinv_n(0:nang+1,0:nang+1,2,2))
      allocate(jac_s(0:nang+1,0:nang+1,2,2))
      allocate(jacinv_s(0:nang+1,0:nang+1,2,2))
      allocate(jac_pp(1:18,0:nang+1,2,2))
      allocate(jac_avg(1:18,0:nang+1,2,2))
      ! Relativistic factor
      allocate(elambda(0:nr+1))
      allocate(enu(0:nr+1))
      ! Geometrical elements
      allocate(vol(0:nr+1, 0:nang+1, 0:nang+1))
      allocate(lr(0:nr+1))
      allocate(lxi(0:nr+1, 0:nang+1, 0:nang+1))
      allocate(leta(0:nr+1, 0:nang+1, 0:nang+1))
      allocate(l_xi(0:nr+1, 0:nang+1, 0:nang+1))
      allocate(l_eta(0:nr+1, 0:nang+1, 0:nang+1))
      allocate(area_xi(0:nr+1, 0:nang+1, 0:nang+1))
      allocate(area_r(0:nr+1, 0:nang+1, 0:nang+1))
      allocate(area_eta(0:nr+1, 0:nang+1, 0:nang+1))
      ! Cartesian and spherical coordinates
      allocate(xc(0:nr+1,0:nang+1,0:nang+1,6))
      allocate(yc(0:nr+1,0:nang+1,0:nang+1,6))
      allocate(zc(0:nr+1,0:nang+1,0:nang+1,6)) 
      allocate(theta(0:nang+1,0:nang+1,6))
      allocate(phi(0:nang+1,0:nang+1,6))
      allocate(theta_meridian(1:2*nang-1))
      allocate(phi_equator(1:4*nang-3))
      allocate(theta_meridian_2PI(1:4*nang-3))
      ! Grid weights
      allocate(edge_w(0:nang+1)) ! Interpolation weights.
      allocate(edge_wt(0:nang+1)) ! Interpolation weights for thermal grid
      allocate(wint(0:nang+1, 0:nang+1)) ! Integration weights to avoid double counting
      ! Structure
      allocate(rtot(0:ncore+nr+1))
      allocate(enu_tot(0:ncore+nr+1))
      allocate(elambda_tot(0:ncore+nr+1))
      allocate(vol_shell(1:ncore))
      allocate(enu_core(1:ncore))
      allocate(rho(1:ncore+nr))
      allocate(nb(1:ncore+nr))
      allocate(xh(1:ncore+nr))
      allocate(ye(1:ncore+nr))
      allocate(yn(1:ncore+nr))
      allocate(yp(1:ncore+nr))
      allocate(ymu(1:ncore+nr))
      allocate(aa(1:ncore+nr))
      allocate(zz(1:ncore+nr))
      allocate(nn(1:ncore+nr))
      allocate(npr(1:ncore+nr))
      allocate(ne(1:ncore+nr))
      allocate(nmu(1:ncore+nr))
      allocate(kFn(1:ncore+nr))
      allocate(kFp(1:ncore+nr))
      allocate(kFe(1:ncore+nr))
      allocate(kFmu(1:ncore+nr))
      allocate(mass(1:ncore+nr))
      allocate(press(1:ncore+nr))
      allocate(effme(1:ncore+nr))
      allocate(effmn(1:ncore+nr))
      allocate(effmp(1:ncore+nr))
      ! Superfluidity (crust and core)
      allocate(tcn(1:ncore+nr))
      allocate(tcp(1:ncore+nr))
      allocate(tccru(1:ncore+nr))
      allocate(gapn_crust(1:ncore+nr))
      allocate(gapn_core(1:ncore+nr))
      allocate(gapp_core(1:ncore+nr))
      ! Temperature in the crust
      allocate(temp(1:nrt, 0:nangt+1, 0:nangt+1, 1:6))
      allocate(tem0(1:nrt, 0:nangt+1, 0:nangt+1, 1:6))
      ! Thermal Fluxes (evalueted in the outermost or rightermost cell interface)
      allocate(flux_r_out(1:nrt, 1:6, 0:nangt, 0:nangt))
      allocate(flux_xi_xip(1:nrt, 1:6, 0:nangt, 0:nangt))
      allocate(flux_eta_etap(1:nrt, 1:6, 0:nangt, 0:nangt))
      ! Envelope 
      allocate(temp_surf(0:nangt+1, 0:nangt+1, 1:6))
      allocate(bb_flux(1:nangt, 1:nangt, 1:6))
      allocate(sfluxb(1:nangt, 1:nangt, 1:6))
      allocate(cfluxb(1:nangt, 1:nangt, 1:6))      
      ! Microphysical quantities for the crust
      allocate(Zimp(0:nr+1))
      allocate(etab(0:nr+1,0:nang+1,0:nang+1,1:6))
      allocate(fh(0:nr+1))
      allocate(omegatau_arr(1:nrt, 0:nangt+1, 0:nangt+1, 1:6))
      allocate(kappa_perp_arr(1:nrt, 0:nangt+1, 0:nangt+1, 1:6))
      allocate(cv_core(1:ncore))
      allocate(cv(1:nrt, 0:nangt+1, 0:nangt+1, 1:6))
      allocate(q_neutrino_core(1:ncore))
      allocate(q_neutrino(1:nrt, 0:nangt+1, 0:nangt+1, 1:6))
      allocate(q_neutrino_der(1:nrt, 0:nangt+1, 0:nangt+1, 1:6))
      ! Contravariant components of the magnetic field
      allocate(br(0:nr+1, 0:nang+1, 0:nang+1, 1:6))
      allocate(bxi(0:nr+1, 0:nang+1, 0:nang+1, 1:6))
      allocate(beta(0:nr+1, 0:nang+1, 0:nang+1, 1:6))
      allocate(bm(0:nr+1, 0:nang+1, 0:nang+1, 1:6))
      allocate(b2(0:nr+1, 0:nang+1, 0:nang+1, 1:6))
      ! Contravariant components of the electric field 
      allocate(er(0:nr+1, 0:nang+1, 0:nang+1, 1:6))
      allocate(exi(0:nr+1, 0:nang+1, 0:nang+1, 1:6))
      allocate(eeta(0:nr+1, 0:nang+1, 0:nang+1, 1:6))
      ! Contravariant components of the electric current
      allocate(jr(0:nr+1, 0:nang+1, 0:nang+1, 1:6))
      allocate(jxi(0:nr+1, 0:nang+1, 0:nang+1, 1:6))
      allocate(jeta(0:nr+1, 0:nang+1, 0:nang+1, 1:6))
      allocate(j2(0:nr+1, 0:nang+1, 0:nang+1, 1:6))        
      ! Associated Legendre polynomial
      allocate(y_lm(0:nang+1, 0:nang+1, 1:6, 0:lmax, -lmax:lmax))
      allocate(dyth_lm(0:nang+1, 0:nang+1, 1:6, 0:lmax, -lmax:lmax))
      allocate(dyphi_lm(0:nang+1, 0:nang+1, 1:6, 0:lmax, -lmax:lmax))
      allocate(frel(1:lmax))
      ! Multipole weights
      allocate(blm(0:lmax, -lmax:lmax))
      allocate(espec_vol(0:lmax, -lmax:lmax))
      allocate(espec_pol(0:lmax, -lmax:lmax))
      allocate(espec_tor(0:lmax, -lmax:lmax))
      ! Phi and Psi scalar functions
      allocate(phi_scalar(0:nr+1,0:nang+1,0:nang+1,1:6))
      allocate(psi_scalar(0:nr+1,0:nang+1,0:nang+1,1:6))
    ! Joule heating calculation
      allocate(q_joule(nrt, nangt, nangt, 6))
    end subroutine allocate_grid

    !---------------------------------------------------------------------------
    !> Grid builder.
    !> Infinitesimal geometrical elements.
    !> @brief Construction the whole grid by initializing the arrays in the 
    !>  cubed-sphere formalism.
    !> We are defining the geometrical elements in one patch one since in the 
    !> other patches we will have the same expressions for the infinitesimal 
    !> geometrical elements
    !> 
    !> Code owners:
    !> Clara Dehman
    !> Daniele Vigano
    !---------------------------------------------------------------------------
    subroutine build_grid()

      implicit none

      ! Subroutine arguments ---------------------------------------------------
      ! None.

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
      integer :: i, j, k!, p
      real*8 :: dang ! Angular step in xi and eta direction, constant
      real*8 :: dr ! Step in radial direction
      real*8 :: da, dat

      if (EoS == "simple") then

        ncore = 0
        if (rmin == 0.) then
          dr = rmax / (dble(nr))
          do i = 0, nr+1
            r(i) = dble(i) * dr
          end do
        else if (rmin < 0.) then
          print*, "<grid.f90> ERROR: rmin negative"
          stop
        else
          dr = (rmax-rmin)/(dble(nr-1))
          do i = 0, nr+1
            r(i) = rmin + dble(i-1) * dr
          end do
        end if
  
        if (r(0) < 0.) then
          print*, "<grid.f90> ERROR: r(0) negative. Check values of rmin and nrt"
          stop
        endif

        elambda = 1d0
        enu = 1d0
        g14 = 1d0
       
      else

        call build_structure()

        if (use_relativistic_grid .eqv. .false.) then
          elambda = 1d0
          enu = 1d0
          enu_core = 1d0
        endif

      endif

      ! We start to evolve the second crustal cells radially
      ievol = 2

      ! The angular variables xi and eta are spanning the range [-Pi/4 : Pi/4] 
      dang = PI / (2d0*dble(nang-1))
      do j = 0, nang+1
       xi(j) = - 0.25d0*PI + (dble(j) - 1) * dang 
      end do
      eta(:) = xi(:)

      da = xi(1) - xi(0)
      pm = datan( dtan(xi)/dtan(0.25d0*PI + da) )

      dat = xi(2) - xi(0)
      pmt = datan( dtan(xi)/dtan(0.25d0*PI + dat))
      !pmt = datan( dtan(xi)/dtan(0.25d0*PI + da))

      X(:) = dtan(xi(:))
      Y(:) = dtan(eta(:))
      C(:) = (1 + X(:)**2)**(0.5d0)
      D(:) = (1 + Y(:)**2)**(0.5d0)
      do j = 0, nang+1
        delta(j, :) = 1 + X(j)**2 + Y(:)**2
      end do
     ! Contravariant components
      lr = 0d0
      lxi = 0d0
      leta = 0d0
     ! Covariant components 
      l_xi = 0d0
      l_eta = 0d0
      area_r = 0d0
      area_xi = 0d0
      area_eta = 0d0
     
      vol = 0d0
      lr(1:nr) = (r(2:nr+1) - r(0:nr-1)) * elambda(1:nr)

      ! Contravariant components of the length elements (appearing in Ronchi)
      do j = 0, nang+1
        do k = 0, nang+1
          lxi(:, j, k) = 2*r(:)*C(j)**2*D(k)*dang/delta(j,k)
          leta(:, j, k) = 2*r(:)*C(j)*D(k)**2*dang/delta(j,k)  
        end do
      end do
       ! Covariant components of the length elements
      do j = 0, nang+1
        do k = 0, nang+1
          l_xi(:, j, k) = 2*r(:)*D(k)/delta(j,k)*(C(j)**2*dang - X(j)*Y(k)*dang)
          l_eta(:, j, k) = 2*r(:)*C(j)/delta(j,k)*(D(k)**2*dang - X(j)*Y(k)*dang) 
        end do
      end do
      
      ! Only the covariant surface element are defined 
      do j = 0, nang+1
        do k = 0, nang+1  
          area_r(:, j, k) = 4*r(:)**2*C(j)**2*D(k)**2*dang*dang/(delta(j,k)**1.5d0)      
          area_xi(1:nr, j, k) = 2*dang*r(1:nr)*(r(2:nr+1) - r(0:nr-1))  &
          &             * elambda(1:nr)*D(k)/(delta(j,k)**0.5d0)
          area_eta(1:nr, j, k) = 2*dang*r(1:nr)*(r(2:nr+1) - r(0:nr-1)) &
          &             * elambda(1:nr)*C(j)/(delta(j,k)**0.5d0)
          vol(1:nr, j, k) = 4*dang*dang*r(1:nr)**2*(r(2:nr+1) - r(0:nr-1)) &
           &             *elambda(1:nr)*C(j)**2*D(k)**2/(delta(j,k)**1.5d0)
        end do
      end do

      ! Define metric (the same for all patches)
      call metric_tensor
      ! Define Jacobian matrixes (direct and inverse) to pass from spherical to CS and viceversa
      call jacobians
      ! Define patch-to-patch Jacobians at all edges -- used for the ghost points
      call jacobian_edges
      ! Define patch-to-patch Jacobians at all interfaces -- used for the average of the field components 
      call jacobian_average

      ! Call CS to Cartesian coordinate transformation,
      ! to assign the (xc,yc,zc) values at each CS point (r,xi,eta,patch)
      call cs_to_cartesian
      ! Call CS to spherical coordinate transformation,
      ! to assign the (theta,phi) values at each CS point (xi,eta,patch)    
      call cs_to_spherical

      call get_wint

      call spherical_harmonics
      
      call edges

    end subroutine build_grid


    !!--------------------------------------------------------------------------
    !> @brief Subroutine metric_tensor
    !!
    !! In this subroutine we are defining our non-orthogonal metric tensor
    !! in the versor basis 
    !! This metric is the same for all the patches
    !! Code owners:
    !!    Clara Dehman
    !!--------------------------------------------------------------------------
    subroutine metric_tensor
      
      implicit none
    
      ! Local variables ------------------------------------------------------------
      integer :: j, k
         ! initialize to zero
      g = 0d0  
     
      ! diagonal terms
      g(:, :, :, 1, 1) = 1d0
      g(:, :, :, 2, 2) = 1d0
      g(:, :, :, 3, 3) = 1d0
       ! off-diagonal terms
      do j = 1, nang
        do k = 1, nang
          g(:, j, k, 2, 3) = - X(j)*Y(k)/(C(j)*D(k)) 
        end do
      end do 
      g(:, :, :, 3, 2) = g(:, :, :, 2, 3) 
  
    end subroutine metric_tensor

!!--------------------------------------------------------------------------
    !> @brief subroutine jacobians
    !!
    !! In this subroutine, we are defining the jacobian matrices and their inverse
    !! introduced in Ronchi's formalism of 1996.  
    !! 
    !! The Jacobian matrices are used to go from the spherical to cubed sphere coordinates.
    !! The inverse of the Jacobian matrices are used to go from cubed sphere to spherical coordinates.
    !! (Axi,Aeta) = jac (Atheta,Aphi);  (Atheta,Aphi) = jacinv (Axi,Aeta)
    !!   Matrix indexing:  | jac(1,1)   jac(1,2) |
    !!                     | jac(2,1)   jac(2,2) |
    !!
    !! jac_eq and jacinv_eq refer to all the four defined equator patches (patch I to patch IV). 
    !! jac_n and jacinv_n refer to the north patch (patch V).
    !! jac_s and jacinv_s refer to the south patch (patch VI). 
    !!
    !> @author
    ! Clara Dehman 
    ! Daniele Viganò
!!--------------------------------------------------------------------------
    subroutine jacobians

      implicit none
  
      integer :: k

      ! equator jacobian, common terms
      jac_eq(:, :, 1, 1) = 0
      jac_eq(:, :, 2, 1) = -1

      jacinv_eq(:, :, 2, 2) = 0
      jacinv_eq(:, :, 1, 2) = -1

      do k = 0, nang+1
        ! equator jacobian terms
        jac_eq(:, k, 2, 2) = X(:)*Y(k)/dsqrt(delta(:,k))
        jac_eq(:, k, 1, 2) = C(:)*D(k)/dsqrt(delta(:,k))

        ! north jacobian matrix
        ! diagonal terms
        jac_n(:, k, 1, 1) = D(k)*X(:)/dsqrt(delta(:,k)-1d0 +1d-50)
        jac_n(:, k, 2, 2) = C(:)*X(:)/dsqrt(delta(:,k)*(delta(:,k)-1d0 +1d-50))
        ! off-diagonal terms
        jac_n(:, k, 1, 2) = - D(k)*Y(k)/dsqrt(delta(:,k)*(delta(:,k)-1d0 +1d-50))
        jac_n(:, k, 2, 1) =  C(:)*Y(k)/dsqrt(delta(:,k)-1d0 +1d-50)
      
        ! south jacobian matrix
        jac_s(:, k, 1, 1) = - D(k)*X(:)/dsqrt(delta(:,k)-1d0 +1d-50)
        jac_s(:, k, 2, 2) = - C(:)*X(:)/dsqrt(delta(:,k)*(delta(:,k)-1d0 +1d-50))
        jac_s(:, k, 1, 2) =   D(k)*Y(k)/dsqrt(delta(:,k)*(delta(:,k)-1d0 +1d-50))
        jac_s(:, k, 2, 1) = - C(:)*Y(k)/dsqrt(delta(:,k)-1d0 +1d-50)
      
        ! inverse of equator jacobian terms
        jacinv_eq(:, k, 1, 1) = X(:)*Y(k)/(C(:)*D(k))
        jacinv_eq(:, k, 2, 1) = dsqrt(delta(:,k))/(C(:)*D(k))
      
        ! inverse of north jacobian matrix
        jacinv_n(:, k, 1, 1) = X(:)/(D(k)*dsqrt(delta(:,k)-1d0 +1d-50))
        jacinv_n(:, k, 2, 2) = X(:)*dsqrt(delta(:,k))/(C(:)*dsqrt(delta(:,k)-1d0 +1d-50))
        jacinv_n(:, k, 1, 2) = Y(k)/(C(:)*dsqrt(delta(:,k)-1d0+1d-50))
        jacinv_n(:, k, 2, 1) = - Y(k)*dsqrt(delta(:,k))/(D(k)*dsqrt(delta(:,k)-1d0 +1d-50))

        ! inverse of south jacobian matrix
        jacinv_s(:, k, 1, 1) = - X(:)/(D(k)*dsqrt(delta(:,k)-1d0) +1d-50 )
        jacinv_s(:, k, 2, 2) = - X(:)*dsqrt(delta(:,k))/(C(:)*dsqrt(delta(:,k)-1d0 +1d-50))
        jacinv_s(:, k, 1, 2) = - Y(k)/(C(:)*dsqrt(delta(:,k)-1d0) +1d-50)
        jacinv_s(:, k, 2, 1) =  Y(k)*dsqrt(delta(:,k))/(D(k)*dsqrt(delta(:,k)-1d0 +1d-50))
      end do

      ! On the axis the Jacobian is singular, and is never used indeed
      jac_n((nang+1)/2,(nang+1)/2,1,1) = 1d0
      jac_s((nang+1)/2,(nang+1)/2,1,1) = 1d0
      jacinv_n((nang+1)/2,(nang+1)/2,1,1) = 1d0
      jacinv_s((nang+1)/2,(nang+1)/2,1,1) = 1d0

      jac_n((nang+1)/2,(nang+1)/2,1,2) = 0d0
      jac_s((nang+1)/2,(nang+1)/2,1,2) = 0d0
      jacinv_n((nang+1)/2,(nang+1)/2,1,2) = 0d0
      jacinv_s((nang+1)/2,(nang+1)/2,1,2) = 0d0

      jac_n((nang+1)/2,(nang+1)/2,2,1) = 0d0
      jac_s((nang+1)/2,(nang+1)/2,2,1) = 0d0
      jacinv_n((nang+1)/2,(nang+1)/2,2,1) = 0d0
      jacinv_s((nang+1)/2,(nang+1)/2,2,1) = 0d0

      jac_n((nang+1)/2,(nang+1)/2,2,2) = 1d0
      jac_s((nang+1)/2,(nang+1)/2,2,2) = 1d0
      jacinv_n((nang+1)/2,(nang+1)/2,2,2) = 1d0
      jacinv_s((nang+1)/2,(nang+1)/2,2,2) = 1d0

    end subroutine jacobians

!!--------------------------------------------------------------------------
    !> @brief subroutine jacinv_edge
    !!
    !! TBD: YOU DON'T NEED TO DO ALL THE CALCULATION, ONLY AT THE EDGE
    !!      SO YOU AVOID ALSO THE /0 PROBLEMS
    !! In this subroutine, we are defining the inverse of the jacobian matrices,
    !! introduced in Ronchi's formalism of 1996, in term of xi_a and eta_a variables.
    !!
    !! The inverse of the Jacobian matrices are used to go from cubed sphere to spherical coordinates.
    !! These jacobians matrices written in term of xi_a and eta_a, which correspond to the interpolated value 
    !! of xi and eta at ghost points in the coordinate system of the adjacent patch, are 
    !! very important to define the jacobian needed to change the vector quantities from one patch to the other.
    !!
    !! Code owners:
    !!    Clara Dehman
!!--------------------------------------------------------------------------
    subroutine jacinv_edge(m,xi_a,eta_a,jac)
     implicit none
      integer :: k
      integer, intent(in) :: m
      real*8, dimension(0:nang+1), intent(in) :: xi_a  ! xi adjacent patch
      real*8, dimension(0:nang+1), intent(in) :: eta_a ! eta adjacent patch
      real*8, dimension (0:nang+1,0:nang+1,2,2), intent(out) :: jac 
      ! adjacent patch
      real*8, dimension(0:nang+1) :: X_a  ! tan(xi_a)
      real*8, dimension(0:nang+1) :: Y_a  ! tan(eta_a)
      real*8, dimension(0:nang+1) :: C_a  ! 1/cos(xi_a)
      real*8, dimension(0:nang+1) :: D_a  ! 1/cos(eta_a)
      real*8, dimension(0:nang+1, 0:nang+1) :: delta_a ! 1 + X_a**2 + Y_a**2
      X_a(:) = dtan(xi_a(:))
      Y_a(:) = dtan(eta_a(:))
      C_a(:) = (1 + X_a(:)**2)**(0.5d0)
      D_a(:) = (1 + Y_a(:)**2)**(0.5d0)

      do k = 0, nang+1
       delta_a(:, k) = 1 + X_a(:)**2 + Y_a(k)**2
      end do
      

    ! jacinv_eq at qm and pm 
      if (m == 1) then 
      jac(:, :, 2, 2) = 0
      jac(:, :, 1, 2) = -1
      do k = 0, nang+1
      jac(:, k, 1, 1) = X_a(:)*Y_a(k)/(C_a(:)*D_a(k))
      jac(:, k, 2, 1) = dsqrt(delta_a(:,k))/(C_a(:)*D_a(k))
      end do
      end if

      ! jacinv_n at qm and pm 
      if (m == 2) then 
      do k = 0, nang+1
      jac(:, k, 1, 1) = X_a(:)/(D_a(k)*dsqrt(delta_a(:,k)-1d0+1d-50))
      jac(:, k, 2, 2) = X_a(:)*dsqrt(delta_a(:,k))/(C_a(:)*dsqrt(delta_a(:,k)-1d0+1d-50))
      jac(:, k, 1, 2) = Y_a(k)/(C_a(:)*dsqrt(delta_a(:,k)-1d0+1d-50))
      jac(:, k, 2, 1) = - Y_a(k)*dsqrt(delta_a(:,k))/(D_a(k)*dsqrt(delta_a(:,k)-1d0+1d-50))
      end do
      end if

      ! jacinv_s at qm and pm 
      if (m == 3) then
      do k = 0, nang+1
      jac(:, k, 1, 1) = - X_a(:)/(D_a(k)*dsqrt(delta_a(:,k)-1d0+1d-50))
      jac(:, k, 2, 2) = - X_a(:)*dsqrt(delta_a(:,k))/(C_a(:)*dsqrt(delta_a(:,k)-1d0+1d-50))
      jac(:, k, 1, 2) = - Y_a(k)/(C_a(:)*dsqrt(delta_a(:,k)-1d0+1d-50))
      jac(:, k, 2, 1) = Y_a(k)*dsqrt(delta_a(:,k))/(D_a(k)*dsqrt(delta_a(:,k)-1d0+1d-50))
      end do
      end if

    end subroutine jacinv_edge

!!--------------------------------------------------------------------------
    !> @brief subroutine jacobian_edges
    !!
    !! In this subroutine, we are determining the 18 (9+9) series of jacobians needed 
    !! to change the cooridnates from one patch to the other (passing by spherical 
    !! coordinates). 
    !!
    !! We go from the adjacent patch to spherical coordinates using the inverse of the 
    !! jacobians, and from spherical coordinates to the coordinate sytem of the original
    !! patch using the jacobian matrices.
    !!
    !!>--------------------------------------------------------------------
    !!> m = 1 refers to jacinv_eq 
    !!> m = 2 refers to jacinv_n 
    !!> m = 3 refers to jacinv_s
    !!>--------------------------------------------------------------------
    !!
    !! Code owners:
    !!    Clara Dehman
    !!    Daniele Viganò
!!--------------------------------------------------------------------------   
    subroutine jacobian_edges 
      implicit none
      integer :: m
      real*8, dimension(0:nang+1) :: xi_adj  ! xi adjacent patch
      real*8, dimension(0:nang+1) :: eta_adj ! eta adjacent patch
      real*8, dimension (0:nang+1,0:nang+1,2,2) :: jacinv 

      !>--------------------------------------------------------------------
      !>-------------------- equator with equator patch --------------------
        xi_adj = xi
        eta_adj = pm 
        m = 1
        call jacinv_edge(m,xi_adj,eta_adj,jacinv)

      ! II -> I, III -> II, IV -> III, I -> IV
      ! e.g., patch II refers to the adjacent patch and patch I to the original patch
        jac_pp(1,:,:,:) = jac_edge(jac_eq(nang+1,:,:,:),jacinv(2,:,:,:))
  
      ! I -> II, II -> III, III -> IV, IV -> I
      ! e.g., patch I refers to the adjacent patch and patch II to the original patch
        jac_pp(2,:,:,:) = jac_edge(jac_eq(0,:,:,:),jacinv(nang-1,:,:,:))
      !>--------------------------------------------------------------------

      !>--------------------------------------------------------------------
      !>-------------------- patch V - VI with patch I - III ---------------

        xi_adj = pm
        eta_adj = eta

      ! jacinv_n
        m = 2
        call jacinv_edge(m,xi_adj,eta_adj,jacinv)  
      ! I -> V
      ! e.g., patch V refers to the adjacent patch and patch I to the original patch 
        jac_pp(3,:,:,:) = jac_edge(jac_eq(:,nang+1,:,:),jacinv(:,2,:,:))
      ! III -> V
      ! e.g., patch V refers to the adjacent patch and patch III to the original patch 
        jac_pp(7,:,:,:) = jac_edge(jac_eq(:,nang+1,:,:),jacinv(nang+1:0:-1,nang-1,:,:))

      ! jacinv_s
        m = 3 
        call jacinv_edge(m,xi_adj,eta_adj,jacinv)  
      ! I -> VI
      ! e.g., patch VI refers to the adjacent patch and patch I to the original patch 
        jac_pp(11,:,:,:) = jac_edge(jac_eq(:,0,:,:),jacinv(:,nang-1,:,:))
      ! III -> VI
      ! e.g., patch VI refers to the adjacent patch and patch III to the original patch
        jac_pp(15,:,:,:) = jac_edge(jac_eq(:,0,:,:),jacinv(nang+1:0:-1,2,:,:))

      ! jacinv_eq
        m = 1
        call jacinv_edge(m,xi_adj,eta_adj,jacinv) 
      ! V -> I
      ! e.g., patch I refers to the adjacent patch and patch V to the original patch 
        jac_pp(4,:,:,:) = jac_edge(jac_n(:,0,:,:),jacinv(:,nang-1,:,:))
      ! V -> III
      ! e.g., patch I refers to the adjacent patch and patch V to the original patch 
        jac_pp(8,:,:,:) = jac_edge(jac_n(nang+1:0:-1,nang+1,:,:),jacinv(:,nang-1,:,:))
      ! VI -> I
      ! e.g., patch I refers to the adjacent patch and patch VI to the original patch 
        jac_pp(12,:,:,:) = jac_edge(jac_s(:,nang+1,:,:),jacinv(:,2,:,:))
      ! VI -> III
      ! e.g., patch I refers to the adjacent patch and patch VI to the original patch 
        jac_pp(16,:,:,:) = jac_edge(jac_s(nang+1:0:-1,0,:,:),jacinv(:,2,:,:))
      !>--------------------------------------------------------------------  

      !>--------------------------------------------------------------------
      !>-------------------- patch V - VI with patch II - IV ---------------

        xi_adj = xi
        eta_adj = pm 

      ! jacinv_n
        m = 2  
        call jacinv_edge(m,xi_adj,eta_adj,jacinv)
      ! II -> V
      ! e.g., patch V refers to the adjacent patch and patch II to the original patch
        jac_pp(5,:,:,:) = jac_edge(jac_eq(:,nang+1,:,:),jacinv(nang-1,:,:,:))
      ! IV -> V
      ! e.g., patch V refers to the adjacent patch and patch IV to the original patch 
        jac_pp(9,:,:,:) = jac_edge(jac_eq(:,nang+1,:,:),jacinv(2,nang+1:0:-1,:,:))
 
      ! jacinv_s
        m = 3 
        call jacinv_edge(m,xi_adj,eta_adj,jacinv)
      ! II -> VI
      ! e.g., patch VI refers to the adjacent patch and patch II to the original patch 
        jac_pp(13,:,:,:) = jac_edge(jac_eq(:,0,:,:),jacinv(nang-1,nang+1:0:-1,:,:))
      ! IV -> VI
      ! e.g., patch VI refers to the adjacent patch and patch IV to the original patch
        jac_pp(17,:,:,:) = jac_edge(jac_eq(:,0,:,:),jacinv(2,:,:,:))

      ! jacinv_eq
        xi_adj = pm
        eta_adj = eta
        m = 1  
        call jacinv_edge(m,xi_adj,eta_adj,jacinv)
      ! V -> II
      ! e.g., patch IV refers to the adjacent patch and patch V to the original patch 
        jac_pp(6,:,:,:) = jac_edge(jac_n(nang+1,:,:,:),jacinv(:,nang-1,:,:))
      ! V -> IV
      ! e.g., patch IV refers to the adjacent patch and patch V to the original patch 
        jac_pp(10,:,:,:) = jac_edge(jac_n(0,nang+1:0:-1,:,:),jacinv(:,nang-1,:,:))
      ! VI -> II
      ! e.g., patch I refers to the adjacent patch and patch VI to the original patch 
        jac_pp(14,:,:,:) = jac_edge(jac_s(nang+1,nang+1:0:-1,:,:),jacinv(:,2,:,:))
      ! VI -> IV
      ! e.g., patch I refers to the adjacent patch and patch VI to the original patch 
        jac_pp(18,:,:,:) = jac_edge(jac_s(0,:,:,:),jacinv(:,2,:,:))
      !>--------------------------------------------------------------------
      end subroutine jacobian_edges
  
  
      function jac_edge(jac_xx,jacinv_xx) result(jac)
  
        implicit none
        ! jac_xx is the jacobian needed to go from spherical to the original patch
        ! jacinv_xx is the jacobian needed to from the adjacent patch to spherical coordinates          
        real*8, dimension(0:nang+1,2,2), intent(in) :: jac_xx, jacinv_xx 
        real*8, dimension(0:nang+1,2,2) :: jac
  
        jac(:,1,1) = jac_xx(:,1,1)*jacinv_xx(:,1,1) + jac_xx(:,1,2)*jacinv_xx(:,2,1)
        jac(:,1,2) = jac_xx(:,1,1)*jacinv_xx(:,1,2) + jac_xx(:,1,2)*jacinv_xx(:,2,2)
        jac(:,2,1) = jac_xx(:,2,1)*jacinv_xx(:,1,1) + jac_xx(:,2,2)*jacinv_xx(:,2,1)
        jac(:,2,2) = jac_xx(:,2,1)*jacinv_xx(:,1,2) + jac_xx(:,2,2)*jacinv_xx(:,2,2)
  
      end function jac_edge

!!--------------------------------------------------------------------------
    !> @brief subroutine fghost
    !!
    !! In this subroutine reads as an input the three vector components
    !! fr, fxi and feta in the coordinate system of the adjacent patch. Then, we 
    !! apply a change of coordinates from adjacent to original patch coordinates.
    !! Later, we interpolate by calling interpolation_edges subroutine, in order to obtain the values 
    !! fr, fxi, and feta at the ghost cells, in the coordinate system of the original patch.
    !!
    !! This subroutine returns the value of fr, fxi and feta at the ghost points and in the coordinate 
    !! system of the original patch.
    !!
    !! f1, f2, and f3 are the field components of the interior points along the set of ghost points 
    !! and they are expressed in the coordinate system of the adjacent patch.
    !! fg1, fg2, and fg3 are the field components at the ghost cells and they are expressed in the coordinate 
    !! system of the original patch
    !!
    !! Code owners:
    !!    Clara Dehman
    !!    Daniele Viganò
!!--------------------------------------------------------------------------
     subroutine fghost(fr,fxi,feta)

      implicit None

      real*8, dimension(0:nr+1,0:nang+1,0:nang+1,6), intent(inout) :: fxi, feta, fr
      real*8, dimension(0:nr+1,0:nang+1) :: f1, f2, f3, fg1, fg2, fg3
      integer pp
 
      pp = 1

      f1 = fxi(:,2,0:nang+1,2)
      f2 = feta(:,2,0:nang+1,2)
      f3 = fr(:,2,0:nang+1,2)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,nang+1,0:nang+1,1) = fg1
      feta(:,nang+1,0:nang+1,1) = fg2
      fr(:,nang+1,0:nang+1,1) = fg3

      f1 = fxi(:,2,0:nang+1,3)
      f2 = feta(:,2,0:nang+1,3)
      f3 = fr(:,2,0:nang+1,3)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,nang+1,0:nang+1,2) = fg1
      feta(:,nang+1,0:nang+1,2) = fg2
      fr(:,nang+1,0:nang+1,2) = fg3

      f1 = fxi(:,2,0:nang+1,4)
      f2 = feta(:,2,0:nang+1,4)
      f3 = fr(:,2,0:nang+1,4)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,nang+1,0:nang+1,3) = fg1
      feta(:,nang+1,0:nang+1,3) = fg2
      fr(:,nang+1,0:nang+1,3) = fg3
      
      f1 = fxi(:,2,0:nang+1,1)
      f2 = feta(:,2,0:nang+1,1)
      f3 = fr(:,2,0:nang+1,1)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,nang+1,0:nang+1,4) = fg1
      feta(:,nang+1,0:nang+1,4) = fg2
      fr(:,nang+1,0:nang+1,4) = fg3

      pp = 2

      f1 = fxi(:,nang-1,0:nang+1,1)
      f2 = feta(:,nang-1,0:nang+1,1)
      f3 = fr(:,nang-1,0:nang+1,1)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0,0:nang+1,2) = fg1
      feta(:,0,0:nang+1,2) = fg2
      fr(:,0,0:nang+1,2)= fg3

      f1 = fxi(:,nang-1,0:nang+1,2)
      f2 = feta(:,nang-1,0:nang+1,2)
      f3 = fr(:,nang-1,0:nang+1,2)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0,0:nang+1,3) = fg1
      feta(:,0,0:nang+1,3) = fg2
      fr(:,0,0:nang+1,3) = fg3

      f1 = fxi(:,nang-1,0:nang+1,3)
      f2 = feta(:,nang-1,0:nang+1,3)
      f3 = fr(:,nang-1,0:nang+1,3)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0,0:nang+1,4) = fg1
      feta(:,0,0:nang+1,4) = fg2
      fr(:,0,0:nang+1,4) = fg3
      
      f1 = fxi(:,nang-1,0:nang+1,4)
      f2 = feta(:,nang-1,0:nang+1,4)
      f3 = fr(:,nang-1,0:nang+1,4)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0,0:nang+1,1) = fg1
      feta(:,0,0:nang+1,1) = fg2
      fr(:,0,0:nang+1,1) = fg3

      pp = 3
      f1 = fxi(:,0:nang+1,2,5)
      f2 = feta(:,0:nang+1,2,5)
      f3= fr(:,0:nang+1,2,5)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0:nang+1,nang+1,1) = fg1
      feta(:,0:nang+1,nang+1,1) = fg2
      fr(:,0:nang+1,nang+1,1) = fg3

      pp = 4
      f1 = fxi(:,0:nang+1,nang-1,1)
      f2 = feta(:,0:nang+1,nang-1,1)
      f3 = fr(:,0:nang+1,nang-1,1)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0:nang+1,0,5) = fg1
      feta(:,0:nang+1,0,5) = fg2
      fr(:,0:nang+1,0,5) = fg3

      pp = 5
      f1 = fxi(:,nang-1,0:nang+1,5)
      f2 = feta(:,nang-1,0:nang+1,5)
      f3 = fr(:,nang-1,0:nang+1,5)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0:nang+1,nang+1,2) = fg1
      feta(:,0:nang+1,nang+1,2) = fg2
      fr(:,0:nang+1,nang+1,2) = fg3

      pp = 6
      f1 = fxi(:,0:nang+1,nang-1,2)
      f2 = feta(:,0:nang+1,nang-1,2)
      f3 = fr(:,0:nang+1,nang-1,2)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,nang+1,0:nang+1,5) = fg1
      feta(:,nang+1,0:nang+1,5) = fg2
      fr(:,nang+1,0:nang+1,5) = fg3

      pp = 7
      f1 = fxi(:,nang+1:0:-1,nang-1,5)
      f2 = feta(:,nang+1:0:-1,nang-1,5) 
      f3 = fr(:,nang+1:0:-1,nang-1,5)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0:nang+1,nang+1,3) = fg1
      feta(:,0:nang+1,nang+1,3) = fg2
      fr(:,0:nang+1,nang+1,3) = fg3

      pp = 8
      f1 = fxi(:,0:nang+1,nang-1,3)
      f2 = feta(:,0:nang+1,nang-1,3)
      f3 = fr(:,0:nang+1,nang-1,3)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,nang+1:0:-1,nang+1,5) = fg1
      feta(:,nang+1:0:-1,nang+1,5) = fg2
      fr(:,nang+1:0:-1,nang+1,5) = fg3

      pp = 9
      f1 = fxi(:,2,nang+1:0:-1,5)
      f2 = feta(:,2,nang+1:0:-1,5) 
      f3 = fr(:,2,nang+1:0:-1,5) 
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0:nang+1,nang+1,4) = fg1
      feta(:,0:nang+1,nang+1,4) = fg2
      fr(:,0:nang+1,nang+1,4) = fg3

      pp = 10 
      f1 = fxi(:,0:nang+1,nang-1,4)
      f2 = feta(:,0:nang+1,nang-1,4)
      f3 = fr(:,0:nang+1,nang-1,4)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0,nang+1:0:-1,5) = fg1
      feta(:,0,nang+1:0:-1,5) = fg2
      fr(:,0,nang+1:0:-1,5) = fg3

      pp = 11
      f1 = fxi(:,0:nang+1,nang-1,6)
      f2 = feta(:,0:nang+1,nang-1,6)
      f3 = fr(:,0:nang+1,nang-1,6)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0:nang+1,0,1) = fg1
      feta(:,0:nang+1,0,1) = fg2
      fr(:,0:nang+1,0,1) = fg3
 
      pp = 12
      f1 = fxi(:,0:nang+1,2,1)
      f2 = feta(:,0:nang+1,2,1)
      f3 = fr(:,0:nang+1,2,1)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0:nang+1,nang+1,6) = fg1
      feta(:,0:nang+1,nang+1,6) = fg2
      fr(:,0:nang+1,nang+1,6) = fg3

      pp = 13 
      f1 = fxi(:,nang-1,nang+1:0:-1,6)
      f2 = feta(:,nang-1,nang+1:0:-1,6) 
      f3 = fr(:,nang-1,nang+1:0:-1,6)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0:nang+1,0,2) = fg1
      feta(:,0:nang+1,0,2) = fg2
      fr(:,0:nang+1,0,2) = fg3
       
      pp = 14
      f1 = fxi(:,0:nang+1,2,2)
      f2 = feta(:,0:nang+1,2,2)
      f3 = fr(:,0:nang+1,2,2)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,nang+1,nang+1:0:-1,6) = fg1
      feta(:,nang+1,nang+1:0:-1,6) = fg2
      fr(:,nang+1,nang+1:0:-1,6) = fg3

      pp = 15
      f1 = fxi(:,nang+1:0:-1,2,6)
      f2 = feta(:,nang+1:0:-1,2,6)
      f3 = fr(:,nang+1:0:-1,2,6)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0:nang+1,0,3) = fg1
      feta(:,0:nang+1,0,3) = fg2
      fr(:,0:nang+1,0,3) = fg3

      pp = 16
      f1 = fxi(:,0:nang+1,2,3)
      f2 = feta(:,0:nang+1,2,3)
      f3 = fr(:,0:nang+1,2,3)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,nang+1:0:-1,0,6) = fg1
      feta(:,nang+1:0:-1,0,6) = fg2
      fr(:,nang+1:0:-1,0,6) =fg3

      pp = 17
      f1 = fxi(:,2,0:nang+1,6)
      f2 = feta(:,2,0:nang+1,6)
      f3 = fr(:,2,0:nang+1,6)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0:nang+1,0,4) = fg1
      feta(:,0:nang+1,0,4) = fg2
      fr(:,0:nang+1,0,4) = fg3

      pp = 18
      f1 = fxi(:,0:nang+1,2,4)
      f2 = feta(:,0:nang+1,2,4)
      f3 = fr(:,0:nang+1,2,4)
      call interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)
      fxi(:,0,0:nang+1,6) = fg1
      feta(:,0,0:nang+1,6) = fg2
      fr(:,0,0:nang+1,6) = fg3

   end subroutine fghost

!!--------------------------------------------------------------------------
    !> @brief subroutine edges
    !!
    !! In this subroutine, we define the weights needed to linearly interpolate
    !! the values of an adjacent patch to have the values at the ghost points. 
    !! The weights are the same for all edges, for both sides.
    !! It consist of a vector, that gives the weight 1-edge_w(i) (from 0 to 1)
    !! related to the point having the same index in the adjacent patch.
    !! Then, there will be a second neighbour, which is the one toward the center of the vector,
    !! for which the weight is (edge_w(i))
    !!
    !! Code owners:
    !!    Daniele Viganò
    !!    Stefano Ascenzi
!!--------------------------------------------------------------------------
    subroutine edges()

      implicit None

      integer i, it
      real*8 da, dat

      da = xi(1) - xi(0)
      dat = xi(2) - xi(0)
    !  pm = datan( dtan(xi)/dtan(0.25d0*PI + da) )

      edge_w = 0d0
      edge_wt = 0d0

      ! Central cell (nang is odd by definition): the weight is 0, since the points coincide
      edge_w(nang/2+1) = 0d0
      edge_wt(nang/2+1) = 0d0

      ! Intermediate points: the weight is between 0 and 1 excluded
      do i = 0, nang/2
        edge_w(i) = (pm(i) - xi(i)) / da
        edge_w(nang+1-i) = (xi(nang+1-i) - pm(nang+1-i)) / da

        ! WARNING: Check better the following
        !edge_wt(i) = (pmt(i) - xi(i)) / dat
        !edge_wt(nang+1-i) = (xi(nang+1-i) - pmt(nang+1-i)) / dat
      enddo

        !new implementation

      do i =1, nang/2
        edge_wt(i) = 0.5*(1.d0 - abs(pmt(i-1) - xi(i+1)) / dat)
        edge_wt(nang+1-i) = 0.5*(1.d0 - abs(xi(nang+1-i-1) - pmt(nang+1-i+1)) / dat)
      end do

      !do i =1, nangt
      !  print*, i, edge_wt(2*i)
      !enddo

      ! the same for the thermal grid

      !do i = 0, nang/2
      !  it = 2*i
      !  edge_wt(i) = (pm(it) - xi(it)) / dat
      !  edge_wt(nangt+1-i) = (xi(nang+1-it) - pm(nang+1-it)) / dat
      !end do
   
    end subroutine edges

!!--------------------------------------------------------------------------
    !> @brief subroutine interpolation_edges
    !!
    !! f1, f2, and f3 are the field components of the interior points along the set of ghost points 
    !! and they are expressed in the coordinate system of the adjacent patch.
    !!
    !! ft1 and ft2 are the the field components of the interior along the set of ghost points, but they are 
    !! expressed in the coordinate system of the original patch.
    !!
    !! fg1, fg2, and fg3 are the field components at the ghost cell and they are expressed in the coordinate 
    !! system of the original patch.
    !!
    !! NB: For the radial component, f3, no change of coordinates from one patch to the other is applied
    !!
    !! Code owners:
    !!    Daniele Viganò
    !!    Clara Dehman
!!--------------------------------------------------------------------------
    subroutine interpolation_edges(pp,f1,f2,f3,fg1,fg2,fg3)

      implicit none
      integer, intent(in) :: pp 
      real*8, dimension(0:nr+1,0:nang+1), intent(in) :: f1, f2, f3
      real*8, dimension(0:nr+1,0:nang+1), intent(out) ::  fg1, fg2, fg3
      real*8, dimension(0:nr+1,0:nang+1) :: ft1, ft2
      real*8, dimension(0:nang+1,2,2) :: jac
      integer j

      ft1 = 0d0
      ft2 = 0d0
      fg1 = 0d0
      fg2 = 0d0
      fg3 = 0d0

    ! interpolation 
      do j = 0, nang/2
      ! Values in the coordinates of the original patch, for which we want the ghost point values
        ft1(:,j) = f1(:,j)*(1d0-edge_w(j)) + f1(:,j+1)*edge_w(j)
        ft2(:,j) = f2(:,j)*(1d0-edge_w(j)) + f2(:,j+1)*edge_w(j)
        fg3(:,j) = f3(:,j)*(1d0-edge_w(j)) + f3(:,j+1)*edge_w(j)
        ! j is replaced by nang+1-j and j+1 by nang-j
        ft1(:,nang+1-j) = f1(:,nang+1-j)*(1d0-edge_w(nang+1-j)) + f1(:,nang-j)*edge_w(nang+1-j)
        ft2(:,nang+1-j) = f2(:,nang+1-j)*(1d0-edge_w(nang+1-j)) + f2(:,nang-j)*edge_w(nang+1-j)
        fg3(:,nang+1-j) = f3(:,nang+1-j)*(1d0-edge_w(nang+1-j)) + f3(:,nang-j)*edge_w(nang+1-j)
      end do

      ! In the middle of the edge
      ft1(:,nang/2+1) = f1(:,nang/2+1)
      ft2(:,nang/2+1) = f2(:,nang/2+1)
      fg3(:,nang/2+1) = f3(:,nang/2+1)

      jac = jac_pp(pp,:,:,:)

      do j = 0, nang + 1
        fg1(:,j) = jac(j,1,1)*ft1(:,j) + jac(j,1,2)*ft2(:,j)
        fg2(:,j) = jac(j,2,1)*ft1(:,j) + jac(j,2,2)*ft2(:,j)
      end do

    end subroutine interpolation_edges

    !!------------------------------------------------------------------------------
    !> @brief This subroutine provide the contravariant component of the curl of a vector
    !> using a finite volume method. 
    !> For the length elements, we are using the contravariant components.
    !> For the surface elements, we are using the convariant components.
    !! DV: Remind me why don't we use the length covariant,
    !!     so that it could be easier, without off-diag terms?
    !!
    !! @param[in]  lr, lxi, leta             contravariant length elements 
    !! @param[in]  area_r, area_xi, area_eta covariant surface elements 
    !! @param[in]   fr                       contravariant radial component
    !! @param[in]   fxi                      contravariant xi component
    !! @param[in]   feta                     contravariant eta component
    !! @param[out]  curlfr                   contravariant radial component of curlf 
    !! @param[out]  curlfxi                  contravariant xi component of curlf
    !! @param[out]  curlfeta                 contravariant eta component of curlf
    !!
    !!
    !! Code owners:
    !!   Clara Dehman
    !!------------------------------------------------------------------------------
    subroutine curl_fnvol(fr,fxi,feta,curlfr,curlfxi,curlfeta,imin)
  
      implicit none
  
      ! Subroutine arguments -------------------------------------------------------
      integer, intent(in) :: imin 
      real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: fr,fxi,feta
      real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(out) :: curlfr, curlfxi, curlfeta

      ! Local variables ------------------------------------------------------------
      ! Auxiliary index for loops
      integer :: i,j,k,p

      curlfr = 0.d0
      curlfxi = 0.d0
      curlfeta = 0.d0
      
      !$OMP Parallel shared(curlfr, curlfxi, curlfeta)

      !$OMP do collapse(2) private(p,j,k)
      do p=1,6
        do k=1,nang
          do j=1,nang
          curlfr(imin:nr,j,k,p) =(leta(imin:nr,j+1,k)*(feta(imin:nr,j+1,k,p) - X(j+1)*Y(k)/(C(j+1)*D(k))*fxi(imin:nr,j+1,k,p)) &
          & - leta(imin:nr,j-1,k)*(feta(imin:nr,j-1,k,p) - X(j-1)*Y(k)/(C(j-1)*D(k))*fxi(imin:nr,j-1,k,p)) & 
          & + lxi(imin:nr,j,k-1)*(fxi(imin:nr,j,k-1,p) - X(j)*Y(k-1)/(C(j)*D(k-1))*feta(imin:nr,j,k-1,p)) &
          & - lxi(imin:nr,j,k+1)*(fxi(imin:nr,j,k+1,p) - X(j)*Y(k+1)/(C(j)*D(k+1))*feta(imin:nr,j,k+1,p)))/area_r(imin:nr,j,k)
          enddo
        enddo
        enddo
        !$OMP end do   
        !$OMP do collapse(2) private(p,i,k)
        do p=1,6
          do k=1,nang
            do i=imin,nr
           curlfxi(i,1:nang,k,p)= (lr(i)*(fr(i,1:nang,k+1,p)-fr(i,1:nang,k-1,p)) &
          & + leta(i-1,1:nang,k)*(feta(i-1,1:nang,k,p)-X(1:nang)*Y(k)/(C(1:nang)*D(k))*fxi(i-1,1:nang,k,p)) &
          & - leta(i+1,1:nang,k)*(feta(i+1,1:nang,k,p)-X(1:nang)*Y(k)/(C(1:nang)*D(k))*fxi(i+1,1:nang,k,p)))/area_xi(i,1:nang,k)
            enddo
          enddo
        enddo
        !$OMP end do 
        !$OMP do collapse(2) private(p,i,j)
        do p=1,6
          do j=1,nang
            do i=imin,nr
            curlfeta(i,j,1:nang,p) = (lr(i)*(fr(i,j-1,1:nang,p)-fr(i,j+1,1:nang,p)) &
            & + lxi(i+1,j,1:nang)*(fxi(i+1,j,1:nang,p)-X(j)*Y(1:nang)/(C(j)*D(1:nang))*feta(i+1,j,1:nang,p)) &
            & - lxi(i-1,j,1:nang)*(fxi(i-1,j,1:nang,p)-X(j)*Y(1:nang)/(C(j)*D(1:nang))*feta(i-1,j,1:nang,p)))/area_eta(i,j,1:nang)
            enddo
          enddo
        enddo
       !$OMP end do 
       !$OMP end Parallel

       ! call cpu_time(t2)
       ! print*, t2-t1

       call edge_average(imin,curlfr,curlfxi,curlfeta)

    end subroutine curl_fnvol

!!--------------------------------------------------------------------------
    !> @brief subroutine jacobian_average
    !!
    !! In this subroutine, we are determining the 18 (9+9) series of jacobians at the 
    !! interfaces needed to average the angular component of the field 
    !!
    !! We go from the adjacent patch to spherical coordinates using the inverse of the 
    !! jacobians, and from spherical coordinates to the coordinate sytem of the original
    !! patch using the jacobian matrices.
    !!
    !! Code owners:
    !!    Clara Dehman
!!--------------------------------------------------------------------------   
    subroutine jacobian_average 
      implicit none

      !>--------------------------------------------------------------------
      !>-------------------- equator with equator patch --------------------
      ! II -> I, III -> II, IV -> III, I -> IV
      ! e.g., from the coordinate system of patch II to the coordinate system of patch I 
        jac_avg(1,:,:,:) = jac_edge(jac_eq(nang,:,:,:),jacinv_eq(1,:,:,:)) 
  
      ! I -> II, II -> III, III -> IV, IV -> I
      ! e.g., from the coordinate system of patch I to the coordinate system of patch II
        jac_avg(2,:,:,:) = jac_edge(jac_eq(1,:,:,:), jacinv_eq(nang,:,:,:))

      !>-------------------- equator with north patch ----------------------
      ! V -> I
      jac_avg(3,:,:,:) = jac_edge(jac_eq(:,nang,:,:),jacinv_n(:,1,:,:))   
      ! I -> V
      jac_avg(4,:,:,:) = jac_edge(jac_n(:,1,:,:), jacinv_eq(:,nang,:,:)) 
      ! V -> II 
      jac_avg(5,:,:,:) = jac_edge(jac_eq(:,nang,:,:),jacinv_n(nang,:,:,:))  
      ! II -> V 
      jac_avg(6,:,:,:) = jac_edge(jac_n(nang,:,:,:), jacinv_eq(:,nang,:,:))
      ! V -> III
      jac_avg(7,:,:,:) = jac_edge(jac_eq(:,nang,:,:),jacinv_n(nang+1:0:-1,nang,:,:))  
      ! III -> V 
      jac_avg(8,:,:,:) = jac_edge(jac_n(nang+1:0:-1,nang,:,:), jacinv_eq(:,nang,:,:))
      ! V -> IV 
      jac_avg(9,:,:,:) = jac_edge(jac_eq(:,nang,:,:),jacinv_n(1,nang+1:0:-1,:,:)) 
      ! IV -> V
      jac_avg(10,:,:,:) = jac_edge(jac_n(1,nang+1:0:-1,:,:), jacinv_eq(:,nang,:,:)) 

      !>-------------------- equator with south patch ----------------------
      ! VI -> I
      jac_avg(11,:,:,:) = jac_edge(jac_eq(:,1,:,:),jacinv_s(:,nang,:,:))
      ! I -> VI
      jac_avg(12,:,:,:) = jac_edge(jac_s(:,nang,:,:), jacinv_eq(:,1,:,:))
      ! VI -> II
      jac_avg(13,:,:,:) = jac_edge(jac_eq(:,1,:,:),jacinv_s(nang,nang+1:0:-1,:,:)) 
      ! II -> VI         
      jac_avg(14,:,:,:) = jac_edge(jac_s(nang,nang+1:0:-1,:,:), jacinv_eq(:,1,:,:))
      ! VI -> III
      jac_avg(15,:,:,:) = jac_edge(jac_eq(:,1,:,:),jacinv_s(nang+1:0:-1,1,:,:))
      ! III -> VI 
      jac_avg(16,:,:,:) = jac_edge(jac_s(nang+1:0:-1,1,:,:), jacinv_eq(:,1,:,:)) 
      ! VI -> IV
      jac_avg(17,:,:,:) = jac_edge(jac_eq(:,1,:,:),jacinv_s(1,:,:,:))   
      ! IV -> VI 
      jac_avg(18,:,:,:) = jac_edge(jac_s(1,:,:,:), jacinv_eq(:,1,:,:))

       end subroutine jacobian_average


  !!------------------------------------------------------------------------------
    !> @brief This subroutine read as an input the contravariant components of a vector 
    !> at the edges between two different patches and it returns the average value of
    !> the two contravariant components at the same point.
    !>
    !> This subroutine is very important to reduce the instability at the edges introduced 
    !> by this cubed-sphere grid 
    !!
    !! @param[inout]   fr                       contravariant radial component of the field
    !! @param[inout]   fxi                      contravariant xi component of the field 
    !! @param[inout]   feta                     contravariant eta component of the field 
    !!
    !! Code owners:
    !!   Clara Dehman
    !!------------------------------------------------------------------------------
  subroutine edge_average(imin,fr,fxi,feta)

   implicit none
  
    ! Subroutine arguments -------------------------------------------------------
      integer, intent(in) :: imin 
      real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(inout) :: fr, fxi, feta
    
    ! Local variables ------------------------------------------------------------
      real*8, dimension(0:nr+1,0:nang+1) :: fr1, fr2, fxi1, fxi2, feta1, feta2
      real*8, dimension(0:nr+1) :: frc1, frc2, frc3, fxic1, fxic2, fxic3, fetac1, fetac2, fetac3
      real*8, dimension(2,2) :: jac1, jac2, jac3, jac4  
      integer :: ed1, ed2 
      
     !----------------------------------------------------------------------------------------
     !----------------------------------  Edges --------------------------------------------
     !----------------------------------------------------------------------------------------
      ed1 = 1 
      ed2 = 2
      ! -------- edge: patch I - II ------------------------- 
      fr1(imin:nr,2:nang-1) = fr(imin:nr,nang,2:nang-1,1)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,nang,2:nang-1,1)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,nang,2:nang-1,1)
      fr2(imin:nr,2:nang-1) = fr(imin:nr,1,2:nang-1,2)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,1,2:nang-1,2)
      feta2(imin:nr,2:nang-1) = feta(imin:nr,1,2:nang-1,2)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,nang,2:nang-1,1) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,nang,2:nang-1,1) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,nang,2:nang-1,1) = feta1(imin:nr,2:nang-1)
      fr(imin:nr,1,2:nang-1,2) = fr2(imin:nr,2:nang-1)
      fxi(imin:nr,1,2:nang-1,2) = fxi2(imin:nr,2:nang-1)
      feta(imin:nr,1,2:nang-1,2) = feta2(imin:nr,2:nang-1)

      ! ---------- edge: patch II - III ---------------------
      fr1(imin:nr,2:nang-1) = fr(imin:nr,nang,2:nang-1,2)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,nang,2:nang-1,2)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,nang,2:nang-1,2)
      fr2(imin:nr,2:nang-1) = fr(imin:nr,1,2:nang-1,3)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,1,2:nang-1,3)
      feta2(imin:nr,2:nang-1) = feta(imin:nr,1,2:nang-1,3)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,nang,2:nang-1,2) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,nang,2:nang-1,2) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,nang,2:nang-1,2) = feta1(imin:nr,2:nang-1)
      fr(imin:nr,1,2:nang-1,3) = fr2(imin:nr,2:nang-1)
      fxi(imin:nr,1,2:nang-1,3) = fxi2(imin:nr,2:nang-1)  
      feta(imin:nr,1,2:nang-1,3) = feta2(imin:nr,2:nang-1)

      ! -------- edge: patch III - IV -----------------------
      fr1(imin:nr,2:nang-1) = fr(imin:nr,nang,2:nang-1,3)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,nang,2:nang-1,3)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,nang,2:nang-1,3)
      fr2(imin:nr,2:nang-1) = fr(imin:nr,1,2:nang-1,4)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,1,2:nang-1,4)
      feta2(imin:nr,2:nang-1) = feta(imin:nr,1,2:nang-1,4)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,nang,2:nang-1,3) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,nang,2:nang-1,3) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,nang,2:nang-1,3) = feta1(imin:nr,2:nang-1)
      fr(imin:nr,1,2:nang-1,4) = fr2(imin:nr,2:nang-1)
      fxi(imin:nr,1,2:nang-1,4) = fxi2(imin:nr,2:nang-1)   
      feta(imin:nr,1,2:nang-1,4)= feta2(imin:nr,2:nang-1)

      ! -------- edge: patch IV - I -------------------------
      fr1(imin:nr,2:nang-1) = fr(imin:nr,nang,2:nang-1,4)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,nang,2:nang-1,4)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,nang,2:nang-1,4)
      fr2(imin:nr,2:nang-1) = fr(imin:nr,1,2:nang-1,1)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,1,2:nang-1,1)
      feta2(imin:nr,2:nang-1) = feta(imin:nr,1,2:nang-1,1)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,nang,2:nang-1,4) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,nang,2:nang-1,4) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,nang,2:nang-1,4) = feta1(imin:nr,2:nang-1)
      fr(imin:nr,1,2:nang-1,1) = fr2(imin:nr,2:nang-1)
      fxi(imin:nr,1,2:nang-1,1) = fxi2(imin:nr,2:nang-1)
      feta(imin:nr,1,2:nang-1,1) = feta2(imin:nr,2:nang-1)
      
      ! -------- edge: patch I - V --------------------------
      ed1 = 3
      ed2 = 4
      fr1(imin:nr,2:nang-1) = fr(imin:nr,2:nang-1,nang,1)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,2:nang-1,nang,1)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,2:nang-1,nang,1)
      fr2(imin:nr,2:nang-1) = fr(imin:nr,2:nang-1,1,5)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,2:nang-1,1,5) 
      feta2(imin:nr,2:nang-1) = feta(imin:nr,2:nang-1,1,5)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,2:nang-1,nang,1) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,2:nang-1,nang,1) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,2:nang-1,nang,1) = feta1(imin:nr,2:nang-1)
      fr(imin:nr,2:nang-1,1,5) = fr2(imin:nr,2:nang-1)
      fxi(imin:nr,2:nang-1,1,5) = fxi2(imin:nr,2:nang-1)
      feta(imin:nr,2:nang-1,1,5) = feta2(imin:nr,2:nang-1)
      
      ! -------- edge: patch II - V -------------------------
      ed1 = 5
      ed2 = 6
      fr1(imin:nr,2:nang-1) = fr(imin:nr,2:nang-1,nang,2)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,2:nang-1,nang,2)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,2:nang-1,nang,2)
      fr2(imin:nr,2:nang-1) = fr(imin:nr,nang,2:nang-1,5)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,nang,2:nang-1,5)  
      feta2(imin:nr,2:nang-1) = feta(imin:nr,nang,2:nang-1,5)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,2:nang-1,nang,2) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,2:nang-1,nang,2) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,2:nang-1,nang,2) = feta1(imin:nr,2:nang-1)
      fr(imin:nr,nang,2:nang-1,5) = fr2(imin:nr,2:nang-1) 
      fxi(imin:nr,nang,2:nang-1,5) = fxi2(imin:nr,2:nang-1)
      feta(imin:nr,nang,2:nang-1,5) = feta2(imin:nr,2:nang-1)

      ! -------- edge: patch III - V ------------------------
      ed1 = 7
      ed2 = 8
      fr1(imin:nr,2:nang-1) = fr(imin:nr,2:nang-1,nang,3)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,2:nang-1,nang,3)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,2:nang-1,nang,3)
      fr2(imin:nr,2:nang-1) = fr(imin:nr,nang-1:2:-1,nang,5)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,nang-1:2:-1,nang,5)
      feta2(imin:nr,2:nang-1) = feta(imin:nr,nang-1:2:-1,nang,5)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,2:nang-1,nang,3) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,2:nang-1,nang,3) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,2:nang-1,nang,3) = feta1(imin:nr,2:nang-1)
      fr(imin:nr,nang-1:2:-1,nang,5) = fr2(imin:nr,2:nang-1)
      fxi(imin:nr,nang-1:2:-1,nang,5) = fxi2(imin:nr,2:nang-1) 
      feta(imin:nr,nang-1:2:-1,nang,5) = feta2(imin:nr,2:nang-1)
      
      ! -------- edge: patch IV - V -------------------------
      ed1 = 9
      ed2 = 10
      fr1(imin:nr,2:nang-1) = fr(imin:nr,2:nang-1,nang,4)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,2:nang-1,nang,4)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,2:nang-1,nang,4)
      fr2(imin:nr,2:nang-1) = fr(imin:nr,1,nang-1:2:-1,5)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,1,nang-1:2:-1,5)
      feta2(imin:nr,2:nang-1) = feta(imin:nr,1,nang-1:2:-1,5)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,2:nang-1,nang,4) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,2:nang-1,nang,4) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,2:nang-1,nang,4) = feta1(imin:nr,2:nang-1)
      fr(imin:nr,1,nang-1:2:-1,5) = fr2(imin:nr,2:nang-1)
      fxi(imin:nr,1,nang-1:2:-1,5) = fxi2(imin:nr,2:nang-1)  
      feta(imin:nr,1,nang-1:2:-1,5) = feta2(imin:nr,2:nang-1)

      ! -------- edge: patch I - VI -------------------------
      ed1 = 11
      ed2 = 12
      fr1(imin:nr,2:nang-1) = fr(imin:nr,2:nang-1,1,1)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,2:nang-1,1,1)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,2:nang-1,1,1)
      fr2(imin:nr,2:nang-1) = fr(imin:nr,2:nang-1,nang,6)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,2:nang-1,nang,6)
      feta2(imin:nr,2:nang-1) = feta(imin:nr,2:nang-1,nang,6)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,2:nang-1,1,1) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,2:nang-1,1,1) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,2:nang-1,1,1) = feta1(imin:nr,2:nang-1)
      fr(imin:nr,2:nang-1,nang,6) = fr2(imin:nr,2:nang-1)
      fxi(imin:nr,2:nang-1,nang,6) = fxi2(imin:nr,2:nang-1)  
      feta(imin:nr,2:nang-1,nang,6) = feta2(imin:nr,2:nang-1)
      
      ! -------- edge: patch II - VI ------------------------
      ed1 = 13
      ed2 = 14
      fr1(imin:nr,2:nang-1) = fr(imin:nr,2:nang-1,1,2)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,2:nang-1,1,2)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,2:nang-1,1,2)
      fr2(imin:nr,2:nang-1) = fr(imin:nr,nang,nang-1:2:-1,6)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,nang,nang-1:2:-1,6)
      feta2(imin:nr,2:nang-1) = feta(imin:nr,nang,nang-1:2:-1,6)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,2:nang-1,1,2) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,2:nang-1,1,2) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,2:nang-1,1,2) = feta1(imin:nr,2:nang-1)
      fr(imin:nr,nang,nang-1:2:-1,6) = fr2(imin:nr,2:nang-1)
      fxi(imin:nr,nang,nang-1:2:-1,6) = fxi2(imin:nr,2:nang-1)  
      feta(imin:nr,nang,nang-1:2:-1,6) = feta2(imin:nr,2:nang-1)
      
      ! -------- edge: patch III - VI -----------------------
      ed1 = 15
      ed2 = 16
      fr1(imin:nr,2:nang-1) = fr(imin:nr,2:nang-1,1,3)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,2:nang-1,1,3)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,2:nang-1,1,3)
      fr2(imin:nr,2:nang-1) = fr(imin:nr,nang-1:2:-1,1,6)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,nang-1:2:-1,1,6)
      feta2(imin:nr,2:nang-1) = feta(imin:nr,nang-1:2:-1,1,6)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,2:nang-1,1,3) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,2:nang-1,1,3) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,2:nang-1,1,3) = feta1(imin:nr,2:nang-1)
      fr(imin:nr,nang-1:2:-1,1,6) = fr2(imin:nr,2:nang-1)
      fxi(imin:nr,nang-1:2:-1,1,6) = fxi2(imin:nr,2:nang-1)  
      feta(imin:nr,nang-1:2:-1,1,6) = feta2(imin:nr,2:nang-1)
    
      ! -------- edge: patch IV - VI ------------------------
      ed1 = 17
      ed2 = 18
      fr1(imin:nr,2:nang-1) = fr(imin:nr,2:nang-1,1,4)
      fxi1(imin:nr,2:nang-1) = fxi(imin:nr,2:nang-1,1,4)
      feta1(imin:nr,2:nang-1) = feta(imin:nr,2:nang-1,1,4) 
      fr2(imin:nr,2:nang-1) = fr(imin:nr,1,2:nang-1,6)
      fxi2(imin:nr,2:nang-1) = fxi(imin:nr,1,2:nang-1,6)
      feta2(imin:nr,2:nang-1) = feta(imin:nr,1,2:nang-1,6)
      call average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)
      fr(imin:nr,2:nang-1,1,4) = fr1(imin:nr,2:nang-1)
      fxi(imin:nr,2:nang-1,1,4) = fxi1(imin:nr,2:nang-1)
      feta(imin:nr,2:nang-1,1,4)  = feta1(imin:nr,2:nang-1)
      fr(imin:nr,1,2:nang-1,6) = fr2(imin:nr,2:nang-1)
      fxi(imin:nr,1,2:nang-1,6) = fxi2(imin:nr,2:nang-1)  
      feta(imin:nr,1,2:nang-1,6) = feta2(imin:nr,2:nang-1)


     !----------------------------------------------------------------------------------------
     !----------------------------------  Corners --------------------------------------------
     !----------------------------------------------------------------------------------------

      ! Corner: patch I - II - V 
      jac1 = jac_avg(1,nang,:,:) ! II -> I
      jac2 = jac_avg(2,nang,:,:) ! I -> II
      jac3 = jac_avg(3,nang,:,:) ! V -> I
      jac4 = jac_avg(4,nang,:,:) ! I -> V
      frc1(imin:nr) = fr(imin:nr,nang,nang,1)
      fxic1(imin:nr) = fxi(imin:nr,nang,nang,1)
      fetac1(imin:nr) = feta(imin:nr,nang,nang,1)
      frc2(imin:nr) = fr(imin:nr,1,nang,2)
      fxic2(imin:nr) = fxi(imin:nr,1,nang,2)
      fetac2(imin:nr) = feta(imin:nr,1,nang,2)
      frc3(imin:nr) = fr(imin:nr,nang,1,5)
      fxic3(imin:nr) = fxi(imin:nr,nang,1,5)
      fetac3(imin:nr) = feta(imin:nr,nang,1,5)
      call average_three(jac1,jac2,jac3,jac4,frc1,fxic1,fetac1,frc2,fxic2,fetac2,frc3,fxic3,fetac3)
      fr(imin:nr,nang,nang,1) = frc1(imin:nr)
      fxi(imin:nr,nang,nang,1) = fxic1(imin:nr)
      feta(imin:nr,nang,nang,1)  = fetac1(imin:nr)
      fr(imin:nr,1,nang,2) = frc2(imin:nr)
      fxi(imin:nr,1,nang,2) = fxic2(imin:nr)  
      feta(imin:nr,1,nang,2) = fetac2(imin:nr)
      fr(imin:nr,nang,1,5) = frc3(imin:nr)
      fxi(imin:nr,nang,1,5) = fxic3(imin:nr)
      feta(imin:nr,nang,1,5)  = fetac3(imin:nr)
     
      ! Corner: patch I - IV - V 
      jac1 = jac_avg(2,nang,:,:) ! IV -> I
      jac2 = jac_avg(1,nang,:,:) ! I -> IV
      jac3 = jac_avg(3,1,:,:)    ! V -> I
      jac4 = jac_avg(4,1,:,:)    ! I -> V
      frc1(imin:nr) = fr(imin:nr,1,nang,1)
      fxic1(imin:nr) = fxi(imin:nr,1,nang,1)
      fetac1(imin:nr) = feta(imin:nr,1,nang,1)
      frc2(imin:nr) = fr(imin:nr,nang,nang,4)
      fxic2(imin:nr) = fxi(imin:nr,nang,nang,4)
      fetac2(imin:nr) = feta(imin:nr,nang,nang,4)
      frc3(imin:nr) = fr(imin:nr,1,1,5)
      fxic3(imin:nr) = fxi(imin:nr,1,1,5)
      fetac3(imin:nr) = feta(imin:nr,1,1,5)
      call average_three(jac1,jac2,jac3,jac4,frc1,fxic1,fetac1,frc2,fxic2,fetac2,frc3,fxic3,fetac3)
      fr(imin:nr,1,nang,1) = frc1(imin:nr)
      fxi(imin:nr,1,nang,1) = fxic1(imin:nr)
      feta(imin:nr,1,nang,1)  = fetac1(imin:nr)
      fr(imin:nr,nang,nang,4) = frc2(imin:nr)
      fxi(imin:nr,nang,nang,4) = fxic2(imin:nr)  
      feta(imin:nr,nang,nang,4) = fetac2(imin:nr)
      fr(imin:nr,1,1,5) = frc3(imin:nr)
      fxi(imin:nr,1,1,5) = fxic3(imin:nr)
      feta(imin:nr,1,1,5)  = fetac3(imin:nr)

      ! Corner: patch II - III - V
      jac1 = jac_avg(1,nang,:,:) ! III -> II
      jac2 = jac_avg(2,nang,:,:) ! II -> III
      jac3 = jac_avg(5,nang,:,:) ! V -> II
      jac4 = jac_avg(6,nang,:,:) ! II -> V
      frc1(imin:nr) = fr(imin:nr,nang,nang,2)
      fxic1(imin:nr) = fxi(imin:nr,nang,nang,2)
      fetac1(imin:nr) = feta(imin:nr,nang,nang,2)
      frc2(imin:nr) = fr(imin:nr,1,nang,3)
      fxic2(imin:nr) = fxi(imin:nr,1,nang,3)
      fetac2(imin:nr) = feta(imin:nr,1,nang,3)
      frc3(imin:nr) = fr(imin:nr,nang,nang,5)
      fxic3(imin:nr) = fxi(imin:nr,nang,nang,5)
      fetac3(imin:nr) = feta(imin:nr,nang,nang,5)
      call average_three(jac1,jac2,jac3,jac4,frc1,fxic1,fetac1,frc2,fxic2,fetac2,frc3,fxic3,fetac3)
      fr(imin:nr,nang,nang,2) = frc1(imin:nr)
      fxi(imin:nr,nang,nang,2) = fxic1(imin:nr)
      feta(imin:nr,nang,nang,2) = fetac1(imin:nr)
      fr(imin:nr,1,nang,3) = frc2(imin:nr)
      fxi(imin:nr,1,nang,3) = fxic2(imin:nr)  
      feta(imin:nr,1,nang,3) = fetac2(imin:nr)
      fr(imin:nr,nang,nang,5) = frc3(imin:nr)
      fxi(imin:nr,nang,nang,5) = fxic3(imin:nr)
      feta(imin:nr,nang,nang,5) = fetac3(imin:nr)
 
      ! Corner: patch III - IV - V
      jac1 = jac_avg(1,nang,:,:) ! IV -> III
      jac2 = jac_avg(2,nang,:,:) ! III -> IV
      jac3 = jac_avg(7,nang,:,:) ! V -> III
      jac4 = jac_avg(8,nang,:,:) ! III -> V
      frc1(imin:nr) = fr(imin:nr,nang,nang,3)
      fxic1(imin:nr) = fxi(imin:nr,nang,nang,3)
      fetac1(imin:nr) = feta(imin:nr,nang,nang,3)
      frc2(imin:nr) = fr(imin:nr,1,nang,4)
      fxic2(imin:nr) = fxi(imin:nr,1,nang,4)
      fetac2(imin:nr) = feta(imin:nr,1,nang,4)
      frc3(imin:nr) = fr(imin:nr,1,nang,5)
      fxic3(imin:nr) = fxi(imin:nr,1,nang,5)
      fetac3(imin:nr) = feta(imin:nr,1,nang,5)
      call average_three(jac1,jac2,jac3,jac4,frc1,fxic1,fetac1,frc2,fxic2,fetac2,frc3,fxic3,fetac3)
      fr(imin:nr,nang,nang,3) = frc1(imin:nr)
      fxi(imin:nr,nang,nang,3) = fxic1(imin:nr)
      feta(imin:nr,nang,nang,3) = fetac1(imin:nr)
      fr(imin:nr,1,nang,4) = frc2(imin:nr)
      fxi(imin:nr,1,nang,4) = fxic2(imin:nr)  
      feta(imin:nr,1,nang,4) = fetac2(imin:nr)
      fr(imin:nr,1,nang,5) = frc3(imin:nr)
      fxi(imin:nr,1,nang,5) = fxic3(imin:nr)
      feta(imin:nr,1,nang,5) = fetac3(imin:nr)

      ! Corner: patch I - II - VI 
      jac1 = jac_avg(1,1,:,:) ! II -> I
      jac2 = jac_avg(2,1,:,:) ! I -> II
      jac3 = jac_avg(11,nang,:,:) ! VI -> I
      jac4 = jac_avg(12,nang,:,:) ! I -> VI
      frc1(imin:nr) = fr(imin:nr,nang,1,1)
      fxic1(imin:nr) = fxi(imin:nr,nang,1,1)
      fetac1(imin:nr) = feta(imin:nr,nang,1,1)
      frc2(imin:nr) = fr(imin:nr,1,1,2)
      fxic2(imin:nr) = fxi(imin:nr,1,1,2)
      fetac2(imin:nr) = feta(imin:nr,1,1,2)
      frc3(imin:nr) = fr(imin:nr,nang,nang,6)
      fxic3(imin:nr) = fxi(imin:nr,nang,nang,6)
      fetac3(imin:nr) = feta(imin:nr,nang,nang,6)
      call average_three(jac1,jac2,jac3,jac4,frc1,fxic1,fetac1,frc2,fxic2,fetac2,frc3,fxic3,fetac3)
      fr(imin:nr,nang,1,1) = frc1(imin:nr)
      fxi(imin:nr,nang,1,1) = fxic1(imin:nr)
      feta(imin:nr,nang,1,1) = fetac1(imin:nr)
      fr(imin:nr,1,1,2) = frc2(imin:nr)
      fxi(imin:nr,1,1,2) = fxic2(imin:nr)  
      feta(imin:nr,1,1,2) = fetac2(imin:nr)
      fr(imin:nr,nang,nang,6) = frc3(imin:nr)
      fxi(imin:nr,nang,nang,6) = fxic3(imin:nr)
      feta(imin:nr,nang,nang,6) = fetac3(imin:nr)

      ! Corner: patch II - III - VI
      jac1 = jac_avg(1,1,:,:) ! III -> II
      jac2 = jac_avg(2,1,:,:) ! II -> III
      jac3 = jac_avg(13,nang,:,:) ! VI -> II
      jac4 = jac_avg(14,nang,:,:) ! II -> VI
      frc1(imin:nr) = fr(imin:nr,nang,1,2)
      fxic1(imin:nr) = fxi(imin:nr,nang,1,2)
      fetac1(imin:nr) = feta(imin:nr,nang,1,2)
      frc2(imin:nr) = fr(imin:nr,1,1,3)
      fxic2(imin:nr) = fxi(imin:nr,1,1,3)
      fetac2(imin:nr) = feta(imin:nr,1,1,3)
      frc3(imin:nr) = fr(imin:nr,nang,1,6)
      fxic3(imin:nr) = fxi(imin:nr,nang,1,6)
      fetac3(imin:nr) = feta(imin:nr,nang,1,6)
      call average_three(jac1,jac2,jac3,jac4,frc1,fxic1,fetac1,frc2,fxic2,fetac2,frc3,fxic3,fetac3)
      fr(imin:nr,nang,1,2) = frc1(imin:nr)
      fxi(imin:nr,nang,1,2)= fxic1(imin:nr)
      feta(imin:nr,nang,1,2) = fetac1(imin:nr)
      fr(imin:nr,1,1,3) = frc2(imin:nr)
      fxi(imin:nr,1,1,3) = fxic2(imin:nr)  
      feta(imin:nr,1,1,3)= fetac2(imin:nr)
      fr(imin:nr,nang,1,6) = frc3(imin:nr)
      fxi(imin:nr,nang,1,6) = fxic3(imin:nr)
      feta(imin:nr,nang,1,6) = fetac3(imin:nr)

      ! Corner: patch I - IV - VI 
      jac1 = jac_avg(2,1,:,:) ! IV -> I
      jac2 = jac_avg(1,1,:,:) ! I -> IV
      jac3 = jac_avg(11,1,:,:) ! VI -> I
      jac4 = jac_avg(12,1,:,:) ! I -> VI
      frc1(imin:nr) = fr(imin:nr,1,1,1)
      fxic1(imin:nr) = fxi(imin:nr,1,1,1)
      fetac1(imin:nr) = feta(imin:nr,1,1,1)
      frc2(imin:nr) = fr(imin:nr,nang,1,4)
      fxic2(imin:nr) = fxi(imin:nr,nang,1,4)
      fetac2(imin:nr) = feta(imin:nr,nang,1,4)
      frc3(imin:nr) = fr(imin:nr,1,nang,6)
      fxic3(imin:nr) = fxi(imin:nr,1,nang,6)
      fetac3(imin:nr) = feta(imin:nr,1,nang,6)
      call average_three(jac1,jac2,jac3,jac4,frc1,fxic1,fetac1,frc2,fxic2,fetac2,frc3,fxic3,fetac3)
      fr(imin:nr,1,1,1) = frc1(imin:nr)
      fxi(imin:nr,1,1,1) = fxic1(imin:nr)
      feta(imin:nr,1,1,1) = fetac1(imin:nr)
      fr(imin:nr,nang,1,4) = frc2(imin:nr)
      fxi(imin:nr,nang,1,4) = fxic2(imin:nr)  
      feta(imin:nr,nang,1,4) = fetac2(imin:nr)
      fr(imin:nr,1,nang,6) = frc3(imin:nr)
      fxi(imin:nr,1,nang,6) = fxic3(imin:nr)
      feta(imin:nr,1,nang,6) = fetac3(imin:nr)

     ! Corner: patch III - IV - VI
      jac1 = jac_avg(1,1,:,:) ! IV -> III
      jac2 = jac_avg(2,1,:,:) ! III -> IV
      jac3 = jac_avg(15,nang,:,:) ! VI -> III
      jac4 = jac_avg(16,nang,:,:) ! III -> VI
      frc1(imin:nr) = fr(imin:nr,nang,1,3)
      fxic1(imin:nr) = fxi(imin:nr,nang,1,3)
      fetac1(imin:nr) = feta(imin:nr,nang,1,3)
      frc2(imin:nr) = fr(imin:nr,1,1,4)
      fxic2(imin:nr) = fxi(imin:nr,1,1,4)
      fetac2(imin:nr) = feta(imin:nr,1,1,4)
      frc3(imin:nr) = fr(imin:nr,1,1,6)
      fxic3(imin:nr) = fxi(imin:nr,1,1,6)
      fetac3(imin:nr) = feta(imin:nr,1,1,6)
      call average_three(jac1,jac2,jac3,jac4,frc1,fxic1,fetac1,frc2,fxic2,fetac2,frc3,fxic3,fetac3)
      fr(imin:nr,nang,1,3) = frc1(imin:nr)
      fxi(imin:nr,nang,1,3) = fxic1(imin:nr)
      feta(imin:nr,nang,1,3) = fetac1(imin:nr)
      fr(imin:nr,1,1,4) = frc2(imin:nr)
      fxi(imin:nr,1,1,4) = fxic2(imin:nr)  
      feta(imin:nr,1,1,4) = fetac2(imin:nr)
      fr(imin:nr,1,1,6) = frc3(imin:nr)
      fxi(imin:nr,1,1,6) = fxic3(imin:nr)
      feta(imin:nr,1,1,6) = fetac3(imin:nr)
  
     end subroutine edge_average

!!------------------------------------------------------------------------------
    !> @brief This subroutine is responsible for doing the average between the
    !> two contravariant components of a vector at the edges belonging to 
    !>  two different patches and it returns the value of the two contravariant 
    !> components at the same point after being averaged
    !>
    !! @param[in]      ed1, ed2              specify the edge between two patches and the
    !!                                       tranformation from one patch to the other or viceversa 
    !! @param[inout]   fr1, fr2              contravariant radial components of the field
    !! @param[inout]   fxi1, fxi2            contravariant xi components of the field 
    !! @param[inout]   feta1, feta2          contravariant eta components of the field 
    !!
    !! Code owners:
    !!   Clara Dehman
    !!------------------------------------------------------------------------------
     subroutine average_two(ed1,ed2,fr1,fxi1,feta1,fr2,fxi2,feta2)

      implicit none
      integer, intent(in) :: ed1, ed2
      real*8, dimension (0:nr+1,0:nang+1), intent(inout) :: fr1,fxi1,feta1,fr2,fxi2,feta2
    
    ! Local variables ------------------------------------------------------------
      real*8, dimension(0:nr+1,0:nang+1) :: fxi_adj, feta_adj
      real*8, dimension(0:nr+1,0:nang+1) :: fr_avg, fxi_avg, feta_avg
      real*8, dimension(0:nang+1,2,2) :: jac1, jac2
      integer :: j
     
      fxi_adj = 0.d0
      feta_adj = 0.d0
      fr_avg = 0.d0
      fxi_avg = 0.d0
      feta_avg = 0.d0

      jac1 = jac_avg(ed1,:,:,:)
      jac2 = jac_avg(ed2,:,:,:)
    
      do j = 2, nang-1
      ! -------- radial component of the field -------------
      fr_avg(:,j) = (fr1(:,j)+fr2(:,j))/2.d0
      fr1(:,j) = fr_avg(:,j)
      fr2(:,j) = fr_avg(:,j)
      ! -------- angular components of the field ------------
      ! changing the field coordinate from one patch to the other 
      fxi_adj(:,j) = jac1(j,1,1)*fxi2(:,j) + jac1(j,1,2)*feta2(:,j)
      feta_adj(:,j) = jac1(j,2,1)*fxi2(:,j) + jac1(j,2,2)*feta2(:,j)
     ! average value in the coordinate system of one of the patches 
      fxi_avg(:,j) = (fxi1(:,j)+fxi_adj(:,j))/2.d0
      feta_avg(:,j) = (feta1(:,j)+feta_adj(:,j))/2.d0
      ! fixing the value of the field components to the average value for the non-tranformed patch
      fxi1(:,j) = fxi_avg(:,j)
      feta1(:,j) = feta_avg(:,j)
      ! fixing the value of the field components to the average value for the tranformed patch
      fxi2(:,j) = jac2(j,1,1)*fxi_avg(:,j) + jac2(j,1,2)*feta_avg(:,j)
      feta2(:,j) = jac2(j,2,1)*fxi_avg(:,j) + jac2(j,2,2)*feta_avg(:,j)
      end do 

      end subroutine average_two


!!------------------------------------------------------------------------------
    !> @brief This subroutine is responsible for doing the average at the corner 
    !> between the three contravariant components of a vector belonging to  
    !> three different patches and it returns the value of the three contravariant 
    !> components at the same point after being averaged
    !>
    !! @param[in]      jac1, jac2            Jococians responsible for the tranformation
    !!                                       from one patch to the other or viceversa 
    !! @param[in]      jac3, jac3            Jococians responsible for the tranformation
    !!                                       from one patch to the other or viceversa
    !! @param[inout]   fr1, fr2, fr3         contravariant radial components of the field
    !! @param[inout]   fxi1, fxi2, fxi3      contravariant xi components of the field 
    !! @param[inout]   feta1, feta2, feta3   contravariant eta components of the field 
    !!
    !! Code owners:
    !!   Clara Dehman
    !!------------------------------------------------------------------------------
     subroutine average_three(jac1,jac2,jac3,jac4,fr1,fxi1,feta1,fr2,fxi2,feta2,fr3,fxi3,feta3)

      implicit none
      real*8, dimension(2,2), intent(in) :: jac1, jac2, jac3, jac4
      real*8, dimension (0:nr+1), intent(inout) :: fr1,fxi1,feta1
      real*8, dimension (0:nr+1), intent(inout) :: fr2,fxi2,feta2
      real*8, dimension (0:nr+1), intent(inout) :: fr3,fxi3,feta3
    ! Local variables ------------------------------------------------------------
      real*8, dimension(0:nr+1) :: fxi_adj1, feta_adj1, fxi_adj2, feta_adj2
      real*8, dimension(0:nr+1) :: fr_avg, fxi_avg, feta_avg
     
      fxi_adj1 = 0.d0
      feta_adj1 = 0.d0
      fxi_adj2 = 0.d0
      feta_adj2 = 0.d0
      fr_avg = 0.d0
      fxi_avg = 0.d0
      feta_avg = 0.d0
    
      ! -------- radial component of the field -------------
      fr_avg = (fr1+fr2+fr3)/3.d0
      fr1 = fr_avg
      fr2 = fr_avg
      fr3 = fr_avg
      ! -------- angular components of the field ------------
      ! changing the field coordinate from one patch to the patch where we will apply the average  
      fxi_adj1(:) = jac1(1,1)*fxi2(:) + jac1(1,2)*feta2(:)
      feta_adj1(:) = jac1(2,1)*fxi2(:) + jac1(2,2)*feta2(:)
      ! changing the field coordinate from the other patch to the patch where we will apply the average  
      fxi_adj2(:) = jac3(1,1)*fxi3(:) + jac3(1,2)*feta3(:)
      feta_adj2(:) = jac3(2,1)*fxi3(:) + jac3(2,2)*feta3(:)
     ! average value in the coordinate system where we are applying the average
      fxi_avg(:) = (fxi1(:)+fxi_adj1(:)+fxi_adj2(:))/3.d0
      feta_avg(:) = (feta1(:)+feta_adj1(:)+feta_adj2(:))/3.d0
      ! fixing the value of the field components to the average value for the non-tranformed patch
      fxi1(:) = fxi_avg(:)
      feta1(:) = feta_avg(:)
      ! fixing the value of the field components to the average value for the first tranformed patch
      fxi2(:) = jac2(1,1)*fxi_avg(:) + jac2(1,2)*feta_avg(:)
      feta2(:) = jac2(2,1)*fxi_avg(:) + jac2(2,2)*feta_avg(:)
      ! fixing the value of the field components to the average value for the second tranformed patch
      fxi3(:) = jac4(1,1)*fxi_avg(:) + jac4(1,2)*feta_avg(:)
      feta3(:) = jac4(2,1)*fxi_avg(:) + jac4(2,2)*feta_avg(:)
      
      end subroutine average_three

!!------------------------------------------------------------------------------
!!------------------------------------------------------------------------------
    !> @brief This subroutines provide the contravariant component of the curl of a vector
    !> using a finite difference method (eq. 24 of Ronchi's formalism )
    !! @param[in]   fr                       contravariant radial component
    !! @param[in]   fxi                      contravariant xi component
    !! @param[in]   feta                     contravariant eta component
    !! @param[out]  curlfr                   contravariant radial component of curlf
    !! @param[out]  curlfxi                  contravariant xi component of curlf
    !! @param[out]  curlfeta                 contravariant eta component of curlf
    !!
    !!
    !! Code owners:
    !!   Clara Dehman
    !!------------------------------------------------------------------------------
  subroutine curl_fndiff(fr,fxi,feta,curlfr,curlfxi,curlfeta,imin)
  
      implicit none
  
      ! Subroutine arguments -------------------------------------------------------
      integer, intent(in) :: imin 
      real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: fr,fxi,feta
      real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(out) :: curlfr, curlfxi, curlfeta
  
      ! Local variables ------------------------------------------------------------
      ! Auxiliary index for loops
      integer :: i,j,k,p
  
      ! ----------------------------------------------------------------------------
      curlfr = 0.d0
      curlfxi = 0.d0
      curlfeta = 0.d0

      do p=1,6
       do i = imin,nr
        do j=1,nang
          do k=1,nang
          curlfr(i,j,k,p)= delta(j,k)**0.5d0/r(i)*(X(j)*Y(k)/(C(j)*D(k))* & 
          & (1/D(k)*(feta(i,j,k+1,p) - feta(i,j,k-1,p))/(eta(k+1) - eta(k-1)) -  &
          & 1/C(j)*(fxi(i,j+1,k,p) - fxi(i,j-1,k,p))/(xi(j+1) - xi(j-1))) - &
          & 1/D(k)*(fxi(i,j,k+1,p) - fxi(i,j,k-1,p))/(eta(k+1) - eta(k-1)) + & 
          & 1/C(j)*(feta(i,j+1,k,p) - feta(i,j-1,k,p))/(xi(j+1) - xi(j-1)))
         enddo
        enddo
       enddo
    
        do i=imin,nr
         do j=1,nang
          do k=1,nang
            curlfxi(i,j,k,p)= 1/r(i)*(X(j)*Y(k)/(delta(j,k)**0.5d0)* &
            & (fxi(i,j,k,p)+ r(i)*(fxi(i+1,j,k,p)-fxi(i-1,j,k,p))/(r(i+1)-r(i-1))) & 
            & - C(j)*D(k)/(delta(j,k)**0.5d0)* & 
            & (feta(i,j,k,p)+ r(i)*(feta(i+1,j,k,p)-feta(i-1,j,k,p))/(r(i+1)-r(i-1))) & 
            & + delta(j,k)**0.5d0/D(k)*(fr(i,j,k+1,p) - fr(i,j,k-1,p))/(eta(k+1) - eta(k-1)))
          enddo
         enddo 
        enddo
  
        do i=imin,nr
          do j=1,nang
           do k=1,nang
            curlfeta(i,j,k,p)= 1/r(i)*(C(j)*D(k)/(delta(j,k)**0.5d0)* &
            & (fxi(i,j,k,p)+ r(i)*(fxi(i+1,j,k,p)-fxi(i-1,j,k,p))/(r(i+1)-r(i-1))) & 
            & - X(j)*Y(k)/(delta(j,k)**0.5d0)* & 
            & (feta(i,j,k,p)+ r(i)*(feta(i+1,j,k,p)-feta(i-1,j,k,p))/(r(i+1)-r(i-1))) & 
            & - delta(j,k)**0.5d0/C(j)*(fr(i,j+1,k,p) - fr(i,j-1,k,p))/(xi(j+1) - xi(j-1)))
          enddo
         enddo
        enddo
      enddo
  end subroutine curl_fndiff
 
!!--------------------------------------------------------------------------
    !> @brief Subroutine Divergence of the magnetic field 
    !! In this subroutine we calculate nabla.B in the cubed-sphere formalism
    !!
    !! Code owners:
    !!    Clara Dehman
    !!--------------------------------------------------------------------------
    subroutine divB_fndiff(fr,fxi,feta,divB)

    implicit none
    integer :: i, j, k, p
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: fr, feta, fxi
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(out) :: divB

  ! initialize to zero
    divB = 0d0     
    do p = 1, 6
     do i = 3, nr - 1  
      do j = 1, nang
      do k = 1, nang
      !  divB(i,j,k,p) = delta(j,k)/(r(i)*D(k)*C(j)**2)*((fxi(i,j+1,k,p)-fxi(i,j-1,k,p))/(xi(j+1)-xi(j-1)) & 
      !  & - fxi(i,j,k,p)/(delta(j,k)**(1/2))*(delta(j+1,k)-delta(j-1,k))/(xi(j+1)-xi(j-1)))  & 
      !  & + delta(j,k)/(r(i)*C(j)*D(k)**2)*((feta(i,j,k+1,p)-feta(i,j,k-1,p))/(eta(k+1)-eta(k-1)) & 
      !  & - feta(i,j,k,p)/(delta(j,k)**(1/2))*(delta(j,k+1)-delta(j,k-1))/(eta(k+1)-eta(k-1)))  & 
      !  & + (fr(i+1,j,k,p)-fr(i-1,j,k,p))/(r(i+1)-r(i-1)) + 2.d0/r(i)*fr(i,j,k,p)
 
       divB(i,j,k,p) = delta(j,k)**(3/2)/(r(i)*D(k)*C(j)**2)*(fxi(i,j+1,k,p)/dsqrt(delta(j+1,k)) & 
       & - fxi(i,j-1,k,p)/dsqrt(delta(j-1,k)))/(xi(j+1)-xi(j-1)) & 
       & + delta(j,k)**(3/2)/(r(i)*C(j)*D(k)**2)*(feta(i,j,k+1,p)/dsqrt(delta(j,k+1)) & 
       & - feta(i,j,k-1,p)/dsqrt(delta(j,k-1)))/(eta(k+1)-eta(k-1)) &
       & + 1.d0/(r(i)**2)*(r(i+1)**2*br(i+1,j,k,p)-r(i-1)**2*br(i-1,j,k,p))/(r(i+1)-r(i-1))
      end do 
      end do 
     end do 
    end do 
    end subroutine divB_fndiff

    !!--------------------------------------------------------------------------
    !> @brief Subroutine Divergence of the magnetic field 
    !> using a finite volume method (divergence theorem). 
    !> For the length elements, we are using the contravariant components.
    !!
    !! @param[in]  lr, lxi, leta             contravariant length elements 
    !! @param[in]   fr                       contravariant radial component
    !! @param[in]   fxi                      contravariant xi component
    !! @param[in]   feta                     contravariant eta component
    !! @param[out]  divB                     div(B)
    !!
    !! Code owners:
    !!   Clara Dehman
    !!--------------------------------------------------------------------------
    subroutine divB_fnvol(fr,fxi,feta,divB)

    implicit none
    integer :: i, j, k, p
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: fr, fxi, feta
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(out) :: divB

  ! initialize to zero
    divB = 0d0     
    do p = 1, 6
     do i = 3, nr - 1
      do j = 1, nang
      do k = 1, nang
       divB(i,j,k,p) = 1.d0/(lxi(i,j,k)*leta(i,j,k))*(1.d0/lr(i)*(br(i+1,j,k,p)*lxi(i+1,j,k)*leta(i+1,j,k) & 
      & - br(i-1,j,k,p)*lxi(i-1,j,k)*leta(i-1,j,k))/(r(i+1)-r(i-1)) & 
      & + C(j)/dsqrt(delta(j,k))*(dsqrt(delta(j+1,k))/C(j+1)*bxi(i,j+1,k,p)*leta(i,j+1,k) & 
      & - dsqrt(delta(j-1,k))/C(j-1)*bxi(i,j-1,k,p)*leta(i,j-1,k))/(xi(j+1) - xi(j-1))  & 
      & + D(k)/dsqrt(delta(j,k))*(dsqrt(delta(j,k+1))/D(k+1)*beta(i,j,k+1,p)*lxi(i,j,k+1) &
      & - dsqrt(delta(j,k-1))/D(k-1)*beta(i,j,k-1,p)*lxi(i,j,k-1))/(eta(k+1) - eta(k-1)))
      end do 
      end do 
     end do 
    end do 
    end subroutine divB_fnvol

!!--------------------------------------------------------------------------
    !> @brief Subroutine dot_prod
    !!
    !! In this subroutine we express the dot product in term of the non-orthogonal
    !! metric tensor
    !! We are working in versor basis
    !!
    !! Code owners:
    !!    Clara Dehman
    !!--------------------------------------------------------------------------
    subroutine dot_prod(ar,fr,axi,fxi,aeta,feta,dprod)

    implicit none
    integer :: p
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: ar, aeta, axi
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: fr, feta, fxi
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(out) :: dprod  

  ! initialize to zero
    dprod = 0d0     
  ! Writing the dot product in term of the metric g
    
    !$OMP Parallel do private(p) shared(dprod, g, ar, fr, axi, fxi, aeta, feta)
    do p = 1, 6
      dprod(:, :, :, p) =  g(:, :, :, 1, 1)*ar(:, :, :, p)*fr(:, :, :, p) & 
           & + g(:, :, :, 2, 2)*axi(:, :, :, p)*fxi(:, :, :, p) &
           & + g(:, :, :, 3, 3)*aeta(:, :, :, p)*feta(:, :, :, p)  &
           & + g(:, :, :, 2, 3)*(axi(:, :, :, p)*feta(:, :, :, p) + aeta(:, :, :, p)*fxi(:, :, :, p))
    end do 
    !$OMP end Parallel do

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
    subroutine crossprod_cont(fr,fxi,feta,gr,gxi,geta,Xprodr,Xprodxi,Xprodeta,imin)

    implicit none
    integer :: i, j, p

    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: fr, fxi, feta
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: gr, gxi, geta
    integer, intent(in) :: imin

    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(out) :: Xprodr, Xprodxi, Xprodeta

    real*8, dimension (0:nang+1) :: tmp1, tmp2, XY, CD
    
  ! initialize to zero
    Xprodr = 0d0
    Xprodxi = 0d0
    Xprodeta = 0d0
 
    !$OMP Parallel do private(j, CD, XY, tmp1, tmp2)
    do p = 1, 6  
      do j = 0, nang+1
        CD(:) = C(j)*D(:)
        XY(:) = X(j)*Y(:)
        do i = imin-1, nr+1
          tmp1 = (feta(i, j, :, p)*gr(i, j, :, p) - fr(i, j, :, p)*geta(i, j, :, p))
          tmp2 = (fr(i, j, :, p)*gxi(i, j, :, p) - fxi(i, j, :, p)*gr(i, j, :, p))
          Xprodr(i, j, :, p) = delta(j, :)**0.5d0/CD(:) &
          & * (fxi(i, j, :, p)*geta(i, j, :, p) - feta(i, j, :, p)*gxi(i, j, :, p))
          Xprodxi(i, j, :, p) = (CD(:)*tmp1 + XY(:)*tmp2)/(delta(j, :)**0.5d0)  
          Xprodeta(i, j, :, p) = (XY(:)*tmp1 + CD(:)*tmp2)/(delta(j, :)**0.5d0)
        end do 
      end do 
    end do 
    !$OMP end Parallel do 

  end subroutine crossprod_cont



!!--------------------------------------------------------------------------
    !> @brief Subroutine crossprod_cont_upwind
    !!
    !! In this subroutine we define the contravariant components of the cross product
    !! in term of the contravariant components of the vectors using the versor basis
    !! for our non-orthogonal metric and considering an upwind scheme
    !! Code owners:
    !!    Clara Dehman
    !!--------------------------------------------------------------------------
    !! Note: Xprodr refers to the contravariant component of the cross product 
    !! in r physical direction
    !!--------------------------------------------------------------------------
  subroutine crossprod_cont_upwind(fr,fxi,feta,gr,gxi,geta,Xprodr,Xprodxi,Xprodeta,imin)

    implicit none
    integer :: i, j, k, p

    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: fr, fxi, feta
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: gr, gxi, geta
    integer, intent(in) :: imin
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(out) :: Xprodr, Xprodxi, Xprodeta

    !--------------------- local variables -------------------------------------
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: grwxi, grweta ! gr with fxi and with feta respectively 
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: gxiwr, getawr ! gxi and geta, with fr respectively 
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: gxiweta, getawxi ! gxi with feta and geta with fxi respectively 

    grwxi = gr ! gr with fxi
    grweta = gr ! gr with feta
    gxiwr = gxi ! gxi with fr
    gxiweta = gxi  ! gxi with feta 
    getawr = geta ! geta with fr 
    getawxi = geta ! geta with fxi 
  
    do p = 1, 6  
      do i = imin+1, nr-1   
        do j = 1, nang
        do k = 1, nang 
        ! condition on the sign of fr
          if (fr(i,j,k,p) > 0.d0) then
            ! gxi with fr  
            gxiwr(i,j,k,p) = gxi(i+1,j,k,p)
            ! geta with fr 
            getawr(i,j,k,p) = geta(i+1,j,k,p)
          elseif (fr(i,j,k,p) < 0.d0) then
             ! gxi with fr  
             gxiwr(i,j,k,p) = gxi(i-1,j,k,p)
             ! geta with fr 
             getawr(i,j,k,p) = geta(i-1,j,k,p)
          endif
        ! condition on the sign of fxi 
          if (fxi(i,j,k,p) > 0.d0) then
           ! geta with fxi 
           getawxi(i,j,k,p)=geta(i,j+1,k,p)
           ! gr with fxi
           grwxi(i,j,k,p)=gr(i,j+1,k,p)
          elseif (fxi(i,j,k,p) < 0.d0) then
           ! geta with fxi 
           getawxi(i,j,k,p)=geta(i,j-1,k,p)
           ! gr with fxi 
           grwxi(i,j,k,p)=gr(i,j-1,k,p)
          endif
        ! condition on the sign of feta
          if (feta(i,j,k,p) > 0.d0) then
            ! gxi with feta 
            gxiweta(i,j,k,p)=gxi(i,j,k+1,p)  
            ! gr with feta 
            grweta(i,j,k,p)=gr(i,j,k+1,p)
          elseif (feta(i,j,k,p) < 0.d0) then
            ! gxi with feta 
            gxiweta(i,j,k,p)=gxi(i,j,k-1,p)
            ! gr with feta 
            grweta(i,j,k,p)=gr(i,j,k-1,p) 
          endif
        end do 
      end do 
    end do 
  end do 

  ! initialize to zero
    Xprodr = 0d0
    Xprodxi = 0d0
    Xprodeta = 0d0
    do p = 1, 6  
      do i = imin-1, nr+1 
        do j = 0, nang+1
          Xprodr(i, j, :, p) = delta(j, :)**0.5d0/(C(j)*D(:)) &
          & * (fxi(i, j, :, p)*getawxi(i, j, :, p) - feta(i, j, :, p)*gxiweta(i, j, :, p)) 
          Xprodxi(i, j, :, p) = (C(j)*D(:)*(feta(i, j, :, p)*grweta(i, j, :, p) - fr(i, j, :, p)*getawr(i, j, :, p)) &
          & + X(j)*Y(:)*(fr(i, j, :, p)*gxiwr(i, j, :, p) - fxi(i, j, :, p)*grwxi(i, j, :, p)))/(delta(j, :)**0.5d0)
          Xprodeta(i, j, :, p) = (X(j)*Y(:)*(feta(i, j, :, p)*grweta(i, j, :, p) - fr(i, j, :, p)*getawr(i, j, :, p)) &
          & + C(j)*D(:)*(fr(i, j, :, p)*gxiwr(i, j, :, p) - fxi(i, j, :, p)*grwxi(i, j, :, p)))/(delta(j, :)**0.5d0)
        end do 
      end do 
    end do 

  end subroutine crossprod_cont_upwind


    !!--------------------------------------------------------------------------
    !> @brief Subroutine crossprod_cov
    !! In this subroutine we define the covariant components of the cross product
    !! in term of the contravariant components of the vectors
    !! using the versor basis for our non-orthogonal metric
    !! Code owners:
    !!    Clara Dehman
    !!--------------------------------------------------------------------------
    !! Note: Xprod_r refers to the covariant component of the cross product 
    !! in r physical direction
    !!--------------------------------------------------------------------------
  subroutine crossprod_cov(fr,fxi,feta,gr,gxi,geta,Xprod_r,Xprod_xi,Xprod_eta,imin)

    implicit none
    integer :: i, j, p
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: fr, fxi, feta
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: gr, gxi, geta
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(out) :: Xprod_r, Xprod_xi, Xprod_eta
    integer, intent(in) :: imin
    
   ! initialize to zero
    Xprod_r = 0.d0
    Xprod_xi = 0.d0
    Xprod_eta = 0.d0

    !$OMP Parallel do private(i,j)
 
    do p = 1, 6
      do i = imin-1, nr+1
        do j = 0, nang+1
          Xprod_r(i, j, :, p) = delta(j, :)**(0.5d0)/(C(j)*D(:)) &
           & * (fxi(i, j, :, p)*geta(i, j, :, p) - feta(i, j, :, p)*gxi(i, j, :, p))
          Xprod_xi(i, j, :, p) = delta(j, :)**(0.5d0)/(C(j)*D(:)) &
           & * (feta(i, j, :, p)*gr(i, j, :, p) - fr(i, j, :, p)*geta(i, j, :, p))
          Xprod_eta(i, j, :, p) = delta(j, :)**(0.5d0)/(C(j)*D(:)) &
           & * (fr(i, j, :, p)*gxi(i, j, :, p) - fxi(i, j, :, p)*gr(i, j, :, p))
        end do 
      end do 
    end do

    !$OMP end Parallel do
 
  end subroutine crossprod_cov

  !!--------------------------------------------------------------------------
  !> @brief Subroutine cs_to_cartesian
  !!
  !! In this subroutine we apply the transformations needed to obtain the 
  !! cartesian coordinates using cubed sphere and radial spherical coordinates
  !!
  !! @param[in]   X, Y, delta              cubed sphere coordinates 
  !! @param[in]   r                        radial coordinate
  !! @param[out]  xc, yc, zc               Cartesian coordinates     
  !!
  !! Code owners:
  !!    Clara Dehman
  !!    Daniele Viganò
  !!--------------------------------------------------------------------------
  subroutine cs_to_cartesian
    implicit none

    integer :: j, k, p
   
    do j = 0, nang+1
      do k = 0, nang+1
      ! patch I 
        p = 1
        xc(:,j,k,p) = r(:)/sqrt(delta(j,k))
        yc(:,j,k,p) = r(:)*X(j)/sqrt(delta(j,k))
        zc(:,j,k,p) = r(:)*Y(k)/sqrt(delta(j,k))
      ! patch II
        p = 2
        xc(:,j,k,p) = - r(:)*X(j)/sqrt(delta(j,k))
        yc(:,j,k,p) = r(:)/sqrt(delta(j,k))
        zc(:,j,k,p) = r(:)*Y(k)/sqrt(delta(j,k))
      ! patch III 
        p = 3
        xc(:,j,k,p) = - r(:)/sqrt(delta(j,k))
        yc(:,j,k,p) = - r(:)*X(j)/sqrt(delta(j,k))
        zc(:,j,k,p) = r(:)*Y(k)/sqrt(delta(j,k))
      ! patch IV
        p = 4
        xc(:,j,k,p) = r(:)*X(j)/sqrt(delta(j,k))
        yc(:,j,k,p) = - r(:)/sqrt(delta(j,k))
        zc(:,j,k,p) = r(:)*Y(k)/sqrt(delta(j,k))
      ! patch V
        p = 5
        xc(:,j,k,p) = - r(:)*Y(k)/sqrt(delta(j,k))
        yc(:,j,k,p) = r(:)*X(j)/sqrt(delta(j,k))
        zc(:,j,k,p) = r(:)/sqrt(delta(j,k)) 
      ! patch VI  
        p = 6 
        xc(:,j,k,p) = r(:)*Y(k)/sqrt(delta(j,k))
        yc(:,j,k,p) = r(:)*X(j)/sqrt(delta(j,k))
        zc(:,j,k,p) = - r(:)/sqrt(delta(j,k))
      end do 
    end do 
 
    xc(:,nang/2,1,1)=xc(:,nang/2,nang,6)
    yc(:,nang/2,1,1)=yc(:,nang/2,nang,6)
    zc(:,nang/2,1,1)=zc(:,nang/2,nang,6)
    
    end subroutine cs_to_cartesian
     !!--------------------------------------------------------------------------
    !> @brief Subroutine cs_to_spherical
    !!
    !! In this subroutine we apply the transformations needed to obtain the 
    !! angular spherical coordinates using cubed sphere coordinates
    !!
    !!
    !! Code owners:
    !!    Clara Dehman
    !!    Daniele Viganò
    !!--------------------------------------------------------------------------
    subroutine cs_to_spherical
      implicit none
    
      real*8, dimension(1:2*nang-1,4) :: theta_mer4
      integer :: j, k, p 


      theta = 0.d0
      phi = 0.d0

      do j = 0, nang+1
        do k = 0, nang+1
          ! patch I 
          p = 1
            if (xi(j) >= 0d0) then
              phi(j,k,p) = xi(j)
            else
              phi(j,k,p) = 2d0*PI + xi(j)
            endif
            theta(j,k,p) = arctan_ratio(C(j),Y(k))
 
          ! patch II
          p = 2
            phi(j,k,p) = xi(j) + 0.5d0*PI
            theta(j,k,p) = arctan_ratio(C(j),Y(k))
  
          ! patch III 
          p = 3
            phi(j,k,p) = xi(j) + PI
            theta(j,k,p) = arctan_ratio(C(j),Y(k))
          
          ! patch IV 
          p = 4
            phi(j,k,p) = xi(j) + 1.5d0*PI
            theta(j,k,p) = arctan_ratio(C(j),Y(k))

          ! patch V 
          p = 5
            phi(j,k,p) = arctan_ratio_2pi(X(j),-Y(k))
            theta(j,k,p) = datan(sqrt(delta(j,k)-1d0))

          ! patch VI
          p = 6
            phi(j,k,p) = arctan_ratio_2pi(X(j),Y(k))
            theta(j,k,p) = PI - datan(sqrt(delta(j,k)-1d0))

        end do 
      end do

      ! For sake of simplicity I call twice the subroutine 
      ! and keep the only one profile of theta and phi 
      ! phi arrives at 2pi, it's periodical
      ! theta goes from 0 to pi
      call get_1d_cuts(theta,(nang+1)/2,theta_mer4,phi_equator)
      theta_meridian(:) = theta_mer4(:,1)
      
      call get_1d_cuts(phi,(nang+1)/2,theta_mer4,phi_equator)
      phi_equator(4*nang - 3) = 2d0*PI

      ! theta_meridian_2PI: 
      theta_meridian_2PI(1:2*nang-1) = theta_meridian(:)
      theta_meridian_2PI(2*nang:4*nang-3) = PI + theta_meridian(2:2*nang-1)

    end subroutine cs_to_spherical

   
    !!--------------------------------------------------------------------------
    !> @brief Subroutine get_1d_cuts
    !!
    !! Here we build the meridional and azimuthal cut of fields,
    !! in order to obtain 1D cuts, used in output mainly.
    !! It can apply to the positions themselves.
    !! The only real 1D profiles possible are: the equatorial one, 
    !! and the four meridians passing through the center of the equatorial patches
    !!
    !! @param[in]   fin            the field that needs to be cut
    !! @param[in]   k              index of the pseudo-azimuthal cut (for equator it's exact, k=(nang+1)/2))
    !! @param[out]  f_mer, f_azi   1D profiles  
    !!
    !! Code owners:
    !!    Daniele Viganò
    !!    Clara Dehman 
    !!--------------------------------------------------------------------------
    subroutine get_1d_cuts(fin,k,f_mer,f_azi)

      implicit None

      real*8, dimension(0:nang+1,0:nang+1,1:6), intent(in) :: fin
      integer, intent(in) :: k
      real*8, dimension(1:2*nang-1,4), intent(out) :: f_mer
      real*8, dimension(1:4*nang-3), intent(out) :: f_azi

      integer c

      c = (nang+1)/2   ! Index of the center

      ! 0 degrees meridian
      f_mer(1:c,1) = fin(c,c:1:-1,5)
      f_mer(c+1:3*c-2,1) = fin(c,nang-1:1:-1,1)
      f_mer(3*c-1:2*nang-1,1) = fin(c,nang-1:c:-1,6)

      ! 90 degrees meridian
      f_mer(1:c,2) = fin(c:nang,c,5)
      f_mer(c+1:3*c-2,2) = fin(c,nang-1:1:-1,2)
      f_mer(3*c-1:2*nang-1,2) = fin(nang-1:c:-1,c,6)
    
      ! 180 degrees meridian
      f_mer(1:c,3) = fin(c,c:nang,5)
      f_mer(c+1:3*c-2,3) = fin(c,nang-1:1:-1,3)
      f_mer(3*c-1:2*nang-1,3) = fin(c,2:c,6)
    
      ! 270 degrees meridian
      f_mer(1:c,4) = fin(c:1:-1,c,5)
      f_mer(c+1:3*c-2,4) = fin(c,nang-1:1:-1,4)
      f_mer(3*c-1:2*nang-1,4) = fin(2:c,c,6)
    
      ! Azimuthal profile
      f_azi(1:c) = fin(c:nang,k,1) 
      f_azi(c+1:(3*nang-1)/2) = fin(2:nang,k,2)
      f_azi((3*nang+1)/2:(5*nang-3)/2) = fin(2:nang,k,3)
      f_azi((5*nang-1)/2:(7*nang-5)/2) = fin(2:nang,k,4)
      f_azi((7*nang-3)/2:4*nang-3) = fin(2:c,k,1) 

    end subroutine get_1d_cuts

    !!--------------------------------------------------------------------------
    !> @brief Subroutine get_1d_cuts
    !!
    !! Here we build the meridional and azimuthal cut of fields,
    !! in order to obtain 1D cuts, used in output mainly.
    !! It can apply to the positions themselves.
    !! The only real 1D profiles possible are: the equatorial one, 
    !! and the four meridians passing through the center of the equatorial patches
    !!
    !! @param[in]   fin            the field that needs to be cut
    !! @param[in]   k              index of the pseudo-azimuthal cut (for equator it's exact, k=(nang+1)/2))
    !! @param[out]  f_mer, f_azi   1D profiles  
    !!
    !! Code owners:
    !!    Stefano Ascenzi
    !!    Daniele Viganò
    !!    Clara Dehman 
    !!--------------------------------------------------------------------------
    subroutine get_1d_cuts_thermal(fin,k,f_mer,f_azi)

      implicit None

      real*8, dimension(0:nangt+1,0:nangt+1,1:6), intent(in) :: fin
      integer, intent(in) :: k
      real*8, dimension(1:2*nangt+1,4), intent(out) :: f_mer
      real*8, dimension(1:4*nangt), intent(out) :: f_azi
       
      integer c

      f_mer = 0d0
      f_azi = 0d0
      c = (nangt+1)/2   ! Index of the center

      ! 0 degrees meridian
      f_mer(1:c,1) = fin(c,c:1:-1,5)
      f_mer(c+1:c+nangt,1) = fin(c,nangt:1:-1,1)
      f_mer(c+1+nangt:2*nangt+1,1) = fin(c,nangt:c:-1,6)
      
      ! 90 degrees meridian
      f_mer(1:c,2) = fin(c:nangt,c,5)
      f_mer(c+1:c+nangt,2) = fin(c,nangt:1:-1,2)
      f_mer(c+1+nangt:2*nangt+1,2) = fin(nangt:c:-1,c,6)

      ! 180 degrees meridian
      f_mer(1:c,3) = fin(c,c:nangt,5)
      f_mer(c+1:c+nangt,3) = fin(c,nangt:1:-1,3)
      f_mer(c+1+nangt:2*nangt+1,3) = fin(c,1:c,6)
    
      ! 270 degrees meridian
      f_mer(1:c,4) = fin(c:1:-1,c,5)
      f_mer(c+1:c+nangt,4) = fin(c,nangt:1:-1,4)
      f_mer(c+1+nangt:2*nangt+1,4) = fin(1:c,c,6)
      
      ! Azimuthal profile
      f_azi(1:c) = fin(c:nangt,k,1) 
      f_azi(c+1:c+nangt) = fin(1:nangt,k,2)
      f_azi(c+1+nangt:c+2*nangt) = fin(1:nangt,k,3)
      f_azi(c+1+2*nangt:c+3*nangt) = fin(1:nangt,k,4)
      f_azi(c+3*nangt+1:4*nangt) = fin(1:c-1,k,1) 
      
    end subroutine get_1d_cuts_thermal



    !!--------------------------------------------------------------------------
    !> @brief Subroutine get_2d_cuts
    !!
    !! Here we build the meridional and azimuthal cut of fields,
    !! in order to obtain 2D cuts, used in output mainly. 
    !! This subroutine is used in order to visualize a the meridional and equatorial 
    !! projections in all the volume at different radii
    !!
    !! for the meridional profiles, we plot longitudinal 0 degrees with longitudinal 180 degrees
    !! and longitudinal 90 degrees with longitudinal 270 degrees
    !!
    !! @param[in]   fin            the field that needs to be cut
    !! @param[in]   k              index of the pseudo-azimuthal cut (for equator it's exact, k=(nang+1)/2))
    !! @param[out]  f_mer, f_azi   2D profiles  
    !!
    !! Code owners:
    !!    Clara Dehman 
    !!--------------------------------------------------------------------------
    subroutine get_2d_cuts(fin,k,tang_comp,f_mer,f_azi)

      implicit None

      real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6), intent(in) :: fin
      integer, intent(in) :: k
      logical, intent(in) :: tang_comp
      real*8, dimension(0:nr+1,1:4*nang-3,2), intent(out) :: f_mer
      real*8, dimension(0:nr+1,1:4*nang-3), intent(out) :: f_azi

      integer c

      c = (nang+1)/2   ! Index of the center
       
      ! 0-180 degrees meridian
       f_mer(:,1:c,1) =  fin(:,c,c:1:-1,5)
       f_mer(:,c+1:(3*nang-1)/2,1) = fin(:,c,nang-1:1:-1,1)
       f_mer(:,(3*nang+1)/2:(5*nang-3)/2,1) = fin(:,c,nang-1:1:-1,6)
       f_mer(:,(5*nang-1)/2:(7*nang-5)/2,1) = fin(:,c,2:nang,3)
       f_mer(:,(7*nang-3)/2:4*nang-3,1) = fin(:,c,nang-1:c:-1,5)
       
      ! 90-270 degrees meridian
       f_mer(:,1:c,2) = fin(:,c:nang,c,5)
       f_mer(:,c+1:(3*nang-1)/2,2) = fin(:,c,nang-1:1:-1,2)
       f_mer(:,(3*nang+1)/2:(5*nang-3)/2,2) = fin(:,nang-1:1:-1,c,6)
       f_mer(:,(5*nang-1)/2:(7*nang-5)/2,2) = fin(:,c,2:nang,4)
       f_mer(:,(7*nang-3)/2:4*nang-3,2) = fin(:,2:c,c,5)

      if (tang_comp) then 
       ! averaging at the axis for visualization, if it's a tangential component
       f_mer(:,1,1:2) = 0.5d0*(f_mer(:,2,1:2) + f_mer(:,4*nang-4,1:2))
       f_mer(:,4*nang-3,1:2) = f_mer(:,1,1:2)
       f_mer(:,2*nang-1,1:2) = 0.5d0*(f_mer(:,2*nang-2,1:2) + f_mer(:,2*nang,1:2)) 
      end if 

      ! Azimuthal profile
      f_azi(:,1:c) = fin(:,c:nang,k,1) 
      f_azi(:,c+1:(3*nang-1)/2) = fin(:,2:nang,k,2)
      f_azi(:,(3*nang+1)/2:(5*nang-3)/2) = fin(:,2:nang,k,3)
      f_azi(:,(5*nang-1)/2:(7*nang-5)/2) = fin(:,2:nang,k,4)
      f_azi(:,(7*nang-3)/2:4*nang-3) = fin(:,2:c,k,1) 

    end subroutine get_2d_cuts

    ! This function is used to define correctly theta [0,pi] in patches I-II-III-IV
    real*8 function arctan_ratio(num,den)

      implicit None
      real*8, intent(in) :: num, den

      if (abs(den) < 1d-10) then
        arctan_ratio = 0.5d0*PI
      else
        if ( (num/den) > 0d0) then
          arctan_ratio = datan(num/den)
        else 
          arctan_ratio = PI + datan(num/den)
        endif 
      endif

    end function

    ! This function is used to define correctly phi [0,2pi] in patches V-VI
    real*8 function arctan_ratio_2pi(num,den)

      implicit None
      real*8, intent(in) :: num, den
  
      if (abs(den) < 1d-10 .and. num > 0d0) then
        arctan_ratio_2pi = 0.5d0*PI
      else if (abs(den) < 1d-10 .and. num < 0d0) then
        arctan_ratio_2pi = 1.5d0*PI
      else
        ! region b in patch V, region a in patch VI
        if ( (num >= 0d0) .and. (den < 0d0) ) then
          arctan_ratio_2pi = PI + datan(num/den)
        ! region a in patch V, region b in patch VI
        elseif ((num >= 0d0) .and. (den > 0d0) ) then
          arctan_ratio_2pi = datan(num/den)
        ! region d in patch V, region c in patch VI
        elseif ((num < 0d0) .and. (den > 0d0) ) then
          arctan_ratio_2pi = 2d0*PI + datan(num/den)
        ! region c in patch V, region d in patch VI
        elseif ((num < 0d0) .and. (den < 0d0) ) then
          arctan_ratio_2pi = PI + datan(num/den)  
        endif
      endif
  
    end function


    !!--------------------------------------------------------------------------
    !> @brief Subroutine f_cs_to_spherical
    !!
    !! In this subroutine we apply the transformations needed for vectors to go from 
    !! cubed sphere coordinate to spherical coordinate 
    !! 
    !! @param[in]   fxi                      contravariant xi component
    !! @param[in]   feta                     contravariant eta component
    !! @param[out]  fth                      contravariant th component
    !! @param[out]  fphi                     contravariant phi component
    !!
    !! Code owners:
    !!    Clara Dehman
    !!    Daniele Viganò
    !!--------------------------------------------------------------------------
    subroutine f_cs_to_spherical(fxi,feta,fth,fphi,imin)

      implicit none
      integer, intent(in) :: imin  ! Minimum index to start the calculation
      real*8, dimension (0:nang+1, 0:nang+1, 2, 2) :: jacinv
      real*8, dimension (0:nr+1, 0:nang+1, 0:nang+1, 6), intent(in) :: fxi, feta
      real*8, dimension (0:nr+1, 0:nang+1, 0:nang+1, 6), intent(out) :: fth, fphi

      integer i, p

      jacinv = 0d0
      fth = 0d0
      fphi = 0d0

      do p = 1 , 6
        if (p <= 4) then
          jacinv = jacinv_eq
        elseif (p == 5) then
          jacinv = jacinv_n
        else
          jacinv = jacinv_s
        endif

        do i = imin, nr + 1
          fth(i,:,:,p) = jacinv(:,:,1,1)*fxi(i,:,:,p) + jacinv(:,:,1,2)*feta(i,:,:,p)
          fphi(i,:,:,p) = jacinv(:,:,2,1)*fxi(i,:,:,p) + jacinv(:,:,2,2)*feta(i,:,:,p)
        end do

      end do 

      fth(:,nang/2+1,nang/2+1,5:6) = 0.d0
      fphi(:,nang/2+1,nang/2+1,5:6) = 0.d0

    end subroutine f_cs_to_spherical

    !!--------------------------------------------------------------------------
    !> @brief Subroutine f_spherical_to_cs
    !!
    !! In this subroutine we apply the transformations needed for vectors to go from 
    !! spherical coordinate to cubed sphere coordinate
    !! 
    !! @param[in]   fth                      contravariant th component
    !! @param[in]   fphi                     contravariant phi component
    !! @param[out]  fxi                      contravariant xi component
    !! @param[out]  feta                     contravariant eta component
    !!
    !! Code owners:
    !!    Clara Dehman
    !!    Daniele Viganò
    !!--------------------------------------------------------------------------
    subroutine f_spherical_to_cs(fth,fphi,fxi,feta,imin)

      implicit none
      integer, intent(in) :: imin  ! Minimum index to start the calculation
      real*8, dimension (0:nang+1,0:nang+1,2,2) :: jac
      real*8, dimension (imin:nr+1, 0:nang+1, 0:nang+1, 6), intent(in) :: fth, fphi
      real*8, dimension (imin:nr+1, 0:nang+1, 0:nang+1, 6), intent(out) :: fxi, feta

      integer i, p

      fxi(:,:,:,:) = 0d0
      feta(:,:,:,:) = 0d0
      jac(:,:,:,:) = 0d0

      do p = 1 , 6
        if (p <= 4) then
          jac = jac_eq
        elseif (p == 5) then
          jac = jac_n
        else
          jac = jac_s
        endif

        do i = imin, nr + 1
          fxi(i,:,:,p) = jac(:,:,1,1)*fth(i,:,:,p) + jac(:,:,1,2)*fphi(i,:,:,p)
          feta(i,:,:,p) = jac(:,:,2,1)*fth(i,:,:,p) + jac(:,:,2,2)*fphi(i,:,:,p)
        end do
      end do

      ! At the axis, we need an average, since the phi components are ill-defined
      ! Note anyway that this is used only for the input, output and BC
      ! We use a higher-order volume-average with the first and second angular neighbours
      i = (nang+1)/2   ! index of eta and xi of the center of the patches 5 and 6

      fxi(:,i,i,5:6) = 1d0/12d0*(fxi(:,i-1,i-1,5:6) + fxi(:,i+1,i+1,5:6)  & 
      & + fxi(:,i+1,i-1,5:6) + fxi(:,i-1,i+1,5:6) + 2d0*(fxi(:,i,i+1,5:6)  &
      & + fxi(:,i,i-1,5:6) + fxi(:,i-1,i,5:6) + fxi(:,i+1,i,5:6)) )

      feta(:,i,i,5:6) = 1d0/12d0*(feta(:,i-1,i-1,5:6) + feta(:,i+1,i+1,5:6)  & 
      & + feta(:,i+1,i-1,5:6) + feta(:,i-1,i+1,5:6) + 2d0*(feta(:,i,i+1,5:6)  &
      & + feta(:,i,i-1,5:6) + feta(:,i-1,i,5:6) + feta(:,i+1,i,5:6)) )

    end subroutine f_spherical_to_cs

    !!--------------------------------------------------------------------------
    !> @brief Subroutine f_spherical_to_cartesian
    !!
    !! In this subroutine we apply the transformations needed for vectors to go from 
    !! spherical coordinate to Cartesian coordinate
    !! 
    !! @param[in]   fth                      r component
    !! @param[in]   fth                      th component
    !! @param[in]   fphi                     phi component
    !! @param[out]  fx                       x component
    !! @param[out]  fy                       y component
    !! @param[out]  fz                       z component
    !!
    !! Code owners:
    !!    Daniele Viganò
    !!--------------------------------------------------------------------------
    subroutine f_spherical_to_cartesian(fr,fth,fphi,fx,fy,fz,imin)

      implicit none
      integer, intent(in) :: imin  ! Minimum index to start the calculation
      real*8, dimension (0:nang+1,0:nang+1,3,3) :: jac
      real*8, dimension (0:nr+1, 0:nang+1, 0:nang+1, 6), intent(in) :: fr, fth, fphi
      real*8, dimension (0:nr+1, 0:nang+1, 0:nang+1, 6), intent(out) :: fx, fy, fz

      integer i, p

      fx(:,:,:,:) = 0d0
      fy(:,:,:,:) = 0d0
      fz(:,:,:,:) = 0d0
      jac(:,:,:,:) = 0d0


      do p = 1 , 6

        do i = imin, nr + 1
          jac(:,:,1,1) = dsin(theta(:,:,p))*dcos(phi(:,:,p))
          jac(:,:,1,2) = dcos(theta(:,:,p))*dcos(phi(:,:,p))
          jac(:,:,1,3) = - dsin(phi(:,:,p))
          jac(:,:,2,1) = dsin(theta(:,:,p))*dsin(phi(:,:,p))
          jac(:,:,2,2) = dcos(theta(:,:,p))*dsin(phi(:,:,p))
          jac(:,:,2,3) = dcos(phi(:,:,p))
          jac(:,:,3,1) = dcos(theta(:,:,p))
          jac(:,:,3,2) = - dsin(theta(:,:,p))

          fx(i,:,:,p) = jac(:,:,1,1)*fr(i,:,:,p) + jac(:,:,1,2)*fth(i,:,:,p) + jac(:,:,1,3)*fphi(i,:,:,p)
          fy(i,:,:,p) = jac(:,:,2,1)*fr(i,:,:,p) + jac(:,:,2,2)*fth(i,:,:,p) + jac(:,:,2,3)*fphi(i,:,:,p)
          fz(i,:,:,p) = jac(:,:,3,1)*fr(i,:,:,p) + jac(:,:,3,2)*fth(i,:,:,p) + jac(:,:,3,3)*fphi(i,:,:,p)
        end do
      
      end do

    end subroutine f_spherical_to_cartesian

    !!--------------------------------------------------------------------------
    !> @brief Subroutine spherical_harmonics
    !!
    !! In this subroutine we calculate the associated legendre polynomials and 
    !! the spherical harmonics, as well as their theta and phi derivatives 
    !! needed to reconstruct the component of the magnetic field at the surface 
    !! of the star and above it (thet 1 ghost cell above the surface). 
    !!            -- Outer magnetic boundary conditions --
    !!
    !! Code owners:
    !!    Clara Dehman
    !!    Daniele Viganò
    !!--------------------------------------------------------------------------
    subroutine spherical_harmonics

      implicit none
      integer l, m
      real*8 :: factorial(0:2*lmax), double_factorial(1:lmax)
      real*8, dimension(0:nang+1, 0:nang+1, 1:6, 0:lmax, 0:lmax) :: p_lm, dpth_lm

    ! integer ln, mn, p, j, k
!      real*8 :: surf_check
   !  real*8, dimension(0:lmax,-lmax:lmax,0:lmax,-lmax:lmax) :: ortho   ! For the orthogonality check
 
      factorial(:) = 0d0
      factorial(0) = 1d0
      factorial(1) = 1d0
      do l = 2, 2*lmax
        factorial(l) = factorial(l-1)*dble(l)
      end do 
 
      double_factorial(:) = 0d0
      double_factorial(1) = 1d0
      ! The double factorial is defined so that: double_factorial(l) := (2l-1)!!
      do l = 2, lmax
        double_factorial(l) = double_factorial(l-1)*dble(2*l-1)
      end do 
     ! ------------------------------------------------------------------------- 
     ! p_lm associated legender polynomial defined from l:0->lmax and m:0->l
     ! Recurrence relations: https://mathworld.wolfram.com/AssociatedLegendrePolynomial.html eq. 7
     
      p_lm(:,:,:,0,0) = 1.d0
      p_lm(:,:,:,1,0) = dcos(theta(:,:,:))
      p_lm(:,:,:,1,1) = - dsin(theta(:,:,:))
 
     ! Used to calculate Y_lm below
      do l = 2, lmax
      do m = 0, l-2 
        p_lm(:,:,:,l,m) = 1.d0/dble(l-m)*( dcos(theta(:,:,:))* &
    &           dble(2*l-1)*p_lm(:,:,:,l-1,m) - dble(l+m-1)*p_lm(:,:,:,l-2,m))
      end do
      ! For the case m=l and m=l-1 we have 
      ! https://en.wikipedia.org/wiki/Associated_Legendre_polynomials , Recurrence relations
      p_lm(:,:,:,l,l-1) = dcos(theta(:,:,:))*dble(2*l-1)*p_lm(:,:,:,l-1,l-1)
      p_lm(:,:,:,l,l) = (-1)**l*double_factorial(l)*dsin(theta(:,:,:))**l

      end do  
 
     ! ------------------------------------------------------------------------- 
     ! Theta derivative of p_lm defined from l:0->lmax and m:0->l
     ! Used to calculate the derivative of Y_lm below
      dpth_lm(:,:,:,0,0) = 0d0
      dyth_lm(:,:,:,0,0) = 0d0
      do l = 1, lmax
        dpth_lm(:,:,:,l,0) = p_lm(:,:,:,l,1)              ! m=0
        do m = 1, l-1
          dpth_lm(:,:,:,l,m) = - 1/2d0*(dble((l+m)*(l-m+1))*p_lm(:,:,:,l,m-1) - p_lm(:,:,:,l,m+1))
        end do 
        dpth_lm(:,:,:,l,l) = - dble(l)*p_lm(:,:,:,l,l-1)  ! m=l
      end do 
  
        do l = 0, lmax
        ! m = 0 case
        y_lm(:,:,:,l,0) = dsqrt(dble(2*l+1)/(4.d0*PI))*p_lm(:,:,:,l,0)
        ! Note: for l=0, dyth_lm = 0
        dyth_lm(:,:,:,l,0) = dsqrt(dble(2*l+1)/(4.d0*PI))*dpth_lm(:,:,:,l,0)
        dyphi_lm(:,:,:,l,0) = 0.d0

        do m = 1, l
         ! y_lm is the real spherical harmonics defined from l:0->lmax and m:-l->l
         ! negative m are storing the sin(m*phi) branch, positive m the cos(m*phi) branch
          y_lm(:,:,:,l,m) = dsqrt(dble(2*l+1)/(2.d0*PI)*factorial(l-m)/factorial(l+m)) &
         &                   *p_lm(:,:,:,l,m)*dcos(dble(m)*phi(:,:,:))
          y_lm(:,:,:,l,-m) = dsqrt(dble(2*l+1)/(2.d0*PI)*factorial(l-m)/factorial(l+m)) &
         &                   *p_lm(:,:,:,l,m)*dsin(dble(m)*phi(:,:,:))
         ! First order theta derivative of y_lm defined from l:0->lmax and m:-l->l
         ! They are used to calculate Btheta at the BC
          dyth_lm(:,:,:,l,m) = dsqrt(dble(2*l+1)/(2.d0*PI)*factorial(l-m)/factorial(l+m)) &
         &                      *dcos(dble(m)*phi(:,:,:))*dpth_lm(:,:,:,l,m)
          dyth_lm(:,:,:,l,-m) = dsqrt(dble(2*l+1)/(2.d0*PI)*factorial(l-m)/factorial(l+m)) & 
         &                      *dsin(dble(m)*phi(:,:,:))*dpth_lm(:,:,:,l,m)
  
          ! First order phi derivative of y_lm defined from l:0->lmax and m:-l->l
          ! They are used to calculate Bphi at the BC
          dyphi_lm(:,:,:,l,m) = - dble(m)*dsqrt(dble(2*l+1)/(2.d0*PI)*factorial(l-m)/factorial(l+m)) & 
     &                            * p_lm(:,:,:,l,m)*dsin(dble(m)*phi(:,:,:))
          dyphi_lm(:,:,:,l,-m) = dble(m)*dsqrt(dble(2*l+1)/(2.d0*PI)*factorial(l-m)/factorial(l+m)) & 
     &                            * p_lm(:,:,:,l,m)*dcos(dble(m)*phi(:,:,:))

        end do 
      end do 

    end subroutine spherical_harmonics


    ! --------------------------------------------------------------------------
    !> Relativistic correction for the vacuum boundary conditions
    !> @brief   It calculates the Barnes hypergeometrical functions, used for the magnetic BC
    
    ! Note: anu_rel_correction: see Eq. (56) of Radler et al. 2001, PRD 64
    ! Note: frel: see f_n factors in Pons, Miralles and Geppert 2009, Eq. (24).

    !> @param[in] eps double_compactness GM/c^2R
    !
    !! Code owners:
    !!    Daniele Viganò
    ! --------------------------------------------------------------------------
    subroutine get_rel_correction(double_compactness)

      ! Modules ----------------------------------------------------------------
      ! None.

      implicit none

      ! Input arguments --------------------------------------------------------
      real*8, intent(in) :: double_compactness  

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
      integer :: i, l
      ! i is the nu of Radler et al. 2001, sec. VI, n is the l
      real*8 :: bb, sum1, sum2, frel_sl 
      ! bb is the b_nu of eq. (58)
      ! frel_sl is eq. (57)
      ! sum1 and sum2 are the two sums appearing in eq. (57)
      ! sum2 is eq. (56) where the sum_nu(a_nu compactness^nu) is equivalent to sum2,
      ! which is the denominator of (58)
      ! the conversion from frel_sl (used for S_l in Radler 2001) 
      ! to frel (used for Phi_l in Pons 2009)
      ! considers that S_l = (1/r)Phi_l

      ! Initialize (used for nonrelativistic case)
      frel = 1d0

      if (double_compactness /= 0.) then

        do l = 1, lmax

          bb = 1d0
          sum1 = bb
          sum2 = bb
          
          ! First 50 terms (more than enough) of the infinite series in nu=0,infinity
          do i = 1, 50
  
            bb = dble((l + i)**2 - 1) * double_compactness * bb / dble((2 * l + i + 1) * i)
            sum1 = sum1 + dble(l + i + 1) * bb / dble(l + 1)
            sum2 = sum2 + bb
  
          end do ! i
  
          frel_sl = sum1 / sum2
          frel(l) = frel_sl + (frel_sl - 1d0) / dble(l)
        
        end do ! l

      endif
      

    end subroutine get_rel_correction


    subroutine get_wint

      implicit none

      integer j, k
  
      wint = 0d0
      do j = 0, nang + 1
        do k = 0, nang + 1
          if ((j .eq. 1 .and. k .eq. 1) .or. (j .eq. 1 .and. k .eq. nang) &
        & .or. (j .eq. nang .and. k .eq. 1) .or. (j .eq. nang .and. k .eq. nang)) then
            ! We are on a corner, triple counted
            wint(j,k) = 1d0/3d0
          else
            if ((j .eq. 1 .or. j .eq. nang) .or. (k .eq. 1 .or. k .eq. nang)) then
              ! We are on an edge, double counted
              wint(j,k) = 0.5d0
            else
              ! Interior points, no corrections
              wint(j,k) = 1
            endif
          endif
        end do
      end do
  
    end subroutine get_wint
  
end module grid

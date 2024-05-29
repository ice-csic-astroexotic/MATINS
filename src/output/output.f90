!-------------------------------------------------------------------------------
! Magneto Thermal 3D
!-------------------------------------------------------------------------------
! Module: Output
!
!> @author
!> Daniele Viganò
!> Clara Dehman
!> Albert Herrando 
!> Stefano Ascenzi
!
!> Output handling.
!> @brief This module is responsible for creating all ygraph output for the simulation,
!>        for the magnetic field evolution.
!
!-------------------------------------------------------------------------------

Module output 

  use utils, only: get_free_unit, real_to_string, int_to_string
  use grid, only: nr, nang, nrt, nangt, ncore, ievol
  use grid, only: r, rtot, theta, phi, theta_meridian, phi_equator, theta_meridian_2PI
  use grid, only: vol, area_r, elambda
  use grid, only: eta, xi, xc, yc, zc
  use grid, only: br, bxi, beta, bm, bpdip
  use grid, only: lmax, blm, espec_vol, espec_pol, espec_tor, phi_scalar, psi_scalar
  use grid, only: jr, jxi, jeta
  use grid, only: er, exi, eeta
  use grid, only: flux_r_out, flux_xi_xip, flux_eta_etap
  use grid, only: en_joule_star_tot, poynting_star_tot, poynting_star_tot_surface, poynting_star_tot_interior
  use grid, only: temp, tem0, T_core, temp_surf, bb_flux
  use grid, only: cv, cv_core, cv_core_tot
  use grid, only: q_joule, q_neutrino, q_neutrino_core, qnu_core_tot
  use grid, only: omegatau_arr, kappa_perp_arr, etab
  use grid, only: f_spherical_to_cartesian, f_cs_to_spherical
  use grid, only: get_1d_cuts, get_1d_cuts_thermal, get_2d_cuts
  use grid, only: enu 
  use grid, only: qnu_mur, qnu_nn, qnu_np, &
      &            qnu_pp, qnu_ep, qnu_cp_con,qnu_cp_cop, qnu_du, &
      &            qnu_ea, qnu_pl, qnu_syn, qnu_cp_cr,qnu_pa
  use constants, only: UNIT_B, UNIT_T, UNIT_R, UNIT_EN, UNIT_TIME, PI
  use constants, only: PI, STEFAN_BOLTZMANN

  implicit None
  character(len=10), parameter :: YG_FORMAT = "(2es12.4)"

  contains

    !---------------------------------------------------------------------------
    !> Outputs in 1D of the magnetic field simulation.
    !
    !> authors
    ! Daniele Viganò
    ! Clara Dehman
    !---------------------------------------------------------------------------
    subroutine output_magnetic_1D(ib,timeMyr)

      use magnetic_analysis, only: analyse_magnetic_field
    
      implicit none

      ! Input parameters -------------------------------------------------------
      real*8, intent(in) :: timeMyr
      integer, intent(in) :: ib
      integer :: unit


      ! Local constants --------------------------------------------------------
      character(len=14), parameter :: ANALYSIS_FORMAT = "(i7,12es12.4)"
      character(len=17), parameter :: BL_FORMAT = "(f5.0,11es15.6)"
      ! ------------------------------------------------------------------------
      real*8, dimension (0:nr+1,0:nang+1,0:nang+1,1:6) :: bth, bphi
      real*8, dimension (0:nr+1,0:nang+1,0:nang+1,1:6) :: jth, jphi
      real*8, dimension (0:nr+1,0:nang+1,0:nang+1,1:6) :: eth, ephi
      real*8 :: b_avg, en_mag_star, divb_L2, j2_star, rey, rey_max, t_hall
      integer :: l, m
      real*8, dimension(1:lmax,-lmax:lmax) :: blm_surf
      real*8, dimension(1:lmax) :: blm_l
      real*8, dimension(-lmax:lmax) :: blm_m, poles
      logical :: isthermalgrid = .false.

      ! Global quantity
      call analyse_magnetic_field(b_avg, en_mag_star, divb_L2, j2_star, rey, rey_max, t_hall)
      unit = get_free_unit()
    
      ! If we just started the simulation, create a new file and write header.
      if (timeMyr == 0.0d0) then
        open(unit=unit, file="out/energy/mag_energy_balance.dat")
        write(unit, "(a80)") 'iter, time[yr], avgB, Emag, Qjtot, Poy_surf, Poy_int, J2, div(B)L2, avgRm, maxRm, t_hall'
        write(6, "(a80)") 'Mag.evo:iter,t[yr],Avg.B,Emag,Joule,Poy_surf,Poy_int,J2,divB,Rm,maxRm,tHall'
        ! If not, just append to the end of the file.
      else
        open(unit=unit, file="out/energy/mag_energy_balance.dat", access="append")
      end if
      write(unit, ANALYSIS_FORMAT) ib, timeMyr*1d6, b_avg, en_mag_star, en_joule_star_tot, poynting_star_tot_surface,  &
    &   poynting_star_tot_interior, j2_star, divb_L2, rey, rey_max, t_hall
      close(unit)
      write(6, ANALYSIS_FORMAT) ib, timeMyr*1d6, b_avg, en_mag_star, en_joule_star_tot, poynting_star_tot_surface,  &
    &   poynting_star_tot_interior, j2_star, divb_L2, rey, rey_max, t_hall
      
      ! Multipoles weights at the surface 
    do l = 1, lmax
      ! blm_surf(l,:) = blm(l,:)**2*(l+1)*(2*l+1)*UNIT_B**2 ! Newtonian limit 
      blm_surf(l,:) = blm(l,:)**2*((l+1)**2/(elambda(nr))**2+l*(l+1))*UNIT_B**2 ! relativistic limit 
      blm_l(l) = sum(blm_surf(l,:))
    end do
    do m = -lmax, lmax
      poles(m) = m
      blm_m(m) = sum(blm_surf(:,m))
    end do
      ! Square of the spherical harmonics weights, integrated over m
      call output_1d_ygraph("out/energy/blm_surface_l.yg", "blm surf int. over m", &
     & timeMyr, poles(1:lmax), blm_l, lmax, BL_FORMAT)
      ! Square of the spherical harmonics weights, integrated over l
      call output_1d_ygraph("out/energy/blm_surface_m.yg", "blm surf int. over l", &
     & timeMyr, poles(-lmax:lmax), blm_m, 2*lmax+1, BL_FORMAT)
  

      ! Component profiles
      call f_cs_to_spherical(bxi,beta,bth,bphi,ievol-2) !ievol-2 is the center of the sphere
      call f_cs_to_spherical(jxi,jeta,jth,jphi,ievol-1)
      call f_cs_to_spherical(exi,eeta,eth,ephi,ievol-1)

      call all_profiles(br, "br", timeMyr, isthermalgrid)
      call all_profiles(bth, "bth", timeMyr, isthermalgrid)
      call all_profiles(bphi, "bphi", timeMyr, isthermalgrid)
      call all_profiles(jr, "jr", timeMyr, isthermalgrid)
      call all_profiles(jth, "jth", timeMyr, isthermalgrid)
      call all_profiles(jphi, "jphi", timeMyr, isthermalgrid)
      call all_profiles(er, "er", timeMyr, isthermalgrid)
      call all_profiles(eth, "eth", timeMyr, isthermalgrid)
      call all_profiles(ephi, "ephi", timeMyr, isthermalgrid)
      call all_profiles(etab, "etab", timeMyr, isthermalgrid)

    end subroutine output_magnetic_1D

    !---------------------------------------------------------------------------
    !> subroutine to write the output for the 2D meridional and equatorial cuts. 
    !> Moreover, we also include the multipoles weight through the volume 
    !
    !> authors
    ! Clara Dehman
    ! Daniele Viganò
    !---------------------------------------------------------------------------
    subroutine output_magnetic_2D(timeMyr)
    
      use magnetic_analysis, only: energy_spectrum

      implicit none

      ! Input parameters -------------------------------------------------------
      real*8, intent(in) :: timeMyr
      integer :: unit
      integer :: i, j, k, p, l, m
      real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: bth, bphi
      real*8, dimension(0:nr+1,1:4*nang-3,2) :: br_merd, bth_merd, bphi_merd, phisc_merd, psisc_merd
      real*8, dimension(0:nr+1,1:4*nang-3) :: br_equator, bth_equator, bphi_equator, phisc_eq, psisc_eq
      real*8, dimension(1:lmax) :: espec_l, espec_axi_l, espec_nonaxi_l
      real*8, dimension(-lmax:lmax) :: espec_m, poles
      real*8 :: Espec, Epol, Etor

      ! Local constants --------------------------------------------------------
      character(len=10), parameter :: ANALYSIS_FORMAT = '(10es12.4)'
      character(len=11), parameter :: BVOL_FORMAT = '(5es14.5)' 
      character(len=11), parameter :: BV_FORMAT = '(300es14.5)' 
      character(len=17), parameter :: BL_FORMAT = '(f5.0,11es15.6)'
      character(len=10), parameter :: BMAP_FORMAT = '(18es15.5)'
    !  character(len=50) :: file_name, time_str

      !-------------------------------------------------------------------------
      ! Magnetic eneryg spectrum 
      !-------------------------------------------------------------------------
      call energy_spectrum(Espec,Epol,Etor)
      
      ! Spectral energy through the volume
      do l = 1, lmax
        espec_l(l) = sum(espec_vol(l,:))
        espec_axi_l(l) = espec_vol(l,0) 
        espec_nonaxi_l(l) = sum(espec_vol(l,1:)) 
      end do
      do m = -lmax, lmax
        poles(m) = m
        espec_m(m) = sum(espec_vol(:,m))
      end do

      ! Square of the spherical harmonics weights, integrated over m
      call output_1d_ygraph("out/energy/espec_volume_l.yg", "blm vol int. over m", &
     & timeMyr, poles(1:lmax), espec_l, lmax, BL_FORMAT)
    ! Square of the spherical harmonics weights, at m=0
      call output_1d_ygraph("out/energy/espec_vol_axi_l.yg", "blm vol for m=0", &
     & timeMyr, poles(1:lmax), espec_axi_l, lmax, BL_FORMAT)
    ! Square of the spherical harmonics weights, integrated over m=1:lmax
      call output_1d_ygraph("out/energy/espec_vol_nonaxi_l.yg", "blm vol int. over m=1:lmax", &
     & timeMyr, poles(1:lmax), espec_nonaxi_l, lmax, BL_FORMAT)
    ! Square of the spherical harmonics weights, integrated over l
      call output_1d_ygraph("out/energy/espec_volume_m.yg", "blm vol int. over l", &
     & timeMyr, poles(-lmax:lmax), espec_m, 2*lmax+1, BL_FORMAT)

     ! Poloidal spectral energy
     do l = 1, lmax
      espec_l(l) = sum(espec_pol(l,:))
      espec_axi_l(l) = espec_pol(l,0) 
      espec_nonaxi_l(l) = sum(espec_pol(l,1:)) 
     end do
     do m = -lmax, lmax
      poles(m) = m
      espec_m(m) = sum(espec_pol(:,m))
     end do

    ! Square of the spherical harmonics weights, integrated over m
     call output_1d_ygraph("out/energy/espec_poloidal_l.yg", "blm pol int. over m", &
    & timeMyr, poles(1:lmax), espec_l, lmax, BL_FORMAT)
    ! Square of the spherical harmonics weights, at m=0
     call output_1d_ygraph("out/energy/espec_pol_axi_l.yg", "blm vol for m=0", &
     & timeMyr, poles(1:lmax), espec_axi_l, lmax, BL_FORMAT)
    ! Square of the spherical harmonics weights, integrated over m=1:lmax
     call output_1d_ygraph("out/energy/espec_pol_nonaxi_l.yg", "blm vol int. over m=1:lmax", &
     & timeMyr, poles(1:lmax), espec_nonaxi_l, lmax, BL_FORMAT)
    ! Square of the spherical harmonics weights, integrated over l
     call output_1d_ygraph("out/energy/espec_poloidal_m.yg", "blm pol int. over l", &
    & timeMyr, poles(-lmax:lmax), espec_m, 2*lmax+1, BL_FORMAT)

     ! toroidal spectral energy 
    do l = 1, lmax
    espec_l(l) = sum(espec_tor(l,:))
    espec_axi_l(l) = espec_tor(l,0) 
    espec_nonaxi_l(l) = sum(espec_tor(l,1:)) 
    end do
    do m = -lmax, lmax
    poles(m) = m
    espec_m(m) = sum(espec_tor(:,m))
    end do

  ! Square of the spherical harmonics weights, integrated over m
    call output_1d_ygraph("out/energy/espec_toroidal_l.yg", "blm tor int. over m", &
   & timeMyr, poles(1:lmax), espec_l, lmax, BL_FORMAT)
  ! Square of the spherical harmonics weights, at m=0
    call output_1d_ygraph("out/energy/espec_tor_axi_l.yg", "blm vol for m=0", &
   & timeMyr, poles(1:lmax), espec_axi_l, lmax, BL_FORMAT)
  ! Square of the spherical harmonics weights, integrated over m=1:lmax
    call output_1d_ygraph("out/energy/espec_tor_nonaxi_l.yg", "blm vol int. over m=1:lmax", &
   & timeMyr, poles(1:lmax), espec_nonaxi_l, lmax, BL_FORMAT)
  ! Square of the spherical harmonics weights, integrated over l
    call output_1d_ygraph("out/energy/espec_toroidal_m.yg", "blm tor int. over l", &
   & timeMyr, poles(-lmax:lmax), espec_m, 2*lmax+1, BL_FORMAT)


   ! output of spectral magnetic energy: toroidal and poloidal parts

    unit = get_free_unit()
   ! If we just started the simulation, create a new file and write header.
    if (timeMyr == 0.0d0) then
      open(unit=unit, file="out/energy/spectral_energy.dat")
      write(unit, *) 'time, Espectral, Epol, Etor'
      ! If not, just append to the end of the file.
    else
      open(unit=unit, file="out/energy/spectral_energy.dat", access="append")
    end if
    write(unit, "(4es12.4)") timeMyr, Espec, Epol, Etor
    close(unit)

    !-------------------------------------------------------------------------
    ! 2D cuts calculations
    !-------------------------------------------------------------------------
      call f_cs_to_spherical(bxi,beta,bth,bphi,ievol-2) 
      br = br*MERGE(0, 1, abs(br) .lt. 1e-20)
      bth = bth*MERGE(0, 1, abs(bth) .lt. 1e-20)
      bphi = bphi*MERGE(0, 1, abs(bphi) .lt. 1e-20)

      call get_2d_cuts(br,(nang+1)/2,.false.,br_merd,br_equator)
      call get_2d_cuts(bth,(nang+1)/2,.true.,bth_merd,bth_equator)    
      call get_2d_cuts(bphi,(nang+1)/2,.true.,bphi_merd,bphi_equator)
      call get_2d_cuts(phi_scalar,(nang+1)/2,.false.,phisc_merd,phisc_eq)
      call get_2d_cuts(psi_scalar,(nang+1)/2,.false.,psisc_merd,psisc_eq) 

    ! 2D meridional Profile: longitudinals 0 + 180 degrees 
      unit = get_free_unit()

      if (timeMyr == 0.0d0) then
        open(unit=unit, file="out/2D/b_merid_l0180_volume.dat")
        write(unit,*) nr+1, 4*nang-3 
        write(unit,BV_FORMAT) r(1:nr+1)
        write(unit,BV_FORMAT) theta_meridian_2PI(1:4*nang-3) 
      else
        open(unit=unit, file="out/2D/b_merid_l0180_volume.dat", access="append")
      end if
       do i = 1, nr+1
       do j = 1, 4*nang-3
        write(unit,BVOL_FORMAT) br_merd(i,j,1),bth_merd(i,j,1),bphi_merd(i,j,1),phisc_merd(i,j,1),psisc_merd(i,j,1)
       end do 
       end do 
      close(unit)

    ! 2D meridional Profile: longitudinals 90 + 270 degrees 
      unit = get_free_unit()

      if (timeMyr == 0.0d0) then
        open(unit=unit, file="out/2D/b_merid_l90270_volume.dat")
        write(unit,*) nr+1, 4*nang-3 
        write(unit,BV_FORMAT) r(1:nr+1)
        write(unit,BV_FORMAT) theta_meridian_2PI(1:4*nang-3) 
      else
        open(unit=unit, file="out/2D/b_merid_l90270_volume.dat", access="append")
      end if
       do i = 1, nr+1
       do j = 1, 4*nang-3
        write(unit,BVOL_FORMAT) br_merd(i,j,2), bth_merd(i,j,2), bphi_merd(i,j,2), phisc_merd(i,j,2), psisc_merd(i,j,2)
       end do 
       end do 
      close(unit)

    ! 2D equatorial profile
      unit = get_free_unit()

      if (timeMyr == 0.0d0) then
        open(unit=unit, file="out/2D/b_equator_volume.dat")
        write(unit, *) nr+1, 4*nang-3 
        write(unit,BV_FORMAT) r(1:nr+1)
        write(unit,BV_FORMAT) phi_equator(1:4*nang-3)
      else
        open(unit=unit, file="out/2D/b_equator_volume.dat", access="append")
      end if
       do i = 1, nr+1
       do j = 1, 4*nang-3
        write(unit,BVOL_FORMAT) br_equator(i,j), bth_equator(i,j), bphi_equator(i,j), phisc_eq(i,j), psisc_eq(i,j) 
       end do 
       end do 
      close(unit)

      ! 2D magnetic map at the surface
       unit = get_free_unit()

      if (timeMyr == 0.0d0) then
        open(unit=unit, file="out/2D/Bmap_surf.dat", action="write")
      else
        open(unit=unit, file="out/2D/Bmap_surf.dat", status="old", access="append")
      end if

      ! Each snapshot starts with "snapshot [time in yr].
      write(unit, *) "snapshot", timeMyr

      do p = 1, 6
      do k = 1, nang
      do j = 1, nang
        write(unit, BMAP_FORMAT) &
            & theta(j,k,p), phi(j,k,p), br(nr,j,k,p), &
            & bth(nr,j,k,p), bphi(nr,j,k,p), br(3*nr/4,j,k,p), &
            & bth(3*nr/4,j,k,p), bphi(3*nr/4,j,k,p), br(nr/2,j,k,p), &
            & bth(nr/2,j,k,p), bphi(nr/2,j,k,p), br(nr/4,j,k,p), &
            & bth(nr/4,j,k,p), bphi(nr/4,j,k,p), br(5,j,k,p), &
            & bth(5,j,k,p), bphi(5,j,k,p)
      end do 
      end do 
      end do 

      ! Each snapshot ends with "snapshot_end".
      write(unit, *) "snapshot_end"
      close(unit)

    end subroutine output_magnetic_2D



    !---------------------------------------------------------------------------
    !> Subroutine to write the thermal 1D profile
    !
    !>authors 
    ! Stefano Ascenzi
    !---------------------------------------------------------------------------
    subroutine output_thermal_1D(timeMyr)

      implicit none 

      real*8, intent(in) :: timeMyr
      logical :: isthermalgrid = .true.
      character(40) :: title
      character(40) :: label
      real*8, dimension(1:2*nangt+1,4) :: f_mer_surf
      real*8, dimension(1:4*nangt) :: f_azi_surf

      ! Make 1D profile from 3D data 
      call all_profiles(temp, "temp", timeMyr, isthermalgrid)
      call all_profiles(kappa_perp_arr, "kappa_perp", timeMyr, isthermalgrid)
      call all_profiles(omegatau_arr, "omegatau", timeMyr, isthermalgrid)
      call all_profiles(cv, "cv", timeMyr, isthermalgrid)
      call all_profiles(q_neutrino/UNIT_TIME, "qnu", timeMyr, isthermalgrid)
      call all_profiles(q_joule/UNIT_TIME, "qj", timeMyr, isthermalgrid)

      label = trim("cv_core")
      title = "out/1D/"//trim(label)
      call output_1d_ygraph(trim(title)//".yg", trim(label), &
      &      timeMyr, rtot(1:ncore), cv_core(1:ncore), ncore, YG_FORMAT)
      label = trim("qnu_core")
      title = "out/1D/"//trim(label)
      call output_1d_ygraph(trim(title)//".yg", trim(label), &
      &      timeMyr, rtot(1:ncore), q_neutrino_core(1:ncore)/UNIT_TIME, ncore, YG_FORMAT)
 

      call get_1d_cuts_thermal(temp_surf(:,:,:),(nangt+1)/2, f_mer_surf,f_azi_surf)

      ! Output of the meridional cuts
      ! Common titles and labels
      label = trim("temp_surf_merid_l")
      title = "out/1D/"//trim(label)

      
      !thermal grid (maxID should be optional and I don't specify it)
      call output_1d_ygraph(trim(title)//"0.yg", trim(label)//"0", &
     &      timeMyr, theta_meridian(1:2*nang-1:2), f_mer_surf(:,1), & 
     &      2*nangt+1, YG_FORMAT)
      call output_1d_ygraph(trim(title)//"90.yg", trim(label)//"90", &
     &      timeMyr, theta_meridian(1:2*nang-1:2), f_mer_surf(:,2), &
     &      2*nangt+1, YG_FORMAT)
      call output_1d_ygraph(trim(title)//"180.yg", trim(label)//"180", &
     &      timeMyr, theta_meridian(1:2*nang-1:2), f_mer_surf(:,3), &
     &      2*nangt+1, YG_FORMAT)
      call output_1d_ygraph(trim(title)//"270.yg", trim(label)//"270", &
     &      timeMyr, theta_meridian(1:2*nang-1:2), f_mer_surf(:,4), &
     &      2*nangt+1, YG_FORMAT)

      ! Equatorial profile
      label = trim("temp_surf_equator")
      title = "out/1D/"//trim(label)//".yg"
      call output_1d_ygraph(trim(title), trim(label), &
    &             timeMyr, phi_equator(2:4*nang-4:2), f_azi_surf, &
    &             4*nangt, YG_FORMAT)

    end subroutine output_thermal_1D

    !----------------------------------------------------------------------------
    !> Subroutine to write the thermal 2D meridional and equatiorial cuts
    !
    !>authors 
    ! Stefano Ascenzi
    !----------------------------------------------------------------------------

    !subroutine output_thermal_2D(timeMyr)
    !
    !  implicit none 
    !
    !  real*8, intent(in) :: timeMyr 
    !
    !  call get_2d_cuts(temp,(nang+1)/2,temp_merd,temp_equator)  
    !
    !end subroutine output_thermal_2D


   
    !---------------------------------------------------------------------------
    !> Module to write the desired 1D profiles
    !
    !> authors
    ! Daniele Viganò
    ! Stefano Ascenzi
    ! Clara Dehman 
    ! 
    !---------------------------------------------------------------------------
    subroutine all_profiles(fin, namef, timeMyr, isthermalgrid)

      implicit None
      real*8, intent(in) :: timeMyr
      logical, intent(in) :: isthermalgrid

      !real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6), intent(in) :: fin
      !real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: f
      real*8, dimension(:,:,:,:), intent(in) :: fin
      real*8, dimension(:,:,:,:), allocatable :: f
      character(*), intent(in) :: namef
      real*8, parameter :: EPS=1d-20
      !real*8, dimension(1:2*nang-1,4) :: f_mer
      !real*8, dimension(1:4*nang-3) :: f_azi
      real*8, dimension(:,:), allocatable :: f_mer
      real*8, dimension(:), allocatable :: f_azi
      integer :: ie, rr, pp
      integer :: ang_dim
      character(50) :: title
      character(50) :: label
      integer, dimension(5) :: ir
      integer, dimension(3) :: ip        
      integer, dimension(3) :: ipos       
      integer :: j, jc !debugging

      ! ANGULAR PROFILES
      if (isthermalgrid) then 
        ! thermal grid
        ang_dim = nangt
        ir(1) = 2
        ir(2) = 3
        ir(3) = nrt/2
        ir(4) = nrt-2
        ir(5) = nrt
        allocate(f(1:nrt,0:nangt+1,0:nangt+1,1:6))
        allocate(f_mer(1:2*nangt+1,4))
        allocate(f_azi(1:4*nangt+1))
      else 

        ! magnetic grid
        ang_dim = nang
  
        ir(1) = 6
        ir(2) = nr/2
        ir(3) = nr-1
        ir(4) = nr
        ir(5) = nr+1

        allocate(f(0:nr+1,0:nang+1,0:nang+1,1:6))
        allocate(f_mer(1:2*nang-1,4))
        allocate(f_azi(1:4*nang-3))
      endif

      ie = (ang_dim+1)/2

      !allocate(fin(0:nr+1,0:nang+1,0:nang+1,1:6))

      ! Get rid of very small values that can create visualization issues
      f = fin*MERGE(0, 1, abs(fin) .lt. EPS)
      
      do rr = 1, 5
        if (isthermalgrid) then
          call get_1d_cuts_thermal(f(ir(rr),:,:,:),ie,f_mer,f_azi)
        else
          call get_1d_cuts(f(ir(rr),:,:,:),ie,f_mer,f_azi)
        endif

        ! Output of the meridional cuts
        ! Common titles and labels
        label = trim(namef//"_merid_r"//int_to_string(ir(rr))//"_l")
        title = "out/1D/"//trim(label)
        
        if (isthermalgrid) then 
          ! thermal grid (maxID should be optional and I don't specify it)
          call output_1d_ygraph(trim(title)//"0.yg", trim(label)//"0", &
         &      timeMyr, theta_meridian(1:2*nang-1:2), f_mer(:,1), & 
         &      2*nangt+1, YG_FORMAT)
          call output_1d_ygraph(trim(title)//"90.yg", trim(label)//"90", &
         &      timeMyr, theta_meridian(1:2*nang-1:2), f_mer(:,2), &
         &      2*nangt+1, YG_FORMAT)
          call output_1d_ygraph(trim(title)//"180.yg", trim(label)//"180", &
         &      timeMyr, theta_meridian(1:2*nang-1:2), f_mer(:,3), &
         &      2*nangt+1, YG_FORMAT)
          call output_1d_ygraph(trim(title)//"270.yg", trim(label)//"270", &
         &      timeMyr, theta_meridian(1:2*nang-1:2), f_mer(:,4), &
         &      2*nangt+1, YG_FORMAT)
        ! Equatorial profile
          label = namef//"_equator_r"//int_to_string(ir(rr))
          title = "out/1D/"//trim(label)//".yg"
          call output_1d_ygraph(trim(title), trim(label), &
      &             timeMyr, phi_equator(2:4*nang-4:2), f_azi, &
      &             4*nangt, YG_FORMAT)
        else
          

          ! magnetic grid
          call output_1d_ygraph(trim(title)//"0.yg", trim(label)//"0", &
         &      timeMyr, theta_meridian(2:2*nang-2), f_mer(2:2*nang-2,1), 2*nang-3, YG_FORMAT)
          call output_1d_ygraph(trim(title)//"90.yg", trim(label)//"90", &
         &      timeMyr, theta_meridian(2:2*nang-2), f_mer(2:2*nang-2,2), 2*nang-3, YG_FORMAT)
          call output_1d_ygraph(trim(title)//"180.yg", trim(label)//"180", &
         &      timeMyr, theta_meridian(2:2*nang-2), f_mer(2:2*nang-2,3), 2*nang-3, YG_FORMAT)
          call output_1d_ygraph(trim(title)//"270.yg", trim(label)//"270", &
         &      timeMyr, theta_meridian(2:2*nang-2), f_mer(2:2*nang-2,4), 2*nang-3, YG_FORMAT)         
          ! Equatorial profile
          label = namef//"_equator_r"//int_to_string(ir(rr))
          title = "out/1D/"//trim(label)//".yg"
          call output_1d_ygraph(trim(title), trim(label), &
          &      timeMyr, phi_equator, f_azi, 4*nang-3, YG_FORMAT)

        endif
      enddo

      ! RADIAL PROFILES
      ! Patches considered in the radial output, arbitrary
      ip(1) = 2
      ip(2) = 5
      ip(3) = 6
      
      ! Positions considered, j=k for simplicity, but can be changed
      ipos(1) = (ang_dim+1)/2
      ipos(2) = ang_dim
      ipos(3) = 1
      
      if (isthermalgrid) then
        ! thermal grid 
        do pp = 1, 3
          do rr = 1, 3
            if (rr == 1) then
              label = "_center"
            elseif (rr == 2) then
              label = "_UpRight"
            elseif (rr == 3) then
              label = "_BotLeft"
            endif
            label = namef//"_radial_patch"//int_to_string(ip(pp))//label
            title = "out/1D/"//trim(label)//".yg"
            call output_1d_ygraph(trim(title), trim(label), &
            &    timeMyr, r(1:nr-1:2), f(1:nrt,ipos(rr),ipos(rr),ip(pp)), &
            &     size(r(1:nr-1:2)), YG_FORMAT)
          enddo
        enddo
      else
        ! magnetic grid
        do pp = 1, 3
          do rr = 1, 3
            if (rr == 1) then
              label = "_center"
            elseif (rr == 2) then
              label = "_UpRight"
            elseif (rr == 3) then
              label = "_BotLeft"
            endif
  
            label = namef//"_radial_patch"//int_to_string(ip(pp))//label
            title = "out/1D/"//trim(label)//".yg"
            call output_1d_ygraph(trim(title), trim(label), &
            &    timeMyr, r(0:nr+1), f(0:nr+1,ipos(rr),ipos(rr),ip(pp)), nr+2, &
            &     YG_FORMAT)
          enddo
        enddo
      endif

    end subroutine all_profiles

    
   !---------------------------------------------------------------------------
    !> Outputs a 1D snapshot of any variable in ygraph format.
    !> @brief
    !> @param[in] file          Output filename.
    !> @param[in] label         Label for the snapshot.
    !> @param[in] timeMyr      Time of the snapshot.
    !> @param[in] positions     Array of positions in 1D axis.
    !> @param[in] values        Array of values for the variable.
    !> @param[in] maxIdx        Maximum index for the positions.
    !> @param[in] format        Format for the output.
    ! 
    !> authors
    !> Jose Pons Botella
    !> Daniele Viganò
    !> Alberto Garcia-Garcia
    !---------------------------------------------------------------------------
  subroutine output_1d_ygraph(file, label, time, positions, values, maxIdx, format)
    
    ! Modules ----------------------------------------------------------------
    ! None.
    
    implicit none
    
    ! Input parameters -------------------------------------------------------
    character(*), intent(in) :: file
    character(*), intent(in) :: label
    real*8, intent(in) :: time
    real*8, intent(in) :: positions(:)
    real*8, intent(in) :: values(:)
    integer, optional, intent(in) :: maxIdx
    character(*), optional, intent(in) :: format
    
    ! Local constants --------------------------------------------------------
    ! None.
    
    ! Local variables --------------------------------------------------------
    integer :: unit
    integer :: i, i_max
    
    ! ------------------------------------------------------------------------

    unit = get_free_unit()
    
    ! If we just started the simulation, create a new file and write header.
    if (time == 0.0d0) then
    open(unit=unit, file=file)
    ! If not, just append to the end of the file.
    else
    open(unit=unit, file=file, access="append")
    end if
    
    write(unit, *) '"Label=', label
    write(unit, *) '"Time=', time*1d6
    
    ! If a maximum number for the positions array is specified, use it, if not
    ! fall back to the actual size of the positions array.
    if (present(maxIdx)) then
    i_max = maxIdx
    else
    i_max = size(positions)
    end if
    
    do i = 1, i_max
    if (present(format)) then
    write(unit, format) positions(i), values(i)
    else
    write(unit, *) positions(i), values(i)
    end if
    end do ! i
    
    write(unit, *)
    
    close(unit)
    
  end subroutine output_1d_ygraph

!-----------------------------------------------------------------------
  !> This subroutine generates checkpoint files
  !
  !> authors 
  ! Stefano Ascenzi
  ! 
  ! NOTE: at the moment the file save only the temperature
!------------------------------------------------------------------------

subroutine output_checkpoint(it, time)
  
  use input_params, only: resume_checkpnt_number
  use grid, only : last_timestep_print

  implicit none

  integer, intent(in) :: it
  real*8, intent(in) :: time
  integer :: u, imin_vis
  integer :: i, j, k, p
  character(len=50) :: file_name, it_str
  logical, save :: firstcall_checkpoint = .TRUE.

!  if (firstcall_checkpoint) then
!    if (resume_checkpnt_number .gt. 0) then 
!      ! if we resume from checkpoint, we start to number 
!      ! the new checkpoints from the one we use to restart
!      timestep_print = resume_checkpnt_number + 1  
!    else
!      timestep_print = 1
!    endif
!    firstcall_checkpoint = .FALSE.   
!  endif

  write(it_str, "(I5)") it

  file_name = "out/3D/checkpoint_"//trim(adjustl(it_str))//".dat"

  print*, '** Writing checkpoint_'//trim(adjustl(it_str))//".dat **"

  imin_vis = 0

  u = get_free_unit()
  open(unit = u, file = file_name, status='replace')

  ! first line is the timestep
  ! second line is the time
  ! third line is the number of last printed timestep
  ! other lines are the temperatures, including ghost cells 

  write(u,*) it
  write(u,*) time 
  write(u,*) last_timestep_print

  do p=1,6
    do i=imin_vis+1,nrt
      do j=0,nangt+1
        do k=0,nangt+1
          write(u,*) temp(i,j,k,p)
        end do
      end do
    end do
  end do

  close(u)

end subroutine output_checkpoint


!---------------------------------------------------------------------------
    !> This subroutine generates the output in vtu format
    !
    !> authors
    ! Albert Herrando
    ! Daniele Viganò
    ! Stefano Ascenzi
    !---------------------------------------------------------------------------
  subroutine output_vtu(time,dt)

    use input_params, only: final_time, resume_checkpnt_number
    use grid, only: last_timestep_print
    
    implicit none
    ! Input parameters -------------------------------------------------------
    real*8, intent(in) :: time, dt
    !logical, intent(in) :: print_T, print_B
    ! ---------------------------------------------
    character(len=70) :: time_str, N_Points_str, N_Cells_str, file_name, file_name_pvd
    character(len=50) :: timestep_str
    integer :: u = 200, upvd = 220
    integer :: N_Points, N_Cells
    integer ::  p, i, j, k, id, imin_vis
    integer, dimension(0:nr,1:nang,1:nang,6) :: id_nodes
    real*8, dimension(0:nr+1, 0:nang+1, 0:nang+1, 1:6) :: fx, fy, fz, fth, fphi     
    logical, save :: FirstCall_output_vtu = .true.
    logical, save :: FirstCall_pvd = .true.
    integer, save :: timestep_print = 0

    if (FirstCall_output_vtu .and. resume_checkpnt_number > 0) then

      ! if we restart from checkpoint we should 1) restart the numbering of output files where 
      ! we stopped in the previous run 2) no need to recreate another pvd file, we just append 
      ! new lines to the one already existing 
      !SA: .pvd file will have a problem anyway because 
      ! the lines written after the checkpoint will stay there. We should find a way to overwrite 
      ! .pvd file starting from the output written after the restarting checkpoint (minor problem) 

      FirstCall_output_vtu = .false. 
      Firstcall_pvd = .false. 

      timestep_print = last_timestep_print + 1 
    endif 


    ! Set this to 0 if you start from the first radial cell
    ! 1 from the second one
    imin_vis = 0
    
    ! Grid cells and points
    N_Points = 6*(nrt+1-imin_vis)*(nangt+1)*(nangt+1)
    N_Cells = 6*(nrt-imin_vis)*nangt*nangt
    write(N_Cells_str, "(I20)") N_Cells
    write(N_Points_str, "(I20)") N_Points
    
    ! Convert (REAL) time  to (STRING) time 
    write(time_str, "(F10.3)") time*1.d6 ! in yrs
    write(timestep_str, "(I3)") timestep_print
   
    
    ! Output files names
    !file_name = "out/3D/output_time_" // trim(adjustl(time_str)) // ".vtu"
    file_name = "out/3D/output_"//trim(adjustl(timestep_str))//".vtu"
    file_name_pvd = "out/3D/animation.pvd"

    !*******************************************
    ! Replacing the . from number to _ to avoid file format problems
    ! THIS IS NEEDED IF WE PRINT THE TIME IN THE FILE NAME
    ! IF WE PRINT THE TIMESTEP IT WILL SUBSTITUTE THE . of .vtu
    ! WITH _, SO BE CAREFULL

    !do i = 1, len(file_name)
    !    if (file_name(i:i) == ".") then
    !        file_name(i:i) = "_"
    !        exit
    !    end if
    !enddo
    ! *********************************************
    
    
    ! Open file
    open(unit = u, file = file_name, status='replace')
    
    ! Write VTU file
    write(u,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(u,'(A)') '  <UnstructuredGrid>'
    write(u,'(A,A,A,A,A)') '    <Piece NumberOfPoints="', trim(adjustl(N_Points_str)), '" NumberOfCells="', & 
                         & trim(adjustl(N_Cells_str)), '">'
    write(u,'(A)') '      <Points>'
    write(u,'(A)') '        <DataArray type="Float32" format="ascii" NumberOfComponents="3">'
    id = 0
    do p=1,6
      do i=imin_vis*2,nr,2
        do j=1,nang,2
          do k=1,nang,2
            id_nodes(i,j,k,p) = id
            write(u,'(f10.5,f10.5,f10.5)') xc(i,j,k,p), yc(i,j,k,p), zc(i,j,k,p)
            id = id + 1
          end do
        end do
      end do
    end do
    write(u,'(A)') '        </DataArray>'
    write(u,'(A)') '      </Points>'
    write(u,'(A)') '      <Cells>'
    write(u,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
    do p=1,6
      do i=imin_vis*2+1,nr-1,2
        do j=2,nang-1,2
          do k=2,nang-1,2
            write(u,*) id_nodes(i-1,j-1,k-1,p), &
                 &                id_nodes(i-1,j+1,k-1,p), &
                 &                id_nodes(i-1,j+1,k+1,p), &
                 &                id_nodes(i-1,j-1,k+1,p), &
                 &                id_nodes(i+1,j-1,k-1,p), &
                 &                id_nodes(i+1,j+1,k-1,p), &
                 &                id_nodes(i+1,j+1,k+1,p), &
                 &                id_nodes(i+1,j-1,k+1,p)
          end do
        end do
      end do
    end do
    write(u,'(A)') '        </DataArray>'
    write(u,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
    do i = 1, N_Cells
      write(u,*) i*8
    enddo
    write(u,'(A)') '        </DataArray>'
    write(u,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
    do i = 1, N_Cells
      write(u,'(A)') '          12'
    enddo
    write(u,'(A)') '        </DataArray>'
    write(u,'(A)') '      </Cells>'
    write(u,'(A)') '      <PointData>'
    write(u,'(A)') '      </PointData>'
    write(u,'(A)') '      <CellData>'
    !if (print_T) then
      write(u,'(A)') '        <DataArray Name="T" NumberOfComponents="1" type="Float32" format="ascii">'
      do p=1,6
        do i=imin_vis+1,nrt
          do j=1,nangt
            do k=1,nangt
              write(u,'(e11.3)') temp(i,j,k,p)
            end do
          end do
        end do
      end do
        write(u,'(A)') '        </DataArray>'
    !else if (print_B) then

      call f_cs_to_spherical(bxi,beta,fth,fphi,ievol-2) !ievol-2 is the center of the sphere
      call f_spherical_to_cartesian(br,fth,fphi,fx,fy,fz,ievol-2)

      write(u,'(A)') '        <DataArray Name="B" NumberOfComponents="3" type="Float32" format="ascii">'
      do p=1,6
        do i=2*imin_vis+1,nr-1,2
          do j=2,nang-1,2
            do k=2,nang-1,2
              write(u,'(e11.3,e11.3,e11.3)') fx(i,j,k,p), fy(i,j,k,p), fz(i,j,k,p)
            end do
          end do
        end do
      end do
      write(u,'(A)') '        </DataArray>'

      call f_cs_to_spherical(jxi,jeta,fth,fphi,ievol-1)
      call f_spherical_to_cartesian(jr,fth,fphi,fx,fy,fz,ievol-1)

      write(u,'(A)') '        <DataArray Name="J" NumberOfComponents="3" type="Float32" format="ascii">'
      do p=1,6
        do i=2*imin_vis+1,nr-1,2
          do j=2,nang-1,2
            do k=2,nang-1,2
              write(u,'(e11.3,e11.3,e11.3)') fx(i,j,k,p), fy(i,j,k,p), fz(i,j,k,p)
            end do
          end do
        end do
      end do
      write(u,'(A)') '        </DataArray>'

    !end if
    write(u,'(A)') '      </CellData>'
    write(u,'(A)') '    </Piece>'
    write(u,'(A)') '  </UnstructuredGrid>'
    write(u,'(A)') '</VTKFile>'
    
    ! Close file
    close(unit=u)

    ! now we write the pvd file PROVA

    if (FirstCall_pvd) then
      open(unit = upvd, file=file_name_pvd, status = 'replace')

      write(upvd, '(A)') '<?xml version="1.0" ?>'
      write(upvd, '(A)') '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
      write(upvd, '(A)') '<Collection>'
      !write(upvd, '(A)') '<DataSet timestep="'//time_str//'"file='//file_name//'/> '
      write(upvd, '(A)') '<DataSet timestep="'//trim(adjustl(time_str)) &
      & //'" group="" part="0" file="output_'//trim(adjustl(timestep_str))//'.vtu"/>'
    else
      open(unit = upvd, file=file_name_pvd, access = 'append')

      !write(upvd, '(A)') '<DataSet timestep="'//time_str//'"file='//file_name//'/> '
      write(upvd, '(A)') '<DataSet timestep="'//trim(adjustl(time_str)) &
      & //'" group="" part="0" file="output_'//trim(adjustl(timestep_str))//'.vtu"/>'
      !write(upvd, '(A)') '<DataSet timestep="output_" // trim(adjustl(timestep_str)) // ".vtu"/'
    endif 

    if (time + dt >= final_time) then
      write(upvd, '(A)') '</Collection>'
      write(upvd, '(A)') '</VTKFile>'
    endif

    close(unit = upvd)

    FirstCall_pvd = .false.
    last_timestep_print = timestep_print
    timestep_print = timestep_print + 1 

  end subroutine

  !---------------------------------------------------------------------------
    !> This subroutine generates the cooling curve output
    !
    !> authors
    ! Stefano Ascenzi
    !---------------------------------------------------------------------------

  subroutine output_temperature_cooling_curve(iter, timeMyr, dt)
    
    implicit none 

    real*8, intent(in) :: timeMyr,dt
    integer, intent(in) :: iter

    logical, save :: Firstcall_out_cooling = .true.
    logical :: Firstcall_cfl = .true.
    real*8 :: lum_em, teff, eint, qnu_crust_tot
    real*8 :: cfl, cfl_serv !courant condition
    real*8, save :: erad = 0.d0, enu_rad = 0.d0
    integer :: j, k, p, i, jc, kc, ic
    integer :: u1, u2, u3, u4, u5 

    character(len=15), parameter :: COOLING_FORMAT = "(i6,13es12.4)"

    lum_em = 0.d0
    qnu_crust_tot = 0.d0
    eint = 0.d0

    do p = 1, 6
      do k = 1, nangt
        do j = 1, nangt
          do i = 1, nrt

            ic = 2*i-1
            jc = 2*j
            kc = 2*k

            qnu_crust_tot = qnu_crust_tot + vol(ic,jc,kc)*q_neutrino(i,j,k,p)*enu(ic)**2

          enddo
        enddo 
      enddo
    enddo

    qnu_crust_tot = qnu_crust_tot*(UNIT_EN/UNIT_TIME) ! convert in cgs units

    lum_em = sum(bb_flux(1:nangt,1:nangt,1:6))*(UNIT_EN/UNIT_TIME) ! convert in cgs units
    teff = (lum_em/(4d0*PI*r(nr)**2*UNIT_R**2*STEFAN_BOLTZMANN))**(1./4.)


    if (Firstcall_out_cooling) then
      u1 = get_free_unit()
      open(unit = u1, file='out/energy/cooling_curve.d', status = 'replace')
      write(u1, "(a95)") "iter, time[yr], BB lum[erg/s], T_core[K], Teff[K], lum.neutrinos[erg/s] (core,crust), Bpdip[G]"
    else
      u1 = get_free_unit()
      open(unit = u1, file='out/energy/cooling_curve.d', access = 'append')
    endif

    write(u1, COOLING_FORMAT) iter, timeMyr*1d6, lum_em, T_core, teff, qnu_core_tot*(UNIT_EN/UNIT_TIME), qnu_crust_tot, bpdip*1e12
    close(unit=u1)


    if (Firstcall_out_cooling) then
      u2 = get_free_unit()
      open(unit = u2, file='out/energy/temperatures.d', status = 'replace')
      write(u2, *) "timeMyr, Tb centers patches, Ts centers patches"
      !Firstcall_out_cooling = .false.
    else
      u2 = get_free_unit()
      open(unit = u2, file='out/energy/temperatures.d', access = 'append')
    endif

    write(u2, COOLING_FORMAT) iter, timeMyr*1d6, temp(nrt,nangt/2,nangt/2,1:6)*UNIT_T, &
   &   temp_surf(nangt/2,nangt/2,1:6)
    close(unit=u2)

    if (Firstcall_out_cooling) then
      u3 = get_free_unit()
      open(unit=u3, file = 'out/energy/core.d', status = 'replace')
      write(u3, *), "time(yr), T_core[K], cv_core [cgs], qnu_core [cgs]"
    !  Firstcall_out_cooling = .false.
    else
      u3 = get_free_unit()
      open(unit=u3, file = 'out/energy/core.d', access = 'append')
    endif
    

    write(u3, "(14es12.4)") timeMyr*1d6, T_core*1d8, cv_core_tot*UNIT_EN/UNIT_T, qnu_core_tot*UNIT_EN/UNIT_TIME
    close(unit=u3)

    if (Firstcall_out_cooling) then
      u4 = get_free_unit()
      open(unit=u4, file = 'out/energy/neutrino_core.d', status = 'replace')
      write(u4, *), "time(yr), mur [cgs], nn, np, pp, ep, cp_con, cp_cop, du, ea, pl, syn, cp_cr, pa"
    !  Firstcall_out_cooling = .false.
    else
      u4 = get_free_unit()
      open(unit=u4, file = 'out/energy/neutrino_core.d', access = 'append')
    endif

    ! convert in cgs 

    qnu_mur = qnu_mur*UNIT_EN/UNIT_TIME
    qnu_nn = qnu_nn*UNIT_EN/UNIT_TIME
    qnu_np = qnu_np*UNIT_EN/UNIT_TIME
    qnu_pp = qnu_pp*UNIT_EN/UNIT_TIME
    qnu_ep = qnu_ep*UNIT_EN/UNIT_TIME
    qnu_cp_con = qnu_cp_con*UNIT_EN/UNIT_TIME
    qnu_cp_cop = qnu_cp_cop*UNIT_EN/UNIT_TIME
    qnu_du = qnu_du*UNIT_EN/UNIT_TIME
    qnu_ea = qnu_ea*UNIT_EN/UNIT_TIME
    qnu_pl = qnu_pl*UNIT_EN/UNIT_TIME
    qnu_syn = qnu_syn*UNIT_EN/UNIT_TIME
    qnu_cp_cr = qnu_cp_cr*UNIT_EN/UNIT_TIME
    qnu_pa = qnu_pa*UNIT_EN/UNIT_TIME

    write(u4, "(14es12.4)") timeMyr*1d6, qnu_mur(1), qnu_nn(1), qnu_np(1), &
    & qnu_pp(1), qnu_ep(1), qnu_cp_con(1), qnu_cp_cop(1), qnu_du(1), qnu_ea(1), qnu_pl(1), qnu_syn(1), & 
    & qnu_cp_cr(1), qnu_pa(1)
    close(unit=u4)

    if (Firstcall_out_cooling) then
      u5 = get_free_unit()
      open(unit=u5, file = 'out/energy/neutrino_crust.d', status = 'replace')
      write(u4, *), "time(yr), mur [cgs], nn, np, pp, ep, cp_con, cp_cop, du, ea, pl, syn, cp_cr, pa"
      Firstcall_out_cooling = .false.
    else
      u5 = get_free_unit()
      open(unit=u5, file = 'out/energy/neutrino_crust.d', access = 'append')
    endif

    write(u5, "(14es12.4)") timeMyr*1d6, qnu_mur(2), qnu_nn(2), qnu_np(2), &
    & qnu_pp(2), qnu_ep(2), qnu_cp_con(2), qnu_cp_cop(2), qnu_du(2), qnu_ea(2), qnu_pl(2), qnu_syn(2), & 
    & qnu_cp_cr(2), qnu_pa(2)
    close(unit=u5)


    write(6,*) "Cooling:iter,Time(yr),Lum_BB,T_core, Tb(N),Ts(N),Teff,Lum.nu core,Lum.nu crust"
    write(6,"(i5,14es12.4)") iter,timeMyr*1d6, lum_em, T_core*UNIT_T, temp(nrt,nangt/2,nangt/2,5)*UNIT_T, &
   &  temp_surf(nangt/2,nangt/2,5), teff, qnu_core_tot*(UNIT_EN/UNIT_TIME), qnu_crust_tot

  end subroutine

!---------------------------------------------------------------------------
    !> This subroutine generates 2D output of the surface temperature
    !> temperature will be printed in separate files for each timestep with 3 columns for
    !> theta, phi, Temp_surf
    !
    !> authors
    ! Stefano Ascenzi
    !---------------------------------------------------------------------------
  subroutine output_surface(time)

    use grid, only: theta, phi

    implicit none 

    real*8, intent(in) :: time

    character(len=50) :: file_name, time_str
    integer :: u = 330
    integer :: j,k,p, jc, kc

    write(time_str, "(F10.3)") time*1.d6

    file_name = "out/2D/tmap_"//trim(adjustl(time_str))//".dat"


    ! Open file
    open(unit = u, file = file_name, status='replace')

    write(u,*) "theta", "phi", "Tmap"

    do p=1,6
      do k=1,nangt
        kc = 2*k
        do j=1,nangt
          jc = 2*j
          write(u,*) theta(jc,kc,p), phi(jc,kc,p), temp_surf(j,k,p)
        end do
      end do
    end do

    close(u)


  end subroutine

end module output 
!-------------------------------------------------------------------------------
! Magneto Thermal 2D
!-------------------------------------------------------------------------------
! Module: Output
!
!> @author
!> Jose Pons Botella
!> Daniele Viganò
!> Alberto Garcia-Garcia
!
!> Output handling.
!> @brief This module is responsible for handling all output for the simulation,
!>        e.g., creating the output files, writing snapshots of the grid for the
!>        magnetic field evolution and for the temperature evolution.
!
!-------------------------------------------------------------------------------
module output

  use constants, only: UNIT_T
  use utils, only: get_free_unit, real_to_string, int_to_string
  use grid, only : aphi, br, bth, bphi, bm, er, eth, ephi, jr, jth, jphi
  use grid, only: nang, np, kmax, lmax, rb, theta, cth
  use grid, only : fh, tem, etab, q_joule_average, q_neutrino, dtb_courant_profile, sfluxb, tss

  contains

    !---------------------------------------------------------------------------
    !> Outputs a resume snapshot.
    !> @brief The minimum required information about the current year, the
    !>        magnetic field, and temperature map is written to a file
    !>        resume.dat so that execution can be resumed from that checkpoint.
    !
    !> @param[in]   timeYear      Timestamp [yr].
    !> @param[in]   period        Period [s].
    !> @param[in]   periodDot     Period derivative (dimensionless).
    !> @param[in]   bpdip         Dipolar poloidal component at polar surface [10**12 G].
    !---------------------------------------------------------------------------
    subroutine output_resume(timeYear, period, periodDot, bpdip)

      ! Modules ----------------------------------------------------------------
      use magnetic_analysis, only: en_joule_star, en_joule_shock_star, poynting_star_tot
      use grid, only: jevol, jmin

      implicit none

      ! Input parameters -------------------------------------------------------
      real*8, intent(in) :: timeYear
      real*8, intent(in) :: period
      real*8, intent(in) :: periodDot
      real*8, intent(in) :: bpdip

      ! Local constants --------------------------------------------------------

      ! Local variables --------------------------------------------------------
      integer :: unit
      integer :: i, j

      ! ------------------------------------------------------------------------

      unit = get_free_unit()

      open(unit=unit, file="resume.dat", action="write")

      write(unit, *) timeYear, period, periodDot, bpdip

      write(unit, *) jevol, jmin

      write(unit, *) en_joule_star, en_joule_shock_star, poynting_star_tot

      do i = 0, nang +1
        do j = 0, np + 2
          write(unit, *) br(i, j), bth(i, j), bphi(i, j), aphi(i, j)
        end do !j
      end do !i

      do i = 1, kmax
        do j = 1, lmax
          write(unit, *) tem(i, j)
        end do !j
      end do !i

      close(unit)

    end subroutine output_resume

    !---------------------------------------------------------------------------
    !> Outputs a complete snapshot of the temperature simulation.
    !> @brief
    !
    !> @param[in]   timeYear      Timestamp [yr].
    !> @param[in]   bindex        Braking index [dimensionless].
    !> @param[in]   bpdip         Surface value of the dipolar component of B [10**12 G].
    !> @param[in]   itert         Counter (number of iterations written in output)
    !> @param[in]   per           Period [s].
    !> @param[in]   pdot          Period derivative (dimensionless).
    !---------------------------------------------------------------------------
    subroutine output_temperature(timeYear, bindex, bpdip, itert, per, pdot)

      ! Modules ----------------------------------------------------------------
      use grid, only: gammac, jcore, shearMaximum, shearModulus

      implicit none

      ! Input parameters -------------------------------------------------------
      integer, intent(in) :: itert
      real*8, intent(in) :: timeYear
      real*8, intent(in) :: bindex
      real*8, intent(in) :: bpdip
      real*8, intent(in) :: per
      real*8, intent(in) :: pdot

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
      ! Meridional index at which we take the slice to show the radial profile.
      integer :: i_slice
      ! Radial index at which we take the slice to show the meridional profile.
      integer :: j_slice
      real*8, dimension(1:np) :: omegatau

      ! ------------------------------------------------------------------------
      i_slice = ( nang + 1) / 3
      j_slice = np - 2

      ! Output cooling curve ---------------------------------------------------
      call output_temperature_cooling_curve(timeYear, itert, bindex, bpdip, per, pdot)

      ! Output resume information.
      call output_resume(timeYear, per, pdot, bpdip)

      ! Output map -------------------------------------------------------------
      call output_temperature_map(timeYear)

      ! Output 1D graphs -------------------------------------------------------

      ! Output the Hall precoefficient for the first time step. It only depends
      ! on the electron density, which does not evolve in our current model.
      if (timeYear == 0.0d0) then
        call output_1d_ygraph("out/fh.yg", &
                            & "fh", &
                            & timeYear, rb(1:), fh(1:), np)
      end if

      ! Magnetic diffusivity, meridional profile.
      call output_1d_ygraph("out/eta_meridional.yg", &
                          & "eta(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(1:), etab(1:, j_slice), nang)
      ! Magnetic diffusivity, radial profile.
      call output_1d_ygraph("out/eta_radial.yg", &
                          & "eta(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:), etab(i_slice, 1:), np)
      ! Magnetization parameter, radial profile.
      omegatau(1:np)=fh(1:np)*bm(i_slice,1:np)/maxval([1d-20,etab(i_slice,1:np)])
      call output_1d_ygraph("out/omegatau_radial.yg", &
                          & "omtau(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:np), omegatau(1:np), np)
      ! Temperature at the pole, radial profile.
      call output_1d_ygraph("out/Tp_radial.yg", &
                          & "Tpole", &
                          & timeYear, rb(3::2), 1.d8 * tem(2, 2:), LMAX-1)
      ! Temperature at the equator, radial profile.
      call output_1d_ygraph("out/Te_radial.yg", &
                          & "Tequa", &
                          & timeYear, rb(3::2), 1.d8 * tem(KMAX/2, 2:), LMAX-1)
      ! Temperature at the bottom, meridional profile.
      call output_1d_ygraph("out/Tb_meridional.yg", &
                          & "Tb", &
                          & timeYear, theta(2::2), &
                          & dlog10(UNIT_T * tem(2:, LMAX)), KMAX-1)
      ! Temperature at the surface, meridional profile.
      call output_1d_ygraph("out/Ts_meridional.yg", &
                          & "Tsur", &
                          & timeYear, theta(2::2), dlog10(tss(2:)), KMAX-1)

      ! Joule dissipated energy, radial profile.
      call output_1d_ygraph("out/qj_radial.yg", &
                          & "qj(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(3::2), q_joule_average((i_slice+1)/2, 2:), LMAX-1)
      ! Joule dissipated energy, meridional profile.
      call output_1d_ygraph("out/qj_meridional.yg", &
                          & "qj(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(2::2), q_joule_average(2:, j_slice/2), KMAX-1)
      ! Local neutrino emissivity, radial profile.
      call output_1d_ygraph("out/qnu_radial.yg", &
                          & "qnu(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(3::2), q_neutrino((i_slice+1)/2, 2:), LMAX-1)
      ! Local neutrino emissivity, meridional profile.
      call output_1d_ygraph("out/qnu_meridional.yg", &
                          & "qnu(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(2::2), q_neutrino(2:, j_slice/2), KMAX-1)

      ! Coulomb parameter, meridional Profile
      call output_1d_ygraph("out/gamma_meridional.yg", &
                           & "Gamma(" // real_to_string(rb(j_slice)) // "th)", &
                           &  timeYear,theta(2::2),gammac(2:kmax,j_slice/2), &
                           &  kmax-1)
      ! Coulomb parameter, radial Profile
      call output_1d_ygraph("out/gamma_equator.yg", &
                           & "Gamma(r,equator)" , &
                           &  timeYear,rb(jcore-1::2), gammac(kmax/2,jcore/2:lmax), &
                           &  lmax-jcore/2+1)

      ! Shear Maximum Meridional Profile
      call output_1d_ygraph("out/shearmaximum_meridional.yg", &
                          & "Shear Maximum(" // real_to_string(rb(j_slice)) // "th)", &
                          &  timeYear,theta(2::2),shearMaximum(2:kmax,j_slice/2), &
                          &  kmax-1)
      ! Shear Maximum Radial Profile
      call output_1d_ygraph("out/shearmaximum_equator.yg", &
                          & "Shear Maximum(r,equator)" , &
                          &  timeYear,rb(jcore-1::2), shearMaximum(kmax/2,jcore/2:lmax), &
                          &  lmax-jcore/2+1)

      !Shear Modulus Meridional Profile
      call output_1d_ygraph("out/shearmodulus_meridional.yg", &
                          & "Shear Modulus(" // real_to_string(rb(j_slice)) // "th)", &
                          &  timeYear,theta(2::2),shearModulus(2:kmax,j_slice/2), &
                          &  kmax-1)
      !Shear Modulus Radial Profile
       call output_1d_ygraph("out/shearmodulus_equator.yg", &
                          & "Shear Modulus(r,equator)" , &
                          &  timeYear,rb(jcore-1::2), shearModulus(kmax/2,jcore/2:lmax), &
                          &  lmax-jcore/2+1)

    end subroutine output_temperature


    !---------------------------------------------------------------------------
    !> Outputs a snapshot of the current cooling curve.
    !> @brief
    !
    !> @param[in]   timeYear      Timestamp in [yr].
    !> @param[in]   itert         Counter (number of iterations in cooling)
    !> @param[in]   bindex        Braking index [dimensionless]
    !---------------------------------------------------------------------------
    subroutine output_temperature_cooling_curve(timeYear, itert, bindex, bpdip, per, pdot)

      ! Modules ----------------------------------------------------------------
      use constants, only: &
        & PI, UNIT_B, UNIT_R, PHFLUX_CONSTANT, STEFAN_BOLTZMANN, T_YEAR

      implicit none

      ! Input parameters -------------------------------------------------------
      real*8, intent(in) :: timeYear
      integer, intent(in) :: itert
      real*8, intent(in) :: bindex
      real*8, intent(in) :: bpdip
      real*8, intent(in) :: per
      real*8, intent(in) :: pdot

      ! Local constants --------------------------------------------------------
      character(len=17), parameter :: COOLING_FORMAT = "(13(1x,es15.7))"
      character(len=72), parameter :: COOLING_CONSOLE_FORMAT = &
        & "(a7, i7, a3, es10.3, 2(a13, 2es10.3), a6, es10.3, a9, es10.3)"

      ! Local variables --------------------------------------------------------
      real*8 :: phflux(KMAX)
      real*8 :: lumerg
      real*8 :: lumph
      real*8 :: teffi
      real*8 :: sd_age
      integer :: unit
      integer :: i

      ! ------------------------------------------------------------------------

      sd_age = 0.d0
      if (pdot /= 0.0d0) then
        sd_age = per / (pdot * T_YEAR)
      end if

      phflux = PHFLUX_CONSTANT * UNIT_R**2 * tss**3 ! photon flux in [ph/km^2/s].
      lumph = 0.d0
      lumerg  = 0.d0

      do i = 2, KMAX
        lumerg = lumerg + sfluxb(i) * &
                & 2.d0 * pi * rb(np)**2 * (cth(2*i-3) - cth(2*i-1))
        lumph = lumph + phflux(i) * &
                & 2.d0 * pi * rb(np)**2 * (cth(2*i-3) - cth(2*i-1))
      end do ! i

      lumerg = lumerg * 1.d40
      teffi = (lumerg / &
              & (STEFAN_BOLTZMANN * 4.d0 * PI * (UNIT_R * rb(np))**2))**(0.25d0)

      unit = get_free_unit()

      if (timeYear == 0.0d0) then
        open(unit=unit, file="out/cool_curve.d")
      else
        open(unit=unit, file="out/cool_curve.d", access="append")
      end if

      write(unit, COOLING_FORMAT) &
        & timeYear, &
        & lumerg, teffi, bpdip * UNIT_B, &
        & per, pdot, sd_age, &
        & tem(2, LMAX) * UNIT_T, &
        & tem(KMAX / 2 + 1, LMAX) * UNIT_T, &
        & tss(2), tss(KMAX / 2 + 1), bindex, lumph

      close(unit)

      write(*, COOLING_CONSOLE_FORMAT) &
        & "COOLING", itert, "t=", timeYear, &
        & "--Tp, Teq[K]=", tem(2, LMAX) * UNIT_T, &
        & tem(KMAX / 2 + 1,LMAX) * UNIT_T, &
        & "--Tsp, Tse[K]=", tss(2), tss(KMAX / 2 + 1), &
        & " Bpdip [G]=", bpdip * UNIT_B, &
        & " Lerg=", lumerg

    end subroutine output_temperature_cooling_curve

    !---------------------------------------------------------------------------
    !> Outputs a snapshot of the provided grid's temperature.
    !> @brief Given a timestamp and the current temperature values of the grid,
    !> this subroutine outputs a snapshot with the polar coordinates for each
    !> cell in the grid with their associated temperature.
    !
    !> @param[in]   timeYear      Timestamp in [yr].
    !> @param[in]   temperature   Temperature values for the grid cells [K].
    !---------------------------------------------------------------------------
    subroutine output_temperature_map(timeYear)

      ! Modules ----------------------------------------------------------------
      use grid, only : shearModulus, shearMaximum

      implicit none

      ! Input parameters -------------------------------------------------------
      real*8, intent(in) :: timeYear

      ! Local constants --------------------------------------------------------
      character(len=19), parameter :: OUTPUT_FORMAT = "(2f10.3,3e10.3)"
      character(len=32), parameter :: TEMPERATURE_FILE = "out/Tmap.dat"

      ! Local variables --------------------------------------------------------
      integer :: unit
      integer :: i, j

      ! ------------------------------------------------------------------------

      unit = get_free_unit()

      ! If we just started the simulation, create a new file and write header.
      if (timeYear == 0.0d0) then
        open(unit=unit, file=TEMPERATURE_FILE, action="write")
        write(unit, *) 'label Temperature'
      ! If not, just append to the end of the file.
      else
        open(unit=unit, file=TEMPERATURE_FILE, status="old", access="append")
      end if

      ! Each snapshot starts with "snapshot [time in yr]".
      write(unit, *) "snapshot", timeYear

      ! Dump each temperature value for each cell in polar coordinates.
      do i = 2, KMAX
        do j = 2, LMAX
          write(unit, OUTPUT_FORMAT) theta(2*i-2), rb(2*j-1), & 
            & 1.d8*tem(i, j), shearModulus(i, j), shearMaximum(i, j)
        end do ! i
      end do ! j

      ! Each snapshot ends with "snapshot_end".
      write(unit, *) "snapshot_end"
      close(unit)

    end subroutine output_temperature_map



    !---------------------------------------------------------------------------
    !> Outputs a complete snapshot of the magnetic field simulation.
    !> @brief
    !
    !> @param[in]   timeYear      Timestamp in [yr].
    !---------------------------------------------------------------------------
    subroutine output_magnetic(iterb,timeYear)

      ! Modules ----------------------------------------------------------------
      use legpol, only: nleg, blout

      implicit none

      ! Input parameters -------------------------------------------------------
      real*8, intent(in) :: timeYear
      integer, intent(in) :: iterb

      ! Local constants --------------------------------------------------------
      character(len=10), parameter :: MAGNETIC_DATA_FORMAT = "(2es12.4)"
      character(len=17), parameter :: BL_FORMAT = "(f3.0,11es15.6)"
      ! Radial index at which we take the slice to show the meridional profile.
      integer :: j_slice
      ! Meridional index at which we take the slice to show the radial profile.
      integer :: i_slice
      ! Loop index for bl_out
      integer :: l

      ! ------------------------------------------------------------------------

      j_slice = np - 2
      i_slice = (nang + 1) / 3

 
      ! Magnetic field map -----------------------------------------------------
      call output_magnetic_map(timeYear)

      ! Magnetic field analysis ------------------------------------------------      
      call output_magnetic_analysis(iterb,timeYear)

      ! Meridional profiles ----------------------------------------------------
      call output_1d_ygraph("outb/a_phi_surface.yg", &
                          & "Aph(" // real_to_string(rb(np)) // ",th)", &
                          & timeYear, theta(1:), aphi(1:, np), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/a_phi_meridional.yg", &
                          & "Aph(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(1:), aphi(1:, j_slice), nang, &
                          & MAGNETIC_DATA_FORMAT)
      ! Magnetic field, meridional profile.
      call output_1d_ygraph("outb/b_r_surface.yg", &
                          & "Br(" // real_to_string(rb(np)) // ",th)", &
                          & timeYear, theta(1:), br(1:, np), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/b_r_meridional.yg", &
                          & "Br(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(1:), br(1:, j_slice), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/b_th_surface.yg", &
                          & "Bth(" // real_to_string(rb(np)) // ",th)", &
                          & timeYear, theta(1:), bth(1:, np), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/b_th_meridional.yg", &
                          & "Bth(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(1:), bth(1:, j_slice), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/b_phi_meridional.yg", &
                          & "Bph(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(1:), bphi(1:, j_slice), nang, &
                          & MAGNETIC_DATA_FORMAT)
      ! Electrical current, meridional profile.
      call output_1d_ygraph("outb/j_r_surface.yg", &
                          & "Jr(" // real_to_string(rb(np)) // ",th)", &
                          & timeYear, theta(1:), jr(1:, np), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/j_r_meridional.yg", &
                          & "Jr(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(1:), jr(1:, j_slice), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/j_th_surface.yg", &
                          & "Jth(" // real_to_string(rb(np)) // ",th)", &
                          & timeYear, theta(1:), jth(1:, np), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/j_th_meridional.yg", &
                          & "Jth(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(1:), jth(1:, j_slice), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/j_phi_surface.yg", &
                          & "Jph(" // real_to_string(rb(np)) // ",th)", &
                          & timeYear, theta(1:),jphi(1:, np), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/j_phi_meridional.yg", &
                          & "Jph(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(1:),jphi(1:, j_slice), nang, &
                          & MAGNETIC_DATA_FORMAT)
      ! Electric field, meridional profile.
      call output_1d_ygraph("outb/e_r_surface.yg", &
                          & "Er(" // real_to_string(rb(np)) // ",th)", &
                          & timeYear, theta(1:), jr(1:, np), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/e_r_meridional.yg", &
                          & "Er(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(1:), er(1:, j_slice), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/e_th_surface.yg", &
                          & "Eth(" // real_to_string(rb(np)) // ",th)", &
                          & timeYear, theta(1:), eth(1:, np), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/e_th_meridional.yg", &
                          & "Eth(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(1:), eth(1:, j_slice), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/e_phi_surface.yg", &
                          & "Eph(" // real_to_string(rb(np)) // ",th)", &
                          & timeYear, theta(1:),ephi(1:, np), nang, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/e_phi_meridional.yg", &
                          & "Eph(" // real_to_string(rb(j_slice)) // ",th)", &
                          & timeYear, theta(1:), ephi(1:, j_slice), nang, &
                          & MAGNETIC_DATA_FORMAT)



      ! Radial profiles --------------------------------------------------------

      ! Potential vector, radial profiles
      call output_1d_ygraph("outb/a_phi_equator.yg", &
                          & "Aphi(r," // real_to_string(theta((nang+1)/2)) // ")", &
                          & timeYear, rb(1:), aphi(((nang+1)/2), 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/a_phi_radial.yg", &
                          & "Aphi(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:), aphi(i_slice, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      ! Magnetic field, radial profile.
      call output_1d_ygraph("outb/b_r_pole.yg", &
                          & "Br(r," // real_to_string(theta(1)) // ")", &
                          & timeYear, rb(1:), br(1, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/b_r_equator.yg", &
                          & "Br(r," // real_to_string(theta((nang+1)/2)) // ")", &
                          & timeYear, rb(1:), br(((nang+1)/2), 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/b_r_radial.yg", &
                          & "Br(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:), br(i_slice, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/b_th_radial.yg", &
                          & "Bth(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:), bth(i_slice, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/b_th_equator.yg", &
                          & "Bth(r," // real_to_string(theta((nang+1)/2)) // ")", &
                          & timeYear, rb(1:), bth(((nang+1)/2), 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/b_phi_radial.yg", &
                          & "Bph(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:), bphi(i_slice, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/b_phi_equator.yg", &
                          & "Bphi(r," // real_to_string(theta((nang+1)/2)) // ")", &
                          & timeYear, rb(1:), bphi(((nang+1)/2), 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      ! Electrical current, radial profile.
      call output_1d_ygraph("outb/j_r_pole.yg", &
                          & "Jr(r," // real_to_string(theta(1)) // ")", &
                          & timeYear, rb(1:), jr(1, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/j_r_equator.yg", &
                          & "Jr(r," // real_to_string(theta((nang+1)/2)) // ")", &
                          & timeYear, rb(1:), jr(((nang+1)/2), 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/j_r_radial.yg", &
                          & "Jr(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:), jr(i_slice, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/j_th_radial.yg", &
                          & "Jth(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:), jth(i_slice, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/j_th_equator.yg", &
                          & "Jth(r," // real_to_string(theta((nang+1)/2)) // ")", &
                          & timeYear, rb(1:), jth(((nang+1)/2), 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/j_phi_radial.yg", &
                          & "Jph(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:), jphi(i_slice, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/j_phi_equator.yg", &
                          & "Jphi(r," // real_to_string(theta((nang+1)/2)) // ")", &
                          & timeYear, rb(1:), jphi(((nang+1)/2), 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      ! Electric field, radial profile.
      call output_1d_ygraph("outb/e_r_pole.yg", &
                          & "Er(r," // real_to_string(theta(1)) // ")", &
                          & timeYear, rb(1:), er(1, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/e_r_equator.yg", &
                          & "Er(r," // real_to_string(theta((nang+1)/2)) // ")", &
                          & timeYear, rb(1:), er(((nang+1)/2), 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/e_r_radial.yg", &
                          & "Er(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:), er(i_slice, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/e_th_radial.yg", &
                          & "Eth(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:), eth(i_slice, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/e_th_equator.yg", &
                          & "Eth(r," // real_to_string(theta((nang+1)/2)) // ")", &
                          & timeYear, rb(1:), eth(((nang+1)/2), 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/e_phi_radial.yg", &
                          & "Eph(r," // real_to_string(theta(i_slice)) // ")", &
                          & timeYear, rb(1:), ephi(i_slice, 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      call output_1d_ygraph("outb/e_phi_equator.yg", &
                          & "Ephi(r," // real_to_string(theta((nang+1)/2)) // ")", &
                          & timeYear, rb(1:), ephi(((nang+1)/2), 1:), np, &
                          & MAGNETIC_DATA_FORMAT)
      ! Multipole weight.
      call output_1d_ygraph("outb/bl_out.yg", "bl (l)", &
                          & timeYear, dble((/ (l, l = 1, nleg) /)) , blout(1:nleg), nleg, &
                          & BL_FORMAT)
      ! Profile of estimation of Courant timestep
      if (maxval(dtb_courant_profile) /= 0) then
      call output_1d_ygraph("outb/dtb_courant_radial.yg", &
                          & "Log(dtb_Cour[year])(r)", &
                          & timeYear, rb(2:np), dlog10(dtb_courant_profile(2:np)), np-1)
      endif



    end subroutine output_magnetic
    !---------------------------------------------------------------------------
    !> Outputs a snapshot of the current magnetic field.
    !> @brief
    !
    !> @param[in]   timeYear      Timestamp in [yr].
    !---------------------------------------------------------------------------
    subroutine output_magnetic_map(timeYear)

      ! Modules ----------------------------------------------------------------
      ! None.

      implicit none

      ! Input parameters -------------------------------------------------------
      real*8, intent(in) :: timeYear

      ! Local constants --------------------------------------------------------
      character(len=10), parameter :: OUTPUT_FORMAT = "(13es15.6)"
      character(len=32), parameter :: MAGNETIC_FILE = "outb/bfield.dat"

      ! Local variables --------------------------------------------------------
      integer :: unit
      integer :: k, l

      ! ------------------------------------------------------------------------

      unit = get_free_unit()

      ! If we just started the simulation, create a new file and write header.
      if (timeYear == 0.0d0) then
        open(unit=unit, file=MAGNETIC_FILE, action="write")
      ! If not, just append to the end of the file.
      else
        open(unit=unit, file=MAGNETIC_FILE, status="old", access="append")
      end if

      ! Each snapshot starts with "snapshot [time in yr].
      write(unit, *) "snapshot", timeYear

      do k = 1, nang
        do l = 1, np
          write(unit, OUTPUT_FORMAT) &
            & theta(k), rb(l), br(k, l), &
            & bth(k, l), bphi(k, l), aphi(k, l), &
            & jr(k, l), jth(k, l), jphi(k, l), &
            & er(k, l), eth(k, l), ephi(k ,l)
        end do ! l
      end do ! k

      ! Each snapshot ends with "snapshot_end".
      write(unit, *) "snapshot_end"
      close(unit)

    end subroutine output_magnetic_map


  !---------------------------------------------------------------------------
  !> Outputs a snapshot of the global magnetic quantities
  !> @brief
  !
  !> @param[in]   timeYear      Timestamp in [yr].
  !---------------------------------------------------------------------------
  subroutine output_magnetic_analysis(iterb,timeYear)
  
    ! Modules ----------------------------------------------------------------
    use magnetic_analysis, only : en_mag_star, en_mag_star_tor, en_mag_magnetosphere
    use magnetic_analysis, only : delta_en_mag_star, delta_en_mag_magnetosphere, helicity_star
    use magnetic_analysis, only : j2_star, en_electric_star, divb_star_l2norm
    use magnetic_analysis, only : q_joule_star, q_joule_shock_star
    use magnetic_analysis, only : en_joule_star, en_joule_shock_star
    use magnetic_analysis, only : poynting_star, poynting_star_tot
  
    implicit none

    ! Input parameters -------------------------------------------------------
    real*8, intent(in) :: timeYear
    integer, intent(in) :: iterb

    ! Local constants --------------------------------------------------------
    character(len=13), parameter :: BTOT_FORMAT = "(i9,11es13.4)"
    character(len=12),  parameter :: BCONS_FORMAT = "(i9,7es13.4)"
    character(len=32), parameter :: BTOT_FILE = "outb/cons_integrated.dat"
    character(len=32), parameter :: BCONS_FILE = "outb/cons_instant.dat"

    ! Local variables --------------------------------------------------------
    integer :: unit
    real*8 :: balance
    ! ------------------------------------------------------------------------

    unit = get_free_unit()

    if (timeYear == 0.0d0) then
      open(unit=unit, file=BTOT_FILE, action="write")
      write(unit, *) "Magnetic conservation: iterb,t,Emag,EmagTor,EmagOut,Helicity,PoyTot,QjTot,QjschockTOT,J2,E2,divB(L2norm)"
    else
      open(unit=unit, file=BTOT_FILE, status="old", access="append")
    end if

    write(unit, BTOT_FORMAT) iterb, timeYear, en_mag_star, en_mag_star_tor,  &
   &  en_mag_magnetosphere, helicity_star, poynting_star_tot, &
   &  en_joule_star, en_joule_shock_star, j2_star, en_electric_star, divb_star_l2norm
    close(unit)

    write(6, *) "MAGNETIC EVOLUTION: iterb, time, Emag, Emag tor, Emag out, divB(L2 norm)"
    write(6, BTOT_FORMAT) iterb, timeYear, en_mag_star, en_mag_star_tor,  &
   &  en_mag_magnetosphere, divb_star_l2norm


    unit = get_free_unit()
    if (timeYear == 0.0d0) then
      open(unit=unit, file=BCONS_FILE, action="write")
      write(unit, *) "Magnetic conservation: iterb,t,dEmag,dEmagOut,Qj,Qjshock,Poynting,dEmag+Poynting-Qj-Qjshock"
    else
      open(unit=unit, file=BCONS_FILE, status="old", access="append")
      balance = poynting_star + delta_en_mag_star - q_joule_star - q_joule_shock_star
      write(unit, BCONS_FORMAT) iterb, timeYear, delta_en_mag_star, delta_en_mag_magnetosphere, &
     &  q_joule_star, q_joule_shock_star, poynting_star, balance
    end if
    close(unit)


    end subroutine output_magnetic_analysis 


    !---------------------------------------------------------------------------
    !> Outputs a 1D snapshot of any variable in ygraph format.
    !> @brief
    !> @param[in] file          Output filename.
    !> @param[in] label         Label for the snapshot.
    !> @param[in] timeYear      Time of the snapshot.
    !> @param[in] positions     Array of positions in 1D axis.
    !> @param[in] values        Array of values for the variable.
    !> @param[in] maxIdx        Maximum index for the positions.
    !> @param[in] format        Format for the output.
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
    write(unit, *) '"Time=', time
    
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

end module output

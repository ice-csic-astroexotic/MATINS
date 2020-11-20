module resume

  contains

    !---------------------------------------------------------------------------
    !> Reads a checkpoint to resume execution.
    !> @brief Reads a checkpoint which contains the minimum required information
    !>        to reinitialize the simulation from a certain time step. All
    !>        values are automatically loaded to the corresponding variables.
    !
    !> @param[in]   timeYear      Timestamp [yr].
    !> @param[in]   period        Period [s].
    !> @param[in]   periodDot     Period derivative (dimensionless).
    !---------------------------------------------------------------------------
    subroutine read_resume(timeYear, period, periodDot, bpdip)

      ! Modules ----------------------------------------------------------------
      use utils, only: get_free_unit
      use magnetic_analysis, only: en_joule_star, en_joule_shock_star
      use magnetic_analysis, only: poynting_star_tot
      use grid, only: br, bth, bphi, aphi, bm, bmed
      use grid, only: tem
      use grid, only: np, nang
      use grid, only: kmax, lmax, jevol, jmin

      implicit none

      ! Input parameters -------------------------------------------------------
      real*8, intent(out) :: timeYear
      real*8, intent(out) :: period
      real*8, intent(out) :: periodDot
      real*8, intent(out) :: bpdip

      ! Local constants --------------------------------------------------------

      ! Local variables --------------------------------------------------------
      integer :: unit
      integer :: i, j

      ! ------------------------------------------------------------------------

      unit = get_free_unit()

      open(unit=unit, file="resume.dat", action="read")

      ! Load scalars.
      read(unit,*) timeYear, period, periodDot, bpdip

      read(unit,*) jevol, jmin

      read(unit, *) en_joule_star, en_joule_shock_star, poynting_star_tot

      write(*, "(a, e17.8)") "RESUMING RUN FROM (TYEAR): ", timeYear

      ! Load magnetic field.
      do i = 0, nang + 1
        do j = 0, np + 2
          read(unit,*) br(i, j), bth(i, j), bphi(i, j), aphi(i, j)
        end do !j
      end do !i

      bm=dsqrt(br**2+bth**2+bphi**2)
      do i = 2, kmax
        do j = 1, lmax
          bmed(i,j) = bm(2*i-2,2*j-1)
        enddo !j
      enddo !i

      ! Load temperature map.
      do i = 1, kmax
        do j = 1, lmax
          read(unit,*) tem(i, j)
        end do !j
      end do !i

      close(unit)

    end subroutine read_resume

  end module resume
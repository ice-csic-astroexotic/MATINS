!-------------------------------------------------------------------------------
!> Compute the Shear modulus and the maximum strength.
!> @brief This subroutine will compute the Shears modulus and the maximum Shear
!>          strength for each cell of the grid following the C.J. Horowitz's
!>          method. A temperature grid is provided as input and the values
!>          for the Shear modulus and maximum strength are computed for each
!>          cell of such grid.
!
!> @param[in] kmax, lmax        Angular and radial dimensions of the arrays.
!> @param[in] tem0              Unredshifted temperature values for each cell.
!> @param[out] shearModulus     Shear modulus for each cell.
!> @param[out] shearMaximum     Maximum Shear strength for each cell.
!
!> @author
!> Jose Pons Botella
!> Alberto Garcia-Garcia
!
!> @note
!> TODO: (calculation following C.J. Horowitz  ADD REFERENCES !)
!-------------------------------------------------------------------------------


subroutine compute_shear()

  ! Modules --------------------------------------------------------------------
  use constants, only: PI, E2, MASS_N
  use grid, only: kmax, lmax
  use grid, only: rho, xh, zz, aa, tem0
  use grid, only: shearMaximum, shearModulus,jcore
  use grid, only: gammac
  implicit none


  ! Local variables ------------------------------------------------------------
  integer k, l, j
  ! TODO: would be nice to add some documentation for this variables.
  real*8 nion, alat
  ! ----------------------------------------------------------------------------

  ! Default values for the fluid phase. If we are not in the solid crust, use
  ! these values for convenience.
  shearModulus = 1.d32
  shearMaximum = 1.d32

  ! If the crust is solid, compute the Shear modules and maximum strength.
  do l = jcore/2+1, lmax

    j = 2 * l - 1 

    if(xh(j) > 0d0) then

      nion = xh(j) * rho(j) / (MASS_N * aa(j)) ! cm**-3
      alat = (3.d0 / (4.d0 * PI * nion))**(1.0 / 3.0)

      do k = 2, kmax

        gammac(k,l) = zz(j)**2 * E2 / alat / (1.38e-8 * tem0(k, l))

        if (gammac(k,l) > 175.d0) then

          shearModulus(k, l) = zz(j)**2 * E2 * nion / alat * &
            &                  (0.1106d0 - 28.7d0 / gammac(k,l)**1.3) ! erg/cm^3
          shearMaximum(k, l) = zz(j)**2 * E2 * nion / alat * &
            &               (0.0195d0 - 1.27d0 / (gammac(k,l) - 71d0))

        end if

      end do

    end if

  end do

end subroutine compute_shear

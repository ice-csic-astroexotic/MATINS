!> @brief TODO:
!!
!! TODO.
!!
!! @param[in]     m         TODO
!! @param[in]     a         TODO
!! @param[in]     b         TODO
!! @param[in]     c         TODO
!! @param[in]     r         TODO
!! @param[in]     gammaIn   TODO
!! @param[in]     gammaOut  TODO
!! @param[in]     uin       TODO
!! @param[in]     uout      TODO
!! @param[in]     i         TODO
!!
!!  Code owners:
!!    Jose A. Pons Botella
!!    Alberto Garcia-Garcia
!!
subroutine getugam(n, a, b, c, r, gammaIn, uIn, gammaOut, uOut)

  ! Module imports -------------------------------------------------------------
  use reals, only: double

  implicit none

  ! Subroutine arguments -------------------------------------------------------
  integer, intent(in) :: n
  real(double), intent(in) :: a(n, n), b(n, n), c(n, n)
  real(double), intent(in) :: r(n)
  real(double), intent(in) :: gammaIn(n, n)
  real(double), intent(in) :: uIn(n)
  real(double), intent(out) :: gammaOut(n, n)
  real(double), intent(out) :: uOut(n)

  ! Local constants ------------------------------------------------------------
  ! None.

  ! Local variables ------------------------------------------------------------
  ! Auxiliary variable for loops.
  integer :: i
  ! TODO: Document
  real(double) :: alpha(n, n)
  ! TODO: Document
  real(double) :: y(n, n)
  ! TODO: Document
  real(double) :: z(n)
  ! TODO: Document
  real(double) :: dd
  ! TODO: Document
  integer :: indices(n)

  ! ----------------------------------------------------------------------------

  ! Check maximum dimension.
  ! TODO.

  ! Initialize alpha.
  alpha = 0_double

  ! Calculate the matrix alpha = b - a * gamma_in.
  alpha = b - matmul(a, gammaIn)

  ! Decompose the alpha matrix just once.
  call ludcmp(n, alpha, indices, dd)

  ! Calculate gammaOut = alpha^-1 * c.
  y = c

  do i = 1, n

    call lubksb(n, alpha, indices, y(1, i))
    gammaOut(:, i) = y(:, i)

  end do ! i

  ! Calculate z = a * uIn.
  z = matmul(a, uIn)

  ! Calculate uOut = alpha^-1 * a * uIn.
  call lubksb(n, alpha, indices, z)
  uOut = -z

  ! Calculate alpha^-1 * r.
  z = r
  call lubksb(n, alpha, indices, z)
  uOut = uOut + z

end subroutine getugam
!> @brief Solves the set of n linear equations Ax = b.
!!
!! This subroutine solves a set of n linear equations of shape Ax = b where a
!! LU decomposition of the A matrix has been provided. It makes use of the LU
!! backsubstitution method.
!!
!! @param[in]     n       Number of linear equations to solve (size of A).
!! @param[in]     a       LU decomposed square matrix (A).
!! @param[in]     indices Vector coming from the LU decomposition that contains
!!                          the row permutations effected by partial pivoting.
!! @param[in/out] b       Vector on the right-hand side of the equation. It gets
!!                          overwritten by the solution vector.
!!
!!  Code owners:
!!    Jose A. Pons Botella
!!    Alberto Garcia-Garcia
!!
subroutine lubksb(n, a, indices, b)

  ! Module imports -------------------------------------------------------------
  use reals, only: double

  implicit none

  ! Subroutine arguments -------------------------------------------------------
  integer, intent(in) :: n
  real(double), intent(in) :: a(n, n)
  integer, intent(in) :: indices(n)
  real(double), intent(inout) :: b(n)

  ! Local constants ------------------------------------------------------------
  ! None.

  ! Local variables ------------------------------------------------------------
  ! Auxiliary loop indices.
  integer :: i, ii, ll
  real(double) :: sum

  ii = 0
  do i = 1, n
    ll = indices(i)
    sum = b(ll)
    b(ll) = b(i)

    if (ii /= 0) then
      sum = sum - dot_product(a(i, ii:i-1), b(ii:i-1))
    else if (sum /= 0.d0) then
      ii = i
    endif

    b(i) = sum
  end do

  do i = n, 1, -1
    b(i) = (b(i) - dot_product(a(i, i+1:n), b(i+1:n))) / a(i, i)
  end do

end subroutine lubksb

!> @brief Block-tridiagonal solver.
!!
!! Solves the block-tridiagonal system:
!!  | b1 c1 0  0  0.....| |u1|    |r1|
!!  | a2 b2 c2 0  0.....| |u2|    |r2|
!!  | 0  a3 b3 c3 0.....| |. |  _ |. |
!!  | ..................| |. |  _ |. |
!!  | ..................| |. |    |. |
!!  | .............an bn| |un|    |. |
!!
!! @param[in]     nx      TODO: What is this?
!! @param[in]     nz      TODO: What is this? Isn't this nx?
!! @param[in]     a       Lower diagonal of the block-tridiagonal matrix.
!! @param[in]     b       Middle diagonal of the block-tridiagonal matrix.
!! @param[in]     c       Upper diagonal of the block-tridiagonal matrix.
!! @param[in]     r       Vector of source values.
!! @param[out]    u       Vector of u-values to solve.
!!
!!  Code owners:
!!    Jose A. Pons Botella
!!    Alberto Garcia-Garcia
!!
subroutine solvtb(nx, nz, a, b, c, r, u)

  ! Module imports -------------------------------------------------------------
  use reals, only: double

  implicit none

  ! Subroutine arguments -------------------------------------------------------
  integer, intent(in) :: nx, nz
  real(double), intent(in) :: a(nx, nz, 3)
  real(double), intent(in) :: b(nx, nz, 3)
  real(double), intent(in) :: c(nx, nz, 3)
  real(double), intent(in) :: r(nx, nz)
  real(double), intent(out) :: u(nx, nz)

  ! Local constants ------------------------------------------------------------
  ! None.

  ! Local variables ------------------------------------------------------------
  ! Auxiliary counters for loops.
  integer :: i, l
  ! TODO: What is this?
  real(double) :: aq(nz, nz), bq(nz, nz), cq(nz, nz)
  ! TODO: What is this?
  real(double) :: rq(nz)
  ! TODO: What is this?
  real(double) :: gamin(nz, nz), gamout(nz, nz), gam(nx, nz, nz)
  ! TODO: What is this?
  real(double) :: uin(nz), uout(nz)

  uin=0.d0
  gamin=0.d0
  aq=0.d0
  bq=0.d0
  cq=0.d0

  do l = 1, nx

    rq = r(l, :)

    do i = 1, nz

      if (i == 1) then
        aq(i, i) =a(l, i, 2)
        aq(i, i+1) =a(l, i ,3)
        bq(i, i) = b(l, i, 2)
        bq(i, i+1) = b(l, i, 3)
        cq(i, i) = c(l, i, 2)
        cq(i, i+1) = c(l, i, 3)
      else if (i == nz) then
        aq(i, i-1) = a(l, i, 1)
        aq(i, i) = a(l, i, 2)
        bq(i, i-1) = b(l, i, 1)
        bq(i, i) = b(l, i, 2)
        cq(i, i-1) = c(l, i, 1)
        cq(i, i) = c(l, i, 2)
      else
        aq(i, i-1) = a(l, i, 1)
        aq(i, i) = a(l, i, 2)
        aq(i, i+1) = a(l, i, 3)
        bq(i, i-1) = b(l, i, 1)
        bq(i, i) = b(l, i, 2)
        bq(i, i+1) = b(l, i, 3)
        cq(i, i-1) = c(l, i, 1)
        cq(i, i) = c(l, i, 2)
        cq(i, i+1) = c(l, i, 3)
      endif
    end do

    call getugam(nz, aq, bq, cq, rq, gamin, uin, gamout, uout)

    ! Fetch output u-values and prepare the input for the next iteration.
    u(l, :) = uout
    uin = uout
    ! Fetch output gamma values and prepare the input for the next iteration.
    gam(l, :, :) = gamout
    gamin = gamout

    ! TODO: Check that we truly need the input values gamin and uin in getugam!

  end do

  ! TODO: What is this for?
  do l = nx-1, 1, -1
    u(l, :) = u(l, :) - matmul(gam(l, :, :), u(l + 1, :))
  end do

end subroutine solvtb
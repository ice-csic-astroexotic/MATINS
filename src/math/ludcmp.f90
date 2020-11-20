!> @brief Replaces a square matrix with the LU decomposition of itself.
!!
!! Given an NxN matrix A, with dimension NP, this routine replaces it by the LU
!! decomposition of a row-wise permutation of itself. In general, it is used in
!! combination with LUBKSB to solve linear equations or invert a matrix.
!!
!! @param[in]     n       Size of the input matrix (n x n).
!! @param[in/out] a       Input matrix.
!! @param[out]    indices Vector that records the row permutation resulting
!!                          from partial pivoting 
!!                          the row permutations effected by partial pivoting.
!! @param[out]    d       Indicates whether the number of row interchanges was
!!                          even (+1) or odd (-1).
!!
!!  Code owners:
!!    Jose A. Pons Botella
!!    Alberto Garcia-Garcia
!!
subroutine ludcmp(n, a, indices, d)

  ! Module imports -------------------------------------------------------------
  use reals, only: double

  implicit none

  ! Subroutine arguments -------------------------------------------------------
  integer, intent(in) :: n
  real(double), intent(inout) :: a(n, n)
  integer, intent(out) :: indices(n)
  real(double), intent(out) :: d

  ! Local constants ------------------------------------------------------------
  ! TODO: Document NMAX
  integer, parameter :: NMAX = 5000
  ! TODO: Document TINY
  real(double), parameter :: TINY = 1.0d-25

  ! Local variables -----------------------------------------------------------
  ! Auxiliary loop indices.
  integer :: i, j, k
  ! TODO: Document
  integer :: imax
  ! TODO: Document
  real(double) :: aamax
  ! TODO: Document
  real(double) :: dum
  ! TODO: Document
  real(double) :: sum
  ! TODO: Document
  real(double) :: vv(NMAX)

  d = 1.d0

  do i = 1, n
    aamax = 0.d0
    do j = 1, n
      if (abs(a(i, j)) > aamax) then
        aamax = abs(a(i, j))
      end if
    end do

    if (aamax == 0.d0) then
      stop 'Singular matrix in ludcmp'
    end if

    vv(i) = 1.d0 / aamax
  end do

  do j = 1, n
    do i = 1, j-1
      sum = a (i, j)
      do k = 1, i-1
        sum = sum - a(i, k) * a(k, j)
      end do
      a(i ,j) = sum
    end do

    aamax = 0.d0

    do i = j, n
      sum = a(i, j)
      do k = 1, j-1
        sum = sum - a(i, k) * a(k, j)
      end do
      a(i, j) = sum

      dum = vv(i) * abs(sum)

      if (dum >= aamax) then
        imax = i
        aamax = dum
      endif
    end do

    if (j /= imax) then
      do k = 1, n
        dum = a(imax, k)
        a(imax, k) = a(j, k)
        a(j, k) = dum
      end do

      d = -d
      vv(imax) = vv(j)
    end if

    indices(j) = imax

    if (a(j, j) == 0.d0) then
      a(j, j) = TINY
    end if

    if (j /= n) then
      dum = 1.d0 / a(j, j)
      do i = j+1, n
        a(i, j) = a(i, j) * dum
      end do
    end if

  end do

end subroutine ludcmp

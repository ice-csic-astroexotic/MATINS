integer function test_math_compute_gausslegendre_weights() result(r)

  use math, only: compute_gausslegendre_weights

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: EPSILON = 1e-4_dp
  integer, parameter :: N = 4

  real(dp) :: x1
  real(dp) :: x2
  real(dp) :: x(N)
  real(dp) :: w(N)

  real(dp) :: x_result(N)
  real(dp) :: w_result(N)
  
  integer :: i, j

  ! Test setup ----------------------------------------------------------------
  x1 = -1.0_dp
  x2 = 1.0_dp
  x = 0.0_dp
  w = 0.0_dp

  ! Expected results ----------------------------------------------------------
  x_result = (/ -0.861136_dp, -0.339981_dp, 0.339981_dp, 0.861136_dp /)
  w_result = (/ 0.347854_dp, 0.652145_dp, 0.652145_dp, 0.347854_dp /)

  ! Test ----------------------------------------------------------------------

  call compute_gausslegendre_weights(n, x1, x2, x, w)

  r = 0

  do i = 1, n
    write(*,*) "x(",i,"): ", x(i)
    if (abs(x(i) - x_result(i)) > epsilon) then
      write(*,*) "Error in x[", i, "] -- expected ", x_result(i), " got ", x(i)
      r = 1
    end if
    write(*,*) "w(",i,"): ", w(i)
    if (abs(w(i) - w_result(i)) > epsilon) then
      write(*,*) "Error in w[", i, "] -- expected ", w_result(i), " got ", w(i)
      r = 1
    end if
  end do

end function
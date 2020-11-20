integer function test_math_odeint() result(r)

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: epsilon = 1e-6_dp

  integer, parameter :: NUM_VARS = 3

  real(dp) :: h1 = -0.53264217542849579_dp
  real(dp) :: x1 = -32.122906915326141_dp
  real(dp) :: x2 = -85.387124458175720_dp
  real(dp) :: y(NUM_VARS) = (/ 77146.901861056540_dp, 0.0_dp, 0.0_dp /)
  real(dp) :: num_ok
  real(dp) :: num_bad

  integer :: i

  real(dp) :: result(NUM_VARS) = (/ 1171100.7665579075_dp, 208652.36643989259_dp, 0.48137774983963327_dp /)

  interface
    subroutine derivs(x, y, dydx)
      real*8 x, y(:), dydx(:)
    end subroutine derivs
  end interface

  call init_eos_tab()
  call odeint(y, NUM_VARS, x1, x2, h1, num_ok, num_bad, derivs)

  r = 0

  do i = 1, size(y)
    if (abs(y(i) - result(i)) > epsilon) then
      write(*,*) "Error in ", i, "-- expected ", result(i), " got ", y(i)
      r = 1
    end if
  end do ! i
  r = 0

end function
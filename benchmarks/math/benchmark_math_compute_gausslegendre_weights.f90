program benchmark_math_compute_gausslegendre_weights

  use reals, only: dp => double
  use math, only: compute_gausslegendre_weights

  implicit none

  integer, parameter :: SIZE = 1000

  real(dp) :: x1
  real(dp) :: x2
  real(dp) :: x(SIZE)
  real(dp) :: w(SIZE)

  real*8 :: start_time, stop_time
  real*8 :: avg_time

  integer :: runs
  integer :: i

  x1 = 1.0_dp
  x2 = -1.0_dp
  x = 0.0_dp
  w = 0.0_dp

  avg_time = 0.0d0
  runs = 100

  do i = 1, runs
    call cpu_time(start_time)
    call compute_gausslegendre_weights(SIZE, x1, x2, x, w)
    call cpu_time(stop_time)
    print *, "Run ", i, " time: ", (stop_time - start_time), " [s]"
    avg_time = avg_time + (stop_time - start_time)
  end do ! i

  avg_time = avg_time / runs

  print *, "Average time:", avg_time

end program benchmark_math_compute_gausslegendre_weights
program benchmark_math_ludcmp

  use reals, only: dp => double
  use math, only: ludcmp

  implicit none

  integer, parameter :: SIZE = 1000

  real(dp) :: a(SIZE, SIZE)
  integer :: indices(SIZE * SIZE)
  real(dp) :: d

  real(dp) :: start_time, stop_time
  real(dp) :: avg_time

  integer :: runs
  integer :: i

  a = 1.0d0
  indices = 0.0d0
  d = 0.0d0
  avg_time = 0.0d0
  runs = 100

  do i = 1, runs
    call cpu_time(start_time)
    call ludcmp(SIZE, a, indices, d)
    call cpu_time(stop_time)
    print *, "Run ", i, " time: ", (stop_time - start_time), " [s]"
    avg_time = avg_time + (stop_time - start_time)
  end do ! i

  avg_time = avg_time / runs

  print *, "Average time:", avg_time

end program benchmark_math_ludcmp
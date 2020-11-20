program benchmark_math_getugam

  use reals, only: dp => double
  use math, only: getugam

  implicit none

  integer, parameter :: SIZE = 1000

  real(dp) :: a(SIZE, SIZE) 
  real(dp) :: b(SIZE, SIZE) 
  real(dp) :: c(SIZE, SIZE)
  real(dp) :: source(SIZE)
  real(dp) :: gamma_in(SIZE, SIZE)
  real(dp) :: u_in(SIZE)
  real(dp) :: gamma_out(SIZE, SIZE)
  real(dp) :: u_out(SIZE)

  real*8 :: start_time, stop_time
  real*8 :: avg_time

  integer :: runs
  integer :: i

  a = 1.0d0
  b = 1.0d0
  c = 1.0d0
  source = 1.0d0
  gamma_in = 1.0d0
  u_in = 1.0d0
  gamma_out = 1.0d0
  u_out = 1.0d0

  avg_time = 0.0d0
  runs = 100

  do i = 1, runs
    call cpu_time(start_time)
    call getugam(SIZE, a, b, c, source, gamma_in, u_in, gamma_out, u_out)
    call cpu_time(stop_time)
    print *, "Run ", i, " time: ", (stop_time - start_time), " [s]"
    avg_time = avg_time + (stop_time - start_time)
  end do ! i

  avg_time = avg_time / runs

  print *, "Average time:", avg_time

end program benchmark_math_getugam
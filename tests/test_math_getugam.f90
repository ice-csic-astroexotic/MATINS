integer function test_math_getugam() result(r)

  use math, only: getugam

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: epsilon = 1e-4_dp
  integer, parameter :: n = 3

  real(dp) :: a(n, n) 
  real(dp) :: b(n, n) 
  real(dp) :: c(n, n)
  real(dp) :: source(n)
  real(dp) :: gamma_in(n, n)
  real(dp) :: u_in(n)
  real(dp) :: gamma_out(n, n)
  real(dp) :: u_out(n)

  real(dp) :: u_result(n)
  
  integer :: i, j

  a = reshape(&
    (/4.0_dp, -2.0_dp, 1.0_dp, &
      -3.0_dp, -1.0_dp, 4.0_dp, &
      1.0_dp, -1.0_dp, 3.0_dp/), &
    (/3,3/))
  b = reshape(&
    (/4.0_dp, -2.0_dp, 1.0_dp, &
      -3.0_dp, -1.0_dp, 4.0_dp, &
      1.0_dp, -1.0_dp, 3.0_dp/), &
    (/3,3/))
  c = reshape(&
    (/4.0_dp, -2.0_dp, 1.0_dp, &
      -3.0_dp, -1.0_dp, 4.0_dp, &
      1.0_dp, -1.0_dp, 3.0_dp/), &
    (/3,3/))

  source = (/1.0_dp, 2.0_dp, 3.0_dp/)
  u_in = (/1.0_dp, 1.0_dp, 1.0_dp/)
  gamma_in = reshape(&
    (/4.0_dp, -2.0_dp, 1.0_dp, &
      -3.0_dp, -1.0_dp, 4.0_dp, &
      1.0_dp, -1.0_dp, 3.0_dp/), &
    (/3,3/))

  u_result = (/-0.83333_dp, -1.36666_dp, 1.56666_dp /)

  call getugam(n, a, b, c, source, gamma_in, u_in, gamma_out, u_out)

  r = 0

  do i = 1, n
    write(*,*) "Uout(",i,"): ", u_out(i)
    if (abs(u_out(i) - u_result(i)) > epsilon) then
      write(*,*) "Error in ", i, " -- expected ", u_result(i), " got ", u_out(i)
      r = 1
    end if
  end do 
end function
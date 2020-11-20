integer function test_math_ludcmp() result(r)

  use math, only: ludcmp

  implicit none

  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: epsilon = 1e-6_dp

  real(dp) :: a(3, 3) = reshape(&
    (/4.0_dp, -2.0_dp, 1.0_dp, &
     -3.0_dp, -1.0_dp, 4.0_dp, &
      1.0_dp, -1.0_dp, 3.0_dp/), &
    (/3,3/))

  integer :: n = 3
  integer :: np = 3
  integer :: indices(3)
  real*8 :: d
  integer :: i, j

  real(dp) :: result(3, 3) = reshape(&
    (/-2.0_dp, -2.0_dp, -0.5_dp, &
      -1.0_dp, -5.0_dp, -0.7_dp, &
      -1.0_dp, -1.0_dp, 1.8_dp/), &
    (/3, 3/))

  call ludcmp(n, a, indices, d)

  r = 0

  do i=1,3
    do j=1,3
      if (abs(a(i,j) - result(i, j)) > epsilon) then
        write(*,*) "Error in ", i, ",", j, "-- expected ", result(i, j), " got ", a(i, j)
        r = 1
      end if
    end do
  end do 
end function
!> @brief Calculates zeros and weights for Gauss-Legendre quadrature.
!!
!! @param[in]     n       Dimension of the xg, wg arrays
!! @param[in]     x1      Lower limit of the range (usually -1)
!! @param[in]     x2      Upper limit of the range (usually +1)
!! @param[out]    x       n zeros in the interval [x1,x2] where weights are computed
!! @param[out]    w       Gauss-Legendre weights
!!
!!  Code owners:
!!    Jose A. Pons
!!------------------------------------------------------------------------------
!! To be used to integrate by Gaussian Quadrature (cf. Numerical Recipes).
!! Gives zeroes and weights of Nth order Legendre polynomials.
!!------------------------------------------------------------------------------
!!
subroutine compute_gausslegendre_weights(n, x1, x2, x, w)

  ! Module imports -------------------------------------------------------------
  use constants, only: PI

  implicit none

  ! Subroutine arguments -------------------------------------------------------
  integer, intent (in) :: n
  real*8, intent (in) :: x1, x2
  real*8, intent (out) :: x(n), w(n)

  ! Local constants ------------------------------------------------------------
  real*8, parameter :: EPS = 3.0E-14

  ! Local variables ------------------------------------------------------------
  ! Auxiliary variables for loops.
  integer :: m, i, j
  ! TODO: What the hell are these?
  real*8 :: xm, xl, z, p1, p2, p3, pp, z1

  ! ----------------------------------------------------------------------------

  ! TODO: This needs a closer look and hopefully some documentation.
  m = (n+1)/2
  xm = 0.5D0*(x2+x1)
  xl = 0.5D0*(x2-x1)

  do i = 1, m

    z1 = 0.D0
    z = cos(PI*(i-.25D0)/(n+.5D0))

    do while (abs(z-z1)>eps)

      p1 = 1.D0
      p2 = 0.D0

      do j = 1, n

        p3 = p2
        p2 = p1
        p1 = ((2.D0*j-1.D0)*z*p2-(j-1.D0)*p3)/j

      end do ! j

      pp = n*(z*p1-p2)/(z*z-1.D0)
      z1 = z
      z = z1 - p1/pp

    end do ! while

    x(i) = xm - xl*z
    x(n+1-i) = xm + xl*z
    w(i) = 2.D0*xl/((1.D0-z*z)*pp*pp)
    w(n+1-i) = w(i)

  end do ! i

end subroutine compute_gausslegendre_weights
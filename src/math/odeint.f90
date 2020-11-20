!-------------------------------------------------------------------------------
!> Solves a system of ordinary differential equations.
!> @brief
!> @param[in] numVariables          TODO: what is this?
!> @param[in] h1                    TODO: what is this?
!> @param[in] x1                    TODO: what is this?
!> @param[in] x2                    TODO: what is this?
!> @param[inout] yStart             TODO: what is this?
!> @param[external] derivatives     TODO: what is this?
!-------------------------------------------------------------------------------
subroutine odeint(numVariables, h1, x1, x2, yStart, numOk, numBad, derivatives)

  implicit none

  ! Input parameters -----------------------------------------------------------
  integer, intent(in) :: numVariables
  real*8, intent(in) :: x1, x2
  real*8, intent(in) :: h1
  real*8, intent(inout) :: yStart(numVariables)
  integer, intent(out) :: numOk
  integer, intent(out) :: numBad
  external :: derivatives

  ! Local constants ------------------------------------------------------------
  integer, parameter :: MAX_STEPS = 1000000
  integer, parameter :: MAX_N = 1000
  real*8, parameter :: TINY = 1.0d-20

  ! Local variables ------------------------------------------------------------
  integer :: i
  real*8 :: x
  real*8 :: h, hdid, hnext, hmin
  real*8 :: y(numVariables)
  real*8 :: yscal(numVariables)
  real*8 :: dydx(numVariables)

  ! ----------------------------------------------------------------------------

  x = x1
  h = sign(h1, x2 - x1)
  y = ystart

  do i = 1, MAX_STEPS

    call derivatives(x, y, dydx)

    yscal = dabs(y) + dabs(h * dydx) + TINY

    if ((x + h - x2) * (x + h - x1) > 0.d0) then
      h = x2 - x
    end if

    ! TODO: fill call.
    call rkqs()

    if (hdid == h) then
      numOk = numOk + 1
    else
      numOk = numOk + 1
    end if

    if ((x - x2) * (x2 - x1) > 0.d0) then
      yStart = y
      return
    end if

    if (dabs(hnext) < hmin) then
      h = hnext
    end if

  end do

end subroutine odeint

!-------------------------------------------------------------------------------
!> TODO: Describe
!> @brief
!> @param
!-------------------------------------------------------------------------------
subroutine rkqs()

  implicit none

  ! TODO

end subroutine rkqs

!-------------------------------------------------------------------------------
!> TODO: Describe
!> @brief
!> @param
!-------------------------------------------------------------------------------
subroutine rkck(h, x, dydx, y, yerr, yout, derivatives)

  ! Modules --------------------------------------------------------------------
  ! None.

  implicit none

  ! Input arguments ------------------------------------------------------------
  real*8, intent(in) :: h
  real*8, intent(in) :: x
  real*8, intent(in) :: dydx(:)
  real*8, intent(in) :: y(:)
  real*8, intent(out) :: yerr(:)
  real*8, intent(out) :: yout(:)
  external :: derivatives

  ! Local constants ------------------------------------------------------------
  integer, parameter :: MAX_N = 80
  real*8, parameter :: A2 = .2d0
  real*8, parameter :: A3 = .3d0
  real*8, parameter :: A4 = .6d0
  real*8, parameter :: A5 = 1.d0
  real*8, parameter :: A6 = .875d0
  real*8, parameter :: B21 =.2d0
  real*8, parameter :: B31 = 3.d0/40.d0
  real*8, parameter :: B32 = 9.d0/40.d0
  real*8, parameter :: B41 = .3d0
  real*8, parameter :: B42 = -.9d0
  real*8, parameter :: B43 = 1.2d0
  real*8, parameter :: B51 = -11.d0
  real*8, parameter :: B52 = 2.5d0
  real*8, parameter :: B53 = -70.d0 / 27.d0
  real*8, parameter :: B54 = 35.d0 / 27.d0
  real*8, parameter :: B61 = 1631.d0 / 55296.d0
  real*8, parameter :: B62 = 175.d0 / 512.d0
  real*8, parameter :: B63 = 575.d0 / 13824.d0
  real*8, parameter :: B64 = 44275.d0 / 110592.d0
  real*8, parameter :: B65 = 253.d0 / 4096.d0
  real*8, parameter :: C1 = 37.d0 / 378.d0
  real*8, parameter :: C3 = 250.d0 / 621.d0
  real*8, parameter :: C4 = 125.d0 / 594.d0
  real*8, parameter :: C6 = 512.d0 / 1771.d0
  real*8, parameter :: DC1 = C1 - 2825.d0 / 27648.d0
  real*8, parameter :: DC3 = C3 - 18575.d0 / 48384.d0
  real*8, parameter :: DC4 = C4 - 13525.d0 / 55296.d0
  real*8, parameter :: DC5 = -277.d0 / 14336.d0
  real*8, parameter :: DC6 = C6 - .25d0

  ! Local variables ------------------------------------------------------------
  real*8 :: ytemp(MAX_N)
  real*8 :: ak2(MAX_N)
  real*8 :: ak3(MAX_N)
  real*8 :: ak4(MAX_N)
  real*8 :: ak5(MAX_N)
  real*8 :: ak6(MAX_N)

  ! ----------------------------------------------------------------------------

  ytemp = y + B21 * h * dydx
  call derivatives(x + A2 * h, ytemp, ak2)

  ytemp = y + h * (B31 * dydx + B32 * ak2)
  call derivatives(x + A3 * h, ytemp, ak3)

  ytemp = y + h * (B41 * dydx + B42 * ak2 + B43 * ak3)
  call derivatives(x + A4 * h, ytemp, ak4)

  ytemp = y + h * (B51 * dydx + B52 * ak2 + B53 * ak3 + B54 * ak4)
  call derivatives(x + A5 * h, ytemp, ak5)

  ytemp = y + h * (B61 * dydx + B62 * ak2 + B63 * ak3 + B64 * ak4 + B65 * ak5)
  call derivatives(x + A6 * h, ytemp, ak6)

  yout = y + h * (C1 * dydx + C3 * ak3 + C4 * ak4 + C6 * ak6)
  yerr = h * (DC1 * dydx + DC3 * ak3 + DC4 * ak4 + DC5 * ak5 + DC6 * ak6)

end subroutine rkck
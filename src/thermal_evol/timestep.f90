!-----------------------------------------------------------------------
!
! Contents: this subroutine prescribe an increasing cooling timestep by hand
! 
!-----------------------------------------------------------------------
subroutine adaptive_cooling_timestep(tyear,dt)
    
  implicit none

  ! Arguments.
  real*8, intent(in) :: tyear
  real*8, intent(out) :: dt

  ! Constants.
    real*8, parameter :: eps = 0.01
  ! None

    if ( tyear <= 9.99d0 + eps  ) then
      dt = 1d-1
    elseif ( tyear <= 9.99d2 + eps ) then
      dt = 1d0
    elseif ( tyear <= 9.99d3 + eps ) then
      dt = 1d1
    elseif ( tyear <= 9.99d4 + eps ) then
      dt = 1d2
    else
      dt = 1d2
    endif

end subroutine adaptive_cooling_timestep

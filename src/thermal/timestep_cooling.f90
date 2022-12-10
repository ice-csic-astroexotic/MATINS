!-----------------------------------------------------------------------
!
! Contents: this subroutine prescribe an increasing cooling timestep by hand
! 
!-----------------------------------------------------------------------
subroutine adaptive_cooling_timestep(time,dt)
    
  use input_params, only: max_dt_cooling
  implicit none

  ! Arguments.
  real*8, intent(in) :: time
  real*8, intent(out) :: dt

  ! Constants.
  real*8, parameter :: eps = 1d-9
  

  ! TBD: optmize the timestep
  if ( time <= 0.999d-6 + eps ) then ! < 1 yrs
    dt = 1d-8
  elseif ( time <= 0.999d-5 + eps ) then ! 1-10 yrs
    dt = 1d-7
  elseif ( time <= 0.999d-4 + eps ) then ! 10-100 yrs
    dt = 1d-6
  elseif ( time <= 0.999d-3 + eps ) then ! 100-1000 kyrs
    dt = 1d-5
  elseif ( time <= 0.999d-2 + eps ) then ! 1000-10000 kyrs 
    dt = 1d-4
  else ! > 10000 kyrs
    dt = 1d-3
  endif

  dt = min(dt,max_dt_cooling)

end subroutine adaptive_cooling_timestep

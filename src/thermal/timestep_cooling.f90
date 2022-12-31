!-----------------------------------------------------------------------
!
! Contents: this subroutine prescribe an increasing cooling timestep
! Authors: Daniele Viganò
! 
!-----------------------------------------------------------------------
subroutine adaptive_cooling_timestep(time,dt)
    
  use input_params, only: min_dt_cooling, max_dt_cooling, eps_dt_cooling
  implicit none

  ! Arguments.
  real*8, intent(in) :: time
  real*8, intent(out) :: dt

  ! A large enough minimum dt is needed at the beginning,
  ! only if one wants to follow well the very first years (usually not)
  dt = max(eps_dt_cooling*time,min_dt_cooling)

  ! A maximum timestep is advised to follow well the steep cooling
  ! at late times
  dt = min(dt,max_dt_cooling)

end subroutine adaptive_cooling_timestep

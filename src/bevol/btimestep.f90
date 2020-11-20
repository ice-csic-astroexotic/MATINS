!-----------------------------------------------------------------------
!!
! Contents:
! 
!-----------------------------------------------------------------------
subroutine adaptive_magnetic_timestep(tyear,tyear_b,tyear_b_print,dt,dtb)
  use grid, only: np, bm, lr, lth, etab, ia_hall, fh, dtb_courant_profile
  use input_params, only : courant_prefactor, dtb0, magnetic_output_dt
    
  implicit none

  ! Arguments.
  real*8, intent(in) :: tyear,tyear_b,tyear_b_print,dt
  real*8, intent(out) :: dtb

  ! Constants.
  real*8, parameter :: eps_dt = 1.d-5 ! This parameter has to be larger than zero to avoid to find dtb=0

  ! Definitions specific to the subroutine.
  integer :: j
  real*8 :: time_left
!-----------------------------------------------------------------------
! Notes:
! Courant time is defined as dt = courant_prefactor*min(dx**2/(fh*B + eta)).
! courant_prefactor is given in the input file (safe value for Hall is 1d-2).
! The maximum is evalued everywhere except at the outer boundary, where
! J is much higher.
! fh*|B| is in km^2/Myr, dr in km.
!-----------------------------------------------------------------------

  ! Initialize the array dtb_courant_profile to some large value.
  dtb_courant_profile = 1d10
  do j=2,np
    dtb_courant_profile(j) = 1d6*minval([lr(j)**2,lth(j)**2])/ &
    &   maxval(ia_hall(j)*fh(j)*bm(:,j) + etab(:,j))
  enddo

  ! The Courant estimation is given by the minimum value of dtb_courant_profile,
  ! multiplied by a Courant prefactor set in the input
  dtb = courant_prefactor*minval(dtb_courant_profile)

  ! Re-adjust the timestep, with a minimum value set in the input to avoid very slow runs
  dtb = maxval( [dtb, dtb0] )

  if (dtb .le. 1d-5) then
    write(*,*) "WARNING btimestep.f90 (ignore if it is sporadic): very small dtb, at tyear_b = ",dtb, tyear_b
  endif

  ! Maximum values are given the time left
  ! to reach the next magnetic output or to reach the next cooling loop
  ! eps_dt avoids that the difference is close to 0 providing a very small timestep
  time_left = tyear + dt + eps_dt - tyear_b
  if (time_left > 0. .and. time_left < dtb) then
    dtb = time_left
  endif

  time_left =  magnetic_output_dt + 0.1*eps_dt - tyear_b_print
  if (time_left > 0. .and. time_left < dtb) then
    dtb = time_left
  endif

  ! Check the value of dtb.
  if (dtb <= 0d0) then
     write(*,*) "ERROR in btimestep.f90: dtb has to be positive, dtb = ",dtb
     stop
  endif

end subroutine adaptive_magnetic_timestep

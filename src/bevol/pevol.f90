!> @brief Rotational evolution
!!
!! This subroutine evolves the spin period and its derivative considering
!! the dipolar value of the surface magnetic field
!!
!! @param[in]  dt          timestep considered (used to calculate the braking index)
!! @param[in]  bpdip       dipolar component of the surface magnetic field at the pole
!! @param[in]  spindown_prefactor         pre-coefficient of the spindown formula (coming from input.f)
!! @param[in]  fchi        coefficient related to the inclination angle
!!                         in the spindown formula (coming from input.f)
!! @param[in/out] per      spin period
!! @param[in/out] pdot     spin period derivative
!! @param[in/out] bindex   braking index, defined as:
!!                         bindex = Omegadotdot Omega / Omegadot^2 = 2 - Pdotdot P / Pdot^2
!!                                = 3 - 2 P Bpdot / Bp Pdot (using P Pdot = k Bp^2).
!!
!!  Code owners:
!!    Daniele Viganò
!!
subroutine period_evol(dtb,bpdip,bpdip_old,per,pdot,bindex)

  ! Modules --------------------------------------------------------------------
  use constants, only: T_YEAR
  use grid, only: spindown_prefactor

  ! Subroutine arguments -------------------------------------------------------
  implicit none
  real*8, intent(in) :: bpdip,bpdip_old,dtb
  real*8, intent(inout) :: per,pdot
  real*8, intent(out)     :: bindex

  ! Local constants ------------------------------------------------------------
  real*8, parameter :: fchi = 1d0   ! 1 for vacuum orthogonal rotator

  ! Local variables ------------------------------------------------------------
  ! Auxiliary variable to store the old period
  real*8                  :: per_old

  bindex = 0.
  per_old=per
  pdot=spindown_prefactor*bpdip**2*fchi/per
  per=dsqrt(per**2+2d0*dtb*T_YEAR*per*pdot)
  if ((per-per_old) /= 0.) bindex=3d0-2d0*per*(bpdip-bpdip_old)/(bpdip*(per-per_old))

  return
end subroutine period_evol

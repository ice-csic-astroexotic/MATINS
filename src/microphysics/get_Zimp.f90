!!-----------------------------------------------------------
!!    Subroutine that provides the impurity parameter Zimp
!!          (input for potekhin2019 routine)  
!!           Note the notation:  Zimp = sqrt(Qimp)
!! -----------------------------------------------------------
!!    The present version interpolates Qimp to
!!      include the effects of the pasta phase smoothly. 
!!       Alternative fits (Jones 2004) are commented out
!! -----------------------------------------------------------
subroutine get_Zimp(nb, Zimp)

  use input_params, only: Qimp, Qpasta
  implicit none
  real*8, intent(in) :: nb
  real*8, intent(out) :: Zimp
  real*8, parameter :: nbpasta=0.04, nbtrans=1.d-5

!  if (nb < nbtrans) then
!    Zimp=dsqrt(Qimp)
!   else if (nb < nbpasta) then
  if (nb < nbpasta) then
    Zimp= dsqrt((nb/2.d-5)**0.45d0)
  else
    Zimp= dsqrt(Qpasta)
  endif
! -----------------------------------------------------------
!     JONES AMORPHOUS CRUST
! Linear fit taking Q from Jones 2004, Tab. I, sp
!   	Zimp=sqrt(6.3+1.3e-13*rho)
! Linear fit taking Q from Jones 2004, Tab. I, p
!   	Zimp=sqrt(4.1+2.4e-13*rho)
! -----------------------------------------------------------
end subroutine get_Zimp

!!-----------------------------------------------------------
!!    Subroutine that provides the impurity parameter Zimp
!!          (input for potekhin2019 routine)  
!!           Note the notation:  Zimp = sqrt(Qimp)
!! -----------------------------------------------------------
!!    The present version interpolates Qimp to
!!      include the effects of the pasta phase smoothly. 
!!       Alternative fits (Jones 2004) are commented out
!! -----------------------------------------------------------
subroutine get_Zimp(rho, Zimp)

  use input_params, only: Qimp, Qpasta
  implicit none
  real*8 rho, Zimp
  real*8, parameter :: rhopasta=8.e13, rhotrans=1.e13

  if (rho < rhotrans) then
    Zimp=dsqrt(Qimp)
  elseif (rho < rhopasta) then
    Zimp=dsqrt(Qimp+(rho-rhotrans)**2/(rhopasta-rhotrans)**2*(Qpasta-Qimp))
  else
    Zimp=dsqrt(Qpasta)
  endif
! -----------------------------------------------------------
!     JONES AMORPHOUS CRUST
! Linear fit taking Q from Jones 2004, Tab. I, sp
!   	Zimp=sqrt(6.3+1.3e-13*rho)
! Linear fit taking Q from Jones 2004, Tab. I, p
!   	Zimp=sqrt(4.1+2.4e-13*rho)
! -----------------------------------------------------------

end subroutine get_Zimp

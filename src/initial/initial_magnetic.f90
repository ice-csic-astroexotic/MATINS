!-------------------------------------------------------------------------------
! initial magnetic field 
!@brief In this subroutine, we aim at defining the initial topology of the magnetic field  
!
!> @author
!>  Clara Dehman
!>  Daniele Viganò
!-------------------------------------------------------------------------------

module initial_magnetic

  use grid, only: nang, nr, r, vol
  use grid, only: br, bxi, beta
  use grid, only: y_lm
  use grid, only: curl_fnvol, fghost, dot_prod
  use magnetic_evolution, only: magnetic_bc
  use input_params, only: bpol_init, btor_init, pol_lm, tor_lm, LMAX_IN

  contains
  
  !---------------------------------------------------------------------------
  !! @brief initial magnetic topology in cubed sphere coordinates 
  !! constructed using the scalar functions
  ! 
  !>> Choice of the radial functions: The requiered condition is the matching with the 
  ! vacuum outside the star. 
  !> Funa function adopted from aguilera et al. 2008, is used for the dipolar poloidal part since 
  ! it is a solution for a dipolar field. This function is responsible on providing a smooth 
  ! matching between the interior radial field and the potential boundary conditions 
  !> A radial function that confine the field to the crust of a neutron star is used for the 
  ! toroidal scalar function psi and for the high multipoles, l > 1, in the poloidal scalar function 
  !
  ! phi_sf is defined up to nr. However, phi_sf is needed at nr+1 to compute ator (axi, aeta) at
  ! nr+1. Or ator at nr+1 is used to compute bxi and beta at nr. But at nr, we impose some averaging
  ! to treat the odd-even decoupling. Therefore, phi_sf and ator at nr+1 are not used.
  !
  ! Note: the phi and psi scalar functions are written as 1/r*sum_{l,m} f(r)*Ylm. 
  ! Then we apply curl(phi vec{r}) and curl(psi vec{r}), with vec{r} = r*er. 
  ! Thus, one can drop the r factor, and what is defined in the code is r*phi and r*psi 
  ! and not psi and phi scalar functions. 
  !
  ! Note: For the poloidal scalar function, we multiply the fr_crustconf by 10 in order to 
  ! increase the weights of the l>1 multipoles. If we remove the factor 10, the contribution 
  ! of the higher order multipoles (l>1) becomes negligible 
  !
  !! Authors:
  !!  Clara Dehman
  !!  Daniele Viganò
  !---------------------------------------------------------------------------   
  subroutine binit()

   implicit none

   ! Internal variables.
   integer i, p, l, m 
   real*8 mu, N1, N2, dummy
   real*8, dimension(0:nang+1, 0:nang+1, 1:6) :: fang_pol, fang_tor, fang_poldip
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: phi_sf, psi_sf
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: brpol, bxipol, betapol   
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: ar, axi, aeta 
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: bxitor, betator 
   real*8, dimension(0:nr+1) :: fr_crustconf

   phi_sf = 0d0 
   psi_sf = 0d0
   fang_pol = 0d0
   fang_tor = 0d0
   fang_poldip = 0d0
   brpol = 0d0
   bxipol = 0d0
   betapol = 0d0
   bxitor = 0d0
   betator = 0d0

! Poloidal and toroidal funtion for l=1 and m=0,+1,-1 
  fang_poldip(:,:,:) = fang_poldip(:,:,:) + pol_lm(1,0)*y_lm(:,:,:,1,0) + &
  & pol_lm(1,1)*y_lm(:,:,:,1,1) + pol_lm(1,-1)*y_lm(:,:,:,1,-1)
  fang_tor(:,:,:) = fang_tor(:,:,:) + tor_lm(1,0)*y_lm(:,:,:,1,0) + &
  & tor_lm(1,1)*y_lm(:,:,:,1,1) + tor_lm(1,-1)*y_lm(:,:,:,1,-1)

   do l=2,LMAX_IN
     fang_pol(:,:,:) = fang_pol(:,:,:) + pol_lm(l,0)*y_lm(:,:,:,l,0) ! symmetric part 
     fang_tor(:,:,:) = fang_tor(:,:,:) + tor_lm(l,0)*y_lm(:,:,:,l,0)
     do m=1,l ! asymmetric part 
       fang_pol(:,:,:) = fang_pol(:,:,:) + pol_lm(l,m)*y_lm(:,:,:,l,m) + pol_lm(l,-m)*y_lm(:,:,:,l,-m)
       fang_tor(:,:,:) = fang_tor(:,:,:) + tor_lm(l,m)*y_lm(:,:,:,l,m) + tor_lm(l,-m)*y_lm(:,:,:,l,-m)
     enddo
   enddo
  
   ! mu is the parameter used for the radial function of the initial poloidal field
   call getmu(r(nr),r(1),mu)
  ! This is the radial function used for the initial toroidal field and high order poloidal field (l>1)
   fr_crustconf(:) = - (r(nr)-r(:))**2*(r(:)-r(1))**2

   do i = 0, nr+1
    ! poloidal and toroidal scalar functions respectively
     phi_sf(i,:,:,:) = funa(mu*r(i),mu*r(nr))*fang_poldip(:,:,:) + fr_crustconf(i)*fang_pol(:,:,:) 
     psi_sf(i,:,:,:) = fr_crustconf(i)*fang_tor(:,:,:)
   enddo

  ! Poloidal magnetic field
  call curl_fnvol(phi_sf,0.*phi_sf,0.*phi_sf,ar,axi,aeta,1)
  call fghost(ar,axi,aeta)
  call curl_fnvol(ar,axi,aeta,brpol,bxipol,betapol,1)

  ! Toroidal magnetic field (brtor = 0 by definition, unused, use ar as a dummy)
  call curl_fnvol(psi_sf,0.*psi_sf,0.*psi_sf,ar,bxitor,betator,1) 

  ! Normalization
  ! Initialization
  N1 = 1d0
  N2 = 1d0

  ! Normalize the poloidal field to the maximum value of |Br| at the surface
  ! (corresponding to the polar value for a dipole)
  dummy = maxval(abs(brpol(nr,1:nang,1:nang,:)))
  if (dummy /= 0.) then
    N1 = bpol_init/dummy
  endif

  ! Normalize the toroidal field to the average root mean square value
  ! ar is used as dummy for the Btor^2 here
  ! TBD: This normalization is not working as expected...
  ar = 0d0
  call dot_prod(0.*bxitor,0.*bxitor,bxitor,bxitor,betator,betator,ar)
  dummy = 0d0
  ar = dsqrt(ar)
  do p=1,6
    dummy = dummy + sum(ar(2:nr-1:2,2:nang-1:2,2:nang-1:2,p)*vol(2:nr-1:2,2:nang-1:2,2:nang-1:2)) & 
    & /sum(vol(2:nr-1:2,2:nang-1:2,2:nang-1:2))
  enddo
  if (dummy /= 0.) then
    N2 = btor_init/dummy         
  endif

  ! Previous normalization
  ! ! Normalization (we can change it, just one choice!)
  ! if (brpol(nr,nang/2+1,nang/2+1,5) /= 0.) then
  !   N1 = bpol_init/brpol(nr,nang/2+1,nang/2+1,5)
  ! endif
  ! if (maxval(abs(bxitor(:,:,:,1:4))) /= 0.) then
  !   N2 = btor_init/maxval(abs(bxitor(:,:,:,1:4)))         
  ! endif

  br = brpol*N1
  bxi = bxipol*N1 + bxitor*N2
  beta = betapol*N1 + betator*N2

  call magnetic_bc(br,bxi,beta)

  call fghost(br,bxi,beta)

  end subroutine binit

  
!!-----------------------------------------------------------------------
!> @brief Definition of the radial function of poloidal field
!!        for crust-confined field, matching with potential solution
!!
!! It uses the Newton-Raphson method to solve the equation 
!! sin(mu*(R_core-R_ns)) - mu*R_core*cos(mu*(R_core-R_ns)) = 0
!! to obtain mu (eq.10 of Aguilera et al. 2008, A&A)
!!
!! Authors:
!!  Daniele Viganò
!-----------------------------------------------------------------------
subroutine getmu(rsurface,rcore,mu)
implicit none
real*8, intent(in) :: rsurface,rcore
real*8, intent(out) :: mu
real*8 dr, dfun, fun
  
dr = rcore - rsurface
! initial guess
mu = 2.d0
fun = 1.d0
do while (dabs(fun) > 1d-5)
  fun=dsin(mu*dr)-mu*rcore*dcos(mu*dr)
  dfun=dr*dcos(mu*dr)-rcore*dcos(mu*dr)+mu*rcore*dr*dsin(mu*dr)
  mu=mu-fun/dfun
enddo
end subroutine getmu


!!-----------------------------------------------------------------------
!> @brief Calculates the Bessel function A(x), linear combination of Bessel functions
!! as defined in equation (8) of Aguilera et al. (2008).
!! where x=mu*r (variable) and xr=mu*r(nr)
!! It allows a smooth matching (function and derivative)
!! with a non-relativistic dipole only
!!
!! Code owners:
!!  Daniele Viganò
!!-----------------------------------------------------------------------
real*8 function funa(x,xr)
implicit none
real*8, intent(in) :: x,xr
real*8 j1,n1

! Spherical Bessel functions j1 and n1 (fist and second kind).
j1=dsin(x)/x**2.d0-dcos(x)/x
n1=-dcos(x)/x**2.d0-dsin(x)/x
funa=x*(j1+dtan(xr)*n1)
return
end function funa


end module initial_magnetic

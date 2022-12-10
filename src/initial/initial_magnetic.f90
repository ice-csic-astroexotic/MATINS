!-------------------------------------------------------------------------------
! initial magnetic field 
!@brief In this subroutine, we aim at defining the initial topology of the magnetic field  
!
!> @author
!>  Clara Dehman
!>  Daniele Viganò
!-------------------------------------------------------------------------------

module initial_magnetic

  use grid, only: nang, nr, ievol, r, nangt
  use grid, only: br, bxi, beta
  use grid, only: lmax
  use grid, only: theta
  use grid, only: y_lm, dyth_lm, dyphi_lm
  use grid, only: curl_fnvol, fghost, f_spherical_to_cs
  use magnetic_evolution, only: magnetic_bc!, magnetic_bc_bessel

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
  !! Code owners:
  !!  Clara Dehman
  !---------------------------------------------------------------------------   
  subroutine binit(bpolmax,btormax)

    implicit none

    real*8, intent(in) :: bpolmax, btormax

  ! Internally variables.
   integer i, j, k, p, l, m 
  ! real*8 funa
   real*8 mu, N1, N2
 !  real*8 frel ! Relativistic factors  used in potentials' derivatives for potential solutions
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: phi_sf, psi_sf   
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: phir, phixi, phieta   
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: psir, psixi, psieta   
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: brpol, bxipol, betapol   
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: ar, axi, aeta 
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: brtor, bxitor, betator 
   real*8, dimension(0:nang+1, 0:nang+1, 1:6, 0:lmax, -lmax:lmax) :: dy2phi_lm
   real*8, dimension(0:nr+1) :: fr_crustconf


  ! In case they are useful 
  ! Second order phi derivative of y_lm defined from l:0->lmax and m:-l->l
  do l = 0, lmax
    do m = -l, l 
    if (m == 0) then 
     dy2phi_lm(:,:,:,l,0) = 0.d0
    else if(m > 0) then
     dy2phi_lm(:,:,:,l,m) = - m**2*y_lm(:,:,:,l,m)
    else 
    dy2phi_lm(:,:,:,l,m) = - abs(m)**2*y_lm(:,:,:,l,m)
    end if 
    end do 
    end do 


   call getmu(r(nr),r(1),mu)

   phi_sf = 0.d0 
   psi_sf = 0.d0 

   do i = 1, nr
    fr_crustconf(i) = - (r(nr)-r(i))**2*(r(i)-r(1))**2
    do p = 1, 6
      do j = 0, nang+1 
      do k = 0, nang+1
        phi_sf(i,j,k,p) = 0.5*funa(mu*r(i),mu*r(nr))*(y_lm(j,k,p,1,0) + &
        &  y_lm(j,k,p,1,1) + y_lm(j,k,p,1,-1) ) +  & 
        &  fr_crustconf(i)*(y_lm(j,k,p,2,1) + y_lm(j,k,p,3,2) + y_lm(j,k,p,3,3)+ y_lm(j,k,p,5,3)) 
        psi_sf(i,j,k,p) = fr_crustconf(i)*(y_lm(j,k,p,2,0) &
        & + y_lm(j,k,p,1,1) + y_lm(j,k,p,2,1) + y_lm(j,k,p,3,2) + y_lm(j,k,p,3,3) + y_lm(j,k,p,5,3))

        !Axysymmetric setup
        !phi_sf(i,j,k,p) = funa(mu*r(i),mu*r(nr))*y_lm(j,k,p,1,0)
        !psi_sf(i,j,k,p)=fr_crustconf(i)*(y_lm(j,k,p,2,0)+y_lm(j,k,p,3,0))
      enddo
      enddo
    enddo
   enddo
   
 ! Poloidal magnetic field
  phir = phi_sf
  phixi = 0.d0 
  phieta = 0.d0  
  call curl_fnvol(phir,phixi,phieta,ar,axi,aeta,1)
  call fghost(ar,axi,aeta)
  call curl_fnvol(ar,axi,aeta,brpol,bxipol,betapol,1)

  ! Toroidal magnetic field 
  psir = psi_sf
  psixi = 0.d0 
  psieta = 0.d0 

  call curl_fnvol(psir,psixi,psieta,brtor,bxitor,betator,1) 

  ! Normalization (we can change it, just one choice!)
  if (brpol(nr,nang/2+1,nang/2+1,5) /= 0.) then
    N1 = bpolmax/brpol(nr,nang/2+1,nang/2+1,5)
  endif
  if (maxval(abs(bxitor(:,:,:,1:4))) /= 0.) then
    N2 = btormax/maxval(abs(bxitor(:,:,:,1:4)))         
  endif

  br = brpol*N1
  bxi = bxipol*N1 + bxitor*N2
  beta = betapol*N1 + betator*N2

  call magnetic_bc(br,bxi,beta)
  call fghost(br,bxi,beta)

  end subroutine binit

  !Same as before but with axisymmetric setup

  subroutine binit_axi(bpolmax,btormax)

    implicit none

    real*8, intent(in) :: bpolmax, btormax

  ! Internally variables.
   integer i, j, k, p, l, m 
  ! real*8 funa
   real*8 mu, N1, N2
 !  real*8 frel ! Relativistic factors  used in potentials' derivatives for potential solutions
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: phi_sf, psi_sf   
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: phir, phixi, phieta   
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: psir, psixi, psieta   
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: brpol, bxipol, betapol   
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: ar, axi, aeta 
   real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: brtor, bxitor, betator 
   real*8, dimension(0:nang+1, 0:nang+1, 1:6, 0:lmax, -lmax:lmax) :: dy2phi_lm
   real*8, dimension(0:nr+1) :: fr_crustconf


  ! In case they are useful 
  ! Second order phi derivative of y_lm defined from l:0->lmax and m:-l->l
  do l = 0, lmax
    do m = -l, l 
    if (m == 0) then 
     dy2phi_lm(:,:,:,l,0) = 0.d0
    else if(m > 0) then
     dy2phi_lm(:,:,:,l,m) = - m**2*y_lm(:,:,:,l,m)
    else 
    dy2phi_lm(:,:,:,l,m) = - abs(m)**2*y_lm(:,:,:,l,m)
    end if 
    end do 
    end do 


   call getmu(r(nr),r(1),mu)

   phi_sf = 0.d0 
   psi_sf = 0.d0 

   do i = 1, nr
    fr_crustconf(i) = - (r(nr)-r(i))**2*(r(i)-r(1))**2
    do p = 1, 6
      do j = 0, nang+1 
      do k = 0, nang+1
        phi_sf(i,j,k,p) = funa(mu*r(i),mu*r(nr))*y_lm(j,k,p,1,0)
        psi_sf(i,j,k,p)=fr_crustconf(i)*(y_lm(j,k,p,2,0)+y_lm(j,k,p,3,0))
      enddo
      enddo
    enddo
   enddo
   
 ! Poloidal magnetic field
  phir = phi_sf
  phixi = 0.d0 
  phieta = 0.d0  
  call curl_fnvol(phir,phixi,phieta,ar,axi,aeta,1)
  call fghost(ar,axi,aeta)
  call curl_fnvol(ar,axi,aeta,brpol,bxipol,betapol,1)

  ! Toroidal magnetic field 
  psir = psi_sf
  psixi = 0.d0 
  psieta = 0.d0 

  call curl_fnvol(psir,psixi,psieta,brtor,bxitor,betator,1) 

  ! Normalization (we can change it, just one choice!)
  if (brpol(nr,nang/2+1,nang/2+1,5) /= 0.) then
    N1 = bpolmax/brpol(nr,nang/2+1,nang/2+1,5)
  endif
  if (maxval(abs(bxitor(:,:,:,1:4))) /= 0.) then
    N2 = btormax/maxval(abs(bxitor(:,:,:,1:4)))         
  endif

  br = brpol*N1
  bxi = bxipol*N1 + bxitor*N2
  beta = betapol*N1 + betator*N2

  call magnetic_bc(br,bxi,beta)
  call fghost(br,bxi,beta)

  end subroutine binit_axi


  !---------------------------------------------------------------------------
  !! @brief initial magnetic topology: Poloidal dipole + Toroidal Quadrupole 
  !! This initial topology is adopted from Aguilera et al 2008
  !
  !! Code owners:
  !!  Clara Dehman
  !--------------------------------------------------------------------------- 
  subroutine DP_TQ(bpolmax,btormax)
   implicit none 

  real*8, intent(in) :: bpolmax, btormax
  real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: bth, bphi

 ! Internally variables.
  integer i, j, k, p
 ! real*8 funa
  real*8 mu
  bth = 0.d0
  bphi = 0.d0
  br = 0.d0

 ! Ingredients for a crustal confined dipolar field 
 call getmu(r(nr),r(1),mu)

 do i = 1, nr
 do p = 1, 6
 do j = 0, nang+1
 do k = 0, nang+1
  br(i,j,k,p) = r(nr)**2.d0*dcos(theta(j,k,p))*funa(mu*r(i),mu*r(nr))/r(i)**2.d0 
  bth(i,j,k,p) = - r(nr)**2.d0*dsin(theta(j,k,p))/(2.d0*r(i)) & 
&     * (funa(mu*r(i+1),mu*r(nr)) - funa(mu*r(i-1),mu*r(nr)))/(r(i+1) - r(i-1))
 ! bphi(i,j,k,p) = mu*r(nr)**2.d0*dsin(theta(j,k,p))*funa(mu*r(i),mu*r(nr))/(2.d0*r(i))
  bphi(i,j,k,p) = r(nr)**2/r(i)*((r(nr)-r(i))**2*(r(i)-r(1))**2*dsin(theta(j,k,p))*dcos(theta(j,k,p)))
 end do
 end do
 end do
 end do 

! normalization: magnetic dipole + toroidal field 
bth = bth*bpolmax/maxval(br(ievol:nr,:,:,:))
br = br*bpolmax/maxval(br(ievol:nr,:,:,:))
bphi = bphi*btormax/maxval(bphi(ievol:nr,:,:,:))
! Spherical to CS
call f_spherical_to_cs(bth,bphi,bxi,beta,0)
! Application of BC
call magnetic_bc(br,bxi,beta)
call fghost(br,bxi,beta)

end subroutine

  !---------------------------------------------------------------------------
  !! @brief initial magnetic topology: Pure Toroidal Quadrupolar Field
  ! Taken from Viganò et al. 2012 
  !
  !! Code owners:
  !!  Clara Dehman
  !--------------------------------------------------------------------------- 
subroutine pureTQ(btormax)
  implicit none 

 real*8, intent(in) :: btormax
 real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: bth, bphi

  ! Internally variables.
 integer i, j, k, p

      bth = 0.d0
      bphi = 0.d0
      br = 0.d0
      do i = 1, nr
      do p = 1, 6
      do j = 0, nang+1
      do k = 0, nang+1
        bphi(:,j,k,p) = r(nr)/r*((r(nr)-r)**2.d0*(r-r(1))**2.d0*dsin(theta(j,k,p))*dcos(theta(j,k,p)))
      end do
      end do
      end do
      end do 
      bphi = bphi*btormax/maxval(bphi(ievol-1:nr,:,:,:))
      ! Spherical to CS
     call f_spherical_to_cs(bth,bphi,bxi,beta,0)
     ! Application of BC
      call magnetic_bc(br,bxi,beta)
      call fghost(br,bxi,beta)
end subroutine



  !---------------------------------------------------------------------------
  !! @brief initial magnetic topology: Bessel test 
  !  Adopted from Viganò et al. 2012 & Pons et al. 2019 
  !  Note: some corrections are implemented to both formalism  
  !
  !! Code owners:
  !!  Clara Dehman
  !  Note: This test is commented out because it need several allocated variables in the grid 
  ! To prevent filling the memory with unused things, I am commenting out this test 
  !--------------------------------------------------------------------------- 
subroutine bessel(bpolmax,btormax)

  implicit none 
  
  real*8, intent(in) :: bpolmax, btormax
 ! real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: bth, bphi

! Internally variables.
 ! integer i, j, k, p

    ! alpha = 2 ! wavenumber in length unit
    ! kbessel = alpha*r(nr) ! wavenumber in dimensionless unit
  !   bth = 0.d0
  !   bphi = 0.d0
  !   br = 0.d0
  ! !  jr_an = 0.d0
  ! !  jth_an = 0.d0
  ! !  jphi_an = 0.d0
     !   do i = 1, nr
          !     do p = 1, 6
          !       do j = 0, nang+1
          !         do k = 0, nang+1
          !       ! ----------- spherical Bessel functions: Pons et al. 2019 ------------------ 
          !       !   br(:,j,k,p) = bpolmax*r(nr)/r*dcos(theta(j,k,p))*(dsin(alpha*r)/(alpha*r)**2.d0 & 
          !       !  &  - dcos(alpha*r)/(alpha*r))     !*dsin(phi(j,k,p))  
          !       !   bth(:,j,k,p) = bpolmax*r(nr)/(2.d0*r)*dsin(theta(j,k,p))*(dsin(alpha*r)/(alpha*r)**2 &
          !       !  &  - dcos(alpha*r)/(alpha*r) - dsin(alpha*r)) !*dsin(phi(j,k,p))
          !       !   bphi(:,j,k,p) = kbessel*bpolmax/2.d0*dsin(theta(j,k,p))*(dsin(alpha*r)/(alpha*r)**2 & 
          !       !  &  - dcos(alpha*r)/(alpha*r)) !*dsin(phi(j,k,p))
          !       ! ----------- spherical Bessel functions: Daniele et al. 2012 ------------------ 
          !       !  br(:,j,k,p) = kbessel*bpolmax/(alpha*r)**2.d0*dcos(theta(j,k,p)) & 
          !       ! &  * (dsin(alpha*r)/(alpha*r) - dcos(alpha*r))  
          !       !   bth(:,j,k,p) = kbessel*bpolmax/(2.d0*(alpha*r)**2.d0)*dsin(theta(j,k,p)) &
          !       ! &  * (dsin(alpha*r)/(alpha*r) - dcos(alpha*r) & 
          !       ! &  - (alpha*r)*dsin(alpha*r)) 
          !       !   bphi(:,j,k,p) = bpolmax*alpha*r(nr)/(2.d0*(alpha*r))*dsin(theta(j,k,p)) &
          !       ! &  * (dsin(alpha*r)/(alpha*r) - dcos(alpha*r))  
          !       !   if (r(1) == 0) then
          !       !    br(0,j,k,p) = kbessel*bpolmax*dcos(theta(j,k,p))/3.d0
          !       !    bth(0,j,k,p) = - kbessel*bpolmax*dsin(theta(j,k,p))/3.d0
          !       !    bphi(0,j,k,p) = 0.d0
          !       !   endif
          !       !  jr_an(:,j,k,p) = alpha*br(:,j,k,p)
          !       !  jth_an(:,j,k,p) = alpha*bth(:,j,k,p)
          !       !  jphi_an(:,j,k,p) = alpha*bphi(:,j,k,p)
  !         end do
  !       end do
  !     end do
  !   end do 
  
  !  brin = br
  !  bxiin = bxi
  !  betain = beta
  
  ! test: Bessel test
  ! initialising the initial analytical solution
  !  bran = br
  !  bxian = bxi
  !  betaan = beta
  !  bthan = bth
  !  bphian = bphi
  
  ! Bessel magnetic boundary condition for the initial B field
  !  call magnetic_bc_bessel(0d0,br,bxi,beta)
  
  end subroutine 
  

  !---------------------------------------------------------------------------
  !! @brief initial magnetic topology: Whistler Test 
  !  Adopted from Justin Notes
  !
  !! Code owners:
  !!  Clara Dehman
  !--------------------------------------------------------------------------- 
subroutine whistler(bpolmax,btormax)
  implicit none 

 real*8, intent(in) :: bpolmax, btormax
 real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: bth, bphi

! Internally variables.
 integer i, j, k, p

      bth = 0.d0
      bphi = 0.d0
      br = 0.d0
      do i = 1, nr
      do p = 1, 6
      do j = 0, nang+1
      do k = 0, nang+1
      bth(i,j,k,p) = (r(nr)-r(i))**2*(r(i)-r(1))**2/(r(i)*dsin(theta(j,k,p)))*(dyphi_lm(j,k,p,1,0) & 
      !   & + dyphi_lm(j,k,p,1,1) + dyphi_lm(j,k,p,2,1)) & 
         & +  1.d-2*(dyphi_lm(j,k,p,13,11) - dyphi_lm(j,k,p,13,-11))) 
      bphi(i,j,k,p) = - (r(nr)-r(i))**2*(r(i)-r(1))**2/r(i)*(dyth_lm(j,k,p,1,0) & 
      !  & + dyth_lm(j,k,p,1,1) + dyth_lm(j,k,p,2,1)) &   
        & +  1.d-2*(dyth_lm(j,k,p,13,11) - dyth_lm(j,k,p,13,-11)))
      end do
      end do
      end do
      end do 
      bth = bth*btormax/maxval(bphi(ievol-1:nr,:,:,:))
      bphi = bphi*btormax/maxval(bphi(ievol-1:nr,:,:,:))
      ! Spherical to CS
     call f_spherical_to_cs(bth,bphi,bxi,beta,0)
     ! Application of BC
      call magnetic_bc(br,bxi,beta)
      call fghost(br,bxi,beta)

end subroutine



! TBD
subroutine hotspot

  use constants, only : PI

  implicit none
  !---------------- Hot Spot Profile Test ----------------------------
  real*8, parameter :: long = 90.d0*PI/180.0, lat = 45.d0*PI/180.0, ang = 20.0*PI/180.0
  real*8, parameter :: T_ext = 1d0
  real*8 :: Ts, deriv_envelope
  real*8 :: condition
  real*8, parameter :: t0 = 1d0, sigma_t = 1d-1, t1 = 1d1
  real*8, parameter :: theta0 = 0.5d0*PI, phi0 = 0.5d0*PI, tilt=0.d0*PI
!!-------------------------------------------------------------------


  integer p,jt,kt,j,k


  do p = 1, 6
    do jt = 1, nangt
      do kt = 1, nangt
        j = 2*jt
        k = 2*kt

      end do
    end do
  end do

  !           ! ------- Hot-Spot IC ------------------------

  !           !condition = dcos(phi(j,k,p))*dsin(theta(j,k,p))*dcos(long)*dsin(lat) + &
  !           !& dsin(phi(j,k,p))*dsin(theta(j,k,p))*dsin(lat)*dsin(long) + &
  !           !& dcos(theta(j,k,p))*dcos(lat)

            
            !if(condition .ge. dcos(ang)) temp(:, jt, kt, p) = T_init
            

  !           ! -------------------------------------------------------------

end subroutine hotspot

! TBD
subroutine test_pa

  implicit none

  integer p,jt,kt,j,k



!
!! --------------- Perez-Azorin 2006 Profile Test -------------------
!
!  real*8, parameter :: kappa_perp = 1.d0, omegatau = 1.d1
!  real*8 :: ang_term
!
!  ! variables for the generalized profile
!  real*8 :: cos_th_prime, beta_ang = 54.735610317245346d0*PI/180., gamma_ang = 45.d0*PI/180.
!    
!! ------------------------------------------------------------------
!
!  ! Parameters and variables for the rotating tilted magnetic field
!  integer, parameter :: nphase = 12, nrot = 2
!  real*8 phi_rot, kx, ky


!  bth = 0.d0
!  bphi = 0.d0
!  br = 0.d0
!  do p = 1, 6
!    do j = 0, nang+1
!      do k = 0, nang+1
!
!       ! --------- Perez-Azorin 2006 (uniform field in z direction) ----------
!    
!       !br(:, j, k, p) = dcos(theta(j, k, p))
!        !bth(:, j, k, p) = -sin(theta(j, k, p))
!
!        ! Generalized PA test
!        !br(:, j, k, p) = dsin(-beta_ang)*dcos(gamma_ang-0.5*PI)*dcos(phi(j, k, p))*dsin(theta(j, k, p)) + &
!        !&                dsin(-beta_ang)*dsin(gamma_ang-0.5*PI)*dsin(phi(j, k, p))*dsin(theta(j, k, p)) + &
!        !&                dcos(-beta_ang)*dcos(theta(j, k, p))
!
!        !bth(:, j, k, p) = dsin(-beta_ang)*dcos(gamma_ang-0.5*PI)*dcos(phi(j, k, p))*dcos(theta(j, k, p)) + &
!        !&                 dsin(-beta_ang)*dsin(gamma_ang-0.5*PI)*dsin(phi(j, k, p))*dcos(theta(j, k, p)) - &
!        !&                 dcos(-beta_ang)*dsin(theta(j, k, p))
!
!        !bphi(:, j, k, p) = - dsin(-beta_ang)*dcos(gamma_ang-0.5*PI)*dsin(phi(j, k, p)) + &
!        !&                    dsin(-beta_ang)*dsin(gamma_ang-0.5*PI)*dcos(phi(j, k, p))
!        ! ---------------------------------------------------------------------
!
!        ! Dipolar Field
!
!        br(:, j, k, p)  = 100.d0*2.d0*dcos(theta(j, k, p))/r(:)/r(:)/r(:)
!        bth(:, j, k, p) = 100.d0*dsin(theta(j, k, p))/r(:)/r(:)/r(:)
!        bphi(:, j, k, p) = 0.d0
!    
!      end do
!    end do
!  end do
!

  !   bth = 0.d0
    !   bphi = 0.d0
    !   br = 0.d0
    !   do p = 1, 6
    !   do j = 0, nang+1
    !   do k = 0, nang+1
    !    br(:, j, k, p) = dcos(theta(j, k, p))
    !    bth(:, j, k, p) = -dsin(theta(j, k, p))
    !   end do
    !   end do
    !   end do
    ! ! CD: I don't see the normalization of Perez-Azorin test 
      ! Spherical to CS - if you need them in cubed sphere coordinates 
    !  call f_spherical_to_cs(bth,bphi,bxi,beta,0)

  do p = 1, 6
    do jt = 1, nangt
      do kt = 1, nangt
        j = 2*jt
        k = 2*kt

            !do it = 1, nrt
                !i = 2*it - 1
                !ang_term = (dsin(theta(j, k, p))**2 + (dcos(theta(j, k, p))**2)/(1.d0 + omegatau*omegatau))
                !temp(it, jt, kt, p) = T_init*exp(-r(i)*r(i)*ang_term/(4*kappa_perp))

                !PA general

                !cos_th_prime = dsin(beta_ang)*dsin(theta(j, k ,p))*(-dcos(gamma_ang)*dcos(phi(j, k, p)) + &
                !&              dsin(gamma_ang)*dsin(phi(j, k, p))) + dcos(beta_ang)*dcos(theta(j, k, p))
                !ang_term = 1.d0 - cos_th_prime*cos_th_prime*omegatau*omegatau/(1.d0 + omegatau*omegatau)
                
                !temp(it, jt, kt, p) = T_init*exp(-r(i)*r(i)*ang_term/(4*kappa_perp))
                
            !enddo

      end do
    end do
  end do


end subroutine test_pa


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
  
!!-----------------------------------------------------------------------
!> @brief Definition of the radial function of poloidal field
!!        for crust-confined field, matching with potential solution
!!
!! It uses the Newton-Raphson method to solve the equation 
!! sin(mu*(R_core-R_ns)) - mu*R_core*cos(mu*(R_core-R_ns)) = 0
!! to obtain mu (eq.10 of Aguilera et al. 2008, A&A)
!!
!! Code owners:
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




end module initial_magnetic



!!----------------------------------------------------------------------------
!> @brief Applies a rotation around z-axis of an angle gamma and around the y 
!! axis of an angle beta
!!----------------------------------------------------------------------------

!subroutine field_rotation(br, bth, bphi, beta_ang, gamma_ang)
!
!  use constants, only: PI
!  use grid, only: nr, nang
!  use grid, only: theta, phi
!
!  real*8, allocatable() :: br_serv(), bth_serv(), bphi_serv()
!
!  br_serv = br
!  bth_serv = bth
!  bphi_serv = bphi 
!
!  br(:,:,:,:) = dsin(-beta_ang)*dcos(gamma_ang-0.5*PI)*dcos(phi(:,:,:))*dsin(theta(:,:,:)) + &
!  &                dsin(-beta_ang)*dsin(gamma_ang-0.5*PI)*dsin(phi(:,:,:))*dsin(theta(:,:,:)) + &
!  &                dcos(-beta_ang)*dcos(theta(:,:,:))
!
!  bth(:,:,:,:) = dsin(-beta_ang)*dcos(gamma_ang-0.5*PI)*dcos(phi(:,:,:))*dcos(theta(:,:,:)) + &
!  &                 dsin(-beta_ang)*dsin(gamma_ang-0.5*PI)*dsin(phi(:,:,:))*dcos(theta(:,:,:)) - &
!  &                 dcos(-beta_ang)*dsin(theta(:,:,:))
!
!  bphi(:,:,:,:) = - dsin(-beta_ang)*dcos(gamma_ang-0.5*PI)*dsin(phi(:,:,:)) + &
!  &                    dsin(-beta_ang)*dsin(gamma_ang-0.5*PI)*dcos(phi(:,:,:))
!end subroutine 

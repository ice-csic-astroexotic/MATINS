module microphysics

use constants, only: PI, UNIT_EN, UNIT_TIME, UNIT_R, UNIT_T, UNIT_B
use grid, only: etab, fh
use grid, only: r, nr, nrt, nangt, rmax, rmin
use grid, only: kappa_perp_arr, omegatau_arr, cv
use grid, only: tem0, bm, cv_core, cv_core_tot, cv_core_tot_der
use grid, only: T_core, enu
use grid, only: q_neutrino, q_neutrino_der, qnu_core_tot

implicit none

contains

 subroutine analytical_microphysics(profile)
   implicit none

   character(len=7), intent(in) :: profile
   integer :: i, j, k, p, ic, jc, kc 
   real*8 :: fh0
   !real*8 :: k_unit_conversion = 3.1375d32
   real*8 :: me = 9.1094d-28 !electron mass in cgs, I should use efficent mass
   real*8 :: me_eff_corr 
   !real*8 :: kb = 1.3807e-16 !Boltzmann constant in cgs 
   real*8 :: e_charge = 4.8032d-10 !electron charge in cgs 
   real*8 :: c_light = 2.99792458d10
   real*8 :: tau0 = 9.9d-19 !cgs units
   real*8 :: tau !9.9d-19 !cgs units, from de Grandis et al. 2021
   real*8 :: mu, cv0, kappa_0 ! Chemical potential, normalization
   real*8 :: tem, t9

   fh0 = 9.79d-7 !fh for a density of n0 = 0.16 fm^-3 (see Aguilera eta al 2008)
 
   !--------------------------------------------------------------------------------------------------  
   ! Define profile for fh and eta

    select case(profile)
      case("uniform")
     !--------------------------------------------------------------------------------------------------
     ! Possible profiles for fh and etab used for the magnetic part
        fh(:) = 1.d1
        etab(:,:,:,:) = 1.d-1
     !--------------------------------------------------------------------------------------------------
      case("Vigan21")
        ! fh[km^2/Myr 1e12G] comes from Viganò et al. 2021 Fig. 1
        fh(:) = 0.011*exp(10.*(r(:)-r(0))**1.8)
        !> case 2:
        ! eta[km^2/Myr] is from the 2D and roughly corresponding to T=3e8 K (with quite low impurities):
        do i = 0,nr+1 
          etab(i,:,:,:) = 0.4*exp(5.*(r(i) - r(0))) !, a factor 10 smaller for T=1e8 k
        enddo
      case("Clara22")
        fh(0)= 0.d0
        etab(0,:,:,:) = 0.d0
        fh(1:nr) = 0.011*exp(10.*(r(1:nr)-r(1))**1.8)
          ! etab[km^2/Myr] providing a similar fit to that in case 2 (from the 2D and roughly corresponding to T=3e8 K),
          ! High impurity close to the inner B.C. which could reflect (pasta layer).
        do i = 1,nr
          etab(i,:,:,:) = 60.d0*(r(i)-r(1))**3.5/(r(nr)-r(1))**3.5 + 30.d0*(r(nr)-r(i))**4/(r(nr)-r(1))**4
        enddo

      end select 

    me_eff_corr = 1.d0 ! Effective mass
    tau = 1.d0 !normalized in unit of 9.9d-19 s (see de Grandis et al. 2021) !3.13887d-32 !9.9d-19 s in Myr
    cv0 = 3.5896874d18*(UNIT_R**3)*UNIT_T/UNIT_EN
    kappa_0 = 5.0154589650d18*UNIT_TIME*UNIT_R*UNIT_T/UNIT_EN 

    do p = 1, 6
      do k = 0, nangt+1
        do j = 0, nangt+1
          do i = 1, nrt 
            ! Define omegatau profile
            ! Units = adimensional 
            ic = 2*i -1
            jc = 2*j
            kc = 2*k

            t9 = 0.1*tem0(i,j,k,p)
            mu = (1.d0 + (1.d0 - r(ic)/rmax)/0.0463)**(4.0/3.0) ! De Grandis fit

            ! DV: Not sure this is ok, since it's not the same magnetic omegatau (tau_e is different)
            !omegatau_arr(i,j,k,p) = fh(ic)*bm(ic,jc,kc,p)/etab(ic,jc,kc,p)
            
            omegatau_arr(i,j,k,p) = (e_charge/me/me_eff_corr/c_light)*bm(ic,jc,kc,p)*UNIT_B*tau*tau0

            ! Define heat capacity per unit volume profile
            ! Units = 10^40 erg/km^3/10^8 K

            !cv(:,:,:,:) = 5.4d-5*((fh(1:nr-1:2)/fh0)**(2.0/3.0))*temp(:,:,:,:) ! this is only the electron contribution
            !cv(i,j,k,p) = 54.0d0*((fh(ic)/fh0)**(2.0/3.0))*tem0(i,j,k,p) !Electron Cv


            cv(i,j,k,p) = cv0*tem0(i,j,k,p)*mu**2!De Grandis like heat capacity
          

            ! Define kappa_perp profile 
            ! Units = 10^40 erg/km/Myr/K

            !kappa_perp(:,:,:,:) = k_unit_conversion*PI*kb*kb*c_light*tau &
            !& *temp(:,:,:,:)/fh(1:nr-1:2)/e_charge/me_eff/(1 + omegatau(:,:,:,:)*omegatau(:,:,:,:))

            !kappa_perp_arr(i,j,k,p) = 1.07288d38*tem0(i,j,k,p) &
            !& *tau/fh(ic)/me_eff_corr/(1d0 + omegatau_arr(i,j,k,p)**2)

            ! The following formula is from De Grandis et al. 2021, eq. 5,
            ! kappa_0 = (pi*kb*c)**2 * n0 * T_0 * tau0/3/mu0
            ! n0 = 2.6d34, T0 = 1.d8, tau0 = 9.9d-19, mu0 = 2.9d-5

            kappa_perp_arr(i,j,k,p) = kappa_0*tem0(i,j,k,p)*tau*mu**2/(1d0 + omegatau_arr(i,j,k,p)**2)
            !kappa_perp_arr(i,j,k,p) = kappa_0/(1d0 + omegatau_arr(i,j,k,p)**2)


            ! neutrino emissivity in the crust

            ! this process emulates a neutrino process with emission at t9=1 
            ! equal to 10^20 erg/cm^3/s, with a dependence with temperature as 
            ! T^8

            q_neutrino(i,j,k,p) = 3.154d8*t9**8
            q_neutrino_der(i,j,k,p) = 8*3.154d8*t9**7

          enddo 
        enddo
      enddo
    enddo

    ! Integration of core quantities 

    ! Electron cv, need to change this
    tem = T_core/enu(1)
    !cv_core_tot = 54.0d0*((fh(1)/fh0)**(2.0/3.0))*tem*    &
    !& 4.d0*PI*rmin**3/3d0

    cv_core_tot = 1.d20*4.d0*PI*rmin**3/3d0*(UNIT_R**3)*UNIT_T/UNIT_EN!54.0d0*((fh(1)/fh0)**(2.0/3.0))*tem*    &
    cv_core_tot_der = 0.d0 !in de Grandis the core heat capacity does not depend on T
    !& 4.d0*PI*rmin**3/3d0
      
 end subroutine analytical_microphysics


 subroutine CoreCooling_Implicit(core_cooling_rhs, core_cooling_rhs_deriv)

  real*8, intent(out) :: core_cooling_rhs, core_cooling_rhs_deriv

  real*8 :: n0, np, fh0, emiss_0, t9_core

  ! Neutrino Cooling
    ! Core 
  !T_core = temp(1,1,1,1)
  n0 = 0.16d39 !reference density of 0.16 fm^-3
  np = 1.d37 !proton density in cgs 
  fh0 = 9.79d-7 !fh for a density of n0 = 0.16 fm^-3 (see Aguilera eta al 2008)
  t9_core = 0.1*T_core

  !Qnu_core = 2.52332d2*((np/n0)**(1.0/3.0))*(T_core**8.0) !MURCA in the core 
  !Qnu_deriv_core = 8.0*Qnu_core/T_core

  ! we now derfine the f(T) and f'(T) for the core 
  ! where f(T) = <e_nu>/<cv>, where <.> denotes the volume integration
  ! WARNING: we took them as constant so we drop the volume dependency 
  ! WARNING: We need to implement a numerical derivation to do that in 
  ! order to be genneral

  ! MURCA Processes (volumes cancels out)
  !core_cooling_rhs = -4.67281481*((np/n0)**(1./3.))*((fh(1)/fh0)**(2.0/3.0))* &
  !& (T_core/enu(1))**7.0 

  !core_cooling_rhs_deriv = -32.7097034*((np/n0)**(1./3.))*((fh(1)/fh0)**(2.0/3.0))* &
  !& (T_core/enu(1))**6.0 


  ! De Grandis (cv core does not depends on T)

  emiss_0 =  1.d21*4.d0*PI*rmin**3/3d0*(UNIT_R**3)*(UNIT_TIME)/UNIT_EN
  cv_core_tot = 1.d20*4.d0*PI*(rmin*UNIT_R)**3/3d0*UNIT_T/UNIT_EN
  cv_core_tot_der = 0.d0 

  qnu_core_tot = emiss_0*(t9_core/enu(1))**8*4*PI*rmin**3/3. !Just for output

  core_cooling_rhs = - (emiss_0/cv_core_tot)*(t9_core/enu(1))**8
  core_cooling_rhs_deriv = -8*(emiss_0/cv_core_tot)*0.1*(t9_core/enu(1))**7 !0.1 because of t9

 end subroutine CoreCooling_Implicit

end module microphysics

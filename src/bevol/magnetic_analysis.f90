module magnetic_analysis

  !> @brief Subroutine for the calculation of magnetic quantities useful for analysis and monitoring
  !!
  !! @param[in]     tyear_b    Time elapsed in years
  !! @param[in]     dtb       Timestep in years
  !!
  !!  Code owners:
  !!    Daniele Viganò
  !!------------------------------------------------------------------------------
  !! en_mag_magnetosphere = 2*pi*int(rns->infty,mu[-1,1]) (Br^2+Bth^2)/8*pi
  !      = rns^3 sum_l b_l^2 (l+1)/2(2l+1)
  ! (to be checked, and only for newtonian)

  ! Energy conservation: Pons et al. 2009, eq. 17
  !
  ! poynting: S=int(Surface) [(cExB)exp(2nu)/4PI] (E in the code contains an enu)
  ! The sign is defined as +div(S) (positive for outgoing flux)
  ! 
  ! cE*enu       [1e12 G*km/Myr = 3.1645569e-2 G km/s ]
  ! E            [3.1645569e-2 G km/s / CLIGHT]
  ! B            [1e12 G]
  ! Surface      [km**2]
  ! one factor exp(nu) is already contained in E
  !!------------------------------------------------------------------------------


  real*8, save :: en_mag_star, en_mag_star_tor, en_mag_magnetosphere
  real*8, save :: delta_en_mag_star, delta_en_mag_magnetosphere, helicity_star
  real*8, save :: j2_star, en_electric_star, divb_star_l2norm
  real*8, save :: q_joule_star, q_joule_shock_star
  real*8, save :: en_joule_star, en_joule_shock_star
  real*8, save :: poynting_star, poynting_star_tot
 
  contains
  !-----------------------------------------------------------------------

    subroutine analyse_magnetic_field(tyear_b,dtb)

      ! Module imports -------------------------------------------------------------
      use constants, only: UNIT_B, UNIT_E, UNIT_R, UNIT_JOULE, CLIGHT
      use constants, only: T_YEAR, PI
      use grid, only: kmax, lmax, nang, np, jevol
      use grid, only: benu
      use grid, only: aphi, br, bth, bphi, er, eth, ephi, j2
      use grid, only: rb, arear, areath, vol, jphi
      use grid, only: q_joule, q_joule_shock
      use legpol, only: nleg, blout
  
      implicit none
  
      ! Subroutine arguments -------------------------------------------------------
      real*8 dtb, tyear_b
        
      ! Local variables ------------------------------------------------------------
      integer :: i,j,k,l
      real*8, save :: en_mag_star_old, en_mag_magnetosphere_old
      real*8 :: facleg, divb_local
  
      ! ----------------------------------------------------------------------------
  
      en_mag_star = 0.d0
      en_mag_star_tor = 0.d0
      en_mag_magnetosphere = 0.d0
      helicity_star = 0.d0
      j2_star = 0.d0
      en_electric_star = 0.d0
      divb_star_l2norm = 0.d0
      q_joule_star = 0.d0
      q_joule_shock_star = 0.d0
      poynting_star = 0.d0
  
      ! Global quantities integrated in the volume.
         
      do i=1,nang
        do j=jevol,np
          ! Multiply by 1/4 since the the i and j run over the full grid (avoid quadruple counting)
          en_mag_star = en_mag_star + UNIT_B**2*UNIT_R**3/( 8d0*PI )*benu(j)*  &
       &              0.25*vol(i,j)*( br(i,j)**2 + bth(i,j)**2 + bphi(i,j)**2)
          en_mag_star_tor = en_mag_star_tor + UNIT_B**2*UNIT_R**3/( 8d0*PI )*benu(j)* &
       &              0.25*vol(i,j)*bphi(i,j)**2
    
          divb_local=( br(i,j+1)*arear(i,j+1) - br(i,j-1)*arear(i,j-1)  &
       &              +  bth(i+1,j)*areath(i+1,j) - bth(i-1,j)*areath(i-1,j) )
  
          divb_star_l2norm = divb_star_l2norm + 0.25d0*(UNIT_B*UNIT_R**2*divb_local)**2
    
          j2_star = j2_star + UNIT_B**2*UNIT_R*0.25d0*vol(i,j)*j2(i,j)/( 16d0*PI**2*benu(j) )
    
          en_electric_star = en_electric_star + UNIT_E**2*UNIT_R**3*0.25d0*vol(i,j)/benu(j)**2* &
       &            ( er(i,j)**2 + eth(i,j)**2 + ephi(i,j)**2 ) / (8.d0*PI)
    
          helicity_star = helicity_star + UNIT_B**2*UNIT_R**4*0.25d0*vol(i,j)*aphi(i,j)*bphi(i,j)
        enddo
        ! Calculation of the poynting flux out of the star (half of area to avoid double counting).
        poynting_star = poynting_star + CLIGHT*UNIT_E*UNIT_B*UNIT_R**2/(4d0*PI)* &
     &      benu(np)*(bphi(i,np)*eth(i,np) - bth(i,np)*ephi(i,np))*(0.5d0*arear(i,np))
    
      ! Calculation of magnetospheric energy according to Appendix A of Akgun et al. 2017
      ! (to be tested, and valid only in the Newtonian case)
!        en_mag_magnetosphere = en_mag_magnetosphere + UNIT_B**2*UNIT_R**3*rb(np)*(0.5d0*arear(i,np))* &
!     &      ( br(i,np)**2 - bth(i,np)**2 - bphi(i,np)**2 ) / (8.*PI)
      enddo

      do k=2,kmax
        i=2*k-2
        do l=2,lmax
          j=2*l-1
          q_joule_star = q_joule_star + 1.d40*vol(i,j)*q_joule(k,l)
          q_joule_shock_star = q_joule_shock_star + 1.d40*vol(i,j)*q_joule_shock(k,l)
        enddo
      enddo

      ! Calculation of the magnetospheric energy using Legendre polynomials
      ! (no relativistic factors are here)
      en_mag_magnetosphere = 0.
      do l=1,nleg
        facleg = (dble(l)+1d0)/(2d0*(2d0*dble(l)+1d0))
        en_mag_magnetosphere = en_mag_magnetosphere + &
     &       UNIT_B**2*UNIT_R**3*facleg*blout(l)**2*rb(np)**3
      enddo

      if (tyear_b == 0.) then
        en_joule_star = 0d0
        en_joule_shock_star = 0d0
        poynting_star_tot = 0d0
        delta_en_mag_magnetosphere = 0d0
        delta_en_mag_star = 0d0
      elseif (dtb /= 0.) then
        ! The following time-integration and time-derivative regard variation from previous loop
        ! Calculation of time-integrated (from 0 to current time) of Joule and Poynting energies
        en_joule_star = en_joule_star + q_joule_star*dtb*T_YEAR
        en_joule_shock_star = en_joule_shock_star + q_joule_shock_star*dtb*T_YEAR
        poynting_star_tot = poynting_star_tot + poynting_star*dtb*T_YEAR
  
        ! Calculation of the rate of internal and magnetospheric energy compared to the previous step
        delta_en_mag_magnetosphere = (en_mag_magnetosphere - en_mag_magnetosphere_old)/(dtb*T_YEAR)
        delta_en_mag_star = (en_mag_star - en_mag_star_old)/(dtb*T_YEAR)
      endif

      en_mag_magnetosphere_old = en_mag_magnetosphere
      en_mag_star_old = en_mag_star
    
      
    end subroutine analyse_magnetic_field


end module magnetic_analysis

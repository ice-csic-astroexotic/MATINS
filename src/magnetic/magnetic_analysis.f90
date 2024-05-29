module magnetic_analysis

!> @brief Subroutine for the calculation of magnetic quantities useful for analysis and monitoring
  !!
  !! @param[in]     tyear_b   Time elapsed for the magnetic field in years
  !! @param[in]     dtb       Timestep in years
  !! 
  !!  Code owners:
  !!    Clara Dehman
  !!    Daniele Viganò
  !!------------------------------------------------------------------------------

  use constants, only: UNIT_B, UNIT_E, UNIT_R, UNIT_EN, CLIGHT, T_YEAR, PI
  use grid, only: C, D, delta, enu, elambda
  use grid, only: nang, nr, nangt, nrt, ievol
  use grid, only: theta, phi, xi, eta, r
  use grid, only: vol, area_r, area_xi, area_eta, wint, lmax, espec_vol, espec_pol, espec_tor
  use grid, only: phi_scalar, psi_scalar
  use grid, only: br, bxi, beta, b2, j2, bm, fh, etab
  use grid, only: y_lm, dyth_lm, dyphi_lm
  use grid, only: br, bxi, beta, er, exi, eeta
  use grid, only: q_joule, f_cs_to_spherical
  use grid, only: en_joule_star_tot, poynting_star_tot, poynting_star_tot_surface, poynting_star_tot_interior 

  implicit none

  contains

  !!------------------------------------------------------------------------------
  !> @brief Subroutine compute_energy_balance
  !!
  !!  Code owners:
  !!    Daniele Viganò
  !!    Clara Dehman
  !!------------------------------------------------------------------------------
  subroutine compute_energy_balance(dtb)
  
  implicit none
  
  ! Subroutine arguments -------------------------------------------------------
    real*8, intent(in) :: dtb
    real*8 :: q_joule_star, poynting_star, poynting_star_surface, poynting_star_interior
  
  ! Local variables ------------------------------------------------------------
    integer :: i,j,k,p
    integer :: it, jt, kt
  
    poynting_star = 0.d0
    poynting_star_surface = 0.d0 
    poynting_star_interior = 0.d0
    q_joule_star = 0.d0

    ! Calculation of the poynting flux out of the star (we divide the area by 1/4 to avoid quadrupole counting). 
    do p= 1,6
      do j = 1, nang
        do k = 1, nang

          poynting_star = poynting_star + delta(j, k)**0.5d0/(C(j)*D(k)) & 
   & * CLIGHT*UNIT_E*UNIT_B*UNIT_R**2/(4d0*PI)*enu(nr)*wint(j,k)*0.25d0* (area_r(nr,j,k) & 
   & * (exi(nr, j, k, p)*beta(nr, j, k, p) - eeta(nr, j, k, p)*bxi(nr, j, k, p) )  &
   & - area_r(ievol,j,k)* (exi(ievol,j,k,p)*beta(ievol,j,k,p) - &
   &   eeta(ievol,j,k,p)*bxi(ievol,j,k,p) ))

   poynting_star_surface = poynting_star_surface + delta(j, k)**0.5d0/(C(j)*D(k)) & 
   & * CLIGHT*UNIT_E*UNIT_B*UNIT_R**2/(4d0*PI)*enu(nr)*wint(j,k)*0.25d0* (area_r(nr,j,k) & 
   & * (exi(nr, j, k, p)*beta(nr, j, k, p) - eeta(nr, j, k, p)*bxi(nr, j, k, p) ) )
   
    poynting_star_interior = poynting_star_interior + delta(j, k)**0.5d0/(C(j)*D(k)) & 
   & * CLIGHT*UNIT_E*UNIT_B*UNIT_R**2/(4d0*PI)*enu(nr)*wint(j,k)*0.25d0* ( &
   & - area_r(ievol,j,k)* (exi(ievol,j,k,p)*beta(ievol,j,k,p) - &
   &   eeta(ievol,j,k,p)*bxi(ievol,j,k,p) ) ) 

        end do 
      end do
    end do 

  ! The factor 1e40 enters in the thermal part in the 2D code, check if it's ok here...
    do p = 1,6 
      do jt = 2, nangt
        j = 2*jt ! corresponds to the center of the thermal cell 
        do kt = 2, nangt
          k = 2*kt
          do it = (ievol+2)/2,nrt
            i=2*it-1  
            q_joule_star = q_joule_star + vol(i,j,k)*q_joule(it,jt,kt,p)*UNIT_EN

          enddo
        enddo
      enddo
    enddo
  
    en_joule_star_tot = en_joule_star_tot + q_joule_star*dtb
    poynting_star_tot = poynting_star_tot + poynting_star*dtb*T_YEAR
    poynting_star_tot_surface = poynting_star_tot_surface + poynting_star_surface*dtb*T_YEAR
    poynting_star_tot_interior = poynting_star_tot_interior + poynting_star_interior*dtb*T_YEAR

  end subroutine compute_energy_balance


  !!------------------------------------------------------------------------------
  !> @brief Subroutine analyse_magnetic_field
  !!
  !!  Code owners:
  !!    Daniele Viganò
  !!    Clara Dehman
  !!------------------------------------------------------------------------------
  subroutine analyse_magnetic_field(b_avg, en_mag_star, divB_L2, j2_star, rey, rey_max, t_hall)
  
    implicit none
  
  ! Subroutine arguments -------------------------------------------------------
    real*8, intent(out) :: b_avg, en_mag_star, divb_L2, j2_star, rey, rey_max, t_hall
  
  ! Local variables ------------------------------------------------------------
    integer :: i, j, k, p
    real*8 :: voltot, divb_local
    real*8, dimension(ievol:nr) :: ratio
  
    en_mag_star = 0.d0
    j2_star = 0.d0
    divB_L2 = 0.d0
    voltot = 0.d0
    rey = 0.d0
    rey_max = 0.d0
    t_hall = 0.d0
    ratio = 0.d0 

    do i = ievol+1, nr - 1 
      do p= 1,6
        do j = 2, nang - 1
          do k = 2, nang - 1 

            if ( ( br(i,j,k,p) /= br(i,j,k,p)) .or. (beta(i,j,k,p) .ne. beta(i,j,k,p))  &
          &       .or. (bxi(i,j,k,p) .ne. bxi(i,j,k,p)) ) then
                print*,"[ERROR] NaN in magnetic field at i,j,k,p = ", i,j,k,p
                print*,"Br, Bxi, Beta = ", br(i,j,k,p), bxi(i,j,k,p), beta(i,j,k,p)
                stop
            endif  
    
           ! Multiply by 1/8=0.125 since the the i, j and k run over the full grid (avoid over counting)
            voltot = voltot + wint(j,k)*0.125d0*vol(i,j,k) 

            en_mag_star = en_mag_star + UNIT_B**2*UNIT_R**3/(8.d0*PI)*enu(i) & 
             & *wint(j,k)*0.125d0*vol(i,j,k)*b2(i,j,k,p) 

           ! we are contracting the contravariant components of the magnetic field with the covariant 
           ! components of the areas in order to get rid of the off-diagonal terms
           ! divb_local includes the integration in the local volume
            divb_local = (br(i+1,j,k,p)*area_r(i+1,j,k) - br(i-1,j,k,p)*area_r(i-1,j,k)  &
             &  + bxi(i,j+1,k,p)*area_xi(i,j+1,k) - bxi(i,j-1,k,p)*area_xi(i,j-1,k) & 
             &  + beta(i,j,k+1,p)*area_eta(i,j,k+1) - beta(i,j,k-1,p)*area_eta(i,j,k-1)) 
           
            divb_L2 = divb_L2 + UNIT_B**2*UNIT_R**4*wint(j,k)*0.125d0*divb_local**2
            j2_star = j2_star + UNIT_B**2*UNIT_R*wint(j,k)*0.125d0*vol(i,j,k)*j2(i,j,k,p)/( 16d0*PI**2*enu(i) )

            rey = rey + wint(j,k)*0.125d0*vol(i,j,k)*bm(i,j,k,p)*fh(i)/etab(i,j,k,p)
            t_hall = t_hall + wint(j,k)*0.125d0*vol(i,j,k)*bm(i,j,k,p)/(fh(i)*j2(i,j,k,p))
            
          end do
        end do 
      end do 
      ratio(i) = maxval(bm(i,2:nang-1,2:nang-1,1:6)/etab(i,2:nang-1,2:nang-1,1:6))*fh(i)
    end do

    b_avg = (en_mag_star*8d0*PI/(UNIT_R**3*voltot))**0.5d0
    rey = rey/voltot
    rey_max = maxval(ratio)
    t_hall = t_hall/voltot

  end subroutine analyse_magnetic_field


 !!------------------------------------------------------------------------------
 !> @brief Subroutine energy spectrum 
 !! The detailed derivation in Dehman et al. 2022.  
 !! 
 !! In this subroutine, we use the reconstructed radial scalar function in order to calculate the 
 !! spectral magnetic energy in the volume of a neutron star 
 !!
 !! Vaccuum magnetic energy in 3D [valid for r>=R]: (used for the spectral magnetic energy at the surface)
 !! espec_vol(l,:) = espec_vol(l,:) + espec_out(l,:)**2*(l+1)*(2*l+1)*r(i)**2*0.5d0*(r(i+1)-r(i-1)) 
 !!
 !!  Code owners:
 !!    Clara Dehman
 !!    Jose Antonio Pons
 !!------------------------------------------------------------------------------ 
  subroutine energy_spectrum(Espec_tot,Epol_tot,Etor_tot)
    implicit None
    integer :: i, l, m
    real*8, dimension(0:nr+1,0:lmax,-lmax:lmax) :: phir, dphir, psir
    real*8, dimension(0:nr+1,0:nang+1,0:nang+1,1:6) :: bth, bphi
    real*8, dimension(0:nr+1) :: Espec, Epol, Etor ! spectral magnetic energies at a given layer 
    real*8, intent(out) :: Espec_tot, Epol_tot, Etor_tot
    
    phi_scalar = 0.d0 
    psi_scalar = 0.d0
    espec_vol = 0.d0
    espec_pol = 0.d0
    espec_tor = 0.d0
    Espec = 0.d0 
    Epol = 0.d0
    Etor = 0.d0
    Espec_tot = 0.d0
    Epol_tot = 0.d0
    Etor_tot = 0.d0

    do i = ievol, nr
      call f_cs_to_spherical(bxi,beta,bth,bphi,i) 
      call get_rad_func(br(i,:,:,:),bth(i,:,:,:),bphi(i,:,:,:),i,phir(i,:,:),dphir(i,:,:),psir(i,:,:)) 
      do l = 1, lmax 
      do m = -l, l 

      ! Phi and Psi scalar functions
        phi_scalar(i,:,:,:) = phi_scalar(i,:,:,:) + 1.d0/r(i)*y_lm(:,:,:,l,m)*phir(i,l,m)
        psi_scalar(i,:,:,:) = psi_scalar(i,:,:,:) + 1.d0/r(i)*y_lm(:,:,:,l,m)*psir(i,l,m)

      ! Magnetic energy spectrum
      Espec(i) = Espec(i) + l*(l+1.d0)*(l*(l+1.d0)/r(i)**4*phir(i,l,m)**2 & 
      &            + dphir(i,l,m)**2/r(i)**2 + psir(i,l,m)**2/r(i)**2)*UNIT_B**2 ! magnetic energy in a layer "i"
 
      Epol(i) = Epol(i) + l*(l+1.d0)*(l*(l+1.d0)/r(i)**4*phir(i,l,m)**2 & 
      &            + dphir(i,l,m)**2/r(i)**2)*UNIT_B**2 ! poloidal energy in a layer "i"
 
      Etor(i) = Etor(i) + l*(l+1.d0)*psir(i,l,m)**2/r(i)**2*UNIT_B**2 ! toroidal energy in a layer "i"
 
      espec_vol(l,m) = espec_vol(l,m) + l*(l+1.d0)*(l*(l+1.d0)/r(i)**2*phir(i,l,m)**2 & 
      &            + dphir(i,l,m)**2 + psir(i,l,m)**2) & 
      &            * UNIT_B**2*0.5d0*(r(i+1)-r(i-1))*UNIT_R**3/(8.d0*PI)*enu(i)*elambda(i) 
 
      espec_pol(l,m) = espec_pol(l,m) + l*(l+1.d0)*(l*(l+1.d0)/r(i)**2*phir(i,l,m)**2 & 
      &            + dphir(i,l,m)**2)*UNIT_B**2*0.5d0*(r(i+1)-r(i-1))*UNIT_R**3/(8.d0*PI)*enu(i)*elambda(i)
 
      espec_tor(l,m) = espec_tor(l,m) + l*(l+1.d0)*psir(i,l,m)**2 & 
      &            * UNIT_B**2*0.5d0*(r(i+1)-r(i-1))*UNIT_R**3/(8.d0*PI)*enu(i)*elambda(i) 
 
      end do  
      end do  
       Epol_tot = Epol_tot + Epol(i)*0.5d0*r(i)**2*(r(i+1)-r(i-1))*UNIT_R**3/(8.d0*PI)*enu(i)*elambda(i)
       Etor_tot = Etor_tot + Etor(i)*0.5d0*r(i)**2*(r(i+1)-r(i-1))*UNIT_R**3/(8.d0*PI)*enu(i)*elambda(i)
       Espec_tot = Espec_tot + Espec(i)*0.5d0*r(i)**2*(r(i+1)-r(i-1))*UNIT_R**3/(8.d0*PI)*enu(i)*elambda(i)
     end do 

  end subroutine energy_spectrum  

 !!------------------------------------------------------------------------------
 !> @brief Subroutine get radial scalar functions in the volume of a neutron star 
 !! The detailed derivation in Dehman et al. 2022. 
 !!   
 !! In this subroutine, we reconstruct the radial scalar function for a given 
 !! field and a spherical harmonics decomposition
 !!
 !!  Code owners:
 !!    Clara Dehman
 !!    Jose Antonio Pons
 !!------------------------------------------------------------------------------ 
  subroutine get_rad_func(brin,bthin,bphiin,i,phir,dphir,psir)

    implicit none

    real*8, dimension(0:nang+1,0:nang+1,1:6), intent(in) :: brin, bthin, bphiin
    integer, intent(in) :: i
    real*8, dimension(0:lmax,-lmax:lmax) :: phir, dphir, psir
    real*8, dimension(0:nang+1,0:nang+1,1:6) :: dphir_yphi, psir_yphi, dphir_yth, psir_yth 
    real*8, dimension(0:nang+1,0:nang+1,1:6,0:lmax,-lmax:lmax) :: dyth
    integer :: l, m, p, c

    phir = 0.d0
    dphir = 0.d0
    psir = 0.d0

    do p = 1, 6    
    do l = 1, lmax
    do m = -l, l
      dyth(1:nang,1:nang,p,l,m) = dyth_lm(1:nang,1:nang,p,l,m)
    end do
    end do
    dphir_yphi(1:nang,1:nang,p) = bphiin(1:nang,1:nang,p)/(dsin(theta(1:nang,1:nang,p)) + 1.d-50)
    psir_yphi(1:nang,1:nang,p) = bthin(1:nang,1:nang,p)/(dsin(theta(1:nang,1:nang,p)) + 1.d-50)  
    dphir_yth(1:nang,1:nang,p) = bthin(1:nang,1:nang,p)
    psir_yth(1:nang,1:nang,p) = bphiin(1:nang,1:nang,p)
    end do

    ! For non-axisymmetric configuration, a non-negligeable contribution comes from the axis
    c = (nang+1)/2

    dphir_yphi(c,c,5:6) = 1d0/12d0*(dphir_yphi(c-1,c-1,5:6) + dphir_yphi(c+1,c+1,5:6)  & 
    & + dphir_yphi(c+1,c-1,5:6) + dphir_yphi(c-1,c+1,5:6) + 2d0*(dphir_yphi(c,c+1,5:6)  &
    & + dphir_yphi(c,c-1,5:6) + dphir_yphi(c-1,c,5:6) + dphir_yphi(c+1,c,5:6)) )

    psir_yphi(c,c,5:6) = 1d0/12d0*(psir_yphi(c-1,c-1,5:6) + psir_yphi(c+1,c+1,5:6)  & 
    & + psir_yphi(c+1,c-1,5:6) + psir_yphi(c-1,c+1,5:6) + 2d0*(psir_yphi(c,c+1,5:6)  &
    & + psir_yphi(c,c-1,5:6) + psir_yphi(c-1,c,5:6) + psir_yphi(c+1,c,5:6)) )

    dphir_yth(c,c,5:6) = 1d0/12d0*(dphir_yth(c-1,c-1,5:6) + dphir_yth(c+1,c+1,5:6)  & 
    & + dphir_yth(c+1,c-1,5:6) + dphir_yth(c-1,c+1,5:6) + 2d0*(dphir_yth(c,c+1,5:6)  &
    & + dphir_yth(c,c-1,5:6) + dphir_yth(c-1,c,5:6) + dphir_yth(c+1,c,5:6)) )

    psir_yth(c,c,5:6) = 1d0/12d0*(psir_yth(c-1,c-1,5:6) + psir_yth(c+1,c+1,5:6)  & 
    & + psir_yth(c+1,c-1,5:6) + psir_yth(c-1,c+1,5:6) + 2d0*(psir_yth(c,c+1,5:6)  &
    & + psir_yth(c,c-1,5:6) + psir_yth(c-1,c,5:6) + psir_yth(c+1,c,5:6)) )

    dyth(c,c,5:6,:,:) = 1d0/12d0*(dyth(c-1,c-1,5:6,:,:) + dyth(c+1,c+1,5:6,:,:)  & 
   & + dyth(c+1,c-1,5:6,:,:) + dyth(c-1,c+1,5:6,:,:) + 2d0*(dyth(c,c+1,5:6,:,:)  &
   & + dyth(c,c-1,5:6,:,:) + dyth(c-1,c,5:6,:,:) + dyth(c+1,c,5:6,:,:)) )

    !$OMP parallel do private(l,m,p) reduction(+: phir, dphir, psir)
    do p = 1,6
      do l = 1, lmax
        do m = -l, l
        
          phir(l,m) = phir(l,m) + 1.d0/(l*(l+1.d0)) * sum(0.25d0*area_r(i,1:nang,1:nang)*wint(1:nang,1:nang) & 
        &      * y_lm(1:nang,1:nang,p,l,m)*brin(1:nang,1:nang,p))
     
          dphir(l,m) = dphir(l,m) + 1.d0/(l*(l+1.d0)*r(i)) * sum(0.25d0*area_r(i,1:nang,1:nang)*wint(1:nang,1:nang)  &
        &      * (dyth(1:nang,1:nang,p,l,m)*dphir_yth(1:nang,1:nang,p) &
        &      + dyphi_lm(1:nang,1:nang,p,l,m)*dphir_yphi(1:nang,1:nang,p) ))

          psir(l,m) = psir(l,m) + 1.d0/(l*(l+1.d0)*r(i)) * sum(0.25d0*area_r(i,1:nang,1:nang)*wint(1:nang,1:nang)  &
        &      * (dyphi_lm(1:nang,1:nang,p,l,m)*psir_yphi(1:nang,1:nang,p) & 
        &      - dyth(1:nang,1:nang,p,l,m)*psir_yth(1:nang,1:nang,p) )) 

        end do
      end do 
    end do
    !$OMP end Parallel do
  end subroutine get_rad_func


end module magnetic_analysis

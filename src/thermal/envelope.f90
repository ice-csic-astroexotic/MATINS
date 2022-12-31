!> @brief This module provides the boundary conditions for the 
!>        temperature evolution (Tb-Ts relation, envelope model)
!>        For the explicit methos, sfluxb is the flux and cfluxb=0.
!>        For the implicit method, it must provide the flux minus the
!>        derivative of the flux w.r.t. the temperature times 
!>        the intensity (sfluxb) and the derivative of the flux (cfluxb)
!>        It contains corrections for magnetic envelopes.
!
subroutine envelope_model()

  use grid, only: nr, nrt, nang, nangt, g14, enu, br, bm, bpdip
  use grid, only: tem0, temp_surf, area_r
  use grid, only: bb_flux, cfluxb, sfluxb 
  use constants, only: STEFAN_BOLTZMANN, UNIT_R, UNIT_EN, UNIT_TIME
  use input_params, only: envelope
      
  implicit none

  ! external functions
  real*8 iron_env_B, light_env, light_mag_env, iron_env_B_PMG2009, iron_env_B_PY2001, iron_env_B_DV2013
  ! internal variables
  integer j,k,jt,kt,p
  real*8 tb9
  !real*8, dimension(1:nangt,1:nangt,1:6) :: temp_surf_incr
  real*8, dimension(0:nangt+1,0:nangt+1,1:6) :: temp_surf_incr
  real*8, parameter :: INCR_TEM = 1.01d0
  real*8 :: sb_constant = 1.78843596d11! Stef-Boltz const in sim. units !

  ! Initialize
  temp_surf = 0.d0
  bb_flux = 0.d0
  cfluxb = 0.d0
  sfluxb = 0.d0
   
  !-----------------------------------------------------------------------
  ! Gudmundsson 1983 (non-magntized)
  !-----------------------------------------------------------------------
  if (envelope == "Gudm") then
    temp_surf(1:nangt,1:nangt,1:6) = 1.d6*(g14**0.455*tem0(nrt,1:nangt,1:nangt,1:6)/1.288d0)**(1./1.82)   

  !-----------------------------------------------------------------------
  ! Accreted non-magnetic envelope.
  !-----------------------------------------------------------------------
  else if (envelope == "Accr") then
    
    do p=1,6
      do kt = 1,nangt
        do jt = 1,nangt
          tb9=0.1d0*tem0(nrt,jt,kt,p)
          temp_surf(jt,kt,p) = light_env(g14,tb9)
          temp_surf_incr(jt,kt,p) = light_env(g14,tb9*INCR_TEM)
        enddo
      enddo
    enddo

  !-----------------------------------------------------------------------
  ! Accreted magnetic envelope.
  !-----------------------------------------------------------------------
  else if (envelope == "Accr_mag") then

    do p=1,6
      do kt = 1,nangt
        do jt = 1,nangt
          tb9=0.1d0*tem0(nrt,jt,kt,p)
          temp_surf(jt,kt,p) = light_mag_env(g14,tb9,br(nr,j,k,p),bm(nr,j,k,p))
          temp_surf_incr(jt,kt,p) = light_mag_env(g14,tb9*INCR_TEM,br(nr,j,k,p),bm(nr,j,k,p))
        enddo
      enddo
    enddo  

  !-----------------------------------------------------------------------
  ! Magnetized envelopes (Potekhin, Pons & Page 2015, Appendix B).
  ! (A surface density of 10^10 g/cm^3 is assumed.)
  !-----------------------------------------------------------------------
  else if (envelope == "Iron_PPP15") then
    do p=1,6
      do kt = 1,nangt
        do jt = 1,nangt
          j = 2*jt
          k = 2*kt
          tb9=0.1d0*tem0(nrt,jt,kt,p)
          temp_surf(jt,kt,p) = iron_env_B(g14,tb9,bpdip,br(nr,j,k,p),bm(nr,j,k,p))
          temp_surf_incr(jt,kt,p) = iron_env_B(g14,tb9*INCR_TEM,bpdip,br(nr,j,k,p),bm(nr,j,k,p))
        enddo
      enddo
    enddo

  !-----------------------------------------------------------------------
  ! Magnetized envelopes (Pons 2009 envelope model).
  !-----------------------------------------------------------------------
    else if (envelope == "Iron_PMG09") then 
    do p=1,6
      do kt = 1,nangt
        do jt = 1,nangt
          j = 2*jt
          k = 2*kt
          tb9=0.1d0*tem0(nrt,jt,kt,p)
          temp_surf(jt,kt,p) = iron_env_B_PMG2009(g14,tb9,br(nr,j,k,p),bm(nr,j,k,p))
          temp_surf_incr(jt,kt,p) = iron_env_B_PMG2009(g14,tb9*INCR_TEM,br(nr,j,k,p),bm(nr,j,k,p))
        enddo
      enddo
    enddo

  !-----------------------------------------------------------------------
  ! Magnetized envelopes (Potekhin 2001 envelope model).
  !-----------------------------------------------------------------------
    else if (envelope == "Iron_PY01") then 
    do p=1,6
      do kt = 1,nangt
        do jt = 1,nangt
          j = 2*jt
          k = 2*kt
          tb9=0.1d0*tem0(nrt,jt,kt,p)
          temp_surf(jt,kt,p) = iron_env_B_PY2001(g14,tb9,br(nr,j,k,p),bm(nr,j,k,p))
          temp_surf_incr(jt,kt,p) = iron_env_B_PY2001(g14,tb9*INCR_TEM,br(nr,j,k,p),bm(nr,j,k,p))
        enddo
      enddo
    enddo

  !-----------------------------------------------------------------------
  ! Magnetized envelopes (Vigano 2013 envelope model).
  !-----------------------------------------------------------------------
  else if (envelope == "Iron_DV13") then 
    do p=1,6
      do kt = 1,nangt
        do jt = 1,nangt
          j = 2*jt
          k = 2*kt
          tb9=0.1d0*tem0(nrt,jt,kt,p)
          temp_surf(jt,kt,p) = iron_env_B_DV2013(g14,tb9,br(nr,j,k,p),bm(nr,j,k,p))
          temp_surf_incr(jt,kt,p) = iron_env_B_DV2013(g14,tb9*INCR_TEM,br(nr,j,k,p),bm(nr,j,k,p))
        enddo
      enddo
    enddo

  else
    print*, "<ERROR> Envelope: invalid name, ",envelope
  endif

  ! Blackbody flux
  ! bb_flux = STEFAN_BOLTZMANN*UNIT_R**2*1d-40*temp_surf**4
  ! Units: 1e40 erg/(Myr*K^4)
  ! SB CONSTANT:[erg/(cm^2*s*K^4)]
  ! bb_flux has to be in units 1e40 erg/Myr

  do p=1,6
    bb_flux(1:nangt, 1:nangt, p) = STEFAN_BOLTZMANN*(temp_surf(1:nangt, 1:nangt, p))**4*   &
&    area_r(nr,2:nang-1:2,2:nang-1:2)*UNIT_R**2*UNIT_TIME/UNIT_EN
  end do
  
  ! For explicit scheme:
  ! sfluxb = bb_flux 
  ! cfluxb(1:nangt, 1:nangt, 1:6) = 0d0

  ! Correction to make the source partially implicit,
  ! by linearizing first order derivative in T
  ! cfluxb: dF/dT
  ! sfluxb: Flux - (dF/dT)*T
  if (envelope == "Gudm") then
    ! For Gudmundsson, it is analytical
    cfluxb(1:nangt, 1:nangt, 1:6) = bb_flux(1:nangt, 1:nangt, 1:6)*enu(nr)*(4d0/1.82d0)/tem0(nrt,1:nangt,1:nangt,1:6)
  else

    ! For the other, we do it numerically
    !cfluxb(:,:,:) = bb_flux*((temp_surf_incr/temp_surf)**4 - 1d0)/((INCR_TEM-1d0)*(0.1*tem0(nr,1:nangt,1:nangt,:)))
    ! HERE the 0.1 at the denominator is not very intuitive. It comes from the fact that the increment 
    !is tb9*INCR_TEM, so it is done in unit of 10^9 and not in unit of 10^8
    
    cfluxb(1:nangt, 1:nangt, 1:6) = bb_flux(1:nangt, 1:nangt, 1:6)*enu(nr)*   &
    & ((temp_surf_incr(1:nangt, 1:nangt, 1:6)/temp_surf(1:nangt, 1:nangt, 1:6))**4  - 1d0)/  &
    &               ((INCR_TEM-1d0)*(0.1*tem0(nrt,1:nangt,1:nangt,1:6)))
  endif
  
  !sfluxb = bb_flux - cfluxb*tem0(nrt,1:nangt,1:nangt,:) !SA: I think it is not correct 


  ! Redshift corrections
  ! It's a bit messy this part, re-check everything.
  ! TBD: check this correction, considering it is used in output for L
  bb_flux(1:nangt, 1:nangt, 1:6) = bb_flux(1:nangt, 1:nangt, 1:6)*enu(nr)**2  

  !sfluxb = sfluxb*enu(nr)**2        ! 10^40 erg/km^2/s. !SA: I think it is not correct 

  
  sfluxb(1:nangt, 1:nangt, 1:6) = bb_flux(1:nangt, 1:nangt, 1:6) - &
  &  cfluxb(1:nangt, 1:nangt, 1:6)*tem0(nrt,1:nangt,1:nangt,1:6) !SA: This should be correct

  temp_surf = temp_surf*enu(nr)


end subroutine envelope_model
  
  
real*8 function iron_env_B(g14,tb9,bpdip,brad,b)

  implicit none
  real*8, intent(in) :: g14, tb9, bpdip, brad, b
  real*8 :: chip, chit, tFe4, t04, t14, tt, tp, aa1, aa2, tmax, costh
    
  t04 = (15.7d0*tb9**1.5d0+1.36d0*tb9)**1.5184d0
  t14 = 1.63d0*bpdip**0.476d0*tb9**(1.348d0/(1d0+1d-2*dsqrt(bpdip)))
  tFe4 = g14*(t14+(1d0+0.15d0*dsqrt(bpdip))*t04)
  tp = 1d6*tFe4**(0.25d0)
  if (bpdip == 0.d0) then
    iron_env_b = tp
  else  
    tmax = 1d6*(5.2d0*g14**0.65d0+0.093d0*dsqrt(g14*bpdip))
    tp = tp*(1d0+(tp/tmax)**4)**(-0.25d0)
    chit = 1d0/(1d0+(1230d0*tb9)**3.35d0*bpdip &
 &          *dsqrt(1d0+2d0*bpdip**2)  &
 &          /(bpdip+450d0*tb9+119d0*bpdip*tb9)**4  &
 &          +6.6d-3*bpdip**2.5d0/(dsqrt(tb9)+2.58d-3*bpdip**2.5d0))
    tt = chit*tp
    if (b /= 0.) then
      costh = dabs(brad/b)
    else
      costh = 1d0
    endif
    aa2 = 10d0*bpdip/(dsqrt(tb9)+0.1d0*bpdip*tb9**(-0.25d0))
    aa1 = aa2*dsqrt(tb9)/3d0
    chip = (1d0+aa1+aa2)*costh**2/(1d0+aa1*costh+aa2*costh**2)
    iron_env_b = tt+(tp-tt)*chip
  endif

end function iron_env_B

!-----------------------------------------------------------------------
! Light magnetised envelopes (Potekhin et al. 2003)
!
!> @author
! Clara Dehman 
!-----------------------------------------------------------------------
real*8 function light_mag_env(g14,tb9,brad,b) 

    implicit none
    real*8, intent(in) :: g14, tb9, b, brad
    real*8 :: tFe4, tacc4, corr, chip, chit, beta, chip1, chip2, alpha, zeta, cosph 

  ! non-magnetised temperature 
    zeta = tb9 - 1.d-3*g14**0.25*dsqrt(7*tb9)
    tFe4 = g14*((7.d0*zeta)**2.25+(0.33d0*zeta)**1.25)
    corr = 0.447d0 + 0.075d0*dlog10(1d9*tb9)/(1d0+(6.2d0*tb9)**4)
    tacc4 = (g14*(18.1d0*tb9)**2.42d0*corr+3.2d0*tb9**1.67d0*tFe4)/(1d0+3.2d0*tb9**1.67d0) 
 
  ! magnetised effect 
    beta = (1.d0 + 0.383d0*tb9**0.367d0)**(-1)
    chit = (1.d0+172.d0*b/(1.d0+155.d0*tb9**2.28))**0.5/(1.d0+383d0*b/(1.d0+94.d0*tb9**1.690))**beta
 
    chip1 = 1.d0 + (4.5d-3+0.055*tb9**2)/(tb9**2+0.0595*tb9**0.328)*(b**0.237)/(1.d0+6.8d-7*b/(tb9**2))**0.113 
    chip2 = 1.d0 + 1.d0/(3.7d0+(163d0+3.4d5*b**(-1.5)*tb9**2)) 
    chip = chip1*chip2**(-1)

    alpha = (2.d0 + chit/chip)**2

    if (b /= 0.) then
      cosph = dabs(brad/b) ! angle between the local magnetic field and the normal 
    else
      cosph = 1.d0
    endif

   ! considering a fully accreted envelope 
    light_mag_env = 1.d6*tacc4**0.25*(chip**alpha*cosph**2+chit**alpha*(1.d0-cosph**2))**(1./alpha)

end function light_mag_env

real*8 function light_env(g14,tb9)

  implicit none
  real*8, intent(in) :: g14, tb9
  real*8 :: tFe4, tacc4, corr

  tFe4 = g14*0.55*(10.*tb9)**2.4/(1.+0.9*tb9)
  corr = 0.447d0+0.075d0*dlog10(1d9*tb9)/(1d0+(6.2d0*tb9)**4)
  tacc4 = (g14*(18.1d0*tb9)**2.42d0*corr+3.2d0*tb9**1.67d0*tFe4)/(1d0+3.2d0*tb9**1.67d0)
  light_env = 1.d6*(tacc4)**(0.25)

end function light_env

real*8 function iron_env_B_DV2013(g14,tb9,brad,b) 

    implicit none

    real*8, intent(in) :: g14, tb9, brad, b 
    real*8 :: tFe4, tacc4, corr, chip, chit, tt, tp, cosph, logb, tpmax, ttmax, ts0

    tFe4 = g14*0.55*(10.*tb9)**2.4/(1.+0.9*tb9)
    corr = 0.447d0+0.075d0*dlog10(1d9*tb9)/(1d0+(6.2d0*tb9)**4)
    tacc4 = (g14*(18.1d0*tb9)**2.42d0*corr+3.2d0*tb9**1.67d0*tFe4)/(1d0+3.2d0*tb9**1.67d0)
    ts0 = 1.d6*(tacc4)**(0.25)

	  logb=dlog10(b)

	  chip=1.d0+b**.5*(0.1*(10.*tb9)**(-0.27) &  
     &	-	0.081*(10.*tb9)**(-0.58)*logb  & 
     &	+	0.0149*(10.*tb9)**(-0.8)*logb**2)

	  chit=1.d0+b**.5*(-0.231-0.514*exp(-(10.*tb9)*0.138) &
     &		+ 0.675*exp(-(10.*tb9)*0.148)*logb &
     &		- 0.204*exp(-(10.*tb9)*0.193)*logb**2 ) &
     &		/(1.+0.0204*b)

	  !dcos2th=br(i,np)**2/bm(i,np)**2
      if (b /= 0.) then
        cosph = dabs(brad/b) ! angle between the local magnetic field and the normal 
      else
        cosph = 1d0
      endif

	  tp = chip*ts0
	  tt = chit*ts0

    if (b >= 0.) then
  ! with the limit due to neutrino emission from the crust
  ! tpmax and ttmax are eq. 32 of Pons et al 2009
     tpmax = 3.6d6*(1.d0+0.02*logb)
     ttmax = 2.8d6/(1.d0+0.6*logb)   
     tp =  (1.d0/tp**4+1.d0/tpmax**4)**(-0.25d0)
     tt = (1.d0/tt**4+1.d0/ttmax**4)**(-0.25d0)
    end if 

    iron_env_B_DV2013 = (tp**4.5*cosph**2+tt**4.5*(1.d0-cosph**2))**(1./4.5)

  end function iron_env_B_DV2013


!-----------------------------------------------------------------------
! Iron magnetised envelopes (Pons et al. 2009)
!
!> @author
! Clara Dehman 
!-----------------------------------------------------------------------
real*8 function iron_env_B_PMG2009(g14,tb9,brad,b) 

    implicit none

    real*8, intent(in) :: g14, tb9, brad, b 
    real*8 :: tacc4, zeta, chip, chit, tt, tp, cosph, logb, tpmax, ttmax, ts0

    zeta = tb9 - 1.d-3*g14**0.25*dsqrt(7.d0*tb9)
    tacc4 = g14*((7.d0*zeta)**2.25+(zeta/3.d0)**1.25)
    ts0 = 1.d6*(tacc4)**(0.25)

	  logb=dlog10(b)

	  chip = 1.d0+0.05d0*b**(0.25d0)/tb9**(0.24d0) 
    
	  chit= dsqrt(1.d0+0.07d0*b*(0.03d0+tb9)**(-0.559d0)) / &
      &          (1.d0+0.9d0*b/(0.03d0+tb9))**(0.4d0)
    
	  !dcos2th=br(i,np)**2/bm(i,np)**2
      if (b /= 0.) then
        cosph = dabs(brad/b) ! angle between the local magnetic field and the normal 
      else
        cosph = 1.d0
      endif

	  tp = chip*ts0
	  tt = chit*ts0

    if (b >= 0.) then
  ! with the limit due to neutrino emission from the crust
  ! tpmax and ttmax are eq. 32 of Pons et al 2009
     tpmax = 3.6d6*(1.d0+0.02*logb)
     ttmax = 2.8d6/(1.d0+0.6*logb)   
     tp =  (1.d0/tp**4+1.d0/tpmax**4)**(-0.25d0)
     tt = (1.d0/tt**4+1.d0/ttmax**4)**(-0.25d0)
    end if 

    iron_env_B_PMG2009 = (tp**4.5*cosph**2+tt**4.5*(1.d0-cosph**2))**(1./4.5)

  end function iron_env_B_PMG2009

!-----------------------------------------------------------------------
! Iron magnetised envelopes (Potekhin et al. 2001)
!
!> @author
! Clara Dehman 
!-----------------------------------------------------------------------
  real*8 function iron_env_B_PY2001(g14,tb9,brad,b) 

    implicit none

    real*8, intent(in) :: g14, tb9, brad, b 
    real*8 :: tacc4, zeta, chip, chit, tt, tp, cosph, logb, tpmax, ttmax, ts0

    zeta = tb9 - 1.d-3*g14**0.25*dsqrt(7.d0*tb9)
    tacc4 = g14*((7.d0*zeta)**2.25+(zeta/3.d0)**1.25)
    ts0 = 1.d6*(tacc4)**(0.25)

	  logb=dlog10(b)

	  chip = 1.d0+0.0492d0*b**(0.292d0)/tb9**(0.24d0)
    
	  chit= dsqrt(1.d0+0.1076d0*b*(0.03d0+tb9)**(-0.559d0)) / &
      &          (1.d0+0.819d0*b/(0.03d0+tb9))**(0.6463d0)
    
	  !dcos2th=br(i,np)**2/bm(i,np)**2
      if (b /= 0.) then
        cosph = dabs(brad/b) ! angle between the local magnetic field and the normal 
      else
        cosph = 1.d0
      endif

	  tp = chip*ts0
	  tt = chit*ts0

    if (b >= 0.) then
  ! with the limit due to neutrino emission from the crust
  ! tpmax and ttmax are eq. 32 of Pons et al 2009
     tpmax = 3.6d6*(1.d0+0.02*logb)
     ttmax = 2.8d6/(1.d0+0.6*logb)   
     tp =  (1.d0/tp**4+1.d0/tpmax**4)**(-0.25d0)
     tt = (1.d0/tt**4+1.d0/ttmax**4)**(-0.25d0)
    end if 

    iron_env_B_PY2001 = (tp**4.5*cosph**2+tt**4.5*(1.d0-cosph**2))**(1./4.5)

  end function iron_env_B_PY2001
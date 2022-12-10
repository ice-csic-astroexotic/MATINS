!-------------------------------------------------------------------------------
!> @brief This subroutine computes the heat capacity. In the core we use 
!> semi-analytical fits for each component (protons, neutrons, electrons)
!> and inthe crust we call EOSMAG (potekhin)
!
!
subroutine compute_heat_capacity

  use grid, only: rtot ! Debug only
  use grid, only: rho, ne, nn, npr, nmu, effmn, effmp
  use grid, only: xh, ye, yn, yp, ymu, aa, zz, tccru, tcn, tcp  
  use grid, only: nrt, nangt, ncore, T_core, enu_core, vol_shell
  use grid, only: tem0, bm, cv, cv_core, cv_core_tot, cv_core_tot_der
  use constants, only: PI, UNIT_R, UNIT_T, UNIT_EN
  implicit none 
  integer i, j, k, p, it, jt, kt, is
  real*8 tem, cv_cgs, taucru, taun, taup
  real*8 cv_cgs1 !cv_cgs calculated at tem*INCR_TEM

  character(len=14), parameter :: CV_FORMAT = "(7es11.3)"
  real*8, parameter :: INCR_TEM = 1.01d0
  real*8, parameter :: UNIT_CV = UNIT_R**3*UNIT_T/UNIT_EN !! erg/K/cm3  ---> 10**40 erg/km^3/10**8 K



  cv = 0d0
  cv_core = 0d0
  cv_core_tot = 0d0
  cv_core_tot_der = 0d0

  do i = 1, ncore

    tem = T_core/enu_core(i)

!      tau = Tcrit/T, it tauX>1 SF corrections are not applied
    taucru = tem/dmax1(tccru(i),0.1d0*tem) 
    taun = tem/dmax1(tcn(i),0.1d0*tem)
    taup = tem/dmax1(tcp(i),0.1d0*tem)
    call cv_ionsenp(tem, rho(i), ne(i), nn(i), npr(i), nmu(i), effmn(i), effmp(i), &
    &             ye(i), yn(i), yp(i), ymu(i), xh(i), aa(i), zz(i), &
    &             taucru, taun, taup, 0d0, cv_cgs)

    ! Recalculate at tem*INCR_TEM for the derivative
    call cv_ionsenp(tem*INCR_TEM, rho(i), ne(i), nn(i), npr(i), nmu(i), effmn(i), effmp(i), &
    &             ye(i), yn(i), yp(i), ymu(i), xh(i), aa(i), zz(i), &
    &             taucru, taun, taup, 0d0, cv_cgs1)

    cv_core(i) = cv_cgs*UNIT_CV
    cv_core_tot = cv_core_tot + vol_shell(i)*cv_core(i)
    cv_core_tot_der = cv_core_tot_der + vol_shell(i)*(cv_cgs1-cv_cgs)*UNIT_CV/((INCR_TEM-1.d0)*tem)
  end do

  do it=1,nrt
   do p=1,6
    do jt=1,nangt
     do kt=1,nangt

      i = it*2 - 1
      is = i + ncore
      j = jt*2
      k = kt*2
      tem = tem0(it,jt,kt,p)    ! in 10**8 K

!      tau = Tcrit/T, it tauX>1 SF corrections are not applied
      taucru = tem/dmax1(tccru(is),0.1d0*tem) 
      taun = tem/dmax1(tcn(is),0.1d0*tem)
      taup = tem/dmax1(tcp(is),0.1d0*tem)
      call cv_ionsenp(tem, rho(is), ne(is), nn(is), npr(is), nmu(is), effmn(is), effmp(is), &
     &             ye(is), yn(is), yp(is), ymu(is), xh(is), aa(is), zz(is), &
     &             taucru, taun, taup, bm(i,j,k,p), cv_cgs)
      cv(it,jt,kt,p) = cv_cgs*UNIT_CV
     enddo
    enddo
   enddo
  enddo

end subroutine compute_heat_capacity



! This function calculates the heat capacity
! Van Riper 1991, ApJ 449-462
subroutine cv_ionsenp(t8, rho, ne, nn, npr, nmu, effmn, effmp, &
 &  ye, yn, yp, ymu, xh, aa, zz, taucru, taun, taup, b12, cv_cgs)

  use constants, only : MASS_N, RHO_NUC
  implicit none
  real*8, intent(in) :: t8, rho, ne, nn, npr, nmu, effmn, effmp
  real*8, intent(in) :: ye, aa, zz, xh, yn, yp, ymu
  real*8, intent(in) :: taucru, taun, taup, b12
  real*8, intent(out) :: cv_cgs
  real*8 cve, cvion, cvn, cvp, cvmu
  real*8 t
  real*8 ccvfn, ccvfp, ccvfe, ccvfmu !,ccvfion
  ! variables connected to routine EOSMAG (Potekhin)
  real*8 rhonuc, temp, gamag, effmp_corr
  real*8, parameter :: UN_B12 = 425.438, UN_T6 = 0.3157746

  t = t8*1.d8
  cvn  = 0.d0
  cvp  = 0.d0
  cvion= 0.d0
  cve  = 0.d0
  cvmu  = 0.d0

  effmp_corr = effmp
  if(rho < 2.854d14) effmp_corr = 1.d0 ! TBD: Why?
  ! Heat capacity
  if (xh == 0.d0) then    ! core
    cvn = ccvfn(t, nn, effmn, xh, taun)
    cvp = ccvfp(t, npr, effmp_corr, taup)
    cve = ccvfe(t, rho, ne, ye)
    cvmu = ccvfmu(t, rho, nmu, ymu)
  else    ! crust
    cvn = ccvfn(t, nn, effmn, xh, taucru)
    cvp = 0.d0
!    cvion = ccvfion(t,rho,xh,ah,zh)
!    cve = ccvfe(t, rho, ne, ye)
    temp = 1.d-6*t/UN_T6
    gamag = b12*UN_B12
    rhonuc = rho*xh
    ! The ion and electron heat capacities are given by eosmag

    call EOSMAG(zz,aa,rhonuc,temp,gamag,cvion,cve)
  endif
  cv_cgs = cvion + cve + cvn + cvp + cvmu

  if (cv_cgs <= 0.d0) stop '<ERROR> heat_capacity.f90: Negative heat capacity'

end subroutine cv_ionsenp

! Neutron heat capacity
real*8 function ccvfn(t, nn, effmn, xh, tau)

  implicit none
  real*8, intent(in) :: t, nn, effmn, xh, tau
  real*8 xnn, xv, v1, v2, Rcv_n, RRcv_B, RRcv_A_cru

  xnn = (6.497229d-14/effmn)*(nn*1d39)**(1.d0/3.d0)
  ccvfn = 4.5507d11*effmn**2*xnn*dsqrt(xnn**2+1.d0)*t 
  
  ! Apply Superfluid corrections 
  Rcv_n = 1d0 
  if ( (tau <= 1d0) .and. (xh == 0d0) ) then  ! SF neutrons in the core
    xv = v2(tau) !3P_2 (nB)
    Rcv_n = RRcv_B(xv)        
  endif
  if ( (tau <= 1d0) .and. (xh > 0d0)) then  ! SF neutrons in the crust
    xv = v1(tau)    !1S_0 (nA)
    Rcv_n = RRcv_A_cru(xv)
  endif
  ccvfn=ccvfn*Rcv_n

end function ccvfn

!   Proton heat capacity
real*8 function ccvfp(t, npr, effmp, taup)

  implicit none
  real*8, intent(in) :: t, npr, effmp, taup
  real*8 xpp, Rcv_p, xv, v1, RRcv_A

  xpp = (6.497229d-14/effmp)*(npr*1d39)**(1.d0/3.d0)
  ccvfp = 4.5507d11*effmp**2*xpp*dsqrt(xpp**2+1.d0)*t 

  !  Apply Superfluid corrections
  Rcv_p=1.d0
  if (taup <= 1d0) then    !1S_0 protons in core
    xv=v1(taup)
    Rcv_p=RRcv_A(xv)
  endif
  ccvfp=ccvfp*Rcv_p

end function ccvfp

!   Electron heat capacity  (old simple formula)
real*8 function ccvfe(t, rho, ne, ye)
  use constants, only : PI
  implicit none
  real*8, intent(in) :: t, rho, ne, ye
  real*8 xe
  !xe = 0.010067d0*(rho*ye)**(1.d0/3.d0)
  xe=(1.0546d-27*(3*PI**2*(ne*1d39))**(1.d0/3.d0))/((9.109d-28)*2.99d10)
  ccvfe = (ne*1d39)*2.298976d-25*t*dsqrt(xe**2+1.d0)/xe**2 !erg/K/cm3
end function ccvfe

! Muon heat capacity  
real*8 function ccvfmu(t, rho, nmu, ymu)
  use constants, only : PI
  implicit none
  real*8, intent(in) :: t, rho, nmu, ymu
  real*8 xmu

  xmu = (1.0546d-27*(3*PI**2*(nmu*1d39))**(1.d0/3.d0))/((1.8835327d-25)*2.99d10) !hkF_mu/m_muc
  if(xmu==0) then
    ccvfmu=0
  else
    ccvfmu = 1.1088d-27*(nmu*1d39)*t*dsqrt(xmu**2+1.d0)/xmu**2 !erg/K/cm3 pi^2 n k^2 T / (mc^2) * sqrt(x2+1)/x2
  end if
end function ccvfmu

!   Ion heat capacity  (old simple formula)
real*8 function ccvfion(t, rho, nb, xh, ah, zh)

  implicit none
  real*8 t,rho,xh,ah,zh,nb
  real*8 nions,tdebye,gamma,fdebye

  nions=nb*1d39*xh/ah
  tdebye=3.48d3*dsqrt(rho)*ZH/AH
  gamma = (2.275d5*ZH**2*(rho*XH/AH)**(1.d0/3.d0))/t
  if (gamma < 1.d0) then  
    ccvfion=nions*4.143d-16/2.d0
  elseif(gamma < 150.d0) then
    ccvfion=nions*4.143d-16*fdebye(t/tdebye) !erg/k/cm3
  else
    ccvfion=nions*(4.143d-16/2.d0)*(1.d0+dlog10(gamma)/dlog10(150.d0))
  end if

end function ccvfion

!   Debye function
real*8 function fdebye(x)
  real*8 x
  if (x <= 0.15d0) then
      fdebye = 77.9273d0*x**3
  elseif (x >= 0.4d0) then
      fdebye = 1.d0 - 1.d0/(20.d0*x**2)
  else
      fdebye = 1.69798d0*x+0.0083073d0
  endif

end function fdebye

! Superfluid/Superconducting suppression factors of the heat capacity
! Levenfish & Yakovlev (1994a)
! proton SF in the core
real*8 function RRcv_A(v)
  implicit none
  real*8 ccc,v,a,zexp1 
     ccc=1.d0/2.2736d0
  ccc=1.d0
  a=0.4186d0+dsqrt(1.007d0**2+(0.5010d0*v)**2)
  zexp1=1.456d0-dsqrt(1.456d0**2+v**2)
  if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff      
  RRcv_A=ccc*a**(2.5d0)*dexp(zexp1)  
end function RRcv_A

! neutron SF in the crust
real*8 function RRcv_A_cru(v) 
  implicit none
  real*8 ccc,v,a,zexp1
     ccc=1.d0/2.3573d0
  ccc=1.d0
  a=0.4186d0+dsqrt(1.007d0**2+(0.5010d0*v)**2)
  zexp1=1.456d0-dsqrt(1.456d0**2+v**2)
  if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff      
  RRcv_A_cru=ccc*a**(2.5d0)*dexp(zexp1)  
end function RRcv_A_cru

! neutron SF in the core
real*8 function RRcv_B(v)
  implicit none
  real*8 ccc,v,a,zexp1
     ccc=1.d0/2.167d0
  ccc=1.d0
  a=0.6893d0+dsqrt(0.790d0**2+(0.2824d0*v)**2)
  zexp1=1.934d0-dsqrt(1.934d0**2+v**2)
  if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff      
  RRcv_B=ccc*a**2*dexp(zexp1)  
end function RRcv_B

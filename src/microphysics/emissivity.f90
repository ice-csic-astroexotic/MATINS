! Compute all neutrino emissivities

subroutine compute_neutrino_emissivity()

  use grid, only: rho, ne, nn, npr, xh, ye, yn, yp, ymu, aa, zz
  use grid, only: kFe, kFn, qnu_core_tot, qnu_core_tot_der
  use grid, only: tccru, tcn, tcp, T_core, enu_core, vol_shell
  use grid, only: enu, vol
  use grid, only: effmn, effmp
  use grid, only: nrt, nangt, ncore
  use grid, only: tem0, bm, q_neutrino_core, q_neutrino, q_neutrino_der
  use grid, only: temp
  use grid, only: qnu_mur, qnu_nn, qnu_np, &
      &            qnu_pp, qnu_ep, qnu_cp_con,qnu_cp_cop, qnu_du, &
      &            qnu_ea, qnu_pl, qnu_syn, qnu_cp_cr,qnu_pa
  use constants, only: PI, UNIT_R, UNIT_TIME, UNIT_EN
  
  implicit none
  integer i,j,k,p,it,jt,kt,is
  real*8 T0, tem_floor, emis, emis1, taucru, taun, taup, qtot, qtotal, &
 & qmur,qnn,qnp,qpp,qep,qcp_con,qcp_cop,qdu,qea,qpl,qsyn,qcp_cr,qpa
  real*8, parameter :: INCR_TEM = 1.01d0
  real*8, parameter :: UNIT_NU = UNIT_R**3*UNIT_TIME/UNIT_EN  ! From cgs to 10**40 erg/(km^3*Myr) 

  character(len=14), parameter :: EMIS_FORMAT = "(15es11.3)"

  q_neutrino = 0d0
  q_neutrino_der = 0d0
  qnu_core_tot = 0d0
  qnu_core_tot_der = 0d0

  ! initialize different neutrino channel to 0
  qnu_mur = 0.0
  qnu_nn = 0.0
  qnu_np = 0.0
  qnu_pp = 0.0
  qnu_ep = 0.0
  qnu_cp_con = 0.0
  qnu_cp_cop = 0.0 
  qnu_du = 0.0
  qnu_ea = 0.0 
  qnu_pl = 0.0 
  qnu_syn = 0.0 
  qnu_cp_cr = 0.0
  qnu_pa = 0.0

  if (T_core > 1.d-2) then 
  ! Cool the core only if its temperature is above the floor 
    do i = 1, ncore

      T0 = T_core/enu_core(i)
      tem_floor = 0.1*T0     ! Physical temperature in 1e8 K
      taucru = T0/dmax1(tccru(i),tem_floor)
      taun   = T0/dmax1(tcn(i),tem_floor)
      taup   = T0/dmax1(tcp(i),tem_floor)

      ! Calculate the emissivities at the point
      emis  = qtot(T0,rho(i),ne(i),nn(i),npr(i), &
      &            kFe(i),kFn(i),effmn(i),effmp(i), &
      &            ye(i),yn(i),yp(i),ymu(i),xh(i),aa(i),zz(i), &
      &            taucru,taun,taup,0d0,qtotal,qmur,qnn,qnp, &
      &            qpp,qep,qcp_con,qcp_cop,qdu,qea,qpl,qsyn,qcp_cr,qpa) !In cgs

      ! Recalculate the emissivities at slightly larger T, to evaluate dQnu/dT
      emis1  = qtot(INCR_TEM*T0,rho(i),ne(i),nn(i),npr(i), &
      &            kFe(i),kFn(i),effmn(i),effmp(i), &
      &            ye(i),yn(i),yp(i),ymu(i),xh(i),aa(i),zz(i), &
      &            taucru,taun,taup,0d0,qtotal,qmur,qnn,qnp, &
      &            qpp,qep,qcp_con,qcp_cop,qdu,qea,qpl,qsyn,qcp_cr,qpa)
      

      q_neutrino_core(i) = emis*UNIT_NU


      qnu_core_tot = qnu_core_tot + vol_shell(i)*q_neutrino_core(i)*enu_core(i)**2

      !qnu_core_tot_der = qnu_core_tot_der +                              &
      !&           vol_shell(i)*(emis1 - emis)*UNIT_NU*enu_core(i)**2/((INCR_TEM-1d0)*T0) ! 10**40 erg * km^3/Myr/10**8 K

      ! new relativistic correction
      qnu_core_tot_der = qnu_core_tot_der +                              &
      &           vol_shell(i)*(emis1 - emis)*UNIT_NU*enu_core(i)/((INCR_TEM-1d0)*T0) ! 10**40 erg * km^3/Myr/10**8 K

      ! calculate all the separate neutrino contribution for the core
      
      qnu_mur(1) = qnu_mur(1) + qmur*vol_shell(i)*UNIT_NU*enu_core(i)**2
      qnu_nn(1) = qnu_nn(1) + qnn*vol_shell(i)*UNIT_NU*enu_core(i)**2
      qnu_np(1) = qnu_np(1) + qnp*vol_shell(i)*UNIT_NU*enu_core(i)**2 
      qnu_pp(1) = qnu_pp(1) + qpp*vol_shell(i)*UNIT_NU*enu_core(i)**2
      qnu_ep(1) = qnu_ep(1) + qep*vol_shell(i)*UNIT_NU*enu_core(i)**2 
      qnu_cp_con(1) = qnu_cp_con(1) + qcp_con*vol_shell(i)*UNIT_NU*enu_core(i)**2
      qnu_cp_cop(1) = qnu_cp_cop(1) + qcp_cop*vol_shell(i)*UNIT_NU*enu_core(i)**2
      qnu_du(1) = qnu_du(1) + qdu*vol_shell(i)*UNIT_NU*enu_core(i)**2
      qnu_ea(1) = qnu_ea(1) + qea*vol_shell(i)*UNIT_NU*enu_core(i)**2 
      qnu_pl(1) = qnu_pl(1) + qpl*vol_shell(i)*UNIT_NU*enu_core(i)**2 
      qnu_syn(1) = qnu_syn(1) + qsyn*vol_shell(i)*UNIT_NU*enu_core(i)**2 
      qnu_cp_cr(1) = qnu_cp_cr(1) + qcp_cr*vol_shell(i)*UNIT_NU*enu_core(i)**2
      qnu_pa(1) = qnu_pa(1) + qpa*vol_shell(i)*UNIT_NU*enu_core(i)**2


    end do
  endif

  do it=1,nrt
   do p=1,6
    do jt=1,nangt
     do kt=1,nangt

      i = it*2 - 1
      is = i + ncore
      j = jt*2
      k = kt*2

      T0 = tem0(it,jt,kt,p)
      tem_floor = 0.1*T0     ! Physical temperature in 1e8 K
      taucru = T0/dmax1(tccru(is),tem_floor)
      taun   = T0/dmax1(tcn(is),tem_floor)
      taup   = T0/dmax1(tcp(is),tem_floor)  

      ! Calculate the emissivities at the point
      emis = qtot(T0,rho(is),ne(is),nn(is),npr(is), &
     &            kFe(is),kFn(is),effmn(is),effmp(is), &
     &            ye(is),yn(is),yp(is),ymu(is),xh(is),aa(is),zz(is), &
     &            taucru,taun,taup,bm(i,j,k,p),qtotal, qmur,qnn,qnp, &
     &            qpp,qep,qcp_con,qcp_cop,qdu,qea,qpl,qsyn,qcp_cr,qpa)

      ! Recalculate the emissivities at slightly larger T, to evaluate dQnu/dT
      emis1 = qtot(T0*INCR_TEM,rho(is),ne(is),nn(is),npr(is), &
      &            kFe(is),kFn(is),effmn(is),effmp(is), &
      &            ye(is),yn(is),yp(is),ymu(is),xh(is),aa(is),zz(is), &
      &            taucru,taun,taup,bm(i,j,k,p),qtotal, qmur,qnn,qnp, &
      &            qpp,qep,qcp_con,qcp_cop,qdu,qea,qpl,qsyn,qcp_cr,qpa)
 
      ! Definition of the elements of matrix used in the cooling scheme
      ! q_neutrino = emis - (d emis/d(T*e^nu))*Te^nu  , given in 10**40 erg/km^3/s
 	    ! q_neutrino_der = d emis/dT (without redshift correction)
      q_neutrino(it,jt,kt,p) = emis*UNIT_NU  ! 10**40 erg/km^3/Myr
      q_neutrino_der(it,jt,kt,p) = (emis1-emis)/((INCR_TEM-1d0)*tem0(it,jt,kt,p))*UNIT_NU  ! 10**40 erg/(km^3*Myr*10**8 K)

      if (temp(it,jt,kt,p)<= 1.d-2) then
        ! no cooling if we are at the floor level
        q_neutrino(it,jt,kt,p) = 0.d0
        q_neutrino_der(it,jt,kt,p) = 0.d0
      endif 

       ! calculate all the separate neutrino contribution for the crust
       ! here the redshift is included 
      qnu_mur(2) = qnu_mur(2) + qmur*vol(i,j,k)*UNIT_NU*enu(i)**2
      qnu_nn(2) = qnu_nn(2) + qnn*vol(i,j,k)*UNIT_NU*enu(i)**2
      qnu_np(2) = qnu_np(2) + qnp*vol(i,j,k)*UNIT_NU*enu(i)**2 
      qnu_pp(2) = qnu_pp(2) + qpp*vol(i,j,k)*UNIT_NU*enu(i)**2
      qnu_ep(2) = qnu_ep(2) + qep*vol(i,j,k)*UNIT_NU*enu(i)**2 
      qnu_cp_con(2) = qnu_cp_con(2) + qcp_con*vol(i,j,k)*UNIT_NU*enu(i)**2
      qnu_cp_cop(2) = qnu_cp_cop(2) + qcp_cop*vol(i,j,k)*UNIT_NU*enu(i)**2
      qnu_du(2) = qnu_du(2) + qdu*vol(i,j,k)*UNIT_NU*enu(i)**2
      qnu_ea(2) = qnu_ea(2) + qea*vol(i,j,k)*UNIT_NU*enu(i)**2 
      qnu_pl(2) = qnu_pl(2) + qpl*vol(i,j,k)*UNIT_NU*enu(i)**2 
      qnu_syn(2) = qnu_syn(2) + qsyn*vol(i,j,k)*UNIT_NU*enu(i)**2 
      qnu_cp_cr(2) = qnu_cp_cr(2) + qcp_cr*vol(i,j,k)*UNIT_NU*enu(i)**2
      qnu_pa(2) = qnu_pa(2) + qpa*vol(i,j,k)*UNIT_NU*enu(i)**2



     enddo
    enddo
   enddo
  enddo

end subroutine compute_neutrino_emissivity

!  This function calls other functions to obtain the neutrino emissivity.
!  Reduction factors due to superfluidity are calculated throughout.
real*8 function qtot(t8,rho,ne,nn,npr,kFe,kFn,effmn,effmp,ye,yn,yp,ymu,xh,aa,zz,taucru,taun,taup,bf, &
&                   qtotal, qmur,qnn,qnp,qpp,qep,qcp_con,qcp_cop,qdu,qea,qpl,qsyn,qcp_cr,qpa)

  implicit none
  real*8, intent(in) :: t8,rho,ne,nn,npr,kFe,kFn
  real*8, intent(in) :: effmn,effmp,ye,yn,yp,ymu,xh,aa,zz,taucru,taun,taup,bf
  real*8, intent(out) :: qtotal, qmur,qnn,qnp,qpp,qep,qcp_con,qcp_cop,qdu,qea,qpl,qsyn,qcp_cr,qpa
  real*8 t,t9
  real*8 qeabrems,qplasma,qsynchro,qcpbf_cr,qpair
  real*8 qmurca,qdurca,qBnn,qBnp,qBpp,qepBrems,qcpbf_con,qcpbf_cop
  integer, parameter :: enable_durca = 0  ! TBD: Pass in the input?

  t = 1.d8*t8            ! t: temperature in K
  t9 = 1.d-1*t8          ! t9: temperature (10**9 K)

  ! Initialization. All Q are in [erg/cm3/sec]
  qmur= 0.d0
  qnn = 0.d0
  qnp = 0.d0
  qpp = 0.d0
  qep = 0.d0
  qcp_con = 0.d0
  qcp_cop = 0.d0
  qdu = 0.d0
  qeA = 0.d0
  qpl = 0.d0 
  qsyn = 0.d0
  qcp_cr = 0.d0
  qpa = 0.d0

  ! Core processes
  if ( yp > 0d0) then
    ! Modified URCA
    qmur = qmurca(t9, nn, npr, ye, yn, yp, ymu, effmn, effmp, taun, taup)
    ! ep Brems  (ep->ep+vv)
    qep = qepbrems(t9, npr, effmp)
    ! np and pp Brems (pp->pp+vv, np->np+vv)
    qnp = qbnp(t9, npr, effmn, effmp, taun, taup)
    qpp = qbpp(t9, npr, effmp, taup)   
    ! Direct Urca  
    if ( ye**(1d0/3d0) + yp**(1d0/3d0) > yn**(1d0/3d0) ) then  
      qdu = qdurca(t9, rho, effmn, effmp, ye, yp, taun, taup)
    end if

    ! Yakovlev et al. (2001) Physics Reports, Volume 354, Issue 1-2, p. 1-155.  (Y01)
    ! see comment pag. 64
    ! The threshold of muonic DURCA is different, but the emissivity is the same of the electronic DURCA
    ! Therefore, if active, you have twice the value in total
    if ( ymu**(1d0/3d0) + yp**(1d0/3d0) > yn**(1d0/3d0) ) qdu = qdu*2.

    ! Cooper Pairs of SF 3P_2 neutrons & SF 1S_0 protons
    if ( taun <= 1d0 )  &
    &    qcp_con = qcpbf_con(t9, nn, effmn, taun)
    if ( taup <= 1d0 )  &
    &    qcp_cop = qcpbf_cop(t9, npr, effmp, taup)
  endif

  ! Core and crust processes
  ! nn Bremsstrahlung (nn->nn+vv)
  qnn = qbnn(t9, nn, kFn, yn, xh, effmn, taucru, taun)
  ! Synchrotron
  qsyn = qsynchro(t, kFe, ye, bf)

  ! Crust processes
  if (xh > 0) then
    ! eA Bremsstrahlung (eA->eA+vv)
    qea  = qeabrems(t, rho, xh)
    ! plasma process
    qpl  = qplasma(t, kFe)
    if ( (taucru <= 1.d0) .and. (xh < 1.d0)) &  ! If in the inner crust and n are SF
 &       qcp_cr = qcpbf_cr(t9, nn, effmn, taucru)
    qpa    = qpair(t, rho, zz, aa)
  endif

  qtot = qmur + qnn + qnp + qpp + qep + qcp_con + qcp_cop + qdu + &
 &       qea + qpl + qsyn + qcp_cr + qpa

  qtotal=qtot

end function qtot


! n-n Bremstrahlung emissivity
! Friman & Maxwell (1979) ApJ 232, 541     (FM79)
! Yakovlev & Levenfish (1995) A&A 297, 717 (YL95)
real*8 function qbnn(t9, nn, kFn, yn, xh, effmn, taucru, taun)

  use constants, only: RHO_TO_N, RHO_NUC
  implicit none
  real*8, intent(in) :: t9, nn, kFn, effmn, yn, xh, taucru, taun
  real*8 rho0n, alphann, betann, nfnu, fmu, xmpi, xu
  real*8 xv, v1, v2, RBnn, RRnn_nA, RRnn_nB

  rho0n = RHO_TO_N*nn/RHO_NUC
  alphann = 0.59d0
  fmu = 1.d0
  if ( (xh > 0d0) .and. (xh < 1d0) ) then ! in the inner crust
    fmu = yn
    xmpi = 139.d0       ! pion mass in MeV
    xu = xmpi/2.d0/(kFn*197.326d0)  ! kFn in MeV
    alphann = 1.d0-3.d0*xu*datan(1.d0/xu)/2.d0+xu**2/2.d0/(1.d0+xu**2)
  endif           !-------------------------------------
  if (xh == 1d0) alphann=0.d0   ! in the outer crust----- 
  betann=0.56d0
  nfnu= 3.d0  ! number of neutrino flavor generated. 
  RBnn=1.d0
  qBnn = 7.4d19*effmn**4*rho0n**(1.d0/3.d0)*alphann*betann*3*T9**8*nfnu*fmu   !(YL95) 
!   qBnn = 7.335d19*effmn**4*(RHO0**(1.d0/3.d0))*T9**8 !(FM79) 

  ! Superfluidity correction (if T is below Tc)
  if ( (taun <= 1d0) .and. (xh == 0.d0) ) then  ! neutrons in the core
    xv = v2(taun)   ! 3P_2 (nB)  
    RBnn = RRnn_nB(xv)
  endif
  if ( (taucru <= 1d0) .and. (xh > 0d0) ) then  ! neutrons in the crust
    xv = v1(taucru)    ! 1S_0 (nA)
    RBnn = RRnn_nA(xv) 
  endif
  qBnn = qBnn*RBnn

end function qbnn

! n-p Bremstrahlung emissivity
! Friman & Maxwell (1979) ApJ 232, 541     (FM79)
! Yakovlev & Levenfish (1995) A&A 297, 717 (YL95)
real*8 function qbnp(t9, npr, effmn, effmp, taun, taup)

  use constants, only: RHO_TO_N, RHO_NUC
  implicit none
  real*8, intent(in) :: t9, npr, effmn, effmp, taun, taup
  real*8 rho0p, alphanp, betanp, nfnu
  real*8 taumax, taumin, xvn, xvp, v1, v2
  real*8 RBnp, RRnp_nB, RRnp_pA, RRnp_BA

  rho0p = RHO_TO_N*npr/RHO_NUC
  alphanp = 1.06d0
  betanp = 0.66d0
  nfnu = 3.d0
  RBnp = 1.d0

! qBnp = 3.15D20*EFFMN**2*EFFMP**2*(RHO0**(2./3.))*T9**8 !(FM79)
  qBnp = 1.5d20*(effmn*effmp)**2*rho0p**(1d0/3d0)*  &
 &       alphanp*betanp*nfnu*T9**8           !(YL95)

! Superfluidity correction
  taumax = max(taun,taup)
  taumin = min(taun,taup)
  if (taumin > 1.d0) go to 2200       ! no SF
  if (taumax <= 1.d0) then            ! n and p SF
    xvn = v2(taun)    ! n in 3P_2 
    xvp = v1(taup)    ! p in 1S_0
    RBnp = RRnp_BA(xvn,xvp) 
  elseif (taumax > 1.d0) then         ! only one component SF
    if (taun <= 1.d0) then            ! only n SF
      xvn = v2(taun)  ! 3P_2 (nB)  
      RBnp = RRnp_nB(xvn)
    elseif (taup <= 1.d0) then        ! only p in 1S_0
      xvp = v1(taup)  ! 1S_0 (pA)
      RBnp = RRnp_pA(xvp)
    endif
  endif

2200  continue

  qBnp = qBnp*RBnp

end function qbnp


! p-p Bremstrahlung emissivity
! Friman & Maxwell (1979) ApJ 232, 541     (FM79)
! Yakovlev & Levenfish (1995) A&A 297, 717 (YL95)
real*8 function qbpp(t9, npr, effmp, taup) 

  use constants, only: RHO_TO_N, RHO_NUC
  implicit none
  real*8, intent(in) :: t9, npr, effmp, taup
  real*8 rho0p, alphapp, betapp, nfnu
  real*8 xv, v1, RBpp, RRpp_pA

  rho0p = RHO_TO_N*npr/RHO_NUC
  alphapp = 0.11d0
  betapp = 0.7d0
  nfnu = 3.d0  
  RBpp = 1.d0
  ! qBpp = 1.709d19*EFFMP**4*(RHO0**(1d0/3d0))*T9**8  !(FM79) 
  qBpp = 7.4d19*effmp**4*(rho0p**(1d0/3d0))*alphapp*betapp*3*nfnu*T9**8  ! (YL95)

  ! Superfluidity correction
  if (taup <= 1.d0) then   ! SF protons
    xv = v1(taup)  ! 1S_0 (pA)
    RBpp = RRpp_pA(xv)
  endif

  qBpp = qBpp*RBpp

end function qbpp


! Modified URCA emissivity
! Friman & Maxwell (1979) ApJ 232, 541     (FM79)
! Yakovlev & Levenfish (1995) A&A 297, 717 (YL95)
! Yakovlev et al. (2001) Physics Reports, Volume 354, Issue 1-2, p. 1-155.  (Y01)
real*8 function qmurca(t9, nn, npr, ye, yn, yp, ymu, effmn, effmp, taun, taup)

  use constants, only: RHO_TO_N, RHO_NUC
  implicit none
  real*8, intent(in) :: t9, nn, npr, ye, yn, yp, ymu, effmn, effmp, taun, taup
  real*8 rho0n, rho0p, alphan, betan, taumax, taumin, xvn, xvp, v1, v2
  real*8 qMn_mu, qMp_mu
  real*8 qMn, qMp, RMn, RMp, RRMn_nB, RRMp_nB, RRMn_pA, RRMp_pA, RRMn_BA, RRMp_BA

  rho0n = RHO_TO_N*nn/RHO_NUC
  rho0p = RHO_TO_N*npr/RHO_NUC  ! TBD UNUSED?
  RMn = 1.d0
  RMp = 1.d0

! rho0 = RHO_TO_N*nb/RHO_NUC
! qmurca = 1.8d21*effmn**3*effmp*(rho0**(2./3.))*T9**8     !FM79
! return       

  alphan = 1.76d0 - 0.63d0/(rho0n+1.d-20)**(2.d0/3.d0) ! OPE FM79
  betan = 0.68d0             ! for a given EoS
  ! neutron branch
  if ( (alphan > 0d0) ) then   
    qMn = 8.55d21*effmn**3*effmp*rho0p**(1d0/3d0)*T9**8*alphan*betan 
  else
    qMn = 0.d0
  endif

  ! proton branch
  if ( (alphan > 0d0) .and. &
 &     (ye**(1d0/3d0) + 3.d0*yp**(1d0/3d0) > yn**(1d0/3d0)) .and. ymu>0d0) then 
    qMp = 8.53d21*effmp**3*effmn*rho0p**(1d0/3d0)*T9**8*alphan*betan &
 &        *max((1.d0-0.25d0*(ye/yp)**(1d0/3d0)),0.d0)
  else
    qMp = 0.d0
  endif

  ! muon branches ! Y01
  if ( (alphan > 0d0) .and. &
 &     (ymu**(1d0/3d0) + 3.d0*yp**(1d0/3d0) > yn**(1d0/3d0)) .and. ymu>0d0) then 
    qMn_mu = 8.1d21*effmn**3*effmp*rho0p**(1d0/3d0)*T9**8*alphan*betan*(ymu/ye)**(1d0/3d0)
    qMp_mu = qMn_mu*((effmp/effmn)**(1d0/3d0))*(ymu/ye)**(1d0/3d0) &
 &   *(ymu**(1d0/3d0)+3.d0*yp**(1d0/3d0)-yn**(1d0/3d0)**2)/(8*ymu**(1d0/3d0)*yp**(1d0/3d0))
  else
    qMp_mu = 0.d0
    qMn_mu = 0.d0
  endif
  qMp = qMp + qMp_mu
  qMn = qMn + qMn_mu

  ! Superfluidity correction
  taumax = max(taun,taup)
  taumin = min(taun,taup)
  if (taumin > 1d0) go to 2400                  ! no SF
  if (taumax <= 1d0) then !2-comp SF
    xvn = v2(taun) !n in 3P_2
    xvp = v1(taup) !p in 1S_0
    RMn = RRMn_BA(xvn,xvp)
    RMp = RRMp_BA(xvn,xvp)
  elseif (taumax > 1d0) then                    ! only one component SF
    if (taun <= 1d0) then      ! only n SF
      xvn = v2(taun) ! 3P_2 (nB)    
      RMn = RRMn_nB(xvn)
      RMp = RRMp_nB(xvn)           
    elseif (taup <= 1.d0) then  !only p SF
      xvp = v1(taup)             !1S_0 (pA)
      RMn = RRMn_pA(xvp)
      RMp = RRMp_pA(xvp)
    endif
  endif

2400   continue
  qMn = qMn*RMn  !neutron branch
  qMp = qMp*RMp  !proton branch
  qMurca = qMn+qMp

end function qmurca

! Direct URCA emissivity
! Lattimer et al. (1991) PRL, 66, 2701
! DUrca e-channel
real*8 function qdurca(t9, rho, effmn, effmp, ye, yp, taun, taup)

  use constants, only: RHO_NUC
  implicit none
  real*8, intent(in) :: t9, rho, effmn, effmp, ye, yp, taun, taup
  real*8 rho0, qdu_e
  real*8 taumax, taumin, xvn, xvp, v1, v2
  real*8 RD, RRD_A, RRD_B, RRD_BA

  rho0 = rho/RHO_NUC

  RD = 1.d0
  qdu_e = 0d0
  qdu_e = 4.0d27*effmn*effmp*t9**6*(ye*rho0)**(1d0/3d0)

  qdurca = qdu_e

  ! Superfluid correction
  taumax = max(taun,taup)
  taumin = min(taun,taup)
  if (taumin > 1d0) go to 2500 ! no SF
  if (taumax <= 1d0) then ! both components SF
    xvn = v2(taun) !n in 3P_2
    xvp = v1(taup) !p in 1S_0
    RD = RRD_BA(xvn,xvp)
  elseif (taumax > 1d0) then ! only one component SF
    if (taun <= 1d0) then     !only n in 3P_2
      xvn = v2(taun)
      RD = RRD_B(xvn)
    elseif (taup <= 1d0) then !only p in 1S_0
      xvp = v1(taup)        
      RD = RRD_A(xvp)
    endif
  endif

2500  continue
  qdurca = qdurca*RD

end function qdurca



! e-N Bremsstrahlung emissivity
! Kaminker et al., A&A (1999), 343, 1009 
! Validity: in the crust, 5e7 < T[K] < 2e9, 1e9 < rho[g/cm3] < 1.4e14 [g/cm3]
real*8 function qeabrems(t, rho, xh)

  use constants, only: RHO_NUC
  implicit none
  real*8, intent(in) :: t, rho, xh
  real*8 rho0, xtau, xr
  real*8 a1, a2, a3, a4, a5, a6, a7, a8, a9

  rho0 = rho/RHO_NUC
  xtau = dlog10(t*1.d-8)
  xr = dlog10(rho*1.d-12)

  ! Range of validity of the approach
  if ( ( (t >= 5d7) .and. (t <= 2d9) ) .and. & 
 &   ( ( rho >= 1d9) .and. (rho <= 1.4d14)) ) then

  a1 = 11.204d0
  a2 = 7.304d0
  a3 = 0.2976d0
  a4 = 0.370d0
  a5 = 0.188d0
  a6 = 0.103d0
  a7 = 0.0547d0
  a8 = 6.77d0
  a9 = 0.228d0
  qeabrems = 10.d0**( a1 + a2*xtau + a3*xr - a4*xtau**2 + a5*xtau*xr &
 &     - a6*xr**2 + a7*xtau**2*xr - a8*dlog(1+a9*rho0) )
  else
     qeabrems = 0.d0
  endif
  
end function qeabrems


! Plasma emissivity
! Yakovlev  et al. Phys.Rep. 354 (2001) 1-155
real*8 function qplasma(t, kFe)

  use constants, only: PI, K_BOLTZMANN, MASS_E_MEV, ALPHA
  implicit none
  real*8, intent(in) :: t, kFe
  real*8 me, Qc, xtr, xr, fp, xip, sumCV2, zexp1

  Qc= 1.203d23
  me = MASS_E_MEV*1d6  ! electron mass in eV/c^2
  xtr = K_BOLTZMANN*t/me   
  xr = kFe*197.326d6/me ! kFe in eV divided by the mass
  fp = dsqrt(4.d0*ALPHA*xr**3/3.d0/PI/(dsqrt(1.d0+xr**2)))/xtr
  zexp1 = -fp
  if ( zexp1 < -1.5d2 ) zexp1 = -1.5d2 ! floor
  xip = xtr**9*(16.23d0*fp**6+4.604d0*fp**(7.5d0))*dexp(zexp1)
  sumcv2 = 0.9248d0
  qplasma = Qc*xip*sumCV2/96.d0/PI**4/ALPHA

end function qplasma

! Synchotron emissivity
! Bezchastnov et al., AAP (1997), 328, 409
real*8 function qsynchro(t, kFe, ye, bf) 

  use constants, only: PI, MASS_E_MEV
  implicit none
  real*8, intent(in) :: t, kFe, ye, bf
  real*8 xr, t9, b13, tb, z, xi
  real*8 cminus2, cplus2, dd1, dd2, alpha1, alpha2, y1, y2
  real*8 ff1, ff2, ff3, ff4, fplus, fminus, sab, ddf1, ddf2, sbc

  xr = kFe*197.326d6/(MASS_E_MEV*1d6) ! kFe in eV divided by the mass in eV
  t9 = t*1d-9
  b13 = bf*1.d-1
  tb = 1.34d9*b13/dsqrt(1.d0+xr**2)
  z = tb/t
  xi = 3.d0/2.d0*z*xr**3

  cminus2 = 0.175d0
  cplus2 = 1.675d0
  dd1 = 44.01d0
  dd2 = 36.97d0
  alpha1 = 3172.d0
  alpha2 = 172.2d0
  y1 = ( (1.d0 + alpha1*xi**(2.d0/3.d0) )**(2.d0/3.d0) - 1.d0 )**(3.d0/2.d0)
  y2 = ( (1.d0 + alpha2*xi**(2.d0/3.d0) )**(2.d0/3.d0) - 1.d0 )**(3.d0/2.d0)
  ff1 = 1.d0 + 3.675d-4*y1
  ff2 = 1.d0 + 2.036d-4*y1 + 7.405d-8*y1**2
  ff3 = 1.d0 + 1.436d-2*y2 + 1.024d-5*y2**2 + 7.647d-8*y2**3
  ff4 = 1.d0 + 3.356d-3*y2 + 1.536d-5*y2**2
  fplus = dd1*ff1**2/ff2**4
  fminus = dd2*ff3/ff4**5
  sab = 27.d0*xi**4*(fplus - cminus2*fminus/cplus2)/PI**2/512.d0/1.0369d0

  ddf1 = 1.d0 + 0.4228d0*z + 0.1014d0*z**2 + 0.006240d0*z**3
  ddf2 = 1.d0 + 0.4535d0*z**(2.d0/3.d0) + 0.03008d0*z - 0.05043d0*z**2 + 0.004314d0*z**3
  sbc = dexp(-z/2.d0)*DDf1/DDf2

  qsynchro = 9.04d14*sab*sbc*b13**2*t9**5

end function qsynchro


! Bremsstrahlung electron-proton emissivity
! Maxwell ApJ (1979) 231, 201
real*8 function qepbrems(t9, npr, effmp)

  use constants, only: RHO_TO_N, RHO_NUC
  implicit none
  real*8, intent(in) :: t9, npr, effmp
  real*8 rho0p

  rho0p = RHO_TO_N*npr/RHO_NUC
  qepbrems = 2.4d17*rho0p**(-2./3.)*t9**8*effmp**2

end function qepbrems


! Pair neutrino emissivity
real*8 function qpair(t, rho, z, a)

  implicit none
  real*8, intent(in) :: t, rho, z, a
  real*8 cv,cvp,x,x2,xi,chi,za,rza,COEF1,COEF2
  real*8 b1,b2,b3,c,f,g,q
  parameter(cv=0.9638, cvp=0.0362)
  parameter(coef1=0.840766, coef2=0.090766)

  if (t < 3.d8) then 
    qpair = 0.d0
    return
  endif

  x = t/5.9302d9
  x2 = x*x
  xi = 1.d0/x
  za = z/a
  rza = rho*za
  chi = 1.d-3*rza**(1./3.)*xi

  if (t < 1.d10) then
    b1 = 9.383d-1
    b2 = -4.141d-1
    b3 = 5.829d-2
    c = 5.5924
  else
    b1 = 1.2383  
    b2 = -0.8141  
    b3 = 0.0     
    c = 4.9924
  endif

  f = (6.002d19 + chi*(2.084d20 + chi*1.872d21))*dexp(-c*chi)/  &
 &      (chi**3 + xi*(b1 + xi*(b2 + xi*b3)) )
  g = 1.0 - x2*(13.04 - x2*(133.5 + x2*(1534. + x2*918.6)))
  q = (1.0050 + 0.3967*dsqrt(x) + 10.7480*x2)**(-1.0)*  &
 &    (1.0+rza/(7.692d7*x**3 + 9.715d6*dsqrt(x)))**(-0.3)

  qpair = (coef1 + coef2*q)*g*dexp(-2.d0*xi)*f

end function qpair


! Cooper pair formation and breaking emissivity in the crust (1S_0)
! Yakovlev et al. Phys.Report 354 (2001), 
! originally in Yakovlev, Kaminker & Levenfish, A&A 343 (1999). 
! Suppresion in the vector channel revised by Steiner, Reddy, PRC 79 (2009)
! implemented in Yakovlev's formula by Page et al. ApJ 707 (2009)
real*8 function qcpbf_cr(t9, nn, effmn, taucru) 

  use constants, only: RHO_TO_N, RHO_NUC
  implicit none
  real*8, intent(in) :: t9, nn, effmn, taucru
  real*8 ffa, v1, xv, Nv, rho0n, pFn_m, cvec_n, cax_n, SSn, aan

  rho0n = RHO_TO_N*nn/RHO_NUC
  pFn_m = 0.353d0*rho0n**(1d0/3d0) ! pF/mc
  Nv = 3.d0  !number of neutrino flavors
  ! Relativistic corrections (aan = 1 for non-relativistic case)
  cvec_n = 1.d0
  cax_n = 1.26d0
  SSn =  4.d0/8.1d1*(pFn_m/effmn)**4  ! suppression due to vector current conservation
  aan = cvec_n**2*SSn + cax_n**2*pFn_m**2*(1.d0 + 1.1d1/4.2d1/effmn**2)
  xv = v1(taucru)   
  qcpbf_cr = 1.17d21*effmn*pFn_m*Nv*aan*ffa(xv)*T9**7 

end function qcpbf_cr

! Cooper pair formation and breaking emissivity in the core (SF neutrons)
real*8 function qcpbf_con(t9, nn, effmn, taun) 

  use constants, only: RHO_TO_N, RHO_NUC
  implicit none
  real*8, intent(in) :: t9, nn, effmn, taun
  real*8 ffb, v2, xvn, nv
  real*8 rho0n, pFn_m, cvec_n, cax_n, SSn, aan

  rho0n = nn*RHO_TO_N/RHO_NUC
  nv = 3.d0
  pFn_m = 0.353d0*rho0n**(1d0/3d0)  ! pF/mc
  ! Relativistic corrections (aan = 1 for non-relativistic)
  cvec_n = 1.d0
  cax_n = 1.26d0
  ! suppression due to vector current conservation (ssn=1 for no suppression)
  ssn =  0.d0
  aan = cvec_n**2*SSn + 2.d0*cax_n**2 
  xvn = v2(taun)   
  qcpbf_con = 1.17d21*effmn*pFn_m*Nv*aan*ffb(xvn)*t9**7

end function qcpbf_con


! Cooper pair formation and breaking emissivity in the core (SF protons)
real*8 function qcpbf_cop(t9, npr, effmp, taup) 

  use constants, only: RHO_TO_N, RHO_NUC
  implicit none
  real*8, intent(in) :: t9, npr, effmp, taup
  real*8 ffa, v1, xvp, nv
  real*8 rho0p, pFp_m, cvec_p, cax_p, SSp, aap, sin2theta
  
  rho0p = npr*RHO_TO_N/RHO_NUC
  pFp_m = 0.353d0*rho0p**(1d0/3d0)
  ! Relativistic corrections (aap = 1 for non-relativistic)
  sin2theta = 0.23d0
  cvec_p = 4.d0*sin2theta - 1.d0
  cax_p = - 1.26d0
  ! suppression due to vector current conservation (ssp=1 for no suppression)
  ssp =  4.d0/8.1d1*(pFp_m/effmp)**4 !suppression due to vector current conservation
  aap = cvec_p**2*SSp + cax_p**2*pFp_m**2*(1.d0 + 1.1d1/4.2d1/effmp**2) 
  xvp = v1(taup)    
  qcpbf_cop = 1.17d21*effmp*pFp_m*Nv*aap*ffa(xvp)*T9**7

end function qcpbf_cop


real*8 function ffa(v) ! CPBF for 1S_0

  implicit none
  real*8, intent(in) :: v
  real*8 AA, BB, zexp1

  AA = 0.602d0*v**2+0.5942d0*v**4+0.288d0*v**6
  BB = 0.5547d0+dsqrt(0.4453d0**2+0.0113d0*v**2)
  zexp1 = 2.245d0-dsqrt(2.245d0**2+4.d0*v**2)
  if (zexp1 < -1.5d2) zexp1=-1.5d2 !cutoff
  ffa = AA*dsqrt(BB)*dexp(zexp1)  

end function ffa


real*8 function ffb(v) ! CPBF for 3P_2

  implicit none
  real*8, intent(in) :: v
  real*8 AA, BB, CC, zexp1

  AA = 1.204d0*v**2+3.733d0*v**4+0.3191d0*v**6
  BB = 0.7591d0+dsqrt(0.2409d0**2+0.3145d0*v**2)
  CC = 1.d0+0.3511d0*v**2
  zexp1 = 0.4616d0-dsqrt(0.4616d0**2+4.d0*v**2)
  if (zexp1 < -1.5d2) zexp1=-1.5d2 !cutoff
  ffb = AA*BB**2/CC*dexp(zexp1)  
  !if (zexp1 < -1.5d2) ffb = 0.d0

end function ffb

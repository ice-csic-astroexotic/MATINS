! ----------------------------------------------------------------
!     INTERFACE TO CALL ALL NEUTRINO EMISSIVITIES ROUTINES
! ----------------------------------------------------------------
      SUBROUTINE compute_neutrino_emissivity()

      use grid, only: rho,xh,ye,yn,yp,aa,zz,tccru,tcn,tcp
      use grid, only: kmax, lmax
      use grid, only: tem0, bmed, q_neutrino, rmc!, emnu
      use constants, only: UNIT_R

      implicit none
      integer k, l, j
      real*8 tem,tem1,emis,emis1,taucru,taun,taup,q3

      do l=2,lmax
        j = 2*l-1
      do k=1,kmax

        tem = tem0(k,l)     ! Physical temperature in 1e8 K

        taucru=tem/dmax1(tccru(j),0.1d0*tem)
        taun=tem/dmax1(tcn(j),0.1d0*tem)
        taup=tem/dmax1(tcp(j),0.1d0*tem)

        emis=q3(tem,rho(j),ye(j),aa(j),zz(j),xh(j),
     &      yn(j),yp(j),taucru,taun,taup,bmed(k,l))
        tem1=tem*1.01d0  ! Slightly larger T to evaluate the derivative
        emis1=q3(tem1,rho(j),ye(j),aa(j),zz(j),xh(j),
     &      yn(j),yp(j),taucru,taun,taup,bmed(k,l))

c ----------------------------------------------------------------
c	Elements of matrix
c    	q_neutrino = emis - (d emis/d(T*e^nu))*Te^nu
c	    rmc = d emis/dT (without redshift correction)
c     q_neutrino is given in 10**40 erg/km^3/s
c ----------------------------------------------------------------

        q_neutrino(k,l)= emis*UNIT_R**3*1.d-40  ! 10**40 erg/km^3/s
        rmc(k,l)= (emis1-emis)/(0.01d0*tem)*UNIT_R**3*1.d-40        ! 10**40 erg/(km^3*s*10**8 K)
      enddo
      enddo

      return
      end

c -----------------------------------------------------------------------------
c  This function calls other functions to obtain the neutrino emissivity.
c  Reduction factors due to superfluidity are calculated throughout.
c -----------------------------------------------------------------------------
      REAL*8 FUNCTION Q3(T8,rho,YE,AH,ZH,XH,YN,YP,taucru,taun,taup,bf)
      IMPLICIT NONE
c     -------------------------------------------------------------------------    
      real*8 t,t8,nbaryon,t9,rho,effmn,effmp
      real*8 ah,zh,xh,ye,yn,yp,taucru,taun,taup,bf
      real*8 fnn,fnp,fye
      real*8   qea,qpl,qsyn,qcp_cr,qpa
     -        ,qeABrems,qplasma,qsynchro,qcpbf_cr,qpair
     -        ,qmur,qnn,qnp,qpp,qep,qcp_co,qdu
     -        ,qmurca,qBnn,qBnp,qBpp,qepBrems,qcpbf_co,qdurca
     -        ,qtot
c
      t = 1.d8*t8            ! t: temperature in K
      t9 = 1.d-1*t8          ! t9: temperature (10**9 K),
      nbaryon=rho/1.66d-24  !nbaryon: baryon number density (cm-3)
c     ----define particle number----------
      fye=ye
      fnn=nbaryon*yn
      fnp=nbaryon*yp
      call eff_mass (fnn,fnp,effmn,effmp)

c     -----Initialization----!all q in [erg/cm3/sec]
      qmur= 0.d0
      qnn = 0.d0
      qnp = 0.d0
      qpp = 0.d0
      qep = 0.d0
      qcp_co = 0.d0
      qdu = 0.d0
      qeA = 0.d0
      qpl = 0.d0 
      qsyn = 0.d0  !here B should be replaced 
      qcp_cr = 0.d0
      qpa = 0.d0
c     -------------------Core processes------------------ 
c            Modified Urca
      qmur = qmurca(t9,rho,effmn,effmp,ye,yn,yp,taun,taup)   
c            NN Brems (nn->nn+vv, pp->pp+vv, np->np+vv)
      qnn = qBnn(t9,rho,effmn,yn,xh,taucru,taun)!also in the crust!!!     
      qnp = qBnp(t9,rho,effmn,effmp,yp,taun,taup)
      qpp = qBpp(t9,rho,effmp,yp,taup)   
c            ep Brems  (ep->ep+vv)
      qep = qepBrems(t9,rho,effmp,yp)
c            Cooper Pairs of SF 3P_2 neutrons & SF 1S_0 protons
      qcp_co= qcpbf_co(t,rho,effmn,effmp,xh,yn,yp,taun,taup)
c            Direct Urca  
c      if(idurca.eq.1) 
      qdu = qdurca(t9,rho,effmn,effmp,fye,yp,taun,taup)    
c     ------------------Crust processes-----------------
c            eA Bremsstrahlung (eA->eA+vv)
      qeA  = qeABrems(t,rho,xh) 
      qpl  = qplasma(t,rho,ye,xh) 
      qsyn = qsynchro(t,rho,ye,bf) !bfield in 1e12 G
      qcp_cr = qcpbf_cr(t,rho,effmn,xh,yn,taucru) 
      qpa    = qpair(t,rho,zh,ah)       

cc      qcp_cr=0.d0
cc      qcp_co=0.d0
c     ----------------------------------------------------------------  
      qtot = qmur+qnn+qnp+qpp+qep+qcp_co+qdu   !Core proc.
     -      +qeA
     -      +qpl
     -      +qsyn
     -     + qcp_cr
     -      +qpa           !Crust proc. 
c      qtot = qmur+qnn+qpp+qnp+qbr+qdu+qpl
      q3 = qtot             

      return
      end

c===================================================================
c     -------------Begin N-N Bremstrahlung----------- 
c===================================================================
      REAL*8 FUNCTION QBNN(T9,RHO,EFFMN,YN,xh,taucru,taun)
c===================================================================
C-----  CALCULATES nn BREMSSTRAHLUNG NEUTRINO EMISSIVITY -----------
C-----  FRIMAN AND MAXWELL (1979) Ap. J. 232, 541  ------------(FM79)
C-----  YAKOVLEV AND LEVENFISH (1995) A&A 297, 717  -----------(YL95)
C-------------------------------------------------------------------
      implicit none
      real*8 t9,rho,effmn,yn,xh,taucru,taun
      real*8 rho0n,e3,alphann,betann,nfnu,fmu,fn,kFn,xmpi,xu
      real*8 xv,v1,v2,RBnn,RRnn_nA,RRnn_nB
c
      rho0n=rho*yn/2.82d14
      e3=1.d0/3.d0
      alphann=0.59d0
      fmu=1.d0
      if ((xh.gt.0.d0).and.(xh.lt.1.d0)) then !in the inner crust----- 
      fmu=yn
      call kFermi(rho,yn,fn,kFn)  !kFn in 1/fm
      kFn=kFn*197.326d0 !kFn in MeV
      xmpi=139.d0  !pion mass in MeV
      xu=xmpi/2.d0/kFn
c      alphann=0.5d0
      alphann=1.d0-3.d0*xu*datan(1.d0/xu)/2.d0+xu**2/2.d0/(1.d0+xu**2)
      endif           !-------------------------------------
      if (xh.eq.1.d0) alphann=0.d0   !in the outer crust----- 
      betann=0.56d0
      nfnu= 3.d0  ! number of neutrino flavor generated. 
      RBnn=1.d0
c     qBnn = 7.335d19*EFFMN**4*(RHO0**e3)*T9**8 (FM79) 
      qBnn = 7.4d19*effmn**4*rho0n**e3*alphann*betann*3*T9**8*nfnu*fmu   !(YL95) 
c      return !to avoid SF
c     ------Superfluidity-----------------------------------
       if ((taun.le.1.d0).and.(xh.eq.0.d0)) then  !only neutrons in the core
cc         xv=v1(taun)   !1S_0 (nA)  !only for tests
cc         RBnn=RRnn_nA(xv)
         xv=v2(taun)   !3P_2 (nB)  
         RBnn=RRnn_nB(xv)
       endif
       if ((taucru.le.1.d0).and.(xh.gt.0.d0)) then   !only neutrons in the crust
         xv=v1(taucru)    !1S_0 (nA)
         RBnn=RRnn_nA(xv) 
       endif
c     -------------------------------------------------------
      qBnn=qBnn*RBnn
      return
      END
c===================================================================
       REAL*8 FUNCTION QBNP(T9,RHO,EFFMN,EFFMP,YP,taun,taup)
c===================================================================
C-----  CALCULATES np BREMSSTRAHLUNG NEUTRINO EMISSIVITY -----------
C-----  FRIMAN AND MAXWELL (1979) Ap. J. 232, 541  ---------- (FM79)
C-----  YAKOVLEV AND LEVENFISH (1995) A&A 297, 717  --------- (YL95)
C-------------------------------------------------------------------
      implicit none
      real*8 t9,rho,effmn,effmp,yp,taun,taup
      real*8 rho0p,e3,alphanp,betanp,nfnu
     -      ,taumax,taumin,xvn,xvp,v1,v2
     -      ,RBnp,RRnp_nB,RRnp_pA,RRnp_BA

      rho0p=rho*yp/2.82d14
      alphanp=1.06d0
      betanp=0.66d0
      nfnu=3.d0
      e3=1.d0/3.d0
      RBnp=1.d0
c
c     qBnp = 3.15D20*EFFMN**2*EFFMP**2*(RHO0**(2./3.))*T9**8 !(FM79)
      qBnp = 1.5d20*(effmn*effmp)**2*rho0p**e3*
     &          alphanp*betanp*nfnu*T9**8           !(YL95)
c     ------------Superfluidity---------------------------------
       taumax= max(taun,taup)
       taumin= min(taun,taup)
       if (taumin.gt.1.d0) go to 2200                  !0-comp SF
c      then, here taumin.le.1.d0, at least 1-comp SF
       if (taumax.le.1.d0) then !2-comp SF
          xvn=v2(taun) !n in 3P_2 
          xvp=v1(taup) !p in 1S_0
          RBnp=RRnp_BA(xvn,xvp) 
       elseif (taumax.gt.1.d0) then                    !1-comp SF
          if (taun.le.1.d0) then     !only n SF
            xvn=v2(taun) !3P_2 (nB)  
            RBnp=RRnp_nB(xvn)
          elseif (taup.le.1.d0) then !only p in 1S_0
            xvp=v1(taup)  !1S_0 (pA)
            RBnp=RRnp_pA(xvp)
          endif
       endif
c     -----------------------------------------------------
2200  continue
      qBnp=qBnp*RBnp
      return
      end

c===================================================================
      REAL*8 FUNCTION QBPP(T9,RHO,EFFMP,YP,taup) 
c===================================================================
C-----  CALCULATES pp BREMSSTRAHLUNG NEUTRINO EMISSIVITY -----------
C-----  FRIMAN AND MAXWELL (1979) Ap. J. 232, 541  ---------- (FM79)
C-----  YAKOVLEV AND LEVENFISH (1995) A&A 297, 717 ---------- (YL95)
c===================================================================

      implicit none
      real*8 t9,rho,effmp,yp,taup
      real*8 rho0p,e3,alphapp,betapp,nfnu
      real*8 xv,v1,RBpp,RRpp_pA
C
      rho0p=rho*yp/2.82d14
      e3=1.d0/3.d0
      alphapp=0.11d0
      betapp=0.7d0
      nfnu= 3.d0  
      RBpp=1.d0
c     qBpp = 1.709d19*EFFMP**4*(RHO0**e3)*T9**8  !(FM79) 
      qBpp = 7.4d19*effmp**4*(rho0p**e3)*alphapp*betapp*3*nfnu*T9**8  ! OK with (YL95)
c     ------------Superfluidity---------------------------------
       if (taup.le.1.d0) then   !only protons
         xv=v1(taup)  !1S_0 (pA)
         RBpp=RRpp_pA(xv)
       endif
C     ----------------------------------------------------------
      qBpp=qBpp*RBpp
      return
      end

c====================================================================
c      -------------end N-N Bremstrahlung----------- 
c====================================================================

c====================================================================
c      -------------begin Urca processes----------- 
c====================================================================
       REAL*8 FUNCTION QMURCA(T9,RHO,EFFMN,EFFMP,YE,YN,YP,taun,taup)
c====================================================================
C----- CALCULATES MODIFIED URCA NEUTRINO EMISSIVITY --------------------
C----- FRIMAN AND MAXWELL (1979) Ap. J. 232, 541  --------------------- (FM79)
C----- Yakovlev & Levenfish (1995) A&A, 297, 717  --------------------- (YL95)
c-----------------------------------------------------------------------
       implicit none
       real*8 t9,rho,effmn,effmp,ye,yn,yp,taun,taup
       real*8 rho0,rho0n,rho0e,rho0p,e3,alphan,betan
     -       ,taumax,taumin,xvn,xvp,v1,v2
       real*8 qMn,qMp,RMn,RMp
     -       ,RRMn_nB,RRMp_nB,RRMn_pA,RRMp_pA
     -       ,RRMn_BA,RRMp_BA
c
       e3=1.d0/3.d0
       rho0 = rho/2.82d14
       rho0n=rho0*yn
       rho0e=rho0*ye
       rho0p=rho0*yp
       RMn=1.d0
       RMp=1.d0
c
c      QMURCA = 1.8d21*effmn**3*effmp*(rho0**(2./3.))*T9**8     !FM79
c      return       
       alphan=1.76d0-0.63d0/(rho0n+1.d-20)**(2.d0/3.d0) !OPE FM79
       betan=0.68d0             ! for a given EoS
c      ----------------------------neutron branch------
       if ((alphan.gt.0.d0).and.(yp.gt.0.d0)) then   
        qMn=8.55d21*effmn**3*effmp*rho0e**e3*T9**8*alphan*betan 
       else
        qMn=0.d0
       endif
c      ----------------------------proton branch-------
       if ((alphan.gt.0.d0).and.(yp.gt.0.d0).and.
     &      (ye**e3+3.d0*yp**e3.gt.yn**e3)) then 
        qMp=8.53d21*effmp**3*effmn*rho0e**e3*T9**8*alphan*betan
     &         *max((1.d0-0.25d0*(ye/yp)**e3),0.d0)
       else
        qMp=0.d0
       endif
c      -------------Superfluidity------------------------
       taumax= max(taun,taup)
       taumin= min(taun,taup)
       if (taumin.gt.1.d0) go to 2400                  !0-comp SF
c      then, here taumin.le.1.d0, at least 1-comp SF
       if (taumax.le.1.d0) then !2-comp SF
          xvn=v2(taun) !n in 3P_2
          xvp=v1(taup) !p in 1S_0
          RMn=RRMn_BA(xvn,xvp)
          RMp=RRMp_BA(xvn,xvp)
       elseif (taumax.gt.1.d0) then                    !1-comp SF
          if (taun.le.1.d0) then      !only n SF
            xvn=v2(taun) !3P_2 (nB)    
            RMn=RRMn_nB(xvn)
            RMp=RRMp_nB(xvn)           
          elseif (taup.le.1.d0) then  !only p SF
            xvp=v1(taup)             !1S_0 (pA)
            RMn=RRMn_pA(xvp)
            RMp=RRMp_pA(xvp)
          endif
       endif
c      ------------------------------------------
2400   continue
       qMn=qMn*RMn  !neutron branch
       qMp=qMp*RMp  !proton branch
       qMurca=qMn+qMp
       return
       end
c====================================================================
       REAL*8 FUNCTION QDURCA(T9,RHO,effmn,effmp,ye,yp,taun,taup)
c====================================================================
C-------------------------------------------------------------------
C----- CALCULATES DIRECT URCA NEUTRINO EMISSIVITY ------------------
C----- Lattimer et al. (1991) Phys. Rev. Lett., 66, 2701  ----------
C-------------------------------------------------------------------
       implicit none
       real*8 t9,rho,ye,yp,effmn,effmp,taup,taun
       real*8 e3,rho0,nbaryon,qdu_e,qdu_mu
     -   ,taumax,taumin,xvn,xvp,v1,v2
     -   ,RD,RRD_A,RRD_B,RRD_BA!,rho_ce,rho_cmu
C       
       e3=1.d0/3.d0            
       nbaryon = rho/1.66d-24   ! u=1.66D-24, atomic mass unit in gr
       rho0 = rho/2.8d14        ! rho_0=2.8D14 nuclear density
c      ---------DU threshold for a given EoS------
!         rho_ce=10.18d14 !for lowd-eos3.tab 
c        rho_ce=26.0d14 !for lowd-eos.tab
c       rho_ce=7.8d14 !for lowd-eos2.tab
c       rho_ce=0.5*1.298D15 !Moderate EoS, see Yakovlev review, p. 61 
c      -------------------------------------------
c       rho_cmu=1.18D15    muons disabled
       RD=1.d0
!        if (rho.gt.rho_ce)  then !e onset 
!           qdu_e = 4.0d27*effmn*effmp*t9**6*(ye*rho0)**e3
!        else
       if (yp.ge.0.11)  then !e onset 
          qdu_e = 4.0d27*effmn*effmp*t9**6*(ye*rho0)**e3
       else
          qdu_e=0.d0
       endif             
c       if (rho.gt.rho_cmu) then  !muons onset (mu_e=mu_muons)
c          qdu_mu = 4.0d27*effmn*effmp*t9**6*(ye*rho0)**e3
c       else
          qdu_mu =0.d0
c       endif
c
       qdurca=qdu_e+qdu_mu
cc       qdurca=0.d0
c       go to 2500   !to avoid SF
c      ------------Superfluidity----------------------------
       taumax= max(taun,taup)
       taumin= min(taun,taup)
       if (taumin.gt.1.d0) go to 2500 !0-comp SF
c      then, here taumin.le.1.d0, at least 1-comp SF
       if (taumax.le.1.d0) then !2-comp SF
          xvn=v2(taun) !n in 3P_2
          xvp=v1(taup) !p in 1S_0
          RD=RRD_BA(xvn,xvp)
       elseif (taumax.gt.1.d0) then !1-comp SF
          if (taun.le.1.d0) then     !only n in 3P_2
           xvn=v2(taun)
           RD=RRD_B(xvn)
          elseif (taup.le.1.d0) then !only p in 1S_0
           xvp=v1(taup)        
           RD=RRD_A(xvp)
          endif
       endif
c      ----------------------------------------------------
2500   continue
       qdurca=qdurca*RD

       return
       end
c====================================================================     
c      -----end Urca processes----------- 
c====================================================================

!====================================================================
!      ----Begin e-N Bremstrahlung processes---- 
!====================================================================
      REAL*8 FUNCTION QeABREMS(t,rho,xh) !new version-called also QABREM
!====================================================================
!       Calculates \nu emissivity from e-N Bremsstrahlung  
!       from Kaminker,Pethick et al., Astron. & Astrophys.(1999)343,1009 
!       Validity: 5d7<t<2d9 [K], 1d9<rho<1.4d14 [g/cm3]
!
      implicit none
      real*8 t,rho,xh,t8,xtau,rho12,rho0,xr,a1,a2,a3,a4,a5,a6,a7,a8,a9

      rho0 = rho/2.8d14 !nuclear sat. dens.
        t8=t*1.d-8
        xtau=dlog10(t8)
        rho12=rho*1.d-12
        xr=dlog10(rho12)

        if (((t.ge.5.d7).and.(t.le.2.d9)).and. !validity of the approach
     -     ((rho.ge.1.d9).and.(rho.le.1.4d14)).and.(xh.gt.0.d0)) then !not in core

         a1=11.204d0
         a2=7.304d0
         a3=0.2976d0
         a4=0.370d0
         a5=0.188d0
         a6=0.103d0
         a7=0.0547d0
         a8=6.77d0
         a9=0.228d0
         QeABREMS=10.d0**(a1+a2*xtau+a3*xr-a4*xtau**2+a5*xtau*xr
     -               -a6*xr**2+a7*xtau**2*xr-a8*dlog(1+a9*rho0))
        else
         QeABREMS=0.d0
        endif
        return
        end
c====================================================================

c====================================================================     
c      ----end e-N Bremstrahlung proceses---- 
c====================================================================
c====================================================================
c
c====================================================================
c====================================================================
c      ----Begin other proceses---- 
c====================================================================
c====================================================================
c====================================================================
      REAL*8 FUNCTION QPLASMA(t,rho,ye,xh) !new version
c====================================================================
C-----------------------------------------------------------------------
C-----  CALCULATES PLASMA NEUTRINO EMISSIVITY --------------------------
C       Yakovlev  et al. Phys.Rep. 354 (2001) 1-155.
C-----------------------------------------------------------------------
        implicit none
        real*8 t,rho,ye,xh
        real*8 pi,kB,me,Qc,fne,kFe,xtr,xr,alpha,fp
     -         ,XIP,sumCV2,zexp1
c
        pi=3.141592653589793d0
        kB=8.617d-5  !(k Boltzman in eV/K)
        me= 0.510998d6 !(me in eV/c^2)
        Qc= 1.203d23
        xtr = kB*T/me   
        call kFermi(rho,ye,fne,kFe)  !kFe in fm-1
        kFe=kFe*197.326d6  !kFe in eV
        xr = kFe/me
        alpha=1.d0/137.d0
        fp = dsqrt(4.d0*alpha*xr**3/3.d0/pi/(dsqrt(1.d0+xr**2)))/xtr
c
        zexp1=-fp
         if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff         
        XIP=xtr**9*(16.23d0*fp**6+4.604d0*fp**(7.5d0))*dexp(zexp1)
        sumCV2=0.9248d0
        if (xh.eq.0.d0) Qc=0.d0  !no plasma if core
c
        QPLASMA=Qc*XIP*sumCV2/96.d0/pi**4/alpha
        RETURN
        END
c===================================================================
      REAL*8 FUNCTION QSYNCHRO(T,rho,ye,Bfield12) 
c====================================================================
C-----------------------------------------------------------------------
C-----  Calculates synchrotron neutrino emissivity
C\bibitem[Bezchastnov et al.(1997)]{1997A&A...328..409B} Bezchastnov, V.~G., Haensel, P., Kaminker, A.~D., \& Yakovlev, D.~G.\ 1997, \aap, 328, 409 
C-----------------------------------------------------------------------
        implicit none
        real*8 pi,T,rho,ye,Bfield12
     -     ,e23,e32,me,fne,kFe,xr,T9,B13,Tb,z,xi
     -     ,cminus2,cplus2,DD1,DD2,alpha1,alpha2,y1,y2
     -     ,ff1,ff2,ff3,ff4,Fplus,Fminus,Sab,DDf1,DDf2,Sbc
        e23=2.d0/3.d0
        e32=3.d0/2.d0
        pi=3.141592653589793d0
        me= 0.510998d6 !(me in eV/c^2)
        call kFermi(rho,ye,fne,kFe)  !kFe in fm-1
        kFe=kFe*197.326d6  !kFe in eV
        xr = kFe/me
c       
        T9=T*1.d-9
c        B13=B*1.d-13
        B13=Bfield12*1.d-1
        TB=1.34d9*B13/dsqrt(1.d0+xr**2)
        z=TB/T
        xi=3.d0/2.d0*z*xr**3
c
        cminus2=0.175d0
        cplus2=1.675d0
        DD1=44.01d0
        DD2=36.97d0
        alpha1=3172.d0
        alpha2=172.2d0
        y1=((1.d0+alpha1*xi**e23)**e23-1.d0)**e32
        y2=((1.d0+alpha2*xi**e23)**e23-1.d0)**e32
        ff1=1.d0+3.675d-4*y1
        ff2=1.d0+2.036d-4*y1+7.405d-8*y1**2
        ff3=1.d0+1.436d-2*y2+1.024d-5*y2**2+7.647d-8*y2**3
        ff4=1.d0+3.356d-3*y2+1.536d-5*y2**2
        Fplus=DD1*ff1**2/ff2**4
        Fminus=DD2*ff3/ff4**5

        Sab=27.d0*xi**4*(Fplus-cminus2*Fminus/cplus2)/
     -             pi**2/512.d0/1.0369d0
c
        DDf1=1.d0+0.4228d0*z+0.1014d0*z**2+0.006240d0*z**3
        DDf2=1.d0+0.4535d0*z**e23
     -          +0.03008d0*z-0.05043d0*z**2+0.004314d0*z**3
        Sbc=dexp(-z/2.d0)*DDf1/DDf2
c
        QSYNCHRO=9.04d14*Sab*Sbc*B13**2*T9**5
        RETURN
        END
c====================================================================
c====================================================================
c      ----end other proceses---- 
c====================================================================
c====================================================================

c=========================================================================
c      Processes not considered - yet D.N.A 02/07
c======================================================
      REAL*8 FUNCTION QepBREMs(T9,RHO,EFFMP,YP)
c======================================================
C-----------------------------------------------------------------------
C-----  Bremsstrahlung electron-proton Ref: Maxwell, 79, ApJ 231, 201 --
C-----------------------------------------------------------------------
C
      implicit none
      real*8 t9,rho,effmp
      real*8 yp,rho0p
C
      rho0p = rho*yp/2.82d14

      IF (YP.GT.0.d0) THEN
        QepBREMs = 2.4d17*RHO0p**(-2./3.)*T9**8*EFFMP**2
      ELSE
        QepBREMs = 0.d0
      ENDIF
C
      RETURN
      END
c====================================================================
      REAL*8 FUNCTION QPAIR(T,RHO,Z,A)
c====================================================================
C-----  CALCULATES PAIR NEUTRINO EMISSIVITY -----------------------------
C-----  ITOH ET AL. (Feb,1996) Ap. J. Suppl. ---------------------------

        implicit none
        real*8 t,rho,z,a
        real*8 cv,cvp,x,x2,xi,chi,za,rza,COEF1,COEF2
        real*8 b1,b2,b3,c,f,g,q
        parameter(cv=0.9638, cvp=0.0362)
        parameter(coef1=0.840766, coef2=0.090766)
C
        if (t.lt.3.d8) then 
          qpair = 0.d0
          return
        endif
C
        x = t/5.9302d9
        x2 = x*x
        xi = 1.d0/x
        za = z/a
        rza = rho*za
        chi = 1.d-3*rza**(1./3.)*xi
C
        if (t.lt.1.d10) then
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
C
        f = (6.002d19 + chi*(2.084d20 + chi*1.872d21))*dexp(-c*chi)/
     -      (chi**3 + xi*(b1 + xi*(b2 + xi*b3)) )
        g = 1.0 - x2*(13.04 - x2*(133.5 + x2*(1534. + x2*918.6)))
        q = (1.0050 + 0.3967*dsqrt(x) + 10.7480*x2)**(-1.0)*
     -      (1.0+rza/(7.692d7*x**3 + 9.715d6*dsqrt(x)))**(-0.3)
C
        QPAIR = (coef1 + coef2*q)*g*dexp(-2.d0*xi)*f
C
        RETURN
      END
c====================================================================

c===================================================================
      REAL*8 FUNCTION QCPBF_cr(T,rho,effmn,xh,yn,taucru) 
c     From Yakovlev et al. Phys.Report 354 (2001), 
c     originally in Yakovlev, Kaminker & Levenfish, A&A 343 (1999). 
c     Suppresion in the vector chanel revised by Steiner, Reddy, PRC 79 (2009)
c     Implemented in Yakovlevs formula by Page et al. ApJ 707 (2009)
c====================================================================
C     Calculates CPBF neutrino emissivity in neutron SF in the crust- 1S_0
C--------------------------------------------------------------------
      implicit none
      real*8 T,rho,effmn,xh,yn,taucru,fmu
     -,e3,T9,rho0,FFA,v1,xv,Nv
     -,rho0n,pFn_m,cvec_n,cax_n,SSn,aan

      e3=1.d0/3.d0
      T9=T*1.d-9
      rho0 = rho/2.82d14
      rho0n=rho0*yn
      pFn_m = 0.353d0*rho0n**e3 !pF/mc
      fmu=yn
      Nv = 3.d0  !number of neutrino flavors
c     ----------relativistic corrections--------------------
      cvec_n = 1.d0
      cax_n = 1.26d0
c      SSn  = 1.d0                       !no suppression
      SSn =  4.d0/8.1d1*(pFn_m/effmn)**4 !suppression due to vector current conservation
      aan = cvec_n**2*SSn +
     -      cax_n**2*pFn_m**2*(1.d0 + 1.1d1/4.2d1/effmn**2) 
c      aan = 1.D0                       !non-relativistic
c     ------------------------------------------------------
      QCPBF_cr=0.d0
      if ((taucru.le.1.d0).and.(xh.gt.0.d0).and.(xh.lt.1.d0)) then   !only SFinner crust
      xv=v1(taucru)   
      QCPBF_cr=1.17d21*effmn*pFn_m*Nv*aan*FFA(xv)*T9**7!*fmu  
      endif
      RETURN
      END
c===================================================================
      REAL*8 FUNCTION QCPBF_co(T,rho,effmn,effmp,xh,yn,yp,taun,taup) 
c====================================================================
C      Calculates CPBF neutrino emissivity in neutron & proton SF in the core
C-----------------------------------------------------------------------
      implicit none
      real*8 T,rho,effmn,effmp,xh,yn,yp,taun,taup
     -,e3,T9,rho0,gA,FFA,FFB,v1,v2,xvn,xvp,Nv
     -,rho0n,pFn_m,cvec_n,cax_n,SSn,aan,QCPBF_con
     -,rho0p,pFp_m,cvec_p,cax_p,SSp,aap,QCPBF_cop,sin2theta
c
      e3=1.d0/3.d0
      T9=T*1.d-9
      rho0 = rho/2.82d14
      gA = 1.26d0 
      Nv = 3.d0

c     ====================Neutrons in 3P_2===============================
      rho0n=rho0*yn
      pFn_m = 0.353d0*rho0n**e3  !pF/mc
c     ----------relativistic corrections--------------------
      cvec_n = 1.d0
      cax_n = gA
c      aan = 1.D0                   !non-relativistic
c      SSn  = 1.d0                  !no suppression
      SSn =  0.d0                   !suppression due to vector current conservation
      aan = cvec_n**2*SSn + 2.d0*cax_n**2 
c     ------------------------------------------------------
      QCPBF_con=0.d0
      if ((taun.le.1.d0).and.(xh.eq.0.d0)) then !only in the SF core
       xvn=v2(taun)   
       QCPBF_con=1.17d21*effmn*pFn_m*Nv*aan*FFB(xvn)*T9**7
      endif


c     =========================Protons in 1S_0=========================
      rho0p=rho0*yp
      pFp_m = 0.353d0*rho0p**e3
c     ----------relativistic corrections--------------------
      sin2theta = 0.23d0
      cvec_p = 4.d0*sin2theta-1.d0
      cax_p = -gA
c      aap = 1.D0                        !non-relativistic
c      SSp  = 1.d0                       !no suppression
      SSp =  4.d0/8.1d1*(pFp_m/effmp)**4 !suppression due to vector current conservation
      aap = cvec_p**2*SSp +
     -      cax_p**2*pFp_m**2*(1.d0 + 1.1d1/4.2d1/effmp**2) 
c     ------------------------------------------------------
      QCPBF_cop=0.d0
      if ((taup.le.1.d0).and.(xh.eq.0.d0))then !only in the SF core
       xvp=v1(taup)    
       QCPBF_cop=1.17d21*effmp*pFp_m*Nv*aap*FFA(xvp)*T9**7
      endif

      QCPBF_co=QCPBF_con+QCPBF_cop
      RETURN
      END

c---------------------------------------------------------------
        Real*8 function FFA(v) !CPBF for 1S_0
      implicit none
        real*8 v,AA,BB, zexp1
c
        AA=0.602d0*v**2+0.5942d0*v**4+0.288d0*v**6
        BB=0.5547d0+dsqrt(0.4453d0**2+0.0113d0*v**2)
        zexp1=2.245d0-dsqrt(2.245d0**2+4.d0*v**2)
        if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff
        FFA=AA*dsqrt(BB)*dexp(zexp1)  
        return
        end
c------------------------------------------------------------------
        Real*8 function FFB(v) ! CPBF for 3P_2
      implicit none
        real*8 v,AA,BB,CC, zexp1
c
        AA=1.204d0*v**2+3.733d0*v**4+0.3191d0*v**6
        BB=0.7591d0+dsqrt(0.2409d0**2+0.3145d0*v**2)
        CC=1.d0+0.3511d0*v**2

        zexp1=0.4616d0-dsqrt(0.4616d0**2+4.d0*v**2)
        if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff
        FFB=AA*BB**2/CC*dexp(zexp1)  
        return
        end

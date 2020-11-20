! PHONONS CONDUCTIVITY (call it to include SF phonons by means of ksfph)

       subroutine cphonon(t,zh,ah,Zimp,ye,rho,xh,yn,taucru,B12,kph_tot)
       implicit none
       real*8 t,rho,ye,ah,zh,xh,yn,Zimp,taucru
       real*8 kii_I,kphph_II,kphimp_II,kph_II,klph
       real*8 ksfph_1,ksfph_2,ksfph_3,ksfph,kph_tot,k_tot_Chu,k_cru
c       real*8 F_kii_I,F_kphph_II
       real*8 kG1,kG2,gammaL,L_ie_Chu,L_ie_cru,L_lph,omtau_th,
     -        mfp_sfph,mfp_sfph_1,mfp_sfph_2,
     -        mfp_sfph_3, mfp_sfph_4,gmix,om,xalpha,pF_M,
     -        Rgmix
       real*8 kb,em,me,mu,clight,hbar,hc,mn,pi,e13
       real*8 rho6,t9,t7,xe,nbaryon
       real*8 fdebye,tdebye,gamma,nions,ni,cvion,ccvfion,alat
       real*8 aprime,a2prime,Tp,q_BZ,wp,cs_Chu,cs,cs_l,cs_t,
     -        cs_Itoh_t,cs_Itoh_l,ktf,k_F,nel
       real*8 mfp,mfpi,gimp,fu,Tu
       real*8 xx,r0,fnn,kFn,effmn,effmp,meffn,vs,ffalpha
       real*8 ccs,llambda
       real*8 dd,eta,aalpha
       real*8 AA,zzeta,a_ni
       real*8 k_tot
       real*8 B12, TEMP, UN_B12, UN_T6, GAMAG, RHONUC, cve
       parameter(UN_B12=425.438, UN_T6=.3157746)
      
c      -------------------------Constants and def.-----------------
       kb=1.381d-16     !erg/K    
       em=1.d0/137.d0
       me=9.11d-28      !gr
       mu=1.66d-24      !gr
       clight=2.9979d10 !cm/s
       hbar=1.0546d-27  !erg*s
       hc=197.326d0     !MeV fm
       mn=939.d0/hc     !neutron mass in 1/fm
       pi=3.141592653589793d0
       e13=1.d0/3.d0
       t9=t*1.d-9
       t7=t*1.d-7
       gammaL = 10.d0
c      -------------------------Specific heat---------------------------
c       cvion=ccvfion(t,rho,xh,ah,zh) 

       TEMP=1.d-6*t/UN_T6
       GAMAG=B12*UN_B12
       RHONUC=RHO*xh
       call EOSMAG(zh,ah,RHONUC,TEMP,GAMAG,cvion,cve)

c      -------------------------Definitions---------------------------
       a2prime=yn*ah/xh
       aprime=ah + a2prime
       xe=1.007d0*((rho/1.d6*ye)**e13) !electron relat. parameter
       alat=(3.d0*ah*mu/rho/xh/4.d0/pi)**e13
c
       nbaryon=rho/mu       !1/cm3
       tdebye=3.48d3*dsqrt(rho)*zh/ah
       nions=nbaryon*xh/ah  !1/cm3
       ni=nions*1.d-39      !1/fm3
       nel =nbaryon*ye      !1/cm3
       k_F = (3*pi**2*nel)**(e13) !1/cm electron
       
       call kFermi(rho,yn,fnn,kFn)  !n kF in 1/fm
       kFn=kFn*hc                   !n kF in MeV (pF)
       call eff_mass (fnn,0.d0,effmn,effmp) 
c       effmn=1.d0  
       meffn=effmn*931.494d0         !meff in MeV, effmn is only the ratio
              
       gamma = 2.275d5*zh**2*(rho*xh/ah)**e13/t
       a_ni = 10.d0                  !fm neutron-ion scattering length
       Tp = 7.832d3*dsqrt(rho*zh**2/ah/aprime) !plasma T in K
       wp = kb*Tp/hbar               !plasma freq. in 1/s 
       q_BZ=(6.d0*pi**2*nions)**e13  !q of Brioulin zone in 1/cm
       Tu = zh**e13*em*Tp/3.d0       !Umklapp T in K
c
c      ============Transv Chug sound speed========================      
       cs_Chu=wp/3.d0/q_BZ 
c       
c      ============Transv Itoh sound speed========================
       cs_Itoh_t=0.7d0*wp/q_BZ 
c
c      ============Long sound speed======================== 
       cs=clight*dsqrt(ye*me*xe/3.d0/mu)
       ktf=0.1d0*(ye*aprime/2.d0)**e13*q_BZ  !aproximate from Itoh
       ktf=2.d0*dsqrt(dsqrt(1.d0+xe**2)/em/pi/xe)*k_F !1/cm
       ktf=2.d0*dsqrt(em/pi)*k_F  !1/cm
       cs_Itoh_l=dsqrt(wp**2/ktf**2) !cm/s
c
c      ============USED sound speed========================
       cs_t=cs_Itoh_t  !cs_Chu or cs_Itoh_t
       cs_l=cs_Itoh_l  !cs or cs_Itoh_l 
       if ((cs_l.ge.clight).or.(cs_t.ge.clight)) pause 'phonon speed>c'
c
c      ============================================================
c      ----Initialization values ----------------------  
       mfp_sfph_1 = 0.d0
       mfp_sfph_2 = 0.d0
       mfp_sfph_3 = 0.d0
c
c      ==========================================================
c      LATTICE PHONONS
c      ==========================================================
c       -----------Chugunov's phonons------------------------
        call cphonon_Chu(t,zh,ah,Zimp,ye,rho,xh,yn,cvion,
     &      k_tot_Chu,L_ie_Chu)    !in cm        
        L_ie_Chu=L_ie_Chu*1.d13         !in fm
c      -----------Reddy's crude estimate---------------------        
        om=3.d0*t9/1.16d1/1.973d2                     !om=3T in 1/fm
        L_ie_cru=2.d0/pi/om                           !Reddy's naive estimate in fm
        k_cru =cvion*cs_l*L_ie_cru*1.d-13/3.d0        !c.g.s 
c      ---------USED L_ph, kph-----------------------
        fu=dexp(-Tu/t)
        L_lph= 1.d0/(fu/L_ie_Chu+(1.d0-fu)/L_ie_cru)  !interpolation between the two 
        klph =cvion*cs_l*L_lph*1.d-13/3.d0            
c      =================================================================
c      SUPERFLUID PHONONS 
c      =================================================================
       IF ((yn.gt.0.d0)                 !only if free neutrons exist
     -  .and.(taucru.le.1.d0)           !and for gap model T< Tcrit 
     -  ) THEN 
c       -----------------------------------------------
        vs=kFn/dsqrt(3.d0)/meffn       !dimensionless SF phonon speed
        r0=2.d0                        !range of strong int. in fm            
        Rgmix= 1.d0/(1.d0+0.35d0*(kFn/hc)*a_ni
     -   +0.55d0*r0*a_ni*(kFn/hc)**2)
        gmix=dsqrt(4.d0*(kFn/hc)*a_ni**2*ni/ah/mn**2)*Rgmix
        llambda= dsqrt(dsqrt(ni*ah*mn)) !in 1/fm
        llambda = llambda*hc            !in MeV
            
c         -----------interpolation of impurities between liquid-solid
         AA = (175.d0-gammaL)/dsqrt(dlog(1.d2))   !for Zimp=0.01
c          AA = (175.d0-gammaL)/dsqrt(dlog(1.d1))  !for Zimp=0.1
          zzeta =  dexp(-((gamma-gammaL)/AA)**2)   !interpolation 
          if (gamma.le.1.d0) zzeta= 1.d0           !liquid
          if (gamma.ge.175.d0) zzeta= 1.d-2        !solid, for Zimp=0.01
c          if (gamma.ge.175.d0) zzeta= 1.d-1       !solid, for Zimp=0.1          
          dd= alat*1.d13/zzeta**e13        !in fm, distance betw. impurities (dd~100.d0)
          xx=dd/r0                         !diluteness parameter 
c       ----CONDUCTIVITIES---------------------------------------------
        IF (gamma.lt.gammaL) THEN   !pure liquid - region I 
          ksfph=0.d0
        ELSE                         !transition liq-sol or pure solid 
c       ................................................        
          IF (gamma.lt.175.d0) THEN  !transition sol-liq        
           ksfph=0.d0
          ELSE                       !pure solid - region II  
c          ============RAYLEIG=====================           
           ksfph_1=1.d20*(vs/0.1d0)**2*(xx/1.d1)**3/(r0/1.d1)**3/t9  
           mfp_sfph_1=6.80d-6*(vs/0.1d0)**4*(xx/1.d1)**3/(r0/1.d1)**3   !mfp in cm
     -                /t9**4     !mfp in cm
c          ============SPLITTING=================== 
           ccs=cs_t/clight      !dimensionless transv cs
           if (vs.gt.ccs) then  !otherwise kinematically forbidden
           ffalpha=1.d0-2.d0*(ccs/vs)**2/3.d0+(ccs/vs)**4/5.d0
           mfp_sfph_2= 863.d-6* (ccs/0.01d0)**7*(vs/0.1d0)   
     -           *(llambda/5.d1)**4/t7**5/ffalpha/gmix**2               !mfp in cm
c          =============ABSORPTION==================
           ccs=cs_l/clight    !dimensionless long cs       
           omtau_th = om*L_lph/ccs !dimensionless, om*L_lph
           xalpha=ccs/vs
           mfp_sfph_3= (vs**2)/(gmix**2)*
     -      (1.d0+(1.d0-xalpha**2)**2*omtau_th**2) 
     -     *L_lph*1.d-13/(xalpha*omtau_th**2)                            !mfp in cm
c           ===========direct decay==========================
           mfp_sfph_4=2.d0*pi*dexp(taucru)/(3.d0*t9/1.16d1/1.973d2)
c           ==========ALTOGETHER SFph==========================
            mfp_sfph=1.d0/
     -      (1.d0/mfp_sfph_1+1.d0/mfp_sfph_2+1.d0/mfp_sfph_3
     -        +1.d0/mfp_sfph_4)    !over mfp
            ksfph=1.5d19*(t7)**3*(0.1/vs)**2*mfp_sfph
           endif
c         ............................................ 
          ENDIF
        ENDIF
       ENDIF
c      ================================================================          
       ksfph=0.d0     !no superfluid phonons
c      =================================================================
       kph_tot=klph+ksfph 
         
      return
      end


c--------------------------------------------------------------------------------
       real*8 function F_kii_I(rho,zh,ah,xh,yn)  !Coulomb liquid as Chugunov
       implicit none
       real*8 rho,zh,ah,xh,yn
       real*8 e13,pi,kb,mu,hbar,nbaryon,nions,alat
       real*8 a2prime,aprime,Tp,wp
       
        e13=1.d0/3.d0
        pi=3.141592653589793d0
        kb=1.381d-16 !erg/K
        mu=1.66d-24 !gr
        hbar=1.0546d-27  !erg*s
        nbaryon=rho/mu  !1/cm3
        nions=nbaryon*xh/ah
        alat=(3.d0*ah*mu/rho/xh/4.d0/pi)**e13
        
        a2prime=yn*ah/xh
        aprime=ah + a2prime
        Tp = 7.832d3*dsqrt(rho*zh**2/aprime/ah) !plasma T in K
        wp = kb*Tp/hbar   !plasma freq. in 1/s 
        
        F_kii_I=0.4d0*kb*wp*nions*alat**2
        
       return
       end

c--------------------------------------------------------------------------------
       real*8 function F_kphph_II(t,rho,zh,ah,xh,ye,yn,xgamma,cs)  
       implicit none
       real*8 t,rho,zh,ah,xh,ye,yn
       real*8 e13,pi,kb,mu,hbar,rho6,nbaryon,nions,alat,tdebye
       real*8 cvion,mfp,xgamma,cs,ccvfion,fdebye
       real*8 Zimp,k_tot

        e13=1.d0/3.d0
        pi=3.141592653589793d0     
        kb=1.381d-16 !erg/K
        mu=1.66d-24 !gr
        nbaryon=rho/mu  !1/cm3
        tdebye=3.48d3*dsqrt(rho)*zh/ah
        nions=nbaryon*xh/ah
        alat=(3.d0*ah*mu/rho/xh/4.d0/pi)**e13

c       ---PH-PH THERMAL CONDUCTIVITY (high T) --------------
        cvion=ccvfion(t,rho,xh,ah,zh) 
        mfp=ah*mu*cs**2*alat/kb/t/2.d0**2    !cm
        F_kphph_II =cvion*cs*mfp/3.d0        !c.g.s.
        return
        end
c--------------------------------------------------------------------------------------


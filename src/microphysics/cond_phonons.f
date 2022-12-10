! PHONONS CONDUCTIVITY (call it to include SF phonons by means of ksfph)
! Previously called CPHONON and CPHONON_CHU, now in the same file

       subroutine cond_phonons(T8,zh,ah,ye,rho,xh,yn,kFn,effmn,
     &              taucru,B12,ksfph)
        
       implicit none
       real*8 t,T8,rho,ye,ah,zh,xh,yn,taucru
       real*8 klph
       real*8 ksfph_1,ksfph,kph_tot,k_tot_Chu,k_cru
       real*8 gammaL,L_ie_Chu,L_ie_cru,L_lph,omtau_th,
     -        mfp_sfph,mfp_sfph_1,mfp_sfph_2,
     -        mfp_sfph_3, mfp_sfph_4,gmix,om,xalpha,
     -        Rgmix
       real*8 kb,em,me,mu,clight,hbar,hc,mn,pi,e13
       real*8 t9,t7,xe,nbaryon
       real*8 tdebye,gamma,nions,ni,cvion,alat
       real*8 aprime,a2prime,Tp,q_BZ,wp,cs_Chu,cs,cs_l,cs_t,
     -        cs_Itoh_t,cs_Itoh_l,ktf,k_F,nel
       real*8 fu,Tu
       real*8 xx,r0,fnn,kFn,effmn,meffn,vs,ffalpha
       real*8 ccs,llambda
       real*8 dd
       real*8 AA,zzeta,a_ni
       real*8 B12, TEMP, UN_B12, UN_T6, GAMAG, RHONUC, cve
       parameter(UN_B12=425.438, UN_T6=.3157746)

      
c      -------------------------Constants and def.-----------------
       t=T8*1.d8
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
       if ((cs_l.ge.clight).or.(cs_t.ge.clight)) then
         write(*,*) "<warning>[cond_phonons] Phonon speed>c"
       endif
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
        call cond_phonons_Chugunov(t,zh,ah,ye,rho,xh,yn,cvion,
     &         k_tot_Chu,L_ie_Chu)    !in cm        
        L_ie_Chu=L_ie_Chu*1.d13         !in fm
c      -----------Reddy's crude estimate---------------------        
        om=3.d0*t9/1.16d1/1.973d2                     !om=3T in 1/fm
        L_ie_cru=2.d0/pi/om                           !Reddy's naive estimate in fm
        k_cru =cvion*cs_l*L_ie_cru*1.d-13/3.d0        !c.g.s 
c      ---------USED L_ph, kph-----------------------
        fu=dexp(-Tu/t)
c        L_lph= 1.d0/(fu/L_ie_Chu+(1.d0-fu)/L_ie_cru)  !interpolation between the two 
        L_lph= L_ie_cru
        klph =cvion*cs_l*L_lph*1.d-13/3.d0            
c      =================================================================
c      SUPERFLUID PHONONS
c      =================================================================
       ksfph=0.d0 

       IF ((yn .gt. 0.d0)                 !only if free neutrons exist
     -  .and.(taucru .le. 1.d0)           !and for gap model T< Tcrit 
     -  ) THEN 
c       -----------------------------------------------
        vs=kFn*hc/dsqrt(3.d0)/meffn       !dimensionless SF phonon speed
        r0=2.d0                        !range of strong int. in fm            
        Rgmix= 1.d0/(1.d0+0.35d0*(kFn)*a_ni
     -   +0.55d0*r0*a_ni*(kFn)**2)
        gmix=dsqrt(4.d0*(kFn)*a_ni**2*ni/ah/mn**2)*Rgmix
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
c       ksfph=0.d0     !no superfluid phonons
c       klph=0.d0     !no lattice phonons
c      =================================================================
       kph_tot=klph+ksfph 
       
      return
      end


      subroutine cond_phonons_Chugunov(t,zh,ah,ye,rho,xh,yn,Cion,
     &       k_tot,L_ie)
c       Chugunov & Haensel 2007
       implicit none
       double precision t,rho,nbaryon
       double precision ye,ah,zh,xh,yn
       double precision gamma,nions
       double precision xi
       double precision k_ii,k_ie,k_tot,kb,nel,k_F,k_Fless
       double precision clight,mu,me,hbar,pi
       double precision cs,q_BZ,Cion,L_ie,v_F,v_F2,eexp,u1,u2
       double precision w_DW,yy,lambda_ie,XEE1,XEE1_1,expintt
       double precision arg1, arg2,FF
       double precision xe
       double precision rho6,k_0,k_ast,beta,wp,Tp,theta,ai,e13
       double precision aprime, a2prime

       external FF, expintt
       kb=1.381d-16 !erg/K
       me=9.11d-28 !gr
       mu=1.66d-24 !gr
       clight=2.9979d10 !cm/s
       hbar=1.0546d-27  !erg*s
       pi=3.141592653589793d0
       rho6= rho/1.0d6 
       e13=1.d0/3.d0

       xe=1.009d0*((rho6*ye)**(e13)) !electron relat. parameter

       nbaryon=rho/1.66d-24   !1/cm3
       nions=nbaryon*xh/ah    !1/cm3
       gamma = 2.275d5*zh**2*(rho*xh/ah)**(e13)/t

      IF (gamma.le.1.d0) THEN  ! a gas
        k_ii=0.d0
        k_ie=0.d0
        k_tot=0.d0
        return
      ELSE    !a solid or liquid
c      alat=(ah*mu/rho)**(1./3.)
c     ---------phonon-phonon Kii ---------------------------
       k_ast= 0.4d0
       beta = e13
       ai= (4.d0*pi*nions/3.d0)**(-e13) !ion sphere radius (cm)
       a2prime=yn*ah/xh    !check with a2prime=yn*ah
       aprime=ah + a2prime
       Tp = 7.832d6*dsqrt(rho6*zh**2/aprime/ah) !plasma T in K
       theta = Tp/t 
       wp = kb*Tp/hbar     !plasma freq. in 1/s
       k_0 = kb*wp*nions*ai**2
       xi=0.7d0    
       k_ii=k_0*
     - dsqrt(k_ast**2+(xi**3*gamma/2.85d0)**2*dexp(2.d0*beta*theta))
c     ---------phonon-electron Kie----------------------------------
       q_BZ=(6*pi**2*nions)**e13  !q of Brioulin zone in 1/cm
c       cs=wp/3.d0/q_BZ            !sound speed-group velocity
       cs=0.7d0*wp/q_BZ  !Itoh
       
       nel =nbaryon*ye  !1/cm3
       k_F = (3*pi**2*nel)**(e13) !1/cm
       k_Fless =2.d0*k_F/q_BZ          !dimensionless y^{-1}
       v_F=clight*xe/dsqrt(1.d0+xe**2) !electron Fermi veloc. 
       v_F2=(v_F/clight)**2            !dimensionless
       u1=2.8d0
       u2=13.d0
       w_DW=1.683d0*dsqrt(xe/ah/zh)
     -  *(0.5d0*u1*dexp(-9.1d0/theta)+u2/theta)
       yy=1.d0/k_Fless 
       arg1=w_DW*yy**2
       arg2=w_DW      
       eexp=(dexp(-arg1)-dexp(-arg2))/arg2
       XEE1=expintt(1,arg1)
       XEE1_1=expintt(1,arg2)
       lambda_ie = 0.5d0*(XEE1-XEE1_1-v_F2*eexp)
       L_ie=(320.d0*ai/(1.d0+xe**2)/lambda_ie)*(26.d0/zh)
     -  *(FF(theta)/0.01d0)*dsqrt(ah*rho6/aprime)   !mean free path  eq.43     
c       Cion=ccvfion(t,rho,xh,ah,zh)
c      
       k_ie=e13*Cion*cs*L_ie    !*1.d3  if  we put a factor 1e3, their results are 
c       similar to our phonons.  
c      -------------Total-------------------------------
       k_tot=1.d0/(1.d0/k_ii+1.d0/k_ie)

      ENDIF

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccc      
      double precision FUNCTION FF(theta)
      implicit none
      double precision theta

      FF=0.014d0+0.03d0/(dexp(theta/5.d0)+1.d0)
      return
      end
cccccccccccccccccccccccccccccccccc
      double precision FUNCTION expintt(n,xx)
      INTEGER n,MAXIT
      double precision x,EPS,FPMIN,EULER
      PARAMETER (MAXIT=100,EPS=1.e-7,FPMIN=1.e-30,EULER=.5772156649)
      INTEGER i,ii,nm1
      double precision a,b,c,d,del,fact,h,psi,xx
      nm1=n-1
      x=xx

      if(n.lt.0.or.x.lt.0..or.(x.eq.0..and.(n.eq.0.or.n.eq.1)))then
        write(*,*) "<warning>[CPHONON_CHU] Bad arguments in expint"
      else if(n.eq.0)then
        expintt=exp(-x)/x
      else if(x.eq.0.)then
        expintt=1./nm1
      else if(x.gt.1.)then
        b=x+n
        c=1./FPMIN
        d=1./b
        h=d
        do 11 i=1,MAXIT
          a=-i*(nm1+i)
          b=b+2.
          d=1./(a*d+b)
          c=b+a/c
          del=c*d
          h=h*del
          if(abs(del-1.).lt.EPS)then
            expintt=h*exp(-x)
            return
          endif
11      continue
        write(*,*) "<warning>Continued fraction failed in expint"
      else
        if(nm1.ne.0)then
          expintt=1./nm1
        else
         expintt=-log(x)-EULER
        endif
        fact=1.
        do 13 i=1,MAXIT
          fact=-fact*x/i
          if(i.ne.nm1)then
            del=-fact/(i-nm1)
          else
            psi=-EULER
            do 12 ii=1,nm1
              psi=psi+1./ii
12          continue
            del=fact*(-log(x)+psi)
          endif
          expintt=expintt+del
          if(abs(del).lt.abs(expintt)*EPS) return
13      continue
        write(*,*) "<warning>[CPHONON_CHU] Series failed in expint"
      endif
      return
      END


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
        

      end

c--------------------------------------------------------------------------------
       real*8 function F_kphph_II(t,rho,zh,ah,xh,cs)  
       implicit none
       real*8 t,rho,zh,ah,xh
       real*8 e13,pi,kb,mu,nbaryon,nions,alat,tdebye
       real*8 cvion,mfp,cs,ccvfion

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


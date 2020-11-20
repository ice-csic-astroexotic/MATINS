c------------------------------------------------------
c     input: rho[g/cm3],y // output: part.num.dens. fnpart[1/cm3], kFermi kF[fm-1]
c------------------------------------------------------
      subroutine kFermi(rho,ypart,fnpart,kFpart) 
      implicit none
      real* 8 e3,pi,nbaryon,fnpart,rho,ypart,kFpart

      e3=1.d0/3.d0
      pi=3.141592653589793d0
c
      nbaryon=rho/1.66d-24
      fnpart=nbaryon*ypart !nbaryon in cm-3   
      kFpart=(3*pi**2*fnpart)**e3 !kFermi in cm-1
      kFpart=kFpart*1.d-13 !in fm-1
      return
      end

c=======================================================
      subroutine gapmodel(imod,isfty,rho,ypart,gapT0,Tc)
c=======================================================
c     Andersson et al, NPA 763 (2005) 212 
c     & Kaminker, Yakovlev, Gnedin, Astron. Astroph. 383 (2002) 1076 
c=======================================================
      implicit none      
      real*8 kF,gapT0,Tc,delta0,k1,k2,k3,k4
     -       ,alpha,alphaTc,rho, ypart, fnpart
      integer imod,isfty

      call kFermi(rho,ypart,fnpart,kF)

      go to (41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56) imod 
c     -----------------------------------------------------     
c     -------------^1S_0 neutrons -------------------------
c     ---------------------------1)  Bare (A - Andersson 05)
41    delta0 = 9.3D0 
      k1 = 0.02d0
      k2 = 0.6d0
      k3 = 1.55d0
      k4 = 0.32d0
       go to 60
c     --------------------------2)    (a - Andersson 05) 
42    delta0 = 68.d0 
      k1 = 0.1d0
      k2 = 4.d0
      k3 = 1.7d0
      k4 = 4.d0
       go to 60
c     --------------------------3) (b - Andersson 05) 
43    delta0 = 4.d0
      k1 = 0.4d0
      k2 = 1.5d0
      k3 = 1.65d0
      k4 = 0.05d0
       go to 60
c     --------------------------4)  (c - Andersson 05)
44    delta0 = 22.d0
      k1 = 0.3d0
      k2 = 0.09d0
      k3 = 1.05d0
      k4 = 4.d0
       go to 60
c     --------------------------5)  (d - Andersson 05)
45    delta0 = 2.9d0
      k1 = 0.3d0
      k2 = 0.017d0
      k3 = 1.3d0
      k4 = 0.07d0
       go to 60
c     ------------------other references--------------------------
c     --------------------------
c46    delta0 = 12.97d0               !6i)  (Gezerlis and Carlson 08 - old fit)
c      k1 = 0.02d0
c      k2 = 1.41d0                    
c      k3 = 1.55d0
c      k4 = 0.59d0                    
c       go to 60
c     --------------------------
46    delta0 = 82.11d0              !6ii)  (Gezerlis and Carlson 08 - new fit)
      k1 = 0.05d0
      k2 = 2.31d0                 
      k3 = 1.63d0
      k4 = 5.81d0                  
       go to 60
c     --------------------------
c47    delta0 = 21.15d0        !7i)  (Ainsworth et al. 1989) - BC09-like gap
c      k1 = 0.2d0
c      k2 = 1.49d0
c      k3 = 1.7d0
c      k4 = 1.61d0
c       go to 60
c     --------------------------
c47    delta0 = 72.7d0         !7ii)  (Schwenk et al. 03)
c      k1 = 0.1d0
c      k2 = 6.2d0                   !2.49d0**2
c      k3 = 1.5d0
c      k4 = 2.79d0                  !1.67d0**2
c       go to 60
47    delta0 = 68.d0         !7iii)  (Ho et al. 2012)
      k1 = 0.1d0
      k2 = 4.d0
      k3 = 1.7d0
      k4 = 4.d0
      go to 60

c     ---------------------------------------------------------------
c     --------------^1S_0 protons------------------------------------
c     ------------from Andersson et al.------------------------------
c     --------------------------
48    delta0 = 120.d0             !8)  (Ho et al. 2012)
      k1 = 0.d0
      k2 = 9.d0
      k3 = 1.3d0
      k4 = 1.8d0
      go to 60

c48    delta0 = 61.d0             !8)  (e - Andersson 05)
c      k1 = 0.d0
c      k2 = 6.d0
c      k3 = 1.1d0
c      k4 = 0.6d0
c       go to 60
c     --------------------------
49    delta0 = 55.d0            !9)  (f - Andersson 05)
      k1 = 0.15d0
      k2 = 4.d0
      k3 = 1.27d0
      k4 = 4.d0
       go to 60
c     --------------------------
50    delta0 = 2.27d0           !10)  (g - Andersson 05)
      k1 = 0.1d0
      k2 = 0.07d0
      k3 = 1.05d0
      k4 = 0.25d0
       go to 60
c     ---------------------------------------------------------------
c     ---------------------^3P_2 neutrons----------------------------
c     --------------------------------------------------------------
51    delta0 =0.068d0            !11)  (Ho et al. 2012, shallow)
      k1 = 1.28d0
      k2 = 0.1d0
      k3 = 2.37d0
      k4 = 0.02d0
       go to 60
c     --------------------------
52    delta0 =0.15d0           !12) (Ho et al. 2012, deep)
      k1 = 2.d0
      k2 = 0.1d0
      k3 = 3.1d0
      k4 = 0.02d0
       go to 60

c51    delta0 = 4.8d0            !11)  (h - Andersson 05)
c      k1 = 1.07d0
c      k2 = 1.8d0
c      k3 = 3.2d0
c      k4 = 2.d0
c       go to 60
c     --------------------------
c52    delta0 =10.2d0           !12) (i - Andersson 05)
c      k1 = 1.09d0
c      k2 = 3.d0
c      k3 = 3.45d0
c      k4 = 2.5d0
c       go to 60
c     --------------------------
53    delta0 = 2.2d0          !13)  (j - Andersson 05)
      k1 = 1.05d0
      k2 = 1.d0
      k3 = 2.82d0
      k4 = 0.6d0
       go to 60
c     --------------------------
54    delta0 = 0.425d0        !14) (k - Andersson 05)
      k1 = 1.1d0
      k2 = 0.5d0
      k3 = 2.7d0
      k4 = 0.5d0
       go to 60
c     --------------------------
55    delta0 = 0.068d0         !15) (l - Andersson 05)
      k1 = 1.28d0
      k2 = 0.1d0
      k3 = 2.37d0
      k4 = 0.02d0
       go to 60
c     --------------------------
56    delta0 = 2.9d0           !16) (m - Andersson 05)
      k1 = 1.21d0
      k2 = 0.5d0
      k3 = 1.62d0
      k4 = 0.5d0
       go to 60
c     --------------------------
60    continue
c     -----------------
      if ((kF.gt.k1).and.(kF.le.k3)) then
      gapT0=delta0*(kF-k1)**2*(kF-k3)**2
     -      /((kF-k1)**2+k2)/((kF-k3)**2+k4) !in MeV
      else
       gapT0=0.d0
      endif
      alpha=alphaTc(isfty)
      Tc=alpha*gapT0*1.d11/8.617342d0 !Tc in K
c
      return
      end
c=======================================================
      real*8 function alphaTc(isfty)
      implicit none
      integer isfty
   
      if (isfty.eq.1) alphaTc=0.5669d0
      if (isfty.eq.2) alphaTc=0.8416d0
      if (isfty.eq.3) alphaTc=0.4926d0

      return
      end

c===============================================================
c     Dimensionless Gap amplitudes-Temperature dependence
c===============================================================
      real*8 function v1(tau) !singlet-pairing SF
      implicit none   
      real*8 tau
      if (tau.le.0.1D0) then 
      v1 = 1/0.5669d0/(tau+1.d-10)
      else
      v1 = dsqrt(1.d0-tau)*(1.456d0-0.157d0/dsqrt(tau)+1.764d0/tau) 
      endif
      RETURN
      END

c===============================================================
      real*8 function v2(tau) !triplet-pairing ^3P_2 (m_j=0)
c===============================================================
      implicit none   
      real*8 tau
      if (tau.le.0.1D0) then
      v2 = 1/0.8416d0/(tau+1.d-10)
      else 
      v2 = dsqrt(1.d0-tau)*(0.7893d0+1.188d0/tau)  !idem v2
      endif
      RETURN
      END

c===============================================================
      real*8 function v3(tau)  !triplet-pairing ^3P_2 (|m_j|=2)   
c===============================================================
      implicit none   
      real*8 tau
      if (tau.le.0.1D0) then
      v3 = 1/0.4926d0/(tau+1.d-10)    
      else
      v3 = dsqrt(1.d0-tau**4)/tau*
     -     (2.03d0-0.4903d0*tau**4+0.1727d0*tau**8) 
      endif
      return
      END
c----------------------------------------------------


c------------------------------------------------------
c--------------DUrca Reduction functions---------------
c------------------------------------------------------
        Real*8 function RRD_A(v) !for 1S_0
      implicit none
        real*8 v,zexp1

        zexp1=3.427d0-dsqrt(3.427d0**2+v**2)
        if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff
        RRD_A=(0.2312d0+dsqrt(0.7688d0**2+(0.1438d0*v)**2))**(5.5d0)*
     &         dexp(zexp1)    
        return
        end
c------------------------------------------------------
        Real*8 function RRD_B(v) !for 3P_2 (m_j=0)
      implicit none
        real*8 v,zexp1

        zexp1=2.701d0-dsqrt(2.701d0**2+v**2)
        if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff
        RRD_B=(0.2546d0+dsqrt(0.7454d0**2+(0.1284d0*v)**2))**(5.d0)*
     &         dexp(zexp1)
        return
        end
c------------------------------------------------------
        Real*8 function RRD_C(v) !for 3P_2 (|m_j|=2)
      implicit none
        real*8 v,zexp1

        zexp1=1.d0-dsqrt(1.d0+(0.4129d0*v)**2)
        if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff        
        RRD_C=(0.5d0+(0.09226d0*v)**2)/
     -        (1.d0+(0.1821d0*v)**2+(0.16736d0*v)**4)+
     &         0.5d0*dexp(zexp1)
        return
        end
c-------------------------------------------------------
c----------Combined superfluidity  3P_2 neutrons (nB) & 1S_0 protons (pA) 
        Real*8 function RRD_BA(xv1,xv2)  !xv1 neutrons , xv2 protons
      implicit none
        real*8 xv1,xv2,xv12,xv22,xv14,xv24,RD_B,RD_A, RRD_B, RRD_A
        xv12=xv1**2
        xv22=xv2**2
        xv14=xv1**4
        xv24=xv2**4
        if ((xv12+xv22).le.2.5d1) then !not strong SF
        RRD_BA=(1.d4-2.839d0*xv24-5.022d0*xv14)/
     &   (1.d4 + 7.57d2*xv22 + 1.494d3*xv12
     &    +2.111d2*xv22*xv12 + 0.4832d0*xv24*xv14)
        else               !strong SF
        RD_B=RRD_B(xv1) !n in 3P_2
        RD_A=RRD_A(xv2) !p in 1S_0
        RRD_BA=min(RD_B,RD_A)
        endif 
        return
        end
c------------------------------------------------------
c--------------MUrca Reduction functions---------------
c------------------------------------------------------
        Real*8 function RRMn_pA(v) !neutron branch for 1S_0 protons
      implicit none
        real*8 v,a,b,zexp1
c

        a=0.1477d0+dsqrt(0.8523d0**2+(0.1175d0*v)**2)
        b=0.1477d0+dsqrt(0.8523d0**2+(0.1297d0*v)**2)
        zexp1=3.4370d0-dsqrt(3.4370d0**2+v**2)
        if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff 
      
        RRMn_pA=0.5d0*(a**(7.5d0)+b**(5.5d0))*dexp(zexp1)  
        return
        end
c------------------------------------------------------
        Real*8 function RRMp_pA(v) !proton branch for 1S_0 protons
      implicit none
        real*8 v,c,zexp1
c
        c=0.2414d0+dsqrt(0.7586d0**2+(0.1318d0*v)**2)
        zexp1=5.339d0-dsqrt(5.339d0**2+(2.d0*v)**2)

c        if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff        
        RRMp_pA=c**7*dexp(zexp1) 
        return
        end
c------------------------------------------------------
c------------RRMn_nA & RRMp_nA obtained from RRMp_pA &  RRMn_pA, respec.
        Real*8 function RRMn_nA(v) !neutron branch for 1S_0 neutrons
      implicit none
        real*8 v, RRMp_pA
c
        RRMn_nA=RRMp_pA(v)
        return
        end
c---------------------------------------------------------
        Real*8 function RRMp_nA(v) !proton branch for 1S_0 neutrons
      implicit none
        real*8 v, RRMn_pA
c
        RRMp_nA=RRMn_pA(v)
        return
        end
c----------------------------------------------------------
c-------Similarity criteria for neutron brance for 3P_2-------
        Real*8 function RRMn_nB(v) !neutron branch for 3P_2 neutrons
      implicit none
        real*8 v, xsi,pi,RRMp_pA,zexp1
c
        if (v.gt.1.5d1) then !limit of strong superfluidity
        xsi= 0.209d0
        pi=3.1415927d0
        zexp1=-2.d0*v
        if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff

        RRMn_nB= 1.20960d0*2.d0*xsi/0.11513d0/pi**8/3.d0/1.7320508d0
     &      *v**6*dexp(zexp1) 
        else
        RRMn_nB=RRMp_pA(v) ! valid only for moderate superfluidity, v=<10
        endif
        return
        end
c---------------------------------------------------------
        Real*8 function RRMp_nB(v) !proton branch for 3P_2 neutrons
      implicit none
        real*8 v,a,b,zexp1

        a=0.1612d0+dsqrt(0.8388d0**2+(0.1117d0*v)**2)
        b=0.1612d0+dsqrt(0.8388d0**2+(0.1274d0*v)**2)
        zexp1=2.398d0-dsqrt(2.398d0**2+v**2)
        if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff

        RRMp_nB=0.5d0*(a**7+b**5)*dexp(zexp1)  
        return
        end
c----------------------------------------------------------
c----------Combined superfluidity  3P_2 neutrons (nB) & 1S_0 protons (pA) 
        Real*8 function RRMn_BA(xvn,xvp)  !neutron branch
      implicit none
        real*8 xvn,xvp,RRD_BA,RRMn_pA,RRD_A

        RRMn_BA=RRD_BA(2.d0*xvn,xvp)*RRMn_pA(xvp)/RRD_A(xvp)
        return
        end
c       --------------------------------------------------
        Real*8 function RRMp_BA(xvn,xvp)  !neutron branch
      implicit none
        real*8 xvn,xvp,RRD_BA,RRMp_nB,RRD_B

        RRMp_BA=RRD_BA(xvn,2.d0*xvp)*RRMp_nB(xvn)/RRD_B(xvn)
        return
        end
c------------------------------------------------------
c--------------Neutrino Bremsstrahlung Reduction functions---------------
c------------------------------------------------------
        Real*8 function RRnp_pA(v) !np for 1S_0 protons
      implicit none
        real*8 v,a,b,zexp1,zexp2

        a=0.9982d0+dsqrt(0.0018d0**2+(0.3815d0*v)**2)
        b=0.3949d0+dsqrt(0.6051d0**2+(0.2666d0*v)**2)
        zexp1=1.306d0-dsqrt(1.306d0**2+v**2)
        zexp2=3.303d0-dsqrt(3.303d0**2+4.d0*v**2)
        if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff
        if (zexp2.lt.-1.5d2) zexp2=-1.5d2

        RRnp_pA=(a*dexp(zexp1)+1.732d0*b**7*dexp(zexp2))/2.732d0  
        return
        end
c------------------------------------------------------
        Real*8 function RRpp_pA(v) !pp for 1S_0 protons
      implicit none
        real*8 v,c,d,zexp1,zexp2

        c=0.1747d0+dsqrt(0.8253d0**2+(0.07933d0*v)**2)
        d=0.7333d0+dsqrt(0.2667d0**2+(0.1678d0*v)**2)
        zexp1=4.228d0-dsqrt(4.228d0**2+4.d0*v**2)
        zexp2=7.762d0-dsqrt(7.762d0**2+9.d0*v**2)
        if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff
        if (zexp2.lt.-1.5d2) zexp2=-1.5d2
c
        RRpp_pA=(c**2*dexp(zexp1)+d**(7.5d0)*dexp(zexp2))/2.d0  
        return
        end
c--------------------------------------------------------
        Real*8 function RRnp_nA(v) !np for 1S_0 neutrons
      implicit none
        real*8 v, RRnp_pA

        RRnp_nA=RRnp_pA(v)
        return
        end
c--------------------------------------------------------
        Real*8 function RRnp_nB(v) !np for 3P_2 neutrons!!!!NOT DEFINED!!!!
      implicit none
        real*8 v, RRnp_nA            ! I define it by similarity criterium

        RRnp_nB=RRnp_nA(v)
        return
        end
c--------------------------------------------------------
        Real*8 function RRnn_nA(v) !nn for 1S_0 neutrons
      implicit none
        real*8 v, RRpp_pA

        RRnn_nA=RRpp_pA(v)
        return
        end
c--------------------------------------------------------
        Real*8 function RRnn_nB(v) !nn for 3P_2 neutrons-
      implicit none
        real*8 v, RRpp_pA

        RRnn_nB=RRpp_pA(v) !Similarity criterium
        return
        end
c--------------------------------------------------------
c----------Combined superfluidity  3P_2 neutrons (nB) & 1S_0 protons (pA) 
        Real*8 function RRnp_BA(xvn,xvp) !xvn for n, xvp for p
      implicit none
        real*8 xvn,xvp,RRD_BA,RRnp_pA,RRD_A

        RRnp_BA=RRD_BA(xvn,xvp)*RRnp_pA(xvp)/RRD_A(xvp)
        return
        end 
c====================================================================
c      ----end functions for SF---- 
c====================================================================


c====================================================================
c    New functions added by Andrea Passamonti, for ambipolar diffusion 
c====================================================================
c       --------------------------------------------------------	
c	Functions for p-n collision Rpn
c       --------------------------------------------------------

c       --------------------------------------------------------------
c       Sp1  (Baiko et al. 2001) eq (29)
        doubleprecision function Sp_1(kn,kp)  !when both are superfluid
        implicit none
        doubleprecision kn,kp

        Sp_1 = 0.8007 * kp / kn**2
     &       * ( 1d0 + 31.28 * kp - 0.0004285 * kp**2
     &         + 26.85 * kn + 0.08012 * kn**2 ) 
     &       * ( 1d0 - 0.5898 * kn + 0.2368 * kn**2 
     &         + 0.5838 * kp**2 + 0.884 * kn * kp 
     &         )**(-1)   !! it is given in mb = 1d-17 cm^2



        return
        end


c       Sp2  (Baiko et al. 2001) eq (29)
        doubleprecision function Sp_2(kn,kp)  !when both are superfluid
        implicit none
        doubleprecision kn,kp

        Sp_2 = 0.3830 * kp**4 / kn**(5.5)
     &       * ( 1d0 + 102.0 * kp + 53.91 * kn ) 
     &       * ( 1d0 - 0.7087 * kn + 0.2537 * kn**2 
     &         + 9.404 * kp**2 - 1.589 * kn * kp 
     &         )**(-1)   !! it is given in mb = 1d-17 cm^2



        return
        end



c       --------------------------------------------------------------
c       Kp1  (Baiko et al. 2001) eq (30)
        doubleprecision function Kp_1(kn,kp)  !when both are superfluid
        implicit none
        doubleprecision kn,kp
        doubleprecision u
        doubleprecision mp_star_mN

        mp_star_mN = 0.8

        u    = kn - 2.126 

        Kp_1 = 1d0/mp_star_mN**2 
     &       * ( 0.04377  + 1.100 * u**2 + 0.1180 * u**3
     &       +   0.1626 * kp + 0.3871 * u * kp - 0.2990 * u**4 )



        return
        end


c       --------------------------------------------------------------
c       Kp2  (Baiko et al. 2001) eq (30)
        doubleprecision function Kp_2(kn,kp)  !when both are superfluid
        implicit none
        doubleprecision kn,kp
        doubleprecision u
        doubleprecision mp_star_mN

        mp_star_mN = 0.8

        u    = kn - 2.116 

        Kp_2 = 1d0/mp_star_mN**2 
     &       * ( 0.0001313 + 1.248 * u**2 + 0.2403 * u**3
     &       +   0.3257 * kp + 0.5536 * u * kp - 0.3237 * u**4
     &       +   0.09786 * u**2 * kp )



        return
        end


c       --------------------------------------------------------------
c       Only neutron superfluid  (Baiko et al. 2001) eq (46)

c       1. Rp1
        doubleprecision function RRp_1n(xvn)  !neutron superfluid
        implicit none
        doubleprecision xvn,term_a

        term_a = 0.4459   
     &         + dsqrt( 0.5541**2 + 0.03016 * xvn**2) 

        RRp_1n = term_a**2 * exp(2.1178 - dsqrt(2.1178**2 + xvn**2))

        return
        end


c       2. Rp2
        doubleprecision function RRp_2n(xvn)  !neutron superfluid
        implicit none
        doubleprecision xvn,term_a

        term_a = 0.801
     &         + dsqrt( 0.199**2 + 0.04645 * xvn**2) 

        RRp_2n = term_a**2 * exp(2.3569 - dsqrt(2.3569**2 + xvn**2))

        return
        end

c       --------------------------------------------------------------
c       Only protons superfluid (Baiko et al. 2001) eq (46)


c	1. Rp1 
        doubleprecision function RRp_1p(xvp)  !proton superfluid
        implicit none
        doubleprecision xvp,term_a,term_b

        term_a = 0.5d0 * ( 0.3695  
     &         + dsqrt( 0.6305**2 + 0.01064 * xvp**2) )

        term_b = 0.5d0 * (1d0 + 0.1917 * xvp**2 )**1.4d0

        RRp_1p = term_a * exp(2.4451 - dsqrt(2.4451**2 + xvp**2))
     &         + term_b * exp(4.6627 - dsqrt(4.6627**2 + 4d0 * xvp**2))


        return
        end


c	2. Rp2 
        doubleprecision function RRp_2p(xvp)  !proton superfluid
        implicit none
        doubleprecision xvp,term_a

        term_a = 0.0436d0 * ( - 3.345d0
     &         + dsqrt( 4.345**2 + 19.55 * xvp**2) )

        RRp_2p = term_a * exp(2.0247 - dsqrt(2.0247**2 + xvp**2))
     &         + 0.0654 * exp(8.992  - dsqrt(8.992**2  + 1.5d0*xvp**2))
     &         + 0.891  * exp(9.627  - dsqrt(9.627**2  + 9d0*xvp**2))

        return
        end


c       --------------------------------------------------------------
c       when both species superfluid (Baiko et al. 2001) eq (47)

c       1. Rp_1
        doubleprecision function RRp_1(xvn,xvp)  !when both are superfluid
        implicit none
        doubleprecision xvn,xvp,xv_min,xv_max
        doubleprecision u_n,u_p,u_min,u_max
        doubleprecision term_a,term_b

        xv_min = min(xvn,xvp)
        xv_max = max(xvn,xvp)

        u_n    = dsqrt( xvn**2    + (1.485)**2 ) - 1.485
        u_p    = dsqrt( xvp**2    + (1.485)**2 ) - 1.485
        u_min  = dsqrt( xv_min**2 + (1.485)**2 ) - 1.485
        u_max  = dsqrt( xv_max**2 + (1.485)**2 ) - 1.485

        term_a = 0.7751 + 0.4823 * u_n + 0.1124 * u_p
     &         + 0.04991  * u_n**2
     &         + 0.08513  * u_n * u_p
     &         + 0.01284  * u_n**2 * u_p

        term_b = 0.2249 + 0.3539 * u_max - 0.2189 * u_min
     &         - 0.6069   * u_n * u_min
     &         + 0.7362   * u_p * u_max

        RRp_1  = term_a * exp(- u_min - u_max) 
     &         + term_b * exp(-2d0 * u_max)  


        return
        end


c       2. Rp_2
        doubleprecision function RRp_2(xvn,xvp)  !when both are superfluid
        implicit none
        doubleprecision xvn,xvp,xv_min,xv_max
        doubleprecision u_n,u_p,u_min,u_max
        doubleprecision term_a,term_b

        xv_min = min(xvn,xvp)
        xv_max = max(xvn,xvp)

        u_n    = dsqrt( xvn**2    + (1.761)**2 ) - 1.761
        u_p    = dsqrt( xvp**2    + (1.761)**2 ) - 1.761
        u_min  = dsqrt( xv_min**2 + (1.761)**2 ) - 1.761
        u_max  = dsqrt( xv_max**2 + (1.761)**2 ) - 1.761

        term_a = 1.1032 + 0.8645 * u_n + 0.2042 * u_p
     &         + 0.07937  * u_n**2
     &         + 0.1451   * u_n * u_p
     &         + 0.01333  * u_n**2 * u_p

        term_b = - 0.1032 - 0.2340 * u_max 
     &         +   0.06152  * u_n  * u_max
     &         +   0.7533   * u_n  * u_min
     &         -   1.007    * u_p  * u_max

        RRp_2  = term_a * exp(- u_min - u_max) 
     &         + term_b * exp(-2d0 * u_max)  

        return
        end



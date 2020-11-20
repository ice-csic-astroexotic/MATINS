c--------------------------------------------------------------------
c   This program calculates the electron and neutron conductivity 
c   of dense matter in the core of neutron stars
c   using the expressions given in
c     Gnedin & Yakovlev 1995 (Nuc. Phys. A 582, 697)
c     Baiko, Haensel & Yakovlev 2001 (A&A 374, 151)
c
c   Update (now commented) 2012:
c   Fluxtube-electron scattering from
c     Ruderman, Zhu & Cheng 1998 (ApJ 492, 267)
c--------------------------------------------------------------------
      subroutine en_cond(B,t8,rho,ye,taun,taup,tcond,eta)
      use constants
      implicit none
      real*8 rho,ye,yn,yp,taun,taup,tcond,econd,ncond, rhon, rhop
      real*8 t8, ne, nb, effme, effmn, effmp, eta
      real*8 Rn, Rp, kfe, q0_kfe, sfn, sfp, Zn, Zp
      real*8 nuee, nuen, nuep, nue, ttaue
      real*8 nunn, nunp, nun, ttaun, RC,xvn,xvp,v1,v2
!     parameter (C=1.2d0) 
      real*8 rfun, zfun, rcfun
      real*8 np2, alpha, n0
      real*8 mp, mel, b
      !real*8 nphi, phi0, nuft, lambda, enfe

      external rfun, zfun, rcfun

      yp = ye         ! adim
      yn = 1.d0-yp    ! adim
      nb = rho/RHO_TO_N    ! from rho[g/cm^3] to fm^-3
      n0 = 0.16d0     ! in fm^-3 Nuclear saturation density
      ne = nb*ye      ! in fm^-3 
      kfe = (3.d0*PI**2*ne)**(1.d0/3.d0)   ! in fm^-1
      np2 = ne
      ALPHA=1./137.d0

! mass in units of c^2/(hbar*c) [fm^-1]
      mel   = MASS_E_MEV/hbarc   ! in fm^-1
      mp    = MASS_N_MEV/hbarc   ! in fm^-1

      effme = dsqrt(1.d0+(kfe/mel)**2)  ! fm^-1*197.32696 MeV fm = m_e [MeV]
! effmn=meffn/mn (normalized meffn), effmp=meffp/mp (normalized meffp)
      rhon = rho/MASS_N*yn
      rhop = rho/MASS_N*yp
      call eff_mass (rhon,rhop,effmn,effmp)

! RUDERMAN FLUXTUBE-SCATTERING RESISTIVITY (added in 2012)
!       LAMBDA = effmp*mp/(4.d0*pi*alpha*np2) 	! in fm (to obtain proton mass, I have to multiply normalized mass meffp*mp
!       PHI0 = 2d-7*1d30			  	! in G cm^2 = 1d30 G fm^2
!       NPHI = b*1d12/phi0			! in fm^-2 (B in units of 1d12 G)
!       ENFE = 0.5d0*kfe**2			! in fm^-2
! energy in units of hbar*c, c=2.9979d23 fm/s
!       NUFT = 3.d0*pi**3*nphi*2.9979d23/
!      &		(64.d0*lambda*enfe*alpha*ne)

c--------------------------------------------------------------------
c       electron conductivity
c--------------------------------------------------------------------
      Zn = 1.d0
      Zp = 1.d0
      Rn = 1.d0
      Rp = 1.d0
      RC = 1.d0

      if (taun.le.1.d0) then  !neutrons SF
            xvn = v2(taun) !3P_2   
            Zn = zfun(xvn)
            Rn = rfun(xvn)
            RC = rcfun(xvn) 
      endif
      if (taup.le.1.d0) then  !protons SF
            xvp = v1(taup) !1S_0  
            Zp = zfun(xvp) 
            Rp = rfun(xvp)
      endif   

      sfn = effmn*(yn*nb*n0/ne**2)**(1.d0/3.d0)*Zn
      sfp = effmp*(yp*nb*n0/ne**2)**(1.d0/3.d0)*Zp

      q0_kfe = dsqrt(0.00929*(1.d0+ 2.83*(sfn+sfp)))

      nuee = 3.58d11/(q0_kfe**3)/((ne/n0)**(1.0/3.0))*t8**2  ! in s^-1
      nuen = 1.15d12/(q0_kfe**3)*effmn**2/(ne/n0)*Rn*t8**2   ! in s^-1
      nuep = 1.15d12/(q0_kfe**3)*effmp**2/(ne/n0)*Rp*t8**2   ! in s^-1
      
      nue  = nuee + nuen + nuep! + nuft
      ttaue = 1.d0/nue*7.763d20  ! electron relaxation time in units of hbar=c=m_e=1
      econd = 1.7d24*1.2*t8*(1.d15/nue)*(ne/n0)**(2.0/3.0)  ! Eq.(65) Gnedin & Yakovlev erg / cm s K

c--------------------------------------------------------------------
C       neutron conductivity
c--------------------------------------------------------------------

c--------------------------------------------------------------------
c The following superfluid and interaction corrections are set to 1,
c because they are very complicated and irrelevant, and the factor
c RC already takes into account SF suppression of n cond.
c--------------------------------------------------------------------

      sfn=1.d0
      sfp=1.d0
  
      nunn = 3.48d15*effmn**3*t8**2*sfn
      nunp = 3.48d15*effmn*effmp**2*t8**2*sfp
 
      nun  = nunn + nunp 
      ttaun = 1.d0/nun*1.d15  ! relaxation time in 10^-15 s
      ncond = 7.2d23*t8/effmn*RC**2*ttaun*(nb/n0*yn)  ! Eq. (52) Baiko-Haensel-Yakovlev 2001

! CONDUCTIVITY (35) Aguilera 2008
! ncond = pi^2*kb^2*n*T*teff/(3*meff)

c--------------------------------------------------------------------
c       total conductivity
c--------------------------------------------------------------------
      tcond = econd + ncond
c--------------------------------------------------------------------
c       resistivity
c--------------------------------------------------------------------
      eta = effme/(1.445e3*ne*ttaue)  ! in Km^2/Myr

      return
      end

c--------------------------------------------------------------------
      real*8 function rfun(y)
c--------------------------------------------------------------------
      real*8 y, y2
      y2 = y*y
      rfun = ( 0.7694 + dsqrt(0.2306**2 + (0.07207*y)**2)
     &       + y2*(27.0+0.1476*y2)*dexp(-dsqrt(4.273**2+y2))
     &       + 0.5051*(dexp(4.273-dsqrt(4.273**2+y2))-1.d0) )
     &      *dexp(1.187-dsqrt(1.187**2+y2))
      return
      end 

c--------------------------------------------------------------------
      real*8 function zfun(y)
c--------------------------------------------------------------------
      real*8 y, y2
      y2 = y*y
      zfun = dsqrt( 0.9443+dsqrt(0.0557**2+(0.1886*y)**2) )
     &      *dexp(1.753-dsqrt(1.753**2+y2))
      return
      end 

c--------------------------------------------------------------------
      real*8 function rcfun(y)
c--------------------------------------------------------------------
      real*8 y, y2
      y2 = y*y
      rcfun = ( 0.647+dsqrt(0.353**2+0.109*y2) )**(1.5)
     &      *dexp(1.39-dsqrt(1.39**2+y2))
      return
      end

      subroutine tensor(x,bm,rho,ye,ah,zh,xh,yn,
     &     taun,taup,ffanis,kt,kh,eta,shearm,sigma_max,
     &	   Qimp,Qpasta)
      
      use constants
      
      implicit none

      include '../decl/dim2.h'
      real*8 x,bm,rho,ye,ah,zh,xh,yn
      real*8 t,rr,zz,aa
      real*8 tau,taun,taup
      real*8 tcond,tcondt,tcondh
      real*8 kpar,ktrans,kt,kh,kfe,ffanis
      real*8 Qimp,Qpasta,Zimp,nel,meff,eta,rsigma
      real*8 mn,nion,alat,shearm,sigma_max
      real*8 qj,qjt,qjh
      real*8 rhopasta,rhotrans

      t=x*1.d8
      rr=rho
      zz=zh
      aa=ah
! -----------------------------------------------------------
!	CORE: en_cond
!	CRUST: potekhin2019 (condconv, condegin)
! -----------------------------------------------------------

      IF (xh.eq.0) then
        call en_cond(bm,t,rr,ye,taun,taup,tcond,tau)
        tcondt = tcond
        ffanis=0.d0
        kpar=tcond
        ktrans=tcondt
        kt=ktrans*1.d13*1.d-40   		! cgs to [10^40 erg/(s*10^8 K*km)]
        kh=0.d0
        nel = (rho/RHO_TO_N)*ye    		! electron density in fm^-3 (1.67d15=1.67d-24 g*(1d13 fm/cm)^3)
        kfe = HBARC*(nel*3.0*PI**2)**(1.d0/3.d0)
        meff = dsqrt(1.d0+(kfe/MASS_E_MEV)**2)      ! effective mass
        eta = meff/(1.445e3*nel*tau)  		! in Km^2/Myr
      ELSE
! -----------------------------------------------------------
!     	Shear modulus and maximum strength
! -----------------------------------------------------------
        call shear(x,rho,ah,zh,xh,shearm,sigma_max)
! -----------------------------------------------------------
! 	impurity parameter (input for Potekhin's routine)
! -----------------------------------------------------------
!	Linear interpolation for Qimp
! -----------------------------------------------------------
      rhopasta=8.e13
      rhotrans=1.e13
      	if (rho .lt. rhopasta .and. rho .gt. rhotrans) then
          	Zimp=dsqrt(Qimp+
     &	(rho-rhotrans)**2/(rhopasta-rhotrans)**2*(Qpasta-Qimp))
	      elseif (rho .ge. 8.e13) then
	          Zimp=dsqrt(Qpasta)
	      else
        	  Zimp=dsqrt(Qimp)
	      endif

!    electron + phonons conductivity (not SF phonons)
        call potekhin2019(t,rr,bm,ah,zh,xh,Zimp,
     &   tcond,tcondt,tcondh,rsigma,qj,qjt,qjh)
! -----------------------------------------------------------
!     JONES AMORPHOUS CRUST
! Linear fit taking Q from Jones 2004, Tab. I, sp
!   	Zimp=sqrt(6.3+1.3e-13*rho)
! Linear fit taking Q from Jones 2004, Tab. I, p
!   	Zimp=sqrt(4.1+2.4e-13*rho)
! -----------------------------------------------------------

! -----------------------------------------------------------
!     JONES AMORPHOUS CRUST
! -----------------------------------------------------------
! Linear fit taking resistivity R from Jones 2004, Tab. II, sp
! 	if (rho .ge. 1.e13) rsigma=1.d0/(8.2d-25 -5.8e-39*rho)
! Linear fit taking resistivity R from Jones 2004, Tab. II, p
! 	if (rho .ge. 1.e13) rsigma=1.d0/(7.3d-25 +4.6e-40*rho)
! -----------------------------------------------------------

! -----------------------------------------------------------
!    phonon conductivity (to include also SF phonon)
!       call cphonon(t,zh,ah,Qimp,ye,rr,xh,yn,taucru,bm,tcondph)
!    radiative (photon) conductivity (unimportant in the crust)
!       call ykurrad(t,ye,rr,bm,tcondp,tcondpt)
! -----------------------------------------------------------
        kpar=tcond    	! +tcondp
        ktrans=tcondt  	! +tcondpt

!         ktrans=dmax1(ktrans,1.d-2*kpar)
        if (bm.gt.0.d0) then
          ffanis=(kpar/ktrans-1.d0)/bm**2
          kh=(tcondh/ktrans)/bm
        else
          ffanis=0.d0
          kh=0.d0
        endif
        kt=ktrans*1.d8*1.d5*1.d-40  		! cgs to [10^40 erg/(s*10^8 K*km)]
        eta=CLIGHT**2/(4.d0*PI*rsigma)*T_YEAR*1.D-4 	! in Km^2/Myr

      ENDIF

      return
      end



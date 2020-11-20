!-------------------------------------------------------------------------------
!> @brief This subroutine computes the heat capacity. In the core we use 
!> semi-analytical fits for each component (protons, neutrons, electrons)
!> and inthe crust we call EOSMAG (potekhin)
!
!> @param[in] kmax       angular dimension of the arrays.
!> @param[in] lmax       radial dimension of the arrays.
!> @param[in] tem0       Temperature values for each cell [10^8 K]
!> @param[in] bmed       Modulus of the magnetic field [10^12 G]
!> @param[out] c_v       Heat capacity at constant volume [10^40 erg/(km^3*10^8 K)]
!
!> @author
!> Jose Pons Botella
!
  subroutine compute_heat_capacity
  use grid, only: rho,xh,ye,yn,yp,aa,zz,tccru,tcn,tcp  
  use grid, only: kmax, lmax
  use grid, only: tem0, bmed, c_v

  implicit none 
  integer k, l, j
  real*8 tem, cv_cgs, taucru,taun,taup

    do l=1,lmax
      j = 2*l-1
      do k=1,kmax
        tem = tem0(k,l)    ! in 10**8 K
!-----------------------------------------------------------------
!      tau = Tcrit/T, it tauX>1 SF corrections are not applied
!-----------------------------------------------------------------
        taucru=tem/dmax1(tccru(j),0.1d0*tem) 
        taun=tem/dmax1(tcn(j),0.1d0*tem)
        taup=tem/dmax1(tcp(j),0.1d0*tem)

        call cv_func(tem,rho(j),ye(j),aa(j),zz(j),xh(j),yn(j),yp(j), &
                   &   taucru,taun,taup,bmed(k,l),cv_cgs)    
        c_v(k,l) = 1.d-17*cv_cgs !! erg/K/cm3  ---> 10**40 erg/km^3/10**8 K
      enddo
    enddo

  return
  end subroutine compute_heat_capacity

!     This function calculates the heat capacity cv
!     Van Riper 1991, ApJ 449-462
!
  subroutine cv_func(t8,rho,YE,AH,ZH,XH,YN,yp,taucru,taun,taup,b12,cv_cgs)
    use constants, only : MASS_N
    implicit none
    real*8 cv_cgs
    real*8 t,t8,rho,YE,AH,ZH,XH,YN,yp,taucru,taup,taun,b12
    real*8 cve,cvion,cvn,cvp,ccvfn,ccvfp,ccvfe!,ccvfion
    real*8 nn,np,effmn,effmp
    !  variables connected to routine EOSMAG (Potekhin)
    real*8 RHONUC, TEMP, UN_B12, UN_T6, GAMAG
    parameter(UN_B12=425.438, UN_T6=.3157746)
    
    t = t8*1.d8
    cvn  = 0.d0
    cvp  = 0.d0
    cvion= 0.d0
    cve  = 0.d0

    nn=rho*yn/MASS_N
    np=rho*yp/MASS_N
    ! Effective masses
    call eff_mass (nn,np,effmn,effmp)
    if(rho < 2.854d14) effmp=1.d0
    ! Heat capacity
      if (xh == 0.d0) then    ! CORE 
        cvn = ccvfn(t,rho,effmn,xh,yn,taun,taucru)
        cvp = ccvfp(t,rho,effmp,yp,taup)
        cvion = 0.d0
        cve = ccvfe(t,rho,ye)
      else    ! CRUST
        cvn = ccvfn(t,rho,effmn,xh,yn,taun,taucru)
        cvp = 0.d0
!       cvion = ccvfion(t,rho,xh,ah,zh)
!       cve = ccvfe(t,rho,ye)
        TEMP=1.d-6*t/UN_T6
        GAMAG=B12*UN_B12
        RHONUC=RHO*xh
        call EOSMAG(zh,ah,RHONUC,TEMP,GAMAG,cvion,cve)
      endif

      cv_cgs=cvion+cve+cvn+cvp

      if(cv_cgs <= 0.d0) stop 'cv_func (cvf.f90): Negative heat capacity'

      return
    end subroutine cv_func

    !--------------------------------------------------------
    !   Neutron heat capacity
    !--------------------------------------------------------
    real*8 function ccvfn(t,rho,effmn,xh,yn,taun,taucru)
      implicit none
      real*8 t,rho,effmn,xh,yn,taun,taucru
      real*8 nn,xnn,xv,v1,v2,Rcv_n,RRcv_B,RRcv_A_cru

      nn=rho*yn/1.66d-24
      xnn = (6.497229d-14/effmn)*nn**(1.d0/3.d0)
      ccvfn = 4.5507d+11*effmn**2*xnn*dsqrt(xnn**2+1.d0)*t 
      Rcv_n=1.d0 
      !  Apply Superfluid corrections 
      if ((taun.le.1.d0).and.(xh.eq.0.d0)) then  !only neutrons in the core
        xv=v2(taun) !3P_2 (nB)
        Rcv_n=RRcv_B(xv)        
      endif
      if ((taucru.le.1.d0).and.(xh.gt.0.d0)) then   !only neutrons in the crust
        xv=v1(taucru)    !1S_0 (nA)
        Rcv_n=RRcv_A_cru(xv) 
      endif
      ccvfn=ccvfn*Rcv_n

      return
    end function ccvfn

    !--------------------------------------------------------
    !   Proton heat capacity
    !--------------------------------------------------------
    real*8 function ccvfp(t,rho,effmp,yp,taup)
      implicit none
      real*8 t,rho,effmp,yp,taup
      real*8 np,xpp,Rcv_p,xv,v1,RRcv_A

      np=rho*yp/1.66d-24
      xpp = (6.497229d-14/effmp)*np**(1.d0/3.d0)
      ccvfp = 4.5507d+11*effmp**2*xpp*dsqrt(xpp**2+1.d0)*t 
      Rcv_p=1.d0
      !  Apply Superfluid corrections
      if (taup.le.1.d0) then    !1S_0 protons in core
        xv=v1(taup)
        Rcv_p=RRcv_A(xv)
      endif
      ccvfp=ccvfp*Rcv_p

      return
    end function ccvfp

    !--------------------------------------------------------
    !   Electron heat capacity  (old simple formula)
    !--------------------------------------------------------
    real*8 function ccvfe(t,rho,ye)
      implicit none
      real*8 t,rho,ye,ne,xe
      ne=rho*ye/1.66d-24
      xe=0.010067d0*(YE*rho)**(1.d0/3.d0)
      ccvfe=ne*2.298976d-25*t*dsqrt(xe**2+1.d0)/xe**2 !erg/K/cm3
      return
    end function ccvfe

    !--------------------------------------------------------
    !   Ion heat capacity  (old simple formula)
    !--------------------------------------------------------
    real*8 function ccvfion(t,rho,xh,ah,zh)
      implicit none
      real*8 t,rho,xh,ah,zh
      real*8 nions,tdebye,gamma,fdebye

      nions=rho*xh/(ah*1.66d-24)
      tdebye=3.48d3*dsqrt(rho)*ZH/AH
      gamma = (2.275d5*ZH**2*(rho*XH/AH)**(1.d0/3.d0))/t

      if(gamma < 1.d0) then  
        ccvfion=nions*4.143d-16/2.d0
      elseif(gamma < 150.d0) then
        ccvfion=nions*4.143d-16*fdebye(t/tdebye) !erg/k/cm3
      else
        ccvfion=nions*(4.143d-16/2.d0)*(1.d0+dlog10(gamma)/dlog10(150.d0))
      end if

      return
    end function ccvfion

    !--------------------------------------------------------
    !   Debye function
    !--------------------------------------------------------
    real*8 function fdebye(x)
      real*8 x
      if (x <= 0.15d0) then
          fdebye = 77.9273d0*x**3
      elseif (x >= 0.4d0) then
          fdebye = 1.d0 - 1.d0/(20.d0*x**2)
      else
          fdebye = 1.69798d0*x+0.0083073d0
      endif
      return
    end function fdebye

    !-----------------------------------------------------------------------
    !   Superfluid/Superconducting suppression factors of the heat capacity
    !       From Levenfish \& Yakovlev (1994a)
    !      proton SFin the core
    !-----------------------------------------------------------------------
    real*8 function RRcv_A(v)
      implicit none
      real*8 ccc,v,a,zexp1 
!        ccc=1.d0/2.2736d0
      ccc=1.d0
      a=0.4186d0+dsqrt(1.007d0**2+(0.5010d0*v)**2)
      zexp1=1.456d0-dsqrt(1.456d0**2+v**2)
      if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff      
      RRcv_A=ccc*a**(2.5d0)*dexp(zexp1)  
      return
    end function RRcv_A
    !-----------------------------------------------------------------------
    !   Superfluid/Superconducting suppression factors of the heat capacity
    !       From Levenfish \& Yakovlev (1994a)
    !      neutron SF in the crust
    !-----------------------------------------------------------------------
    real*8 function RRcv_A_cru(v) 
      implicit none
      real*8 ccc,v,a,zexp1
!        ccc=1.d0/2.3573d0
      ccc=1.d0
      a=0.4186d0+dsqrt(1.007d0**2+(0.5010d0*v)**2)
      zexp1=1.456d0-dsqrt(1.456d0**2+v**2)
      if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff      
      RRcv_A_cru=ccc*a**(2.5d0)*dexp(zexp1)  
      return
    end function RRcv_A_cru
    !-----------------------------------------------------------------------
    !   Superfluid/Superconducting suppression factors of the heat capacity
    !       From Levenfish \& Yakovlev (1994a)
    !      neutron SF in the core
    !-----------------------------------------------------------------------
    real*8 function RRcv_B(v)
      implicit none
      real*8 ccc,v,a,zexp1
!        ccc=1.d0/2.167d0
      ccc=1.d0
      a=0.6893d0+dsqrt(0.790d0**2+(0.2824d0*v)**2)
      zexp1=1.934d0-dsqrt(1.934d0**2+v**2)
      if (zexp1.lt.-1.5d2) zexp1=-1.5d2 !cutoff      
      RRcv_B=ccc*a**2*dexp(zexp1)  
      return
    end function RRcv_B

  

!> @author
!> Jose A. Pons 
!
!> @brief This module provides the boundary conditions for the 
!>        temperature evolution (Tb-Ts relation, envelope model)
!>        For the explicit methos, sfluxb is the flux and cfluxb=0.
!>        For the implicit method, it must provide the flux minus the
!>        derivative of the flux w.r.t. the temperature times 
!>        the intensity (sfluxb) and the derivative of the flux (cfluxb)
!>        It contains corrections for magnetic envelopes.
!
  subroutine fluxesb(bpdip)

  use grid, only: kmax, lmax, np, g14, benu, br, bm
  use grid, only: tem0, cfluxb, sfluxb, tss
  use constants, only: STEFAN_BOLTZMANN, UNIT_R
  use input_params, only: ienv
      
  implicit none

  ! Dipolar magnetic field at the surface
  real*8, intent(in) :: bpdip
  ! internal variables
  integer k,i
  real*8, dimension(kmax) :: x
  real*8 :: chip,chit,tb9,tFe4,tacc4,corr
  real*8 :: t04, t14, tt, tp, aa1, aa2, tmax, costh

  x(:)=tem0(:,lmax)

  do k=1,kmax
        
    tb9=0.1d0*x(k)
 
    !-----------------------------------------------------------------------
    ! Gudmundsson 1983 (non-magntized)
    !-----------------------------------------------------------------------
    if (ienv == 2 .or. ienv == 3) then
      tss(k)=1.d6*(g14**0.455*x(k)/1.288d0)**(1./1.82)

    !-----------------------------------------------------------------------
    ! Accreted non-magnetic envelope.
    !-----------------------------------------------------------------------
    else if (ienv == 1) then
      corr=0.447d0+0.075d0*dlog10(1d9*tb9)/(1d0+(6.2d0*tb9)**4)
      tacc4=(g14*(18.1d0*tb9)**2.42d0*corr+3.2d0*tb9**1.67d0*tFe4)/(1d0+3.2d0*tb9**1.67d0)
      tss(k)=1.d6*(tacc4)**(0.25)
    !-----------------------------------------------------------------------
    ! Magnetized envelopes (Potekhin, Pons & Page 2015, Appendix B).
    ! (A surface density of 10^10 g/cm^3 is assumed.)
    !-----------------------------------------------------------------------
    else 
      t04=(15.7d0*tb9**1.5d0+1.36d0*tb9)**1.5184d0
      t14=1.63d0*bpdip**0.476d0*tb9**(1.348d0/(1d0+1d-2*dsqrt(bpdip)))
      tFe4=g14*(t14+(1d0+0.15d0*dsqrt(bpdip))*t04)
      tp=1d6*tFe4**(0.25d0)

      if (bpdip.eq.0.d0) then
        tss(k) = tp
      else  
        tmax=1d6*(5.2d0*g14**0.65d0+0.093d0*dsqrt(g14*bpdip))
        tp=tp*(1d0+(tp/tmax)**4)**(-0.25d0)
        chit=1d0/(1d0+(1230d0*tb9)**3.35d0*bpdip &
     &          *dsqrt(1d0+2d0*bpdip**2)  &
     &          /(bpdip+450d0*tb9+119d0*bpdip*tb9)**4  &
     &          +6.6d-3*bpdip**2.5d0/(dsqrt(tb9)+2.58d-3*bpdip**2.5d0))
        tt=chit*tp

        i=2*k-2 
        if (bm(i,np) /= 0.) then
          costh = dabs(br(i,np)/bm(i,np))
        else
          costh = 1d0
        endif
        aa2=10d0*bpdip/(dsqrt(tb9)+0.1d0*bpdip*tb9**(-0.25d0))
        aa1=aa2*dsqrt(tb9)/3d0
        chip=(1d0+aa1+aa2)*costh**2/(1d0+aa1*costh+aa2*costh**2)
        tss(k)=tt+(tp-tt)*chip
      endif

    endif

  enddo

  if (ienv == 3) then
    ! Implicit (only for Gudmundsson), by linearizing first order derivative wrt T
    ! In this case sfluxb contained already the derivative term
    sfluxb=(STEFAN_BOLTZMANN*tss**4*UNIT_R**2*1d-40)*(1d0-4d0/1.82d0)  
    cfluxb=4d0*(STEFAN_BOLTZMANN*tss**4*UNIT_R**2*1d-40)/x/1.82d0
    ! Flux - (dF/dT)*T.
    cfluxb=cfluxb*benu(np)
  else
    ! Explicit source 
    sfluxb = STEFAN_BOLTZMANN*tss**4*UNIT_R**2*1d-40
    cfluxb=0d0
  endif

  ! Redshift corrections
  sfluxb=sfluxb*benu(np)**2        ! 10^40 erg/km^2/s.
  tss=tss*benu(np)


  end subroutine fluxesb 

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   ----------------------   MAIN block   -----------------------      *
*       This is auxiliary INTERFACE      *
*  MODIFIED AND ADAPTED FOR OUR PURPOSES IN Nov 2019 (J.A. PONS)
*     Calculations are performed in the subroutine CONDUCT19 (below),    *
*     where most of internal quantities are in the relativistic units  *
*         (\hbar = m_e = c = 1)                                        *
*     Subroutine CONDCONV converts the conductivities to CGS units.    *
* In this MAIN program we also demonstrate the convertions to SI units *
* for electrical conductivity and to cm^2/g for thermal opacity        *
* Calculations are performed in CONDEGIN for the inner  crust          *  
*   ---------------------------------------------------------------    *
      subroutine potekhin2019(T8,RHO,B12,CMI,Zion,XH,Zimp,ye,
     & yn,taucru,CKAPPA,CKAPPAT,CKAPPAH,SIGMA,QJ,QJT,QJH,CRKAPI,Ksfph)

      implicit none
      real*8 T8,RHO, B12,CMI,CMI1,Zion,XH,Zimp,T6,B!,RHOlg,Tlg
      real*8 SIGMA,CKAPPA,QJ,SIGMAT,CKAPPAT,QJT,SIGMAH,CKAPPAH,QJH
      real*8 AUM, AUD, TEMP, DENSI, UNISIG, UNIKAP
      real*8 RSIGMA,RTSIGMA,RHSIGMA,RKAPPA,RTKAPPA,RHKAPPA, RKAPI
      real*8 CRKAPi, ksfph
      real*8 ye,yn,taucru
* NOTATIONS:
*        AUM - atomic mass unit divided by the electron mass
*        AUD - relativistic unit of density in g/cm^3
* conv.th.conductivity [erg cm^{-1} s^{-1} K^{-1}] into opacity [cm^2/g]

      T6= 1.d2*T8
      CMI1=CMI/XH
      QJ=0.
      QJH=0.
      QJT=0.

       Rkapi=0.d0
       Ksfph=0.d0

      UNISIG=7.763d20
      UNIKAP=2.778d15
      TEMP=T6/5930. ! Temperature in mc^2
      AUM=1822.9
      AUD=15819.4
      DENSI=XH*RHO/(AUD*AUM*CMI) ! number density of ions

      if (xh.eq.1) then

* DESIGNED for outer crust (xh=1)
* 'conduct08.f' of Potekhin: CONDCONV + CONDUCT
* ("Magnetic conductivity tensors" in the website; 'condeg08.f' does not take into account thermal averaging and finite size nuclei)
* Thermal averaging (important in outer, partially deg.): YES
* Finite size nuclei (important in inner crust): YES, added by Pons
* Quantizing effect (important in outer crust, large field): YES, see MAGNET=0/1 in CONDUCT

      B=1.d12*B12
      call CONDCONV(T6,RHO,B,Zion,CMI,Zimp, ! input
     *  SIGMA,RKAPPA,QJ,SIGMAT,RTKAPPA,QJT,SIGMAH,RHKAPPA,QJH) ! output

* If you want to allow for ion thermal conduction in the approximation
*  of Chugunov & Haensel (2007), then uncomment the next line:
c-------- Phonon thermal conductivity: outer crust -----------------------------------------------
      call CONDI(TEMP,DENSI,Zion,CMI,RKAPi)
      RKAPPA=RKAPPA+RKAPi
      RTKAPPA=RTKAPPA+RKAPi
*   -------  CONVERSION TO ORDINARY PHYSICAL (CGSE) UNITS:   --------
      CKAPPA=RKAPPA*UNIKAP ! KAPPA in erg/(K cm s)
c------- added for visualization purposes -----------------------------
      CRKAPi=RKAPi*UNIKAP
      CKAPPAT=RTKAPPA*UNIKAP
      CKAPPAH=RHKAPPA*UNIKAP


      else

* DESIGNED for inner crust (xh<1)
* 'condegin08.f' of Potekhin: CONDEGIN 
* ("For inner crust" in the website; 'condesc08.f' contains form factor, but for a fixed composition.
* Differences come from the composition itself (Eos) rather than form factors.)
* Thermal averaging (important in outer, partially deg.): NO
* Finite size nuclei (important in inner crust): YES
* Quantizing effect (important in outer crust, large field): YES, see conditions in CONDEGIN

c      AUM=1822.9
c      AUD=15819.4
c      DENSI=XH*RHO/(AUD*AUM*CMI) ! number density of ions
      B=B12/44.14 ! B is the magnetic field in relativistic units
c      TEMP=T6/5930. ! Temperature in mc^2
*  Call for the central subroutine which calculates the transport
*  coefficients SIGMA,KAPPA,Q (in Relativistic units):
      call CONDEGIN(TEMP,DENSI,B,Zion,CMI,CMI1,Zimp,
     & RSIGMA,RTSIGMA,RHSIGMA,RKAPPA,RTKAPPA,RHKAPPA)

c      Tlg=6.d0+dlog10(T6)
c      RHOlg=dlog10(RHO)
c      call CONDEGsc(Tlg,RHOlg,B12,Zimp,
c     & RSIGMA,RTSIGMA,RHSIGMA,RKAPPA,RTKAPPA,RHKAPPA)

c----- Allow for ion/lattice phonon thermal conduction ------------------
      call CONDIN(TEMP,DENSI,Zion,CMI,CMI1,RKAPi)
c      call CONDIsc(TEMP,DENSI,Zion,CMI,CMI1,RKAPi)
c----- Allow for the SF phonon thermal conductivity (cgs unit) ----------
      call cphonon(T8,Zion,CMI,ye,RHO,xh,yn,taucru,B12,ksfph)
c------------------------------------------------------------------------ 
      RKAPPA=RKAPPA+RKAPi
      RTKAPPA=RTKAPPA+RKAPi
*   -------  CONVERSION TO ORDINARY PHYSICAL (CGSE) UNITS:   --------  *
c------- added for visualization purposes -------------------------------
      CRKAPi=RKAPi*UNIKAP
c------------------------------------------------------------------------
      SIGMA=RSIGMA*UNISIG ! SIGMA in s^{-1}
      CKAPPA=RKAPPA*UNIKAP ! KAPPA in erg/(K cm s)
      CKappa=Ckappa+Ksfph  
      SIGMAT=RTSIGMA*UNISIG
      CKAPPAT=RTKAPPA*UNIKAP
      CKAPPAT=CKAPPAT+Ksfph
      SIGMAH=RHSIGMA*UNISIG
      CKAPPAH=RHKAPPA*UNIKAP

      endif

      if (qj .ne. qj .or. qjh.ne. qjh .or. qjt.ne. qjt) then
        print*,'potekhin2012.f error:',t6,rho,qj,qjh,qjt,xh
        stop
      endif

      return
      end


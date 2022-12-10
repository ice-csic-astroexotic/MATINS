* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   CALCULATION OF ELECTRON TRANSPORT COEFFICIENTS IN MAGNETIC FIELDS  *
*           for the case of very strong electron degeneracy            *
*      (thermal averaging is completely neglected in this version).    *
*   For theoretical background and references see CONDEGIN:            *
*           http://www.ioffe.rssi.ru/astro/conduct/                    *
* Difference from generic CONDEGEN - allowance for nuclear form factor *
*      Remarks and suggestions are welcome. Please send them to        *
*       Alexander Potekhin <palex@astro.ioffe.ru>                      *
* Last updates: 18.02.2008 - ion thermal conduction included           *
*              06.10.2008 - ion thermal conduction CORRECTED           *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
* * * * * * *   Block cond_crust_electrons  * * * * * * * * * * * * * **
*  This subroutine calculates the electron electric conductivity tensor*
*                   in strongly degenerate matter                      *
*   --------------------------------------------------  Version 12.11.07
! PREVIOUSLY CALLED CONDEGIN
      subroutine cond_inner_crust_electrons(TEMP,DENSI,B,
     & Zion,CMI,CMI1,Zimp,
     & RSIGMA,RTSIGMA, RHSIGMA, RKAPPA, RTKAPPA, RHKAPPA)
* Input: TEMP - temperature, DENSI - number density of ions, B - magn.f.
*        Zion and CMI - ion charge and mass numbers,
*        CMI1 - number of nucleons per nucleus, Zimp -impurity parameter
*          (eff.Z) - that is, Z_{imp}^2 = < n_j (Z-Z_j)^2 > / n,
*         where Z_j, n_j are partial charges and densities of impurities
* Output: RSIGMA, RTSIGMA, RHSIGMA - 
*         longitudinal, transverse and off-diagonal conductivities
*         (all in the relativistic units: \hbar = m_e = c = 1)
      implicit double precision (A-H), double precision (O-Z)
      save
      data PI/3.14159265d0/
      data BOHR/137.036/

      iquant = 0
*   -------------------   RESTRICTIONS:   --------------------------   *
      if (TEMP.le.0..or.DENSI.le.0..or.B.lt.0..or.Zion.le.0..or.
     &        CMI.le.0.) stop '<ERROR>:Non-positive input parameter'
      if (CMI1.lt.CMI)   stop '<ERROR>: Incorrect CMI1'
      if (Zion.lt..5)    stop '<ERROR>: Too small ion charge'
      if (CMI.lt.1.)     stop '<ERROR>: Too small ion mass'
      if (DENSI.gt.1.d6) stop '<ERROR>: Too high density'
*   -----------------   PLASMA PARAMETERS   ------------------------   *
      DENS=DENSI*Zion ! DENS - number density of electrons
      SPHERION=(.75/PI/DENSI)**.3333333 ! Ion sphere radius
      GAMMA=Zion**2/BOHR/TEMP/SPHERION ! Ion coupling parameter
      XSR=(3.*PI**2*DENS)**.3333333 ! special relativity parameter 
*   XSR equals the non-magnetic Fermi momentum
      EF0=dsqrt(1.+XSR**2) ! non-mag.Fermi energy
      Q2e0=4./PI/BOHR*XSR*EF0 ! non-magn.e-screening at strong degen.
      CST=XSR**3/.75/PI*Zion/BOHR**2 ! =4\pi n_i(Ze^2)^2
      if (CMI.eq.CMI1) then ! outer envelope
         xnuc=.00155*(CMI/Zion)**.33333*XSR ! i.e. r_nuc=1.15 A^{1/3} fm
      else ! inner envelope
         xnuc=.00247*XSR ! i.e. r_nuc=1.83 Z^{1/3} fm (Itoh&Kohyama'83)
      endif ! xnuc=r_nuc/a_i - nucleus size parameter
*   -------------------  Chemical potential   ----------------------   *
*   ---------------- QUANTIZING EFFECTS ----------------------------   *
      if (XSR**2.gt.4.d2*B .or. iquant .le. 1) then ! non-quantizing case
         PCL=XSR
         Q2e=Q2e0
      else ! quantizing magn.field
         PM0=XSR**3/1.5/B ! p_F in strongly quantizing field
        if (PM0**2.le.2.*B) then ! strongly quantizing case
           PCL=PM0
           Q2e=B/PI/BOHR/PCL ! e-screening
        else ! weakly quantizing case - find p_F by iteration
           Pmax=PM0
           Pmin=XSR/2.
   22      PCL=(Pmax+Pmin)/2.
           NL=PCL**2/2./B
           SN=PCL
           SM=1./PCL
          do N=1,NL
             PN=dsqrt(PCL**2-2.*B*N) ! =p_n
             SN=SN+2.*PN
             SM=SM+2./PN
          enddo
           D=SN*B/2./PI**2 ! estimate n_e
          if (D.lt.DENS) then ! increase PCL
             Pmin=PCL
          else ! decrease PCL
             Pmax=PCL
          endif
          if (dabs(D-DENS).gt.1.d-4*DENS) goto 22 ! next iteration
           Q2e=B/PI/BOHR*SM
        endif
      endif
*   -------------------   Relaxation times    ----------------------   *
      call COULIN(PCL,XSR,GAMMA,B,Zion,CMI,Q2e,xnuc,iquant,
     *   CLeff,CLlong,CLtran,SN,THtoEL)
      E=dsqrt(1.+PCL**2) ! magn.Fermi energy
      TAU=PCL**3/E/4./PI/DENSI/(Zion/BOHR)**2/CLlong/SN
      GYROM=B/E ! gyrofrequency
      TAUt0=PCL**3/E/CST/CLtran*SN
      if (Zimp.gt.0.) call COUL99I(PCL,XSR,GAMMA,B,Q2e,xnuc, ! incl.impurity
     *   CLeffI,CLlongI,CLtranI,SN,iquant)
      TAUlong=TAU*CLlong*Zion**2/(CLlong*Zion**2+CLlongI*Zimp**2)
      TAUt=TAUt0*CLtran*Zion**2/(CLtran*Zion**2+CLtranI*Zimp**2)
      TAUtran=TAUt/(1.+(TAUt*GYROM)**2)
      TAUhall=TAUt*GYROM*TAUtran
* Modification of thermal conductivity: inclusion of THtoEL (21.11.99)
      TAUlongT=TAU*CLlong*Zion**2/
     /  (CLlong*Zion**2*THtoEL+CLlongI*Zimp**2)
      TAUtT=TAUt0*CLtran*Zion**2/(CLtran*Zion**2*THtoEL+CLtranI*Zimp**2)
      TAUtranT=TAUtT/(1.+(TAUt*GYROM)**2)
      TAUhallT=TAUtT*GYROM*TAUtranT
*   ----------------------------------------------------------------   *
*   Longitudinal transport coefficients:
      C=SN*PCL**3/3./PI**2/E/BOHR ! common factor

      RSIGMA=C*TAUlong
      RTSIGMA=C*TAUtran
      RHSIGMA=C*TAUhall
* Find thermal conductivity from the Wiedemann-Franz law:
      CTH=C*PI**2*TEMP/3.*BOHR
      if (B.eq.0.) then ! corrected 12.11.07
         call TAUEESY(XSR,TEMP,TAUEE) ! eff.e-e relax.time
         EECOR=TAUEE/(TAUlongT+TAUEE)
      else
         EECOR=1.
      endif
      RKAPPA=CTH*TAUlongT*EECOR
      RTKAPPA=CTH*TAUtranT*EECOR
      RHKAPPA=CTH*TAUhallT*EECOR
      return
      end

*   ----------------------------------------------------
      subroutine cond_inner_crust_ions(TEMP,DENSI,Zion,CMI,CMI1,RKAPi)
* ion thermal conductivities in the inner crust, Chugunov & Haensel'07
* Input: TEMP - temperature, DENSI - number density of ions
*        Zion and CMI - ion charge and mass numbers
*        CMI1 - number of nucleons per nucleus
* Output: RKAPi - ion thermal conduction
* All quantities are in the relativistic units: \hbar = m_e = c = 1.
*                                                       Version 06.10.08
      implicit double precision (A-H), double precision (O-Z)
      save
      data Uminus1/2.8/,Uminus2/13./,AUM/1822.9/,BOHR/137.036/
* Dimensional quantities are in the relativistic units (m_e=\hbar=c=1)
*        Uminus1,Uminus2 - dimensionless frequency moments of phonons
*        AUM - atomic mass unit divided by the electron mass
*        BOHR - radius of the first Bohr orbit in the rel.units
      parameter(PI=3.14159265d0,EPS=1.d-8,TINY=1.d-99)
*   ----------------------   Preliminaries   -----------------------   *
      DENS=DENSI*Zion ! number density of electrons
      XSR=(3.d0*PI**2*DENS)**.3333333 ! x_r - density parameter
      PCL=XSR ! classical Fermi momentum; equality due to B=0.
      VCL=PCL/dsqrt(1.d0+PCL**2)
      SPHERION=(.75/PI/DENSI)**.3333333 ! Ion sphere radius
      QBZ=PCL*(2./Zion)**.3333333 ! q_{BZ}
      GAMI=Zion**2/(BOHR*SPHERION*TEMP)
      TRP=Zion/GAMI*dsqrt(CMI*AUM*SPHERION/3.d0/BOHR) ! =T/T_pi
      OMPI=TEMP/TRP ! ion plasma frequency
      if (GAMI.lt.TINY**.2) then
        write(*,'(a)') '<ERROR>:GAMI too low, kappa too large'
        stop
      endif
*   ---------------------- reduced ion heat capacity:
      call HLfit8(1.d0/TRP,F,U,CV,S,1)
*   ---------------------- ion-electron:
      A0=1.683*sqrt(PCL/CMI/Zion) ! zero-vibr.param.(Baiko&Yakovlev95)
      WDW=A0*(.5*Uminus1*exp(-9.1*TRP)+Uminus2*TRP) ! DW factor (B&Y'95)
      if (CMI.eq.CMI1) then ! outer envelope
         xnuc=.00155*(CMI/Zion)**.33333*XSR ! i.e. r_nuc=1.15 A^{1/3} fm
      else ! inner envelope
         xnuc=.00247*XSR ! i.e. r_nuc=1.83 Z^{1/3} fm (Itoh&Kohyama'83)
      endif ! xnuc=r_nuc/a_i - nucleus size parameter
         
      W=WDW+43.*xnuc**2
      F=.014+.03/(1.d0+dexp(.2/TRP)) ! CH'07, Eq.(45)
      Y=QBZ/(2.*PCL)
      WY2=W*Y**2
      if (W.gt.EPS) then
         EXPW=dexp(-W)
         EXPWY2=dexp(-WY2)
         CL=(EXPINT(WY2,0)*EXPWY2-EXPINT(W,0)*EXPW-
     -     VCL*(EXPWY2-EXPW)/W)/2.
      else ! small-W asymptote
         CL=dlog(1.d0/Y)-VCL*(1.d0-Y**2)/2.
      endif
      FLEN=8.32d5*SPHERION/(1.d0+PCL**2)/CL*F/Zion*
     *  dsqrt((XSR/1.00884)**3*CMI/Zion) ! Eq.(43)
      CS=OMPI/(3.*QBZ) ! sound speed
      RKAPie=CV*DENSI*CS*FLEN/3. ! Eq.(42) of CH'07
*   ---------------------- ion-ion:
      if (TRP.lt.1.d-4) then
         write(*,'(a)') '<ERROR>:T too low, kappa too large'
         RKAPi=RKAPie
         stop
      else
         RKAP0=OMPI*DENSI*SPHERION**2
         RKAPii=RKAP0*dsqrt(1./GAMI**5/
     /     dlog(2.d0+.57735/dsqrt(GAMI)**3)**2+
     +     0.16+(GAMI/77.)**2*dexp(.666667/TRP)) ! Eq.(26) of CH'07
*   ---------------------- total ion:
         RKAPi=1./(1./RKAPii+1./RKAPie)
      endif
      return
      end
      
*  ================   EFFECTIVE COULOMB LOGARITHM  ===================  *
      subroutine COULIN(PCL,XSR,GAMMA,B,Zion,CMI,Q2e,xnuc,iquant,
     *   CLeff,CLlong,CLtran,SN,THtoEL)
*                                                        Version 24.02.00
*   Input: PCL - non-magnetic electron momentum \equiv \sqrt(E^2-1),
*          Q2e - squared electron screening wavenumber (IN RELATIV.UNITS)
*          XSR = p_F/mc - relativity (density) parameter,
*          GAMMA - Coulomb coupling parameter of ions,
*          B - magnetic field,
*          Zion - mean charge of the ion,
*          CMI - mean atomic weight,
*          xnuc - ion radius divided by Wigner-Seitz cell radius
*   Output: CLlong, CLtran - eff.Coulomb log.,
*           SN = N_e(E)/N_0(E) = (3/2)(eB\hbar/c)\sum_{ns} p_n/p_0^3
      implicit double precision (A-H), double precision (O-Z)
      integer iquant
      save
      data Uminus1/2.78/,Uminus2/12.973/,AUM/1822.9/,BOHR/137.036/
* Dimensional quantities are in the relativistic units (m_e=\hbar=c=1)
*        Uminus1,Uminus2 - dimensionless frequency moments of phonons
*        AUM - atomic mass unit divided by the electron mass
*        BOHR - radius of the first Bohr orbit in the rel.units
      data PI/3.14159265/
*   ----------------------   Preliminaries   ------------------------   *
      DENS=XSR**3/3./PI**2 ! number density of electrons
      DENSI=DENS/Zion ! number density of ions (rel.)
      SPHERION=(.75/PI/DENSI)**.3333333 ! Ion sphere radius
      Q2icl=3.*GAMMA/SPHERION**2 ! squared Debye screening momentum
      ECL=sqrt(1.+PCL**2) ! Energy
      VCL=PCL/ECL ! Velocity
      PM2=(2.*PCL)**2 ! squared max.momentum transfer
      TRP=Zion/GAMMA*sqrt(CMI*AUM*SPHERION/3./BOHR) ! =T/T_p
      BORNCOR=VCL*Zion*PI/BOHR ! first non-Born correction
*   ---------------------   Non-magnetic fit   ----------------------   *
      C=(1.+.06*GAMMA)*dexp(-dsqrt(GAMMA))
      Q2s=(Q2icl*C+Q2e)*dexp(-BORNCOR) ! eff.scr.wavenumber in rel.un.
      XS=Q2s/PM2 ! eff.screening param.
      R2W=Uminus2/Q2icl*(1.+.3333*BORNCOR)
      XW=R2W*PM2 ! eff. Debye-Waller param.
** Modification WITH FINITE SIZES OF NUCLEI; xnuc=r_{nuc}/a_i
      XW1=14.7327*xnuc**2 ! =4(9\pi/4)^{2/3} x_{nucl}^2 =coeff.at q^2
      XW1=XW1*(1.+.3333*BORNCOR)*(1.+Zion/13.*dsqrt(xnuc))
      CL=COULAN2(XS,XW,VCL,XW1)
      A0=1.683*sqrt(PCL/CMI/Zion) ! zero-vibr.param.(Baiko&Yakovlev95)
      VIBRCOR=exp(-A0/4.*Uminus1*exp(-9.1*TRP)) ! corr.for zero-vibr.
      T0=.19/Zion**.16667 ! critical T/T_p parameter
      G0=TRP/sqrt(TRP**2+T0**2)*(1.+(Zion/125.)**2) ! mod.10.01.99
      GW=G0*VIBRCOR
      CLeff=CL*GW ! 1st FIT (for non-magnetic electrical conductivity)
      G2=TRP/sqrt(.0081+TRP**2)**3
      THtoEL=1.+G2/G0*(1.+BORNCOR*VCL**3)*.0105*(1.-1./Zion)*
     *  (1.d0+xnuc**2*dsqrt(2.d0*Zion))
      TRU=TRP*3.*VCL*BOHR/Zion**.3333333 ! T/Tu
      if (TRU.lt.20.) then ! correction for dying-out umklapp processes
         CLhigh=CLeff
         EU=dexp(-1.d0/TRU)
         CLlowK=50.*sqrt(XSR/CMI)/Zion*TRP**3 ! low-T lim.for kappa
         CLlowS=CLlowK/VCL/BOHR/.75*TRP**2 ! low-T lim.for sigma
         CLeff=CLhigh*EU+CLlowS*(1.-EU)
         THtoEL=(CLhigh*THtoEL*EU+CLlowK*(1.-EU))/CLeff
      endif
      if (PCL**2.gt.4.d2*B .or. iquant .eq. 0) then ! Non-magnetic case
         CLlong=CLeff
         CLtran=CLeff
         SN=1.d0
         goto 50
      endif
*   -----------------------   Magnetic fit   ------------------------   *
      ENU=PCL**2/2.d0/B
      NL=ENU
      SN=0.
      do N=0,NL
         PB=dsqrt(ENU-N) ! =p_n/sqrt(2b)
         SN=SN+PB
        if (N.ne.0) SN=SN+PB
      enddo
      SN=SN*1.5d0*B*dsqrt(2.d0*B)/PCL**3
      if (ENU.le.1.d0) then ! Exact calculation     
         Xis=Q2s/2./B ! Screening parameter, scaled magnetically
         ZETA=R2W*2.*B ! magn.scaled exponential coefficient
         Xi=2.*PCL**2/B
         Xsum=Xi+Xis
         Q2M=(EXPINT(Xsum,1)-
     -     dexp(-ZETA*Xi)*EXPINT((1.+ZETA)*Xsum,1))/Xsum
         CLlong=(PCL*VCL/B)**2*Q2M/1.5*GW
         QtranM=(1.+Xsum)*EXPINT(Xsum,0)-1.-dexp(-ZETA*Xi)*
     *      ((1.+(1.+ZETA)*Xsum)*EXPINT((1.+ZETA)*Xsum,0)-1.)
         QtranP=(1.+Xis)*EXPINT(Xis,0)-1.-
     -   ((1.+(1.+ZETA)*Xis)*EXPINT((1.+ZETA)*Xis,0)-1.)
         Q=(ECL**2*QtranP+QtranM)*B/PCL**2 ! Q(E,b)
         CLtran=.375*Q/ECL**2*GW
      else
*   Preliminaries:
         DNU=ENU-NL
         XS1=(dsqrt(XS)+1./(2.+XW/2.))**2
         PN=dsqrt(2.*B*DNU)
         SQB=dsqrt(B)
         X=dmax1(PN/SQB,1.d-10)
*   Longitudinal:
        if (XW.lt..01) then
           EXW=1.
        elseif (XW.gt.50.) then
           EXW=1./XW
        else
           EXW=(1.d0-dexp(-XW))/XW
        endif
         A1=(30.-15.*EXW-(15.-6.*EXW)*VCL**2)/
     /     (30.-10.*EXW-(20.-5.*EXW)*VCL**2)
         Q1=.25*VCL**2/(1.-.667*VCL**2)
         DLT=SQB/PCL*(A1/X-sqrt(X)*(1.5-.5*EXW+Q1)+
     +     (1.-EXW+.75*VCL**2)/(1.+VCL**2)*(X-sqrt(X))/NL)
         Y1=1./(1.+DLT)
         CL0=dlog(1.d0+1.d0/XS1)
         P2=CL0*(.07+.2*EXW)
         Y2=1.5*CL0*(X**3-X/3.)/(NL+.75/(1.+2.*B)**2*X**2)+P2*X
         PY=1.+.06*CL0**2/NL**2
         DT=dsqrt(PY*Y1**2+Y2**2) ! ratio of relax.times
         CLlong=CLeff/DT
*   Transverse:
         DB=1./(1.+.5/B)
         CL1=XS1*CL0
         P1=.8*(1.+CL1)+.2*CL0
         P2=1.42-.1*DB+sqrt(CL1)/3.
         P3=(.68-.13*DB)*CL1**.165
         P4=(.52-.1*DB)*sqrt(sqrt(CL1))
         DLT=SQB/PCL*(P1/X**2*SQB/PCL+P3*alog(NL+0.)/X-
     -     (P2+P4*alog(NL+0.))*sqrt(X))
         CLtran=CLeff*(1.+DLT)
      endif
   50 return
      end

      function COULAN2(XS,XW0,V,XW1)
*  ------    Analytic expression for Coulomb logarithm - Version 23.05.00
*   XS=(q_s/2p)^2, where p - momentum, q_s - eff.scr.momentum
*   XW0=u_{-2} (2p/\hbar q_D)^2, where u_{-2}=13, q_D^2=3\Gamma/a_i^2
*   V=p/(mc)
*   XW1=s1*(2p/\hbar)^2, s1 \approx r_{nuc}^2
      implicit double precision (A-H), double precision (O-Z)
      save
      data EPS/1.d-2/,EPS1/1.D-3/,EULER/0.5772156649 d0/
      if(XS.lt.0..or.XW0.lt.0..or.V.lt.0..or.XW1.lt.0.)stop'COULAN2'
      do I=0,1
        if (I.eq.0) then
           XW=XW0+XW1
           B=XS*XW
        else ! to do the 2nd term
          XW=XW1
          B=XS*XW
        endif
        if (I.eq.0.or.KEY.eq.2) then ! 23.05.00: for KEY=2 re-check
* Check applicability of asymptotes:
          if (XW.lt.EPS) then
             KEY=1
             goto 50
          endif
          if (XW.gt.1./EPS.and.B.gt.1./EPS) then
             KEY=2
          elseif (XS.lt.EPS1.and.B.lt.EPS1/(1.+XW)) then
             KEY=3
          else
             KEY=4
          endif
        endif
   50   continue
         EA=dexp(-XW)
         E1=1.-EA

        if (KEY.ne.1) E2=(XW-E1)/XW
        if (KEY.eq.1) then
           CL0=dlog((XS+1.)/XS)
           CL1=.5*XW*(2.-1./(XS+1.)-2.*XS*CL0)
           CL2=.5*XW*(1.5-3.*XS-1./(XS+1.)+3.*XS**2*CL0)
        elseif (KEY.eq.2) then
           CL0=dlog(1.d0+1./XS)
           CL1=(CL0-1.d0/(1.+XS))/2.
           CL2=(2.*XS+1.)/(2.*XS+2.)-XS*CL0
        elseif (KEY.eq.3) then
           CL1=.5*(EA*EXPINT(XW,0)+dlog(XW)+EULER)
           CL2=.5*E2
        elseif (KEY.eq.4) then
           CL0=dlog((XS+1.)/XS)
           EL=EXPINT(B,0)-EXPINT(B+XW,0)*EA
           CL1=.5*(CL0+XS/(XS+1.)*E1-(1.+B)*EL)
           CL2=.5*(E2-XS*XS/(1.+XS)*E1-2.*XS*CL0+XS*(2.+B)*EL)
        else
           stop'<ERROR>[COULAN2 cond_inner_crust_ei.f]: invalid KEY'
        endif
        if (I.eq.0) then ! 1st term calculated
           COULAN2=CL1-V**2*CL2
          if (XW1.lt.EPS1) return ! don't calculate the 2nd term
        else ! 2nd term calculated
           COULAN2=COULAN2-(CL1-V**2*CL2)
        endif
      enddo
      return
      end

      subroutine COUL99I(PCL,XSR,GAMMA,B,Q2e,xnuc,
     *   CLeff,CLlong,CLtran,SN,iquant) ! IMPURITY
*                                                       Version 18.11.99
*  This is a simplified version of COUL01 for the el.-impurity scattering
*   Input: XSR = p_F/mc - relativity (density) parameter,
*          PCL - non-magnetic electron momentum \equiv \sqrt(E^2-1),
*          GAMMA - Coulomb coupling parameter of ions,
*          B - magnetic field,
*          Q2e - squared electron screening wavenumber
*    (ALL IN THE RELATIVISTIC UNITS) 
*   Output: CLlong, CLtran - eff.Coulomb log.,
*           SN = N_e(E)/N_0(E) = (3/2)(eB\hbar/c)\sum_{ns} p_n/p_0^3
      implicit double precision (A-H), double precision (O-Z)
      integer iquant
      save

* Dimensional quantities are in the relativistic units (m_e=\hbar=c=1)
*        BOHR - radius of the first Bohr orbit in the rel.units
      data PI/3.14159265/,XW/1.d99/ ! XW=infinity
*   ----------------------   Preliminaries   -----------------------   *
      DENS=XSR**3/3./PI**2 ! number density of electrons
      ECL=sqrt(1.+PCL**2) ! Energy
      VCL=PCL/ECL ! Velocity
      PM2=(2.*PCL)**2 ! squared max.momentum transfer
*   ---------------------   Non-magnetic fit   ---------------------   *
      C=(1.+.06*GAMMA)*dexp(-dsqrt(GAMMA))
      Q2s=Q2e ! eff.scr.wavenumber in rel.un.
      XS=Q2s/PM2 ! eff.screening param.
** Modification WITH FINITE SIZES OF NUCLEI; xnuc=r_{nuc}/a_i
      XW1=14.7327*xnuc**2 ! =4(9\pi/4)^{2/3} x_{nucl}^2 =coeff.at q^2
      XW1=XW1*(1.+.3333*BORNCOR)*(1.+Zion/13.*dsqrt(xnuc))
      CLeff=COULAN2(XS,XW,VCL,XW1)
c      CLeff=COULAN(XS,XW,VCL) ! 1st FIT (for non-magn.el.conductivity)


      if (PCL**2.ge.4.d2*B .or. iquant.eq.0) then ! Non-magnetic case
         CLlong=CLeff
         CLtran=CLeff
         SN=1.d0
         goto 50
      endif
*   -----------------------   Magnetic fit   -----------------------   *
      ENU=PCL**2/2.d0/B
      NL=ENU
      SN=0.
      do N=0,NL
         PB=dsqrt(ENU-N) ! =p_n/sqrt(2b)
         SN=SN+PB
        if (N.ne.0) SN=SN+PB
      enddo
      SN=SN*1.5d0*B*dsqrt(2.d0*B)/PCL**3


      if (ENU.le.1.d0) then ! Exact calculation     
         Xis=Q2s/2./B ! Screening parameter, scaled magnetically
         Xi=2.*PCL**2/B
         Xsum=Xi+Xis
         Q2M=EXPINT(Xsum,1)
         CLlong=(PCL*VCL/B)**2*Q2M/1.5
         QtranM=(1.+Xsum)*EXPINT(Xsum,0)-1.
         QtranP=(1.+Xis)*EXPINT(Xis,0)-1.
         Q=(ECL**2*QtranP+QtranM)*B/PCL**2 ! Q(E,b)
         CLtran=.375*Q/ECL**2

      else
*   Preliminaries:
         DNU=ENU-NL
         XS1=(dsqrt(XS)+1./(2.+XW/2.))**2
         PN=dsqrt(2.*B*DNU)
         SQB=dsqrt(B)
         X=dmax1(PN/SQB,1.d-10)
*   Longitudinal:
        if (XW.lt..01) then
           EXW=1.
        elseif (XW.gt.50.) then
           EXW=1./XW
        else
           EXW=(1.d0-dexp(-XW))/XW
        endif
         A1=(30.-15.*EXW-(15.-6.*EXW)*VCL**2)/
     /     (30.-10.*EXW-(20.-5.*EXW)*VCL**2)
         Q1=.25*VCL**2/(1.-.667*VCL**2)
         DLT=SQB/PCL*(A1/X-sqrt(X)*(1.5-.5*EXW+Q1)+
     +     (1.-EXW+.75*VCL**2)/(1.+VCL**2)*(X-sqrt(X))/NL)
         Y1=1./(1.+DLT)
         CL0=dlog(1.d0+1.d0/XS1)
         P2=CL0*(.07+.2*EXW)
         Y2=1.5*CL0*(X**3-X/3.)/(NL+.75/(1.+2.*B)**2*X**2)+P2*X
         PY=1.+.06*CL0**2/NL**2
         DT=dsqrt(PY*Y1**2+Y2**2) ! ratio of relax.times
         CLlong=CLeff/DT
*   Transverse:
         DB=1./(1.+.5/B)
         CL1=XS1*CL0
         P1=.8*(1.+CL1)+.2*CL0
         P2=1.42-.1*DB+sqrt(CL1)/3.
         P3=(.68-.13*DB)*CL1**.165
         P4=(.52-.1*DB)*sqrt(sqrt(CL1))
         DLT=SQB/PCL*(P1/X**2*SQB/PCL+P3*alog(NL+0.)/X-
     -     (P2+P4*alog(NL+0.))*sqrt(X))
         CLtran=CLeff*(1.+DLT)
      endif
   50 return
      end

      function COULAN(XS,XW,V)
*  ------   Analytic expression for Coulomb logarithm - Version 18.11.99
*   XS=(q_s/2p)^2, where p - momentum, q_s - eff.scr.momentum
*   XW=u_{-2} (2p/\hbar q_D), 
*         where u_{-2}=13, q_D^2=3\Gamma/a_i^2
*   V=p/(mc)
      implicit double precision (A-H), double precision (O-Z)
      save
      data EPS/1.d-2/,EULER/0.5772156649 d0/
      if (XS.lt.0..or.XW.lt.0..or.V.lt.0.) stop'COULAN:range'
      EA=dexp(-XW)
      E1=1.-EA
      B=XS*XW

      if (XW.gt.EPS) then
         E2=(XW-E1)/XW
      else
         E2=XW/2.-XW**2/6.
         CL0=dlog((XS+1.)/XS)
         CL1=.5*XW*(2.-1./(XS+1.)-2.*XS*CL0)
         CL2=.5*XW*(1.5-3.*XS-1./(XS+1.)+3.*XS**2*CL0)
         goto 50
      endif
      if (XW.gt.1./EPS.and.B.gt.1./EPS) then
         CL0=dlog(1.d0+1./XS)
         CL1=(CL0-1.d0/(1.+XS))/2.
         CL2=(2.*XS+1.)/(2.*XS+2.)-XS*CL0
      elseif (XS.lt.EPS.and.B.lt.EPS) then
         CL1=.5*(EA*EXPINT(XW,0)+dlog(XW)+EULER)
         CL2=.5*E2
      else
         CL0=dlog((XS+1.)/XS)
         EL=EXPINT(B,0)-EXPINT(B+XW,0)*EA
         CL1=.5*(CL0+XS/(XS+1.)*E1-(1.+B)*EL)
         CL2=.5*(E2-XS*XS/(1.+XS)*E1-2.*XS*CL0+XS*(2.+B)*EL)

      endif

   50 continue
      COULAN=CL1-V**2*CL2
      return
      end


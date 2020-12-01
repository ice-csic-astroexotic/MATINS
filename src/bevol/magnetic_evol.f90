!-------------------------------------------------------------------------------
! Magneto Thermal 2D
!-------------------------------------------------------------------------------
! Module: Magnetic Evolution
!
!> @author
!> Daniele Viganò
!> Clara Dehman
!


! TBD: put compute_dB and compute_E instead of the pol and tor ones




!> @brief Magnetic field time advance with finite-different methods
!>        It contains:
!>        subroutine magnetic_evol (called in main)
!>          It calls one of the following subroutines:
!>              euler_alternate
!>              euler
!>              AB4_alternate
!>              AB4
!>              RK4 (which calls subroutine RK4_substep)
!>              RK4_vecpot (which calls subroutine RK4_vecpot_substep)
!>        subroutine curl_phi
!>        subroutine curl_pol
!>        subroutine add_hyperresistivity_bphi
!>        subroutine compute_Epol
!>        subroutine compute_Etor
!>        subroutine compute_Etor_upwind
!>        subroutine compute_dBpol
!>        subroutine compute_dBtor
!>        subroutine compute_dBtor_burgers
!>        subroutine magnetic_bc
!>        subroutine magnetic_bc_vecpot
!>        subroutine potential_magnetic_legendre
!>        subroutine magnetic_bc_flat
!>        subroutine compute_joule
!
! Notes:
!
!   Formalism and units used in all the following subroutines
!
!     Physical current (G):
!       J = c/(4*pi*e^nu) curl [e^nu B]
!       Our numerical J is e^nu*4*pi/c times physical J:
!       Jnum = curl [e^nu B]
!
!     Physical E-field (G):
!       E = 4 pi eta /c^2 J + 1/4pi e n_e B x J
!         = eta e^(-nu) curl [e^nu B] + tauh e^(-nu) curl [e^nu B] x B
!       with tauh = c/4pi e n_e
!       Our electric field is c*e^nu times the physical E-field:
!       E = eta Jnum + tauh Jnum x B
!         = eta curl[e^nu B] + tauh curl[e^nu B] x B
!
!     Physical temporal advance of B:
!       dB/dt = - c*curl[e^nu E]= - curl[Enum]
!
!-------------------------------------------------------------------------------
module magnetic_evolution

  ! Module imports -------------------------------------------------------------
  use input_params, only: etor_scheme, magnetic_advance_method
  use input_params, only: enable_burgers_correction, coeff_hyper, bgeom
  use grid, only: jevol, jmin
  use grid, only: np, nxi, neta

  implicit none

  ! Procedure pointer to the toroidal electric field evolution method.
  ! Either centered or upwind methods are implemented. It will be chosen at
  ! initialization of the magnetic evolution routine but by default it is set to
  ! the centered scheme.
  procedure(), pointer, save :: compute_etor_p => compute_Etor

  contains

    !---------------------------------------------------------------------------
    !> Magnetic field evolution.
    !> @brief It calls the time advance method.
    !
    !> @param[in] dtb_myr Timestep in Myr
    !> @param[in] iterb   Iteration number (used only for AB4 method)
    !
    !>  Code owners:
    !>    Daniele Viganò
    !>    Clara Dehman
    !---------------------------------------------------------------------------
    subroutine magnetic_evol(dtb_myr,iterb)

      ! Modules ----------------------------------------------------------------
      use grid, only: kmax, lmax
      use grid, only: br, bth, bphi, bm, bmed, benu
      use grid, only: jr, jth, jphi, er, eth, ephi

      implicit none

      ! Subroutine arguments ---------------------------------------------------
      real*8 :: dtb_myr
      integer :: iterb
 
      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------

      ! Index for loops.
      integer :: i, k, l

      ! Auxiliary fields for intermediate calculations and increments.
      real*8, dimension (0:nang+1, 0:np+2) :: xr, xth, xphi
      real*8, dimension (0:nang+1, 0:np+2) :: dbr, dbth, dbphi

      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------

      ! Select scheme for Etor computation.
      select case(etor_scheme)

        case ("ETOR_C")
          compute_etor_p => compute_Etor
        case ("ETOR_U")
          compute_etor_p => compute_Etor_upwind
        case default
          write(*,*) "<error>", &
                  & "[BEVOL]", &
                  & "Invalid Etor scheme: ", &
                  & etor_scheme
          stop

      end select

      ! Choose method for magnetic field evolution and evolve it.
      select case (magnetic_advance_method)

        case ("EUL")
          call euler(dtb_myr,dbr,dbth,dbphi)
        case ("EULA")
          call euler_alternate(dtb_myr,dbr,dbth,dbphi)
        case ("AB4")
          call AB4(dtb_myr,iterb)
        case ("AB4A")
          call AB4_alternate(dtb_myr,iterb)
        case ("RK4")
          call RK4(dtb_myr)
        case ("RK4V")
          call RK4_vecpot(dtb_myr)
        case ("SPEC")
          call magnetic_evol_spectral(dtb_myr)
        case default
          write(*,*) "<error>", &
                   & "[BEVOL]", &
                   & "Invalid magnetic advance method: ", &
                   & magnetic_advance_method
          stop

      end select

      ! Calculation of new electric currents
      xr   = 0d0
      xth  = 0d0
      xphi = 0d0
      do i=0,nang+1
        xr(i,jmin-1:) = br(i,jmin-1:)*benu(jmin-1:)
        xth(i,jmin-1:) = bth(i,jmin-1:)*benu(jmin-1:)
        xphi(i,jmin-1:) = bphi(i,jmin-1:)*benu(jmin-1:)
      enddo
      call curl_pol(xr,xth,jmin,jphi)
      call curl_phi(xphi,jmin,jr,jth)

      ! Calculation of new electric fields
      call compute_Epol(br,bth,bphi,jr,jth,jphi,er,eth)
      call compute_etor_p(br,bth,jr,jth,jphi,ephi)

      ! Calculation of new modulus and its average at the center of the thermal cell
      bm = dsqrt( br**2 + bth**2 + bphi**2 )
      do k=2,kmax
        do l=1,lmax
          bmed(k,l) = 0.25d0*bm(2*k-2,2*l-1) + &
    &       0.125*(bm(2*k-1,2*l-1) + bm(2*k-3,2*l-1) + bm(2*k-2,2*l) + bm(2*k-2,2*l-2)) + &
    &       0.0625*(bm(2*k-1,2*l-2) + bm(2*k-1,2*l) + bm(2*k-3,2*l-2) + bm(2*k-3,2*l-2)) 
        enddo
      enddo
        
    end subroutine magnetic_evol




    !---------------------------------------------------------------------------
    !> @brief Euler time advance.
    !
    !> @param[in]     dtb_myr   Timestep in Myr
    !> @param[out]    dbr       increment in br   in UNIT_B
    !> @param[out]    dbth      increment in bth  in UNIT_B
    !> @param[out]    dbphi     increment in bphi in UNIT_B
    !
    !> Code owners:
    !>    Daniele Viganò
    !---------------------------------------------------------------------------
    subroutine euler(dtb_myr,dbr,dbth,dbphi)

      ! Modules ----------------------------------------------------------------
      use grid, only: aphi, br, bth, bphi, benu
      use grid, only: jr, jth, jphi
      use grid, only: er, eth, ephi

      implicit none
      
      ! Subroutine arguments ---------------------------------------------------
      real*8, intent(in) :: dtb_myr
      real*8, dimension (0:nang+1,0:np+2), intent(out) :: dbr,dbth,dbphi

      ! Local variables --------------------------------------------------------

      ! Auxiliary variable for loops.
      integer :: i

      ! Auxiliary fields for intermediate calculations and increments.
      real*8, dimension (0:nang+1, 0:np+2) :: xr, xth, xphi

      ! ------------------------------------------------------------------------

      call compute_dBpol(dtb_myr, ephi, dbr, dbth)
      if (enable_burgers_correction) then
        call compute_dBtor_burgers(dtb_myr,br,bth,bphi,dbphi)
      else
        call compute_dBtor(dtb_myr,er,eth,dbphi)
      endif
      bphi = bphi + dbphi
      br  = br  + dbr
      bth = bth + dbth

      if (coeff_hyper > 0) then
        call add_hyperresistivity_bphi(dtb_myr, bphi)
      endif

      call magnetic_bc(aphi, br, bth, bphi)

    end subroutine euler


    !---------------------------------------------------------------------------
    !> @brief Subroutine for the Runge-Kutta 4th-order time-advance method
    !!
    !! @param[in]     dtb_myr   Timestep in Myr
    !!
    !!  Code owners:
    !!    Daniele Viganò
    !!    Clara Dehman
    !!------------------------------------------------------------------------------
    !! Currently, no hyper-resistivity is applied, but it could be called easily if needed.
    !!------------------------------------------------------------------------------
    subroutine RK4(dtb_myr)
  
      ! Module imports -------------------------------------------------------------
      use grid, only: aphi, br, bth, bphi
      use grid, only: jr, jth, jphi, er, eth, ephi
  
      implicit none
      
      ! Subroutine arguments -------------------------------------------------------
      real*8, intent(in) :: dtb_myr
  
      ! Local variables ------------------------------------------------------------
      ! Auxiliary fields for intermediate calculations and increments.
      real*8, dimension (0:nang+1, 0:np+2) :: dbr1,dbth1,dbphi1
      real*8, dimension (0:nang+1, 0:np+2) :: dbr2,dbth2,dbphi2
      real*8, dimension (0:nang+1, 0:np+2) :: dbr3,dbth3,dbphi3
      real*8, dimension (0:nang+1, 0:np+2) :: dbr4,dbth4,dbphi4
      real*8, dimension (0:nang+1, 0:np+2) :: brint,bthint,bphiint
  
      ! ----------------------------------------------------------------------------
  
      ! In the first substep the values of B, J and E coming from the previous magnetic loop are used
      call compute_Epol(br,bth,bphi,jr,jth,jphi,er,eth)
  
      call compute_etor_p(br,bth,jr,jth,jphi,ephi)
  
      call compute_dBpol(dtb_myr,ephi,dbr1,dbth1)
  
      if (enable_burgers_correction) then
        call compute_dBtor_burgers(dtb_myr,br,bth,bphi,dbphi1)
      else
        call compute_dBtor(dtb_myr,er,eth,dbphi1)
      endif
  
      brint   = br   + 0.5d0*dbr1
      bthint  = bth  + 0.5d0*dbth1
      bphiint = bphi + 0.5d0*dbphi1
      call magnetic_bc(aphi,brint,bthint,bphiint)
  
      ! In the following three sub-steps J and E are recalculated using Bint
      call RK4_substep(dtb_myr, brint,bthint,bphiint,dbr2,dbth2,dbphi2)
      brint   = br   + 0.5d0*dbr2
      bthint  = bth  + 0.5d0*dbth2
      bphiint = bphi + 0.5d0*dbphi2
      call magnetic_bc(aphi,brint,bthint,bphiint)
  
      call RK4_substep(dtb_myr, brint,bthint,bphiint,dbr3,dbth3,dbphi3)
      brint   = br   + 0.5d0*dbr3
      bthint  = bth  + 0.5d0*dbth3
      bphiint = bphi + 0.5d0*dbphi3
      call magnetic_bc(aphi,brint,bthint,bphiint)
  
      call RK4_substep(dtb_myr, brint,bthint,bphiint,dbr4,dbth4,dbphi4)
      br   = br   + (1./6.*dbr1   + 1./3.*dbr2   + 1./3.*dbr3   + 1./6.*dbr4   )
      bth  = bth  + (1./6.*dbth1  + 1./3.*dbth2  + 1./3.*dbth3  + 1./6.*dbth4  )
      bphi = bphi + (1./6.*dbphi1 + 1./3.*dbphi2 + 1./3.*dbphi3 + 1./6.*dbphi4 )
  
      if (coeff_hyper > 0) then
        call add_hyperresistivity_bphi(dtb_myr, bphi)
      endif
  
      call magnetic_bc(aphi,br,bth,bphi)
  
    end subroutine RK4
  
  
    !---------------------------------------------------------------------------
    !> @brief Subroutine for the substeps of the Runge-Kutta time-advance method
    !!
    !! @param[in]     dtb_myr   Timestep in Myr
    !! @param[in]     brint     br   in UNIT_B
    !! @param[in]     bthint    bth  in UNIT_B
    !! @param[in]     bphiint   bphi in UNIT_B
    !! @param[out]    dbr       increment in br   in UNIT_B
    !! @param[out]    dbth      increment in bth  in UNIT_B
    !! @param[out]    dbphi     increment in bphi in UNIT_B
    !!
    !! Code owners:
    !!    Daniele Viganò
    !!    Clara Dehman
    !!------------------------------------------------------------------------------
    !! Important: do NOT use the B field from grid, since they are in general different
    !!------------------------------------------------------------------------------
    subroutine RK4_substep(dtb_myr,brint,bthint,bphiint,dbr,dbth,dbphi)
  
      ! Module imports -------------------------------------------------------------
      use grid, only: jr, jth, jphi, er, eth, ephi, benu
  
      implicit none
      
      ! Subroutine arguments -------------------------------------------------------
      real*8, intent(in) :: dtb_myr
  
      ! Local variables ------------------------------------------------------------
      ! Auxiliary fields for intermediate calculations and increments.
      integer :: i
      real*8, dimension (0:nang+1, 0:np+2), intent(in) :: brint, bthint, bphiint
      real*8, dimension (0:nang+1, 0:np+2), intent(out) :: dbr, dbth, dbphi
      real*8, dimension (0:nang+1,0:np+2) :: xr, xth, xphi
  
      ! ----------------------------------------------------------------------------
  
      xr   = 0d0
      xth  = 0d0
      xphi = 0d0
      do i=0,nang+1
        xr(i,jmin-1:) = brint(i,jmin-1:)*benu(jmin-1:)
        xth(i,jmin-1:) = bthint(i,jmin-1:)*benu(jmin-1:)
        xphi(i,jmin-1:) = bphiint(i,jmin-1:)*benu(jmin-1:)
      enddo
      call curl_pol(xr,xth,jmin,jphi) 
      call curl_phi(xphi,jmin,jr,jth)
  
      call compute_Epol(brint,bthint,bphiint,jr,jth,jphi,er,eth)
  
      call compute_etor_p(brint,bthint,jr,jth,jphi,ephi)
  
      call compute_dBpol(dtb_myr,ephi,dbr,dbth)
      
      if (enable_burgers_correction) then
        call compute_dBtor_burgers(dtb_myr,brint,bthint,bphiint,dbphi)
      else
        call compute_dBtor(dtb_myr,er,eth,dbphi)
      endif
  
    end subroutine RK4_substep
  
  
  
    !!------------------------------------------------------------------------------
    !> @brief Subroutine for the substeps of the Adams-Bashforth 4th-order time-advance method
    !!
    !! @param[in]     dtb_myr   Timestep in Myr
    !! @param[in]     iterb     Number of iteration
    !!
    !! Code owners:
    !!    Clara Dehman
    !!    Daniele Viganò
    !!------------------------------------------------------------------------------
    subroutine AB4(dtb_myr,iterb)
  
      ! Module imports -------------------------------------------------------------
      use grid, only: aphi, br, bth, bphi
      use grid, only: er, eth, ephi
  
      implicit none
      
      ! Subroutine arguments -------------------------------------------------------
      real*8, intent(in) :: dtb_myr
      integer, intent(in) :: iterb
  
      ! Local variables ------------------------------------------------------------
      ! Auxiliary fields for intermediate calculations and previous increments.
      real*8, dimension (0:nang+1, 0:np+2) :: dbr0, dbth0, dbphi0
      real*8, dimension (:, :), allocatable, save :: dbr1, dbth1, dbphi1
      real*8, dimension (:, :), allocatable, save :: dbr2, dbth2, dbphi2
      real*8, dimension (:, :), allocatable, save :: dbr3, dbth3, dbphi3
   
      ! First three iteration with Euler alternate
      if (iterb <= 3) then 
        call euler_alternate(dtb_myr,dbr0,dbth0,dbphi0)
        if (iterb == 0) then
          allocate(dbr1(0:nang + 1, 0:np + 2))
          allocate(dbth1(0:nang + 1, 0:np + 2))
          allocate(dbphi1(0:nang + 1, 0:np + 2))
          allocate(dbr2(0:nang + 1, 0:np + 2))
          allocate(dbth2(0:nang + 1, 0:np + 2))
          allocate(dbphi2(0:nang + 1, 0:np + 2))
          allocate(dbr3(0:nang + 1, 0:np + 2))
          allocate(dbth3(0:nang + 1, 0:np + 2))
          allocate(dbphi3(0:nang + 1, 0:np + 2))
  
          dbr3=dbr0  
          dbth3=dbth0
          dbphi3=dbphi0
        end if 
        if (iterb == 1) then
          dbr2=dbr0  
          dbth2=dbth0
          dbphi2=dbphi0
        end if 
        if (iterb == 2) then
          dbr1=dbr0  
          dbth1=dbth0
          dbphi1=dbphi0
        endif
      end if 
  
      ! From the third iteration it is possible to use the 4th order Adams-Bashforth method
      if (iterb >= 3) then
  
        call compute_dBpol(dtb_myr,ephi,dbr0,dbth0)
        if (enable_burgers_correction) then
          call compute_dBtor_burgers(dtb_myr,br,bth,bphi,dbphi0)
        else
          call compute_dBtor(dtb_myr,er,eth,dbphi0)
        endif
  
        br = br+ (55./24.*dbr0-59./24.*dbr1+37./24.*dbr2-9./24.*dbr3)
        dbr3=dbr2
        dbr2=dbr1
        dbr1=dbr0
    
        bth = bth + (55./24.*dbth0-59./24.*dbth1+37./24.*dbth2-9./24.*dbth3)
        dbth3=dbth2
        dbth2=dbth1
        dbth1=dbth0
  
        bphi = bphi + (55./24.*dbphi0-59./24.*dbphi1+37./24.*dbphi2-9./24.*dbphi3)
        dbphi3=dbphi2
        dbphi2=dbphi1
        dbphi1=dbphi0
    
        call magnetic_bc(aphi,br,bth,bphi)
  
      end if 
  
    end subroutine AB4

  
  
    !!------------------------------------------------------------------------------
    !> @brief Subroutine for poloidal magnetic field increment, given a Ephi
    !!
    !! @param[in]     dtb_myr   Timestep in Myr
    !! @param[in]     ephi      toroidal electric field
    !! @param[out]    dbr       increment in br   in UNIT_B
    !! @param[out]    dbth      increment in bth  in UNIT_B
    !!
    !! Code owners:
    !!    Daniele Viganò
    !!    Clara Dehman
    !!------------------------------------------------------------------------------
    subroutine compute_dBpol(dtb_myr,ephi,dbr,dbth)
  
      ! Module imports -------------------------------------------------------------
      use grid, only: lphi, arear, areath
  
      implicit none
  
      ! Subroutine arguments -------------------------------------------------------
      real*8, intent(in) :: dtb_myr
      real*8, dimension (0:nang+1, 0:np+2), intent(in) :: ephi
      real*8, dimension (0:nang+1, 0:np+2), intent(out) :: dbr,dbth
  
      ! Local variables ------------------------------------------------------------
      ! Auxiliary index for loops
      integer :: i,j
  
      ! ----------------------------------------------------------------------------
      dbr = 0.d0
      dbth = 0.d0
      do i=2,nang-1
        dbr(i,jevol:np) = -dtb_myr*(ephi(i+1,jevol:np)*lphi(i+1,jevol:np) &
     &                           -ephi(i-1,jevol:np)*lphi(i-1,jevol:np))/arear(i,jevol:np)
      enddo
      dbr(1,jevol:np) = -dtb_myr*(ephi(2,jevol:np)*lphi(2,jevol:np))/arear(1,jevol:np)
      dbr(nang,jevol:np) = -dtb_myr*(-ephi(nang-1,jevol:np)*lphi(nang-1,jevol:np))/ &
     &                      arear(nang,jevol:np)
      dbr(0,:) = dbr(2,:)
      dbr(nang+1,:) = dbr(nang-1,:)
  
      do j=jevol+1,np-1
        dbth(2:nang-1,j) = -dtb_myr*(ephi(2:nang-1,j-1)*lphi(2:nang-1,j-1) &
     &                          -ephi(2:nang-1,j+1)*lphi(2:nang-1,j+1))/areath(2:nang-1,j)
      enddo
      dbth(0,:) = -dbth(2,:) 
      dbth(nang+1,:) = -dbth(nang-1,:)
  
    end subroutine compute_dBpol
  
  
    !!------------------------------------------------------------------------------
    !> @brief Subroutine for toroidal magnetic field increment, given a Er and Eth
    !!
    !! @param[in]     dtb_myr   Timestep in Myr
    !! @param[in]     er,eth    Poloidal electric field
    !! @param[out]    dbphi     increment in bphi in UNIT_B
    !!
    !! Code owners:
    !!    Daniele Viganò
    !!    Clara Dehman
    !!------------------------------------------------------------------------------
    subroutine compute_dBtor(dtb_myr,er,eth,dbphi)
  
      ! Module imports -------------------------------------------------------------
      use grid, only: lr, lth, areaphi
  
      implicit none
  
      ! Subroutine arguments -------------------------------------------------------
      real*8, intent(in) :: dtb_myr
      real*8, dimension (0:nang+1, 0:np+2), intent(inout) :: er,eth
      real*8, dimension (0:nang+1, 0:np+2), intent(out) :: dbphi
  
      ! Local variables ------------------------------------------------------------
      ! Auxiliary index for loops
      integer :: i,j
  
      ! ----------------------------------------------------------------------------
  
      dbphi=0.d0
  
      do j=jevol+1,np-1
        do i=2,nang-1
          dbphi(i,j) = dtb_myr*((er(i+1,j)-er(i-1,j))*lr(j)  &
   &        -(eth(i,j+1)*lth(j+1)-eth(i,j-1)*lth(j-1)))/areaphi(j)
        enddo
      enddo
      dbphi(0,:) = -dbphi(2,:)
      dbphi(nang+1,:) = -dbphi(nang-1,:)
  
    end subroutine compute_dBtor

    
    !!------------------------------------------------------------------------------
    !> @brief Subroutine for poloidal electric field
    !!
    !! @param[in]     dtb_myr             Timestep in Myr
    !! @param[in]     brin,bthin,bphiin   magnetic field in UNIT_B
    !! @param[in]     jrin,jthin,jphiin   Electrical currents in UNIT_B/UNIT_R
    !! @param[out]    erin,ethin          poloidal electric field
    !!
    !! Code owners:
    !!    Daniele Viganò
    !!    Clara Dehman
    !!------------------------------------------------------------------------------
    !! Do not put the B ad J fields from grid module, they can be different
    !!------------------------------------------------------------------------------
    subroutine compute_Epol(brin,bthin,bphiin,jrin,jthin,jphiin,erin,ethin)
  
      ! Module imports -------------------------------------------------------------
      use grid, only: fh, ia_hall, etab
  
      implicit none
  
      ! Subroutine arguments -------------------------------------------------------
      real*8, dimension (0:nang+1,0:np+2), intent(in) :: brin,bthin,bphiin,jrin,jthin,jphiin
      real*8, dimension (0:nang+1,0:np+2), intent(out) :: erin,ethin
   
      ! Local variables ------------------------------------------------------------
      ! Auxiliary index for loops
      integer :: i
  
      ! ----------------------------------------------------------------------------
  
      erin=0d0
      ethin=0d0
  
      do i=1,nang
        erin(i,jevol:np-1)=etab(i,jevol:np-1)*jrin(i,jevol:np-1) &
     &        -ia_hall(jevol:np-1)*fh(jevol:np-1)* &
     &        (jphiin(i,jevol:np-1)*bthin(i,jevol:np-1)-jthin(i,jevol:np-1)*bphiin(i,jevol:np-1))
  
        ethin(i,jevol:np)=etab(i,jevol:np)*jthin(i,jevol:np) &
     &   +ia_hall(jevol:np)*fh(jevol:np)*(jphiin(i,jevol:np)*  &
     &   brin(i,jevol:np)-jrin(i,jevol:np)*bphiin(i,jevol:np))
      enddo
      ethin(0,:) = -ethin(2,:)
      ethin(nang+1,:) = -ethin(nang-1,:)
      erin(0,:) = erin(2,:)
      erin(nang+1,:) = erin(nang-1,:)
  
    end subroutine compute_Epol
  
  
    !!------------------------------------------------------------------------------
    !> @brief Subroutine for toroidal electric field
    !!
    !! @param[in]     dtb_myr             Timestep in Myr
    !! @param[in]     brin,bthin          magnetic field in UNIT_B
    !! @param[in]     jrin,jthin,jphiin   Electrical currents in UNIT_B/UNIT_R
    !! @param[out]    ephiin              toroidal electric field
    !!
    !! Code owners:
    !!    Daniele Viganò
    !!    Clara Dehman
    !!------------------------------------------------------------------------------
    !! Do not put the B ad J fields from grid module, they can be different
    !!------------------------------------------------------------------------------
    subroutine compute_Etor(brin,bthin,jrin,jthin,jphiin,ephiin)
  
      ! Module imports -------------------------------------------------------------
      use grid, only: fh, ia_hall, etab
  
      implicit none
  
      ! Subroutine arguments -------------------------------------------------------
      real*8, dimension (0:nang+1,0:np+2), intent(in) :: brin,bthin,jrin,jthin,jphiin
      real*8, dimension (0:nang+1,0:np+2), intent(out) :: ephiin
   
      ! Local variables ------------------------------------------------------------
      ! Auxiliary index for loops
      integer :: i
  
      ! ----------------------------------------------------------------------------
  
      ephiin=0d0
      do i=1,nang
        ephiin(i,jevol:np)=etab(i,jevol:np)*jphiin(i,jevol:np)  &
       &        -ia_hall(jevol:np)*fh(jevol:np)*(-jrin(i,jevol:np)*bthin(i,jevol:np)+jthin(i,jevol:np)*brin(i,jevol:np))
      enddo
      ephiin(0,:) = -ephiin(2,:)
      ephiin(nang+1,:) = -ephiin(nang-1,:)
  
    end subroutine compute_Etor

    

    !!------------------------------------------------------------------------------
    !> @brief This subroutines provide the contravatiant component of the curl of a vector. 
    !> For the geometrical elements, we are using the covariant components 
    !!
    !! @param[in]  l_r, l_xi, l_eta          covariant length elements 
    !! @param[in]  area_r, area_xi, area_eta covariant surface elements 
    !! @param[in]   fr                       contravariant radial component
    !! @param[in]   fxi                      contravariant xi component
    !! @param[in]   feta                     contravariant eta component
    !! @param[out]  curlfr                   contravariant radial component of curlf
    !! @param[out]  curlfxi                  contravariant xi component of curlf
    !! @param[out]  curlfeta                 contravariant eta component of curlf
    !!
    !! Code owners:
    !!   Clara Dehman
    !!------------------------------------------------------------------------------
    subroutine curl_cont(fr,fxi,feta,jminin,curlfr,curlfxi,curlfeta)
  
      ! Module imports -------------------------------------------------------------
      use grid, only: l_r, l_xi, l_eta, area_r, area_xi, area_eta
  
      implicit none
  
      ! Subroutine arguments -------------------------------------------------------
      integer, intent(in) :: jminin
      real*8, dimension (0:np+2,0:nxi+1,0:neta+1), intent(in) :: fr,fxi,feta
      real*8, dimension (0:np+2,0:nxi+1,0:neta+1), intent(out) :: curlfr, curlfxi, curlfeta
  
      ! Local variables ------------------------------------------------------------
      ! Auxiliary index for loops
      integer :: i,j,k
  
      ! ----------------------------------------------------------------------------
  
      curlfr = 0.d0
      curlfxi = 0.d0
      curlfeta = 0.d0

      do j=1,nxi
        do k=1,neta
        curlfr(jminin:np,j,k)=(l_xi(jminin:np,j,k-1)*fxi(jminin:np,j,k-1) &
       & - l_xi(jminin:np,j,k+1)*fxi(jminin:np,j,k+1) &
       & + l_eta(jminin:np,j+1,k)*feta(jminin:np,j+1,k) &
       & - l_eta(jminin:np,j-1,k)*feta(jminin:np,j-1,k))/area_r(jminin:np,j,k)
        enddo
      enddo
  
      do i=jminin,np
        do k=1,neta
        curlfxi(i,1:nxi,k)=(l_r(i,1:nxi,k+1)*fr(i,1:nxi,k+1) &
       & - l_r(i,1:nxi,k-1)*fr(i,1:nxi,k-1) &
       & + l_eta(i-1,1:nxi,k)*feta(i-1,1:nxi,k) &
       & - l_eta(i+1,1:nxi,k)*feta(i+1,1:nxi,k))/area_xi(i,1:nxi,k)
        enddo
      enddo

      do i=jminin,np
        do j=1,nxi
        curlfeta(i,j,1:neta)=(l_xi(i+1,j,1:neta)*fxi(i+1,j,1:neta) &
       & - l_xi(i-1,j,1:neta)*fxi(i-1,j,1:neta) &
       & + l_r(i,j-1,1:neta)*fr(i,j-1,1:neta) &
       & - l_r(i,j+1,1:neta)*fr(i,j+1,1:neta))/area_eta(i,j,1:neta)
        enddo
      enddo
  
    end subroutine curl_cont
  
  
    !!-----------------------------------------------------------------------
    !> @brief Subroutine for the Boundary Conditions
    !!
    !! @param[inout]     aphiin    Aphi in UNIT_B*UNIT_R
    !! @param[inout]     brin      br   in UNIT_B
    !! @param[inout]     bthin     bth  in UNIT_B
    !! @param[inout]     bphiin    bphi in UNIT_B
    !!
    !!  Code owners:
    !!    Daniele Viganò
    !!    Clara Dehman
    !!
    !!-----------------------------------------------------------------------
    !! Axis at i=1,nang.
    !! Inner boundary at j=jevol, and surface at j=np.
    !! BC at the axis: Fth=Fphi=0, where F is any vector (we are in axial symmetry)
    !!-----------------------------------------------------------------------
    subroutine magnetic_bc(aphiin,brin,bthin,bphiin)
  
      ! Module imports -------------------------------------------------------------
      use grid, only: lphi, cth, arear, jcore
      use legpol, only: blout, getbl
  
      implicit none 
  
      ! Subroutine arguments -------------------------------------------------------
      real*8, dimension (0:nang+1, 0:np+2), intent(inout) :: aphiin, brin, bthin, bphiin
  
      ! Local variables ------------------------------------------------------------
      ! Auxiliary variable for loops
      integer :: i,j
  
      ! ----------------------------------------------------------------------------
  
      ! Spectral decomposition of the radial magnetic field at the surface.
      call getbl(nang-1,cth(1:nang-1),brin(1:nang-1,np))
  
      ! Retrieve the vector potential aphi from the poloidal field
      aphiin = 0d0
      do j=jmin,np
        do i=2,nang/2+1
          aphiin(i,j)=(aphiin(i-2,j)*lphi(i-2,j)+brin(i-1,j)*arear(i-1,j))/lphi(i,j)
          aphiin(nang+1-i,j)=(aphiin(nang-i+3,j)*lphi(nang-i+3,j)  &
       &      - brin(nang+2-i,j)*arear(nang+2-i,j))/lphi(nang+1-i,j)
        enddo
      enddo
      aphiin(0,:)=-aphiin(2,:)
      aphiin(nang+1,:)=-aphiin(nang-1,:)
  
      ! Boundary condition for aphi and bphi at the surface.
      ! Spectral vacuum BC.  Other options now disabled (to be implemented)
      call potential_magnetic_legendre(blout, aphiin, bphiin)
  
      ! Boundary condition for aphi and bphi in the core
      if (bgeom <= 1) then
        aphiin(:,0:jcore) = 0.d0
        bphiin(:,0:jcore) = 0.d0
        aphiin(:,jcore+1) = 0.5d0*aphiin(:,jcore+2)
        bphiin(:,jcore+1) = 0.5d0*bphiin(:,jcore+2)
      endif
  
      ! retreive poloidal magnetic field components from Aphi.
      call curl_phi(aphiin,jmin,brin,bthin)
  
      ! Boundary condition for aphi and bphi on the axis
      bphiin(1,:)=0
      bphiin(nang,:)=0
      bphiin(0,:)=-bphiin(2,:)
      bphiin(nang+1,:)=-bphiin(nang-1,:)
    
    end subroutine magnetic_bc  
  
  
    !!-----------------------------------------------------------------------
    !> @brief Subroutine for the Boundary Conditions when evolving Aphi
    !!
    !! @param[inout]     aphiin    Aphi in UNIT_B*UNIT_R
    !! @param[inout]     brin      br   in UNIT_B
    !! @param[inout]     bthin     bth  in UNIT_B
    !! @param[inout]     bphiin    bphi in UNIT_B
    !!
    !!  Code owners:
    !!    Daniele Viganò
    !!    Clara Dehman
    !!
    !!-----------------------------------------------------------------------
    !! Axis at i=1,nang.
    !! Inner boundary at j=jevol, and surface at j=np.
    !!-----------------------------------------------------------------------
    subroutine magnetic_bc_vecpot(aphiin,brin,bthin,bphiin)
  
      ! Module imports -------------------------------------------------------------
      use grid, only: cth, jcore
      use legpol, only: blout, getbl
  
      implicit none 
  
      ! Subroutine arguments -------------------------------------------------------
      real*8, dimension (0:nang+1, 0:np+2), intent(inout) :: aphiin, brin, bthin, bphiin
  
      ! Local variables ------------------------------------------------------------
      ! None
  
      ! ----------------------------------------------------------------------------
  
      ! Boundary condition for aphi and bphi at the surface.
      ! Spectral vacuum BC.  Other options now disabled (to be implemented)
      call getbl(nang-1,cth(1:nang-1),brin(1:nang-1,np))
      call potential_magnetic_legendre(blout, aphiin, bphiin)
  
      ! Boundary condition for aphi and bphi in the core
      if (bgeom <= 1) then
        aphiin(:,0:jcore) = 0.d0
        bphiin(:,0:jcore) = 0.d0
      endif
  
      ! retreive poloidal magnetic field components from Aphi.
      call curl_phi(aphiin,jmin,brin,bthin)
  
      ! Boundary condition for aphi and bphi on the axis
      bphiin(1,:)=0
      bphiin(nang,:)=0
      bphiin(0,:)=-bphiin(2,:)
      bphiin(nang+1,:)=-bphiin(nang-1,:)
    
    end subroutine magnetic_bc_vecpot
  
    
    !!-----------------------------------------------------------------------  
    !> @brief Subroutine for the potential solution outside the star
    !!
    !! @param[in]     blout     Magnetic multipole weights
    !! @param[in]     aphiin    Azimuthal component of vector potential, in UNIT_B*UNIT_R
    !! @param[in]     bthin     bth in UNIT_B
    !! @param[in]     bphiin    bphi in UNIT_B
    !!
    !!  Code owners:
    !!    Daniele Viganò
    !!------------------------------------------------------------------------------
    !! Important: here, do NOT use the B field from grid, since they are in general different
    !!------------------------------------------------------------------------------
    subroutine potential_magnetic_legendre(blout, aphiin, bphiin)
  
      ! Module imports -------------------------------------------------------------
      use legpol, only: nleg, dpln
      use grid, only: rb
  
      implicit none
  
      ! Subroutine arguments -------------------------------------------------------
      real*8, dimension (0:nleg), intent(in) :: blout
      real*8, dimension (0:nang+1, 0:np+2), intent(inout) :: aphiin, bphiin
  
      ! Local variables ------------------------------------------------------------
      ! Auxiliary variable for loops
      integer :: n
      ! Auxiliary variable for sums
      real*8, dimension (1:nang) :: sum1, sum2
  
      ! ----------------------------------------------------------------------------
  
      bphiin(:,np-1)=0.5d0*bphiin(:,np-2)
      bphiin(:,np:np+2)=0d0
  
      sum1 = 0d0
      sum2 = 0d0
      do n=1,nleg
        sum1 = sum1 + (blout(n)/n)*dpln(1:nang,n)*rb(np)**(n+2)/rb(np+1)**(n+1)
        sum2 = sum2 + (blout(n)/n)*dpln(1:nang,n)*rb(np)**(n+2)/rb(np+2)**(n+1)
      enddo
      aphiin(1:nang,np+1)=sum1
      aphiin(1:nang,np+2)=sum2
      aphiin(0,:)=-aphiin(2,:)
      aphiin(nang+1,:)=-aphiin(nang-1,:)
      
    end subroutine potential_magnetic_legendre
  
    
    
  
    !!------------------------------------------------------------------------------
    !> @brief Subroutine for the Joule heating calculation
    !!
    !!  Code owners:
    !!    Daniele Viganò
    !!------------------------------------------------------------------------------
    !! The subroutine calculates also the Joule estimated from the shock
    !! See Viganò et al. 2012 paper for details and notation
    !! The Joule rate is defined as
    !! Q_joule = e^2nu*eta*J^2
    !!------------------------------------------------------------------------------  
    subroutine compute_joule()
  
      ! Module imports -------------------------------------------------------------
      use constants, only: UNIT_JOULE
      use grid, only: kmax, lmax, lr, lth, jcore
      use grid, only: bphi
      use grid, only: jr, jth, jphi
      use grid, only: etab, benu
      use grid, only: q_joule, q_joule_shock, j2, lamr, lamth
  
      implicit none
  
      ! Local variables ------------------------------------------------------------
      ! Auxiliary variable for loops
      integer :: i,j,k,l
      ! Auxiliary variables for shock correction
      real*8 jump

      ! ----------------------------------------------------------------------------

      j2 = 0d0
      ! In the crust-confined case, 
      ! the supercurrents arising from the discontinuity at the interface are neglected.
      ! They are located at jevol-2,jevol-1,jevol
      j2(:,jevol+1:) = jr(:,jevol+1:)**2 + jth(:,jevol+1:)**2 + jphi(:,jevol+1:)**2

      q_joule=0d0
      q_joule_shock=0d0
      do k=2,kmax
        i=2*k-2

        ! Supercurrent contribution, with the eta of the core, for crust-confined cases.
        ! It would be relevant only if associated to etab of the crust
        if (jevol == jcore+1) then
          q_joule_shock(k,jevol/2)= UNIT_JOULE*etab(i,jcore)* &
            &     ( jr(i,jcore)**2 + jth(i,jcore)**2 + jphi(i,jcore)**2 )
        endif
        ! The first cell where Joule is to be calculated is at l=(jevol+2)/2, which
        ! corresponds to the center j=jevol+1 if jevol is even (jcore is odd)
        !                           j=jevol   if jevol is odd  (jcore is even)
        do l= (jevol+2)/2,lmax
          j=2*l-1
        ! Considering the weighted contributions of the 9 points
        ! (center, edges, corners) corresponding to the thermal cell
          q_joule(k,l)= - UNIT_JOULE*( 0.25d0*etab(i,j)*j2(i,j) + &
          &           0.125d0*(etab(i+1,j)*j2(i+1,j) + etab(i-1,j)*j2(i-1,j) + &
          &                    etab(i,j+1)*j2(i,j+1) + etab(1,j-1)*j2(i,j-1) ) + & 
          &          0.0625d0*(etab(i+1,j+1)*j2(i+1,j+1) + etab(i-1,j+1)*j2(i-1,j+1) + &
          &                    etab(i+1,j-1)*j2(i+1,j-1) + etab(i-1,j-1)*j2(i-1,j-1) ) )
          ! Shock correction calculation for Bphi starts
          jump=0.5d0*(bphi(i+1,j)-bphi(i-1,j))
          if(i /= 2 .and. i /= nang .and. jump*lamth(j) < 0d0) then
            q_joule_shock(k,l) = UNIT_JOULE*2d0/3d0*lamth(j)*jump**3*benu(j)/lth(j)
          endif
          jump=0.5d0*(bphi(i,j+1)-bphi(i,j-1))
          if(jump*lamr(i,j) < 0d0) then
            q_joule_shock(k,l)= q_joule_shock(k,l) + &
          &          UNIT_JOULE*2d0/3d0*lamr(i,j)*jump**3*benu(j)/lr(j)
          endif
  
        enddo
      enddo
  
    end subroutine compute_joule
  
end module magnetic_evolution

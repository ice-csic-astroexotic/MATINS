!-------------------------------------------------------------------------------
! Magneto Thermal 3D
!-------------------------------------------------------------------------------
! 
!> @author
!> Daniele Viganò
!> Jose Pons Botella
!
!> Neutron star structure.
!-----------------------------------------------------------


subroutine build_structure()
  ! Modules ----------------------------------------------------------------
  use constants, only : PI, UNIT_R, UNIT_B, T_YEAR, STATC_TO_C, E_CHARGE
  use grid

  implicit none
  ! Subroutine arguments ---------------------------------------------------
  ! None.
  ! Local constants --------------------------------------------------------
  integer, parameter :: nrh = 5000    ! Number of points used in the homogeneous radial structure
  ! Local variables --------------------------------------------------------
  real*8 rh(0:nrh), dr_temp(0:nrh), pressh(0:nrh), massh(0:nrh), elamh(0:nrh), enuh(0:nrh)
  real*8 r_cut, r_core, arg
  ! variable radial grid size
  real*8 dr
  ! supefluid gaps and critical temperatures
  real*8 gapT0, Tc
  ! Auxiliary variables
  integer i, l, i_cut, i_core
  real*8 rho_0, nb_0, zz_0, aa_0, yn_0, xh_0, ymu_0
  real*8 alfa

  ! ------------------------------------------------------------------------
  ! Load the EoS table
  call init_eos_tab()
  ! Solve the TOV equations for a very fine uniform radial spacing
  ! It provides the profilies of quantities "...h", the surface gravity, 
  ! and the last radius of the domain (crust-envelope interface definition)
  call solve_tov_structure(nrh, rh, i_cut, pressh, massh, enuh, elamh)
  r_cut  = rh(i_cut)

  ! Look for the location of the core
  r_core = 0
  i = 1
  do while (r_core == 0)
    ! I use for comfort the index 0 of the array that I will use later to store the variables
    call geteost(pressh(i), rho_0, nb_0, zz_0, aa_0, yn_0, xh_0, ymu_0)
    if (xh_0 .gt. 0) then
      r_core = rh(i - 1)
    endif
    i = i+1
  end do

  moment_inertia = 0d0
  fh = 0d0
  ! Define the crust and core gaps.
  tccru = 0.d0
  tcn = 0.d0
  tcp = 0.d0
  Zimp = 0d0

  ! The first nr points are in the core, used only for calculating the microsphysical quantities

  dr = (r_cut-r_core)/dble(nr)
  do i = 0, nr+1
    r(i) = r_core + dble(i) * dr
    rtot(ncore + i) = r(i)
  end do

  dr = r_core/dble(ncore)
  do i = 0, ncore-1
    rtot(i) = dble(i) * dr
  end do

  do i = ncore + nr, 1, -1 
    ! Interpolate the quantities found by the TOV solver in the new radial grid
    dr_temp = dabs(rtot(i) - rh)
    l = MINLOC(dr_temp,DIM=1)
    !write(*,*), "l ", l
    alfa = (rtot(i)-rh(l))/(rh(l+1)-rh(l-1))
    press(i) = pressh(l) + alfa*( pressh(l+1) - pressh(l-1))
    mass(i)  = massh(l)   + alfa*( massh(l+1)  - massh(l-1))
    elambda_tot(i) = elamh(l)  + alfa*( elamh(l+1)  - elamh(l-1) )
    enu_tot(i)  = enuh(l)   + alfa*( enuh(l+1)   - enuh(l-1)  )
    ! Retrieve density and composition
    call geteost(press(i), rho(i), nb(i), zz(i), aa(i), yn(i), xh(i), ymu(i))

    if (i == ncore + 1) then
      xh(i) = xh(i+1)
      aa(i) = aa(i+1)
      zz(i) = zz(i+1)
      ymu(i) = ymu(i+1)
    endif

    ! Calculate the moment of inertia with the given profile
    moment_inertia = moment_inertia + rho(i)*4d0*PI*rtot(i)**4*(rtot(i)-rtot(i-1))*UNIT_R**5
    if (xh(i) == 0d0) then
      ! Charge neutrality in the core.
      yp(i) = 1.d0-yn(i)
      ye(i) = yp(i) - ymu(i)
    else
      ! Zero proton fraction in the crust.
      yp(i) = 0d0
      ye(i) = xh(i) * zz(i) / aa(i)
    end if
    nn(i) = nb(i)*yn(i)
    npr(i) = nb(i)*yp(i)
    ne(i) = nb(i)*ye(i)
    nmu(i) = nb(i)*ymu(i)
    ! Fermi momenta in fm**-1 for n, p, e
    kFn(i) = (3*PI**2*nn(i))**(1.d0/3.d0)
    kFp(i) = (3*PI**2*npr(i))**(1.d0/3.d0)
    kFe(i) = (3*PI**2*ne(i))**(1.d0/3.d0)
    kFmu(i) = (3*PI**2*nmu(i))**(1.d0/3.d0)

    ! Superfluid gaps
    ! For neutron 1s0 crust.
    if (superfluid_n_crust /= '0') then
      if(xh(i) > 0.d0)then
        call gapmodel(superfluid_n_crust, 1, kFn(i), gapT0, Tc)
        tccru(i)=Tc*1.d-8
        gapn_crust(i) = gapT0
      end if
    end if
    ! For neutron 3p2 core.
    if (superfluid_n_core /= '0') then
      if(xh(i) == 0.d0) then
        call gapmodel(superfluid_n_core, 2, kFn(i), gapT0, Tc)
        tcn(i) = Tc*1.d-8
        gapn_core(i) = gapT0
      end if
    end if
    ! For proton 1s0 core.
    if (superfluid_p_core /= '0') then
      if(xh(i) == 0.d0)then
        call gapmodel(superfluid_p_core, 1, kFp(i), gapT0, Tc)
        tcp(i)=Tc*1.d-8
        gapp_core(i) = gapT0
      end if
    end if

    if ( i >= ncore) then

      call get_Zimp(nb(i),Zimp(i-ncore))

    ! Hall factor fh = c/(4*PI*e*n_e)
    ! we need it in code units [km**2/(Myr*10**12 G)]
    ! 
    ! e/c in cgs is 1.602e-19/STATC_TO_C (coming from unit conversion)
    ! n_e in the code is in [fm**-3]
    ! Hall factor = c / 4*pi*e*n_e
    !             = c*1.d12 / (4*pi*1.6d-19*1 statC*(rho*Ye/m_b)
    !             = 3d10*1.67d-27 / (4*pi*1.6d-19*3d10*rho*Ye)
    !             = 1d8*1.673 / (1.602*4*pi)
    ! Conversion from cgs to code units: k = 1d6*T_YEAR*UNIT_B/UNIT_R**2
    ! fh(i)[km**2/(Myr*1e12 G)] = 1.568d-5 / ne[fm**-3]
      fh(i-ncore) = 1d6*T_YEAR*UNIT_B*STATC_TO_C/(UNIT_R**2) / (4d0*PI*ne(i)*E_CHARGE*1d39)
    endif

  end do ! i


  vol_shell(1:ncore) = 2d0*PI/3d0*(rtot(2:ncore+1)**3 - rtot(0:ncore-1)**3)/elambda_tot(1:ncore)
  ! The ghost (center, above surface) values of the metric factor are defined manually
  enu(0:nr) = enu_tot(ncore:ncore+nr)
  enu_core(1:ncore) = enu_tot(1:ncore)
  elambda(0:nr) = elambda_tot(ncore:ncore+nr)
  enu(nr+1) = enu(nr)
  elambda(nr+1) = elambda(nr)

  call effective_mass()

  ! Output structure grid and gaps to files, write screen output information.
  call output_structure()

end subroutine build_structure


!-----------------------------------------------------------
! @brief The subroutine defines the derivatives in the first integrations
!-----------------------------------------------------------
!! @param[in]  nrh      Maximum number of points in the radius
!! @param[out] radius   homogeneous radii of the output profiles
!! @param[out] i_cut    index of output profiles identifying the p_cut (crust/envelope location)
!! @param[out] pressh   profile of pressure
!! @param[out] massh    profile of enclosed mass
!! @param[out] enuh     profile of timelike element of metric e**nu
!! @param[out] elamh    profile of radial stretch metric factor e**lambda
!---------------------------------------------------------------
! --------------------------------------------------------------------------
!> Main structure calculation routine.
!> @brief Main interface to calculate the structure. This subroutine:
!>          1) Defines the radial grid (irregular, finer in the crust)
!>          2) Calls the TOV solver to calculate the star structure
!>          3) Calls the EOS to obtain needed variables (composition)
!>          4) Calculates and stores the supercon/superfluid gaps (uses fits,
!              there are different models implemented)
!---------------------------------------------------------------------------
subroutine solve_tov_structure(nrh, radius, i_cut, pressh, massh, enuh, elamh)

  ! Modules ----------------------------------------------------------------
  use input_params, only: use_relativistic_grid, p_central, p_cut, EoS
  use constants, only: PI, c2dg, c4dg, UNIT_R
  use grid, only: g14, get_rel_correction
  
  implicit none

  ! Subroutine arguments ---------------------------------------------------
  integer, intent(in) :: nrh
  integer, intent(out) :: i_cut
  real*8, intent(out) :: radius(0:nrh), pressh(0:nrh), massh(0:nrh), enuh(0:nrh), elamh(0:nrh)

  ! Local constants --------------------------------------------------------
  ! Maxumum value allowed for the star radius in km
  real*8, parameter :: rstar_max  = 15d0 

  real*8, parameter :: p_min = 3d22   ! Minimum pressure in cgs (very low)
                                      ! used for the first integration
  ! Local variables --------------------------------------------------------
  external derivs, derivs_tov_eqs         ! What external means?
  integer i          ! Cycling index
  integer nok, nbad  ! Used by the default subroutine odeint
  real*8 y(3)        ! Variables to be integrated by odeint
  real*8 x1, x2      ! Range of integration at each step in odeint
  real*8 h1, hh      ! Increments used in the integrations
  real*8 pc, rhoc    ! Pressure and density in geometric units
  real*8 nu0double, pcgs           ! 2*nu (in the metric) and pressure in cgs units
  real*8 rhocgs, nb, z, a, xn, xh, ymu  ! Energy density, baryon density, atomic and mass number,
                                   ! fraction of free n and p given by EoS
  real*8 schw_radius_ratio  ! Ratio R_schwarschild/R = 2*G*M/c**2*R
                            ! i.e. twice the compactness
  
  ! ------------------------------------------------------------------------

  pcgs = p_central
  write(*,*)
  write(*,"(a,es10.3)") '<info>[STRUCTURE] Central pressure [cgs]:', pcgs
  
  ! Get the central density
  call geteost(pcgs,rhocgs,nb,z,a,xn,xh,ymu) 
  pc = pcgs/c4dg      ! cm**-2
  rhoc = rhocgs/c2dg  ! cm**-2

!----------------------------------------------------
! Integrate once down to very low pressure
! to obtain precise values of M, R and nu(R)
! The variable is ln(P), so the TOV equations become:
! dr/d(ln(p)) = p/(dp/dr)
! dm/d(ln(p)) = 4*PI*rho*r**2 dr/d(ln(p)) = 4*PI*rho*r**2*p/(dp/dr)
! d(2nu)/d(ln(p)) = - 2.d0*p/(rho+p)
! y = [r, m, 2*nu - C], C constant of integration (see below)
!----------------------------------------------------

  hh = 1.d-2*pc
  x1 = dlog(pc-hh)    ! consider as first position a radius just outside the center, to avoid r=0
  x2 = dlog(p_min/c4dg)
  h1 = 1.d-2*(x2-x1)

  ! Starting values for the three variables to be integrated:
  y(1) = dsqrt((hh)/((rhoc+pc)*(rhoc/3.d0+pc)*2.d0*PI))  ! This is r(dp) as the inverse of dp(r)
                                                         ! in the next iteration, y(1), see below
  y(2) = 4.d0*PI*rhoc*y(1)**3/3.d0   ! mass assuming a ball of constant density
  y(3) = 0.d0                        ! 2*nu = y3(r) + C, choosing y3(0) = 0 implies C = 2*nu(0)

  ! Integration
  call odeint(y,3,x1,x2,h1,nok,nbad,derivs)

  ! C is found by imposing that e**2*nu(R) = sqrt(1-2M/R), where M = Y2(R) and R = Y1(R)
  ! Therefore, 2nu(0) = log(sqrt(1-2M/R)) - Y3(R)
  nu0double = - y(3) + dlog(1.d0-2.d0*y(2)/y(1))
  write(*, '(a,f8.3)') "<info>[STRUCTURE] Star gravitational Mass [Msun]: ", y(2)/1.4766d5
  write(*, '(a,f8.3)') "<info>[STRUCTURE] Star Radius [km]: ", y(1)/UNIT_R

  !----------------------------------------------------
  !  Integration of the TOV equations in geometrical units (length in cm)
  !  dP/dr  = -(rho + P)(m + 4 PI P r**3)/(r**2 - 2m*r)
  !  dm/dr  = 4PI r**2 rho
  !  dnu/dr = -(2/(rho + P))*dP/dr
  !  y = [P, m, 2*nu]  (the metric element is: e**(2nu)*(c*dt)**2 )
  !----------------------------------------------------
    
  ! Define the radius, homogeneous and up to a maximum value defined above
  ! This radius has nrh points but it is more extended than needed
  ! It is in fact used to get the profile p(r) down to p_cut
  ! Ideally, it would be better to refine r for lower p,
  ! but with enough points it doesn't matter.

  radius(0) = 0.d0
  do i = 1, nrh
    radius(i) = radius(i-1) + rstar_max/nrh
  end do
  x1 = 1.d-2*radius(1)*UNIT_R    ! The first radial point considered cannot be zero,
                                 ! but is just 1% of the first radius, in cm
  ! Initializing values of variables
  y(1) = pc - ( (rhoc+pc)*(rhoc/3.d0+pc)*2.d0*PI )*x1**2  ! in the limit of small r dp/dr propto r
                                   ! then integrate in r and the result is proportional to r**2/2
  y(2) = 4.d0*PI*rhoc*x1**3/3.d0   ! Enclosed mass in the limit of small r
  y(3) = nu0double                 ! Here the central value of 2nu comes from before
  
  ! Initializing pressure, enclosed mass, metric factors and index
  pressh = 0d0
  massh  = 0d0
  enuh   = 0d0
  elamh  = 0d0
  i = 0
  ! Integration until the pressure falls below the pressure p_cut (input),
  ! which will correspond to the last numerical point in the grid (grid.f90)
  ! The factor 0.01 is to ensure that we can also calculate the ghost cells
  do while (pcgs .ge. 0.01*p_cut)
    i = i + 1
    x2 = radius(i) * UNIT_R
    h1 = 1.d-2 * (x2 - x1) ! h1 is the increment in radius, taken as
                           ! 1% of the radius difference at each step
    ! in odeint, y is both input and output
    ! geteost is used inside derivs_tov_eqs to get the density
    call odeint(y,3,x1,x2,h1,nok,nbad,derivs_tov_eqs)
    x1 = x2
    pcgs = y(1)*c4dg       ! pressure in cgs
    pressh(i) = pcgs
    massh(i)= y(2)/1.4766d5     ! From geometric units to Solar masses (1.4766 km)
    enuh(i) = dexp(y(3)/2.d0)
    elamh(i) = 1.d0/dsqrt(1.d0-2.d0*y(2)/x2)
    if (pressh(i) .gt. p_cut) then 
      i_cut = i ! Index corresponding to the last pressh above pcut
    end if
  end do

  ! Set to one the Relativistic factors e^lam and e^nu (newtonian limit).
  if (use_relativistic_grid .eqv. .true.) then
    write(*,"(a)") "<info>[STRUCTURE] TOV solved, relativistic grid"
    ! Nothing specifically done, handled by default.
  else
    write(*,"(a)") "<info>[STRUCTURE] Newtonian grid"
    enuh = 1d0
    elamh = 1d0
  end if

  g14 = 1.d-14*y(2)/x2*(3.d10**2/x2)/enuh(i)   ! Gravity at the surface in units of 10**14 cm/s**2
  schw_radius_ratio = 1d0 - 1d0/elamh(i)**2    ! Ratio between Schwarschild and star radius: 2*G*M/c**2*R      
  call get_rel_correction(schw_radius_ratio)   ! Relativistic corrections used in magnetic field

end subroutine solve_tov_structure


!-----------------------------------------------------------
! @brief The subroutine defines the derivatives in the first integrations
!! @param[in]  x      ln(p), indipendent variable
!! @param[in]  y      sets of variables to be integrated
!! @param[out] dydx   derivatives for each y variable
!---------------------------------------------------------------
subroutine derivs(x,y,dydx)

  ! Modules ---------------------------------------------------
  use constants, only: PI, c2dg, c4dg

  implicit none
  
  ! Subroutine arguments ---------------------------------------
  real*8, intent(in) :: x, y(3)
  real*8, intent(out) :: dydx(3)

  ! Local variables --------------------------------------------
  real*8 p, rho, z, a, xn, xh, ymu
  real*8 dpdr, r, mass, nb
  real*8 pcgs, rhocgs
  ! ------------------------------------------------------------

  r = y(1)
  mass = y(2)
  p = dexp(x)
  
  pcgs = p*c4dg
  call geteost(pcgs,rhocgs,nb,z,a,xn,xh,ymu)
  rho = rhocgs/c2dg
  dpdr = - (rho+p)*(mass+4.d0*PI*r**3*p)/(r*r-2.d0*mass*r)
  dydx(1)= p/dpdr                    ! radius
  dydx(2)= 4.d0*PI*r*r*rho*(p/dpdr)  ! mass
  dydx(3) = - 2.d0*p/(rho+p)         ! 2*nu

end subroutine derivs


!-----------------------------------------------------------
! @brief The subroutine defines the derivatives in TOV equations
!-----------------------------------------------------------
!! @param[in]  r      radius, indipendent variable
!! @param[in]  y      sets of variables to be integrated
!! @param[out] dydx   derivatives for each y variable
!---------------------------------------------------------------
subroutine derivs_tov_eqs(r,y,dydx)

  ! Modules ---------------------------------------------------
  use constants, only: PI, c2dg, c4dg

  implicit none
  
  ! Subroutine arguments ---------------------------------------
  real*8, intent(in) :: r, y(3)
  real*8, intent(out) :: dydx(3)

  ! Local variables --------------------------------------------
  real*8 p, rho, z, a, xn, xh, ymu
  real*8 dpdr, mass, nb
  real*8 pcgs, rhocgs
  ! ------------------------------------------------------------

  p=y(1)
  mass=y(2)
  pcgs = p*c4dg
  ! Retrieve the density from the EoS table
  call geteost(pcgs,rhocgs,nb,z,a,xn,xh,ymu)
  rho = rhocgs/c2dg
  dpdr = -(rho+p)*(mass+4.d0*PI*r**3*p)/(r*r-2.d0*mass*r)
  dydx(1)= dpdr                     ! radius
  dydx(2)= 4.d0*PI*r*r*rho          ! mass
  dydx(3) = - 2.d0*dpdr/(rho+p)      ! 2*nu

end subroutine derivs_tov_eqs

!----------------------------------------------------------------
! @brief The subroutine retrieve the variables from the EoS table
!-----------------------------------------------------------
!! @param[in]  p       pressure
!! @param[out]  rho    mass density retrieved from the table
!! @param[out]  nb     baryon density retrieved from the table
!! @param[out]  z      atomic number of ions retrieved from the table (crust) or proton fraction (core)
!! @param[out]  a      mass number of ions retrieved from the table (crust)
!! @param[out]  xn     fraction of free neutrons (inner crust and core)
!! @param[out]  xh     fraction of nucleons clustered in ions (crust)
!! @param[out]  ymu    fraction of muons (core)
!----------------------------------------------------------------
subroutine geteost(p,rho,nb,z,a,xn,xh,ymu)

  ! Modules ----------------------------------------------------
  ! None
  use input_params, only: EoS

  implicit none

  ! Subroutine arguments ---------------------------------------
  real*8, intent(in) :: p
  real*8, intent(out) :: rho, nb, z, a, xn, xh, ymu

  ! Local variables --------------------------------------------
  integer, parameter :: nmax=2000
  real*8, dimension(nmax), save :: pt, rhot, nbt, zt, at, xnt, ymut
  ! Variables to store and saved from the table values
  integer, save :: neos  ! Number of lines read in the table
  real*8, dimension(nmax) :: p_temp  ! Array to locate the closest value of pressure
  real*8 alpha           ! Dummy variable used in the interpolation
  integer i              ! Loop index

  ! Pre-conditions ---------------------------------------------
  ! If the pressure is too high or too low, stop
  if( (p<pt(2)) .or. (p>pt(neos-1)) ) then
    write(*,'(a,es12.4,a,2es12.4)') '<ERROR>[ns_structure.f90] The required pressure ', p, &
    & ' is out of EoS table range: ',pt(2),pt(neos-1)
    stop
  end if
  !----------------------------------------------------------------

  ! Locate the input value of pressure in the array,
  ! using the MINLOC intrinsic function that returns the location of the
  ! minimum value in an array
  p_temp = dabs(pt - p)
  i = MINLOC(p_temp,DIM=1)

  ! Linear interpolation in log10(p) for energy and baryon densities
  ! p = k*nb**adiabatic_index, approximating that rho is almost proportional to nb
  ! alpha = 10.d0**(log10(p/pt(i))/gammat(i))
  ! rho = rhot(i)*alpha
  ! nb = nbt(i)*alpha

  ! Interpolation of xn
  ! Order 0: closest value of xn in the table:
  ! xn = xnt(i)

  ! Order 1: linear interpolation
  ! alpha = (rho-rhot(i))/(rhot(i+1)-rhot(i-1))
  ! xn = xnt(i) + alpha*(xnt(i+1)-xnt(i-1))

  ! Order 2: parabolic interpolation of rho, nb and xn:
  ! It is the best option to avoid artificial steps in the profiles
  alpha = (p-pt(i)) * (p-pt(i+1)) / ( (pt(i-1)-pt(i)) * (pt(i-1)-pt(i+1)) )
  rho = rhot(i-1)*alpha
  nb = nbt(i-1)*alpha
  xn = xnt(i-1)*alpha
  ymu = ymut(i-1)*alpha
  alpha = (p-pt(i+1))*(p-pt(i-1))/((pt(i)-pt(i+1))*(pt(i)-pt(i-1)))
  rho = rho + rhot(i)*alpha
  nb = nb + nbt(i)*alpha
  xn = xn + xnt(i)*alpha
  ymu = ymu + ymut(i)*alpha
  alpha = (p-pt(i-1))*(p-pt(i))/((pt(i+1)-pt(i-1))*(pt(i+1)-pt(i)))
  rho = rho + rhot(i+1)*alpha
  nb = nb + nbt(i+1)*alpha
  xn = xn + xnt(i+1)*alpha
  ymu = ymu + ymut(i+1)*alpha
  
  if ((xn >= 1.d0) .or. (xn < 0.d0)) then
  ! Put by hand limits on the allowed fraction of free neutrons
    xn = xnt(i)
    xn = dmax1(xn,0.d0)
  end if

  if ((ymu >= 1.d0) .or. (ymu < 0.d0)) then
  ! Put by hand limits on the allowed fraction of muons
      ymu = ymut(i)
      ymu = dmax1(ymu,0.d0)
  end if
  
  ! In the crust, with ions + free neutrons + electrons
  if (at(i)>0.d0) then
    xh = 1.d0 - xn
    z = zt(i) 
    a = at(i)
    ymu = 0.d0
  ! In the core, with protons + neutrons + electrons + muons
  ! xn + yp = 1
  else
    xh = 0.d0
    z = 0.d0
    a = 1.d0
  endif

  return

  ! If init_eos_tab is called, it enters here to read and save the table:
  entry init_eos_tab()

  ! Read the EOS table
  ! EOS_DH.tab is a manually copied table from Tables of Douchin & Heansel 2001, A&A)
  select case (EoS)
  !  case ("DH")
  !    open(10,file='in/EOS_DH.tab')
  !  case ("BSk21")
  !    open(10,file='in/EOS_BSk21.tab')
    case ("BSk22")
      open(10,file='in/EOS_BSk22.tab')
    case ("BSk24")
      open(10,file='in/EOS_BSk24.tab')
    case ("BSk25")
      open(10,file='in/EOS_BSk25.tab')
    case ("BSk26")
      open(10,file='in/EOS_BSk26.tab')
  !  case ("DH_1_u")
  !    open(10,file='in/EOS_DH_1_u.tab')
    case ("SLy4")
      open(10,file='in/Compose/EOS_RG(SLy4).tab')
  !  case ("CMF6")
  !    open(10,file='in/Compose/EOS_DS(CMF)-6.tab')
  !  case ("SLy2")
  !    open(10,file='in/Compose/EOS_RG(SLy2).tab')
  !  case ("SKa")
  !    open(10,file='in/Compose/EOS_RG(SKa).tab')
  !  case ("CMF2")
  !    open(10,file='in/Compose/EOS_DS(CMF)-2.tab')
  !  case ("SKb")
  !    open(10,file='in/Compose/EOS_RG(SKb).tab')
  !  case ("SkMp")
  !    open(10,file='in/Compose/EOS_RG(SkMp).tab')
  !  case ("TNTYST")
  !    open(10,file='in/Compose/EOS_TNTYST.tab')
    case ("SLy4mu")
      open(10,file='in/Compose/EOS_RG(SLy4)mu.tab')
    case ("SKamu")
      open(10,file='in/Compose/EOS_RG(SKa)mu.tab')
    case ("SkMpmu")
      open(10,file='in/Compose/EOS_RG(SkMp)mu.tab')
    case ("SLy2mu")
      open(10,file='in/Compose/EOS_RG(SLy2)mu.tab')
    case ("SKbmu")
      open(10,file='in/Compose/EOS_RG(SKb)mu.tab')
    case ("CMF2mu")
      open(10,file='in/Compose/EOS_DS(CMF)-2mu.tab')
    case ("CMF6mu")
      open(10,file='in/Compose/EOS_DS(CMF)-6mu.tab')
    case default
      write(*,*) "<ERROR>[ns_structure]: Invalid EoS name ", EoS
      stop
  end select
  ! Reading: baryon density [fm**-3], energy density [g/cm**3], pressure [erg/cm**3],
  !          ion atomic number, ion mass number, free neutrons fraction, fraction of muons

  do i = 1, nmax
    read(10,*,end=20) nbt(i), rhot(i), pt(i), zt(i), at(i), xnt(i), ymut(i)
  enddo
  20 neos = i - 1
  close(10)
  return

end subroutine geteost


    ! --------------------------------------------------------------------------
    !> Output struture information to files and screen.
    !> @brief Outputs structure grid and gaps information to files and also
    !>        useful information (surface density, core radius...) to screen.
    !---------------------------------------------------------------------------
subroutine output_structure()

  ! Modules ----------------------------------------------------------------
  use utils, only: get_free_unit
  use grid
  use input_params
  use constants, only: T8_TO_MEV, GRAV_CONSTANT, MSUN, CLIGHT, UNIT_R

  implicit none

  ! Subroutine arguments ---------------------------------------------------
  ! None.

  ! Local constants --------------------------------------------------------
  character(len=35), parameter :: STRUCTURE_FORMAT = "(f8.4,3es11.3,2f8.4)"
  character(len=35), parameter :: COMPOSITION_FORMAT = "(f8.4,2es11.3,4f7.3,3f8.3,7es11.3)"
  character(len=14), parameter :: GAPS_FORMAT = "(8es11.3)"

  ! Local variables --------------------------------------------------------
  integer :: j, l
  integer :: unit
  ! Dimensionless parameter GM/c^2R
  real*8 compactness

  ! ------------------------------------------------------------------------

  compactness = GRAV_CONSTANT*MSUN/(CLIGHT**2*UNIT_R)*mass(nr)/r(nr)

  unit = get_free_unit()

  open(unit, file="out/1D/structure.d")
  write(unit, '(2a)') "r[km], en.density[g/cm**3], ", &
 &                   "mass[Msun], p[dyn/cm**2], e**nu, e**lambda"
  do j = 1, ncore+nr
    write(unit, STRUCTURE_FORMAT) &
      & rtot(j), rho(j), mass(j), press(j), enu_tot(j), elambda_tot(j)
  end do ! j
  close(unit)

  open(unit, file="out/1D/composition.d")
  write(unit, '(5a)') "r[km], en.density[g/cm**3], n_b[fm**(-3)], ", &
 &                   "Ye, Yp, Yn, Ymu, Xh, A, Z, ", &
 &                   "kFe[1/fm], kFn[1/fm], kFp[1/fm], kFmu[1/fm] ", &
 &                   "rel.eff.me, rel.eff.mn, rel.eff.mp"
  do j = 1, ncore + nr
    write(unit, COMPOSITION_FORMAT) &
      & rtot(j), rho(j), nb(j), ye(j), yp(j), yn(j), ymu(j), xh(j), &
      & aa(j), zz(j), kFe(j), kFn(j), kFp(j), kFmu(j), &
      & effme(j), effmn(j), effmp(j)
  end do ! j
  close(unit)

  open(unit, file="out/1D/sf_gaps.d")
  write(unit, '(3a)') "en.density[g/cm**3], nb[fm**(-3)], ", &
&                     "Tc_n_crust[MeV], Tc_n_core[MeV], Tc_p_core[MeV]" , &
&                     "gap_n_crust[MeV], gap_n_core[MeV], gap_p_core[MeV]"
  do j = 1, ncore + nr
    write(unit, GAPS_FORMAT) &
      & rho(j), nb(j), T8_TO_MEV*tccru(j), T8_TO_MEV*tcn(j), T8_TO_MEV*tcp(j), &
      & gapn_crust(j), gapn_core(j), gapp_core(j)
  end do
  close(unit)

  open(unit, file="out/1D/fh.yg")
  write(unit, '(3a)') "Hall_prefactor[km^2/(Myr*10^12 G)]"
  do j = 1, nr
    write(unit, COMPOSITION_FORMAT) r(j), fh(j)
  end do
  close(unit)

  ! Screen output.
  write(*, "(a,1pe11.3)") "<info>[STRUCTURE] Crust/envelope density [g/cm^3]:", rho(ncore+nr)
  write(*, "(a,f8.3)")    "<info>[STRUCTURE] Compactness GM/Rc^2 :", compactness
  write(*, "(a,f8.3)")    "<info>[STRUCTURE] Surface gravity [10^14 cm/s^2] :", g14
  write(*, "(a,f8.3)")    "<info>[STRUCTURE] Surface redshift (1+z) :", 1.d0/enu(nr)
  write(*, "(a,1pe11.3)") "<info>[STRUCTURE] Moment of Inertia [g*cm^2]:", moment_inertia
  write(*, *)

end subroutine output_structure


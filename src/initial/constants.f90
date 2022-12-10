!> @brief Physical Constants of common use in the code
!! 
!!------------------------------------------------------------------------------ 
module constants

  implicit none

  real*8, parameter :: PI = dacos(-1.d0)
  real*8, parameter :: MASS_P = 1.6726219d-24         ! proton mass [g]
  real*8, parameter :: MASS_N = 1.674927471d-24       ! neutron mass [g]
  real*8, parameter :: MASS_N_MEV = 939.565560d0      ! neutron mass [MeV]
  real*8, parameter :: MASS_E_MEV = 0.510998d0        ! electron mass [MeV]
  real*8, parameter :: E_CHARGE = 1.602e-19           ! electron charge [C]
  real*8, parameter :: STATC_TO_C = 10                ! conversion of charge (e/c) from statC to cgs
  real*8, parameter :: E2 = 2.3070771d-19             ! elementary charge squared [erg*cm]
  real*8, parameter :: CLIGHT = 2.99792458d10         ! speed of light [cm/s]
  real*8, parameter :: OMEGA_P_12 = 9.6043d15         ! proton Larmor angular frequency at 10**12 G [rad / s], omega_p = eB/m
  real*8, parameter :: RHO_NUC = 2.82d14              ! nuclear saturation energy density [g cm^-3]
  real*8, parameter :: N0_NUC = 0.16d0                ! nuclear saturation baryon density [fm^-3]
  real*8, parameter :: MSUN = 1.98847d30              ! Solar mass [kg]
  real*8, parameter :: GRAV_CONSTANT = 6.674d-5       ! Gravitational constant [cm^3/(kg*s^2)]
  real*8, parameter :: HBARC = 197.3269               ! hbar*c in [MeV*fm]
  real*8, parameter :: T_YEAR = 3.1556926d7           ! conversion years to seconds
  real*8, parameter :: T8_TO_MEV = 8.617328d-3        ! Conversion from 10**8 K to MeV [MeV/10**8 K]
  real*8, parameter :: RHO_TO_N = MASS_N*1.d39        ! conversion from [gr/cm^3] to [particles/fm^3]
  real*8, parameter :: K_BOLTZMANN = 8.617d-5         ! Boltzman constant [eV/K]
  real*8, parameter :: STEFAN_BOLTZMANN = 5.670373d-5 ! black-body constant [erg/(cm^2*s*K^4)]
  real*8, parameter :: PHFLUX_CONSTANT = 1.5205163d11 ! constant for photon flux [ph/(cm^2*s)]
  real*8, parameter :: ALPHA = 1d0/137d0             ! fine structure constant
  ! Integrated photon flux formula:
  ! N = 4pi*Apery_const*k_b^3/(h^3*c^2) T^3 =
  !   = 1.5205163e11 T[K]^3 photons cm^-2 s^-1 =
  !   = 1.5205163e21 T[K]^3 photons km^-2 s^-1 =
  ! Apery's constant = 1.202057  (Riemann zeta function with argument 3)

  ! Conversion factors to geometrized units (G=c=1)
  real*8, parameter :: c2dg=1.347459039d028          !    c2/G=1.347459039d28 g cm-1
  real*8, parameter :: c4dg=1.211035789d49           !    c4/G=1.211035789d49 g cm s-2

  ! Units for the numerical values used in the code
  real*8, parameter :: UNIT_T = 1.d8                 ! Units of temperature [K]
  real*8, parameter :: UNIT_B = 1.d12                ! Units of magnetic field [G]
  real*8, parameter :: UNIT_R = 1.d5                 ! Units of lengths [cm]
  real*8, parameter :: UNIT_E = UNIT_B*UNIT_R/(T_YEAR*CLIGHT) ! Units of electric field [statVolt = G]
  real*8, parameter :: UNIT_EN = 1.d40              ! Units of the energy [erg/s] used in the thermal part
  real*8, parameter :: UNIT_JOULE = UNIT_B**2*UNIT_R**3/(UNIT_EN*4d0*PI*T_YEAR)  ! Joule heating density units [10^40 erg/(km^3*s)]
  real*8, parameter :: UNIT_TIME = T_YEAR*1.d6         ! Units of time [s]

!-----------------------------------------------------------------------
! UNITS FOR MAGNETIC ANALYSIS AND JOULE HEATING:
!-----------------------------------------------------------------------
! Numerically we use:
!
! Jnum = e^nu 4pi J / c
! Q_joule = - (4*pi/c^2)*eta*J^2*e^(2nu)
!     = - eta*Jnum^2/(4*pi) = - UNIT_JOULE*eta*Jnum^2/(4*pi)
! with UNIT_JOULE = UNIT_J^2*UNIT_ETA*UNIT_R**3/(4*pi) / [10^40] = 2.522e(-16)
!  
! UNIT_J   [1e12 G/km]
! UNIT_ETA [km^2/Myr]
! UNIT_R^3 [km^3] to pass from G^2/s = erg/cm^3/s to erg/km^3/s
! UNIT_JOULE [1e40 erg/km^3 s]
!
!-----------------------------------------------------------------------
  
end module constants

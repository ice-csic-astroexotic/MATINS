!-------------------------------------------------------------------------------
! Magneto Thermal 2D
!-------------------------------------------------------------------------------
! Main program
!
!> @author
!> Daniele Viganò
!> José A. Pons
!> 
!-----------------------------------------------------------------------

program mt2d

  ! Modules ----------------------------------------------------------------
  use input_params
  use input_params, only: angular_dimension, radial_dimension
  use input_params, only: read_input_file
  use input_params, only: enable_breaking
  use input_params, only: resume_run
  use initial_magnetic, only: binit
  use resume, only: read_resume
  use grid, only: allocate_grid
  use grid, only: set_grid_size, build_structure, build_grid
  use grid, only: kmax, lmax, benu, spindown_prefactor
  use grid, only: tem, tem0, q_joule, q_joule_shock, q_joule_average
  use grid, only: sfluxb, tss
  use microphysics, only: compute_transport
  use magnetic_evolution, only: magnetic_evol, compute_joule
  use magnetic_stress, only: breaking
  use magnetic_analysis, only: analyse_magnetic_field
  use thermal_evol, only: tevol
  use legpol, only: allocate_legendre_polynomials
  use legpol, only: compute_legendre_polynomials, getbl_init, blout
  use output, only: output_temperature, output_magnetic, output_resume
  
  implicit none

  ! Local constants --------------------------------------------------------
  real*8, parameter :: eps_dt_print = 1.d-5    ! Tolerance for printing time [years]

  ! Local variables --------------------------------------------------------
  integer :: l
  integer :: itert, iterb
  real*8 :: dt, dtb, tyear_print, tyear_b_print, tyear, tyear_b
  real*8 :: pdot, bindex
  real*8 :: bpdip, bpdip_old

  ! Read input and initialize values
  call read_input_file()
  ! Initialize magnetic field normalization (bpdip).

  ! Set the size of the angular and radial dimensions of the grid read from
  ! input so that the grid can be allocated with the proper dimensions.
  call set_grid_size(angular_dimension, radial_dimension)

  ! Allocate dynamic arrays throughout the whole code.
  call allocate_grid()
  call allocate_legendre_polynomials()

  ! Generate the radial grid and compute the star structure
  call build_structure

  ! Generate the angular grid and geometry arrays (length, surface, volume).
  call build_grid

  ! Initialize the Legendre polynomials.
  call compute_legendre_polynomials()

  ! Initialize GETBL.
  ! Needed for the normalization and for the boundary conditions.
  call getbl_init

  ! Resume checkpoint if requested.
  if (resume_run .eqv. .true.) then
    call read_resume(tyear, per, pdot, bpdip)
    tyear_b=tyear
  else

    ! Initialize the magnetic field.
    call binit
    tyear=0d0
    tyear_b=0d0
    tem=tinit
    tss=tinit
    bpdip = bpolin
    pdot=spindown_prefactor*bpdip**2/per
  endif

  ! Initialize time, redshifted temperature (tem=T*e^nu), fluxes and iteration numbers
  itert=0
  iterb=0
  tyear_print=0d0
  sfluxb=0d0


  ! Radiative conductivities (usually not used):
  !	call yakitl

  write(*,*) "Starting evolution, up to years: ",final_time
  write(*,*)

  ! While loop is used because it allows to set both the physical total time
  ! regardless of the fixed/adaptive timestep used
  do while ( tyear <= final_time )

    if (itevol .eqv. .false.) then
      dt = 1.01*final_time  ! If T does not change, give only one dt step
    else
      ! Option 1: adaptive cooling timestep, increasing at 10,100,1000 years
      call adaptive_cooling_timestep(tyear,dt)
      ! Option 2: constant dt to be set here
      ! dt = 1d0
    endif

    ! Recover the physical temperature.
    do l=1,lmax
      tem0(1:kmax,l) = tem(1:kmax,l)/benu(2*l-1)
    enddo

    ! Shear modulus and maximum strength
    call compute_shear

    ! Microphysics: diffusion coefficients, anisotropy factors, magnetic diffusivity
    call compute_transport

    ! Microphysics: specific heat, neutrino emissivities (in 1e40 erg/s)
    call compute_heat_capacity
    call compute_neutrino_emissivity

    ! Surface boundary conditions for heat transport
    call fluxesb(bpdip)

    if (itevol .eqv. .false. .or. &
        & tyear == 0d0 .or. &
        & (tyear_print + eps_dt_print) >= cooling_output_dt) then
      call output_temperature(tyear, bindex, bpdip, itert, per, pdot)
      tyear_print=0d0
    endif

    ! Initialize the Joule heating
    q_joule_average = 0d0

    ! No magnetic field evolution:
    if(ibevol .eqv. .false.)then

      call period_evol(dt,bpdip,bpdip,per,pdot,bindex)

    ! Magnetic field evolution:
    else

      do while ( tyear_b <= tyear + dt )

        ! IT IS IMPORTANT TO MAINTAIN THIS ORDER OF CALLS TO PRESERVE CONSISTENCY IN OUTPUT AND ANALYSIS
        ! Compute the Joule power for the fields entering the loop, used in the cooling advance
        call compute_joule

        ! Consider the average over dt of q_joule (important if dtb<<dt)
        q_joule_average = q_joule_average + (q_joule+q_joule_shock)*dtb/dt

        ! Analyse fields. Some quantities (integrals and rates) regard the 
        ! time variation dtb used in the previous iteration
        ! (which led fields to the values entering the loop)
        ! In the case of time=0 (first loop), these quantities are ignored
        call analyse_magnetic_field(tyear_b,dtb)

        ! Activate magnetic stress if needed.
        if (enable_breaking) then
          call breaking(dtb,tyear_b)
        end if

        ! Option 1: adaptive timestep (set courant_prefactor > 0)
        ! Option 2: fixed timestep from input (set courant_prefactor=0 in input)
        if (courant_prefactor > 0) then 
          call adaptive_magnetic_timestep(tyear,tyear_b,tyear_b_print,dt,dtb)
        else
          dtb = dtb0
        endif

        ! Write output of current fields before evolving them
        ! eps_dt_print is the tolerance given
        if ( (tyear_b_print + eps_dt_print) >= magnetic_output_dt .or. tyear_b == 0d0 ) then
          call output_magnetic(iterb,tyear_b)
          tyear_b_print = 0.d0
        endif

        ! Magnetic field evolution.
        ! dtb in input and in main is in YEARS but in magnetic_evol in Myrs.
          call magnetic_evol(1.d-6*dtb,iterb)

        ! Dipolar surface field updated
        bpdip_old = bpdip
        bpdip = - 2.d0*blout(1)
        call period_evol(dtb,bpdip,bpdip_old,per,pdot,bindex)
        if(tyear.eq.0d0) bindex=0d0

        tyear_b_print = tyear_b_print + dtb
        tyear_b = tyear_b + dtb
        iterb = iterb + 1

      enddo

    endif

    if (itevol .eqv. .false.) then
    ! No temperature evolution case
      tem = tinit
    else
    ! Advance temperature one timestep dt
      call tevol(dt)
    endif

    ! Update times.
    itert = itert + 1
    tyear = tyear + dt 
    tyear_print = tyear_print + dt

    ! End of time loop
    enddo

  stop

end program

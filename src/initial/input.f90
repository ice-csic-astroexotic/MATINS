module input_params

  ! General parameters.

  ! Total time [yr].
  real*8, save :: final_time
  ! Cooling and magnetic output dt [yr].
  real*8, save :: cooling_output_dt,magnetic_output_dt
  ! Grid dimensions kmax and lmax.
  integer, save :: angular_dimension, radial_dimension

  ! Magnetic parameters.

  ! Initial poloidal and toroidal amplitudes [10^12 G].
  real*8, save :: bpolin, btorin
  ! Initial magnetic field topology.
  integer, save :: bgeom
  ! Enable magnetic field evolution.
  logical, save :: ibevol
  ! Minimum magnetic timestep [yr].
  real*8, save :: dtb0
  ! Courant pre-factor.
  real*8, save :: courant_prefactor
  ! Intial spin period.
  real*8, save :: per
  ! Pre-coefficient of toroidal field hyper-resistivity.
  real*8, save :: coeff_hyper
  ! Turn on/off Hall effect and/or ambipolar diffusion completely. 
  logical, save :: enable_hall_effect

  ! Temperature parameters.

  ! Initial temperature [10^8 K].
  real*8, save :: tinit
  ! Evolution of temperature.
  logical, save :: itevol
  ! Impurity parameters in the crust and pasta phase.
  real*8, save :: Qimp, Qpasta
  ! Envelope model.
  real*8, save :: ienv
  ! Gap models.  Set any parameter to 0 to deactivate that gap.
  !   SF n in the crust: (1-7)
  !   SF n in the core: (11-16)
  !   SF p in the core:  (8-10)
  integer, save :: superfluid_n_crust
  integer, save :: superfluid_n_core
  integer, save :: superfluid_p_core

  ! Advanced parameters.

  ! Flag to Specify if we are in the Newtonian (false) or relativistic (true)
  ! case for the grid. In the Newtonian case, all metric factors are one.
  logical, save :: use_relativistic_grid
  ! Magnetic advance method.
  character(len=4), save :: magnetic_advance_method
  ! For the Etor, you can choose either the centered difference, or the upwind.
  character(len=6), save :: etor_scheme
  ! Activate burgers correction in magnetic evolution.
  logical, save :: enable_burgers_correction
  ! Number of initial multipoles for phi and psi.
  integer, save :: n_initial_multipoles_phi
  integer, save :: n_initial_multipoles_psi
  ! Arrays of initial multipole value for phi and psi.
  real*8, dimension(:), allocatable, save :: initial_multipoles_phi
  real*8, dimension(:), allocatable, save :: initial_multipoles_psi
  ! Enable calls to breaking subroutine (or not, to save time).
  logical, save :: enable_breaking
  ! Resume run from file or not.
  logical, save :: resume_run

  contains

    !> @brief This Subroutine reads different flags and parameters from a file.
    subroutine read_input_file()

      implicit none
  
      open(1,file="in/input.dat")
      ! Simulation time and frequencies of snapshots.
      read(1,*)
      read(1,*) final_time
      ! Resume run.
      read(1,*) resume_run
      ! Checkpointing times.
      read(1,*) cooling_output_dt
      read(1,*) magnetic_output_dt
      ! Grid dimensions.
      read(1,*) angular_dimension
      read(1,*) radial_dimension
      read(1,*)
      ! Magnetic block.
      read(1,*) ibevol
      read(1,*) bpolin
      read(1,*) btorin
      read(1,*) bgeom
      read(1,*) n_initial_multipoles_phi
      allocate(initial_multipoles_phi(n_initial_multipoles_phi))
      read(1,*) initial_multipoles_phi
      read(1,*) n_initial_multipoles_psi
      allocate(initial_multipoles_psi(n_initial_multipoles_psi))
      read(1,*) initial_multipoles_psi
      read(1,*) enable_breaking
      ! Temperature block.
      read(1,*)
      read(1,*) itevol
      read(1,*) ienv
      read(1,*) Qimp
      read(1,*) Qpasta
      ! Advanced parameters.
      read(1,*)
      read(1,*) use_relativistic_grid
      read(1,*) enable_hall_effect
      read(1,*) tinit
      read(1,*) superfluid_n_crust
      read(1,*) superfluid_n_core
      read(1,*) superfluid_p_core
      read(1,*) per
      read(1,*) 
      ! Numerical method choices for magnetic field.
      read(1,*) magnetic_advance_method
      read(1,*) courant_prefactor
      read(1,*) dtb0
      read(1,*) etor_scheme
      read(1,*) enable_burgers_correction
      read(1,*) coeff_hyper

      close(1)

    end subroutine read_input_file

end module input_params

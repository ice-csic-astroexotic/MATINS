module input_params

  ! General parameters.

  ! Total time [yr].
  real*8, save :: final_time
  ! Cooling and magnetic output dt [yr].
  real*8, save :: tau_output,tau_output3D
  ! Grid dimensions nangt, nrt
  integer, save :: angular_dimension, radial_dimension, radial_dimension_core

  ! Magnetic parameters.
  integer, parameter :: LMAX_IN = 10
  ! Initial poloidal and toroidal normalizations [10^12 G].
  real*8, save :: bpol_init, btor_init
  ! Initial magnetic field topology.
  real*8, dimension(1:LMAX_IN,-LMAX_IN:LMAX_IN), save :: pol_lm, tor_lm
  ! Enable magnetic field evolution.
  logical, save :: ibevol
  ! Time advance method for the magnetic evolution.
  character(len=3), save :: time_advance
  ! Minimum magnetic timestep [yr].
  real*8, save :: dtb0
  ! Courant pre-factor.
  real*8, save :: courant_prefactor

  ! Structure parameters.
  real*8, save :: p_central, p_cut
  character(len=6), save :: EoS
  character(len=7), save :: profile
  ! Temperature parameters.

  ! Initial temperature [10^8 K]
  real*8, save :: T_init
  ! Minimum, maximum timestep for cooling [yr], and parameter used for the progressive timestep.
  real*8, save :: max_dt_cooling, min_dt_cooling, eps_dt_cooling
  ! Evolution of temperature.
  logical, save :: itevol
  ! Impurity parameters in the crust and pasta phase.
  real*8, save :: Qimp, Qpasta
  ! Envelope model.
  character(len=10), save :: envelope
  ! Gap models.  Set any parameter to '0' to deactivate that gap.
  character(len=6), save :: superfluid_n_crust
  character(len=6), save :: superfluid_n_core
  character(len=6), save :: superfluid_p_core

  ! Advanced parameters.

  ! Flag to Specify if we are in the Newtonian (false) or relativistic (true)
  ! case for the grid. In the Newtonian case, all metric factors are one.
  logical, save :: use_relativistic_grid
  ! Magnetic advance method.
  character(len=4), save :: magnetic_advance_method
  ! For the Etor, you can choose either the centered difference, or the upwind.
  character(len=6), save :: e_scheme
  integer, save :: checkpnt_freq, resume_checkpnt_number
  integer, save :: num_threads

  contains

    !> @brief This Subroutine reads different flags and parameters from a file.
    subroutine read_input_file()

      implicit none

      character :: dummy
      integer l

      open(1,file="in/input.dat")
      ! Simulation time and frequencies of snapshots.
      read(1,*)
      read(1,*) final_time
      ! Resume run options.
      read(1,*) checkpnt_freq
      read(1,*) resume_checkpnt_number
      ! Frequency of 1D and 2/3D output
      read(1,*) tau_output
      read(1,*) tau_output3D
      ! Grid dimensions.
      read(1,*) angular_dimension
      read(1,*) radial_dimension
      read(1,*) radial_dimension_core
      read(1,*)
      ! Magnetic block.
      read(1,*) ibevol
      read(1,*) bpol_init
      read(1,*) btor_init
      ! Temperature block
      read(1,*)
      read(1,*) itevol
      read(1,*) envelope
      read(1,*) T_init
      ! Structure block
      read(1,*)
      read(1,*) EoS
      read(1,*) p_central
      read(1,*) p_cut
      read(1,*) profile
      read(1,*) use_relativistic_grid
      read(1,*) Qimp
      read(1,*) Qpasta
      ! Superfluidity block
      read(1,*)
      read(1,*) superfluid_n_crust
      read(1,*) superfluid_n_core
      read(1,*) superfluid_p_core
      ! Numerical method choices for magnetic field.
      read(1,*) 
      read(1,*) time_advance
      read(1,*) courant_prefactor
      read(1,*) dtb0
      read(1,*) e_scheme
      ! Cooling timestep.
      read(1,*) 
      read(1,*) min_dt_cooling
      read(1,*) max_dt_cooling
      read(1,*) eps_dt_cooling

      ! Magnetic topology: initial poloidal and toroidal multipoles
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      pol_lm = 0d0
      do l=1,LMAX_IN
        read(1,*) dummy, pol_lm(l,-l:l)
      end do
      ! Poloidal field
      read(1,*)
      read(1,*)
      tor_lm = 0d0
      do l=1,LMAX_IN
        read(1,*) dummy, tor_lm(l,-l:l)
      end do

      close(1)


      if (mod(angular_dimension,2) == 0.) then
        write(6,*) "<ERROR>[input]: the angular dimension has to be an odd number, >=3."
        stop
      endif        

      if (EoS .eq. "simple" .and. profile .eq. "realist") then
        write(6,*) "<ERROR>[input]: You need a realistic EoS to have realistic profiles."
        stop
      endif

    end subroutine read_input_file

end module input_params

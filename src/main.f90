program main

  use input_params
  use grid, only: allocate_grid, build_grid
  use grid, only: nrt, tem0, temp, enu
  use grid, only: nangt, T_core !for debugging
  use microphysics, only: analytical_microphysics
  use conductivities, only: compute_conductivities
  use magnetic_evolution, only: magnetic_evol, dtb_adaptive
  use thermal_evolution, only: tevol
  use output, only: output_magnetic_1D, output_magnetic_2D, output_surface
  use output, only: output_thermal_1D, output_vtu, output_temperature_cooling_curve
  use OMP_LIB

  implicit none

  real*8 :: dt, dtb
  real*8 :: time, time_b
  real*8 :: tau_output_serv, tau_output3D_serv
  integer :: nradial, nangular
  integer :: ib, it, ic, nout, nout3D, noutT, nout3DT

  !Read input and initialize values
  call read_input_file()

  time = 0.d0 !1.0 this is for Perez-Azorin test!3.d-1
  time_b = 0.d0

  ! Definition of the grid, metric and all common quantities
  call allocate_grid()

  ! Generate all the grid and metric
  ! The structure of the star is called within it
  call build_grid()
  print*,"Grid built"

  ! Definition of initial magnetic field and temperature in spherical coordinates
  call initial_condition()
  print*,"Initial conditions set"

  ib = 0
  it = 0
  if (courant_prefactor .le. 0.) then
    dtb = dtb0
  endif


 ! Thermal evolution
  do while (time < final_time)

    if (time == 0.d0) then
      do ic = 1,nrt
        tem0(ic,:,:,:) = temp(ic,:,:,:)/enu(2*ic-1)
      enddo 
    endif

     if (itevol .eqv. .false.) then
    !  dt = 1.d-5*final_time 
      dt = 1.01d0*final_time 
    ! Analytical evolution of temperature (from Yakovlev, see Dehman 2022)
      ! if (time /= 0.d0) then
      !   do ic = 1,nrt
      !     tem0(ic,:,:,:) = 6.8695d0/((time*1.d6)**(1.d0/6.d0)*enu(2*ic-1)) 
      !   enddo
      ! end if 
    else
      call adaptive_cooling_timestep(time,dt)
    endif


    if (profile .eq. "realist") then
      ! Conductivities: diffusion coefficients, anisotropy factors, magnetic diffusivity
      ! Heat capacity, neutrino emissivities (in 1e40 erg/s)
      call compute_conductivities
      call compute_heat_capacity
      call compute_neutrino_emissivity

      ! Surface boundary conditions for heat transport
      call envelope_model
    
    ! Call microphysics analytical formula if they are prescribed analytically
    else
      call analytical_microphysics(profile)
      call envelope_model
    endif

    if (itevol .eqv. .true.) then

      ! TBD: choose a synchronized way for the output, and less frequent
      tau_output_serv = min(10*dt, tau_output)
      tau_output3D_serv = min(10*dt,tau_output3D)

      noutT = floor(tau_output_serv/dt)
      nout3DT = floor(tau_output3D_serv/dt)
      if (noutT == 0) noutT = 1
      if (nout3DT == 0) nout3DT = 1

      if (modulo(it,noutT) == 0) then
        call output_thermal_1D(time*1.d6)
        call output_temperature_cooling_curve(it, time, dt)
      end if

      if (modulo(it,nout3DT) == 0) then
        call output_vtu(time,tau_output_serv)
        call output_surface(time)
      end if

    endif

    if (ibevol .eqv. .true.) then
   
      if (courant_prefactor .gt. 0.) then
        call dtb_adaptive(dtb)
        dtb = maxval([dtb*courant_prefactor,dtb0])
      endif
 
      nout = floor(tau_output/dtb)
      nout3D = floor(tau_output3D/dtb)
      if (nout == 0) nout = 1
      if (nout3D == 0) nout3D = 1

      do while (time_b <= time + dt)

        if (courant_prefactor .gt. 0.) then
          call dtb_adaptive(dtb)
          dtb = maxval([dtb*courant_prefactor,dtb0])
        endif

        if (modulo(ib,nout) == 0) then
          call output_magnetic_1D(ib,time_b)
        end if

        if (modulo(ib,nout3D) == 0) then
          call output_magnetic_2D(time_b)
          call output_vtu(time_b, dt)
        end if
    
        call magnetic_evol(dtb, time_b)
    
        time_b = time_b + dtb
        ib = ib + 1
    
      end do 

    else ! In case of no B evolution, print the initial conditions

      if (time == 0.) then
        call output_magnetic_1D(ib,time_b)
        call output_magnetic_2D(time_b)
      end if

    endif
  
    time = time + dt
    it = it + 1

    if (itevol .eqv. .true.) then
      call tevol(dt)
    endif
  
  enddo

  print*
  print*,"DONE"

  stop

end program 

!-----------------------------------------------------------------------
! 2D magneto-thermal evolution of neutron stars.
! Based on v3.7.01 (22/07/2013), written by Daniele Vigano (DV).
!
! TBD: bm, bmed, eta, etab are duplicated
! Contents:
!-----------------------------------------------------------------------
      program mt2d

      use input_params
      use input_params, only: angular_dimension, radial_dimension
      use input_params, only: read_input_file
      use grid, only: kmax, lmax, nang, np, fh, benu, belam, jcore, g14, spindown_prefactor
      use grid, only: allocate_grid
      use grid, only: set_grid_size, build_structure, build_grid
      use grid, only: tem0, q_neutrino, q_joule, q_joule_shock, q_joule_average, rmc, tdc, emnu
      use grid, only: br, bth, bphi, aphi, bm
      use grid, only: bmed, j2
      use grid, only: cfluxb, sfluxb, tss
      use grid, only: eta, omegatau
      use microphysics, only : difdrv, etab
<<<<<<< HEAD:src/main.f
      use magnetic_evolution, only: magnetic_evol
      use magnetic_stress, only: breaking
=======
      use magnetic_evolution, only: magnetic_evol, q_joule_calc
      use magnetic_analysis, only: calculate_magnetic_analysis
>>>>>>> 49bee1ffe50a22ce8ff2d367e28275b90191f8c5:src/obsolete/main.f
      use thermal_evol, only: tevol, tem
      use legpol, only: allocate_legendre_polynomials, get_frel
      use legpol, only: compute_legendre_polynomials, getbl_init, blout
      use output, only: output_temperature, output_magnetic
      
      implicit none

      integer i,j,k,l
      integer icount, ndt,iterb
      integer its,nits
      integer, parameter :: nwrite=100
      real*8 tpart, tbpart, tyear, tbyear, cooling_output_dtad
      real*8 dtb_previous,pdot,bindex

      real*8 bpdip,bpdip_old
      real*8 eps

      !-------------------------------------------------------------------------
      ! Read input and initialize values
      !-------------------------------------------------------------------------
      call read_input_file()

      ! Initialize magnetic field normalization (bpdip).
      bpdip = bpolin

      ! Set the size of the angular and radial dimensions of the grid read from
      ! input so that the grid can be allocated with the proper dimensions.
      call set_grid_size(angular_dimension, radial_dimension)

      !-------------------------------------------------------------------------
      ! Allocate dynamic arrays throughout the whole code.
      !-------------------------------------------------------------------------
      call allocate_grid()
      call allocate_legendre_polynomials()

!-----------------------------------------------------------------------
! Generate the radial grid and compute the star structure
!-----------------------------------------------------------------------
      call build_structure
      !-----------------------------------------------------------------------
      ! Initialize the frel factors (boundary conditions correction for vacuum, NOT USED now).
      !-----------------------------------------------------------------------
      eps=1d0-1d0/belam(np+1)**2    !! 2*G*M/c**2*R
      call get_frel(eps)

!-----------------------------------------------------------------------
! Generate the angular grid and geometry arrays (length, surface, volume).
!-----------------------------------------------------------------------
      call build_grid
      print*,"Thermal grid: kmax, lmax: ",kmax,lmax
      print*,"Magnetic grid: nang, np: ",nang,np
      pdot=spindown_prefactor*bpdip**2*fchi/per

!-----------------------------------------------------------------------
! Initialize the Legendre polynomials.
!-----------------------------------------------------------------------
      call compute_legendre_polynomials()
!-----------------------------------------------------------------------
! Write output.   This is obsolete ! If we are not going to use IDL, remove it !
!-----------------------------------------------------------------------
! Array dimensions.
! Needed by some IDL macros. (See load.pro macro file.)
      open(2,file="outb/info_run.dat")
      write(2,*)lmax,kmax
      write(2,*)jcore/2
      write(2,*)np,nang
      close(2)
!-----------------------------------------------------------------------      
! Initialize time, redshifted temperature (tem=T*e^nu), etc.
!-----------------------------------------------------------------------
      tyear=0d0
      tpart=0d0
      tbyear=0d0
      tem=tinit
      tss=tinit
      sfluxb=0d0
!-----------------------------------------------------------------------
! Only needed for radiative conductivities.
!-----------------------------------------------------------------------
cc	call YAKITL
!-----------------------------------------------------------------------
! Initialize GETBL.
!-----------------------------------------------------------------------
! Needed for the normalization and for the boundary conditions.
      call getbl_init
!-----------------------------------------------------------------------
! Initialize the magnetic field.
!-----------------------------------------------------------------------
      if(resume == 0)then
        call binit

      else
        write(*,*)  'Restart not implemented yet !!!'
        stop
! Resume execution from a previous run.  
! Data is read from the file out/lastmaps.dat   
! TODO: Modify output to write this file !!
        open(2,file="out/lastmaps.dat")
        read(2,*)tyear,per,pdot
        write(*,'(a,e17.8)')"MAIN: Resume from (tyear):",tyear
        tbyear=tyear
        do i=0,nang+1
          do j=0,np+2
            read(2,*)br(i,j),bth(i,j),bphi(i,j),aphi(i,j)
          enddo
        enddo
        do k=1,kmax
          do l=1,lmax
            read(2,*) tem(k,l)
          enddo
        enddo
        close(2)
      endif
! Time advance: dt and dtb are in years (cooling advance).
      icount=0
      iterb=0

!-----------------------------------------------------------------------
! Start of time loop.
!-----------------------------------------------------------------------
! final_time and cooling_output_dt are read from the input file.
      write(*,'(a)')"MAIN: Starting time evolution..................."
cc	write(*,*)
cc	write(*,'(a,1pe10.3,a10,1pe10.3,a27,2i6)')
cc     &	"MAIN: Time (tyear)",tyear,"of (final_time)",final_time,
cc     &	". Counters (icount, iterb):",icount,iterb

      do while(tyear.lt.final_time)

! Fix the time step (overriding the above calculations).
! If dt and dtb are not sufficiently small, you may get numerical
! instabilities, resulting in a crash.
        if (itevol.eq.0) then
          dt=final_time  ! If T does not change, give only one dt step
        endif  
        cooling_output_dtad=cooling_output_dt
! -----------------------------------------------------------        
! Construct the tem0 array   (unredshifted temperature).
! -----------------------------------------------------------
        do l=1,lmax
          tem0(1:kmax,l)=tem(1:kmax,l)/benu(2*l-1)
        enddo
! -----------------------------------------------------------
!     	Shear modulus and maximum strength
! -----------------------------------------------------------
        call compute_shear(tem0)
! -----------------------------------------------------------    
! DIFFUSION COEFFICIENTS, ANISOTROPY FACTORS, MAGNETIC DIFFUSIVITY
! -----------------------------------------------------------
        call difdrv(tem0, br, bth, bphi, bm)
      
! -----------------------------------------------------------
!   bmed, eta (for tevol) and Omega_B*tau (only used for output purposes)      
! -----------------------------------------------------------
      bmed=0d0
      eta=0d0
      do k=2,kmax
       do l=1,lmax
        j = 2*l-1
        i = 2*k-2
        bmed(k,l)=bm(2*k-2,2*l-1)
        eta(k,l)=etab(i,j)
        omegatau(k,l)=fh(j)*bm(i,j)/etab(i,j)   
       enddo
      enddo
           
! MICROPHYSICS:  specific heat, emissivity (in 1e40 erg/s)
      
      call cvtot(kmax,lmax,tem0,bmed,tdc)
      call emiss(kmax,lmax,tem0,bmed,q_neutrino,rmc,emnu)

! SURFACE B.C. COEFFICIENTS (FLUXES AND EFFECTIVE TEMPERATURES)

      call fluxesb(tem0,sfluxb,cfluxb,tss,ienv)

      if(itevol.eq.0.or.tyear.eq.0d0.or.tpart+1d-4.ge.cooling_output_dtad)then
        call output_temperature(tyear,
     & tem, etab, omegatau,
     & tss, q_joule, emnu,
     & bindex, bpdip, icount, per, pdot, sfluxb)
        icount=icount+1
        tpart=0d0
      endif

! No magnetic field evolution:
      if(ibevol.eq.0)then
        q_joule=0d0
        call period_evol(dt,bpdip,bpdip_old,per,pdot,bindex)
! Magnetic field evolution:
      else
! Initialize the value dtb_previous, the old value of dtb, used by magnetic_evol.
        dtb_previous=dtb
! Determine the step size dtb and the number of iterations nits.
cc	      call BTIMESTEP(j2,dtb,nits)

! Magnetic field update.
        do while(tbyear.le.tyear+dt)

          q_joule_average = 0d0
          call calculate_magnetic_analysis(iterb,dtb_previous)

          !  Write output
          if(tbpart .ge. (magnetic_output_dt - 1d-6) .or. (tbyear .eq. 0d0))then
            call output_magnetic(iterb,tbyear)
            tbpart=0.d0
          endif

<<<<<<< HEAD:src/main.f
! REAL-TIME CRUSTAL FRACTURES (VERY HEAVY OUTPUTS)
!        call breaking(dtb,tbyear)

        tbyear=tbyear+dtb
        iterb=iterb+1
        dtbold=dtb


=======
          ! Magnetic field evolution.
          ! - dtb in main is in YEARS but in magnetic_evol in Myrs.
          call magnetic_evol(1.d-6*dtb,iterb)
          call q_joule_calc

          ! Consider the average over dt of q_joule (important if dtb<<dt)
          q_joule_average = q_joule_average + q_joule*dtb/dt
  
          ! Dipolar surface field updated
          bpdip = -2.d0*blout(1)  
          call period_evol(dt,bpdip,bpdip_old,per,pdot,bindex)
          if(tyear.eq.0d0) bindex=0d0
  
          tbpart=tbpart+dtb
          tbyear=tbyear+dtb
          iterb=iterb+1
          dtb_previous=dtb
  
          ! Magnetic stresses
cc		 call BREAKING(dtb,tbyear)
>>>>>>> 49bee1ffe50a22ce8ff2d367e28275b90191f8c5:src/obsolete/main.f
        enddo
      endif

      if (itevol .eq. 0) then
      ! No temperature evolution case
        tem=tinit
      else
      ! Advance temperature one timestep dt
        call tevol(dt,cfluxb,sfluxb,rmc,tdc)
      endif

      ! Update times.
      tyear=tyear+dt 
      tpart=tpart+dt

      ! End of time loop
      enddo

      stop
      end program

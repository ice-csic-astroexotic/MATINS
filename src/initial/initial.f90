!-------------------------------------------------------------------------------
! initial conditions
!@brief In this subroutine, we are defining the initial conditions for the magnetic 
! the thermal parts. 
! 
! To define the initial magnetic topology, we use select case
!>          It calls one of the following subroutines:
!>              bessel 
!>              test_pa
!>              hotspot
!>              pureTQ
!>              DP+TQ
!>              ScalarFunc
!>              Whister
! and they return the initial topology of the magnetic field in cubed-sphere coordinates
!
!
!> @Authors:
!>  Clara Dehman
!>  Stefano Ascenzi
!>  Daniele Viganò
!-------------------------------------------------------------------------------


subroutine initial_condition(it, time_read)

  use input_params
  use constants, only: PI
  use grid, only: nr, nang, nrt, nangt
  use grid, only: r, theta, phi, tem0, temp, temp_surf, T_core
  use grid, only: br, beta, bxi, b2, bm, lmax, q_joule
  use grid, only: jr, jeta, jxi, j2
  use grid, only: er, exi, eeta
  use grid, only: ievol, enu, etab
  use grid, only: f_spherical_to_cs
  use grid, only: curl_fnvol, fghost
  use grid, only: dot_prod
  use initial_magnetic
  use magnetic_evolution, only: compute_E, compute_joule

  implicit none

  integer, intent(out) :: it
  real*8, intent(out) :: time_read
  real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: xr, xxi, xeta, bth, bphi
  

  integer p, j, k, ic
  
  ! Temperature definition
  temp = T_init     ! Uniform temperature
 ! tem0 = T_init
  do ic = 1,nrt
    tem0(ic,:,:,:) = temp(ic,:,:,:)/enu(2*ic-1)
  enddo 
  T_core = T_init   ! Temperature of the core
  q_joule = 0.d0

  if (resume_checkpnt_number .gt. 0) then 

    if (ibevol) then
    ! Delete this once the checkpoint for magnetic field is written 
      print*, 'WARNING: At the moment the resume from checkpoint only contain the thermal evolution: B is set as the initial one'
    endif
    call set_from_checkpoint(resume_checkpnt_number, it, time_read)

  else 
    time_read = 0.d0
    it = 0
  endif  


  ! Set the initial magnetic topology.
  if (bpol_init == 0. .and. btor_init == 0.) then
    br = 0d0
    bxi = 0d0
    beta = 0d0
    b2 = 0d0
    bm = 0d0
    jr = 0d0
    jxi = 0d0
    jeta = 0d0
    er = 0d0
    exi = 0d0
    eeta = 0d0    
    j2 = 0d0

  else

    call binit()
  
   ! Calculation of electric currents
    xr   = 0d0
    xxi = 0d0
    xeta  = 0d0
    do p = 1, 6
      do j = 0, nang+1
        do k = 0, nang+1
          xr(ievol-1:,j,k,p) = br(ievol-1:,j,k,p)*enu(ievol-1:)
          xxi(ievol-1:, j,k,p) = bxi(ievol-1:,j,k,p)*enu(ievol-1:)
          xeta(ievol-1:,j,k,p) = beta(ievol-1:,j,k,p)*enu(ievol-1:)
        enddo
      enddo
    enddo
  
    ! Calculation of electrical currents
    call curl_fnvol(xr,xxi,xeta,jr,jxi,jeta,ievol)
    call fghost(jr,jxi,jeta)
  
    ! Calculation of electric fields
    call compute_E(br,bxi,beta,jr,jxi,jeta,er,exi,eeta)
    
    ! Calculation of magnetic field intensity
    call dot_prod(br,br,bxi,bxi,beta,beta,b2)
    bm = sqrt(b2)

    call envelope_model()

    ! Calculate J**2 used for the Joule dissipation
    call dot_prod(jr,jr,jxi,jxi,jeta,jeta,j2)
    call compute_joule

  endif

  
end subroutine initial_condition

!-------------------------------------------------------------
! set from checkpoint
!> @brief: This routine read the checkpoint file and fill the variable 
!> arreys
!
!> checkpoint_number[in]  ID number of the checkpoint file from which we want to restart 
!> it [out]               iteration number at the checkpoint
!> time_read [out]        time at the checkpoint
!
! 
!> @author:
!> Stefano Ascenzi
!-------------------------------------------------------------

subroutine set_from_checkpoint(checkpoint_number, it, time_read)

  use grid, only: nrt, nangt
  use grid, only: T_core, temp, last_timestep_print

  implicit none 

  integer, intent(out) :: it
  real*8, intent(out) :: time_read
  integer, intent(in) :: checkpoint_number
  character(len = 50)  :: filename, number_str
  integer :: i,j,k,p, imin_vis 
  imin_vis = 0

  write(number_str, "(I3)") checkpoint_number

  filename = "out/3D/checkpoint_"//trim(adjustl(number_str))//".dat"


  open (unit = 230, file = filename, status = 'old')
  
  read(230,*) it
  read(230,*) time_read
  read(230,*) last_timestep_print
  do p=1,6
    do i=imin_vis+1,nrt
      do j=0,nangt+1
        do k=0,nangt+1
          read(230,*) temp(i,j,k,p)
        end do
      end do
    end do
  end do

  close(230)

  T_core = temp(1,1,1,1)
  
end subroutine set_from_checkpoint
!-------------------------------------------------------------------------------
! Magneto Evolution 3D 
!-------------------------------------------------------------------------------
! Module: Magnetic Evolution
!
!> @author
!>  Clara Dehman
!>  Daniele Viganò
!
!> @brief Magnetic field time advance with finite-different methods
!>        It contains:
!>        subroutine magnetic_evol (called in main)
!>        The contained subroutines:
!>        subroutine euler
!>        subroutine RK4
!>        subroutine RK_substep
!>        subroutine curl_fnvol
!>        subroutine curl_fndiff
!>        subroutine compute_E
!>        subroutine compute_dB
!>        subroutine magnetic_bc_bessel
!>        subroutine magnetic_bc
!>        subroutine dtb_adaptive
!>        subroutine 
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
 
  ! Module imports 
  use input_params, only: e_scheme, time_advance
  use grid, only: nang, nr, ievol, etab, fh, nangt, nrt
  use grid, only: curl_fnvol, curl_fndiff, crossprod_cont, crossprod_cont_upwind
  use grid, only: f_cs_to_spherical, f_spherical_to_cs
  use grid, only: edge_average
  use grid, only: r, theta, phi, xi, eta
  use grid, only: fghost, dot_prod
  use grid, only: br, bxi, beta, bm, b2, enu, elambda
  use grid, only: jr, jeta, jxi, j2
  use grid, only: er, exi, eeta
  use grid, only: y_lm, dyth_lm, dyphi_lm, blm, bpdip, lmax
  use grid, only: area_r, vol
  use grid, only: q_joule
  use magnetic_analysis, only: compute_energy_balance, wint
  use OMP_LIB 

  implicit none

  procedure(), pointer, save :: crossprod_p

  contains

    subroutine magnetic_evol(dtb, t)
    
    implicit none

    ! Subroutine arguments ---------------------------------------------------
      real*8 :: dtb, t

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------

      ! Index for loops.
      integer :: j, k, p

!      integer :: l,m

      ! Auxiliary fields for intermediate calculations and increments.
      real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: xr, xxi, xeta

      !-------------------------------------------------------------------------
      ! Time-advance method for magnetic field evolution 
      select case (time_advance)
      case("Eul")
        call euler(dtb) 
      case ("RK2") 
        call RK2(dtb)
      case ("RK3") 
        call RK3(dtb)
      case ("RK4")
        call RK4(dtb)
      end select

      ! Calculation of new electric currents
      xr = 0d0
      xxi = 0d0
      xeta = 0d0
  
     !$OMP Parallel do private(p,j,k) collapse(3)
      do p = 1, 6
        do j = 0, nang+1
         do k = 0, nang+1
           xr(ievol-1:,j,k,p) = br(ievol-1:,j,k,p)*enu(ievol-1:)
           xxi(ievol-1:,j,k,p) = bxi(ievol-1:,j,k,p)*enu(ievol-1:)
           xeta(ievol-1:,j,k,p) = beta(ievol-1:,j,k,p)*enu(ievol-1:)
          enddo
        enddo
      enddo
     !$OMP end Parallel do
    
      call curl_fnvol(xr,xxi,xeta,jr,jxi,jeta,ievol)
      call fghost(jr,jxi,jeta)

      ! Calculation of new electric fields
      call compute_E(br,bxi,beta,jr,jxi,jeta,er,exi,eeta) 

      ! Calculation of magnetic field intensity
      call dot_prod(br,br,bxi,bxi,beta,beta,b2)
      bm = sqrt(b2)
      ! Calculate J**2 used for the Joule dissipation
      call dot_prod(jr,jr,jxi,jxi,jeta,jeta,j2)
      call compute_joule

      ! It uses the new Joule and Poynting flux related to the current timestep
      call compute_energy_balance(dtb)

  end subroutine magnetic_evol


  !---------------------------------------------------------------------------
  !! @brief Euler time advance.
  !
  !! @param[in]     dtb       Timestep in Myr
  !! @param[out]    dbr       increment in br   in UNIT_B
  !! @param[out]    dbxi      increment in bxi  in UNIT_B
  !! @param[out]    dbeta     increment in beta in UNIT_B 
  !
  !! Code owners:
  !!  Clara Dehman
  !---------------------------------------------------------------------------
  subroutine euler(dtb) 
  
    implicit none
      
    ! Subroutine arguments ---------------------------------------------------
    real*8, intent(in) :: dtb
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: dbr,dbxi,dbeta
    ! Local variables --------------------------------------------------------
    ! None
    ! ------------------------------------------------------------------------
    call compute_dB(dtb,er,exi,eeta,dbr,dbxi,dbeta)
    br  = br  + dbr
    bxi = bxi + dbxi
    beta = beta + dbeta
    call magnetic_bc(br,bxi,beta)
    call fghost(br,bxi,beta)
    !   call magnetic_bc_bessel(t+dtb,br,bxi,beta)

  end subroutine euler



  !---------------------------------------------------------------------------
  !> @brief Subroutine for the Runge-Kutta 2nd-order time-advance method
  !!
  !! @param[in]     dtb   Timestep in Myr
  !!
  !!--------------------------------------------------------------------------
  subroutine RK2(dtb)

    implicit none
    
    ! Subroutine arguments ---------------------------------------------------
    real*8, intent(in) :: dtb

    ! Local variables --------------------------------------------------------
    ! Auxiliary fields for intermediate calculations and increments.
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: dbr,dbxi,dbeta
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: brint,bxiint,betaint
    ! ------------------------------------------------------------------------

    call compute_dB(dtb,er,exi,eeta,dbr,dbxi,dbeta)
    brint   = br   + dbr
    bxiint  = bxi  + dbxi
    betaint = beta + dbeta
    call magnetic_bc(brint,bxiint,betaint)
    call fghost(brint,bxiint,betaint)

    call RK_substep(dtb,brint,bxiint,betaint,dbr,dbxi,dbeta)
    br   = 0.5d0*(brint + br + dbr)
    bxi  = 0.5d0*(bxiint + bxi + dbxi)
    beta = 0.5d0*(betaint + beta + dbeta)
    call magnetic_bc(br,bxi,beta)
    call fghost(br,bxi,beta)

  end subroutine RK2


  !---------------------------------------------------------------------------
  !> @brief Subroutine for the Runge-Kutta 3rd-order time-advance method
  !!
  !! @param[in]     dtb   Timestep in Myr
  !!
  !!--------------------------------------------------------------------------
  subroutine RK3(dtb)

    implicit none
    
    ! Subroutine arguments ---------------------------------------------------
    real*8, intent(in) :: dtb

    ! Local variables --------------------------------------------------------
    ! Auxiliary fields for intermediate calculations and increments.
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: dbr,dbxi,dbeta
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: brint,bxiint,betaint
    ! ------------------------------------------------------------------------
    ! In the first substep the values of B, J and E coming from the previous magnetic loop are used
    call compute_dB(dtb,er,exi,eeta,dbr,dbxi,dbeta)
    brint   = br   + 0.5d0*dbr
    bxiint  = bxi  + 0.5d0*dbxi
    betaint = beta + 0.5d0*dbeta
    call magnetic_bc(brint,bxiint,betaint)
    call fghost(brint,bxiint,betaint)

    call RK_substep(dtb,brint,bxiint,betaint,dbr,dbxi,dbeta)
    brint   = brint + 0.5d0*dbr
    bxiint  = bxiint + 0.5d0*dbxi
    betaint = betaint + 0.5d0*dbeta
    call magnetic_bc(brint,bxiint,betaint)
    call fghost(brint,bxiint,betaint)

    call RK_substep(dtb,brint,bxiint,betaint,dbr,dbxi,dbeta)
    br   = 1./3.*br + 2./3.*(brint + 0.5d0*dbr)
    bxi   = 1./3.*bxi + 2./3.*(bxiint + 0.5d0*dbxi)
    beta   = 1./3.*beta + 2./3.*(betaint + 0.5d0*dbeta)
    call magnetic_bc(br,bxi,beta)
    call fghost(br,bxi,beta)

  end subroutine RK3




  !---------------------------------------------------------------------------
  !> @brief Subroutine for the Runge-Kutta 4th-order time-advance method
  !!
  !! @param[in]     dtb   Timestep in Myr
  !!
  !!  Code owners:
  !!    Clara Dehman
  !!--------------------------------------------------------------------------
  subroutine RK4(dtb)

    implicit none
    
    ! Subroutine arguments ---------------------------------------------------
    real*8, intent(in) :: dtb

    ! Local variables --------------------------------------------------------
    ! Auxiliary fields for intermediate calculations and increments.
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: dbr1,dbxi1,dbeta1
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: dbr2,dbxi2,dbeta2
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: dbr3,dbxi3,dbeta3
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: dbr4,dbxi4,dbeta4
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: brint,bxiint,betaint
    ! ------------------------------------------------------------------------
    ! In the first substep the values of B, J and E coming from the previous magnetic loop are used
    call compute_dB(dtb,er,exi,eeta,dbr1,dbxi1,dbeta1)

    brint   = br   + 0.5d0*dbr1
    bxiint  = bxi  + 0.5d0*dbxi1
    betaint = beta + 0.5d0*dbeta1
    call magnetic_bc(brint,bxiint,betaint)
    call fghost(brint,bxiint,betaint)
   ! call magnetic_bc_bessel(t+0.25d0*dtb,brint,bxiint,betaint)

    ! In the following three sub-steps J and E are recalculated using Bint
    call RK_substep(dtb,brint,bxiint,betaint,dbr2,dbxi2,dbeta2)
    brint   = br   + 0.5d0*dbr2
    bxiint  = bxi  + 0.5d0*dbxi2
    betaint = beta + 0.5d0*dbeta2
    call magnetic_bc(brint,bxiint,betaint)
    call fghost(brint,bxiint,betaint)
   ! call magnetic_bc_bessel(t+0.5d0*dtb,brint,bxiint,betaint)
    call RK_substep(dtb,brint,bxiint,betaint,dbr3,dbxi3,dbeta3)
    brint   = br   + 0.5d0*dbr3
    bxiint  = bxi  + 0.5d0*dbxi3
    betaint = beta + 0.5d0*dbeta3
    call magnetic_bc(brint,bxiint,betaint)
    call fghost(brint,bxiint,betaint)
    !call magnetic_bc_bessel(t+0.75d0*dtb,brint,bxiint,betaint)
    call RK_substep(dtb,brint,bxiint,betaint,dbr4,dbxi4,dbeta4)
    br   = br   + (1./6.*dbr1   + 1./3.*dbr2   + 1./3.*dbr3   + 1./6.*dbr4   )
    bxi  = bxi  + (1./6.*dbxi1  + 1./3.*dbxi2  + 1./3.*dbxi3  + 1./6.*dbxi4  )
    beta = beta + (1./6.*dbeta1 + 1./3.*dbeta2 + 1./3.*dbeta3 + 1./6.*dbeta4 )
    call magnetic_bc(br,bxi,beta)
    call fghost(br,bxi,beta)
 !   call magnetic_bc_bessel(t+dtb,br,bxi,beta)
  end subroutine RK4

  !!------------------------------------------------------------------------
  !> @brief Subroutine for the substeps of the Runge-Kutta time-advance method
  !!
  !! @param[in]     dtb   Timestep in Myr
  !! @param[in]     brint     br  
  !! @param[in]     bxiint    bth  
  !! @param[in]     betaint   bphi 
  !! @param[out]    dbr       increment in br  
  !! @param[out]    dbxi      increment in bxi
  !! @param[out]    dbeta     increment in beta 
  !!
  !! Code owners:
  !!    Clara Dehman
  !!------------------------------------------------------------------------------
  subroutine RK_substep(dtb,brint,bxiint,betaint,dbr,dbxi,dbeta)

    ! Module imports -------------------------------------------------------------

    implicit none
    
    ! Subroutine arguments -------------------------------------------------------
    real*8, intent(in) :: dtb

    ! Local variables ------------------------------------------------------------
    ! Auxiliary fields for intermediate calculations and increments.
    integer :: j,k,p
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: brint, bxiint, betaint
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(out) :: dbr, dbxi, dbeta
    real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: xr, xxi, xeta

    ! ----------------------------------------------------------------------------    
    xr   = 0d0
    xxi = 0d0
    xeta  = 0d0
    ! TBD: OPENMP?
    do p = 1, 6
      do j = 0, nang+1
        do k = 0, nang+1
          xr(ievol-1:,j,k,p) = brint(ievol-1:,j,k,p)*enu(ievol-1:)
          xxi(ievol-1:,j,k,p) = bxiint(ievol-1:,j,k,p)*enu(ievol-1:)
          xeta(ievol-1:,j,k,p) = betaint(ievol-1:,j,k,p)*enu(ievol-1:)
        enddo
      enddo
    enddo
    call curl_fnvol(xr,xxi,xeta,jr,jxi,jeta,ievol)
  !  call curl_fndiff(xr,xxi,xeta,jr,jxi,jeta,ievol)
    call fghost(jr,jxi,jeta)
    ! Calculation of new electric fields
    call compute_E(brint,bxiint,betaint,jr,jxi,jeta,er,exi,eeta)
      
    call compute_dB(dtb,er,exi,eeta,dbr,dbxi,dbeta)
    
  end subroutine RK_substep
    
!!------------------------------------------------------------------------------
!> @brief Subroutine for magnetic field increment, given the three components
!! er, exi, eeta of the electric field in the cubed sphere coordinate system 
!!
!! @param[in]     dtb       Timestep in Myr
!! @param[in]     er        radial component of the electric field
!! @param[in]     exi       xi component of the electric field
!! @param[in]     eeta      eta component of the electric field
!! @param[out]    dbr       increment in br in UNIT_B
!! @param[out]    dbxi      increment in bxi in UNIT_B
!! @param[out]    dbeta     increment in beta in UNIT_B
!!
!! Code owners:
!!   Clara Dehman
!!------------------------------------------------------------------------------
  subroutine compute_dB(dtb,er,exi,eeta,dbr,dbxi,dbeta)
    
        implicit none
    
        ! Subroutine arguments -------------------------------------------------------
        real*8, intent(in) :: dtb
        real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(in) :: er, exi, eeta
        real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(out) :: dbr,dbxi,dbeta
        real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6) :: curler, curlexi, curleeta
        integer :: i, j, k, p
        ! ----------------------------------------------------------------------------
        dbr = 0.d0
        dbxi = 0.d0
        dbeta = 0.d0

        call curl_fnvol(er,exi,eeta,curler,curlexi,curleeta,ievol-1)
      !  call curl_fndiff(er,exi,eeta,curler,curlexi,curleeta,ievol-1)

        !$OMP Parallel do private(i,j,k)  collapse(4)
        do p=1, 6
        do k=0, nang+1
        do j=0, nang+1
        do i=0, nr+1
          dbr(i,j,k,p) =  - dtb*curler(i,j,k,p)
          dbxi(i,j,k,p) =  - dtb*curlexi(i,j,k,p)
          dbeta(i,j,k,p) =  - dtb*curleeta(i,j,k,p)
        end do
        end do
        end do
        end do
        !$OMP end Parallel do

    end subroutine compute_dB


  !!------------------------------------------------------------------------
  !> @brief Subroutine for electric field
  !!
  !! @param[in]     brin,bxiin,betain   magnetic field in UNIT_B
  !! @param[in]     jrin,jxiin,jetain   Electrical currents in UNIT_B/UNIT_R
  !! @param[out]    erin,exiin,eetain   Electric field
  !!
  !! Code owners:
  !!  Clara Dehman
  !!------------------------------------------------------------------------
    subroutine compute_E(brin,bxiin,betain,jrin,jxiin,jetain,erin,exiin,eetain)
    
        implicit none
    
        ! Subroutine arguments -------------------------------------------------------
        real*8, dimension (0:nr+1,0:nang+1,0:nang+1, 6), intent(in) :: brin,bxiin,betain,jrin,jxiin,jetain
        real*8, dimension (0:nr+1,0:nang+1,0:nang+1, 6), intent(out) :: erin,exiin, eetain
        real*8, dimension (0:nr+1, 0:nang+1, 0:nang+1, 6) :: ehallr, ehallxi, ehalleta 
        
        ! Local variables ------------------------------------------------------------
        ! Auxiliary index for loops
         integer :: j, k, p
     
      ! ----------------------------------------------------------------------------
        erin=0d0
        exiin=0d0
        eetain=0d0


        select case(e_scheme)

          case ("Center")
            crossprod_p => crossprod_cont
          case ("Upwind")
            crossprod_p => crossprod_cont_upwind
          case default
            write(*,*) "<error>", &
                    & "[BEVOL]", &
                    & "Invalid scheme, choose eithr Center or Upwind: ", &
                    & e_scheme
          stop

        end select

    
       call crossprod_p(jrin,jxiin,jetain,brin,bxiin,betain,ehallr,ehallxi,ehalleta,ievol)
             
      !$OMP Parallel do private(p,j,k)   collapse(3)
       do p = 1, 6        
        do k=0,nang+1
        do j=0,nang+1
      ! Ohmic term  + Hall term
        erin(:,j,k,p) = etab(:,j,k,p) * jrin(:,j,k,p) + fh(:)*ehallr(:,j,k,p)
        exiin(:,j,k,p) = etab(:,j,k,p) * jxiin(:,j,k,p) + fh(:)*ehallxi(:,j,k,p)
        eetain(:,j,k,p) = etab(:,j,k,p) * jetain(:,j,k,p) + fh(:)*ehalleta(:,j,k,p)
        end do 
        end do 
       end do   
      !$OMP end Parallel do

        exiin(0:1,:,:,:) = 0.d0
        eetain(0:1,:,:,:) = 0.d0
        exiin(2,:,:,:) = 0.5d0*exiin(3,:,:,:)
        eetain(2,:,:,:) = 0.5d0*eetain(3,:,:,:)
     
       call edge_average(ievol,erin,exiin,eetain)

  end subroutine compute_E


  !!-----------------------------------------------------------------------
  !> @brief Subroutine for the Boundary Conditions considering an analytical 
  !! solution for the bessel test mainly 
  !! 
  !! Code owners:
  !!  Clara Dehman
  !!  Daniele Viganò
  !!-----------------------------------------------------------------------
   subroutine magnetic_bc_bessel(time,brout,bxiout,betaout) 

    implicit none 
      
    ! Subroutine arguments -------------------------------------------------------
       real*8, intent(in) :: time
       real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(out) :: brout, bxiout, betaout

       brout(0,:,:,:) = 0 !decay(time)*brin(0:1,:,:,:) 
      ! brout(nr+1,:,:,:) = 0 !decay(time)*brin(nr:nr+1,:,:,:)
       bxiout(0,:,:,:) = 0 !decay(time)*bxiin(0:1,:,:,:)
      ! bxiout(nr+1,:,:,:) = 0 !decay(time)*bxiin(nr:nr+1,:,:,:)
       betaout(0,:,:,:) = 0 !decay(time)*betain(0:1,:,:,:)
      ! betaout(nr+1,:,:,:) = 0 !decay(time)*betain(nr:nr+1,:,:,:)
       brout(1,:,:,:) = 0.5d0*brout(2,:,:,:)
      ! brout(nr,:,:,:) = 0.5d0*brout(nr-1,:,:,:)
       bxiout(1,:,:,:) = 0.5d0*bxiout(2,:,:,:)
      ! bxiout(nr,:,:,:) = 0.5d0*bxiout(nr-1,:,:,:)
       betaout(1,:,:,:) = 0.5d0*betaout(2,:,:,:)
      ! betaout(nr,:,:,:) = 0.5d0*betaout(nr-1,:,:,:)
    
    end subroutine magnetic_bc_bessel  


  !-----------------------------------------------------------------------
  !> @brief This function calculates the decaying time of Bessel Function
  ! 
  ! Code owners:
  !  Daniele Viganò
  !-----------------------------------------------------------------------
  ! real*8 function decay(time)
  !   implicit none
  !   real*8, intent(in) :: time

  !   decay = dexp(- etab*alpha**2*(time))
  ! end function decay

 !!-----------------------------------------------------------------------
  !> @brief In this subroutine, we calculate the weights using the spectral 
  !! decomposition of the radial component of the magnetic field at the surface,
  !! br(i=nr). Once the weights are calculated, we reconstruct the radial 
  !! component of the magnetic field above the surface, Br for r > R. Moreover, 
  !! these weights are also used to reconstruct the theta and phi components of 
  !! magnetic field, Bth and Bphi respectively, for r >= R. 
  !! Our formalism consists of considering a potential outer boundary conditions. 
  !! 
  !! Code owners:
  !!  Clara Dehman
  !!  Daniele Viganò 
  !!-----------------------------------------------------------------------
    subroutine magnetic_bc(br_bc,bxi_bc,beta_bc)
    implicit none 
     integer i, j, k, p, l, m
     real*8, dimension (0:nr+1,0:nang+1,0:nang+1,6), intent(inout) :: br_bc, bxi_bc, beta_bc  
     real*8, dimension (nr+1:nr+1,0:nang+1,0:nang+1,6) :: br_temp, bth_temp, bphi_temp, bxi_temp, beta_temp 
     real*8, dimension(0:lmax,-lmax:lmax) :: blm_out
     real*8, dimension(0:nang+1,0:nang+1,6) :: sin_theta
     real*8 :: r_const, rl
    ! ------------------------------------------------------------------------- 
    ! Inner Boundary Conditions

    !$OMP Parallel do private(j,k,p) collapse(3)
    do p = 1, 6  
      do k = 0, nang+1
        do j = 0, nang+1
           br_bc(0,j,k,p) = 0.d0
           bxi_bc(0,j,k,p) = 0.d0
           beta_bc(0,j,k,p) = 0.d0 
          ! solving the inner odd-even decoupling of the magnetic field 
           bxi_bc(1,j,k,p) = r(2)/r(1)*bxi_bc(2,j,k,p) 
           beta_bc(1,j,k,p) = r(2)/r(1)*beta_bc(2,j,k,p) 
         end do
       end do
     end do
     !$OMP end Parallel do 
    ! ------------------------------------------------------------------------- 
    ! Outer Boundary Conditions  

     blm = 0.d0

     call get_blm(br_bc(nr,:,:,:),blm_out)
     blm = blm_out
     bpdip = dsqrt(sum(blm(1,:)**2))/elambda(nr)

    ! Reconstructing the magnetic field components for r >= R 
     br_temp = 0.d0 
     bth_temp = 0.d0 
     bphi_temp = 0.d0 

    !$OMP Parallel private(beta_temp, bxi_temp, i)

     !$OMP do private(i,j,p) collapse(3)
     do p = 1, 6
       do j = 0, nang+1
         do i = 0, nang+1
           sin_theta(i,j,p) = dsin(theta(i,j,p)) + 1d-50
         end do
       end do
     end do
     !$OMP end do 
     i = nr+1


     !$OMP do reduction(+: br_temp, bth_temp, bphi_temp) private(m, l, p, k, j, rl) collapse(2)
     do l = 1, lmax
       do m = 0, lmax
         rl = (r(nr)/r(i))**(l+2)
         if (m == 0) then
          br_temp(i,:,:,:) = br_temp(i,:,:,:) + blm(l,m)*(l+1)*y_lm(:,:,:,l,m)*rl/elambda(i)
          bth_temp(i,:,:,:) =  bth_temp(i,:,:,:) - blm(l,m)*dyth_lm(:,:,:,l,m)*rl
          bphi_temp(i,:,:,:) = bphi_temp(i,:,:,:) - blm(l,m)* &
                        dyphi_lm(:,:,:,l,m)/(sin_theta(:,:,:))*rl
        else 
          br_temp(i,:,:,:) = br_temp(i,:,:,:) + (blm(l,m)*y_lm(:,:,:,l,m)+blm(l,-m)*y_lm(:,:,:,l,-m))*(l+1)*rl/elambda(i)
          bth_temp(i,:,:,:) =  bth_temp(i,:,:,:) - (blm(l,m)*dyth_lm(:,:,:,l,m)+blm(l,-m)*dyth_lm(:,:,:,l,-m))*rl
          bphi_temp(i,:,:,:) = bphi_temp(i,:,:,:) - (blm(l,m)*dyphi_lm(:,:,:,l,m)+blm(l,-m)*dyphi_lm(:,:,:,l,-m)) & 
          &  /(sin_theta(:,:,:))*rl
        end if 

       end do 
     end do
     !$OMP end do 

     ! Manually set that Bphi at the pole is 0 (since it's ill-defined)
     bphi_temp(i,nang/2+1,nang/2+1,5:6) = 0.d0

!    method 1:      
     bxi_temp = 0.d0
     beta_temp = 0.d0
     call f_spherical_to_cs(bth_temp,bphi_temp,bxi_temp,beta_temp,nr+1)
    !$OMP do private(j,k,p) collapse(3) 
    do p = 1, 6  
      do k = 0, nang+1
        do j = 0, nang+1
         br_bc(nr+1,j,k,p) = br_temp(nr+1,j,k,p)
         bxi_bc(nr+1,j,k,p) = bxi_temp(nr+1,j,k,p)
         beta_bc(nr+1,j,k,p) = beta_temp(nr+1,j,k,p)

         ! Solving the outer odd-even decoupling for the angular components of the magnetic field 
         bxi_bc(nr,j,k,p) = 0.5d0*(bxi_bc(nr-1,j,k,p) + bxi_bc(nr+1,j,k,p))
         beta_bc(nr,j,k,p) = 0.5d0*(beta_bc(nr-1,j,k,p) + beta_bc(nr+1,j,k,p))

         end do
       end do
     end do
!     method 2: 
!
!     call f_cs_to_spherical(bxi_bc,beta_bc,bth_bc,bphi_bc,nr-1)
!     bth_temp(nr,:,:,:) = 0.5d0*(bth_bc(nr-1,:,:,:) + bth_temp(nr+1,:,:,:))
!     bphi_temp(nr,:,:,:) = 0.5d0*(bphi_bc(nr-1,:,:,:) + bphi_temp(nr+1,:,:,:))
!
!    ! the output to be used as outer magnetic boundary conditions 
!     bxi_temp = 0.d0
!     beta_temp = 0.d0
!     br_bc(nr+1,:,:,:) = br_temp(nr+1,:,:,:)
!     call f_spherical_to_cs(bth_temp,bphi_temp,bxi_temp,beta_temp,nr)
!     bxi_bc(nr:nr+1,:,:,:) = bxi_temp(nr:nr+1,:,:,:)
!     beta_bc(nr:nr+1,:,:,:) = beta_temp(nr:nr+1,:,:,:)
     !$OMP end do 
          !$OMP end Parallel
    end subroutine magnetic_bc  


    subroutine get_blm(brin,blm_out)

      implicit none
  
      real*8, dimension(0:nang+1,0:nang+1,1:6), intent(in) :: brin
      real*8, dimension(0:lmax,-lmax:lmax) :: blm_out
  
      integer :: l, m, p
      blm_out = 0.d0
      !$OMP Parallel do reduction(+: blm_out) private(m, l, p) collapse(2)
      do p = 1, 6
        do l = 0, lmax
          do m = -l, l
          
            blm_out(l,m) = blm_out(l,m) + elambda(nr)/((l+1.d0)*(r(nr)*r(nr))) * & 
            &       0.25d0*sum(area_r(nr,1:nang,1:nang)*wint(1:nang,1:nang)*     &
            &           y_lm(1:nang,1:nang,p,l,m)*brin(1:nang,1:nang,p))       
    
          ! !Summing the even grid points 
          !  blm_out(l,m) = blm_out(l,m) + 1.d0/((l+1.d0)*(r(nr)*r(nr))) * & 
          !  &       sum(area_r(nr,2:nang-1:2,2:nang-1:2)*y_lm(2:nang-1:2,2:nang-1:2,p,l,m) &
          !  &             *brin(2:nang-1:2,2:nang-1:2,p))       
  
  ! HIGHER ORDER (TO BE TESTED, IT GIVES MUCH MORE RESIDUALS...)
  !          blm_out(l,m) = blm_out(l,m) + 1.d0/((l+1.d0)*(r(nr))**2.d0) * & 
  !   &           sum(area_r(nr,2:nang-1:2,2:nang-1:2)*y_lm(2:nang-1:2,2:nang-1:2,p,l,m)*( &
  !   &            0.25d0*brin(2:nang-1:2,2:nang-1:2,p)   &
  !   &         +  0.125d0*brin(1:nang-2:2,2:nang-1:2,p)  &
  !   &         +  0.125d0*brin(3:nang:2,2:nang-1:2,p)    &
  !   &         +  0.125d0*brin(2:nang-1:2,1:nang-2:2,p)  &
  !   &         +  0.125d0*brin(2:nang-1:2,3:nang:2,p)    &
  !   &         +  0.0625d0*brin(1:nang-2:2,1:nang-2:2,p) &
  !   &         +  0.0625d0*brin(2:nang:2,1:nang-2:2,p)   &
  !   &         +  0.0625d0*brin(2:nang:2,2:nang:2,p)     &
  !   &         +  0.0625d0*brin(1:nang-2:2,1:nang-2:2,p) ))
          end do
        end do 
      end do
       !$OMP end Parallel do
    end subroutine get_blm

    !!------------------------------------------------------------------------------
    !> @brief Subroutine for the Joule heating calculation
    !!
    !!  Code owners:
    !!    Clara Dehman
    !!    Daniele Viganò
    !!------------------------------------------------------------------------------
    !! See Viganò et al. 2012 paper for details and notation
    !! The Joule rate is defined as
    !! Q_joule = e^2nu*eta*J^2
    !!------------------------------------------------------------------------------  
    subroutine compute_joule()
  
      ! Module imports -------------------------------------------------------------
      use constants, only: UNIT_JOULE
      use grid, only: nangt, nrt
      use grid, only: etab, enu
      use grid, only: q_joule, j2
  
      implicit none
  
      ! Local variables ------------------------------------------------------------
      ! Auxiliary variable for loops
      integer :: i, j, k, p
      integer :: it, jt, kt
    
      ! ----------------------------------------------------------------------------

      q_joule = 0.d0

      ! In the crust-confined case, 
      ! the supercurrents arising from the discontinuity at the interface are neglected.
      ! In what we are considering as inner BC, they are located at ievol-2 

      do p = 1,6 
        do jt = 1, nangt
          j = 2*jt ! corresponds to the center of the thermal cell 
          do kt = 1, nangt
            k = 2*kt
            ! The first cell where Joule is to be calculated is at it=(ievol+2)/2, which
            ! corresponds to the centre of the thermal cell since ievol=2 
            do it = (ievol+2)/2,nrt
              i=2*it-1
              ! Zero-order approximation: center of the thermal cell 
              ! q_joule(it,jt,kt,p) = - UNIT_JOULE*etab(i,j,k,p)*j2(i,j,k,p)

              ! Refined version by considering the centers, face centers, edges and corners:
              q_joule(it,jt,kt,p) = - UNIT_JOULE*etab(i,j,k,p)*( (1./8.)*j2(i,j,k,p) + &
          &           (1./16.)*(j2(i+1,j,k,p) + j2(i-1,j,k,p) + j2(i,j+1,k,p) + j2(i,j-1,k,p) + &
          &                      j2(i,j,k+1,p) + j2(i,j,k-1,p) ) + & 
          &           (1./32.)*(j2(i+1,j+1,k,p) + j2(i+1,j-1,k,p) + j2(i+1,j,k+1,p) + j2(i+1,j,k-1,p) + &
          &                     j2(i-1,j+1,k,p) + j2(i-1,j-1,k,p) + j2(i-1,j,k+1,p) + j2(i-1,j,k-1,p) + &
          &                     j2(i,j+1,k+1,p) + j2(i,j+1,k-1,p) + j2(i,j-1,k+1,p) + j2(i,j-1,k-1,p) ) + & 
          &           (1./64.)*(j2(i+1,j+1,k+1,p) + j2(i+1,j+1,k-1,p) + j2(i+1,j-1,k+1,p) + &
          &                     j2(i+1,j-1,k-1,p) + j2(i-1,j+1,k+1,p) + j2(i-1,j+1,k-1,p) + &
          &                     j2(i-1,j-1,k+1,p) + j2(i-1,j-1,k-1,p) ) )
  
            end do
          end do
        end do
      end do 
  
    end subroutine compute_joule

!!------------------------------------------------------------------------------
    !> @brief Subroutine dtb_adaptive: it calculates the needed timestep
    !!
!!------------------------------------------------------------------------------

    subroutine dtb_adaptive(dtb_ad)

      use grid, only: lr, lxi

      implicit none

      real*8, intent(out) :: dtb_ad

      dtb_ad = minval([lr(ievol:nr)**2,lxi(ievol:nr,1:nang,1:nang)**2])/ &
     &        (maxval(fh(ievol:nr))*maxval(bm(ievol:nr,1:nang,1:nang,:)) + maxval(etab(ievol:nr,1:nang,1:nang,:)))

    ! Alternative:
    !   real*8, dimension(0:nr+1) :: dtb_local
    !   integer i 

    end subroutine dtb_adaptive

end module magnetic_evolution

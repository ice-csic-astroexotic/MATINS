!-------------------------------------------------------------------------------
! Module: thermal_evolution
!
!> @author
! Stefano Ascenzi 
! Daniele Viganò
!
!> @brief This module is responsible for handling the grid and structure.
!
!-------------------------------------------------------------------------------

Module thermal_evolution

  use constants, only: PI, UNIT_R, UNIT_EN, UNIT_TIME
  use grid, only: nrt, nangt
  use grid, only: r, xi, eta
  use grid, only: area_r, area_xi, area_eta, vol
  use grid, only: C, D, X, Y, delta
  use grid, only: enu, elambda
  use grid, only: edge_wt, pmt
  use grid, only: br, bxi, beta, bm! magnetic field
  use grid, only: temp
  use grid, only: rmax
  use grid, only: theta, phi ! DEBUG

  !use mycrophysics, only: core_cooling_rhs, core_cooling_rhs_deriv

  ! mycrophysical quantities

  use grid, only: cv, kappa_perp_arr, omegatau_arr
  use grid, only: q_neutrino, q_neutrino_der, q_joule
  use grid, only: T_core

  ! envelope quantities

  use grid, only: cfluxb, sfluxb

  ! diagnostic quantities

  use grid, only: flux_r_out, flux_xi_xip, flux_eta_etap

  ! WARNING: when you call cv from grid remember that the order of indices is different 
  ! with respect to the one used here (not a big problem though, we will get a segmentation
  ! fault)
 contains


subroutine tevol(dt)



 ! Define all the matrix elements interior to the border
 ! Assume by now constant scalars for cv and kappa (isotropical), with no sources and no fluxes
 ! Here you need to use the things defined in grid, assume you have them (Daniele will call them properly)

  implicit none

  real*8, intent(in) :: dt
  real*8, parameter :: c_v = 1d0!, kappa_perp = 1d0, omegatau = 1.d1
  real*8, parameter :: T_ext = 1d0, T_int = 100d0
  real*8, parameter :: alpha = 8.d0, T_0 = 20.d0 !deprecated
  real*8, parameter :: sb_constant = 1.78843596d11! Stef-Boltz const in sim. units !1.d-5 
  character(len=20) :: name, name_kappa, name_flux

  integer i, p, j, k, ic, jc, kc, n
  integer n_jinf, n_jsup, n_kinf, n_ksup
  integer kprime
  integer i_print ! choose the radial layer to print in gnuplot
  integer bandwd ! half-width of the matrix central band
  integer INFO   ! 0 if the linear system solver worked succesfully
  integer n1_jinf, n2_jinf, n1_jsup, n2_jsup
  integer n1_kinf, n2_kinf, n1_ksup, n2_ksup
  integer, dimension(1:4) :: n1_corner, n2_corner
  integer, dimension(1:6*nrt*(nangt+2)*(nangt+2)) :: IPIV !vector output of the linear sistem solver
  integer, save :: counter = 1
  real*8 :: precoeff
  real*8 :: identity1, identity2 !variables to check if the matrix elements are correct
  real*8 :: T_surface, Flux_surface, Flux_core
  real*8 :: kappa_perp, omegatau
  real*8 :: omegatau_out, omegatau_xip, omegatau_etap
  real*8 :: kappa_perp_out, kappa_perp_xip, kappa_perp_etap
  real*8 :: deriv_envelope ! dT_s/dT
  !real*8 :: A_kappa, B_kappa, G_kappa
  !real*8 :: E_kappa, F_kappa, H_kappa, I_kappa, J_kappa, K_kappa
  !real*8 :: L_kappa, M_kappa, N_kappa, O_kappa, P_kappa
  !real*8 :: value_inf, value_sup
!  real*8 :: zxixi_j0, zxir_j0, zxieta_j0 !maybe it's better to definie this element in another way, like including it in the vector with index 0
!  real*8, dimension(0:1,0:2,0:2) :: value_boundary
  real*8, dimension(:, :, :, :), allocatable :: kappa_rr_out, kappa_xixi_xip, kappa_etaeta_etap
  real*8, dimension(:, :, :, :), allocatable :: kappa_rxi_out, kappa_reta_out
  real*8, dimension(:, :, :, :), allocatable :: kappa_xir_xip, kappa_xieta_xip
  real*8, dimension(:, :, :, :), allocatable :: kappa_etar_etap, kappa_etaxi_etap
  real*8, dimension(:, :, :, :), allocatable :: zrr, zxixi, zetaeta
  real*8, dimension(:, :, :, :), allocatable :: zrxi, zreta
  real*8, dimension(:, :, :, :), allocatable :: zxir, zxieta
  real*8, dimension(:, :, :, :), allocatable :: zetar, zetaxi
  real*8, dimension(:, :, :, :), allocatable :: h
  !real*8, dimension(:, :, :, :), allocatable :: cv
  real*8, dimension(:, :), allocatable, save :: matrix
  real*8, dimension(:,:), allocatable :: matrix_solve
  real*8, dimension(:), allocatable :: source


  ! ***** VARIABLES FOR TESTs (delete later) *******
  
  !real*8 :: pl_index
  real*8 :: T_analit !Analitic temperature profile
  real*8, save :: time = 0.d0
  real*8 :: ang_term !For the Anisotropic Perez-Azorin test
  real*8 :: cos_th_prime, beta_ang = 54.735610317245346d0*PI/180., gamma_ang = 45.d0*PI/180.
  real*8 :: temp_serv1, temp_serv2, temp_serv3, temp_serv4, temp_serv5 !service variable for temperature 

  logical, save :: FirstCall = .TRUE.
  logical :: CoolingOff = .FALSE. !When CoolingOff is .true. the cooling is off
  logical :: CalculateFlux = .FALSE. !When is true the netflux in each cell of a given radial layer is calc. and printed

  real*8, dimension(:, :, :, :), allocatable :: flux_r_anl_out, flux_xi_anl_xip, flux_eta_anl_etap
  !real*8, dimension(:, :, :), allocatable :: netflux_r, netflux_xi, netflux_eta !we calculate the fluxes at only one radial layer 

  ! temp_analt is calculated at t+dt at the end of the routine, we save the value because it can be used at the next
  ! time step to calculate the flux.
  real*8, dimension(:, :, :, :), allocatable, save :: temp_anlt
  !real*8, save, dimension(:), allocatable :: temp
 ! real*8 :: lat, long, ang
 real*8 :: maxbm


  !if (FirstCall) then
  !  allocate(temp(0: 6*nrt*(nangt+2)*(nangt+2)-1))
  !endif

  
  !pl_index = log(T_ext/T_int)/log(r(1)/r(2*nrt-1)) !this is the index of a PL temperature profile
  !*************************************************


  ! ******** Variables for the source term ******

   real*8 :: Q, nu_loss, nu_loss_deriv
   real*8 :: lat = 0.0, long = 0.0, ang = 20.0*PI/180.0
   real*8 :: condition
  ! *********************************************


!  OPEN(unit = 100, file='matrix_input.dat',status='replace')
  !OPEN(unit = 101, file='matrix_output.dat',status='unknown')
  !OPEN(unit = 102, file='source_input.dat',status='unknown')
  !OPEN(unit = 103, file='source_output.dat',status='unknown')


! allocate the variable h. The index p is not used now, but we keep it.
  allocate(h(1:nrt, 1:6, 1:nangt, 1:nangt))
  !allocate(cv(1:nrt, 1:6, 1:nangt, 1:nangt))

  ! By now we play with a constant heat capacity
  !cv = c_v

! allocate the component of the conductivity tensor
    
  ! BE CAREFUL: basically we are excluding one interface here
  ! that I assume will enter in the boundary condition. THINK ABOUT IT
  allocate(kappa_rr_out(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(kappa_xixi_xip(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(kappa_etaeta_etap(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(kappa_rxi_out(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(kappa_reta_out(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(kappa_xir_xip(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(kappa_xieta_xip(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(kappa_etar_etap(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(kappa_etaxi_etap(1:nrt, 1:6, 0:nangt, 0:nangt))

! Allocate the auxiliary variables
! WARNING: check the dimensions are correct

  allocate(zrr(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(zxixi(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(zetaeta(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(zrxi(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(zreta(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(zxir(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(zxieta(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(zetar(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(zetaxi(1:nrt, 1:6, 0:nangt, 0:nangt))

  !set the half-width of the matrix central band (diagonal excluded)

  bandwd = (nangt +2)*(6*(nangt+2)+1)

! Allocate the matrix and source.
! (Reduced) matrix dimension MxN
! M = 3*bandwd+1
! N = 6*(nrt+1)*(nangt+2)*(nangt+2)

  !allocate(matrix(0: 6*(nrt+1)*(nangt+2)*(nangt+2)-1, 0:6*(nrt+1)*(nangt+2)*(nangt+2)-1)) ! full matrix
  if (FirstCall) then
    allocate(matrix(0: 3*bandwd, 0:6*nrt*(nangt+2)*(nangt+2)-1))  ! reduced matrix
  endif
  allocate(matrix_solve(0: 3*bandwd, 0:6*nrt*(nangt+2)*(nangt+2)-1))  ! reduced matrix_solve
  allocate(source(0: 6*nrt*(nangt+2)*(nangt+2)-1))

  ! Allocation Variables for TESTs ********************************************


  allocate(flux_r_anl_out(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(flux_xi_anl_xip(1:nrt, 1:6, 0:nangt, 0:nangt))
  allocate(flux_eta_anl_etap(1:nrt, 1:6, 0:nangt, 0:nangt))

!  allocate(netflux_r(1:6, 1:nangt, 1:nangt))
!  allocate(netflux_xi(1:6, 1:nangt, 1:nangt))
!  allocate(netflux_eta(1:6, 1:nangt, 1:nangt))

  if (FirstCall) then
    allocate(temp_anlt(1:nrt, 0:nangt+1, 0:nangt+1, 1:6))
  endif

  ! ****************************************************

! Initialize the coeffficients of the auxiliary variables to zero

  zrr = 0.d0
  zxixi = 0.d0
  zetaeta = 0.d0
  zrxi = 0.d0
  zreta = 0.d0
  zxir = 0.d0
  zxieta = 0.d0
  zetar = 0.d0
  zetaxi = 0.d0

! Initialize Flux_core 

  Flux_core = 0.d0

  ! Initialize Fluxes for diagnostic

  flux_r_out(:,:,:,:) = 0.d0
  flux_xi_xip(:,:,:,:) = 0.d0
  flux_eta_etap(:,:,:,:) = 0.d0

! Here we calculate the auxiliary variables for each thermal cell
! indexs i, j, k refers to the thermal grid

  maxbm = maxval(bm)

  !$OMP Parallel do reduction(+: Flux_core) private(p,i,j,k, ic, kc, jc, omegatau_out, kappa_perp_out, & 
  !$OMP & omegatau_xip, omegatau_etap, kappa_perp_xip, kappa_perp_etap)
  do k = 0, nangt
    do j = 0, nangt
      do p = 1, 6
        do i = 1, nrt

          ic = 2*i-1 ! center in the magnetic grid
          kc = 2*k   ! center in the magnetic grid
          jc = 2*j   ! center in the magnetic grid

        !write(*,*) i, j, k
        !write(*,*) "Br =", br(ic, jc, kc, p), "Bxi =", bxi(ic, jc, kc, p), &
        !&           "Beta=", beta(ic, jc, kc, p), "Bm=", bm(ic, jc, kc, p)

          if ((k /= 0) .and. (j /= 0)) then
            h(i, p, j, k) = dt/(cv(i,j,k,p)*vol(ic, jc, kc))
            ! units h = dt/(cv*vol) = Myr*1e8 K / 10**40 erg
          endif

          ! Average


          if (i == nrt) then
            omegatau_out = omegatau_arr(i,j,k,p)
            kappa_perp_out = kappa_perp_arr(i,j,k,p)
          else if (i == 1) then
            omegatau_out = omegatau_arr(i+1,j,k,p)
            kappa_perp_out = kappa_perp_arr(i+1,j,k,p)
          else
            omegatau_out = 0.5d0*(omegatau_arr(i, j, k, p) + omegatau_arr(i+1, j, k, p))
            if (omegatau_out /= omegatau_out) then
              print*,i,j,k,p,omegatau_out,omegatau_arr(i,j,k,p),omegatau_arr(i+1,j,k,p)
              stop
            endif
            kappa_perp_out = kappa_perp_arr(i,j,k,p)+ &
            & (kappa_perp_arr(i+1,j,k,p) - kappa_perp_arr(i,j,k,p))*(r(ic+1)-r(ic))/(r(ic+2)-r(ic))

          endif 

          omegatau_xip = 0.5d0*(omegatau_arr(i, j, k, p) + omegatau_arr(i, j+1, k, p))
          omegatau_etap = 0.5d0*(omegatau_arr(i, j, k, p) + omegatau_arr(i, j, k+1, p))

          ! Linear interpolation to obtain the value at the borders

          kappa_perp_xip = kappa_perp_arr(i,j,k,p) + &
          & (kappa_perp_arr(i,j+1,k,p) - kappa_perp_arr(i,j,k,p))*(xi(jc+1)-xi(jc))/(xi(jc+2)-xi(jc))

          kappa_perp_etap = kappa_perp_arr(i,j,k,p) + &
          & (kappa_perp_arr(i,j,k+1,p) - kappa_perp_arr(i,j,k,p))*(eta(kc+1)-eta(kc))/(eta(kc+2)-eta(kc))


          ! TEST --------------------

         ! if ((p==1 .and. j==nangt) .or. (p==2 .and. j==0))then
         !   omegatau_xip = 0.25d0*(omegatau_arr(i,nangt, k, 1) + omegatau_arr(i, nangt+1, k, 1) &
         !   & + omegatau_arr(i,0, k, 2) + omegatau_arr(i, 1, k, 2))
          !
!
!            kappa_perp_xip = 0.5d0*(kappa_perp_arr(i,nangt,k,1) + &
!            & (kappa_perp_arr(i,nangt+1,k,1) - &
!            &  kappa_perp_arr(i,nangt,k,1))*(xi(2*nangt+1)-xi(2*nangt))/(xi(2*nangt+2)-xi(2*nangt)) + &
!            & kappa_perp_arr(i,0,k,2) + &
!            & (kappa_perp_arr(i,1,k,2) - kappa_perp_arr(i,0,k,2))*(xi(1)-xi(0))/(xi(2)-xi(0)))
!
!          elseif((p==2 .and. j==nangt) .or. (p==3 .and. j==0)) then
!
!            omegatau_xip = 0.25d0*(omegatau_arr(i,nangt, k, 2) + omegatau_arr(i, nangt+1, k, 2) &
!            & + omegatau_arr(i,0, k, 3) + omegatau_arr(i, 1, k, 3))
!
!
!            kappa_perp_xip = 0.5d0*(kappa_perp_arr(i,nangt,k,2) + &
!            & (kappa_perp_arr(i,nangt+1,k,2) - &
!            &  kappa_perp_arr(i,nangt,k,2))*(xi(2*nangt+1)-xi(2*nangt))/(xi(2*nangt+2)-xi(2*nangt)) + &
!            & kappa_perp_arr(i,0,k,3) + &
!            & (kappa_perp_arr(i,1,k,3) - kappa_perp_arr(i,0,k,3))*(xi(1)-xi(0))/(xi(2)-xi(0)))
!
!          elseif((p==3 .and. j == nangt) .or. (p==4 .and. j==0)) then
!
!            omegatau_xip = 0.25d0*(omegatau_arr(i,nangt, k, 3) + omegatau_arr(i, nangt+1, k, 3) &
!            & + omegatau_arr(i,0, k, 4) + omegatau_arr(i, 1, k, 4))
!
!
!            kappa_perp_xip = 0.5d0*(kappa_perp_arr(i,nangt,k,3) + &
!            & (kappa_perp_arr(i,nangt+1,k,3) - &
!            &  kappa_perp_arr(i,nangt,k,3))*(xi(2*nangt+1)-xi(2*nangt))/(xi(2*nangt+2)-xi(2*nangt)) + &
!            & kappa_perp_arr(i,0,k,4) + &
!            & (kappa_perp_arr(i,1,k,4) - kappa_perp_arr(i,0,k,4))*(xi(1)-xi(0))/(xi(2)-xi(0)))
!
!          elseif((p==4 .and. j == nangt) .or. (p==1 .and. j==0)) then 
!
!            omegatau_xip = 0.25d0*(omegatau_arr(i,nangt, k, 4) + omegatau_arr(i, nangt+1, k, 4) &
!            & + omegatau_arr(i,0, k, 1) + omegatau_arr(i, 1, k, 1))
!
!
!            kappa_perp_xip = 0.5d0*(kappa_perp_arr(i,nangt,k,4) + &
!            & (kappa_perp_arr(i,nangt+1,k,4) - &
!            &  kappa_perp_arr(i,nangt,k,4))*(xi(2*nangt+1)-xi(2*nangt))/(xi(2*nangt+2)-xi(2*nangt)) + &
!            & kappa_perp_arr(i,0,k,1) + &
!            & (kappa_perp_arr(i,1,k,1) - kappa_perp_arr(i,0,k,1))*(xi(1)-xi(0))/(xi(2)-xi(0)))  
!
!          elseif((p==1 .and. k==nangt) .or. (p==5 .and. k==0)) then    
!            omegatau_etap = 0.25d0*(omegatau_arr(i,j,nangt,1) + omegatau_arr(i,j,nangt+1,1) &
!            & + omegatau_arr(i,j,0,5) + omegatau_arr(i,j,1,5))
!
!            kappa_perp_xip = 0.5d0*(kappa_perp_arr(i,j,nangt,1) + &
!            & (kappa_perp_arr(i,j,nangt+1,1) - &
!            &  kappa_perp_arr(i,j,nangt,1))*(eta(2*nangt+1)-eta(2*nangt))/(eta(2*nangt+2)-eta(2*nangt)) + &
!            & kappa_perp_arr(i,j,0,5) + &
!            & (kappa_perp_arr(i,j,1,5) - kappa_perp_arr(i,j,0,5))*(eta(1)-eta(0))/(eta(2)-eta(0)))
!
!          elseif((p==6 .and. k==nangt) .or. (p==1 .and. k==0)) then    
!            omegatau_etap = 0.25d0*(omegatau_arr(i,j,nangt,6) + omegatau_arr(i,j,nangt+1,6) &
!            & + omegatau_arr(i,j,0,1) + omegatau_arr(i,j,1,1))
!
!            kappa_perp_xip = 0.5d0*(kappa_perp_arr(i,j,nangt,6) + &
!            & (kappa_perp_arr(i,j,nangt+1,6) - &
!            &  kappa_perp_arr(i,j,nangt,6))*(eta(2*nangt+1)-eta(2*nangt))/(eta(2*nangt+2)-eta(2*nangt)) + &
!            & kappa_perp_arr(i,j,0,1) + &
!            & (kappa_perp_arr(i,j,1,1) - kappa_perp_arr(i,j,0,1))*(eta(1)-eta(0))/(eta(2)-eta(0)))
!
!          endif

          !--------------------------------

          if (maxval(bm) == 0d0) then
            kappa_rr_out(i, p, j, k) = A_kappa(ic+1, kappa_perp_out)
            kappa_xixi_xip(i, p, j, k) = B_kappa(ic, kc, kappa_perp_xip) 
            kappa_etaeta_etap(i, p, j, k) = G_kappa(ic, jc, kappa_perp_etap)
            kappa_xieta_xip(i, p, j, k) = E_kappa(ic, jc+1, kc, kappa_perp_xip)
            kappa_etaxi_etap(i, p, j, k) = F_kappa(ic, jc, kc+1, kappa_perp_etap)

            !--------- TEST ------------------------
            !if (i==nrt-2) then
            !  if ((p==1 .and. j==nangt) .or. (p==2 .and. j==0))then
            !   print*, p, k, 'B:', B_kappa(ic, kc, kappa_perp_xip),  &
            !    & 'E:', E_kappa(ic, jc+1, kc, kappa_perp_xip)
            !  endif
            !endif
          !--------------------------------- delate after


          else

            if (omegatau_out /= omegatau_out) then
              print*,"CK2",i,j,k,p,kappa_perp_out,omegatau_out,omegatau_arr(i,j,k,p),omegatau_arr(i+1,j,k,p)
              stop
            endif

          kappa_rr_out(i, p, j, k) = A_kappa(ic+1, kappa_perp_out) + &
          &            H_kappa(ic+1, p ,jc, kc, omegatau_out, kappa_perp_out)*br(ic+1, jc, kc, p)/bm(ic+1, jc, kc, p)

          kappa_xixi_xip(i, p, j, k) = B_kappa(ic, kc, kappa_perp_xip) + & 
          &              I_kappa(ic, p, jc+1, kc, omegatau_xip, kappa_perp_xip)*     &
          &              bxi(ic, jc+1, kc, p)/bm(ic, jc+1, kc, p)


          ! ------------TEST----------
          !if (i==nrt-2) then
          !  if ((p==1 .and. j==nangt) .or. (p==2 .and. j==0))then
          !    print*, p, k, 'B:', B_kappa(ic, kc, kappa_perp_xip), 'I:', I_kappa(ic, p, jc+1, kc, omegatau_xip, kappa_perp_xip), &
          !    & 'E:', E_kappa(ic, jc+1, kc, kappa_perp_xip), 'N:', &
          !    & N_kappa(ic, p, jc+1, kc, omegatau_xip, kappa_perp_xip), 'J:', &
          !    & J_kappa(ic, p, jc+1, kc, omegatau_xip, kappa_perp_xip), 'H:', &
          !    & H_kappa(ic, p, jc+1, kc, omegatau_xip, kappa_perp_xip), 'M:', &
          !    & M_kappa(ic, p, jc+1, kc, omegatau_xip, kappa_perp_xip)
          !  endif
          !endif

          ! --------------------------------- delate after

          kappa_etaeta_etap(i, p, j, k) = G_kappa(ic, jc, kappa_perp_etap) + & 
          &              J_kappa(ic, p, jc, kc+1, omegatau_etap, kappa_perp_etap)*   &
          &              beta(ic, jc, kc+1, p)/ bm(ic, jc, kc+1, p)

          kappa_rxi_out(i, p, j, k) = I_kappa(ic+1, p, jc, kc, omegatau_out, kappa_perp_out)* &
          &             br(ic+1, jc, kc, p)/bm(ic+1, jc, kc, p) +                             &
          &             K_kappa(ic+1, p, jc, kc, omegatau_out, kappa_perp_out) !0.d0

          kappa_reta_out(i, p, j, k) = J_kappa(ic+1, p, jc, kc, omegatau_out, kappa_perp_out)* &
          &             br(ic+1, jc, kc, p)/bm(ic+1, jc, kc, p) +                              &
          &             L_kappa(ic+1, p, jc, kc, omegatau_out, kappa_perp_out)!0.d0

          kappa_xir_xip(i, p, j, k) = H_kappa(ic, p, jc+1, kc, omegatau_xip, kappa_perp_xip)*  &
          &              bxi(ic, jc+1, kc, p)/bm(ic, jc+1, kc, p) +                            &
          &              M_kappa(ic, p, jc+1, kc, omegatau_xip, kappa_perp_xip)!0.d0

          kappa_xieta_xip(i, p, j, k) = E_kappa(ic, jc+1, kc, kappa_perp_xip) -                &
          &              N_kappa(ic, p, jc+1, kc, omegatau_xip, kappa_perp_xip) +              &
          &              J_kappa(ic, p, jc+1, kc, omegatau_xip, kappa_perp_xip)*               &
          &              bxi(ic, jc+1, kc, p)/bm(ic, jc+1, kc, p)!0.d0

          kappa_etar_etap(i, p, j, k) = H_kappa(ic, p, jc, kc+1, omegatau_etap, kappa_perp_etap)* &
          &              beta(ic, jc, kc+1, p)/bm(ic, jc, kc+1, p) +                              &
          &              O_kappa(ic, p, jc, kc+1, omegatau_etap, kappa_perp_etap)!0.d0

          kappa_etaxi_etap(i, p, j, k) = F_kappa(ic, jc, kc+1, kappa_perp_etap) +            &
          &              P_kappa(ic, p, jc, kc+1, omegatau_etap, kappa_perp_etap) +          &
          &              I_kappa(ic, p, jc, kc+1, omegatau_etap, kappa_perp_etap)*           &
          &              beta(ic, jc, kc+1, p)/bm(ic, jc, kc+1, p)!0.d0
          endif
          !----------------------------------------

          ! The areas here are the covariant components
          if (i .ne. nrt) then
            zrr(i, p, j, k) = enu(ic+1)*kappa_rr_out(i, p, j, k)*area_r(ic+1, jc, kc)/(r(ic+2) - r(ic))
          else
            zrr(i, p, j, k) = 0.d0 !null surface flux
          endif
          zxixi(i, p, j, k) = enu(ic)*kappa_xixi_xip(i, p, j, k)*area_xi(ic, jc+1, kc)/(xi(jc+2) - xi(jc))
          zetaeta(i, p, j, k) = enu(ic)*kappa_etaeta_etap(i, p, j, k)*area_eta(ic, jc, kc+1)/(eta(kc+2) - eta(kc))
          
          if (k .ne. 0) then
            !These terms with k = 0 are never used
            if (i .ne. nrt) then 
              zreta(i, p, j, k) = 0.25*enu(ic+1)*kappa_reta_out(i, p, j, k)*area_r(ic+1, jc, kc)/(eta(kc+1) - eta(kc-1))
            else 
              zreta(i, p, j, k) = 0.d0 !null surface flux
            endif
            zxieta(i, p, j, k) = 0.25*enu(ic)*kappa_xieta_xip(i, p, j, k)*area_xi(ic, jc+1, kc)/(eta(kc+1) - eta(kc-1))
          endif
          
          if (j .ne. 0) then
            !These terms with j = 0 are never used
            if (i .ne. nrt) then
              zrxi(i, p, j, k) = 0.25*enu(ic+1)*kappa_rxi_out(i, p, j, k)*area_r(ic+1, jc, kc)/(xi(jc+1) - xi(jc-1))
            else
              zrxi(i, p, j, k) = 0.d0 !null surface flux
            endif
            zetaxi(i, p, j, k) = 0.25*enu(ic)*kappa_etaxi_etap(i, p, j, k)*area_eta(ic, jc, kc+1)/(xi(jc+1) - xi(jc-1))
          endif
        
          if (i .ne. nrt) then
            zxir(i, p, j, k) = 0.25*enu(ic)*kappa_xir_xip(i, p, j, k)*area_xi(ic, jc+1, kc)/(r(ic+1) - r(ic-1))
            zetar(i, p, j, k) = 0.25*enu(ic)*kappa_etar_etap(i, p, j, k)*area_eta(ic, jc, kc+1)/(r(ic+1) - r(ic -1))
          else ! if we are at the last radial layer we use the center of the cell to estimate the derivative
          ! DV: Why? SA: because we don't have a i+1 cell
            ! NOTE: I set them to 0 just for debugging
            zxir(i, p, j, k) = 0.d0!0.25*enu(ic)*kappa_xir_xip(i, p, j, k)*area_xi(ic, jc+1, kc)/(r(ic) - r(ic-1))
            zetar(i, p, j, k) = 0.d0!0.25*enu(ic)*kappa_etar_etap(i, p, j, k)*area_eta(ic, jc, kc+1)/(r(ic) - r(ic -1))
          endif
          

          !Compute the flux between the core and the crust 
          if (i==1)then
            if ((j .ne. 0) .and. (k .ne. 0)) then ! exclude the ghost cells 
              Flux_core = Flux_core - zrr(i,p,j,k)*(temp(2,j,k,p)-T_core) - &
              &           zrxi(i,p,j,k)*(temp(2,j+1,k,p)-temp(2,j-1,k,p)) - &
              &           zreta(i,p,j,k)*(temp(2,j,k+1,p)-temp(2,j,k-1,p)) 
              !Flux_core = Flux_core + 1 
              !print*, Flux_core
            endif
          endif
          
          
         !Compute the fluxes for test


          if (i .ne. 1 .and. i .ne. nrt) then
            if (j .ne. 0 .and. k .ne.0) then ! exclude radial flux in ghost cell
              flux_r_out(i,p,j,k) = - zrr(i,p,j,k)*(temp(i+1,j,k,p) - temp(i,j,k,p))    &
              &                 - zrxi(i,p,j,k)*(temp(i+1, j+1, k, p) + temp(i, j+1, k, p) - temp(i+1, j-1, k, p) &
              &                 - temp(i,j-1,k,p))  &
              &                 - zreta(i,p,j,k)*(temp(i+1,j,k+1,p) + temp(i,j,k+1,p)  &
              &                 - temp(i+1,j,k-1,p) - temp(i,j,k-1,p))
            endif
          
            if (k .ne. 0) then ! exclude k=0 ghost layer
              !flux_xi_xip(i,p,j,k) = - zxir(i,p,j,k)*(temp(i+1,j+1,k,p) + temp(i+1,j,k,p)  &
              !&                       - temp(i-1,j+1,k,p) - temp(i-1,j,k,p))                  
              flux_xi_xip(i,p,j,k) = - zxir(i,p,j,k)*(temp(i+1,j+1,k,p) + temp(i+1,j,k,p)  &
              &                       - temp(i-1,j+1,k,p) - temp(i-1,j,k,p))                  &
              &                       - zxixi(i,p,j,k)*(temp(i,j+1,k,p) - temp(i,j,k,p))      &
              &                       - zxieta(i,p,j,k)*(temp(i,j+1,k+1,p) + temp(i,j,k+1,p)  &
              &                       - temp(i,j+1,k-1,p) - temp(i,j,k-1,p))
            endif

            if (j .ne. 0) then ! exclude j=0 ghost layer
              flux_eta_etap(i,p,j,k) = - zetar(i,p,j,k)*(temp(i+1,j,k+1,p) + temp(i+1,j,k,p)  &
              &                       - temp(i-1,j,k+1,p) - temp(i-1,j,k,p))  &
              &                       - zetaxi(i,p,j,k)*(temp(i,j+1,k+1,p) + temp(i,j+1,k,p)  &
              &                       - temp(i,j-1,k+1,p) - temp(i,j-1,k,p))                  &
              &                       - zetaeta(i,p,j,k)*(temp(i,j,k+1,p) - temp(i,j,k,p))
            endif
          endif

          if (CalculateFlux) then 
        
            if ((i .ne. nrt) .and. (j .ne. nangt) .and. (k .ne. nangt) .and. (i > 1) .and. (j > 1) .and. (k > 1)) then

              ! ******************    START DEPRECATED  *****************************
              ! Areas Missing
              !flux_r_out(i, p, j, k) = - kappa_rr_out(i, p, j, k)*(temp(i+1, j, k, p) - temp(i, j, k, p))/(r(ic+2)- r(ic)) + &
              !&                 - kappa_rxi_out(i, p, j, k)*(temp(i+1, j+1, k, p) + temp(i, j+1, k, p) - temp(i+1, j-1, k, p) + &
              !&                      - temp(i, j-1, k, p))/(4.*(xi(jc+1)-xi(jc-1))) + &
              !&                      - kappa_reta_out(i, p, j, k)*(temp(i+1, j, k+1, p) + temp(i, j, k+1, p) + &
              !&                      - temp(i+1, j, k-1, p) - temp(i, j, k-1, p))/(4.*(eta(kc+1)-eta(kc-1)))
              !
              !flux_xi_xip(i, p, j, k) = - kappa_xir_xip(i, p, j, k)*(temp(i+1, j+1, k, p) + temp(i+1, j, k, p) + &
              !&                       - temp(i-1, j+1, k, p) - temp(i-1, j, k, p))/(4.0*(r(ic+1)-r(ic-1))) + &
              !&                       - kappa_xixi_xip(i, p, j, k)*(temp(i, j+1, k, p) - temp(i, j, k, p))/(xi(jc+2)- xi(jc)) + &
              !&                       - kappa_xieta_xip(i, p, j, k)*(temp(i, j+1, k+1, p) + temp(i, j, k+1, p) + &
              !&                       - temp(i, j+1, k-1, p) - temp(i, j, k-1, p))/(4.0*(eta(kc+1)-eta(kc-1)))
              !
              !flux_eta_etap(i, p, j, k) = - kappa_etar_etap(i, p, j, k)*(temp(i+1, j, k+1, p) + temp(i+1, j, k, p) + &
              !&                       - temp(i-1, j, k+1, p) - temp(i-1, j, k, p))/(4.0*(r(ic+1)-r(ic-1))) + &
              !&                       - kappa_etaxi_etap(i, p, j, k)*(temp(i, j+1, k+1, p) + temp(i, j+1, k, p) + &
              !&                       - temp(i, j-1, k+1, p) - temp(i, j-1, k, p))/(4.0*(xi(jc+1)-xi(jc-1))) + &
              !&                       - kappa_etaeta_etap(i, p, j, k)*(temp(i, j, k+1, p) - temp(i, j, k, p))/(eta(kc+2)- eta(jc))

              ! ---- Now the following is for Perez-Azorin Test -------

              if (p == 5) then
                flux_xi_anl_xip(i,p,j,k) = - omegatau*dsin(theta(jc+1, kc, p))*D(kc)*Y(kc)/dsqrt(delta(jc+1, k)*(delta(jc+1, k)-1.))
                flux_eta_anl_etap(i,p,j,k) = omegatau*dsin(theta(jc, kc, p))*C(jc)*X(jc)/dsqrt(delta(jc, k+1)*(delta(jc, k+1)-1.))
              elseif(p == 6) then
                flux_xi_anl_xip(i,p,j,k) = omegatau*dsin(theta(jc+1, kc, p))*D(kc)*Y(kc)/dsqrt(delta(jc+1, k)*(delta(jc+1, k)-1.))
                flux_eta_anl_etap(i,p,j,k) = - omegatau*dsin(theta(jc, kc, p))*C(jc)*X(jc)/dsqrt(delta(jc, k+1)*(delta(jc, k+1)-1.))
              else
                  flux_xi_anl_xip(i,p,j,k) = omegatau*dsin(theta(jc+1, kc, p))*C(jc+1)*D(kc)/dsqrt(delta(jc+1, k))
                  flux_eta_anl_etap(i,p,j,k) = omegatau*dsin(theta(jc, kc, p))*X(jc)*Y(kc+1)/dsqrt(delta(jc, k+1))
              endif

              ! the temp_anlt is not assigned at the first time step!!!
              flux_r_anl_out(i,p,j,k) = - 0.25*(temp_anlt(i+1, j, k, p) + temp_anlt(i, j, k, p))*r(ic+1)/time
              flux_xi_anl_xip(i,p,j,k) = - 0.25*flux_xi_anl_xip(i,p,j,k)*(temp_anlt(i, j+1, k, p) + &
              &                            temp_anlt(i, j, k, p))*r(ic)/time
              flux_eta_anl_etap(i,p,j,k) = - 0.25*flux_eta_anl_etap(i,p,j,k)*(temp_anlt(i, j, k+1, p) + &
              &                             temp_anlt(i, j, k, p))*r(ic)/time

            else
              ! I set the flux at the patch boundary to 0, later I will change it

              flux_r_out(i, p, j, k) = 0.d0
              flux_xi_xip(i, p, j, k) = 0.d0
              flux_eta_etap(i, p, j, k) = 0.d0
            
              flux_r_anl_out(i, p, j, k) = 0.d0
              flux_xi_anl_xip(i, p, j, k) = 0.d0
              flux_eta_anl_etap(i, p, j, k) = 0.d0
            endif

          endif 

          !if (i == 2 .and. (p .ne. 5) .and. (p .ne. 6) .and. k ==4 .and. (i*j>0) .and. (i < nrt) .and. (j < nangt)) then
            !write(*,*) 'j =', j, 'p =', p, 'r =', r(ic), 'time =', time, 'T =', 0.5*(temp(i, j+1, k, p) + temp(i, j, k, p))
            !write(*,*) 'Flux_xi =', flux_xi_xip(i,p,j,k), 'Flux_anl_xi =', flux_xi_anl_xip(i,p,j,k)
            !write(*,*) 'error =', 100*(flux_xi_xip(i,p,j,k)/flux_xi_anl_xip(i,p,j,k) -1.0)
          !endif
          
        enddo
      enddo
    enddo
  enddo
  !$OMP end Parallel do



  ! ----------- Test for the fluxes ---------------------
  !print*, 'Flux test'
  !print*, 'k', 'Flux'
  !do k =1, nangt
  !!  print*, k, (zxixi(nrt-2, 1, nangt,k)-zxixi(nrt-2, 2, 0,k))/zxixi(nrt-2, 2, 0,k), &
  !!  & (kappa_xixi_xip(nrt-2, 1, nangt,k)-kappa_xixi_xip(nrt-2, 2, 0,k))/(kappa_xixi_xip(nrt-2, 2, 0,k))
  !!  print*, k, -zxir(28,1,nangt,k)*(temp(29,nangt+1,k,1)+temp(29,nangt,k,1)-temp(27,nangt+1,k,1)-temp(27,nangt,k,1)) - &
  !!  & zxixi(28,1,nangt,k)*(temp(28,nangt+1,k,1)-temp(28,nangt,k,1)) - &
  !!  & zxieta(28,1,nangt,k)*(temp(28,nangt+1,k+1,1)+temp(28,nangt,k+1,1)-temp(28,nangt+1,k,1)-temp(28,nangt,k,1)), &
  !!  & -zxir(28,2,0,k)*(temp(29,1,k,2)+temp(29,0,k,2)-temp(27,1,k,2)-temp(27,0,k,2)) - &
  !!  & zxixi(28,2,0,k)*(temp(28,1,k,2)-temp(28,0,k,2)) - &
  !!  & zxieta(28,2,0,k)*(temp(28,1,k+1,2)+temp(28,0,k+1,2)-temp(28,1,k,2)-temp(28,0,k,2))
  !  if (abs(flux_xi_xip(28,2,0,k))>0.d0) then 
  !    print*, k, flux_xi_xip(28,1,nangt,k),  flux_xi_xip(28,2,0,k), &
  !     &       (flux_xi_xip(28,1,nangt,k) - flux_xi_xip(28,2,0,k))/flux_xi_xip(28,2,0,k)
  !  else 
  !    print*, k, flux_xi_xip(28,1,nangt,k),  flux_xi_xip(28,2,0,k), &
  !    & kappa_xixi_xip(28,1,nangt,k)*(temp(28,nangt+1,k,1) - temp(28,nangt,k,1)) , &
  !    &kappa_xixi_xip(28,2,0,k)*(temp(28,1,k,2) - temp(28,0,k,2)) , &
  !    & kappa_xieta_xip(28,1,nangt,k)*(temp(28,nangt+1,k+1,1) + temp(28,nangt,k+1,1)  &
  !    &                       - temp(28,nangt+1,k-1,1) - temp(28,nangt,k-1,1)), &
  !    & kappa_xieta_xip(28,2,0,k)*(temp(28,1,k+1,2) + temp(28,0,k+1,2)  &
  !    &                       - temp(28,1,k-1,2) - temp(28,0,k-1,2))
  !  endif
  !enddo

  ! -----------------------------------------------------




  if (FirstCall) then
    matrix = 0.d0 ! initialize all the elements to 0
  endif
  source = 0.d0 ! Initialize the source to 0 each timestep 


  call envelope_model() ! we calculate the surface temperature, the flux and flux derivative


  ! Start filling the matrix in the interior's points, excluding the radial extremes, 1 and nrt+1
  !$OMP Parallel do private(p,i,j,k, ic, kc, jc, precoeff, nu_loss, nu_loss_deriv, Q, n)
  do i = 2, nrt !before nrt-1, now we include the last cell
    do p = 1, 6
      do k = 1, nangt
        do j = 1, nangt
          precoeff = dt/c_v ! we need to substitute c_v with an array

        
          ic = 2*i-1
          jc = 2*j
          kc = 2*k
              
          ! Here we have to fill the reduced matrix AB (called "matrix" here)
          ! given the full matrix A the rule in LAPACK is the following:
          ! AB(2*bandwd+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
          ! we have to consider however that here our matrix start from index 0 and not 1, so the previous
          ! formula is modified in the following:
          ! AB(2*bandwd+i-j,j) = A(i,j)
          !
          ! Here we do not start from a full matrix but we fill the reduced matrix directly
          ! NOTE: we could do both to double check: fill directly and fill the full matrix
          !       and than convert to the reduced one.



          ! ------------------------------- Sources ---------------------------------
          nu_loss = - q_neutrino(i,j,k,p) - q_joule(i,j,k,p) !- (temp(i, j, k, p)/T_0)**alpha
          nu_loss_deriv = - q_neutrino_der(i,j,k,p) !alpha*nu_loss/temp(i, j, k, p) !derivative with respect to temperature
          Q = (nu_loss*enu(ic) - nu_loss_deriv*temp(i, j, k, p))*enu(ic)

          if (CoolingOff) then
            !If CoolingOff is set true the cooling is switched off
            Q = 0.d0
            nu_loss_deriv = 0.d0
          endif
          ! -------------------------------------------------------------------------


          n = n_matrix(i, p, j, k)

          ! Diagonal elements
        
          matrix(2*bandwd,n) = 1.d0 + h(i, p, j, k)*(zrr(i,p, j, k) + zrr(i-1, p, j, k) +   &
         &                                    zxixi(i, p, j, k) + zxixi(i, p, j-1, k)   +   &
         &                                    zetaeta(i, p, j,k) + zetaeta(i, p, j, k-1))   &
         &       - nu_loss_deriv*enu(ic)*dt/cv(i,j,k,p) ! source term

         

         if(matrix(2*bandwd,n) > 1.d50) then ! FOR DEBUGGING
            write(*,*) "Tevol debug", i, j, k, p, matrix(2*bandwd,n), h(i,p,j,k), cv(i,j,k,p)
         endif
         
         if (i == nrt) then ! if we are at the surface we add the linearization of the BB term
            matrix(2*bandwd,n) = matrix(2*bandwd,n) + &
            &  h(i, p, j, k)*cfluxb(j,k,p) + & !linearization of the stiff term 
            &  h(i, p, j, k)*(-zxir(i, p, j, k) + zxir(i, p, j-1,k) - zetar(i, p, j, k) + & !transverse flux at the border
            & zetar(i, p, j, k-1))

         endif 

          

        ! condition = dcos(phi(jc,kc,p))*dsin(theta(jc,kc,p))*dcos(long)*dsin(lat) + &
        ! & dsin(phi(jc,kc,p))*dsin(theta(jc,kc,p))*dsin(lat)*dsin(long) + &
         !& dcos(theta(jc,kc,p))*dcos(lat)
         
        
         !if ((condition .ge. dcos(ang)) .and. (r(ic) .le. rmax/2.0)) then
         !   Q = 100.0*exp(-(time+dt)/1.5)
         !endif

          ! Source for the internal patches
          source(n) = temp(i, j, k, p) + Q*dt/cv(i,j,k,p)

          ! First neighbours
          matrix(2*bandwd+n-n_matrix(i-1, p, j, k), n_matrix(i-1, p, j, k)) = h(i, p, j, k)*(-zrr(i-1, p, j, k) + &
         &                              zxir(i, p, j, k) - &
         &                              zxir(i, p, j-1, k) + zetar(i, p, j, k) - zetar(i, p, j, k-1))           !i-1 element

          if (i .ne. nrt) then
            matrix(2*bandwd+n-n_matrix(i+1, p, j, k), n_matrix(i+1, p, j, k)) = h(i, p, j, k)*(-zrr(i, p, j, k) - &
            &                              zxir(i, p, j, k) + &
            &                              zxir(i, p, j-1, k) - zetar(i, p, j, k) + zetar(i, p, j, k-1))           !i+1 element
          endif

          matrix(2*bandwd+n-n_matrix(i, p, j-1, k), n_matrix(i, p, j-1, k)) = h(i, p, j, k)*(-zxixi(i, p, j-1, k) + &
         &                               zrxi(i, p, j, k) - &
         &                               zrxi(i-1, p, j, k) + zetaxi(i, p, j, k) - zetaxi(i, p, j, k-1))        !j-1 element

        

          matrix(2*bandwd+n-n_matrix(i, p, j, k-1), n_matrix(i, p, j, k-1)) = h(i, p, j, k)*(- zetaeta(i, p, j, k-1) + &
         &                                   zreta(i, p, j, k) - &
         &                                   zreta(i-1, p, j, k) + zxieta(i, p, j, k) - zxieta(i, p, j-1, k))   !k-1 element

          matrix(2*bandwd+n-n_matrix(i, p, j, k+1), n_matrix(i, p, j, k+1)) = h(i, p, j, k)*(-zetaeta(i, p, j, k) - &
         &                                     zreta(i, p, j, k) + &
         &                                     zreta(i-1, p, j, k) - zxieta(i, p, j, k) + zxieta(i, p, j-1, k)) !k+1 element

          matrix(2*bandwd+n-n_matrix(i, p, j+1, k), n_matrix(i, p, j+1, k)) = h(i, p, j, k)*(-zxixi(i, p, j, k) - &
         &                                 zrxi(i, p, j, k) + &
         &                                 zrxi(i-1, p, j, k) - zetaxi(i, p, j, k) + zetaxi(i, p, j, k-1))      !j+1 element

         if (i == nrt) then ! correct the first neighbours for transverse flux at the outer border
            matrix(2*bandwd+n-n_matrix(i, p, j-1, k), n_matrix(i, p, j-1, k)) =                                     &
            & matrix(2*bandwd+n-n_matrix(i, p, j-1, k), n_matrix(i, p, j-1, k)) + h(i, p, j, k)*zxir(i, p, j-1, k)

            matrix(2*bandwd+n-n_matrix(i, p, j+1, k), n_matrix(i, p, j+1, k)) =                                     &
            & matrix(2*bandwd+n-n_matrix(i, p, j+1, k), n_matrix(i, p, j+1, k)) - h(i, p, j, k)*zxir(i, p, j, k)

            matrix(2*bandwd+n-n_matrix(i, p, j, k-1), n_matrix(i, p, j, k-1)) =                                     &
            & matrix(2*bandwd+n-n_matrix(i, p, j, k-1), n_matrix(i, p, j, k-1)) + h(i, p, j, k)*zetar(i, p, j, k-1)

            matrix(2*bandwd+n-n_matrix(i, p, j, k+1), n_matrix(i, p, j, k+1)) =                                     &
            & matrix(2*bandwd+n-n_matrix(i, p, j, k+1), n_matrix(i, p, j, k+1)) - h(i, p, j, k)*zetar(i, p, j, k)
          endif


          ! Second neighbours
          matrix(2*bandwd+n-n_matrix(i, p, j-1, k-1), n_matrix(i, p, j-1, k-1)) = h(i, p, j, k)*(-zxieta(i, p, j-1, k) - &
          &                                             zetaxi(i, p, j, k-1))  !j-1, k-1 element
          matrix(2*bandwd+n-n_matrix(i, p, j-1, k+1), n_matrix(i, p, j-1, k+1)) = h(i, p, j, k)*(zetaxi(i, p, j, k) + &
          &                                             zxieta(i, p, j-1, k))     !j-1, k+1 element
          matrix(2*bandwd+n-n_matrix(i-1, p, j-1, k), n_matrix(i-1, p, j-1, k)) = h(i, p, j, k)*(-zrxi(i-1, p, j, k) - &
          &                                             zxir(i, p, j-1, k))      !j-1, i-1 element
          if (i .ne. nrt) then
            matrix(2*bandwd+n-n_matrix(i+1, p, j-1, k), n_matrix(i+1, p, j-1, k)) = h(i, p, j, k)*(zrxi(i, p, j, k) + &
            &                                             zxir(i, p, j-1, k))         !j-1, i+1 element
          endif
          matrix(2*bandwd+n-n_matrix(i, p, j+1, k-1), n_matrix(i, p, j+1, k-1)) = h(i, p, j, k)*(zxieta(i, p, j, k) + &
          &                                             zetaxi(i, p, j, k-1))     !j+1, k-1 element
          matrix(2*bandwd+n-n_matrix(i, p, j+1, k+1), n_matrix(i, p, j+1, k+1)) = h(i, p, j, k)*(-zxieta(i, p, j, k) - &
          &                                             zetaxi(i, p, j, k))      !j+1, k+1 element
          matrix(2*bandwd+n-n_matrix(i-1, p, j+1, k), n_matrix(i-1, p, j+1, k)) = h(i, p, j, k)*(zxir(i, p, j, k) + &
          &                                             zrxi(i-1, p, j, k))         !j+1, i-1 element
          if (i .ne. nrt) then
            matrix(2*bandwd+n-n_matrix(i+1, p, j+1, k), n_matrix(i+1, p, j+1, k)) = h(i, p, j, k)*(-zrxi(i, p, j, k) - &
            &                                             zxir(i, p, j, k))          !j+1, i+1 element
          endif
          matrix(2*bandwd+n-n_matrix(i-1, p, j, k-1), n_matrix(i-1, p, j, k-1)) = h(i, p, j, k)*(-zreta(i-1, p, j, k) - &
          &                                             zetar(i, p, j, k-1))    !k-1, i-1 element
          if (i .ne. nrt) then
            matrix(2*bandwd+n-n_matrix(i+1, p, j, k-1), n_matrix(i+1, p, j, k-1)) = h(i, p, j, k)*(zreta(i, p, j, k) + &
            &                                             zetar(i, p, j, k-1))       !k-1, i+1 element
          endif 
          matrix(2*bandwd+n-n_matrix(i-1, p, j, k+1), n_matrix(i-1, p, j, k+1)) = h(i, p, j, k)*(zetar(i, p, j, k) + &
          &                                             zreta(i-1, p, j, k))       !k+1, i-1 element
          if (i .ne. nrt) then
            matrix(2*bandwd+n-n_matrix(i+1, p, j, k+1), n_matrix(i+1, p, j, k+1)) = h(i, p, j, k)*(-zreta(i, p, j, k) - &
            &                                             zetar(i, p, j, k))        !k+1, i+1 element
          endif

          ! Check the identities

          !identity1 = matrix(2*bandwd,n) +                                                      &
          !&          matrix(2*bandwd+n-n_matrix(i-1, p, j, k), n_matrix(i-1, p, j, k)) +       &
          !&          matrix(2*bandwd+n-n_matrix(i+1, p, j, k), n_matrix(i+1, p, j, k)) +       &
          !&          matrix(2*bandwd+n-n_matrix(i, p, j-1, k), n_matrix(i, p, j-1, k)) +       &
          !&          matrix(2*bandwd+n-n_matrix(i, p, j+1, k), n_matrix(i, p, j+1, k)) +       &
          !&          matrix(2*bandwd+n-n_matrix(i, p, j, k-1), n_matrix(i, p, j, k-1)) +       &
          !&          matrix(2*bandwd+n-n_matrix(i, p, j, k+1), n_matrix(i, p, j, k+1))

          !identity2 = matrix(2*bandwd+n-n_matrix(i, p, j-1, k-1), n_matrix(i, p, j-1, k-1)) +   &
          !&          matrix(2*bandwd+n-n_matrix(i, p, j-1, k+1), n_matrix(i, p, j-1, k+1)) +   &
          !&          matrix(2*bandwd+n-n_matrix(i-1, p, j-1, k), n_matrix(i-1, p, j-1, k)) +   &
          !&          matrix(2*bandwd+n-n_matrix(i+1, p, j-1, k), n_matrix(i+1, p, j-1, k)) +   &
          !&          matrix(2*bandwd+n-n_matrix(i, p, j+1, k-1), n_matrix(i, p, j+1, k-1)) +   &
          !&          matrix(2*bandwd+n-n_matrix(i, p, j+1, k+1), n_matrix(i, p, j+1, k+1)) +   &
          !&          matrix(2*bandwd+n-n_matrix(i-1, p, j+1, k), n_matrix(i-1, p, j+1, k)) +   &
          !&          matrix(2*bandwd+n-n_matrix(i+1, p, j+1, k), n_matrix(i+1, p, j+1, k)) +   &
          !&          matrix(2*bandwd+n-n_matrix(i-1, p, j, k-1), n_matrix(i-1, p, j, k-1)) +   &
          !&          matrix(2*bandwd+n-n_matrix(i+1, p, j, k-1), n_matrix(i+1, p, j, k-1)) +   &
          !&          matrix(2*bandwd+n-n_matrix(i-1, p, j, k+1), n_matrix(i-1, p, j, k+1)) +   &
          !&          matrix(2*bandwd+n-n_matrix(i+1, p, j, k+1), n_matrix(i+1, p, j, k+1))

          !if (dabs(identity1 - 1.d0) > 1.e-5 .or. dabs(identity2) > 1.e-5) then
          !  write(*,*) '************************ WARNING *******************************************'
          !  write(*,*) 'One or both the identities are not satisfied'
          !  write(*,*) '(p, i, j, k) =', '(',p, i, j, k,')'
          !  write(*,*) 'First identity =', identity1
          !  write(*,*) 'Second identity =', identity2
          !endif

        enddo
      enddo
    enddo
  enddo
  !$OMP end Parallel do

  !--------------------------------------- Boundary Conditions --------------------------------------------------
  !   -------------------------------- Angular Boundary Conditions -------------------------------------------
  !
  ! We cycle over i and p.
  ! We the index of the temperature that we have to fill for the border j,k = 0 (j,k = nangt+1) and
  ! we call them n_jinf, n_jsup, n_kinf, n_ksup, where inf and sup stand for 0 and nangt+1.
  ! We identify the index of the two coupled elements as n1_jinf, n2_jinf etc. where inf and sup
  ! stands for the border j,k = 0, or j,k = nangt+1, while the index 1 indicates that we are considering
  ! the temperature with the same coordinate of the cell in the border and the index 2 the temperature
  ! just one cell towards the centre (identified as kprime).
  ! The values of the indexes depend on the patch, see the notes.

  
  if (FirstCall) then
    !do i =2, nrt !before i=2, nrt-1 !check I included the last layer
    !do i =2, nrt ! MAYBE USE THIS FOR PEREZ-AZORIN TEST(exclud first layer or maybe also the last) 
    do i=1,nrt!i =2, nrt
      do p = 1, 6
        do k = 1, nangt   ! CHECK if the interval is right. Before it was k = 2, nangt-1
        
          if (k <= nangt/2) then
            kprime = k + 1
          else
            kprime = k - 1
          endif

          ! We determine the matrix row we have to fill number
          n_jinf = n_matrix(i, p, 0, k)        ! For j = 0
          n_jsup = n_matrix(i, p, nangt+1, k)  ! For j = nangt+1
          n_kinf = n_matrix(i, p, k, 0)        ! For k = 0
          n_ksup = n_matrix(i, p, k, nangt+1)  ! For k = nangt+1

          source(n_jinf) = 0.d0 !source always 0 in the ghost cells
          source(n_jsup) = 0.d0
          source(n_kinf) = 0.d0
          source(n_ksup) = 0.d0
            
          if (p == 1) then
            n1_jinf = n_matrix(i, 4, nangt, k)
            n2_jinf = n_matrix(i, 4, nangt, kprime)
            n1_jsup = n_matrix(i, 2, 1, k)
            n2_jsup = n_matrix(i, 2, 1, kprime)
            n1_kinf = n_matrix(i, 6, k, nangt)
            n2_kinf = n_matrix(i, 6, kprime, nangt)
            n1_ksup = n_matrix(i, 5, k, 1)
            n2_ksup = n_matrix(i, 5, kprime, 1)

          elseif (p == 2) then
            n1_jinf = n_matrix(i, 1, nangt, k)
            n2_jinf = n_matrix(i, 1, nangt, kprime)
            n1_jsup = n_matrix(i, 3, 1, k)
            n2_jsup = n_matrix(i, 3, 1, kprime)
            n1_kinf = n_matrix(i, 6, nangt, nangt+1-k)
            n2_kinf = n_matrix(i, 6, nangt, nangt+1-kprime)
            n1_ksup = n_matrix(i, 5, nangt, k)
            n2_ksup = n_matrix(i, 5, nangt, kprime)
 
          elseif (p == 3) then
            n1_jinf = n_matrix(i, 2, nangt, k)
            n2_jinf = n_matrix(i, 2, nangt, kprime)
            n1_jsup = n_matrix(i, 4, 1, k)
            n2_jsup = n_matrix(i, 4, 1, kprime)
            n1_kinf = n_matrix(i, 6, nangt+1-k, 1)
            n2_kinf = n_matrix(i, 6, nangt+1-kprime, 1)
            n1_ksup = n_matrix(i, 5, nangt+1-k, nangt)
            n2_ksup = n_matrix(i, 5, nangt+1-kprime, nangt)

          elseif (p == 4) then
            n1_jinf = n_matrix(i, 3, nangt, k)
            n2_jinf = n_matrix(i, 3, nangt, kprime)
            n1_jsup = n_matrix(i, 1, 1, k)
            n2_jsup = n_matrix(i, 1, 1, kprime)
            n1_kinf = n_matrix(i, 6, 1, k)
            n2_kinf = n_matrix(i, 6, 1, kprime)
            n1_ksup = n_matrix(i, 5, 1, nangt+1-k)
            n2_ksup = n_matrix(i, 5, 1, nangt+1-kprime)
      
          elseif (p == 5) then
            n1_jinf = n_matrix(i, 4, nangt+1-k, nangt)
            n2_jinf = n_matrix(i, 4, nangt+1-kprime, nangt)
            n1_jsup = n_matrix(i, 2, k, nangt)
            n2_jsup = n_matrix(i, 2, kprime, nangt)
            n1_kinf = n_matrix(i, 1, k, nangt)
            n2_kinf = n_matrix(i, 1, kprime, nangt)
            n1_ksup = n_matrix(i, 3, nangt+1-k, nangt)
            n2_ksup = n_matrix(i, 3, nangt+1-kprime, nangt)

          else
            n1_jinf = n_matrix(i, 4, k, 1)
            n2_jinf = n_matrix(i, 4, kprime, 1)
            n1_jsup = n_matrix(i, 2, nangt+1-k, 1)
            n2_jsup = n_matrix(i, 2, nangt+1-kprime, 1)
            n1_kinf = n_matrix(i, 3, nangt+1-k, 1)
            n2_kinf = n_matrix(i, 3, nangt+1-kprime, 1)
            n1_ksup = n_matrix(i, 1, k, 1)
            n2_ksup = n_matrix(i, 1, kprime, 1)

          endif  ! End of patch cases

          ! We assign the value to the matrix elements
          kc = 2*k
          ! jt = 0
          matrix(2*bandwd, n_jinf) = 1.d0
          matrix(2*bandwd+n_jinf-n1_jinf, n1_jinf) = - (1.d0 - edge_wt(kc))!- (1.d0 - 0.5*edge_wt(kc)) !- (1.d0 - 0.5*edge_w(kc))
          matrix(2*bandwd+n_jinf-n2_jinf, n2_jinf) = - edge_wt(kc)!- 0.5*edge_wt(kc) !- 0.5*edge_w(kc)
          ! jt = nangt + 1
          matrix(2*bandwd, n_jsup) = 1.d0
          matrix(2*bandwd+n_jsup-n1_jsup, n1_jsup) = - (1.d0 - edge_wt(kc))!- (1.d0 - 0.5*edge_wt(kc)) !- (1.d0 - 0.5*edge_w(kc))
          matrix(2*bandwd+n_jsup-n2_jsup, n2_jsup) = - edge_wt(kc)!- 0.5*edge_wt(kc) !- 0.5*edge_w(kc)
          ! kt = 0
          matrix(2*bandwd, n_kinf) = 1.d0
          matrix(2*bandwd+n_kinf-n1_kinf, n1_kinf) = - (1.d0 - edge_wt(kc))!- (1.d0 - 0.5*edge_wt(kc)) !- (1.d0 - 0.5*edge_w(kc))
          matrix(2*bandwd+n_kinf-n2_kinf, n2_kinf) = - edge_wt(kc)!- 0.5*edge_wt(kc) !- 0.5*edge_w(kc)
          ! kt = nangt + 1
          matrix(2*bandwd, n_ksup) = 1.d0
          matrix(2*bandwd+n_ksup-n1_ksup, n1_ksup) = - (1.d0 - edge_wt(kc))!- (1.d0 - 0.5*edge_wt(kc)) !- (1.d0 - 0.5*edge_w(kc))
          matrix(2*bandwd+n_ksup-n2_ksup, n2_ksup) = - edge_wt(kc)!- 0.5*edge_wt(kc) !- 0.5*edge_w(kc)

          !write(*,*) n_jinf, n1_jinf, n2_jinf
          !write(*,*) n_jsup, n1_jsup, n2_jsup
          !write(*,*) matrix(2*bandwd, n_jinf), matrix(2*bandwd+n_jinf-n1_jinf, n1_jinf), matrix(2*bandwd+n_jinf-n2_jinf, n2_jinf)
          !write(*,*) matrix(2*bandwd, n_jsup), matrix(2*bandwd+n_jsup-n1_jsup, n1_jsup), matrix(2*bandwd+n_jsup-n2_jsup, n2_jsup)
          !write(*,*) "*****************************************************"

        enddo  ! End of k loop



        ! Interpolation constraints at the four corners
        ! Hereafter the indexes of the corners are:
        ! (1) for jt=0, kt=0
        ! (2) for jt=nangt+1, kt=0
        ! (3) for jt=0, kt=nangt+1
        ! (4) for jt=nangt+1, kt=nangt+1
        if (p==1) then
          n1_corner(1) = n_matrix(i, 4, nangt, 1)
          n2_corner(1) = n_matrix(i, 6, 1, nangt)
          n1_corner(2) = n_matrix(i, 2, 1, 1)
          n2_corner(2) = n_matrix(i, 6, nangt, nangt)
          n1_corner(3) = n_matrix(i, 4, nangt, nangt)
          n2_corner(3) = n_matrix(i, 5, 1, 1)
          n1_corner(4) = n_matrix(i, 2, 1, nangt)
          n2_corner(4) = n_matrix(i, 5, nangt, 1)

        elseif(p==2) then
          n1_corner(1) = n_matrix(i, 1, nangt, 1)
          n2_corner(1) = n_matrix(i, 6, nangt, nangt)
          n1_corner(2) = n_matrix(i, 3, 1, 1)
          n2_corner(2) = n_matrix(i, 6, nangt, 1)
          n1_corner(3) = n_matrix(i, 1, nangt, nangt)
          n2_corner(3) = n_matrix(i, 5, nangt, 1)
          n1_corner(4) = n_matrix(i, 3, 1, nangt)
          n2_corner(4) = n_matrix(i, 5, nangt, nangt)

        elseif(p==3) then
          n1_corner(1) = n_matrix(i, 2, nangt, 1)
          n2_corner(1) = n_matrix(i, 6, nangt, 1)
          n1_corner(2) = n_matrix(i, 4, 1, 1)
          n2_corner(2) = n_matrix(i, 6, 1, 1)
          n1_corner(3) = n_matrix(i, 2, nangt, nangt)
          n2_corner(3) = n_matrix(i, 5, nangt, nangt)
          n1_corner(4) = n_matrix(i, 4, 1, nangt)
          n2_corner(4) = n_matrix(i, 5, 1, nangt)!n_matrix(i, 5, nangt, nangt)<-error

        elseif(p==4) then
          n1_corner(1) = n_matrix(i, 3, nangt, 1)
          n2_corner(1) = n_matrix(i, 6, 1, 1)
          n1_corner(2) = n_matrix(i, 1, 1, 1)
          n2_corner(2) = n_matrix(i, 6, 1, nangt)
          n1_corner(3) = n_matrix(i, 3, nangt, nangt)
          n2_corner(3) = n_matrix(i, 5, 1, nangt)
          n1_corner(4) = n_matrix(i, 1, 1, nangt)
          n2_corner(4) = n_matrix(i, 5, 1, 1)

        elseif(p==5) then
          n1_corner(1) = n_matrix(i, 4, nangt, nangt)
          n2_corner(1) = n_matrix(i, 1, 1, nangt)
          n1_corner(2) = n_matrix(i, 2, 1, nangt)
          n2_corner(2) = n_matrix(i, 1, nangt, nangt)
          n1_corner(3) = n_matrix(i, 4, 1, nangt)
          n2_corner(3) = n_matrix(i, 3, nangt, nangt)
          n1_corner(4) = n_matrix(i, 2, nangt, nangt)
          n2_corner(4) = n_matrix(i, 3, 1, nangt)

        else
          n1_corner(1) = n_matrix(i, 4, 1, 1)
          n2_corner(1) = n_matrix(i, 3, nangt, 1)
          n1_corner(2) = n_matrix(i, 2, nangt, 1)
          n2_corner(2) = n_matrix(i, 3, 1, 1)
          n1_corner(3) = n_matrix(i, 4, nangt, 1)
          n2_corner(3) = n_matrix(i, 1, 1, 1)
          n1_corner(4) = n_matrix(i, 2, 1, 1)
          n2_corner(4) = n_matrix(i, 1, nangt, 1)
        endif

        matrix(2*bandwd, n_matrix(i, p, 0, 0)) = 1.d0
        matrix(2*bandwd, n_matrix(i, p, nangt+1, 0)) = 1.d0
        matrix(2*bandwd, n_matrix(i, p, 0, nangt+1)) = 1.d0
        matrix(2*bandwd, n_matrix(i, p, nangt+1, nangt+1)) = 1.d0
     
        ! All the neighbours have weigth 1/2 by construction (see notes)
        matrix(2*bandwd+n_matrix(i, p, 0, 0)-n1_corner(1), n1_corner(1)) = -0.5d0
        matrix(2*bandwd+n_matrix(i, p, 0, 0)-n2_corner(1), n2_corner(1)) = -0.5d0
        matrix(2*bandwd+n_matrix(i, p, nangt+1, 0)-n1_corner(2), n1_corner(2)) = -0.5d0
        matrix(2*bandwd+n_matrix(i, p, nangt+1, 0)-n2_corner(2), n2_corner(2)) = -0.5d0
        matrix(2*bandwd+n_matrix(i, p, 0, nangt+1)-n1_corner(3), n1_corner(3)) = -0.5d0
        matrix(2*bandwd+n_matrix(i, p, 0, nangt+1)-n2_corner(3), n2_corner(3)) = -0.5d0
        matrix(2*bandwd+n_matrix(i, p, nangt+1, nangt+1)-n1_corner(4), n1_corner(4)) = -0.5d0
        matrix(2*bandwd+n_matrix(i, p, nangt+1, nangt+1)-n2_corner(4), n2_corner(4)) = -0.5d0

        !source corner (source is 0 in the ghost-cell!)

        source(n_matrix(i, p, 0, 0)) = 0.d0
        source(n_matrix(i, p, nangt+1, 0)) = 0.d0
        source(n_matrix(i, p, 0, nangt+1)) = 0.d0
        source(n_matrix(i, p, nangt+1, nangt+1)) = 0.d0
      
      
      enddo
    enddo
  endif
  

  ! Radial Boundary Conditions


  call core_cooling(T_core, dt, Flux_core) !This routine evolve the temperature in the core
  !print*, time*1.d6, T_core, Flux_core, temp(2,2,2,2), zrr(1,2,2,2)*(temp(2,2,2,2) - T_core)

  !$OMP Parallel do private(p,j,k, kc, jc, n_jinf, n_jsup)  collapse(2)
  do p = 1, 6
    do k = 0, nangt+1
      do j = 0, nangt+1

        jc = 2*j
        kc = 2*k

        n_jinf = n_matrix(1, p, j, k)
        n_jsup = n_matrix(nrt, p, j, k)

        matrix(2*bandwd, n_jinf) = 1.d0


        ! --------------------  Realistic BC  --------------------------
        ! Here we consider the cooling of the core and the BB cooling from the surface.
        ! we also call the envelope model
     
        
        if (k .ne. 0 .and. j .ne. 0 .and. j .ne.  nangt+1 .and. k .ne. nangt+1) then

          ! units h = dt/(cv*vol) = Myr*1e8 K / 10**40 erg
          ! units source: Qnu*dt/cv : 1e8 K
          ! units of sfluxb have to be then: 1e40 erg/Myr

          !if (k == 2 .and. j == 2 .and. p ==2 ) then !DEBUGGING
          !  print*, time, h(nrt,p,j,k), sfluxb(j,k,p)*UNIT_EN/UNIT_TIME, h(nrt,p,j,k)*sfluxb(j,k,p)
          !endif 

          source(n_jsup) = source(n_jsup) - h(nrt,p,j,k)*sfluxb(j,k,p)
          source(n_jinf) = T_core !inner cooling core

        else
            source(n_jsup) = 0.d0
            source(n_jinf) = 0.d0
        endif 

        

        ! ------ Perez-Azorin 2006 test ----------
        ! set omegatau = 0 for isotropic case

        !ang_term = dsin(theta(jc, kc, p))**2 + (dcos(theta(jc, kc, p))**2)/(1.d0 + omegatau*omegatau)
        
        ! Generalized PA test
        !cos_th_prime = dsin(beta_ang)*dsin(theta(jc, kc ,p))*(-dcos(gamma_ang)*dcos(phi(jc, kc, p)) + &
        !&              dsin(gamma_ang)*dsin(phi(jc, kc, p))) + dcos(beta_ang)*dcos(theta(jc, kc, p))
        !ang_term = 1.d0 - cos_th_prime*cos_th_prime*omegatau*omegatau/(1.d0 + omegatau*omegatau)

        !source(n_jinf) = T_int*((1./((time+dt)))**(3./2.))*exp(-r(1)*r(1)*ang_term/(4.*kappa_perp*(time+dt)))!Perez Azorin
        !source(n_jsup) = T_int*((1./((time+dt)))**(3./2.))*exp(-r(2*nrt-1)*r(2*nrt-1)*ang_term/(4.*kappa_perp*(time+dt)))!Perez Azorin

        ! --- Hot Spot Test -------
        !source(n_jinf) = T_ext
        !source(n_jsup) = T_ext
        ! ----------------------------

       if (k == 0 .or. j == 0 .or. j ==  nangt+1 .or. k == nangt+1) then
            source(n_jinf) = 0.d0
            source(n_jsup) = 0.d0
        endif

      enddo
    enddo
  enddo
  !$OMP end Parallel do

  !--------------------------------------------------------------------------------------------------------------
  ! End Boundary Conditions -------------------------------------------------------------------------------------

  ! end filling the matrix and source
  
  ! call inversion(...., x) ! Public library? No priority now
  !call dgbsv(N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB, INFO)
  ! N    = number of columns
  ! KL   = lower band dimension
  ! KU   = upper band dimension
  ! NRHS = number of column of the matrix B
  ! AB   = reduced matrix
  ! LDAB = reduced matrix row number
  ! IPIV = pivot indices that define the permutation matrix P
  ! B    = source matrix (right hand side)
  ! LDB  = lead dimension of B
  ! INFO = success/not success, 0 if succesfull

!  do i = 0, 3*bandwd
!    write(100,*) matrix(i,:)
!  enddo

  !write(*,*) "Source Before"
  !if (FirstCall) then
  !  do i = 0,6*nrt*(nangt+2)*(nangt+2)-1
  !      write(102, *) source(i)
        !write(*,*) i, source(i)
  !  enddo
  !endif

  

  !if (FirstCall) then
  !  do p = 1, 6
  !     write (name,'(a12,i1,a2)') "out/temp_ini", p, ".d"
  !     open(unit = 120, file=trim(name), status = 'replace')
  !     
  ! !open(unit=115,file=trim(name),status='replace')
  !
  !     do j = 1, nangt
  !         do k =1, nangt
  !             jc = 2*j
  !             kc = 2*k
  !             !ang_term = (dsin(theta(jc, kc, p))**2 + (dcos(theta(jc, kc, p))**2)/(1.d0 + omegatau*omegatau))
  !             
  !
  !             ! generlized PA solution
  !              cos_th_prime = dsin(beta_ang)*dsin(theta(j, k ,p))*(-dcos(gamma_ang)*dcos(phi(j, k, p)) + &
  !              &              dsin(gamma_ang)*dsin(phi(j, k, p))) + dcos(beta_ang)*dcos(theta(j, k, p))
  !              ang_term = 1.d0 - cos_th_prime*cos_th_prime*omegatau*omegatau/(1.d0 + omegatau*omegatau)
  !             T_analit = T_int*((1./(time+dt))**(3./2.))!*exp(-r(2)*r(2)*ang_term/(4.*kappa_perp*(time+dt))) !Perez-Azorin
  !                !T_analit = 0.d0
  !
  !             write(120, *) theta(jc, kc, p), phi(jc, kc, p), (source(n_matrix(2,p, j, k)) - T_analit)/T_analit, &
  !             &     source(n_matrix(2,p, j, k))
  !         enddo
  !     enddo
  !     close(120)
  ! enddo
  !endif
  ! Open a file for the L2 measure

  ! *************    PEREZ-AZORIN TEST ----- WRITE L2 Error    ****************
  !write (name,'(a10)') "out/L2_T.d"
  !if (FirstCall) then
  !  open(unit = 130, file=trim(name), status = 'replace')

    ! we write the analitic temperature at the first timestep 
    ! WARNING The following loop seems pretty useless and it's repeted after matrix 
    ! inversion. I comment it but think about it
    !do p = 1, 6
    !  !do k = 1, nangt ! do not include ghost cells
    !  do k = 0, nangt+1 ! here we also save ghost cells 
    !      !do j = 1, nangt
    !      do j = 0, nangt+1
    !          do i =1, nrt
    !              n = n_matrix(i, p, j, k)
    !!              !temp(i, j, k, p) = source(n)
    !              if (abs(source(n))>= 1.d50) then             
    !                write(*,*) i, j, k, p, source(n), temp(i, j, k, p)
    !              endif
    !!
    !!              ! we save analitic temperature
    !!              ! Perez-Azorin 2006 Isotropic
    !!              
    !!              ic = 2*i - 1
    !!              jc = 2*j
    !!              kc = 2*k
    !!
    !!              !ang_term = dsin(theta(jc, kc, p))**2 + (dcos(theta(jc, kc, p))**2)/(1.d0 + omegatau*omegatau) !standard PA
    !!              cos_th_prime = dsin(beta_ang)*dsin(theta(jc, kc ,p))*(-dcos(gamma_ang)*dcos(phi(jc, kc, p)) + &
    !!              &              dsin(gamma_ang)*dsin(phi(jc, kc, p))) + dcos(beta_ang)*dcos(theta(jc, kc, p))
    !!              ang_term = 1.d0 - cos_th_prime*cos_th_prime*omegatau*omegatau/(1.d0 + omegatau*omegatau) !generalized PA
    !!
    !!              temp_anlt(i, j, k, p) = T_int*((1./(time))**(3./2.))*exp(-r(ic)*r(ic)*ang_term/(4.*kappa_perp*(time)))
    !          enddo
    !      enddo
    !  enddo
    !enddo

    ! at the first time step, write the L2 before matrix is inverted. This is L2 at time 0
  !  write(130, *) time, sum((temp(:,1:nangt,1:nangt,:)-temp_anlt(:,1:nangt,1:nangt,:))**2)/sum(temp_anlt(:,1:nangt,1:nangt,:)**2) 
  !else
  !  open(unit = 130, file=trim(name), access = 'append')
  !endif
  ! ****** PEREZ_AZORIN TEST - L2 *********



  !************************ MATRIX INVERSION *********************************

  ! we do not have to pass matrix in the inversion routine because it will be modified 
  ! while we want to preserve the sparse elements. We define an appropriate array
  ! matrix_solve for that


  matrix_solve = matrix

 call dgbsv(6*nrt*(nangt+2)*(nangt+2), bandwd, bandwd, 1, matrix_solve, 3*bandwd+1, &
 &          IPIV, source, 6*nrt*(nangt+2)*(nangt+2), INFO)



  !*************************************************************************

 ! We assign to the temperature vector the values resulting from the linear system solving

 ! We assign to the temperature array the value stored in source
! This loop seems pretty slow...
!$OMP Parallel do private(p,i,j,k, n)
do i =1, nrt
 do p = 1, 6
    !do k = 1, nangt ! do not include ghost cells
    do k = 0, nangt+1 ! here we also save ghost cells 
        !do j = 1, nangt
        do j = 0, nangt+1
            
                n = n_matrix(i, p, j, k)
                temp(i, j, k, p) = source(n)

                if(temp(i,j,k,p)<1.d-2) temp(i,j,k,p)=1.d-2 ! we impose a floor to the temperature

                !if (k == 0 .or. j == 0 .or. j ==  nangt+1 .or. k == nangt+1) then
                !  write(*,*) i, j, k, p, temp(i,j,k,p)
                !endif 

                if (temp(i,j,k,p)>1d2 .or. temp(i,j,k,p)<0.d0) then 
                  write(*,*) i, j, k, p, 'temp:', temp(i,j,k,p)
                endif 


                ! we save analitic temperature
                ! Perez-Azorin 2006 
                !
                !ic = 2*i - 1
                !jc = 2*j
                !kc = 2*k

                !ang_term = dsin(theta(jc, kc, p))**2 + (dcos(theta(jc, kc, p))**2)/(1.d0 + omegatau*omegatau) !standard PA
                !cos_th_prime = dsin(beta_ang)*dsin(theta(jc, kc ,p))*(-dcos(gamma_ang)*dcos(phi(jc, kc, p)) + &
                !&              dsin(gamma_ang)*dsin(phi(jc, kc, p))) + dcos(beta_ang)*dcos(theta(jc, kc, p))
                !ang_term = 1.d0 - cos_th_prime*cos_th_prime*omegatau*omegatau/(1.d0 + omegatau*omegatau) !generalized PA

                !temp_anlt(i, j, k, p) = T_int*((1./(time+dt))**(3./2.))*exp(-r(ic)*r(ic)*ang_term/(4.*kappa_perp*(time+dt)))
            enddo
        enddo
    enddo
  enddo
  !$OMP end Parallel do

! *********** DEGUGGING ghost cells and flux ***************

 !print*, 'Cells'
 !do k = 1,2*nangt+1
 ! print*, k, eta(k), pmt(k)
 !enddo 

 !print*, 'PESI'
 !do k =1, nangt
 ! kc = 2*k
 ! if (k<(nangt+1)/2) then
 !   print*, k, 0.5*edge_wt(kc), 1.0 - 0.5*(1.0 + dabs(pmt(kc-1) - eta(kc+1))/(eta(2)-eta(0))), dabs(pmt(kc-1) - eta(kc+1))
 ! elseif (k == (nangt+1)/2) then
 !   print*, k, 0.5*edge_wt(kc), 0.d0, dabs(pmt(kc-1) + eta(kc+1))
 ! else
 !   print*, k, 0.5*edge_wt(kc), 1.0 - 0.5*(1.0 + dabs(pmt(kc+1) - eta(kc-1))/(eta(2)-eta(0))), dabs(pmt(kc+1) - eta(kc-1))
 ! endif
 !enddo 
 !stop

  !do k=1, nangt
  !  print*, eta(2*k), edge_wt(2*k), pmt(2*k)
  !enddo 
  !do p = 1, 6
  !  do k = 1, nangt
  !    do j = 1, nangt
  !      kc = 2*k
  !      i = nrt-1
  !      if (p == 1 .and. k==8) then
  !        !write(*,*) p, j, flux_r_out(i,p,j,k)-flux_r_out(i-1,p,j,k) + flux_xi_xip(i,p,j,k) - &
  !        !& flux_xi_xip(i,p,j-1,k) + flux_eta_etap(i,p,j,k) - flux_eta_etap(i,p,j,k-1)
  !        !write(*,*) p, j, (temp(i+1,j,k+1,p) + temp(i+1,j,k,p)  &
  !        !    &  - temp(i-1,j,k+1,p) - temp(i-1,j,k,p)), temp(i,j,k+1,p) - temp(i,j,k,p), &
  !        !    &  (temp(i,j+1,k+1,p) + temp(i,j+1,k,p)  &
  !        !    &  - temp(i,j-1,k+1,p) - temp(i,j-1,k,p))
  !        !write(*,*) p, j, flux_xi_xip(i,j,k,p) - flux_xi_xip(i,j-1,k,p)  
  !        !write(*,*) p, j, flux_eta_etap(i,j,k,p) - flux_eta_etap(i,j,k-1,p)
  !        !write(*,*) p, j, flux_r_out(i-1,j,k,p)-flux_r_out(i-1,j,k,p), zreta(i,p,j,k)
  !        !write(*,*) p, j, temp(i,j,k+1,p), temp(i,j,k,p),temp(i,j,k+1,p) - temp(i,j,k,p)
  !        write(*,*) p, j, (temp(i+1,j,k+1,p) + temp(i+1,j,k,p)  &
  !            &  - temp(i-1,j,k+1,p) - temp(i-1,j,k,p))
  !      elseif (p == 4 .and. k == 8) then
  !         !write(*,*) p, j, flux_r_out(i,p,j,k)-flux_r_out(i-1,p,j,k) + flux_xi_xip(i,p,j,k) - &
  !        !& flux_xi_xip(i,p,j-1,k) + flux_eta_etap(i,p,j,k) - flux_eta_etap(i,p,j,k-1)
  !        !write(*,*) p, j, flux_xi_xip(i,j,k,p) - flux_xi_xip(i,j-1,k,p) 
  !        !write(*,*) p, j, flux_eta_etap(i,j,k,p) - flux_eta_etap(i,j,k-1,p)
  !        !write(*,*) p, j, flux_r_out(i-1,j,k,p)-flux_r_out(i-1,j,k,p), zreta(i,p,j,k)
  !        !write(*,*) p, j, (temp(i+1,j,k+1,p) + temp(i+1,j,k,p)  &
  !        !    &  - temp(i-1,j,k+1,p) - temp(i-1,j,k,p)), temp(i,j,k+1,p) - temp(i,j,k,p), &
  !        !    &  (temp(i,j+1,k+1,p) + temp(i,j+1,k,p)  &
  !        !    &  - temp(i,j-1,k+1,p) - temp(i,j-1,k,p))
  !        !write(*,*) p, j, temp(i,j,k+1,p), temp(i,j,k,p), temp(i,j,k+1,p) - temp(i,j,k,p)
  !        write(*,*) p, j, (temp(i+1,j,k+1,p) + temp(i+1,j,k,p)  &
  !            &  - temp(i-1,j,k+1,p) - temp(i-1,j,k,p))
  !      endif
  !    enddo
  !  enddo
  !enddo
  ! *************************************************

   !Now we write the fluxes for a diagnostic 
!  i_print = 5
!  
!  if (CalculateFlux) then
!
!    do p = 1, 6
!      do k = 1, nangt
!        do j = 1, nangt
!          ic = 2*i_print -1 
!          jc = 2*j
!          kc = 2*k
!          ! here we calculate the net_fluxs along different directions
!        
!          netflux_r(p, j, k) = -zrr(i_print, p, j, k)*(temp(i_print+1, j, k, p) - temp(i_print, j, k, p)) &
!          & - zrxi(i_print, p, j, k)*(temp(i_print+1, j+1, k, p) + temp(i_print, j+1, k, p)               &
!          & - temp(i_print+1, j-1, k, p) - temp(i_print, j-1, k, p))                                      &
!          & - zreta(i_print, p, j, k)*(temp(i_print+1, j, k+1, p) + temp(i_print, j, k+1, p)              &
!          & - temp(i_print+1, j, k-1, p) - temp(i_print, j, k-1, p))                                      &
!          & + zrr(i_print-1, p, j, k)*(temp(i_print, j, k, p) - temp(i_print-1, j, k, p))                 &
!          & + zrxi(i_print-1, p, j, k)*(temp(i_print, j+1, k, p) + temp(i_print-1, j, k, p)               &
!          & - temp(i_print, j-1, k, p) - temp(i_print-1, j-1, k, p))
!
!
!          netflux_xi(p, j, k) = -zxir(i_print, p, j, k)*(temp(i_print+1, j+1, k, p) + temp(i_print+1, j, k, p)  &
!          & - temp(i_print-1, j+1, k, p) - temp(i_print-1, j, k, p))                                            &
!          & - zxixi(i_print, p, j, k)*(temp(i_print, j+1, k, p) - temp(i_print, j, k, p))                       &
!          & - zxieta(i_print, p, j, k)*(temp(i_print, j+1, k+1, p) + temp(i_print, j, k+1, p)                   &
!          & - temp(i_print, j+1, k-1, p) - temp(i_print, j, k-1, p))                                            &
!          & + zxir(i_print, p, j-1, k)*(temp(i_print+1, j, k, p) + temp(i_print+1, j-1, k, p)                   &
!          & - temp(i_print-1, j, k, p) - temp(i_print-1, j-1, k, p))                                            &
!          & + zxixi(i_print, p, j-1, k)*(temp(i_print, j, k, p) - temp(i_print, j-1, k, p))                     &
!          & + zxieta(i_print, p, j-1, k)*(temp(i_print, j, k+1, p) + temp(i_print, j-1, k+1, p)                 &
!          & - temp(i_print, j, k-1, p) - temp(i_print, j-1, k-1, p)) 
!
!          netflux_eta(p, j, k) = - zetar(i_print, p, j, k)*(temp(i_print+1, j, k+1, p) + temp(i_print+1, j, k, p)  &
!          & - temp(i_print-1, j, k+1, p) - temp(i_print-1, j, k, p))                                               &
!          & - zetaxi(i_print, p, j, k)*(temp(i_print, j+1, k+1, p) + temp(i_print, j+1, k, p)                      &
!          & - temp(i_print, j-1, k+1, p) - temp(i_print, j-1, k, p))                                               &
!          & - zetaeta(i_print, p, j, k)*(temp(i_print, j, k+1, p) - temp(i_print, j, k, p))                        &
!          & + zetar(i_print, p, j, k-1)*(temp(i_print+1, j, k, p) + temp(i_print+1, j, k-1, p)                     &
!          & - temp(i_print-1, j, k, p) - temp(i_print-1, j, k-1, p))                                               &
!          & + zetaxi(i_print, p, j, k-1)*(temp(i_print, j+1, k, p) + temp(i_print, j+1, k-1, p)                    &
!          & - temp(i_print, j-1, k, p) - temp(i_print, j-1, k-1, p))                                               &
!          & + zetaeta(i_print, p, j, k-1)*(temp(i_print, j, k, p) - temp(i_print, j, k-1, p)) 
!        
!        enddo
!      enddo
!    enddo
!  endif  
  !write(*,*) "L2_Temp =", sum((temp(:,:,:,:)-temp_anlt(:,:,:,:))**2)/sum(temp_anlt(:,:,:,:)**2)

  

  ! PEREZ AZORIN - Write the L2 measure on file 
  !
  !write(130, *) time+dt, sum((temp(:,1:nangt,1:nangt,:)-temp_anlt(:,1:nangt,1:nangt,:))**2)/sum(temp_anlt(:,1:nangt,1:nangt,:)**2)
  !close(130)
  ! PEREZ AZORIN

  !do i = 0, 3*bandwd
  !  write(100,*) matrix(i,:)
  !enddo
 

!  write(*,*) "INFO = 0 --> succesfull"
!  write(*,*) "INFO = -i, the i-th argument had an illegal value"
!  write(*,*) "INFO = i, U(i,i) is exactly zero. U singular"
!  write(*,*) "INFO =", INFO



  if (INFO .ne. 0) then
    write(*,*) "WARNING:: Matrix inversion not succesful"
    write(*,*) "Simulation Step =", counter
    write(*,*) "Time =", time
    write(*,*) "INFO=", INFO
  endif
 
  !source vector

!  write(*,*) "Source After"

  !do i = 0,6*nrt*(nangt+2)*(nangt+2)-1
  !  write(103, *) source(i)
!    write(*,*) i, source(i)
  !enddo

  do p = 1, 6

!    close(100)
  !close(101)
  !close(102)
  !close(103)
!:
!: ! Print the temperature profile in real coordinate
!:    write (name,'(a13,i1,a2)') "out/temp_prof", p, ".d"
!:    write (name_kappa,'(a14,i1,a2)') "out/kappa_prof", p, ".d"
!:    write (name_flux,'(a13,i1,a2)') "out/flux_prof", p, ".d"
!:    if (FirstCall) then
!:        open(unit = 115, file=trim(name), status = 'replace')
!:        if (CalculateFlux) then
!:          open(unit = 215, file=trim(name_kappa), status = 'replace')
!:          open(unit = 315, file=trim(name_flux), status = 'replace')
!:        endif
!:    else
!:        open(unit = 115, file=trim(name), access = 'append')
!:        if (CalculateFlux) then
!:          open(unit = 215, file=trim(name_kappa), access = 'append')
!:          open(unit = 315, file=trim(name_flux), access = 'append')
!:        endif
!:
!:
!:        ! We leave two blank line
!:        write(115, *)
!:        write(115, *)
!:        if (CalculateFlux) then
!:          write(215, *)
!:          write(215, *)
!:          write(315, *)
!:          write(315, *)
!:        endif
!:
!:    endif
    !open(unit=115,file=trim(name),status='replace')

    do j = 1, nangt 
        do k =1, nangt
            jc = 2*j
            kc = 2*k

            ! PEREZ AZORIN ANALYTIC PROFILE 
            !cos_th_prime = dsin(beta_ang)*dsin(theta(jc, kc ,p))*(-dcos(gamma_ang)*dcos(phi(jc, kc, p)) + &
            !&              dsin(gamma_ang)*dsin(phi(jc, kc, p))) + dcos(beta_ang)*dcos(theta(jc, kc, p))
            !ang_term = 1.d0 - cos_th_prime*cos_th_prime*omegatau*omegatau/(1.d0 + omegatau*omegatau) !generalized PA

            !T_analit = T_int*((1./(time+dt))**(3./2.))*exp(-r(i_print)*r(i_print)*ang_term/(4.*kappa_perp*(time+dt)))
            !T_analit = 0.d0

            !1) theta 2) phi 3) rel_error_t 4) T 5) L2_T 6) F_r 7) F_xi 8) F_eta 9) F_r_anlt 10) F_xi_anlt 11) F_eta_anlt
            !12) rel_error_Fr 13) rel_error_Fxi 14) rel_error_Feta 6) br 7) bxi 8) beta 9) T_anlt 10) T- T_anlt
            !write(115, *) theta(jc, kc, p), phi(jc, kc, p),         &
            !&   (temp(i_print, j, k, p) - temp_anlt(i_print, j, k, p))/temp_anlt(i_print, j, k, p), &
            !&   temp(i_print, j, k, p), &
            !&   ((temp(i_print, j, k, p) - temp_anlt(i_print, j, k, p))**2.0)/(temp_anlt(i_print, j, k, p)**2.0), &
            !!  &(source(n_matrix(i_print,p, j, k)) - temp_anlt(i_print, j, k, p))/temp_anlt(i_print, j, k, p), &
            !!&  source(n_matrix(i_print,p, j, k)), &
            !!&  (source(n_matrix(i_print,p, j, k)) - T_analit)**2.0/(T_analit**2.0), &
            !!& flux_r_out(i_print, p, j, k), flux_xi_xip(i_print, p, j, k), flux_eta_etap(i_print, p, j, k),   &
            !!& flux_r_anl_out(i_print, p, j, k), flux_xi_anl_xip(i_print, p, j, k), flux_eta_anl_etap(i_print, p, j, k),  &
            !!& flux_r_out(i_print, p, j, k)/flux_r_anl_out(i_print, p, j, k) - 1.d0,            &
            !!& flux_xi_xip(i_print, p, j, k)/flux_xi_anl_xip(i_print, p, j, k) - 1.d0,           &
            !!& flux_eta_etap(i_print, p, j, k)/flux_eta_anl_etap(i_print, p, j, k) - 1.d0, &
            !& br(i_print, jc, kc, p), bxi(i_print, jc, kc, p), beta(i_print, jc, kc, p), &
            !& temp_anlt(i_print, j, k, p), temp(i_print, j, k, p) - temp_anlt(i_print, j, k, p)
            
            if (CalculateFlux) then
              if (k .ne. 0 .and. j .ne. 0 .and. k .ne. nangt+1 .and. j .ne. nangt+1) then

                ! 1) theta 2)phi 3) kappa_rr 4) kappa_rxi 5) kappa_reta 6) kappa_xir 7) kappa_xixi 8) kappa_xieta
                ! 9) kappa_etar 10) kappa_etaxi 11) kappa_etaeta

                write(215, *) theta(jc, kc, p), phi(jc, kc, p), kappa_rr_out(i_print, p, j, k), &
                & kappa_rxi_out(i_print, p, j, k), kappa_reta_out(i_print, p, j, k), kappa_xir_xip(i_print, p, j, k), &
                & kappa_xixi_xip(i_print, p, j, k), kappa_xieta_xip(i_print, p, j, k), kappa_xieta_xip(i_print, p, j, k), &
                & kappa_etar_etap(i_print, p, j, k), kappa_etaxi_etap(i_print, p, j, k), kappa_etaeta_etap(i_print, p, j, k)
            
                ! 1) theta 2) phi 3) netflux_r 4) netflux_xi 5) netflux_eta 6) netflux_angular 7) netflux_total
!                write(315, *) theta(jc, kc, p), phi(jc, kc, p), netflux_r(p, j, k), netflux_xi(p, j, k), netflux_eta(p, j, k), &
!               &  netflux_xi(p, j, k) +  netflux_eta(p, j, k),  netflux_r(p, j, k) +  netflux_xi(p, j, k) +  netflux_eta(p, j, k)
              endif
            endif
        enddo
    enddo
    close(115)
    if (CalculateFlux) then
      close(215)
      close(315)
    endif
 enddo

  ! tem(1:nrt+1,1:6,0:nangt+1,0:nangt+1) = x

  counter = counter + 1 ! variable to count the timesteps in termal evolution
  time = time + dt !For Debugging (or not?)
  if (FirstCall) then
    FirstCall = .FALSE.
  endif
  contains

  integer function n_matrix(i, p, j, k)
    implicit none
    integer, intent(in) :: i, p, j, k
 
    n_matrix = j + (nangt+2)*(k + (nangt+2)*((p - 1) + 6*(i-1)))  ! change i ---> i-1
  end function n_matrix

  !----- Functions used to build the conductivity tensors. We leave the kappa_perp term outside the definition --
  !----- In general I keep here only the geometrical terms and put all the other terms (magnetic field, omega*tau) outside

  real*8 function A_kappa(i, kappa)
    implicit none
    integer, intent(in) :: i
    real*8, intent(in) :: kappa

    A_kappa = kappa/elambda(i) !check if one or 2 elambda

  end function A_kappa

  real*8 function B_kappa(i, k, kappa)
    implicit none
    integer, intent(in) :: i,k
    real*8, intent(in) :: kappa

    B_kappa = kappa*D(k)/r(i)
  end function B_kappa

  real*8 function G_kappa(i, j, kappa)
    implicit none
    integer, intent(in) :: i, j
    real*8, intent(in) :: kappa

    G_kappa = kappa*C(j)/r(i)
  end function G_kappa

  real*8 function E_kappa(i, j, k, kappa)
    implicit none
    integer, intent(in) :: i, j, k
    real*8, intent(in) :: kappa
    
    E_kappa = kappa*X(j)*Y(k)/r(i)/D(k)
  end function E_kappa

  real*8 function F_kappa(i, j, k, kappa)
    implicit none
    integer, intent(in) :: i, j, k
    real*8, intent(in) :: kappa
    
    F_kappa = kappa*X(j)*Y(k)/r(i)/C(j)
  end function F_kappa

  real*8 function H_kappa(i, p, j, k, omegatau, kappa)
    implicit none
    integer, intent(in) :: i, j, k, p
    real*8, intent(in) :: omegatau
    real*8, intent(in) :: kappa

    H_kappa = kappa*omegatau*omegatau*br(i, j, k, p)/bm(i, j, k, p)/elambda(i)
    if (H_kappa /= H_kappa) then
    print*,bm(i,j,k,p),elambda(i),i,j,k,p,H_kappa,omegatau,kappa,"HK"
   
    stop
    endif
    !H_kappa = 0.d0
  end function H_kappa

  real*8 function I_kappa(i, p, j, k, omegatau, kappa)
    implicit none
    integer, intent(in) :: i, j, k, p
    real*8, intent(in) :: omegatau
    real*8, intent(in) :: kappa

    I_kappa = (D(k) - X(j)*X(j)*Y(k)*Y(k)/C(j)/C(j)/D(k))/r(i)
    I_kappa = kappa*I_kappa*omegatau*omegatau*bxi(i, j, k, p)/bm(i, j, k, p)
    !I_kappa = 0.d0
  end function I_kappa

  real*8 function J_kappa(i, p, j, k, omegatau, kappa)
    implicit none
    integer, intent(in) :: i, j, k, p
    real*8, intent(in) :: omegatau
    real*8, intent(in) :: kappa

    J_kappa = (C(j) - X(j)*X(j)*Y(k)*Y(k)/D(k)/D(k)/C(j))/r(i)
    J_kappa = kappa*J_kappa*omegatau*omegatau*beta(i, j, k, p)/bm(i, j, k, p)
    !J_kappa = 0.d0
  end function J_kappa

  real*8 function K_kappa(i, p, j, k, omegatau, kappa)
    implicit none
    integer, intent(in) :: i, j, k, p
    real*8, intent(in) :: omegatau
    real*8, intent(in) :: kappa
    real*8 :: sqrtdelta
  
    sqrtdelta = dsqrt(1 + X(j)**2 + Y(k)**2)

    K_kappa = (bxi(i, j, k, p)*X(j)*Y(k)/C(j) - D(k)*beta(i, j, k, p))/bm(i, j, k, p)
    K_kappa = kappa*K_kappa*omegatau*sqrtdelta/C(j)/D(k)/r(i)
    !K_kappa = 0.d0
  end function K_kappa

  real*8 function L_kappa(i, p, j, k, omegatau, kappa)
    implicit none
    integer, intent(in) :: i, p, j, k
    real*8, intent(in) :: omegatau
    real*8, intent(in) :: kappa
    real*8 :: sqrtdelta

    sqrtdelta = dsqrt(1 + X(j)**2 + Y(k)**2)

    L_kappa = (bxi(i, j, k, p)*C(j) - beta(i, j, k, p)*X(j)*Y(k)/D(k))/bm(i, j, k, p)
    L_kappa = kappa*L_kappa*omegatau*sqrtdelta/C(j)/D(k)/r(i)
    !L_kappa = 0.d0
  end function L_kappa

  real*8 function M_kappa(i, p, j, k, omegatau, kappa)
    implicit none
    integer, intent(in) :: i, p, j, k
    real*8, intent(in) :: omegatau
    real*8, intent(in) :: kappa
    real*8 :: sqrtdelta

    sqrtdelta = dsqrt(1 + X(j)**2 + Y(k)**2)

    M_kappa = (C(j)*D(k)*beta(i, j, k, p) - X(j)*Y(k)*bxi(i, j, k, p))/bm(i, j, k, p)
    M_kappa = kappa*M_kappa*omegatau/sqrtdelta/elambda(i)
    !M_kappa = 0.d0
  end function M_kappa

  real*8 function N_kappa(i, p, j, k, omegatau, kappa)
    implicit none
    integer, intent(in) :: i, p, j, k
    real*8, intent(in) :: omegatau
    real*8, intent(in) :: kappa
    real*8 :: sqrtdelta

    sqrtdelta = dsqrt(1 + X(j)**2 + Y(k)**2)

    N_kappa = kappa*omegatau*sqrtdelta*br(i, j, k, p)/bm(i, j, k, p)/D(k)/r(i)
    !N_kappa = 0.d0
  end function N_kappa

  real*8 function O_kappa(i, p, j, k, omegatau, kappa)
    implicit none
    integer, intent(in) :: i, p, j, k
    real*8, intent(in) :: omegatau
    real*8, intent(in) :: kappa
    real*8 :: sqrtdelta

    sqrtdelta = dsqrt(1 + X(j)**2 + Y(k)**2)

    O_kappa = (X(j)*Y(k)*beta(i, j, k, p) - C(j)*D(k)*bxi(i, j, k, p))/bm(i, j, k, p)
    O_kappa = kappa*O_kappa*omegatau/sqrtdelta/elambda(i)
    !O_kappa = 0.d0
 end function O_kappa

  real*8 function P_kappa(i, p, j, k, omegatau, kappa)
    implicit none
    integer, intent(in) :: i, p, j, k
    real*8, intent(in) :: omegatau
    real*8, intent(in) :: kappa
    real*8 :: sqrtdelta

    sqrtdelta = dsqrt(1 + X(j)**2 + Y(k)**2)

    P_kappa = kappa*omegatau*sqrtdelta*br(i, j, k, p)/bm(i, j, k, p)/C(j)/r(i)
    !P_kappa = 0.d0
  end function P_kappa
!---------------------------------------------------------------

end subroutine tevol



subroutine core_cooling(T_core, dt, flux)
  
  use input_params, only: profile
  use grid, only: cv_core_tot, qnu_core_tot
  use grid, only: cv_core_tot_der, qnu_core_tot_der
  use grid, only: enu
  use microphysics, only: CoreCooling_Implicit

  implicit none

  real*8, intent(in) :: dt
  real*8, intent(in) :: flux
  real*8, intent(inout) :: T_core
  real*8, dimension(1:4) :: integrated_cv, integrated_loss
  ! DV: Why do they have 4 values? Even for RK4, you don't need them, you only need k1,k2,k3,k4...
  real*8 :: T_0 = 20.d0, alpha = 8.d0
  real*8 :: k1, k2, k3, k4 !Runge-Kutta Parameters
  real*8 :: deriv_rhs
  real*8 :: ddt, t_prime
  real*8 :: core_cooling_rhs, core_cooling_rhs_deriv

  integer :: i, N

  !Perform the integral here-----
  !------------------------------
  

  N = 10000
  t_prime = 0.0
  ddt = dt/N
!  write(*,*) 'T_core before cooling:', T_core

  ! Explicit Euler
  !do i = 1, N
  !  integrated_cv(1) = 1.d0
  !  integrated_loss(1) = - (T_core/T_0)**alpha
  !  T_core = T_core + ddt*integrated_loss(1)/integrated_cv(1)
  !enddo
  !4th order RK
  
  !do i = 1, N
 
  !  t_prime = t_prime + ddt

 !   integrated_cv(1) = 1.d0
  !  integrated_loss(1) = - (T_core/T_0)**alpha
  !  k1 = integrated_loss(1)/integrated_cv(1)

  !  integrated_cv(2) = 1.d0
  !  integrated_loss(2) = - ((T_core + 0.5*k1*ddt)/T_0)**alpha
  !  k2 = integrated_loss(2)/integrated_cv(2)

  !  integrated_cv(3) = 1.d0
  !  integrated_loss(3) = - ((T_core + 0.5*k2*ddt)/T_0)**alpha
  !  k3 = integrated_loss(3)/integrated_cv(3)

  !  integrated_cv(4) = 1.d0
  !  integrated_loss(4) = - ((T_core + k3*ddt)/T_0)**alpha
  !  k4 = integrated_loss(4)/integrated_cv(4)

  !  T_core = T_core + (k1 + 2*k2 + 2*k3 + k4)*ddt/6.d0

  !enddo

  ! Implicit Euler
 
  !deriv_rhs = alpha*integrated_loss(1)/T_core  ! this is valid only for this cv =1

  !T_core = (T_core + dt*(integrated_loss(1)/integrated_cv(1) - deriv_rhs*T_core))/(1.d0 - deriv_rhs*dt)

  ! now we call a routine to calculate f(T) and f'(T)
  if (profile .eq. 'realist') then ! realistic mycrophysics 
    core_cooling_rhs = - qnu_core_tot/cv_core_tot
    core_cooling_rhs_deriv = - (qnu_core_tot_der*cv_core_tot - &
    &                         qnu_core_tot*cv_core_tot_der)/cv_core_tot**2 
  else 
    !semianalytic
    call CoreCooling_Implicit(core_cooling_rhs, core_cooling_rhs_deriv)
  endif

  ! now we calculate T in the next timestep 
  !T_core = (T_core + dt*(core_cooling_rhs - core_cooling_rhs_deriv*enu(1)*T_core - flux/cv_core_tot))/ &
  !&        (1.d0 - core_cooling_rhs_deriv*enu(1)*dt)

  !New core cooling (here cv_core_tot in the flux term is calculayed implicitly)

  !correct for flux term
  core_cooling_rhs = core_cooling_rhs - flux/cv_core_tot
  core_cooling_rhs_deriv = core_cooling_rhs_deriv +flux*cv_core_tot_der/cv_core_tot**2

  T_core = (T_core + dt*(core_cooling_rhs - core_cooling_rhs_deriv*enu(1)*T_core))/ &
  &        (1.d0 - core_cooling_rhs_deriv*enu(1)*dt)


end subroutine core_cooling

end module thermal_evolution

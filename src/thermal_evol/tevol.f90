!!-------------------------------------------------------------------------------
!! Module: Thermal evolution
!!
!!> @author
!!> Jose Pons Botella
!!> Daniele Viganò
!!> Alberto García-García
!!
!!> 
!!> @brief This module is responsible for the temperature evolution.
!!
!!-------------------------------------------------------------------------------
module thermal_evol

  use constants, only: PI, T_YEAR
  use grid, only: kmax, lmax, tem, rmc, c_v
  use grid, only: rb, theta, belam, benu
  use grid, only: q_joule_average, q_neutrino, cfluxb, sfluxb
  use microphysics, only: dfc2, dfc3, anis2, anis3
  use math, only: solvtb

  contains

      subroutine tevol(dty)
!!
!!   This routine generates the block-tridiagonal system,
!!   where a, b, c are submatrices with non-zero elements 
!!   only in the diagonal and its lower and upper elements
!!
!!    | b1 c1 0  0  0.....| |x1|    |r1|
!!    | a2 b2 c2 0  0.....| |x2|    |r2|
!!    | 0  a3 b3 c3 0.....| |. |  _ |. |
!!    | ..................| |. |  _ |. |
!!    | ..................| |. |    |. |
!!    | .............an bn| |xn|    |rn|
!!
!!     And solves the system to obtain x (here, updated temperature)
!!
!!     Note that the temperature is solved on a half-resolution grid,
!!      the cell centers correspond to rc(l) = rb(2*l-1) (odd)
!!                                 and zc(k) = theta(2*i-2)  (even)
!!
!   ********************VARIABLE DEFINITIONS***************************
!
!     ARK  - K-FACE CELL AREAS
!     ARL  - L-FACE CELL AREAS
!     CVOL - CELL VOLUMES
!     DFC1 - DIFFUSION COEFFICIENT FOR LOWER LEFT REGION OF CELL K,L
!     DFC2 - DIFFUSION COEFFICIENT FOR UPPER LEFT REGION OF CELL K,L
!     DFC3 - DIFFUSION COEFFICIENT FOR UPPER RIGHT REGION OF CELL K,L
!     DFC4 - DIFFUSION COEFFICIENT FOR LOWER RIGHT REGION OF CELL K,L
!     c_v  - TIME-DERIVATIVE COEFFICIENTS
!     RMC  - REMOVAL COEFFICIENTS
!
!     DT   - TIME STEP
!
!     KMAX - MAXIMUM K-INDEX OF REAL CELLS AND EDGES
!     LMAX - MAXIMUM L-INDEX OF REAL CELLS AND EDGES
!     KD   - K-DIMENSIONING PARAMETER
!     LD   - L-DIMENSIONING PARAMETER
!
!   **********************END DEFINITIONS******************************
!
      implicit none

      integer k, l, ik, il

!   input variables

      real*8 dty, dt, ark, arl, dsource, source
      real*8, dimension(kmax,lmax) :: crr,czz,cxz,cxr
      real*8 dan1rr,dan2rr,dan3rr,dan4rr,dan1rz,dan2rz,dan3rz,dan4rz
      real*8 dan1zz,dan2zz,dan3zz,dan4zz,d1,d2,d3,d4,xaef
      real*8 dan1zr,dan2zr,dan3zr,dan4zr
!      real*8 Theat

      real*8, dimension(lmax-1,kmax-1,3) :: a, b, c
      real*8, dimension(lmax-1,kmax-1) :: x, rr
      real*8 z(kmax+1), r(lmax+1)

        dt = dty*T_YEAR    ! change timestep from years to seconds

!     initialize a,b,c,x,rr to zero
        x=0.d0
        rr=0.d0
        a=0.d0
        b=0.d0
        c=0.d0

!--------------------------------------------------------
!     LOCAL ARRAYS OF RADIUS (R) AND ANGLES (Z) AT THE THERMAL GRID INTERFACES
!     WE SOLVE THE TEMPERATURE IN A SMALLER GRID 
!        (HALF THE POINTS IN BOTH DIRECTIONS)
!     FOR CLARITY, WE KEEP THE DUPLICATED GRID ONLY HERE
!     TODO: we should use directly what is in grid, using rb an theta below
!--------------------------------------------------------
       do l=1,lmax+1
         r(l) = rb(2*l)
       enddo
       do k=1,kmax   
         z(k)=theta(2*k-1)
       enddo
       z(kmax+1) = 2*z(kmax)-z(kmax-1)
!--------------------------------------------------------
!     CALCULATE MATRIX ELEMENTS FOR CELL-CENTERED EQUATIONS
!--------------------------------------------------------
      do k=2,kmax
        do l=2,lmax

          d1 = dfc3(k-1,l)*benu(2*l-1)
          d2 = dfc2(k,l)*benu(2*l)
          d3 = dfc3(k,l)*benu(2*l-1)
          d4 = dfc2(k,l-1)*benu(2*l-2)

          dan1rr=d1*anis3(1,1,k-1,l)/belam(2*l-1)
          dan2rr=d2*anis2(1,1,k,l)/belam(2*l)
          dan3rr=d3*anis3(1,1,k,l)/belam(2*l-1)
          dan4rr=d4*anis2(1,1,k,l-1)/belam(2*l-2)

          dan1rz=d1*anis3(1,2,k-1,l)
          dan2rz=d2*anis2(1,2,k,l)
          dan3rz=d3*anis3(1,2,k,l)
          dan4rz=d4*anis2(1,2,k,l-1)

          dan1zr=d1*anis3(2,1,k-1,l)/belam(2*l-1)
          dan2zr=d2*anis2(2,1,k,l)/belam(2*l)
          dan3zr=d3*anis3(2,1,k,l)/belam(2*l-1)
          dan4zr=d4*anis2(2,1,k,l-1)/belam(2*l-2)

          dan1zz=d1*anis3(2,2,k-1,l)
          dan2zz=d2*anis2(2,2,k,l)
          dan3zz=d3*anis3(2,2,k,l)
          dan4zz=d4*anis2(2,2,k,l-1)

!---------------------------------------------------
!     CALCULATE NORMAL AREA ELEMENTS
!     TODO: in principle these areas are already calculated.
!          In any case, we don't need to do it every time...
!---------------------------------------------------

          ark = 1.5d0*(r(l)+r(l-1))/(r(l)**2+r(l)*r(l-1)+r(l-1)**2) &
     &         *sin(z(k))/(cos(z(k-1))-cos(z(k)))

          arl = 3.d0*r(l)**2/(r(l)**3-r(l-1)**3)/belam(2*l-1)
         
          crr(k,l)=dan2rr*arl/(0.5d0*(r(l+1)-r(l-1))) 
          czz(k,l)=dan3zz*ark/(0.25d0*(r(l)+r(l-1))*(z(k+1)-z(k-1)))
          cxz(k,l)=dan3zr*ark/4.d0/(r(l)-r(l-1))
          cxr(k,l)=dan2rz*arl/4.d0/r(l)/(z(k)-z(k-1))
       enddo
      enddo
!----------------------------------------------------------------
!      Angular Boundary Conditions: 
!         - the theta area elements at the axis are zero --> cxz=czz=0
!         - symmetry imposed across the axis for the r_coefs 
!----------------------------------------------------------------
      cxz(1,:) = 0.d0
      czz(1,:) = 0.d0
      
!----------------------------------------------------------------
!      Radial Boundary Conditions: 
!         - just copy in the innermost cell
!         - set to zero the radial flux coefficients in the outermost point
!           because it will be included later using the envelope model.
!----------------------------------------------------------------
      crr(:,1)=crr(:,2)
      cxr(:,1)=cxr(:,2)
      crr(:,lmax)=0.d0
      cxr(:,lmax)=0.d0

!----------------------------------------------------------------------
!     Calculate matrix elements
!     There are only mkax-1, and lmax-1 centered values of temperature
!     Redimension and shift everything by one
!----------------------------------------------------------------------
      do ik=1,kmax-1
        do il=1,lmax-1
         k = ik+1
         l = il+1
        
         ! Q_joule already includes the benu^2
         xaef = dt/c_v(k,l) 
         ! In dsource, there is a benu instead of benu^2, because the derivative of emissivity 
         ! is calculated as a function of physical T, while the linearization is in terms of
         ! the redshifted one.
         dsource = xaef*rmc(k,l)*benu(2*l-1) 
         source = xaef*(q_neutrino(k,l)*benu(2*l-1)**2 + q_joule_average(k,l))

         a(il,ik,1) = (-cxz(k-1,l)-cxr(k,l-1))*xaef
         b(il,ik,1) = (-czz(k-1,l)-cxr(k,l-1)+cxr(k,l))*xaef
         c(il,ik,1) = (cxz(k-1,l)+cxr(k,l))*xaef
         a(il,ik,2) = (-cxz(k-1,l)+cxz(k,l)-crr(k,l-1))*xaef
         c(il,ik,2) = (cxz(k-1,l)-cxz(k,l)-crr(k,l))*xaef 
         a(il,ik,3) = (cxz(k,l)+cxr(k,l-1))*xaef
         b(il,ik,3) =  (-czz(k,l)+cxr(k,l-1)-cxr(k,l))*xaef
         c(il,ik,3) = (-cxz(k,l)-cxr(k,l))*xaef

         b(il,ik,2) = 1.d0-(a(il,ik,1)+b(il,ik,1)+c(il,ik,1)+a(il,ik,2)  &
     &                   + c(il,ik,2)+a(il,ik,3)+b(il,ik,3)+c(il,ik,3))  &
     &                   + dsource
         
! Old definition if depsilon/dT is included in the source
!         rr(il,ik) =tem(k,l) - source 

! Definition of vector with depsilon/dT not included in source
         rr(il,ik) = tem(k,l)*(1d0 + dsource) - source 

!------------------------------------------------------------------------------
       enddo
      enddo 
      
!----------------------------------------------------------------------
!             Angular boundary conditions
!----------------------------------------------------------------------
!----------------------
!    Theta = 0
!----------------------
         k = 1

         a(:,k,2) = a(:,k,2) + a(:,k,1)
         b(:,k,2) = b(:,k,2) + b(:,k,1)
         c(:,k,2) = c(:,k,2) + c(:,k,1)
         a(:,k,1) = 0.d0
         b(:,k,1) = 0.d0
         c(:,k,1) = 0.d0
!----------------------
!    Theta = pi
!----------------------
         k = kmax-1

         a(:,k,2) = a(:,k,2) + a(:,k,3)
         b(:,k,2) = b(:,k,2) + b(:,k,3)
         c(:,k,2) = c(:,k,2) + c(:,k,3)
         a(:,k,3) = 0.d0
         b(:,k,3) = 0.d0
         c(:,k,3) = 0.d0
!----------------------------------------------------------------------
!             Inner core boundary conditions
!----------------------------------------------------------------------
      l=1
!
!    Fixed T core boundary conditions
!   Comment the next two lines for cooling calculations
!
!        Tcore = tem(k,1) 
!        rr(l,k) =rr(l,k) - Tcore*(a(l,k,1)+a(l,k,2)+a(l,k,3))
!
!      Comment the next three lines for fixed Tcore B.C.
!
         b(l,:,1) = b(l,:,1) + a(l,:,1)
         b(l,:,2) = b(l,:,2) + a(l,:,2)
         b(l,:,3) = b(l,:,3) + a(l,:,3)

         a(l,:,1:3) = 0.d0
!----------------------------------------------------------------------
!             Surface boundary conditions
!----------------------------------------------------------------------
      il=lmax-1
      l = lmax
! ****************************************************
! HEATING BC!!!!!!!
!      do k=1,kmax/10
!      	Theat = 1.e1
! 	rr(il,ik) = rr(il,ik) - Theat*(c(il,ik,1)+c(il,ik,2)+c(il,ik,3))
!       c(il,ik,1:3) = 0.d0
!      enddo
! ***************************************************** 
      do ik=1,kmax-1
        k = ik+1

        arl = 3.d0*r(l)**2/(r(l)**3-r(l-1)**3)/belam(2*l-1)
        xaef = dt/c_v(k,l) 

        b(il,ik,2) = b(il,ik,2) + xaef*arl*cfluxb(k)
        rr(il,ik) = rr(il,ik) - xaef*arl*sfluxb(k)
      enddo

!--------------------------------------------------
!     DIFFUSION MATRIX ELEMENTS CALCULATION COMPLETE
!         SOLVE THE BLOCK TRIDIAGONAL MATRIX
!--------------------------------------------------
      call solvtb(lmax-1,kmax-1,a,b,c,rr,x)

! Update temperatures.
      do k=2,kmax
      do l=2,lmax
         tem(k,l)=x(l-1,k-1)
       enddo
       enddo
       tem(1,:)=tem(2,:)

!  Set a minimum temperature in case some transient problem appears locally
       do k=1,kmax
         do l=1,lmax
           tem(k,l)=dmax1(tem(k,l),1d-2)
         enddo
       enddo

    end subroutine tevol

end module thermal_evol


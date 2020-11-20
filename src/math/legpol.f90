!> @brief TODO: Pending but this will probably be refactored.
!!
!! TODO.
module legpol

  ! Modules --------------------------------------------------------------------
  use grid, only: nang, nleg, ng

  ! Module constants -----------------------------------------------------------
  ! None

  ! Module variables -----------------------------------------------------------
  real*8, dimension(:, :), allocatable, save :: pln     ! Legendre polynomials function of theta
  real*8, dimension(:, :), allocatable, save :: dpln    ! theta-derivatives of Legendre polynomials
  real*8, dimension(:, :), allocatable, save :: plnx    ! Legendre polynomials at zeroes positions
  real*8, dimension(:), allocatable, save :: xg         ! Locations of zeros for Gauss-Legendre quadrature
  real*8, dimension(:), allocatable, save :: wg         ! Weights of each zero for Gauss-Legendre quadrature
  real*8, dimension(:), allocatable, save :: frel       ! Relativistic factors  used in potentials' derivatives for potential solutions
  real*8, dimension(:), allocatable, save :: blout      ! Weights representing the multipoles

  contains

    subroutine allocate_legendre_polynomials()

      ! Modules ----------------------------------------------------------------
      ! None.

      implicit none

      ! Input arguments --------------------------------------------------------
      ! None.

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
      ! None.

      ! ------------------------------------------------------------------------

      allocate(pln(nang, 0:nleg))
      allocate(dpln(nang, 0:nleg))
      allocate(plnx(ng, 0:nleg))
      allocate(xg(ng))
      allocate(wg(ng))
      allocate(frel(1:nleg))
      allocate(blout(0:nleg))

    end subroutine allocate_legendre_polynomials

    ! --------------------------------------------------------------------------
    !> Compute Legendre polynomials.
    !> @brief It calculates the Legendre polynomials using recurrence relations
    ! --------------------------------------------------------------------------
    subroutine compute_legendre_polynomials()

      ! Modules ----------------------------------------------------------------
      use grid, only: cth, sth

      implicit none

      ! Input arguments --------------------------------------------------------
      ! None.

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
      integer :: i, n

      ! ------------------------------------------------------------------------

      ! n=0 and n=1.
      pln(:, 0) = 1d0
      dpln(:, 0) = 0d0
      pln(:, 1) = cth(1:nang)
      dpln(:, 1) = -sth(1:nang)

      ! The rest of n.  
      do n = 1, nleg-1
        pln(1:nang, n+1) = (dble(2 * n + 1) * cth(1:nang) *pln(:,n) &
          &                 - dble(n) * pln(1:nang, n-1)) / dble(n+1)
      end do ! n

      do i = 1, nang
        do n = 2, nleg
          if(sth(i) == 0d0) then
            dpln(i, n) = 0d0
          else
            dpln(i, n) = (n * cth(i) * pln(i, n) - n * pln(i, n-1)) / sth(i)
          end if
        end do ! n
      end do ! i

    end subroutine compute_legendre_polynomials

    ! --------------------------------------------------------------------------
    !> Relativistic correction for the vacuum boundary conditions
    !> @brief   It calculates the Barnes hypergeometrical functions, used for the magnetic BC
    
    ! Note: anu_rel_correction: see Eq. (56) of Radler et al. 2001, PRD 64
    ! Note: frel: see f_n factors in Pons, Miralles and Geppert 2009, Eq. (24).

    !> @param[in] eps double_compactness GM/c^2R
    !
    ! --------------------------------------------------------------------------
    subroutine get_rel_correction(double_compactness)

      ! Modules ----------------------------------------------------------------
      ! None.

      implicit none

      ! Input arguments --------------------------------------------------------
      real*8, intent(in) :: double_compactness  

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
      integer :: i, n
      ! i is the nu of Radler et al. 2001, sec. VI, n is the l
      real*8 :: bb, sum1, sum2, frel_sl 
      ! bb is the b_nu of eq. (58)
      ! frel_sl is eq. (57)
      ! sum1 and sum2 are the two sums appearing in eq. (57)
      ! sum2 is eq. (56) where the sum_nu(a_nu compactness^nu) is equivalent to sum2,
      ! which is the denominator of (58)
      ! the conversion from frel_sl (used for S_l in Radler 2001) 
      ! to frel (used for Phi_l in Pons 2009)
      ! considers that S_l = (1/r)Phi_l

      ! Initialize (used for nonrelativistic case)
      frel = 1d0

      if (double_compactness /= 0.) then

        do n = 1, nleg

          bb = 1d0
          sum1 = bb
          sum2 = bb
          
          ! First 50 terms (more than enough) of the infinite series in nu=0,infinity
          do i = 1, 50
  
            bb = dble((n + i)**2 - 1) * double_compactness * bb / dble((2 * n + i + 1) * i)
            sum1 = sum1 + dble(n + i + 1) * bb / dble(n + 1)
            sum2 = sum2 + bb
  
          end do ! i
  
          frel_sl = sum1 / sum2
          frel(n) = frel_sl + (frel_sl - 1d0) / dble(n)
        
        end do ! n

      endif

    end subroutine get_rel_correction

    ! --------------------------------------------------------------------------
    !> @brief compute the Legendre polynomials at the location of the zeros.
    ! --------------------------------------------------------------------------
    subroutine getbl_init()

      ! Modules ----------------------------------------------------------------
      use math, only: compute_gausslegendre_weights

      implicit none

      ! Input arguments --------------------------------------------------------
      ! None.

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
      integer :: i, l

      ! ------------------------------------------------------------------------
    
      call compute_gausslegendre_weights(ng, -1d0, 1d0, xg, wg)
  
      do i = 1, ng

        plnx(i, 0) = 1d0
        plnx(i, 1) = xg(i)

        do l = 1, nleg-1
          plnx(i, l+1) = (dble(2 * l + 1) * xg(i) * plnx(i, l) &
            &             - dble(l) * plnx(i, l - 1)) / dble(l + 1)
        end do ! l

      end do ! i

    end subroutine getbl_init

    ! --------------------------------------------------------------------------
    !> @brief Computes the Multipole decomposition (bl coefs.) given the angles 
    !>        (th) and the function (bzout) on a grid (dimension m) by numerically
    !>        evaluating an integral (Gauss-Legendre quadrature)
    !> @param[in] m     dimension of the function
    !> @param[in] cth   cos(theta)
    !> @param[in] fun   function to be weighted
    ! Note: sum = int(mu [-1,1]) P_l(mu)*Brout(rns,mu) = - 2b_l (l+1)/(2l+1).
    ! Note: The normalization follows b_l = - l*a_l (cf. Vigano thesis, eqs. 2.93 and 2.94).
    ! --------------------------------------------------------------------------
    subroutine getbl(m, cth, fun)

      ! Modules ----------------------------------------------------------------
      ! None.

      implicit none

      ! Input arguments --------------------------------------------------------
      integer, intent(in) :: m
      real*8, intent(in) :: cth(m)
      real*8, intent(in) :: fun(m)

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables  -------------------------------------------------------
      integer :: i, l
      ! Sum is the integral, funx is the function interpolated at zeroes location
      real*8 :: sum
      real*8 :: funx

      ! ------------------------------------------------------------------------

      blout = 0d0
      do l = 0, nleg

        sum = 0d0

        do i = 1, ng

          call fun_interp(m, cth, fun, xg(i), funx)
          sum = sum + wg(i) * plnx(i,l) * funx

        end do ! i
      
        if(dabs(sum) < 1d-10) then
            sum=0d0
        end if

        blout(l) = -dble(2 * l + 1) / (2d0 * dble(l + 1)) * sum

      end do ! l

    end subroutine getbl

    ! --------------------------------------------------------------------------
    !> Interpolation function with a linear interpolation of the cosine.
    !> @brief Given a 1-D array theta and a function f with dimension m
    !>        returns the value of the function funx=f(csx)
    !> @param[in] m Dimension of the array.
    !> @param[in] cth Angular grid (cth = cos (theta)).
    !> @param[in] fun Function to be interpolated.
    !> @param[in] csx Input value (cos(x)).
    !> @param[out] funx Output value of the function at (cos(x)).
    !----------------------------------------------------------------------- ---
    subroutine fun_interp(m, cth, fun, csx, funx)

      ! Modules ----------------------------------------------------------------
      ! None.

      implicit none

      ! Input arguments --------------------------------------------------------
      integer, intent(in) :: m
      real*8, intent(in) :: cth(m)
      real*8, intent(in) :: fun(m)
      real*8, intent(in) :: csx
      real*8, intent(out) :: funx

      ! Local constants --------------------------------------------------------
      ! None.

      ! Local variables --------------------------------------------------------
      integer i1, i2, i
      logical :: located

      ! ------------------------------------------------------------------------

      ! Locate the position in the grid .
      ! TODO: Is there an f90 system function to locate ?
      i1 = 1
      i2 = 2
      located = .false.

      do i = 2, m-1

        if(csx >= cth(i)) then
          i1 = i - 1
          i2 = i
          located = .true.
          exit
        end if

      end do ! i

      if (located .eqv. .false.) then
        i1 = m - 1
        i2 = m
      end if

      ! interpolate (in cos(x))
      funx = fun(i1) + (fun(i2) - fun(i1)) / (cth(i2) - cth(i1)) * (csx - cth(i1))

    end subroutine fun_interp

end module legpol
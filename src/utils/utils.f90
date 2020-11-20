!-------------------------------------------------------------------------------
! Magneto Thermal 2D
!-------------------------------------------------------------------------------
! Module: Utils
!
!> @author
!> Jose Pons Botella
!> Daniele Viganò
!> Alberto Garcia-Garcia
!
!> Handful of utility functions and subroutines.
!> @brief This module contains a bunch of useful routines which serve various
!>        purposes and are too general to be included in any other module.
!-------------------------------------------------------------------------------
module utils

  contains

    !---------------------------------------------------------------------------
    !> Converts an integer into a string.
    !> @brief
    !> @return The string representation in (i0) format.
    !---------------------------------------------------------------------------
    function int_to_string(input) result(str)

      ! Modules ----------------------------------------------------------------
      ! None.

      implicit none

      ! Input parameters -------------------------------------------------------
      integer, intent(in) :: input

      ! Local constants --------------------------------------------------------
      character(len = 4), parameter :: INT_STR_FORMAT = "(i0)"

      ! Local variables --------------------------------------------------------
      character(:), allocatable :: str
      character(range(input)+2) :: tmp

      ! ------------------------------------------------------------------------

      write(tmp, INT_STR_FORMAT) input
      str = trim(tmp)

    end function int_to_string

    !---------------------------------------------------------------------------
    !> Converts a real into a string.
    !> @brief
    !> @return The string representation in (f5.2) format.
    !---------------------------------------------------------------------------
    function real_to_string(input) result(str)

      ! Modules ----------------------------------------------------------------
      ! None.

      implicit none

      ! Input parameters -------------------------------------------------------
      character(:), allocatable :: str

      ! Local constants --------------------------------------------------------
      character(len = 6), parameter :: REAL_STR_FORMAT = "(f5.2)"

      ! Local variables --------------------------------------------------------
      real*8, intent(in) :: input
      character(range(input)+2) :: tmp

      ! ------------------------------------------------------------------------

      write(tmp, REAL_STR_FORMAT) input
      str = trim(tmp)

    end function real_to_string

    !---------------------------------------------------------------------------
    !> Gets an available unit to write.
    !> @brief
    !> @return The integer id for a free unit to associate to a file.
    !---------------------------------------------------------------------------
    function get_free_unit() result(free_unit)

      ! Modules ----------------------------------------------------------------
      ! None.

      implicit none

      ! Input parameters -------------------------------------------------------
      ! None.

      ! Local constants --------------------------------------------------------
      integer, parameter :: UNIT_MIN = 10, UNIT_MAX = 1000

      ! Local variables --------------------------------------------------------
      integer :: free_unit
      logical :: opened
      integer :: unit

      ! ------------------------------------------------------------------------

      free_unit = -1

      do unit = UNIT_MIN, UNIT_MAX
        inquire(unit=unit, opened=opened)
        if (.not. opened) then
          free_unit = unit
          return 
        end if
      end do ! unit

    end function get_free_unit

end module utils
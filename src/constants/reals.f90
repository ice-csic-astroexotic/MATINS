!> @brief Floating point precision types.
!!
!! In Fortran90 and later, intrinsic types such as real have a kind attribute
!! which guarantees a specific precision or range. REAL*8 and counterparts
!! should no longer be used [2]. double precision is also no longer needed and
!! can be thought as a real kind. The intrinsic function select_real_kind(p, r)
!! should be used as a portable way of getting guaranteed types of p significant
!! digits of precision and an exponent range of at least r.
!!
!! [1] Chin, Worth, and Greenough. Thoughts on Using the Features of Fortran 95.
!!
module reals

  implicit none

  integer, parameter :: float = selected_real_kind(6, 37)
  integer, parameter :: double = selected_real_kind(15, 307)
  integer, parameter :: quad = selected_real_kind(33, 4931)

end module reals
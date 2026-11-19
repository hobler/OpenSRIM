module Vector_mod
!
!  Defines a (3-d) vector data type with vector operations
!
!  Used Modules:
!  ============
!
   use Const_mod, only: INT4, REAL8
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
   save
!
!  Public Procedures and Operators:
!  ===============================
!
   public :: operator(.dot.)                 ! r=v1.v2
   public :: operator(.cross.)               ! r=v1 x v2
   public :: abs                             ! r=|v|
   public :: abs2                            ! r=|v|**2

   interface operator(.dot.)
      module procedure dotProduct
   end interface

   interface operator(.cross.)
      module procedure crossProduct
   end interface

   interface abs
      module procedure absVector
   end interface

   interface abs2
      module procedure abs2Vector
   end interface
!
!=======================================================================
!
contains
!
!=======================================================================
!
   function dotProduct( v1, v2 ) result( r )
!           ==========  i   i            o
!
!     Take dot product of two vectors
!
!  Arguments:
!
      real(REAL8), intent(in) :: v1(3), v2(3)
      real(REAL8) :: r
!
!  Start of Executable Code:
!  ------------------------
!
      r = sum( v1*v2 )

   end function dotProduct
!
!=======================================================================
!
   function crossProduct( v1, v2 ) result( v3 )
!           ============  i   i            o
!
!     Take cross product of two vectors
!
!  Arguments:
!
      real(REAL8), intent(in) :: v1(3), v2(3)
      real(REAL8) :: v3(3)
!
!  Start of Executable Code:
!  ------------------------
!
      v3(1) = v1(2) * v2(3) - v1(3) * v2(2)
      v3(2) = v1(3) * v2(1) - v1(1) * v2(3)
      v3(3) = v1(1) * v2(2) - v1(2) * v2(1)

   end function crossProduct
!
!=======================================================================
!
   function absVector( v ) result( r )
!           =========  i           o
!
!     Calculate the length of a vector.
!
!  Arguments:
!
      real(REAL8), intent(in) :: v(3)
      real(REAL8) :: r
!
!  Start of Executable Code:
!  ------------------------
!
      r = sqrt( sum( v*v ) )

   end function absVector
!
!=======================================================================
!
   function abs2Vector( v ) result( r )
!           ==========  i           o
!
!     Calculate the square of the length of a vector.
!
!  Arguments:
!
      real(REAL8), intent(in) :: v(3)
      real(REAL8) :: r
!
!  Start of Executable Code:
!  ------------------------
!
      r = sum( v*v )

   end function abs2Vector
!
!=======================================================================
!
end module Vector_mod

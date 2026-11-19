module Const_mod
!
!  Contains parameters to be used in all parts of the program including kind numbers
!
!  Used Modules:
!  ============
!
!  None.
!
!  Default Type and Access:
!  =======================
!
   implicit none
   public
!
!  Version of Code:
!
   character(10), parameter :: VERSION = '2024-03+'
!
!  Kind Numbers:
!  ============
!
!  Finally, the kind numbers should be defined like this:
!
   integer, parameter :: INT1 = selected_int_kind( r=2 )    ! usually 1-byte integer (< 128)
   integer, parameter :: INT2 = selected_int_kind( r=4 )    ! usually 2-byte integer (< 32768)
!  integer, parameter :: INT4 = selected_int_kind( r=9 )    ! usually 4-byte integer (< 2.15e9)
   integer, parameter :: INT8 = selected_int_kind( r=12 )   ! usually 8-byte integer (< 9.22e18)

!  integer, parameter :: REAL4 = selected_real_kind (p=6,r=36)   ! usually < 1.23456e37
!  integer, parameter :: REAL8 = selected_real_kind (p=13,r=100) ! usually < 1.23456789012345e308
!
!  For consistency with old data types [integer(INT4) == integer,
!  real(REAL4) == real, real(REAL8) == double precision], however,
!  we need to define, until all type declarations have been updated:

   integer, parameter :: INT4 = kind(1)
   integer, parameter :: REAL4 = kind(1.0)
   integer, parameter :: REAL8 = kind(1.0d0)
!
!  Numerical Constants:
!  ===================

   real(REAL8), parameter :: EPSMB = 1.e-07_REAL8
   real(REAL8), parameter :: EPSMA = EPSMB**2                  ! 1.e-14_REAL8
   real(REAL8), parameter :: SMALL = 1.e-20_REAL8
   real(REAL8), parameter :: BIG = 1.e+30_REAL8
   real(REAL8), parameter :: NEGBIG = - BIG
   real(REAL8), parameter :: UNDEFINED = -8.97518e+11_REAL8     ! some negative number

   real(REAL8), parameter :: ZERO  = 0
   real(REAL8), parameter :: ONE   = 1
   real(REAL8), parameter :: TWO   = 2
   real(REAL8), parameter :: THREE = 3
   real(REAL8), parameter :: FOUR  = 4
   real(REAL8), parameter :: FIVE  = 5
   real(REAL8), parameter :: SIX   = 6
   real(REAL8), parameter :: EIGHT = 8
   real(REAL8), parameter :: NINE = 9
   real(REAL8), parameter :: TEN   = 10

   real(REAL8), parameter :: HALF  = ONE/TWO
   real(REAL8), parameter :: THIRD = ONE/THREE
   real(REAL8), parameter :: QUART = ONE/FOUR
   real(REAL8), parameter :: SIXTH = ONE/SIX
   real(REAL8), parameter :: EIGHTH = ONE/EIGHT
   real(REAL8), parameter :: TENTH = ONE/TEN
!
!  Mathematical Constants:
!  ======================

   real(REAL8), parameter :: PI = 3.1415926535897932384_REAL8
   real(REAL8), parameter :: EULER = 2.71828182846_REAL8
   real(REAL8), parameter :: SQRT2 = 1.41421356237_REAL8
   real(REAL8), parameter :: SQRT3 = 1.732050807568877_REAL8
!
!  Enumerations:
!  ============

   type State_type                              ! projectile state
!       ==========
      integer(INT4) :: SUBTHRESHOLD = 0
      integer(INT4) :: STARTING = 1
      integer(INT4) :: MOVING = 2
      integer(INT4) :: STOPPED = 3
      integer(INT4) :: LEAVING = 4
      integer(INT4) :: REENTERED = 5
      integer(INT4) :: ERROR = -1
   end type State_type

   type(State_type), parameter :: STATE = State_type()

   type Grid_type                               ! cell grid
!       =========
      integer(INT4) :: regular = 1
      integer(INT4) :: rectilinear = 2
      integer(INT4) :: curvilinear = 3
   end type Grid_type

   type RelaxMethod_type                        ! numerical relaxation method
!       ================
      integer(INT4) :: euler = 1
      integer(INT4) :: lagrange = 2
   end type RelaxMethod_type

   type RelaxModel_type                         ! physical relaxation model
!       ===============
      integer(INT4) :: none = 1
      integer(INT4) :: full = 2
      integer(INT4) :: viscoelastic = 3
   end type RelaxModel_type

   type BndType_type                       ! boundary condition as defined by lshyperoctree
!       ============
      integer(INT4) :: unknown = 0
      integer(INT4) :: periodic = 1
      integer(INT4) :: reflect = 2
      integer(INT4) :: extend = 3
      integer(INT4) :: free = 4
   end type BndType_type

   type Enum_type
!       =========
      type(Grid_type) :: grid = Grid_type()
      type(RelaxMethod_type) :: relaxmethod = RelaxMethod_type()
      type(RelaxModel_type) :: relaxmodel = RelaxModel_type()
      type(BndType_type) :: bndtype = BndType_type()
   end type Enum_type

   type(Enum_type), parameter :: ENUM = Enum_type()
!
!=======================================================================
!
end module Const_mod

module Array_mod
!
!  Defines array operations
!
!  Used Modules:
!  ============
!
   use Const_mod, only: INT4, REAL8,  &
                        ZERO
   use Types_mod, only: Proj_type, Ipos_type
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
   save
!
!  Public Procedures
!  =================
!
   public :: alloc                           ! allocate array and check
!            =====
   public :: changeShape                     ! change the shape of an array, copying existing
!            ===========                     !   elements if possible and filling empty
!                                            !   positions with ZEROs
   public :: addReflect                      ! Extend array by the reflection of itself
!            ==========
   public :: check_allocation                ! check if allocation was successful
!            ================

   interface alloc
      module procedure alloc1dLog
      module procedure alloc1dLogLb
      module procedure alloc2dLog
      module procedure alloc1dChar
      module procedure alloc1dCharLb
      module procedure alloc1dInt4
      module procedure alloc1dInt4Lb
      module procedure alloc2dInt4
      module procedure alloc2dInt4Lb
      module procedure alloc1dReal8
      module procedure alloc1dReal8Lb
      module procedure alloc1dReal8Src
      module procedure alloc2dReal8
      module procedure alloc2dReal8Lb
      module procedure alloc3dReal8
      module procedure alloc3dReal8Lb
      module procedure alloc4dReal8
      module procedure alloc1dIpos
      module procedure alloc2dIpos
      module procedure alloc1dProj
   end interface

   interface changeShape
      module procedure changeLogShape1d
      module procedure changeLogShape2d
      module procedure changeInt4Shape1d
      module procedure changeInt4Shape2d
      module procedure changeReal8Shape1d
      module procedure changeReal8Shape2d
      module procedure changeReal8Shape3d
      module procedure changeReal8Shape4d
      module procedure changeIposShape1d
      module procedure changeIposShape2d
      module procedure changeProjShape1d
   end interface

   interface addReflect
      module procedure addReflect1dLog
      module procedure addReflect2dLog
      module procedure addReflect1dReal8
      module procedure addReflect2dReal8
      module procedure addReflect3dReal8
      module procedure addReflect4dReal8
   end interface
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine alloc1dLog( array, n, array_name, proc_name )
!             ==========   io    i      i           i
!
!     Allocate a 1D logical array with size n and check allocation
!
!  Arguments
!
      logical, allocatable, intent(inout) :: array(:) ! array to be allocated
      integer(INT4) :: n                              ! array size
      character(*), intent(in) :: array_name          ! array name in the calling procedure
      character(*), intent(in) :: proc_name           ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc1dLog
!
!=======================================================================
!
   subroutine alloc1dLogLb( array, n1, n2, array_name, proc_name )
!             ============   io    i   i       i           i
!
!     Allocate a 1D logical array with lower bound n1 and upper bound n2 and check allocation
!
!  Arguments
!
      logical, allocatable, intent(inout) :: array(:) ! array to be allocated
      integer(INT4) :: n1                             ! lower bound of array
      integer(INT4) :: n2                             ! upper bound of array
      character(*), intent(in) :: array_name          ! array name in the calling procedure
      character(*), intent(in) :: proc_name           ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n1:n2 ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc1dLogLb
!
!=======================================================================
!
   subroutine alloc2dLog( array, n, array_name, proc_name )
!             ==========   io    i      i           i
!
!     Allocate a 2D logical array with shape n and check allocation
!
!  Arguments
!
      logical, allocatable, intent(inout) :: array(:,:)  ! array to be allocated
      integer(INT4), intent(in) :: n(2)                  ! array shape
      character(*), intent(in) :: array_name             ! array name in the calling
                                                         !    procedure
      character(*), intent(in) :: proc_name              ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n(1), n(2) ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc2dLog
!
!=======================================================================
!
   subroutine alloc1dChar( array, n, array_name, proc_name )
!             ===========   io    i      i           i
!
!     Allocate a 1D character array with size n and check allocation
!
!  Arguments
!
      character(*), allocatable, intent(inout) :: array(:)  ! array to be allocated
      integer(INT4) :: n                                    ! array size
      character(*), intent(in) :: array_name                ! array name in calling procedure
      character(*), intent(in) :: proc_name                 ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc1dChar
!
!=======================================================================
!
   subroutine alloc1dCharLb( array, n1, n2, array_name, proc_name )
!             =============   io    i   i       i           i
!
!     Allocate a 1D character array with lower bound n1 and upper bound n2 and check allocation
!
!  Arguments
!
      character(*), allocatable, intent(inout) :: array(:)  ! array to be allocated
      integer(INT4) :: n1                                   ! lower bound of array
      integer(INT4) :: n2                                   ! upper bound of array
      character(*), intent(in) :: array_name                ! array name in calling procedure
      character(*), intent(in) :: proc_name                 ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n1:n2 ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc1dCharLb
!
!=======================================================================
!
   subroutine alloc1dInt4( array, n, array_name, proc_name )
!             ===========   io    i      i           i
!
!     Allocate a 1D Int4 array with size n and check allocation
!
!  Arguments
!
      integer(INT4), allocatable, intent(inout) :: array(:) ! array to be allocated
      integer(INT4) :: n                                    ! array size
      character(*), intent(in) :: array_name                ! array name in the calling
                                                            !    procedure
      character(*), intent(in) :: proc_name                 ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc1dInt4
!
!=======================================================================
!
   subroutine alloc1dInt4Lb( array, n1, n2, array_name, proc_name )
!             =============   io    i   i       i           i
!
!     Allocate a 1D Int4 array with lower bound n1 and upper bound n2 and check allocation
!
!  Arguments
!
      integer(INT4), allocatable, intent(inout) :: array(:) ! array to be allocated
      integer(INT4) :: n1                                   ! lower bound of array
      integer(INT4) :: n2                                   ! upper bound of array
      character(*), intent(in) :: array_name                ! array name in the calling
                                                            !    procedure
      character(*), intent(in) :: proc_name                 ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n1:n2 ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc1dInt4Lb
!
!=======================================================================
!
   subroutine alloc2dInt4( array, n, array_name, proc_name )
!             ===========   io    i      i           i
!
!     Allocate a 2D Int4 array with shape n and check allocation
!
!  Arguments
!
      integer(INT4), allocatable, intent(inout) :: array(:,:)  ! array to be allocated
      integer(INT4), intent(in) :: n(2)                        ! array shape
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n(1), n(2) ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc2dInt4
!
!=======================================================================
!
   subroutine alloc2dInt4Lb( array, n1, n2, array_name, proc_name )
!             ==============   io    i   i       i           i
!
!     Allocate a 2D INT4 array with lower bounds n1(:) and upper bounds n2(:),
!     and check allocation
!
!  Arguments
!
      integer(INT4), allocatable, intent(inout) :: array(:,:)  ! array to be allocated
      integer(INT4) :: n1(2)                                   ! lower bounds
      integer(INT4) :: n2(2)                                   ! upper bounds
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n1(1):n2(1), n1(2):n2(2) ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc2dInt4Lb
!
!=======================================================================
!
   subroutine alloc1dReal8( array, n, array_name, proc_name )
!             ============   io    i      i           i
!
!     Allocate a 1D Real8 array with size n and check allocation
!
!  Arguments
!
      real(REAL8), allocatable, intent(inout) :: array(:)   ! array to be allocated
      integer(INT4) :: n                                    ! array size
      character(*), intent(in) :: array_name                ! array name in the calling
                                                            !    procedure
      character(*), intent(in) :: proc_name                 ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc1dReal8
!
!=======================================================================
!
   subroutine alloc1dReal8Lb( array, n1, n2, array_name, proc_name )
!             ==============   io    i   i       i           i
!
!     Allocate a 1D Real8 array with lower bound n1 and upper bound n2 and check allocation
!
!  Arguments
!
      real(REAL8), allocatable, intent(inout) :: array(:)   ! array to be allocated
      integer(INT4) :: n1                                   ! lower bound of array
      integer(INT4) :: n2                                   ! upper bound of array
      character(*), intent(in) :: array_name                ! array name in the calling
                                                            !    procedure
      character(*), intent(in) :: proc_name                 ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n1:n2 ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc1dReal8Lb
!
!=======================================================================
!
   subroutine alloc1dReal8Src( array, array_src, array_name, proc_name )
!             ===============   io        i          i           i
!
!     Allocate a 1D Real8 array from a source and check allocation
!
!     NOTE: array's lower bounds will always be 1
!
!  Arguments
!
      real(REAL8), allocatable, intent(inout) :: array(:)   ! array to be allocated
      real(REAL8), intent(in) :: array_src(:)               ! source array
      character(*), intent(in) :: array_name                ! array name in the calling
                                                            !    procedure
      character(*), intent(in) :: proc_name                 ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array, source=array_src, stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc1dReal8Src
!
!=======================================================================
!
   subroutine alloc2dReal8( array, n, array_name, proc_name )
!             ============   io    i      i           i
!
!     Allocate a 2D Real8 array with shape n and check allocation
!
!  Arguments
!
      real(REAL8), allocatable, intent(inout) :: array(:,:) ! array to be allocated
      integer(INT4) :: n(2)                                 ! array shape
      character(*), intent(in) :: array_name                ! array name in the calling
                                                            !    procedure
      character(*), intent(in) :: proc_name                 ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n(1), n(2) ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc2dReal8
!
!=======================================================================
!
   subroutine alloc2dReal8Lb( array, n1, n2, array_name, proc_name )
!             ==============   io    i   i       i           i
!
!     Allocate a 2D Real8 array with lower bounds n1(:) and upper bounds n2(:),
!     and check allocation
!
!  Arguments
!
      real(REAL8), allocatable, intent(inout) :: array(:,:) ! array to be allocated
      integer(INT4) :: n1(2)                                ! lower bounds
      integer(INT4) :: n2(2)                                ! upper bounds
      character(*), intent(in) :: array_name                ! array name in the calling
                                                            !    procedure
      character(*), intent(in) :: proc_name                 ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n1(1):n2(1), n1(2):n2(2) ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc2dReal8Lb
!
!=======================================================================
!
   subroutine alloc3dReal8( array, n, array_name, proc_name )
!             ============   io    i      i           i
!
!     Allocate a 3D Real8 array with shape n and check allocation
!
!  Arguments
!
      real(REAL8), allocatable, intent(inout) :: array(:,:,:)  ! array to be allocated
      integer(INT4) :: n(3)                                    ! array shape
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n(1), n(2), n(3) ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc3dReal8
!
!=======================================================================
!
   subroutine alloc3dReal8Lb( array, n1, n2, array_name, proc_name )
!             ==============   io    i   i       i           i
!
!     Allocate a 3D Real8 array with lower bounds n1(:) and upper bounds n2(:),
!     and check allocation
!
!  Arguments
!
      real(REAL8), allocatable, intent(inout) :: array(:,:,:)  ! array to be allocated
      integer(INT4) :: n1(3)                                   ! lower bounds
      integer(INT4) :: n2(3)                                   ! upper bounds
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n1(1):n2(1), n1(2):n2(2), n1(3):n2(3) ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc3dReal8Lb
!
!=======================================================================
!
   subroutine alloc4dReal8( array, n, array_name, proc_name )
!             ============   io    i      i           i
!
!     Allocate a 4D Real8 array with shape n and check allocation
!
!  Arguments
!
      real(REAL8), allocatable, intent(inout) :: array(:,:,:,:)! array to be allocated
      integer(INT4), intent(in) :: n(4)                        ! array shape
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n(1), n(2), n(3), n(4) ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc4dReal8
!
!=======================================================================
!
   subroutine alloc1dIpos( array, n, array_name, proc_name )
!             ===========   io    i      i           i
!
!     Allocate a 1D Ipos array with size n and check allocation
!
!  Arguments
!
      type(Ipos_type), allocatable, intent(inout) :: array(:)  ! array to be allocated
      integer(INT4) :: n                                       ! array size
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc1dIpos
!
!=======================================================================
!
   subroutine alloc2dIpos( array, n, array_name, proc_name )
!             ===========   io    i      i           i
!
!     Allocate a 2D Int4 array with shape n and check allocation
!
!  Arguments
!
      type(Ipos_type), allocatable, intent(inout) :: array(:,:)! array to be allocated
      integer(INT4) :: n(2)                                    ! array shape
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n(1), n(2) ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc2dIpos
!
!=======================================================================
!
   subroutine alloc1dProj( array, n, array_name, proc_name )
!             ===========   io    i      i           i
!
!     Allocate a 1D Proj_type array with size n and check allocation
!
!  Arguments
!
      type(Proj_type), allocatable, intent(inout) :: array(:)  ! array to be allocated
      integer(INT4) :: n                                       ! array size
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
!
!  Local variables
!
      integer(INT4) :: istat                 ! status returned by the allocate statement
      character(256) :: errmsg               ! error message returned by the allocate statement
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( array( n ), stat=istat, errmsg=errmsg )

      call check_allocation( proc_name, array_name, istat, errmsg )
!          ================      i          i         i      i

   end subroutine alloc1dProj
!
!=======================================================================
!
   subroutine changeLogShape1d( array, nnew, array_name, proc_name, shift_, move_right_ )
!             ================   io     i        i           i        i          i
!
!     Change the shape of a 1D array, copying existing elements if possible and filling
!     empty positions with .false.
!
!  Arguments:
!
      logical, allocatable, intent(inout) :: array(:) ! array
      integer(INT4), intent(in) :: nnew               ! new dimension
      character(*), intent(in) :: array_name          ! array name in the calling procedure
      character(*), intent(in) :: proc_name           ! name of calling procedure
      integer(INT4), optional, intent(in) :: shift_   ! shift of data
      logical, optional, intent(in) :: move_right_    ! if true, right-align old
!                                                     !    elements in new array
!  Local Variables:
!
      integer(INT4) :: shift, iold, jold, nold, inew, jnew
      logical, allocatable :: tmp(:)
!
!  Start of Executable Code:
!  ------------------------
!
      nold = size( array )

      if( present( move_right_ ) ) then
         if( present( shift_ ) ) then
            error stop 'changeLogShape1d: move_right_ and shift must not both be present'
         else
            if( move_right_ ) then
               shift = nnew - nold
            else
               shift = 0
            end if
         end if
      else
         if( present( shift_ ) ) then
            shift = shift_
         else
            shift = 0
         end if
      end if

      iold = max( 1, 1-shift )
      inew = max( 1+shift, 1)
      jold = min( nold, nnew-shift )
      jnew = min( nold+shift, nnew )

      call alloc( tmp, nnew, proc_name, array_name )
!          =====  io    i        i          i
      tmp(:) = .false.
      tmp(inew:jnew) = array(iold:jold)

      call move_alloc( tmp, array )

   end subroutine changeLogShape1d
!
!=======================================================================
!
   subroutine changelogShape2d( array, nnew, array_name, proc_name, shift_, move_right_ )
!             =================   io     i        i           i        i          i
!
!     Change the shape of a 2D array, copying existing elements if possible and filling
!     empty positions with .false.
!
!  Arguments:
!
      logical, allocatable, intent(inout) :: array(:,:)  ! array
      integer(INT4), intent(in) :: nnew(2)               ! new dimensions
      character(*), intent(in) :: array_name             ! array name in the calling procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
      integer(INT4), optional, intent(in) :: shift_(2)         ! shift of data
      logical, optional, intent(in) :: move_right_(2)          ! where true, right-align old
!                                                              !    elements in new array
!  Local Variables:
!
      integer(INT4) :: shift(2), iold(2), jold(2), nold(2), inew(2), jnew(2)
      logical, allocatable :: tmp(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      nold(1) = size( array, 1 )
      nold(2) = size( array, 2 )

      if( present( move_right_ ) ) then
         if( present( shift_ ) ) then
            error stop 'changeLogShape2d: move_right_ and shift must not both be present'
         else
            where( move_right_ )
               shift(:) = nnew(:) - nold(:)
            elsewhere
               shift(:) = 0
            end where
         end if
      else
         if( present( shift_ ) ) then
            shift(:) = shift_(:)
         else
            shift(:) = 0
         end if
      end if

      iold(:) = max( 1, 1-shift(:) )
      inew(:) = max( 1+shift(:), 1)
      jold(:) = min( nold(:), nnew(:)-shift(:) )
      jnew(:) = min( nold(:)+shift(:), nnew(:) )

      call alloc( tmp, [nnew(1), nnew(2)], proc_name, array_name )
!          =====  io           i               i          i
      tmp(:,:) = .false.
      tmp(inew(1):jnew(1),inew(2):jnew(2)) = array(iold(1):jold(1),iold(2):jold(2))

      call move_alloc( tmp, array )

   end subroutine changeLogShape2d
!
!=======================================================================
!
   subroutine changeInt4Shape1d( array, nnew, array_name, proc_name, shift_, move_right_ )
!             =================   io     i        i           i        i          i
!
!     Change the shape of a 1D array, copying existing elements if possible and filling
!     empty positions with ZEROs
!
!  Arguments:
!
      integer(INT4), allocatable, intent(inout) :: array(:) ! array
      integer(INT4), intent(in) :: nnew                     ! new dimension
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
      integer(INT4), optional, intent(in) :: shift_         ! shift of data
      logical, optional, intent(in) :: move_right_          ! if true, right-align old
!                                                           !    elements in new array
!  Local Variables:
!
      integer(INT4) :: shift, iold, jold, nold, inew, jnew
      integer(INT4), allocatable :: tmp(:)
!
!  Start of Executable Code:
!  ------------------------
!
      nold = size( array )

      if( present( move_right_ ) ) then
         if( present( shift_ ) ) then
            error stop 'changeInt4Shape1d: move_right_ and shift must not both be present'
         else
            if( move_right_ ) then
               shift = nnew - nold
            else
               shift = 0
            end if
         end if
      else
         if( present( shift_ ) ) then
            shift = shift_
         else
            shift = 0
         end if
      end if

      iold = max( 1, 1-shift )
      inew = max( 1+shift, 1)
      jold = min( nold, nnew-shift )
      jnew = min( nold+shift, nnew )

      call alloc( tmp, nnew, proc_name, array_name )
!          =====  io    i        i          i
      tmp(:) = 0
      tmp(inew:jnew) = array(iold:jold)

      call move_alloc( tmp, array )

   end subroutine changeInt4Shape1d
!
!=======================================================================
!
   subroutine changeInt4Shape2d( array, nnew, array_name, proc_name, shift_, move_right_ )
!             =================   io     i        i           i        i          i
!
!     Change the shape of a 2D array, copying existing elements if possible and filling
!     empty positions with ZEROs
!
!  Arguments:
!
      integer(INT4), allocatable, intent(inout) :: array(:,:)  ! array
      integer(INT4), intent(in) :: nnew(2)                     ! new dimensions
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
      integer(INT4), optional, intent(in) :: shift_(2)         ! shift of data
      logical, optional, intent(in) :: move_right_(2)          ! where true, right-align old
!                                                              !    elements in new array
!  Local Variables:
!
      integer(INT4) :: shift(2), iold(2), jold(2), nold(2), inew(2), jnew(2)
      integer(INT4), allocatable :: tmp(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      nold(1) = size( array, 1 )
      nold(2) = size( array, 2 )

      if( present( move_right_ ) ) then
         if( present( shift_ ) ) then
            error stop 'changeInt4Shape2d: move_right_ and shift must not both be present'
         else
            where( move_right_ )
               shift(:) = nnew(:) - nold(:)
            elsewhere
               shift(:) = 0
            end where
         end if
      else
         if( present( shift_ ) ) then
            shift(:) = shift_(:)
         else
            shift(:) = 0
         end if
      end if

      iold(:) = max( 1, 1-shift(:) )
      inew(:) = max( 1+shift(:), 1)
      jold(:) = min( nold(:), nnew(:)-shift(:) )
      jnew(:) = min( nold(:)+shift(:), nnew(:) )

      call alloc( tmp, [nnew(1), nnew(2)], proc_name, array_name )
!          =====  io           i               i          i
      tmp(:,:) = 0
      tmp(inew(1):jnew(1),inew(2):jnew(2)) = array(iold(1):jold(1),iold(2):jold(2))

      call move_alloc( tmp, array )

   end subroutine changeInt4Shape2d
!
!=======================================================================
!
   subroutine changeReal8Shape1d( array, nnew, array_name, proc_name, shift_, move_right_ )
!             ==================   io     i        i           i        i          i
!
!     Change the shape of a 1D array, copying existing elements if possible and filling
!     empty positions with ZEROs
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: array(:)   ! array
      integer(INT4), intent(in) :: nnew                     ! new dimension
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
      integer(INT4), optional, intent(in) :: shift_         ! shift of data
      logical, optional, intent(in) :: move_right_          ! if true, right-align old
!                                                           !    elements in new array
!  Local Variables:
!
      integer(INT4) :: shift, iold, jold, nold, inew, jnew
      real(REAL8), allocatable :: tmp(:)
!
!  Start of Executable Code:
!  ------------------------
!
      nold = size( array )

      if( present( move_right_ ) ) then
         if( present( shift_ ) ) then
            error stop 'changeReal8Shape1d: move_right_ and shift must not both be present'
         else
            if( move_right_ ) then
               shift = nnew - nold
            else
               shift = 0
            end if
         end if
      else
         if( present( shift_ ) ) then
            shift = shift_
         else
            shift = 0
         end if
      end if

      iold = max( 1, 1-shift )
      inew = max( 1+shift, 1)
      jold = min( nold, nnew-shift )
      jnew = min( nold+shift, nnew )

      call alloc( tmp, nnew, proc_name, array_name )
!          =====  io    i        i          i
      tmp(:) = ZERO
      tmp(inew:jnew) = array(iold:jold)

      call move_alloc( tmp, array )

   end subroutine changeReal8Shape1d
!
!=======================================================================
!
   subroutine changeReal8Shape2d( array, nnew, array_name, proc_name, shift_, move_right_ )
!             ==================   io     i        i           i        i          i
!
!     Change the shape of a 2D array, copying existing elements if possible and filling
!     empty positions with ZEROs
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: array(:,:) ! array
      integer(INT4), intent(in) :: nnew(2)                  ! new dimensions
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
      integer(INT4), optional, intent(in) :: shift_(2)      ! shift of data
      logical, optional, intent(in) :: move_right_(2)       ! where true, right-align old
!                                                           !    elements in new array
!  Local Variables:
!
      integer(INT4) :: shift(2), iold(2), jold(2), nold(2), inew(2), jnew(2)
      real(REAL8), allocatable :: tmp(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      nold(1) = size( array, 1 )
      nold(2) = size( array, 2 )

      if( present( move_right_ ) ) then
         if( present( shift_ ) ) then
            error stop 'changeReal8Shape2d: move_right_ and shift must not both be present'
         else
            where( move_right_ )
               shift(:) = nnew(:) - nold(:)
            elsewhere
               shift(:) = 0
            end where
         end if
      else
         if( present( shift_ ) ) then
            shift(:) = shift_(:)
         else
            shift(:) = 0
         end if
      end if

      iold(:) = max( 1, 1-shift(:) )
      inew(:) = max( 1+shift(:), 1)
      jold(:) = min( nold(:), nnew(:)-shift(:) )
      jnew(:) = min( nold(:)+shift(:), nnew(:) )

      call alloc( tmp, [nnew(1), nnew(2)], proc_name, array_name )
!          =====  io           i               i          i
      tmp(:,:) = ZERO
      tmp(inew(1):jnew(1),inew(2):jnew(2)) = array(iold(1):jold(1),iold(2):jold(2))

      call move_alloc( tmp, array )

   end subroutine changeReal8Shape2d
!
!=======================================================================
!
   subroutine changeReal8Shape3d( array, nnew, array_name, proc_name, shift_, move_right_ )
!             ==================   io     i        i           i        i          i
!
!     Change the shape of a 3D array, copying existing elements if possible and filling
!     empty positions with ZEROs
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: array(:,:,:)  ! array
      integer(INT4), intent(in) :: nnew(3)                     ! new dimensions
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
      integer(INT4), optional, intent(in) :: shift_(3)         ! shift of data
      logical, optional, intent(in) :: move_right_(3)          ! where true, right-align old
!                                                              !    elements in new array
!  Local Variables:
!
      integer(INT4) :: shift(3), iold(3), jold(3), nold(3), inew(3), jnew(3)
      real(REAL8), allocatable :: tmp(:,:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      nold(1) = size( array, 1 )
      nold(2) = size( array, 2 )
      nold(3) = size( array, 3 )

      if( present( move_right_ ) ) then
         if( present( shift_ ) ) then
            error stop 'changeReal8Shape3d: move_right_ and shift must not both be present'
         else
            where( move_right_ )
               shift(:) = nnew(:) - nold(:)
            elsewhere
               shift(:) = 0
            end where
         end if
      else
         if( present( shift_ ) ) then
            shift(:) = shift_(:)
         else
            shift(:) = 0
         end if
      end if

      iold(:) = max( 1, 1-shift(:) )
      inew(:) = max( 1+shift(:), 1)
      jold(:) = min( nold(:), nnew(:)-shift(:) )
      jnew(:) = min( nold(:)+shift(:), nnew(:) )

      call alloc( tmp, [nnew(1), nnew(2), nnew(3)], proc_name, array_name )
!          =====  io               i                    i          i
      tmp(:,:,:) = ZERO
      tmp(inew(1):jnew(1),inew(2):jnew(2),inew(3):jnew(3))  &
            = array(iold(1):jold(1),iold(2):jold(2),iold(3):jold(3))

      call move_alloc( tmp, array )

   end subroutine changeReal8Shape3d
!
!=======================================================================
!
   subroutine changeReal8Shape4d( array, nnew, array_name, proc_name, shift_, move_right_ )
!             ==================   io     i        i           i        i          i
!
!     Change the shape of a 4D array, copying existing elements if possible and filling
!     empty positions with ZEROs
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: array(:,:,:,:)! array
      integer(INT4), intent(in) :: nnew(4)                     ! new dimensions
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
      integer(INT4), optional, intent(in) :: shift_(4)         ! shift of data
      logical, optional, intent(in) :: move_right_(4)          ! where true, right-align old
!                                                              !    elements in new array
!  Local Variables:
!
      integer(INT4) :: shift(4), iold(4), jold(4), nold(4), inew(4), jnew(4)
      real(REAL8), allocatable :: tmp(:,:,:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      nold(1) = size( array, 1 )
      nold(2) = size( array, 2 )
      nold(3) = size( array, 3 )
      nold(4) = size( array, 4 )

      if( present( move_right_ ) ) then
         if( present( shift_ ) ) then
            error stop 'changeReal8Shape4d: move_right_ and shift must not both be present'
         else
            where( move_right_ )
               shift(:) = nnew(:) - nold(:)
            elsewhere
               shift(:) = 0
            end where
         end if
      else
         if( present( shift_ ) ) then
            shift(:) = shift_(:)
         else
            shift(:) = 0
         end if
      end if

      iold(:) = max( 1, 1-shift(:) )
      inew(:) = max( 1+shift(:), 1)
      jold(:) = min( nold(:), nnew(:)-shift(:) )
      jnew(:) = min( nold(:)+shift(:), nnew(:) )

      call alloc( tmp, [nnew(1), nnew(2), nnew(3), nnew(4)], proc_name, array_name )
!          =====  io                    i                        i          i
      tmp(:,:,:,:) = ZERO
      tmp(inew(1):jnew(1),inew(2):jnew(2),inew(3):jnew(3),inew(4):jnew(4))  &
            = array(iold(1):jold(1),iold(2):jold(2),iold(3):jold(3),iold(4):jold(4))

      call move_alloc( tmp, array )

   end subroutine changeReal8Shape4d
!
!=======================================================================
!
   subroutine changeIposShape1d( array, nnew, array_name, proc_name, shift_, move_right_ )
!             =================   io     i        i           i        i          i
!
!     Change the shape of a 1D array, copying existing elements if possible and leaving
!     new elements undefined
!
!  Arguments:
!
      type(Ipos_type), allocatable, intent(inout) :: array(:)  ! array
      integer(INT4), intent(in) :: nnew                        ! new dimension
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
      integer(INT4), optional, intent(in) :: shift_            ! shift of data
      logical, optional, intent(in) :: move_right_             ! if true, right-align old
!                                                              !    elements in new array
!  Local Variables:
!
      integer(INT4) :: shift, iold, jold, nold, inew, jnew
      type(Ipos_type), allocatable :: tmp(:)
!
!  Start of Executable Code:
!  ------------------------
!
      nold = size( array )

      if( present( move_right_ ) ) then
         if( present( shift_ ) ) then
            error stop 'changeIposShape1d: move_right_ and shift must not both be present'
         else
            if( move_right_ ) then
               shift = nnew - nold
            else
               shift = 0
            end if
         end if
      else
         if( present( shift_ ) ) then
            shift = shift_
         else
            shift = 0
         end if
      end if

      iold = max( 1, 1-shift )
      inew = max( 1+shift, 1)
      jold = min( nold, nnew-shift )
      jnew = min( nold+shift, nnew )

      call alloc( tmp, nnew, proc_name, array_name )
!          =====  io    i        i          i
      tmp(inew:jnew) = array(iold:jold)

      call move_alloc( tmp, array )
!
   end subroutine changeIposShape1d
!
!=======================================================================
!
   subroutine changeIposShape2d( array, nnew, array_name, proc_name, shift_, move_right_ )
!             =================   io     i        i           i        i          i
!
!     Change the shape of a 2D array, copying existing elements if possible and leaving
!     new elements undefined
!
!  Arguments:
!
      type(Ipos_type), allocatable, intent(inout) :: array(:,:)! array
      integer(INT4), intent(in) :: nnew(2)                     ! new dimensions
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
      integer(INT4), optional, intent(in) :: shift_(2)         ! shift of data
      logical, optional, intent(in) :: move_right_(2)          ! where true, right-align old
!                                                              !    elements in new array
!  Local Variables:
!
      integer(INT4) :: shift(2), iold(2), jold(2), nold(2), inew(2), jnew(2)
      type(Ipos_type), allocatable :: tmp(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      nold(1) = size( array, 1 )
      nold(2) = size( array, 2 )

      if( present( move_right_ ) ) then
         if( present( shift_ ) ) then
            error stop 'changeIposShape2d: move_right_ and shift must not both be present'
         else
            where( move_right_ )
               shift(:) = nnew(:) - nold(:)
            elsewhere
               shift(:) = 0
            end where
         end if
      else
         if( present( shift_ ) ) then
            shift(:) = shift_(:)
         else
            shift(:) = 0
         end if
      end if

      iold(:) = max( 1, 1-shift(:) )
      inew(:) = max( 1+shift(:), 1)
      jold(:) = min( nold(:), nnew(:)-shift(:) )
      jnew(:) = min( nold(:)+shift(:), nnew(:) )

      call alloc( tmp, [nnew(1), nnew(2)], proc_name, array_name )
!          =====  io           i               i          i
      tmp(inew(1):jnew(1),inew(2):jnew(2)) = array(iold(1):jold(1),iold(2):jold(2))

      call move_alloc( tmp, array )

   end subroutine changeIposShape2d
!
!=======================================================================
!
   subroutine changeProjShape1d( array, nnew, array_name, proc_name, shift_, move_right_ )
!             =================   io     i        i           i        i          i
!
!     Change the shape of a 1D array, copying existing elements if possible and leaving
!     empty positions undefined
!
!  Arguments:
!
      type(Proj_type), allocatable, intent(inout) :: array(:)  ! array
      integer(INT4), intent(in) :: nnew                        ! new dimension
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
      integer(INT4), optional, intent(in) :: shift_            ! shift of data
      logical, optional, intent(in) :: move_right_             ! if true, right-align old
!                                                              !    elements in new array
!  Local Variables:
!
      integer(INT4) :: shift, iold, jold, nold, inew, jnew
      type(Proj_type), allocatable :: tmp(:)
!
!  Start of Executable Code:
!  ------------------------
!
      nold = size( array )

      if( present( move_right_ ) ) then
         if( present( shift_ ) ) then
            error stop 'changeProjShape1d: move_right_ and shift must not both be present'
         else
            if( move_right_ ) then
               shift = nnew - nold
            else
               shift = 0
            end if
         end if
      else
         if( present( shift_ ) ) then
            shift = shift_
         else
            shift = 0
         end if
      end if

      iold = max( 1, 1-shift )
      inew = max( 1+shift, 1)
      jold = min( nold, nnew-shift )
      jnew = min( nold+shift, nnew )

      call alloc( tmp, nnew, proc_name, array_name )
!          =====  io    i        i          i
      tmp(inew:jnew) = array(iold:jold)

      call move_alloc( tmp, array )

   end subroutine changeProjShape1d
!
!=======================================================================
!
   subroutine addReflect1dLog( array, direction, reflectbnd, array_name, proc_name )
!             ===============   io        i          i           i           i
!
!     Extend array by the reflection of itself in given direction.
!
!  Arguments:
!
      logical, allocatable, intent(inout) :: array(:) ! array
      integer(INT4), intent(in) :: direction          ! flag indicating reflection about
                                                      !    +1: upper boundary
                                                      !    -1: lower boundary
                                                      !     0: no reflection
      logical, intent(in) :: reflectbnd               ! flag indicating whether the boundary
                                                      !    is reflected
      character(*), intent(in) :: array_name          ! array name in the calling procedure
      character(*), intent(in) :: proc_name           ! name of calling procedure elements in
!                                                     !    new array
!  Local Variables:
!
      integer(INT4) :: nold, nnew
!
!  Start of Executable Code:
!  ------------------------
!
      if( direction == 0 )  return

      nold = size( array )                      ! old and new size of array
      if( reflectbnd ) then
         nnew = 2*nold
      else
         nnew = 2*nold - 1
      end if

      call changeShape( array, nnew, array_name, proc_name, move_right_=(direction<0) )
!          ===========   io     i        i           i                i

      if( direction > 0 ) then                  ! reflect about upper boundary
         array(nnew-nold+1:nnew) = array(1:nold)
      else                                      ! reflect about lower boundary
         array(1:nold) = array(nnew-nold+1:nnew)
      end if

   end subroutine addReflect1dLog
!
!=======================================================================
!
   subroutine addReflect2dLog( array, direction, reflectbnd, array_name, proc_name )
!             ===============   io        i          i           i           i
!
!     Extend array by the reflection of itself in given direction.
!
!  Arguments:
!
      logical, allocatable, intent(inout) :: array(:,:)  ! array
      integer(INT4), intent(in) :: direction(2)          ! flag indicating reflection about
                                                         !    +1: upper boundary
                                                         !    -1: lower boundary
                                                         !     0: no reflection
      logical, intent(in) :: reflectbnd(2)               ! flag indicating whether the
                                                         !    boundary is reflected
      character(*), intent(in) :: array_name             ! array name in the calling
                                                         !    procedure
      character(*), intent(in) :: proc_name              ! name of calling procedure
!                                                        !    elements in new array
!  Local Variables:
!
      integer(INT4) :: nold(2), nnew(2)
!
!  Start of Executable Code:
!  ------------------------
!
      nold(:) = shape( array )
      where( direction(:) /= 0 )
         where( reflectbnd(:) )
            nnew(:) = 2*nold(:)
         else where
            nnew(:) = 2*nold(:) - 1
         end where
      else where
         nnew(:) = nold(:)
      end where

      call changeShape( array, nnew, array_name, proc_name, move_right_=(direction<0) )
!          ===========   io     i        i           i                i

      if( direction(1) > 0 ) then               ! reflect about upper boundary along 1st axis
         array(nnew(1)-nold(1)+1:nnew(1), 1:nold(2)) = array(1:nold(1), 1:nold(2))
      else if( direction(1) < 0 ) then          ! reflect about lower boundary along 1st axis
         array(1:nold(1), 1:nold(2)) = array(nnew(1)-nold(1)+1:nnew(1), 1:nold(2))
      end if

      if( direction(2) > 0 ) then               ! reflect about upper boundary along 2nd axis
         array(:, nnew(2)-nold(2)+1:nnew(2)) = array(:, 1:nold(2))
      else if( direction(2) < 0 ) then          ! reflect about lower boundary along 2nd axis
         array(:, 1:nold(2)) = array(:, nnew(2)-nold(2)+1:nnew(2))
      end if

   end subroutine addReflect2dLog
!
!=======================================================================
!
   subroutine addReflect1dReal8( array, direction, reflectbnd, array_name, proc_name )
!             =================   io        i          i           i           i
!
!     Extend array by the reflection of itself in given direction.
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: array(:)   ! array
      integer(INT4), intent(in) :: direction                ! flag indicating reflection about
                                                            !    +1: upper boundary
                                                            !    -1: lower boundary
                                                            !     0: no reflection
      logical, intent(in) :: reflectbnd                     ! flag indicating whether the
                                                            !    boundary is reflected
      character(*), intent(in) :: array_name                ! array name in the calling
                                                            !    procedure
      character(*), intent(in) :: proc_name                 ! name of calling procedure
!                                                           !    elements in new array
!  Local Variables:
!
      integer(INT4) :: nold, nnew
!
!  Start of Executable Code:
!  ------------------------
!
      if( direction == 0 )  return

      nold = size( array )                      ! old and new size of array
      if( reflectbnd ) then
         nnew = 2*nold
      else
         nnew = 2*nold - 1
      end if

      call changeShape( array, nnew, array_name, proc_name, move_right_=(direction<0) )
!          ===========   io     i        i           i                i

      if( direction > 0 ) then                  ! reflect about upper boundary
         array(nnew-nold+1:nnew) = array(1:nold)
      else                                      ! reflect about lower boundary
         array(1:nold) = array(nnew-nold+1:nnew)
      end if

   end subroutine addReflect1dReal8
!
!=======================================================================
!
   subroutine addReflect2dReal8( array, direction, reflectbnd, array_name, proc_name )
!             =================   io        i          i           i           i
!
!     Extend array by the reflection of itself in given direction.
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: array(:,:) ! array
      integer(INT4), intent(in) :: direction(2)             ! flag indicating reflection about
                                                            !    +1: upper boundary
                                                            !    -1: lower boundary
                                                            !     0: no reflection
      logical, intent(in) :: reflectbnd(2)                  ! flag indicating whether the
                                                            !    boundary is reflected
      character(*), intent(in) :: array_name                ! array name in the calling
                                                            !    procedure
      character(*), intent(in) :: proc_name                 ! name of calling procedure
!                                                           !    elements in new array
!  Local Variables:
!
      integer(INT4) :: nold(2), nnew(2)
!
!  Start of Executable Code:
!  ------------------------
!
      nold(:) = shape( array )
      where( direction(:) /= 0 )
         where( reflectbnd(:) )
            nnew(:) = 2*nold(:)
         else where
            nnew(:) = 2*nold(:) - 1
         end where
      else where
         nnew(:) = nold(:)
      end where

      call changeShape( array, nnew, array_name, proc_name, move_right_=(direction<0) )
!          ===========   io     i        i           i                i

      if( direction(1) > 0 ) then               ! reflect about upper boundary along 1st axis
         array(nnew(1)-nold(1)+1:nnew(1), 1:nold(2)) = array(1:nold(1), 1:nold(2))
      else if( direction(1) < 0 ) then          ! reflect about lower boundary along 1st axis
         array(1:nold(1), 1:nold(2)) = array(nnew(1)-nold(1)+1:nnew(1), 1:nold(2))
      end if

      if( direction(2) > 0 ) then               ! reflect about upper boundary along 2nd axis
         array(:, nnew(2)-nold(2)+1:nnew(2)) = array(:, 1:nold(2))
      else if( direction(2) < 0 ) then          ! reflect about lower boundary along 2nd axis
         array(:, 1:nold(2)) = array(:, nnew(2)-nold(2)+1:nnew(2))
      end if

   end subroutine addReflect2dReal8
!
!=======================================================================
!
   subroutine addReflect3dReal8( array, direction, reflectbnd, array_name, proc_name )
!             =================   io        i          i           i           i
!
!     Extend array by the reflection of itself in given direction.
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: array(:,:,:)  ! array
      integer(INT4), intent(in) :: direction(3)                ! flag indicating reflection
                                                               !    +1: about upper boundary
                                                               !    -1: about lower boundary
                                                               !     0: no reflection
      logical, intent(in) :: reflectbnd(3)                     ! flag indicating whether the
                                                               !    boundary is reflected
      character(*), intent(in) :: array_name                   ! array name in the calling
                                                               !    procedure
      character(*), intent(in) :: proc_name                    ! name of calling procedure
!                                                              !    elements in new array
!  Local Variables:
!
      integer(INT4) :: nold(3), nnew(3)
!
!  Start of Executable Code:
!  ------------------------
!
      nold(:) = shape( array )
      where( direction(:) /= 0 )
         where( reflectbnd(:) )
            nnew(:) = 2*nold(:)
         else where
            nnew(:) = 2*nold(:) - 1
         end where
      else where
         nnew(:) = nold(:)
      end where

      call changeShape( array, nnew, array_name, proc_name, move_right_=(direction<0) )
!          ===========   io     i        i           i                i

      if( direction(1) > 0 ) then               ! reflect about upper boundary along 1st axis
         array(nnew(1)-nold(1)+1:nnew(1), 1:nold(2), 1:nold(3)) =  &
               array(1:nold(1), 1:nold(2), 1:nold(3))
      else if( direction(1) < 0 ) then          ! reflect about lower boundary along 1st axis
         array(1:nold(1), 1:nold(2), 1:nold(3)) =  &
               array(nnew(1)-nold(1)+1:nnew(1), 1:nold(2), 1:nold(3))
      end if

      if( direction(2) > 0 ) then               ! reflect about upper boundary along 2nd axis
         array(:, nnew(2)-nold(2)+1:nnew(2), 1:nold(3)) =  &
               array(:, 1:nold(2), 1:nold(3))
      else if( direction(2) < 0 ) then          ! reflect about lower boundary along 2nd axis
         array(:, 1:nold(2), 1:nold(3)) =  &
               array(:, nnew(2)-nold(2)+1:nnew(2), 1:nold(3))
      end if

      if( direction(3) > 0 ) then               ! reflect about upper boundary along 2nd axis
         array(:, :, nnew(3)-nold(3)+1:nnew(3)) =  &
               array(:, :, 1:nold(3))
      else if( direction(3) < 0 ) then          ! reflect about lower boundary along 2nd axis
         array(:, :, 1:nold(3)) =  &
               array(:, :, nnew(3)-nold(3)+1:nnew(3))
      end if

   end subroutine addReflect3dReal8
!
!=======================================================================
!
   subroutine addReflect4dReal8( array, direction, reflectbnd, array_name, proc_name )
!             =================   io        i          i           i           i
!
!     Extend array by the reflection of itself in given direction.
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: array(:,:,:,:)   ! array
      integer(INT4), intent(in) :: direction(4)                   ! flag indicating reflection
                                                                  !    +1: about upper boundary
                                                                  !    -1: about lower boundary
                                                                  !     0: no reflection
      logical, intent(in) :: reflectbnd(4)                        ! flag indicating whether the
                                                                  !    boundary is reflected
      character(*), intent(in) :: array_name                      ! array name in the calling
                                                                  !    procedure
      character(*), intent(in) :: proc_name                       ! name of calling procedure
!                                                                 !    elements in new array
!  Local Variables:
!
      integer(INT4) :: nold(4), nnew(4)
!
!  Start of Executable Code:
!  ------------------------
!
      nold(:) = shape( array )
      where( direction(:) /= 0 )
         where( reflectbnd(:) )
            nnew(:) = 2*nold(:)
         else where
            nnew(:) = 2*nold(:) - 1
         end where
      else where
         nnew(:) = nold(:)
      end where

      call changeShape( array, nnew, array_name, proc_name, move_right_=(direction<0) )
!          ===========   io     i        i           i                i

      if( direction(1) > 0 ) then               ! reflect about upper boundary along 1st axis
         array(nnew(1)-nold(1)+1:nnew(1), 1:nold(2), 1:nold(3), 1:nold(4)) =  &
               array(1:nold(1), 1:nold(2), 1:nold(3), 1:nold(4))
      else if( direction(1) < 0 ) then          ! reflect about lower boundary along 1st axis
         array(1:nold(1), 1:nold(2), 1:nold(3), 1:nold(4)) =  &
               array(nnew(1)-nold(1)+1:nnew(1), 1:nold(2), 1:nold(3), 1:nold(4))
      end if

      if( direction(2) > 0 ) then               ! reflect about upper boundary along 2nd axis
         array(:, nnew(2)-nold(2)+1:nnew(2), 1:nold(3), 1:nold(4)) =  &
               array(:, 1:nold(2), 1:nold(3), 1:nold(4))
      else if( direction(2) < 0 ) then          ! reflect about lower boundary along 2nd axis
         array(:, 1:nold(2), 1:nold(3), 1:nold(4)) =  &
               array(:, nnew(2)-nold(2)+1:nnew(2), 1:nold(3), 1:nold(4))
      end if

      if( direction(3) > 0 ) then               ! reflect about upper boundary along 2nd axis
         array(:, :, nnew(3)-nold(3)+1:nnew(3), 1:nold(4)) =  &
               array(:, :, 1:nold(3), 1:nold(4))
      else if( direction(3) < 0 ) then          ! reflect about lower boundary along 2nd axis
         array(:, :, 1:nold(3), 1:nold(4)) =  &
               array(:, :, nnew(3)-nold(3)+1:nnew(3), 1:nold(4))
      end if

      if( direction(4) > 0 ) then               ! reflect about upper boundary along 2nd axis
         array(:, :, :, nnew(4)-nold(4)+1:nnew(4)) =  &
               array(:, :, :, 1:nold(4))
      else if( direction(4) < 0 ) then          ! reflect about lower boundary along 2nd axis
         array(:, :, :,  1:nold(4)) =  &
               array(:, :,  :, nnew(4)-nold(4)+1:nnew(4))
      end if

   end subroutine addReflect4dReal8
!
!=======================================================================
!
   subroutine check_allocation( proc_name, array_name, istat, errmsg )
!             ================      i          i         i      i
!
!     Check success of allocation and abort program if allocation was not successful
!
!     To be used immeadiately after an allocation statement
!
!  Arguments:
!
      character(*), intent(in) :: proc_name  ! name of procedure where allocation was issued
      character(*), intent(in) :: array_name ! name of array attempted to be allocated
      integer(INT4), intent(in) :: istat     ! status returned by the allocate statement
      character(*), intent(in) :: errmsg     ! error message returned by the allocate statement
!
!  Local Variables:
!
!    None
!
!  Start of Executable Code:
!  ------------------------
!
      if( istat /= 0 ) then
         print *, proc_name // ': Error allocating ' // trim( array_name )
         print *, 'Status =', istat
         print *, trim( errmsg )
         error stop
      end if

   end subroutine check_allocation
!
!=======================================================================
!
end module Array_mod

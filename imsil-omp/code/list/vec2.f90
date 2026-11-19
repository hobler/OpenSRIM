module List_Vec2_mod
!
!  Defines a doubly linked list of REAL8 vectors of size 2, using a sentinel 
!  to mark the beginning and the end of the list.
!
!  The list has defined a current position.
!  Reference given in list_int4_mod.f90
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
!  Public Procedures:
!  =================
!
!  Initialization
!
   public :: newVec2List                     ! Initializes the list.
!            ===========
!
!  Managing current position
!
   public :: setPosVec2List                ! Sets current position
!            ==============
   public :: getPosVec2List                ! Gets current position
!            ==============      
   public :: rewindVec2List                ! Moves current position before first value.
!            ==============
   public :: windVec2List                  ! Moves current position to the last value.
!            ============
   public :: iterateVec2List               ! Moves current position forward.
!            ===============
!
!  Managing list contents
!
   public :: getValVec2List                ! Gets value at current position.
!            ==============
   public :: insertVec2List                ! Inserts value after current position and
!            ==============                !   moves current position to new value.
   public :: removeVec2List                ! Removes value at current position and moves
!            ==============                !   current position before removed value.
   public :: clearVec2List                 ! Clears the whole list.
!            =============
!
!  Getting information
!
   public :: sizeOfVec2List                ! Returns the number of values in the list.
!            ==============
   public :: dumpVec2List                  ! Dumps the list to standard output.
!            ============
!
!  Public Data Types:
!  =================
!
   type, public :: Vec2Node                     ! list node
!                  ========
      private
      real(REAL8) :: value(2)
      type(Vec2Node), pointer :: next
      type(Vec2Node), pointer :: prev
   end type Vec2Node

   type, public :: Vec2List                     ! list for Vec2D numbers
!                  ========
      private
      type(Vec2Node), pointer :: nil => null( ) ! sentinel:
                                                !   nil%next => first node
                                                !   nil%prev => last node
                                                !   nil%value ... undefined
      type(Vec2Node), pointer :: pos => null( ) ! pointer to current position
      integer(INT4) :: size = 0
   end type Vec2List

contains
!
!=========================================================================================
!  subroutine for initialization
!=========================================================================================
!
   subroutine newVec2List( list )
!             ===========   io
!
!     Initializes the list. Must be called before the list can be used.
!
!  Arguments:
!
      type(Vec2List), intent(inout) :: list     ! empty list
!
!  Start of Executable Code:
!  ------------------------
!
      if( associated( list%nil ) )  error stop 'newVec2List: list already initialized'
      
      allocate( list%nil )
      list%nil%next => list%nil
      list%nil%prev => list%nil

      list%pos => list%nil

      list%size = 0
      
   end subroutine newVec2List
!
!=========================================================================================
!  subroutines for managing the current position
!=========================================================================================
!
   subroutine setPosVec2List( list, pos )
!             ==============   io    i
!
!     Put current position at pos.
!
!     WARNING: Be careful and make sure that pos is a node of list. This cannot easily
!              be checked within this subroutine!
!
!  Arguments:
!
      type(Vec2List), intent(inout) :: list     ! list
      type(Vec2Node), pointer :: pos            ! position
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      list%pos => pos

   end subroutine setPosVec2List
!
!=========================================================================================
!
   subroutine getPosVec2List( list, pos )
!             ================    i    o
!
!     Get current position
!
!  Arguments:
!
      type(Vec2List), intent(in) :: list        ! list
      type(Vec2Node), pointer :: pos            ! position
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      pos => list%pos

   end subroutine getPosVec2List
!
!=========================================================================================
!
   subroutine rewindVec2List( list )
!             ==============   io
!
!     Put current position at sentinel (nil)
!
!  Arguments:
!
      type(Vec2List), intent(inout) :: list     ! list
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      list%pos => list%nil

   end subroutine rewindVec2List
!
!=========================================================================================
!
   subroutine windVec2List( list )
!             ==============   io
!
!     Moves current position to the last node
!
!  Arguments:
!
      type(Vec2List), intent(inout) :: list     ! list
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      list%pos => list%nil%prev

   end subroutine windVec2List
!
!=========================================================================================
!
   subroutine iterateVec2List( list, isAtEnd )
!             ==============   io      o
!
!     Iterate current position and read value
!
!  Arguments:
!
      type(Vec2List), intent(inout) :: list           ! list
      logical, intent(out) :: isAtEnd                 ! flag indicating current position is at
!                                                     !   end of list or list is empty
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      isAtEnd = associated( list%pos%next, list%nil )
      if( isAtEnd ) return

      list%pos => list%pos%next

   end subroutine iterateVec2List
!
!=========================================================================================
!  subroutines for managing the list contents
!=========================================================================================
!
   subroutine getValVec2List( list, value, isAtBeg_ )
!             ==============   io     o       o
!
!     Return value at current position
!
!  Arguments:
!
      type(Vec2List), intent(inout) :: list     ! list
      real(REAL8), intent(out) :: value(2)      ! value read
      logical, intent(out), optional :: isAtBeg_! flag indicating current position is before
!                                               !   first value of list or list is empty
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( isAtBeg_ ) ) then
         isAtBeg_ = associated( list%pos, list%nil )
         if( isAtBeg_ ) return
      else
         if( associated( list%pos, list%nil ) )  &
               error stop 'getValVec2List: current position is before first value'
      end if
      
      value = list%pos%value

   end subroutine getValVec2List
!
!=========================================================================================
!
   subroutine insertVec2List( list, value )
!             =============   io     i
!
!     Inserts value after current value and make new value current
!
!  Arguments:
!
      type(Vec2List), intent(inout) :: list     ! list 
      real(REAL8), intent(in) :: value(2)       ! value to be inserted
!
!  Local Variables:
!
      type(Vec2Node), pointer :: p             ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( p )
      p%value = value
      p%next => list%pos%next
      p%prev => list%pos
      p%next%prev => p
      p%prev%next => p

      list%pos => p

      list%size = list%size + 1

   end subroutine insertVec2List
!
!=========================================================================================
!
   subroutine removeVec2List( list, isAtBeg_ )
!             ==============   io      o
!
!     Remove current element from the list and move current position before removed value
!
!  Arguments:
!
      type(Vec2List), intent(inout) :: list     ! list
      logical, intent(out), optional :: isAtBeg_ ! flag indicating current position is before
!                                               !   first value of list or list is empty
!
!  Local Variables:
!
      type(Vec2Node), pointer :: p
!
!  Start of Executable code:
!  ------------------------
!
      if( present( isAtBeg_ ) ) then
         isAtBeg_ = associated( list%pos, list%nil )
         if( isAtBeg_ ) return
      else
         if( associated( list%pos, list%nil ) )  &
               error stop 'removeVec2List: current position is before first value'
      end if
      
      p => list%pos
      p%next%prev => p%prev
      p%prev%next => p%next
      
      list%pos => p%prev
      
      deallocate( p )

      list%size = list%size - 1

   end subroutine removeVec2List
!
!=========================================================================================
!
   subroutine clearVec2List( list )
!             ===============   io
!
!     Clears the whole list
!
!  Arguments:
!
      type(Vec2List), intent(inout) :: list
!
!  Local Variables:
!
      logical :: isAtEnd
!
!  Start of Executable code:
!  ------------------------
!
      call rewindVec2List( list )
!          ================   io
      do
         call iterateVec2List( list, isAtEnd )
!             =================   io      o
         if( isAtEnd )  exit
!
         call removeVec2List( list )
!             ================   io
      end do
!
   end subroutine clearVec2List
!
!=========================================================================================
!  subroutines for getting information about the list
!=========================================================================================
!
   integer(INT4) function sizeOfVec2List( list )
!                         ==============   i  
!
!     Returns the size of the list, i.e. the number of items
!
!  Arguments:
!  
      type(Vec2List), intent(in) :: list      ! list
!
!  Start of Executable Code:
!  ------------------------
!
      sizeOfVec2List = list%size

   end function sizeOfVec2List
!
!=========================================================================================
!
   subroutine dumpVec2List( list )
!             ============   io
!
!     Dumps the whole list to standard output
!
!  Arguments:
!
      type(Vec2List), intent(inout) :: list      ! list
!
!  Local Veriables:
!
      real(REAL8) :: value(2)
      logical :: isAtEnd
      type(Vec2Node), pointer :: p
!
!  Start of Executable Code:
!  ------------------------
!
      p => list%pos                            ! remember current position

      call rewindVec2List( list )
!          ==============   io
      print *,  'List contents:'
      do
         call iterateVec2List( list, isAtEnd )
!             ===============   io     o
         if( isAtEnd )  exit
!
         call getValVec2List( list, value )
!             ==============   io     o         
         print *,  value
      end do
      print *,  'End of list contents.'

      list%pos => p                            ! restore current position

   end subroutine dumpVec2List
!
!=========================================================================================
!
end module List_Vec2_mod

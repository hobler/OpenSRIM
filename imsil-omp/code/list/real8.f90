module List_REAL8_mod
!
!  Defines a doubly linked list of data type REAL4, using a sentinel to mark
!  the beginning and the end of the list.
!
!  Ref: T.H.Cormen et al., Introduction to Algorithms, 2nd Ed., MIT Press, 2001, pp.204ff
!
!  The list has defined a current position.
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
   public :: newReal8List                   ! Initializes the list.
!            ============
!
!  Managing current position
!
   public :: setPosReal8List                ! Sets current position
!            ===============
   public :: getposReal8List                ! Gets current position
!            ===============      
   public :: rewindReal8List                ! Moves current position before first value.
!            ===============
   public :: windReal8List                  ! Moves current position to the last value.
!            =============
   public :: iterateReal8List               ! Moves current position forward.
!            ===============
!
!  Managing list contents
!
   public :: getValReal8List                ! Gets value at current position.
!            ===============
   public :: insertReal8List                ! Inserts value after current position and
!            ===============                !   moves current position to new value.
   public :: removeReal8List                ! Removes value at current position and moves
!            ===============                !   current position before removed value.
   public :: clearReal8List                 ! Clears the whole list.
!            ==============
!
!  Getting information
!
   public :: sizeOfReal8List                ! Returns the number of values in the list.
!            ===============
   public :: dumpReal8List                  ! Dumps the list to standard output.
!            =============
!
!  Public Data Types:
!  =================
!
   type, public :: REAL8Node                    ! list node
!                  =========
      private
      real(REAL8) :: value
      type(REAL8Node), pointer :: next
      type(REAL8Node), pointer :: prev
   end type REAL8Node

   type, public :: REAL8List                     ! list for REAL8 numbers
!                  =========
      private
      type(REAL8Node), pointer :: nil => null( ) ! sentinel:
                                                 !   nil%next => first node
                                                 !   nil%prev => last node
                                                 !   nil%value ... undefined
      type(REAL8Node), pointer :: pos => null( ) ! pointer to current position
      integer(INT4) :: size = 0
   end type REAL8List
!
contains
!
!=======================================================================
!  subroutine for initialization
!=======================================================================
!
   subroutine newReal8List( list )
!             ============   io
!
!     Initializes the list. Must be called before the list can be used.
!
!  Arguments:
!
      type(REAL8List), intent(inout) :: list     ! empty list
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( associated( list%nil ) )  error stop 'newReal8List: list already initialized'
      
      allocate( list%nil )
      list%nil%next => list%nil
      list%nil%prev => list%nil

      list%pos => list%nil

      list%size = 0
      
   end subroutine newReal8List
!
!=======================================================================
!  subroutines for managing the current position
!=======================================================================
!
   subroutine setPosReal8List( list, pos )
!             ===============   io    i
!
!     Put current position at pos.
!
!     WARNING: Be careful and make sure that pos is a node of list. This cannot easily
!              be checked within this subroutine!
!
!  Arguments:
!
      type(REAL8List), intent(inout) :: list     ! list
      type(REAL8Node), pointer :: pos            ! position
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      list%pos => pos

   end subroutine setPosReal8List
!
!=======================================================================
!
   subroutine getposReal8List( list, pos )
!             ===============    i    o
!
!     Get current position
!
!  Arguments:
!
      type(REAL8List), intent(in) :: list        ! list
      type(REAL8Node), pointer :: pos            ! position
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      pos => list%pos

   end subroutine getposReal8List
!
!=======================================================================
!
   subroutine rewindReal8List( list )
!             ===============   io
!
!     Put current position at sentinel (nil)
!
!  Arguments:
!
      type(REAL8List), intent(inout) :: list     ! list
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      list%pos => list%nil

   end subroutine rewindReal8List
!
!=======================================================================
!
   subroutine windReal8List( list )
!             =============   io
!
!     Moves current position to the last node
!
!  Arguments:
!
      type(REAL8List), intent(inout) :: list     ! list
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      list%pos => list%nil%prev

   end subroutine windReal8List
!
!=======================================================================
!
   subroutine iterateReal8List( list, isAtEnd )
!             ================   io      o
!
!     Iterate current position and read value
!
!  Arguments:
!
      type(REAL8List), intent(inout) :: list           ! list
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

   end subroutine iterateReal8List
!
!=======================================================================
!  subroutines for managing the list contents
!=======================================================================
!
   subroutine getValReal8List( list, value, isAtBeg_ )
!             ===============   io     o       o
!
!     Return value at current position
!
!  Arguments:
!
      type(REAL8List), intent(inout) :: list    ! list
      real(REAL8), intent(out) :: value         ! value read
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
               error stop 'getValReal8List: current position is before first value'
      end if
      
      value = list%pos%value

   end subroutine getValReal8List
!
!=======================================================================
!
   subroutine insertReal8List( list, value )
!             ===============   io     i
!
!     Inserts value after current value and make new value current
!
!  Arguments:
!
      type(REAL8List), intent(inout) :: list     ! list 
      real(REAL8), intent(in) :: value        ! value to be inserted
!
!  Local Variables:
!
      type(REAL8Node), pointer :: p             ! temporary pointer
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

   end subroutine insertReal8List
!
!=======================================================================
!
   subroutine removeReal8List( list, isAtBeg_ )
!             ===============   io      o
!
!     Remove current element from the list and move current position before removed value
!
!  Arguments:
!
      type(REAL8List), intent(inout) :: list    ! list
      logical, intent(out), optional :: isAtBeg_! flag indicating current position is before
!                                               !   first value of list or list is empty
!
!  Local Variables:
!
      type(REAL8Node), pointer :: p
!
!  Start of Executable code:
!  ------------------------
!
      if( present( isAtBeg_ ) ) then
         isAtBeg_ = associated( list%pos, list%nil )
         if( isAtBeg_ ) return
      else
         if( associated( list%pos, list%nil ) )  &
               error stop 'removeReal8List: current position is before first value'
      end if
      
      p => list%pos
      p%next%prev => p%prev
      p%prev%next => p%next
      
      list%pos => p%prev
      
      deallocate( p )

      list%size = list%size - 1

   end subroutine removeReal8List
!
!=======================================================================
!
   subroutine clearReal8List( list )
!             ==============   io
!
!     Clears the whole list
!
!  Arguments:
!
      type(REAL8List), intent(inout) :: list
!
!  Local Variables:
!
      logical :: isAtEnd
!
!  Start of Executable code:
!  ------------------------
!
      call rewindReal8List( list )
!          ===============   io
      do
         call iterateReal8List( list, isAtEnd )
!             ================   io      o
         if( isAtEnd )   exit

         call removeReal8List( list )
!             ===============   io
      end do

   end subroutine clearReal8List
!
!=======================================================================
!  subroutines for getting information about the list
!=======================================================================
!
   integer(INT4) function sizeOfReal8List( list )
!                         ===============   i  
!
!     returns the size of the list, i.e. the number of items
!
!  Arguments:
!  
      type(REAL8List), intent(in) :: list      ! list
!
!  Start of Executable Code:
!  ------------------------
!
      sizeOfReal8List = list%size

   end function sizeOfReal8List
!
!=======================================================================
!
   subroutine dumpReal8List( list )
!             =============   io
!
!     Dumps the whole list to standard output
!
!  Arguments:
!
      type(REAL8List), intent(inout) :: list      ! list
!
!  Local Veriables:
!
      real(REAL8) :: value
      logical :: isAtEnd
      type(REAL8Node), pointer :: p
!
!  Start of Executable Code:
!  ------------------------
!
      p => list%pos                            ! remember current position

      call rewindReal8List( list )
!          ===============   io
      print *,  'List contents:'
      do
         call iterateReal8List( list, isAtEnd )
!             ================   io     o
         if( isAtEnd )   exit
!
         call getValReal8List( list, value )
!             ===============   io     o         
         print *,  value
      end do
      print *,  'End of list contents.'

      list%pos => p                            ! restore current position

   end subroutine dumpReal8List
!
!========================================================================================
!
end module List_REAL8_mod

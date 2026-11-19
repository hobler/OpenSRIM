!  To get working code for other list types replace
!
!     INT4    by    the name that should appear in the Node and List type names
!     INTEGER(INT4)    by    the Fortran type declaration statement
!
!  E.g., to define an Integer of kind INT4 list with node type Int4Node and 
!  list type Int4List, replace
!
!     INT4    by    Int4
!     INTEGER(INT4)    by    integer(INT4)
!
!  To define a list of derived data type MyType with node type MyTypeNode and 
!  list type MyTypeList, replace
!
!     INT4    by    MyType
!     INTEGER(INT4)    by    type(MyType)
!
module List_INT4_mod
!
!  Defines a doubly linked list of data type INT4, using a sentinel to mark 
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
   public :: newInt4List                   ! Initializes the list.
!            ===========
!
!  Managing current position
!
   public :: setPosInt4List                ! Sets current position
!            ==============
   public :: getPosInt4List                ! Gets current position
!            ==============      
   public :: rewindInt4List                ! Moves current position before first value.
!            ==============
   public :: windInt4List                  ! Moves current position to the last value.
!            ============
   public :: iterateInt4List               ! Moves current position forward.
!            ===============
!
!  Managing list contents
!
   public :: getvalInt4List                ! Gets value at current position.
!            ==============
   public :: insertInt4List                ! Inserts value after current position and
!            ==============                !   moves current position to new value.
   public :: removeInt4List                ! Removes value at current position and moves
!            ==============                !   current position before removed value.
   public :: clearInt4List                 ! Clears the whole list.
!            =============
!
!  Getting information
!
   public :: sizeOfInt4List                ! Returns the number of values in the list.
!            ==============
   public :: dumpInt4List                  ! Dumps the list to standard output.
!            ============
!
!  Public Data Types (with private components):
!  =================
!
   type, public :: INT4Node                     ! list node
!                  ========
      private
      integer(INT4) :: value
      type(INT4Node), pointer :: next
      type(INT4Node), pointer :: prev
   end type INT4Node

   type, public :: INT4List                     ! list for INT4 numbers
!                  ========
      private
      type(INT4Node), pointer :: nil => null( ) ! sentinel:
                                                !   nil%next => first node
                                                !   nil%prev => last node
                                                !   nil%value ... undefined
      type(INT4Node), pointer :: pos => null( ) ! pointer to current position
      integer(INT4) :: size = 0
   end type INT4List

contains
!
!=======================================================================
!  subroutine for initialization
!=======================================================================
!
   subroutine newInt4List( list )
!             ===========   io
!
!     Initializes the list. Must be called before the list can be used.
!
!  Arguments:
!
      type(INT4List), intent(inout) :: list     ! empty list
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( associated( list%nil ) )  error stop 'newInt4List: list already initialized'
      
      allocate( list%nil )
      list%nil%next => list%nil
      list%nil%prev => list%nil

      list%pos => list%nil

      list%size = 0
      
   end subroutine newInt4List 
!
!=======================================================================
!  subroutines for managing the current position
!=======================================================================
!
   subroutine setPosInt4List( list, pos )
!             ==============   io    i
!
!     Put current position at pos.
!
!     WARNING: Be careful and make sure that pos is a node of list. This cannot easily
!              be checked within this subroutine!
!
!  Arguments:
!
      type(INT4List), intent(inout) :: list     ! list
      type(INT4Node), pointer :: pos            ! position
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      list%pos => pos

   end subroutine setPosInt4List
!
!=======================================================================
!
   subroutine getPosInt4List( list, pos )
!             ==============    i    o
!
!     Get current position
!
!  Arguments:
!
      type(INT4List), intent(in) :: list        ! list
      type(INT4Node), pointer :: pos            ! position
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      pos => list%pos

   end subroutine getPosInt4List
!
!=======================================================================
!
   subroutine rewindInt4List( list )
!             ==============   io
!
!     Put current position at sentinel (nil)
!
!  Arguments:
!
      type(INT4List), intent(inout) :: list     ! list
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      list%pos => list%nil

   end subroutine rewindInt4List
!
!=======================================================================
!
   subroutine windInt4List( list )
!             ============   io
!
!     Moves current position to the last node
!
!  Arguments:
!
      type(INT4List), intent(inout) :: list     ! list
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      list%pos => list%nil%prev

   end subroutine windInt4List
!
!=======================================================================
!
   subroutine iterateInt4List( list, isAtEnd )
!             ===============   io      o
!
!     Iterate current position and read value
!
!  Arguments:
!
      type(INT4List), intent(inout) :: list           ! list
      logical, intent(out) :: isAtEnd                 ! flag indicating current position is at
!                                                     !   end of list or list is empty
!
!  Local Variables:
!
!    None
!
!  Start of Executable Code:
!  ------------------------
!
      isAtEnd = associated( list%pos%next, list%nil )
      if( isAtEnd ) return

      list%pos => list%pos%next

   end subroutine iterateInt4List
!
!=======================================================================
!  subroutines for managing the list contents
!=======================================================================
!
   subroutine getValInt4List( list, value, isAtBeg_ )
!             ==============   io     o       o
!
!     Return value at current position
!
!  Arguments:
!
      type(INT4List), intent(inout) :: list     ! list
      integer(INT4), intent(out) :: value       ! value read
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
               error stop 'getValInt4List: current position is before first value'
      end if
      
      value = list%pos%value

   end subroutine getValInt4List
!
!=======================================================================
!
   subroutine insertInt4List( list, value )
!             ==============   io     i
!
!     Inserts value after current value and make new value current
!
!  Arguments:
!
      type(INT4List), intent(inout) :: list     ! list 
      integer(INT4), intent(in) :: value        ! value to be inserted
!
!  Local Variables:
!
      type(INT4Node), pointer :: p             ! temporary pointer
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

   end subroutine insertInt4List
!
!=======================================================================
!
   subroutine removeInt4List( list, isAtBeg_ )
!             ==============   io      o
!
!     remove current element from the list and move current position before removed value
!
!  Arguments:
!
      type(INT4List), intent(inout) :: list     ! list
      logical, intent(out), optional :: isAtBeg_! flag indicating current position is before
!                                               !   first value of list or list is empty
!
!  Local Variables:
!
      type(INT4Node), pointer :: p
!
!  Start of Executable code:
!  ------------------------
!
      if( present( isAtBeg_ ) ) then
         isAtBeg_ = associated( list%pos, list%nil )
         if( isAtBeg_ ) return
      else
         if( associated( list%pos, list%nil ) )  &
               error stop 'removeInt4List: current position is before first value'
      end if
      
      p => list%pos
      p%next%prev => p%prev
      p%prev%next => p%next
      
      list%pos => p%prev
      
      deallocate( p )

      list%size = list%size - 1

   end subroutine removeInt4List
!
!=======================================================================
!
   subroutine clearInt4List( list )
!             =============   io
!
!     clears the whole list
!
!  Arguments:
!
      type(INT4List), intent(inout) :: list
!
!  Local Variables:
!
      logical :: isAtEnd
!
!  Start of Executable code:
!  ------------------------
!
      call rewindInt4List( list )
!          ==============   io
      do
         call iterateInt4List( list, isAtEnd )
!             ===============   io      o
         if( isAtEnd )  exit
!
         call removeInt4List( list )
!             ==============   io
      end do

   end subroutine clearInt4List
!
!=======================================================================
!  subroutines for getting information about the list
!=======================================================================
!
   integer(INT4) function sizeOfInt4List( list )
!                         ==============   i  
!
!     returns the size of the list, i.e. the number of items
!
!  Arguments:
!  
      type(INT4List), intent(in) :: list      ! list
!
!  Start of Executable Code:
!  ------------------------
!
      sizeOfInt4List = list%size

   end function sizeOfInt4List
!
!=======================================================================
!
   subroutine dumpInt4List( list )
!             ============   io
!
!     Dumps the whole list to standard output
!
!  Arguments:
!
      type(INT4List), intent(inout) :: list      ! list
!
!  Local Veriables:
!
      integer(INT4) :: value
      logical :: isAtEnd
      type(INT4Node), pointer :: p
!
!  Start of Executable Code:
!  ------------------------
!
      p => list%pos                            ! remember current position

      call rewindInt4List( list )
!          ==============   io
      print *,  'List contents:'
      do
         call iterateInt4List( list, isAtEnd )
!             ===============   io     o
         if( isAtEnd )  exit

         call getValInt4List( list, value )
!             ==============   io     o         
         print *,  value
      end do
      print *,  'End of list contents.'

      list%pos => p                            ! restore current position

   end subroutine dumpInt4List
!
!========================================================================================
!
end module List_INT4_mod

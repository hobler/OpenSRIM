module Stack_Int4_mod
!
!  Defines a stack of Int4 data type
!
!  Used Modules:
!  ============
!
   use Const_mod, only: INT4, REAL8
   use Array_mod, only: check_allocation
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
   public :: newInt4Stack                    ! Creates a new empty stack.
!            ============
   public :: sizeOfInt4Stack                 ! Returns the number of items.
!            ===============
   public :: pushInt4Stack                   ! Pushes a value to the stack.
!            =============
   public :: popInt4Stack                    ! Pops a value from the stack.
!            ============
   public :: dumpInt4Stack                   ! Dumps the stack to standard output.
!            =============
   public :: clearInt4Stack                  ! Clears the whole stack.
!            ==============
   public :: int4Stack2Array                 ! Convert stack to array, clearing the stack
!            ===============
!
!  Public Data Type:
!  ================
!
   type, private :: Int4Node                 ! node for Int4 stack
!                   ========
      integer(INT4) :: value
      type(Int4Node), pointer :: next
   end type Int4Node

   type, public :: Int4Stack                 ! stack for Int4 numbers
!                  =========
      private
      type(Int4Node), pointer :: head
      integer(INT4) :: size
   end type Int4Stack
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine newInt4Stack( stack )
!             ============    o
!
!     Creates a new empty stack    
!
!  Arguments:
!
      type(Int4Stack), intent(out) :: stack     ! new empty stack
!
!  Start of Executable Code:
!  ------------------------
!
      nullify( stack%head )

      stack%size = 0

   end subroutine newInt4Stack
!
!=======================================================================
!
   integer(INT4) function sizeOfInt4Stack( stack )
!                         ===============    i
!
!     Returns the size of the stack, i.e. the number of items
!
!  Arguments:
!  
      type(Int4Stack), intent(in) :: stack      ! stack
!
!  Start of Executable Code:
!  ------------------------
!
      sizeOfInt4Stack = stack%size

   end function sizeOfInt4Stack
!
!=======================================================================
!
   subroutine pushInt4Stack( value, stack )
!             =============    i     io
!
!     Pushes a value to the top of the stack
!
!  Arguments:
!
      integer(INT4), intent(in) :: value        ! value to be pushed
      type(Int4Stack), intent(inout) :: stack   ! stack 
!
!  Local Variables:
!
      type(Int4Node), pointer :: p             ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      p => stack%head
      allocate( stack%head )
      stack%head%value = value
      stack%head%next => p

      stack%size = stack%size + 1

   end subroutine pushInt4Stack
!
!=======================================================================
!
   subroutine popInt4Stack( stack, value, isEmpty )
!             ============   io      o       o
!  
!     Pops value from the stack
!
!  Arguments:
!
      integer(INT4), intent(out) :: value       ! value popped
      type(Int4Stack), intent(inout) :: stack   ! stack
      logical, intent(out) :: isEmpty           ! flag indicating empty stack
!
!  Local Variables:
!
      type(Int4Node), pointer :: p              ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      isEmpty = .not. associated( stack%head )
      if( isEmpty ) return

      value = stack%head%value

      p => stack%head%next
      deallocate( stack%head )
      stack%head => p

      stack%size = stack%size - 1

   end subroutine popInt4Stack
!
!=======================================================================
!
   subroutine dumpInt4Stack( stack )
!             =============    i
!  
!     Dumps the whole stack to standard output
!
!  Arguments:
!
      type(Int4Stack), intent(in) :: stack      ! stack
!
!  Local Veriables:
!
      type(Int4Node), pointer :: p              ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      p => stack%head
      print *,  'Contents of stack:'
      do
         if( .not. associated( p ) )  exit
         print *,  p%value
         p => p%next
      enddo
      print *,  'End of contents of stack'

   end subroutine dumpInt4Stack
!
!=======================================================================
!
   subroutine clearInt4Stack( stack )
!             ==============   io
!  
!     Clears the whole stack
!
!  Arguments:
!
      type(Int4Stack), intent(inout) :: stack
!
!  Local Variables:
!
      integer(INT4) :: value
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      do
         call popInt4Stack( stack, value, isEmpty)
!             ============    io     o       o
         if( isEmpty )  exit
      enddo

   end subroutine clearInt4Stack
!
!=======================================================================
!
   subroutine int4Stack2Array( stack, array )
!             ===============   io      o
!
!     Convert stack to array, clearing the stack
!
!  Arguments:
!
      type(Int4Stack), intent(inout) :: stack
      integer(INT4), allocatable :: array(:)
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: i, n, istat, value
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      n = sizeOfInt4Stack( stack )
!         ===============    i
      allocate( array( n ), stat=istat, errmsg=errmsg )
      call check_allocation( 'int4Stack2Array', 'array', istat, errmsg )

      do i = n, 1, -1
         call popInt4Stack( stack, value, isEmpty)
!             ============   io      o      o
         array(i) = value
      enddo

   end subroutine int4Stack2Array
!
!=======================================================================
!
end module Stack_Int4_mod

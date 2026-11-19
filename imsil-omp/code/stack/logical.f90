module Stack_Logical_mod
!
!  Defines a stack of logical data type
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
!
!  Public Procedures:
!  =================
!
   public :: newLogicalStack                 ! Creates a new empty stack.
!            ===============
   public :: sizeOfLogicalStack              ! Returns the number of items.  
!            ==================
   public :: pushLogicalStack                ! Pushes a value to the stack.
!            ================
   public :: popLogicalStack                 ! Pops a value from the stack.
!            ===============
   public :: dumpLogicalStack                ! Dumps the stack to standard output
!            ================
   public :: clearLogicalStack               ! Clears the whole stack.
!            =================
   public :: logicalStack2Array              ! Convert stack to array, clearing the stack
!            ==================
!
!  Public Data Type:
!  ================
!
   type, private :: LogicalNode              ! node for Logical stack
!                   ===========
      logical :: value
      type(LogicalNode), pointer :: next
   end type LogicalNode

   type, public :: LogicalStack              ! stack for logical values
!                  ============
      private
      type(LogicalNode), pointer :: head
      integer(INT4) :: size
   end type LogicalStack
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine newLogicalStack( stack )
!             ===============    o
!
!     Creates a new empty stack    
!
!  Arguments:
!
      type(LogicalStack), intent(out) :: stack  ! new empty stack
!
!  Start of Executable Code:
!  ------------------------
!
      nullify( stack%head )

      stack%size = 0

   end subroutine newLogicalStack
!
!=======================================================================
!
   integer(INT4) function sizeOfLogicalStack( stack )
!                         ==================    i
!  
!     Returns the size of the stack, i.e. the number of items
!
!  Arguments:
!  
      type(LogicalStack), intent(in) :: stack   ! stack
!
!  Start of Executable Code:
!  ------------------------
!
      sizeOfLogicalStack = stack%size

   end function sizeOfLogicalStack
!
!=======================================================================
!
   subroutine pushLogicalStack( value, stack )
!             ================    i     io
!  
!     Pushes a value to the top of the stack
!
!  Arguments:
!
      logical, intent(in) :: value                 ! value to be pushed
      type(LogicalStack), intent(inout) :: stack   ! stack 
!
!  Local Variables:
!
      type(LogicalNode), pointer :: p           ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      p => stack%head
      allocate( stack%head )
      stack%head%value = value
      stack%head%next => p

      stack%size = stack%size + 1

   end subroutine pushLogicalStack
!
!=======================================================================
!
   subroutine popLogicalStack( stack, value, isEmpty )
!             ===============    io     o       o
!  
!     Pops value from the stack
!
!  Arguments:
!
      logical, intent(out) :: value                ! value popped
      type(LogicalStack), intent(inout) :: stack   ! stack
      logical, intent(out) :: isEmpty              ! flag indicating empty stack
!
!  Local Variables:
!
      type(LogicalNode), pointer :: p              ! temporary pointer
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

   end subroutine popLogicalStack
!
!=======================================================================
!
   subroutine dumpLogicalStack( stack )
!             ================    i
!  
!     Dumps the whole stack to standard output
!
!  Arguments:
!
      type(LogicalStack), intent(in) :: stack   ! stack
!
!  Local Veriables:
!
      type(LogicalNode), pointer :: p           ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      p => stack%head
      print *, 'Contents of stack:'
      do
         if( .not. associated( p ) )  exit
         print *,  p%value
         p => p%next
      enddo
      print *,  'End of contents of stack'

   end subroutine dumpLogicalStack
!
!=======================================================================
!
   subroutine clearLogicalStack( stack )
!             =================    io
!  
!     Clears the whole stack
!
!  Arguments:
!
      type(LogicalStack), intent(inout) :: stack
!
!  Local Variables:
!
      logical :: value
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      do
         call popLogicalStack( stack, value, isEmpty )
!             ===============    io     o       o
         if( isEmpty )  exit
      enddo

   end subroutine clearLogicalStack
!
!=======================================================================
!
   subroutine logicalStack2Array( stack, array )
!             ==================   io      o
!
!     Convert stack to array, clearing the stack
!
!  Arguments:
!
      type(LogicalStack), intent(inout) :: stack
      logical, allocatable :: array(:)
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: i, n, istat
      logical :: value
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      n = sizeOfLogicalStack( stack )
!         ==================    i
      allocate( array( n ), stat=istat, errmsg=errmsg )
      call check_allocation( 'logicalStack2Array', 'array', istat, errmsg )

      do i = n, 1, -1
         call popLogicalStack( stack, value, isEmpty)
!             ===============   io      o      o
         array(i) = value
      enddo

   end subroutine logicalStack2Array
!
!=======================================================================
!
end module Stack_Logical_mod

module Stack_Real8_mod
!
!  Defines a stack of Real8 data type
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
   public :: newReal8Stack                   ! Creates a new empty stack.
!            =============
   public :: sizeOfReal8Stack                ! Returns the number of items.
!            ================
   public :: pushReal8Stack                  ! Pushes a value to the stack.
!            ==============
   public :: popReal8Stack                   ! Pops a value from the stack.
!            =============
   public :: dumpReal8Stack                  ! Dumps the stack to standard output.
!            ==============
   public :: clearReal8Stack                 ! Clears the whole stack.
!            ===============
   public :: real8Stack2Array                ! Convert stack to array, clearing the stack
!            ================
!
!  Public Data Types:
!  =================
!
   type, private :: Real8Node                ! node for real8 stack
!                   =========
      real(REAL8) :: value
      type(Real8Node), pointer :: next
   end type Real8Node

   type, public :: Real8Stack                ! stack for real8 numbers
!                  ==========
      private
      type(Real8Node), pointer :: head
      integer(INT4) :: size
   end type Real8Stack
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine newReal8Stack( stack )
!             =============    o
!
!     Creates a new empty stack    
!
!  Arguments:
!
      type(Real8Stack), intent(out) :: stack    ! new empty stack
!
!  Start of Executable Code:
!  ------------------------
!
      nullify( stack%head )

      stack%size = 0

   end subroutine newReal8Stack
!
!=======================================================================
!
   integer(INT4) function sizeOfReal8Stack( stack )
!                         ================    i
!  
!     Returns the size of the stack, i.e. the number of items
!
!  Arguments:
!  
      type(Real8Stack), intent(in) :: stack     ! stack
!
!  Start of Executable Code:
!  ------------------------
!
      sizeOfReal8Stack = stack%size

   end function sizeOfReal8Stack
!
!=======================================================================
!
   subroutine pushReal8Stack( value, stack )
!             ==============    i     io
!  
!     pushes a value to the top of the stack
!
!  Arguments:
!
      real(REAL8), intent(in) :: value          ! value to be pushed
      type(Real8Stack), intent(inout) :: stack  ! stack 
!
!  Local Variables:
!
      type(Real8Node), pointer :: p             ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      p => stack%head
      allocate( stack%head )
      stack%head%value = value
      stack%head%next => p

      stack%size = stack%size + 1

   end subroutine pushReal8Stack
!
!=======================================================================
!
   subroutine popReal8Stack( stack, value, isEmpty )
!             =============   io      o       o
!  
!     Pops value from the stack
!
!  Arguments:
!
      real(REAL8), intent(out) :: value         ! value popped
      type(Real8Stack), intent(inout) :: stack  ! stack
      logical, intent(out) :: isEmpty           ! flag indicating empty stack
!
!  Local Variables:
!
      type(Real8Node), pointer :: p             ! temporary pointer
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

   end subroutine popReal8Stack
!
!=======================================================================
!
   subroutine dumpReal8Stack( stack )
!             ==============    i
!  
!     Dumps the whole stack to standard output
!
!  Arguments:
!
      type(Real8Stack), intent(in) :: stack     ! stack
!
!  Local Veriables:
!
      type(Real8Node), pointer :: p             ! temporary pointer
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
      print *, 'End of contents of stack'

   end subroutine dumpReal8Stack
!
!=======================================================================
!
   subroutine clearReal8Stack( stack )
!             ===============   io
!  
!     Clears the whole stack
!
!  Arguments:
!
      type(Real8Stack), intent(inout) :: stack
!
!  Local Variables:
!
      real(REAL8) :: value
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      do
         call popReal8Stack( stack, value, isEmpty )
!             =============   io      o       o
         if( isEmpty )  exit
      enddo

   end subroutine clearReal8Stack
!
!=======================================================================
!
   subroutine real8Stack2Array( stack, array )
!             ================   io      o
!
!     Convert stack to array, clearing the stack
!
!  Arguments:
!
      type(Real8Stack), intent(inout) :: stack
      real(REAL8), allocatable :: array(:)
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: i, n, istat
      real(REAL8) :: value
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      n = sizeOfReal8Stack( stack )
!         ================    i
      allocate( array( n ), stat=istat, errmsg=errmsg )
      call check_allocation( 'real8Stack2Array', 'array', istat, errmsg )

      do i = n, 1, -1
         call popReal8Stack( stack, value, isEmpty )
!             =============   io      o      o
         array(i) = value
      enddo

   end subroutine real8Stack2Array
!
!=======================================================================
!
end module Stack_Real8_mod

module Stack_Vec2_mod
!
!  Defines a stack of Vec2 data type
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
   public :: newVec2Stack                    ! Creates a new empty stack.
!            ============
   public :: sizeOfVec2Stack                 ! Returns the number of items.
!            ===============
   public :: pushVec2Stack                   ! Pushes a value to the stack.
!            =============
   public :: popVec2Stack                    ! Pops a value from the stack.
!            ============
   public :: dumpVec2Stack                   ! Dumps the stack to standard output.
!            =============
   public :: clearVec2Stack                  ! Clears the whole stack.
!            ==============
   public :: vec2Stack2Array                 ! Convert stack to array, clearing the stack
!            ===============
!
!  Public Data Types:
!  =================
!
   type, private :: Vec2Node                 ! node for Vec2 stack
!                   ========
      real(REAL8) :: value(2)
      type(Vec2Node), pointer :: next
   end type Vec2Node

   type, public :: Vec2Stack                 ! stack of Vec2s
!                  =========
      private
      type(Vec2Node), pointer :: head
      integer(INT4) :: size
   end type Vec2Stack
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine newVec2Stack( stack )
!             ============    o
!
!     Creates a new empty stack    
!
!  Arguments:
!
      type(Vec2Stack), intent(out) :: stack    ! new empty stack
!
!  Start of Executable Code:
!  ------------------------
!
      nullify( stack%head )

      stack%size = 0

   end subroutine newVec2Stack
!
!=======================================================================
!
   integer(INT4) function sizeOfVec2Stack( stack )
!                         ===============    i
!  
!     Returns the size of the stack, i.e. the number of items
!
!  Arguments:
!  
      type(Vec2Stack), intent(in) :: stack      ! stack
!
!  Start of Executable Code:
!  ------------------------
!
      sizeOfVec2Stack = stack%size

   end function sizeOfVec2Stack
!
!=======================================================================
!
   subroutine pushVec2Stack( value, stack )
!             =============    i     io
!  
!     Pushes a value to the top of the stack
!
!  Arguments:
!
      real(REAL8), intent(in) :: value(2)       ! value to be pushed
      type(Vec2Stack), intent(inout) :: stack   ! stack 
!
!  Local Variables:
!
      type(Vec2Node), pointer :: p              ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      p => stack%head
      allocate( stack%head )
      stack%head%value(:) = value(:)
      stack%head%next => p

      stack%size = stack%size + 1

   end subroutine pushVec2Stack
!
!=======================================================================
!
   subroutine popVec2Stack( stack, value, isEmpty )
!             ============   io      o       o
!  
!     Pops value from the stack
!
!  Arguments:
!
      real(REAL8), intent(out) :: value(2)       ! value popped
      type(Vec2Stack), intent(inout) :: stack  ! stack
      logical, intent(out) :: isEmpty            ! flag indicating empty stack
!
!  Local Variables:
!
      type(Vec2Node), pointer :: p             ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      isEmpty = .not. associated( stack%head )
      if( isEmpty ) return

      value(:) = stack%head%value(:)

      p => stack%head%next
      deallocate( stack%head )
      stack%head => p

      stack%size = stack%size - 1

   end subroutine popVec2Stack
!
!=======================================================================
!
   subroutine dumpVec2Stack( stack )
!             =============    i
!  
!     Dumps the whole stack to standard output
!
!  Arguments:
!
      type(Vec2Stack), intent(in) :: stack     ! stack
!
!  Local Variables:
!
      type(Vec2Node), pointer :: p             ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      p => stack%head
      print *,  'Contents of stack:'
      do
         if( .not. associated( p ) )  exit
         print *,  p%value(:)
         p => p%next
      enddo
      print *,  'End of contents of stack'

   end subroutine dumpVec2Stack
!
!=======================================================================
!
   subroutine clearVec2Stack( stack )
!             ==============   io
!  
!     Clears the whole stack
!
!  Arguments:
!
      type(Vec2Stack), intent(inout) :: stack
!
!  Local Variables:
!
      real(REAL8) :: value(2)
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      do
         call popVec2Stack( stack, value, isEmpty )
!             ============    io     o       o
         if( isEmpty ) exit
      enddo

   end subroutine clearVec2Stack
!
!=======================================================================
!
   subroutine vec2Stack2Array( stack, array )
!             ===============   io      o
!
!     Convert stack to array, clearing the stack
!
!  Arguments:
!
      type(Vec2Stack), intent(inout) :: stack
      real(REAL8), allocatable :: array(:,:)
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: i, n, istat
      real(REAL8) :: value(2)
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      n = sizeOfVec2Stack( stack )
!         ===============    i
      allocate( array( 2, n ), stat=istat, errmsg=errmsg )
      call check_allocation( 'vec2StackArray', 'array', istat, errmsg )

      do i = n, 1, -1
         call popVec2Stack( stack, value, isEmpty)
!             ============   io      o      o
         array(:,i) = value(:)
      enddo

   end subroutine vec2Stack2Array
!
!=======================================================================
!
end module Stack_Vec2_mod

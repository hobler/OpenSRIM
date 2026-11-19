module Stack_Proj_mod
!
!  Defines a stack of Proj data type
!
!  Used Modules:
!  ============
!
   use Const_mod, only: INT4, REAL8
   use Array_mod, only: check_allocation
   use Types_mod, only: Proj_type
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
   public :: newProjStack                    ! Creates a new empty stack.
!            ============
   public :: sizeOfProjStack                 ! Returns the number of items.
!            ===============
   public :: pushProjStack                   ! Pushes a value to the stack.
!            =============
   public :: popProjStack                    ! Pops a value from the stack.
!            ============
   public :: dumpProjStack                   ! Dumps the stack to standard output.
!            =============
   public :: clearProjStack                  ! Clears the whole stack.
!            ==============
   public :: ProjStack2Array                 ! Convert stack to array, clearing the stack
!            ===============
!
!  Public Data Types:
!  =================
!
   type, private :: ProjNode                 ! node for Proj stack
!                   ========
      type(Proj_type) :: value
      type(ProjNode), pointer :: next
   end type ProjNode

   type, public :: ProjStack                 ! stack of Proj_type
!                  =========
      private
      type(ProjNode), pointer :: head
      integer(INT4) :: size
   end type ProjStack
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine newProjStack( stack )
!             ============    o
!
!     Creates a new empty stack    
!
!  Arguments:
!
      type(ProjStack), intent(out) :: stack     ! new empty stack
!
!  Start of Executable Code:
!  ------------------------
!
      nullify( stack%head )

      stack%size = 0

   end subroutine newProjStack
!
!=======================================================================
!
   integer(INT4) function sizeOfProjStack( stack )
!                         ===============    i
!
!     Returns the size of the stack, i.e. the number of items
!
!  Arguments:
!  
      type(ProjStack), intent(in) :: stack      ! stack
!
!  Start of Executable Code:
!  ------------------------
!
      sizeOfProjStack = stack%size

   end function sizeOfProjStack
!
!=======================================================================
!
   subroutine pushProjStack( value, stack )
!             =============    i     io
!
!     Pushes a value to the top of the stack
!
!  Arguments:
!
      type(Proj_type), intent(in) :: value      ! value to be pushed
      type(ProjStack), intent(inout) :: stack   ! stack
!
!  Local Variables:
!
      type(ProjNode), pointer :: p             ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      p => stack%head
      allocate( stack%head )
      stack%head%value = value
      stack%head%next => p

      stack%size = stack%size + 1

   end subroutine pushProjStack
!
!=======================================================================
!
   subroutine popProjStack( stack, value, isEmpty )
!             ============   io      o       o
!  
!     Pops value from the stack
!
!  Arguments:
!
      type(Proj_type), intent(out) :: value       ! value popped
      type(ProjStack), intent(inout) :: stack   ! stack
      logical, intent(out) :: isEmpty           ! flag indicating empty stack
!
!  Local Variables:
!
      type(ProjNode), pointer :: p              ! temporary pointer
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

   end subroutine popProjStack
!
!=======================================================================
!
   subroutine dumpProjStack( stack )
!             =============    i
!  
!     Dumps the whole stack to standard output
!
!  Arguments:
!
      type(ProjStack), intent(in) :: stack      ! stack
!
!  Local Veriables:
!
      type(ProjNode), pointer :: p              ! temporary pointer
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

   end subroutine dumpProjStack
!
!=======================================================================
!
   subroutine clearProjStack( stack )
!             ==============   io
!  
!     Clears the whole stack
!
!  Arguments:
!
      type(ProjStack), intent(inout) :: stack
!
!  Local Variables:
!
      type(Proj_type) :: value
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      do
         call popProjStack( stack, value, isEmpty)
!             ============    io     o       o
         if( isEmpty )  exit
      enddo

   end subroutine clearProjStack
!
!=======================================================================
!
   subroutine ProjStack2Array( stack, array )
!             ===============   io      o
!
!     Convert stack to array, clearing the stack
!
!  Arguments:
!
      type(ProjStack), intent(inout) :: stack
      type(Proj_type), allocatable :: array(:)
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: i, n, istat
      type(Proj_type) :: value
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      n = sizeOfProjStack( stack )
!         ===============    i
      allocate( array( n ), stat=istat, errmsg=errmsg )
      call check_allocation( 'ProjStack2Array', 'array', istat, errmsg )

      do i = n, 1, -1
         call popProjStack( stack, value, isEmpty)
!             ============   io      o      o
         array(i) = value
      enddo

   end subroutine ProjStack2Array
!
!=======================================================================
!
end module Stack_Proj_mod

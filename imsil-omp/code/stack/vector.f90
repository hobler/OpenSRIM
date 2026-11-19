module Stack_Vector_mod
!
!  Defines a stack of Vector data type
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
   public :: newVectorStack                  ! Creates a new empty stack.
!            ==============
   public :: sizeOfVectorStack               ! Returns the number of items.
!            =================
   public :: pushVectorStack                 ! Pushes a value to the stack.
!            ===============
   public :: popVectorStack                  ! Pops a value from the stack.
!            ==============
   public :: dumpVectorStack                 ! Dumps the stack to standard output.
!            ===============
   public :: clearVectorStack                ! Clears the whole stack.
!            ================
   public :: vectorStack2Array               ! Convert stack to array, clearing the stack
!            =================
!
!  Public Data Types:
!  =================
!
   type, private :: VectorNode               ! node for vector stack
!                   ==========
      real(REAL8) :: value(3)
      type(VectorNode), pointer :: next
   end type VectorNode

   type, public :: VectorStack               ! stack of Vectors
!                  ===========
      private
      type(VectorNode), pointer :: head
      integer(INT4) :: size
   end type VectorStack
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine newVectorStack( stack )
!             ==============    o
!
!     Creates a new empty stack    
!
!  Arguments:
!
      type(VectorStack), intent(out) :: stack    ! new empty stack
!
!  Start of Executable Code:
!  ------------------------
!
      nullify( stack%head )

      stack%size = 0

   end subroutine newVectorStack
!
!=======================================================================
!
   integer(INT4) function sizeOfVectorStack( stack )
!                         =================    i
!  
!     Returns the size of the stack, i.e. the number of items
!
!  Arguments:
!  
      type(VectorStack), intent(in) :: stack    ! stack
!
!  Start of Executable Code:
!  ------------------------
!
      sizeOfVectorStack = stack%size

   end function sizeOfVectorStack
!
!=======================================================================
!
   subroutine pushVectorStack( value, stack )
!             ===============    i     io
!  
!     Pushes a value to the top of the stack
!
!  Arguments:
!
      real(REAL8), intent(in) :: value(3)        ! value to be pushed
      type(VectorStack), intent(inout) :: stack  ! stack 
!
!  Local Variables:
!
      type(VectorNode), pointer :: p               ! temporary pointer
!
!  Start of Executable Code:
!  ------------------------
!
      p => stack%head
      allocate( stack%head )
      stack%head%value(:) = value(:)
      stack%head%next => p

      stack%size = stack%size + 1

   end subroutine pushVectorStack
!
!=======================================================================
!
   subroutine popVectorStack( stack, value, isEmpty )
!             ==============   io      o       o
!  
!     Pops value from the stack
!
!  Arguments:
!
      real(REAL8), intent(out) :: value(3)       ! value popped
      type(VectorStack), intent(inout) :: stack  ! stack
      logical, intent(out) :: isEmpty            ! flag indicating empty stack
!
!  Local Variables:
!
      type(VectorNode), pointer :: p             ! temporary pointer
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

   end subroutine popVectorStack
!
!=======================================================================
!
   subroutine dumpVectorStack( stack )
!             ===============    i
!  
!     Dumps the whole stack to standard output
!
!  Arguments:
!
      type(VectorStack), intent(in) :: stack     ! stack
!
!  Local Variables:
!
      type(VectorNode), pointer :: p             ! temporary pointer
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

   end subroutine dumpVectorStack
!
!=======================================================================
!
   subroutine clearVectorStack( stack )
!             ================   io
!  
!     Clears the whole stack
!
!  Arguments:
!
      type(VectorStack), intent(inout) :: stack
!
!  Local Variables:
!
      real(REAL8) :: value(3)
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      do
         call popVectorStack( stack, value, isEmpty )
!             ==============    io     o       o
         if( isEmpty ) exit
      enddo

   end subroutine clearVectorStack
!
!=======================================================================
!
   subroutine vectorStack2Array( stack, array )
!             =================   io      o
!
!     Convert stack to array, clearing the stack
!
!  Arguments:
!
      type(VectorStack), intent(inout) :: stack
      real(REAL8), allocatable :: array(:,:)
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: i, n, istat
      real(REAL8) :: value(3)
      logical :: isEmpty
!
!  Start of Executable code:
!  ------------------------
!
      n = sizeOfVectorStack( stack )
!         =================    i
      allocate( array( 3, n ), stat=istat, errmsg=errmsg )
      call check_allocation( 'vectorStack2Array', 'array', istat, errmsg )

      do i = n, 1, -1
         call popVectorStack( stack, value, isEmpty)
!             ==============   io      o      o
         array(:,i) = value(:)
      enddo

   end subroutine vectorStack2Array
!
!=======================================================================
!
end module Stack_Vector_mod

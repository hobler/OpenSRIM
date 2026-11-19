module Stack_mod
!
!  Defines a generic stack data type ("first in, last out")
!
!  Used Modules:
!  ============
!
   use Stack_Int4_mod
   use Stack_Logical_mod
   use Stack_Real8_mod
   use Stack_Vec2_mod
   use Stack_Vector_mod
   use Stack_Proj_mod
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
   public :: newStack                        ! Creates a new empty stack.
!            ========
   interface newStack
      module procedure newInt4Stack          ! stack of Int4
      module procedure newLogicalStack       ! stack of logical
      module procedure newReal8Stack         ! stack of Real8
      module procedure newVec2Stack          ! stack of Vec2
      module procedure newVectorStack        ! stack of Vector
      module procedure newProjStack          ! stack of Proj_type
   end interface

   public :: sizeOfStack                     ! Returns the number of items.
!            ===========
   interface sizeOfStack
      module procedure sizeOfInt4Stack       ! stack of Int4
      module procedure sizeOfLogicalStack    ! stack of logical
      module procedure sizeOfReal8Stack      ! stack of Real8
      module procedure sizeOfVec2Stack       ! stack of Vec2
      module procedure sizeOfVectorStack     ! stack of Vector
      module procedure sizeOfProjStack       ! stack of Proj_type
   end interface

   public :: pushStack                       ! Pushes a value to the stack.
!            =========
   interface pushStack
      module procedure pushInt4Stack         ! stack of Int4
      module procedure pushLogicalStack      ! stack of logical
      module procedure pushReal8Stack        ! stack of Real8
      module procedure pushVec2Stack         ! stack of Vec2
      module procedure pushVectorStack       ! stack of Vector
      module procedure pushProjStack         ! stack of Proj_type
   end interface

   public :: popStack                        ! Pops a value from the stack.
!            ========
   interface popStack
      module procedure popInt4Stack          ! stack of Int4
      module procedure popLogicalStack       ! stack of logical
      module procedure popReal8Stack         ! stack of Real8
      module procedure popVec2Stack          ! stack of Vec2
      module procedure popVectorStack        ! stack of Vector
      module procedure popProjStack          ! stack of Proj_type
   end interface

   public :: dumpStack                       ! Dumps the stack to standard output.
!            =========
   interface dumpStack
      module procedure dumpInt4Stack         ! stack of Int4
      module procedure dumpLogicalStack      ! stack of logical
      module procedure dumpReal8Stack        ! stack of Real8
      module procedure dumpVec2Stack         ! stack of Vec2
      module procedure dumpVectorStack       ! stack of Vector
      module procedure dumpProjStack         ! stack of Proj_type
   end interface

   public :: clearStack                      ! Clears the whole stack.
!            ==========
   interface clearStack
      module procedure clearInt4Stack        ! stack of Int4
      module procedure clearLogicalStack     ! stack of logical
      module procedure clearReal8Stack       ! stack of Real8
      module procedure clearVec2Stack        ! stack of Vec2
      module procedure clearVectorStack      ! stack of Vector
      module procedure clearProjStack        ! stack of Proj_type
   end interface

   public :: stack2Array                     ! Convert stack to array, clearing the stack
!            ===========
   interface stack2Array
      module procedure int4Stack2Array       ! stack of Int4
      module procedure logicalStack2Array    ! stack of logical
      module procedure real8Stack2Array      ! stack of Real8
      module procedure vec2Stack2Array       ! stack of Vec2
      module procedure vectorStack2Array     ! stack of Vector
      module procedure ProjStack2Array       ! stack of Proj_type
   end interface
!
!  Public Data Types:
!  =================
!
   public :: Int4Stack
   public :: LogicalStack
   public :: Real8Stack
   public :: Vec2Stack
   public :: VectorStack
   public :: ProjStack
!
!=======================================================================
!
end module Stack_mod

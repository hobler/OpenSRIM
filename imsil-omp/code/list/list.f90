module List_mod
!
!  Defines a generic list data type (add lists of other types as they become available)
!
!  Used Modules:
!  ============
!
   use List_Int4_mod
   use List_Vec2_mod
   use List_Real8_mod
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
   public :: newList                         ! Initializes the list.
!            =======
!
!  Managing current position
!
   public :: setPosList                      ! Sets current position.
!            ==========      
   public :: getPosList                      ! Gets current position.
!            ==========      
   public :: rewindList                      ! Moves current position before first value.
!            ==========
   public :: windList                        ! Moves current position to last value.
!            ========
   public :: iterateList                     ! Moves current position forward.
!            ===========
!
!  Managing list contents
!
   public :: getvalList                      ! Gets value at current position.
!            ==========
   public :: insertList                      ! Inserts value after current position and
!            ==========                      !   moves current position to new value.
   public :: removeList                      ! Remove value at current position and moves 
!            ==========                      !   current position before removed value.
   public :: clearList                       ! Clears the whole list.
!            =========
!
!  Getting information
!
   public :: sizeOfList                      ! Returns the number of values it the list.
!            ==========
   public :: dumpList                        ! Dumps the list to standard output.
!            ========
!
!  Generic procedure definition
!
   interface newList
      module procedure newInt4List         ! list of Int4
      module procedure newVec2List         ! list of Vec2
      module procedure newReal8List        ! list of Real8
   end interface

   interface setPosList
      module procedure setPosInt4List      ! list of Int4
      module procedure setPosVec2List      ! list of Vec2
      module procedure setPosReal8List     ! list of Real8
   end interface

   interface getPosList
      module procedure getPosInt4List      ! list of Int4
      module procedure getPosVec2List      ! list of Vec2
      module procedure getPosReal8List     ! list of Real8
   end interface

   interface rewindList
      module procedure rewindInt4List      ! list of Int4
      module procedure rewindVec2List      ! list of Vec2
      module procedure rewindReal8List     ! list of Real8
   end interface

   interface windList
      module procedure windInt4List        ! list of Int4
      module procedure windVec2List        ! list of Vec2
      module procedure windReal8List       ! list of Real8
   end interface

   interface iterateList
      module procedure iterateInt4List     ! list of Int4
      module procedure iterateVec2List     ! list of Vec2
      module procedure iterateReal8List    ! list of Real8
   end interface

   interface getValList
      module procedure getValInt4List      ! list of Int4
      module procedure getValVec2List      ! list of Vec2
      module procedure getValReal8List     ! list of Real8
   end interface

   interface insertList
      module procedure insertInt4List      ! list of Int4
      module procedure insertVec2List      ! list of Vec2
      module procedure insertReal8List     ! list of Real8
   end interface

   interface removeList
      module procedure removeInt4List      ! list of Int4
      module procedure removeVec2List      ! list of Vec2
      module procedure removeReal8List     ! list of Real8
   end interface

   interface clearList
      module procedure clearInt4List       ! list of Int4
      module procedure clearVec2List       ! list of Vec2
      module procedure clearReal8List      ! list of Real8
   end interface

   interface sizeOfList
      module procedure sizeOfInt4List      ! list of Int4
      module procedure sizeOfVec2List      ! list of Vec2
      module procedure sizeOfReal8List     ! list of Real8
   end interface

   interface dumpList
      module procedure dumpInt4List        ! list of Int4
      module procedure dumpVec2List        ! list of Vec2
      module procedure dumpReal8List       ! list of Real8
   end interface
!
!  Public Data Types:
!  =================
!
   public :: Int4List                        ! Int4 list type
   public :: Real8List                       ! Real8 list type
   public :: Vec2List                        ! 2D vector list type
!
!=======================================================================
!
end module List_mod

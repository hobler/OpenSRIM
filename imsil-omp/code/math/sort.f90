module Math_Sort_mod
!
!  Quick sort algorithm.
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8
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
   public :: quickSort                ! Quick sort 
!            =========
   interface quickSort
      module procedure quickSortReal
      module procedure quickSortInt
   end interface
!
!  Private Procedures:
!  ==================
!
   private :: swop                    ! swap two numbers
!             ====
   interface swop
      module procedure swopReal
      module procedure swopInt
   end interface

contains
!
!===========================================================================================
!
   recursive subroutine quickSortReal( x, y, l, r )
!                       =============  io io i  i
!
!     The subroutine partitions into two vectors the vector to be sorted, then recursively
!     calls itself once for each partition.  This repeats until each partition is sorted.
!
!     X is the sorted version of the original X.
!     Y is integer array where each value Y(n) holds the position in the original X of 
!     the nth value in the returned X.
!
!     For example, if the original vectors are
!        X= (0.2, 0.3, 0.1),
!        Y= (1, 2, 3)
!     then the returned vectors are
!        X= (0.1, 0.2, 0.3)
!        Y= (3, 1, 2)
!     Thus, X is sorted, and Y gives the positions in the original X of each element of the 
!     returned X.
!
!  Arguments
!
      real(REAL8), intent(inout) :: x(:)     !  array to be sorted.
      integer(INT4), intent(inout) :: y(:)   !  [i, i=l:r] 
      integer(INT4), intent(in) :: l         !  lbound(x) 
      integer(INT4), intent(in) :: r         !  ubound(x) 
!
!  Local variables
!  
      integer :: l1, r1
!
!  Start of executable code
!  -------------------------
!
      if( l < r ) then
         l1 = l
         r1 = r

         do
            do while( l1 < r .and. x(l1) <= x(l) )
               l1 = l1 + 1
            end do

            do while( l < r1 .and. x(r1) >= x(l))
               r1 = r1 - 1
            end do

            if( l1 < r1 ) then
               call swop( x(l1), x(r1) )
!                   ====

               call swop( y(l1), y(r1) )
!                   ====
            end if
!
            if( l1 >= r1 )  exit        !Partition has been found.

         end do
!
!  Partition with x(l) at r1
!
         call swop( x(l), x(r1) )
!             ====  io     io
         call swop( y(l), y(r1) )
!             ====   io     io
!
!  Quick sort the left partition.
!
         call quickSort( x, y, l, r1-1 )
!             ========= io  io i   i
!
!  Quick sort the right partition.
!
         call quickSort( x, y, r1+1, r )
!             ========= io io   i    i
      end if

   end subroutine quickSortReal
!
!==========================================================================================
!
   recursive subroutine quickSortInt( x, y, l, r )
!                       ============  io io i  i
!
!     Sorts the integer values
!
!  Calling arguments
!
      integer(INT4), intent(inout) :: x(:)
      integer(INT4), intent(inout) :: y(:)
      integer(INT4), intent(in) :: l
      integer(INT4), intent(in) :: r
!
!  Local variables
!  
      integer :: l1, r1
!
!  Start of executable code
!  -------------------------
!
      if( l < r ) then
         l1 = l
         r1 = r

         do
            do while( l1 < r .and. x(l1) <= x(l) )
               l1 = l1 + 1
            end do

            do while( l < r1 .and. x(r1) >= x(l))
               r1 = r1 - 1
            end do

            if( l1 < r1 ) then
               call swop( x(l1), x(r1) )
!                   ====
               call swopInt( y(l1), y(r1) )
!                   =======
            end if

            if( l1 >= r1 )  exit !Partition has been found.

         end do
!
!  Partition with x(l) at r1
!
         call swop( x(l), x(r1) )
!             ====  io     io
!
         call swopInt( y(l), y(r1) )
!             =======   io     io
!
!  Quick sort the left partition.
!
         call quickSort( x, y, l, r1-1 )
!             ========= io  io i   i
!
!  Quick sort the right partition.
!
         call quickSort( x, y, r1+1, r )
!             ========= io io   i    i

      end if

   end subroutine quickSortInt
!
!==========================================================================================
!
   subroutine swopReal( a, b )
!             ========  io io
!
!      Subroutine SWOP swaps two real values.
!
!  Calling arguments
!
      real(REAL8), intent(inout) :: a,b
!
!  Local variables
!  
      real(REAL8)  :: temp
!
!  Start of executable code
!  -------------------------
!
      temp = a
      a = b
      b = temp

   end subroutine swopReal
!
!==========================================================================================
!
   subroutine swopInt( a, b )
!             =======  io io
!
!      Swaps two integer values.
!
!  Calling arguments
!
      integer(INT4), intent(inout) :: a, b
!
!  Local variables
!  
      integer(INT4) :: temp
!
!  Start of executable code
!  -------------------------
!
      temp = a
      a = b
      b = temp

   end subroutine swopInt
!
!==========================================================================================
!
end module Math_Sort_mod

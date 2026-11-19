module Math_BinSearch_mod
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
   public :: binSearch                ! Binary search
!            =========
!
contains
!
!===========================================================================================
!
   subroutine binSearch( n, array, val, i )
!             =========  i    i     i   o
!
!     Does a binary search in an array sorted in ascending order.
!     If value is outside the range of array, the first or last interval is returned
!
!  Arguments:
!
      integer(INT4), intent(in) :: n            ! dimension of array
      real(REAL8), intent(in) :: array(n)       ! the array to be searched in
      real(REAL8), intent(in) :: val            ! the value to be searched for
      integer(INT4), intent(out) :: i           ! the index of the left boundary of the
                                                !    interval containing val (1<=i<n)
      integer(INT4) :: il, im, iu
!
!  Start of Executable Code:
!  ------------------------
!
      il=1
      iu=n

      do
         im=(il+iu)/2
         if( val < array(im) ) then
            iu=im
         else
            il=im
         endif
         if( iu-il <= 1 )  exit
      enddo

      i=il

   end subroutine binSearch
!
!=======================================================================
!
end module Math_BinSearch_mod

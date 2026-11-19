module Math_IntMath_mod
!
!  Contains algorithms of integer mathematics.
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, INT8, REAL8
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
   public :: index1to2              ! convert index of 1d array to indices of 2d array
!            =========
   public :: index2to1              ! convert indices of 2d array to index of 1d array
!            =========
   public :: index1to3              ! convert index of 1d array to indices of 3d array
!            =========
   public :: index3to1              ! convert indices of 3d array to index of 1d array
!            =========
   public :: ggT                    ! greates common factor
!            ===

   interface index1to2
      module procedure index1to2_int4        ! 1D index is INT4
      module procedure index1to2_int8        ! 1D index is INT8
   end interface index1to2

   interface index2to1
      module procedure index2to1_int4        ! 1D index is INT4
      module procedure index2to1_int8        ! 1D index is INT8
   end interface index2to1

   interface index1to3
      module procedure index1to3_scalar     ! indices and extents are passed as scalars
      module procedure index1to3_scalar_int8! indices and extents are passed as scalars
      module procedure index1to3_array      ! indices and extents are passed as arrays
   end interface index1to3

   interface index3to1
      module procedure index3to1_scalar     ! indices and extents are passed as scalars
      module procedure index3to1_scalar_int8! indices and extents are passed as scalars
      module procedure index3to1_array      ! indices and extents are passed as arrays
   end interface index3to1

contains
!
!===========================================================================================
!
   subroutine index1to2_int4( i, m, j, k )
!             ==============  i  i  o  o
!
!     Converts the one-based index of a 1D array to the one-based indices of the corresponding
!     2D array:
!                                a1(i) = a2(j,k)
!
!     The shape of a1 is (m*..), while the shape of a2 is (m,..)
!
!  Arguments:
!
      integer(INT4), intent(in) :: i            ! index of 1D array
      integer(INT4), intent(in) :: m            ! first dimension of a2
      integer(INT4), intent(out) :: j           ! first index of 2D array
      integer(INT4), intent(out) :: k           ! second index of 2D array
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      j = mod( i-1, m ) + 1
      k = (i-1) / m + 1

   end subroutine index1to2_int4
!
!===========================================================================================
!
   subroutine index1to2_int8( i, m, j, k )
!             ==============  i  i  o  o
!
!     Converts the one-based index of a 1D array to the one-based indices of the corresponding
!     2D array:
!                                a1(i) = a2(j,k)
!
!     The shape of a1 is (m*..), while the shape of a2 is (m,..)
!
!  Arguments:
!
      integer(INT8), intent(in) :: i            ! index of 1D array
      integer(INT4), intent(in) :: m            ! first dimension of a2
      integer(INT4), intent(out) :: j           ! first index of 2D array
      integer(INT4), intent(out) :: k           ! second index of 2D array
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      j = mod( i-1, int( m, INT8 ) ) + 1
      k = (i-1) / m + 1

   end subroutine index1to2_int8
!
!===========================================================================================
!
   subroutine index2to1_int4( j, k, m, i )
!             ==============  i  i  i  o
!
!     Converts the one-based indices of a 2D array to the one-based index of the corresponding
!     1D array:
!                                a1(i) = a2(j,k)
!
!     The shape of a1 is (m*..), while the shape of a2 is (m,..)
!
!  Arguments:
!
      integer(INT4), intent(in) :: j            ! first index of 2D array
      integer(INT4), intent(in) :: k            ! second index of 2D array
      integer(INT4), intent(in) :: m            ! first dimension of a2
      integer(INT4), intent(out) :: i           ! index of 1D array
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      i = (k-1) * m + j

   end subroutine index2to1_int4
!
!===========================================================================================
!
   subroutine index2to1_int8( j, k, m, i )
!             ==============  i  i  i  o
!
!     Converts the one-based indices of a 2D array to the one-based index of the corresponding
!     1D array:
!                                a1(i) = a2(j,k)
!
!     The shape of a1 is (m*..), while the shape of a2 is (m,..)
!
!  Arguments:
!
      integer(INT4), intent(in) :: j            ! first index of 2D array
      integer(INT4), intent(in) :: k            ! second index of 2D array
      integer(INT4), intent(in) :: m            ! first dimension of a2
      integer(INT8), intent(out) :: i           ! index of 1D array
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      i = (k-1) * m + j

   end subroutine index2to1_int8
!
!===========================================================================================
!
   subroutine index1to3_scalar( i, m, n, j, k, l )
!             ================  i  i  i  o  o  o
!
!     Converts the one-based index of a 1D array to the one-based indices of the corresponding
!     3D array:
!                                a1(i) = a2(j,k,l)
!
!     The shape of a1 is (m*n*..), while the shape of a2 is (m,n,..)
!
!  Arguments:
!
      integer(INT4), intent(in) :: i            ! index of 1D array
      integer(INT4), intent(in) :: m            ! first dimension of a2
      integer(INT4), intent(in) :: n            ! second dimension of a2
      integer(INT4), intent(out) :: j           ! first index of 3D array
      integer(INT4), intent(out) :: k           ! second index of 3D array
      integer(INT4), intent(out) :: l           ! third index of 3D array
!
!  Local Variables:
!
      integer(INT4) :: ii
!
!  Start of Executable Code:
!  ------------------------
!
      call index1to2( i, m*n, ii, l )
!          =========  i   i   o   o
      call index1to2( ii, m, j, k )
!          =========  i   i  o  o

   end subroutine index1to3_scalar
!
!===========================================================================================
!
   subroutine index1to3_scalar_int8( i, m, n, j, k, l )
!             =====================  i  i  i  o  o  o
!
!     Converts the one-based index of a 1D array to the one-based indices of the corresponding
!     3D array:
!                                a1(i) = a2(j,k,l)
!
!     The shape of a1 is (m*n*..), while the shape of a2 is (m,n,..)
!
!  Arguments:
!
      integer(INT8), intent(in) :: i            ! index of 1D array
      integer(INT4), intent(in) :: m            ! first dimension of a2
      integer(INT4), intent(in) :: n            ! second dimension of a2
      integer(INT4), intent(out) :: j           ! first index of 3D array
      integer(INT4), intent(out) :: k           ! second index of 3D array
      integer(INT4), intent(out) :: l           ! third index of 3D array
!
!  Local Variables:
!
      integer(INT4) :: ii
!
!  Start of Executable Code:
!  ------------------------
!
      call index1to2( i, m*n, ii, l )
!          =========  i   i   o   o
      call index1to2( ii, m, j, k )
!          =========  i   i  o  o

   end subroutine index1to3_scalar_int8
!
!===========================================================================================
!
   subroutine index3to1_scalar( j, k, l, m, n, i )
!             ================  i  i  i  i  i  o
!
!     Converts the one-based indices of a 3D array to the one-based index of the corresponding
!     1D array:
!                                a1(i) = a2(j,k,l)
!
!     The shape of a1 is (m*n*..), while the shape of a2 is (m,n,..)
!
!  Arguments:
!
      integer(INT4), intent(in) :: j            ! first index of 3D array
      integer(INT4), intent(in) :: k            ! second index of 3D array
      integer(INT4), intent(in) :: l            ! third index of 3D array
      integer(INT4), intent(in) :: m            ! first dimension of a2
      integer(INT4), intent(in) :: n            ! second dimension of a2
      integer(INT4), intent(out) :: i           ! index of 1D array
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      i = (l-1) * (m*n) + (k-1) * m + j

   end subroutine index3to1_scalar
!
!===========================================================================================
!
   subroutine index3to1_scalar_int8( j, k, l, m, n, i )
!             =====================  i  i  i  i  i  o
!
!     Converts the one-based indices of a 3D array to the one-based index of the corresponding
!     1D array:
!                                a1(i) = a2(j,k,l)
!
!     The shape of a1 is (m*n*..), while the shape of a2 is (m,n,..)
!
!  Arguments:
!
      integer(INT4), intent(in) :: j            ! first index of 3D array
      integer(INT4), intent(in) :: k            ! second index of 3D array
      integer(INT4), intent(in) :: l            ! third index of 3D array
      integer(INT4), intent(in) :: m            ! first dimension of a2
      integer(INT4), intent(in) :: n            ! second dimension of a2
      integer(INT8), intent(out) :: i           ! index of 1D array
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      i = (l-1) * (m*n) + (k-1) * m + j

   end subroutine index3to1_scalar_int8
!
!===========================================================================================
!
   subroutine index1to3_array( i, shape, j )
!             ===============  i    i    o
!
!     Converts the one-based index of a 1D array to the one-based indices of the corresponding
!     3D array:
!                                a1(i) = a2(j(1),j(2),j(3))
!
!     The shape of a1 is (shape(1)*shape(2)*shape(3)), while the shape of a2 is
!     (shape(1),shape(2),shape(3))
!
!  Arguments:
!
      integer(INT4), intent(in) :: i            ! index of 1D array
      integer(INT4), intent(in) :: shape(3)     ! dimensions of a2
      integer(INT4), intent(out) :: j(3)        ! indices of 3D array
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call index1to3( i, shape(1), shape(2), j(1), j(2), j(3) )
!          =========  i     i         i        o    o     o

   end subroutine index1to3_array
!
!===========================================================================================
!
   subroutine index3to1_array( j, shape, i )
!             ===============  i    i    o
!
!     Converts the one-based indices of a 3D array to the one-based index of the corresponding
!     1D array:
!                                a1(i) = a2(j(1),j(2),j(3))
!
!     The shape of a1 is (shape(1)*shape(2)*shape(3)), while the shape of a2 is
!     (shape(1),shape(2),shape(3))
!
!  Arguments:
!
      integer(INT4), intent(in) :: j(3)         ! indices of 3D array
      integer(INT4), intent(in) :: shape(3)     ! dimensions of a2
      integer(INT4), intent(out) :: i           ! index of 1D array
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call index3to1( j(1), j(2), j(3), shape(1), shape(2), i )
!          =========   i     i     i       i         i      o

   end subroutine index3to1_array
!
!===========================================================================================
!
   integer(INT4) function ggT( i, j )
!                         ===  i  i
!
!     Calculates the greates common factor of i and j
!     If value is outside the range of array, the first or last interval is returned
!
!  Arguments:
!
      integer(INT4), intent(in) :: i            ! first number
      integer(INT4), intent(in) :: j            ! second number
!
!  Local Variables:
!
      integer(INT4) :: dividend, divisor, remainder
!
!  Start of Executable Code:
!  ------------------------
!
      dividend = i
      divisor = j

      do
         remainder = mod( dividend, divisor )
         if( remainder == 0 )  exit
         dividend = divisor
         divisor = remainder
      end do

      ggT = divisor

   end function ggT
!
!=======================================================================
!
end module Math_IntMath_mod

module Math_Solve_mod
!
!  Solve linear systems of equations
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
   save
!
!  Public Procedures
!  =================
!
   public :: tridiag                         ! initialize the table
!            =======
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine tridiag( al, am, au, b, x )
!             =======  i   i   i   i  o
!
!     Solve the tridiagonal system of linear equations
!
!        / am(1) au(1)                         \   / x(1) \   / b(1) \
!        | al(2) am(2) au(2)                   |   | x(2) |   | b(2) |
!        |       al(3) am(3) au(3)             | * | x(3) | = | b(3) |
!        |               .     .     .         |   |  .   |   |  .   |
!        |                     .     .     .   |   |  .   |   |  .   |
!        \                         al(n) am(n) /   \ x(n) /   \ b(n) /
!
!     This is done by Gaussian elimination, i.e., by bringing the system of the equations
!     into the form
!
!        /   1   c(1)                       \   / x(1) \   / d(1) \
!        |        1   c(2)                  |   | x(2) |   | d(2) |
!        |             1   c(3)             | * | x(3) | = | d(3) |
!        |                  .    .          |   |  .   |   |  .   |
!        |                       .   c(n-1) |   |  .   |   |  .   |
!        \                             1    /   \ x(n) /   \ d(n) /
!
!     (forward elimination) and then solving for x(i), i=n..1 (backward substitution)
!
!  Arguments:
!
      real(REAL8), intent(in) :: al(:)          ! lower diagonal
      real(REAL8), intent(in) :: am(:)          ! main diagonal
      real(REAL8), intent(in) :: au(:)          ! upper diagonal
      real(REAL8), intent(in) :: b(:)           ! right-hand side
      real(REAL8), intent(out) :: x(:)          ! solution vector
!
!  Local Variables:
!
      integer(INT4) :: i, n
      real(REAL8), allocatable :: c(:), d(:), div
!
!  Start of Executable Code:
!  ------------------------
!
      n = size( b )
      allocate( c( n ), d( n ) )
!
!  Forward elimination
!
      c(1) = au(1) / am(1)
      d(1) = b(1) / am(1)
      do i = 2, n
         div = am(i) - c(i-1) * al(i)
         c(i) = au(i) / div
         d(i) = ( b(i) - d(i-1) * al(i) ) / div
      end do
!
!  Backward substitution
!
      x(n) = d(n)
      do i = n-1, 1, -1
         x(i) = d(i) - c(i) * x(i+1)
      end do

      deallocate( c, d )

   end subroutine tridiag
!
!=======================================================================
!
end module Math_Solve_mod

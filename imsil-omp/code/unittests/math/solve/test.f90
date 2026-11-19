module Test_Solve_mod
!
!     Tests the Solve_mod module
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO
   use Math_Solve_mod
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
   public :: testTridiag                        ! test tridiagonal solver
!            ===========
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine testTridiag()
!
!     Test the tridiagonal solver
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: i, j, n
      real(REAL8) :: result
      real(REAL8), allocatable :: al(:), am(:), au(:), b(:), x(:)
!
!  Start of Executable Code:
!  ------------------------
!
      print *, 'Enter number of equations'
      read *, n
      allocate( al(n), am(n), au(n), b(n), x(n) )

      do
         call random_number( al )
!             =============  o
         call random_number( am )
!             =============  o
         call random_number( au )
!             =============  o
         call random_number( b )
!             =============  o
         al = 2*al-1
         am = 2*am-1
         au = 2*au-1
         b = 2*b-1
         call tridiag( al, am, au, b, x )
!             =======  i   i   i   i  o
         do i = 1, n
            write( *, '(a)', advance='no' )  '|'
            do j = 1, i-2
               write( *, '(i5)', advance='no' )  0
            end do
            if( i /= 1 )  write( *, '(i5)', advance='no' )  nint(1000*al(i))
            write( *, '(i5)', advance='no' )                nint(1000*am(i))
            if( i /= n )  write( *, '(i5)', advance='no' )  nint(1000*au(i))
            do j = i+2, n
               write( *, '(i5)', advance='no' )  0
            end do
            if( i == (n+1)/2 ) then
               write( *, '(a)', advance='no' )  ' | * |'
            else
               write( *, '(a)', advance='no' )  ' |   |'
            end if
            write( *, '(es10.3)', advance='no' )  x(i)
            if( i == (n+1)/2 ) then
               write( *, '(a)', advance='no' ) ' | = |'
            else
               write( *, '(a)', advance='no' ) ' |   |'
            end if
            if( i == 1 ) then
               result = am(i)*x(i) + au(i)*x(i+1)
            else if( i == n ) then
               result = al(i)*x(i-1) + am(i)*x(i)
            else
               result = al(i)*x(i-1) + am(i)*x(i) + au(i)*x(i+1)
            end if
            write( *, '(i5)', advance='no' )  nint(1000*result)
            if( i == (n+1)/2 ) then
               write( *, '(a)', advance='no' ) ' | = |'
            else
               write( *, '(a)', advance='no' ) ' |   |'
            end if
            write( *, '(i5)', advance='no' )  nint(1000*b(i))
            write( *, '(a)' )  ' |'
         end do

         print *, 'Enter <CR> to continue'
         read '(a)'

      end do

   end subroutine testTridiag
!
!=======================================================================
!
end module Test_Solve_mod
!
!#######################################################################
!
program test

   use const_mod, only: INT4
   use Test_Solve_mod

   integer(INT4) :: choice

   do
      print *, ' '
      print *, 'Choose among the following tests:'
      print *, ' '
      print *, '   1 ... test tridiagonal equation solver'
      print *, ' '
      read *, choice

      select case( choice )
         case( 1 )
            call testTridiag()
         case default
            print *, 'Invalid answer.'
      end select
   end do

end program test

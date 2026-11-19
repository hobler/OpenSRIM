module Random_Haas_mod
!
!  Provides the uniform random number generator according to
!
!  A.Haas, 'The multiple prime random number generator', ACM Trans. Math.
!  Software 13 (4) 368-381, 1987.
!
!  Used Modules:
!  ============
!
!$ use omp_lib,         only: omp_get_thread_num
   use Const_mod,       only: REAL8, INT4
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
   public :: ransetHaas             ! set the seed of the uniform random number generator
!            ==========
   public :: rangetHaas             ! get the seed of the uniform random number generator
!            ==========
   public :: rannumHaas             ! generate random number uniformly distributed in [0,1]
!            ==========
!
!  Private Data:
!  ============
!
   integer(INT4) :: m = 6103        ! seed of Haas uniformly distributed random number
   integer(INT4) :: ia = 15678      ! generator (own choice)
   integer(INT4) :: ib = 67638      ! Changing irand has a limited effect on the
   integer(INT4) :: irand = 5421    ! sequence of random numbers (see Reference)

!   integer(INT4) :: m = 971     ! seed of Haas uniformly distributed random number
!   integer(INT4) :: ia = 11113  ! generator that was used in the original paper
!   integer(INT4) :: ib = 104322 ! (without shift)
!   integer(INT4) :: irand = 4181
!$omp threadprivate( m, ia, ib, irand )
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine ransetHaas( mrand, iarand, ibrand, jrand )
!             ==========    i      i       i       i
!
!     Set the seed of the random number generator
!
!     The values of m, ia, ib, and irand are mapped into the intervals
!                    1000 < m < 9000
!                  20000 < ia < 90000
!                 130000 < ib < 220000
!                       irand < 10000,
!     respectively.
!
!     Changing only irand results in an identical series of
!     random numbers after a relatively small number of numbers.
!
!     This subroutine must be called from a parallel region. It sets different seed for each
!     thread.
!
!  Arguments:
!
      integer(INT4), intent(in) :: mrand      ! seed
      integer(INT4), intent(in) :: iarand     !
      integer(INT4), intent(in) :: ibrand     !
      integer(INT4), intent(in) :: jrand      !
!
!  Local Variables:
!
      integer(INT4) :: ithread, mrand_, iarand_, ibrand_, jrand_
!
!  Start of Executable Code:
!  ------------------------
!
! Modify the seed to get different seeds on different threads
!
      ithread = 1
!$    ithread = 1 + omp_get_thread_num( )
      mrand_ = modulo( mrand + (ithread-2), 7999 ) + 1
      iarand_ = modulo( iarand + (ithread-2), 69999 ) + 1
      ibrand_ = modulo( ibrand + (ithread-2), 89999 ) + 1
      jrand_ = modulo( jrand + (ithread-2), 9999 ) + 1

      if( mrand_ < 1 .or. mrand_ > 8000)  goto 931
      if( iarand_ < 1 .or. iarand_ > 70000 )  goto 932
      if( ibrand_ < 1 .or. ibrand_ > 90000 ) goto 933
      if( jrand_ < 1 .or. jrand_ > 10000 )  goto 934

      m = 1000 + mrand_
      ia = 20000 + iarand_
      ib = 130000 + ibrand_
      irand = jrand_

      return

931   error stop 'ransetHaas: 0<mrand<=8000 required'
932   error stop 'ransetHaas: 0<iarand<=70000 required'
933   error stop 'ransetHaas: 0<ibrand<=90000 required'
934   error stop 'ransetHaas: 0<irand<=10000 required'

   end subroutine ransetHaas
!
!=======================================================================
!
   subroutine rangetHaas( mrand, iarand, ibrand, jrand )
!             ==========    o      o       o       o
!
!     Get the seed of the random number generator
!
!  Arguments:
!
      integer(INT4), intent(out) :: mrand    ! seed
      integer(INT4), intent(out) :: iarand   !
      integer(INT4), intent(out) :: ibrand   !
      integer(INT4), intent(out) :: jrand    !
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      mrand = m
      iarand = ia
      ibrand = ib
      jrand = irand

   end subroutine rangetHaas
!
!=======================================================================
!
   subroutine rannumHaas( rand )
!             ==========   o
!
!     Generate uniformly distributed random numbers in [0,1]
!
!  Arguments:
!
      real(REAL8), intent(out) :: rand          ! random number
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      m = m + 7
      if( m >= 9973 )  m = m - 9871
      ia = ia + 1907
      if( ia >= 99991 )  ia = ia - 89989
      ib = ib + 73939
      if( ib >= 224729 )  ib = ib - 96233
      irand = mod( irand*m+ia+ib, 100000 ) / 10
      rand = 10000 * irand

      m = m + 7
      if( m >= 9973 )  m = m - 9871
      ia = ia + 1907
      if( ia >= 99991 )  ia = ia - 89989
      ib = ib + 73939
      if( ib >= 224729 )  ib = ib - 96233
      irand = mod( irand*m+ia+ib, 100000 ) / 10
      rand = (rand+irand+0.5) / 100000000.

   end subroutine rannumHaas
!
!=======================================================================
!
end module Random_Haas_mod

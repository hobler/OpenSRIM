module Random_Knuth_mod
!
!  Knuth's subtractive random number generator Ranfib as described in
!
!  W.H. Press et al., Numerical Recipes. The Art of Scientific Computing, 3rd edition,
!  Cambridge University Press, 2007; p.354
!
!  It provides a uniform random number in the range [0,1). It provides different random
!  number sequences on different threads.
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: REAL8, INT4,  &
                              ZERO, ONE
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
   public :: rannumKnuth            ! generate random number uniformly distributed in [0,1]
!            ===========
   public :: initKnuth              ! initialize the random number generator
!            =========
!
!  Private Data:
!  ============
!
   integer(INT4) :: inext1, inext2  ! indices of next random numbers to be used
   real(REAL8) :: rantab(0:54)      ! table of random numbers
!$omp threadprivate( inext1, inext2, rantab )
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initKnuth( )
!             =========
!
!     Initialize the random number generator
!
!     This subroutine must be called from a parallel region
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      inext1 = 0
      inext2 = 31
      call random_number( rantab )
!          =============    o

   end subroutine initKnuth
!
!=======================================================================
!
   subroutine rannumKnuth( rand )
!             ===========   o
!
!     Generate a uniformly distributed random number in [0,1]
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
      inext1 = mod( inext1+1, 55 )
      inext2 = mod( inext2+1, 55 )
      rand = rantab(inext1) - rantab(inext2)
      if( rand < ZERO )  rand = rand + ONE
      rantab(inext1) = rand

   end subroutine rannumKnuth
!
!=======================================================================
!
end module Random_Knuth_mod

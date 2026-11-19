module Random_Fortran_mod
!
!  Wrapper for setting and getting the seed of the Fortran default random number generator
!  random_seed()
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
   public :: ransetFortran          ! set the seed of the uniform random number generator
!            =============
   public :: rangetFortran          ! get the seed of the uniform random number generator
!            =============
   public :: rannumFortran          ! generate uniformly distributed random numbers in [0,1]
!            =============
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine ransetFortran( seed1, seed2, seed3, seed4 )
!             =============    i      i     i      i
!
!     Initialize the random number generator
!
!     It is assumed that it is sufficient to initialize the random number generator on one
!     thread and that, nevertheless, different random numbers are generated on different
!     threads by random_number().
!
!  Arguments:
!
      integer(INT4), intent(in) :: seed1  ! currently corresponding to IMSIL parameter m
      integer(INT4), intent(in) :: seed2  ! currently corresponding to IMSIL parameter ia
      integer(INT4), intent(in) :: seed3  ! currently corresponding to IMSIL parameter ib
      integer(INT4), intent(in) :: seed4  ! currently corresponding to IMSIL parameter irand
!
!  Local Variables:
!
      integer :: nseed
      integer, allocatable :: seed(:)
      integer :: topos(4)
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
!
!  Check validity of the seed (to avoid different seeds to produce identical results)
!
      if( seed1 >= 2**14 .or. seed2 >= 2**17 .or. seed3 >= 2**18 .or. seed4 >= 2**14 )  &
            error stop 'ransetFortran: Seed to large'
      if( seed1 <=0 .or. seed2 <= 0 .or. seed3 <= 0 .or. seed4 <= 0 )  &
            error stop 'ransetFortran: Seed must be positive'
!
!  Allocate the seed array
!
      call random_seed( size=nseed )
!          ===========       o
      allocate( seed( nseed ) )
      seed(:) = 0
!
!  Cheat the compiler so it does not complain about topos + len > bit_size.
!  If this does not work, comment mvbits below that are not used
!
      if( nseed == 1 ) then
         topos(:) = [ 0, 14, 31, 49 ]
      else if( nseed == 2 ) then
         topos(:) = [ 0, 14, 0, 18 ]
      else if( nseed == 3 ) then
         topos(:) = [ 0, 14, 0, 0 ]
      else
         topos(:) = 0
      end if
!
!  Convert the seed
!
      if( nseed == 1 ) then
         call mvbits( int(seed1), 0, 14, seed(1), topos(1) )
         call mvbits( int(seed2), 0, 17, seed(1), topos(2) )
         call mvbits( int(seed3), 0, 18, seed(1), topos(3) )
         call mvbits( int(seed4), 0, 14, seed(1), topos(4) )
      else if( nseed == 2 ) then
         call mvbits( int(seed1), 0, 14, seed(1), topos(1) )
         call mvbits( int(seed2), 0, 17, seed(1), topos(2) )
         call mvbits( int(seed3), 0, 18, seed(2), topos(3) )
         call mvbits( int(seed4), 0, 14, seed(2), topos(4) )
      else if( nseed == 3 ) then
         call mvbits( int(seed1), 0, 14, seed(1), topos(1) )
         call mvbits( int(seed2), 0, 17, seed(1), topos(2) )
         call mvbits( int(seed3), 0, 18, seed(2), topos(3) )
         call mvbits( int(seed4), 0, 14, seed(3), topos(4) )
      else
         call mvbits( int(seed1), 0, 14, seed(1), topos(1) )
         call mvbits( int(seed2), 0, 17, seed(2), topos(2) )
         call mvbits( int(seed3), 0, 18, seed(3), topos(3) )
         call mvbits( int(seed4), 0, 14, seed(4), topos(4) )
      end if
!
!  Set the seed
!
      call random_seed( put=seed(:) )
!          ===========        i
!$omp end master
!$omp barrier

   end subroutine ransetFortran
!
!=======================================================================
!
   subroutine rangetFortran( seed1, seed2, seed3, seed4 )
!             =============    o      o     o      o
!
!     Initialize the random number generator
!
!  Arguments:
!
      integer(INT4), intent(out) :: seed1 ! currently corresponding to IMSIL parameter m
      integer(INT4), intent(out) :: seed2 ! currently corresponding to IMSIL parameter ia
      integer(INT4), intent(out) :: seed3 ! currently corresponding to IMSIL parameter ib
      integer(INT4), intent(out) :: seed4 ! currently corresponding to IMSIL parameter irand
!
!  Local Variables:
!
      integer :: nseed
      integer, allocatable :: seed(:)
      integer :: pos(4)
!
!  Start of Executable Code:
!  ------------------------
!
!  Allocate the seed array
!
      call random_seed( size=nseed )
!          ===========       o
      allocate( seed( nseed ) )
!
!  Get the seed
!
      call random_seed( get=seed(:) )
!          ===========        i
!
!  Cheat the compiler so it does not complain about topos + len > bit_size.
!  If this does not work, comment mvbits below that are not used
!
      if( nseed == 1 ) then
         pos(:) = [ 0, 14, 31, 49 ]
      else if( nseed == 2 ) then
         pos(:) = [ 0, 14, 0, 18 ]
      else if( nseed == 3 ) then
         pos(:) = [ 0, 14, 0, 0 ]
      else
         pos(:) = 0
      end if
!
!  Convert the seed
!
      seed1 = 0
      seed2 = 0
      seed3 = 0
      seed4 = 0
      if( nseed == 1 ) then
         seed1 = ibits( seed(1), pos(1), 14 )
         seed2 = ibits( seed(1), pos(2), 17 )
         seed3 = ibits( seed(1), pos(3), 18 )
         seed4 = ibits( seed(1), pos(4), 14 )
      else if( nseed == 2 ) then
         seed1 = ibits( seed(1), pos(1), 14 )
         seed2 = ibits( seed(1), pos(2), 17 )
         seed3 = ibits( seed(2), pos(3), 18 )
         seed4 = ibits( seed(2), pos(4), 14 )
      else if( nseed == 3 ) then
         seed1 = ibits( seed(1), pos(1), 14 )
         seed2 = ibits( seed(1), pos(2), 17 )
         seed3 = ibits( seed(2), pos(3), 18 )
         seed4 = ibits( seed(3), pos(4), 14 )
      else
         seed1 = ibits( seed(1), pos(1), 14 )
         seed2 = ibits( seed(2), pos(2), 17 )
         seed3 = ibits( seed(3), pos(3), 18 )
         seed4 = ibits( seed(4), pos(4), 14 )
         if( nseed > 4 ) then
            print *, 'WARNING: Seed incomplete to continue random number sequence'
         end if
      end if

   end subroutine rangetFortran
!
!=======================================================================
!
   subroutine rannumFortran( rand )
!             =============   o
!
!     Generate uniformly distributed random numbers in [0,1] using the default Fortran
!     random number generator.
!
!     This layer is necessary as random_number may not be used as a procedure target.
!
!  Arguments:
!
      real(REAL8), intent(out) :: rand    ! random number
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call random_number( rand )
!          =============   o

   end subroutine rannumFortran
!
!=======================================================================
!
end module Random_Fortran_mod

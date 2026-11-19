module Random_mod
!
!  Provide random number generators
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: REAL8, INT4, &
                                 ZERO, ONE, TWO, THIRD

   use Param_Setup_mod,    only: rng, mrand, iarand, ibrand, irand, randskip

   use Random_Haas_mod,    only: rannumHaas, ransetHaas, rangetHaas
   use Random_Knuth_mod,   only: rannumKnuth, initKnuth
   use Random_Fortran_mod, only: rannumFortran, ransetFortran, rangetFortran
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
   public :: initRannum             ! initialize the rannum procedure pointer
!            ==========

   procedure( rannumHaas ), pointer, public :: rannum => null( )
!                                              ======
                                    ! generate random number uniformly distributed in [0,1]

   procedure( rangetHaas ), pointer, public :: ranget
!                                              ======
                                    ! get the seed of the uniform random number generator

   procedure( ransetHaas ), pointer, public :: ranset
!                                              ======
                                    ! set the seed of the uniform random number generator

   public :: ranno                  ! generate random number normal distributed with mean 
!            =====                  !    value 0 and standard deviation 1
   public :: rnnget                 ! get the status of the normal random number generator
!            ======
   public :: rnnset                 ! set the status of the normal random number generator
!            ======

   public :: ranint                 ! generate discrete random numbers with given
!            ======                 !    probabilities
   public :: ranpo                  ! generate Poisson distributed random number 
!            =====
   public :: randir                 ! generate random direction (point on unit sphere)
!            ======
   public :: ranrot                 ! generate a random rotation matrix (3 orthogonal vectors
!            ======                 !    from origin to points on unit sphere)
!
!  Private Data:
!  ============
!
   real(REAL8) :: rand2             ! second normal random number
   logical :: first = .true.        ! indicates that a new pair of normal random numbers is
                                    !    to be calculated
!$omp threadprivate( rand2, first )
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initRannum( )
!             ==========
!
!     Initialize the rannum procedure pointer
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: i
      real(REAL8) :: rand
!
!  Start of Executable Code:
!  ------------------------
!
!  Select the random number generator
!
!$omp master
      select case( rng )
         case( 'Haas' )
            rannum => rannumHaas
            ranset => ransetHaas
            ranget => rangetHaas
         case( 'Fortran' )
            rannum => rannumFortran
            ranset => ransetFortran
            ranget => rangetFortran
         case( 'Knuth' )
            rannum => rannumKnuth
            ranset => ransetFortran
            ranget => rangetFortran
         case default
      end select
!$omp end master
!$omp barrier
!
!  Set the seed
!
      call ranset( mrand, iarand, ibrand, irand )
!          ======    i      i       i       i
      if( rng == 'Knuth' )  call initKnuth( )
!                                =========
!
!  Skip the first randskip random numbers
!
      do i = 1, randskip
         call rannum( rand )
!             ======   o
      end do

   end subroutine initRannum
!
!=======================================================================
!
   subroutine ranno( rand )
!             =====   o
!
!     Generate normal distributed random numbers with mean value 0 and
!     standard deviation 1     
!
!  Arguments:
!
      real(REAL8), intent(out) :: rand
!
!  Local Variables:
!
      real(REAL8) :: u1, u2, v1, v2, s, fac
!
!  Start of Executable Code:
!  ------------------------
!
      if( first ) then
10       continue
         call rannum( u1 )
!             ======  o
         call rannum( u2 )
!             ======  o
         v1 = two*u1 - ONE
         v2 = two*u2 - ONE
         s = v1**2 + v2**2
         if( s >= ONE .or. s == ZERO )  goto 10
         
         fac = sqrt( -TWO*log(s)/s )
         rand = v1 * fac
         rand2 = v2 * fac
      else
         rand = rand2
      endif

      first = .not. first      
!
   end subroutine ranno
!
!=======================================================================
!
   subroutine rnnget( firstget, randget )
!             ======     i         i
!
!     Get the status of the normal random number generator
!
!  Arguments:
!
      logical, intent(out) :: firstget ! flag indicating that new random numbers
                                       !    need to be generated at next call
      real(REAL8) :: randget           ! second random number to be returned at
                                       !    next call if first=.false.
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      firstget = first
      randget = rand2

   end subroutine rnnget
!
!=======================================================================
!
   subroutine rnnset( randset_ )
!             ======     i
!
!     Get the status of the normal random number generator
!
!  Arguments:
!
      real(REAL8), intent(in), optional :: randset_! second random number to be returned at
                                                   !    next call of ranno if not specified,
                                                   !    a new number will be calculated
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( randset_ ) ) then
         first = .false.
         rand2 = randset_
      else
         first = .true.
      end if

   end subroutine rnnset
!
!=======================================================================
!
   subroutine ranint( prob, krand )
!             ======   i      o
!  
!     Generate a random number (krand) with given probabilities (prob)
!
!     sum( prob(:) ) is assumed to equal 1.
!
!  Arguments:
!
      real(REAL8), intent(in) :: prob(:)        ! probabilities of numbers 1:size( prob )
      integer(INT4), intent(out) :: krand       ! random number
!
!  Local Variables:
!
      real(REAL8) :: rand
!
!  Start of Executable Code:
!  ------------------------
!
      if( size( prob ) == 1 ) then

         krand = 1

      else

         call rannum( rand )
!             ======   o

         do krand = 1, size( prob ) - 1
            rand = rand - prob(krand)
            if( rand < ZERO )  exit
         end do

      end if

   end subroutine ranint
!
!=======================================================================
!
   subroutine ranpo( mean, krand )
!             =====   i      o
!
!     Generate Poisson distributed random number with mean value mean
!
!     In order to avoid possibly infinite computation times, the upper end
!     of the distribution function is cut off. 
!
!     The mean execution time is approximately proportional to mean.
!
!  Arguments:
!
      real(REAL8), intent(in) :: mean           ! mean value
      integer(INT4), intent(out) :: krand       ! random number
!
!  Local Variables:
!
      integer(INT4), parameter :: ACC = 1e-7    ! cumulative probability of distribution cut off
      real(REAL8) :: rand, probk, probks
!
!  Start of Executable Code:
!  ------------------------
!
      call rannum( rand )
!          ======   o

      rand = min( rand, ONE-ACC )

      krand = 0
      probk = exp(-mean)
      probks = probk

 10   if( probks > rand ) return
      krand = krand+1
      probk = probk*mean/krand
      probks = probks + probk
      goto 10

   end subroutine ranpo
!
!=======================================================================
!
   subroutine randir( dir )
!             ======   o
!
!     Generate a random direction (point on the unit sphere)
!
!     Algorithm due to G. Marsaglia, Ann. Math. Stat. 43 (2) 645, 1972,
!     also described in W. Press et al., Numerical Recipes 3rd Ed., Cambridge University Press,
!     New York, 2007, pp.1128-1131.
!
!  Arguments:
!
      real(REAL8), intent(out) :: dir(3)        ! direction vector
!
!  Local Variables:
!
      real(REAL8) :: u1, u2, rr, s
!
!  Start of Executable Code:
!  ------------------------
!
!  Get point inside unit circle
!
 10   call rannum( u1 )
!          ======  o
      call rannum( u2 )
!          ======  o
      rr = u1**2 + u2**2
      if( rr > ONE ) goto 10
!
!  Get point on unit sphere
!
      s = sqrt(ONE - rr)
      dir(1) = TWO*u1*s
      dir(2) = TWO*u2*s
      dir(3) = ONE - TWO*rr

   end subroutine randir
!
!=======================================================================
!
   subroutine ranrot( dirx, diry, dirz )
!             ======   o     o     o
!
!     Generate a random rotation matrix (three orthogonal vectors from the origin to points on
!     the unit sphere)
!
!     Algorithm due to K. Shoemake, Proc. SIGGRAPH'85, 245, 1985.
!     also described in W. Press et al., Numerical Recipes 3rd Ed., Cambridge University Press,
!     New York, 2007, pp.1128-1131.
!
!  Arguments:
!
      real(REAL8), intent(out) :: dirx(3)        ! first direction vector
      real(REAL8), intent(out) :: diry(3)        ! second direction vector
      real(REAL8), intent(out) :: dirz(3)        ! third direction vector
!
!  Local Variables:
!
      real(REAL8) :: u1, u2, u3, u4, rr12, rr34, s, x1, x2, x3, x4
      real(REAL8) :: x11, x12, x13, x14, x22, x23, x24, x33, x34
!
!  Start of Executable Code:
!  ------------------------
!
!  Get two points inside unit circle
!
 10   call rannum( u1 )
!          ======  o
      call rannum( u2 )
!          ======  o
      rr12 = u1**2 + u2**2
      if( rr12 > ONE ) goto 10

 20   call rannum( u3 )
!          ======  o
      call rannum( u4 )
!          ======  o
      rr34 = u3**2 + u4**2
      if( rr34 > ONE ) goto 20
!
!  Get point on 4D unit sphere
!
      s = sqrt((ONE - rr12) / rr34)
      x1 = u1
      x2 = u2
      x3 = u3*s
      x4 = u4*s
!
!  Get products needed for the matrix
!
      x11 = x1*x1
      x12 = x1*x2
      x13 = x1*x3
      x14 = x1*x4
      x22 = x2*x2
      x23 = x2*x3
      x24 = x2*x4
      x33 = x3*x3
      x34 = x3*x4
!
!  Get matrix
!
      dirx(:) = [ ONE - TWO*(x22+x33), TWO*(x12-x34), TWO*(x13+x24) ]
      diry(:) = [ TWO*(x12+x34), ONE - TWO*(x11+x33), TWO*(x23-x14) ]
      dirz(:) = [ TWO*(x13-x24), TWO*(x23+x14), ONE - TWO*(x11+x22) ]

   end subroutine ranrot
!
!=======================================================================
!
end module Random_mod

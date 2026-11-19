module Stat_Mom1d_mod
!
!  Provide data type and manage 1d moments
!
!  The algorithm for the uncertainties of the moments takes into account that values from the
!  same cascade may be correlated while values from different cascades are uncorrelated.
!  Also, values may have weights.
!
!  Used Modules:
!  ============
!
   use Const_mod, only : INT4, REAL8,  &
                         ZERO, ONE, TWO, THREE, FOUR, SIX, EIGHT, NINE, HALF
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
   public :: initMom1d                       ! Initialize moments
!            =========
   public :: scoreMom1d                      ! Accumulate weighted power of values in scasc
!            ==========
   public :: updateMom1d                     ! Update s and ss using scasc
!            ===========
   public :: mergeIntoMasterMom1d            ! Merge moments of other threads into master
!            ====================
   public :: mergeMom1d                      ! Merge two Mom1d instances before evaluation
!            ==========
   public :: evalMom1d                       ! Evaluate mean..kurt
!            =========
   public :: writeMom1d                      ! Write moments to output file
!            ==========
   interface writeMom1d
      module procedure writeMom1d_scalar
      module procedure writeMom1d_array
   end interface
!
!  Private Procedures:
!  ==================
!
   private :: calcFormat            ! Calculate format for writing values
!             ==========
!
!  Public Data Types:
!  =================
!
   type, public :: Mom1d_type                ! Array of moments
!                  ==========      
      !private   ! does not allow writing whole structure to file with single write (gfortran)
                 ! needed for writing backup file
      integer(INT4) :: nscore                ! number of scores within cascade
      integer(INT4) :: ncasc                 ! number of cascades with scores
      real(REAL8) :: scasc(0:4)              ! sum within cascade
      real(REAL8) :: s(0:4)                  ! cumulated sum of scasc
      real(REAL8) :: ss(0:4,0:4)             ! cumulated sum of products of scasc
      real(REAL8) :: mean                    ! mean value
      real(REAL8) :: std                     ! standard deviation
      real(REAL8) :: skew                    ! skewness
      real(REAL8) :: kurt                    ! kurtosis
      real(REAL8) :: mean_err                ! standard deviation of mean value
      real(REAL8) :: std_err                 ! standard deviation of standard deviation
      real(REAL8) :: skew_err                ! standard deviation of skewness
      real(REAL8) :: kurt_err                ! standard deviation of kurtosis
   end type Mom1d_type
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initMom1d( mom1d )
!             =========    o
!
!     Initialize all sums with zero
!
!  Arguments:
!
      type(Mom1d_type), intent(out) :: mom1d    ! moments array
!
!  Local Variables: none
!
!
!  Start of Executable Code:
!  ------------------------
!
      mom1d%nscore = 0
      mom1d%ncasc = 0
      mom1d%scasc(:) = ZERO
      mom1d%s(:) = ZERO
      mom1d%ss(:,:) = ZERO

   end subroutine initMom1d
!
!=======================================================================
!
   subroutine scoreMom1d( mom1d, value, weight )
!             ==========   io      i      i
!
!     Accumulate the weighted powers of the values in the intra-cascade sums scasc
!
!  Arguments:
!
      type(Mom1d_type), intent(inout) :: mom1d  ! moments array
      real(REAL8), intent(in) :: value          ! value
      real(REAL8), intent(in) :: weight         ! weight
!
!  Local Variables:
!
      integer(INT4) :: i
!
!  Start of Executable Code:
!  ------------------------
!
      mom1d%nscore = mom1d%nscore + 1

      mom1d%scasc(0) = mom1d%scasc(0) + weight
      do i = 1, 4
         mom1d%scasc(i) = mom1d%scasc(i) + weight * value**i
      end do

   end subroutine scoreMom1d
!
!=======================================================================
!
   subroutine updateMom1d( mom1d )
!             ===========   io
!
!     Add contributions of scasc to the sums s and ss
!
!  Arguments:
!
      type(Mom1d_type), intent(inout) :: mom1d  ! moments array
!
!  Local Variables:
!
      integer(INT4) :: k1, k2
!
!  Start of Executable Code:
!  ------------------------
!
      if( mom1d%nscore == 0 )  return

      mom1d%ncasc = mom1d%ncasc + 1

      mom1d%s(:) = mom1d%s(:) + mom1d%scasc(:)
      do k1 = 0, 4
         do k2 = 0, k1
            mom1d%ss(k1,k2) = mom1d%ss(k1,k2) + mom1d%scasc(k1)*mom1d%scasc(k2)
         end do
      end do
!
!  Reset nscore and scasc
!
      mom1d%nscore = 0
      mom1d%scasc(:) = ZERO
      
   end subroutine updateMom1d
!
!=======================================================================
!
   subroutine mergeIntoMasterMom1d( mom )
!             ====================  io
!
!     Merge mom of other threads into master thread
!
!  Arguments:
!
      type(Mom1d_type), intent(inout) :: mom    ! moments to be merged
!
!  Local Variables:
!
      type(Mom1d_type), save :: mom_tmp         ! shared temporary array
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      call initMom1d( mom_tmp )
!          =========     o
!$omp end master

!$omp barrier
!$omp critical( critical_mom1d_mergeIntoMasterMom1d )
      call mergeMom1d( mom, mom_tmp )
!          ==========   i     io
!$omp end critical( critical_mom1d_mergeIntoMasterMom1d )

!$omp barrier
!$omp master
      mom = mom_tmp
!$omp end master

   end subroutine mergeIntoMasterMom1d
!
!===========================================================================================
!
   subroutine mergeMom1d( mom1, mom2 )
!             ==========   i     io
!
!     Merge mom1 into mom2
!
!  Arguments:
!
      type(Mom1d_type), intent(in) :: mom1      ! moments array
      type(Mom1d_type), intent(inout) :: mom2   ! moments array
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      mom2%nscore = mom2%nscore + mom1%nscore
      mom2%ncasc = mom2%ncasc + mom1%ncasc
      mom2%scasc(:) = mom2%scasc(:) + mom1%scasc(:)
      mom2%s(:) = mom2%s(:) + mom1%s(:)
      mom2%ss(:,:) = mom2%ss(:,:) + mom1%ss(:,:)

   end subroutine mergeMom1d
!
!=======================================================================
!
   subroutine evalMom1d( mom1d )
!             =========    io
!  
!     Evaluate mean, std, skew, kurtosis, and their standard deviations
!
!  Arguments:
!
      type(Mom1d_type), intent(inout) :: mom1d  ! moments array
!
!  Local Variables:
!
      integer(INT4) :: i, k, k1, k2, l, l1, l2
      integer(INT4), allocatable :: kk1(:), kk2(:)
      real(REAL8) :: kbinom(0:4,0:4)
      real(REAL8) :: n, mom(0:4), cenmom(2:4)
      real(REAL8) :: var1, var2, var3, var4, var_skew, var_kurt
      real(REAL8) :: cov(2:4,1:4)
!
!  Start of Executable Code:
!  ------------------------
!
      mom1d%mean = ZERO
      mom1d%std = ZERO
      mom1d%skew = ZERO
      mom1d%kurt = ZERO
      mom1d%mean_err = ZERO
      mom1d%std_err = ZERO
      mom1d%skew_err = ZERO
      mom1d%kurt_err = ZERO

      if( mom1d%ncasc == 0 )  return
!
!  Get raw moments (divide sums by sum of weights)
!
      n = mom1d%s(0)
      mom(:) = mom1d%s(:) / n
!
!  Mean value and error
!
      mom1d%mean = mom(1)
      var1 = mom1d%ss(1,1)/n**2 - mom(1)**2 / mom1d%ncasc
      mom1d%mean_err = sqrt( max( ZERO, var1 ) )
!
!  Set up binomial coefficients
!
      do k = 0, 4
         kbinom(0,k) = 1
         kbinom(k,k) = 1
         do l = 1, k-1
            kbinom(l,k) = kbinom(l-1,k-1) + kbinom(l,k-1)
         end do
      end do
!
!  Get central momemts
!
      do k = 2, 4
         cenmom(k) = mom(k)
         do l = 1, k
            cenmom(k) = cenmom(k) + kbinom(l,k) * mom(k-l) * (-mom(1))**l
         end do
      end do

      if( cenmom(2) <= ZERO )  return
!
!  Make ss symmetric
!
      do k1 = 0, 4
         do k2 = k1+1, 4
            mom1d%ss(k1,k2) = mom1d%ss(k2,k1)
         end do
      end do
!
!  Covariances
!
      do k = 2, 4
         cov(k,1) = ZERO
         do l = 0, k
            cov(k,1) = cov(k,1) + kbinom(l,k) * mom1d%ss(k-l,1) * (-mom(1))**l
         end do
         cov(k,1) = cov(k,1) / n**2
         cov(k,1) = cov(k,1) - cenmom(k) * mom(1) / mom1d%ncasc
      end do

      kk1 = [ 2, 3, 4, 3, 4 ]
      kk2 = [ 2, 2, 2, 3, 4 ]
      do i = 1, size( kk1 )
         k1 = kk1(i)
         k2 = kk2(i)
         cov(k1,k2) = ZERO
         do l1 = 0, k1
            do l2 = 0, k2
               cov(k1,k2) = cov(k1,k2) + kbinom(l1,k1) * kbinom(l2,k2)  &
                            * mom1d%ss(k1-l1,k2-l2) * (-mom(1))**(l1+l2)
            end do
         end do
         cov(k1,k2) = cov(k1,k2) / n**2
         cov(k1,k2) = cov(k1,k2) - cenmom(k1) * cenmom(k2) / mom1d%ncasc
      end do
!
!  Standard deviation and its error
!
      mom1d%std = sqrt( max( ZERO, cenmom(2) ) )
      if( mom1d%std == ZERO )  return
      var2 = cov(2,2)                           ! cov relates to moments with fixed mean
      mom1d%std_err = sqrt(var2) / (TWO*mom1d%std)
!
!  Skewness and its error
!
      mom1d%skew = cenmom(3) / mom1d%std**3
      var3 = cov(3,3) + (THREE*cenmom(2))**2*var1 - SIX*cenmom(2)*cov(3,1)
      var_skew = ( var3 + (THREE*HALF*cenmom(3)/cenmom(2))**2*var2  &
                   - THREE*cenmom(3)/cenmom(2)*cov(3,2) + NINE*cenmom(3)*cov(2,1) )  &
                 / cenmom(2)**3
      mom1d%skew_err = sqrt( max( ZERO, var_skew ) )
!
!  Kurtosis and its error
!
      mom1d%kurt = cenmom(4) / mom1d%std**4
      var4 = cov(4,4) + (FOUR*cenmom(3))**2*var1 - EIGHT*cenmom(3)*cov(4,1)
      var_kurt = ( var4 + (TWO*cenmom(4)/cenmom(2))**2*var2  &
                   - FOUR*cenmom(4)/cenmom(2)*cov(4,2)  &
                   + 16*cenmom(3)*cenmom(4)/cenmom(2)*cov(2,1) )  &
                 / cenmom(2)**4
      mom1d%kurt_err = sqrt( max(ZERO, var_kurt ) )

   end subroutine evalMom1d
!
!=======================================================================
!
   subroutine writeMom1d_scalar( mom1d, text, units, nmax )
!             =================    io    i      i     i
!
!     Evaluate the moments and confidence intervals and write them to
!     the output file
!
!     Note: this subroutine is presently unused and therefore untested
!
!  Arguments:
!
      type(Mom1d_type), intent(inout) :: mom1d  ! Moments array
      character(*), intent(in) :: text          ! Title text
      character(*), intent(in) :: units         ! Units
      integer(INT4), intent(in) :: nmax         ! Maximum order of moment to be calculated
!
!  Local Variables:
!
      character(8) :: fmt_val, fmt_err
!
!  Start of Executable Code:
!  ------------------------
!
!  Evaluate moments
!
      call evalMom1d( mom1d )
!          =========    io
!
!  Write moments to the output file
!
      write( 21,'(/t5,a)' ) trim( text )
      call calcFormat( mom1d%mean, mom1d%mean_err, fmt_val, fmt_err )
!          ==========     i             i             o        o
      write( 21,'(t8,a,' // fmt_val // ',a,' // fmt_err // ',a,a)' )  &
           'mean value                :', mom1d%mean, ' +/-', mom1d%mean_err, ' ', units

      if( mom1d%ncasc == 0 ) return

      if( nmax > 1) then
         call calcFormat( mom1d%std, mom1d%std_err, fmt_val, fmt_err )
!             ==========    i             i            o        o
         write( 21,'(t8,a,' // fmt_val // ',a,' // fmt_err // ',a,a)' )  &
              'standard deviation        :', mom1d%std, ' +/-', mom1d%std_err, ' ', units
      end if

      if( nmax > 2) then
         call calcFormat( mom1d%skew, mom1d%skew_err, fmt_val, fmt_err )
!             ==========      i             i            o        o
         write( 21,'(t8,a,' // fmt_val // ',a,' // fmt_err // ',a,a)' )  &
              'skewness                  :', mom1d%skew, ' +/-', mom1d%skew_err
      end if

      if( nmax > 3) then
         call calcFormat( mom1d%kurt, mom1d%kurt_err, fmt_val, fmt_err )
!             ==========      i             i            o        o
         write( 21,'(t8,a,' // fmt_val // ',a,' // fmt_err // ',a,a)' )  &
              'kurtosis                  :', mom1d%kurt, ' +/-', mom1d%kurt_err
      end if

   end subroutine writeMom1d_scalar
!
!=======================================================================
!
   subroutine writeMom1d_array( mom1d, title, header, nmax )
!             ================   io      i      i      i
!
!     Evaluate the moments and confidence intervals and write them to
!     the output file
!
!  Arguments:
!
      type(Mom1d_type), intent(inout) :: mom1d(:)  ! Moments array
      character(*), intent(in) :: title            ! Title text
      character(*), intent(in) :: header(:)        ! Column headers
      integer(INT4), intent(in) :: nmax            ! Maximum order of moment to be calculated
!
!  Local Variables:
!
      character(8) :: fmt_val(size(mom1d)), fmt_err(size( mom1d ))
      integer(INT4) :: i, n
!
!
!  Start of Executable Code:
!  ------------------------
!
!      if( mom1d%ncasc == 0 ) return
!
!  Evaluate moments
!
      n = size(mom1d)
      do i = 1, n
         call evalMom1d( mom1d(i) )
!             =========    io
      end do
!
!  Write moments to the output file
!
      write( 21,'(/t5,a)' ) trim( title )

      write( 21, '(t5,a)', advance='no' ) '+------------+'
      write( 21, '(100a)' ) ('------------+', i=1,n)

      write( 21, '(t5,a)', advance='no' ) '|            |'
      write( 21, '(100(a,a))' )  (header(i), '|', i=1,size( header ))

      write( 21, '(t5,a)', advance='no' ) '+------------+'
      write( 21, '(100(a))' ) ('------------+', i=1,n)
!
!  Mean value
!
      do i = 1, n
         call calcFormat( mom1d(i)%mean, mom1d(i)%mean_err, fmt_val(i), fmt_err(i) )
!             ==========        i                i              o           o
      end do
      !print *, 'mom1d(1)%mean, mom1d(1)%mean_err, fmt_val(1), fmt_err(1) = ', &
      !          mom1d(1)%mean, mom1d(1)%mean_err, fmt_val(1), fmt_err(1)

      write( 21, '(t5,a)', advance='no' ) '| mean value |'
      do i = 1, n
         write( 21, '(' // fmt_val(i) // ',a)', advance='no' ) mom1d(i)%mean, ' |'
      end do
      write( 21, '()' )

      write( 21, '(t5,a)', advance='no' ) '|        +/- |'
      do i = 1, n
         write( 21, '(' // fmt_err(i) // ',a)', advance='no' ) mom1d(i)%mean_err, ' |'
      end do
      write( 21, '()' )
!
!  Standard deviation
!
      if( nmax > 1 ) then
         do i = 1, n
            call calcFormat( mom1d(i)%std, mom1d(i)%std_err, fmt_val(i), fmt_err(i) )
!                ==========       i               i              o           o
         end do
         !print *, 'mom1d(1)%std, mom1d(1)%std_err, fmt_val(1), fmt_err(1) = ', &
         !          mom1d(1)%std, mom1d(1)%std_err, fmt_val(1), fmt_err(1)

         write( 21, '(t5,a)', advance='no' ) '| std.dev.   |'
         do i = 1, n
            write( 21, '(' // fmt_val(i) // ',a)', advance='no' ) mom1d(i)%std, ' |'
         end do
         write( 21, '()' )

         write( 21, '(t5,a)', advance='no' ) '|        +/- |'
         do i = 1, n
            write( 21, '(' // fmt_err(i) // ',a)', advance='no' ) mom1d(i)%std_err, ' |'
         end do
         write( 21, '()' )
      end if
!
!  Skewness
!
      if( .not.all( mom1d(:)%std == ZERO ) .and. nmax > 2 ) then
         do i = 1, n
            call calcFormat( mom1d(i)%skew, mom1d(i)%skew_err, fmt_val(i), fmt_err(i) )
!                ==========        i               i               o           o
         end do
         !print *, 'mom1d(1)%skew, mom1d(1)%skew_err, fmt_val(1), fmt_err(1) = ', &
         !          mom1d(1)%skew, mom1d(1)%skew_err, fmt_val(1), fmt_err(1)

         write( 21, '(t5,a)', advance='no' ) '| skewness   |'
         do i = 1, n
            write( 21, '(' // fmt_val(i) // ',a)', advance='no' ) mom1d(i)%skew, ' |'
         end do
         write( 21, '()' )

         write( 21, '(t5,a)', advance='no' ) '|        +/- |'
         do i = 1, n
            write( 21, '(' // fmt_err(i) // ',a)', advance='no' ) mom1d(i)%skew_err, ' |'
         end do
         write( 21, '()' )
      end if
!
!  Kurtosis
!
      if( .not.all( mom1d(:)%std == ZERO ) .and. nmax > 3 ) then
         do i = 1, n
            call calcFormat( mom1d(i)%kurt, mom1d(i)%kurt_err, fmt_val(i), fmt_err(i) )
!                ==========        i               i               o           o
         end do
         !print *, 'mom1d(1)%kurt, mom1d(1)%kurt_err, fmt_val(1), fmt_err(1) = ', &
         !          mom1d(1)%kurt, mom1d(1)%kurt_err, fmt_val(1), fmt_err(1)

         write( 21, '(t5,a)', advance='no' ) '| kurtosis   |'
         do i = 1, n
            write( 21, '(' // fmt_val(i) // ',a)', advance='no' ) mom1d(i)%kurt, ' |'
         end do
         write( 21, '()' )

         write( 21, '(t5,a)', advance='no' ) '|        +/- |'
         do i = 1, n
            write( 21, '(' // fmt_err(i) // ',a)', advance='no' ) mom1d(i)%kurt_err, ' |'
         end do
         write( 21, '()' )
      end if

      write( 21, '(t5,a)', advance='no' ) '+------------+'
      write( 21, '(100(a))' ) ('------------+', i=1,n)

   end subroutine writeMom1d_array
!
!=======================================================================
!
   subroutine calcFormat( value, error, fmt_val, fmt_err )
!             ==========    i      i       o        o
!
!     Calculate the format for writing values and errors to the output file
!
!  Arguments:
!
      real(REAL8), intent(in) :: value
      real(REAL8), intent(in) :: error
      character(8), intent(out) :: fmt_val
      character(8), intent(out) :: fmt_err
!
!  Local Variables:
!
      integer(INT4) :: n, v, v1, v_err
!
!
!  Start of Executable Code:
!  ------------------------
!
      if( value == ZERO ) then                  ! v1 ... digits before decimal point
         v1 = 0                                 !        (if negative: zeros after decimal
      else                                      !         point before first signifivant digit)
         v1 = ceiling(log10(abs(value)))
      end if
      v = max( 0, v1 )                          ! v ... digits before decimal point including
      if( value < ZERO )  v = v + 1             !       sign
      if( error == ZERO ) then
         n = min( 9, max( 0, 4-v1 ) )           ! n ... digits after decimal point so that
      else                                      !       error has two significant digits
         n = min( 9, max( 0, ceiling(-log10(error))+1 ) )
      end if

      if( error == ZERO ) then
         v_err = 0
      else
         v_err = max( 0, ceiling(log10(abs(error))) )
      end if

      if( n+v+1 <= 11 ) then
         write( fmt_val, '(a,i1)' )  'f11.', n
         if( v_err <= v ) then
            fmt_err = fmt_val
         else
            fmt_err = 'es11.3e2'
         end if
      else if( v < 10 ) then
         write( fmt_val, '(a,i1)' )  'f11.', 10-v
         fmt_err = 'es11.3e2'
      else
         fmt_val = 'es11.3e2'
         fmt_err = 'es11.3e2'
      end if

      !print *, value, error, fmt_val, fmt_err

      !if( value > -99999.9999d0 .and. value < 999999.9999d0) then
      !   fmt_val = 'f11.4'
      !else
      !   fmt_val = 'es11.3e2'
      !end if

      !if( error > -99999.9999d0 .and. error < 999999.9999d0) then
      !   fmt_err = 'f11.4'
      !else
      !   fmt_err = 'es11.3e2'
      !end if

   end subroutine calcFormat
!
!=======================================================================
!
end module Stat_Mom1d_mod

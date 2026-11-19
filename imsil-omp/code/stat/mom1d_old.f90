module Stat_Mom1d_Old_mod
!
!  Provide data type and manage 1d moments
!
!  The formulas for the standard deviations of the central moments is taken from
!  E.L. Lehmann: Elements of Large-Sample Theory, Springer (1999), p.298. They are also given
!  in M. Fisz, Wahrscheinlichkeitsrechnung und Mathematische Statistik, VEB Deutscher Verlag
!  der Wissnschaften (1980), p.434. It is further assumed that the uncertainty of the second
!  moment is negligible for uncertainties of the skewness and the kurtosis.
!
!  Used Modules:
!  ============
!
   use Const_mod, only : INT4, REAL8,  &
                         ZERO, ONE, TWO
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
   public :: initMom1d              ! Initialize moments
!            =========
   public :: incrMom1d              ! Increase moments
!            =========
   public :: evalMom1d              ! Evaluate mean-kurt
!            =========
   public :: writeMom1d             ! Write moments to output file
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
   type, public :: Mom1d_type                   ! Array of moments
!                  ==========      
      private
      real(REAL8) :: count                      ! number of entries
      real(REAL8) :: mom(0:8)                   ! 0th to 8th moment
      real(REAL8) :: mean                       ! mean value
      real(REAL8) :: std                        ! standard deviation
      real(REAL8) :: skew                       ! skewness
      real(REAL8) :: kurt                       ! kurtosis
      real(REAL8) :: mean_err                   ! standard deviation of mean value
      real(REAL8) :: std_err                    ! standard deviation of standard deviation
      real(REAL8) :: skew_err                   ! standard deviation of skewness
      real(REAL8) :: kurt_err                   ! standard deviation of kurtosis
      logical :: evaluated                      ! flag indicating mean-kurt have been evaluated
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
!     Initialize the moments with zero
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
      mom1d%count = ZERO

      mom1d%mom(:) = ZERO

      mom1d%evaluated = .false.

   end subroutine initMom1d
!
!=======================================================================
!
   subroutine incrMom1d( mom1d, value, weight )
!             =========    io     i      i
!
!     Increase the (unnormalized) moments
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
      mom1d%count = mom1d%count + ONE

      mom1d%mom(0) = mom1d%mom(0) + weight

      do, i=1,8
         mom1d%mom(i) = mom1d%mom(i) + weight * value**i
      end do
      
   end subroutine incrMom1d
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
      integer(INT4) :: i, j
      real(REAL8) :: kbinom(0:8,0:8)
      real(REAL8) :: sum
      real(REAL8) :: mom(0:8), cenmom(0:8), cenmom_err(2:4)
!
!  Start of Executable Code:
!  ------------------------
!
      mom1d%evaluated = .true.

      mom1d%mean = ZERO
      mom1d%std = ZERO
      mom1d%skew = ZERO
      mom1d%kurt = ZERO
      mom1d%mean_err = ZERO
      mom1d%std_err = ZERO
      mom1d%skew_err = ZERO
      mom1d%kurt_err = ZERO

      if( mom1d%count == ZERO )  return

!  Set up binomial coefficients

      do, j=0,8
         kbinom(0,j)=1
         kbinom(j,j)=1
         do, i=1,j-1
            kbinom(i,j)=kbinom(i-1,j-1)+kbinom(i,j-1)
         end do
      end do
!
!  Normalize moments (divide by sum of weights)
!
      sum = mom1d%mom(0)
      do, i=0,8
         mom(i) = mom1d%mom(i) / sum
      end do
!
!  Convert moments into central momemts
!
      do, i=0,8
         cenmom(i) = kbinom(0,i) * mom(i)
         do, j=1,i
            cenmom(i) = cenmom(i) + kbinom(j,i) * mom(i-j) * (-mom(1))**j
         end do
      end do
!
!  Calculate standard deviations of central moments (2:4)
!
      do, i=2,4
         cenmom_err(i) = sqrt( max( ZERO , cenmom(2*i) - 2*i*cenmom(i-1)*cenmom(i+1) -  &
            cenmom(i)**2 + i**2*cenmom(2)*cenmom(i-1)**2) / mom1d%count )
      end do
!
!  Convert moments into more convenient quantities
!
      mom1d%mean = mom(1)
      mom1d%mean_err = sqrt( max( ZERO, cenmom(2) ) / mom1d%count )

      if( cenmom(2) <= ZERO )  return

      mom1d%std = sqrt(cenmom(2))

      if( mom1d%std == ZERO )  return

      mom1d%std_err = cenmom_err(2) / (TWO*mom1d%std)

      mom1d%skew = cenmom(3) / mom1d%std**3
      mom1d%skew_err = cenmom_err(3) / mom1d%std**3

      mom1d%kurt = cenmom(4) / mom1d%std**4
      mom1d%kurt_err = cenmom_err(4) / mom1d%std**4

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

      if( mom1d%count == ZERO ) return

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
!             =================    io    i      i      i
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
      character(8) :: fmt_val(size( mom1d )), fmt_err(size( mom1d ))
      integer(INT4) :: i, n
!
!
!  Start of Executable Code:
!  ------------------------
!
!      if( mom1d%count == ZERO ) return
!
!  Evaluate moments
!
      n = size( mom1d )
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
      integer(INT4) :: n, v, v1
!
!
!  Start of Executable Code:
!  ------------------------
!
      if( value == ZERO ) then
         v1 = 0
      else
         v1 = ceiling(log10(abs(value)))
      end if
      v = max( 0, v1 )
      if( value < ZERO )  v = v + 1
      if( error == ZERO ) then
         n = max( 0, 4-v1 )
      else
         n = min( 9, max( 0, ceiling(-log10(error))+1 ) )
      end if

      if( n+v+1 <= 11 ) then
         write( fmt_val, '(a,i1)' )  'f11.', n
         fmt_err = fmt_val
      else if( v < 11 ) then
         write( fmt_val, '(a,i1)' )  'f11.', 11-v
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

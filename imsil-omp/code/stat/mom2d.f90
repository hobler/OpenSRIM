module Stat_Mom2d_mod
!
!  Provide data type and manage 2d moments
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
   public :: initMom2d                       ! Initialize moments
!            =========
   public :: incrMom2d                       ! Increase moments
!            =========
   public :: mergeIntoMasterMom2d            ! Merge moments of other threads into master
!            ====================
   public :: mergeMom2d                      ! Merge two Mom2d instances before evaluation
!            ==========
   public :: writeMom2d                      ! Evaluate moments and write to output file
!            ==========
   interface writeMom2d
      module procedure writeMom2d_scalar
      module procedure writeMom2d_array
   end interface
!
!  Public Data Types:
!  =================
!
   type, public :: Mom2d_type                ! Array of moments
!                  ==========
      !private   ! does not allow writing whole structure to file with single write (gfortran)
                 ! needed for writing backup file
      real(REAL8) :: nscore                  ! number of entries
!      real(REAL8) :: mom(0:2,0:2)            ! 0th to 2nd moment
      real(REAL8) :: mom(0:4,0:4)            ! 0th to 4th mixed moments
   end type Mom2d_type
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initMom2d( mom2d )
!             =========    o
!
!     Initialize the moments with zero
!
!  Arguments:
!
      type(Mom2d_type), intent(out) :: mom2d    ! moments array
!
!  Local Variables: none
!
!
!  Start of Executable Code:
!  ------------------------
!
      mom2d%nscore = ZERO

      mom2d%mom(:,:) = ZERO

   end subroutine initMom2d
!
!=======================================================================
!
   subroutine incrMom2d( mom2d, xval, yval, weight )
!             =========   io     i     i      i
!
!     Increase the (unnormalized) moments
!
!  Arguments:
!
      type(Mom2d_type), intent(inout) :: mom2d  ! moments array
      real(REAL8), intent(in) :: xval           ! x-value
      real(REAL8), intent(in) :: yval           ! y-value
      real(REAL8), intent(in) :: weight         ! weight
!
!  Local Variables:
!
      integer(INT4) :: i, j, k
!
!  Start of Executable Code:
!  ------------------------
!
      mom2d%nscore = mom2d%nscore + ONE

      mom2d%mom(0,0) = mom2d%mom(0,0) + weight

      do, i=1,4
!      do, i=1,2
         mom2d%mom(i,0) = mom2d%mom(i,0) + weight * xval**i
         mom2d%mom(0,i) = mom2d%mom(0,i) + weight * yval**i
         do, j=1,i-1
            k = i - j
            mom2d%mom(j,k) = mom2d%mom(j,k) + weight * xval**j * yval**k
         end do
      end do

   end subroutine incrMom2d
!
!=======================================================================
!
   subroutine mergeIntoMasterMom2d( mom )
!             ====================  io
!
!     Merge mom of other threads into master thread
!
!  Arguments:
!
      type(Mom2d_type), intent(inout) :: mom    ! moments to be merged
!
!  Local Variables:
!
      type(Mom2d_type), save :: mom_tmp         ! shared temporary array
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      call initMom2d( mom_tmp )
!          =========     o
!$omp end master

!$omp barrier
!$omp critical( critical_mom2d_mergeIntoMasterMom2d )
      call mergeMom2d( mom, mom_tmp )
!          ==========   i     io
!$omp end critical( critical_mom2d_mergeIntoMasterMom2d )

!$omp barrier
!$omp master
      mom = mom_tmp
!$omp end master

   end subroutine mergeIntoMasterMom2d
!
!===========================================================================================
!
   subroutine mergeMom2d( mom1, mom2 )
!             ==========   i     io
!
!     Merge mom1d1 into mom1d2
!
!  Arguments:
!
      type(Mom2d_type), intent(in) :: mom1      ! moments array
      type(Mom2d_type), intent(inout) :: mom2   ! moments array
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      mom2%nscore = mom2%nscore + mom1%nscore
      mom2%mom(:,:) = mom2%mom(:,:) + mom1%mom(:,:)

   end subroutine mergeMom2d
!
!=======================================================================
!
   subroutine evalMom2d( mom2d )
!             =========    io
!  
!     Evaluate mean value, standard deviations, skewnesses, and kurtoses
!
!  Arguments:
!
      type(Mom2d_type), intent(inout) :: mom2d  ! moments array
!
!  Local Variables:
!
      integer(INT4) :: i, j, k, l, m
      real(REAL8) :: kbinom(0:4,0:4)
!      real(REAL8) :: kbinom(0:2,0:2)
      real(REAL8) :: sum
      real(REAL8) :: cenmom(0:4,0:4)
!      real(REAL8) :: cenmom(0:2,0:2)
!
!  Start of Executable Code:
!  ------------------------
!
      if( mom2d%nscore == ZERO )  return
!
!  Set up binomial coefficients
!
      do, j=0,4
!      do, j=0,2
         kbinom(0,j) = 1
         kbinom(j,j) = 1
         do, i=1,j-1
            kbinom(i,j) = kbinom(i-1,j-1) + kbinom(i,j-1)
         end do
      end do
!
! normalize moments (divide by sum of weights)
!
      sum = mom2d%mom(0,0)
      do, j=0,4
!      do, j=0,2
         do, i=0,j
            k = j - i
            mom2d%mom(i,k) = mom2d%mom(i,k) / sum
         end do
      end do
!
!  Convert moments into central momemts
!
      do, j=0,4
!      do, j=0,2
         do, i=0,j
            k = j - i
            cenmom(i,k) = mom2d%mom(i,k)
            do, l=1,i
               cenmom(i,k) = cenmom(i,k) + kbinom(l,i) * mom2d%mom(i-l,k) * (-mom2d%mom(1,0))**l
            end do
            do, m=1,k
               cenmom(i,k) = cenmom(i,k) + kbinom(m,k) * mom2d%mom(i,k-m) * (-mom2d%mom(0,1))**m
            end do
            do, l=1,i
               do, m=1,k
                  cenmom(i,k) = cenmom(i,k) + kbinom(l,i) * kbinom(m,k) * mom2d%mom(i-l,k-m) * &
                                              (-mom2d%mom(1,0))**l * (-mom2d%mom(0,1))**m
               end do
            end do
         end do
      end do
!
!  Convert moments into more convenient quantities
!
      if( cenmom(2,0) <= ZERO .or. cenmom(0,2) <= ZERO ) return

      mom2d%mom(2,0) = sqrt( max( ZERO, cenmom(2,0) ) )
      mom2d%mom(0,2) = sqrt( max( ZERO, cenmom(0,2) ) )
      mom2d%mom(1,1) = cenmom(1,1) / (mom2d%mom(2,0) * mom2d%mom(0,2))

      do, j=3,4
         mom2d%mom(j,0) = cenmom(j,0) / mom2d%mom(2,0)**j
         mom2d%mom(0,j) = cenmom(0,j) / mom2d%mom(0,2)**j
         do i = 1, j-1
            k = j - i
            mom2d%mom(i,k) = cenmom(i,k) / (mom2d%mom(2,0)**i * mom2d%mom(0,2)**k)
         end do
      end do

   end subroutine evalMom2d
!
!=======================================================================
!
   subroutine writeMom2d_scalar( mom2d, text, units, nmax )
!             =================   io     i      i     i
!  
!     Evaluate the moments and confidence intervals and write them to
!     the output file
!
!  Arguments:
!
      type(Mom2d_type), intent(inout) :: mom2d  ! Moments array
      character(*), intent(in) :: text          ! Title text
      character(*), intent(in) :: units         ! Units
      integer(INT4), intent(in) :: nmax         ! maximum order of moments
!
!  Local Variables: none
!
!
!  Start of Executable Code:
!  ------------------------
!
      if( mom2d%nscore == ZERO )  return
!
!  Evaluate moments
!
      call evalMom2d( mom2d )
!          =========    io
!
!  Write moments to the output file
!
      write( 21,'(/t5,a)' ) trim( text )

      if( nmax > 1 )  write( 21,'(t8,a,f14.6)' )  &
                           'correlation coefficient   : ', mom2d%mom(1,1)
      if( nmax > 2 )  write( 21,'(t8,a,f12.4)' )  &
                           'skewness(12)              : ', mom2d%mom(1,2)
      if( nmax > 2 )  write( 21,'(t8,a,f12.4)' )  &
                           'skewness(21)              : ', mom2d%mom(2,1)
      if( nmax > 3 )  write( 21,'(t8,a,f12.4)' )  &
                           'kurtosis(13)              : ', mom2d%mom(1,3)
      if( nmax > 3 )  write( 21,'(t8,a,f12.4)' )  &
                           'kurtosis(22)              : ', mom2d%mom(2,2)
      if( nmax > 3 )  write( 21,'(t8,a,f12.4)' )  &
                           'kurtosis(31)              : ', mom2d%mom(3,1)

   end subroutine writeMom2d_scalar
!
!=======================================================================
!
   subroutine writeMom2d_array( mom2d, title, header, nmax )
!             ================   io      i      i      i
!
!     Evaluate the moments and confidence intervals and write them to
!     the output file
!
!  Arguments:
!
      type(Mom2d_type), intent(inout) :: mom2d(:)  ! Moments array
      character(*), intent(in) :: title            ! Title text
      character(*), intent(in) :: header(:)        ! Column headers
      integer(INT4), intent(in) :: nmax         ! maximum order of moments
!
!  Local Variables:
!
      integer(INT4) :: i, n
!
!  Start of Executable Code:
!  ------------------------
!
!      if( mom2d%nscore == ZERO )  return
!
!  Evaluate moments
!
      n = size( mom2d )
      do i = 1, n
         call evalMom2d( mom2d(i) )
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
!  Correlation coefficient
!
      if( nmax > 1 ) then
         write( 21, '(t5,a)', advance='no' ) '| corr.coef. |'
         if( minval(mom2d%mom(1,1)) > -99999.9999d0 .and.  &
             maxval(mom2d%mom(1,1)) < 999999.9999d0) then
            write( 21, '(100(f11.4,a))' ) (mom2d(i)%mom(1,1), ' |', i=1,n)
         else
            write( 21, '(100(es11.3e2,a))' ) (mom2d(i)%mom(1,1), ' |', i=1,n)
         end if
      end if
!
!  Skewnesses
!
      if( nmax > 2 ) then
         write( 21, '(t5,a)', advance='no' ) '| skewness12 |'
         if( minval(mom2d%mom(1,2)) > -99999.9999d0 .and.  &
             maxval(mom2d%mom(1,2)) < 999999.9999d0) then
            write( 21, '(100(f11.4,a))' ) (mom2d(i)%mom(1,2), ' |', i=1,n)
         else
            write( 21, '(100(es11.3e2,a))' ) (mom2d(i)%mom(1,2), ' |', i=1,n)
         end if

         write( 21, '(t5,a)', advance='no' ) '| skewness21 |'
         if( minval(mom2d%mom(2,1)) > -99999.9999d0 .and.  &
             maxval(mom2d%mom(2,1)) < 999999.9999d0) then
            write( 21, '(100(f11.4,a))' ) (mom2d(i)%mom(2,1), ' |', i=1,n)
         else
            write( 21, '(100(es11.3e2,a))' ) (mom2d(i)%mom(2,1), ' |', i=1,n)
         end if
      end if
!
!  Kurtoses
!
      if( nmax > 3 ) then
         write( 21, '(t5,a)', advance='no' ) '| kurtosis13 |'
         if( minval(mom2d%mom(1,3)) > -99999.9999d0 .and.  &
             maxval(mom2d%mom(1,3)) < 999999.9999d0) then
            write( 21, '(100(f11.4,a))' ) (mom2d(i)%mom(1,3), ' |', i=1,n)
         else
            write( 21, '(100(es11.3e2,a))' ) (mom2d(i)%mom(1,3), ' |', i=1,n)
         end if

         write( 21, '(t5,a)', advance='no' ) '| kurtosis22 |'
         if( minval(mom2d%mom(2,2)) > -99999.9999d0 .and.  &
             maxval(mom2d%mom(2,2)) < 999999.9999d0) then
            write( 21, '(100(f11.4,a))' ) (mom2d(i)%mom(2,2), ' |', i=1,n)
         else
            write( 21, '(100(es11.3e2,a))' ) (mom2d(i)%mom(2,2), ' |', i=1,n)
         end if

         write( 21, '(t5,a)', advance='no' ) '| kurtosis31 |'
         if( minval(mom2d%mom(3,1)) > -99999.9999d0 .and.  &
             maxval(mom2d%mom(3,1)) < 999999.9999d0) then
            write( 21, '(100(f11.4,a))' ) (mom2d(i)%mom(3,1), ' |', i=1,n)
         else
            write( 21, '(100(es11.3e2,a))' ) (mom2d(i)%mom(3,1), ' |', i=1,n)
         end if

      end if

      write( 21, '(t5,a)', advance='no' ) '+------------+'
      write( 21, '(100(a))' ) ('------------+', i=1,n)

   end subroutine writeMom2d_array
!
!=======================================================================
!
end module Stat_Mom2d_mod

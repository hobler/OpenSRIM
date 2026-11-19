module Stat_His1d_mod
!
!  Provides data type and manage 1d histograms. 
!  Arrays of histograms are supported. The histograms share a common grid.
!
!  Used Modules:
!  ============
!
!$ use omp_lib,               only: omp_get_thread_num
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, HALF, ONE, QUART,  &
                                    EPSMA, EPSMB, SMALL, BIG
   use Array_mod,             only: check_allocation

   use Undefine_mod,          only: defined
   use Open_mod,              only: openFile

   use Param_Setup_mod,       only: nelem, reflect, shiftto

   use Stat_His_Adapt_mod,    only: adaptHis1d
   use Stat_His_Extend_mod,   only: extendHis1d
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
   public :: initHis1d                       ! Initialize histograms
!            =========
   public :: scoreHis1d                      ! Enter score into histogram
!            ==========
   public :: getHis1d                        ! Get value at position pos
!            ========
   public :: mergeIntoMasterHis1d            ! Merge histograms of other threads into master
!            ====================
   public :: mergeHis1d                      ! Merge one histogram into another
!            ==========
   public :: broadcastHis1d                  ! Copy histogram of master thread to other threads
!            ==============
   public :: writeHis1d                      ! Output histograms
!            ==========
   interface writeHis1d
      module procedure writeHis1dnn          ! Version without factor and hmin
      module procedure writeHis1dns          ! Version without factor with scalar hmin
      module procedure writeHis1dss          ! Version with scalar factor and hmin
      module procedure writeHis1dsa          ! Version with scalar factor and array hmin
      module procedure writeHis1dan          ! Version with array factor without hmin
   end interface
   public :: readHis1d                       ! Read and allocate his
!            =========
   public :: saveHis1d                       ! Save to backup file
!            =========
   public :: loadHis1d                       ! Load from backup file
!            =========
!
!  Private Procedure:
!  =================
!
   private :: writeHis1daa                   ! Version with array factor and hmin
!             ============
   private :: sumHis1d                       ! Integrate his (currently unused)
!             ========
!
!  Public Data Types:
!  =================
!
   type, public :: His1d_type                ! histogram data type
!                  ==========
      real(REAL8), allocatable :: hist(:,:)  ! (size(indhis),nbox) histograms
      integer(INT4), allocatable :: indhis(:)! indicator of the meaning of the histograms
      integer(INT4) :: nbox                  ! number of boxes (last dim.of his)
      integer(INT4) :: imin                  ! index of first non-empty box
      integer(INT4) :: imax                  ! index of last non-empty box
      real(REAL8) :: hbeg                    ! beginning of first box
      real(REAL8) :: wbox                    ! box width
      logical :: fixedsize                   ! indicates whether nbox is fixed
   end type His1d_type
!
!==========================================================================================
!
contains
!
!==========================================================================================
!
   subroutine initHis1d( indhis, nbox, wbox, hbeg, his )
!             =========    i      i     i     i     o
!
!     Allocate memory and initialize histogram with zero. The first 
!     dimension of hist is set to the dimension of indhis. The second
!     dimension of hist is set to nbox.
!
!  Arguments:
!
      integer(INT4), intent(in) :: indhis(:)    ! indicator of histogram meaning
      integer(INT4), intent(in) :: nbox         ! histogram dimension
      real(REAL8), intent(in) :: wbox           ! initial box width
      real(REAL8), intent(in) :: hbeg           ! beginning of a histogram box
      type(His1d_type), intent(out) :: his      ! histogram structure
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: istat
!
!  Start of Executable Code:
!  ------------------------
!
      his%fixedsize = ( nbox /= 0 )

      if( allocated( his%hist ) )  deallocate( his%hist )
      allocate( his%hist( size( indhis ), nbox ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initHis1d', 'his%hist', istat, errmsg )
      his%hist(:,:) = ZERO

      if( allocated( his%indhis ) )  deallocate( his%indhis )
      allocate( his%indhis( size( indhis ) ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initHis1d', 'his%indhis', istat, errmsg )
      his%indhis(:) = indhis(:)

      his%nbox = nbox
      his%imin = nbox + 1
      his%imax = 0
      his%hbeg = hbeg
      his%wbox = wbox

   end subroutine initHis1d
!
!========================================================================================
!
   subroutine scoreHis1d( ielem, pos, his, hincr, power_ )
!             ==========    i     i   io     i      i
!
!     Enter score into the histogram, and adjust histogram grid if necessary.
!     If power_ is specified, pos**(1/power_) will be scored. This can be reversed
!     when writing the histogram, see writeHis1d. power_ must not be specified
!     together with negative pos.
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem           ! index of his(ielem,ibox)
      real(REAL8), intent(in) :: pos               ! position determining ibox
      type(His1d_type), intent(inout) :: his       ! histogram structure
      real(REAL8), intent(in) :: hincr             ! increment of his(ielem,ibox) times the
                                                   !    box width
      integer(INT4), intent(in), optional :: power_! power to take of the abscissa values
!
!  Local Variables:
!
      integer(INT4) :: ibox
      real(REAL8) :: pos1
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( power_ ) ) then
         if( power_ == 1 ) then
            pos1 = pos
         else
            pos1 = pos ** (ONE/power_)
         end if
      else
         pos1 = pos
      end if
!
!  Determine box index (ibox) and adjust histogram grid if necessary
!
      if( his%fixedsize ) then

         call adaptHis1d( pos1, ibox, his%hist, his%hbeg, his%wbox, his%imin, his%imax )
!             ==========   i     o       io        io        io        io        io

      else

         call extendHis1d( pos1, ibox, his%hist, his%hbeg, his%wbox, his%imin, his%imax )
!             ===========   i     o       io        io        i         io        io
         his%nbox = size( his%hist, 2 )

      end if
!
!  Score value in histogram
!
      his%hist(ielem,ibox) = his%hist(ielem,ibox) + hincr / his%wbox

   end subroutine scoreHis1d
!
!===========================================================================================
!
   subroutine getHis1d( ielem, pos, his, hval )
!             ========    i     i    i    o
!
!     Get the histogram value at a specified position
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem        ! index of his(ielem,ibox)
      real(REAL8), intent(in) :: pos            ! position determining ibox
      type(His1d_type), intent(in) :: his       ! histogram structure
      real(REAL8), intent(out) :: hval          ! value of his(ielem,ibox)
!
!  Local Variables:
!
      integer(INT4) :: ibox
!
!  Start of Executable Code:
!  ------------------------
!
      ibox = ceiling( (pos - his%hbeg) / his%wbox )

      if( ibox < his%imin .or. ibox > his%imax ) then
         hval = ZERO
      else
         hval = his%hist(ielem,ibox)
      end if

   end subroutine getHis1d
!
!=======================================================================
!
   subroutine mergeIntoMasterHis1d( his )
!             ====================  io
!
!     Merge his of other threads into master thread, while leaving his of other threads
!     untouched
!
!     This subroutine must be called from a parallel region
!
!  Arguments:
!
      type(His1d_type), intent(inout) :: his    ! histogram to be merged
!
!  Local Variables:
!
      integer(INT4) :: ithread
      type(His1d_type), save :: his_tmp         ! shared temporary histogram
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      his_tmp = his
!$omp end master

!$omp barrier
!$omp critical( critical_his1d_mergeIntoMasterHis1d )
      ithread = 1
!$    ithread = 1 + omp_get_thread_num( )
      if( ithread > 1 ) then
         call mergeHis1d( his, his_tmp )
!             ==========   i     io
      end if
!$omp end critical( critical_his1d_mergeIntoMasterHis1d )

!$omp barrier
!$omp master
      his = his_tmp
!$omp end master

   end subroutine mergeIntoMasterHis1d
!
!===========================================================================================
!
   subroutine mergeHis1d( his1, his2, ielem1_, ielem2_ )
!             ==========   i     io      i        i
!
!     Merge his1 into his2 (add his1 to his2)
!
!     If ielem1_ and ielem2_ are present, only his2(ielem2_,:) is merged into his2(ielem1_,:),
!     otherwise all histograms are merged.
!
!  Arguments:
!
      type(His1d_type), intent(in) :: his1            ! histogram to be merged
      type(His1d_type), intent(inout) :: his2         ! histogram to be merged into
      integer(INT4), intent(in), optional :: ielem1_  ! atom index in his1
      integer(INT4), intent(in), optional :: ielem2_  ! atom index in his2
!
!  Local Variables:
!
      integer(INT4) :: i1, i2, imin1, imin2, imax1, imax2
      integer(INT4) :: ielem1min, ielem1max, ielem2min, ielem2max
      real(REAL8) :: eps, xold
      real(REAL8), allocatable :: x1(:), x2(:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Nothing to do if his1 is empty
!
      if( his1%imin > his1%imax )  return
!
!  Box boundaries
!
      imin1 = his1%imin - 1
      imax1 = his1%imax
      allocate( x1( imin1:imax1 ) )
      x1(:) = [ (his1%hbeg + i1*his1%wbox, i1=imin1,imax1) ]
!
!  Make sure his1 completely fits into his2
!
      eps = max( EPSMA*abs(x1(imin1)), SMALL )
      call scoreHis1d( 1, x1(imin1)+eps, his2, ZERO )
!          ==========  i         i        io    i

      eps = max( EPSMA*abs(x1(imax1)), SMALL )
      call scoreHis1d( 1, x1(imax1)-eps, his2, ZERO )
!          ==========  i         i        io    i
      imin2 = his2%imin - 1
      imax2 = his2%imax
      allocate( x2( imin2:imax2 ) )
      x2(:) = [ (his2%hbeg + i2*his2%wbox, i2=imin2,imax2) ]
!
!  Merge his1 into his2
!
      if( present( ielem1_ ) ) then
         ielem1min = ielem1_
         ielem1max = ielem1_
         ielem2min = ielem2_
         ielem2max = ielem2_
      else
         ielem1min = 1
         ielem1max = size( his1%hist, 1 )
         ielem2min = 1
         ielem2max = size( his2%hist, 1 )
      end if

      xold = x1(imin1)
      i1 = imin1 + 1
      do i2 = imin2+1, imax2
         if( x2(i2) > xold )  exit
      end do
      i2 = min( i2, imax2 )

      do
         if( x1(i1) < x2(i2) ) then
            his2%hist(ielem2min:ielem2max,i2) = his2%hist(ielem2min:ielem2max,i2)  &
                  + his1%hist(ielem1min:ielem1max,i1) * (x1(i1)-xold) / his2%wbox
            xold = x1(i1)
            i1 = i1 + 1
            if( i1 > imax1 )  exit
         else
            his2%hist(ielem2min:ielem2max,i2) = his2%hist(ielem2min:ielem2max,i2)  &
                  + his1%hist(ielem1min:ielem1max,i1) * (x2(i2)-xold) / his2%wbox
            xold = x2(i2)
            i2 = i2 + 1
            if( i2 > imax2 )  exit
         end if
      end do

   end subroutine mergeHis1d
!
!===========================================================================================
!
   subroutine broadcastHis1d( his )
!             ==============  io
!
!     Spread his of master thread to other threads
!
!     This subroutine must be called from within a parallel region
!
!  Arguments:
!
      type(His1d_type), intent(inout) :: his    ! histogram to be spread
!
!  Local Variables:
!
      type(His1d_type), save :: his_tmp         ! shared temporary histogram
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      his_tmp = his
!$omp end master
!$omp barrier
      his = his_tmp

   end subroutine broadcastHis1d
!
!===========================================================================================
!
   subroutine writeHis1dnn( fname, fmode, his, lsum_, power_ )
!             ============    i      i     i     i      i
!
!     Writes the histogram to a specified file. This is the version without factor and
!     cutoff values.
!
!  Arguments:
!
      character(*), intent(in) :: fname            ! File name
      character(*), intent(in) :: fmode            ! indicates whether output should be
                                                   !    appended or overwrite file
      type(His1d_type), intent(in) :: his          ! Histogram structure
      logical, intent(in), optional :: lsum_       ! indicates whether the sum of histogram
                                                   ! values should be written to the file
      integer(INT4), intent(in), optional :: power_! power to take of the abscissa values
!
!  Local Variables:
!
      real(REAL8) :: factor(size( his%hist, 1 ))
      real(REAL8) :: hmin(size( his%hist, 1 ))
      real(REAL8) :: hmax(size( his%hist, 1 ))
      real(REAL8) :: hcrit(size( his%hist, 1 ))
!
!  Start of Executable Code:
!  ------------------------
!
      factor(:) = ONE
      hmin(:) = -BIG
      hmax(:) = BIG
      hcrit(:) = BIG

      call writeHis1daa( fname, fmode, factor, his, hmin, hmax, hcrit, lsum_, power_ )
!          ============    i      i      i      i    i     i      i      i      i

   end subroutine writeHis1dnn
!
!===========================================================================================
!
   subroutine writeHis1dns( fname, fmode, his, hmin, lsum_, power_ )
!             ============    i      i     i    i      i      i
!
!     Writes the histogram to a specified file. This is the version without factor and
!     with scalar cutoff values.
!
!  Arguments:
!
      character(*), intent(in) :: fname            ! File name
      character(*), intent(in) :: fmode            ! indicates whether output should be
                                                   !    appended or overwrite file
      type(His1d_type), intent(in) :: his          ! Histogram structure
      real(REAL8), intent(in) :: hmin              ! Lower cut-off value times wbox
      logical, intent(in), optional :: lsum_       ! indicates whether the sum of histogram
                                                   ! values should be written to the file
      integer(INT4), intent(in), optional :: power_! power to take of the abscissa values
!
!  Local Variables:
!
      real(REAL8) :: factor(size( his%hist, 1 ))
      real(REAL8) :: hmin_(size( his%hist, 1 ))
      real(REAL8) :: hmax(size( his%hist, 1 ))
      real(REAL8) :: hcrit(size( his%hist, 1 ))
!
!  Start of Executable Code:
!  ------------------------
!
      factor(:) = ONE
      hmin_(:) = hmin
      hmax(:) = BIG
      hcrit(:) = BIG

      call writeHis1daa( fname, fmode, factor, his, hmin_, hmax, hcrit, lsum_, power_ )
!          ============    i      i      i      i     i      i       i     i      i

   end subroutine writeHis1dns
!
!===========================================================================================
!
   subroutine writeHis1dss( fname, fmode, factor, his, hmin_, hmax_, hcrit_, lsum_, power_ )
!             ============    i      i      i      i     i     i       i       i      i
!
!     Writes the histogram to a specified file. This is the version with scalar factor and
!     with scalar cutoff values.
!
!  Arguments:
!
      character(*), intent(in) :: fname            ! File name
      character(*), intent(in) :: fmode            ! indicates whether output should be
                                                   !    appended or overwrite file
      real(REAL8), intent(in) :: factor            ! factor to multiply the histogram with
      type(His1d_type), intent(in) :: his          ! Histogram structure
      real(REAL8), intent(in) :: hmin_             ! Lower cut-off value times wbox
      real(REAL8), intent(in), optional :: hmax_   ! Upper cut-off value
      real(REAL8), intent(in), optional :: hcrit_  ! Critical value
      logical, intent(in), optional :: lsum_       ! indicates whether the sum of histogram
                                                   ! values should be written to the file
      integer(INT4), intent(in), optional :: power_! power to take of the abscissa values
!
!  Local Variables:
!
      real(REAL8) :: factor_(size( his%hist, 1 ))
      real(REAL8) :: hmin(size( his%hist, 1 ))
      real(REAL8) :: hmax(size( his%hist, 1 ))
      real(REAL8) :: hcrit(size( his%hist, 1 ))
!
!  Start of Executable Code:
!  ------------------------
!
      factor_(:) = factor
      hmin(:) = hmin_
      if( present( hmax_ ) ) then
         hmax(:) = hmax_
      else
         hmax(:) = BIG
      end if
      if( present( hcrit_ ) ) then
         hcrit(:) = hcrit_
      else
         hcrit(:) = BIG
      end if

      call writeHis1daa( fname, fmode, factor_, his, hmin, hmax, hcrit, lsum_, power_ )
!          ============    i      i       i      i    i     i      i      i      i

   end subroutine writeHis1dss
!
!===========================================================================================
!
   subroutine writeHis1dsa( fname, fmode, factor, his, hmin, hmax_, hcrit_, lsum_, power_ )
!             ============    i      i      i      i    i      i      i       i      i
!
!     Writes the histogram to a specified file. This is the version with scalar factor and
!     with array cutoff values.
!
!  Arguments:
!
      character(*), intent(in) :: fname               ! File name
      character(*), intent(in) :: fmode               ! indicates whether output should be
                                                      !    appended or overwrite file
      real(REAL8), intent(in) :: factor               ! factor to multiply the histogram with
      type(His1d_type), intent(in) :: his             ! Histogram structure
      real(REAL8), intent(in) :: hmin(:)              ! Lower cut-off value times wbox
      real(REAL8), intent(in), optional :: hmax_(:)   ! Upper cut-off value
      real(REAL8), intent(in), optional :: hcrit_(:)  ! Critical value
      logical, intent(in), optional :: lsum_          ! indicates whether the sum of histogram
                                                      ! values should be written to the file
      integer(INT4), intent(in), optional :: power_   ! power to take of the abscissa values
!
!  Local Variables:
!
      real(REAL8) :: factor_(size( his%hist, 1 ))
      real(REAL8) ::    hmax(size( his%hist, 1 ))
      real(REAL8) ::   hcrit(size( his%hist, 1 ))
!
!  Start of Executable Code:
!  ------------------------
!
      factor_(:) = factor
      if( present( hmax_ ) ) then
         hmax(:) = hmax_(:)
      else
         hmax(:) = BIG
      end if
      if( present( hcrit_ ) ) then
         hcrit(:) = hcrit_(:)
      else
         hcrit(:) = BIG
      end if

      call writeHis1daa( fname, fmode, factor_, his, hmin, hmax, hcrit, lsum_, power_ )
!          ============    i      i       i      i    i     i      i      i      i

   end subroutine writeHis1dsa
!
!===========================================================================================
!
   subroutine writeHis1dan( fname, fmode, factor, his, lsum_, power_ )
!             ============    i      i      i      i     i      i
!
!     Writes the histogram to a specified file. This is the version with scalar factor and
!     with scalar cutoff values.
!
!  Arguments:
!
      character(*), intent(in) :: fname            ! File name
      character(*), intent(in) :: fmode            ! indicates whether output should be
                                                   !    appended or overwrite file
      real(REAL8), intent(in) :: factor(:)         ! factor to multiply the histogram with
      type(His1d_type), intent(in) :: his          ! Histogram structure
      logical, intent(in), optional :: lsum_       ! indicates whether the sum of histogram
                                                   ! values should be written to the file
      integer(INT4), intent(in), optional :: power_! power to take of the abscissa values
!
!  Local Variables:
!
      real(REAL8) ::  hmin_(size( his%hist, 1 ))
      real(REAL8) ::  hmax_(size( his%hist, 1 ))
      real(REAL8) :: hcrit_(size( his%hist, 1 ))
!
!  Start of Executable Code:
!  ------------------------
!
      hmin_(:) = - BIG
      hmax_(:) = BIG
      hcrit_(:) = BIG

      call writeHis1daa( fname, fmode, factor, his, hmin_, hmax_, hcrit_, lsum_, power_ )
!          ============    i      i      i      i     i      i      i       i      i

   end subroutine writeHis1dan
!
!===========================================================================================
!
   subroutine writeHis1daa( fname, fmode, factor, his, hmin, hmax, hcrit, lsum_, power_ )
!             ============    i      i      i      i    i     i      i      i      i
!  
!     Writes the histogram to a specified file. Histogram values are limited 
!     to the range [hmin,hmax]. Histogram values larger than hcrit are set
!     to hmax. This is the version with array factor and array cutoff values.
!
!  Arguments:
!
      character(*), intent(in) :: fname            ! File name
      character(*), intent(in) :: fmode            ! indicates whether output should be
                                                   !    appended or overwrite file
      real(REAL8), intent(in) :: factor(:)         ! factors to multiply the histogram with
      type(His1d_type), intent(in) :: his          ! Histogram structure
      real(REAL8), intent(in) :: hmin(:)           ! Lower cut-off value times wbox
      real(REAL8), intent(in) :: hmax(:)           ! Upper cut-off value
      real(REAL8), intent(in) :: hcrit(:)          ! Critical value (see Purpose)
      logical, intent(in), optional :: lsum_        ! indicates whether the sum of histogram
                                                   ! values should be written to the file
      integer(INT4), intent(in), optional :: power_ ! power to take of the abscissa values
!
!  Local Variables:
!
      character(20) :: fmt
      logical :: lsum
      integer(INT4) :: power
      integer(INT4) :: na, ibox, nx, v
      real(REAL8) :: xmax
      real(REAL8), allocatable :: x(:), hist(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( lsum_ ) ) then
         lsum = lsum_
      else
         lsum = .false.
      end if

      if( present( power_ ) ) then
         power = power_
      else
         power = 1
      end if
!
!  Do not write file if histogram is empty
!
      if( his%imin > his%imax ) return
!
!  Copy to more convenient arrays
!
      allocate( x( his%imax-his%imin+2 ) )
      x(:) = [ ( his%hbeg+ibox*his%wbox, ibox=his%imin-1,his%imax ) ]
      if( power /= 1 )  x(:) = x(:)**power
      allocate( hist( size( his%hist, 1 ), his%imax-his%imin+1 ) )
      hist(:,:) = his%hist(:,his%imin:his%imax)
      hist(:,:) = hist(:,:) * his%wbox       ! convert to integral content

      do ibox = 1, size( hist, 2 )
!
!  Cut off values lower than hmin
!
         hist(:,ibox) = max( hist(:,ibox), hmin(:) )
!
!  Convert to density
!
         hist(:,ibox) = hist(:,ibox) / (x(ibox+1)-x(ibox))
!
!  Set damage concentration to hmax if it exceeds hcrit
!
         where( hist(:,ibox) > hcrit(:) )  hist(:,ibox) = hmax(:)

      end do
!
!  Open file
!
      call openFile( 22, fname, 'formatted', fmode )
!          ========  i     i         i         i
!
!  Write header
!
      na = size( hist, 1 )
      nx = size( x )

      if( lsum ) then
         write( 22, '( i7, i13 )' ) 2, 2*nx
         write( 22, '( i7, i13 )' ) 0, 0
      else
         write( 22, '( i7, 100i14 )' ) 1+na, 2*nx
         write( 22, '( i7, 100i14 )' ) 0, his%indhis(:)
      end if
!
!  Format for writing the lines
!
      fmt = '(f'
      xmax = maxval( x(:) )
      if( any( x(:) < ZERO ) )  xmax = max( xmax, 10*maxval( abs(x(:)), mask=(x(:)<ZERO) ) )
      v = max( 2, ceiling(log10(xmax))+1 )
      write( fmt(3:), '(i2)' )  v+4
      fmt = trim(fmt) // '.3, '
      if( .not.lsum )  fmt = trim(fmt) // '100'

      if( any( hist(:,:) < ZERO ) .or. any( hmin(:) < ZERO ) ) then
         fmt = trim(fmt) // 'es15.6'
      else
         fmt = trim(fmt) // 'es14.6'
      end if
      if( any( abs( hist(:,:) ) <= 1d-100 .and. hist(:,:) /= ZERO ) .or.  &
          any( abs( hmin(:) ) <= 1d-100 .and. hmin(:) /= ZERO ) ) then
         fmt = trim(fmt) // 'e3'
      end if
      fmt = trim(fmt) // ')'
!
!  Write first point
!
      if( lsum ) then
         write( 22, fmt )  x(1), sum( factor(:)*hmin(:)/his%wbox )
      else
         write( 22, fmt )  x(1), factor(:)*hmin(:)/his%wbox
      end if
!
!  Write histogram values
!
      do ibox = 1, size( hist(:,:), 2 )
         if( lsum ) then
            write( 22, fmt ) x(ibox), sum( factor(:)*hist(:,ibox) )
            write( 22, fmt ) x(ibox+1), sum( factor(:)*hist(:,ibox) )
         else
            write( 22, fmt ) x(ibox), factor(:)*hist(:,ibox)
            write( 22, fmt ) x(ibox+1), factor(:)*hist(:,ibox)
         end if
      enddo
!
!  Write last point
!
      if( lsum ) then
         write( 22, fmt )  x(size( x )), sum( factor(:)*hmin(:)/his%wbox )
      else
         write( 22, fmt )  x(size( x )), factor(:)*hmin(:)/his%wbox
      end if
!
!  Deallocate auxiliary array and close the file
!
      deallocate( x, hist )

      close( 22 )

   end subroutine writeHis1daa
!
!==========================================================================================
!
   subroutine readHis1d( fname, his, ielem1_, ielem2_, factor_ )
!             =========    i    io     i       i       i
!  
!     Reads the histogram from the specifiled file and adds it to his. 
!
!     Note: his%hist will be noisy unless the original box width of his is
!           2^n times (n=0,1,2,...) the histogram box width on the file
!
!  Arguments:
!
      character(*), intent(in) :: fname               ! File name
      type(His1d_type), intent(inout) :: his          ! Histogram structure
      integer(INT4), intent(in), optional :: ielem1_  ! his file column number
      integer(INT4), intent(in), optional :: ielem2_  ! atom index in his (1,..nelem)
      real(REAL8), intent(in), optional :: factor_    ! factor to multiply the histogram with
!
!  Local Variables:
!
      integer(INT4) :: ncol, nrow, nelem1, idummy, i, ibox
      real(REAL8) :: z1, z2, factor
      type(His1d_type) :: his1
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( factor_ ) ) then
         factor = factor_
      else
         factor = ONE
      end if
!
!  Open file
!
      call openFile( 22, fname, 'formatted', 'read' )
!          ========  i     i         i          i
!
!  Read first line and allocate his1%hist
!
      read( 22, *, err=900, end=900 )  ncol, nrow

      nelem1 = ncol - 1
      his1%nbox = ( nrow - 2 ) / 2
      his1%imin = 1
      his1%imax = his1%nbox
      allocate( his1%hist( nelem1, his1%nbox ) )

      if( present( ielem1_ ) ) then
         if( ielem1_ > nelem1 )  error stop 'readHis1d: requested file column number too big'
         if( .not.present( ielem2_ ) )  &
               error stop 'readHis1d: either none or both of ielem1, ielem2 must be present'
      else
         if( nelem1 /= nelem )  error stop 'readHis1d: inconsistent number of columns on file'
      end if
!
!  Second and third line can be ignored
!
      read( 22, *, err=900, end=900 )  ( idummy, i=1,ncol )
      read( 22, *, err=900, end=900 )  his1%hbeg, his1%hist(:,1)
!
!  Read histogram values
!
      do ibox = 1, his1%nbox
         read( 22, *, err=900, end=900 ) z1, his1%hist(:,ibox)  ! for a legal histogram file,
         read( 22, *, err=900, end=900 ) z2, his1%hist(:,ibox)  ! hisf is same at both points

         if( ibox == 1 ) then
            his1%wbox = z2 - z1
         else
            if( abs( z2-z1 - his1%wbox ) > EPSMB * abs( his1%wbox ) )  &
                  error stop 'readHis1d: Histogram boxes must be same size'
         end if
      end do

      his1%hist(:,:) = factor * his1%hist(:,:)
!
!  Shift and/or reflect histogram if desired
!
      if( defined( shiftto ) )  his1%hbeg = shiftto
      if( reflect ) then
         his1%hbeg = - (his1%hbeg + his1%nbox*his1%wbox)
         his1%hist(:,:) = his1%hist(:,::-1)
      end if
!
!  Merge read histogram into his
!
      call mergeHis1d( his1, his, ielem1_, ielem2_ )
!          ==========   i    io     i       i

      return

900   error stop 'readHis1d: error during reading histogram file'

   end subroutine readHis1d
!
!=======================================================================
!
   subroutine sumHis1d( ielem, pos1, pos2, his, sval )
!             ========    i     i     i     i    o
!  
!     Integrate the histogram within given bounds
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem        ! index of his(ielem,ibox)
      real(REAL8), intent(in) :: pos1           ! position determining lower point
      real(REAL8), intent(in) :: pos2           ! position determining upper point
      type(His1d_type), intent(in) :: his       ! histogram structure
      real(REAL8), intent(out) :: sval          ! Int((pos1,pos2){his(ielem,x)}dx) [cm^-3]
!
!  Local Variables:
!
      integer(INT4) :: imin, imax
      real(REAL8) :: dr1, dr2, val1, val2
!
!  Start of Executable Code:
!  ------------------------
!
      if( ielem > size( his%indhis, 1 ) .or. ielem < 1 )  &
            error stop 'sumHis1d: index out of bounds'

      if( pos1 >= pos2 ) then
         sval = ZERO
         return
      end if

      imin = nint( ( pos1 - his%hbeg ) / his%wbox + HALF )
      imax = nint( ( pos2 - his%hbeg ) / his%wbox + HALF )
                                                ! shring indices if out of his bounds
      imin = max( his%imin - 1, imin )
      imax = min( his%imax + 1, imax )

      call getHis1d( ielem, pos1, his, val1 )
!          ========    i     i     i    o
      call getHis1d( ielem, pos2, his, val2 )
!          ========    i     i     i    o
!
!  Count lower/upper part of depth which may not be completely covered by imin/imax
!
      dr1 = imin * his%wbox - pos1
      dr2 = pos2 - ( imax - 1 ) * his%wbox

      sval = ZERO

      if( imin == imax ) then
         sval = ( dr2 + dr1 - his%wbox ) * val1
      else 
         sval = sum( his%hist(ielem,imin + 1: imax - 1) ) * his%wbox
         sval = sval + dr1 * val1 + dr2 * val2  
      end if

   end subroutine sumHis1d
!
!==========================================================================================
!
   subroutine saveHis1d( his )
!             =========   i
!  
!     Saves the his structure to the file with unit number 22
!
!  Arguments:
!
      type(His1d_type), intent(in) :: his       ! histogram structure
!
!  Local Variables:
!
      integer(INT4) :: nhis
!
!  Start of Executable Code:
!  ------------------------
!
      nhis = size( his%indhis(:) )
      write( 22 )  nhis, his%nbox
      write( 22 )  his%fixedsize
      write( 22 )  his%hbeg, his%wbox, his%imin, his%imax
      write( 22 )  his%hist(:,:), his%indhis(:)

   end subroutine saveHis1d
!
!=======================================================================
!
   subroutine loadHis1d( his )
!             =========  io
!  
!     Load the his structure to the file with unit number 22
!
!  Arguments:
!
      type(His1d_type), intent(inout) :: his    ! histogram structure
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: nhis
      integer(INT4) :: istat
!
!  Start of Executable Code:
!  ------------------------
!
      read( 22, iostat=istat )  nhis, his%nbox
      read( 22, iostat=istat )  his%fixedsize
      read( 22, iostat=istat )  his%hbeg, his%wbox, his%imin, his%imax

      if( allocated( his%hist ) ) then
         deallocate( his%hist )
         deallocate( his%indhis )
      end if

      allocate( his%hist( nhis, his%nbox ), stat=istat, errmsg=errmsg )
      call check_allocation( 'loadHis1d', 'his%hist', istat, errmsg )

      allocate( his%indhis( nhis ), stat=istat, errmsg=errmsg )
      call check_allocation( 'loadHis1d', 'his%indhis', istat, errmsg )

      read( 22, iostat=istat )  his%hist(:,:), his%indhis(:)
      if( istat /= 0 )  error stop 'error during reading backup file'      
     
   end subroutine loadHis1d
!
!=======================================================================
!
end module Stat_His1d_mod

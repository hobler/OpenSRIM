module Stat_His3d_mod
!
!  Provides data type and manage 3d histograms. 
!  Arrays of histograms are supported. The histograms share a common grid.
!
!  Used Modules:
!  ============
!
!$ use omp_lib,               only: omp_get_thread_num
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, HALF, ONE,  &
                                    EPSMA, SMALL, BIG
   use Array_mod,             only: check_allocation

   use Open_mod,              only: openFile
   use Stat_His_adapt_mod,    only: adaptHis3d
   use Stat_His_Extend_mod,   only: extendHis3d
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
   public :: initHis3d                       ! Initialize histograms
!            =========
   public :: scoreHis3d                      ! Enter score into histogram
!            ==========
   public :: getHis3d                        ! Get value at position pos
!            ========
   public :: mergeIntoMasterHis3d            ! Merge histograms of other threads into master
!            ====================
   public :: mergeHis3d                      ! Merge one histogram into another
!            ==========
   public :: broadcastHis3d                  ! Copy histogram of master thread to other threads
!            ==============
   public :: writeHis3d                      ! Output histograms
!            ==========
   interface writeHis3d
      module procedure writeHis3dScalar
      module procedure writeHis3dArray
   end interface
   public :: saveHis3d                       ! Save to backup file
!            =========
   public :: loadHis3d                       ! Load from backup file
!            =========
!
!  Public Data Types:
!  =================
!
   type, public :: His3d_type                   ! histogram data type
!                  ==========
      real(REAL8), allocatable :: hist(:,:,:,:) ! (size(indhis),nbox(1),nbox(2),nbox(3))
                                                !   contains the histograms
      integer(INT4), allocatable :: indhis(:)   ! indicator of the meaning of the histograms
      integer(INT4) :: nbox(3)                  ! number of boxes (last dim.of his)
      integer(INT4) :: imin(3)                  ! index of first non-empty box
      integer(INT4) :: imax(3)                  ! index of last non-empty box
      real(REAL8) :: hbeg(3)                    ! beginning of first box
      real(REAL8) :: wbox(3)                    ! box width
      logical :: fixedsize                      ! indicates whether nbox is fixed
   end type His3d_type
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initHis3d( indhis, nbox, wbox, hbeg, his)
!             =========    i      i     i     i     o
!
!     Allocate memory and initialize histogram with zero. The first 
!     dimension of hist is set to the dimensions of indhis. The second
!     to fourth dimension of hist are set to the components of nbox.
!
!  Arguments:
!
      integer(INT4), intent(in) :: indhis(:)    ! indicator of histogram meaning
      integer(INT4), intent(in) :: nbox(3)      ! histogram dimension
      real(REAL8), intent(in) :: wbox(3)        ! initial box width
      real(REAL8), intent(in) :: hbeg(3)        ! corner of a histogram box
      type(His3d_type), intent(out) :: his      ! histogram structure
!
!  Local Variables:
!
      character(256) :: errmsg
      integer(INT4) :: istat
!
!  Start of Executable Code:
!  ------------------------
!
      his%fixedsize = all( nbox(:) /= 0 )

      if( allocated( his%hist ) )  deallocate( his%hist )
      allocate( his%hist( size( indhis, 1 ), nbox(1), nbox(2), nbox(3) ), stat=istat,  &
                errmsg=errmsg )
      call check_allocation( 'initHis3d', 'his%hist', istat, errmsg )
      his%hist(:,:,:,:) = ZERO

      if( allocated( his%indhis ) )  deallocate( his%indhis )
      allocate( his%indhis( size( indhis, 1 ) ), stat=istat, errmsg=errmsg )
      call check_allocation( 'initHis3d', 'his%indhis', istat, errmsg )
      his%indhis(:) = indhis(:)

      his%nbox(:) = nbox(:)
      his%imin(:) = nbox(:) + 1
      his%imax(:) = 0
      his%hbeg(:) = hbeg(:)
      his%wbox(:) = wbox(:)

   end subroutine initHis3d
!
!=======================================================================
!
   subroutine scoreHis3d( ia, pos, his, hincr )
!             ==========   i   i   io     i
!
!     Adjust histogram grid if necessary and enter score into the histogram
!
!  Arguments:
!
      integer(INT4), intent(in) :: ia           ! index of his(ia,ibox)
      real(REAL8), intent(in) :: pos(3)         ! position determining ibox
      type(His3d_type), intent(inout) :: his    ! histogram structure
      real(REAL8), intent(in) :: hincr          ! increment of his(ia,ibox) times the
!                                               !    product of the box widths
!  Local Variables:
!
      integer(INT4) :: ibox(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  Adjust histogram grid if necessary and determine box index (ibox)
!
      if( his%fixedsize ) then

         call adaptHis3d( pos(:), ibox(:), &
!             ==========    i       o
                          his%hist, his%hbeg(:), his%wbox(:), his%imin(:), his%imax(:) )
!                           io         io           io           io           io

      else

         call extendHis3d( pos(:), ibox(:),  &
!             ===========    i       o
                           his%hist, his%hbeg(:), his%wbox(:), his%imin(:), his%imax(:) )
!                            io         io           i            io           io
         his%nbox(:) = [ size( his%hist, 2 ), size( his%hist, 3 ), size( his%hist, 4 ) ]

      end if
!
!  Score value in histogram
!
      his%hist(ia,ibox(1),ibox(2),ibox(3)) = his%hist(ia,ibox(1),ibox(2),ibox(3))  &
            + hincr / ( his%wbox(1) * his%wbox(2) * his%wbox(3) )

   end subroutine scoreHis3d
!
!=======================================================================
!
   subroutine getHis3d( ia, pos, his, hval )
!             ========   i   i    i    o
!
!     Get the histogram value at a specified position
!
!  Arguments:
!
      integer(INT4), intent(in) :: ia           ! index of his(ia,iboxx,iboxz)
      real(REAL8), intent(in) :: pos(3)         ! position determining ibox
      type(His3d_type), intent(in) :: his       ! histogram structure
      real(REAL8), intent(out) :: hval          ! value of his(ia,iboxx,iboxz)
!
!  Local Variables:
!
      integer(INT4) :: iboxx, iboxy, iboxz
!
!  Start of Executable Code:
!  ------------------------
!
      iboxx = nint( ( pos(1)-his%hbeg(1) ) / his%wbox(1)+HALF )
      iboxy = nint( ( pos(2)-his%hbeg(2) ) / his%wbox(2)+HALF )
      iboxz = nint( ( pos(3)-his%hbeg(3) ) / his%wbox(3)+HALF )

      if( iboxx < his%imin(1) .or. iboxx > his%imax(1) .or. &
          iboxy < his%imin(2) .or. iboxy > his%imax(2) .or. &
          iboxz < his%imin(3) .or. iboxz > his%imax(3) ) then
         hval = ZERO
      else
         hval = his%hist( ia, iboxx, iboxy, iboxz )
      end if

   end subroutine getHis3d
!
!=======================================================================
!
   subroutine mergeIntoMasterHis3d( his )
!             ====================  io
!
!     Merge his of other threads into master thread, while leaving his of other threads
!     untouched
!
!     This subroutine must be called from a parallel region
!
!  Arguments:
!
      type(His3d_type), intent(inout) :: his    ! histogram to be merged
!
!  Local Variables:
!
      integer(INT4) :: ithread
      type(His3d_type), save :: his_tmp         ! shared temporary histogram
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      his_tmp = his
!$omp end master

!$omp barrier
!$omp critical( critical_his3d_mergeIntoMasterHis3d )
      ithread = 1
!$    ithread = 1 + omp_get_thread_num( )
      if( ithread > 1 ) then
         call mergeHis3d( his, his_tmp )
!             ==========   i     io
      end if
!$omp end critical( critical_his3d_mergeIntoMasterHis3d )

!$omp barrier
!$omp master
      his = his_tmp
!$omp end master

   end subroutine mergeIntoMasterHis3d
!
!===========================================================================================
!
   subroutine mergeHis3d( his1, his2 )
!             ==========   i     io
!
!     Merge his1 into his2 (add his1 to his2)
!
!  Arguments:
!
      type(His3d_type), intent(in) :: his1      ! histogram to be merged from
      type(His3d_type), intent(inout) :: his2   ! histogram to be merged into
!
!  Local Variables:
!
      integer(INT4) :: ix1, ix2, iy1, iy2, iz1, iz2
      integer(INT4) :: ixmin1, ixmin2, iymin1, iymin2, izmin1, izmin2
      integer(INT4) :: ixmax1, ixmax2, iymax1, iymax2, izmax1, izmax2
      real(REAL8) :: epsx, epsy, epsz, xold, yold, zold, dx, dy, dz
      real(REAL8), allocatable :: x1(:), x2(:), y1(:), y2(:), z1(:), z2(:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Nothing to do if his1 is empty
!
      if( any( his1%imin(:) > his1%imax(:) ) )  return
!
!  Box boundaries of his1
!
      ixmin1 = his1%imin(1) - 1
      ixmax1 = his1%imax(1)
      iymin1 = his1%imin(2) - 1
      iymax1 = his1%imax(2)
      izmin1 = his1%imin(3) - 1
      izmax1 = his1%imax(3)
      allocate( x1( ixmin1:ixmax1 ) )
      allocate( y1( iymin1:iymax1 ) )
      allocate( z1( izmin1:izmax1 ) )
      x1(:) = [ (his1%hbeg(1) + ix1*his1%wbox(1), ix1=ixmin1,ixmax1) ]
      y1(:) = [ (his1%hbeg(2) + iy1*his1%wbox(2), iy1=iymin1,iymax1) ]
      z1(:) = [ (his1%hbeg(3) + iz1*his1%wbox(3), iz1=izmin1,izmax1) ]
!
! Box boundaries of his2. Adjust so his1 completely fits into his2
!
      epsx = max( EPSMA*abs(x1(ixmin1)), SMALL )
      epsy = max( EPSMA*abs(y1(iymin1)), SMALL )
      epsz = max( EPSMA*abs(z1(izmin1)), SMALL )
      call scoreHis3d( 1, [ x1(ixmin1)+epsx, y1(iymin1)+epsy, z1(izmin1)+epsz ], his2, ZERO )
!          ==========  i                            i                             io    i

      epsx = max( EPSMA*abs(x1(ixmax1)), SMALL )
      epsy = max( EPSMA*abs(y1(iymax1)), SMALL )
      epsz = max( EPSMA*abs(z1(izmax1)), SMALL )
      call scoreHis3d( 1, [ x1(ixmax1)-epsx, y1(iymax1)-epsy, z1(izmax1)-epsz ], his2, ZERO )
!          ==========  i                            i                             io    i

      ixmin2 = his2%imin(1) - 1
      ixmax2 = his2%imax(1)
      iymin2 = his2%imin(2) - 1
      iymax2 = his2%imax(2)
      izmin2 = his2%imin(3) - 1
      izmax2 = his2%imax(3)
      allocate( x2( ixmin2:ixmax2 ) )
      allocate( y2( iymin2:iymax2 ) )
      allocate( z2( izmin2:izmax2 ) )
      x2(:) = [ (his2%hbeg(1) + ix2*his2%wbox(1), ix2=ixmin2,ixmax2) ]
      y2(:) = [ (his2%hbeg(2) + iy2*his2%wbox(2), iy2=iymin2,iymax2) ]
      z2(:) = [ (his2%hbeg(3) + iz2*his2%wbox(3), iz2=izmin2,izmax2) ]
!
!  Merge his1 into his2
!
      xold = x1(ixmin1)
      ix1 = ixmin1 + 1
      do ix2 = ixmin2+1, ixmax2
         if( x2(ix2) > xold )  exit
      end do
      ix2 = min( ix2, ixmax2 )

      do                                        ! loop over boxes in x direction

         yold = y1(iymin1)
         iy1 = iymin1 + 1
         do iy2 = iymin2+1, iymax2
            if( y2(iy2) > yold )  exit
         end do
         iy2 = min( iy2, iymax2 )

         do                                     ! loop over boxes in y direction

            zold = z1(izmin1)
            iz1 = izmin1 + 1
            do iz2 = izmin2+1, izmax2
               if( z2(iz2) > zold )  exit
            end do
            iz2 = min( iz2, izmax2 )

            do                                  ! loop over boxes in z direction

               dx = min( x1(ix1)-xold, x2(ix2)-xold )
               dy = min( y1(iy1)-yold, y2(iy2)-yold )
               dz = min( z1(iz1)-zold, z2(iz2)-zold )
               his2%hist(:,ix2,iy2,iz2) = his2%hist(:,ix2,iy2,iz2)  &
                     + his1%hist(:,ix1,iy1,iz1) * dx * dy * dz  &
                        / (his2%wbox(1) * his2%wbox(2) * his2%wbox(3))

               if( z1(iz1) < z2(iz2) ) then
                  zold = z1(iz1)
                  iz1 = iz1 + 1
                  if( iz1 > izmax1 )  exit
               else
                  zold = z2(iz2)
                  iz2 = iz2 + 1
                  if( iz2 > izmax2 )  exit
               end if

            end do

            if( y1(iy1) < y2(iy2) ) then
               yold = y1(iy1)
               iy1 = iy1 + 1
               if( iy1 > iymax1 )  exit
            else
               yold = y2(iy2)
               iy2 = iy2 + 1
               if( iy2 > iymax2 )  exit
            end if

         end do

         if( x1(ix1) < x2(ix2) ) then
            xold = x1(ix1)
            ix1 = ix1 + 1
            if( ix1 > ixmax1 )  exit
         else
            xold = x2(ix2)
            ix2 = ix2 + 1
            if( ix2 > ixmax2 )  exit
         end if

      end do

   end subroutine mergeHis3d
!
!===========================================================================================
!
   subroutine broadcastHis3d( his )
!             =====================  io
!
!     Spread his of master thread to other threads
!
!     This subroutine must be called from within a parallel region
!
!  Arguments:
!
      type(His3d_type), intent(inout) :: his    ! histogram to be merged
!
!  Local Variables:
!
      type(His3d_type), save :: his_tmp         ! shared temporary histogram
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      his_tmp = his
!$omp end master
!$omp barrier
      his = his_tmp

   end subroutine broadcastHis3d
!
!===========================================================================================
!
   subroutine writeHis3dScalar( fname, fmode, his, hmin_, hmax_, hcrit_, factor_ )
!             ================    i      i     i     i      i       i       i
!  
!     Writes the histogram to a specified file. Histogram values are limited 
!     to the range [hmin_,hmax_]. Histogram values larger than hcrit_ are set
!     to hmax_. hmin_, hmax_, and hcrit_ apply to all histograms.
!
!  Arguments:
!
      character(*), intent(in) :: fname            ! File name
      character(*), intent(in) :: fmode            ! indicates whether output should be
                                                   !    appended or overwrite file
      type(His3d_type), intent(in) :: his          ! Histogram structure
      real(REAL8), intent(in), optional :: hmin_   ! Lower cut-off value times wbox
      real(REAL8), intent(in), optional :: hmax_   ! Upper cut-off value
      real(REAL8), intent(in), optional :: hcrit_  ! Critical value (see Purpose)
      real(REAL8), intent(in), optional :: factor_ ! Factor to multiply the histogram with
!
!  Local Variables:
!
      real(REAL8) :: hmin(size( his%hist , 1 ))
      real(REAL8) :: hmax(size( his%hist , 1 ))
      real(REAL8) :: hcrit(size( his%hist , 1 ))
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( hmin_ ) ) then
         hmin(:) = hmin_
      else
         hmin(:) = -BIG
      end if

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

      call writeHis3dRaw( fname, fmode, his, hmin, hmax, hcrit, factor_ )
!          =============    i      i     i    i     i      i       i

   end subroutine writeHis3dScalar
!
!=======================================================================
!
   subroutine writeHis3dArray( fname, fmode, his, hmin, hmax_, hcrit_, factor_ )
!             ===============    i       i    i    i      i      i        i
!  
!     Writes the histogram to a specified file. Histogram values are limited 
!     to the range [hmin,hmax_]. Histogram values larger than hcrit_ are set
!     to hmax_. hmin, hmax_, and hcrit_ apply to all histograms.
!
!  Arguments:
!
      character(*), intent(in) :: fname               ! File name
      character(*), intent(in) :: fmode               ! indicates whether output should be
                                                      !    appended or overwrite file
      type(His3d_type), intent(in) :: his             ! Histogram structure
      real(REAL8), intent(in) :: hmin(:)              ! Lower cut-off value times wbox
      real(REAL8), intent(in), optional :: hmax_(:)    ! Upper cut-off value
      real(REAL8), intent(in), optional :: hcrit_(:)   ! Critical value (see Purpose)
      real(REAL8), intent(in), optional :: factor_     ! Factor to multiply the histogram with
!
!  Local Variables:
!
      real(REAL8) :: hmax(size( his%hist, 1 ))
      real(REAL8) :: hcrit(size( his%hist, 1 ))
!
!  Start of Executable Code:
!  ------------------------
!
      if( size( hmin ) /= size( his%hist, 1 ) )  &
            error stop 'writeHis3dArray: incompatible size of hmin'

      if( present( hmax_ ) ) then
         if( size( hmax_ ) /= size( his%hist, 1 ) )  &
               error stop 'writeHis3dArray: incompatible size of hmax'
         hmax(:) = hmax_(:)
      else
         hmax(:) = BIG
      end if

      if( present( hcrit_ ) ) then
         if( size( hcrit_ ) /= size( his%hist, 1 ) )  &
               error stop 'writeHis3dArray: incompatible size of hcrit'
         hcrit(:) = hcrit_(:)
      else
         hcrit(:) = BIG
      end if

      call writeHis3dRaw( fname, fmode, his, hmin, hmax, hcrit, factor_ )
!          =============    i      i     i    i     i      i       i

   end subroutine writeHis3dArray
!
!=======================================================================
!
   subroutine writeHis3dRaw( fname, fmode, his, hmin, hmax, hcrit, factor_ )
!             =============    i      i     i    i     i      i      i
!  
!     Writes the histogram to a specified file. Histogram values are limited 
!     to the range [hmin,hmax]. Histogram values larger than hcrit are set
!     to hmax.
!
!  Arguments:
!
      character(*), intent(in) :: fname         ! File name
      character(*), intent(in) :: fmode         ! indicates whether output should be
                                                !    appended or overwrite file
      type(His3d_type), intent(in) :: his       ! Histogram structure
      real(REAL8), intent(in) :: hmin(:)        ! Lower cut-off value times wbox
      real(REAL8), intent(in) :: hmax(:)        ! Upper cut-off value
      real(REAL8), intent(in) :: hcrit(:)       ! Critical value (see Purpose)
      real(REAL8), intent(in), optional :: factor_  ! Factor to multiply the histogram with
!
!  Local Variables:
!
      character(80) :: text
      character(256) :: errmsg
!      logical, allocatable :: mask(:,:,:)
      integer(INT4) :: ia, i, iboxx, iboxy, iboxz, istat, i1, i2, i3
      real(REAL8) :: factor
      real(REAL8), allocatable :: hist(:,:,:)
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
!  Do not write file if histogram is empty
!
      if( his%imin(1) > his%imax(1) )  return
!
!  Open file
!
      call openFile( 22, fname, 'unformatted', fmode )
!
!  Write header
!
      text = 'Atom numbers = '
      do ia = 1, size( his%hist, 1 ) - 1
         write( text(len_trim(text)+1:), '(i3,a2)' )  his%indhis(ia), ', '
      end do
      write( text(len_trim(text)+1:), '(i3)' )  his%indhis(ia)

      write( 22 )  text
!
!  Write dimensions and position vectors
!
      write( 22 )  size( his%hist, 1 ),  &
                   his%imax(1)-his%imin(1)+1,  &
                   his%imax(2)-his%imin(2)+1,  &
                   his%imax(3)-his%imin(3)+1
      write( 22 )  ( his%hbeg(1)+( i-HALF )*his%wbox(1), i=his%imin(1),his%imax(1) )
      write( 22 )  ( his%hbeg(2)+( i-HALF )*his%wbox(2), i=his%imin(2),his%imax(2) )
      write( 22 )  ( his%hbeg(3)+( i-HALF )*his%wbox(3), i=his%imin(3),his%imax(3) )
!
!  Write histogram
!
      allocate( hist( his%imin(1):his%imax(1), his%imin(2):his%imax(2),  &
                      his%imin(3):his%imax(3) ), stat=istat, errmsg=errmsg )
      call check_allocation( 'writeHis3dRaw', 'hist', istat, errmsg )

      do ia = 1, size( his%hist, 1 )
!
!  Set damage concentration to hmax if it exceeds hcrit
!
         hist(:,:,:) = his%hist(ia,his%imin(1):his%imax(1), his%imin(2):his%imax(2),  &
                                  his%imin(3):his%imax(3))
         forall( i1=his%imin(1):his%imax(1), i2=his%imin(2):his%imax(2),  &
                 i3=his%imin(3):his%imax(3),  &
                 his%hist(ia,i1,i2,i3) > hcrit(ia) )
            hist(i1,i2,i3) = hmax(ia)
         end forall
!
!  Cut off values lower than hmin/(wbox(1)*wbox(2)*wbox(3))
!
#ifdef GTS
!         hist(:,:,:) = max( hist(:,:,:), ZERO )
#else
         hist(:,:,:) = max( hist(:,:,:), hmin(ia)/(his%wbox(1)*his%wbox(2)*his%wbox(3)) )   ! good for superposition ?
#endif
!
!  Write histogram values
!
         write( 22 )  ( ( ( factor * hist(iboxx,iboxy,iboxz),  &
                            iboxx=his%imin(1),his%imax(1) ),  &
                          iboxy=his%imin(2),his%imax(2) ),  &
                        iboxz=his%imin(3),his%imax(3) )

      end do
!
!  Deallocate auxiliary array and close the file
!
      deallocate( hist )

      close( 22 )
      
   end subroutine writeHis3dRaw
!
!=======================================================================
!
   subroutine saveHis3d( his )
!             =========   i
!  
!     saves the his structure to the file with unit number 22
!
!  Arguments:
!
      type(His3d_type), intent(in) :: his       ! histogram structure
!
!  Local Variables:
!
      integer(INT4) :: nhis
!
!  Start of Executable Code:
!  ------------------------
!
      nhis = size( his%indhis(:) )
      write( 22 )  nhis, his%nbox(:)
      write( 22 )  his%fixedsize
      write( 22 )  his%hbeg(:), his%wbox(:), his%imin(:), his%imax(:)
      write( 22 )  his%hist(:,:,:,:), his%indhis(:)

   end subroutine saveHis3d
!
!=======================================================================
!
   subroutine loadHis3d( his )
!             =========   o
!  
!     Load the his structure to the file with unit number 22
!
!  Arguments:
!
      type(His3d_type), intent(out) :: his      ! histogram structure
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
      read( 22, iostat=istat )  nhis, his%nbox(:)
      read( 22, iostat=istat )  his%fixedsize
      read( 22, iostat=istat )  his%hbeg(:), his%wbox(:), his%imin(:), his%imax(:)

      if( allocated( his%hist ) ) then
         deallocate( his%hist )
         deallocate( his%indhis )
      end if

      allocate( his%hist( nhis, his%nbox(1), his%nbox(2), his%nbox(3) ), stat=istat,  &
                errmsg=errmsg )
      call check_allocation( 'loadHis3d', 'his%hist', istat, errmsg )

      allocate( his%indhis( nhis ), stat=istat, errmsg=errmsg )
      call check_allocation( 'loadHis3d', 'his%indhis', istat, errmsg )

      read( 22, iostat=istat )  his%hist(:,:,:,:), his%indhis(:)

      if( istat /= 0 )  error stop 'error during reading backup file'      
     
   end subroutine loadHis3d
!
!=======================================================================
!
end module Stat_His3d_mod

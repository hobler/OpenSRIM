module Eval_RRobin_mod
!
!  Statistical evaluation for COSIRES'94 Round Robin simulation
!  (still useful as a benchmark)
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, PI

   use Param_Ion_mod,   only: ion
   use Param_Xtal_mod,  only: xtal

   use Open_mod,        only: openFile
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
   public :: rob10                  ! initialize statistic quantities
!            =====
   public :: rob25                  ! score statistic quantities
!            =====
   public :: rob30                  ! evaluate statistic quantities and 
!            =====                  !    write results to a file
!
!  Private Procedures:
!  ==================
!
   private :: mergeIntoMasterInt1d           ! Merge integer array into master thread
!             ====================
   private :: mergeIntoMasterInt             ! Merge integer variable into master thread
!             ==================
   private :: mergeIntoMasterReal            ! Merge real variable into master thread
!             ===================
!
!  Private Data:
!  ============
!
   integer(INT4), parameter :: NE = 100
   integer(INT4), parameter :: NA = 30

   integer(INT4) :: n(NE), m(NA)
   integer(INT4) :: ntrans
   real(REAL8) :: eloss, angle2
!$omp threadprivate( n, m, ntrans, eloss, angle2 )
!
contains
!
!=======================================================================
!
   subroutine rob10( )
!             =====
!
!     rob10 - initialize statistic quantities for COSIRES'94 Round Robin
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
      n(:) = 0
      m(:) = 0

      ntrans = 0
      eloss = ZERO
      angle2 = ZERO

   end subroutine rob10
!
!=======================================================================
!
   subroutine rob25( e, pos, dir, lnot )
!             =====  i   i    i    o
!
!     Score statistic quantities of transmitted ions for COSIRES'94
!     Round Robin
!
!  Arguments:
!
      real(REAL8), intent(in) :: e              ! energy of the transmitted ion [eV]
      real(REAL8), intent(in) :: pos(3)         ! position of the ion
      real(REAL8), intent(in) :: dir(3)         ! direction of motion of the ion
      logical, intent(out) :: lnot              ! flag indicating that the ion was not counted
                                                !    as transmitted since it was emitted from 
                                                !    the side walls
!
!  Local Variables:
!
      integer(INT4) :: ia, ie
      real(REAL8) :: angle, fac, xp, yp, xmax
!
!  Start of Executable Code:
!  ------------------------
!
      angle = 180/PI * acos( dir(3) )
      fac = (10.86 - pos(3)) / dir(3)
      xp = pos(1) + fac*dir(1)
      yp = pos(2) + fac*dir(2)

      xmax = 13.575

      if( angle <= 80 .and.  &
          xp >= -xmax .and. xp <= xmax .and. yp >= -xmax .and. yp <= xmax ) then

         lnot = .false.
                                                ! another transmitted ion...
                                                !    (notice: the lateral limits must be 
                                                !     the same as in selncc)
         ntrans = ntrans + 1
                                                ! evaluate angle statistics
         ia = min( NA, max( 1, int( NA*angle/90 + 1 ) ) )
         m(ia) = m(ia) + 1
         angle2 = angle2 + angle**2
                                                ! evaluate energy statistics      
         ie = min( NE, max( 1, int( NE*e/ion%e + 1 ) ) )
         n(ie) = n(ie) + 1
         eloss = eloss + (ion%e - e)
      else
         lnot = .true.
      endif   

   end subroutine rob25
!
!=======================================================================
!
   subroutine rob30( )
!             =====
!
!     Write output of statistics for COSIRES'94 Round Robin
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ie, ia
!
!  Start of Executable Code:
!  ------------------------
!
!  Gather statistics of other threads into master thread
!
      call mergeIntoMasterInt1d( n )
!          ====================  io
      call mergeIntoMasterInt1d( m )
!          ====================  io
      call mergeIntoMasterInt( ntrans )
!          ==================    io
      call mergeIntoMasterReal( eloss )
!          ===================   io
      call mergeIntoMasterReal( angle2 )
!          ===================    io

!$omp master
      call openFile( 22, 'DAT', 'formatted', 'write' )
!          ========  i     i        i           i

      write( 22,'(2i4)' ) 2, NE
      write( 22,'(2i4)' ) 0, 0
      do ie=1,NE
         write( 22,'(i5,i8)' ) ie, n(ie)
      end do
      write( 22,'(a,e15.5)' ) '*', eloss / max( ntrans, 1 )
      write( 22,'(2i4)' ) 2,NA
      write( 22,'(2i4)' ) 0,0
      do ia=1,NA
         write( 22,'(i5,i8)' ) ia, m(ia)
      end do
      write( 22,'(a,e15.5)' ) '*', angle2 / max( ntrans, 1 )
      write( 22,'(a,i5)' ) '*', ntrans

      close (22)
!$omp end master
      
   end subroutine rob30
!
!=======================================================================
!
   subroutine mergeIntoMasterInt1d( i )
!             ====================  io
!
!     Merge integer array of other threads into array of master thread
!
!     To be called from within a parallel region
!
!  Arguments:
!
      integer(INT4), intent(inout) :: i(:)   ! array to be merged
!
!  Local Variables:
!
      integer(INT4), allocatable, save :: i_tmp(:) ! shared temporary array
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      allocate( i_tmp, source=i )
      i_tmp(:) = 0
!$omp end master

!$omp barrier
!$omp critical( critical_rrobin_mergeIntoMasterInt1d )
      i_tmp(:) = i_tmp(:) + i(:)
!$omp end critical( critical_rrobin_mergeIntoMasterInt1d )

!$omp barrier
!$omp master
      i(:) = i_tmp(:)
      deallocate( i_tmp )
!$omp end master

   end subroutine mergeIntoMasterInt1d
!
!=======================================================================
!
   subroutine mergeIntoMasterInt( i )
!             ==================  io
!
!     Merge integer variables of other threads into array of master thread
!
!     To be called from within a parallel region
!
!  Arguments:
!
      integer(INT4), intent(inout) :: i      ! integer to be merged
!
!  Local Variables:
!
      integer(INT4), allocatable, save :: i_tmp    ! shared temporary variable
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      i_tmp = 0
!$omp end master

!$omp barrier
!$omp critical( critical_rrobin_mergeIntoMasterInt )
      i_tmp = i_tmp + i
!$omp end critical( critical_rrobin_mergeIntoMasterInt )

!$omp barrier
!$omp master
      i = i_tmp
!$omp end master

   end subroutine mergeIntoMasterInt
!
!=======================================================================
!
   subroutine mergeIntoMasterReal( r )
!             ===================  io
!
!     Merge real variable of other threads into array of master thread
!
!     This subroutine must be called from within a parallel region
!
!  Arguments:
!
      real(REAL8), intent(inout) :: r      ! variable to be merged
!
!  Local Variables:
!
      real(REAL8), allocatable, save :: r_tmp ! shared temporary variable
!
!  Start of Executable Code:
!  ------------------------
!
!$omp master
      r_tmp = 0
!$omp end master

!$omp barrier
!$omp critical( critical_rrobin_mergeIntoMasterReal )
      r_tmp = r_tmp + r
!$omp end critical( critical_rrobin_mergeIntoMasterReal )

!$omp barrier
!$omp master
      r = r_tmp
!$omp end master

   end subroutine mergeIntoMasterReal
!
!=======================================================================
!
end module Eval_RRobin_mod

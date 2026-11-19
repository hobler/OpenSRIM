module Time_mod
!
!  Manage getting and formatting the CPU time
!
!  Used Modules:
!  ============
!
   use Const_mod,    only: INT4, REAL8,  &
                           ZERO
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
   public :: initCPUtime                     ! set the CPU time to 0 or a given time
!            ===========
   public :: getCPUtime            
!            ==========
   interface getCPUtime
      module procedure getCPUtime_r          ! get CPU time in seconds
      module procedure getCPUtime_c          ! get CPU time in "hh:mm:ss"
   end interface
   public :: initElapsedTime                 ! initialize the elapsed time
!            ===============
   public :: getElapsedTime                  ! get the elapsed time in "hh:mm:ss"
!            ==============
!
!  Private Data:
!  ============
!
   integer(INT4) :: itimbeg
   real(REAL8) :: timbeg
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initCPUtime( time_init_ )
!             ===========      i
!
!     Set the CPU time to 0 or time_init_
!
!  Arguments:
!
      real(REAL8), intent(in), optional :: time_init_
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call cpu_time( timbeg )                    ! Fortran95 function !
!          ========    o

      if( present( time_init_ ) )  timbeg = timbeg - time_init_

   end subroutine initCPUtime
!
!=======================================================================
!
   subroutine getCPUtime_c( ctime )
!             ============    o
!
!     Get the CPU time in the format "hhh:mm:ss"
!
!  Arguments:
!
      character(9), intent(out) :: ctime
!
!  Local Variables:
!
      integer(INT4) :: itimh, itimm, itims
      real(REAL8) :: time
!
!  Start of Executable Code:
!  ------------------------
!
      call getCPUtime_r( time )
!          ============   o
      time = time - timbeg
      itimh = time / 3600
      itimm = (time-itimh*3600) / 60
      itims = time - itimh*3600 - itimm*60

      write( ctime,'(i3.2,a,i2.2,a,i2.2)' ) itimh, ':', itimm, ':', itims

   end subroutine getCPUtime_c
!
!=======================================================================
!
   subroutine getCPUtime_r( time )
!             ============   o
!
!     Get the CPU time in seconds
!
!  Arguments:
!
      real(REAL8), intent(out) :: time
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call cpu_time( time )                     ! Fortran95 function
!          ========   o

   end subroutine getCPUtime_r
!
!=======================================================================
!
   subroutine initElapsedTime( )
!             ===============
!
!     Set the CPU time to 0 or time_init
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
      call system_clock( itimbeg )                    ! Fortran95 function !
!          ============     o

   end subroutine initElapsedTime
!
!=======================================================================
!
   subroutine getElapsedTime( ctime )
!             ==============    o
!
!     Get the CPU time in the format "hhh:mm:ss"
!
!  Arguments:
!
      character(9), intent(out) :: ctime
!
!  Local Variables:
!
      integer(INT4) :: itimend, itime, count_rate, count_max
      integer(INT4) :: itimh, itimm, itims
      real(REAL8) :: time
!
!  Start of Executable Code:
!  ------------------------
!
      call system_clock( itimend, count_rate, count_max )
!          ============     o         o           o
      itime = modulo( itimend - itimbeg, count_max )
      time = real( itime ) / real( count_rate )
      itimh = time / 3600
      itimm = (time-itimh*3600) / 60
      itims = time - itimh*3600 - itimm*60

      write( ctime,'(i3.2,a,i2.2,a,i2.2)' ) itimh, ':', itimm, ':', itims

   end subroutine getElapsedTime
!
!=======================================================================
!
end module Time_mod

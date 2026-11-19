module Eval_Tra_mod
!
!  Records the trajectories
!
!  Lines are written to the file delayed to eliminate duplicate points on the file when
!  a trajectory starts or ends
!
!  Code requires single-threaded execution
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 BIG, ZERO, ONE, PI
   use Vector_mod,         only: operator(.dot.)
   use Types_mod,          only: Proj_type

   use Param_Setup_mod,    only: niontra, nthread
   use Param_Ion_mod,      only: ion
   use Param_Mat_mod,      only: mat, coord_xtal, irxtal
   use Param_Out_mod,      only: out

   use Geom_CoordSys_mod,  only: CoordSys_type, initCoordSys, toCoordSys
   use Ion_Proj_mod,       only: getXStart
   use Open_mod,           only: openFile
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
   public :: tra10                              ! initialization
!            =====
   public :: tra12                              ! score at start of ion trajectory
!            =====
   public :: tra2x                              ! score projectile
!            =====
   public :: tra30                              ! output of statistical quantities
!            =====
!
!  Private Data:
!  ============
!
   character(78) :: line                        ! line to be written to TRA file
   type(CoordSys_type) :: coord_statis
!
contains
!
!=======================================================================
!
   subroutine tra10( )
!             =====
!
!     Initialize the trajectory recording
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
      if( nthread > 1 )  stop 'tra10: Trajectory statistic requires sequential execution'
!
!  Initialize line to be written to TRA file
!
      line = ''
!
!  Initialize reference coordinate system for statistics
!
      select case( out%mcoord )
         case( 3 )
            call initCoordSys( coord_statis, rotate_=ion%rotate, tiltx_=ion%tiltx,  &
!                ============       o                i                  i
                               tilty_=ion%tilty )
!                                      i
         case( 4 )
            call initCoordSys( coord_statis, out%dira)
!                ============       o          i
         case( 5 )
            if( .not. mat(irxtal)%lpoly )  coord_statis = coord_xtal
      end select
!
!  Open the trajectory file
!
      call openFile( 23, 'TRA', 'formatted', 'write' )
!          ========  i     i        i           i

      return

   end subroutine tra10
!
!=======================================================================
!
   subroutine tra12( proj )
!             =====   i
!
!     Score statistic quantities at the beginning of the simulation of a 
!     new ion atom
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! starting point of trajectory
!
!  Local Variables:
!
      real(REAL8) :: posdum(3)
      real(REAL8) :: dir_s(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  Initialize the statistics coordinate system
!
      select case( out%mcoord )
         case( 2 )
            call getXStart( posdum, dir_s )
!                =========    o       o
            call initCoordSys( coord_statis, dir_s )
!                ============       o          i
         case( 5 )
            if( mat(irxtal)%lpoly )  coord_statis = coord_xtal
      end select
!
!  Save starting point or write it to the trajectory file
! 
      if( out%ltra(-1) .or. out%ltra(1) )  call tra2x( proj, 1 )
!                                               =====   i    i

   end subroutine tra12
!
!=======================================================================
!
   subroutine tra2x( proj, iflag )
!             =====   i      i
!
!     Store trajectory points of projectile stack
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      integer(INT4), intent(in) :: iflag        ! flag indicating position within trajectory:
                                                !    iflag=0: no trajectory
                                                !    iflag=-1,1: start of trajectory
                                                !    iflag=2: along a trajectory
                                                !    iflag=3: end of trajectory inside target
                                                !    iflag=4: backscattered projectile
                                                !    iflag=5: transmitted projectile
!
!  Local Variables:
!
      logical :: lstart, lend
      integer(INT4) :: ielem, ig
      real(REAL8) :: x, y, z, dirx, diry, dirz, e, angle
      real(REAL8) :: pos_s(3),  dir_s(3)
!
!  Start of Executable Code:
!  ------------------------
!
      if( mod( proj%iion, niontra ) > 0 )  return

      if( out%mcoord==1 ) then
         if( iflag >= 0 ) then
            pos_s(:) = proj%pos(:)
            dir_s(:) = proj%dir(:)
         else
            pos_s(:) = proj%origin%pos(:)
            dir_s(:) = proj%origin%dir(:)
         end if
      else
         if( iflag >= 0 ) then
            call toCoordSys( coord_statis, proj%pos, pos_s )
!                ==========       i           i        o
            call toCoordSys( coord_statis, proj%dir, dir_s )
!                ==========       i           i        o
         else
            call toCoordSys( coord_statis, proj%origin%pos, pos_s )
!                ==========       i               i           o
            call toCoordSys( coord_statis, proj%origin%dir, dir_s )
!                ==========       i               i           o
         end if
      end if

      x = pos_s(1)
      y = pos_s(2)
      z = pos_s(3)
      dirx = dir_s(1)
      diry = dir_s(2)
      dirz = dir_s(3)

      e = proj%e
      if( proj%virtual ) then
         ielem = - proj%ielem
      else
         ielem = proj%ielem
      end if
      ig = proj%ig
!
!  In case of backscattered and transmitted projectiles, count them only if they are within
!  the analyzer
!
      select case( iflag )

         case( 4 )

            angle = 180/PI * acos( -dirz )
            if( angle < out%amina .or. angle > out%amaxa .or. &
                max( ZERO, e ) < out%emina .or. e > out%emaxa )  return

         case( 5 )

            angle = 180/PI * acos( dirz )
            if( angle < out%amina .or. angle > out%amaxa .or. &
                max( ZERO, e ) < out%emina .or. e > out%emaxa )  return

      end select
!
!  Write point to the trajectory file, avoiding "along the trajectory" points duplicating
!  "beginning of trajectory" or "end of trajectory" points
!
      lstart = (line(77:78) == ' 1' .and. iflag == 2)
      lend = (line(77:78) == ' 2' .and. iflag == 3)
      if( lstart ) then
         write( 23, '(a)' )  line
         line = ''
      else
         if( .not.lend .and. line /= '' )  write( 23, '(a)' )  line
         write( line, '(3f11.3, 3f8.4, f12.2, 3i3)' )  &
               x, y, z, dirx, diry, dirz, e, ielem, ig, iflag
       end if

   end subroutine tra2x
!
!=======================================================================
!
   subroutine tra30( )
!             =====
!
!     Trajectory file related actions at the end of the simulation
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
!  Write last line
!
      write( 23, '(a)' )  line
!
!  Close trajectory file
! 
      close( 23 )

   end subroutine tra30
!
!=======================================================================
!
end module Eval_Tra_mod

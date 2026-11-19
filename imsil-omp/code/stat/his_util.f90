module Stat_His_Util_mod
!
!  Utility programs for histograms
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, HALF

   use Param_Setup_mod,       only: ndim, nelemIon, nelem, nr, ldyn
   use Param_Geom_mod,        only: geom
   use Param_Dam_mod,         only: dam
   use Param_Mat_mod,         only: irxtal

   use Stat_His1d_mod,        only: His1d_type
   use Stat_His2d_mod,        only: His2d_type
   use Stat_His3d_mod,        only: His3d_type
   use Target_Surface2d_mod,  only: getDistSurface2d
   use Target_Surface3d_mod,  only: getDistSurface3d
!   use Target_mod,            only: getDomainIndex
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
   public :: getAmoCPDD                      ! get amorphous zones
!            ==========

   public :: deleteInAir                     ! delete concentrations in the air
!            ===========
   interface deleteInAir
      module procedure deleteInAir1d
      module procedure deleteInAir2d
      module procedure deleteInAir3d
   end interface
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine deleteInAir1d( his1 )
!             =============   io
!
!     Remove all histogram boxes that are completely in the air
!
!  Arguments:
!
      type(His1d_type), intent(inout) :: his1   ! histogram
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( ndim /= 1 ) return
      if( geom%lcell )  return

      if( geom%posif(0) > his1%hbeg ) then
         his1%imin = max( his1%imin, int( ( geom%posif(0) - his1%hbeg ) / his1%wbox ) )
         his1%hist(:,1:his1%imin-1) = ZERO
      end if

      if( geom%posif(nr) < his1%hbeg + his1%nbox*his1%wbox ) then
         his1%imax = min( his1%imax, int( ( geom%posif(nr) - his1%hbeg ) / his1%wbox ) )
         his1%hist(:, his1%imax+1:his1%nbox ) = ZERO
      end if

   end subroutine deleteInAir1d
!
!=======================================================================
!
   subroutine deleteInAir2d( his2 )
!             =============   io
!
!     Remove all histogram boxes that are completely in the air
!
!  Arguments:
!
      type(His2d_type), intent(inout) :: his2   ! histogram
!
!  Local Variables:
!
!      logical :: outside
      integer(INT4) :: iboxx, iboxz
      real(REAL8) :: sdf, pos(3)
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcell )  return

      do iboxz = his2%imin(2), his2%imax(2)
         do iboxx = his2%imin(1), his2%imax(1)
            pos(:) = [ his2%hbeg(1) + (iboxx-0.5) * his2%wbox(1), ZERO, &
                       his2%hbeg(2) + (iboxz-0.5) * his2%wbox(2) ]
!            call getDomainIndex( pos(:), outside_=outside )
!!                ==============    i                o
!            if( outside ) then
            call getDistSurface2d( pos(:), sdf )
!                ================    i      o
            if( sdf < ZERO ) then
               his2%hist(:,iboxx,iboxz) = ZERO
            end if
         end do
      end do

   end subroutine deleteInAir2d
!
!=======================================================================
!
   subroutine deleteInAir3d( his3 )
!             =============   io
!
!     Remove all histogram boxes that are completely in the air
!
!  Arguments:
!
      type(His3d_type), intent(inout) :: his3   ! histogram
!
!  Local Variables:
!
!      logical :: outside
      integer(INT4) :: iboxx, iboxy, iboxz
      real(REAL8) :: sdf, pos(3)
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcell )  return

      do iboxz = his3%imin(3), his3%imax(3)
         do iboxy = his3%imin(2), his3%imax(2)
            do iboxx = his3%imin(1), his3%imax(1)
               pos(:) = [ his3%hbeg(1) + (iboxx-0.5) * his3%wbox(1), &
                          his3%hbeg(2) + (iboxy-0.5) * his3%wbox(2), &
                          his3%hbeg(3) + (iboxz-0.5) * his3%wbox(3) ]
!               call getDomainIndex( pos(:), outside_=outside )
!!                   ==============    i            o
!               if( outside )  his3%hist(:,iboxx,iboxy,iboxz) = ZERO
               call getDistSurface3d( pos(:), sdf )
!                   ================    i      o
               if( sdf < ZERO )  his3%hist(:,iboxx,iboxy,iboxz) = ZERO
            end do
         end do
      end do

   end subroutine deleteInAir3d
!
!=======================================================================
!
   subroutine getAmoCPDD( his1, namo, zamo_left, zamo_right )
!             ==========   i     o        o          o
!
!     Get the amorphous zones
!
!  Arguments:
!
      type(His1d_type), intent(in) :: his1      ! 1d histogram containing ions and interstitials
      integer(INT4), intent(out) :: namo        ! number of amorphous regions
      real(REAL8), intent(out) :: zamo_left(:)  ! left boundary of amorphous region
      real(REAL8), intent(out) :: zamo_right(:) ! right boundary of amorphous region
!
!  Local Variables:
!
      logical :: amorph
      integer(INT4) :: ibox
      real(REAL8) :: zl, zu, h, hl, ha, sumhis, sumhisl
!
!  Start of Executable Code:
!  ------------------------
!
      namo = 0
      amorph = .false.

      do ibox = his1%imin, his1%imax

         zl = his1%hbeg + (ibox-1)*his1%wbox
         zu = his1%hbeg + ibox*his1%wbox

         if( ldyn ) then
            sumhis = sum( his1%hist(nelemIon+1:nelemIon+nelem,ibox) )
            if( ibox > 1 )  sumhisl = sum( his1%hist(nelemIon+1:nelemIon+nelem,ibox-1) )
         else
            sumhis = sum( his1%hist(nelemIon+1:nelem,ibox) )
            if( ibox > 1 )  sumhisl = sum( his1%hist(nelemIon+1:nelem,ibox-1) )
         end if
                                                !
         if( amorph ) then                      ! If amorphous, check whether target gets
                                                !    crystalline again
            if( sumhis > dam%damamo ) then
               if( ibox < his1%imax )  cycle               ! remain amorphous
               zamo_right(namo) = zu
            else                                              ! last box amorphous
               h = log(sumhis)                              ! get crystalline
               hl = log(sumhisl)
               ha = log(dam%damamo)
               zamo_right(namo) = zl + HALF*his1%wbox*(2*ha-h-hl)/(h-hl)
            end if

            zamo_right(namo) = max( min( zamo_right(namo), geom%posif(irxtal) ), &
                                    geom%posif(irxtal-1) )

            if( zamo_right(namo) == zamo_left(namo) )  namo = namo - 1

            amorph = .false.
                                                !
         else                                   ! If crystalline, check whether target gets
                                                !    amorphous (again)
            if( sumhis > dam%damamo ) then
               namo = namo + 1                  ! gets amorphous
               if( namo > size( zamo_left ) )  &
                     error stop 'getAmoCPDD: too many amorphous regions'
               if( ibox == his1%imin ) then
                  zamo_left(namo) = zl
               else
                  h = log(sumhis)
                  hl = log(sumhisl)
                  ha = log(dam%damamo)
                  zamo_left(namo) = zl + HALF*his1%wbox*(2*ha-h-hl)/(h-hl)
               end if

               zamo_left(namo) = max( min( zamo_left(namo), geom%posif(irxtal) ),  &
                                      geom%posif(irxtal-1) )

               amorph = .true.

            end if

         end if

      end do

   end subroutine getAmoCPDD
!
!=======================================================================
!
end module Stat_His_Util_mod

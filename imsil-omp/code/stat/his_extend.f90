module Stat_His_Extend_mod
!
!  Contains subroutines that extends histogram as to accomodate new scores without increasing
!  the histogram box widths. THESE SUBROUTINES ALSO RETURN THE INDEX OF A SCORE.
!
!  Used Modules:
!  ============
!
   use Const_mod, only: INT4, REAL8,  &
                        EPSMA, ZERO, ONE, TWO, HALF, THIRD
   use Array_mod, only: changeShape
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
   public :: extendHis1d                        ! determine 1d histogram index and extend
!            ===========                        !    histogram if necessary
   public :: extendHis2d                        ! determine 1d histogram index and extend
!            ===========                        !    histogram if necessary
   public :: extendHis3d                        ! determine 1d histogram index and extend
!            ===========                        !    histogram if necessary
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine extendHis1d( value, ibox,  &
!             ===========    i     o
                           hist, hbeg, wbox, imin, imax )
!                           io    io    i     io    io
!
!     Determines box indices in a 1d histogram and, if necessary, rearranges the histogram so
!     that value is contained in it.
!
!  Arguments:
!
      real(REAL8), intent(in) :: value             ! value to be scored
      integer(INT4), intent(out) :: ibox           ! box index for value
      real(REAL8), allocatable, intent(inout) :: hist(:,:)  ! (nhist,nbox) histograms
      real(REAL8), intent(inout):: hbeg            ! value corresponding to the lower
                                                   !    boundary of the histogram
      real(REAL8), intent(in) :: wbox              ! box width of the histogram
      integer(INT4), intent(inout) :: imin         ! index of first non-empty histogram box
      integer(INT4), intent(inout) :: imax         ! index of last non-empty histogram box
!
!  Local Variables:
!
      logical move_right(2)
      integer(INT4) :: nbox, mbox, newshape(2)
!
!  Start of Executable Code:
!  ------------------------
!
      ibox = ceiling( (value - hbeg) / wbox )

      nbox = size( hist, 2 )

      if( nbox == 0 ) then                      ! first score
         hbeg = hbeg + (ibox-1)*wbox
         nbox = 1
         move_right(2) = .false.
         ibox = 1
      else if( ibox < 1 ) then                  ! score below lower histogram limit
         mbox = 1 - ibox
         hbeg = hbeg - mbox*wbox
         nbox = nbox + mbox
         move_right(2) = .true.
         ibox = 1
      else if( ibox > nbox ) then               ! score above upper histogram limit
         move_right(2) = .false.
         nbox = ibox
      end if

      if( nbox /= size( hist, 2 ) ) then
         newshape(1) = size( hist, 1 )
         newshape(2) = nbox
         move_right(1) = .false.
         call changeShape( hist, newshape(:), 'hist', 'extendHis1d', move_right_=move_right(:) )
!             ===========   io        i         i           i                    i
      end if

      imin = 1
      imax = nbox

   end subroutine extendHis1d
!
!=======================================================================
!
   subroutine extendHis2d( value, ibox,  &
!             ===========    i     o
                           hist, hbeg, wbox, imin, imax )
!                           io    io    i     io    io
!
!     Determines box indices in a 2d histogram and, if necessary, rearranges the histogram so
!     that value is contained in it.
!
!  Arguments:
!
      real(REAL8), intent(in) :: value(2)          ! values to be scored
      integer(INT4), intent(out) :: ibox(2)        ! box indices for value
      real(REAL8), allocatable, intent(inout) :: hist(:,:,:)   ! (nhist,nbox1,nbox2) histograms
      real(REAL8), intent(inout):: hbeg(2)         ! values corresponding to the lower
                                                   !    boundaries of the histogram
      real(REAL8), intent(in) :: wbox(2)           ! box widths of the histogram
      integer(INT4), intent(inout) :: imin(2)      ! indices of first non-empty histogram boxes
      integer(INT4), intent(inout) :: imax(2)      ! indices of last non-empty histogram boxes
!
!  Local Variables:
!
      logical :: move_right(3)
      integer(INT4) :: nbox(2), mbox(2), oldshape(3), newshape(3)
!
!  Start of Executable Code:
!  ------------------------
!
      ibox(:) = ceiling( (value(:) - hbeg(:)) / wbox(:) )

      nbox(:) = [ size( hist, 2 ), size( hist, 3 ) ]
!
      where( nbox(:) == 0 )                     ! first score
         hbeg(:) = hbeg(:) + (ibox(:)-1)*wbox(:)
         nbox(:) = 1
         move_right(2:) = .false.
         ibox(:) = 1
      elsewhere( ibox(:) < 1 )                  ! score below lower histogram limit
         mbox(:) = 1 - ibox(:)
         hbeg(:) = hbeg(:) - mbox(:)*wbox(:)
         nbox(:) = nbox(:) + mbox(:)
         move_right(2:) = .true.
         ibox(:) = 1
      elsewhere( ibox(:) > nbox(:) )            ! score above upper histogram limit
         move_right(2:) = .false.
         nbox(:) = ibox(:)
      end where

      oldshape(:) = shape(hist)
      if( any( nbox(:) /= oldshape(2:) ) ) then
         newshape(1) = oldshape(1)
         newshape(2:) = nbox(:)
         move_right(1) = .false.
         call changeShape( hist, newshape(:), 'hist', 'extendHis2d', move_right_=move_right(:) )
!             ===========   io        i         i           i                    i
      end if

      imin(:) = 1
      imax(:) = nbox(:)

   end subroutine extendHis2d
!
!=======================================================================
!
   subroutine extendHis3d( value, ibox,  &
!             ===========    i     o
                           hist, hbeg, wbox, imin, imax)
!                           io    io    i     io    io
!
!     Determines box indices in a 3d histogram and, if necessary, rearranges the histogram so
!     that value is contained in it.
!
!  Arguments:
!
      real(REAL8), intent(in) :: value(3)          ! values to be scored
      integer(INT4), intent(out) :: ibox(3)        ! box indices for value
      real(REAL8), allocatable, intent(inout) :: hist(:,:,:,:)
                                                   ! (nhist,nbox1,nbox2,nbox3) histograms
      real(REAL8), intent(inout):: hbeg(3)         ! values corresponding to the lower
                                                   !    boundaries of the histogram
      real(REAL8), intent(in) :: wbox(3)           ! box widths of the histogram
      integer(INT4), intent(inout) :: imin(3)      ! indices of first non-empty histogram boxes
      integer(INT4), intent(inout) :: imax(3)      ! indices of last non-empty histogram boxes
!
!  Local Variables:
!
      logical :: move_right(4)
      integer(INT4) :: nbox(3), mbox(3), oldshape(4), newshape(4)
!
!  Start of Executable Code:
!  ------------------------
!
      ibox(:) = ceiling( (value(:) - hbeg(:)) / wbox(:) )

      nbox(:) = [ size( hist, 2 ), size( hist, 3 ), size( hist, 4 ) ]
!
      where( nbox(:) == 0 )                     ! first score
         hbeg(:) = hbeg(:) + (ibox(:)-1)*wbox(:)
         nbox(:) = 1
         move_right(2:) = .false.
         ibox(:) = 1
      elsewhere( ibox(:) < 1 )                  ! score below lower histogram limit
         mbox(:) = 1 - ibox(:)
         hbeg(:) = hbeg(:) - mbox(:)*wbox(:)
         nbox(:) = nbox(:) + mbox(:)
         move_right(2:) = .true.
         ibox(:) = 1
      elsewhere( ibox(:) > nbox(:) )            ! score above upper histogram limit
         move_right(2:) = .false.
         nbox(:) = ibox(:)
      end where

      oldshape(:) = shape(hist)
      if( any( nbox(:) /= oldshape(2:) ) ) then
         newshape(1) = oldshape(1)
         newshape(2:) = nbox(:)
         move_right(1) = .false.
         call changeShape( hist, newshape(:), 'hist', 'extendHis3d', move_right_=move_right(:) )
!             ===========   io        i         i           i                    i
      end if

      imin(:) = 1
      imax(:) = nbox(:)

   end subroutine extendHis3d
!
!=======================================================================
!
 end module Stat_His_Extend_mod

module Stat_His_Adapt_mod
!
!  Contains subroutines that adapt histogram as to accomodate new scores without increasing
!  the number of histogram boxes. THESE SUBROUTINES ALSO RETURN THE INDEX OF A SCORE.
!
!  Used Modules:
!  ============
!
   use Const_mod, only: INT4, REAL8,  &
                        EPSMA, ZERO, ONE, TWO, HALF, THIRD
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
   public :: adaptHis1d                         ! determine 1d histogram index and adapt
!            ==========                         !    histogram if necessary
   public :: adaptHis2d                         ! determine 1d histogram index and adapt
!            ==========                         !    histogram if necessary
   public :: adaptHis3d                         ! determine 1d histogram index and adapt
!            ==========                         !    histogram if necessary
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine adaptHis1d( value, ibox,  &
!             ==========    i     o
                          hist, hbeg, wbox, imin, imax )
!                          io    io    io    io    io
!
!     Determines box indices in a 1d histogram and, if necessary, rearranges the histogram so
!     that value is contained in it.
!
!  Remarks
!
!     1) The dimensions of hist remain fixed, i.e., the number of histogram boxes is unchanged
!        on exit
!     2) imax(1)=0 on entry indicates the histogram to be empty
!     3) If value is outside the histogram, the histogram is tried to be shifted as to contain
!        the new value. If this is not possible without loss of data, the histogram is
!        coarsened by increasing the box widths by a factor of 2 until both old data and the
!        new value are contained in the histogram.
!     4) Coarsening is done in a way to guarantee that the box boundary closest to the
!        origin is maintained. For the first score, it is guaranteed that hbeg remains a
!        box boundary
!     5) hist is NOT incremented (and cannot be, since no increment is passed).
!
!  Arguments:
!
      real(REAL8), intent(in) :: value             ! value to be scored
      integer(INT4), intent(out) :: ibox           ! box index for value
      real(REAL8), intent(inout) :: hist(:,:)      ! (nhist,nbox) histograms
      real(REAL8), intent(inout):: hbeg            ! value corresponding to the lower
                                                   !    boundary of the histogram
      real(REAL8), intent(inout) :: wbox           ! box width of the histogram
      integer(INT4), intent(inout) :: imin         ! index of first non-empty histogram box
      integer(INT4), intent(inout) :: imax         ! index of last non-empty histogram box
!
!  Local Variables:
!
      logical coarsen
      integer(INT4) :: ishift, iodd, jbox
      integer(INT4) :: nbox, i1, i2
      real(REAL8) :: hist_(size( hist, 1 )), xxxdeb
!
!  Start of Executable Code:
!  ------------------------
!
      nbox = size( hist, 2 )
!
!  For first score, move histogram such that value is approximately in center of histogram
!
      if( imax == 0 )  &
            hbeg = hbeg + nint( (value-hbeg)/wbox - nbox/2 ) * wbox
!
!  Coarsening loop
!
      do

         ibox = ceiling( (value - hbeg) / wbox )
!
!  Test whether shifting or coarsening of histogram is necessary
!
         if( ibox < 1 ) then
            ishift = 1 - ibox
            coarsen = ( ishift > nbox - imax )
         else if( ibox > nbox ) then
            ishift = ibox - nbox
            coarsen = ( ishift > imin - 1 )
         else
            ishift = 0
            coarsen = .false.
         end if

         if( .not.coarsen )  exit
!
!  Double box size if necessary
!
         ! iodd is needed to maintain origin as box boundary
         iodd = abs( mod( nint( hbeg/wbox ), 2 ) )
         ! update histogram values
         do i1 = 1, nbox                        ! loop over new histogram boxes
            jbox = 2*i1 - iodd - 1              ! first of the two indices in the old histogram
            hist_(:) = ZERO
            do i2 = max(imin,jbox), min(imax,jbox+1)
               hist_(:) = hist_(:) + HALF*hist(:,i2)
            end do
            hist(:,i1) = hist_(:)
         end do
         ! adjust control variables
         hbeg = hbeg - iodd*wbox
         wbox = TWO*wbox
         imin = (imin + 1 + iodd) / 2
         imax = (imax + 1 + iodd) / 2

      end do                                    ! end of coarsening loop
!
!  Shift histogram if necessary
!
      if( ibox < 1 ) then
!
!  Shift histogram to the left
!
         ! modify ishift so occupied boxes are in the middle of hist
         ishift = (ishift + nbox - imax) / 2
         ! update histogram values
         do i1 = imax, imin, -1                 ! loop over old histogram boxes
            i2 = i1 + ishift                    ! index in new histogram
            hist(:,i2) = hist(:,i1)
            hist(:,i1) = ZERO                   ! i1<i2, therfore hist(:,i1) will be overwritten
         end do                                 !    later until i2 is exhausted
         ! adjust control variables
         hbeg = hbeg - ishift*wbox
         imin = imin + ishift
         imax = imax + ishift

         ibox = ibox + ishift

      else if( ibox > nbox ) then
!
!  Shift histogram to the right
!
         ! modify ishift so occupied boxes are in the middle of hist
         ishift = (ishift + imin - 1) / 2
         ! update histogram
         do i1 = imin, imax                     ! loop over old histogram boxes
            i2 = i1 - ishift                    ! index in new histogram
            hist(:,i2) = hist(:,i1)
            hist(:,i1) = ZERO                   ! i1>i2, therfore hist(:,i1) will be overwritten
         end do                                 !    later until i2 is exhausted
         ! adjust control variables
         hbeg = hbeg + ishift*wbox
         imin = imin - ishift
         imax = imax - ishift

         ibox = ibox - ishift

      end if

! for debug purposes
      if( nint(hbeg/wbox)+ibox == 0 ) then
         xxxdeb=ZERO
      endif
! end debug

!
!  Adjust minimum and maximum occupied-box index
!
      imin = min( imin, ibox )
      imax = max( imax, ibox )

   end subroutine adaptHis1d
!
!=======================================================================
!
   subroutine adaptHis2d( value, ibox,  &
!             ==========    i     o
                          hist, hbeg, wbox, imin, imax)
!                          io    io    io    io    io
!
!     Determines box indices in a 2d histogram and, if necessary, rearranges the histogram so
!     that value is contained in it.
!
!  Remarks
!
!     1) The dimensions of hist remain fixed, i.e., the number of histogram boxes is unchanged
!        on exit
!     2) imax(1)=0 on entry indicates the histogram to be empty
!     3) If value is outside the histogram, the histogram is tried to be shifted as to contain
!        the new value. If this is not possible without loss of data, the histogram is
!        coarsened by increasing the box widths by a factor of 2 until both old data and the
!        new value are contained in the histogram.
!     4) Coarsening is done in a way to guarantee that the box boundary closest to the
!        origin is maintained. For the first score, it is guaranteed that hbeg remains a
!        box boundary
!     5) hist is NOT incremented (and cannot be, since no increment is passed).
!
!  Arguments:
!
      real(REAL8), intent(in) :: value(2)          ! values to be scored
      integer(INT4), intent(out) :: ibox(2)        ! box indices for value
      real(REAL8), intent(inout) :: hist(:,:,:)    ! (nhist,nbox1,nbox2) histograms
      real(REAL8), intent(inout):: hbeg(2)         ! values corresponding to the lower
                                                   !    boundaries of the histogram
      real(REAL8), intent(inout) :: wbox(2)        ! box widths of the histogram
      integer(INT4), intent(inout) :: imin(2)      ! indices of first non-empty histogram boxes
      integer(INT4), intent(inout) :: imax(2)      ! indices of last non-empty histogram boxes
!
!  Local Variables:
!
      logical coarsen
      integer(INT4) :: i, ishift, iodd, j1, j2, jbox, jj
      integer(INT4) :: nbox(2), i1(2), i2(2), n1(2), n2(2)
      real(REAL8) :: hist_( size( hist, 1 )), xxxdeb
!
!  Start of Executable Code:
!  ------------------------
!
      nbox(:) = [ size( hist, 2 ), size( hist, 3 ) ]
!
!  For first score, move histogram such that value is approximately in center of histogram
!
      if( imax(1) == 0 )  &
            hbeg(:) = hbeg(:) + nint( (value(:)-hbeg(:))/wbox(:) - nbox(:)/2 ) * wbox(:)

! Loop over dimensions
! --------------------

      do i = 1, 2
!
!  Coarsening loop
!
         do

            ibox(i) = ceiling( (value(i) - hbeg(i)) / wbox(i) )
!
!  Test whether shifting or coarsening of histogram is necessary
!
            if( ibox(i) < 1 ) then
               ishift = 1 - ibox(i)
               coarsen = ( imax(i) + ishift > nbox(i) )
            else if( ibox(i) > nbox(i) ) then
               ishift = ibox(i) - nbox(i)
               coarsen = ( imin(i) - ishift < 1 )
            else
               ishift = 0
               coarsen = .false.
            end if

            if( .not.coarsen )  exit
!
!  Double box size if necessary
!
            ! iodd is needed to maintain origin as box boundary
            iodd = abs( mod( nint( hbeg(i)/wbox(i) ), 2 ) )
            ! n1, n2 are the indices to be updated (currently: all)
            n1(:) = 1
            n2(:) = nbox(:)
            ! update histogram values
            do j1 = n1(1), n2(1)
               i1(1) = j1                       ! index in new histogram
               i2(1) = j1                       ! index in old histogram
               do j2 = n1(2), n2(2)
                  i1(2) = j2
                  i2(2) = j2
                  jbox = 2*i1(i) - iodd - 1
                  hist_(:) = ZERO
                  do jj = max(imin(i),jbox), min(imax(i),jbox+1)
                     i2(i) = jj
                     hist_(:) = hist_(:) + HALF*hist(:,i2(1),i2(2))
                  end do
                  hist(:,i1(1),i1(2)) = hist_(:)
               end do
            end do
            ! adjust control variables
            hbeg(i) = hbeg(i) - iodd*wbox(i)
            wbox(i) = TWO*wbox(i)
            imin(i) = (imin(i) + 1 + iodd) / 2
            imax(i) = (imax(i) + 1 + iodd) / 2

         end do                                 ! end of coarsening loop
!
!  Shift histogram if necessary
!
         if( ibox(i) < 1 ) then
!
!  Shift histogram to the left
!
            ! modify ishift so occupied boxes are in the middle of hist
            ishift = (ishift + nbox(i) - imax(i)) / 2
            ! n1, n2 are the indices to be updated
            n1(:) = 1
            n2(:) = nbox(:)
            n1(i) = imin(i)
            n2(i) = imax(i)
            ! update histogram values
            do j1 = n2(1), n1(1), -1
               i1(1) = j1                       ! index in new histogram
               i2(1) = j1                       ! index in old histogram
               do j2 = n2(2), n1(2), -1
                  i1(2) = j2
                  i2(2) = j2
                  i2(i) = i1(i) + ishift
                  hist(:,i2(1),i2(2)) = hist(:,i1(1),i1(2))
                  hist(:,i1(1),i1(2)) = ZERO
               end do
            end do
            ! adjust control variables
            hbeg(i) = hbeg(i) - ishift*wbox(i)
            imin(i) = imin(i) + ishift
            imax(i) = imax(i) + ishift

            ibox(i) = ibox(i) + ishift

         else if( ibox(i) > nbox(i) ) then
!
!  Shift histogram to the right
!
            ! modify ishift so occupied boxes are in the middle of hist
            ishift = (ishift + imin(i) - 1) / 2
            ! n1, n2 are the indices to be updated
            n1(:) = 1
            n2(:) = nbox(:)
            n1(i) = imin(i)
            n2(i) = imax(i)
            ! update histogram
            do j1 = n1(1), n2(1)
               i1(1) = j1
               i2(1) = j1
               do j2 = n1(2), n2(2)
                  i1(2) = j2
                  i2(2) = j2
                  i2(i) = i1(i) - ishift
                  hist(:,i2(1),i2(2)) = hist(:,i1(1),i1(2))
                  hist(:,i1(1),i1(2)) = ZERO
               end do
            end do
            ! adjust control variables
            hbeg(i) = hbeg(i) + ishift*wbox(i)
            imin(i) = imin(i) - ishift
            imax(i) = imax(i) - ishift

            ibox(i) = ibox(i) - ishift

         end if

! for debug purposes
         if( nint(hbeg(i)/wbox(i))+ibox(i)==0 ) then
            xxxdeb=ZERO
         endif
! end debug

!
!  Adjust minimum and maximum occupied-box index
!
         imin(i) = min( imin(i), ibox(i) )
         imax(i) = max( imax(i), ibox(i) )

      end do

   end subroutine adaptHis2d
!
!=======================================================================
!
   subroutine adaptHis3d( value, ibox,  &
!             ==========    i     o
                          hist, hbeg, wbox, imin, imax)
!                          io    io    io    io    io
!
!     Determines box indices in a 3d histogram and, if necessary, rearranges the histogram so
!     that value is contained in it.
!
!  Remarks
!
!     1) The dimensions of hist remain fixed, i.e., the number of histogram boxes is unchanged
!        on exit
!     2) imax(1)=0 on entry indicates the histogram to be empty
!     3) If value is outside the histogram, the histogram is tried to be shifted as to contain
!        the new value. If this is not possible without loss of data, the histogram is
!        coarsened by increasing the box widths by a factor of 2 until both old data and the
!        new value are contained in the histogram.
!     4) Coarsening is done in a way to guarantee that the box boundary closest to the
!        origin is maintained. For the first score, it is guaranteed that hbeg remains a
!        box boundary
!     5) hist is NOT incremented (and cannot be, since no increment is passed).
!
!  Arguments:
!
      real(REAL8), intent(in) :: value(3)          ! values to be scored
      integer(INT4), intent(out) :: ibox(3)        ! box indices for value
      real(REAL8), intent(inout) :: hist(:,:,:,:)  ! (nhist,nbox1,nbox2,nbox3) histograms
      real(REAL8), intent(inout):: hbeg(3)         ! values corresponding to the lower
                                                   !    boundaries of the histogram
      real(REAL8), intent(inout) :: wbox(3)        ! box widths of the histogram
      integer(INT4), intent(inout) :: imin(3)      ! indices of first non-empty histogram boxes
      integer(INT4), intent(inout) :: imax(3)      ! indices of last non-empty histogram boxes
!
!  Local Variables:
!
      logical coarsen
      integer(INT4) :: i, ishift, iodd, j1, j2, j3, jbox, jj
      integer(INT4) :: nbox(3), i1(3), i2(3), n1(3), n2(3)
      real(REAL8) :: hist_(size( hist, 1 )), xxxdeb
!
!  Start of Executable Code:
!  ------------------------
!
      nbox(:) = [ size( hist, 2 ), size( hist, 3 ), size( hist, 4 ) ]
!
!  For first score, move histogram such that value is approximately in center of histogram
!
      if( imax(1) == 0 )  &
            hbeg(:) = hbeg(:) + nint( (value(:)-hbeg(:))/wbox(:) - nbox(:)/2 ) * wbox(:)

! Loop over dimensions
! --------------------

      do i = 1, 3
!
!  Coarsening loop
!
         do

            ibox(i) = ceiling( (value(i) - hbeg(i)) / wbox(i) )
!
!  Test whether shifting or coarsening of histogram is necessary
!
            if( ibox(i) < 1 ) then
               ishift = 1 - ibox(i)
               coarsen = ( imax(i) + ishift > nbox(i) )
            else if( ibox(i) > nbox(i) ) then
               ishift = ibox(i) - nbox(i)
               coarsen = ( imin(i) - ishift < 1 )
            else
               ishift = 0
               coarsen = .false.
            end if

            if( .not.coarsen )  exit
!
!  Double box size if necessary
!
            ! iodd is needed to maintain origin as box boundary
            iodd = abs( mod( nint( hbeg(i)/wbox(i) ), 2 ) )
            ! n1, n2 are the indices to be updated (currently: all)
            n1(:) = 1
            n2(:) = nbox(:)
            ! update histogram values
            do j1 = n1(1), n2(1)
               i1(1) = j1                       ! index in new histogram
               i2(1) = j1                       ! index in old histogram
               do j2 = n1(2), n2(2)
                  i1(2) = j2
                  i2(2) = j2
                  do j3 = n1(3), n2(3)
                     i1(3) = j3
                     i2(3) = j3
                     jbox = 2*i1(i) - iodd - 1
                     hist_(:) = ZERO
                     do jj = max(imin(i),jbox), min(imax(i),jbox+1)
                        i2(i) = jj
                        hist_(:) = hist_(:) + HALF*hist(:,i2(1),i2(2),i2(3))
                     end do
                     hist(:,i1(1),i1(2),i1(3)) = hist_(:)
                  end do
               end do
            end do
            ! adjust control variables
            hbeg(i) = hbeg(i) - iodd*wbox(i)
            wbox(i) = TWO*wbox(i)
            imin(i) = (imin(i) + 1 + iodd) / 2
            imax(i) = (imax(i) + 1 + iodd) / 2

         end do                                 ! end of coarsening loop
!
!  Shift histogram if necessary
!
         if( ibox(i) < 1 ) then
!
!  Shift histogram to the left
!
            ! modify ishift so occupied boxes are in the middle of hist
            ishift = (ishift + nbox(i) - imax(i)) / 2
            ! n1, n2 are the indices to be updated
            n1(:) = 1
            n2(:) = nbox(:)
            n1(i) = imin(i)
            n2(i) = imax(i)
            ! update histogram values
            do j1 = n2(1), n1(1), -1
               i1(1) = j1                       ! index in new histogram
               i2(1) = j1                       ! index in old histogram
               do j2 = n2(2), n1(2), -1
                  i1(2) = j2
                  i2(2) = j2
                  do j3 = n2(3), n1(3), -1
                     i1(3) = j3
                     i2(3) = j3
                     i2(i) = i1(i) + ishift
                     hist(:,i2(1),i2(2),i2(3)) = hist(:,i1(1),i1(2),i1(3))
                     hist(:,i1(1),i1(2),i1(3)) = ZERO
                  end do
               end do
            end do
            ! adjust control variables
            hbeg(i) = hbeg(i) - ishift*wbox(i)
            imin(i) = imin(i) + ishift
            imax(i) = imax(i) + ishift

            ibox(i) = ibox(i) + ishift

         else if( ibox(i) > nbox(i) ) then
!
!  Shift histogram to the right
!
            ! modify ishift so occupied boxes are in the middle of hist
            ishift = (ishift + imin(i) - 1) / 2
            ! n1, n2 are the indices to be updated
            n1(:) = 1
            n2(:) = nbox(:)
            n1(i) = imin(i)
            n2(i) = imax(i)
            ! update histogram
            do j1 = n1(1), n2(1)
               i1(1) = j1
               i2(1) = j1
               do j2 = n1(2), n2(2)
                  i1(2) = j2
                  i2(2) = j2
                  do j3 = n1(3), n2(3)
                     i1(3) = j3
                     i2(3) = j3
                     i2(i) = i1(i) - ishift
                     hist(:,i2(1),i2(2),i2(3)) = hist(:,i1(1),i1(2),i1(3))
                     hist(:,i1(1),i1(2),i1(3)) = ZERO
                  end do
               end do
            end do
            ! adjust control variables
            hbeg(i) = hbeg(i) + ishift*wbox(i)
            imin(i) = imin(i) - ishift
            imax(i) = imax(i) - ishift

            ibox(i) = ibox(i) - ishift

         end if

! for debug purposes
         if( nint(hbeg(i)/wbox(i))+ibox(i)==0 ) then
            xxxdeb=ZERO
         endif
! end debug

!
!  Adjust minimum and maximum occupied-box index
!
         imin(i) = min( imin(i), ibox(i) )
         imax(i) = max( imax(i), ibox(i) )

      end do

   end subroutine adaptHis3d
!
!=======================================================================
!
 end module Stat_His_Adapt_mod

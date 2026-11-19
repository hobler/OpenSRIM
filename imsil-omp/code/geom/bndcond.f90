module Geom_BndCond_mod
!
!  Manages the boundary conditions
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: REAL8, INT4,  &
                              ZERO, ONE, TWO,  &
                              ENUM

   use Param_Geom_mod,  only: geom
   use Param_Out_mod,   only: out
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
   public :: applyBndCond1d                  ! apply boundary condition to z position
!            ==============
   public :: applyBndCond2d                  ! apply boundary condition to x and z position
!            ==============
   public :: applyBndCond3d                  ! apply boundary condition to x, y, and z position
!            ==============
!
!  Private Procedures:
!  ==================
!
   private :: applyBndCond                   ! apply boundary condition in one direction
!             ============
!
contains
!
!=======================================================================
!
   subroutine applyBndCond1d( z, lgeom_, allow_score_ )
!             ============== io   (i)        (o)
!
!     Apply the boundary conditions along the z coordinate. Depending on lgeom_, the mapping
!     is done into the geometry or score domain.
!
!     If lgeom_ is .false. on entry, allow_score_ indicates on exit whether a score is allowed.
!     allow_score_ must be present if lgeom_ is .false.
!
!  Arguments:
!
      real(REAL8), intent(inout) :: z                 ! On entry: position before mapping
                                                      ! On exit: position after mapping
      logical, optional, intent(in) :: lgeom_         ! (default: .true.) indicates mapping into
                                                      !    geometry or score domain
      logical, optional, intent(out) :: allow_score_  ! indicates whether score is allowed
!
!  Local Variables:
!
      logical :: lgeom, allow_scorez
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( lgeom_ ) ) then
         lgeom = lgeom_
      else
         lgeom = .true.
      end if

      call applyBndCond( z, geom%bndtypez(:), geom%zlim(:), out%lextendz(:), lgeom,  &
!          ============ io          i              i              i            i
                         allow_scorez )
!                             o

      if( present( allow_score_ ) )  allow_score_ = allow_scorez

   end subroutine applyBndCond1d
!
!=======================================================================
!
   subroutine applyBndCond2d( x, z, lgeom_, allow_score_ )
!             ============== io io   (i)        (o)
!
!     Apply the boundary conditions along the x and z coordinates. Depending on lgeom_, the
!     mapping is done into the geometry or score domain.
!
!     If lgeom_ is .false. on entry, allow_score_ indicates on exit whether a score is allowed.
!     allow_score_ must be present if lgeom_ is .false.
!
!  Arguments:
!
      real(REAL8), intent(inout) :: x                 ! On entry: position before mapping
                                                      ! On exit: position after mapping
      real(REAL8), intent(inout) :: z                 ! On entry: position before mapping
                                                      ! On exit: position after mapping
      logical, optional, intent(in) :: lgeom_         ! (default: .true.) indicates mapping into
                                                      !    geometry or score domain
      logical, optional, intent(out) :: allow_score_  ! indicates whether score is allowed
!
!  Local Variables:
!
      logical :: lgeom, allow_scorex, allow_scorez
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( lgeom_ ) ) then
         lgeom = lgeom_
      else
         lgeom = .true.
      end if

      call applyBndCond( x, geom%bndtypex(:), geom%xlim(:), out%lextendx(:), lgeom,  &
!          ============ io          i              i              i            i
                         allow_scorex )
!                             o
      call applyBndCond( z, geom%bndtypez(:), geom%zlim(:), out%lextendz(:), lgeom,  &
!          ============ io          i              i              i            i
                         allow_scorez )
!                             o

      if( present( allow_score_ ) )  allow_score_ = allow_scorex .and. allow_scorez

   end subroutine applyBndCond2d
!
!=======================================================================
!
   subroutine applyBndCond3d( x, y, z, lgeom_, allow_score_ )
!             ============== io io io   (i)        (o)
!
!     Apply the boundary conditions along the x, y, and z coordinates. The target of the mapping
!     may be the geometry or score domain.
!
!     - The absence of the optional argument allow_score_ indicates that the target is the
!       geometry domain. The result is always valid.
!
!     - The presence of the optional argument allow_score_ indicates that the target is the
!       score domain. The returned value of allow_score_ indicates whether the returned value
!       of z is valid.
!
!  Arguments:
!
      real(REAL8), intent(inout) :: x                 ! On entry: position before mapping
                                                      ! On exit: position after mapping
      real(REAL8), intent(inout) :: y                 ! On entry: position before mapping
                                                      ! On exit: position after mapping
      real(REAL8), intent(inout) :: z                 ! On entry: position before mapping
                                                      ! On exit: position after mapping
      logical, optional, intent(in) :: lgeom_         ! (default: .true.) indicates mapping into
                                                      !    geometry or score domain
      logical, optional, intent(out) :: allow_score_  ! indicates whether score is allowed
!
!  Local Variables:
!
      logical :: lgeom, allow_scorex, allow_scorey, allow_scorez
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( lgeom_ ) ) then
         lgeom = lgeom_
      else
         lgeom = .true.
      end if

      call applyBndCond( x, geom%bndtypex(:), geom%xlim(:), out%lextendx(:), lgeom,  &
!          ============ io          i              i              i            i
                         allow_scorex )
!                             o
      call applyBndCond( y, geom%bndtypey(:), geom%ylim(:), out%lextendy(:), lgeom,  &
!          ============ io          i              i              i            i
                         allow_scorey )
!                             o
      call applyBndCond( z, geom%bndtypez(:), geom%zlim(:), out%lextendz(:), lgeom,  &
!          ============ io          i              i              i            i
                         allow_scorez )
!                             o

      if( present( allow_score_ ) )  &
            allow_score_ = allow_scorex .and. allow_scorey .and. allow_scorez

   end subroutine applyBndCond3d
!
!=======================================================================
!
   subroutine applyBndCond( x, bndtype, xlim, lextend, lgeom, allow_score )
!             ============ io     i      i       i       i         o
!
!     Apply the boundary conditions to map the position x into the geometry or score domain.
!
!     There are a few special cases:
!
!     - At reflective boundaries, first the boundary type of the opposite boundary becomes
!       effective, and then the mapping at the reflective boundary is done.
!
!     - At extending boundaries and with mapping into the score domain, the mapping is labeled
!       invalid if the point is mapped to the domain boundary.
!
!     - No mapping is done at free boundaries
!
!     No mapping is done at free boundaries, so x may end up outside the geometry domain at a
!     free boundary.
!
!  Arguments:
!
      real(REAL8), intent(inout) :: x              ! on entry: coordinate to be mapped
                                                   ! on exit: mapped coordinate
      integer(INT4), intent(in) :: bndtype(2)      ! boundary type
      real(REAL8), intent(in) :: xlim(3)           ! xlim(1): left boundary of domain
                                                   ! xlim(2): right boundary of domain
                                                   ! xlim(3): center of domain (used only
                                                   !    at reflective boundaries
      logical, intent(in) :: lextend(2)            ! flag indicating that no mapping should be
                                                   !    done if target is the score domain
      logical, intent(in) :: lgeom                 ! indicates mapping into geometry or score
                                                   !    domain
      logical, intent(out) :: allow_score          ! flag indicating that returned value
                                                   !    of x may be used in scores.
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      allow_score = .true.
!
!  Treat periodic and extending boundaries
!
      if( x < xlim(1) ) then
         if( bndtype(1) == ENUM%bndtype%periodic .or.  &
             ( bndtype(1) == ENUM%bndtype%reflect .and.  &
               bndtype(2) == ENUM%bndtype%periodic ) ) then
            if( lgeom .or. .not.lextend(1) )  &
                  x = xlim(1) + modulo( x-xlim(1), xlim(2)-xlim(1) )
         else if( bndtype(1) == ENUM%bndtype%extend .or.  &
             ( bndtype(1) == ENUM%bndtype%reflect .and.  &
               bndtype(2) == ENUM%bndtype%extend ) ) then
            if( lgeom ) then
               x = xlim(1)
            else
               allow_score = lextend(1)
            end if
         end if
      else if( x > xlim(2) ) then
         if( bndtype(2) == ENUM%bndtype%periodic .or.  &
             ( bndtype(2) == ENUM%bndtype%reflect .and.  &
               bndtype(1) == ENUM%bndtype%periodic ) ) then
            if( lgeom .or. .not.lextend(2) )  &
                  x = xlim(1) + modulo( x-xlim(1), xlim(2)-xlim(1) )
         else if( bndtype(2) == ENUM%bndtype%extend .or.  &
             ( bndtype(2) == ENUM%bndtype%reflect .and.  &
               bndtype(1) == ENUM%bndtype%extend ) ) then
            if( lgeom ) then
               x = xlim(2)
            else
               allow_score = lextend(2)
            end if
         end if
      end if
!
!  Treat reflective boundaries (only one can be reflective)
!  Note that for cell-based geometry, reflected boundaries have been eliminated
!
      if( .not. geom%lcell ) then
         if( bndtype(1) == ENUM%bndtype%reflect ) then
            if( x < xlim(3) )  x = xlim(1) + xlim(2) - x
         else if( bndtype(2) == ENUM%bndtype%reflect ) then
            if( x > xlim(3) )  x = xlim(1) + xlim(2) - x
         end if
      end if

   end subroutine applyBndCond
!
!=======================================================================
!
   subroutine applyBndCondOld( x, bndtype, xlim, lextend, allow_score_ )
!             =============== io     i      i       i          io
!
!     Apply the boundary conditions to map the position x into the geometry or score domain.
!
!     The absence of the optional argument allow_score_ indicates that the mapping is to be
!     performed into the geometry domain. The result is always valid.
!
!     The presence of the optional argument allow_score_ indicates that the mapping is to be
!     performed into the score domain. The returned value of allow_score_ indicates whether the
!     returned value of z is valid.

!     There are a few special cases:
!
!     - At reflective boundaries, first the boundary type of the opposite boundary becomes
!       effective, and then the mapping at the reflective boundary is done.
!
!     - At extending boundaries and with mapping into the score domain, the mapping is labeled
!       invalid if the point is mapped to the domain boundary.
!
!     - No mapping is done at free boundaries
!
!     No mapping is done at free boundaries, so x may end up outside the geometry domain at a
!     free boundary.
!
!  Arguments:
!
      real(REAL8), intent(inout) :: x                    ! on entry: coordinate to be mapped
                                                         ! on exit: mapped coordinate
      integer(INT4), intent(in) :: bndtype(2)            ! boundary type
      real(REAL8), intent(in) :: xlim(3)                 ! xlim(1): left boundary of domain
                                                         ! xlim(2): right boundary of domain
                                                         ! xlim(3): center of domain (used only
                                                         !    at reflective boundaries
      logical, intent(in) :: lextend(2)                  ! flag indicating that no mapping
                                                         !    should be done if target is the
                                                         !    score domain
      logical, optional, intent(inout) :: allow_score_   ! flag indicating that returned value
                                                         !    of x may be used in scores.
                                                         !    Never set to .true. in this
                                                         !    subroutine.
                                                         !    Presence of allow_score_ indicates
                                                         !    target is the score domain
!
!  Local Variables:
!
      logical :: lscore, lgeom, allow_map(2)
!
!  Start of Executable Code:
!  ------------------------
!
      lscore = present( allow_score_ )
      lgeom = .not.lscore
      allow_map(:) = lgeom .or. .not.lextend(:)
!
!  Treat periodic and extending boundaries
!
      if( x < xlim(1) .and. allow_map(1) ) then
         if( bndtype(1) == ENUM%bndtype%periodic .or.  &
             ( bndtype(1) == ENUM%bndtype%reflect .and.  &
               bndtype(2) == ENUM%bndtype%periodic ) ) then
            x = xlim(1) + modulo( x-xlim(1), xlim(2)-xlim(1) )
         else if( bndtype(1) == ENUM%bndtype%extend .or.  &
             ( bndtype(1) == ENUM%bndtype%reflect .and.  &
               bndtype(2) == ENUM%bndtype%extend ) ) then
            if( lscore ) then
               allow_score_ = .false.
               return
            end if
            x = xlim(1)
         end if
      else if( x > xlim(2) .and. allow_map(2) ) then
         if( bndtype(2) == ENUM%bndtype%periodic .or.  &
             ( bndtype(2) == ENUM%bndtype%reflect .and.  &
               bndtype(1) == ENUM%bndtype%periodic ) ) then
            x = xlim(1) + modulo( x-xlim(1), xlim(2)-xlim(1) )
         else if( bndtype(2) == ENUM%bndtype%extend .or.  &
             ( bndtype(2) == ENUM%bndtype%reflect .and.  &
               bndtype(1) == ENUM%bndtype%extend ) ) then
            if( lscore ) then
               allow_score_ = .false.
               return
            end if
            x = xlim(2)
         end if
      end if
!
!  Treat reflective boundaries (only one can be reflective)
!
      if( bndtype(1) == ENUM%bndtype%reflect ) then
         if( x < xlim(3) .and. allow_map(1) )  x = xlim(1) + xlim(2) - x
      else if( bndtype(2) == ENUM%bndtype%reflect ) then
         if( x > xlim(3) .and. allow_map(2) )  x = xlim(1) + xlim(2) - x
      end if

   end subroutine applyBndCondOld
!
!=======================================================================
!
end module Geom_BndCond_mod

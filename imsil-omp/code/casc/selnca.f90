module Casc_Selnca_mod
!
!  Select the next collision partner in amorphous material
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, ONE,  &
                                 STATE

   use Types_mod,          only: Proj_type
   use Random_mod,         only: ranint

   use Param_Setup_mod,    only: nelem
   use Param_Sn_mod,       only: sn

   use Random_mod,         only: rannum
   use Sn_FFP_mod,         only: evalFFP
   use Sn_FFPcyl_mod,      only: setFFPcylRadius, isInOldFFPcyl
   use Sn_RestFFP_mod,     only: getRestFFP, setRestFFP
   use Sn_Ndira_mod,       only: ndira
   use Xtal_Ipos_mod,      only: initIpos
   use Xtal_Tvib_mod,      only: tvib2
   use Target_mod,         only: getDomainIndex
   use Target_Surface_mod, only: getDistSurface
   use Target_mod,         only: getDomainProperty, getDomainDensities
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
   public :: selnca                 ! select next collision partner in amorphous target
!            ======
!
contains
!
!=======================================================================
!
   subroutine selnca( proj, id, posc, ffp, pmax, lrecoil, recoil, p )
!             ======   i    i    o    o     o       o       o     o
!
!     Select the next collision partner in amorphous material
!
!     The next collision partner is within a circle perpendicular
!     to the direction of motion at the end of the free flight path
!
!     Collisions with atoms outside the layer are not considered.
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      integer(INT4), intent(in) :: id           ! index of domain where the collision partner
                                                !    is to be selected
      real(REAL8), intent(out) :: posc(3)       ! collision point (where the impact parameter
                                                !   is perpendicular to direction of motion)
                                                !   or end of free flight path
      real(REAL8), intent(out) :: ffp           ! free flight path (= distance between
                                                !   proj%pos and posc)
      real(REAL8), intent(out) :: pmax          ! maximum impact parameter
      logical, intent(out) :: lrecoil           ! indicates whether collision partner (recoil)
                                                !   has been found inside the region
      type(Proj_type), intent(out) :: recoil    ! recoil
      real(REAL8), intent(out) :: p             ! impact parameter (= distance between
                                                !   recoil%pos and posc)
!
!  Local Variables:
!
      real(REAL8), parameter :: RANMIN = 1e-10

      integer(INT4) :: ielemr, nelemr, ielem2(nelem)
      real(REAL8) :: ffpmax, rand, cospsi, sinpsi, dir1(3)
      real(real8) :: prob(nelem), vib(nelem)
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine free flight path
!
      call getRestFFP( id, proj%ig, ffp, pmax )
!          ==========  i      i      o    o

      if( ffp <= ZERO ) then                    ! if no valid ffp in list, determine new one

         call evalFFP( proj, id, ffp, pmax )    ! mean ffp
!             =======   i    i    o    o
         if( sn%lstffp ) then                                  ! statistical ffp
            call rannum( rand )
!                ======   o
            rand = max( rand, RANMIN )
            ffp = - ffp * log(rand)
         else                                                  ! random first ffp
            if( proj%ig == 1 .and. all( proj%pos(:) == proj%origin%pos(:) ) ) then
               call rannum( rand )
!                   ======   o
               ffp = ffp * rand
            end if
         end if

         call setRestFFP( id, proj%ig, ffp, pmax )    ! store new value in list
!             ==========  i      i      i    i

      end if
!
!  Set the radius of the ffp cylinder
!
      call setFFPcylRadius( pmax, proj%ig, id )
!          ===============   i       i     i
!
!  Effectively move forward projectile by proj%ffpmin (.ne.0 due to time integral)
!
      ffp = ffp + proj%ffpmin
!
!  If free flight path exceeds its upper limit, set it to its limit and abort collision
!  (RestFFP will be updated in selncp)
!
      ffpmax = min( sn%ffpmax, max( ZERO, proj%distsurf )+sn%pmax )

      if( ffp > ffpmax ) then
         ffp = ffpmax
         posc(:) = proj%pos(:) + ffp * proj%dir(:)
         lrecoil = .false.
         return
      end if
!
!  Determine collision point
!
      posc(:) = proj%pos(:) + ffp * proj%dir(:)
!
!  Calculate the impact parameter 
!
      call rannum( rand )
!          ======   o
      p = pmax * sqrt(rand)
!
!  Determine the position of the collision partner
!
      cospsi = ZERO
      sinpsi = ONE
      call ndira( proj%dir(:), cospsi, sinpsi, dir1(:) )
!          =====      i          i       i       o
      recoil%pos(:) = posc(:) + p * dir1(:)
!
!  Check whether the collision partner is inside the layer
!
      recoil%id = id    ! guess
      call getDomainIndex( recoil%pos(:), recoil%id )
!          ==============        i            io
      call getDistSurface( recoil%pos(:), recoil%distsurf )
!          ==============       i                 o

      if( recoil%id /= id .or. recoil%distsurf < ZERO ) then
         lrecoil = .false.
         return
      end if
!
!  Check whether recoil is within the old ffp-cylinder. If so, abort collision
!
      if( sn%chkffpcyl ) then
         lrecoil = .not.isInOldFFPcyl( recoil%pos(:), proj%ig, id )
!                       =============        i           i     i
         if( .not.lrecoil ) return
      else
         lrecoil = .true.
      end if
!
!  Get target atom index and vibration amplitude
!
      call getDomainDensities( recoil%id, nelemr, ielem2(:), prob_=prob(:) )
!          ==================      i        o        o           o
      call ranint( prob(1:nelemr), ielemr )
!          ======        i           o
      recoil%ielem = ielem2(ielemr)
!
!  Initialize remaining properties of recoil
!  (recoil%e and recoil%dir(:) will be initialized in scatter)
!
      recoil%ig = proj%ig + 1
      recoil%iion = proj%iion
      recoil%nreflect = 0
      recoil%reentered = proj%reentered

      recoil%origin%internal = .true.
      recoil%origin%amorphous = .true.
      recoil%origin%lattice = .false.
      call initIpos( recoil%origin%ipos )
!          ========           o
      recoil%origin%isa = -1
      recoil%origin%pos(:) = recoil%pos(:)
      recoil%ffpmin = ZERO
      recoil%emin = recoil%e

      recoil%virtualparent = proj%virtual
      recoil%virtual = proj%virtual             ! virtual may be modified later
      recoil%ededuct = ZERO
      recoil%state = STATE%SUBTHRESHOLD
!
!  Take into account thermal vibration for recoil too
!
      call getDomainProperty( recoil%id, vib_=vib(:) )
!          =================      i          o
      call tvib2( vib(recoil%ielem), proj%dir(:), recoil%origin%pos(:) )
!          =====         i               i                  io

   end subroutine selnca
!
!=======================================================================
!
end module Casc_Selnca_mod

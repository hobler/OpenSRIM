module Casc_Selncp_mod
!
!  Select the next collision partner and/or determine free flight path
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8
   use Random_mod,      only: rannum

   use Types_mod,       only: Proj_type
   use Param_Setup_mod, only: nr
   use Param_Geom_mod,  only: geom
   use Param_Sn_mod,    only: sn
   use Param_Out_mod,   only: out

   use Target_mod,      only: getOtherNearDomainIndices
   use Target_mod,      only: getDomainProperty
   use Casc_Selncc_mod, only: selncc
   use Casc_Selnca_mod, only: selnca
   use Sn_RestFFP_mod,  only: updateRestFFP
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
   public :: selncp                          ! select next collision partner
!            ======
!
contains
!
!=======================================================================
!
   subroutine selncp( proj, posc, ffp, pmaxse, lrecoil, recoil, p )
!             ======   i     o     o     o        o       o     o
!
!     Select next collision partner, taking into account potential collision
!     partners from several regions. It is possible that no collision partner
!     is found.
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      real(REAL8), intent(out) :: posc(3)       ! position where the collision takes place
                                                !    (= where the impact parameter intersects
                                                !    the direction of motion) or 
                                                !    end of free flight path
      real(REAL8), intent(out) :: ffp           ! free flight path (= distance between
                                                !    proj%pos and posc)
      real(REAL8), intent(out) :: pmaxse        ! maximum impact parameter for Se
      logical, intent(out) :: lrecoil           ! flag indicating whether a collision partner
                                                !    (= recoil) has been found
      type(Proj_type), intent(out) :: recoil    ! collision partner (recoil) if any
      real(REAL8), intent(out) :: p             ! impact parameter (= distance between posc
                                                !    and recoil%pos)
!
!  Local Variables:
!
      logical, save :: first = .true.
      logical :: lrecoil2, lxtal
      integer(INT4) :: id, idn, ndn !, nregions
      integer(INT4), allocatable, save :: idnear(:)   ! save needed to save allocation status
      real(REAL8) :: rand
      real(REAL8) :: posc2(3), ffp2, p2, pmax2
      type(Proj_type) :: recoil2
!$omp threadprivate( first, idnear )
!
!  Start of Executable Code:
!  ------------------------
!
!  Allocate memory for near regions (may be reallocated in getOtherNearDomainIndices)
!
      if( first ) then
         first =.false.
         if( geom%lcell ) then
            allocate( idnear( 1 ) )
         else
            allocate( idnear( max( 1, nr-1 ) ) )
         end if
      end if
!
!  We need a ffp to start with ... take that of the nearest domain
!
!  Determine collision partner in domain proj%id. Note that it is possible that no collision
!  partner in that domain is found (lrecoil=.false.), and that recoil and p are not set.
!
      call getDomainProperty( proj%id, lxtal_=lxtal )
!          =================     i          o

      if( lxtal ) then

         call selncc( proj, proj%id, posc, ffp, pmaxse, lrecoil, recoil, p )
!             ======   i       i      o     o     o        o       o     o
      else

         call selnca( proj, proj%id, posc, ffp, pmaxse, lrecoil, recoil, p )
!             ======   i       i      o     o     o        o       o     o
      end if
!
!  Get other nearby regions, i.e. regions with potential collision partners with ffp smaller 
!  than the above
!
      call getOtherNearDomainIndices( proj, posc, ffp, sn%pmax, ndn, idnear )
!          =========================   i     i     i      i      o     o
!
!  Check the nearby regions for additional collision partners and choose smallest free 
!  flight path. Note that pmaxse must correspond to proj%id, so it must not be updated
!
      do idn = 1, ndn

         id = idnear(idn)

         call getDomainProperty( id, lxtal_=lxtal )
!             =================  i        o

         if( lxtal ) then

            call selncc( proj, id, posc2, ffp2, pmax2, lrecoil2, recoil2, p2 )
!                ======   i    i     o     o      o       o         o     o
         else

            call selnca( proj, id, posc2, ffp2, pmax2, lrecoil2, recoil2, p2 )
!                ======   i    i     o     o      o       o         o     o
         end if

         if( ffp2 < ffp ) then
            posc(:) = posc2(:)
            ffp = ffp2
            lrecoil = lrecoil2
            recoil = recoil2
            p = p2
         end if

      end do
!
!  Consider surface roughness layer
!
!  Note: this also guarrantees that no recoil can start outside the target, even in case of
!        a dynamic simulation with a regular grid, where surface cells can extend into the
!        vacuum
!
      if( lrecoil ) then
         if( recoil%distsurf < geom%wrough ) then
            call rannum( rand )
!                ======   o
            if( recoil%distsurf < rand*geom%wrough )  lrecoil = .false.
         end if
      end if
!
!  Update old free flight path list
!
      call updateRestFFP( proj%ig, ffp )
!          =============     i      i

   end subroutine selncp
!
!=======================================================================
!
end module Casc_Selncp_mod

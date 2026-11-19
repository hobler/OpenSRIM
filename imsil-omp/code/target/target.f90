module Target_mod
!
!  Provide procedures for geometry and material related operations on the target
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO

   use Types_mod,             only: Proj_type

   use Param_Setup_mod,       only: ndim, nelem
   use Param_Geom_mod,        only: geom
   use Param_Ion_mod,         only: ion
   use Param_Mat_mod,         only: mat
   use Param_Out_mod,         only: out

   use Target_Region_mod,     only: initRegion, getRegionIndex, getOtherNearRegionIndices
   use Target_Cell_mod,       only: initCell, getCellIndex, getOtherNearCellIndices,  &
                                    getProbCell, getDensCell, getVibCell,  &
                                    getXtalCell, getEdCell, getEsurfCell
   use Geom_Polygon_old_mod,  only: writePolygons
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
   public :: initTarget                ! initialize the target
!            ==========
   public :: getDomainIndex            ! determine domain nearest to given point and whether
!            ==============            !    point is outside region
   public :: getOtherNearDomainIndices ! determine domains other than nearest domain that
!            ========================= !    intersect with the free flight path cylinder
   public :: getDomainDensities        ! get list of densities and related properties of a
!            ==================        !    domain
   public :: getDomainProperty         ! get basic properties of a domain
!            =================
!
!  Program History:
!  ===============
!
!   8-May-2013  G.Hobler  creation
!  13-May-2013  G.Hobler  include procedures from target_geom_mod and target_mat_mod
!
contains
!
!=======================================================================
!
   subroutine initTarget( )
!             ==========
!
!     Initialize the target
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: np
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcell ) then
         call initCell( )
!		        ========
      else
         call initRegion( )
!             ==========
      end if
!
!  Write the surface polygon
!
      if( ndim == 2 .and. .not.geom%usesdf .and. out%lsurf ) then
         np = geom%region(0)%np + 1
         call writePolygons( [np], geom%region(0)%px(1:np), geom%region(0)%pz(1:np), ZERO )
!             =============    i           i                        i                 i
      end if

   end subroutine initTarget
!
!=======================================================================
!
   subroutine getDomainIndex( pos, id_, outside_ )
!             ==============   i   io      o
!
!     Determine the nearest region to the point pos
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)            ! position to be checked
      integer(INT4), optional, intent(inout) :: id_! domain index
      logical, optional, intent(out) :: outside_   ! flag indicating whether pos(:) is outside
                                                   !    all regions (thus id_ being only the
                                                   !    region nearest to id_)
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
     if( geom%lcell ) then
        call getCellIndex( pos, id_, outside_ )
!            ============   i   io     o
      else
         call getRegionIndex( pos, id_, outside_ )
!             ==============   i   io     o
      end if

   end subroutine getDomainIndex
!
!=======================================================================
!
   subroutine getOtherNearDomainIndices( proj, posc, ffp, pmax, ndn, idnear )
!             =========================   i     i     i    i     o     o
!
!     Determine regions other than proj%id that are within the free flight path cylinder
!     (= cylinder with height ffp and radius pmax)
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      real(REAL8), intent(in) :: posc(3)        ! position of end of free flight path
      real(REAL8), intent(in) :: ffp            ! free flight path (= height of cylinder)
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter (= radius of cyl.)
      integer(INT4), intent(out) :: ndn         ! number of domains found
      integer(INT4), allocatable, intent(inout) :: idnear(:)
                                                ! indices of domains found
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcell ) then
         call getOtherNearCellIndices( proj, posc, ffp, pmax, ndn, idnear )
!             =======================   i     i     i    i     o     o
      else
         call getOtherNearRegionIndices( proj, posc, ffp, pmax, ndn, idnear )
!             =========================   i     i     i    i     o     o
      end if

   end subroutine getOtherNearDomainIndices
!
!=======================================================================
!
   subroutine getDomainDensities( id, nelemr, ielem, prob_, dens_ )
!             ==================  i     o       o      o      o
!
!     Get a list of densities and related properties of a domain.
!
!  Arguments:
!
      integer(INT4), intent(in) :: id                 ! domain index
      integer(INT4), intent(out) :: nelemr            ! number of elements in domain
      integer(INT4), intent(out) :: ielem(:)          ! index of elements in domain
      real(REAL8), optional, intent(out) :: prob_(:)  ! probability of encountering each
                                                      !    element in domain
      real(REAL8), optional, intent(out) :: dens_(:)  ! atomic density of each element in
                                                      !    domain
!
!  Local Variables:
!
      integer(INT4) :: i
      real(REAL8) :: prob(nelem), dens_tot
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%lcell ) then
         call getProbCell( id, prob )
!             ===========  i    o
         nelemr = count( prob > ZERO )
         ielem(1:nelemr) = pack( [ (i, i=1,nelem) ], prob>ZERO )
      else
         nelemr = mat(id)%nelemr
         ielem(1:nelemr) = mat(id)%ielem(1:nelemr)
      end if

      if( present( prob_ ) ) then
         if( geom%lcell ) then
            prob_(1:nelemr) = pack( prob, prob>ZERO )
         else
            prob_(1:nelemr) = mat(id)%prob(1:nelemr)
         end if
      end if

      if( present( dens_ ) ) then
         if( geom%lcell ) then
            call getDensCell( id, dens_tot )
!                ===========  i       o
            dens_(1:nelemr) = dens_tot * pack( prob, prob>ZERO )
         else
            dens_(1:nelemr) = mat(id)%dens * mat(id)%prob(1:nelemr)
         end if
      end if

   end subroutine getDomainDensities
!
!=======================================================================
!
   subroutine getDomainProperty( id, dens_, dens_tot_, vib_, lxtal_, ed_, esurf_ )
!             =================  i     o        o       o      o      o     o
!
!     Get basic properties of the target.
!
!  Arguments:
!
      integer(INT4), intent(in) :: id                    ! domain index
      real(REAL8), optional, intent(out) :: dens_(:)     ! (nelem) atomic density of each
                                                         !    element in domain
      real(REAL8), optional, intent(out) :: dens_tot_    ! total atomic density, =sum(dens_(:))
      real(REAL8), optional, intent(out) :: vib_(:)      ! (nelem) thermal vibration amplitude
      logical, optional, intent(out) :: lxtal_           ! flag indicating whether target is
                                                         !    crystalline
      real(REAL8), optional, intent(out) :: ed_(:)       ! (nelem) displacement energy (eV)
      real(REAL8), optional, intent(out) :: esurf_(:)    ! (nelem) surface binding energy (eV)
!
!  Local Variables:
!
      integer(INT4) :: nelemr, ielem(nelem)
      real(REAL8) :: dens_tot_l, prob(nelem)
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( dens_ ) ) then
         if( geom%lcell ) then
            call getProbCell( id, prob )
!                ===========  i    o
            call getDensCell( id, dens_tot_l )
!                ===========  i      o
         else
            nelemr = mat(id)%nelemr
            ielem(1:nelemr) = mat(id)%ielem(1:nelemr)
            prob(:) = ZERO
            prob(ielem(1:nelemr)) = mat(id)%prob(1:nelemr)
            dens_tot_l = mat(id)%dens
         end if
         dens_(:) = dens_tot_l * prob(:)
      end if

      if( present( dens_tot_ ) ) then
         if( geom%lcell ) then
            call getDensCell( id, dens_tot_ )
!                ===========  i       o
         else
            dens_tot_ = mat(id)%dens
         end if
      end if

      if( present( vib_ ) ) then
         if( geom%lcell ) then
            call getVibCell( id, vib_(:) )
!                ==========  i     o
         else
            vib_(:) = mat(id)%vib(:)
         end if
      end if

      if( present( lxtal_ ) ) then
         if( geom%lcell ) then
            call getXtalCell( id, lxtal_ )
!                ===========  i     o
         else
            lxtal_ = mat(id)%lxtal
         end if
      end if

      if( present( ed_ ) ) then
         if( geom%lcell ) then
            call getEdCell( id, ed_(:) )
!                =========  i   o
         else
            ed_(:) = mat(id)%ed(:)
         end if
      end if

      if( present( esurf_ ) ) then
         if( geom%lcell ) then
            call getEsurfCell( id, esurf_(:) )
!                ============  i     o
         else
            esurf_(:) = mat(id)%esurf(:)
         end if
      end if

   end subroutine getDomainProperty
!
!=======================================================================
!
end module Target_mod

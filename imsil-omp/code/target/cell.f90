module Target_Cell_mod
!
!  Holds cell structure for dynamic targets and provides subroutines operating on it.
!
!  Used Modules:
!  ============
!
   use Const_mod,                only: INT4, REAL8,  &
                                       ZERO, BIG,  &
                                       ENUM

   use Types_mod,                only: Proj_type

   use Param_Setup_mod,          only: ndim, ldyn, cellfile, cellrec
   use Param_Geom_mod,           only: geom
   use Param_Out_mod,            only: out

   use Target_Cell1d_NU_mod,     only: initCell1dNU, readAtomsFromCell1dNU, getCellIndex1dNU,  &
                                       getOtherNearCellIndices1dNU, scoreCell1dNU,  &
                                       broadcastCell1dNU, mergeCell1dNU,  &
                                       relaxCell1dNU, writeCell1dNU, writeDamCell1dNU,  &
                                       getProbCell1dNU, getDensCell1dNU, getVibCell1dNU,  &
                                       getXtalCell1dNU, getEdCell1dNU, getEsurfCell1dNU,  &
                                       readBndFromCell1dNU
   use Target_Cell_Reg_mod,      only: initCellReg, readAtomsFromCellReg,  &
                                       scoreCellReg, getCellIndexReg,  &
                                       getOtherNearCellIndicesReg,  &
                                       getProbCellReg, getDensCellReg, getVibCellReg,  &
                                       getXtalCellReg, getEdCellReg, getEsurfCellReg,  &
                                       broadcastCellReg, mergeCellReg,  &
                                       relaxCellReg, writeCellReg, writeDamCellReg
   use Target_Cell2d_Old_mod,    only: initCell2dOld, readAtomsFromCell2dOld,  &
                                       getCellIndex2dOld, getOtherNearCellIndices2dOld,  &
                                       scoreCell2dOld, relaxCell2dOld,  &
                                       readCell2dOld, writeCell2dOld,  &
                                       getProbCell2dOld, getDensCell2dOld, getVibCell2dOld,  &
                                       getXtalCell2dOld, getEdCell2dOld, getEsurfCell2dOld
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
   public :: initCell                  ! initialize the cell target
!            ========
   public :: readAtomsFromCell         ! read atomic numbers from cell file
!            =================
   public :: readBndFromCell           ! read the geometry boundaries from cell file
!            ===============
   public :: readCell                  ! read the cell structure from file (presently unused!)
!            ========
   public :: scoreCell                 ! score a count in the target cell
!            =========
   public :: getCellIndex              ! determine the cell nearest to point pos
!            ============
   public :: getOtherNearCellIndices   ! determine cells other than nearest cell that
!            =======================   !    intersect with the free flight path cylinder
   public :: getProbCell               ! get the probabilities to encounter each of the atoms
!            ===========               !    in cell ic
   public :: getDensCell               ! get the total atomic density in cell ic
!            ===========
   public :: getVibCell                ! get the thermal vibration amplitude in cell ic
!            ==========
   public :: getXtalCell               ! get the crystallinity flag of the cell ic
!            ===========
   public :: getEdCell                 ! get the displacement energy in cell ic
!            =========
   public :: getEsurfCell              ! get the surface binding energy in cell ic
!            ============
   public :: broadcastCell             ! copy cells from master thread to other threads
!            =============
   public :: mergeCell                 ! merge the contents of cell_tmp into cell
!            =========
   public :: relaxCell                 ! redistribute the cell densities so that stress relaxes
!            =========
   public :: writeCell                 ! write the cell structure to file
!            =========
   public :: writeDamCell              ! write the damage to file
!            ============
!
contains
!
!=======================================================================
!
   subroutine initCell( )
!             ========
!
!     Initialize the target
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
      if( geom%grid == ENUM%grid%regular ) then
         call initCellReg( )
!             ===========
      else
         if( ndim == 1 ) then
            call initCell1dNU( )
!                ============
         else
            call initCell2dOld( )
!                =============
         end if
      end if

!$omp parallel
      call broadcastCell( )
!          =============
!$omp barrier
!$omp end parallel

      if( ldyn .and. out%lcell )  call writeCell( 0 )
!                                      =========  i

   end subroutine initCell
!
!=======================================================================
!
   subroutine readAtomsFromCell( )
!             =================
!
!     Read the atomic numbers from record cellrec of the CELL file cellfile.
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
      if( geom%grid == ENUM%grid%regular ) then
         call readAtomsFromCellReg( )
!             ====================
      else
         if( ndim == 1 ) then
            call readAtomsFromCell1dNU( )
!                =====================
         else
            call readAtomsFromCell2dOld( )   ! not yet implemented!
!                ======================
         end if
      end if

   end subroutine readAtomsFromCell
!
!=======================================================================
!
   subroutine readBndFromCell( )
!             ===============
!
!     Read the atomic numbers from the CELL file cellfile.
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
      if( ndim == 1 ) then
         call readBndFromCell1dNU( geom%zlim(1), geom%zlim(2) )
!             ===================       o             o
      else
         error stop 'readBndFromCell: readBndFromCell2dOld not yet implemented'
!            call readSurfaceFromCell2dOld( cellfile, cellrec )   ! not yet implemented!
!!                ========================     i         i
      end if

   end subroutine readBndFromCell
!
!=======================================================================
!
   subroutine readCell( )
!             ========
!
!     Read the target cell structure from record cellrec of cellfile
!
!     NOTE: This is presently not called from anywhere !
!
!  Arguments:
!
!     none
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         continue                               ! is now called from initCell1dReg
      else
         if( ndim == 1 ) then
            continue                            ! is now called from initCell1dReg
         else
            call readCell2dOld( )
!                =============
         end if
      end if

   end subroutine readCell
!
!=======================================================================
!
   subroutine scoreCell( ielem, ic, w, ddami, ddamv, edam )
!             =========    i    i   i    i      i     i
!
!     Score a count in the cell histogram
!
!     Cell contents are forced not to get less than zero.
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem        ! index of atom species
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(in) :: w              ! weight (can also be negative)
      real(REAL8), intent(in) :: ddami          ! interstitial increment
      real(REAL8), intent(in) :: ddamv          ! vacancy increment
      real(REAL8), intent(in) :: edam           ! damage energy
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call scoreCellReg( ielem, ic, w, ddami, ddamv, edam )
!             ============    i    i   i    i      i     i
      else
         if( ndim == 1 ) then
            call scoreCell1dNU( ielem, ic, w, ddami, ddamv, edam )
!                =============    i    i   i    i      i     i
         else
            call scoreCell2dOld( ielem, ic, w, edam )
!                ==============    i    i   i   i
         end if
      end if

   end subroutine scoreCell
!
!=======================================================================
!
   subroutine getCellIndex( pos, ic_, outside_ )
!             ============   i    o      o
!
!     Determine the cell nearest to point pos
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)            ! position to be checked
      integer(INT4), optional, intent(out) :: ic_  ! cell index
      logical, optional, intent(out) :: outside_   ! flag indicating whether pos(:) is outside
                                                   !    all cells (thus ic_ being only the
                                                   !    cell nearest to ic_)
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call getCellIndexReg( pos, ic_, outside_ )
!             ===============   i   io      o
      else
         if( ndim == 1 ) then
            call getCellIndex1dNU( pos, ic_, outside_ )
!                ================   i   io      o
         else
            call getCellIndex2dOld( pos, ic_, outside_ )
!                =================   i   io      o
         end if
      end if

   end subroutine getCellIndex
!
!=======================================================================
!
   subroutine getOtherNearCellIndices( proj, posc, ffp, pmax, ncn, icnear )
!             =======================   i     i     i    i     o     o
!
!     Determine cells other than proj%id that are within the free flight path cylinder
!     (= cylinder with height ffp and radius pmax)
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile
      real(REAL8), intent(in) :: posc(3)        ! position of end of free flight path
      real(REAL8), intent(in) :: ffp            ! free flight path (= height of cylinder)
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter (= radius of cyl.)
      integer(INT4), intent(out) :: ncn         ! number of cells found
      integer(INT4), allocatable, intent(inout) :: icnear(:)
                                                ! index of cells found
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call getOtherNearCellIndicesReg( proj, posc, ffp, pmax, ncn, icnear )
!             ==========================   i     i     i    i    io     o
      else
         if( ndim == 1 ) then
            call getOtherNearCellIndices1dNU( proj, posc, ffp, pmax, ncn, icnear )
!                ===========================   i     i     i    i    io     o
         else
            call getOtherNearCellIndices2dOld( proj, posc, ffp, pmax, ncn, icnear )
!                ============================   i     i     i    i    io     o
         end if
      end if

   end subroutine getOtherNearCellIndices
!
!=======================================================================
!
   subroutine getProbCell( ic, prob )
!             ===========  i    o
!
!     Get the probabilities to encounter each of the atoms in cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(out) :: prob(:)       ! (nelem) probabilities
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call getProbCellReg( ic, prob )
!             ==============  i    o
      else
         if( ndim == 1 ) then
            call getProbCell1dNU( ic, prob )
!                ===============  i    o
         else
            call getProbCell2dOld( ic, prob )
!                ================  i    o
         end if
      end if

   end subroutine getProbCell
!
!=======================================================================
!
   subroutine getDensCell( ic, dens )
!             ===========  i    o
!
!     Get the total atomic density in cell ic
!
!     NOTE: In case of surface roughness this is the density before reduction due to roughness
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(out) :: dens          ! density
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call getDensCellReg( ic, dens )
!                ================  i    o
      else
         if( ndim == 1 ) then
            call getDensCell1dNU( ic, dens )
!                ===============  i    o
         else
            call getDensCell2dOld( ic, dens )
!                ================  i    o
         end if
      end if

   end subroutine getDensCell
!
!=======================================================================
!
   subroutine getVibCell( ic, vib )
!             ==========  i    o
!
!     Get the thermal vibration amplitude in cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(out) :: vib(:)        ! thermal vibration amplitudes
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call getVibCellReg( ic, vib(:) )
!             =============  i     o
      else
      if( ndim == 1 ) then
            call getVibCell1dNU( ic, vib(:) )
!                ==============  i     o
         else
            call getVibCell2dOld( ic, vib(:) )
!                ===============  i     o
         end if
      end if

   end subroutine getVibCell
!
!=======================================================================
!
   subroutine getXtalCell( ic, lxtal )
!             ===========  i     o
!
!     Get the crystallinity flag of the cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      logical, intent(out) :: lxtal             ! crystallinity flag
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call getXtalCellReg( ic, lxtal )
!             ==============  i     o
      else
         if( ndim == 1 ) then
            call getXtalCell1dNU( ic, lxtal )
!                =============  i     o
         else
            call getXtalCell2dOld( ic, lxtal )
!                ================  i     o
         end if
      end if

   end subroutine getXtalCell
!
!=======================================================================
!
   subroutine getEdCell( ic, ed )
!             =========  i   o
!
!     Get the displacement energy in cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(out) :: ed(:)         ! displaceent energy
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call getEdCellReg( ic, ed(:) )
!             ============  i   o
      else
         if( ndim == 1 ) then
            call getEdCell1dNU( ic, ed(:) )
!                =============  i   o
         else
            call getEdCell2dOld( ic, ed(:) )
!                ==============  i   o
         end if
      end if

   end subroutine getEdCell
!
!=======================================================================
!
   subroutine getEsurfCell( ic, esurf )
!             ============  i     o
!
!     Get the surface binding energy in cell ic
!
!  Arguments:
!
      integer(INT4), intent(in) :: ic           ! cell index
      real(REAL8), intent(out) :: esurf(:)      ! surface binding energy
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call getEsurfCellReg( ic, esurf(:) )
!             ===============  i     o
      else
      if( ndim == 1 ) then
            call getEsurfCell1dNU( ic, esurf(:) )
!                ================  i     o
         else
            call getEsurfCell2dOld( ic, esurf(:) )
!                =================  i     o
         end if
      end if

   end subroutine getEsurfCell
!
!=======================================================================
!
   subroutine broadcastCell( )
!             =============
!
!     Copy cells from master thread to other threads
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call broadcastCellReg( )
!             ================
      else
         if( ndim == 1 ) then
            call broadcastCell1dNU( )
!                =================
         else
            error stop 'broadcastCell: broadcastCell2d not yet implemented'
!            call broadcastCell2d( )
!!                ===============
         end if
      end if

   end subroutine broadcastCell
!
!=======================================================================
!
   subroutine mergeCell( )
!             =========
!
!     Merge the contents of cell_tmp into cell
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call mergeCellReg( )
!             ============
      else
         if( ndim == 1 ) then
            call mergeCell1dNU( )
!                =============
         else
            error stop 'mergeCell: mergeCell2d not yet implemented'
!            call mergeCell2d( )
!!                ===========
         end if
      end if

   end subroutine mergeCell
!
!=======================================================================
!
   subroutine relaxCell( iion )
!             =========   i
!
!     Relax the cells
!
!  Arguments:
!
      integer(INT4), intent(in) :: iion         ! index of ion
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code
!  ------------------------
!
!$omp master
      if( geom%grid == ENUM%grid%regular ) then
         call relaxCellReg( )
!             ============
      else
         if( ndim == 1 ) then
            call relaxCell1dNU( )
!                =============
         else
            call relaxCell2dOld( iion )
!                ==============    i
         end if
      end if
!$omp end master
!$omp barrier

      call broadcastCell( )
!          =============

   end subroutine relaxCell
!
!=======================================================================
!
   subroutine writeCell( iion )
!             =========   i
!
!     Write the cell structure of the target to file
!
!  Arguments:
!
      integer(INT4), intent(in) :: iion         ! index of ion
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( geom%grid == ENUM%grid%regular ) then
         call writeCellReg( )
!             ============
      else
         if( ndim == 1 ) then
            call writeCell1dNU( )
!                =============
         else
            call writeCell2dOld( iion )
!                ==============    i
         end if
      end if
!
   end subroutine writeCell
!
!=======================================================================
!
   subroutine writeDamCell( )
!             ============
!
!     Write damage to file
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
      if( geom%grid == ENUM%grid%regular ) then
         call writeDamCellReg( )
!             ===============
      else
         if( ndim == 1 ) then
            call writeDamCell1dNU( )
!                ================
         else
            return   ! Cell2d case not yet implemented'
!            call writeCell2dOld( )
!                ==============
         end if
      end if

   end subroutine writeDamCell
!
!=======================================================================
!
end module Target_Cell_mod

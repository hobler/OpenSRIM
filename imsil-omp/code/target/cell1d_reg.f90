module Target_Cell1d_Reg_mod
!
!  Holds 1d cell structureon a regular grid and provides subroutines operating on it.
!  A regular grid is maintained even upon relaxation.
!
!  Note: avoid ldyn=t with lspher=t
!
!  Used Modules:
!  ============
!
!$ use omp_lib,                     only: omp_get_thread_num
   use Const_mod,                   only: INT4, REAL8,  &
                                          ZERO, ONE, HALF, BIG, EPSMA,  &
                                          ENUM
   use Array_mod,                   only: alloc, changeShape, addReflect
   use Open_mod,                    only: openFile
   use Stack_mod,                   only: real8Stack, newStack, pushStack, popStack, sizeOfStack
   use String_mod,                  only: str
   use Types_mod,                   only: Proj_type
   use Undefine_mod,                only: defined

   use Param_Setup_mod,             only: nr, nelem, nelemIon, nelemTarg,  &
                                          ldyn, usecell, cellfile, cellrec, usehis, hisfile,  &
                                          reflect, shiftto
   use Param_Geom_mod,              only: geom
   use Param_Elem_mod,              only: elem
   use Param_Mat_mod,               only: mat, irxtal
   use Param_Dam_mod,               only: dam
   use Param_Out_mod,               only: out

   use Geom_CurviCoord_mod,         only: toCurviCoord1d
   use Math_Solve_mod,              only: tridiag

#ifdef AMAT
   use Cinterface_Mechanics1d_mod,  only: relaxViscoElasticEuler1d, initSurface1dFromFF
#endif
   use Stat_His1d_mod,              only: His1d_type, writeHis1d
   use Target_Cell1d_NU_mod,        only: readCell1dNUraw
   use Target_Surface1d_mod,        only: updateSurface1d
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
   public :: readAtomsFromCell1dReg          ! read atomic numbers from cell file
!            ======================
   public :: initCell1dReg                   ! initialize grid and cells
!            =============
   public :: initHis1dReg                    ! initialize damage
!            ============
   public :: getNearestCell1dReg             ! get nearest non-empty cell index
!            ===================
   public :: getOtherNearCellIndices1dReg    ! determine cells other than nearest cell that
!            ============================    !    intersect with the free flight path cylinder
   public :: fixNegDensCell1dReg             ! fix negative densities by drawing material from
!            ===================             !    neighboring cells
   public :: diffuseCell1dReg                ! diffuse the atom species
!            ================
   public :: relaxCell1dReg                  ! relax the cell densities
!            ==============
   public :: relaxCell1dRegFull              ! fully relax the cell densities
!            ==================
   public :: writeCell1dReg                  ! write target cell contents to CELL file
!            ==============
   public :: writeDamCell1dReg               ! write damage to output file
!            =================
!
!  Private Procedures:
!  ==================
!
   private :: initCell1dRegFromGeom          ! initialize grid and cell from region geometry
!             =====================
   private :: readCell1dReg                  ! read target cell contents from file
!             =============
   private :: readHis1dReg                   ! read damage from file
!             ============

!  Public Data:
!  ===========
!
   real(REAL8), parameter :: HISMIN = 0.0001
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine readAtomsFromCell1dReg( )
!             ======================
!
!     Read the atomic numbers from record cellrec of the CELL file cellfile.
!     The CELL file has the same format as the HIS file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: idummy
      integer(INT4) :: i, ielem
      real(REAL8), allocatable :: z(:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Open file
!
      call openFile( 22, cellfile, 'formatted', 'read' )
!          ========  i      i           i         i
!
!  First line can be ignored
!
      read( 22, *, err=901, end=901 )  ( idummy, i=1,2 )
!
!  Read the atomic numbers and determine the atom masses
!
      if( ldyn ) then
         call alloc( z, nelem, 'z', 'readAtomsFromCell1dReg' )
!             =====  io   i     i               i
         read( 22, *, err=900, end=900 )  idummy, z(:)
         do ielem = 1, nelemIon
            if( z(ielem) /= elem(ielem)%z )  &
               error stop 'readAtomsFromCell1dReg: Inconsistent ion atomic number'
         end do
         do ielem = nelemIon+1, nelem
            elem(ielem)%z = z(ielem)
         end do
         deallocate( z )
      else
         read( 22, *, err=900, end=900 )  idummy, elem(nelemIon+1:nelem)%z
      end if

      close( 22 )

      return

900   error stop 'readAtomsFromCell1dReg: error during reading CELL file'
901   error stop 'readAtomsFromCell1dReg: no data on CELL file'

   end subroutine readAtomsFromCell1dReg
!
!=======================================================================
!
   subroutine initCell1dReg( z, dens, ffgeom, lxtal, pressure, stress )
!             =============  io  io     io     io       io       io
!
!     Read the contents of the target cells from the cellfile file (same format as HIS file)
!
!     At the target boundaries, the densities are averaged over the cell.
!     Geometric fill fractions (ffgeom) are determined
!
!     TODO: call readCell1dNUraw
!     TODO: eliminate shiftto and reflect or consider it in initSurface1d
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: z(:)          ! positions of grid points
      real(REAL8), allocatable, intent(inout) :: dens(:,:)     ! (nelem,nc) densities of atoms
      real(REAL8), allocatable, intent(inout) :: ffgeom(:)     ! geometric fill fraction
      logical, allocatable, intent(inout) :: lxtal(:)          ! crystallinity flag
      real(REAL8), allocatable, intent(inout) :: pressure(:)   ! pressure
      real(REAL8), allocatable, intent(inout) :: stress(:,:,:) ! zz component of stress tensor
!
!  Local Variables:
!
      integer(INT4) :: direction
!
!  Start of Executable Code:
!  ------------------------
!
!  Initialize the arrays
!
      if( usecell ) then
         call readCell1dReg( z, dens, ffgeom, lxtal, pressure, stress )
!             =============  io  io     io     io       io       io
      else
         call initCell1dRegFromGeom( z, dens, ffgeom, lxtal, pressure, stress )
!             =====================  io  io     io     io       io       io
      end if
!
!  Apply reflective boundary conditions (if any)
!
      if( geom%bndtypez(1) == ENUM%bndtype%reflect ) then
         direction = -1
         geom%bndtypez(1) = geom%bndtypez(2)
      else if( geom%bndtypez(2) == ENUM%bndtype%reflect ) then
         direction = +1
         geom%bndtypez(2) = geom%bndtypez(1)
      else
         direction = 0
      end if

      call addReflect( z, direction, .false., 'z', 'initCell1dReg' )
!          ========== io      i         i      i          i
      call addReflect( dens, [0, direction], [.false., .true.], 'dens', 'initCell1dReg' )
!          ==========   io         i                i             i            i
      call addReflect( ffgeom, direction, .true., 'ffgeom', 'initCell1dReg' )
!          ==========    io        i        i        i             i
      call addReflect( lxtal, direction, .true., 'lxtal', 'initCell1dReg' )
!          ==========   io        i        i        i            i
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         call addReflect( pressure, direction, .true., 'pressure', 'initCell1dReg' )
!             ==========     io         i        i         i              i
         call addReflect( stress, [0, 0, direction], [.false., .false., .true.], 'stress', &
!             ==========    io          i                         i                i
                         'initCell1dReg' )
!                               i
      end if

   end subroutine initCell1dReg
!
!=======================================================================
!
   subroutine readCell1dReg( z, dens, ffgeom, lxtal, pressure, stress )
!             =============  io  io     io     io       io       io
!
!     Read the contents of the target cells from the cellfile file (same format as HIS file)
!
!     At the target boundaries, the densities are averaged over the cell.
!     Geometric fill fractions (ffgeom) are determined, and the surface is initialized
!
!     TODO: call readCell1dNUraw
!     TODO: eliminate shiftto and reflect or consider it in initSurface1d
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: z(:)          ! positions of grid points
      real(REAL8), allocatable, intent(inout) :: dens(:,:)     ! (nelem,nc) densities of atoms
      real(REAL8), allocatable, intent(inout) :: ffgeom(:)     ! geometric fill fraction
      logical, allocatable, intent(inout) :: lxtal(:)          ! crystallinity flag
      real(REAL8), allocatable, intent(inout) :: pressure(:)   ! pressure
      real(REAL8), allocatable, intent(inout) :: stress(:,:,:) ! zz component of stress tensor
!
!  Local Variables:
!
      logical :: regular
      integer(INT4) :: ncol, nrow, irec, nrec, nz, nz_tmp
      integer(INT4) :: iz, jz, izmin, izmax, i, ielem1, idummy
      real(REAL8) :: rdummy, wz, dwz, offset, shift, zminr, zmaxr, zmin, zmax
      real(REAL8), allocatable :: ff(:)                        ! fill fraction
      real(REAL8), allocatable :: z_tmp(:), dens_tmp(:,:), pressure_tmp(:), stress_tmp(:,:,:)
!
!  Start of Executable Code:
!  ------------------------
!
!  Open file
!
      call openFile( 22, cellfile, 'formatted', 'read' )
!          ========  i      i           i         i

      nrec = cellrec
      if( nrec == 0 )  nrec = huge(0_INT4)
      nz = 0

      do irec = 1, nrec
!
!  Read first line and allocate memory
!
         read( 22, *, err=800, end=800 )  ncol, nrow

         if( ncol-1 /= nelemTarg )  &
            error stop 'readCell1dReg: inconsistent number of columns on file'

         nz = ( nrow - 2 ) / 2

         if( allocated( z ) )  deallocate( z )
         call alloc( z, 0, nz, 'z', 'readCell1dReg' )
!             =====  io i  i    i          i
         if( allocated( dens ) )  deallocate( dens )
         call alloc( dens, [nelem, nz], 'dens', 'readCell1dReg' )
!             =====   io       i          i            i
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            if( allocated( pressure ) )  deallocate( pressure )
            call alloc( pressure, nz, 'pressure', 'readCell1dReg' )
!                =====     io     i       i              i
            if( allocated( stress ) )  deallocate( stress )
            call alloc( stress, [1,1,nz], 'stress', 'readCell1dReg' )
!                =====    io       i         i             i
         end if
!
!  Second and third line can be ignored
!
         read( 22, *, err=900, end=900 )  ( idummy, i=1,ncol )
         read( 22, *, err=900, end=900 )  ( rdummy, i=1,ncol )
!
!  Read cell densities
!
         if( ldyn ) then
            ielem1 = 1
         else                                   ! Is ldyn=f possible?
            dens(1:nelemIon,:) = ZERO
            ielem1 = nelemIon+1
         end if

         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            do iz = 1, nz
               read( 22, *, err=900, end=900 )  z(iz-1), dens(ielem1:nelem,iz),  &
                                                pressure(iz), stress(1,1,iz)
               read( 22, *, err=900, end=900 )  z(iz), dens(ielem1:nelem,iz),  &
                                                pressure(iz), stress(1,1,iz)
            end do
         else
            do iz = 1, nz
               read( 22, *, err=900, end=900 )  z(iz-1), dens(ielem1:nelem,iz)
               read( 22, *, err=900, end=900 )  z(iz), dens(ielem1:nelem,iz)
            end do
         end if

         dens(ielem1:nelem,:) = dens(ielem1:nelem,:) / 1e24
!
!  Ignore last line
!
         read( 22, *, err=900, end=900 )  ( rdummy, i=1,ncol )

      end do

800   if( irec == 1 )  goto 901

      close( 22 )
!
!  Check cell widths
!
      wz = z(1) - z(0)
      dwz = 0.001 * wz
      regular = .true.
      do iz = 2, nz
         if( abs( z(iz)-z(iz-1) - wz ) > dwz ) then
            regular = .false.
            exit
         end if
      end do
!
!  If the grid on the cell file is regular, take it as is
!  If it is not regular, construct a new one as if it were created from regions
!
      if( .not.regular ) then
!
!  grid
!
         if( irxtal > 0 ) then
            offset = -100 * EPSMA * z(nz)  ! avoid cell boundaries at lattice planes
         else
            offset = ZERO
         end if

         if( out%nbox == 0 ) then
            wz = out%wbox

            izmin = floor( ( z(0) - offset ) / wz )
            izmax = ceiling( ( z(nz) - offset ) / wz )
            nz_tmp = izmax - izmin

            call alloc( z_tmp, 0, nz_tmp, 'z_tmp', 'readCell1dReg' )
!                =====    o    i    i        i             i
            z_tmp(:) = [ (offset + iz*wz, iz=izmin,izmax) ]
         else if( out%nbox == 1 ) then
            error stop 'readCell1dReg: requires out%nbox /= 1'
         else
            nz_tmp = out%nbox

            wz = (z(nz) - z(0)) / (nz_tmp - 1)

            call alloc( z_tmp, 0, nz_tmp, 'z_tmp', 'readCell1dReg' )
!                =====   io    i    i        i             i
            z_tmp(:) = [ (z(0) + iz*wz, iz=0,nz_tmp) ]

            shift = modulo( z_tmp(0)-offset, wz )
            z_tmp(:) = z_tmp(:) - shift
         end if
!
!  cells
!
         call alloc( dens_tmp, [nelem, nz_tmp], 'dens_tmp', 'readCell1dReg' )
!             =====     io            i             i              i
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            call alloc( pressure_tmp, nz_tmp, 'pressure_tmp', 'readCell1dReg' )
!                =====       io         i           i                i
            call alloc( stress_tmp, [1,1,nz_tmp], 'stress_tmp', 'readCell1dReg' )
!                =====      io           i             i               i
         end if

         dens_tmp(:,:) = ZERO
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            pressure_tmp(:) = ZERO
            stress_tmp(1,1,:) = ZERO
         end if

         do jz = 1, nz
            zminr = ( z(jz-1) - z_tmp(0) ) / wz
            zmaxr = ( z(jz) - z_tmp(0) ) / wz

            izmin = max( 1, ceiling( zminr ) )
            izmax = max( 1, ceiling( zmaxr ) )

            do iz = izmin, izmax
               dens_tmp(:,iz) = dens_tmp(:,iz) + dens(:,jz)
            end do

            if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
               pressure_tmp(izmin:izmax) = pressure_tmp(izmin:izmax) + pressure(jz)
               stress_tmp(1,1,izmin:izmax) = stress_tmp(1,1,izmin:izmax) + stress(1,1,jz)
            end if

            dens_tmp(:,izmin) = dens_tmp(:,izmin) - dens(:,jz) * ( zminr - (izmin-1) )
            dens_tmp(:,izmax) = dens_tmp(:,izmax) - dens(:,jz) * ( izmax - zmaxr )

            if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
               pressure_tmp(izmin) = pressure_tmp(izmin) - pressure(jz) * ( zminr - (izmin-1) )
               stress_tmp(1,1,izmin) = stress_tmp(1,1,izmin)  &
                                       - stress(1,1,jz) * ( zminr - (izmin-1) )
               pressure_tmp(izmax) = pressure_tmp(izmax) - pressure(jz) * ( izmax - zmaxr )
               stress_tmp(1,1,izmax) = stress_tmp(1,1,izmax)  &
                                       - stress(1,1,jz) * ( izmax - zmaxr )
            end if

         end do

         nz = nz_tmp
         call move_alloc( z_tmp, z )
         call move_alloc( dens_tmp, dens )
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            call move_alloc( pressure_tmp, pressure )
            call move_alloc( stress_tmp, stress )
         end if

      end if
!
!  Shift and reflect the distribution if desired
!
      if( defined( shiftto ) )  &
         z(:) = shiftto + z(:) - z(0)

      if( reflect ) then
         z(:) = - z(nz:0:-1)
         dens(:,:) = dens(:,nz:1:-1)
         if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
            pressure(:) = pressure(nz:1:-1)
            stress(1,1,:) = stress(1,1,nz:1:-1)
         end if
      end if
!
!  Determine fill fractions (neglecting surface effect)
!
      call alloc( ff, nz, 'ff', 'readCell1dReg' )
!          =====  io  i    i           i
      ff(:) = ONE

      call alloc( ffgeom, nz, 'ffgeom', 'readCell1dReg' )
!          =====    io    i      i              i
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         ffgeom(:) = ff(:) * (ONE - pressure(:)/dam%modbulk)
      else
         ffgeom(:) = ff(:)
      end if
!
!  Correct for surface effects and determine target boundaries
!
      do iz = 1, nz
         if( sum( dens(:,iz) ) > ZERO )  exit
         ff(iz) = ZERO
         ffgeom(iz) = ZERO
      end do
      if( iz > nz )  goto 902
      ff(iz) = sum( dens(:,iz) / elem(:)%dens )
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         ffgeom(:) = ff(:) * (ONE - pressure(:)/dam%modbulk)
      else
         ffgeom(:) = ff(:)
      end if
      zmin = z(iz) - ffgeom(iz) * (z(iz)-z(iz-1))

      do iz = nz, 1, -1
         if( sum( dens(:,iz) ) > ZERO )  exit
         ff(iz) = ZERO
         ffgeom(iz) = ZERO
      end do
      if( iz < 1 )  goto 902
      ff(iz) = sum( dens(:,iz) / elem(:)%dens )
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         ffgeom(:) = ff(:) * (ONE - pressure(:)/dam%modbulk)
      else
         ffgeom(:) = ff(:)
      end if
      zmax = z(iz-1) + ffgeom(iz) * (z(iz)-z(iz-1))
!
!  Initialize the surface
!
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
#ifdef AMAT
         call initSurface1dFromFF( z(:), ffgeom(:), dens(:,:), 0_INT4 )
!             ===================   i        i          i        i
#else
         error stop 'ViscoElasticEuler1d model requires separate license'
#endif
      else
         call updateSurface1d( zmin, zmax )
!             ===============   i     i
      end if
!
!  CELL file does not contain info on crystallinity, therefore set to amorphous
!
      if( allocated( lxtal ) )  deallocate( lxtal )
      call alloc( lxtal, nz, 'lxtal', 'readCell1dReg' )
!          =====   io    i      i            i
      lxtal(:) = .false.

      return

900   print *, 'iz=', iz
      error stop 'readCell1dReg: error during reading CELL file'
901   error stop 'readCell1dReg: no data on CELL file'
902   error stop 'readCell1dReg: Target seems to be empty'


   end subroutine readCell1dReg
!
!=======================================================================
!
   subroutine initCell1dRegFromGeom( z,  &
!             =====================  io
                                     dens, ffgeom, lxtal, pressure, stress )
!                                     io     io     io       io       io
!
!     Initialize the cell structure of the target from a 1d region definition
!
!     The target extends from geom%posif(0) to geom%posif(nr). On entry, geom%zlim(:) has
!     already been initialized to include a vacuum surface layer at free boundaries and the
!     reflected region at reflected boundaries.
!
!     This procedure does NOT perform the reflection at a reflective boundary.
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: z(:)          ! positions of grid points
      real(REAL8), allocatable, intent(inout) :: dens(:,:)     ! (nelem,nc) densities of atoms
      real(REAL8), allocatable, intent(inout) :: ffgeom(:)     ! geometric fill fraction
      logical, allocatable, intent(inout) :: lxtal(:)          ! crystallinity flag
      real(REAL8), allocatable, intent(inout) :: pressure(:)   ! pressure
      real(REAL8), allocatable, intent(inout) :: stress(:,:,:) ! zz component of stress tensor
!
!  Local Variables:
!
      integer(INT4) :: nz, iz, izmin, izmax
      integer(INT4) :: ir, ielem, ielemr

      real(REAL8) :: zmin, zmax, zminr, zmaxr
      real(REAL8) :: ddens
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine (unreflected) grid from geom%zlim(:) and target geometry
!
      if( geom%bndtypez(1) == ENUM%bndtype%reflect ) then
         zmin = geom%zlim(1)
         zmax = geom%zlim(3)
      else if( geom%bndtypez(2) == ENUM%bndtype%reflect ) then
         zmin = geom%zlim(3)
         zmax = geom%zlim(2)
      else
         zmin = geom%zlim(1)
         zmax = geom%zlim(2)
      end if

      izmin = - ceiling( (geom%posif(0)-zmin) / geom%wbox )
      izmax = ceiling( (zmax-geom%posif(0)) / geom%wbox )
      nz = izmax - izmin

      if( allocated( z ) )  deallocate( z )
      call alloc( z, 0, nz, 'z', 'initCell1dRegFromGeom' )
!          =====  o  i        i         i              i
      z(:) = [ ( geom%posif(0)+iz*geom%wbox, iz=izmin,izmax ) ]
!
!  Determine fill fractions
!
      call alloc( ffgeom, nz, 'ffgeom', 'initCell1dRegFromGeom' )
!          =====    io    i      i                 i
      ffgeom(:) = ( max( zmin, min( zmax, z(1:nz) ) )  &
                    - max( zmin, min( zmax, z(0:nz-1) ) ) ) / geom%wbox
!
!  Cell densities
!
      call alloc( dens, [nelem, nz], 'dens', 'initCell1dRegFromGeom' )
!          =====   io        i         i                i
      dens(:,:) = ZERO

      do ir = 1, nr
         zminr = (geom%posif(ir-1)-z(0)) / geom%wbox
         zmaxr = (geom%posif(ir)-z(0)) / geom%wbox
         izmin = max( 1, ceiling( zminr ) )
         izmax = max( 1, ceiling( zmaxr ) )
         do ielemr = 1, mat(ir)%nelemr
            ielem = mat(ir)%ielem(ielemr)
            ddens = mat(ir)%dens * mat(ir)%prob(ielemr)
            dens(ielem,izmin:izmax) = dens(ielem,izmin:izmax) + ddens
            dens(ielem,izmin) = dens(ielem,izmin) - ddens * ( zminr - (izmin-1) )
            dens(ielem,izmax) = dens(ielem,izmax) - ddens * ( izmax - zmaxr )
         end do
      end do
!
!  Pressure and stress
!
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         call alloc( pressure, nz, 'pressure', ' initCell1dRegFromGeom' )
!             =====     io     i       i                   i
         pressure(:) = ZERO
         call alloc( stress, [1,1,nz], 'stress', ' initCell1dRegFromGeom' )
!             =====    io       i         i                  i
         stress(1,1,:) = ZERO
      end if
!
!  Crystal flag
!
      call alloc( lxtal, nz, 'lxtal', ' initCell1dRegFromGeom' )
!          =====   io    i     i                  i
      lxtal(:) = .true.
      do ir = 1, nr
         if( mat(ir)%lxtal )  cycle
         zminr = (geom%posif(ir-1)-z(0)) / geom%wbox
         zmaxr = (geom%posif(ir)-z(0)) / geom%wbox
         izmin = max( 1, ceiling( zminr ) )
         izmax = max( 1, ceiling( zmaxr ) )
         lxtal(izmin:izmax) = .false.
      end do

   end subroutine initCell1dRegFromGeom
!
!=======================================================================
!
!   subroutine initCell1dRegFromGeom( z,  &
!!             =====================  io
!                                     dens, ffgeom, lxtal, pressure, stress )
!!                                     io     io     io       io       io
!!
!!     Initialize the cell structure of the target from a 1d region definition
!!
!!     The target extends from geom%posif(0) to geom%posif(nr). On entry, geom%zlim(:) has
!!     already been initialized to include a vacuum surface layer at free boundaries and the
!!     reflected region at reflected boundaries. geom%posif(:) still is as defined in INP.
!!
!!  Arguments:
!!
!      real(REAL8), allocatable, intent(inout) :: z(:)          ! positions of grid points
!      real(REAL8), allocatable, intent(inout) :: dens(:,:)     ! (nelem,nc) densities of atoms
!      real(REAL8), allocatable, intent(inout) :: ffgeom(:)     ! geometric fill fraction
!      logical, allocatable, intent(inout) :: lxtal(:)          ! crystallinity flag
!      real(REAL8), allocatable, intent(inout) :: pressure(:)   ! pressure
!      real(REAL8), allocatable, intent(inout) :: stress(:,:,:) ! zz component of stress tensor
!!
!!  Local Variables:
!!
!      integer(INT4) :: nzadd(2), nz, iz, izmin, izmax, ir, irr
!      integer(INT4) :: ielem, ielemr
!
!      integer(INT4), allocatable :: jr(:)       ! region indizes corresponding to posifr
!      real(REAL8), allocatable :: posifr(:)     ! region interfaces including reflection
!      real(REAL8) :: wz, zmin, zmax
!      real(REAL8) :: offset, zminr, zmaxr
!      real(REAL8) :: ddens, shift
!
!  Start of Executable Code:
!  ------------------------
!
!  Reflect geometry if there is a reflecting boundary
!
!     if( geom%bndtypez(1) == ENUM%bndtype%reflect ) then
!        nrr = 2*nr
!        call alloc( jr, nrr, 'jr', 'initCell1dRegFromGeom' )
!             =====   o    i     i               i
!        jr(:) = [ (ir, ir=nr,1,-1), (ir, ir=1,nr) ]
!        call alloc( posifr, 0, nrr, 'posifr', 'initCell1dRegFromGeom' )
!             =====     o    i   i    i               i
!        posifr(:) = [ (geom%posif(ir), ir=nr,1,-1), (geom%posif(ir), ir=0,nr) ]
!     else if( geom%bndtypez(2) == ENUM%bndtype%reflect ) then
!        nrr = 2*nr
!        call alloc( jr, nrr, 'jr', 'initCell1dRegFromGeom' )
!             =====   o    i     i               i
!        jr(:) = [ (ir, ir=1,nr), (ir, ir=nr,1,-1) ]
!        call alloc( posifr, 0, nrr, 'posifr', 'initCell1dRegFromGeom' )
!             =====     o    i   i    i               i
!        posifr(:) = [ (geom%posif(ir), ir=0,nr), (geom%posif(ir), ir=nr-1,0,-1) ]
!     else
!        nrr = nr
!        call alloc( jr, nr, 'jr', 'initCell1dRegFromGeom' )
!             =====   o   i     i               i
!        jr(:) = [ (ir, ir=1,nr) ]
!        call alloc( posifr, 0, 2*nr+1, 'posifr', 'initCell1dRegFromGeom' )
!             =====     o    i    i      i               i
!        posifr(:) = posif(:)
!     end if
!
!  Determine grid from target geometry
!
!     zmin = posifr(0)
!     zmax = posifr(nrr)
!     izmin = - ceiling( (zmin-geom%zlim(1)) / geom%wbox )
!     izmax = ceiling( (geom%zlim(2)-geom%posif(0)) / geom%wbox )
!     nz = izmax - izmin + 1
!     if( allocated( z ) )  deallocate( z )
!     call alloc( z, 0, nz, 'z', 'initCell1dRegFromGeom' )
!          =====  o  i        i         i              i
!     z(:) = [ ( posifr(0)+iz*geom%wbox, iz=izmin,izmax ) ]
!
!  Determine fill fractions
!
!     call alloc( ffgeom, nz, 'ffgeom', 'initCell1dRegFromGeom' )
!          =====    io    i      i                 i
!     ffgeom(:) = ( max( zmin, min( zmax, z(1:nz) ) )  &
!                   - max( zmin, min( zmax, z(0:nz-1) ) ) ) / geom%wbox
!
!  Cell densities
!
!     call alloc( dens, [nelem, nz], 'dens', 'initCell1dRegFromGeom' )
!          =====   io        i         i                i
!     dens(:,:) = ZERO
!     do irr = 1, nrr
!        zminr = (posifr(irr-1)-z(0)) / geom%wbox
!        zmaxr = (posifr(irr)-z(0)) / geom%wbox
!        izmin = max( 1, ceiling( zminr ) )
!        izmax = max( 1, ceiling( zmaxr ) )
!        ir = jr(irr)
!        do ielemr = 1, mat(ir)%nelemr
!           ielem = mat(ir)%ielem(ielemr)
!           ddens = mat(ir)%dens * mat(ir)%prob(ielemr)
!           dens(ielem,izmin:izmax) = dens(ielem,izmin:izmax) + ddens
!           dens(ielem,izmin) = dens(ielem,izmin) - ddens * ( zminr - (izmin-1) )
!           dens(ielem,izmax) = dens(ielem,izmax) - ddens * ( izmax - zmaxr )
!        end do
!     end do
!
!  Pressure and stress
!
!     if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
!        call alloc( pressure, nz, 'pressure', ' initCell1dRegFromGeom' )
!             =====     io     i       i                   i
!        pressure(:) = ZERO
!        call alloc( stress, [1,1,nz], 'stress', ' initCell1dRegFromGeom' )
!             =====    io       i         i                  i
!        stress(1,1,:) = ZERO
!     end if
!
!  Crystal flag
!
!     call alloc( lxtal, nz, 'lxtal', ' initCell1dRegFromGeom' )
!          =====   io    i     i                  i
!     lxtal(:) = .true.
!     do irr = 1, nr
!        ir = jr(irr)
!        if( mat(ir)%lxtal )  cycle
!        zminr = (posifr(irr-1)-z(0)) / geom%wbox
!        zmaxr = (posifr(irr)-z(0)) / geom%wbox
!        izmin = max( 1, ceiling( zminr ) )
!        izmax = max( 1, ceiling( zmaxr ) )
!        lxtal(izmin:izmax) = .false.
!     end do
!
!  end subroutine initCell1dRegFromGeom
!
!=======================================================================
!
   subroutine initHis1dReg( z, damdens )
!             ============  i    io
!
!     Initialize the damage concentrations either by setting them zero or by reading
!     them from the hisfile file
!
!  Arguments:
!
      real(REAL8), intent(in) :: z(0:)                         ! positions of grid points
      real(REAL8), allocatable, intent(inout) :: damdens(:,:)  ! (nelemIn,nc) damage densities
!
!  Local Variables:
!
      integer(INT4) :: nz, nelemIn
!
!  Start of Executable Code:
!  ------------------------
!
      nelemIn = 2*nelem        ! maximum atom index (note: ielem>nelem ... vacancies)
      nz = ubound( z, 1 )

      if( usehis ) then
         call readHis1dReg( z, damdens )
!             ============  i    io
      else
         call alloc( damdens, [nelemIn, nz], 'damdens', 'initCell1dReg' )
!             =====     io           i           i              i
         damdens(:,:) = ZERO
      end if

   end subroutine initHis1dReg
!
!=======================================================================
!
   subroutine getNearestCell1dReg( iz, denstot, jz )
!             ===================  i      i     o
!
!     Determine the cell the point pos is located in or closest to.
!
!  Arguments:
!
      integer(INT4), intent(in) :: iz           ! index of cell to which the nearest non-empty
                                                !    cell is sought
      real(REAL8), intent(in) :: denstot(:)     ! total densities
      integer(INT4), intent(out) :: jz          ! index of non-empty cell nearest to iz
!
!  Local Variables:
!
      integer(INT4) :: nz, i
!
!  Start of Executable Code:
!  ------------------------
!
      nz = size( denstot )

      if( denstot(iz) == ZERO ) then
         do i = 1, nz
            jz = iz + i
            if( jz <= nz ) then
               if( denstot(jz) /= ZERO )  exit
            end if
            jz = iz - i
            if( jz >= 1 ) then
               if( denstot(jz) /= ZERO )  exit
            end if
         end do
      else
         jz = iz
      end if

   end subroutine getNearestCell1dReg
!
!=======================================================================
!
   subroutine getOtherNearCellIndices1dReg( z, zmin2, zmax1,  &
!             ============================  i    i      i
                                            proj, posc, ffp, pmax, ncn, icnear )
!                                            i     i     i    i     o     o
!
!     Determine cells other than proj%id that are within the free flight path cylinder
!     (= cylinder with height ffp and radius pmax)
!
!     Only cells in the range (1:nz) are returned.
!
!  Arguments:
!
      real(REAL8), intent(in) :: z(0:)          ! positions of grid points
      real(REAL8), intent(in) :: zmin2          ! smallest z larger than z(0)
      real(REAL8), intent(in) :: zmax1          ! largest z smaller than z(nz)
      type(Proj_type), intent(in) :: proj       ! projectile
      real(REAL8), intent(in) :: posc(3)        ! position of end of free flight path
      real(REAL8), intent(in) :: ffp            ! free flight path (= height of cylinder)
      real(REAL8), intent(in) :: pmax           ! maximum impact parameter (= radius of cyl.)
      integer(INT4), intent(out) :: ncn         ! number of regions found
      integer(INT4), allocatable, intent(inout) :: icnear(:)   
                                                ! index of regions found except proj%id
!
!  Local Variables:
!
      integer(INT4) :: nz, i, iz, izmin, izmax
      real(REAL8) :: wz, dz, z1, z2, zmin, zmax, distother, zcen, ffpcen(3)
!
!  Start of Executable Code:
!  ------------------------
!
      nz = ubound( z, 1 )
      wz = z(1) - z(0)
!
      iz = proj%id                              ! cell to be excluded from the search
!
!  First do a rough check whether other cells are far from the projectile
!
      if( proj%pos(3) < z(1) ) then
         distother = z(1) - proj%pos(3)
      else if( proj%pos(3) > z(nz-1) ) then
         distother = proj%pos(3) - z(nz-1)
      else
         distother = min( z(iz) - proj%pos(3), proj%pos(3) - z(iz-1) )
      end if

      if( distother**2 > ffp**2 + pmax**2 ) then
         ncn = 0
         return
      end if
!
!  Determine minimum and maximum z-values of the ffp cylinder or the sphere circumscribing
!  the ffp cylinder
!
      if( geom%lcart ) then
         dz = pmax * sqrt( proj%dir(1)**2 + proj%dir(2)**2 )
         z1 = proj%pos(3) - sign( dz, proj%dir(3) )
         z2 = posc(3) + sign( dz, proj%dir(3) )
         zmin = min( z1, z2 )
         zmax = max( z1, z2 )
      else
         ffpcen(:) = HALF * ( proj%pos(:) + posc(:) )
         call toCurviCoord1d( ffpcen(:), zcen )
!             ==============      i       o
         dz = sqrt( (HALF*ffp)**2 + pmax**2 )
         zmin = zcen - dz
         zmax = zcen + dz
      end if

      zmin = max( zmin2, min( zmax1, zmin ) )          ! don't apply bounds to izmin,izmax
      zmax = max( zmin2, min( zmax1, zmax ) )          !    because of integer overflow
      izmin = ceiling( ( zmin - z(0) ) / wz )          ! izmin,izmax are in [1,nz]
      izmax = ceiling( ( zmax - z(0) ) / wz )
!
!  Make sure icnear array is big enough
!
      if( size( icnear ) < nz-1 )  &
            call changeShape( icnear, nz-1, 'icnear', 'getOtherNearCellIndices1dReg' )
!                ===========    io     i       i                     i
!
!  Set icnear array
!
      if( iz < izmin .or. iz > izmax ) then
         print *, 'izmin, izmax, iz = ', izmin, izmax, iz
         error stop 'getOtherNearCellInidices1dReg: cell index error'
         ncn = izmax - izmin + 1
         icnear(1:ncn) = [ ( i, i=izmin,izmax ) ]        ! ???
      else
         ncn = izmax - izmin
         icnear(1:ncn) = [ ( i, i=izmin,iz-1 ), ( i, i=iz+1,izmax ) ]
      end if

   end subroutine getOtherNearCellIndices1dReg
!
!=======================================================================
!
   subroutine fixNegDensCell1dReg( dens )
!             ===================   io
!
!     Eliminate negative cell densities by drawing material from neighboring cells.
!     Upon exit, densities less than zero on entry, are set to zero.
!
!     OpenMP: This subroutine must be called from the master thread
!
!  Arguments:
!
      real(REAL8), intent(inout) :: dens(:,:)   ! (nelem,nc) densities of atoms
!
!  Local Variables:
!
      logical :: mask(2)
      integer(INT4) :: ielem, nelem_, iz, nz, k, jz(2)
      real(REAL8) :: sum_dens
!
!  Start of Executable Code:
!  ------------------------
!
      nelem_ = size( dens, 1 )
      nz = size( dens, 2 )
      if( nelem_ /= nelem )  error stop 'fixNegDensCell1dReg: Inconsistent nelem'

      do iz = 1, nz
         do ielem = 1, nelem
            if( dens(ielem,iz) >= ZERO )  cycle
            print '(a)', 'fixNegDensitiesCell1dReg: dens(' //  &
                         str( ielem, 'i3' ) // ',' // str( iz, 'i9' ) // ')=' //  &
                         str( dens(ielem,iz), '1es10.3' )
            do k = 1, nz
               jz(:) = [ iz-k, iz+k ]
               mask(:) = ( jz(:) >= 1 .and. jz(:) <= nz )
               if( .not.any( mask(:) ) ) then
                  print *, 'WARNING: Negative density not completely eliminated'
                  exit
               end if

               ! avoid redistributing negative densities
               where( mask(:) )  mask(:) = ( dens(ielem,jz(:)) > ZERO )
               sum_dens = sum( dens(ielem,jz(:)), mask(:) )

               if( sum_dens > -dens(ielem,iz) ) then
                  where( mask(:) )  &
                     dens(ielem,jz(:)) = dens(ielem,jz(:)) * (ONE + dens(ielem,iz)/sum_dens)
                  dens(ielem,iz) = ZERO
                  exit
               else
                  where( mask(:) )  dens(ielem,jz(:)) = ZERO
                  dens(ielem,iz) = dens(ielem,iz) + sum_dens
               end if
            end do
         end do
      end do

   end subroutine fixNegDensCell1dReg
!
!=======================================================================
!
   subroutine diffuseCell1dReg( z, dens, dpa )
!             ================  i   io    i
!
!     Diffuse dens(posif) with diffusivity induced by dpa
!     The concentrations of the diffusing species in the first and last material cell are
!     set to zero.
!
!  Arguments:
!
      real(REAL8), intent(in) :: z(0:)          ! (0:nz) cell boundaries
      real(REAL8), intent(inout) :: dens(:,:)   ! (nelem,nz) cell concentrations
      real(REAL8), intent(in) :: dpa(:)         ! (nz) dpa's
!
!  Local Variables:
!
      integer(INT4) :: mz, nz, ielem, izmin, izmax
      real(REAL8) :: wz, fac
      real(REAL8), allocatable :: a(:), b(:), al(:), am(:), au(:)
!
!  Start of Executable Code:
!  ------------------------
!
      nz = ubound( z, 1 )
      wz = z(1) - z(0)
!
!  Determine first and last cell with material
!
      do izmin = 1, nz
         if( sum( dens(:,izmin) ) > ZERO )  exit
      end do
      do izmax = nz, 1, -1
         if( sum( dens(:,izmax) ) > ZERO )  exit
      end do

      mz = izmax - izmin + 1
      call alloc( a, mz, 'a', 'diffuseCell1dReg' )
!          =====  o  i    i           i
      call alloc( b, mz, 'b', 'diffuseCell1dReg' )
!          =====  o  i    i           i
      call alloc( al, mz, 'al', 'diffuseCell1dReg' )
!          =====  o   i    i            i
      call alloc( am, mz, 'am', 'diffuseCell1dReg' )
!          =====  o   i    i            i
      call alloc( au, mz, 'au', 'diffuseCell1dReg' )
!          =====  o   i    i            i

      do ielem = 1, nelem

         if( dam%diffusivity(ielem) == ZERO )  cycle
!
!  Calculate coefficients of current densities: J(i) = a(j)*C(i) - b(j)*C(i+1)
!  J(i) is the current density from cell i into cell i+1, where j=i-izmin+1
!
         fac = dam%diffusivity(ielem) / wz
         a(1:mz-1) = fac * ( (ONE-dam%fracdest(ielem))*dpa(izmin:izmax-1) +  &
                                    dam%fracdest(ielem)*dpa(izmin+1:izmax) )
         b(1:mz-1) = fac * ( (ONE-dam%fracdest(ielem))*dpa(izmin+1:izmax) +  &
                                    dam%fracdest(ielem)*dpa(izmin:izmax-1) )
         a(mz) = ZERO
         b(mz) = ZERO                              ! unused
!
!  Calculate coefficients of density changes dC(i):
!
!     wz*dC(i) = -(J(i)-J(i-1)) = al(i)*C(i-1) + am(i)*C(i) + au(i)*C(i+1)
!
         al(1) = ZERO                              ! unused
         am(1) = ZERO                              ! Dirichlet boundary condition
         au(1) = ZERO                              ! Dirichlet boundary condition
         al(2:mz-1) = a(1:mz-2)
         am(2:mz-1) = - a(2:mz-1) - b(1:mz-2)
         au(2:mz-1) = b(2:mz-1)
         al(mz) = ZERO                             ! Dirichlet boundary condition
         am(mz) = ZERO                             ! Dirichlet boundary condition
         au(mz) = ZERO                             ! unused
!
!  Modify coefficients so they describe the linear system obtained by setting dC(i)=C(i)-C0(i)
!  on the left-hand side, where C0 are the initial densities, and leaving C(i) unknowns on
!  the right-hand side (Backward Euler discretization). The resulting equations are
!
!     al(i)*C(i-1) + am(i)*C(i) + au(i)*C(i+1) = b(i)
!
         al(:) = -al(:)
         am(:) = wz - am(:)
         au(:) = -au(:)
         b(1) = ZERO
         b(2:mz-1) = wz * dens(ielem,izmin+1:izmax-1)
         b(mz) = ZERO
!
!  Solve the equations
!
         a(:) = dens(ielem,izmin:izmax)
         call tridiag( al, am, au, b, a )
!             =======  i   i   i   i  o
         dens(ielem,izmin:izmax) = a(:)

      end do

   end subroutine diffuseCell1dReg
!
!=======================================================================
!
   subroutine relaxCell1dReg( z,  &
!             ==============  io
                              dens, damdens, dens0, ffgeom, lxtal, dpa, pressure, stress,  &
!                              io     io      io      io     io    io      io       io
                              istep )
!                               i
!
!     Redistribute the cell densities so that stress relaxes.
!     Either a viscoelastic or a full relaxation model is employed.
!
!     This subroutine should not be called too often in order to avoid artificial diffusion,
!     but often enough to keep the surface sufficiently accurate.
!
!     OpenMP: This subroutine must be called from the master thread
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: z(:)          ! positions of grid points
      real(REAL8), allocatable, intent(inout) :: dens(:,:)     ! (nelem,nc) densities of atoms
                                                               !    on entry: before relax
                                                               !    on exit: after relax
      real(REAL8), allocatable, intent(inout) :: damdens(:,:)  ! (nelemIn,nc) damage densities
      real(REAL8), allocatable, intent(inout) :: dens0(:,:)    ! (nelem,nc) densities of atoms
                                                               !    after previous relax
      real(REAL8), allocatable, intent(inout) :: ffgeom(:)     ! geometric fill fraction:
                                                               !    on entry: after last relax
                                                               !    on exit: after this relax
      logical, allocatable, intent(inout) :: lxtal(:)          ! crystallinity flag
      real(REAL8), allocatable, intent(inout) :: dpa(:)        ! displacements per atom
      real(REAL8), allocatable, intent(inout) :: pressure(:)   ! pressure
      real(REAL8), allocatable, intent(inout) :: stress(:,:,:) ! zz component of stress tensor
      integer(INT4), intent(in) :: istep                       ! step counter
!
!  Local Variables:
!
      integer(INT4) :: iz, izmin, izmax, nz, nz2, nelemIn, nzadd
      real(REAL8) :: wz
      real(REAL8), allocatable :: z_tmp(:)
!
!  Start of Executable Code:
!  ------------------------
!
      nelemIn = 2*nelem
!
!  Add empty frame to allow swelling
!
      nz = ubound( z, 1 )
      wz = z(1) - z(0)

      nzadd = 1               ! layers to be added (one or two)

      call alloc( z_tmp, 0, nz+nzadd, 'z_tmp', 'relaxCell1dReg' )
!          =====   io    i   i       i             i
      if( nzadd == 2 ) then
         z_tmp(0:) = [ z(0)-wz, z(:), z(nz)+wz ]
      else
         z_tmp(0:) = [ z(0)-wz, z(:) ]
      end if

      call changeShape( ffgeom, nz+nzadd, 'ffgeom', 'relaxCell1dReg', shift_=1 )
!          ===========    io       i         i               i           i
      call changeShape( lxtal, nz+nzadd, 'lxtal', 'relaxCell1dReg', shift_=1 )
!          ===========   io       i         i              i            i
      call changeShape( dpa, nz+nzadd, 'dpa', 'relaxCell1dReg', shift_=1 )
!          ===========  io      i        i            i            i
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         call changeShape( pressure, nz+nzadd, 'pressure', 'relaxCell1dReg', shift_=1 )
!             ===========     io        i          i               i            i
         call changeShape( stress, [1,1,nz+nzadd], 'stress', 'relaxCell1dReg', shift_=[0,0,1] )
!             ===========    io          i            i               i            i
      end if

      call changeShape( dens, [nelem, nz+nzadd], 'dens', 'relaxCell1dReg', shift_=[0,1] )
!          ===========   io           i            i             i              i
      !dens(:,1) = dens(:,2)
      !dens(:,nz+2) = dens(:,nz+1)

      if( out%lhis ) then
         call changeShape( damdens, [nelemIn, nz+nzadd], 'damdens', 'relaxCell1dReg', &
!             ===========     io             i               i              i
                           shift_=[0,1] )
!                               i
         !damdens(:,1) = damdens(:,2)
         !damdens(:,nz+2) = damdens(:,nz+1)
      end if

      call changeShape( dens0, [nelem, nz+nzadd], 'dens0', 'relaxCell1dReg', shift_=[0,1] )
!          ===========   io           i              i             i              i
      !dens0(:,1) = dens0(:,2)
      !dens0(:,nz+2) = dens0(:,nz+1)
!
!  Relax the fill fractions
!
#ifdef AMAT
      call relaxViscoElasticEuler1d( z_tmp(:), dens(:,:), damdens(:,:), dens0(:,:), &
!          ========================    io        io           io           io
                                     ffgeom(:), dpa(:),  &
!                                       io       io
                                     pressure(:), stress(1,1,:), istep )
!                                        io            io          i
#else
      error stop 'ViscoElasticEuler1d model requires separate license'
#endif
!
!  Remove empty cells on both sides
!
      nz2 = nz + nzadd

      do iz = 1, nz2
         if( any( dens(:,iz) > ZERO ) )  exit
      end do
      izmin = iz

      do iz = nz2, 1, -1
         if( any( dens(:,iz) > ZERO ) )  exit
      end do
      izmax = iz

      nz = izmax - izmin + 1

      deallocate( z )
      call alloc( z, 0, nz, 'z', 'relaxCell1dReg' )
!          =====  io i  i    i           i
      z(0:) = z_tmp(izmin-1:izmax)

      call changeShape( ffgeom, nz, 'ffgeom', 'relaxCell1dReg', shift_=1-izmin )
!          ===========    io    i      i               i              i
      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         call changeShape( pressure, nz, 'pressure', 'relaxCell1dReg', shift_=1-izmin )
!             ===========     io     i       i                i              i
         call changeShape( stress, [1,1,nz], 'stress', 'relaxCell1dReg', shift_=[0,0,1-izmin] )
!             ===========    io       i         i               i              i
      end if

      call changeShape( lxtal, nz, 'lxtal', 'relaxCell1dReg', shift_=1-izmin )
!          ===========   io    i      i              i              i
      call changeShape( dens, [nelem, nz], 'dens', 'relaxCell1dReg', shift_=[0,1-izmin] )
!          ===========   io        i         i              i               i
      call changeShape( dens0, [nelem, nz], 'dens0', 'relaxCell1dReg', shift_=[0,1-izmin] )
!          ===========   io         i          i              i               i
      if( out%lhis )  &
            call changeShape( damdens, [nelemIn, nz], 'damdens', 'relaxCell1dReg',  &
!                ===========     io            i          i              i
                              shift_=[0,1-izmin] )
!                                    i
!
!  Update surface positions
!
#ifdef AMAT
      call initSurface1dFromFF( z(:), ffgeom(:), dens(:,:), istep )
!          ===================   i        i        i          i
#else
         error stop 'ViscoElasticEuler1d model requires separate license'
#endif

   end subroutine relaxCell1dReg
!
!=======================================================================
!
   subroutine relaxCell1dRegFull( z, zmin1, zmax2,  &
!             ==================  io   i      i
                                  dens, damdens, ffgeom, lxtal )
!                                  io     io       io     io
!
!     Redistribute the cell densities so that stress relaxes.
!     This is done by shuffling cell contents to the left and right of posfix.
!
!     Cells are removed so that no empty cells remain (this is not possible in 2-d!).
!     Currently no cells are added (no swelling!).
!
!     Initializes the surface.
!
!     NOTE: ff does NOT contain useful data on exit, but must be declared inout as it is
!     modified in the subroutine.
!
!  Arguments:
!
      real(REAL8), allocatable, intent(inout) :: z(:)          ! positions of grid points
      real(REAL8), intent(in) :: zmin1                         ! largest z smaller than z(0)
      real(REAL8), intent(in) :: zmax2                         ! smallest z larger than z(nz)
      real(REAL8), allocatable, intent(inout) :: dens(:,:)     ! (nelem,nc) densities of atoms
      real(REAL8), allocatable, intent(inout) :: damdens(:,:)  ! (nelemIn,nc) damage densities
      real(REAL8), allocatable, intent(inout) :: ffgeom(:)     ! geometric fill fraction:
                                                               !    on entry: after last relax
                                                               !    on exit: after this relax
      logical, allocatable, intent(inout) :: lxtal(:)          ! crystallinity flag
!
!  Local Variables:
!
      logical :: outside, empty, full_new, isEmpty
      logical, allocatable :: lxtalnew(:)
      integer(INT4) :: nz, iz, izfix, jz, jzmin, jzmax, nelemIn, ielem
      real(REAL8) :: wz, zmin, zmax, z_tmp
      real(REAL8) :: densfix_left(nelem), densfix_right(nelem)
      real(REAL8) :: ddens(nelem), ddam(2*nelem)
      real(REAL8) :: posfix, z0, fac, ff_tot, ff_des, fffix_left, fffix_right, dff
      real(REAL8), allocatable :: densnew(:,:), damnew(:,:), edamnew(:), ff(:), ffnew(:)
      type(Real8Stack) :: zsurfStack   ! stack of surface positions
!
!  Start of Executable Code:
!  ------------------------
!
      nelemIn = 2*nelem

      nz = ubound( z, 1 )
      wz = z(1) - z(0)

      call alloc( ff, nz, 'ff', 'relaxCell1dRegFull' )
!          =====  io  i    i             i
      ff(:) = ZERO
      do ielem = 1, nelem
         ff(:) = ff(:) + dens(ielem,:) / elem(ielem)%dens
      end do
!
!  Determine position posfix which remains fixed. We will update the grid starting at the cell
!  containing posfix
!
      if( geom%posfix > z(0) .and. geom%posfix < z(nz) ) then
         posfix = geom%posfix
      else
         posfix = HALF * (z(0) + z(nz))         ! Is this necessary?
      end if
                                       ! don't apply bounds to iz because of integer overflow
      z_tmp = max( zmin1, min( zmax2, posfix ) )
      izfix = ceiling( ( z_tmp - z(0) ) / wz )
      outside = ( izfix < 1 .or. izfix > nz )
      izfix = max( 1, min( nz, izfix ) )
!
!  Divide ff and dens of izfix into contributions left and right of posfix
!  fac is the desired contribution to the fill fraction of the space of izfix left of posfix
!  ONE-fac is the corresponding contribution of the space of izfix left of posfix
!  For all other cells the desired fill fraction is ONE
!
      if( outside ) then
         if( izfix == 1 ) then
            fac = ZERO
         else if( izfix == nz ) then
            fac = ONE
         else
            stop 'relaxCell1dRegFull: Inconsistent value of izfix'
         end if
      else
         fac = ( posfix - z(izfix-1) ) / wz
      end if

      densfix_left(:) = dens(:,izfix) * fac
      fffix_left = ff(izfix) * fac

      densfix_right(:) = dens(:,izfix) * (ONE - fac)
      fffix_right = ff(izfix) * (ONE - fac)
!
!  Prepare arrays for new fill fractions and cell contents
!
      ff_tot = sum( ff(1:izfix-1) ) + (fffix_left - fac)
      jzmin = izfix - ceiling( ff_tot ) - 1     ! safety margin of 1
      ff_tot = sum( ff(izfix+1:nz) ) + (fffix_right - (ONE-fac))
      jzmax = izfix + ceiling( ff_tot ) + 1     ! safety margin of 1

      call alloc( ffnew, jzmin, jzmax, 'ffnew', 'relaxCell1dRegFull' )
!          =====   io      i      i       i               i
      ffnew(:) = ZERO                           ! fill fraction in the current new cell
      call alloc( densnew, [1, jzmin], [nelem, jzmax], 'densnew', 'relaxCell1dRegFull' )
!          =====    io         i             i             i                i
      densnew(:,:) = ZERO                       ! new cell densities
      if( out%lhis ) then
         call alloc( damnew, [1, jzmin], [nelemIn, jzmax], 'damnew', 'relaxCell1dRegFull' )
!             =====    io        i              i             i                i
         damnew(:,:) = ZERO                     ! new cell damage
      end if
      call alloc( edamnew, jzmin, jzmax, 'edamnew', 'relaxCell1dRegFull' )
!          =====    io       i      i        i                i
      edamnew(:) = ZERO                          ! new cell damage energy
      call alloc( lxtalnew, jzmin, jzmax, 'lxtalnew', 'relaxCell1dRegFull' )
!          =====     io       i      i        i                 i
      lxtalnew(:) = .false.                     ! new crystallinity flag
!
!  We empty the old density array dens and fill the new density array densnew so that the new
!  fill fractions equal unity
!
!  Shuffle cell contents to the left of posfix from old grid to new grid
!
      ff(izfix) = fffix_left
      dens(:,izfix) = densfix_left(:)

      iz = izfix                                ! grid index of source cell
      jz = izfix                                ! grid index of target cell

      ff_des = fac                              ! desired fill new factor of jz

      do
         empty = (ff(iz) <= ff_des - ffnew(jz))   ! cell iz contents fit into cell jz
         full_new = (ff(iz) >= ff_des - ffnew(jz))! cell jz full after shuffling from cell iz
                                                  !    MIND the "=" in both cases
         if( empty ) then
            dff = ff(iz)
            ddens(:) = dens(:,iz)
            if( out%lhis )  ddam(:) = damdens(:,iz)
         else
            dff = ff_des - ffnew(jz)
            ddens(:) = dens(:,iz) * dff / ff(iz)  ! ddens and ddam are proportional to dff
            if( out%lhis )  ddam(:) = damdens(:,iz) * dff / ff(iz)
         end if

         ffnew(jz) = ffnew(jz) + dff
         ff(iz) = ff(iz) - dff
         densnew(:,jz) = densnew(:,jz) + ddens(:)
         dens(:,iz) = dens(:,iz) - ddens(:)
         if( out%lhis ) then
            damnew(:,jz) = damnew(:,jz) + ddam(:)
            damdens(:,iz) = damdens(:,iz) - ddam(:)
         end if
         lxtalnew(jz) = lxtal(iz)          ! dirty

         if( full_new ) then
            jz = jz - 1
            if( jz < jzmin )  error stop 'relaxCell1dRegFull: WARNING: jz too small'
            ff_des = ONE
         end if
         if( empty )  iz = iz - 1
         if( iz < 1 )  exit
      end do
      jzmin = jzmin + 1                         ! remove safety margin
!
!  Shuffle cell contents to the right of posfix from old grid to new grid
!
      ff(izfix) = fffix_right
      dens(:,izfix) = densfix_right(:)

      iz = izfix
      jz = izfix

      do
         empty = (ff(iz) <= ONE - ffnew(jz))    ! cell iz contents fit into cell jz
         full_new = (ff(iz) >= ONE - ffnew(jz)) ! cell jz full after shuffling from cell iz
                                                !    MIND the "=" in both cases
         if( empty ) then
            dff = ff(iz)
            ddens(:) = dens(:,iz)
            if( out%lhis )  ddam(:) = damdens(:,iz)
         else
            dff = ONE - ffnew(jz)
            ddens(:) = dens(:,iz) * dff / ff(iz)  ! ddens is proportional to dff
            if( out%lhis )  ddam(:) = damdens(:,iz) * dff / ff(iz)
         end if

         ffnew(jz) = ffnew(jz) + dff
         ff(iz) = ff(iz) - dff
         densnew(:,jz) = densnew(:,jz) + ddens(:)
         dens(:,iz) = dens(:,iz) - ddens(:)
         if( out%lhis ) then
            damnew(:,jz) = damnew(:,jz) + ddam(:)
            damdens(:,iz) = damdens(:,iz) - ddam(:)
         end if
         lxtalnew(jz) = lxtal(iz)          ! dirty

         if( full_new ) then
            jz = jz + 1
            if( jz > jzmax )  error stop 'relaxCell1dRegFull: WARNING: jz too large'
         end if
         if( empty )  iz = iz + 1
         if( iz > nz )  exit
      end do
      jzmax = jzmax - 1                         ! remove safety margin
!
!  Deallocate old grid and cells and create a new grid
!
      z0 = z(0)
      nz = jzmax - jzmin + 1
      deallocate( z )
      call alloc( z, 0, nz, 'z', 'relaxCell1dRegFull' )
!          =====  io i  i    i            i
      z(:) = [( z0 + jz*wz, jz=jzmin-1,jzmax )]
!
!  Allocate and set cell contents
!
      deallocate( dens )
      call alloc( dens, [nelem, nz], 'dens', 'relaxCell1dRegFull' )
!          =====   io        i         i              i
      dens(:,:) = densnew(:,jzmin:jzmax)

      if( out%lhis ) then
         deallocate( damdens )
         call alloc( damdens, [nelemIn, nz], 'damdens', 'relaxCell1dRegFull' )
!             =====     io           i           i                i
         damdens(:,:) = damnew(:,jzmin:jzmax)
      end if

      deallocate( lxtal )
      call alloc( lxtal, nz, 'lxtal', 'relaxCell1dRegFull' )
!          =====    o    i      i              i
      lxtal(:) = lxtalnew(jzmin:jzmax)
!
!  Determine fill fractions
!
      deallocate( ffgeom )
      call alloc( ffgeom, nz, 'ffgeom', 'relaxCell1dRegFull' )
!          =====     o    i      i                i
      do iz = 1, nz
         ffgeom(iz) = sum( dens(:,iz) / elem(:)%dens )
      end do
!
!  Determine surfaces
!
      call newStack( zsurfStack )
!          ========      o
      do iz = 1, nz
         if( ffgeom(iz) == ZERO .or.  &
             ( abs(ffgeom(iz) - ONE) < EPSMA .and. iz /= 1 .and. iz /= nz ) )  cycle
         if( iz == 1 ) then
            call pushStack( z(iz)-wz*ffgeom(iz), zsurfStack )
!                =========         i                 io
         else if( iz == nz ) then
            call pushStack( z(iz-1)+wz*ffgeom(iz), zsurfStack )
!                =========           i                 io
         else
            if( abs( ffgeom(iz+1) - ffgeom(iz-1) ) < EPSMA ) then
               call pushStack( z(iz-1)+HALF*wz*ffgeom(iz), zsurfStack )
!                   =========               i                  io
               call pushStack( z(iz)-HALF*wz*ffgeom(iz), zsurfStack )
!                   =========              i                 io
            else if( ffgeom(iz+1) > ffgeom(iz-1) ) then
               call pushStack( z(iz)-wz*ffgeom(iz), zsurfStack )
!                   =========          i                io
            else
               call pushStack( z(iz-1)+wz*ffgeom(iz), zsurfStack )
!                   =========           i                 io
            end if
         end if
      end do
!
!  Reinitialize the surface (for now, only two surfaces are allowed)
!
      if( sizeOfStack( zsurfStack ) /= 2 ) then
!         ===========      i
         print *, 'stack size=', sizeOfStack( zsurfStack )
         isEmpty = .false.
         do while( .not.isEmpty )
            call popStack( zsurfStack, zmin, isEmpty )
!                ========      io       o       o
            print *, 'zsurf=', zmin
         end do
         error stop 'relaxCell1dRegFull: more than 2 surfaces'
      end if

      call popStack( zsurfStack, zmax, isEmpty )
!          ========      io       o       o
      call popStack( zsurfStack, zmin, isEmpty )
!          ========      io       o       o
!
!  Initialize the Surface
!
      call updateSurface1d( zmin, zmax )
!          ===============   i     i

   end subroutine relaxCell1dRegFull
!
!=======================================================================
!
   subroutine writeCell1dReg( z, dens, pressure, stress )
!             ==============  i   i       i        i
!
!     Write the contents of the target cells to the CELL file (same format as HIS file)
!
!  Arguments:
!
      real(REAL8), intent(in) :: z(0:)                      ! positions of grid points
      real(REAL8), intent(in) :: dens(:,:)                  ! (nelem,nz) densities
      real(REAL8), intent(in), allocatable :: pressure(:)   ! (nz) pressure
      real(REAL8), intent(in), allocatable :: stress(:,:,:) ! (1,1,nz) zz component of stress
!                                                           !    tensor
!  Local Variables:
!
      type(His1d_type) :: his
      integer(INT4) :: nz, iz
      real(REAL8) :: wz
!
      character(6), save :: fmode = 'write'
!
!  Start of Executable Code:
!  ------------------------
!
      nz = ubound( z, 1 )
      wz = z(1) - z(0)

      his%nbox = nz
      his%imin = 1
      his%imax = nz
      his%hbeg = z(0)
      his%wbox = wz

      if( dam%relaxmodel == ENUM%relaxmodel%viscoelastic ) then
         call alloc( his%hist, [nelem+2, nz], 'his%hist', 'writeCell1dReg' )
!             =====     io           i             i               i
         do iz = 1, nz
            his%hist(1:nelem,:) = dens(:,:)
            his%hist(nelem+1,:) = pressure(:) / 1e24     ! compensate for scaling
            his%hist(nelem+2,:) = stress(1,1,:) / 1e24   ! compensate for scaling
         end do

         call alloc( his%indhis, nelem+2, 'his%indhis', 'writeCell1dReg' )
!             =====      io         i           i               i
         his%indhis(1:nelem) = elem(1:nelem)%z
         his%indhis(nelem+1:nelem+2) = 0
      else
         call alloc( his%hist, [nelem, nz], 'his%hist', 'writeCell1dReg' )
!             =====     io          i           i               i
         his%hist(:,:) = dens(:,:)

         call alloc( his%indhis, nelem, 'his%indhis', 'writeCell1dReg' )
!             =====      io        i          i               i
         his%indhis(1:nelem) = elem(1:nelem)%z
      end if
!
      call writeHis1d( 'CELL', fmode, 1e24_REAL8, his, -BIG )
!          ==========    i       i        i        i     i

      deallocate( his%hist )
      deallocate( his%indhis )

      fmode = 'append'

   end subroutine writeCell1dReg
!
!=======================================================================
!
   subroutine writeDamCell1dReg( z, damdens )
!             =================  i     i
!
!     Write the damage of the target cells to the HIS file
!
!  Arguments:
!
      real(REAL8), intent(in) :: z(0:)          ! positions of grid points
      real(REAL8), intent(in) :: damdens(:,:)   ! (nelemIn,nz) damage densities
!
!  Local Variables:
!
      integer(INT4) :: nz, iz, ielem, nelemIn
      real(REAL8) :: wz
!
!  Start of Executable Code:
!  ------------------------
!
      nelemIn = 2*nelem

      nz = ubound( z, 1 )
      wz = z(1) - z(0)

      call openFile( 22, 'HIS', 'formatted', 'write' )
!          ========  i     i         i          i
!
!  Write header
!
      write( 22, * )  1+nelemIn, 2*nz+2
      write( 22, * )  0, elem(1:nelem)%z, -elem(1:nelem)%z
!
!  Write first point
!
      write( 22, '(f11.3, t13, 100es13.6)' )  &
           z(0), (1e24*HISMIN*out%fachis(1)/wz, ielem=1,nelemIn)
!
!  Write histogram values
!
      do iz = 1, nz
         write( 22, '(f11.3, t13, 100es13.6)' )  &
            z(iz-1),  &
            ( 1e24 * max( HISMIN*out%fachis(1)/wz, damdens(ielem,iz) ), ielem=1,nelemIn )
         write( 22, '(f11.3, t13, 100es13.6)' )  &
            z(iz),  &
            ( 1e24 * max( HISMIN*out%fachis(1)/wz, damdens(ielem,iz) ), ielem=1,nelemIn )
      end do
!
!  Write last point
!
      write( 22, '(f11.3, t13, 100es13.6)' )  &
           z(nz), (1e24 * HISMIN*out%fachis(1)/wz, ielem=1,nelemIn)

      close( 22 )

   end subroutine writeDamCell1dReg
!
!=======================================================================
!
   subroutine readHis1dReg( z, damdens )
!             ============  i    io
!
!     Read the damage concentrations from the hisfile file (HIS file format)
!
!  Arguments:
!
      real(REAL8), intent(in) :: z(0:)                         ! positions of grid points
      real(REAL8), allocatable, intent(inout) :: damdens(:,:)  ! (nelemIn,nc) damage densities
!
!  Local Variables:
!
      integer(INT4) :: ncol, nrow, idummy, nz, nz_tmp, nelemIn
      integer(INT4) :: iz, i
      integer(INT4) :: jz, izmin, izmax
      real(REAL8) :: zminr, zmaxr
      real(REAL8) :: wz, rdummy
      real(REAL8), allocatable :: z_tmp(:), damdens_tmp(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
      nelemIn = 2*nelem

      nz = ubound( z, 1 )
      wz = z(1) - z(0)
!
!  Open file
!
      call openFile( 22, hisfile, 'formatted', 'read' )
!          ========  i     i         i          i
!
!  Read first line and allocate memory
!
      read( 22, *, err=900, end=900 ) ncol, nrow

      if( ncol-1 /= nelemIn )  error stop 'readHis1dNU: inconsistent number of columns on file'

      nz_tmp = ( nrow - 2 ) / 2
      call alloc( z_tmp, 0, nz_tmp, 'damdens_tmp', 'readHis1dReg' )
!          =====   io    i    i           i              i
      call alloc( damdens_tmp, [nelemIn, nz_tmp], 'damdens_tmp', 'readHis1dReg' )
!          =====       io              i                i              i
!
!  Second and third line can be ignored
!
      read( 22, *, err=900, end=900 )  ( idummy, i=1,ncol )
      read( 22, *, err=900, end=900 )  ( rdummy, i=1,ncol )
!
!  Read histogram values
!
      do iz = 1, nz_tmp
         read( 22, *, err=900, end=900 )  z_tmp(iz-1), damdens_tmp(:,iz)
         read( 22, *, err=900, end=900 )  z_tmp(iz), damdens_tmp(:,iz)
      end do
      damdens_tmp(:,:) = damdens_tmp(:,:) / 1e24
!
!  Ignore last line
!
      close( 22 )
!
!  Shift and reflect the distribution if desired
!
      if( defined( shiftto ) )  &
         z_tmp(:) = shiftto + z_tmp(:) - z_tmp(0)

      if( reflect ) then
         z_tmp(:) = - z_tmp(nz_tmp:0:-1)
         damdens_tmp(:,:) = damdens_tmp(:,nz:1:-1)
      end if
!
!  Interpolate to grid
!
      if( allocated( damdens ) )  deallocate( damdens )
      call alloc( damdens, [nelemIn, nz], 'damdens', 'readHis1dReg' )
!          =====     io           i           i            i
      damdens(:,:) = ZERO

      do jz = 1, nz_tmp
         zminr = ( z_tmp(jz-1) - z(0) ) / wz
         zmaxr = ( z_tmp(jz) - z(0) ) / wz

         izmin = max( 1, ceiling( zminr ) )
         izmax = max( 1, ceiling( zmaxr ) )

         do iz = izmin, izmax
            damdens(:,iz) = damdens(:,iz) + damdens_tmp(:,jz)
         end do

         damdens(:,izmin) = damdens(:,izmin) - damdens_tmp(:,jz) * ( zminr - (izmin-1) )
         damdens(:,izmax) = damdens(:,izmax) - damdens_tmp(:,jz) * ( izmax - zmaxr )
      end do

      return

900   error stop 'readHis1dReg: error during reading histogram file'

   end subroutine readHis1dReg
!
!=======================================================================
!
end module Target_Cell1d_Reg_mod

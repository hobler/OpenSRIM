module Param_Mat_mod
!
!  Holds data defining the materials. The material is defined separately in each region. 
!  For region 0 (vacuum) only mat(0)%lxtal=f will be defined.
!
!  Uses public data (mat structure array) and private namelist (&material)
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, ONE, FOUR, EIGHT, THIRD, QUART, PI, BIG, EPSMB
   use Undefine_mod,       only: undefine, defined
   use String_mod,         only: str
   use Vector_mod,         only: operator(.dot.), operator(.cross.), abs
   use Open_mod,           only: openFile

   use Param_Setup_mod,    only: nr, nelemIon, nelem, ndim, nthread, temp, ldyn, usecell
   use Param_Geom_mod,     only: geom
   use Param_Elem_mod,     only: elem, anaMolecule
   use Param_Ion_mod,      only: ion
   use Param_Xtal_mod,     only: xtal_struct => xtal

   use Geom_CoordSys_mod,  only: CoordSys_type, initCoordSys, invertCoordSys, fromCoordSys
   use Math_Debye_mod,     only: debye1
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
   public :: readMat                ! read &material record from input file
!            =======                !    and initialize the mat structure
   public :: initMat                ! initialize the mat structures
!            =======                !
   public :: writeMat               ! write data derived from &material
!            ========               !    parameters to output file
!
!  Public Data:
!  ===========
!
   type :: Mat_type                          ! Material data type
!          ========

      logical :: l100                        ! unrotated crystal flag
      logical :: lpoly                       ! polycrystalline flag
      logical :: lsput                       ! sputtering flag
      logical :: lxtal                       ! crystal flag

      character(80) :: name                  ! (chemical) name of initial material
      character(80) :: xtal                  ! crystallographic system (unused)

      integer(INT4), allocatable :: ielem(:) ! (nelemr) atom indices
      integer(INT4) :: nelemr                ! number of atom species in region

      real(REAL8) :: dens                    ! atomic density [A-3]
      real(REAL8), allocatable :: ed(:)      ! (nelem) displacement energy [eV]
                                             !    will be initialized outside this module
      real(REAL8), allocatable :: esurf(:)   ! (nelem) surface binding energy [eV]
                                             !    will be initialized outside this module
      real(REAL8), allocatable :: prob(:)    ! (nelemr) probability to encounter atom in region
      real(REAL8) :: rotx                    ! rotation angle of x axis [deg]
      real(REAL8), allocatable :: syield(:)  ! (nelem) sputtering yield
      real(REAL8), allocatable :: tdebye(:)  ! (nelem) Debye temperature
      real(REAL8) :: xscale                  ! scaling length for some crystal parameters
      real(REAL8), allocatable :: vib(:)     ! (nelem) thermal vibration amplitude [A]
      real(REAL8) :: vsurf(3)                ! a vector parallel to the surface
      real(REAL8) :: wafer(3)                ! surface normal

   end type Mat_type

   type(Mat_type), allocatable, public :: mat(:)   ! (nr) materials structure

   integer(INT4), public :: irxtal                 ! should finally be eliminated
   type(CoordSys_type), public :: coord_xtal       ! crystal coordinate system
!$omp  threadprivate( coord_xtal )
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine readMat( )
!             =======
!
!     Read the &material records from the input file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: region                   ! region to which input corresponds
      integer(INT4) :: atom1                    ! projectile atom index
      integer(INT4) :: atom2                    ! target atom index

      logical :: lpoly                          ! polycrystalline flag

      character(80) :: name                     ! material
      character(80) :: xtal                     ! crystallographic system

      real(REAL8) :: density                    ! atomic density of the material [cm-3]
      real(REAL8) :: rotx                       ! rotation angle of x axis [deg]
      real(REAL8) :: syield                     ! sputtering yield
      real(REAL8) :: tdebye                     ! Debye temperature
      real(REAL8) :: vib                        ! thermal vibration amplitude [A]
      real(REAL8) :: vsurf(3)                   ! crystallographic direction along x
      real(REAL8) :: wafer(3)                   ! crystallographic direction along z

      namelist /material/ region, atom1, atom2,  &
                          lpoly,  &
                          name, xtal,  &
                          density, rotx, syield, tdebye, vib, vsurf, wafer

      logical :: lr(nr)
      character(2) :: ele
      character(256) :: msg
      integer(INT4) :: ios, ir, ielem, ielem1, nelemr, iz, natom2
      integer(INT4) :: ielem1min, ielem1max, ielem2min, ielem2max
!
!  Start of Executable Code:
!  ------------------------
!
!  Ignore &material if geometry and material information is read from a cell file
!  Note: irxtal must be set
!
      if( usecell ) then
         irxtal = 0
         return
      end if
!
!  Allocate memory for target regions
!
      allocate( mat( 0:nr+1 ) )

      do ir = 1, nr
         allocate( mat(ir)%ielem( nelem ) )       ! don't know nelemr in advance.
         allocate( mat(ir)%prob( nelem ) )        !    allocate with dimension nelem
         allocate( mat(ir)%tdebye( nelem ) )      !    since nelem >= nelemr
         allocate( mat(ir)%vib( nelem ) )
         allocate( mat(ir)%syield( nelem ) )
      end do
!
!  Set flag indicating definition of region ir
!
      lr(:) = .false.
!
!  Set default values of structure parameters
!
      mat(0)%lxtal = .false.                    ! should finally be eliminated
      mat(nr+1)%lxtal = .false.

      do ir = 1, nr
         mat(ir)%lpoly = .false.
         mat(ir)%lxtal = .false.

         call undefine( mat(ir)%name )

         call undefine( mat(ir)%dens )
         call undefine( mat(ir)%prob(:) )
         call undefine( mat(ir)%rotx )
         call undefine( mat(ir)%syield(:) )
         call undefine( mat(ir)%tdebye(:) )
         call undefine( mat(ir)%vib(:) )
         call undefine( mat(ir)%vsurf(:) )
         call undefine( mat(ir)%wafer(:) )

      end do
!
!  Loop over reading records
!  -------------------------
!
      call openFile( 11, 'INP', 'formatted', 'read' )
!          ========  i     i         i         i

      do
!
!  Set default values of namelist variables
!
         lpoly = .false.

         call undefine( region )
         call undefine( atom1 )
         call undefine( atom2 )

         call undefine( density )
         call undefine( name )
         call undefine( rotx )
         call undefine( syield )
         call undefine( tdebye )
         call undefine( vib )
         call undefine( vsurf )
         call undefine( wafer )
         call undefine( xtal ) 
!
!  Read &material record from the input file
!
         read( 11, nml=material, iostat=ios, iomsg=msg )
         if( ios < 0 )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 901               ! error during reading target record
         end if
!
!  Check if region is specified
!
         if( .not.defined( region ) )  goto 905
         ir = region
         if( ir < 1 .or. ir > nr )  goto 906

         lr(ir) = .true.
!
!  Check parameters on consistency with index parameter
!
         if( defined( atom1 ) .or. defined( atom2 ) ) then
!            =======               =======
            if( defined( name ) )  goto 915
            if( defined( density ) )  goto 910
            if( defined( rotx ) )  goto 914
            if( any( defined( vsurf ) ) )  goto 913
            if( any( defined( wafer ) ) )  goto 912
            if( defined( xtal ) )  goto  911
         end if            

         if( defined( atom1 ) ) then
!            =======
            if( defined( tdebye ) )  goto 922
            if( defined( vib ) )  goto 921
         end if            

         if( defined( atom2 ) ) then
!            =======
            if( defined( syield ) )  goto 923
         end if
!
!  Define range of atoms to which parameters apply
!
         if( defined( atom1 ) ) then
!            =======
            ielem1min = atom1
            ielem1max = atom1
         else
            ielem1min = 1
            ielem1max = nelem
         end if
         ielem1 = ielem1min

         if( defined( atom2 ) ) then
!            =======
            ielem2min = atom2
            ielem2max = atom2
         else
            ielem2min = 1
            ielem2max = nelem
         end if
!
!  Copy parameters to mat structure
!
         if( defined( name ) )  mat(ir)%name = name

         if( defined( xtal ) ) then
            if( xtal == 'on' .or. xtal == 'yes' ) then
               mat(ir)%lxtal = .true.
            else if( xtal == 'off' .or. xtal == 'no' ) then
               mat(ir)%lxtal = .false.
            else
               goto 950
            end if
         end if

         mat(ir)%lpoly = lpoly

         if( defined( density ) )  mat(ir)%dens = density / 1e24     ! 1/cm3 -> 1/A^3
         if( defined( rotx ) )  mat(ir)%rotx = rotx
         if( defined( syield ) )  mat(ir)%syield(ielem1min:ielem1max) = syield
         if( defined( tdebye ) )  mat(ir)%tdebye(ielem2min:ielem2max) = tdebye
         if( defined( vib ) )  mat(ir)%vib(ielem2min:ielem2max) = vib
         if( any( defined( vsurf(:) ) ) )  mat(ir)%vsurf(:) = vsurf(:)
         if( any( defined( wafer(:) ) ) )  mat(ir)%wafer(:) = wafer(:)

      end do
!
!  End of loop over reading records
!  --------------------------------
!
      close( 11 )
!
!  Check if all regions have been defined
!
      if( .not.all( lr(:) ) )  goto 900
!
!  Extract atoms and probabilities from material name
!  If they are not yet in the atoms list (elem(*)%name), insert them
!
      do ir = 1, nr

         if( .not.defined( mat(ir)%name ) )  goto 960
         name = mat(ir)%name                    ! name will be eaten away by anaMolecule
         nelemr = 0

         do                                     ! loop over atoms in region

            call anaMolecule( name, ele, iz, natom2)
!                ===========   io    o   o     o
            if( iz == 0 )  exit

            do ielem = nelemIon+1, nelem           ! check if atom already in list
               if( ele == elem(ielem)%name )  exit
            end do

            if( ielem > nelem ) then               ! new entry in atom list
               do ielem = nelemIon+1, nelem
                  if( .not.defined( elem(ielem)%name ) ) exit
               end do
               if( ielem > nelem )  goto 963
               elem(ielem)%name = ele
            end if

            nelemr = nelemr + 1                 ! new entry in region atom list
            mat(ir)%ielem(nelemr) = ielem
            mat(ir)%prob(nelemr) = natom2

         end do

         mat(ir)%nelemr = nelemr
         mat(ir)%prob(1:nelemr) = mat(ir)%prob(1:nelemr) / sum( mat(ir)%prob(1:nelemr) )

      end do

      return
!
!  Error messages
!
900   error stop '&material: not all regions defined'
901   error stop '&material: error during reading record'

905   error stop '&material: region must be specified'
906   error stop '&material: region out of range'

910   error stop '&material: (atom1 or atom2) and density must not be on same record'
911   error stop '&material: (atom1 or atom2) and xtal must not be on same record'
912   error stop '&material: (atom1 or atom2) and wafer must not be on same record'
913   error stop '&material: (atom1 or atom2) and vsurf must not be on same record'
914   error stop '&material: (atom1 or atom2) and rotx must not be on same record'
915   error stop '&material: (atom1 or atom2) and name must not be on same record'

921   error stop '&material: atom1 and vib must not be on same record'
922   error stop '&material: atom1 and tdebye must not be on same record'
923   error stop '&material: atom2 and syield must not be on same record'

950   error stop '&material: invalid xtal'

960   error stop '&material: name must be specified'
963   error stop '&material: maximum number of atoms (natom) exceeded'

   end subroutine readMat
!
!=======================================================================
!
   subroutine initMat( )
!             =======
!
!     Initialize the mat structure
!
!  Arguments:
!
!     None
!
!  Program History:
!
!      7-Jan-2003  G.Hobler  creation
!
!  Local Variables:
!
      integer(INT4) :: ir, ielemr, ielem, nelemr
      integer(INT4) :: isite, n
      real(REAL8) :: x, u1, prob, wsurf(3), cosx, sinx
      real(REAL8) :: alat(3), alat_scal
      type(CoordSys_type) :: coord_xtal_shared
!
!  Start of Executable Code:
!  ------------------------
!
!  Nothing to do if geometry is defined by cells
!
      if( usecell )  return
!
!  Coordinate system (may be overwritten if a crystalline region exists)
!
      call initCoordSys( coord_xtal_shared )
!          ============         o

      do ir = 1, nr
!
!  Checks that have to be done before assigning default values
!
         if( any( defined( mat(ir)%vib ) .and. defined( mat(ir)%tdebye(:) ) ) )  goto 9503

         if( any( defined( mat(ir)%vsurf(:) ) ) .and.  &
             .not.all( defined( mat(ir)%vsurf(:) ) ) )  goto 956
         if( any( defined( mat(ir)%wafer(:) ) ) .and.  &
             .not.all( defined( mat(ir)%wafer(:) ) ) )  goto 952
         if( defined( mat(ir)%wafer(1) ) .and. .not.defined( mat(ir)%vsurf(1) ) )  goto 953
!
!  Set default values
!
         if( .not.defined( mat(ir)%rotx ) )  mat(ir)%rotx = ZERO

         do ielem = 1, nelem
            if( .not.defined( mat(ir)%syield(ielem) ) )  mat(ir)%syield(ielem) = ZERO

            if( .not.defined( mat(ir)%tdebye(ielem) ) )  &
                  mat(ir)%tdebye(ielem) = elem(ielem)%tdebye

            if( .not.defined( mat(ir)%vib(ielem) ) ) then
               if( mat(ir)%tdebye(ielem) == ZERO ) then
                  mat(ir)%vib(ielem) = 0.083
               else
                  x = mat(ir)%tdebye(ielem) / temp
                  u1 = 12.1 * sqrt( (debye1( x )/x + QUART) /  &
!                                    ======  i
                        (mat(ir)%tdebye(ielem) * elem(ielem)%m) )

                  if( u1 > BIG ) then
                     mat(ir)%vib(ielem) = 0.083
                  else
                     mat(ir)%vib(ielem) = u1
                  end if
               end if
            end if
         end do

         if( .not.defined( mat(ir)%dens ) ) then
            select case( mat(ir)%name )
!               case( 'SIO2', 'SiO2', 'sio2' )
!                  mat(ir)%dens = 7.e22
!               case( 'SI3N4', 'Si3N4', 'si3n4' )
!                  mat(ir)%dens = 9.e22
!               case( 'GAAS', 'GaAs', 'gaas' )
!                  mat(ir)%dens = 4.4e-2      ! = 4.4e22 cm-3
!               case( 'PMMA', 'pmma', 'C5H8O2' )
!                  mat(ir)%name = 'C5H8O2'
!                  mat(ir)%dens = 8.52e-2     ! = 8.52e22 cm-3
!               case( 'WSI2', 'WSi2', 'wsi2' )
!                  mat(ir)%dens = 7.611e-2    ! = 7.611e22 cm-3; H.-H. Vuong, private commun.
               case default
                  nelemr = mat(ir)%nelemr
                  mat(ir)%dens = ONE /  &
                        sum( mat(ir)%prob(1:nelemr) / elem(mat(ir)%ielem(1:nelemr))%dens )
            end select
         end if

         if( .not.all( defined( mat(ir)%vsurf(:) ) ) )  mat(ir)%vsurf(:) = [ ONE, ONE, ZERO ]
         if( .not.all( defined( mat(ir)%wafer(:) ) ) )  mat(ir)%wafer(:) = [ ZERO, ZERO, ONE ]
!
!  Normalize vsurf and wafer (needed for checking orthogonality)
!
         mat(ir)%wafer(:) = mat(ir)%wafer(:) / abs( mat(ir)%wafer(:) )
!                                              ===
         mat(ir)%vsurf(:) = mat(ir)%vsurf(:) / abs( mat(ir)%vsurf(:) )
!                                              ===
!
!  Check validity of the parameters
!
         if( any( mat(ir)%tdebye(:) < ZERO ) )  goto 9624
         if( any( mat(ir)%vib(:) < ZERO ) )  goto 9622
         if( any( mat(ir)%syield(:) > ZERO ) ) then
            if( ion%ranrot )  goto 957
            if( ndim > 1 .or. ldyn .or. usecell )  goto 958
         end if
         if( any( mat(ir)%syield(:) < ZERO ) )  goto 951
         if( mat(ir)%lpoly .and..not. mat(ir)%lxtal )  goto 999
         if( ldyn .and. mat(ir)%lxtal )  goto 1000
         if( abs( mat(ir)%wafer(:).dot.mat(ir)%vsurf(:) ) > 0.001 )  goto 955

         if( mat(ir)%lxtal .and. xtal_struct%uc%name /= ' ' ) then
            do ielemr = 1, mat(ir)%nelemr       ! check if elements of region are compatible
               n = 0                            !    with elements in xtal%uc
               do isite = 1, xtal_struct%uc%nsite
                  if( xtal_struct%uc%ielemr(isite) == ielemr )  n = n + 1
               end do
               prob = float( n ) / float( xtal_struct%uc%nsite )
               if( prob < (ONE-EPSMB)*mat(ir)%prob(ielemr) .or.  &
                   prob > (ONE+EPSMB)*mat(ir)%prob(ielemr) )  goto 970
            end do
         end if
!
!  Calculate derived quantities
!
         mat(ir)%lsput = any( mat(ir)%syield(:) /= ZERO )
         if( ndim > 1 .and. mat(ir)%lsput )  goto 973
                                                ! crystal unit vectors
         if( mat(ir)%rotx /= ZERO ) then
            wsurf(:) = mat(ir)%wafer(:) .cross. mat(ir)%vsurf(:)
            cosx = cos( mat(ir)%rotx*PI/180 )
            sinx = sin( mat(ir)%rotx*PI/180 )
            mat(ir)%vsurf(:) = mat(ir)%vsurf(:)*cosx + wsurf(:)*sinx
         end if

         if( mat(ir)%lxtal ) then
            if( xtal_struct%uc%name == ' ' ) then        ! quarter lattice constant
               mat(ir)%xscale = (EIGHT / mat(ir)%dens) ** THIRD / FOUR
               alat(:) = mat(ir)%xscale
            else                                         ! lattice constants
               alat_scal =  &
                     ( xtal_struct%uc%nsite /  &
                       ( mat(ir)%dens * product( xtal_struct%uc%alat(:) ) ) ) ** THIRD
               alat(:) = alat_scal * xtal_struct%uc%alat(:)
               mat(ir)%xscale = minval(alat(:))
            end if
                                                ! coordinate system
            call initCoordSys( coord_xtal_shared, mat(ir)%wafer, mat(ir)%vsurf, ONE/alat(:) )
!                ============         o                 i              i            i
            call invertCoordSys( coord_xtal_shared )
!                ==============         io
         else
            mat(ir)%xscale = ONE
         end if

      end do

!$omp parallel
      coord_xtal = coord_xtal_shared
!$omp end parallel

      do irxtal = 1, nr                         ! first (and only) crystalline region
         if( mat(irxtal)%lxtal )  exit
      end do
      if( irxtal > nr )  irxtal = 0
!
!  Modify rotation angle (rotx's are added !?!, assume only one crystalline region)
!
      do ir = 1, nr
         if( .not. mat(ir)%lxtal )  cycle
         ion%rotate = ion%rotate - mat(ir)%rotx
      end do
!
!  Checks for lrobin=t
!
      if( geom%lrobin ) then
         if( .not.mat(1)%lxtal )  goto 971
         if( mat(irxtal)%wafer(1)/=ZERO .or. mat(irxtal)%wafer(2)/=ZERO .or.  &
             mat(irxtal)%wafer(3)/=ONE .or.  &
             mat(irxtal)%vsurf(1)/=ONE .or. mat(irxtal)%vsurf(2)/=ZERO .or.  &
             mat(irxtal)%vsurf(3)/=ZERO )  goto 972

      end if
!
!  Initialize ion%irxtal
!
      ion%irxtal = irxtal

      return

9503  error stop '&material: vib and tdebye must not be specified in the same time'
951   error stop '&material: syield>=0 required'
952   error stop '&material: all or no components of wafer must be specified'
953   error stop '&material: vsurf obligatory if wafer is specified'
!954   error stop '&material: wafer not allowed in amorphous region'
955   error stop '&material: vsurf not orthogonal to wafer'
956   error stop '&material: all or no components of vsurf must be specified'
957   error stop '&material: syield not allowed when ranrot=t'
958   error stop '&material: syield requires ndim=1, ldyn=f, usecell=f'

9622  error stop '&material: vib>=0 required'
9624  error stop '&material: tdebye>=0 required'
963   error stop '&material: maximum number of atoms natom exceeded'

970   error stop '&material: elements in region not compatible with elements in crystal'
971   error stop '&material: lxtal=t required if lrobin=t'
972   error stop '&material: wafer=0,0,1 vsurf=1,0,0 required if lrobin=t'
973   error stop '&material: sputtering not allowed for ndim>1'

999   error stop '&material: lpoly=t requires lxtal=t'
1000  error stop '&material: crystalline regions not allowed with ldyn=t'

   end subroutine initMat
!
!=======================================================================
!
   subroutine writeMat( )
!             ========
!
!     Write data derived from &material parameters to the output file
!
!  Arguments:
!
!     None
!
!  Program History:
!
!      7-Jan-2003  G.Hobler  creation
!
!  Local Variables:
!
      integer(INT4) :: ir, i, ielemr, ielem, ielemmin, ielemmax
      real(REAL8) :: a(3), b(3), c(3), dir(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  Open the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i

      if( usecell ) then
         write( 21, '(t2,a)' )  &
               '&material: information read from cell file'
         write( 21, '()' )
         close( 21 )
         return
      end if

      do ir = 1, nr

         if( ir == 1 ) then
            write( 21, '(t2,a)', advance='no' ) '&material: '
         else
            write( 21, '(t2,a)', advance='no' ) '           '
         end if
         write( 21, '(a,i2,a)', advance='no' )  &
               'Region ' // str( ir ) // ' (' // str( mat(ir)%name ) // '): '

         if( mat(ir)%lxtal .and. mat(ir)%lpoly ) then
            write( 21, '(a)', advance='no' )  &
                  'polycrystalline, '
         else if( mat(ir)%lxtal ) then
            write( 21, '(a)', advance='no' )  &
                  'crystalline, '
         else
            write( 21, '(a)', advance='no' )  &
                  'amorphous, '
         end if
         write( 21, '(a)' )  &
               'density=' // str( 1e24*mat(ir)%dens, 'es8.2e2' )

         if( mat(ir)%lxtal ) then

            if( xtal_struct%uc%name == ' ' ) then
               write( 21, '(t13,a)', advance='no' )  &
                     '   crystal system=default diamond, '
            else
               write( 21, '(t13,a)', advance='no' )  &
                     '   crystal system=' // str( xtal_struct%uc%name ) // ', '
            end if
            write( 21, '(a)' )  &
                  'xscale=' // str( mat(ir)%xscale, 'f99.3' )

            write( 21, '(t13,a)', advance='no' )  &
                  '   Lattice vector a: '
            dir = [ 1 , 0 , 0 ]
            call fromCoordSys( coord_xtal, dir, a)
!                ============         i        i   o
            write( 21, '(a)' )  &
                  '(' // str( a(1), 'f99.3' ) //  &
                  ',' // str( a(2), 'f99.3' ) //  &
                  ',' // str( a(3), 'f99.3' ) // ')'

            write( 21, '(t13,a)', advance='no' )  &
                  '   Lattice vector b: '
            dir = [ 0 , 1 , 0 ]
            call fromCoordSys( coord_xtal, dir, b)
!                ============         i        i   o
            write( 21, '(a)' )  &
                  '(' // str( b(1), 'f99.3' ) //  &
                  ',' // str( b(2), 'f99.3' ) //  &
                  ',' // str( b(3), 'f99.3' ) // ')'

            write( 21, '(t13,a)', advance='no' )  &
                  '   Lattice vector c: '
            dir = [ 0 , 0 , 1 ]
            call fromCoordSys( coord_xtal, dir, c)
!                ============         i        i   o
            write( 21, '(a)' )  &
                  '(' // str( c(1), 'f99.3' ) //  &
                  ',' // str( c(2), 'f99.3' ) //  &
                  ',' // str( c(3), 'f99.3' ) // ')'

            do ielemr = 1, mat(ir)%nelemr
               if( ielemr == 1 ) then
                  write( 21, '(t13,a)', advance='no' )  &
                        '   Thermal vibration amplitude='
               else
                  write( 21, '(t13,a)', advance='no' )  &
                        '                               '
               end if
               write( 21, '(a)' )  &
                     str( mat(ir)%vib(mat(ir)%ielem(ielemr)) , 'f99.4' ) //  &
                     ' (' // str( elem(mat(ir)%ielem(ielemr))%name ) // ')'
            end do
         end if

         if( any( mat(ir)%syield > ZERO ) ) then
            write( 21, '(t13,a,(t28,5(f7.4:'','' )))' ) 'sputter yield=', mat(ir)%syield(:)
         end if

      end do

      write( 21, '(/t13,a)' ) 'Atoms in the regions:'

      do i = 0, (nelem-1)/8

         ielemmin = 8*i + 1
         ielemmax = min( 8*i + 8, nelem )

         write( 21, '(t13,a)', advance='no' ) '+--------+'
         write( 21, '(8a)' ) &
               ('------+', ielem=ielemmin,ielemmax)

         write( 21, '(t13,a, t22,a)', advance='no' ) '|', '|'
         write( 21, '(8(a,i2,a))' ) &
               ('Atom', ielem, '|', ielem=ielemmin,ielemmax)

         write( 21, '(t13,a, t22,a)', advance='no' ) '|', '|'
         write( 21, '(8(a,a2,a))' ) &
               (' (', elem(ielem)%name, ') |', ielem=ielemmin,ielemmax)

         write( 21, '(t13,a)', advance='no' ) '+--------+'
         write( 21, '(8a)' ) &
               ('------+', ielem=ielemmin,ielemmax)

         do ir = 1, nr
            write( 21, '(t13,a,i2,a)', advance='no' ) '|Region', ir, '|'
            do ielem = ielemmin, ielemmax
               do ielemr = 1, mat(ir)%nelemr
                  if( ielem == mat(ir)%ielem(ielemr) )  exit
               end do
               if( ielemr <= mat(ir)%nelemr ) then
                  write( 21, '(f5.1,a)', advance='no' ) &
                     100*mat(ir)%prob(ielemr), '%|'
               else
                  write( 21, '(f5.1,a)', advance='no' ) &
                     ZERO, '%|'
               end if
            end do
            write( 21, '()' )
         end do

         write( 21, '(t13,a)', advance='no' ) '+--------+'
         write( 21, '(8a)' ) &
               ('------+', ielem=ielemmin,ielemmax)
      end do

      write( 21, '()' )

      close( 21 )

   end subroutine writeMat
!
!=======================================================================
!
end module Param_Mat_mod

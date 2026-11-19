module Param_Geom_mod
!
!  Holds data about the geometry of the target
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, ONE, TWO, SIX, TEN, HALF, EPSMB, SMALL,  &
                                    BIG, NEGBIG,  &
                                    ENUM
   use Array_mod,             only: alloc, check_allocation
   use Open_mod,              only: openFile
   use Stack_mod,             only: LogicalStack, Int4Stack, Real8Stack,  &
                                    newStack, sizeOfStack, popStack, pushStack, sizeOfStack,  &
                                    stack2Array
   use String_mod,            only: str, to_lower
   use Undefine_mod,          only: undefine, defined, definedlength

   use Geom_Polygon_mod,      only: Polygon_type, initPolygon

   use Param_Setup_mod,       only: ndim, nr, ldyn, usecell
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
   public :: readGeom               ! read &geometry record from input file
!            ========
   public :: initGeom               ! initialize the geom structure
!            ========
   public :: writeGeom              ! write data derived from &geometry
!            =========              !    record to the output file
!
!  Private Procedures:
!  ==================
!
   private :: updateRegionPolygon   ! update region polygon taking boundary types into account
!             ===================
!
!  Public Data:
!  ===========
!
   type :: Geom_type                         ! geometry data type
!          =========
!
!  Parameters defining the geometry coordinate system
!
      logical :: lcart                       ! flag indicating Cartesian geometry coordinate
                                             !    system (default)
      logical :: lcylx                       ! flag indicating cylindrical geometry coordinate
                                             !    system with axis parallel x
      logical :: lcyly                       ! flag indicating cylindrical geometry coordinate
                                             !    system with axis parallel y
      logical :: lcylz                       ! flag indicating cylindrical geometry coordinate
                                             !    system with axis parallel z
      logical :: lrobin                      ! flag indicating special geometry for 1994
                                             !    round robin simulation
      logical :: lspher                      ! flag indicating spherical geometry coordinate
                                             !    system
      real(real8) :: xcen                    ! x-coordinate of geometry coordinate system center
      real(real8) :: ycen                    ! y-coordinate of geometry coordinate system center
      real(real8) :: zcen                    ! z-coordinate of geometry coordinate system center
!
!  Other parameters
!
      logical :: freex(2)                    ! flag indicating whether left/right surface
                                             !    is free (ldyn=t, ndim=2, old mode)
      logical :: freez(2)                    ! flag indicating whether top/bottom surface
                                             !    is free (ldyn=t, ndim=2, old mode)
      logical :: lcell                       ! cell-based geometry definition flag
      logical :: usehyperoctree              ! whether or not to use the hyperoctree variant of
                                             !    the sdf variants (ndim=2,3)
      logical :: usesdf                      ! flag indicating sdffiles shall be used rather
                                             !    than polygons (ndim=2,3)

      character(80) :: geomfile              ! Name of file containing further &GEOMETRY records
      character(80), allocatable :: sdffile(:) ! (0:nr) Name of an signed distance function file

      integer(INT4) :: bndtypex(2)           ! boundary type for left and right boundary
      integer(INT4) :: bndtypey(2)           ! boundary type for front and back boundary
      integer(INT4) :: bndtypez(2)           ! boundary type for upper and lower boundary
                                             !    (values consistent with ENUM%bndtype)
                                             !    (finalized in initCell*d if geom%lcell=t)
      integer(INT4) :: grid                  ! grid type (1..Cartesian, 2..curvilinear)
      integer(INT4) :: nbox                  ! number of cell boxes in z direction
      integer(INT4) :: nbox2(2)              ! number of mesh or cell boxes in x and z direction

      real(REAL8) :: fracupd                 ! minimum fractional cell width change triggering
                                             !    mesh update (ldyn=t, ndim=1)
      real(REAL8) :: posfix                  ! position which remains fixed during
                                             !    relaxation (ldyn=t and ndim=1)
      real(REAL8), allocatable :: posif(:)   ! (0:nr) interface positions (ndim=1)
      real(REAL8) :: wbox                    ! width of cell boxes in z direction
      real(REAL8) :: wbox2(2)                ! width of mesh boxes in x and z direction
      real(REAL8) :: wrough                  ! thickness of surface roughness layer
      real(REAL8) :: xlim(3)                 ! *lim(1:2): limits of x/y/z coordinates in
      real(REAL8) :: ylim(3)                 !            geometry domain except at reflective
      real(REAL8) :: zlim(3)                 !            boundary where the reflected
                                             !            geometry is included
                                             ! *lim(3): position of reflective boundary if any
                                             !    (finalized in initSurface*)
      real(REAL8) :: xper(2)                 ! x/y/z coordinates of the period of periodic
      real(REAL8) :: yper(2)                 !    structures. Relevant only if geometry is
      real(REAL8) :: zper(2)                 !    defined by positions or polygons and if
                                             !    periodic boundaries are defined
                                             !    Obsolete (kept until eliminated everywhere)

      type(Polygon_type), allocatable :: region(:) ! (0:nr) region polygon (ndim=2)

   end type Geom_type

   type(Geom_type), public :: geom           ! geometry structure
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine readGeom( )
!             ========
!
!     Read the &geometry record from the input file INP and initialize the geom structure.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      logical :: freex(2)
      logical :: freez(2)
      logical :: lcylx
      logical :: lcyly
      logical :: lcylz
      logical :: lperx                          ! obsolete
      logical :: lrobin
      logical :: lspher
#ifdef GTS
      logical :: usehyperoctree
#endif

      character(80) :: bndtypex(2)
      character(80) :: bndtypey(2)
      character(80) :: bndtypez(2)
      character(80) :: filgeom                  ! obsolete
      character(80) :: geomfile
      character(80) :: grid
      character(80) :: sdffile

      integer(INT4) :: nbox
      integer(INT4) :: nbox2(2)
      integer(INT4) :: nionupd                  ! obsolete
      integer(INT4) :: point
      integer(INT4) :: points(21)
      integer(INT4) :: region

      real(REAL8) :: center(3)
      real(REAL8) :: fracupd
      real(REAL8) :: pos(2)
      real(REAL8) :: posfix
      real(REAL8) :: posif(20)
      real(REAL8) :: wbox
      real(REAL8) :: wbox2(2)
      real(REAL8) :: wrough
      real(REAL8) :: xlim(2)
      !real(REAL8) :: ylim(2)                   ! useless without 3D region-based geometry
      real(REAL8) :: zlim(2)
      real(REAL8) :: xper(2)                    ! obsolete
      real(REAL8) :: yper(2)                    ! obsolete
      real(REAL8) :: xscale                     ! scaling factor for x coordinates
      real(REAL8) :: zscale                     ! scaling factor for z coordinates

      namelist /geometry/ freex, freez, lcylx, lcyly, lcylz, lrobin, lspher,  &
                          bndtypex, bndtypey, bndtypez, geomfile, grid, sdffile,  &
                          nbox, nbox2, point, points, region,  &
                          center, fracupd, pos, posfix, posif, wbox, wbox2, wrough,  &
                          xlim, xscale, zscale  !, ylim
#ifdef GTS
      namelist /geometry/ usehyperoctree
#endif
      ! obsolete
      namelist /geometry/ lperx, filgeom, nionupd, xper, yper

      logical :: isEmpty
      integer(INT4) :: np, ip, jp
      type(Int4Stack) :: ipStack
      type(Real8Stack) :: xStack, zStack

      character(256) :: msg
      integer(INT4) :: ir, iir, irmin, jr, iip, npm, npold, ip0, i, incr, ios, ninp
      real(REAL8) :: eps, dummy

      integer(INT4), allocatable :: ipg(:,:), itmp(:,:)       ! (max(npr),0:nr)
      integer(INT4), allocatable :: neighr(:,:), neighp(:,:)  ! (max(npr),0:nr)
      real(REAL8), allocatable :: px(:), pz(:)
      type(Int4Stack), allocatable :: ipgStack(:)
!
!  Start of Executable Code:
!  ------------------------
!
      call alloc( geom%posif, 0, nr, 'geom%posif', 'readGeom' )
!          =====      o       i  i        i            i
      call alloc( geom%sdffile, 0, nr, 'geom%sdffile', 'readGeom' )
!          =====      o         i  i          i            i
!
!  Set default values of structure variables
!
      geom%freex = [ .false., .false. ]
      geom%freez = [ .true., .false. ]
      geom%lcylx = .false.
      geom%lcyly = .false.
      geom%lcylz = .false.
      geom%lrobin = .false.
      geom%lspher = .false.
!#ifdef GTS
      geom%usehyperoctree = .false.
!#else
      geom%usehyperoctree = .false.
!#endif

      call undefine( geom%geomfile )
      call undefine( geom%sdffile(:) )

      call undefine( geom%bndtypex(:) )
      call undefine( geom%bndtypey(:) )
      call undefine( geom%bndtypez(:) )
      call undefine( geom%grid )

      call undefine( geom%nbox )
      call undefine( geom%nbox2(:) )

      call undefine( geom%fracupd )
      call undefine( geom%posfix )
      call undefine( geom%posif(:) )
      call undefine( geom%wbox )
      call undefine( geom%wbox2(:) )
      call undefine( geom%wrough )
      call undefine( geom%xlim(:) )
      call undefine( geom%ylim(:) )
      call undefine( geom%zlim(:) )
      call undefine( geom%xcen )
      call undefine( geom%ycen )
      call undefine( geom%zcen )
!
!  First loop over reading records from input file
!  -----------------------------------------------
!
      call openFile( 11, 'INP', 'formatted', 'read' )
!          ========  i     i         i         i
!
!  Set default values of xscale, zscale (they have no geom% counterparts)
!
      ir = -1
      xscale = ONE
      zscale = ONE

      do                                        ! loop over all records
!
!  Set default values of namelist variables
!
         freex = geom%freex
         freez = geom%freez
         lcylx = geom%lcylx
         lcyly = geom%lcyly
         lcylz = geom%lcylz
         lperx = .false.                        ! obsolete
         lrobin = geom%lrobin
         lspher = geom%lspher
#ifdef GTS
         usehyperoctree = geom%usehyperoctree
#endif

         call undefine( bndtypex(:) )
         call undefine( bndtypey(:) )
         call undefine( bndtypez(:) )
         call undefine( filgeom )               ! obsolete
         call undefine( geomfile )
         call undefine( grid )
         call undefine( sdffile )

         call undefine( nbox )
         call undefine( nbox2(:) )
         call undefine( nionupd )               ! obsolete
         call undefine( point )
         call undefine( points(:) )
         call undefine( region )

         call undefine( center )
         call undefine( fracupd )
         call undefine( pos(:) )
         call undefine( posfix )
         call undefine( posif(:) )
         call undefine( wbox )
         call undefine( wbox2(:) )
         call undefine( wrough )
         call undefine( xlim(:) )
         !call undefine( ylim(:) )
         call undefine( zlim(:) )
         call undefine( xper(:) )               ! obsolete
         call undefine( yper(:) )               ! obsolete
!
!  Read &geometry record from the input file
!
         read( 11, nml=geometry, iostat=ios, iomsg=msg )
         if( ios < 0 )  exit
         if( ios > 0 ) then
            write( *, '(/a)' )  trim( msg )
            goto 901
         end if
!
!  Check obsolete parameters
!
         if( defined( filgeom ) )  goto 909
         if( defined( nionupd ) )  goto 980
         if( lperx )  goto 981
         if( any( defined( xper(:) ) ) )  goto 982
         if( any( defined( yper(:) ) ) )  goto 983
!
!  Issue warnings for parameters that will be ignored
!  TODO: move to initGeom
!
         if( .not.lcylx .and..not.lcyly .and..not.lcylz ) then
            if( any( defined( center ) ) )  &
                 write( *, * )  '(WARNING) &geometry: center ignored unless lcyl=t'
         end if
         if( .not.ldyn .and. ( defined( fracupd ) .or. defined( posfix ) ) )  &
              print *,  '(WARNING) &geometry: fracupd, posfix ignored when ldyn=f'
         if( usecell ) then
            if( any( defined( posif(:) ) ) )  &
                 print *,  '(WARNING) &geometry: posif ignored when usecell=t'
            if( defined( point ) .or. any( defined( pos(:) ) ) .or. defined( region ) .or.  &
                any( defined( points ) ) ) then
                  write( *, * )  '(WARNING) &geometry: point, pos, region, points ignored'
                  write( *, * )  '                     if usecell=t'
            end if
         end if
!
!  Check consistency with ndim
!
         if( ndim /= 2 .and. defined( point ) )  goto 970
         if( ndim == 1 .and. defined( region ) )  goto 971
         if( ndim /= 1 ) then
            if( lrobin )  goto 972
            if( lspher )  goto 973
            if( defined( fracupd ) )  goto 974
            if( defined( posfix ) )  goto 975
            if( any( defined( posif ) ) )  goto 976
         end if
!
!  Check index parameters
!
         if( defined( point ) ) then
            if( point < 1 )  goto 910
         end if

         if( defined( region ) ) then
            if( region < 0 .or. region > nr )  goto 920
         end if
!
!  Check parameters on consistency with index parameters
!
         if( defined( point ) ) then
            if( .not.all( defined( pos(:) ) ) )  goto 911
            if( defined( region ) )  goto 917
            if( any( defined( points ) ) )  goto 918
            if( defined( nbox ) )  goto 9190
            if( any( defined( nbox2(:) ) ) )  goto 9191
            if( defined( wbox ) )  goto 9192
            if( any( defined( wbox2(:) ) ) )  goto 9193
            if( defined( wrough ) )  goto 954
            if( any( defined( xlim(:) ) ) )  goto 952
            !if( any( defined( ylim(:) ) ) )  goto 953
            if( any( defined( zlim(:) ) ) )  goto 9530
            if( defined( geomfile ) )  goto 956
            if( defined( sdffile ) )  goto 959
         else if( defined( region ) ) then
            if( .not.any( defined( points ) ) .and. .not.defined( sdffile ) )  goto 923
            if( any( defined( pos(:) ) ) )  goto 926
            if( defined( nbox ) )  goto 9270
            if( any( defined( nbox2(:) ) ) )  goto 9271
            if( defined( wbox ) )  goto 9272
            if( any( defined( wbox2(:) ) ) )  goto 9273
            if( defined( wrough ) )  goto 954
            if( any( defined( xlim(:) ) ) )  goto 952
            !if( any( defined( ylim(:) ) ) )  goto 953
            if( any( defined( zlim(:) ) ) )  goto 9530
            if( defined( geomfile ) )  goto 955
         end if

         if( defined( sdffile ) ) then
            if( .not.defined( region ) )  goto 922
         end if
!
!  Copy parameters to geom structure
!
         geom%freex = freex
         geom%freez = freez
         geom%lcylx = lcylx
         geom%lcyly = lcyly
         geom%lcylz = lcylz
         geom%lrobin = lrobin
         geom%lspher = lspher
#ifdef GTS
         geom%usehyperoctree = usehyperoctree
#endif

         do i = 1, 2
            if( defined( bndtypex(i) ) ) then
               if( to_lower( bndtypex(i) ) == 'free' ) then
!                  ========
                  geom%bndtypex(i) = ENUM%bndtype%free
               else if( to_lower( bndtypex(i) ) == 'extending' ) then
!                       ========
                  geom%bndtypex(i) = ENUM%bndtype%extend
               else if( to_lower( bndtypex(i) ) == 'periodic' ) then
!                       ========
                  geom%bndtypex(i) = ENUM%bndtype%periodic
               else if( to_lower( bndtypex(i) ) == 'reflective' ) then
!                       ========
                  geom%bndtypex(i) = ENUM%bndtype%reflect
               else
                  goto 956
               end if
            end if
         end do
         do i = 1, 2
            if( defined( bndtypey(i) ) ) then
               if( to_lower( bndtypey(i) ) == 'free' ) then
!                  ========
                  geom%bndtypey(i) = ENUM%bndtype%free
               else if( to_lower( bndtypey(i) ) == 'extending' ) then
!                       ========
                  geom%bndtypey(i) = ENUM%bndtype%extend
               else if( to_lower( bndtypey(i) ) == 'periodic' ) then
!                       ========
                  geom%bndtypey(i) = ENUM%bndtype%periodic
               else if( to_lower( bndtypey(i) ) == 'reflective' ) then
!                       ========
                  geom%bndtypey(i) = ENUM%bndtype%reflect
               else
                  goto 9561
               end if
            end if
         end do
         do i = 1, 2
            if( defined( bndtypez(i) ) ) then
               if( to_lower( bndtypez(i) ) == 'free' ) then
!                  ========
                  geom%bndtypez(i) = ENUM%bndtype%free
               else if( to_lower( bndtypez(i) ) == 'extending' ) then
!                       ========
                  geom%bndtypez(i) = ENUM%bndtype%extend
               else if( to_lower( bndtypez(i) ) == 'periodic' ) then
!                       ========
                  geom%bndtypez(i) = ENUM%bndtype%periodic
               else if( to_lower( bndtypez(i) ) == 'reflective' ) then
!                       ========
                  geom%bndtypez(i) = ENUM%bndtype%reflect
               else
                  goto 957
               end if
            end if
         end do
         if( defined( grid ) ) then
            if( to_lower( grid ) == 'cartesian' .or. to_lower( grid ) == 'regular' ) then
!               ========                             ========
               geom%grid = ENUM%grid%regular
            else if( to_lower( grid ) == 'rectilinear' ) then
!                    ========
               if( ndim /= 1 )  goto 9581
               geom%grid = ENUM%grid%rectilinear
            else if( to_lower( grid ) == 'curvilinear' ) then
!                    ========
               if( ndim /= 2 )  goto 9582
               geom%grid = ENUM%grid%curvilinear
            else
               goto 958
            end if
         end if

         if( defined( geomfile ) )  geom%geomfile = geomfile
         if( defined( sdffile ) )  geom%sdffile(region) = sdffile

         if( defined( nbox ) )  geom%nbox = nbox
         if( defined( nbox2(1) ) )  geom%nbox2(1) = nbox2(1)
         if( defined( nbox2(2) ) )  geom%nbox2(2) = nbox2(2)

         if( defined( center(1) ) )  geom%xcen = center(1)
         if( defined( center(2) ) )  geom%ycen = center(2)
         if( defined( center(3) ) )  geom%zcen = center(3)
         if( defined( fracupd ) )  geom%fracupd = fracupd
         if( defined( posfix ) )  geom%posfix = posfix
         if( defined( wbox2(1) ) )  geom%wbox = wbox
         if( defined( wbox2(1) ) )  geom%wbox2(1) = wbox2(1)
         if( defined( wbox2(2) ) )  geom%wbox2(2) = wbox2(2)
         if( defined( wrough ) )  geom%wrough = wrough
         if( defined( xlim(1) ) )  geom%xlim(1) = xlim(1)
         if( defined( xlim(2) ) )  geom%xlim(2) = xlim(2)
         !if( defined( ylim(1) ) )  geom%ylim(1) = ylim(1)
         !if( defined( ylim(2) ) )  geom%ylim(2) = ylim(2)
         if( defined( zlim(1) ) )  geom%zlim(1) = zlim(1)
         if( defined( zlim(2) ) )  geom%zlim(2) = zlim(2)

         if( .not.usecell ) then                ! Append positions to geom%posif
            do i = 1, definedlength( posif )
               ir = ir + 1
               if( ir > nr )  goto 906          ! too many posif's
               geom%posif(ir) = posif(i)
            end do
         end if

      end do
!
!  End of first loop over reading records from the input file
!  ----------------------------------------------------------
!
!  Check if enough posif values have been specified
!
      if( ndim == 1 .and. nr > 1 .and. .not.usecell ) then
         if( ir < nr )  goto 905
      end if
!
!  Check if no or all sdffiles are defined
!
      if( all( defined( geom%sdffile(:) ) ) ) then
         geom%usesdf = .true.
      else if( .not.any( defined( geom%sdffile(:) ) ) ) then
         geom%usesdf = .false.
      else
         goto 960
      end if
!
!  In the remainder of this subroutine the polygons and the internal mesh are defined
!  We can skip this if ndim/=2, usecell=t, or usesdf=f
!
      if( ndim /= 2 .or. usecell .or. geom%usesdf ) then
         close( 11 )
         return
      end if
!
!  Open and check integrity of geometry file if specified
!
      if( defined( geom%geomfile ) ) then

         call openFile( 111, trim(geom%geomfile), 'formatted', 'read' )
!             ========   i            i                i         i

         do
            read( 111, nml=geometry, iostat=ios, iomsg=msg )
            if( ios < 0 )  exit
            if( ios > 0 ) then
               write( *, '(/a)' )  trim( msg )
               goto 901
            end if
         end do

      end if
!
!  Read points and put them on stacks
!
      call newStack( ipStack )
!          ========     o
      call newStack( xStack )
!          ========    o
      call newStack( zStack )
!          ========    o
!
!  Second loop over reading records from input and geometry file: Read the point definitions
!  -----------------------------------------------------------------------------------------
!
      rewind 11
      if( defined( geom%geomfile ) ) then
         ninp = 111                             ! start reading geomfile so any point redefined
      else                                      ! in INP file overwrites point in geomfile
         ninp = 11
      end if
      rewind ninp

      np = 0

      do

         call undefine( point )
!
!  Read &geometry record from input or geometry file
!
         read( ninp, nml=geometry, iostat=ios, iomsg=msg )
         if( ios < 0 ) then
            if( ninp == 111 ) then
               ninp = 11
               rewind 11
               cycle
            else
               exit
            end if
         end if
         if( ios > 0 ) then
            write( *, '(/a)' )  trim( msg )
            goto 901
         end if

         if( .not.defined( point ) )  cycle

         np = max( np, point )
         call pushStack( point, ipStack )
!             =========    i      io
         call pushStack( pos(1)*xscale, xStack )
!             =========        i          io
         call pushStack( pos(2)*zscale, zStack )
!             =========        i          io
      end do
!
!  End of second loop over reading records
!  ---------------------------------------
!
!  Transfer points from stacks to arrays
!
      allocate( px( np ) )
      allocate( pz( np ) )
      call undefine( px(:) )

      do
         call popStack( ipStack, ip, isEmpty )
!             ========    io     o      o
         if( isEmpty )  exit
         if( ip > np )  goto 915
         if( defined( px(ip) ) ) then
            print *, '(WARNING): &geometry: duplicate point definition. ' // &
                     'Earlier definition overridden.'
            call popStack( xStack, dummy, isEmpty )
!                ========    io      o       o
            call popStack( zStack, dummy, isEmpty )
!                ========    io      o       o
            cycle
         end if
         call popStack( xStack, px(ip), isEmpty )
!             ========    io      o        o
         call popStack( zStack, pz(ip), isEmpty )
!             ========    io      o        o
      end do
!
!  Make sure coordinates intended to be equal are indeed equal
!  (numbers read from input file are converted from default real to real(REAL8))
!
      eps = TEN**( -precision( 1. ) )
      do ip = 1, np-1
         do jp = ip+1, np
            if( abs( px(ip) - px(jp) ) < eps * max( abs( px(ip) ), abs( px(jp) ) ) )  &
                  px(jp) = px(ip)
            if( abs( pz(ip) - pz(jp) ) < eps * max( abs( pz(ip) ), abs( pz(jp) ) ) )  &
                  pz(jp) = pz(ip)
         end do
      end do
!
!  Read the definition of the regions and put them on stack
!
      allocate( ipgStack( 0:nr ) )
      do ir = 0, nr
         call newStack( ipgStack(ir) )
!             ========       o
      end do
!
!  Third loop over reading records from input and geometry file: Read the region definitions
!  -----------------------------------------------------------------------------------------
!
      rewind 11
      if( defined( geom%geomfile ) ) then
         ninp = 111                             ! start reading geomfile so any point redefined
      else                                      ! in INP file overwrites point in geomfile
         ninp = 11
      end if
      rewind ninp

      do

         call undefine( region )
         call undefine( points(:) )
!
!  Read &geometry record from input or geometry file
!
         read( ninp, nml=geometry, iostat=ios, iomsg=msg )
         if( ios < 0 ) then
            if( ninp == 111 ) then
               ninp = 11
               rewind 11
               cycle
            else
               exit
            end if
         end if
         if( ios > 0 ) then
            write( *, '(/a)' )  trim( msg )
            goto 901
         end if

         if( .not.defined( region ) )  cycle
         ir = region
!
!  Extract point indices from points and put them on stack
!
         if( defined( points(size( points )) ) )  goto 941
         do i = 1, definedlength( points(:) )
            if( points(i) < 1 .or. points(i) > np )  goto 925
            call pushStack( points(i), ipgStack(ir) )
!                =========      i           io
         end do

      end do
!
!  End of third loop over reading records from input and geometry file
!  -------------------------------------------------------------------
!
      close( 11 )
      if( defined( geom%geomfile ) )  close( 111 )
!
!  Transfer points from stack to arrays and initialize global point index arrays
!
      allocate( geom%region( 0:nr ) )

      do ir = 0, nr
         geom%region(ir)%np = sizeOfStack( ipgStack(ir) )
!                             ===========       i
         if( ir > 0 .and. geom%region(ir)%np == 0 )  goto 932
      end do

      do ir = 0, nr
         np = geom%region(ir)%np
         allocate( geom%region(ir)%px( np+1 ) )
         allocate( geom%region(ir)%pz( np+1 ) )
      end do

      npm = maxval(geom%region(:)%np)
      allocate( ipg( npm+1, 0:nr ) )

      do ir = 0, nr
         np = geom%region(ir)%np
         do ip = np, 1, -1
            call popStack( ipgStack(ir), ipg(ip,ir), isEmpty )
!                ========      io            o          o
         end do
         geom%region(ir)%px(1:np) = px(ipg(1:np,ir))
         geom%region(ir)%pz(1:np) = pz(ipg(1:np,ir))
                                                ! close the polygon
         geom%region(ir)%px(np+1) = geom%region(ir)%px(1)
         geom%region(ir)%pz(np+1) = geom%region(ir)%pz(1)
         ipg(np+1,ir) = ipg(1,ir)
      end do
!
!  Setup neighbor list
!
!  Segment (ip,ip+1) of region ir is the same as
!  segment (neighp(ip,ir),neighp(ip+1,ir)) of region neighr(ip,ir).
!
      if( geom%region(0)%np == 0 ) then
         irmin = 1
      else
         irmin = 0
      end if

      allocate( neighr( npm, irmin:nr ) )
      allocate( neighp( npm, irmin:nr ) )

      neighr(:,:) = -1
      neighp(:,:) = -1

      do ir = 1, nr                                ! loop over untreated points
         do ip = 1, geom%region(ir)%np
            iir_loop: do iir = irmin, ir-1         ! loop over already treated points
               do iip = 1, geom%region(iir)%np
                  if( ( ipg(ip,ir) == ipg(iip,iir) .and.  &
                        ipg(ip+1,ir) == ipg(iip+1,iir) ) .or.  &
                      ( ipg(ip,ir) == ipg(iip+1,iir) .and.  &
                        ipg(ip+1,ir) == ipg(iip,iir) ) ) then ! found=.true.
                     neighr(iip,iir) = ir
                     neighr(ip,ir) = iir
                     neighp(iip,iir) = ip
                     neighp(ip,ir) = iip
                     exit iir_loop
                  end if
               end do
            end do iir_loop
         end do
      end do

      if( geom%region(0)%np == 0 ) then
!
!  Determine region 0 (circumference)
!
         ir_loop: do ir = 1, nr                 ! Find first point on circumference
            do ip = 1, geom%region(ir)%np
               if( neighr(ip,ir) == -1 )  exit ir_loop
            end do
         end do ir_loop
         if( ir == nr+1 )  goto 933
         incr = +1
         ip0 = 1
         neighr(ip,ir) = 0
         neighp(ip,ir) = ip0
         call pushStack( ipg(ip,ir), ipgStack(0) )
!             =========      i           io

         do
            ip = ip + incr
            ip = modulo( ip-1, geom%region(ir)%np ) + 1
            if( neighr(ip,ir) == 0 ) then       ! circumference closed
               exit
            else if( neighr(ip,ir) == -1 ) then ! new segment on circumference
               ip0 = ip0 + 1
               neighr(ip,ir) = 0
               neighp(ip,ir) = ip0
               if( incr == +1 ) then
                  call pushStack( ipg(ip,ir), ipgStack(0) )
!                      =========      i           io
               else
                  call pushStack( ipg(ip+1,ir), ipgStack(0) )
!                      =========       i            io
               end if
            else                                ! segment goes into interior
               if( ipg(ip,ir) == ipg(neighp(ip,ir),neighr(ip,ir)) ) then
                  jr = neighr(ip,ir)
                  ip = neighp(ip,ir)
                  ir = jr
                  incr = -incr
               else if( ipg(ip,ir) == ipg(neighp(ip,ir)+1,neighr(ip,ir)) ) then
                  jr = neighr(ip,ir)
                  ip = neighp(ip,ir)
                  ir = jr
               else
                  goto 934
               end if
            end if
         end do
                                                ! convert stacks to arrays
         geom%region(0)%np = sizeOfStack( ipgStack(0) )
!                            ===========      i
         np = geom%region(0)%np
         npold = size( ipg, dim=1 )
         if( np+1 > npold ) then
            allocate( itmp( np+1, 0:nr ) )
            itmp(:,:) = -1
            itmp(1:npold,1:nr) = ipg(:,1:nr)
            call move_alloc( itmp, ipg )
         end if
         deallocate( geom%region(0)%px, geom%region(0)%pz )
         allocate( geom%region(0)%px( np+1 ) )
         allocate( geom%region(0)%pz( np+1 ) )

         do ip = np, 1, -1
            call popStack( ipgStack(0), ipg(ip,0), isEmpty )
!                ========      io          o         o
         end do
         geom%region(0)%px(1:np) = px(ipg(1:np,0))
         geom%region(0)%pz(1:np) = pz(ipg(1:np,0))
                                                ! close the polygon
         geom%region(0)%px(np+1) = geom%region(0)%px(1)
         geom%region(0)%pz(np+1) = geom%region(0)%pz(1)
         ipg(np+1,0) = ipg(1,0)

      end if
!
!  Check validity of geometry definition
!
      do ir = irmin, nr
         np = geom%region(ir)%np
         if( any( neighr(1:np,ir) < 0 ) )  goto 945
         if( any( neighr(1:np,ir) > nr ) )  goto 945
      end do

      return

901   error stop '&geometry: error during reading namelist'
905   error stop '&geometry: too few posif values'
906   error stop '&geometry: too many posif values'

909   error stop 'geometry: filgeom is obsolete. Use geomfile instead'

910   error stop '&geometry: point>0 required'
911   error stop '&geometry: pos must be specified on point record'
915   error stop '&geometry: point indizes must be contiguous'
917   error stop '&geometry: point and region must not be on same record'
918   error stop '&geometry: point and points must not be on same record'
9190  error stop '&geometry: point and nbox must not be on same record'
9191  error stop '&geometry: point and nbox2 must not be on same record'
9192  error stop '&geometry: point and wbox must not be on same record'
9193  error stop '&geometry: point and wbox2 must not be on same record'

920   error stop '&geometry: 0<=region<nr required'
922   error stop '&geometry: sdffile requires region to be specified'
923   error stop '&geometry: points or sdffile must be specified on region record'
924   error stop '&geometry: error during reading points'
925   error stop '&geometry: illegal point index referenced in region record'
926   error stop '&geometry: pos and region must not be on same record'
9270  error stop '&geometry: region and nbox must not be on same record'
9271  error stop '&geometry: region and nbox2 must not be on same record'
9272  error stop '&geometry: region and wbox must not be on same record'
9273  error stop '&geometry: region and wbox2 must not be on same record'

932   error stop '&geometry: undefined region'
933   error stop '&geometry: could not find any point of the circumference'
934   error stop 'readGeom: internal error in finding the circumference'

941   error stop '&geometry: Too many values of points'
944   error stop '&geometry: points variable erroneous'
945   error stop '&geometry: illegal geometry (holes?)'

952   error stop '&geometry: xlim and region or point must not be on same record'
!953   error stop '&geometry: ylim and region or point must not be on same record'
9530  error stop '&geometry: zlim and region or point must not be on same record'
954   error stop '&geometry: wrough and region or point must not be on same record'
955   error stop '&geometry: geomfile and region or point must not be on same record'
956   error stop '&geometry: bndtypex must be ''free'', ''extending'', ''periodic'', ' //  &
                 'or ''reflective'''
9561  error stop '&geometry: bndtypey must be ''free'', ''extending'', ''periodic'', ' //  &
                 'or ''reflective'''
957   error stop '&geometry: bndtypez must be ''free'', ''extending'', ''periodic'', ' //  &
                 'or ''reflective'''
958   error stop "&geometry: grid must be 'Cartesian', 'rectilinear', or 'curvilinear'"
9581  error stop "&geometry: grid='rectilinear' requires ndim=1"
9582  error stop "&geometry: grid='curvilinear' requires ndim=2"
959   error stop '&geometry: sdffile and point must not be on same record'
960   error stop '&geometry: sdffile must be specified for all or no regions'

970   error stop '&geometry: point must only be specified for ndim=2'
971   error stop '&geometry: region must only be specified for ndim>1'
972   error stop '&geometry: lrobin must only be specified for ndim=1'
973   error stop '&geometry: lspher must only be specified for ndim=1'
974   error stop '&geometry: fracupd must only be specified for ndim=1'
975   error stop '&geometry: posfix must only be specified for ndim=1'
976   error stop '&geometry: posif must only be specified for ndim=1'

980   error stop '&geometry: nionupd is obsolete. Use niongeom of &setup record instead'
981   error stop '&geometry: lperx is obsolete. Use bndtypex instead'
982   error stop '&geometry: xper is obsolete. Use xlim instead'
983   error stop '&geometry: yper is obsolete.'    ! Use ylim instead'

   end subroutine readGeom
!
!=======================================================================
!
   subroutine initGeom( )
!             ========
!
!     Initialize the geom structure
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ir, ip, np, nflag
      real(REAL8) :: x1, x2, z1, z2, eps
      real(REAL8) :: xmin, xmax, zmin, zmax
      real(REAL8), allocatable :: px(:), pz(:)
!
!  Start of Executable Code:
!  ------------------------
!
      geom%lcart = .not.any( [ geom%lcylx, geom%lcyly, geom%lcylz, geom%lspher ] )
      geom%lcell = ( ldyn .or. usecell )
!
!  Set default values
!
      if( .not.defined( geom%fracupd ) ) then
         if( ndim == 1 ) then
            geom%fracupd = 0.1
         else
            geom%fracupd = ZERO
         end if
      end if
      if( .not.defined( geom%wrough ) )  geom%wrough = ZERO
      if( .not.defined( geom%xcen ) )  geom%xcen = ZERO
      if( .not.defined( geom%ycen ) )  geom%ycen = ZERO
      if( .not.defined( geom%zcen ) )  geom%zcen = ZERO

      if( .not.usecell ) then
         if( nr == 1 .and. .not.defined( geom%posif(0) ) ) then
            geom%posif(0) = -0.001
            geom%posif(1) = BIG
         end if

         if( .not.defined( geom%posfix ) ) then
            geom%posfix = HALF * ( geom%posif(0) + geom%posif(nr) )
         end if
      end if

      if( .not.defined( geom%bndtypex(1) ) )  geom%bndtypex(1) = ENUM%bndtype%free
      if( .not.defined( geom%bndtypex(2) ) )  geom%bndtypex(2) = geom%bndtypex(1)
      if( .not.defined( geom%bndtypey(1) ) )  geom%bndtypey(1) = ENUM%bndtype%free
      if( .not.defined( geom%bndtypey(2) ) )  geom%bndtypey(2) = geom%bndtypey(1)
      if( .not.defined( geom%bndtypez(1) ) )  geom%bndtypez(1) = ENUM%bndtype%free
      if( .not.defined( geom%bndtypez(2) ) )  geom%bndtypez(2) = geom%bndtypez(1)
      if( .not.defined( geom%grid ) )  geom%grid = ENUM%grid%regular

      ! NOTE: geom%nbox, geom%wbox, geom%nbox2, geom%wbox2 might have been defined in readOut

      ! geom%nbox and geom%wbox are only needed for ndim=1, ldyn=t, usecell=f
      if( .not.defined( geom%nbox ) )  geom%nbox = 32
      if( .not.defined( geom%wbox ) )  geom%wbox = ONE

      ! geom%nbox2 and geomwbox2 are only needed for ndim=2, usesdf=f
      ! The determine the internal mesh. If ldyn=t, they also determine the cell size and
      ! To small cells cause much noise.
      if( ldyn ) then
         if( .not.defined( geom%nbox2(1) ) )  geom%nbox2(1) = 32
      else
         if( .not.defined( geom%nbox2(1) ) )  geom%nbox2(1) = 1000
      end if
      if( .not.defined( geom%nbox2(2) ) )  geom%nbox2(2) = geom%nbox2(1)
      if( .not.defined( geom%wbox2(1) ) )  geom%wbox2(1) = ONE
      if( .not.defined( geom%wbox2(2) ) )  geom%wbox2(2) = geom%wbox2(1)
         ! NOTE: If ldyn=t and out%wbox2 is defined, geom%wbox2 was set to out%wbox2 in
         !       readOut
!
!  Set default values, snap to target limits, and check validity of xlim, zlim   ! , ylim
!
      !geom%ylim(1:2) = [ NEGBIG, BIG ]          ! 3D geometry can only be defined by SDF files

      if( usecell .or. geom%usesdf ) then

         geom%xlim(1:2) = [ NEGBIG, BIG ]

         if( geom%lcart ) then
            geom%zlim(1:2) = [ NEGBIG, BIG ]
         else
            geom%zlim(1:2) = [ ZERO, BIG ]
         end if

      else                                      ! this implies ndim<3

         eps = TEN**( -precision( 1. ) )

         if( ndim == 2 ) then
            xmin = minval( geom%region(0)%px(:) )
            if( defined( geom%xlim(1) ) ) then
               if( abs( geom%xlim(1) - xmin ) <  &
                   eps * max( abs( geom%xlim(1) ), abs( xmin ) ) )  &
                     geom%xlim(1) = xmin
               if( geom%xlim(1) > xmin )  goto 901
            else
               geom%xlim(1) = xmin
            end if

            xmax = maxval( geom%region(0)%px(:) )
            if( defined( geom%xlim(2) ) ) then
               if( abs( geom%xlim(2) - xmax ) <  &
                   eps * max( abs( geom%xlim(2) ), abs( xmax ) ) )  &
                     geom%xlim(2) = xmax
               if( geom%xlim(2) < xmax )  goto 902
            else
               geom%xlim(2) = xmax
            end if
         end if

         if( ndim == 1 ) then
            zmin = geom%posif(0)
         else
            zmin = minval( geom%region(0)%pz(:) )
         end if
         if( defined( geom%zlim(1) ) ) then
            if( abs( geom%zlim(1) - zmin ) < eps * max( abs( geom%zlim(1) ), abs( zmin ) ) )  &
                  geom%zlim(1) = zmin
            if( geom%zlim(1) > zmin )  goto 903
         else
            geom%zlim(1) = zmin
         end if

         if( ndim == 1 ) then
            zmax = geom%posif(nr)
         else
            zmax = maxval( geom%region(0)%pz(:) )
         end if
         if( defined( geom%zlim(2) ) ) then
            if( abs( geom%zlim(2) - zmax ) < eps * max( abs( geom%zlim(2) ), abs( zmax ) ) )  &
                  geom%zlim(2) = zmax
            if( geom%zlim(2) < zmax )  goto 904
         else
            geom%zlim(2) = zmax
         end if

         if( .not.geom%lcart ) then
            if( geom%zlim(1) < ZERO )  goto 905
         end if

      end if
!
!  Check validity of the parameters
!
      if( ndim == 1 .and. .not.usecell .and. .not.geom%usesdf ) then
         if( any( geom%posif(1:nr) <= geom%posif(0:nr-1) ) )  goto 900
      end if

      if( geom%lcell .and. .not.geom%lcart )  goto 913
      if( geom%lcell .and. geom%lrobin )  goto 914

      nflag = count( [ geom%lcylx, geom%lcyly, geom%lcylz, geom%lspher, geom%lrobin ] )
      if( nflag > 1 )  goto 911
      if( nflag > 0 .and. ldyn )  goto 912
      if( ndim > 1 .and. geom%lspher )  goto 964
      if( ndim > 2 .and. any( [ geom%lcylx, geom%lcyly, geom%lcylz ] ) )  goto 965

      if( geom%fracupd < ZERO )  goto 954
      if( geom%wrough < ZERO ) goto 955

      if( ( geom%bndtypex(1) == ENUM%bndtype%periodic .and.  &
            geom%bndtypex(2) == ENUM%bndtype%free ) .or.  &
          ( geom%bndtypex(2) == ENUM%bndtype%periodic .and.  &
            geom%bndtypex(1) == ENUM%bndtype%free ) )  goto 920
      if( ( geom%bndtypey(1) == ENUM%bndtype%periodic .and.  &
            geom%bndtypey(2) == ENUM%bndtype%free ) .or.  &
          ( geom%bndtypey(2) == ENUM%bndtype%periodic .and.  &
            geom%bndtypey(1) == ENUM%bndtype%free ) )  goto 921
      if( ( geom%bndtypez(1) == ENUM%bndtype%periodic .and.  &
            geom%bndtypez(2) == ENUM%bndtype%free ) .or.  &
          ( geom%bndtypez(2) == ENUM%bndtype%periodic .and.  &
            geom%bndtypez(1) == ENUM%bndtype%free ) )  goto 922

      if( ( geom%bndtypex(1) == ENUM%bndtype%periodic .and.  &
            geom%bndtypex(2) == ENUM%bndtype%extend ) .or.  &
          ( geom%bndtypex(2) == ENUM%bndtype%periodic .and.  &
            geom%bndtypex(1) == ENUM%bndtype%extend ) )  goto 923
      if( ( geom%bndtypey(1) == ENUM%bndtype%periodic .and.  &
            geom%bndtypey(2) == ENUM%bndtype%extend ) .or.  &
          ( geom%bndtypey(2) == ENUM%bndtype%periodic .and.  &
            geom%bndtypey(1) == ENUM%bndtype%extend ) )  goto 924
      if( ( geom%bndtypez(1) == ENUM%bndtype%periodic .and.  &
            geom%bndtypez(2) == ENUM%bndtype%extend ) .or.  &
          ( geom%bndtypez(2) == ENUM%bndtype%periodic .and.  &
            geom%bndtypez(1) == ENUM%bndtype%extend ) )  goto 925

      if( all( geom%bndtypex(:) == ENUM%bndtype%reflect ) )  goto 926
      if( all( geom%bndtypey(:) == ENUM%bndtype%reflect ) )  goto 927
      if( all( geom%bndtypez(:) == ENUM%bndtype%reflect ) )  goto 928

      if( any( geom%bndtypex(:) == ENUM%bndtype%reflect ) .and. usecell )  goto 929
      if( any( geom%bndtypey(:) == ENUM%bndtype%reflect ) .and. usecell )  goto 930
      if( any( geom%bndtypez(:) == ENUM%bndtype%reflect ) .and. usecell )  goto 931

      if( ndim >= 2 .and. geom%xlim(2) <= geom%xlim(1) )  goto 961
      !if( ndim == 3 .and. geom%ylim(2) <= geom%ylim(1) )  goto 962
      if( geom%zlim(2) <= geom%zlim(1) )  goto 963

      if( geom%nbox < 0 )  goto 9280
      if( geom%wbox <= 0 )  goto 9283

      if( geom%nbox2(1) < 0 )  goto 9281
      if( geom%nbox2(2) < 0 )  goto 9281
      if( geom%wbox2(1) <= 0 )  goto 9282
      if( geom%wbox2(2) <= 0 )  goto 9282
!
!  Modify the internal 1-D mesh so that an integer number of mesh boxes, not exceeding
!  geom%nbox, fits between geom%zlim(:)
!
      if( ndim == 1 ) then
         if( geom%zlim(2) /= BIG ) then
            if( geom%nbox == 0 ) then
               geom%nbox = max( 1, nint( (geom%zlim(2)-geom%zlim(1))/geom%wbox ) )
            else
               geom%nbox = min( geom%nbox,  &
                                max( 1, nint( (geom%zlim(2)-geom%zlim(1))/geom%wbox ) ) )
            end if
            geom%wbox = (geom%zlim(2)-geom%zlim(1)) / geom%nbox
         end if
      end if
!
!  Note that at this point we cannot include the vacuum layer in geom%zlim(:) at free
!  boundaries because sn%pmax is not yet available
!
!
!  The remainder of this subroutine applies only to 2D geometry defined by polygons.
!  We can skip this if ndim/=2, usecell=t, or usesdf=t
!
      if( ndim /= 2 .or. usecell .or. geom%usesdf )  return
!     -----------------------------------------------------
!
!  Initialize region polygons
!
      do ir = 0, nr
         call initPolygon( geom%region(ir) )
!             ===========        io
      end do
!
!  Adapt surface polygon taking the boundary conditions into account
!
!  NOTE: We currently cannot do this with the region polygons because the procedure mcdis4,
!        which is called for region polygons, does not take virtual segments into account.
!        This may lead to some loss in performance.
!
      call updateRegionPolygon( geom%region(0) )
!          ===================       io
!
!  Modify the internal mesh so that an integer number of mesh boxes, not exceeding geom%nbox2,
!  fits between geom%xlim(:) and geom%zlim(:)
!
      if( geom%xlim(2) /= BIG ) then
         if( geom%nbox2(1) == 0 ) then
            geom%nbox2(1) = max( 1, nint( (geom%xlim(2)-geom%xlim(1))/geom%wbox2(1) ) )
         else
            geom%nbox2(1) = min( geom%nbox2(1),  &
                                 max( 1, nint( (geom%xlim(2)-geom%xlim(1))/geom%wbox2(1) ) ) )
         end if
         geom%wbox2(1) = (geom%xlim(2)-geom%xlim(1)) / geom%nbox2(1)
      end if

      if( geom%zlim(2) /= BIG ) then
         if( geom%nbox2(2) == 0 ) then
            geom%nbox2(2) = max( 1, nint( (geom%zlim(2)-geom%zlim(1))/geom%wbox2(2) ) )
         else
            geom%nbox2(2) = min( geom%nbox2(2),  &
                                 max( 1, nint( (geom%zlim(2)-geom%zlim(1))/geom%wbox2(2) ) ) )
         end if
         geom%wbox2(2) = (geom%zlim(2)-geom%zlim(1)) / geom%nbox2(2)
      end if
!
!  In case of dynamic simulation and artificial surface, check whether the two
!  left/right/front/bottom most points are parallel to vertical/horizontal lines
!
!  This is a leftover of the old 2-D dynamic mode. It should be adapted to new
!  requirements
!
      if( ldyn ) then   ! maybe this can be applied for non dynamic simulation?

         ! not sure if this works when also lperx=t
         np = geom%region(0)%np
         allocate( px( np ), pz( np ) )
         px(1:np) = geom%region(0)%px(1:np)
         pz(1:np) = geom%region(0)%pz(1:np)

         if( .not.geom%freex(1) ) then
            ip = minloc(px(1:np),dim=1)
            x1 = px(ip)
            if( ip == 1 ) then                  ! special treatment for the first point
               x2 = min( px(2), px(np) )
            else                                ! find which of the neighboring points has minx
               x2 = min( px(ip-1), px(ip+1) )
            end if
            if( abs( x1 - x2 ) > max( SMALL, EPSMB * abs(x1) ) )  goto 950
         end if

         if( .not.geom%freex(2) ) then
            ip = maxloc( px(1:np), dim=1 )
            x1 = px(ip)
            if( ip == 1 ) then                  ! special treatment for the first point
               x2 = max( px(np), px(2)  )
            else                                ! find which of the neighboring points has minx
               x2 = max( px(ip-1), px(ip+1) )
            end if
            if( abs( x1 - x2 ) > max( SMALL, EPSMB * abs(x1) ) )  goto 950
         end if

         if( .not.geom%freez(1) ) then
            ip = minloc(pz(1:np),dim=1)
            z1 = pz(ip)
            if( ip == 1 ) then                  ! special treatment for the first point
               z2 = min( pz(2), pz(np) )
            else                                ! find which of the neighboring points has minz
               z2 = min( pz(ip-1), pz(ip+1) )
            end if
            if( abs( z1 - z2 ) > max( SMALL, EPSMB * abs(z1) ) )  goto 950
         end if

         if( .not.geom%freez(2) ) then
            ip = maxloc( pz(1:np), dim=1 )
            z1 = pz(ip)
            if( ip == 1 ) then                  ! special treatment for the first point
               z2 = max( pz(np), pz(2)  )
            else                                ! find which of the neighboring points has minz
               z2 = max( pz(ip-1), pz(ip+1) )
            end if
            if( abs( z1 - z2 ) > max( SMALL, EPSMB * abs(z1) ) )  goto 950
         end if

         deallocate( px, pz )

      end if

      return

900   error stop '&geometry: posif values must be increasing'
901   error stop '&geometry: xlim(1)<=min(pos(1)) required.'
902   error stop '&geometry: xlim(2)>=max(pos(1)) required.'
903   error stop '&geometry: zlim(1)<=min(posif or pos(2)) required.'
904   error stop '&geometry: zlim(2)>=max(posif or pos(2)) required.'
905   error stop '&geometry: zlim(1)>=0 required for lspher=t, lcylx=t, lcyly=t, or lcylz=t'

911   error stop '&geometry: only one of lcylx=t, lcyly=t, lcylz=t, lspher=t, lrobin=t allowed'
912   error stop '&geometry: lcylx=t, lcyly=t, lcylz=t, lspher=t, and lrobin=t require ldyn=f'
913   error stop '&geometry: lcylx=t, lcyly=t, lcylz=t, lspher=t require ldyn=f, usecell=f'
914   error stop '&geometry: lrobin=t requires ldyn=f, usecell=f'

920   error stop '&geometry: bndtypex: periodic and free boundaries must not be combined'
921   error stop '&geometry: bndtypey: periodic and free boundaries must not be combined'
922   error stop '&geometry: bndtypez: periodic and free boundaries must not be combined'
923   error stop '&geometry: bndtypex: periodic and extending boundaries must not be combined'
924   error stop '&geometry: bndtypey: periodic and extending boundaries must not be combined'
925   error stop '&geometry: bndtypez: periodic and extending boundaries must not be combined'
926   error stop '&geometry: Not both bndtypex must be reflective'
927   error stop '&geometry: Not both bndtypey must be reflective'
928   error stop '&geometry: Not both bndtypez must be reflective'
929   error stop '&geometry: bndtypex reflective not allowed for usecell=t'
930   error stop '&geometry: bndtypey reflective not allowed for usecell=t'
931   error stop '&geometry: bndtypez reflective not allowed for usecell=t'

9280  error stop '&geometry: nbox>=0 required'
9281  error stop '&geometry: nbox2>=0 required'
9282  error stop '&geometry: wbox>0 required'
9283  error stop '&geometry: wbox2>0 required'

950   error stop '&geometry: line segment closest to artificial boundary must be parallel to it'
954   error stop '&geometry: fracupd>=0 required'
955   error stop '&geometry: wrough >= 0 required'

961   error stop '&geometry: xlim(2)>xlim(1) required'
!962   error stop '&geometry: ylim(2)>ylim(1) required'
963   error stop '&geometry: zlim(2)>zlim(1) required'
964   error stop '&geometry: lspher=t requires ndim=1'
965   error stop '&geometry: lcylx=t, lcyly=t, lcylz=t require ndim=1 or 2'

970   error stop '&geometry: Region completely at negative radii'

   end subroutine initGeom
!
!=======================================================================
!
   subroutine writeGeom( )
!             =========
!
!     Write parameters derived from the &geometry record to the output file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ir, ip, i, ilast, lenmax
!
!  Start of Executable Code:
!  ------------------------
!
!  Open the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i
!
      write( 21, '(t2,a)' )  &
            '&geometry: ' // str( ndim ) // 'd geometry'

      if( geom%wrough > ZERO )  &
         write( 21, '(t13,a,f8.3)' )  &
               'Surface roughness layer of thickness wrough=' // str( geom%wrough, 'f99.3' )

      if( geom%lcylx )  &
         write( 21, '(t13,a)' )  &
               'lcylx=t: Cylindrical symmetry, center=(' // str( geom%xcen, 'f99.1' ) //  &
               ',' // str( geom%ycen, 'f99.1' ) // ',' // str( geom%zcen, 'f99.1' ) // ')'
      if( geom%lcyly )  &
         write( 21, '(t13,a)' )  &
               'lcyly=t: Cylindrical symmetry, center=(' // str( geom%xcen, 'f99.1' ) //  &
               ',' // str( geom%ycen, 'f99.1' ) // ',' // str( geom%zcen, 'f99.1' ) // ')'
      if( geom%lcylz )  &
         write( 21, '(t13,a)' )  &
               'lcylz=t: Cylindrical symmetry, center=(' // str( geom%xcen, 'f99.1' ) //  &
               ',' // str( geom%ycen, 'f99.1' ) // ',' // str( geom%zcen, 'f99.1' ) // ')'
      if( geom%lspher )  &
         write( 21, '(t13,a)' )  &
               'lspher=t: Spherical layers, center=(' // str( geom%xcen, 'f99.1' ) //  &
               ',' // str( geom%ycen, 'f99.1' ) //  &
               ',' // str( geom%zcen, 'f99.1' ) // ')'
      if( geom%lrobin )  &
         write( 21, '(t13,a)' )  &
               'lrobin=t: Layer laterally bounded'

      if( ldyn .and. .not. usecell ) then
         if( ndim == 1 ) then
            write( 21, '(t13,a)' )  &
                  'Maximum number of cells: nbox=' // str( geom%nbox )
            write( 21, '(t13,a)' ) &
                  'Cell width: wbox=' // str( geom%wbox, 'f99.3' )
         else if( ndim == 2 ) then
            write( 21, '(t13,a)' )  &
                  'Maximum number of cell boxes: nbox2=' //  &
                  str( geom%nbox2(1) ) // ',' // str( geom%nbox2(2) )
            write( 21, '(t13,a)' ) &
                  'Cell widths: wbox2=' //  &
                  str( geom%wbox2(1), 'f99.3' ) // ',' // str( geom%wbox2(2), 'f99.3' )
         end if
      else
         if( ndim == 2 .and. .not.geom%usesdf ) then
            write( 21, '(t13,a)' )  &
                  'Maximum number of boxes for internal mesh: nbox2=' //  &
                  str( geom%nbox2(1) ) //  &
                  ',' // str( geom%nbox2(2) )
            write( 21, '(t13,a)' ) &
                  'Box widths for internal mesh: wbox2=' //  &
                  str( geom%wbox2(1), 'f99.3' ) //  &
                  ',' // str( geom%wbox2(2), 'f99.3' )
         end if
      end if

      if( ndim > 1 ) then
         write( 21, '(t13,a)', advance='no' )  &
               'Left/right (x) boundary type: '
         if( geom%bndtypex(1) == ENUM%bndtype%free ) then
            write( 21, '(a)', advance='no' )  &
                  'free/'
         else if( geom%bndtypex(1) == ENUM%bndtype%extend ) then
            write( 21, '(a)', advance='no' )  &
                  'extend/'
         else if( geom%bndtypex(1) == ENUM%bndtype%periodic ) then
            write( 21, '(a)', advance='no' )  &
                  'periodic/'
         else if( geom%bndtypex(1) == ENUM%bndtype%reflect ) then
            write( 21, '(a)', advance='no' )  &
                  'reflective/'
         end if
         if( geom%bndtypex(2) == ENUM%bndtype%free ) then
            write( 21, '(a)', advance='no' )  &
                  'free'
         else if( geom%bndtypex(2) == ENUM%bndtype%extend ) then
            write( 21, '(a)', advance='no' )  &
                  'extending'
         else if( geom%bndtypex(2) == ENUM%bndtype%periodic ) then
            write( 21, '(a)', advance='no' )  &
                  'periodic'
         else if( geom%bndtypex(2) == ENUM%bndtype%reflect ) then
            write( 21, '(a)', advance='no' )  &
                  'reflective'
         end if
         if( any( geom%bndtypex(:) == ENUM%bndtype%periodic ) ) then
            write( 21, '(a)' )  &
               ', xlim=' // str( geom%xlim(1), 'f99.3' ) // ',' // str( geom%xlim(2), 'f99.3' )
         else
            write( 21, '()' )
         end if
      end if

      if( ndim == 3 ) then
         write( 21, '(t13,a)', advance='no' )  &
               'Front/back (y) boundary type: '
         if( geom%bndtypey(1) == ENUM%bndtype%free ) then
            write( 21, '(a)', advance='no' )  &
                  'free/'
         else if( geom%bndtypey(1) == ENUM%bndtype%extend ) then
            write( 21, '(a)', advance='no' )  &
                  'extend/'
         else if( geom%bndtypey(1) == ENUM%bndtype%periodic ) then
            write( 21, '(a)', advance='no' )  &
                  'periodic/'
         else if( geom%bndtypey(1) == ENUM%bndtype%reflect ) then
            write( 21, '(a)', advance='no' )  &
                  'reflective/'
         end if
         if( geom%bndtypey(2) == ENUM%bndtype%free ) then
            write( 21, '(a)', advance='no' )  &
                  'free'
         else if( geom%bndtypey(2) == ENUM%bndtype%extend ) then
            write( 21, '(a)', advance='no' )  &
                  'extending'
         else if( geom%bndtypey(2) == ENUM%bndtype%periodic ) then
            write( 21, '(a)', advance='no' )  &
                  'periodic'
         else if( geom%bndtypey(2) == ENUM%bndtype%reflect ) then
            write( 21, '(a)', advance='no' )  &
                  'reflective'
         end if
         if( any( geom%bndtypey(:) == ENUM%bndtype%periodic ) ) then
            write( 21, '(a)' )  &
               ', ylim=' // str( geom%ylim(1), 'f99.3' ) // ',' // str( geom%ylim(2), 'f99.3' )
         else
            write( 21, '()' )
         end if
      end if

      write( 21, '(t13,a)', advance='no' )  &
            'Upper/lower (z) boundary type: '
      if( geom%bndtypez(1) == ENUM%bndtype%free ) then
         write( 21, '(a)', advance='no' )  &
               'free/'
      else if( geom%bndtypez(1) == ENUM%bndtype%extend ) then
         write( 21, '(a)', advance='no' )  &
               'extend/'
      else if( geom%bndtypez(1) == ENUM%bndtype%periodic ) then
         write( 21, '(a)', advance='no' )  &
               'periodic/'
      else if( geom%bndtypez(1) == ENUM%bndtype%reflect ) then
         write( 21, '(a)', advance='no' )  &
               'reflective/'
      end if
      if( geom%bndtypez(2) == ENUM%bndtype%free ) then
         write( 21, '(a)' )  &
               'free'
      else if( geom%bndtypez(2) == ENUM%bndtype%extend ) then
         write( 21, '(a)' )  &
               'extending'
      else if( geom%bndtypez(2) == ENUM%bndtype%periodic ) then
         write( 21, '(a)' )  &
               'periodic'
      else if( geom%bndtypez(2) == ENUM%bndtype%reflect ) then
         write( 21, '(a)' )  &
               'reflective'
      end if
      if( any( geom%bndtypez(:) == ENUM%bndtype%periodic ) ) then
         write( 21, '(a)' )  &
            ', zlim=' // str( geom%zlim(1), 'f99.3' ) // ',' // str( geom%zlim(2), 'f99.3' )
      else
         write( 21, '()' )
      end if

      if( ldyn ) then
         if( ndim == 1 ) then
            write( 21, '(t13,a)' )  &
                  'Cell update if fracupd=' // str( geom%fracupd, 'f99.4' ) // ' is exceeded'
            write( 21, '(t13,a)' )  &
                  '   such that z=' // str( geom%posfix, 'f99.1' ) // ' remains fixed'
         end if

         if( geom%grid == ENUM%grid%regular ) then
            write( 21, '(/t13,a)' )  &
                  'Regular cell grid'
         else if( geom%grid == ENUM%grid%rectilinear ) then
            write( 21, '(/t13,a)' )  &
                  'Rectilinear cell grid'
         else if( geom%grid == ENUM%grid%curvilinear ) then
            write( 21, '(/t13,a)' )  &
                  'Curvilinear cell grid'
         end if
      end if
      write( 21, '()' )

      if( usecell ) then
         write( 21, '(t13,a)' )  'Geometry defined by cell file'
      else if( ndim == 1 ) then
         do ir = 1, nr
            write( 21, '(t13,a)' )  &
                  'Region ' // str( ir ) //  &
                  ': between z=' // str( geom%posif(ir-1), 'f99.1' ) //  &
                  ' and z=' // str( geom%posif(ir), 'f99.1' )
         end do
      else if( geom%usesdf ) then
         write( 21, '(t13,a)' )  'Geometry defined by signed distance function files:'
         do ir = 1, nr
            write( 21, '(t13,a)' ) 'Region ' // str( ir ) // ': ' // trim( geom%sdffile(ir) )
         end do
      else if( defined( geom%geomfile ) ) then
         write( 21, '(t13,a)' )  'Geometry file: ' // trim( geom%geomfile )
      else if( ndim == 2 ) then
         do ir = 1, nr
            write( 21, '(t13,a)' ) 'Region ' // str( ir ) // ': bounded by'

            lenmax = 1
            do i = 1, geom%region(ir)%np
               lenmax = max( lenmax, len( str( geom%region(ir)%px(i), 'f99.1' ) ) )
               lenmax = max( lenmax, len( str( geom%region(ir)%pz(i), 'f99.1' ) ) )
            end do
            do ip = 1, geom%region(ir)%np, 5
               ilast = min( ip+4, geom%region(ir)%np )
               write( 21, '(t13,a,5(f' // str( lenmax ) // '.1:a))' ) &
                  '   x(' // str( ip ) //  &
                  ':' // str( ilast ) // ')=',  &
                  (geom%region(ir)%px(i), ',', i=ip,ilast-1), geom%region(ir)%px(ilast)
               write( 21, '(t13,a,5(f' // str( lenmax ) // '.1:a))' ) &
                  '   z(' // str( ip ) //  &
                  ':' // str( ilast ) // ')=',  &
                  (geom%region(ir)%pz(i), ',', i=ip,ilast-1), geom%region(ir)%pz(ilast)
            end do
         end do
      else
         error stop 'writeGeom: Internal error'
      end if

      if( ldyn .and. ndim == 2 ) then
         write( 21, '(t13,a)' )  &
               'Left/Right boundary is free: ' // str( geom%freex(1) ) //  &
               '/' // str( geom%freex(2) )
         write( 21, '(t13,a)' )  &
               'Top/Bottom boundary is free: ' // str( geom%freez(1) ) //  &
               '/' // str( geom%freez(2) )
      end if

      write( 21, '()' )

      close( 21 )

   end subroutine writeGeom
!
!=======================================================================
!
   subroutine updateRegionPolygon( region )
!             ===================    io
!
!     Update the region polygon taking the boundary types into account
!
!  Arguments:
!
      type(Polygon_type), intent(inout) :: region  ! region polygon
!
!  Local Variables:
!
      logical, allocatable :: skip(:)
      integer(INT4) :: i, imax, ip, np
      real(REAL8) :: lambda, lambdamax
      real(REAL8), allocatable :: px(:), pz(:), pzleft(:), pzright(:), pzother(:)
      real(REAL8), allocatable :: pxx(:), pzz(:)
      type(LogicalStack) :: virtualStack
      type(Real8Stack) :: pxStack, pzStack, pzleftStack, pzrightStack
!
!  Start of Executable Code:
!  ------------------------
!
!  Abbreviations
!
      np = geom%region(0)%np
      allocate( px( np+1 ), pz( np+1 ) )
      px(:) = region%px(:)
      pz(:) = region%pz(:)
!
!  Simplify polygon by combining successive segments which are parallel to the same axis
!
      allocate( pxx(0:np+1), pzz(0:np+1) )
      pxx(0) = px(np)
      pzz(0) = pz(np)
      pxx(1:np+1) = px(:)
      pzz(1:np+1) = pz(:)
      allocate( skip( np+1 ) )
      skip(:) = .false.

      do ip = 1, np
         skip(ip) = ( ( pxx(ip-1) == pxx(ip) .and. pxx(ip) == pxx(ip+1) ) .or.  &
                      ( pzz(ip-1) == pzz(ip) .and. pzz(ip) == pzz(ip+1) ) )
      end do

      if( count( skip(:) ) > 0 ) then
         deallocate( px, pz )
         px = pack( pxx(1:np+1), mask=.not.skip(1:np+1) )
         pz = pack( pzz(1:np+1), mask=.not.skip(1:np+1) )
         np = size( px ) - 1
         px(np+1) = px(1)
         pz(np+1) = pz(1)
      end if

      deallocate( pxx, pzz )
!
!  If both x boundaries are periodic, add points on boundaries, as necessary, to make sure
!  points on boundaries all have a counterpart on the other boundary
!
      if( all( geom%bndtypex(:) == ENUM%bndtype%periodic ) ) then
!
!  ... Determine surface segments on left and right boundary
!
         call newStack( pzleftStack )
!             ========       o
         call newStack( pzrightStack )
!             ========       o
         do ip = 1, np
            if( px(ip) == geom%xlim(1) .and. px(ip+1) == geom%xlim(1) ) then
               call pushStack( pz(ip), pzleftStack )
!                   =========    i         io
               call pushStack( pz(ip+1), pzleftStack )
!                   =========     i          io
            else if( px(ip) == geom%xlim(2) .and. px(ip+1) == geom%xlim(2) ) then
               call pushStack( pz(ip), pzrightStack )
!                   =========    i          io
               call pushStack( pz(ip+1), pzrightStack )
!                   =========     i           io
            end if
         end do
         call stack2Array( pzleftStack, pzleft )
!             ===========      io         o
         call stack2Array( pzrightStack, pzright )
!             ===========       io          o
!
!  ... Cycle segments so that z-values are monotonous
!
         do i = 3, size( pzleft ), 2
            lambda = ( pzleft(i) - pzleft(i-2) ) / ( pzleft(i-1) - pzleft(i-2) )
            if( lambda < 0 )  exit
         end do
         pzleft = cshift( pzleft, i-1 )

         do i = 3, size( pzright ), 2
            lambda = ( pzright(i) - pzright(i-2) ) / ( pzright(i-1) - pzright(i-2) )
            if( lambda < 0 )  exit
         end do
         pzright = cshift( pzright, i-1 )
!
!  ... Subdivide segments as necessary
!
         call newStack( pxStack )
!             ========     o
         call newStack( pzStack )
!             ========     o
         call newStack( virtualStack )
!             ========       o

         do ip = 1, np

            call pushStack( px(ip), pxStack )
!                =========    i       io
            call pushStack( pz(ip), pzStack )
!                =========    i       io

            ! Cycle if not any segment is on a boundary
            if( px(ip) == geom%xlim(1) .and. px(ip+1) == geom%xlim(1) ) then
               pzother = pzright(:)
            else if( px(ip) == geom%xlim(2) .and. px(ip+1) == geom%xlim(2) ) then
               pzother = pzleft(:)
            else
               call pushStack( .false., virtualStack )
!                   =========     i          io
               cycle
            end if

            ! Determine pzother(imax:imax+1) interval which pz(ip) is in or immediately precedes
            imax = size( pzother ) + 1
            lambdamax = -BIG
            do i = 1, size( pzother )-1, 2
               lambda = ( pz(ip) - pzother(i) ) / ( pzother(i+1) - pzother(i) )
               if( lambda > lambdamax .and. lambda < ONE ) then
                  imax = i
                  lambdamax = lambda
               end if
            end do

            if( lambdamax <= ZERO .or. lambdamax == -BIG ) then
                                                ! pz(ip) not in any pzother(i:i+1)
               call pushStack( .false., virtualStack )
!                   =========     i          io
               if( imax == 1 )  cycle
               lambda = ( pzother(imax-1) - pz(ip) ) / ( pz(ip+1) - pz(ip) )
               if( lambda >= ONE .or. lambda <= ZERO )  cycle  ! .or. lambda <= ZERO ?

               ! pzother(imax-1) inside pz(ip:ip+1)
               call pushStack( px(ip), pxStack )
!                   =========    i       io
               call pushStack( pzother(imax-1), pzStack )
!                   =========        i           io
               imax = imax - 2
            end if
            call pushStack( .true., virtualStack )
!                =========    i          io

            do i = imax, 1, -2
               lambda = ( pzother(i) - pz(ip) ) / ( pz(ip+1) - pz(ip) )
               if( lambda >= ONE )  exit

               ! pzother(i) inside pz(ip:ip+1)
               call pushStack( px(ip), pxStack )
!                   =========    i       io
               call pushStack( pzother(i), pzStack )
!                   =========     i         io
               call pushStack( .false., virtualStack )
!                   =========     i          io
               if( i == 1 )  exit
               lambda = ( pzother(i-1) - pz(ip) ) / ( pz(ip+1) - pz(ip) )
               if( lambda >= ONE )  exit

               !  pzother(i-1) inside pz(ip:ip+1)
               call pushStack( px(ip), pxStack )
!                   =========    i       io
               call pushStack( pzother(i-1), pzStack )
!                   =========       i         io
               call pushStack( .true., virtualStack )
!                   =========    i          io
            end do

         end do

         call pushStack( ZERO, pxStack )        ! placeholder for last=first point
!             =========   i      io
         call pushStack( ZERO, pzStack )
!             =========   i      io
         call pushStack( .false., virtualStack )
!             =========     i          io

         region%np = sizeOfStack( pxStack ) - 1
!                    ===========    i
         if( region%np < 0 )  goto 980

         if( allocated( region%px ) )  deallocate( region%px )
         if( allocated( region%pz ) )  deallocate( region%pz )
         if( allocated( region%virtual ) )  deallocate( region%virtual )
         call stack2Array( pxStack, region%px )
!             ===========    io         o
         call stack2Array( pzStack, region%pz )
!             ===========    io         o
         call stack2Array( virtualStack, region%virtual )
!             ===========       io              o
         region%px(np+1) = region%px(1)
         region%pz(np+1) = region%pz(1)

      else

         if( allocated( region%virtual ) )  deallocate( region%virtual )
         call alloc( region%virtual, np+1, 'virtual', 'updateRegionPolygon' )
!             =====        o          i        i              i
         do ip = 1, np
            region%virtual(ip) =  &
                  ( px(ip) == geom%xlim(1) .and. px(ip+1) == geom%xlim(1) .and.  &
                    geom%bndtypex(1) /= ENUM%bndtype%free ) .or.  &
                  ( px(ip) == geom%xlim(2) .and. px(ip+1) == geom%xlim(2) .and.  &
                    geom%bndtypex(2) /= ENUM%bndtype%free )
         end do

      end if

      return

960   error stop '&geometry: Whole target outside histogram boundaries'
970   error stop '&geometry: Whole target completely at negative radii'
980   error stop 'updateRegionPolygon: Internal error'

   end subroutine updateRegionPolygon
!
!=======================================================================
!
end module Param_Geom_mod

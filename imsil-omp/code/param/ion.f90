module Param_Ion_mod
!
!  Hold data about the ions and manages their initialization and output. 
!  The ions may be atoms to be implanted, molecules, or recoils starting 
!  from the interior of the target.
!
!  Nomenclature: Ion = Atom/molecule/recoil to be implanted or started, 
!                      before interaction with the target
!
!                Projectile = Particle moving in the target 
!
!  Uses public data (ion structure) and private namelist (&ions)
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: REAL8, INT4,  &
                              ZERO, ONE, TWO, EIGHT, PI,  &
                              ENUM
   use String_mod,      only: str
   use Undefine_mod,    only: undefine, defined
   use Stack_mod,       only: Int4Stack, newStack, popStack, pushStack
   use Open_mod,        only: openFile

   use Param_Setup_mod, only: nelem, nelemIon, nelemTarg, ndim, nion, nion1, nthread,  &
                              ldyn, usebk, useion, bkfile
   use Param_Elem_mod,  only: elem, anaMolecule
   use Param_Geom_mod,  only: geom

   use Element_mod,     only: getIonMass  ! to be eliminated after splitting into readIon and initIon
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
   public :: readIon                         ! read &ions record from input file
!            =======
   public :: initIon                         ! initialize ion structure
!            =======
   public :: writeIon                        ! write data derived from &ions record
!            ========                        !    to the output file
!
!  Public Data:
!  ===========
!
   type :: Ion_type
!          ========
!
      logical :: ldet                        ! same as dam%ldet (ugly)
      logical :: liinit                      ! internal starting point flag
      logical :: llinit                      ! start as lattice atom flag
      logical :: lpinit                      ! plane implantation flag
      logical :: lrinit                      ! circle implantation flag
      logical :: ltinit                      ! start as target atom flag
      logical :: ranrot                      ! random rotation angle
      logical :: reenter                     ! reenter flag for ions
      logical :: usetra                      ! use trajectory file to initialize ions

      character(4) :: doseunits              ! units of dose
      character(80) :: name                  ! chemical name of the ion (e.g. 'BF2')
      character(80) :: trafile               ! file containing initial conditions for each ion

      integer(INT4) :: irxtal                ! same as irxtal from mat structure (ugly)
      integer(INT4), allocatable :: moddiv(:)! (nbeam) divergence model
      integer(INT4), allocatable :: natom(:) ! (nelemIon) => number of atoms of kind ielem in
                                             !    the ion
      integer(INT4) :: nbeam                 ! number of beams

      real(REAL8) :: dinit                   ! square edge length
      real(REAL8), allocatable :: diverg(:,:)! (2,nbeam) beam divergence in beam x and y
                                             !    direction [deg]
      real(REAL8) :: dose                    ! implant dose [1/A^2]
      real(REAL8) :: ddose                   ! dose increment per ion [1/A^2]
      real(REAL8) :: dose1                   ! 0 or read from backup file
      real(REAL8) :: dxinit                  ! xinit(2) - xinit(1)
      real(REAL8) :: dyinit                  ! yinit(2) - yinit(1)
      real(REAL8) :: dzinit                  ! zinit(2) - zinit(1)
      real(REAL8) :: e                       ! ion energy [eV]
      real(REAL8) :: m                       ! ion molecule mass [amu]
      real(REAL8) :: rinit                   ! radius/distance for circle/plane implant
      real(REAL8) :: rotate                  ! rotation angle [deg]
      real(REAL8) :: sigma(2)                ! standard deviation of beam broadening in x and y
                                             !    direction [A]
      real(REAL8) :: tiltx                   ! tilt angle about x axis [deg]
      real(REAL8) :: tilty                   ! tilt angle about y axis [deg]
      real(REAL8), allocatable :: weight(:)  ! (nbeam) weight of the beam
      real(REAL8) :: xinit(2)                ! x coordinate(s) of initial point (window)
      real(REAL8) :: yinit(2)                ! y coordinate(s) of initial point (window)
      real(REAL8) :: zinit(2)                ! z coordinate(s) of initial point (window)

   end type Ion_type

   type(Ion_type), public :: ion                     ! ion structure
!
!=======================================================================
!
contains
!
!==============================================================================
!
   subroutine readIon( )
!             =======
!
!     Read the &ions record from the input file.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: atom1                    ! atom index
      integer(INT4) :: beam                     ! beam index

      logical :: liinit                         ! internal starting point flag
      logical :: llinit                         ! start as lattice atom flag
      logical :: lpinit                         ! plane implantation flag
      logical :: lrinit                         ! circle implantation flag
      logical :: ltinit                         ! start as target atom flag
      logical :: ranrot                         ! random rotation angle flag
      logical :: reenter                        ! reenter flag for ions
      logical :: usetra                         ! use trajectory file for ion ininitialization
      logical :: volatile                       ! volatility flag

      character(4) :: doseunits                 ! units of dose
      character(80) :: filtra                   ! obsolete
      character(80) :: trafile                  ! file containing initial conditions for each
                                                ! ion (TRA format)
      character(80) :: name                     ! chemical name of ion species

      integer(INT4) :: moddiv                   ! beam divergence model
      integer(INT4) :: nion_setup               ! number of simulated ions from &setup

      real(REAL8) :: dinit                      ! edge length of square where ions start if
                                                !    lpinit=t
      real(REAL8) :: diverg(2)                  ! beam divergence in x and y direction [deg]
      real(REAL8) :: dose                       ! ion dose [ions/cm2]
      real(REAL8) :: energy                     ! ion energy [keV]
      real(REAL8) :: fwhm(2)                    ! full width at half maximum of Gaussian beam
                                                !    in x and y direction [A]
      real(REAL8) :: mass                       ! atom mass of ion [amu]
      real(REAL8) :: rinit                      ! radius (distance) for circle (plane) if
                                                !    lrinit=t (lpinit=t)
      real(REAL8) :: rotate                     ! rotation angle of ion beam [deg]
      real(REAL8) :: tilt                       ! tilt angle of ion beam [deg]
      real(REAL8) :: tiltx                      ! tilt angle about x axis [deg]
      real(REAL8) :: tilty                      ! tilt angle about y axis [deg]
      real(REAL8) :: weight                     ! weight of the beam
      real(REAL8) :: xinit(2)                   ! x coordinate(s) of initial point (window)
      real(REAL8) :: yinit(2)                   ! y coordinate(s) of initial point (window)
      real(REAL8) :: zinit(2)                   ! z coordinate(s) of initial point (window)
      real(REAL8) :: xwin(2)                    ! same as xinit (kept for compatibility)
      real(REAL8) :: ywin(2)                    ! same as yinit (kept for compatibility)
      real(REAL8) :: zwin(2)                    ! same as zinit (kept for compatibility)


      namelist /ions/ atom1, beam,  &
                      liinit, llinit, lpinit, lrinit, ltinit, ranrot, reenter, usetra,  &
                      volatile,  &
                      doseunits, trafile, name,  &
                      moddiv,  &
                      dinit, diverg, dose, energy, fwhm, mass, rinit, rotate, tilt, tiltx,  &
                      tilty, weight, xinit, yinit, zinit, xwin, ywin, zwin
      ! obsolete
      namelist /ions/ filtra,  &
                      nion

      logical, allocatable :: setbeam(:)
      logical :: isEmpty
      character(2) :: ele
      character(80) :: molec                    ! temporary variable holding name
      character(256) :: msg
      integer(INT4) :: ielem, ielem1min, ielem1max, j1, iz, ios, natom1
      integer(INT4) :: ibeam, jbeam, ibeammin, ibeammax
      type(Int4Stack) :: beamStack
!
!  Start of Executable Code:
!  ------------------------
!
      call openFile( 11, 'INP', 'formatted', 'read' )
!          ========  i     i         i         i
!
!  Determine number of beams and number of ion atom species
!  --------------------------------------------------------
!
      ion%nbeam = 0
      call newStack( beamStack )
!          ========      o

      call undefine( name )

      do

         call undefine( beam )
!
!  Read &ions record from the input file
!
         read( 11, nml=ions, iostat=ios, iomsg=msg )
         if( ios < 0 )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 900
         end if
!
!  Handle beam index
!
         if( defined( beam ) ) then
            if( beam < 1 )  goto 906
            call pushStack( beam, beamStack )
!                =========   i       io
            ion%nbeam = max( ion%nbeam, beam )
         end if

      end do

      if( ion%nbeam == 0 ) then
         ion%nbeam = 1
         call pushStack( 1, beamStack )
!             =========  i     io
      end if
!
!   Check that all beam indices are specified up to the maximum index
!
      allocate( setbeam( ion%nbeam ) )
      setbeam(:) = .false.
      do jbeam = 1, ion%nbeam
         call popStack( beamStack, ibeam, isEmpty )
!             ========     io        o      o
         if( isEmpty )  exit
         setbeam(ibeam) = .true.
      end do

      if( any( .not.setbeam(:) ) )  goto 907
      deallocate( setbeam )
!
!  Analyze name for number of ion atom species
!
      if( defined( name ) ) then
         molec = name
         do j1 = 1, nelem
            call anaMolecule( molec, ele, iz, natom1 )
!                ===========   io     o   o     o
            if( iz == 0 )  exit
            nelemIon = j1
         end do
         if( molec /= ' ' )  goto 904
      else
         goto 908
      end if
!
!  Allocate memory and set default values
!  --------------------------------------
!
!  Allocate memory for natom and beam dependent parameters
!
      allocate( ion%natom( nelemIon ) )
      ion%natom(:) = 0
      allocate( ion%moddiv( ion%nbeam ) )

      allocate( ion%diverg( 2, ion%nbeam ) )
      allocate( ion%weight( ion%nbeam ) )

!
!  Set default values of structure variables
!
      ion%liinit = .false.
      ion%llinit = .false.
      ion%lpinit = .false.
      ion%lrinit = .false.
      ion%ltinit = .false.
      ion%ranrot = .false.
      ion%reenter = ( .not.geom%lcart .or. ndim > 1 )
      ion%usetra = .false.

      call undefine( ion%doseunits )
      call undefine( ion%trafile )
      call undefine( ion%name )

      call undefine( ion%moddiv(:) )

      call undefine( ion%dinit )
      call undefine( ion%diverg(:,:) )
      call undefine( ion%dose )
      call undefine( ion%e )
      call undefine( ion%rinit )
      call undefine( ion%rotate )
      call undefine( ion%tiltx )
      call undefine( ion%tilty )
      call undefine( ion%sigma(1:2) )
      call undefine( ion%weight(:) )
      call undefine( ion%xinit(1:2) )
      call undefine( ion%yinit(1:2) )
      call undefine( ion%zinit(1:2) )
!
!  Remember nion from &setup record
!
      nion_setup = nion
!
!  Loop over reading records
!  -------------------------
!
      rewind 11

      do
!
!  Set default values of namelist variables
!
         call undefine( atom1 )
         call undefine( beam )

         reenter = ion%reenter
         liinit = ion%liinit
         llinit = ion%llinit
         lpinit = ion%lpinit
         lrinit = ion%lrinit
         ltinit = ion%ltinit
         ranrot = ion%ranrot
         usetra = ion%usetra
         volatile = .false.

         call undefine( doseunits )
         call undefine( filtra )                ! obsolete
         call undefine( name )
         call undefine( trafile )

         call undefine( moddiv )
         call undefine( nion )

         call undefine( dinit )
         call undefine( diverg(1:2) )
         call undefine( dose )
         call undefine( energy )
         call undefine( fwhm(1:2) )
         call undefine( mass )
         call undefine( rinit )
         call undefine( rotate )
         call undefine( tilt )
         call undefine( tiltx )
         call undefine( tilty )
         call undefine( weight )
         call undefine( xinit(1:2) )
         call undefine( yinit(1:2) )
         call undefine( zinit(1:2) )
         call undefine( xwin(1:2) )
         call undefine( ywin(1:2) )
         call undefine( zwin(1:2) )
!
!  Read &ions record from the input file
!
         read( 11, nml=ions, iostat=ios, iomsg=msg )
         if( ios < 0 )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 900
         end if
!
!  Check obsolete parameters
!
         if( defined( filtra ) )  goto 930
         if( defined( nion ) )  goto 931
!
!  Check index parameters
!
         if( defined( atom1 ) ) then
            if( atom1 < 1 .or. atom1 > nelemIon )  goto 905
         end if
!
!  Check parameters on consistency with index parameters
!
         if( defined( atom1 ) ) then
                                                ! Implant parameters
            if( defined( doseunits ) )  goto 909
            if( defined( name ) )  goto 910
            if( defined( energy ) )  goto 911
            if( defined( tilt ) .or. defined( tiltx ) .or. defined( tilty ) )  goto 912
            if( defined( rotate ) )  goto 913
            if( any( defined( diverg(:) ) ) )  goto 9140
            if( defined( moddiv ) )  goto 9141
            if( any( defined( fwhm(:) ) ) )  goto 915
            if( defined( dose ) )  goto 916
                                                ! starting points of trajectories
            if( any( defined( xinit(:) ) ) )  goto 918
            if( any( defined( yinit(:) ) ) )  goto 919
            if( any( defined( zinit(:) ) ) )  goto 9100
            if( any( defined( xwin(:) ) ) )  goto 9101
            if( any( defined( ywin(:) ) ) )  goto 9102
            if( any( defined( zwin(:) ) ) )  goto 9103
            if( defined( trafile ) )  goto 9104
                                                ! undocumented modes
            if( defined( rinit ) )  goto 9105
            if( defined( dinit ) )  goto 9106
         end if

         if( defined( beam ) ) then
            if( defined( atom1 ) )  goto 9200
            if( defined( doseunits ) )  goto 9201
            if( defined( trafile ) )  goto 9202
            if( defined( name ) )  goto 9203
            if( defined( dinit ) )  goto 9205
            if( defined( dose ) )  goto 9206
            if( defined( energy ) )  goto 9207
            if( any( defined( fwhm(:) ) ) )  goto 9208
            if( defined( mass ) )  goto 9209
            if( defined( rinit ) )  goto 9210
            if( defined( rotate ) )  goto 9211
            if( defined( tilt ) .or. defined( tiltx ) .or. defined( tilty ) )  goto 9212
            if( any( defined( xinit(:) ) ) )  goto 9213
            if( any( defined( yinit(:) ) ) )  goto 9214
            if( any( defined( zinit(:) )) )  goto 9215
            if( any( defined( xwin(:) ) ) )  goto 9216
            if( any( defined( ywin(:) ) ) )  goto 9217
            if( any( defined( zwin(:) ) ) )  goto 9218
         end if

         if( defined( tilt ) .and. defined( tilty ) )  goto 9219
!
!  Define range of index variables to which parameters apply
!
         if( defined( atom1 ) ) then
            ielem1min = atom1
            ielem1max = atom1
         else
            ielem1min = 1
            ielem1max = nelemIon
         end if
         ielem = ielem1min

         if( defined( beam ) ) then
            ibeammin = beam
            ibeammax = beam
         else
            ibeammin = 1
            ibeammax = ion%nbeam
         end if
!
!  Analyze the ion name and initialize the elements
!
         if( defined( name ) ) then
            ion%name = name
                                                ! Get element information from ion species name
            molec = name
            do j1 = 1, nelemIon
               call anaMolecule( molec, ele, iz, natom1 )
!                   ===========   io     o   o     o
               if( iz == 0 )  exit
               ion%natom(j1) = natom1
               elem(j1)%name = ele
               elem(j1)%z = iz
               call getIonMass( iz, elem(j1)%m )
!                   ==========  i      o
            end do
         end if

         if( defined( mass ) ) then
            if( mass <= ZERO )  goto 980
            elem(ielem1min:ielem1max)%m = mass
         end if

         if( volatile )  elem(ielem1min:ielem1max)%volatile = .true.
!
!  Copy parameters to the ion structure
!
         ion%liinit = liinit
         ion%llinit = llinit
         ion%lpinit = lpinit
         ion%lrinit = lrinit
         ion%ltinit = ltinit
         ion%ranrot = ranrot
         ion%reenter = reenter
         ion%usetra = usetra

         if( defined( doseunits ) )  ion%doseunits = adjustl( doseunits )
         if( defined( trafile ) )  ion%trafile = trafile

         if( defined( moddiv ) )  ion%moddiv(ibeammin:ibeammax) = moddiv

         if( defined( dinit ) )  ion%dinit = dinit
         if( defined( diverg(1) ) )  ion%diverg(1,ibeammin:ibeammax) = diverg(1)
         if( defined( diverg(2) ) )  ion%diverg(2,ibeammin:ibeammax) = diverg(2)
         if( defined( dose ) )  ion%dose = dose          ! conversion to 1/A^2 or 1/A later
         if( defined( energy ) )  ion%e = 1000 * energy  ! convert from keV to eV
         if( defined( rotate ) )  ion%rotate = rotate
         if( defined( rinit ) )  ion%rinit = rinit
         if( defined( fwhm(1) ) )  ion%sigma(1) = fwhm(1) / sqrt( EIGHT * log(TWO) )
         if( defined( fwhm(2) ) )  ion%sigma(2) = fwhm(2) / sqrt( EIGHT * log(TWO) )
         if( defined( tilt ) )  ion%tilty = tilt
         if( defined( tiltx ) )  ion%tiltx = tiltx
         if( defined( tilty ) )  ion%tilty = tilty
         if( defined( weight ) )  ion%weight(ibeammin:ibeammax) = weight

         if( defined( xwin(1) ) )  ion%xinit(1) = xwin(1)   ! compatibility
         if( defined( xwin(2) ) )  ion%xinit(2) = xwin(2)
         if( defined( ywin(1) ) )  ion%yinit(1) = ywin(1)
         if( defined( ywin(2) ) )  ion%yinit(2) = ywin(2)
         if( defined( zwin(1) ) )  ion%zinit(1) = zwin(1)
         if( defined( zwin(2) ) )  ion%zinit(2) = zwin(2)
         if( defined( xinit(1) ) )  ion%xinit(1) = xinit(1)
         if( defined( xinit(2) ) )  ion%xinit(2) = xinit(2)
         if( defined( yinit(1) ) )  ion%yinit(1) = yinit(1)
         if( defined( yinit(2) ) )  ion%yinit(2) = yinit(2)
         if( defined( zinit(1) ) )  ion%zinit(1) = zinit(1)
         if( defined( zinit(2) ) )  ion%zinit(2) = zinit(2)

      end do
!
!  End of loop over reading records
!  --------------------------------
!
      close( 11 )
!
!  Restore nion from the &setup record
!
      nion = nion_setup
!
!  Set parameters in other records
!
      if( ldyn ) then
         nelemTarg = nelem
      else
         nelemTarg = nelem - nelemIon
      end if

      return
!
!  Error messages
!
900   error stop '&ions: error during reading namelist'

904   stop '&ions: illegal name'
905   stop '&ions: 1<=atom1<=number of ion atom species required'
906   stop '&ions: beam>=1 required'
907   stop '&ions: beam indices must be contiguous'
908   stop '&ions: name must be specified'

909   error stop '&ions: atom1 and doseunits must not be on same record'
910   error stop '&ions: atom1 and name must not be on same record'
911   error stop '&ions: atom1 and energy must not be on same record'
912   error stop '&ions: atom1 and tilt, tiltx, or tilty must not be on same record'
913   error stop '&ions: atom1 and rotate must not be on same record'
9140  error stop '&ions: atom1 and diverg must not be on same record'
9141  error stop '&ions: atom1 and moddiv must not be on same record'
915   error stop '&ions: atom1 and fwhm must not be on same record'
916   error stop '&ions: atom1 and dose must not be on same record'
918   error stop '&ions: atom1 and xinit must not be on same record'
919   error stop '&ions: atom1 and yinit must not be on same record'
9100  error stop '&ions: atom1 and zinit must not be on same record'
9101  error stop '&ions: atom1 and xwin must not be on same record'
9102  error stop '&ions: atom1 and ywin must not be on same record'
9103  error stop '&ions: atom1 and zwin must not be on same record'
9104  error stop '&ions: atom1 and trafile must not be on same record'
9105  error stop '&ions: atom1 and rinit must not be on same record'
9106  error stop '&ions: atom1 and dinit must not be on same record'

9200  error stop '&ions: beam and atom1 must not be on same record'
9201  error stop '&ions: beam and doseunits must not be on same record'
9202  error stop '&ions: beam and trafile must not be on same record'
9203  error stop '&ions: beam and name must not be on same record'
9205  error stop '&ions: beam and dinit must not be on same record'
9206  error stop '&ions: beam and dose must not be on same record'
9207  error stop '&ions: beam and energy must not be on same record'
9208  error stop '&ions: beam and fwhm must not be on same record'
9209  error stop '&ions: beam and mass must not be on same record'
9210  error stop '&ions: beam and rinit must not be on same record'
9211  error stop '&ions: beam and rotate must not be on same record'
9212  error stop '&ions: beam and tilt, tiltx or tilty must not be on same record'
9213  error stop '&ions: beam and xinit must not be on same record'
9214  error stop '&ions: beam and yinit must not be on same record'
9215  error stop '&ions: beam and zinit must not be on same record'
9216  error stop '&ions: beam and xwin must not be on same record'
9217  error stop '&ions: beam and ywin must not be on same record'
9218  error stop '&ions: beam and zwin must not be on same record'
9219  error stop '&ions: tilt and tilty must not be on same record'

930   error stop '&ions: filtra is obsolete. Use trafile instead'
931   error stop '&ions: nion must be specified on &setup record'

980   error stop '&ions: m>0 required'

   end subroutine readIon
!
!==============================================================================
!
   subroutine initIon( )
!             =======
!
!     Initialize the ion structure.
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: idummy
!
!  Start of Executable Code:
!  ------------------------
!
!  Check that has to be done before assigning default values
!
      if( ion%ranrot .and. defined( ion%rotate ) )  goto 953
!
!  Set default values and check obligatory parameters
!
      if( .not.defined( ion%doseunits ) )  ion%doseunits = 'cm-2'
      if( .not.defined( ion%name ) )  goto 957

      where( .not.defined( ion%moddiv(:) ) )  ion%moddiv(:) = 1

      if( .not.defined( ion%dinit ) )  ion%dinit = 5.43
      where( .not.defined( ion%diverg(1,:) ) )  ion%diverg(1,:) = ZERO
      where( .not.defined( ion%diverg(2,:) ) )  ion%diverg(2,:) = ion%diverg(1,:)
      if( defined( ion%dose ) ) then
         if( ion%doseunits == 'cm-2' ) then
            ion%dose = ion%dose / 1e16
         else if( ion%doseunits == 'cm-1' ) then
            ion%dose = ion%dose / 1e8
         end if
      else
         ion%dose = ONE                         ! 1/A^2=1e16/cm2, 1/A=1e8/cm
      end if
      if( .not.defined( ion%e ) )  goto 958
      if( .not.defined( ion%rinit ) )  ion%rinit = ZERO
      if( .not.defined( ion%rotate ) )  ion%rotate = ZERO
      if( .not.defined( ion%tiltx ) )  ion%tiltx = ZERO
      if( .not.defined( ion%tilty ) )  ion%tilty = ZERO
      if( .not.defined( ion%sigma(1) ) )  ion%sigma(1) = ZERO
      if( .not.defined( ion%sigma(2) ) )  ion%sigma(2) = ion%sigma(1)
      where( .not.defined( ion%weight(:) ) )  ion%weight(:) = ONE
!
!  Check validity of the parameters
!
      if( ion%lrinit .and. ion%lpinit )  goto 972
      if( ion%lrinit .or. ion%lpinit ) then
         if( ndim /= 1 )  goto 973
         if( any( ion%diverg(:,:) /= ZERO ) )  goto 974
      end if
      if( ion%lpinit .and. ion%tilty == ZERO )  goto 975
      if( .not.ion%liinit .and. .not.  &
          ( ( ndim > 1 .and. any( geom%bndtypex(:) == ENUM%bndtype%free ) ) .or.  &
            ( ndim > 2 .and. any( geom%bndtypey(:) == ENUM%bndtype%free ) ) .or.  &
            any( geom%bndtypez(:) == ENUM%bndtype%free ) ) )  goto 976

      if( ion%usetra .and. .not. defined( ion%trafile ) )  goto 990
      if( ion%usetra .and. nthread > 1 )  goto 991

      if( ion%doseunits /= 'cm-2' .and. ion%doseunits /= 'cm-1' .and.  &
          ion%doseunits /= '1' )  goto 952

      if( any( ion%moddiv(:) < 1 .or. ion%moddiv(:) > 2 ) )  goto 9541

      if( any( ion%diverg(:,:) < ZERO ) )  goto 9540
      if( ion%dose <= ZERO )  goto 956
      if( ion%e <= ZERO )  goto 951
      if( any( ion%sigma(:) < ZERO ) )  goto 955
      if( any( ion%weight(:) < ZERO ) )  goto 959
      if( all( ion%weight(:) == ZERO ) )  goto 9591
      if( ion%xinit(1) > ion%xinit(2) )  goto 970
      if( ion%yinit(1) > ion%yinit(2) )  goto 971
      if( ion%zinit(1) > ion%zinit(2) )  goto 9711
!
!  Calculate derived quantities
!
      ! rotate may later be modified if rotx/=0 in &xtal

      if( usebk .and. useion ) then
         call openFile( 22, bkfile, 'unformatted', 'read' )
!             ========  i     i           i          i
         !read( 22, err=930, end=930 ) dummy, ion%dose1
         rewind 22
         read( 22 ) idummy, ion%dose1
         if( ion%dose - ion%dose1 < - 1e-10*ion%dose )  &!goto 921
               error stop '&ions: dose must be larger than dose on backup file'
         close( 22 )
      else
         ion%dose1 = ZERO
      end if

      if( nion > nion1 ) then
         ion%ddose = ( ion%dose - ion%dose1 ) / ( nion - nion1 )
      else                                      ! nion=nion1 possible if read from backup file
         ion%ddose = ZERO                       !    but then ion%ddose is not needed
      end if

      ion%weight(:) = ion%weight(:) / sum( ion%weight(:) )

      if( .not.defined( ion%xinit(1) ) )  ion%xinit(1) = ZERO
      if( .not.defined( ion%xinit(2) ) )  ion%xinit(2) = ion%xinit(1)
      if( .not.defined( ion%yinit(1) ) )  ion%yinit(1) = ZERO
      if( .not.defined( ion%yinit(2) ) )  ion%yinit(2) = ion%yinit(1)
      if( .not.defined( ion%zinit(1) ) )  ion%zinit(1) = ZERO
      if( .not.defined( ion%zinit(2) ) )  ion%zinit(2) = ion%zinit(1)

      ion%dxinit = ion%xinit(2)-ion%xinit(1)
      ion%dyinit = ion%yinit(2)-ion%yinit(1)
      ion%dzinit = ion%zinit(2)-ion%zinit(1)

      ion%m = sum( ion%natom(1:nelemIon) * elem(1:nelemIon)%m )

      return

!921   error stop '&ions: dose must be larger than dose on backup file'
!930   error stop '&ions: error during reading backup file'

951   error stop '&ions: energy>0 required'
952   error stop '&ions: invalid doseunits'
953   error stop '&ions: ranrot and rotate must not be both specified'
9540  error stop '&ions: diverg>=0 required'
9541  error stop '&ions: 1<=moddiv<=2 required'
955   error stop '&ions: fwhm>=0 required'
956   error stop '&ions: dose>0 required'
957   error stop '&ions: name must be specified'
958   error stop '&ions: energy must be specified'
959   error stop '&ions: weight>=0 required'
9591  error stop '&ions: at least one weight>0 required'

970   error stop '&ions: xinit(2)>=xinit(1) required'
971   error stop '&ions: yinit(2)>=yinit(1) required'
9711  error stop '&ions: zinit(2)>=zinit(1) required'
972   error stop '&ions: lrinit=lpinit=t is not allowed'
973   error stop '&ions: ndim=1 required for lrinit=t or lpinit=t'
974   error stop '&ions: diverg=0 required for lrinit=t or lpinit=t'
975   error stop '&ions: lpinit=t is not allowed for tilt=0'
976   error stop '&ions: liinit=f requires at least one free boundary'

990   error stop '&setup: trafile required if usetra=t'
991   error stop '&setup: usetra=t not allowed for parallel execution'

   end subroutine initIon
!
!=======================================================================
!
   subroutine writeIon( )
!             ========
!
!     Write data derived from the &ions record to the output file
!
!  Arguments:
!
      integer(INT4) :: i, ielem, ielemmin, ielemmax, ibeam
!
!  Local Variables:
!
      character(9) :: mdname(2) = [ 'isotropic', 'Gaussian ' ]

      integer(INT4) :: nequal

      !real(REAL8) :: rotate, tilt
!
!  Start of Executable Code:
!  ------------------------
!
!  Open the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i

      if( ion%usetra ) then
!
!  Cannot report initial conditions if they are read from the file
!
         write( 21, '(t2,a,t13,a)' )  &
               '&ions:',  &
               'Initial conditions are read from the file'
         write( 21, '(t16,a)' )  trim( ion%trafile )
         if( ion%irxtal /= 0 ) then
            if( ion%llinit ) then
               write( 21, '(t13,a)' )  &
                     'llinit=t: Initial positions are moved to the nearest lattice site'
            else
               write( 21, '(t13,a)' )  &
                     'llinit=f: Initial positions are taken as provided'
            end if
         end if

      else
!
!  Ion species
!
         write( 21, '(t2,a,t13,a,a,a)' )  &
               '&ions:',  &
               'Ion species: ''' // str( ion%name ) // ''', atoms per ion:'

         do i = 0, (nelemIon-1)/8

            ielemmin = 8*i + 1
            ielemmax = min( 8*i + 8, nelemIon )

            write( 21, '(t13,a)', advance='no' ) '+--------+'
            write( 21, '(8a)' )  &
                  ('------+', ielem=ielemmin,ielemmax)

            write( 21, '(t13,a, t22,a)', advance='no' ) '|', '|'
            write( 21, '(8(a,i2,a))' )  &
                  ('Atom', ielem, '|', ielem=ielemmin,ielemmax)

            write( 21, '(t13,a, t22,a)', advance='no' ) '|', '|'
            write( 21, '(8(a,a2,a))' )  &
                  (' (', elem(ielem)%name, ') |', ielem=ielemmin,ielemmax)

            write( 21, '(t13,a)', advance='no' ) '+--------+'
            write( 21, '(8a)' )  &
                  ('------+', ielem=ielemmin,ielemmax)

            write( 21, '(t13,a)', advance='no' ) '|  Ion   |'
            write( 21, '(8(i4,a))' )  &
                  (ion%natom(ielem), '  |', ielem=ielemmin,ielemmax)

            write( 21, '(t13,a)', advance='no' ) '+--------+'
            write( 21, '(8a)' )  &
                  ('------+', ielem=ielemmin,ielemmax)
            write( 21, '(a)' )  ' '
         end do
!
!  Dose and energy
!
         if( ion%doseunits == 'cm-2' ) then
            write ( 21, '(t13,a)' )  &
               'dose=' // str( 1e16*ion%dose, 'es8.2e2' ) // ' cm-2'
         else if( ion%doseunits == 'cm-1' ) then
            write ( 21, '(t13,a)' )  &
               'dose=' // str( 1e8*ion%dose, 'es8.2e2' ) // ' cm-1'
         else
            write ( 21, '(t13,a)' )  &
               'dose=' // str( ion%dose, 'es8.2e2' )
         end if

         write( 21, '(t13,a)' )  &
               'energy=' // str( nint(ion%e) )
!
!  Initial directions
!
         if( ion%lrinit ) then

            write( 21, '(t13,a)' )  &
                   'tilty=' // str( ion%tilty, 'f99.2' ) //  &
                   'rotate parallel to vector from reference'
            write( 21, '(t13,a)' )  &
                   '   starting point to starting point'

         else if( ion%lpinit ) then

            write( 21, '(t13,a)' )  &
                  'tilty=' // str( ion%tilty, 'f99.2' ) //  &
                  ', rotate=' // str( ion%rotate, 'f99.2' )

         else

            write( 21, '(t13,a)', advance='no' )  &
                  'tiltx=' // str( ion%tiltx, 'f99.2' )
            write( 21, '(a)', advance='no' )  &
                  ', tilty=' // str( ion%tilty, 'f99.2' )
            if( ion%ranrot ) then
               write( 21, '(a)', advance='no' )  &
                     ', rotate=random'
            else
               write( 21, '(a)', advance='no' )  &
                     ', rotate=' // str( ion%rotate, 'f99.2' )
            end if

            if( ion%nbeam == 1 ) then
               write( 21, '(a)', advance='no' )  &
                     ', diverg=' // str( ion%diverg(1,1), 'f99.2' )
               if( ion%diverg(2,1) /= ion%diverg(1,1) )  &
                     write( 21, '(a,f6.2)', advance='no' )  &
                           ',' // str( ion%diverg(2,1), 'f99.2' )
               if( any( ion%diverg(:,1) > ZERO ) ) then
                  write( 21, '(a)' )  &
                        ', moddiv=' // mdname(ion%moddiv(1))
               else
                  write( 21, '(a)' )  ''
               end if
            else
               write( 21, '(a)' )
               do ibeam = 1, ion%nbeam
                  write( 21, '(t13,a)', advance='no' )  &
                        'beam=' // str( ibeam, 'i99' ) // ': '
                  write( 21, '(a)', advance='no' )  &
                        'weight=' // str( ion%weight(ibeam), 'f99.4' )
                  write( 21, '(a)', advance='no' )  &
                        ', diverg=' // str( ion%diverg(1,ibeam), 'f99.2' )
                  if( ion%diverg(2,ibeam) /= ion%diverg(1,ibeam) )  &
                        write( 21, '(a,f6.2)', advance='no' )  &
                              ',' // str( ion%diverg(2,ibeam), 'f99.2' )
                  if( any( ion%diverg(:,ibeam) > ZERO ) ) then
                        write( 21, '(t13,a)' )  &
                              ', moddiv=' // mdname(ion%moddiv(ibeam))
                  else
                     write( 21, '(a)' )  ''
                  end if
               end do
            end if

         end if
!
!  Initial positions
!
         nequal = count( [ ion%xinit(1) == ion%xinit(2),  &
                           ion%yinit(1) == ion%yinit(2),  &
                           ion%zinit(1) == ion%zinit(2) ] )

         if( nequal == 3 ) then
            write( 21, '(/t13,a)' )  &
                  'Reference starting point: (' // str( ion%xinit(1), 'f99.3' ) //  &
                  ',' // str( ion%yinit(1), 'f99.3' ) //  &
                  ',' // str( ion%zinit(1), 'f99.3' ) // ')'
         else
            write( 21, '(/t13,a)', advance='no' )  &
                  'The reference starting points are randomly chosen '
            select case( nequal )
               case( 2 )
                  write( 21, '(a)' )  &
                        'along the line'
               case( 1 )
                  write( 21, '(a)' )  &
                        'in the rectangle'
               case( 0 )
                  write( 21, '(a)' )  &
                        'in the cuboid'
            end select
            write( 21, '(t13,a)' )  &
                  '   spanned by ('//str( ion%xinit(1), 'f99.1' ) //  &
                  ',' // str( ion%yinit(1), 'f99.1' ) //  &
                  ',' // str( ion%zinit(1), 'f99.1' ) //  &
                  ') and (' // str( ion%xinit(2), 'f99.1' ) //  &
                  ',' // str( ion%yinit(2), 'f99.1' ) //  &
                  ',' // str( ion%zinit(2), 'f99.1' ) // ')'
         end if

         if( ion%lrinit ) then

            write( 21, '(t13,a)' )  &
                  'lrinit=t: The starting points are randomly chosen along a circle'
            write( 21, '(t13,a)' )  &
                   '   parallel to the xy plane around the reference point.'
            write( 21, '(t13,a)' )  &
                   '   The circle radius is rinit=' // str( ion%rinit, 'f99.2' )

         else if( ion%lpinit ) then

            write( 21, '(t13,a)' )  &
                  'lpinit=t: The starting points are randomly chosen within an'
            write( 21, '(t13,a)' )  &
                  '   interval of length dinit='//str( ion%dinit , 'f99.3' ) //  &
                  'along a line which is'
            write( 21, '(t13,a)' )  &
                  '   a distance rinit='//str( ion%rinit, 'f99.3' ) //  &
                  'from the reference point,'
            write( 21, '(t13,a)' )  &
                  '   which is parallel to the xy plane, and which is at an angle of'
            write( 21, '(t13,a)' )  &
                  '   rotate (see below) with respect to the x axis'

         else

            if( any( ion%sigma(:) > ZERO ) ) then
               write( 21, '(t13,a)' )  &
                     'The starting points are Gaussian distributed around the reference point'
               write( 21, '(t13,a)', advance='no' )  &
                     '   with a standard deviation(s) of ' // str( ion%sigma(1), 'f99.2' )
               if( ion%sigma(2) /= ion%sigma(1) ) then
                  write( 21, '(a)' )  &
                        ', ' // str( ion%sigma(2), 'f99.2' )
               else
                  write( 21, '(a)' ) ''
               end if
            end if

            if( ion%irxtal /= 0 ) then
               if( ion%llinit ) then
                  write( 21, '(t13,a)' )  &
                        'The starting points are moved to the nearest lattice site'
               else if( ion%ldet ) then
                  write( 21, '(t13,a)' )  &
                        'The starting points are randomly moved perpendicular to the'
                  write( 21, '(t13,a)' )  &
                        '   column used for deterministic defects'
               else
                  write( 21, '(t13,a)' )  &
                        'The starting point is moved to a random position within a unit'
                  write( 21, '(t13,a)' )  &
                        '   cell centered around the reference point'
               end if
            end if

         end if
      end if

      if( ion%liinit ) then
         write( 21, '(t13,a)' )  &
               'liinit=t: The ions may start inside the target'
      else
         write( 21, '(t13,a)' )  &
               'liinit=f: The starting points are forward or backward projected'
         write( 21, '(t13,a)' )  &
               '   to the surface along the initial direction'
      end if

      if( ndim > 1 .or. .not.geom%lcart ) then
         if( ion%reenter ) then
            write( 21, '(t13,a)' )  &
                 'Ions are followed through vacuum regions'
         else
            write( 21, '(t13,a)' )  &
                 'Ions are not followed through vacuum regions'
         end if
      end if

      write( 21, '()' )

      close( 21 )

   end subroutine writeIon
!
!=======================================================================
!
end module Param_Ion_mod

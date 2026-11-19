module Param_Sn_mod
!
!  Holds data on nuclear scattering and free flight path
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, ONE, TEN, HALF, THIRD, PI, BIG

   use String_mod,         only: str, capitalize, to_lower
   use Undefine_mod,       only: undefine, defined
   use Open_mod,           only: openFile
   use Param_WriteTab_mod, only: writeTab

   use Param_Setup_mod,    only: nelem, nr, lsput, usecell
   use Param_Geom_mod,     only: geom
   use Param_Elem_mod,     only: elem
   use Param_Ion_mod,      only: ion
   use Param_Xtal_mod,     only: xtal
   use Param_Mat_mod,      only: mat
   use Param_Dam_mod,      only: dam
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
   public :: readSn                 ! read &snpar record from input file
!            ======                 !
   public :: initSn                 ! initialize the sn structure
!            ======                 !
   public :: writeSn                ! write data derived from &snpar record
!            =======                !    to the output file
!
!  Public Data:
!  ===========
!
   type :: Sn_type
!          =======
!  Scattering parameters
!
      logical :: chkffpcyl                   ! check ffp cylinder flag
      logical :: lapsisef                    ! apsis cutoff energy flag
      logical :: linel                       ! inelastic scattering flag
      logical :: lmd                         ! molecular dynamics flag
      logical :: lmehl                       ! Gauss-Mehler quadrature flag
      logical :: lnegx1                      ! enforce negative time integral shift flag
      logical :: lpmaxtab                    ! flag indicating a table of pmax values is
                                             !    to be written to the 'PMAX' file
      logical :: lposffp                     ! enforce positive ffp flag
      logical :: lpsif                       ! trajectory cutoff angle flag
      logical, allocatable :: lquad(:,:)     ! (nelem,nelem) numerical quadrature flag
      logical :: lstffp                      ! statistical ffp flag

      character(8), allocatable :: quad(:,:)          ! (nelem,nelem) method of scattering angle
                                                      !    quadrature
      character(8), allocatable :: quadtime(:,:)      ! (nelem,nelem) method of time integral
                                                      !    quadrature
      character(80), allocatable :: scatfile(:,:)     ! (nelem,nelem) scattering table files
      character(80), allocatable :: scrcoeffile(:,:)  ! (nelem,nelem) screening function
                                                      !    coefficients table files
      character(256), allocatable :: scrfile(:,:)      ! (nelem,nelem) screening function table
                                                      !    files
      character(11), allocatable :: scrfun(:,:)       ! (nelem,nelem) screening functions
      character(8), allocatable :: scrlen(:,:)        ! (nelem,nelem) screening lengths

      integer(INT4), allocatable :: nabsc(:,:)        ! (nelem,nelem) number of abscissas for
                                                      !    scattering angle quadrature
      integer(INT4), allocatable :: nabsctime(:,:)    ! (nelem,nelem) number of abscissas for
                                                      !    time integral quadrature
      integer(INT4) :: nmehl                 ! number of abscissas for Gauss-Mehler quadrature
      integer(INT4) :: reflectmax            ! maximum number of surface reflections

      real(REAL8) :: cospsf                  ! sin(trajectory cutoff angle)   !!! rename: sinpsf
      real(REAL8) :: denffp                  ! minimum energy transfer in amorphous material
                                             !    [eV]
      real(REAL8) :: dirl2f                  ! sin(trajectory cutoff angle)^2 !!! eliminate ?
      real(REAL8) :: ef                      ! cutoff energy for ion trajectories [eV]
      real(REAL8) :: eminpsi                 ! minimum energy for the scattering angle 
                                             !    criterion [eV]
      real(REAL8) :: eminpsi2                ! minimum energy for the scattering angle 
                                             !    criterion [eV] outside target
      real(REAL8) :: etime                   ! energy below which the time integral is
                                             !    considered [eV]
      real(REAL8) :: ffpmax                  ! maximum free flight path [A]
                                             !    (default value will not be set for usesdf=t
                                             !    until initTarget)
      real(REAL8) :: pmax                    ! maximum of any used pmaxmax, pmaxc(:)
      real(REAL8), allocatable :: pmaxc(:)   ! (nr) maximum impact parameter in crystalline
                                             !    material [A]
      real(REAL8) :: pmaxmax                 ! maximum maximum impact parameter [A]
      real(REAL8) :: pmaxmin                 ! minimum maximum impact parameter [A]
      real(REAL8) :: pmaxpsi                 ! maximum maximum impact parameter for the
                                             !    scattering angle criterion [A]
      real(REAL8) :: psiffp                  ! minimum scattering angle in amorphous materials
                                             !    [deg]
      real(REAL8) :: psiffp2                 ! minimum scattering angle outside target close
                                             !    to amorphous material [deg]
   end type Sn_type

   type(Sn_type), public :: sn                       ! sn structure
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine readSn( )
!             ======
!
!     Read the &snpar record from the input file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: atom1                    ! atom index of projectile
      integer(INT4) :: atom2                    ! atom index of target atom

      logical :: chkffpcyl                      ! check ffp cylinder flag
      logical :: lapsisef                       ! apsis cutoff energy flag
      logical :: linel                          ! inelastic scattering flag
      logical :: lmd                            ! molecular dynamics flag
      logical :: lnegx1                         ! enforce negative time integral shift flag
      logical :: lmehl                          ! Gauss-Mehler quadrature flag
      logical :: lposffp                        ! enforce positive ffp flag
      logical :: lpmaxtab                       ! pmax write table file
      logical :: lstffp                         ! statistical ffp flag

      character(80) :: coeffil                  ! obsolete
      character(80) :: coeffile                 ! obsolete
      character(80) :: scrcoeffile              ! name of screening function coefficient file
      character(8) :: quad                      ! method of scattering angle quadrature
      character(8) :: quadtime                  ! method of time integral quadrature
      character(80) :: scafil                   ! obsolete
      character(80) :: scatfile                 ! name of scattering table file
      character(80) :: scrfil                   ! obsolete
      character(256) :: scrfile                  ! name of screening function table file
      character(11) :: scrfun                   ! screening function
      character(8) :: scrlen                    ! screening length

      integer(INT4) :: nabsc                    ! number of abscissas for scattering angle
                                                !    quadrature
      integer(INT4) :: nabsctime                ! number of abscissas for time integral
                                                !    quadrature
      integer(INT4) :: nmehl                    ! Number of abscissas in Gauss-Mehler
                                                !    quadrature
      integer(INT4) :: reflectmax               ! Maximum number of surface reflections

      real(REAL8) :: denffp                     ! minimum energy transfer [eV]
      real(REAL8) :: ef                         ! cutoff energy for ion trajectory
      real(REAL8) :: eminpsi                    ! lower bound to energy for application of
                                                !    scattering angle criterion
      real(REAL8) :: eminpsi2                   ! lower bound to energy for application of
                                                !    scattering angle criterion outside target
      real(REAL8) :: etime                      ! maximum energy for time integral [keV]
      real(REAL8) :: ffpmax                     ! maximum free flight path [A]
      real(REAL8) :: pmaxmax                    ! maximum maximum impact parameter in
                                                !    amorphous regions [A]
      real(REAL8) :: pmaxmin                    ! minimum maximum impact parameter in
                                                !    amorphous regions [A]
                                                !    calculated from scattering angle criterion
      real(REAL8) :: pmaxpsi                    ! maximum maximum impact parameter for the
                                                !    scattering angle criterion [A]
      real(REAL8) :: psif                       ! cutoff angle for ion trajectory
      real(REAL8) :: psiffp                     ! minimum scattering angle [deg]
      real(REAL8) :: psiffp2                    ! minimum scattering angle outside target [deg]

      namelist /snpar/ atom1, atom2,  &
                       lapsisef, chkffpcyl, linel, lmd, lnegx1, lmehl, lposffp,  &
                       lpmaxtab, lstffp,  &
                       quad, quadtime, scatfile, scrcoeffile, scrfile, scrfun, scrlen,  &
                       nabsc, nabsctime, nmehl, reflectmax,  &
                       denffp, ef, eminpsi, eminpsi2, etime, ffpmax,  &
                       pmaxmax, pmaxmin, pmaxpsi, psif, psiffp, psiffp2
      ! obsolete
      namelist /snpar/ coeffil, coeffile, scafil, scrfil

      character(256) :: msg
      integer(INT4) :: ios, ielem1, ielem1min, ielem1max, n
      integer(INT4) :: ielem2, ielem2min, ielem2max
!
!  Start of Executable Code:
!  ------------------------
!
!  Allocate memory
!
      allocate( sn%nabsc( nelem, nelem ) )
      allocate( sn%nabsctime( nelem, nelem ) )

      if( .not.usecell )  allocate( sn%pmaxc( nr ) )

      allocate( sn%quad( nelem, nelem ) )
      allocate( sn%quadtime( nelem, nelem ) )
      allocate( sn%scatfile( nelem, nelem ) )
      allocate( sn%scrcoeffile( nelem, nelem ) )
      allocate( sn%scrfile( nelem, nelem ) )
      allocate( sn%scrfun( nelem, nelem ) )
      allocate( sn%scrlen( nelem, nelem ) )
!
!  Set default values of structure variables
!
      sn%chkffpcyl = .false.
      sn%lapsisef = .false.
      sn%linel = .false.
      sn%lmd = .false.
      sn%lnegx1 = .true.
      sn%lmehl = .false.
      sn%lposffp = .false.
      sn%lpmaxtab = .false.
      sn%lstffp = .true.

      call undefine( sn%quad(:,:) )
      call undefine( sn%quadtime(:,:) )
      call undefine( sn%scatfile(:,:) )
      call undefine( sn%scrcoeffile(:,:) )
      call undefine( sn%scrfile(:,:) )
      call undefine( sn%scrfun(:,:) )
      call undefine( sn%scrlen(:,:) )

      call undefine( sn%nabsc )
      call undefine( sn%nabsctime )
      call undefine( sn%nmehl )
      call undefine( sn%reflectmax )

      call undefine( sn%cospsf )
      call undefine( sn%denffp )
      call undefine( sn%ef )
      call undefine( sn%eminpsi )
      call undefine( sn%eminpsi2 )
      call undefine( sn%etime )
      call undefine( sn%ffpmax )
      call undefine( sn%pmaxmax )
      call undefine( sn%pmaxmin )
      call undefine( sn%pmaxpsi )
      call undefine( sn%psiffp )
      call undefine( sn%psiffp2 )
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
         call undefine( atom1 )
         call undefine( atom2 )

         chkffpcyl = sn%chkffpcyl
         lapsisef = sn%lapsisef
         linel = sn%linel
         lmd = sn%lmd
         lnegx1 = sn%lnegx1
         lmehl = sn%lmehl
         lposffp = sn%lposffp
         lpmaxtab = sn%lpmaxtab
         lstffp = sn%lstffp

         call undefine( coeffil )               ! obsolete
         call undefine( coeffile )              ! obsolete
         call undefine( quad )
         call undefine( quadtime )
         call undefine( scafil )                ! obsolete
         call undefine( scatfile )
         call undefine( scrcoeffile )
         call undefine( scrfil )                ! obsolete
         call undefine( scrfile )
         call undefine( scrfun )
         call undefine( scrlen )

         call undefine( nabsc )
         call undefine( nabsctime )
         call undefine( nmehl )
         call undefine( reflectmax )

         call undefine( denffp )
         call undefine( ef )
         call undefine( eminpsi )
         call undefine( eminpsi2 )
         call undefine( etime )
         call undefine( ffpmax )
         call undefine( pmaxmax )
         call undefine( pmaxmin )
         call undefine( pmaxpsi )
         call undefine( psif )
         call undefine( psiffp )
         call undefine( psiffp2 )
!
!  Read &snpar record from the input file
!
         read( 11, nml=snpar, iostat=ios, iomsg=msg )
         if( ios < 0  )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 900
         end if
!
!  Check obsolete parameters
!
         if( defined( scafil ) )  goto 930
         if( defined( scrfil ) )  goto 931
         if( defined( coeffil ) )  goto 932
         if( defined( coeffile ) )  goto 933
!
!  Check index parameters
!
         if( defined( atom1 ) ) then
            if( atom1 < 1 .or. atom1 > nelem )  goto 906
         end if

         if( defined( atom2 ) ) then
            if( atom2 < 1 .or. atom2 > nelem )  goto 907
         end if
!
!  Check parameters on consistency with index parameters
!
         if( defined( atom1 ) .or. defined( atom2 ) ) then
            if( defined( nmehl ) )  goto 920
            if( defined( reflectmax ) )  goto 9201
            if( defined( etime ) )  goto 921
            if( defined( psiffp ) )  goto 922
            if( defined( psiffp2 ) )  goto 922
            if( defined( denffp ) )  goto 923
            if( defined( pmaxmax ) )  goto 9240
            if( defined( pmaxpsi ) )  goto 92410
            if( defined( eminpsi ) )  goto 92411
            if( defined( eminpsi2 ) )  goto 92412
            if( defined( pmaxmin ) )  goto 9242
            if( defined( ffpmax ) )  goto 9243
            if( defined( ef ) )  goto 925
            if( defined( psif ) )  goto 926
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
         ielem2 = ielem2min
!
!  Make sure no more than one of scatfile, scrcoeffile, scrfile, and scrfun is specified
!
         n = 0
         if( defined( scatfile ) )  n = n + 1
         if( defined( scrcoeffile ) )  n = n + 1
         if( defined( scrfile ) )  n = n + 1
         if( defined( scrfun ) )  n = n + 1
         if( n > 1 )  goto 940
!
!  Copy parameters to sn structure
!
         sn%chkffpcyl = chkffpcyl
         sn%lapsisef = lapsisef
         sn%linel = linel
         sn%lmd = lmd
         sn%lnegx1 = lnegx1
         sn%lmehl = lmehl
         sn%lposffp = lposffp
         sn%lpmaxtab = lpmaxtab
         sn%lstffp = lstffp

         if( defined( quad ) )  sn%quad(ielem2min:ielem2max,ielem1min:ielem1max) = quad
         if( defined( quadtime ) )  &
               sn%quadtime(ielem2min:ielem2max,ielem1min:ielem1max) = quadtime

         if( defined( scatfile ) ) then
            sn%scatfile(ielem2min:ielem2max,ielem1min:ielem1max) = scatfile
            call undefine( sn%scrcoeffile(ielem2min:ielem2max,ielem1min:ielem1max) )
            call undefine( sn%scrfile(ielem2min:ielem2max,ielem1min:ielem1max) )
            call undefine( sn%scrfun(ielem2min:ielem2max,ielem1min:ielem1max) )
         else if( defined( scrcoeffile ) ) then
            sn%scrcoeffile(ielem2min:ielem2max,ielem1min:ielem1max) = scrcoeffile
            call undefine( sn%scatfile(ielem2min:ielem2max,ielem1min:ielem1max) )
            call undefine( sn%scrfile(ielem2min:ielem2max,ielem1min:ielem1max) )
            call undefine( sn%scrfun(ielem2min:ielem2max,ielem1min:ielem1max) )
         else if( defined( scrfile ) ) then
            sn%scrfile(ielem2min:ielem2max,ielem1min:ielem1max) = scrfile
            call undefine( sn%scatfile(ielem2min:ielem2max,ielem1min:ielem1max) )
            call undefine( sn%scrcoeffile(ielem2min:ielem2max,ielem1min:ielem1max) )
            call undefine( sn%scrfun(ielem2min:ielem2max,ielem1min:ielem1max) )
         else if( defined( scrfun ) ) then
            sn%scrfun(ielem2min:ielem2max,ielem1min:ielem1max) = scrfun
            call undefine( sn%scrcoeffile(ielem2min:ielem2max,ielem1min:ielem1max) )
            call undefine( sn%scatfile(ielem2min:ielem2max,ielem1min:ielem1max) )
            call undefine( sn%scrfile(ielem2min:ielem2max,ielem1min:ielem1max) )
         end if

         if( defined( scrlen ) )  sn%scrlen(ielem2min:ielem2max,ielem1min:ielem1max) = scrlen

         if( defined( nabsc ) )  sn%nabsc(ielem2min:ielem2max,ielem1min:ielem1max) = nabsc
         if( defined( nmehl ) )  sn%nmehl = nmehl
         if( defined( nabsctime ) ) then
            sn%nabsctime(ielem2min:ielem2max,ielem1min:ielem1max) = nabsctime
            error stop 'nabsctime not yet implemented'
         end if

         if( defined( denffp ) )  sn%denffp = denffp
         if( defined( reflectmax ) )  sn%reflectmax = reflectmax

         if( defined( psif ) )  sn%cospsf = sin( psif*PI/180 )
         if( defined( ef ) )  sn%ef = ef
         if( defined( eminpsi ) )  sn%eminpsi = eminpsi
         if( defined( eminpsi2 ) )  sn%eminpsi2 = eminpsi2
         if( defined( etime ) )  sn%etime = 1000*etime
         if( defined( ffpmax ) )  sn%ffpmax = ffpmax
         if( defined( pmaxmax ) )  sn%pmaxmax = pmaxmax
         if( defined( pmaxmin ) )  sn%pmaxmin = pmaxmin
         if( defined( pmaxpsi ) )  sn%pmaxpsi = pmaxpsi
         if( defined( psiffp ) )  sn%psiffp = psiffp
         if( defined( psiffp2 ) )  sn%psiffp2 = psiffp2
!
!  Some checks that cannot be done in initSn
!
         if( defined( psif ) ) then
            if( psif < 0 .or. psif > 90 )  goto 962
            if( (ion%lrinit .or. ion%lpinit) .and. psif == ZERO )  goto 963
         end if

      end do
!
!  End of loop over reading records
!  --------------------------------
!
      close( 11 )

      return
!
!  Error messages
!
900   error stop '&snpar: error during reading namelist'

904   error stop '&snpar: region not allowed with usecell=T'
905   error stop '&snpar: 1<=region<=nr required'
906   error stop '&snpar: 1<=atom1<=natom required'
907   error stop '&snpar: 1<=atom2<=natom required'

920   error stop '&snpar: atom1 or atom2 and nmehl must not be on same record'
9201  error stop '&snpar: atom1 or atom2 and reflectmax must not be on same record'
921   error stop '&snpar: atom1 or atom2 and etime must not be on same record'
922   error stop '&snpar: atom1 or atom2 and psiffp or psiffp2 must not be on same record'
923   error stop '&snpar: atom1 or atom2 and denffp must not be on same rec.'
9240  error stop '&snpar: atom1 or atom2 and pmaxmax must not be on same record'
92410 error stop '&snpar: atom1 or atom2 and pmaxpsi must not be on same record'
92411 error stop '&snpar: atom1 or atom2 and eminpsi must not be on same record'
92412 error stop '&snpar: atom1 or atom2 and eminpsi2 must not be on same record'
9242  error stop '&snpar: atom1 or atom2 and pmaxmin must not be on same record'
9243  error stop '&snpar: atom1 or atom2 and ffpmax must not be on same record'
925   error stop '&snpar: atom1 or atom2 and ef must not be on same record'
926   error stop '&snpar: atom1 or atom2 and psif must not be on same record'

930   error stop '&snpar: scafil is obsolete. Use scatfile instead.'
931   error stop '&snpar: scrfil is obsolete. Use scrfile instead.'
932   error stop '&snpar: coeffil is obsolete. Use scrcoeffile instead.'
933   error stop '&snpar: coeffile is obsolete. Use scrcoeffile instead.'

940   error stop '&snpar: only one of scatfile, scrcoeffile, scrfile, scrfun must be specified'
962   error stop '&snpar: psif=0 or 0<psif<180 required'
963   error stop '&snpar: psif=0 is not allowed for lrinit=t or lpinit=t'

   end subroutine readSn
!
!=======================================================================
!
   subroutine initSn( )
!             =======
!
!     Initialize the sn structure
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(256) :: line
      integer(INT4) :: i, ielem1, ielem2, ielemmin, ielemmax, ir
!
!  Start of Executable Code:
!  ------------------------
!
!  Check that has to be done before assigning default values
!
      if( defined( sn%pmaxpsi ) .and. defined( sn%eminpsi ) )  goto 971
!
!  Set default values (except scrlen)
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            if( .not.defined( sn%quad(ielem2,ielem1) ) )  &
                  sn%quad(ielem2,ielem1) = 'Legendre'
            if( .not.defined( sn%quadtime(ielem2,ielem1) ) )  &
                  sn%quadtime(ielem2,ielem1) = 'no'

            if( .not.defined( sn%scatfile(ielem2,ielem1) )  &
                .and. .not.defined( sn%scrfile(ielem2,ielem1) )  &
                .and. .not.defined( sn%scrcoeffile(ielem2,ielem1) )  &
                .and. .not.defined( sn%scrfun(ielem2,ielem1) ) ) then
               if( (elem(ielem1)%z == 5 .and. elem(ielem2)%z == 8) .or.  &
                   (elem(ielem1)%z == 8 .and. elem(ielem2)%z == 5) ) then
                  sn%scatfile(ielem2,ielem1) = 'BO'
               else if( (elem(ielem1)%z == 5 .and. elem(ielem2)%z == 14) .or.  &
                        (elem(ielem1)%z == 14 .and. elem(ielem2)%z == 5) ) then
                  sn%scatfile(ielem2,ielem1) = 'BSI'
               else if( elem(ielem1)%z == 14 .and. elem(ielem2)%z == 14 ) then
                  sn%scatfile(ielem2,ielem1) = 'SISI'
               else if( (elem(ielem1)%z == 15 .and. elem(ielem2)%z == 14) .or.  &
                        (elem(ielem1)%z == 14 .and. elem(ielem2)%z == 15) ) then
                  sn%scatfile(ielem2,ielem1) = 'PSI'
               else if( (elem(ielem1)%z == 33 .and. elem(ielem2)%z == 14) .or.  &
                        (elem(ielem1)%z == 14 .and. elem(ielem2)%z == 33) ) then
                  sn%scatfile(ielem2,ielem1) = 'ASSI'
               else if( (elem(ielem1)%z == 68 .and. elem(ielem2)%z == 14) .or.  &
                        (elem(ielem1)%z == 14 .and. elem(ielem2)%z == 68) ) then
                  sn%scatfile(ielem2,ielem1) = 'ERSI'
               else
                  sn%scatfile(ielem2,ielem1) = 'ZBLuniv'
               end if
            end if
         end do
      end do

      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            if( .not.defined( sn%nabsc(ielem2,ielem1) ) )  sn%nabsc(ielem2,ielem1) = 4
            if( .not.defined( sn%nabsctime(ielem2,ielem1) ) )  sn%nabsctime(ielem2,ielem1) = 4
         end do
      end do
      if( .not.defined( sn%nmehl ) )  sn%nmehl = 32
      if( .not.defined( sn%reflectmax ) )  sn%reflectmax = huge( sn%reflectmax )

      if( .not.defined( sn%cospsf ) )  sn%cospsf = ZERO
      if( .not.defined( sn%denffp ) )  sn%denffp = 15.
      if( .not.defined( sn%ef ) )  sn%ef = 10.
      if( .not.defined( sn%etime ) )  sn%etime = 5000.

      if( .not.defined( sn%ffpmax ) ) then
         if( geom%usehyperoctree ) then
            sn%ffpmax = 30.                     ! finite value since it determines cell distmin
         else
            sn%ffpmax = BIG
         end if
      end if

      if( .not.defined( sn%pmaxmax ) )  sn%pmaxmax = 4.                           ! Angstrom
      if( .not.defined( sn%pmaxmin ) )  sn%pmaxmin = 0.004                        ! Angstrom

      if( .not.defined( sn%pmaxpsi ) .and. .not.defined( sn%eminpsi ) )  &
            sn%pmaxpsi = 1.53 / min( ONE, sqrt( 0.15 + 0.85*ion%e/10000 ) )
      if( defined( sn%pmaxpsi ) )  &
            sn%pmaxpsi = max( sn%pmaxmin, min( sn%pmaxmax, sn%pmaxpsi ) )
      if( lsput .and. .not.defined( sn%psiffp2 ) )  &
            sn%psiffp2 = 0.01
      if( .not.defined( sn%eminpsi2 ) ) then
         if( lsput ) then
            sn%eminpsi2 = HALF * ion%e
         else
            sn%eminpsi2 = sn%eminpsi
         end if
      end if
!
!     Note: sn%eminpsi, sn%eminpsi2, and sn%psiffp2 may remain undefined. This has special
!           meaning in Sn_FFP_mod
!
      if( .not.defined( sn%psiffp ) )  sn%psiffp = 1.
!
!  Check validity of the parameters and convert case of character-type parameters
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            if( sn%scatfile(ielem2,ielem1) == 'SCATTAB'  &
                .or. sn%scatfile(ielem2,ielem1) == 'SCATTAB_big' ) goto 930

            sn%quad(ielem2,ielem1) = to_lower( sn%quad(ielem2,ielem1) )
!                                    ========            i
            sn%quadtime(ielem2,ielem1) = to_lower( sn%quadtime(ielem2,ielem1) )
!                                        ========            i

            if( sn%quad(ielem2,ielem1) == 'legendre' ) then
               sn%quad(ielem2,ielem1) = 'Legendre'
            else if( sn%quad(ielem2,ielem1) == 'mehler' ) then
               sn%quad(ielem2,ielem1) = 'Mehler'
            else
               goto 940
            end if

            if( sn%quadtime(ielem2,ielem1) == 'no' ) then
               sn%quadtime(ielem2,ielem1) = 'no'
            else if( sn%quadtime(ielem2,ielem1) == 'legendre' ) then
               sn%quadtime(ielem2,ielem1) = 'Legendre'
            else if( sn%quadtime(ielem2,ielem1) == 'mehler' ) then
               sn%quadtime(ielem2,ielem1) = 'Mehler'
            else
               goto 941
            end if
            if( sn%quadtime(ielem2,ielem1) /= 'no' )  goto 942

            if( defined( sn%scrfun(ielem2,ielem1) ) ) then
               sn%scrfun(ielem2,ielem1) = to_lower( sn%scrfun(ielem2,ielem1) )
!                                         ========            i
               if( sn%scrfun(ielem2,ielem1) == 'zbluniv' ) then
                  sn%scrfun(ielem2,ielem1) = 'ZBLuniv'
               else if( sn%scrfun(ielem2,ielem1) == 'krc' ) then
                  sn%scrfun(ielem2,ielem1) = 'KrC'
               else if( sn%scrfun(ielem2,ielem1) == 'moliere' ) then
                  sn%scrfun(ielem2,ielem1) = 'Moliere'
               else if( sn%scrfun(ielem2,ielem1) == 'zblspec' ) then
                  sn%scrfun(ielem2,ielem1) = 'ZBLspec'
               else if( sn%scrfun(ielem2,ielem1) == 'nlh' ) then
                  sn%scrfun(ielem2,ielem1) = 'NLH'
               else if( sn%scrfun(ielem2,ielem1) == 'nlhlin' ) then
                  sn%scrfun(ielem2,ielem1) = 'NLHlin'
               else
                  goto 945
               end if
            end if
         end do
      end do
      if( any( sn%quad(:,:) == 'Legendre' .and. sn%nabsc(:,:) > 64 ) )  goto 9431
      if( any( sn%nabsc(:,:) <= 0 ) )  goto 943
      if( any( sn%quadtime(:,:) == 'Legendre' .and. sn%nabsctime(:,:) > 64 ) )  goto 9441
      if( any( sn%nabsctime(:,:) <= 0 ) )  goto 944

      if( sn%nmehl <= 0 )  goto 950
      if( sn%reflectmax < 0 )  goto 9501

      if( sn%denffp <= ZERO )  goto 953
      if( sn%lmd .and. dam%ldam )  goto 958
      if( dam%lflow .and. (dam%eminflow < ZERO .or. dam%eminflow > sn%ef) )  goto 960
      if( dam%lflow ) then
         if( dam%eminflow > ion%e )  goto 9601
      else
         if( sn%ef > ion%e )  goto 9602
      end if
      if( sn%ef <= ZERO )  goto 961
      if( defined( sn%eminpsi ) .and. sn%eminpsi <= ZERO )  goto 9544
      if( defined( sn%eminpsi2 ) .and. sn%eminpsi2 <= ZERO )  goto 9544
      if( sn%etime < ZERO )  goto 9510
      if( defined( sn%ffpmax ) .and. sn%ffpmax <= ZERO )  goto 9543

      if( sn%pmaxmax <= ZERO )  goto 9540
      if( sn%pmaxpsi <= ZERO )  goto 9541
      if( sn%psiffp <= ZERO .or. sn%psiffp >= 180 )  goto 952
      if( defined(sn%psiffp2 ) .and. ( sn%psiffp2 <= ZERO .or. sn%psiffp2 >= 180 ) )  goto 952
      if( sn%pmaxmin <= ZERO )  goto 9542
      if( sn%pmaxmax < sn%pmaxmin )  goto 970
!
!  Replace ZBLspec, NLH, and NLHlin screening functions by screening coefficients file
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            if( defined( sn%scrfun(ielem2,ielem1) ) ) then
               if( sn%scrfun(ielem2,ielem1) == 'ZBLspec' ) then
                  call undefine( sn%scrfun(ielem2,ielem1) )
                  sn%scrcoeffile(ielem2,ielem1) = 'ZBLspec_coeffs.dat'
               else if( sn%scrfun(ielem2,ielem1) == 'NLH' ) then
                  call undefine( sn%scrfun(ielem2,ielem1) )
                  sn%scrcoeffile(ielem2,ielem1) = 'NLH_coeffs.dat'
               else if( sn%scrfun(ielem2,ielem1) == 'NLHlin' ) then
                  call undefine( sn%scrfun(ielem2,ielem1) )
                  sn%scrcoeffile(ielem2,ielem1) = 'NLHlin_coeffs.dat'
               end if
            end if
         end do
      end do
!
!  Modify paths
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
            if( defined( sn%scatfile(ielem2,ielem1) ) ) then
               sn%scatfile(ielem2,ielem1) = 'scattab/' // sn%scatfile(ielem2,ielem1)
            else if( defined( sn%scrfile(ielem2,ielem1) ) ) then
               if( elem(ielem1)%z <= elem(ielem2)%z ) then
                  ielemmin = ielem1
                  ielemmax = ielem2
               else
                  ielemmin = ielem2
                  ielemmax = ielem1
               end if
               if( to_lower( sn%scrfile(ielem2,ielem1) ) == 'zblspec' ) then
!                  ========            i
                  sn%scrfile(ielem2,ielem1) = 'ZBLspec_pos/screening_zblspec_' //  &
                        str( elem(ielemmin)%z ) // '_' // str( elem(ielemmax)%z ) // '.dat'
               else if( to_lower( sn%scrfile(ielem2,ielem1) ) == 'zbluniv' ) then
!                       ========            i
                  sn%scrfile(ielem2,ielem1) = 'ZBLuniv/screening_zbluniv_' //  &
                        str( elem(ielemmin)%z ) // '_' // str( elem(ielemmax)%z ) // '.dat'
               else if( to_lower( sn%scrfile(ielem2,ielem1) ) == 'dmol' ) then
!                       ========            i
                  sn%scrfile(ielem2,ielem1) = 'DMol_pos/screening_dmol_' //  &
                        str( elem(ielemmin)%z ) // '_' // str( elem(ielemmax)%z ) // '.dat'
               end if
               sn%scrfile(ielem2,ielem1) = 'screentab/' // sn%scrfile(ielem2,ielem1)
            else if( defined( sn%scrcoeffile(ielem2,ielem1) ) ) then
               sn%scrcoeffile(ielem2,ielem1) = 'screencoefs/' // sn%scrcoeffile(ielem2,ielem1)
            end if
         end do
      end do
!
!  Handle scrlen
!
      do ielem1 = 1, nelem
         do ielem2 = 1, nelem
                                                ! Default values
            if( .not.defined( sn%scrlen(ielem2,ielem1) ) ) then
               if( defined( sn%scrcoeffile(ielem2,ielem1) ) ) then
                  sn%scrlen(ielem2,ielem1) = 'zbl'
                  call openFile( 11, sn%scrcoeffile(ielem2,ielem1), 'formatted', 'read' )
!                      ========  i              i                        i         i
                  do                           ! Read scrlen from scrcoeffile
                     read( 11, '(a)' )  line
                     if( line(1:1) /= '#' )  exit
                     i = index( line, 'scrlen' )
                     if( i == 0 )  cycle
                     line = adjustl( line(i+6:) )
                     if( line(1:1) /= '=' )  cycle
                     line = adjustl( line(2:) )
                     i = scan( line, ', ' )
                     sn%scrlen(ielem2,ielem1) = to_lower( line(:i-1) )
!                                               ========      i
                     exit
                  end do
                  close( 11 )
               else if( defined( sn%scrfile(ielem2,ielem1) ) ) then
                  sn%scrlen(ielem2,ielem1) = '1'
               else if( defined( sn%scrfun(ielem2,ielem1) ) ) then
                  if( to_lower( sn%scrfun(ielem2,ielem1) ) == 'krc'  &
!                     ========            i
                      .or. to_lower( sn%scrfun(ielem2,ielem1) ) == 'moliere' ) then
!                          ========            i
                     sn%scrlen(ielem2,ielem1) = 'firsov'
                  else
                     sn%scrlen(ielem2,ielem1) = 'zbl'
                  end if
               else
                  sn%scrlen(ielem2,ielem1) = 'zbl'
               end if
            end if
                                                ! Adjust case of screening length names
            if( sn%scrlen(ielem2,ielem1) /= '1' ) then
               sn%scrlen(ielem2,ielem1) = to_lower( sn%scrlen(ielem2,ielem1) )
               if( sn%scrlen(ielem2,ielem1) == 'zbl' ) then
                  sn%scrlen(ielem2,ielem1) = 'ZBL'
               else if( sn%scrlen(ielem2,ielem1) == 'firsov' ) then
                  sn%scrlen(ielem2,ielem1) = 'Firsov'
               else if( sn%scrlen(ielem2,ielem1) == 'lindhard' ) then
                  sn%scrlen(ielem2,ielem1) = 'Lindhard'
               else
                  print *, 'scrlen=', sn%scrlen(ielem2,ielem1)
                  goto 946
               end if
            end if
         end do
      end do
!
!  Calculate derived quantities
!
      sn%dirl2f = sn%cospsf**2
      sn%lpsif = (sn%cospsf /= ZERO)

      allocate( sn%lquad( nelem, nelem ) )
      sn%lquad(:,:) = .not. defined( sn%scatfile(:,:) )
      if( sn%lmehl )  sn%lquad(:,:) = .true.

      if( .not.usecell ) then
         do ir = 1, nr
            sn%pmaxc(ir) = xtal%pmaxc * mat(ir)%xscale
         end do
      end if
!
!     sn%pmax will be initialized in initFFP
!
      return

930   error stop '&snpar: Scattering table SCATTAB(_big) has been renamed to ZBLuniv(_big)'

940   error stop '&snpar: invalid quad'
941   error stop '&snpar: invalid quadtime'
942   error stop '&snpar: quadtime must be ''no'''
943   error stop '&snpar: nabsc>0 required'
9431  error stop '&snpar: 0<nabsc<65 required for quad=''Legendre'''
944   error stop '&snpar: nabsctime>0 required'
9441  error stop '&snpar: 0<nabsctime<65 required for quadtime=''Legendre'''
945   error stop '&snpar: invalid scrfun'
946   error stop '&snpar: invalid scrlen'

950   error stop '&snpar: nmehl>0 required'
9501  error stop '&snpar: reflectmax>=0 required'
9510  error stop '&snpar: etime>0 or =0 required'
952   error stop '&snpar: 0<psiffp<180 and 0<psiffp2<180 required'
953   error stop '&snpar: denffp>0 required'
9540  error stop '&snpar: pmaxmax>0 required'
9541  error stop '&snpar: pmaxpsi>0 required'
9542  error stop '&snpar: pmaxmin>0 required'
9543  error stop '&snpar: ffpmax>0 required'
9544  error stop '&snpar: eminpsi>0, eminpsi2>0 required'

958   error stop '&snpar: damage may not be recorded for md simulation'
960   error stop '&damage,&snpar: 0<=eminflow<=ef required'
9601  error stop '&damage,&snpar: e0>eminflow required'
9602  error stop '&snpar: e0>ef required'
961   error stop '&snpar: ef>0 required'

970   error stop '&snpar: pmaxmin <= pmaxmax required'
971   error stop '&snpar: pmaxpsi and eminpsi must not be on same record'

   end subroutine initSn
!
!=======================================================================
!
   subroutine writeSn( )
!             =======
!
!     Write parameters derived from the &snpar record to the output file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(80) :: text
      character(80), allocatable :: carr(:,:)
      integer(INT4) :: ielem1, ielem2
!
!  Start of Executable Code:
!  ------------------------
!
      allocate( carr( nelem, nelem ) )
!
!  Open the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i

      write( 21, '(t2,a)', advance='no' ) '&snpar:    '
!
!  Simulation mode
!
      if( sn%lmd ) then
         write( 21, '(a)' )  &
               'Molecular dynamics simulation'
      else
         write( 21, '(a)' )  &
               'Binary collision simulation'
      end if

      if( .not.sn%lmd ) then
!
!  Calculation of the scattering quantities
!
         if( sn%lmehl ) then
            write( 21, '(t16,a,i2,a)' )  &
                  'using ', sn%nmehl, '-point Gauss-Mehler quadrature'
            if( sn%etime > sn%ef ) &
                  write( 21, '(t16,a/)' )  'Calculation of time integral: table'
         else
            text = 'Calculation of scattering angles'
            if( any( sn%lquad(:,:) ) )  text = trim(text) // ' (number of abscissae)'
            do ielem1 = 1, nelem
               do ielem2 = 1, nelem
                  if( sn%lquad(ielem2,ielem1) ) then
                     carr(ielem2,ielem1) = capitalize( trim(sn%quad(ielem2,ielem1)) ) //  &
                           '(' // str( sn%nabsc(ielem2,ielem1) ) // ')'
                  else
                     carr(ielem2,ielem1) = 'table'
                  end if
               end do
            end do
            write( 21, '()' )
            call writeTab( carr(:,:), text )
!                ========      i       i

            if( sn%etime > sn%ef ) then
               text = 'Calculation of time integral'
               if( any( sn%quadtime(:,:) /= 'no' ) )  &
                     text = trim(text) // ' (number of abscissae)'
               carr(:,:) = sn%quadtime(:,:)
               do ielem1 = 1, nelem
                  do ielem2 = 1, nelem
                     if( sn%quadtime(ielem2,ielem1) /= 'no' ) then
                        carr(ielem2,ielem1) = capitalize( trim(carr(ielem2,ielem1)) ) //  &
                              '(' // str( sn%nabsctime(ielem2,ielem1) ) // ')'
                     else
                        carr(ielem2,ielem1) = 'table'
                     end if
                  end do
               end do
               call writeTab( carr(:,:), text )
!                   ========      i       i
            end if
         end if

         if( sn%etime > sn%ef ) then
            if( sn%etime < ion%e )  &
                  write( 21, '(t13,a/)' )  &
                        'Time integral is considered only for e<' // str( sn%etime, 'f99.0' )
         else
            write( 21, '(t13,a/)' )  'Time integral is not considered'
         end if

         call writeTab( sn%scrfun(:,:), 'Screening function' )
!             ========       i                    i
         call writeTab( sn%scatfile(:,:), 'Scattering table file', 8 )
!             ========       i                      i              i
         call writeTab( sn%scrcoeffile(:,:), 'Screening function coefficients file', 12 )
!             ========         i                                i                    i
         call writeTab( sn%scrfile(:,:), 'Screening function file', 10 )
!             ========       i                          i                       i
         call writeTab( sn%scrlen(:,:), 'Screening length' )
!             ========       i                   i
!
!  Scattering parameters
!
         if( sn%chkffpcyl ) then
            write( 21, '(t13,a)' )  &
                  'chkffpcyl=t: collision partners may not be in last ffp cylinder'
         else
            write( 21, '(t13,a)' )  &
                  'chkffpcyl=f: collision partners may be in last ffp cylinder'
         end if

         if( sn%lnegx1 ) then
            write( 21, '(t13,a)' )  &
                  'lnegx1=t: negative (or zero) time integral shift is enforced'
         else
            write( 21, '(t13,a)' )  &
                  'lnegx1=f: negative (or zero) time integral shift is not enforced'
         end if

         if( sn%lposffp ) then
            write( 21, '(t13,a)' )  &
                  'lposffp=t: positive (or zero) ffps are enforced'
         else
            write( 21, '(t13,a)' )  &
                  'lposffp=f: positive (or zero) ffps are not enforced'
         end if

         if( sn%lapsisef ) then
            write( 21, '(t13,a)' )  &
                  'lapsisef=t: Apsis energy is checked against ef for cutoff'
         else
            write( 21, '(t13,a)' )  &
                  'lapsisef=f: Energy after collision is checked against ef for cutoff'
         end if

         if( sn%linel ) then
            write( 21, '(t13,a)' )  &
                  'linel=t: electronic stopping considered in collision kinetics'
         else
            write( 21, '(t13,a)' )  &
                  'linel=f: electronic stopping neglected for collision kinetics'
         end if
!
!  Free flight path parameters
!
         write( 21, '(/t13,a)' )  &
               'Free flight path parameters in amorphous regions: '

         write( 21, '(t13,a)' )  &
               '   (pmax)max=' // str( sn%pmaxmax, 'f99.3' ) //  &
               ', (pmax)min=' // str( sn%pmaxmin, 'f99.3' )
         if( defined( sn%pmaxpsi ) ) then
            write( 21, '(t13,a)', advance='no' )  &
                  '   (pmax)psi=' // str( sn%pmaxpsi, 'f99.3' )
         else
            write( 21, '(t13,a)', advance='no' )  &
                  '   (emin)psi=' // str( nint(sn%eminpsi) )
         end if
         if( defined( sn%eminpsi2 ) )  &
            write( 21, '(a)', advance='no' )  &
                  ', (emin)psi2=' // str( nint(sn%eminpsi2) )
         write( 21, '()' )

         write( 21, '(t13,a)', advance='no' )  &
               '   denffp=' // str( sn%denffp, 'f99.1' ) //  &
               ', psiffp=' // str( sn%psiffp, 'f99.2' )
         if( defined( sn%psiffp2 ) )  &
            write( 21, '(a)', advance='no' )  &
                  ', psiffp2=' // str( sn%psiffp2, 'f99.4' )
         write( 21, '()' )

         if( sn%lstffp ) then
            write( 21, '(t13,a)' )  &
                  '   The free flight paths are chosen statistically'
         else
            write( 21, '(t13,a)' )  &
                  '   The free flight paths are chosen deterministically'
         end if
         if( sn%ffpmax == BIG ) then
            write( 21, '(t13,a)' )  &
                  '   ffpmax is not set'
         else
            write( 21, '(t13,a)' )  &
                  '   ffpmax=' // str( sn%ffpmax, 'f99.3' )
         endif

         if( sn%lpmaxtab )  &
            write( 21, '(t13,a)' )  &
                  'Tables of pmax values are written to file PMAX'

      end if
!
!  Cutoff parameters
!
      write( 21, '(/t13,a)', advance='no' )  &
            'Cutoff: energy=' // str( sn%ef, 'f99.2' )
      if( sn%lpsif ) then
         write( 21, '(a,f6.2)' )  &
               ', angle=' // str( asin( sqrt(sn%dirl2f) )*180/PI, 'f99.1' )
      else
         write( 21, '()' )
      end if
      if( sn%reflectmax < huge( sn%reflectmax ) )  &
         write( 21, '(t13,a)' )  &
               'Maximum number of surface reflections: ' // str( sn%reflectmax )

      close( 21 )

   end subroutine writeSn
!
!=======================================================================
!
end module Param_Sn_mod

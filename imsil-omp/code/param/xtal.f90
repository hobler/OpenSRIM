module Param_Xtal_mod
!
!  Holds data on crystal structures
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: REAL8, INT4,  &
                              ZERO, ONE, TWO, THREE, EIGHT, HALF, QUART, PI
   use String_mod,      only: str
   use Undefine_mod,    only: undefine, defined
   use Open_mod,        only: openFile

   use Param_Setup_mod, only: usecell
   use Param_Ion_mod,   only: ion
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
   public :: readXtal               ! read &crystal record from input file
!            ========
   public :: initXtal               ! initialize the xtal structure
!            ========
   public :: writeXtal              ! write data derived from &crystal
!            =========              !    record to the output file
!
!  Public Data:
!  ===========
!
!  Note: some variables to be eliminated
!
!  Note: There seems to be a compiler error that does not allow to allocate components of a
!        structure outside the subroutine where it has been used for the first time.
!        Therefore name had to be defined outside uc as uc_name.
!
   type :: UnitCell_type                     ! Unit cell data type
!          =============
      character(80) :: name = ' '            ! name of crystallographic system
                                             !    blank means default diamond lattice
      integer(INT4) :: nsite                 ! number of lattice sites in the unit cell
      integer(INT4), allocatable :: ielemr(:)! (nsite) chemical element index (per region)
                                             !    of atom at lattice site. E.g. SiC:
                                             !    ielemr=1 <=> Si, ielemr=2 <=> C
      real(REAL8) :: alat(3)                 ! lattice constants [A] (will be scaled to give
                                             !    right density)
      real(REAL8), allocatable :: pos(:,:)   ! (3,nsite) position of lattice site in units 
                                             !    of the unit cell size (range of pos: [0,1[)
   end type UnitCell_type

   type :: Xtal_type                         ! Crystal data type
!          =========
      logical :: lvib2                       ! 2-d lattice vibration flag

      integer(INT4) :: nrad2                 ! square of the search radius [xscale]
      integer(INT4) :: nloldm                ! maximum number of old lattice sites

      real(REAL8) :: cpsolm                  ! cos(psoldm)
      real(REAL8) :: ffoldm                  ! maximum ffp for old lattice sites [xscale]
      real(REAL8) :: ffpmax                  ! maximum free flight path [xscale]
      real(REAL8) :: pmaxc                   ! maximum impact parameter [xscale]
      real(REAL8) :: rrad                    ! sqrt(nrad2+1) [xscale]
      real(REAL8) :: rrad2                   ! rrad**2

      character(80) :: uc_name               ! temporary, see note above
      type(UnitCell_type) :: uc              ! definition of the unit cell

   end type Xtal_type

   type(Xtal_type), public :: xtal
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine readXtal( )
!             ========
!
!     Read the &crystal record from the input file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      logical :: lvib2
      character(80) :: name
      integer(INT4) :: nloldm, nrad2
      real(REAL8) :: ffpmax, pmax, ffoldm, psoldm

      namelist /crystal/ lvib2,  &
                         name,  &
                         nloldm, nrad2,  &
                         ffoldm, ffpmax, pmax, psoldm

      character(256) :: msg
      integer(INT4) :: ios
!
!  Start of Executable Code:
!  ------------------------
!
!  Set default values of structure variables
!
      xtal%lvib2 = .false.

      call undefine( xtal%uc_name )

      call undefine( xtal%nloldm )
      call undefine( xtal%nrad2 )

      call undefine( xtal%cpsolm )
      call undefine( xtal%ffoldm )
      call undefine( xtal%ffpmax )
      call undefine( xtal%pmaxc )
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
         lvib2 = xtal%lvib2

         call undefine( name )

         call undefine( nloldm )
         call undefine( nrad2 ) 

         call undefine( ffoldm )
         call undefine( ffpmax )
         call undefine( pmax )
         call undefine( psoldm )
!
!  Read parameters from the crystal namlist
!
         read( 11, nml=crystal, iostat=ios, iomsg=msg )
         if( ios < 0 )  exit
         if( ios > 0 ) then
            print '(/a)',  trim( msg )
            goto 900
         end if
!
!  Copy parameters to the xtal structure
!
         xtal%lvib2 = lvib2

         if( defined( name ) )  xtal%uc_name = name

         if( defined( nloldm ) )  xtal%nloldm = nloldm
         if( defined( nrad2 ) )  xtal%nrad2 = nrad2

         if( defined( psoldm ) )  xtal%cpsolm = cos( psoldm*PI/180 )
         if( defined( ffoldm ) )  xtal%ffoldm = ffoldm
         if( defined( ffpmax ) )  xtal%ffpmax = ffpmax
         if( defined( pmax ) )  xtal%pmaxc = pmax
!
!  Check that cannot be done in initXtal
!
         if( defined( psoldm ) ) then
            if( psoldm < ZERO )  goto 962
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
900   error stop '&crystal: error during reading namelist'

962   error stop '&crystal: psoldm>0 required'

   end subroutine readXtal
!
!=======================================================================
!
   subroutine initXtal( )
!             ========
!
!     Initialize the xtal structure
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: isite
      real(REAL8) :: fac(3)
!
!  Start of Executable Code:
!  ------------------------
!
!  Set default values
!
      if( .not.defined( xtal%uc_name ) )  xtal%uc_name = ' '

      if( .not.defined( xtal%nloldm ) )  xtal%nloldm = 6

      if( .not.defined( xtal%cpsolm ) )  xtal%cpsolm = cos(PI/6) ! psoldm = 30
      if( .not.defined( xtal%ffoldm ) )  xtal%ffoldm = 0.50

      if( xtal%uc_name == ' ' ) then            ! default diamond lattice

         if( .not.defined( xtal%ffpmax ) )  xtal%ffpmax = ONE

         if( .not.defined( xtal%nrad2 ) ) then
            if( .not.defined( xtal%pmaxc ) ) then
               xtal%pmaxc = 1.72
               if( ion%e < 10000 ) xtal%pmaxc = xtal%pmaxc / sqrt( 0.3 + 0.7*ion%e/10000 )
            end if
            xtal%nrad2 = xtal%pmaxc**2 + xtal%ffpmax**2
         end if

         if( .not.defined( xtal%pmaxc ) )  &
               xtal%pmaxc = sqrt( xtal%nrad2 + 1 - xtal%ffpmax**2 )

      else                                      ! arbitrary orthogonal lattice

         if( .not.defined( xtal%ffpmax ) )  xtal%ffpmax = QUART

         if( defined( xtal%nrad2 ) )  goto 970

         if( .not.defined( xtal%pmaxc ) ) then
            if( xtal%uc_name == 'diamond' ) then
               xtal%pmaxc=0.43
            else
               goto 971
            end if
         end if
      end if
!
!  Check validity of the parameters
!
      if( xtal%ffpmax < ZERO )  goto 951
      if( xtal%uc_name == ' ' .and. xtal%nrad2 <= 0 )  goto 953
      if( xtal%pmaxc <= ZERO )  goto 954
      if( xtal%nloldm < 0 )  goto 960
      if( xtal%ffoldm < ZERO )  goto 961
!
!  Calculate derived quantities
!
      xtal%uc%name = xtal%uc_name
      if( xtal%uc%name == ' ' ) then            ! default diamond lattice
         xtal%rrad2 = xtal%nrad2 + ONE
      else
         xtal%rrad2 = xtal%pmaxc**2 + xtal%ffpmax**2
      end if
      xtal%rrad = sqrt(xtal%rrad2)

      select case( xtal%uc%name )

         case( ' ' )

            continue                            ! default diamond lattice: nothing to do

         case( 'sc' )

            xtal%uc%alat( 1:3 ) = ONE           ! these parameters are for a simple
            xtal%uc%nsite = 1                   !    cubic lattice
            allocate( xtal%uc%pos( 3, xtal%uc%nsite ) )
            xtal%uc%pos( 1:3, 1 ) = [ ZERO, ZERO, ZERO ]
            allocate( xtal%uc%ielemr( xtal%uc%nsite ) )
            xtal%uc%ielemr( 1:1 ) = 1

         case( 'bcc' )

            xtal%uc%alat( 1:3 ) = ONE           ! these parameters are for a body-centered
            xtal%uc%nsite = 2                   !     cubic lattice
            allocate( xtal%uc%pos( 3, xtal%uc%nsite ) )
            xtal%uc%pos( 1:3, 1 ) = [ ZERO, ZERO, ZERO ]
            xtal%uc%pos( 1:3, 2 ) = [ HALF, HALF, HALF ]
            allocate( xtal%uc%ielemr( xtal%uc%nsite ) )
            xtal%uc%ielemr( 1:2 ) = 1

         case( 'fcc' )

            xtal%uc%alat( 1:3 ) = ONE           ! these parameters are for a face-centered
            xtal%uc%nsite = 4                   !    cubic lattice
            allocate( xtal%uc%pos( 3, xtal%uc%nsite ) )
            xtal%uc%pos( 1:3, 1 ) = [ ZERO, ZERO, ZERO ]
            xtal%uc%pos( 1:3, 2 ) = [ HALF, HALF, ZERO ]
            xtal%uc%pos( 1:3, 3 ) = [ HALF, ZERO, HALF ]
            xtal%uc%pos( 1:3, 4 ) = [ ZERO, HALF, HALF ]
            allocate( xtal%uc%ielemr( xtal%uc%nsite ) )
            xtal%uc%ielemr( 1:4 ) = 1

         case( 'diamond' )

            xtal%uc%alat( 1:3 ) = ONE           ! these parameters are for a diamond lattice
            xtal%uc%nsite = 8                   !    defined via the general cubic lattice
            allocate( xtal%uc%pos( 3, xtal%uc%nsite ) )
            xtal%uc%pos( 1:3, 1 ) = [ ZERO, ZERO, ZERO ]
            xtal%uc%pos( 1:3, 2 ) = [ HALF, HALF, ZERO ]
            xtal%uc%pos( 1:3, 3 ) = [ HALF, ZERO, HALF ]
            xtal%uc%pos( 1:3, 4 ) = [ ZERO, HALF, HALF ]
            xtal%uc%pos( 1:3, 5:8 ) = xtal%uc%pos( 1:3, 1:4 ) + QUART
            allocate( xtal%uc%ielemr( xtal%uc%nsite ) )
            xtal%uc%ielemr( 1:8 ) = 1

         case( 'diamond2' )

            xtal%uc%alat( 1:3 ) = [ ONE, ONE, TWO ]  ! like 'diamond, but with a twice as large
            xtal%uc%nsite = 16                       !    unit cell (meant for testing)
            allocate( xtal%uc%pos( 3, xtal%uc%nsite ) )
            xtal%uc%pos( 1:3, 1 ) = [ ZERO, ZERO, ZERO ]
            xtal%uc%pos( 1:3, 2 ) = [ HALF, HALF, ZERO ]
            xtal%uc%pos( 1:3, 3 ) = [ HALF, ZERO, HALF ]
            xtal%uc%pos( 1:3, 4 ) = [ ZERO, HALF, HALF ]
            xtal%uc%pos( 1:3, 5:8 ) = xtal%uc%pos( 1:3, 1:4 ) + QUART
            xtal%uc%pos( 1:3, 9:16 ) = xtal%uc%pos( 1:3, 1:8 )
            xtal%uc%pos( 3, 9:16 ) = xtal%uc%pos( 3, 9:16 ) + ONE
            xtal%uc%pos( 3, 1:16 ) = xtal%uc%pos( 3, 1:16 ) / TWO
            allocate( xtal%uc%ielemr( xtal%uc%nsite ) )
            xtal%uc%ielemr( 1: 16 ) = 1

         case( 'zincblende', '3C', '3c' )

            xtal%uc%alat( 1:3 ) = ONE
            xtal%uc%nsite = 8
            allocate( xtal%uc%pos( 3, xtal%uc%nsite ) )
            xtal%uc%pos( 1:3, 1 ) = [ ZERO, ZERO, ZERO ]
            xtal%uc%pos( 1:3, 2 ) = [ HALF, HALF, ZERO ]
            xtal%uc%pos( 1:3, 3 ) = [ HALF, ZERO, HALF ]
            xtal%uc%pos( 1:3, 4 ) = [ ZERO, HALF, HALF ]
            xtal%uc%pos( 1:3, 5:8 ) = xtal%uc%pos( 1:3, 1:4 ) + QUART
            allocate( xtal%uc%ielemr( xtal%uc%nsite ) )
            xtal%uc%ielemr( 1:4 ) = 1
            xtal%uc%ielemr( 5:8 ) = 2

         case( 'wurtzite', 'wurzite', '2H', '2h' )

            xtal%uc%alat( 1:3 ) = [ ONE, sqrt(THREE), sqrt(EIGHT/THREE) ]
            xtal%uc%nsite = 8
            allocate( xtal%uc%pos( 3, xtal%uc%nsite ) )
            fac( 1:3 ) = [ ONE/2, ONE/6, ONE/8 ]
            xtal%uc%pos( 1:3, 1 ) = [ 0, 5, 0 ] * fac
            xtal%uc%pos( 1:3, 2 ) = [ 1, 2, 0 ] * fac
            xtal%uc%pos( 1:3, 3 ) = [ 0, 1, 1 ] * fac
            xtal%uc%pos( 1:3, 4 ) = [ 1, 4, 1 ] * fac
            xtal%uc%pos( 1:3, 5 ) = [ 0, 1, 4 ] * fac
            xtal%uc%pos( 1:3, 6 ) = [ 1, 4, 4 ] * fac
            xtal%uc%pos( 1:3, 7 ) = [ 0, 5, 5 ] * fac
            xtal%uc%pos( 1:3, 8 ) = [ 1, 2, 5 ] * fac
            allocate( xtal%uc%ielemr( xtal%uc%nsite ) )
            xtal%uc%ielemr( 1:2 ) = 1
            xtal%uc%ielemr( 3:4 ) = 2
            xtal%uc%ielemr( 5:6 ) = 1
            xtal%uc%ielemr( 7:8 ) = 2

         case( '4H', '4h' )

            xtal%uc%alat( 1:3 ) = [ ONE, sqrt(THREE), sqrt(32/THREE) ]
            xtal%uc%nsite = 16
            allocate( xtal%uc%pos( 3, xtal%uc%nsite ) )
            fac( 1:3 ) = [ ONE/2, ONE/6, ONE/16 ]
            xtal%uc%pos( 1:3, 1 ) = [ 0, 0, 0 ] * fac     ! A
            xtal%uc%pos( 1:3, 2 ) = [ 1, 3, 0 ] * fac     ! A
            xtal%uc%pos( 1:3, 3 ) = [ 0, 0, 3 ] * fac     ! A
            xtal%uc%pos( 1:3, 4 ) = [ 1, 3, 3 ] * fac     ! A
            xtal%uc%pos( 1:3, 5 ) = [ 1, 1, 4 ] * fac     ! B
            xtal%uc%pos( 1:3, 6 ) = [ 0, 4, 4 ] * fac     ! B
            xtal%uc%pos( 1:3, 7 ) = [ 1, 1, 7 ] * fac     ! B
            xtal%uc%pos( 1:3, 8 ) = [ 0, 4, 7 ] * fac     ! B
            xtal%uc%pos( 1:3, 9 ) = [ 0, 0, 8 ] * fac     ! A
            xtal%uc%pos( 1:3, 10 ) = [ 1, 3, 8 ] * fac    ! A
            xtal%uc%pos( 1:3, 11 ) = [ 0, 0, 11 ] * fac   ! A
            xtal%uc%pos( 1:3, 12 ) = [ 1, 3, 11 ] * fac   ! A
            xtal%uc%pos( 1:3, 13 ) = [ 0, 2, 12 ] * fac   ! C
            xtal%uc%pos( 1:3, 14 ) = [ 1, 5, 12 ] * fac   ! C
            xtal%uc%pos( 1:3, 15 ) = [ 0, 2, 15 ] * fac   ! C
            xtal%uc%pos( 1:3, 16 ) = [ 1, 5, 15 ] * fac   ! C
!            xtal%uc%pos( 1:3, 1 ) = [ 0, 5, 0 ] * fac
!            xtal%uc%pos( 1:3, 2 ) = [ 1, 2, 0 ] * fac
!            xtal%uc%pos( 1:3, 3 ) = [ 0, 3, 1 ] * fac
!            xtal%uc%pos( 1:3, 4 ) = [ 1, 0, 1 ] * fac
!            xtal%uc%pos( 1:3, 5 ) = [ 0, 3, 4 ] * fac
!            xtal%uc%pos( 1:3, 6 ) = [ 1, 0, 4 ] * fac
!            xtal%uc%pos( 1:3, 7 ) = [ 0, 1, 5 ] * fac
!            xtal%uc%pos( 1:3, 8 ) = [ 1, 4, 5 ] * fac
!            xtal%uc%pos( 1:3, 9 ) = [ 0, 1, 8 ] * fac
!            xtal%uc%pos( 1:3, 10 ) = [ 1, 4, 8 ] * fac
!            xtal%uc%pos( 1:3, 11 ) = [ 0, 0, 9 ] * fac
!            xtal%uc%pos( 1:3, 12 ) = [ 1, 3, 9 ] * fac
!            xtal%uc%pos( 1:3, 13 ) = [ 0, 0, 12 ] * fac
!            xtal%uc%pos( 1:3, 14 ) = [ 1, 3, 12 ] * fac
!            xtal%uc%pos( 1:3, 15 ) = [ 0, 5, 13 ] * fac
!            xtal%uc%pos( 1:3, 16 ) = [ 1, 2, 13 ] * fac
            allocate( xtal%uc%ielemr( xtal%uc%nsite ) )
            xtal%uc%ielemr( 1:2 ) = 1
            xtal%uc%ielemr( 3:4 ) = 2
            xtal%uc%ielemr( 5:6 ) = 1
            xtal%uc%ielemr( 7:8 ) = 2
            xtal%uc%ielemr( 9:10 ) = 1
            xtal%uc%ielemr( 11:12 ) = 2
            xtal%uc%ielemr( 13:14 ) = 1
            xtal%uc%ielemr( 15:16 ) = 2

         case( '6H', '6h' )

            xtal%uc%alat( 1:3 ) = [ ONE, sqrt(THREE), sqrt(24d0) ]
            xtal%uc%nsite = 24
            allocate( xtal%uc%pos( 3, xtal%uc%nsite ) )
            fac( 1:3 ) = [ ONE/2, ONE/6, ONE/24 ]
            xtal%uc%pos( 1:3, 1 ) = [ 0, 0, 0 ] * fac     ! A
            xtal%uc%pos( 1:3, 2 ) = [ 1, 3, 0 ] * fac     ! A
            xtal%uc%pos( 1:3, 3 ) = [ 0, 0, 3 ] * fac     ! A
            xtal%uc%pos( 1:3, 4 ) = [ 1, 3, 3 ] * fac     ! A
            xtal%uc%pos( 1:3, 5 ) = [ 0, 4, 4 ] * fac     ! B
            xtal%uc%pos( 1:3, 6 ) = [ 1, 1, 4 ] * fac     ! B
            xtal%uc%pos( 1:3, 7 ) = [ 0, 4, 7 ] * fac     ! B
            xtal%uc%pos( 1:3, 8 ) = [ 1, 1, 7 ] * fac     ! B
            xtal%uc%pos( 1:3, 9 ) = [ 0, 2, 8 ] * fac     ! C
            xtal%uc%pos( 1:3, 10 ) = [ 1, 5, 8 ] * fac    ! C
            xtal%uc%pos( 1:3, 11 ) = [ 0, 2, 11 ] * fac   ! C
            xtal%uc%pos( 1:3, 12 ) = [ 1, 5, 11 ] * fac   ! C
            xtal%uc%pos( 1:3, 13 ) = [ 0, 0, 12 ] * fac   ! A
            xtal%uc%pos( 1:3, 14 ) = [ 1, 3, 12 ] * fac   ! A
            xtal%uc%pos( 1:3, 15 ) = [ 0, 0, 15 ] * fac   ! A
            xtal%uc%pos( 1:3, 16 ) = [ 1, 3, 15 ] * fac   ! A
            xtal%uc%pos( 1:3, 17 ) = [ 0, 2, 16 ] * fac   ! C
            xtal%uc%pos( 1:3, 18 ) = [ 1, 5, 16 ] * fac   ! C
            xtal%uc%pos( 1:3, 19 ) = [ 0, 2, 19 ] * fac   ! C
            xtal%uc%pos( 1:3, 20 ) = [ 1, 5, 19 ] * fac   ! C
            xtal%uc%pos( 1:3, 21 ) = [ 0, 4, 20 ] * fac   ! B
            xtal%uc%pos( 1:3, 22 ) = [ 1, 1, 20 ] * fac   ! B
            xtal%uc%pos( 1:3, 23 ) = [ 0, 4, 23 ] * fac   ! B
            xtal%uc%pos( 1:3, 24 ) = [ 1, 1, 23 ] * fac   ! B
            allocate( xtal%uc%ielemr( xtal%uc%nsite ) )
            xtal%uc%ielemr( 1:2 ) = 1
            xtal%uc%ielemr( 3:4 ) = 2
            xtal%uc%ielemr( 5:6 ) = 1
            xtal%uc%ielemr( 7:8 ) = 2
            xtal%uc%ielemr( 9:10 ) = 1
            xtal%uc%ielemr( 11:12 ) = 2
            xtal%uc%ielemr( 13:14 ) = 1
            xtal%uc%ielemr( 15:16 ) = 2
            xtal%uc%ielemr( 17:18 ) = 1
            xtal%uc%ielemr( 19:20 ) = 2
            xtal%uc%ielemr( 21:22 ) = 1
            xtal%uc%ielemr( 23:24 ) = 2

         case default

            call openFile( 129, trim(xtal%uc%name), 'formatted', 'read' )

            read( 129, *, err=990, end=990 )  xtal%uc%alat(1:3)
            read( 129, *, err=990, end=990 )  xtal%uc%nsite
            allocate( xtal%uc%pos( 3, xtal%uc%nsite ) )
            allocate( xtal%uc%ielemr( xtal%uc%nsite ) )
            do isite = 1, xtal%uc%nsite
               read( 129, *, err=990, end=990 )  xtal%uc%pos(1:3,isite), xtal%uc%ielemr(isite)
            end do

      end select ! can we check the correct number and abundance of elements? -> param_mat_mod !

      return

951   error stop '&crystal: ffpmax>0 required'
953   error stop '&crystal: nrad2>0 required'
954   error stop '&crystal: pmax>0 required'

960   error stop '&crystal: nloldm>0 required'
961   error stop '&crystal: ffoldm>0 required'

970   error stop '&crystal: nrad2 and name must not both be specified'
971   error stop '&crystal: no default values of pmax available for lattice other than diamond'

990   error stop '&crystal: file specified by name has wrong format'

   end subroutine initXtal
!
!=======================================================================
!
   subroutine writeXtal( )
!             =========
!
!     Write parameters derived from the &crystal record to the output file.
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
      if( usecell )  return                     ! No crystal defined if geometry and material
!                                               !    information is read from a cell file
!  Open the output file
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i

      if( xtal%uc%name == ' ' ) then
         write( 21, '(t2,a)' )  &
               '&crystal:  default diamond lattice (units: [xscale])'
      else
         write( 21, '(/t2,a)' )  &
               '&crystal:  ' // str( xtal%uc%name ) // ' lattice (units: [xscale])'
      end if

      write( 21, '(t13,a)' )  &
            'pmax=' // str( xtal%pmaxc, 'f99.3' ) //  &
            ', search distance=' // str( xtal%rrad, 'f99.3' ) //  &
            ', ffpmax=' // str( xtal%ffpmax, 'f99.3' )

      write( 21, '(t13,a)' )  &
            'psoldm=' // str( 180./PI*acos( xtal%cpsolm ), 'f99.1' ) //  &
            ', ffoldm=' // str( xtal%ffoldm, 'f99.3' ) //  &
             ', nloldm=' // str( xtal%nloldm )

      write( 21, '()' )

      close( 21 )

   end subroutine writeXtal
!
!=======================================================================
!
end module Param_Xtal_mod

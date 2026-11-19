module Param_mod
!
!  Read input data from input file INP and write information (except simulation results) 
!  to the output file OUT
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: VERSION,  &
                                    INT4, REAL8,  &
                                    ZERO

   use Open_mod,              only: openFile
   use Time_mod,              only: getCPUtime, initCPUtime, getElapsedTime, initElapsedTime

   use Param_Setup_mod,       only: nr, usebk, usegeom,  &
                                    readSetup, initSetup, writeSetup,  &
                                    ldyn, usecell, cellfile, cellrec
   use Param_Geom_mod,        only: geom, readGeom, initGeom, writeGeom
   use Param_Elem_mod,        only: elem, readAtoms, initAtoms, writeAtoms
   use Param_Ion_mod,         only: readIon, initIon, writeIon
   use Param_Xtal_mod,        only: readXtal, initXtal, writeXtal
   use Param_Mat_mod,         only: mat, irxtal, readMat, initMat, writeMat
   use Param_Dam_mod,         only: dam, readDam, initDam, writeDam
   use Param_Sn_mod,          only: readSn, initSn, writeSn
   use Param_Se_mod,          only: readSe, initSe, writeSe
   use Param_Out_mod,         only: out, readOut, initOut, writeOut

   use Target_Cell_mod,       only: readAtomsFromCell, readBndFromCell
! 
!      use Surf_mod
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
   public :: readParams             ! read parameters from the input file
!            ==========
   public :: initParams             ! initialize parameters in Param_*_mod
!            ==========
   public :: writeParams            ! write input prameters and some derived parameters to
!            ===========            !   the output file
!
!  Private Procedures:
!  ==================
!
   private :: writeHeader           ! write header of output file
!             ===========
!
!  Private Data:
!  ============
!
!  None
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine readParams( )
!             ==========
!
!     Read the input file
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
!  Initialize the CPU time
!
      call initCPUtime( )
!          ===========
      call initElapsedTime( )
!          ===============
!
!  Read the records of the input file and initializes the
!  model parameter structures
!
      call readSetup( )
!          =========

      if( .not.(usebk .and. usegeom) )  call readGeom( )
!                                            ========

      call readAtoms( )
!          =========

      call readIon( )
!          =======

      call readXtal( )
!          ========

      call readMat( )
!          =======

      if( usecell )  call readAtomsFromCell( )
!                         =================

      call readDam( )
!          =======

      call readSn( )
!          ======

      call readSe( )
!          ======

      call readOut( )
!          =======

   end subroutine readParams
!
!=======================================================================
!
   subroutine initParams( )
!             ==========
!
!     Initialize parameters in Param_*_mod, using the parameters from the input file
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
!  Initialize input parameters
!
      call initSetup( )
!          =========

      if( .not.(usebk .and. usegeom) )  call initGeom( )
!                                            ========

      call initAtoms( )
!          =========

      if( usecell )  call readBndFromCell( )    ! initialize geom%zlim(:) etc.
!                         ===============

      call initIon( )
!          =======

      call initXtal( )
!          ========

      call initMat( )
!          =======

      call initDam( )
!          =======

      call initSn( )
!          ======

      call initSe( )
!          ======

      call initOut( )
!          =======

   end subroutine initParams
!
!=======================================================================
!
   subroutine writeParams( )
!             ===========
!
!     Report input and derived parameters to the output file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(9) :: ctime
!
!  Start of Executable Code:
!  ------------------------
!
!  Write the header of the output file
!
      call writeHeader( )
!          ===========
!
!  Write the parameters derived from the input data to the output file
!
      call writeSetup( )
!          ==========

      call writeGeom( )
!          =========

      call writeAtoms( )
!          ==========

      call writeIon( )
!          ========

      if( irxtal > 0 )  call writeXtal( )
!                            =========

      call writeMat( )
!          ========

      call writeDam( )
!          ========

      call writeSn( )
!          =======

      call writeSe( )
!          =======

      call writeOut( )
!          ========
!
!  Write CPU and elapsed time for initialization
!
      call openFile( 21, 'OUT', 'formatted', 'append' )
!          ========  i     i         i          i

      call getCPUtime( ctime )
!          ==========    o
      write( 21, '(t2,a,a,a)' ) 'Initialization time: ', ctime, ' (hh:mm:ss)'

      call getElapsedTime( ctime )
!          ==============    o
      write( 21, '(t2,a,a,a)' ) 'Elapsed time       : ', ctime, ' (hh:mm:ss)'

      write( 21, '(a)' ) ' '
      close( 21 )

   end subroutine writeParams
!
!=======================================================================
!
   subroutine writeHeader( )
!             ===========
!  
!     Write header to output file and copy header of input file to output
!     file. The header includes a list of atoms with their properties and 
!     abundance in ions and target materials.
!
!  Arguments:
!
!     None
!
!  Program History:
!
!
!  Namelist:
!
      character(256) :: hostname, workpath, tablepath
      namelist / path / hostname, workpath, tablepath
!
!  Local Variables:
!
      logical expath
      character(3) :: cmonth(12) = [ 'Jan','Feb','Mar','Apr','May','Jun', &
                                      'Jul','Aug','Sep','Oct','Nov','Dec' ]
      character(20) :: cdate             
      character(40) :: hostname32
      integer :: idat(8), istat, nh

      character(80) :: line
!
!  Start of Executable Code:
!  ------------------------
!
!  Open input and output file
!
      call openFile( 11, 'INP', 'formatted', 'read' )
!          ========  i     i         i          i
      call openFile( 21, 'OUT', 'formatted', 'write' )
!          ========  i     i         i          i
!
!  Check if path file exists, and read paths
!
      hostname = ' '
      inquire( file='PATH', exist=expath )
      if( expath ) then
         open( 42, file='PATH', form='formatted', action='read', position='rewind' )
         read( 42, nml=path, iostat=istat )
         close( 42 )
      end if
      nh = len_trim(hostname)
      hostname32 = ' '
      hostname32(max(1,15-nh/2):) = hostname
!
!  Date and time
!
      call date_and_time( values=idat )
!          =============       o
      write( cdate, '(i2,a,a,a,i4,i3,a,i2.2,a,i2.2)' ) &
            idat(3), '-', cmonth(idat(2)), '-', idat(1), idat(5), ':', idat(6), ':', idat(7)
!
!  Write header of output file
!
      write( 21, '(t2,a,t29,a,t60,a20)' ) &
            'IMSIL (Version ' // trim(VERSION) // ')', hostname32, cdate
!
!  Copy the header of the input file to the output file
!
      write( 21, '(t2,a)' ) ' '
      do
         read( 11, '(a)', iostat=istat ) line
         if( istat /= 0  )  exit
         line = adjustl( line )                  ! remove leading blanks
!               =======   i
         if( line(1:1) == '&'  )  exit
         write( 21, '(t2,a)' ) line(1:78)
      end do
!
!  State units of input parameters
!
      write( 21, '(/t2,a)' ) &
            'INPUT PARAMETERS:'

      write( 21, '(t2,a)' ) &
            '(units: [energy]=eV, [angle]=deg, [length]=A, [conc.]=cm-3)'

      write( 21, '()' )

      close( 21 )

   end subroutine writeHeader
!
!=======================================================================
!
end module Param_mod

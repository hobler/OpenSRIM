module Open_mod
!
!  Open files
!
!  Used Modules:
!  ============
!
   use Const_mod, only: INT4, REAL8
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
   public :: openFile               ! opens a file
!            ========
!
contains
!
!=======================================================================
!
   subroutine openFile( nunit, fname, fform, fmode )
!             ========    i      i      i      i
!
!     Opens the file fname with unit number nunit, format fform, and action
!     fmode. According to the value of fmode three cases are distinguished:
!
!     fmode='read':   file must exist and will be positioned at the beginning.
!     fmode='write':  file may or may not exist and will be positioned at the
!                     beginning if it already exists, i.e. it will be overwritten.
!     fmode='append': file may or may not exist and will be positioned at the
!                     end if it already exists, i.e. new data will be appended.
!
!     If the file 'PATH' exists and it contains a namelist record of the form
!
!        &path workpath='...' tablepath='...' /
!
!     "read" files are searched for in the local directory, then in workpath,
!     and finally in tablepath. "write" files are searched for and written to
!     workpath only.
!
!     fform may be 'formatted' or 'unformatted'
!
!  Arguments:
!
      integer(INT4), intent(in) :: nunit        ! file unit number
      character(*), intent(in) :: fname         ! file name
      character(*), intent(in) :: fform         ! file format (formatted/unformatted)
      character(*), intent(in) :: fmode         ! file action mode (read/write/append)
!
!  Namelist:
!
      character(256) :: hostname, workpath, inppath, tablepath
      namelist / path / hostname, workpath, inppath, tablepath
!
!  Local Variables:
!
      logical :: expath, ex
      character(10) :: fact, fstat, fpos
      character(256) :: pfname
      integer(INT4) :: istat
!
!  Start of Executable Code:
!  ------------------------
!
!  Set status, action, and position
!
      select case( fmode )

         case( 'append' )

            fact = 'write'
            fstat = 'unknown'
            fpos = 'append'

         case( 'write' )

            fact = 'write'
            fstat = 'replace'
            fpos = 'rewind'

         case( 'read' )

            fact = 'read'
            fstat = 'old'
            fpos = 'rewind'

         case default

            error stop 'openFile: illegal value of fmode'

      end select
!
!  Check if path file exists, and read paths
!
      workpath = ' '
      inppath = ' '
      tablepath = ' '
      inquire( file='PATH', exist=expath )
      if( expath ) then
         open( 42, file='PATH', form='formatted', action='read', position='rewind' )
         read( 42, nml=path, iostat=istat )
         close( 42 )
      end if
!
!  For a "read" file try workpath first, then inppath, finally tablepath
!
      if( fact=='read' ) then
         pfname = trim( workpath ) // trim( fname )
!         pfname = fname
         inquire( file=pfname, exist=ex )
         if( .not.ex ) then
            if( fname == 'INP' )  error stop 'INP file not found'
            pfname = trim( inppath ) // trim( fname )
            inquire( file=pfname, exist=ex )
            if( .not.ex ) then
               pfname = trim( tablepath ) // trim( fname )
               inquire( file=pfname, exist=ex )
               if( .not.ex ) then
                  !call tracebackqq()
                  print *,  'openFile: file ''' // trim( fname ) // ''' does not exist'
                  error stop
               end if
            end if
         end if
!
!  For a "write" file use workpath
!
      else

         pfname = trim( workpath ) // trim( fname )
         inquire( file=pfname, exist=ex )
         
      end if
!
!  Open the file
!
      open( nunit, file=trim(pfname), form=fform, action=fact, status=fstat, position=fpos,  &
            iostat=istat )

      if( istat /= 0 ) then
         !call tracebackqq()
         print *,  'openFile: error during opening file ''' // trim(pfname) // ''''
         error stop
      end if

   end subroutine openFile
!
!=======================================================================
!
end module Open_mod

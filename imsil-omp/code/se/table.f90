module Se_Table_mod
!
!  Handle electronic stopping power defined in file. The file may be one of the following:
!
!  1) se%semodel='table': a list of (energy,dE/dx) values. Units are ([eV],[eV/A]).
!  2) se%semodel='srim': a SRIM-2013 stopping power file containing (energy,se(1:92)) tuples
!     where se(1:92) denotes the stopping power of the elemental targets with Z2=1:92 for a
!     Z1 specific to the file. Units are ([eV],[eV A^2]).
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8, ZERO 
   use Open_mod,        only: openFile

   use Param_Setup_mod, only: nelem
   use Param_Elem_mod,  only: elem
   use Param_Se_mod,    only: se
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
   public :: initSeTable         ! read se from se%sefile
!            ===========
   public :: evalSeTable         ! evaluate dE/dx for given energy
!            ===========
!
!  Private Data:
!  ============
!
   logical :: init = .false.
   integer(INT4), allocatable :: nse(:,:)       ! (NELEM,NELEM)
                                                ! number of valid e/se pairs
   real(REAL8), allocatable :: etab(:,:,:)      ! (NELEM,NELEM,NPOINTS)
   real(REAL8), allocatable :: setab(:,:,:)     ! (NELEM,NELEM,NPOINTS)
!
!=============================================================================================
!
contains
!
!=============================================================================================
!
   subroutine initSeTable( )
!             ===========
!
!     Read the electronic stopping defined in the file
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(80) :: text
      integer(INT4) :: ios, ise, icom, ielem1, ielem2, i
      integer(INT4), allocatable :: ncom(:,:)
      real(REAL8) :: dummy
!
!  Start of Executable Code:
!  ------------------------
!
      init = .true.
!
!  Allocate memory
!
      allocate( nse( 1:nelem, 1:nelem) )
      allocate( ncom( 1:nelem, 1:nelem) )
      nse(:,:) = 0
      ncom(:,:) = 0

      do ielem2 = 1, nelem
         do ielem1 = 1, nelem
            if( se%semodel(ielem2,ielem1) /= 'table' .and.  &
                se%semodel(ielem2,ielem1) /= 'srim' ) cycle

            call openFile( 22, se%sefile(ielem2,ielem1), 'formatted', 'read' )
!                ========  i            i                     i         i

            ise = 0
            icom = 0
!
! Determine number of e/se points
!
            do
               read( 22, *, iostat = ios ) text
               if( ios < 0  )  exit
               if( ios > 0 )  goto 901
               if( text(1:1) == '*' .or. text(1:1) == '#' ) then
                  if( ise > 0 )  goto 905
                  icom = icom + 1
               else
                  ise = ise + 1
               end if
            end do

            if( ise == 0 )  goto 902
            nse(ielem2,ielem1) = ise
            ncom(ielem2,ielem1) = icom

            close( 22 )

         end do
      end do

      allocate( etab( 1:nelem, 1:nelem, 0:maxval(nse(:,:)) ) )
      allocate( setab( 1:nelem, 1:nelem, 0:maxval(nse(:,:)) ) )
      etab(:,:,:) = ZERO
      setab(:,:,:) = ZERO

      do ielem2 = 1, nelem
         do ielem1 = 1, nelem
            if( se%semodel(ielem2,ielem1) /= 'table' .and.  &
                se%semodel(ielem2,ielem1) /= 'srim' ) cycle

            call openFile( 22, se%sefile(ielem2,ielem1), 'formatted', 'read' )
!                ========  i            i                     i         i
!
!  Zero point. Necessary for interpolation
!
            etab(ielem2,ielem1,0) = ZERO
            setab(ielem2,ielem1,0) = ZERO
!
!  Read the comment lines
!
            do icom = 1, ncom(ielem2,ielem1)
               read( 22, *, iostat = ios ) text
               if( ios < 0 )  exit
               if( ios > 0 )  goto 901
               if( text(1:1) /= '*' .and. text(1:1) /= '#' )  goto 904
            end do
!
!  Read the valid data
!
            do ise = 1, nse(ielem2,ielem1)
               if( se%semodel(ielem2,ielem1) == 'table' ) then
                  read( 22, *, iostat=ios )  &
                        etab(ielem2,ielem1,ise), setab(ielem2,ielem1,ise)
               else
                  read( 22, *, iostat=ios )  &
                        etab(ielem2,ielem1,ise),  &
                        (dummy, i=1,elem(ielem2)%z-1),  &
                        setab(ielem2,ielem1,ise),  &
                        (dummy, i=elem(ielem2)%z+1,92)
               end if
               if( ios < 0 )  exit
               if( ios > 0 )  goto 901
            end do
            close( 22 )

         end do
      end do
!
!  Error messages
!
      return

901   error stop 'initSeTable: error during reading sefile'
902   error stop 'initSeTable: no valid data points in sefile'
904   error stop 'initSeTable: error reading comment'
905   error stop 'initSeTable: comments in sefile can only precede data'

   end subroutine initSeTable
!
!===========================================================================================
!        
   subroutine evalSeTable( e, ielem1, ielem2, se_val )
!             ===========  i      i      i      o
!
!     Interpolate the electronic stopping from the data read from the file
!
!  Arguments:
!
      real(REAL8), intent(in) :: e              ! ion energy [eV]
      integer(INT4), intent(in) :: ielem1       ! chemical element index of the ion
      integer(INT4), intent(in) :: ielem2       ! index of the atom
      real(REAL8), intent(out) :: se_val        ! stopping power [eV A^2]
!
!  Local Variables:
!
      integer(INT4) :: im, imin, imax, ie
      real(REAL8) :: a, b
!
!  Start of Executable Code:
!  ------------------------
!
      if( .not.init )  goto 901
      if( se%semodel(ielem2,ielem1) /= 'table' .and.  &
          se%semodel(ielem2,ielem1) /= 'srim' )  goto 902
!
!  Set se to zero if energy is beyond maximum tabulated value
!
      if( e > etab(ielem2,ielem1,nse(ielem2,ielem1)) ) then
         se_val = ZERO
         return
      end if
!
!  Find the energy box by binary search
!
      imin = 0
      imax = nse( ielem2, ielem1 )

111   im = ( imin + imax ) / 2
      if( e > etab( ielem2, ielem1, im ) ) then
         imin = im
      else
         imax = im
      end if
      if( imax - imin > 1 )  goto 111
      ie = imin
!
!  Interpolate se which corresponds to eProj
!
      a = ( setab(ielem2,ielem1,ie) - setab(ielem2,ielem1,ie+1) ) /  &
           ( etab(ielem2,ielem1,ie) - etab(ielem2,ielem1,ie+1) )
      b = ( etab(ielem2,ielem1,ie) * setab(ielem2,ielem1,ie+1) -  &
             setab(ielem2,ielem1,ie) *  &
           etab(ielem2,ielem1,ie+1) ) / ( etab( ielem2, ielem1, ie ) -  &
             etab(ielem2,ielem1,ie+1) )

      se_val = a * e + b  ! e = a * z + b

      if( se%semodel(ielem2,ielem1) == 'table' )  se_val = se_val / elem(ielem2)%dens

      return

901   stop 'evalSeTable: not initialized'
902   stop 'evalSeTable: model 7 or 8 not defined for this (atom2,atom1)'

   end subroutine evalSeTable
!
!============================================================================================
!
end module Se_Table_mod

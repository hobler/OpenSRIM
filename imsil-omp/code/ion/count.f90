module Ion_Count_mod
!
!  Contains the ion and projectile counter
!
!  Terms: Ion =        Atom/molecule/recoil to be implanted or started, 
!                      before interaction with the target
!
!         Projectile = Particle moving in the target 
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8
   use Open_mod,           only: openFile
   use Math_Intmath_mod,   only: index1to2, index2to1

   use Param_Setup_mod,    only: nelemIon, nion, nion1, nionchunk, nthread
   use Param_Ion_mod,      only: ion
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
   public :: initCounter            ! initialize the projectiles counters
!            ===========
   public :: advanceCounter         ! advance the ion/projectile counter
!            ==============         !    and return some info
   public :: setChunkCounter        ! set the ion chunk counter
!            ===============
   public :: lastProjectile         ! determine if this is the last projectile,
!            ==============         !    optionally out of a group
   public :: getIonIndex            ! get the index of the last ion of the most recent chunk
!            ===========
!
!  Private data:
!  ============
!
   integer(INT4) :: ichunk             ! ion chunk counter
   integer(INT4) :: iionchunk          ! ion counter within current chunk
   integer(INT4) :: ielem = 1          ! ion atom element index (e.g.: index of B, F in
                                       !    molecule)
   integer(INT4) :: iatomi = 0         ! ion atom index corresponding to atoms of the
                                       !    same element in the ion
                                       !    (e.g. ielem=1 (B): iatomi=1..1
                                       !          ielem=2 (F): iatomi=1..2)
contains
!
!==============================================================================
!
   subroutine initCounter( )
!             ===========
!
!     Initialize the projectiles counters to one before the first projectile
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
      if( ion%usetra ) then
         call openFile( 123, trim(ion%trafile), 'formatted', 'read' )
!             ========   i          i                i         i
      end if

      ielem = nelemIon
      iatomi = ion%natom(ielem)

   end subroutine initCounter
!
!=======================================================================
!
   subroutine advanceCounter( ielem1, newmolecule, done )
!             ==============    o         o         o
!
!     Advance the counter and return atom index and flag indicating whether
!     this is a new molecule
!
!  Arguments:
!
      integer(INT4), intent(out) :: ielem1      ! atom index of projectile (undefined if done=t)
      logical, intent(out) :: newmolecule       ! flag indicating whether projectile is from a
                                                !    new molecule (undefined if done=t)
      logical, intent(out) :: done              ! flag indicating no more ions to simulate
                                                !    in current chunk
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      done = .false.

      if( ion%usetra ) then                     ! no molecules can be identified on a tra file

         iionchunk = iionchunk + 1
         ielem1 = 1
         newmolecule = .true.

      else

         if( iatomi == ion%natom(ielem) ) then
            if( ielem == nelemIon ) then
               done = ( iionchunk == nionchunk/nthread )
               if( done )  return
               iionchunk = iionchunk + 1
               ielem = 1
               iatomi = 1
               newmolecule = .true.
            else
               ielem = ielem + 1
               iatomi = 1
               newmolecule = .false.
            end if
         else
            iatomi = iatomi + 1
            newmolecule = .false.
         end if

         ielem1 = ielem
! debug
!      print *, 'ichunk, iionchunk, ielem, iatomi=', ichunk, iionchunk, ielem, iatomi
! end-debug

      end if

   end subroutine advanceCounter
!
!=======================================================================
!
   subroutine setChunkCounter( ichunk_, iionchunk_ )
!             ===============     i         i
!
!     Set the chunk counter
!
!  Arguments:
!
      integer(INT4), intent(in) :: ichunk_      ! chunk counter
      integer(INT4), intent(in) :: iionchunk_   ! first ion index within chunk minus 1
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      ichunk = ichunk_
      iionchunk = iionchunk_

   end subroutine setChunkCounter
!
!=======================================================================
!
   function lastProjectile( ngroup_ )  result( last )
!           ==============     i                o
!
!     Determine if the last ion of the most recent chunk was the last one out of a group of
!     ngroup projectile molecules AND the current ion is the last atom of an ion molecule
!
!     Result is false if ngroup<=0
!
!     If usetra=t and ngroup_ is not given, last=.true. for the last projectile on the TRA file
!     (unit 123)
!
!  Arguments:
!
      integer(INT4), intent(in), optional :: ngroup_  ! size of group
      logical :: last                                 ! flag indicating projectile is last in
                                                      !    group
!
!  Local Variables:
!
      integer(INT4) :: iion, ngroup, i, n, iflag, idummy
      real(REAL8) :: rdummy
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( ngroup_ ) ) then
         ngroup = ngroup_
      else
         ngroup = nion
      end if

      if( ngroup <= 0 ) then
         last = .false.
         return
      end if

      iion = ichunk * nionchunk

      if( ion%usetra ) then

         if( present( ngroup_ ) ) then
            last = ( mod( iion, ngroup ) == 0 )
         else                                   ! last ion if no further lines with iflag/=1,2
            last = .true.
            n = 0
            do
               n = n + 1
               read( 123, *, end=100 )  (rdummy, i=1,7), idummy, idummy, iflag
               if( iflag /= 1 .and. iflag /= 2 )  exit
            end do
            last = .false.
100         do i = 1, n                         ! restore original file position
               backspace 123
            enddo
         end if

      else

         last = ( mod( iion, ngroup ) == 0  &
                  .and. ielem == nelemIon .and. iatomi == ion%natom(ielem) )
         if( ngroup /= 1 )  &
               print *, 'lastProjectile: ngroup=', ngroup, ', ichunk=', ichunk,  &
                        ', iion=', iion, ', last=', last

      end if

   end function lastProjectile
!
!=======================================================================
!
   subroutine getIonIndex( iion )
!             ===========   o
!
!     Get the index of the ion. In multiprocessing (nthread>1), it is assumed that all
!     threads have processed the same number of ions, and the current thread processes
!     the last one of the ions currently running in parallel. Note that this is not
!     strictly the case.
!
!  Arguments:
!
      integer(INT4), intent(out) :: iion                 ! index of ion
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      call index2to1( iionchunk, ichunk, nionchunk/nthread, iion )
!          =========      i         i           i            o
      iion = iion * nthread

   end subroutine getIonIndex
!
!=======================================================================
!
end module Ion_Count_mod

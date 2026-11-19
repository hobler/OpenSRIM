module Ion_Count_mod_old
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
   use iso_fortran_env,    only: lock_type

   use Const_mod,          only: INT4, REAL8
   use Open_mod,           only: openFile

   use Param_Setup_mod,    only: nelemIon, nion, nion1
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
   public :: advanceCounter_single  ! single-image version of advanceCounter
!            =====================
   public :: advanceCounter         ! advances the ion/projectile counter
!            ==============         !    and returns some info
   public :: lastProjectile         ! determine if this is the last projectile,
!            ==============         !    optionally out of a group
   public :: getIonIndex            ! get the index of the most recent ion
!            ===========
!
!  Private data:
!  ============
!
   type(lock_type) :: iion_lock[*]     ! lock for iion_global

   integer(INT4) :: iion_global[*] = 0 ! global ion counter (e.g.: of BF2 ions)
   ! The following variables are per image:
   integer(INT4) :: iion               ! ion counter
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

      iion_global = nion1

      iion = nion1
      ielem = nelemIon
      iatomi = ion%natom(ielem)

   end subroutine initCounter
!
!=======================================================================
!
   subroutine advanceCounter_single( ielem1, newmolecule, done )
!             =====================    o         o         o
!
!     Advance the counter and return atom index and flag indicating whether
!     this is a new molecule
!
!     This is the version for single-image execution. It may enhance efficiency for
!     single image execution significantly, in particular if cascade simulation is
!     not done or cheap.
!
!  Arguments:
!
      integer(INT4), intent(out) :: ielem1      ! atom index of projectile
      logical, intent(out) :: newmolecule       ! flag indicating whether projectile is from a
                                                !    new molecule
      logical, intent(out) :: done              ! flag indicating no more ions to simulate
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      done = .false.

      if( ion%usetra ) then

         iion_global = iion_global + 1
         ielem1 = 1
         newmolecule = .true.

      else

         if( iatomi == ion%natom(ielem) ) then
            if( ielem == nelemIon ) then
               done = ( iion_global == nion )
               iion_global = iion_global + 1
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
!      print *, 'iion_global, ielem, iatomi=', iion, ielem, iatomi
! end-debug

      end if

      iion = iion_global

   end subroutine advanceCounter_single
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
      integer(INT4), intent(out) :: ielem1      ! atom index of projectile
      logical, intent(out) :: newmolecule       ! flag indicating whether projectile is from a
                                                !    new molecule
      logical, intent(out) :: done              ! flag indicating no more ions to simulate
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      lock( iion_lock[1] )

         done = .false.

         if( ion%usetra ) then

            iion_global[1] = iion_global[1] + 1
            ielem1 = 1
            newmolecule = .true.

         else

            if( iatomi == ion%natom(ielem) ) then
               if( ielem == nelemIon ) then
                  done = ( iion_global[1] >= nion )
                  iion_global[1] = iion_global[1] + 1
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
            !print *,  'Image=', this_image(), ' iion=', iion_global[1]
!      print *, 'iion_global, ielem, iatomi=', iion, ielem, iatomi
! end-debug

         end if

         iion = iion_global[1]

      unlock( iion_lock[1] )

   end subroutine advanceCounter
!
!=======================================================================
!
   function lastProjectile( ngroup_, k_ )  result( last )
!           ==============     i     i              o
!
!     Determine if this is the last atom of one of the k_ last projectile molecules out of a
!     group of ngroup_ projectile molecules
!
!     Result is false if ngroup=0
!
!     If usetra=t and ngroup is not given, last=t for the last projectile on the TRA file
!     (unit 123)
!
!  Arguments:
!
      integer(INT4), intent(in), optional :: ngroup_  ! size of group
      integer(INT4), intent(in), optional :: k_       ! the k_ last last-molecules atoms will be
                                                      !    considered "last"
      logical :: last                                 ! flag indicating one of k_ last
                                                      !    projectiles
!
!  Local Variables:
!
      integer(INT4) :: ngroup, k, i, n, iflag, idummy
      real(REAL8) :: rdummy
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( ngroup_ ) ) then
         if( ngroup_ == 0 ) then
            last = .false.
            return
         end if
      end if

      if( present( k_ ) ) then
         k = k_
      else
         k = 1
      end if

      if( ion%usetra ) then

         if( present( ngroup_ ) ) then
            last = ( mod( iion+k-1, ngroup_ ) < k )
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

         if( present( ngroup_ ) ) then
            ngroup = ngroup_
         else
            ngroup = nion
         end if

         last = ( mod( iion+k-1, ngroup ) < k  &
                  .and. ielem == nelemIon .and. iatomi == ion%natom(ielem) )

      end if

   end function lastProjectile
!
!=======================================================================
!
   subroutine getIonIndex( iion_ )
!             ===========    o
!
!     Get the index of the most recent ion
!
!  Arguments:
!
      integer(INT4), intent(out) :: iion_       ! index of most recent ion
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      iion_ = iion

   end subroutine getIonIndex
!
!=======================================================================
!
end module Ion_Count_mod_old

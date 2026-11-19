module Xtal_Atoms_Stat_mod
!
!  Provide an atom list for given lattice sites. Defects are 
!
!     -) generated statistically from the damage histograms
!     -) at random positions
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO
   use Types_mod,          only: Ipos_type, Site_type, Atom_type, NATOMM
   use Geom_CoordSys_mod,  only: fromCoordSys, getVolume

   use Param_Setup_mod,    only: nelem, ldamdyn
   use Param_Xtal_mod,     only: xtal
   use Param_Mat_mod,      only: mat, coord_xtal
   use Param_Dam_mod,      only: dam

   use Random_mod,         only: rannum, ranpo
   use Xtal_Ipos_mod,      only: getPos, operator( == )
   use Xtal_Lattice_mod,   only: getNearestSite
   use Eval_Stat_mod,      only: getDamage
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
   public :: getAtomsStat                    ! get a list of target atoms for given
!            ============                    !    lattice sites
!
!=========================================================================================
!
contains
!
!==========================================================================================
!
   subroutine getAtomsStat( pos, ir, nsite, site, natom, atom, cubeo, cube )
!             ============   i   i     io    io     o     o      i     i   
!
!     Get a list of target atoms for given lattice sites and initialize all components of site
!     and atom for cross-linking.
!
!     The lattice site list may be modified due to random interstitials and vacancies
!     if lamo2=.false.
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)         ! position where the damage histogram is
                                                !    evaluated (target coordinate system)
      integer(INT4), intent(in) :: ir           ! region index
      integer(INT4), intent(inout) :: nsite     ! number of lattice sites
      type(Site_type), intent(inout) :: site(:) ! lattice sites for which atoms are requested
      integer(INT4), intent(out) :: natom       ! number of lattice site atoms
      type(Atom_type), intent(out) :: atom(:)   ! lattice site atoms
      real(REAL8), intent(in) :: cubeo(3,2)     ! old cube boundaries
      real(REAL8), intent(in) :: cube(3,2)      ! new cube boundaries
!
!  Local Variables:
!
      integer(INT4) :: isite, ielemr, ielem, k, ia, isa
      type(Ipos_type) :: ipos
      real(REAL8) :: probv(nelem), histv(nelem), histi(nelem)
      real(REAL8) :: rand, vol, rmean
      real(REAL8) :: posai(3), posa(3), rand1, rand2, rand3
!
!  Start of Executable Code:
!  ------------------------
!
!  First set up a perfect lattice at the given lattice sites
!
      natom = 0

      do isite = 1, nsite
         natom = natom + 1
         if( natom > size( atom ) )  error stop 'getAtomsStat: size of atom array too small.'

         site(isite)%nsa = 1
         site(isite)%iatom(1) = natom
         site(isite)%latticeAtom = .true.

         atom(natom)%isite = isite
         atom(natom)%isa = 1
         call getPos( site(isite)%ipos, ir, atom(natom)%pos(:) )
!             ======        i           i          o
         atom(natom)%lvib = .false.
         atom(natom)%ielem = site(isite)%ielem    ! default species at site
      end do
!   
!  If no random interstitials are to be taken into account, we are done
!
      if( .not.ldamdyn .or. dam%lamo2 ) return
!
!  Calculate probabilities of random interstitials and vacancies
!  (for the moment we neglect atoms recoiled from other regions)
!
      probv(:) = ZERO

      do ielemr = 1, mat(ir)%nelemr
         ielem = mat(ir)%ielem(ielemr)

         call getDamage( pos(:), ielem, histi(ielemr), histv(ielemr) )
!             =========    i       i         o              o

         if( mat(ir)%prob(ielemr) > 0 )  &
               probv(ielem) = histv(ielemr) / ( mat(ir)%prob(ielemr) * mat(ir)%dens )
      end do
!
!  Create vacancies randomly
!
      do isite = 1, nsite
         call rannum( rand )
!             ======   o
         if( rand < probv(site(isite)%ielem) ) then   ! vacancy
            site(isite)%latticeAtom = .false.
            site(isite)%nsa = 0
         end if
      end do
!
!  Pack atom list to exclude vacancies (vacancy sites are not eliminated)
!
      natom = 0

      do isite = 1, nsite
         do isa = 1, site(isite)%nsa
            natom = natom + 1
            if( natom > size( atom ) ) error stop 'getAtomsStat: size of atom too small.'
            atom(natom) = atom(site(isite)%iatom(isa))
            site(isite)%iatom(isa) = natom
         end do
      end do                  !  now natom is the number of atoms in the list with vacancies !
!
!  Create random interstitials within the new search cube but outside the old search cube
!
      call getVolume( coord_xtal, vol )
!          =========        i         o
      vol = vol * product( cube(:,2)-cube(:,1) )

      do ielemr = 1, mat(ir)%nelemr
         rmean = histi(ielemr) * vol

         call ranpo( rmean, k )                 ! number of interstitials is Poisson
!             =====    i    o                   ! distributed

         do ia = 1, k                         

            call rannum( rand1 )              
!                ======    o
            call rannum( rand2 )
!                ======    o
            call rannum( rand3 )
!                ======    o
            posai(:) = cube(:,1) + [rand1,rand2,rand3] * (cube(:,2) - cube(:,1))

            if( all( posai(:) > cubeo(:,1) ) .and. all( posai(:) < cubeo(:,2) ) )  cycle
                                                ! reject if interstitial is inside old cube

            call fromCoordSys( coord_xtal, posai(:), posa(:) )
!                ============        i          i         o
!
!  Find nearest lattice site position and determine if it is already in the list
!
            call getNearestSite( posa(:), ir, ipos )
!                ==============      i    i    o

            do isite = 1, nsite
               if( site(isite)%ipos == ipos )  exit
            end do                                 ! if not in list, isite=nsite+1
!
!  Check whether there is enough space to store this atom
!
            if( natom >= size( atom ) )  return        ! not enough space for additional atoms

            if( isite == nsite + 1 ) then
               if( isite > size( site ) )  return      ! not enough sites for additional site
            else
               if( site(isite)%nsa == NATOMM )  cycle  ! not enough space at the existing site
            end if
!
!  If the site does not exist yet, create a new one
!
            if( isite == nsite + 1) then
               nsite = isite
               site(isite)%ipos = ipos
               if( xtal%uc%name == ' ' ) then
                  site(isite)%ielem = mat(ir)%ielem(1)
               else
                  site(isite)%ielem = mat(ir)%ielem(xtal%uc%ielemr(site(isite)%ipos%iucs))
               end if
               site(isite)%nsa = 0               ! atom will be assigned below
            end if
!
!  Assign atom to site
! 
            natom = natom + 1

            site(isite)%nsa = site(isite)%nsa + 1
            site(isite)%latticeAtom = .false.
            site(isite)%iatom(site(isite)%nsa) = natom

            atom(natom)%isite = isite
            atom(natom)%isa = site(isite)%nsa
            atom(natom)%pos(:) = posa(:)
            atom(natom)%lvib = .false.
            atom(natom)%ielem = mat(ir)%ielem(ielemr)

         end do

      end do

   end subroutine getAtomsStat
!
!=============================================================================================
!
end module Xtal_Atoms_Stat_mod

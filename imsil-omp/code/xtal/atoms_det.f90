module Xtal_Atoms_Det_mod
!
!  Provide an atom list for given lattice sites. Defects are
!
!     -) deterministic, i.e. their types and positions do not change during the simulation
!     -) of certain types defined by the positions of their atoms relative to the lattice
!        sites 
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO
   use Types_mod,             only: Site_type, Atom_type, LSDefect_type

   use Param_Setup_mod,       only: nelem
   use Param_Mat_mod,         only: coord_xtal

   use Eval_Det_mod,          only: det60
   use Geom_CoordSys_mod,     only: fromCoordSys
   use Random_mod,            only: ranint
   use Target_mod,            only: getDomainDensities
   use Xtal_Ipos_mod,         only: getPos
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
   public :: getAtomsDet                   ! get a list of target atoms for given
!            ===========                   !   lattice sites
!
!===========================================================================================
!
contains
!
!============================================================================================
!
   subroutine getAtomsDet( ir, site, natom, atom )
!             ===========  i    io     o     o
!
!     Get a list of target atoms for given lattice sites and initialize all component
!     for cross-linking
!
!  Arguments:
!
      integer(INT4), intent(in) :: ir           ! region index
      type(Site_type), intent(inout) :: site(:) ! lattice sites
                                                !    site(:)%ipos is (input), all other
                                                !    components of site(:) are (output)
      integer(INT4), intent(out) :: natom       ! number of lattice site atoms
      type(Atom_type), intent(out) :: atom(:)   ! lattice site atoms
!
!  Local Variables:
!
      logical :: isDefect

      integer(INT4) :: nsite, isite, iatom, isa, ielemr, nelemr, ielem(nelem)
      integer(INT4) :: jpos(3)
      real(REAL8) :: posi(3), prob(nelem)
      type(LSDefect_type) :: lsdef
!
!  Start of Executable Code:
!  ------------------------
!
      nsite = size( site )

      iatom = 0
!
!  Loop over all lattice sites
!
      do isite = 1, nsite

         jpos(:) = site(isite)%ipos%ipos(:)

         call det60( jpos(:), ir, lsdef, isDefect )   ! check whether there is a defect and
!             =====    i      i     o       o         !    if so what

         if( isDefect ) then                    ! there is a defect

            site(isite)%latticeAtom = .false.

            site(isite)%nsa = lsdef%nsa

            do isa = 1, site(isite)%nsa
               iatom = iatom + 1
               if( iatom > size( atom ) )  error stop 'getAtomsDet: size of atom too small.'
               site(isite)%iatom(isa) = iatom

               atom(iatom)%isite = isite
               atom(iatom)%isa = isa
               posi(:) = site(isite)%ipos%ipos(:) + lsdef%dpos(:,isa)
               call fromCoordSys( coord_xtal, posi(:), atom(iatom)%pos(:) )
!                   ============        i          i              o
               atom(iatom)%ielem = lsdef%ielem(isa)
               atom(iatom)%lvib = .false.
            end do

         else                                   ! no defect, only a lattice atom

            site(isite)%latticeAtom = .true.

            iatom = iatom + 1
            if( iatom > size( atom ) )  error stop 'getAtomsDet: size of atom too small.'

            site(isite)%nsa = 1
            site(isite)%iatom(1) = iatom

            atom(iatom)%isite = isite
            atom(iatom)%isa = 1
            call getPos( site(isite)%ipos, ir, atom(iatom)%pos(:) )
!                ======          i                     o
            atom(iatom)%lvib = .false.

            call getDomainDensities( ir, nelemr, ielem(:), prob_=prob(:) )
!                ==================  i     o       o           o
            call ranint( prob(1:nelemr), ielemr )
!                ======        i           o
            atom(iatom)%ielem = ielem( ielemr )

         end if

      end do

      natom = iatom

   end subroutine getAtomsDet
!
!=========================================================================================
!
end module Xtal_Atoms_Det_mod

module Xtal_Atoms_mod
!
!  Provide a list of target atoms around a given position
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    HALF, ONE, ZERO, TWO
   use Types_mod,             only: Ipos_type, Site_type, Atom_type

   use Param_Dam_mod,         only: dam
   use Param_Mat_mod,         only: mat, coord_xtal
   use Param_Xtal_mod,        only: xtal

   use Geom_CoordSys_mod,     only: toCoordSys
   use Xtal_Ipos_mod,         only: operator( == )
   use Xtal_Lattice_mod,      only: maxSizeLattice, getSites
   use Xtal_Atoms_Det_mod,    only: getAtomsDet
   use Xtal_Atoms_Stat_mod,   only: getAtomsStat
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
   public :: getAtoms                        ! get a list of targets atoms and their
!            ========                        !   corresponding lattice sites. If cell 
!                                            !   overlaps with previous cell, atoms remain
!                                            !   the same in the overlap area.
!
!======================================================================================
!
contains
!
!======================================================================================
!
   subroutine getAtoms( pos, ir, nsite, site, natom, atom )
!             ========   i   i    io     io    io     io
!
!     Get a list of target atoms and their corresponding lattice sites around a given
!     position
!
!     -) If no new unit cube is entered, nothing is done
!     -) If a new unit cube is entered with partly identical sites to
!        the old one, the occupation of the sites belonging to both
!        cubes will be the same
!     -) Displacements due to lattice vibrations are not considered at
!        the time of generation of an atom in order to save the cpu time 
!        required to displace atoms never involved in any search procedure
!
!  Arguments:
!
      real(REAL8), intent(in) :: pos(3)         ! position (target coordinate system)
      integer(INT4), intent(in) :: ir           ! region index
      integer(INT4), intent(inout) :: nsite     ! number of lattice sites
      type(Site_type), intent(inout) :: site(:) ! lattice sites
      integer(INT4), intent(inout) :: natom     ! number of target atoms
      type(Atom_type), intent(inout) :: atom(:) ! target atoms
!                                                
!  Local Variables:
!
      logical, save :: first = .true.
      integer(INT4), save :: icube(3)           ! Note: needs attention when extending to more
      real(REAL8), save :: cube(3,2)            !       than one crystalline region !

      logical :: newcube
      integer(INT4) :: icubeo(3)
      real(REAL8) :: cubeo(3,2)

      integer(INT4) :: nsite_new
      type(Ipos_type), allocatable, save :: newsite_ipos(:)
      logical, allocatable, save :: newsite_flag(:)

      integer(INT4) :: isite, isite_new, jsite, nsite_old, iatom, jatom, natom_new, isa

      real(REAL8) :: posi(3)
!$omp threadprivate( first, icube, cube, newsite_ipos, newsite_flag )
! debug
!      integer(INT4), save :: counter = 0
!      character(6) :: fname
! end-debug
!
!  Start of Executable Code:
!  ------------------------
!
!  Convert position into the crystal coordinate system
!
      call toCoordSys( coord_xtal, pos(:), posi(:) )
!          ==========        i          i        o
!
!  Determine whether the unit cube is the same as before
!                         
      if( first ) then
         if( xtal%uc%name == ' ' ) then
            icube(:) = nint( HALF*posi(:) - HALF )
         else
            icube(:) = nint( posi(:) - HALF )
         end if
         if( .not. allocated( newsite_ipos ) )  &
               allocate( newsite_ipos( maxsizeLattice() ) )
!                                      ==============
         if( .not. allocated( newsite_flag ) )  &
               allocate( newsite_flag( maxsizeLattice() ) )
!                                      ==============
         first = .false.
      else      
         icubeo(:) = icube(:)
         if( xtal%uc%name == ' ' ) then
            icube(:) = nint( HALF*posi(:) - HALF )
         else
            icube(:) = nint( posi(:) - HALF )
         end if
         newcube = any( icube(:) /= icubeo(:) )
                                                ! if this is not a new cube, nothing has to
         if( .not.newcube )  return             !    be done

         cubeo(:,:) = cube(:,:)                 ! save old search cube boundaries (needed
                                                !    for random interstitials)
      end if
!
!  Get new lattice site list
!
      call getSites( posi(:), nsite_new, newsite_ipos(:), cube(:,:) )
!          ========    i          o            o              o
!
!  Label lattice sites that are identical in both lists
!
      site(1:nsite)%flag = .false.
      newsite_flag(1:nsite_new) = .false.

      do isite = 1, nsite
         do isite_new = 1, nsite_new
            if( site(isite)%ipos == newsite_ipos(isite_new) ) then
               site(isite)%flag = .true.
               newsite_flag(isite_new) = .true.
               exit
            end if
         end do
      end do
!
!  Pack site
!
      jsite = 0

      do isite = 1, nsite
         if( site(isite)%flag ) then
            jsite = jsite + 1
            site(jsite) = site(isite)
            do isa = 1, site(isite)%nsa
               atom(site(isite)%iatom(isa))%isite = jsite
            end do
         end if
      end do

      nsite = jsite
      nsite_old = nsite
!
!  Pack atom
!
      jatom = 0

      do isite = 1, nsite
         do isa = 1, site(isite)%nsa
            jatom = jatom + 1
            atom(jatom) = atom(site(isite)%iatom(isa))
            site(isite)%iatom(isa) = jatom
         end do
      end do

      natom = jatom
!
!  Append new lattice sites to site
!
      jsite = nsite

      do isite_new = 1, nsite_new
         if( newsite_flag(isite_new) )  cycle
         jsite = jsite + 1
         if( jsite > size( site ) ) then
            print *,  jsite, size( site ), 'getAtoms: size of site too small'
            error stop
         end if
         site(jsite)%ipos = newsite_ipos(isite_new)
         site(jsite)%flag = .false.
                                                ! index of element normally on this site
         if( xtal%uc%name == ' ' ) then
            site(jsite)%ielem = mat(ir)%ielem(1)
         else
            site(jsite)%ielem = mat(ir)%ielem(xtal%uc%ielemr(site(jsite)%ipos%iucs))
         end if
      end do

      nsite = jsite
!
!  Determine lattice site atoms at the new lattice sites
!
      if( dam%ldet ) then
         call getAtomsDet( ir, site(nsite_old+1:nsite), natom_new, atom(natom+1:) )
!             ===========  i           io                     o             o
      else
         nsite_new = nsite - nsite_old
         call getAtomsStat( pos(:), ir, nsite_new, site(nsite_old+1:),  &
!             ============    i      i     io              io
                            natom_new, atom(natom+1:), cubeo(:,:), cube(:,:) )
!                               o            o             i           i
         nsite = nsite_old + nsite_new
      end if
!                                               ! shift pointers
      do isite = nsite_old+1, nsite
         do isa = 1, site(isite)%nsa
            site(isite)%iatom(isa) = site(isite)%iatom(isa) + natom
         end do
      end do

      do iatom = natom+1, natom+natom_new
         atom(iatom)%isite = atom(iatom)%isite + nsite_old
      end do

      natom = natom + natom_new
! debug
!      if( counter < 5 ) then
!         counter = counter + 1
!         write( fname, '(a, i1)' ) 'SITES', counter
!         open( 487, file=fname, form='formatted', status='replace' )
!         do isite = 1, nsite
!            write( 487, '(5i3)' )  site(isite)%ipos%iuc(:)-icube(:), site(isite)%ipos%iucs,  &
!                                   site(isite)%ielem
!         end do
!         close( 487 )
!      end if
! end-debug
!
   end subroutine getAtoms
!
!=======================================================================
!
end module Xtal_Atoms_mod

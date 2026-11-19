module Casc_Selncc_mod
!
!  Select the next collision partner in crystalline material
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 EPSMA, SMALL, ZERO, ONE, THREE, HALF,  &
                                 STATE
   use Vector_mod,         only: operator(.dot.), operator(.cross.), abs, abs2
   use Types_mod,          only: Proj_type, Site_type, Atom_type, NATOMM

   use Param_Setup_mod,    only: ldamdyn
   use Param_Geom_mod,     only: geom
   use Param_Xtal_mod,     only: xtal
   use Param_Mat_mod,      only: mat
   use Param_Dam_mod,      only: dam
   use Param_Sn_mod,       only: sn

   use Target_mod,         only: getDomainIndex
   use Target_Surface_mod, only: getDistSurface
   use Random_mod,         only: rannum
   use Sn_Ndira_mod,       only: ndira
   use Eval_Stat_mod,      only: getDamage
   use Xtal_Ipos_mod,      only: getPos
   use Xtal_Lattice_mod,   only: maxSizeLattice
   use Xtal_Atoms_mod,     only: getAtoms
   use Xtal_Tvib_mod,      only: tvib2, tvib3, addFakeVib
   use Xtal_OldSites_mod,  only: testOldSites
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
   public :: selncc                 ! select next collision partner in
!            ======                 !   crystalline material
!
!  Private Data:
!  ============
!
   logical :: first = .true.
                                                ! arrays for lattice atoms or first atom
                                                !    of split configuration

   integer(INT4) :: nsite                       ! number of lattice sites
   integer(INT4) :: natom                       ! number of lattice site atoms

   type(Site_type), allocatable :: site(:)      ! (NLSMAX) lattice sites
   type(Atom_type), allocatable :: atom(:)      ! (NATOMM*NLSMAX) lattice site atoms
!$omp threadprivate( first, nsite, natom, site, atom )
!
contains
!
!=======================================================================
!
   subroutine selncc( proj, ir, posc, ffp, pmax, lrecoil, recoil, p )
!             ======   i    i    o     o    o       o       o     o
!
!     Select next collision partner in crystalline material
!
!     1) the next collision partner must
!        a) be among the stored lattice sites, i.e. those lattice sites
!           which are closer than sqrt(nrad2+1) to the unit-cube 
!           containing x.
!        b) be within the minimum and the maximum free flight path.
!        c) lead to an impact parameter smaller than the maximum
!           impact parameter.
!        d) not be identical to previous collision partners stored in
!           the old site list (see subroutine oslist)
!     2) thermal vibrations are considered by a random displacement of 
!        the collision partner from its lattice site in the plane
!        perpendicular to the flight direction (lamo2=t). The probability
!        density function is Gaussian.
!        thermal vibrations are not considered during the search
!        procedure.
!
!  Arguments:
!
      type(Proj_type), intent(in) :: proj       ! projectile properties
      integer(INT4), intent(in) :: ir           ! index of region where the collision partner
                                                !    is to be selected
      real(REAL8), intent(out) :: posc(3)       ! position at which collision takes place 
                                                !    (= where the impact parameter is 
                                                !    perpendicular to the direction of 
                                                !    motion) or end of free flight path
      real(REAL8), intent(out) :: ffp           ! free flight path (= distance between x and
                                                !    xc)
      real(REAL8), intent(out) :: pmax          ! maximum impact parameter
      logical, intent(out) :: lrecoil           ! indicates whether a collision partner
                                                !    (recoil) has been found
      type(Proj_type), intent(out) :: recoil    ! recoil
      real(REAL8), intent(out) :: p             ! impact parameter (= distance between xr and
!                                               !    xc)
!  Local Variables:
!
      logical :: found, lamo2

      integer(INT4) :: ielem, ielemr, iatom, jatom
      integer(INT4) :: n
      integer(INT4) :: irl

      real(REAL8) :: ffpmin, ffpmax, ffplim
      real(REAL8) :: vib
      real(REAL8) :: probi
      real(REAL8) :: histi, histv, rand, ffpadd, p2l, ffpl, ffp1
      real(REAL8) :: cospsi, sinpsi
      real(REAL8) :: dir1(3)
      real(REAL8) :: postry(3), poscl(3)
      real(REAL8) :: x, y, distsurf
!
!  Start of Executable Code:
!  ------------------------
!
      lamo2=.false.

      if( first ) then
         n = maxSizeLattice( )
!            ==============
         if( .not. dam%ldet .and. .not.lamo2 ) then   ! actual number of sites varies
            n = n + THREE * sqrt( real(n) )           !    statistically => increase space
         end if

         allocate( site( n ) )
         allocate( atom( NATOMM * n ) )
         first = .false.
      end if

      if( xtal%uc%name == ' ' ) then
         ffplim = 4 * mat(ir)%xscale
      else
         ffplim = mat(ir)%xscale
      end if

      pmax = sn%pmaxc(ir)
      ffpmin = proj%ffpmin + EPSMA
      ffpmax = min( xtal%ffpmax*mat(ir)%xscale, max( ZERO, proj%distsurf )+sn%pmax )
!
!  Some auxiliary variables in order to handle the case when no collision
!  partner is found in the first search            
!
      ffpadd = ZERO
      postry(:) = proj%pos(:)
!
!  Loop over attempts to find a collision partner
!  ==============================================
!
      Search_loop: do               
!
!  Update the lattice and interstitial site list
!
         call getAtoms( postry(:), ir, nsite, site(:), natom, atom(:) )
!             ========    i        i    io     io       io     io
!
!debug
!         open( 113, file='POS' )
!         do iatom=1,natom
!            write( 113, '(3f6.2)' ) atom(iatom)%pos(:)
!         end do
!         close( 113 )
!         error stop
!end debug
! 
!  Loop over atoms
!  ---------------
!
!  Search for site with appropriate ffp, p, and which are not
!  identical to one of the sites in the old site list
!
         lrecoil = .false.
         ffp = ffpmax

         Atoms_loop: do iatom = 1, natom
                                                !
                                                ! modify coordinates due to thermal 
                                                !    vibrations if necessary
                                                !
            if( .not. (xtal%lvib2 .or. atom(iatom)%lvib) ) then
                                                ! make a rough check if atom may be within pmax
               p2l = abs2( proj%dir(:) .cross. (atom(iatom)%pos(:) - postry(:)) )
!                    ====
               ielem = atom(iatom)%ielem
               vib = mat(ir)%vib( ielem )
               call addFakeVib( vib, proj%e, ielem )
!                   ==========  io     i       i
               if( p2l > (pmax + 2*vib)**2 )  cycle Atoms_loop

               call tvib3( vib, atom(iatom)%pos(:) )
!                   =====   i          io
               atom(iatom)%lvib = .true.
            end if
                                                !
                                                ! check ffpl
                                                !
            ffpl = (atom(iatom)%pos(:) - postry(:)) .dot. proj%dir(:)
            if( ffpl < ffpmin )  cycle Atoms_loop
            if( ffpl > ffp )  cycle Atoms_loop
            poscl(:) = postry(:) + ffpl*proj%dir(:)
                                                !
                                                ! check p2l
                                                !
            p2l = abs2( poscl(:) - atom(iatom)%pos(:) )
!                 ====
            if( p2l > pmax**2 )  cycle Atoms_loop
                                                ! 
                                                ! check old site list
                                                !
            call testOldSites( site(atom(iatom)%isite)%ipos, atom(iatom)%isa, proj%ig, found )
!                ============                i                       i           i       o
            if( found )  cycle Atoms_loop
                                                !
                                                ! check if atom is in appropriate region
                                                !
            if( geom%lrobin ) then              ! for round robin the target is laterally 
                                                !    limited (NOTICE: the boundary must be 
                                                !    the same as is rrobin)
               if( xtal%uc%name == ' ' ) then
                  if( site( atom(iatom)%isite )%ipos%ipos(1) < -10 .or.  &
                     site( atom(iatom)%isite )%ipos%ipos(1) > 10 .or.  &
                     site( atom(iatom)%isite )%ipos%ipos(2) < -10 .or.  &
                     site( atom(iatom)%isite )%ipos%ipos(2) > 10 )  cycle Atoms_loop
               else if( xtal%uc%name == 'diamond' .or. xtal%uc%name == 'diamond2' ) then
                  x = site( atom(iatom)%isite )%ipos%iuc(1) +  &
                     xtal%uc%pos(1,site( atom(iatom)%isite )%ipos%iucs)
                  y = site( atom(iatom)%isite )%ipos%iuc(2) +  &
                     xtal%uc%pos(2,site( atom(iatom)%isite )%ipos%iucs)
                  if( x < -2.625 .or. x > 2.625 .or. y < -2.625 .or. y > 2.625 )  &
                        cycle Atoms_loop
               else
                  error stop 'selncc: unknown crystal name with lrobin=t'
               end if
            end if
            irl = ir
            call getDomainIndex( atom(iatom)%pos(:), irl )
!                ==============         i            io
            call getDistSurface( atom(iatom)%pos(:), distsurf )
!                ==============         i               o
            if( irl /= ir .or. distsurf < ZERO )  cycle Atoms_loop
                                                !
                                                ! if all tests have been passed, preliminarily
                                                !    accept atom as collision partner
                                                !
            lrecoil = .true.
            jatom = iatom
            ffp = ffpl
            posc(:) = poscl(:)
            recoil%id = irl

         end do Atoms_loop
!
!  End of loop over atoms
!  ----------------------
!
         if( lrecoil ) then

            exit Search_loop

         else
!
!  If no collision partner has been found, repeat search with new ion position
!  unless free flight path is already larger than ffplim
!
            ffp1 = ffpmax - EPSMA
            postry(:) = postry(:) + ffp1*proj%dir(:)
            ffpadd = ffpadd + ffp1
            ffpmin = ffpmin - ffp1
            if( ffpadd <= ffplim ) then
               cycle Search_loop
            else
               ffp = ffpadd
               posc(:) = proj%pos(:) + ffp*proj%dir(:)
               return
            end if

         end if

      end do Search_loop
!
!  End of loop over attempts to find a collision partner
!  =====================================================
!
!  ( here we arrive only if a recoil has been found )
!
!  Copy remaining atom properties to recoil
!
      recoil%origin%ipos = site( atom(jatom)%isite )%ipos
      call getPos( recoil%origin%ipos, ir, recoil%origin%pos(:) )
!          ======          i           i          o
      recoil%origin%internal = .true.
      recoil%origin%amorphous = .false.
      recoil%origin%lattice = site( atom(jatom)%isite )%latticeAtom
      recoil%origin%isa = atom(jatom)%isa

      recoil%pos(:) = atom(jatom)%pos(:)
      recoil%ielem = atom(jatom)%ielem
!
!  Adjust free flight path
!
      ffp = ffp + ffpadd
!
!  Calculate recoil coordinates taking 2-d amorphization model or thermal 
!  vibrations into account
!
      lamo2 = .false.

      if( ldamdyn .and. recoil%origin%lattice .and. dam%lamo2 .and. .not. dam%ldet ) then
                                                ! determine if 2d random interstitial is hit
         call rannum( rand )                    !   (lamo2) and which atom it is (ielem)
!             ======   o
         do ielemr = 1, mat(ir)%nelemr
            ielem = mat(ir)%ielem(ielemr)

            call getDamage( proj%pos(:), ielem, histi, histv )
!                =========      i          i      o      o
            probi = histi / mat(ir)%dens

            if( rand < probi )  exit
            rand = rand - probi
         end do

         lamo2 = ( ielemr <= mat(ir)%nelemr )

      end if

      if( lamo2 ) then
         recoil%ielem = ielem                   ! overwrite atom species index
                              
         call rannum( rand )                    ! get position of 2d random interstitial
!             ======   o

         p = sqrt(rand) * pmax
         cospsi = ZERO
         sinpsi = ONE
         call ndira( proj%dir(:), cospsi, sinpsi, dir1(:) )
!             =====      i          i       i        o
         recoil%pos(:) = posc(:) + p*dir1(:)
         recoil%origin%pos(:) = recoil%pos(:)

         vib = mat(ir)%vib(recoil%ielem)
         call addFakeVib( vib, proj%e, recoil%ielem )
!             ==========  io     i          i
         call tvib2( vib, proj%dir(:), recoil%pos(:) )
!             =====   i        i            io
                                                ! 2d random interstitial does not count 
         recoil%origin%lattice = .false.        !    as lattice site

      else
                                                ! apply 2d lattice vibrations
         if( xtal%lvib2 ) then

            vib = mat(ir)%vib(recoil%ielem)
            call addFakeVib( vib, proj%e, recoil%ielem )
!                ==========  io     i          i
            call tvib2( vib, proj%dir(:), recoil%pos(:) )
!                =====   i        i            io
         end if
         p = abs( recoil%pos(:) - posc(:) )
!            ===

      end if
!
!  Check if recoil has been moved out of target by the lattice vibrations
!  Note: lrecoil would be set to zero anyway in selncp if distsurf<0
!
      call getDistSurface( recoil%pos(:), recoil%distsurf )
!          ==============       i                 o
      if( recoil%distsurf < ZERO ) then
         lrecoil = .false.
         return
      end if
!
!  Initialize remaining properties of recoil
!  (recoil%e and recoil%dir(:) will be initialized in scatter)
!
      recoil%ffpmin = ZERO
      recoil%emin = recoil%e

      recoil%ig = proj%ig + 1
      recoil%iion = proj%iion
      recoil%nreflect = 0
      recoil%reentered = proj%reentered

      recoil%virtualparent = proj%virtual
      recoil%virtual = proj%virtual             ! virtual may be modified later
      recoil%ededuct = ZERO
      recoil%state = STATE%SUBTHRESHOLD

   end subroutine selncc
!
!=======================================================================
!
end module Casc_Selncc_mod

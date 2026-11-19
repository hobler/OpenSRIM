module Ion_Proj_mod
!
!  Manage data about the ions. The ions may be atoms to be implanted, molecules, or recoils 
!  starting from the interior of the target.
!
!  Terms: Ion =        Atom/molecule/recoil to be implanted or started, 
!                      before interaction with the target
!         Projectile = Particle moving in the target 
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: REAL8, INT4,  &
                                    ZERO, ONE, TWO, HALF, PI, SMALL, BIG,  &
                                    STATE
   use Undefine_mod,          only: defined
   use Types_mod,             only: Proj_type, Ipos_type

   use Param_Setup_mod,       only: nelemIon, nelem, nthread
   use Param_Geom_mod,        only: geom
   use Param_Elem_mod,        only: elem
   use Param_Ion_mod,         only: ion
   use Param_Xtal_mod,        only: xtal
   use Param_Mat_mod,         only: mat, coord_xtal, irxtal
   use Param_Dam_mod,         only: dam
   use Param_Sn_mod,          only: sn

   use Casc_Surface_mod,      only: enterTarget
   use Eval_Det_mod,          only: getCol
   use Geom_CoordSys_mod,     only: CoordSys_type, initCoordSys, toCoordSys, fromCoordSys
   use Random_mod,            only: rannum, ranno, ranint
   use Sn_FFPcyl_mod,         only: newFFPcyl
   use Target_mod,            only: getDomainIndex
   use Target_Surface_mod,    only: getDistSurface
   use Target_Surface_mod,    only: getEntranceGrid
   use Xtal_Ipos_mod,         only: getPos
   use Xtal_Lattice_mod,      only: getNearestSite
   use Xtal_Tvib_mod,         only: tvib3
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
   public :: getProjectile          ! get new projectile from ions
!            =============          !
   public :: getXStart              ! get the most recently initialized projectile (could
!            =========              !    eventually be renamed to "getRecentProjectile")
!
!  Private data:
!  ============
!
   type(Proj_type) :: proj_start
   logical :: proj_valid
!$omp threadprivate( proj_start, proj_valid )

contains
!
!==============================================================================
!
   subroutine getProjectile( ielem1, newmolecule, proj, found, done )
!             =============    i          i        o      o     o
!
!     Initialize a projectile
!
!     proj%id is not initialized if no projectile is found
!
!  Arguments:
!
      integer(INT4), intent(in) :: ielem1       ! element index
      logical, intent(in) :: newmolecule        ! flag indicating the projectile belongs to a
                                                !    new molecule
      type(Proj_type), intent(out) :: proj      ! projectile
      logical, intent(out) :: found             ! flag indicating that valid projectile has been
                                                !    found
      logical, intent(out) :: done              ! flag indicating no more ions to simulate
!
!  Local Variables:
!
      integer(INT4) :: ibeam
      integer(INT4) :: ielem2
      integer(INT4) :: idummy, iflag
      integer(INT4) :: icolmin(2), icolmax(2)
      integer(INT4) :: istat

      type(Ipos_type) :: ipos

      real(REAL8) :: dirx, diry
      real(REAL8) :: posi(3)
      real(REAL8) :: rand, rand1, rand2, rand3
      real(REAL8) :: rotate
      real(REAL8) :: psi, sinpsi, cospsi
      real(REAL8) :: sinphi, cosphi
      real(REAL8) :: wcol(2)

      type(CoordSys_type) :: coord_beam
!
!  Start of Executable Code:
!  ------------------------
!
      done = .false.
!
!  If this is another atom of the same ion molecule, take initial conditions from last ion atom
!  (overwrite ielem and e below...)
!
      if( .not.newmolecule ) then
         found = proj_valid
         proj = proj_start
      end if
!
!  Initialize atom index and energy
!
      if( ion%ltinit .and. any( elem(nelemIon+1:nelem)%z == elem(ielem1)%z ) ) then
         do ielem2 = nelemIon+1, nelem
            if( elem(ielem2)%z == elem(ielem1)%z )  exit
         end do
         proj%ielem = ielem2                    ! ion==target: use target
      else
         proj%ielem = ielem1
      end if

      proj%e = ion%e * elem( proj%ielem )%m / ion%m
!
!  If this is another atom of the same ion molecule, we are more or less done
!
      if( newmolecule ) then
!
!  Some default values that may be modified later
!
         proj%origin%amorphous = .false.
         proj%origin%lattice = .false.
!
!  Read starting point from the trajectory file (overwriting proj%e and proj%ielem)
!
         if( ion%usetra ) then

            do
               read( 123, *, iostat=istat )  proj%pos(1), proj%pos(2), proj%pos(3),  &
                                             proj%dir(1), proj%dir(2), proj%dir(3),  &
                                             proj%e, proj%ielem, idummy, iflag
               if( istat < 0 ) then
                  done = .true.
                  return
               end if
               if( proj%ielem > nelem )  &
                     error stop 'getProjectile: element index exceeds maximum number of elements'
               ! Skip projectiles with iflag=1 and iflag=2
               select case( iflag )
                  case( 1 )                     ! update old sites with started recoils
                     if( irxtal /= 0 ) then     !    ignoring the recoil generation index
                        call getNearestSite( proj%pos(:), irxtal, ipos )
!                            ==============      i          i      o
                        call getPos( ipos, irxtal, proj%pos(:) )
!                            ======   i      i        o
                        proj%origin%lattice = .true.
                        proj%origin%ipos = ipos
                        proj%origin%isa = 1
                     end if
                     cycle
                  case( 2 )
                     cycle
                  case default
                     exit
               end select
            end do
            !
            !  Move starting point to the closest lattice site if llinit=t
            !
            if( irxtal /= 0 .and. ion%llinit ) then   ! position of the closest lattice atom

               call getNearestSite( proj%pos(:), irxtal, ipos )
!                   ==============      i          i      o
               call getPos( ipos, irxtal, proj%pos(:) )
!                   ======   i      i        o
               proj%origin%lattice = .true.
               proj%origin%ipos = ipos
               proj%origin%isa = 1

               if( proj%ielem > nelemIon ) then
                  call tvib3( mat(irxtal)%vib(proj%ielem), proj%pos(:) )
!                      =====              i                   io
               end if

            end if
!
!  Select the starting point along a circle
!
         else if( ion%lrinit ) then

            call rannum( rand )
!                ======   o
            rotate = 360*rand
            proj%dir(1) = sin(ion%tilty*PI/180) * cos(rotate*PI/180)
            proj%dir(2) = sin(ion%tilty*PI/180) * sin(rotate*PI/180)
            proj%dir(3) = cos(ion%tilty*PI/180)
            proj%pos(1) = ion%xinit(1) + ion%rinit * cos(rotate*PI/180)
            proj%pos(2) = ion%yinit(1) + ion%rinit * sin(rotate*PI/180)
            proj%pos(3) = ion%zinit(1)
!
!  Select the starting point along a line along rotate
!
         elseif( ion%lpinit ) then

            call rannum( rand )
!                ======   o
            proj%dir(1) = sin(ion%tilty*PI/180) * cos(ion%rotate*PI/180)
            proj%dir(2) = sin(ion%tilty*PI/180) * sin(ion%rotate*PI/180)
            proj%dir(3) = cos(ion%tilty*PI/180)
            proj%pos(1) = ion%xinit(1) - ion%rinit * sin(ion%rotate*PI/180)  &
                          + rand * ion%dinit * cos(ion%rotate*PI/180)
            proj%pos(2) = ion%yinit(1) + ion%rinit * cos(ion%rotate*PI/180)  &
                          + rand * ion%dinit * sin(ion%rotate*PI/180)
            proj%pos(3) = ion%zinit(1)
!
!  All other cases
!
         else
!
!  Beam coordinate system
!
            if( ion%ranrot ) then
               call rannum( rand )
!                   ======   o
               rotate = 360*rand
            else
               rotate = ion%rotate
            end if
            call initCoordSys( coord_beam, rotate_=rotate, tiltx_=ion%tiltx, tilty_=ion%tilty )
!                ============      o              i                 i             i
!
!  Initial direction
!
            call ranint( ion%weight(:), ibeam )
!                ======        i          o

            if( all( ion%diverg(:,ibeam) == ZERO ) ) then

               proj%dir(:) = [ ZERO, ZERO, ONE ]

            else                                ! consider beam divergence

               do
                  select case( ion%moddiv(ibeam) ) ! neglecting the curvature of the sphere
                     case( 1 )                     ! isotropic within ellipse
                        do
                           call rannum( rand1 )
!                               ======    o
                           call rannum( rand2 )
!                               ======    o
                           rand1 = TWO*rand1 - ONE
                           rand2 = TWO*rand2 - ONE
                           if( rand1**2 + rand2**2 < ONE ) exit
                        end do
                        dirx = ion%diverg(1,ibeam)*PI/180 * rand1
                        diry = ion%diverg(2,ibeam)*PI/180 * rand2
                     case( 2 )                  ! Gaussian
                        call ranno( rand1 )
!                            =====    o
                        call ranno( rand2 )
!                            =====    o
                        dirx = ion%diverg(1,ibeam)*PI/180 * rand1
                        diry = ion%diverg(2,ibeam)*PI/180 * rand2
                     case default
                        error stop 'unknown beam divergence model in getProjectile'
                  end select

                  psi = sqrt( dirx**2 + diry**2 )  ! consider sphere curvature
                  sinpsi = sin( psi )
                  call rannum( rand )
!                      ======   o
                  if( psi < PI .and. rand*psi < sinpsi )  exit
               end do

               cospsi = cos( psi )
               if( psi > SMALL ) then
                  cosphi = dirx / psi
                  sinphi = diry / psi
               else
                  cosphi = ONE
                  sinphi = ZERO
               end if

               proj%dir(:) = [ sinpsi*cosphi, sinpsi*sinphi, cospsi ]

            end if

            call fromCoordSys( coord_beam, proj%dir(:), proj%dir(:) )
!                ============      i           i            o
!
!  Starting point
!
            !
            !  Beam profile
            !
            if( any( ion%sigma(:) > ZERO ) ) then
               call ranno( rand1 )
!                   =====    o
               call ranno( rand2 )
!                   =====    o
               proj%pos(:) = [ ion%sigma(1)*rand1, ion%sigma(2)*rand2, ZERO ]

               call fromCoordSys( coord_beam, proj%pos(:), proj%pos(:) )
!                   ============      i           i            o
            else
               proj%pos(:) = [ ZERO, ZERO, ZERO ]
            end if
            !
            !  Add to focal point
            !
            call rannum( rand1 )
!                ======    o
            call rannum( rand2 )
!                ======    o
            call rannum( rand3 )
!                ======    o
            proj%pos(1) = (ion%xinit(1) + rand1*ion%dxinit) + proj%pos(1)
            proj%pos(2) = (ion%yinit(1) + rand2*ion%dyinit) + proj%pos(2)
            proj%pos(3) = (ion%zinit(1) + rand3*ion%dzinit) + proj%pos(3)
            !
            !  modify starting point in a crystal as requested
            !
            if( irxtal /= 0 ) then

               if( ion%llinit ) then            ! position of the closest lattice atom

                  call getNearestSite( proj%pos(:), irxtal, ipos )
!                      ==============      i          i      o
                  call getPos( ipos, irxtal, proj%pos(:) )
!                      ======   i      i        o
                  proj%origin%lattice = .true.
                  proj%origin%ipos = ipos
                  proj%origin%isa = 1

                  if( proj%ielem > nelemIon ) then
                     call tvib3( mat(irxtal)%vib(proj%ielem), proj%pos(:) )
!                         =====              i                   io
                  end if

               else
                                                ! rotate point to the crystal coordinate system
                  call toCoordSys( coord_xtal, proj%pos, posi )
!                      ==========      i          i       o

                  if( dam%ldet ) then           ! need to distribute ions over intersection
                                                !    of the column with the surface
                     call getCol( icolmin, icolmax )
!                         ======     o        o
                     wcol(1:2) = icolmax(1:2) - icolmin(1:2)

                     call rannum( rand1 )
!                         ======    o
                     call rannum( rand2 )
!                         ======    o
                     posi(1) = posi(1) + rand2*wcol(1)
                     posi(2) = posi(2) + rand1*wcol(2)
                     posi(3) = posi(3) -  &
                               rand1 * wcol(2) * mat(irxtal)%wafer(2) / mat(irxtal)%wafer(3) -  &
                               rand2 * wcol(1) * mat(irxtal)%wafer(1) / mat(irxtal)%wafer(3)

                  else                          ! position inside the unit cube centered around
                                                !    the origin
                     call rannum( rand1 )
!                         ======    o
                     call rannum( rand2 )
!                         ======    o
                     call rannum( rand3 )
!                         ======    o
                     if( xtal%uc%name == ' ' ) then
                        posi(1) = posi(1) + (4*rand1-2)
                        posi(2) = posi(2) + (4*rand2-2)
                        posi(3) = posi(3) + (4*rand3-2)
                     else
                        posi(1) = posi(1) + (rand1-0.5)
                        posi(2) = posi(2) + (rand2-0.5)
                        posi(3) = posi(3) + (rand3-0.5)
                     end if

                  end if
                                                ! rotate point to the target coordinate system
                  call fromCoordSys( coord_xtal, posi, proj%pos )
!                      ============          i           i       o
               end if

            end if

         end if
!
!  Perform projection to the surface (if desired) and determine region
!
         call getEntranceGrid( proj, found, onlyforward_=ion%liinit )
!             ===============   io     o              i

         if( found ) then
            call getDistSurface( proj%pos(:), proj%distsurf )  ! we need distsurf in enterTarget
!                ==============      i              o
            if( proj%distsurf < - sn%pmax ) then
               call enterTarget( proj, found )
!                   ===========   io     o
            else
               proj%id = 1    ! guess
               call getDomainIndex( proj%pos(:), proj%id )
!                   ==============      i           io
            end if
         end if
!
!  Initialize the first ffp cylinder
!
         call newFFPcyl( proj%pos(:), ig=1 )
!             =========      i         i

      end if
!
!  Complete the projectile structure
!
      proj%ig = 1
      proj%nreflect = 0
      proj%reentered = .false.

      proj%origin%e = proj%e
      proj%origin%pos(:) = proj%pos(:)
      proj%origin%dir(:) = proj%dir(:)
      proj%origin%internal = ion%liinit

      proj%ffpmin = ZERO

      proj%posold(:) = proj%pos(:)
      proj%emin = proj%e

      proj%virtual = .false.                 ! projectile real if generated here
      proj%virtualparent = .false.           ! no parent
      proj%state = STATE%MOVING              ! state is moving
      proj%ededuct = ZERO

      proj_start = proj
      proj_valid = found

   end subroutine getProjectile
!
!=======================================================================
!
   subroutine getXstart( pos, dir )
!             =========   o    o
!
!     Get the starting point and direction most recently calculated
!
!  Arguments:
!
      real(REAL8), intent(out) :: pos(3)        ! starting point of trajectory
      real(REAL8), intent(out) :: dir(3)        ! initial direction of motion
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      pos(:) = proj_start%pos(:)
      dir(:) = proj_start%dir(:)

   end subroutine getXStart
!
!=======================================================================
!
end module Ion_Proj_mod

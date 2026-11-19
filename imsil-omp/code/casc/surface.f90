module Casc_Surface_mod
!
!  Treats trajectories that cross the surface (apply the surface binding energy)
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, ONE, SMALL, EPSMA, BIG
   use Types_mod,          only: Proj_type
   use Vector_mod,         only: operator(.dot.), abs

   use Param_Setup_mod,    only: nelem, lsput
   use Param_Sn_mod,       only: sn
   use Param_Dam_mod,      only: dam

   use Target_mod,         only: getDomainIndex
   use Target_Surface_mod, only: getDistSurface, getEntrancePoint, getExitPoint
   use Target_mod,         only: getDomainProperty
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
   public :: enterTarget                     ! get the entrance point into the target
!            ===========
   public :: exitTarget                      ! apply the surface binding energy
!            ==========
!
!  Private Procedures:
!  ==================
!
   private :: refractEntrance                ! change the direction of motion due to 
!             ===============                !    refraction/reflection at the surface
!                                            !    potential on entrance into target
   private :: refractExit                    ! change the direction of motion due to 
!             ===========                    !    refraction/reflection at the surface
!                                            !    potential on exit from target
!
contains
!
!===============================================================================
!
   subroutine enterTarget( proj, found )
!             ===========   io     o
!
!     Determine entrance point into the target and perform refraction 
!
!     NOTE: it is assumed that the projectile is inside the grid covering the target
!
!  Arguments:
!
      type(Proj_type), intent(inout) :: proj    ! projectile
      logical, intent(out) :: found             ! flag indicating whether an entrance point
                                                !    has been found
!
!  Local Variables:
!
      logical :: outside
      real(REAL8) :: dist, ffp, esurf(nelem), dirsurf(3), pos(3), posold(3)
!
!  Start of Executable Code:
!  ------------------------
!
      pos(:) = proj%pos(:)                      ! avoid changing proj%pos until entrance is
                                                !    guarranteed
      dist = proj%distsurf

      if( dist > - sn%pmax*(ONE-EPSMA) ) then
         error stop 'enterTarget: entered closer than pmax'
      end if
!
!  Step forward until distance from surface is smaller than pmax
!
      do
         posold(:) = pos(:)
         ffp = min( sn%ffpmax, abs( dist ) )
         pos(:) = pos(:) + ffp * proj%dir(:)

         call getDistSurface( pos(:), dist, outside )
!             ==============    i      o       o
         if( outside ) then
            found = .false.                     ! no entrance point found
            return
         else if( dist > - sn%pmax ) then
            exit
         end if
      end do
!
!  Now entrance is guarranteed and dist=-sn%pmax is between pos and posold
!
      proj%pos(:) = pos(:)                      ! remember pos and posold
      proj%posold(:) = posold(:)

      call getEntrancePoint( posold(:), proj%pos(:), found )
!          ================    io          io          o

      if( .not.found ) then
         print *, 'enterTarget: Could not determine entrance point'
         return
      end if
!
!  Update projectile region and distance from surface
!
      proj%id = 1                                    ! guess
      call getDomainIndex( proj%pos(:), proj%id )
!          ==============      i          io
      call getDistSurface( proj%pos(:), proj%distsurf, dir_=dirsurf(:) )
!          ==============      i              o               o
!
!  Consider the surface binding energy
!
      if( lsput ) then
         call getDomainProperty( proj%id, esurf_=esurf(:) )
!             =================     i          o
         call refractEntrance( proj, dirsurf(:), esurf(proj%ielem) )
!             ===============   io      i              i
      end if

   end subroutine enterTarget
!
!=======================================================================
!
   subroutine exitTarget( proj, leaving, error )
!             ==========   io      o       o
!
!     Determine the exit point and perform refraction / reflection by applying a planar
!     surface potential
!
!     It is assumed that the present position is outside -sn%pmax and the old position
!     was inside. If this condition is not fulfilled, the error flag is set.
!
!  Arguments:
!
      type(Proj_type), intent(inout) :: proj    ! projectile
      logical, intent(out) :: leaving           ! flag indicating projectile is leaving
                                                !    the target (=.false. if error=.true.)
      logical, intent(out) :: error             ! flag indicating that both present and last
                                                !    position are outside target
!
!  Local variables:
!
      logical :: found, reflect
      real(REAL8) :: esurf(nelem), dirsurf(3), posold(3)
!
!  Start of executable code:
!  ------------------------
!
      error = .false.
      leaving = .true.
!
!  Put the ion back to the boundary of the interface layer. If not found continue simulation
!  outside the target until grid boundary is reached.
!
      posold(:) = proj%posold(:)                ! avoid overwriting of proj%posold

      call getExitPoint( posold(:), proj%pos(:), found )
!          ============    io           io         o

      if( .not.found ) then
         error = .true.
         leaving = .false.
         return
      end if
!
!  Update projectile region and distance from surface
!
      proj%id = 1                                    !
      call getDomainIndex( proj%pos(:), proj%id )    ! inserted in analogy to enterTarget
!          ==============      i          io         !
      call getDistSurface( proj%pos(:), proj%distsurf, dir_=dirsurf(:) )
!          ==============      i              o              o
!
!  Consider the surface binding energy
!
      if( lsput ) then
         call getDomainProperty( proj%id, esurf_=esurf(:) )
!             =================     i          o
         call refractExit( proj, dirsurf(:), esurf(proj%ielem), reflect )
!             ===========   io      i              i               o

         if( reflect ) then                     ! if reflected, projectile must be inside
            proj%pos(:) = posold(:)
            call getDomainIndex( proj%pos(:), proj%id )
!                ==============      i           io
            call getDistSurface( proj%pos(:), proj%distsurf )
!                ==============      i              o
         end if
      else
         reflect = .false.
      end if

      leaving = .not.reflect                    ! projectile is leaving if it is not reflected

   end subroutine exitTarget
!
!===============================================================================
!
   subroutine refractEntrance( proj, dirsurf, esurf )
!             ===============   io      i       i
!
!     Change the direction of motion due to refraction at the surface potential
!
!     NOTE: The potential is always higher on that side of the surface where the projectile
!           comes from, i.e. the projectile is always assumed to move downhill.
!           THEREFORE THIS PROCEDURE MUST ONLY BE USED FOR PROJECTILES ENTERING THE TARGET.
!
!  Arguments:
!
      type(Proj_type), intent(inout) :: proj    ! projectile
      real(REAL8), intent(in) :: dirsurf(3)     ! surface normal
      real(REAL8), intent(in) :: esurf          ! surface binding energy [eV]
!
!  Local Variables:
!
      real(REAL8) :: cosalfa, cosbeta, sinbeta, dirpar(3), dirperp(3)
      real(REAL8) :: eperp, dirlen
!
!  Start of Executable Code:
!  ------------------------
!
!  Isotropic surface potential
!
      if( abs(dam%ksurf) > EPSMA ) then
         proj%e = proj%e + esurf
         return
      end if
!
!  Planar surface potential
!
      cosalfa = proj%dir(:) .dot. dirsurf(:)    ! alfa = angle between incidence direction
!                                               !        and surface normal
      dirperp(:) = dirsurf(:) * cosalfa         ! dirperp = component of incidence direction
!                                               !           perpendicular to surface (inwards
                                                !           target)
      dirpar(:) = proj%dir(:) - dirperp(:)      ! dirpar = component of incidence direction
                                                !          parallel to the surface

      eperp = proj%e * cosalfa**2
!
!  Refraction
!
      dirlen = abs(dirpar(:))

      if( dirlen > SMALL ) then                 ! beta = angle between exit direction and
                                                !        negative surface normal
         cosbeta = sqrt( (eperp+esurf)/(proj%e+esurf) )
         sinbeta = sqrt( (proj%e-eperp)/(proj%e+esurf) )

         dirpar(:) = dirpar(:) / dirlen
         proj%dir(:) = - dirsurf(:) * cosbeta + dirpar(:) * sinbeta

      else
         ! do nothing
      end if

      proj%e = proj%e + esurf

   end subroutine refractEntrance
!
!=======================================================================
!
   subroutine refractExit( proj, dirsurf, esurf, reflect )
!             ===========   io      i       i       o
!
!     Change the direction of motion due to refraction/reflection at the surface potential
!
!     NOTE: The potential is always higher on that side of the surface where the projectile
!           moves towards, i.e. the projectile is always assumed to move uphill.
!           THEREFORE THIS PROCEDURE MUST ONLY BE USED FOR PROJECTILES LEAVING THE TARGET.
!
!  Arguments:
!
      type(Proj_type), intent(inout) :: proj    ! projectile
      real(REAL8), intent(in) :: dirsurf(3)     ! surface normal
      real(REAL8), intent(in) :: esurf          ! surface binding energy [eV]
      logical, intent(out) :: reflect           ! flag indicating reflection
!
!  Local Variables:
!
      real(REAL8) :: cosalfa, cosbeta, sinbeta, dirpar(3), dirperp(3)
      real(REAL8) :: dirsurfp(3)
      real(REAL8) :: eperp, dirlen
!
!  Start of Executable Code:
!  ------------------------
!
!  Isotropic surface potential
!
      if( abs( dam%ksurf - ONE ) < EPSMA ) then
         reflect = (proj%e <= esurf)
         if( reflect ) then
            proj%dir(:) = - proj%dir(:)
            proj%nreflect = proj%nreflect + 1
         else
            proj%e = max( proj%e - esurf, 0.01_REAL8 )
         end if
         return
      end if
!
!  All other cases
!
!  Calculate surface potential normal
!
      if( abs(dam%ksurf) < EPSMA ) then         ! planar surface potential

         dirsurfp(:) = dirsurf(:)

      else                                      ! ellipsoidal/hyperboloidal surface potential

         cosalfa = proj%dir(:) .dot. dirsurf(:) ! alfa = angle between incidence direction
                                                !        and surface normal
         dirperp(:) = dirsurf(:) * cosalfa      ! dirperp = component of incidence direction
                                                !           perpendicular to surface
         dirpar(:) = proj%dir(:) - dirperp(:)   ! dirpar = component of incidence direction
                                                !          parallel to the surface
         dirsurfp(:) = dirperp(:) + dam%ksurf * dirpar(:)
         dirlen = abs(dirsurfp(:))
         if( dirlen > SMALL )  dirsurfp(:) = dirsurfp(:) / dirlen

      end if
!
!  Components perpendicular and normal to the surface potential
!
      cosalfa = proj%dir(:) .dot. dirsurfp(:)   ! alfa = angle between incidence direction
!                                               !        and surface potential normal

      if( cosalfa <= SMALL ) then               ! for alfa>90 deg, fall back to isotropic model
         reflect = (proj%e <= esurf)            ! (can happen in hyperbolic case, i.e. ksurf<0)
         if( reflect ) then                     ! dirlen<=SMALL is included in this case
            proj%dir(:) = - proj%dir(:)
            proj%nreflect = proj%nreflect + 1
         else
            proj%e = proj%e - esurf
         end if
         return
      end if

      dirperp(:) = dirsurfp(:) * cosalfa        ! dirperp = component of incidence direction
!                                               !           perpendicular to surface potential
      dirpar(:) = proj%dir(:) - dirperp(:)      ! dirpar = component of incidence direction
                                                !          parallel to the surface potential
!
!  Reflection/refracion from the surface potential, detemining the new projectile direction
!
      eperp = proj%e * cosalfa**2

      if( eperp <= esurf ) then
!
!  Reflection
!
         proj%dir(:) = - dirperp(:) + dirpar(:)

      else
!
!  Refraction
!
         dirlen = abs(dirpar(:))
!
         if( dirlen > SMALL ) then              ! beta = angle between exit direction and
                                                !        surface potential normal
            cosbeta = sqrt( (eperp-esurf)/(proj%e-esurf) )
            sinbeta = sqrt( (proj%e-eperp)/(proj%e-esurf) )

            dirpar(:) = dirpar(:) / dirlen
            proj%dir(:) = dirsurfp(:) * cosbeta + dirpar(:) * sinbeta

         else
            ! do nothing
         end if

      end if
!
!  Actual reflection/refraction relative to the surface
!
      reflect = ( (proj%dir(:) .dot. dirsurf(:)) <= ZERO )

      if( reflect ) then
         proj%nreflect = proj%nreflect + 1
      else
         proj%e = max( proj%e - esurf, 0.01_REAL8 )
      end if

   end subroutine refractExit
!
!=======================================================================
!
end module Casc_Surface_mod

module Target_Mechanics2d_mod
!
!  Provides functions that relax the target in response to the changes in cell contents
!
!  ATTENTION: This module is replaced by a stub on branch master. When merging master (or
!             any branch derived from master) into the dyn2d version, the dyn2d version of
!             the module must be kept.
!
!  Note: Cartesian = vertices have integer coordinates
!        regular = elements are congruent rectangles/rectilinearparallelepipeds
!        rectilinear = same as regular except that grid spacing may vary
!        curvilinear/structured = same combinatorial structure as regular grid
!        (according to Wikipedia)
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, ONE, TWO, THREE, FOUR, HALF, QUART, BIG, EPSMB

   use Array_mod,       only: check_allocation
   use Open_mod,        only: openFile

   use Param_Setup_mod, only: nionout
   use Param_Geom_mod,  only: geom
   use Param_Dam_mod,   only: dam

   !use Ion_Count_mod,   only: lastProjectile
!
#ifdef __INTEL_COMPILER
   use f95_precision,   only: wp => dp     !uncomment gesv!
   use lapack95,        only: gesv
#else
   use Const_mod,       only: wp => REAL8
#endif
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
   !integer(INT4), parameter :: wp = REAL8
!
!  Public Procedures:
!  =================
!
   public :: relaxAreas2d                       ! relax 2D grid such that cells have desired
!            ============                       !    area (version avoiding small cells)
   public :: relaxAreas2d_simple                ! relax 2D grid such that cells have desired
!            ===================                !    area
   public :: relaxViscoElastic2d                ! relax 2D grid according to a viscoelastic
!            ===================                !    model
!
!  Private Procedures:
!  ==================
!
   private :: dumpGrid                           ! dump grid to file
!             ========
!
contains
!
!=======================================================================
!
   subroutine relaxAreas2d( x, z, area_rlx, mask, movable )
!             ============  io io    io      i       i
!
!     Relax 2d grid such that cell areas have desired area area_rlx
!
!     Grid line indices run from 0:n
!     Cell indices run from 1:n
!     Cell i is between grid lines i-1 and i
!
!  Arguments:
!
      real(REAL8), intent(inout) :: x(0:,0:)       ! x-coordinates of vertices
      real(REAL8), intent(inout) :: z(0:,0:)       ! z-coordinates of vertices
      real(REAL8), intent(inout) :: area_rlx(:,:)  ! on entry: desired relaxed cell areas
                                                   ! on exit: obtained relaxed cell areas
      logical, intent(in) :: mask(:,:)             ! indicates cells to be considered
      logical, intent(in) :: movable(0:,0:)        ! indicates points to be movable
!
!  Local Variables:
!
      integer(INT4) :: nx, nz, nx0, nz0, nx1, nz1
      real(REAL8) :: maxerror, maxerror_old, u, denom, fac, maxarea
      real(REAL8) :: area(size( area_rlx, 1 ),size( area_rlx, 2 ))
      real(REAL8) :: f(size( area_rlx, 1 ),size( area_rlx, 2 ))
      real(REAL8) :: error(size( area_rlx, 1 ),size( area_rlx, 2 ))
      real(REAL8) :: xold(0:size( area_rlx, 1 ),0:size( area_rlx, 2 ))
      real(REAL8) :: zold(0:size( area_rlx, 1 ),0:size( area_rlx, 2 ))
      real(REAL8) :: dudx(0:size( area_rlx, 1 ),0:size( area_rlx, 2 ))
      real(REAL8) :: dudz(0:size(area_rlx,1),0:size(area_rlx,2))
!
      real(REAL8), parameter :: ERROR_MAX = 0.01    ! stop iteration when error < ERROR_MAX *
                                                    !    maxarea everywhere
      real(REAL8), parameter :: ERROR_CHANGE = 0.0001 ! stop iteration when relative change in
                                                    !    error < ERROR_CHANGE
      real(REAL8), parameter :: AREA_MIN = 0.1      ! unmasked cells should not become smaller
                                                    !    than area_min, masked cells head
                                                    !    towards (area_rlx+area_min)/2 if
                                                    !    area_rlx < area_min
                                                    !    (area_min = AREA_MIN*maxarea)
!
!  Start of Executable Code:
!  ------------------------
!
      maxerror_old = ZERO
!
      nx = size( area_rlx, 1 )
      nz = size( area_rlx, 2 )
!
      nx0 = nx - 1
      nz0 = nz - 1
      nx1 = nx
      nz1 = nz
!
      area(:,:) = HALF * (  &
         ( x(1:nx1,0:nz0) - x(0:nx0,0:nz0) ) * ( z(0:nx0,1:nz1) - z(0:nx0,0:nz0) ) -  &
         ( x(0:nx0,1:nz1) - x(0:nx0,0:nz0) ) * ( z(1:nx1,0:nz0) - z(0:nx0,0:nz0) ) +  &
         ( x(1:nx1,1:nz1) - x(0:nx0,1:nz1) ) * ( z(1:nx1,1:nz1) - z(1:nx1,0:nz0) ) -  &
         ( x(1:nx1,1:nz1) - x(1:nx1,0:nz0) ) * ( z(1:nx1,1:nz1) - z(0:nx0,1:nz1) ) )
      maxarea = maxval( area(:,:) )
      print *, 'maxarea = ', maxarea
!
      where( .not.mask(:,:) )  area_rlx(:,:) = max( area_rlx(:,:), AREA_MIN*maxarea )
!
      do
!
         call dumpGrid( x, z )
!             ========  i  i
!
!  Check error in area
!
         where( mask(:,:) .and. area_rlx(:,:) > AREA_MIN*maxarea )
            error(:,:) = abs( area(:,:) - area_rlx(:,:) )
         elsewhere( mask(:,:) )
            error(:,:) = abs( area(:,:) - HALF*(area_rlx(:,:)+AREA_MIN*maxarea) )
         elsewhere( area > area_rlx(:,:) )
            error(:,:) = ZERO
         elsewhere
            error(:,:) = abs( area(:,:) - area_rlx(:,:) )
         end where

         maxerror = maxval( error(:,:) )
         write( *, '(a,1p,e10.2,a,2i4)' )  &
            'relaxAreas2d: Max.Error=', maxerror, ' at ix,iz=', maxloc( error(:,:) )
!
         if( maxerror < ERROR_MAX * maxarea .or.  &
             abs( maxerror-maxerror_old ) < ERROR_CHANGE*max( maxerror, maxerror_old ) )  exit
         maxerror_old = maxerror
!
!  Calculate energy function
!
         u = sum( (area(:,:)-area_rlx(:,:))**2,  &
                  mask=(mask(:,:).or.(area(:,:)<area_rlx(:,:))) ) +  &
             sum( (area(:,:)-AREA_MIN*maxarea)**2,  &
                  mask=(mask(:,:).and.(area(:,:)<AREA_MIN*maxarea)) )
!
!  Calculate half of the derivative of energy wrt area
!
         f(:,:) = ZERO
         where( mask(:,:) .or. area(:,:) < area_rlx(:,:) )  f(:,:) = area(:,:) - area_rlx(:,:)
         where( mask(:,:) .and. area(:,:) < AREA_MIN*maxarea )  &
            f(:,:) = f(:,:) + area(:,:) - AREA_MIN*maxarea
!
!  Calculate derivative of energy wrt coordinates
!
         dudx(:,:) = ZERO
         dudz(:,:) = ZERO
!
         dudx(0:nx0,0:nz0) = dudx(0:nx0,0:nz0) + f(:,:) * ( z(1:nx1,0:nz0) - z(0:nx0,1:nz1) )
         dudx(0:nx0,1:nz1) = dudx(0:nx0,1:nz1) + f(:,:) * ( z(0:nx0,0:nz0) - z(1:nx1,1:nz1) )
         dudx(1:nx1,0:nz0) = dudx(1:nx1,0:nz0) + f(:,:) * ( z(1:nx1,1:nz1) - z(0:nx0,0:nz0) )
         dudx(1:nx1,1:nz1) = dudx(1:nx1,1:nz1) + f(:,:) * ( z(0:nx0,1:nz1) - z(1:nx1,0:nz0) )
!
         dudz(0:nx0,0:nz0) = dudz(0:nx0,0:nz0) + f(:,:) * ( x(0:nx0,1:nz1) - x(1:nx1,0:nz0) )
         dudz(0:nx0,1:nz1) = dudz(0:nx0,1:nz1) + f(:,:) * ( x(1:nx1,1:nz1) - x(0:nx0,0:nz0) )
         dudz(1:nx1,0:nz0) = dudz(1:nx1,0:nz0) + f(:,:) * ( x(0:nx0,0:nz0) - x(1:nx1,1:nz1) )
         dudz(1:nx1,1:nz1) = dudz(1:nx1,1:nz1) + f(:,:) * ( x(1:nx1,0:nz0) - x(0:nx0,1:nz1) )
!
         denom = sum( dudx(:,:)**2 + dudz(:,:)**2 )
!
!  Calculate relaxed positions, but reduce increment if any area becomes negative
!
         fac = ONE
         xold(:,:) = x(:,:)
         zold(:,:) = z(:,:)
!
         do
            where( movable(:,:) )
               x(:,:) = xold(:,:) - fac * dudx(:,:) * u / denom
               z(:,:) = zold(:,:) - fac * dudz(:,:) * u / denom
            end where
!
           area(:,:) = HALF * (  &
              ( x(1:nx1,0:nz0) - x(0:nx0,0:nz0) ) * ( z(0:nx0,1:nz1) - z(0:nx0,0:nz0) ) -  &
              ( x(0:nx0,1:nz1) - x(0:nx0,0:nz0) ) * ( z(1:nx1,0:nz0) - z(0:nx0,0:nz0) ) +  &
              ( x(1:nx1,1:nz1) - x(0:nx0,1:nz1) ) * ( z(1:nx1,1:nz1) - z(1:nx1,0:nz0) ) -  &
              ( x(1:nx1,1:nz1) - x(1:nx1,0:nz0) ) * ( z(1:nx1,1:nz1) - z(0:nx0,1:nz1) ) )
!
            if( any( area(:,:) < ZERO .and. mask(:,:) ) ) then
               write( *, * )  'relaxAreas2d: negative area, reducing increment...'
               fac = HALF * fac
            else
               exit
            end if
         end do
!
      end do
!
      call dumpGrid( x, z )
!          ========  i  i
!
      area_rlx(:,:) = area(:,:)
!
   end subroutine relaxAreas2d
!
!=======================================================================
!
   subroutine relaxAreas2d_simple( x, z, area_rlx, mask, movable )
!             ===================  io io    io      i       i
!
!     Relax 2d grid such that cell areas have desired area area_rlx
!
!     Grid line indices run from 0:n
!     Cell indices run from 1:n
!     Cell i is between grid lines i-1 and i
!
!  Arguments:
!
      real(REAL8), intent(inout) :: x(0:,0:)       ! x-coordinates of vertices
      real(REAL8), intent(inout) :: z(0:,0:)       ! z-coordinates of vertices
      real(REAL8), intent(inout) :: area_rlx(:,:)  ! on entry: desired relaxed cell areas
                                                   ! on exit: obtained relaxed cell areas
      logical, intent(in) :: mask(:,:)             ! indicates cells to be considered
      logical, intent(in) :: movable(0:,0:)        ! indicates points to be movable
!
!  Local Variables:
!
      integer(INT4) :: nx, nz, nx0, nz0, nx1, nz1
      real(REAL8) :: rmserror, rmserror_old, u, denom, fac, maxarea
      real(REAL8) :: area(size( area_rlx, 1 ),size( area_rlx, 2 ))
      real(REAL8) :: f(size( area_rlx, 1 ),size( area_rlx, 2 ))
      real(REAL8) :: error(size( area_rlx, 1 ),size( area_rlx, 2 ))
      real(REAL8) :: xold(0:size( area_rlx, 1 ),0:size( area_rlx, 2 ))
      real(REAL8) :: zold(0:size( area_rlx, 1 ),0:size( area_rlx, 2 ))
      real(REAL8) :: dudx(0:size( area_rlx, 1 ),0:size( area_rlx, 2 ))
      real(REAL8) :: dudz(0:size( area_rlx, 1 ),0:size( area_rlx, 2 ))
!
      real(REAL8), parameter :: ERROR_MAX = 0.0001   ! stop iteration when error < ERROR_MAX *
                                                    !    maxarea everywhere
      real(REAL8), parameter :: ERROR_CHANGE = 0.000001 ! stop iteration when relative change in
                                                    !    error < ERROR_CHANGE
      real(REAL8), parameter :: AREA_MIN = 0.1      ! unmasked cells should not become smaller
                                                    !    than area_min, masked cells head
                                                    !    towards (area_rlx+area_min)/2 if
                                                    !    area_rlx < area_min
                                                    !    (area_min = AREA_MIN*maxarea)
!
!  Start of Executable Code:
!  ------------------------
!
      rmserror_old = BIG
!
      nx = size( area_rlx, 1 )
      nz = size( area_rlx, 2 )
!
      nx0 = nx - 1
      nz0 = nz - 1
      nx1 = nx
      nz1 = nz
!
      area(:,:) = HALF * (  &
         ( x(1:nx1,0:nz0) - x(0:nx0,0:nz0) ) * ( z(0:nx0,1:nz1) - z(0:nx0,0:nz0) ) -  &
         ( x(0:nx0,1:nz1) - x(0:nx0,0:nz0) ) * ( z(1:nx1,0:nz0) - z(0:nx0,0:nz0) ) +  &
         ( x(1:nx1,1:nz1) - x(0:nx0,1:nz1) ) * ( z(1:nx1,1:nz1) - z(1:nx1,0:nz0) ) -  &
         ( x(1:nx1,1:nz1) - x(1:nx1,0:nz0) ) * ( z(1:nx1,1:nz1) - z(0:nx0,1:nz1) ) )
      maxarea = maxval( area(:,:) )
      print *, 'maxarea = ', maxarea
!
      error(:,:) = ZERO
      !where( .not.mask(:,:) )  area_rlx(:,:) = max( area_rlx(:,:), AREA_MIN*maxarea )
!
!  Check error in area
!
      where( mask(:,:) )
         error(:,:) = abs( area(:,:) - area_rlx(:,:) )
      end where
!
      rmserror = sqrt( sum( error(:,:), mask=mask(:,:) ) / count( mask(:,:) ) )
      write( *, '(a,1p,e10.2,a,2i4)' )  &
         'relaxAreas2d: RMS Error=', rmserror, ', maximum at ix,iz=', maxloc( error(:,:) )
!
!  Interate
!
      do
!
         call dumpGrid( x, z )
!             ========  i  i
!
!  Calculate energy function
!
         u = sum( (area(:,:)-area_rlx(:,:))**2,  &
                  mask=mask(:,:) )
         !u = sum( (area(:,:)-area_rlx(:,:))**2,  &
         !         mask=(mask(:,:).or.(area(:,:)<area_rlx(:,:))) ) +  &
         !    sum( (area(:,:)-AREA_MIN*maxarea)**2,  &
         !         mask=(mask(:,:).and.(area(:,:)<AREA_MIN*maxarea)) )
!
!  Calculate half of the derivative of energy wrt area
!
         f(:,:) = ZERO
         where( mask(:,:) )  f(:,:) = area(:,:) - area_rlx(:,:)
         !where( mask(:,:) .or. area(:,:) < area_rlx(:,:) )  f(:,:) = area(:,:) - area_rlx(:,:)
         !where( mask(:,:) .and. area(:,:) < AREA_MIN*maxarea )  &
         !   f(:,:) = f(:,:) + area(:,:) - AREA_MIN*maxarea
!
!  Calculate derivative of energy wrt coordinates
!
         dudx(:,:) = ZERO
         dudz(:,:) = ZERO
!
         dudx(0:nx0,0:nz0) = dudx(0:nx0,0:nz0) + f(:,:) * ( z(1:nx1,0:nz0) - z(0:nx0,1:nz1) )
         dudx(0:nx0,1:nz1) = dudx(0:nx0,1:nz1) + f(:,:) * ( z(0:nx0,0:nz0) - z(1:nx1,1:nz1) )
         dudx(1:nx1,0:nz0) = dudx(1:nx1,0:nz0) + f(:,:) * ( z(1:nx1,1:nz1) - z(0:nx0,0:nz0) )
         dudx(1:nx1,1:nz1) = dudx(1:nx1,1:nz1) + f(:,:) * ( z(0:nx0,1:nz1) - z(1:nx1,0:nz0) )
!
         dudz(0:nx0,0:nz0) = dudz(0:nx0,0:nz0) + f(:,:) * ( x(0:nx0,1:nz1) - x(1:nx1,0:nz0) )
         dudz(0:nx0,1:nz1) = dudz(0:nx0,1:nz1) + f(:,:) * ( x(1:nx1,1:nz1) - x(0:nx0,0:nz0) )
         dudz(1:nx1,0:nz0) = dudz(1:nx1,0:nz0) + f(:,:) * ( x(0:nx0,0:nz0) - x(1:nx1,1:nz1) )
         dudz(1:nx1,1:nz1) = dudz(1:nx1,1:nz1) + f(:,:) * ( x(1:nx1,0:nz0) - x(0:nx0,1:nz1) )
!
         denom = sum( dudx(:,:)**2 + dudz(:,:)**2 )
!
!  Calculate relaxed positions, but reduce increment if any area becomes negative
!
         fac = ONE
         xold(:,:) = x(:,:)
         zold(:,:) = z(:,:)
!
         do
            where( movable(:,:) )
               x(:,:) = xold(:,:) - fac * dudx(:,:) * u / denom
               z(:,:) = zold(:,:) - fac * dudz(:,:) * u / denom
            end where
!
           area(:,:) = HALF * (  &
              ( x(1:nx1,0:nz0) - x(0:nx0,0:nz0) ) * ( z(0:nx0,1:nz1) - z(0:nx0,0:nz0) ) -  &
              ( x(0:nx0,1:nz1) - x(0:nx0,0:nz0) ) * ( z(1:nx1,0:nz0) - z(0:nx0,0:nz0) ) +  &
              ( x(1:nx1,1:nz1) - x(0:nx0,1:nz1) ) * ( z(1:nx1,1:nz1) - z(1:nx1,0:nz0) ) -  &
              ( x(1:nx1,1:nz1) - x(1:nx1,0:nz0) ) * ( z(1:nx1,1:nz1) - z(0:nx0,1:nz1) ) )
!
            if( any( area(:,:) < ZERO .and. mask(:,:) ) ) then
               write( *, * )  'relaxAreas2d: negative area, reducing increment...'
               fac = HALF * fac
               cycle
            end if
!
            where( mask(:,:) )
               error(:,:) = abs( area(:,:) - area_rlx(:,:) )
            end where
!
            rmserror = sqrt( sum( error(:,:), mask=mask(:,:) ) / count( mask(:,:) ) )
            if( rmserror > rmserror_old ) then
               !write( *, * )  'relaxAreas2d: increasing error, reducing increment...'
               fac = HALF * fac
               cycle
            end if
!
            exit
!
         end do
!
         if( rmserror < ERROR_MAX * maxarea .or.  &
             abs( rmserror-rmserror_old ) < ERROR_CHANGE*max( rmserror, rmserror_old ) )  exit
!
         rmserror_old = rmserror
!
      end do
!
      write( *, '(a,1p,e10.2,a,2i4)' )  &
         'relaxAreas2d: RMS Error=', rmserror, ', maximum at ix,iz=', maxloc( error(:,:) )
!
      call dumpGrid( x, z )
!          ========  i  i
!
      area_rlx(:,:) = area(:,:)
!
   end subroutine relaxAreas2d_simple
!
!=======================================================================
!
   subroutine relaxViscoElastic2d( x, y, mask, movable, ff0, ff, dpa, stress )
!             ===================  io io  i       i      i   i    i     io
!
!     Relax 2d grid such that cell areas head towards the fill fraction ff using a
!     viscoelastic model
!
!     Grid line indices run from 0:n
!     Cell indices run from 1:n
!     Cell i is between grid lines i-1 and i
!
!     NOTE: This subroutine uses large automatic arrays. Use the -heap-arrays option of
!           the Intel Fortran compiler or say "ulimit -s unlimited" under Linux
!
!     NOTE: The suffices _l, _r, _b, _t denote the half-midlines of a cell:
!
!                           01|             |11
!                           --+------+------+--
!                             |      |t     |
!                             |  l   |   r  |
!                             +------+------+
!                             |      |      |
!                             |      |b     |
!                           --+------+------+--
!                           00|      |      |10
!
!  Arguments:
!
      real(REAL8), intent(inout) :: x(0:,0:)          ! (0:nx,0:ny) x-coordinates of nodes
      real(REAL8), intent(inout) :: y(0:,0:)          ! (0:nx,0:ny) y-coordinates of nodes
      logical, intent(in) :: mask(:,:)                ! (nx,ny) cells to be relaxed
      logical, intent(in) :: movable(0:,0:)           ! (0:nx,0:ny) movable nodes
      real(REAL8), intent(in) :: ff0(:,:)             ! (nx,ny) fill fractions at beginning of
                                                      !    relaxation
      real(REAL8), intent(in) :: ff(:,:)              ! (nx,ny) fill fractions at end of
                                                      !    relaxation, if there were no
                                                      !    volume change
      real(REAL8), intent(in) :: dpa(:,:)             ! (nx,ny) displacements per atom
      real(REAL8), intent(inout) :: stress(:,:,:)     ! (13,nx,ny) pressure plus deviatoric
                                                      !     stress components in 4 cell parts
                                                      !     (sxxdev, syydev, sxydev; l, r, b, t)
!
!  Local Variables:
!
      integer(INT4), parameter :: nt = 1      ! number of time steps
!
      integer(INT4) :: ix, iy, nx, ny, iu, nb, info, i, ixx, iyy, iuu, istat, it
      real(REAL8) :: dx_dxi, dx_deta, dy_dxi, dy_deta, d2x_dxideta, d2y_dxideta
      real(REAL8) :: detxy_c, detxy_d_xi, detxy_d_eta
      real(REAL8) :: ucoef_a, ucoef_b_xi, ucoef_b_eta
      real(REAL8) :: dux_dx_uxcoef(0:1,0:1),  &
                     dux_dx_uxcoef_l(0:1,0:1), dux_dx_uxcoef_r(0:1,0:1),  &
                     dux_dx_uxcoef_b(0:1,0:1), dux_dx_uxcoef_t(0:1,0:1)
      real(REAL8) :: duy_dx_uycoef(0:1,0:1),  &
                     duy_dx_uycoef_l(0:1,0:1), duy_dx_uycoef_r(0:1,0:1),  &
                     duy_dx_uycoef_b(0:1,0:1), duy_dx_uycoef_t(0:1,0:1)
      real(REAL8) :: dux_dy_uxcoef(0:1,0:1),  &
                     dux_dy_uxcoef_l(0:1,0:1), dux_dy_uxcoef_r(0:1,0:1),  &
                     dux_dy_uxcoef_b(0:1,0:1), dux_dy_uxcoef_t(0:1,0:1)
      real(REAL8) :: duy_dy_uycoef(0:1,0:1),  &
                     duy_dy_uycoef_l(0:1,0:1), duy_dy_uycoef_r(0:1,0:1),  &
                     duy_dy_uycoef_b(0:1,0:1), duy_dy_uycoef_t(0:1,0:1)
      real(REAL8) :: exx_uxcoef(0:1,0:1), exx_uxcoef_l(0:1,0:1), exx_uxcoef_r(0:1,0:1),  &
                                          exx_uxcoef_b(0:1,0:1), exx_uxcoef_t(0:1,0:1)
      real(REAL8) :: eyy_uycoef(0:1,0:1), eyy_uycoef_l(0:1,0:1), eyy_uycoef_r(0:1,0:1),  &
                                          eyy_uycoef_b(0:1,0:1), eyy_uycoef_t(0:1,0:1)
      real(REAL8) :: exy_uxcoef(0:1,0:1), exy_uxcoef_l(0:1,0:1), exy_uxcoef_r(0:1,0:1),  &
                                          exy_uxcoef_b(0:1,0:1), exy_uxcoef_t(0:1,0:1)
      real(REAL8) :: exy_uycoef(0:1,0:1), exy_uycoef_l(0:1,0:1), exy_uycoef_r(0:1,0:1),  &
                                          exy_uycoef_b(0:1,0:1), exy_uycoef_t(0:1,0:1)
      real(REAL8) :: dt_eta, geff, time, dtime, ffi
      real(REAL8) :: sxxdev_0_l, sxxdev_0_r, sxxdev_0_b, sxxdev_0_t
      real(REAL8) :: syydev_0_l, syydev_0_r, syydev_0_b, syydev_0_t
      real(REAL8) :: sxydev_0_l, sxydev_0_r, sxydev_0_b, sxydev_0_t
      real(REAL8) :: sxxdev_const_l, sxxdev_const_r, sxxdev_const_b, sxxdev_const_t
      real(REAL8) :: sxxdev_uxcoef_l(0:1,0:1), sxxdev_uxcoef_r(0:1,0:1),  &
                     sxxdev_uxcoef_b(0:1,0:1), sxxdev_uxcoef_t(0:1,0:1)
      real(REAL8) :: sxxdev_uycoef_l(0:1,0:1), sxxdev_uycoef_r(0:1,0:1),  &
                     sxxdev_uycoef_b(0:1,0:1), sxxdev_uycoef_t(0:1,0:1)
      real(REAL8) :: syydev_const_l, syydev_const_r, syydev_const_b, syydev_const_t
      real(REAL8) :: syydev_uxcoef_l(0:1,0:1), syydev_uxcoef_r(0:1,0:1),  &
                     syydev_uxcoef_b(0:1,0:1), syydev_uxcoef_t(0:1,0:1)
      real(REAL8) :: syydev_uycoef_l(0:1,0:1), syydev_uycoef_r(0:1,0:1),  &
                     syydev_uycoef_b(0:1,0:1), syydev_uycoef_t(0:1,0:1)
      real(REAL8) :: p_const, p_uxcoef(0:1,0:1), p_uycoef(0:1,0:1)
      real(REAL8) :: sxx_const_l, sxx_const_r, sxx_const_b, sxx_const_t
      real(REAL8) :: sxx_uxcoef_l(0:1,0:1), sxx_uxcoef_r(0:1,0:1),  &
                     sxx_uxcoef_b(0:1,0:1), sxx_uxcoef_t(0:1,0:1)
      real(REAL8) :: sxx_uycoef_l(0:1,0:1), sxx_uycoef_r(0:1,0:1),  &
                     sxx_uycoef_b(0:1,0:1), sxx_uycoef_t(0:1,0:1)
      real(REAL8) :: syy_const_l, syy_const_r, syy_const_b, syy_const_t
      real(REAL8) :: syy_uxcoef_l(0:1,0:1), syy_uxcoef_r(0:1,0:1),  &
                     syy_uxcoef_b(0:1,0:1), syy_uxcoef_t(0:1,0:1)
      real(REAL8) :: syy_uycoef_l(0:1,0:1), syy_uycoef_r(0:1,0:1),  &
                     syy_uycoef_b(0:1,0:1), syy_uycoef_t(0:1,0:1)
      real(REAL8) :: sxy_const_l, sxy_const_r, sxy_const_b, sxy_const_t
      real(REAL8) :: sxy_uxcoef_l(0:1,0:1), sxy_uxcoef_r(0:1,0:1),  &
                     sxy_uxcoef_b(0:1,0:1), sxy_uxcoef_t(0:1,0:1)
      real(REAL8) :: sxy_uycoef_l(0:1,0:1), sxy_uycoef_r(0:1,0:1),  &
                     sxy_uycoef_b(0:1,0:1), sxy_uycoef_t(0:1,0:1)
      real(REAL8) :: F_const_l(2), F_const_r(2), F_const_b(2), F_const_t(2)
      real(REAL8) :: F_coef_l(2,2,0:1,0:1), F_coef_r(2,2,0:1,0:1),  &
                     F_coef_b(2,2,0:1,0:1), F_coef_t(2,2,0:1,0:1)
!
      logical, allocatable :: maskbb(:,:,:)
      real(REAL8), allocatable :: AA(:,:,:,:,:,:)
      real(REAL8), allocatable :: bb(:,:,:)
      real(wp), allocatable :: A(:,:), b(:)
      real(REAL8), allocatable :: u(:,:,:)
      real(REAL8), allocatable :: vfac(:,:)
      real(REAL8) :: exx, eyy, p
      real(REAL8) :: sxxdev_l, sxxdev_r, sxxdev_b, sxxdev_t,  &
                     syydev_l, syydev_r, syydev_b, syydev_t,  &
                     sxydev_l, sxydev_r, sxydev_b, sxydev_t

      character(256) :: msg
!
!  Start of Executable Code:
!  ------------------------
!
      !print *, 'Entered relaxViscoElastic2d'
!
      nx = ubound(x,1)
      ny = ubound(x,2)
!
!  Allocate work arrays (NOTE: They are too big to be automatic arrays at time of writing)
!
      print *, nx, ny
      allocate( AA( 2, 0:nx, 0:ny, 2, 0:nx, 0:ny ), stat=istat, errmsg=msg )
         call check_allocation('relaxViscoElastic2d', 'AA', istat, msg )
      allocate( bb( 2, 0:nx, 0:ny ), stat=istat, errmsg=msg )
         call check_allocation('relaxViscoElastic2d', 'bb', istat, msg )
      allocate( maskbb( 2, 0:nx, 0:ny ), stat=istat, errmsg=msg )
         call check_allocation('relaxViscoElastic2d', 'maskbb', istat, msg )
      allocate( u( 2, 0:nx, 0:ny ), stat=istat, errmsg=msg )
         call check_allocation('relaxViscoElastic2d', 'u', istat, msg )
      allocate( vfac( nx, ny ),  stat=istat, errmsg=msg )
         call check_allocation('relaxViscoElastic2d', 'vfac', istat, msg )
!
      !print *, 'Successfully allocated work arrays'
!
      vfac(:,:) = ONE
!
      do it = 1, nt                             ! loop over time steps
!
         dtime = ONE / nt
         time = it * dtime                      ! Note: 0<time<=1
!
         maskbb(:,:,:) = .false.
         AA(:,:,:,:,:,:) = ZERO
         bb(:,:,:) = ZERO
!
         do iy = 1, ny                          ! loop over cells
            do ix = 1, nx
!
               if( .not. mask(ix,iy) )  cycle
!
!  Calculate cell geometry dependent parameters and strain and stress
!
               call calcStrainAndStress( )
!                   ===================
!
!  Calculate the forces:
!
!  first index: force components (1: Fx, 2: Fy)
!  second index: displacement components (1: ux, 2: uy)
!
!  (dx_deta,dy_deta) and (dx_dxi,dy_dxi) are the midline vectors
!
!  We omit a factor HALF since it cancels in the equations
!
               F_const_b(1) = sxx_const_b * dy_deta - sxy_const_b * dx_deta
               F_const_t(1) = sxx_const_t * dy_deta - sxy_const_t * dx_deta
               F_coef_b(1,1,:,:) = sxx_uxcoef_b(:,:) * dy_deta - sxy_uxcoef_b(:,:) * dx_deta
               F_coef_t(1,1,:,:) = sxx_uxcoef_t(:,:) * dy_deta - sxy_uxcoef_t(:,:) * dx_deta
               F_coef_b(1,2,:,:) = sxx_uycoef_b(:,:) * dy_deta - sxy_uycoef_b(:,:) * dx_deta
               F_coef_t(1,2,:,:) = sxx_uycoef_t(:,:) * dy_deta - sxy_uycoef_t(:,:) * dx_deta
!
               F_const_b(2) = sxy_const_b * dy_deta - syy_const_b * dx_deta
               F_const_t(2) = sxy_const_t * dy_deta - syy_const_t * dx_deta
               F_coef_b(2,1,:,:) = sxy_uxcoef_b(:,:) * dy_deta - syy_uxcoef_b(:,:) * dx_deta
               F_coef_t(2,1,:,:) = sxy_uxcoef_t(:,:) * dy_deta - syy_uxcoef_t(:,:) * dx_deta
               F_coef_b(2,2,:,:) = sxy_uycoef_b(:,:) * dy_deta - syy_uycoef_b(:,:) * dx_deta
               F_coef_t(2,2,:,:) = sxy_uycoef_t(:,:) * dy_deta - syy_uycoef_t(:,:) * dx_deta
!
               F_const_l(1) = - sxx_const_l * dy_dxi + sxy_const_l * dx_dxi
               F_const_r(1) = - sxx_const_r * dy_dxi + sxy_const_r * dx_dxi
               F_coef_l(1,1,:,:) = - sxx_uxcoef_l(:,:) * dy_dxi + sxy_uxcoef_l(:,:) * dx_dxi
               F_coef_r(1,1,:,:) = - sxx_uxcoef_r(:,:) * dy_dxi + sxy_uxcoef_r(:,:) * dx_dxi
               F_coef_l(1,2,:,:) = - sxx_uycoef_l(:,:) * dy_dxi + sxy_uycoef_l(:,:) * dx_dxi
               F_coef_r(1,2,:,:) = - sxx_uycoef_r(:,:) * dy_dxi + sxy_uycoef_r(:,:) * dx_dxi
!
               F_const_l(2) = - sxy_const_l * dy_dxi + syy_const_l * dx_dxi
               F_const_r(2) = - sxy_const_r * dy_dxi + syy_const_r * dx_dxi
               F_coef_l(2,1,:,:) = - sxy_uxcoef_l(:,:) * dy_dxi + syy_uxcoef_l(:,:) * dx_dxi
               F_coef_r(2,1,:,:) = - sxy_uxcoef_r(:,:) * dy_dxi + syy_uxcoef_r(:,:) * dx_dxi
               F_coef_l(2,2,:,:) = - sxy_uycoef_l(:,:) * dy_dxi + syy_uycoef_l(:,:) * dx_dxi
               F_coef_r(2,2,:,:) = - sxy_uycoef_r(:,:) * dy_dxi + syy_uycoef_r(:,:) * dx_dxi
!
!  Set up matrix and rhs for the system of linear equations
!
               AA(:,ix-1,iy-1,:,ix-1:ix,iy-1:iy) = AA(:,ix-1,iy-1,:,ix-1:ix,iy-1:iy)  &
                                                   + F_coef_b(:,:,:,:) + F_coef_l(:,:,:,:)
               AA(:,ix-1,iy  ,:,ix-1:ix,iy-1:iy) = AA(:,ix-1,iy  ,:,ix-1:ix,iy-1:iy) &
                                                   + F_coef_t(:,:,:,:) - F_coef_l(:,:,:,:)
               AA(:,ix  ,iy-1,:,ix-1:ix,iy-1:iy) = AA(:,ix  ,iy-1,:,ix-1:ix,iy-1:iy)  &
                                                   - F_coef_b(:,:,:,:) + F_coef_r(:,:,:,:)
               AA(:,ix  ,iy  ,:,ix-1:ix,iy-1:iy) = AA(:,ix  ,iy  ,:,ix-1:ix,iy-1:iy)  &
                                                   - F_coef_t(:,:,:,:) - F_coef_r(:,:,:,:)
!
!                                               ! opposite sign compared to A, since b is on rhs
               bb(:,ix-1,iy-1) = bb(:,ix-1,iy-1) - F_const_b(:) - F_const_l(:)
               bb(:,ix-1,iy  ) = bb(:,ix-1,iy  ) - F_const_t(:) + F_const_l(:)
               bb(:,ix  ,iy-1) = bb(:,ix  ,iy-1) + F_const_b(:) - F_const_r(:)
               bb(:,ix  ,iy  ) = bb(:,ix  ,iy  ) + F_const_t(:) + F_const_r(:)
!
!  Set nodes (ix-1:ix,iy-1:iy) movable
!
               maskbb(1,ix-1:ix,iy-1:iy) = .true.
!
            end do
         end do
!
!  Consider unmovable points
!
         maskbb(1,:,:) = maskbb(1,:,:) .and. movable(:,:)
         maskbb(2,:,:) = maskbb(1,:,:)
!
!  Eliminate equations for non-movable nodes
!
         nb = count( maskbb )
         print *, 'nb=', nb
!
         allocate( b( nb ), stat=istat )
         if( istat /= 0 )  stop 'relaxViscoElastic2d: Error during allocating b'
         !print *, 'relaxViscoElastic2d: b successfully allocated'
         allocate( A( nb, nb ), stat=istat )
         if( istat /= 0 )  stop 'relaxViscoElastic2d: Error during allocating A'
         !print *, 'relaxViscoElastic2d: A successfully allocated'
!
         i = 0                                  ! pack() doesn't seem to work with 6d arrays
         do iy = 0, ny
            do ix = 0, nx
               do iu = 1, 2
                  do iyy = 0, ny
                     do ixx = 0, nx
                        do iuu = 1, 2
                           if( maskbb(iu,ix,iy) .and. maskbb(iuu,ixx,iyy) ) then
                              i = i + 1
                              A( mod(i-1,nb)+1, (i-1)/nb+1 ) = AA(iuu,ixx,iyy,iu,ix,iy)
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do
         !print *, 'i=', i
         !print *, 'size(A)=', size( A )
         b = pack( bb, maskbb )
!
! debug
!
!         i = 0
!outerb:  do iy = 0, ny
!            do ix = 0, nx
!               do iu = 1, 2
!                  if( maskbb(iu,ix,iy) ) then
!                     i=i+1
!                     if( i > 10 )  exit outerb
!                     print '(a,3i4,2(a,1pe12.3))', 'iu,ix,iy=',  &
!                        iu,ix,iy, ' bb=', bb(iu,ix,iy),  &
!                        ' b=', b(i)
!                  end if
!               end do
!            end do
!         end do outerb
!
         deallocate( AA, bb )
!
!  Replace the last two equations by demanding the average displacement to be zero
!
!         A(nb-1,:) = ONE
!         A(nb,:) = ONE
!         b(nb) = ZERO
!
!  Solve the linear system of equations
!
         !print *, 'Reached LAPACK subroutine'
         !allocate( A_save, source=A )
         !allocate( b_save, source=b )
!
#ifdef __INTEL_COMPILER
         call gesv( A, b, info=info )
!             ====  io io       o
#else
         error stop 'relaxViscoElastic2d: MKL library procedure gesv not available'
#endif
!
         !do i = 1, size( b )
         !   print *, A_save(162,i), b(i), A_save(162,i)*b(i)
         !   !print *, sum(A_save(i,:)*b(:)), b_save(i)
         !end do
         !print  *, sum(A_save(162,:)*b(:)), b_save(162)
         !print *, info
!         stop
         !print *, 'Returned from LAPACK subroutine'
!
         if( info == 0 ) then
            ! nothing to do
         else if( info == -1 ) then
            print *, 'relaxViscoElastic2d: Matrix A had an illegal value'
         else if( info == -2 ) then
            print *, 'relaxViscoElastic2d: Vector b had an illegal value'
         else if( info > 0 ) then
            print *, 'relaxViscoElastic2d: LU decomposition produced a singular matrix at i=', &
                     info
         else
            print *, 'relaxViscoElastic2d: Error code: ', i
         end if
!
!  Update the nodes
!
         i = 0
         do iy = 0, ny
            do ix = 0, nx
               do iu = 1, 2
                  if( maskbb(iu,ix,iy) ) then
                     i = i + 1
                     u(iu,ix,iy) = b(i)
                  else
                     u(iu,ix,iy) = ZERO
                  end if
               end do
            end do
         end do
!
         deallocate( A, b )
!
!  Debug
!
         !write(12,*) 'ux:'
         !do iy=0,ny
         !   write( 12, '(20(1pe10.3))' )  u(1,:,iy)
         !end do
         !write(12,*) 'uy:'
         !do iy=0,ny
         !   write( 12, '(20(1pe10.3))' )  u(2,:,iy)
         !end do
         !!print *, 'nb=', nb, i
         !!u(:,:,:) = unpack( b, maskbb, field=ZERO )
         !!print *, 'lbound(u)=', lbound(u)
         !!print *, 'ubound(u)=', ubound(u)
         !print *, 'x:'
         !do iy=0,3
         !   print *, x(0:3,iy)
         !end do
         !print *, 'y:'
         !do iy=0,3
         !   print *, y(0:3,iy)
         !end do
         !print *, 'ux:'
         !do iy=0,3
         !   print *, u(1,0:3,iy)
         !end do
         !print *, 'uy:'
         !do iy=0,3
         !   print *, u(2,0:3,iy)
         !end do
!
         print *, 'Maximum ux: ', maxval(u(1,:,:)), ' at ix,iz=', maxloc(u(1,:,:))-(/1,1/)
         print *, 'Minimum ux: ', minval(u(1,:,:)), ' at ix,iz=', minloc(u(1,:,:))-(/1,1/)
         print *, 'Maximum uz: ', maxval(u(2,:,:)), ' at ix,iz=', maxloc(u(2,:,:))-(/1,1/)
         print *, 'Minimum uz: ', minval(u(2,:,:)), ' at ix,iz=', minloc(u(2,:,:))-(/1,1/)
!
!  Calculate the deviatoric stress
!
         do iy = 1, ny
            do ix = 1, nx
!
               if( .not. mask(ix,iy) )  cycle
!
               call calcStrainAndStress( )
!                   ===================
!
               exx = sum( exx_uxcoef(:,:) * u(1,ix-1:ix,iy-1:iy) )
               eyy = sum( eyy_uycoef(:,:) * u(2,ix-1:ix,iy-1:iy) )
!
               vfac(ix,iy) = vfac(ix,iy) * ( ONE + exx + eyy )
!
               p = dam%modbulk * ( ONE - vfac(ix,iy) / ff(ix,iy) )

               sxxdev_l = sxxdev_const_l + sum( sxxdev_uxcoef_l(:,:) * u(1,ix-1:ix,iy-1:iy) )  &
                                         + sum( sxxdev_uycoef_l(:,:) * u(2,ix-1:ix,iy-1:iy) )
               sxxdev_r = sxxdev_const_r + sum( sxxdev_uxcoef_r(:,:) * u(1,ix-1:ix,iy-1:iy) )  &
                                         + sum( sxxdev_uycoef_r(:,:) * u(2,ix-1:ix,iy-1:iy) )
               sxxdev_b = sxxdev_const_b + sum( sxxdev_uxcoef_b(:,:) * u(1,ix-1:ix,iy-1:iy) )  &
                                         + sum( sxxdev_uycoef_b(:,:) * u(2,ix-1:ix,iy-1:iy) )
               sxxdev_t = sxxdev_const_t + sum( sxxdev_uxcoef_t(:,:) * u(1,ix-1:ix,iy-1:iy) )  &
                                         + sum( sxxdev_uycoef_t(:,:) * u(2,ix-1:ix,iy-1:iy) )
!
               syydev_l = syydev_const_l + sum( syydev_uxcoef_l(:,:) * u(1,ix-1:ix,iy-1:iy) )  &
                                         + sum( syydev_uycoef_l(:,:) * u(2,ix-1:ix,iy-1:iy) )
               syydev_r = syydev_const_r + sum( syydev_uxcoef_r(:,:) * u(1,ix-1:ix,iy-1:iy) )  &
                                         + sum( syydev_uycoef_r(:,:) * u(2,ix-1:ix,iy-1:iy) )
               syydev_b = syydev_const_b + sum( syydev_uxcoef_b(:,:) * u(1,ix-1:ix,iy-1:iy) )  &
                                         + sum( syydev_uycoef_b(:,:) * u(2,ix-1:ix,iy-1:iy) )
               syydev_t = syydev_const_t + sum( syydev_uxcoef_t(:,:) * u(1,ix-1:ix,iy-1:iy) )  &
                                         + sum( syydev_uycoef_t(:,:) * u(2,ix-1:ix,iy-1:iy) )
!
               sxydev_l = sxy_const_l + sum( sxy_uxcoef_l(:,:) * u(1,ix-1:ix,iy-1:iy) ) +  &
                                        sum( sxy_uycoef_l(:,:) * u(2,ix-1:ix,iy-1:iy) )
               sxydev_r = sxy_const_r + sum( sxy_uxcoef_r(:,:) * u(1,ix-1:ix,iy-1:iy) ) +  &
                                        sum( sxy_uycoef_r(:,:) * u(2,ix-1:ix,iy-1:iy) )
               sxydev_b = sxy_const_b + sum( sxy_uxcoef_b(:,:) * u(1,ix-1:ix,iy-1:iy) ) +  &
                                        sum( sxy_uycoef_b(:,:) * u(2,ix-1:ix,iy-1:iy) )
               sxydev_t = sxy_const_t + sum( sxy_uxcoef_t(:,:) * u(1,ix-1:ix,iy-1:iy) ) +  &
                                        sum( sxy_uycoef_t(:,:) * u(2,ix-1:ix,iy-1:iy) )
!
               stress(:,ix,iy) = (/ p, sxxdev_l, syydev_l, sxydev_l,  &
                                       sxxdev_r, syydev_r, sxydev_r,  &
                                       sxxdev_b, syydev_b, sxydev_b,  &
                                       sxxdev_t, syydev_t, sxydev_t /)
!
!  Check consistency of pressure
!
               ! check still missing
!
            end do
         end do
!
!  Update the coordinates
!
         x(:,:) = x(:,:) + u(1,:,:)
         y(:,:) = y(:,:) + u(2,:,:)
!
         !print *, 'x:'
         !do iy=0,3
         !   print *, x(0:3,iy)
         !end do
         !print *, 'y:'
         !do iy=0,3
         !   print *, y(0:3,iy)
         !end do
!
      end do
!
      !print *, 'Maximum sxxdev: ', maxval(stressdev(1,:,:)), ' at ix,iz=', maxloc(stressdev(1,:,:))
      !print *, 'Minimum sxxdev: ', minval(stressdev(1,:,:)), ' at ix,iz=', minloc(stressdev(1,:,:))
      !print *, 'Maximum szzdev: ', maxval(stressdev(2,:,:)), ' at ix,iz=', maxloc(stressdev(2,:,:))
      !print *, 'Minimum szzdev: ', minval(stressdev(2,:,:)), ' at ix,iz=', minloc(stressdev(2,:,:))
      !print *, 'Maximum sxzdev: ', maxval(stressdev(3,:,:)), ' at ix,iz=', maxloc(stressdev(3,:,:))
      !print *, 'Minimum sxzdev: ', minval(stressdev(3,:,:)), ' at ix,iz=', minloc(stressdev(3,:,:))
!
      deallocate( maskbb, u, vfac, stat=istat, errmsg=msg )
      if( istat /= 0 )  goto 900
!
      return
!
900   print *, 'relaxViscoElastic2d: ' // trim( msg )
      stop
!
   contains
!
!-----------------------------------------------------------------------
!
      subroutine calcStrainAndStress( )
!                ===================
!
!  Cell geometry dependent parameters
!
100      dx_dxi  = HALF * ( x(ix,iy) + x(ix,iy-1) - x(ix-1,iy) - x(ix-1,iy-1) )
         dx_deta = HALF * ( x(ix,iy) - x(ix,iy-1) + x(ix-1,iy) - x(ix-1,iy-1) )
         dy_dxi  = HALF * ( y(ix,iy) + y(ix,iy-1) - y(ix-1,iy) - y(ix-1,iy-1) )
         dy_deta = HALF * ( y(ix,iy) - y(ix,iy-1) + y(ix-1,iy) - y(ix-1,iy-1) )
!
         detxy_c = dx_dxi * dy_deta - dx_deta * dy_dxi
!
         d2x_dxideta = x(ix,iy) - x(ix,iy-1) - x(ix-1,iy) + x(ix-1,iy-1)
         d2y_dxideta = y(ix,iy) - y(ix,iy-1) - y(ix-1,iy) + y(ix-1,iy-1)
!
         detxy_d_xi = dx_dxi * d2y_dxideta - d2x_dxideta * dy_dxi
         detxy_d_eta = d2x_dxideta * dy_deta - dx_deta * d2y_dxideta
!
!  The midlines must not cross the edges
!
         if( abs(detxy_d_xi) > TWO * abs(detxy_c) * ( ONE - EPSMB ) .or.  &
             abs(detxy_d_eta) > TWO * abs(detxy_c) * ( ONE - EPSMB ) ) then
            if( ( x(ix-1,iy) - x(ix-1,iy-1) ) * ( y(ix,iy-1) - y(ix-1,iy) ) - &
                ( y(ix-1,iy) - y(ix-1,iy-1) ) * ( x(ix,iy-1) - x(ix-1,iy) ) > 0 ) then
               print '(a,2(1pe15.6))', 'WARNING! calcStrainAndStress: moved point',  &
                                      x(ix-1,iy-1), y(ix-1,iy-1)
                  x(ix-1,iy-1) = HALF * ( x(ix-1,iy) + x(ix,iy-1) )
                  y(ix-1,iy-1) = HALF * ( y(ix-1,iy) + y(ix,iy-1) )
               print '(a,2(1pe15.6))', 'to', x(ix-1,iy-1), y(ix-1,iy-1)
               goto 100
            else if( ( x(ix-1,iy-1) - x(ix,iy-1) ) * ( y(ix,iy) - y(ix-1,iy-1) ) - &
                     ( y(ix-1,iy-1) - y(ix,iy-1) ) * ( x(ix,iy) - x(ix-1,iy-1) ) > 0 ) then
               print '(a,2(1pe15.6))', 'WARNING! calcStrainAndStress: moved point',  &
                                      x(ix,iy-1), y(ix,iy-1)
                  x(ix,iy-1) = HALF * ( x(ix-1,iy-1) + x(ix,iy) )
                  y(ix,iy-1) = HALF * ( y(ix-1,iy-1) + y(ix,iy) )
               print '(a,2(1pe15.6))', 'to', x(ix,iy-1), y(ix,iy-1)
               goto 100
            else if( ( x(ix,iy-1) - x(ix,iy) ) * ( y(ix-1,iy) - y(ix,iy-1) ) - &
                     ( y(ix,iy-1) - y(ix,iy) ) * ( x(ix-1,iy) - x(ix,iy-1) ) > 0 ) then
               print '(a,2(1pe15.6))', 'WARNING! calcStrainAndStress: moved point',  &
                                      x(ix,iy), y(ix,iy)
                  x(ix,iy) = HALF * ( x(ix,iy-1) + x(ix-1,iy) )
                  y(ix,iy) = HALF * ( y(ix,iy-1) + y(ix-1,iy) )
               print '(a,2(1pe15.6))', 'to', x(ix,iy), y(ix,iy)
               goto 100
            else if( ( x(ix,iy) - x(ix-1,iy) ) * ( y(ix-1,iy-1) - y(ix,iy) ) - &
                     ( y(ix,iy) - y(ix-1,iy) ) * ( x(ix-1,iy-1) - x(ix,iy) ) > 0 ) then
               print '(a,2(1pe15.6))', 'WARNING! calcStrainAndStress: moved point',  &
                                      x(ix-1,iy), y(ix-1,iy)
                  x(ix-1,iy) = HALF * ( x(ix,iy) + x(ix-1,iy-1) )
                  y(ix-1,iy) = HALF * ( y(ix,iy) + y(ix-1,iy-1) )
               print '(a,2(1pe15.6))', 'to', x(ix-1,iy), y(ix-1,iy)
               goto 100
            else                                ! ther should have been a reentrant corner
               stop 'calcStrainAndStress: internal error'
            end if
         end if
!
!  Spatial derivatives of displacements (ux,uy)
!
! XXX_uxcoef = coefficient of ux(ix,iz) in expression for XXX
! XXX_uycoef = coefficient of uy(ix,iz) in expression for XXX
!
!  Values calculated are for the total cell (no suffix) or
!  mean values over the half-midlines (_l, _r, _b, _t)
!
         !dux_dx_uxcoef(1,1) = HALF * ( dy_deta - dy_dxi ) / detxy_c
         ucoef_a = HALF * ( y(ix-1,iy) - y(ix,iy-1) )    ! D(u,y), (1,1)
         ucoef_b_xi = y(ix-1,iy-1) - y(ix,iy-1)
         ucoef_b_eta = y(ix-1,iy) - y(ix-1,iy-1)
         dux_dx_uxcoef(1,1) = ucoef_a / detxy_c
         dux_dx_uxcoef_l(1,1) = meanNeg( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dx_uxcoef_r(1,1) = meanPos( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dx_uxcoef_b(1,1) = meanNeg( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
         dux_dx_uxcoef_t(1,1) = meanPos( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
!
         !dux_dx_uxcoef(1,0) = HALF * ( dy_deta + dy_dxi ) / detxy_c
         ucoef_a = HALF * ( y(ix,iy) - y(ix-1,iy-1) )    ! D(u,y), (1,0)
         ucoef_b_xi = y(ix,iy) - y(ix-1,iy)
         ucoef_b_eta = y(ix-1,iy-1) - y(ix-1,iy)
         dux_dx_uxcoef(1,0) = ucoef_a / detxy_c
         dux_dx_uxcoef_l(1,0) = meanNeg( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dx_uxcoef_r(1,0) = meanPos( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dx_uxcoef_b(1,0) = meanNeg( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
         dux_dx_uxcoef_t(1,0) = meanPos( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
!
         !dux_dx_uxcoef(0,1) = - dux_dx_uxcoef(1,0)
         ucoef_a = HALF * ( y(ix-1,iy-1) - y(ix,iy) )    ! D(u,y), (0,1)
         ucoef_b_xi = y(ix,iy-1) - y(ix-1,iy-1)
         ucoef_b_eta = y(ix,iy-1) - y(ix,iy)
         dux_dx_uxcoef(0,1) = ucoef_a / detxy_c
         dux_dx_uxcoef_l(0,1) = meanNeg( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dx_uxcoef_r(0,1) = meanPos( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dx_uxcoef_b(0,1) = meanNeg( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
         dux_dx_uxcoef_t(0,1) = meanPos( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
!
         !dux_dx_uxcoef(0,0) = - dux_dx_uxcoef(1,1)
         ucoef_a = HALF * ( y(ix,iy-1) - y(ix-1,iy) )    ! D(u,y), (0,0)
         ucoef_b_xi = y(ix-1,iy) - y(ix,iy)
         ucoef_b_eta = y(ix,iy) - y(ix,iy-1)
         dux_dx_uxcoef(0,0) = ucoef_a / detxy_c
         dux_dx_uxcoef_l(0,0) = meanNeg( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dx_uxcoef_r(0,0) = meanPos( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dx_uxcoef_b(0,0) = meanNeg( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
         dux_dx_uxcoef_t(0,0) = meanPos( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
!
         duy_dx_uycoef(:,:) = dux_dx_uxcoef(:,:)
         duy_dx_uycoef_l(:,:) = dux_dx_uxcoef_l(:,:)
         duy_dx_uycoef_r(:,:) = dux_dx_uxcoef_r(:,:)
         duy_dx_uycoef_b(:,:) = dux_dx_uxcoef_b(:,:)
         duy_dx_uycoef_t(:,:) = dux_dx_uxcoef_t(:,:)
!
         !dux_dy_uxcoef(1,1) = HALF * ( dx_dxi - dx_deta ) / detxy_c
         ucoef_a = HALF * ( x(ix,iy-1) - x(ix-1,iy) )    ! D(x,u), (1,1)
         ucoef_b_xi = x(ix,iy-1) - x(ix-1,iy-1)
         ucoef_b_eta = x(ix-1,iy-1) - x(ix-1,iy)
         dux_dy_uxcoef(1,1) = ucoef_a / detxy_c
         dux_dy_uxcoef_l(1,1) = meanNeg( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dy_uxcoef_r(1,1) = meanPos( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dy_uxcoef_b(1,1) = meanNeg( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
         dux_dy_uxcoef_t(1,1) = meanPos( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
!
         !dux_dy_uxcoef(1,0) = - HALF * ( dx_dxi + dx_deta ) / detxy_c
         ucoef_a = HALF * ( x(ix-1,iy-1) - x(ix,iy) )    ! D(x,u), (1,0)
         ucoef_b_xi = x(ix-1,iy) - x(ix,iy)
         ucoef_b_eta = x(ix-1,iy) - x(ix-1,iy-1)
         dux_dy_uxcoef(1,0) = ucoef_a / detxy_c
         dux_dy_uxcoef_l(1,0) = meanNeg( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dy_uxcoef_r(1,0) = meanPos( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dy_uxcoef_b(1,0) = meanNeg( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
         dux_dy_uxcoef_t(1,0) = meanPos( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
!
         !dux_dy_uxcoef(0,1) = - dux_dy_uxcoef(1,0)
         ucoef_a = HALF * ( x(ix,iy) - x(ix-1,iy-1) )    ! D(x,u), (0,1)
         ucoef_b_xi = x(ix-1,iy-1) - x(ix,iy-1)
         ucoef_b_eta = x(ix,iy) - x(ix,iy-1)
         dux_dy_uxcoef(0,1) = ucoef_a / detxy_c
         dux_dy_uxcoef_l(0,1) = meanNeg( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dy_uxcoef_r(0,1) = meanPos( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dy_uxcoef_b(0,1) = meanNeg( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
         dux_dy_uxcoef_t(0,1) = meanPos( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
!
         !dux_dy_uxcoef(0,0) = - dux_dy_uxcoef(1,1)
         ucoef_a = HALF * ( x(ix-1,iy) - x(ix,iy-1) )    ! D(x,u), (0,0)
         ucoef_b_xi = x(ix,iy) - x(ix-1,iy)
         ucoef_b_eta = x(ix,iy-1) - x(ix,iy)
         dux_dy_uxcoef(0,0) = ucoef_a / detxy_c
         dux_dy_uxcoef_l(0,0) = meanNeg( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dy_uxcoef_r(0,0) = meanPos( ucoef_a, ucoef_b_xi, detxy_c, detxy_d_xi )
!                               =======     i         i          i         i
         dux_dy_uxcoef_b(0,0) = meanNeg( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
         dux_dy_uxcoef_t(0,0) = meanPos( ucoef_a, ucoef_b_eta, detxy_c, detxy_d_eta )
!                               =======     i         i          i         i
!
         duy_dy_uycoef(:,:) = dux_dy_uxcoef(:,:)
         duy_dy_uycoef_l(:,:) = dux_dy_uxcoef_l(:,:)
         duy_dy_uycoef_r(:,:) = dux_dy_uxcoef_r(:,:)
         duy_dy_uycoef_b(:,:) = dux_dy_uxcoef_b(:,:)
         duy_dy_uycoef_t(:,:) = dux_dy_uxcoef_t(:,:)
!
!  Strain tensor elements:  exx = dux_dx,  eyy = duy_dy,  exy = HALF*(dux_dy+duy_dx)
!  Quantities calculated are in the cell center (so suffix) or
!  mean values over the half-midlines (_l, _r, _b, _t)
!
         exx_uxcoef(:,:) = dux_dx_uxcoef(:,:)
         exx_uxcoef_l(:,:) = dux_dx_uxcoef_l(:,:)
         exx_uxcoef_r(:,:) = dux_dx_uxcoef_r(:,:)
         exx_uxcoef_b(:,:) = dux_dx_uxcoef_b(:,:)
         exx_uxcoef_t(:,:) = dux_dx_uxcoef_t(:,:)
!
         eyy_uycoef(:,:) = duy_dy_uycoef(:,:)
         eyy_uycoef_l(:,:) = duy_dy_uycoef_l(:,:)
         eyy_uycoef_r(:,:) = duy_dy_uycoef_r(:,:)
         eyy_uycoef_b(:,:) = duy_dy_uycoef_b(:,:)
         eyy_uycoef_t(:,:) = duy_dy_uycoef_t(:,:)
!
         exy_uxcoef(:,:) = HALF * dux_dy_uxcoef(:,:)
         exy_uxcoef_l(:,:) = HALF * dux_dy_uxcoef_l(:,:)
         exy_uxcoef_r(:,:) = HALF * dux_dy_uxcoef_r(:,:)
         exy_uxcoef_b(:,:) = HALF * dux_dy_uxcoef_b(:,:)
         exy_uxcoef_t(:,:) = HALF * dux_dy_uxcoef_t(:,:)
!
         exy_uycoef(:,:) = HALF * duy_dx_uycoef(:,:)
         exy_uycoef_l(:,:) = HALF * duy_dx_uycoef_l(:,:)
         exy_uycoef_r(:,:) = HALF * duy_dx_uycoef_r(:,:)
         exy_uycoef_b(:,:) = HALF * duy_dx_uycoef_b(:,:)
         exy_uycoef_t(:,:) = HALF * duy_dx_uycoef_t(:,:)
!
!  Physical parameters (viscosity)
!
         dt_eta = dtime * dpa(ix,iy) / dam%viscosity    ! dpa is passed to this procedure
         geff = ONE / ( ONE / dam%modshear + dt_eta )   ! effective shear modulus
!
!  Stress tensor elements
!
         sxxdev_0_l = stress(2,ix,iy)
         syydev_0_l = stress(3,ix,iy)
         sxydev_0_l = stress(4,ix,iy)
         sxxdev_0_r = stress(5,ix,iy)
         syydev_0_r = stress(6,ix,iy)
         sxydev_0_r = stress(7,ix,iy)
         sxxdev_0_b = stress(8,ix,iy)
         syydev_0_b = stress(9,ix,iy)
         sxydev_0_b = stress(10,ix,iy)
         sxxdev_0_t = stress(11,ix,iy)
         syydev_0_t = stress(12,ix,iy)
         sxydev_0_t = stress(13,ix,iy)
!
         sxxdev_const_l = geff / dam%modshear * sxxdev_0_l
         sxxdev_const_r = geff / dam%modshear * sxxdev_0_r
         sxxdev_const_b = geff / dam%modshear * sxxdev_0_b
         sxxdev_const_t = geff / dam%modshear * sxxdev_0_t
         sxxdev_uxcoef_l(:,:) = FOUR/THREE * geff * exx_uxcoef_l(:,:)
         sxxdev_uxcoef_r(:,:) = FOUR/THREE * geff * exx_uxcoef_r(:,:)
         sxxdev_uxcoef_b(:,:) = FOUR/THREE * geff * exx_uxcoef_b(:,:)
         sxxdev_uxcoef_t(:,:) = FOUR/THREE * geff * exx_uxcoef_t(:,:)
         sxxdev_uycoef_l(:,:) = -TWO/THREE * geff * eyy_uycoef_l(:,:)
         sxxdev_uycoef_r(:,:) = -TWO/THREE * geff * eyy_uycoef_r(:,:)
         sxxdev_uycoef_b(:,:) = -TWO/THREE * geff * eyy_uycoef_b(:,:)
         sxxdev_uycoef_t(:,:) = -TWO/THREE * geff * eyy_uycoef_t(:,:)
!
         syydev_const_l = geff / dam%modshear * syydev_0_l
         syydev_const_r = geff / dam%modshear * syydev_0_r
         syydev_const_b = geff / dam%modshear * syydev_0_b
         syydev_const_t = geff / dam%modshear * syydev_0_t
         syydev_uxcoef_l(:,:) = -TWO/THREE * geff * exx_uxcoef_l(:,:)
         syydev_uxcoef_r(:,:) = -TWO/THREE * geff * exx_uxcoef_r(:,:)
         syydev_uxcoef_b(:,:) = -TWO/THREE * geff * exx_uxcoef_b(:,:)
         syydev_uxcoef_t(:,:) = -TWO/THREE * geff * exx_uxcoef_t(:,:)
         syydev_uycoef_l(:,:) = FOUR/THREE * geff * eyy_uycoef_l(:,:)
         syydev_uycoef_r(:,:) = FOUR/THREE * geff * eyy_uycoef_r(:,:)
         syydev_uycoef_b(:,:) = FOUR/THREE * geff * eyy_uycoef_b(:,:)
         syydev_uycoef_t(:,:) = FOUR/THREE * geff * eyy_uycoef_t(:,:)
!
         sxy_const_l = geff / dam%modshear * sxydev_0_l
         sxy_const_r = geff / dam%modshear * sxydev_0_r
         sxy_const_b = geff / dam%modshear * sxydev_0_b
         sxy_const_t = geff / dam%modshear * sxydev_0_t
         sxy_uxcoef_l(:,:) = TWO * geff * exy_uxcoef_l(:,:)
         sxy_uxcoef_r(:,:) = TWO * geff * exy_uxcoef_r(:,:)
         sxy_uxcoef_b(:,:) = TWO * geff * exy_uxcoef_b(:,:)
         sxy_uxcoef_t(:,:) = TWO * geff * exy_uxcoef_t(:,:)
         sxy_uycoef_l(:,:) = TWO * geff * exy_uycoef_l(:,:)
         sxy_uycoef_r(:,:) = TWO * geff * exy_uycoef_r(:,:)
         sxy_uycoef_b(:,:) = TWO * geff * exy_uycoef_b(:,:)
         sxy_uycoef_t(:,:) = TWO * geff * exy_uycoef_t(:,:)
!
         ffi = ff0(ix,iy) + time * ( ff(ix,iy) - ff0(ix,iy) )
         p_const = dam%modbulk * ( ONE - vfac(ix,iy) / ffi )
         p_uxcoef(:,:) = - dam%modbulk * vfac(ix,iy) / ffi * exx_uxcoef(:,:)
         p_uycoef(:,:) = - dam%modbulk * vfac(ix,iy) / ffi * eyy_uycoef(:,:)
!
         sxx_const_l = sxxdev_const_l - p_const
         sxx_const_r = sxxdev_const_r - p_const
         sxx_const_b = sxxdev_const_b - p_const
         sxx_const_t = sxxdev_const_t - p_const
         sxx_uxcoef_l(:,:) = sxxdev_uxcoef_l(:,:) - p_uxcoef(:,:)
         sxx_uxcoef_r(:,:) = sxxdev_uxcoef_r(:,:) - p_uxcoef(:,:)
         sxx_uxcoef_b(:,:) = sxxdev_uxcoef_b(:,:) - p_uxcoef(:,:)
         sxx_uxcoef_t(:,:) = sxxdev_uxcoef_t(:,:) - p_uxcoef(:,:)
         sxx_uycoef_l(:,:) = sxxdev_uycoef_l(:,:) - p_uycoef(:,:)
         sxx_uycoef_r(:,:) = sxxdev_uycoef_r(:,:) - p_uycoef(:,:)
         sxx_uycoef_b(:,:) = sxxdev_uycoef_b(:,:) - p_uycoef(:,:)
         sxx_uycoef_t(:,:) = sxxdev_uycoef_t(:,:) - p_uycoef(:,:)
!
         syy_const_l = syydev_const_l - p_const
         syy_const_r = syydev_const_r - p_const
         syy_const_b = syydev_const_b - p_const
         syy_const_t = syydev_const_t - p_const
         syy_uxcoef_l(:,:) = syydev_uxcoef_l(:,:) - p_uxcoef(:,:)
         syy_uxcoef_r(:,:) = syydev_uxcoef_r(:,:) - p_uxcoef(:,:)
         syy_uxcoef_b(:,:) = syydev_uxcoef_b(:,:) - p_uxcoef(:,:)
         syy_uxcoef_t(:,:) = syydev_uxcoef_t(:,:) - p_uxcoef(:,:)
         syy_uycoef_l(:,:) = syydev_uycoef_l(:,:) - p_uycoef(:,:)
         syy_uycoef_r(:,:) = syydev_uycoef_r(:,:) - p_uycoef(:,:)
         syy_uycoef_b(:,:) = syydev_uycoef_b(:,:) - p_uycoef(:,:)
         syy_uycoef_t(:,:) = syydev_uycoef_t(:,:) - p_uycoef(:,:)

      end subroutine calcStrainAndStress
!
!-----------------------------------------------------------------------
!
      real(REAL8) function meanNeg( a, b, c, d )
!                          =======  i  i  i  i
!
!  Calculate  2 * Integral( (a+bx)/(c+dx) ) dx  from -1/2 to 0
!
         real(REAL8), intent(in) :: a, b, c, d
!
         if( abs(d) < EPSMB*abs(c) ) then
            meanNeg = ( a - QUART*b ) / c
         else
            meanNeg = (b/d) - TWO * ( (a/d) - (b/d)*(c/d) ) * log( ONE - HALF/(c/d) )
         end if
!
      end function meanNeg
!
!-----------------------------------------------------------------------
!
      real(REAL8) function meanPos( a, b, c, d )
!                          =======  i  i  i  i
!
!  Calculate  2 * Integral( (a+bx)/(c+dx) ) dx  from 0 to 1/2
!
         real(REAL8), intent(in) :: a, b, c, d
!
         if( abs(d) < EPSMB*abs(c) ) then
            meanPos = ( a + QUART*b ) / c
         else
            meanPos = (b/d) + TWO * ( (a/d) - (b/d)*(c/d) ) * log( ONE + HALF/(c/d) )
         end if
!
      end function meanPos
!
   end subroutine relaxViscoElastic2d
!
!=======================================================================
!
   subroutine dumpGrid( x, z )
!             ========  i  i
!
!     Dump grid to file
!
!  Arguments:
!
      real(REAL8), intent(in) :: x(:,:)         ! x-coordinates
      real(REAL8), intent(in) :: z(:,:)         ! z-coordinates
!
!  Local Variables:
!
      integer(INT4) :: ix, iz, nx, nz, ix1, ix2, ix3, iz1, iz2, iz3
      character(6), save :: fmode = 'write'
!
!  Start of Executable Code:
!  ------------------------
!
      !if( .not. lastProjectile( nionout ) .and. .not. lastProjectile() )  &
      !   return
!
      call openFile( 96, 'GRID', 'formatted', fmode )
!          ========  i      i          i         i
      fmode = 'append'
!
      nx = size( x, 1 )
      nz = size( x, 2 )
!
      write( 96, * ) 2, 2*nx*nz
      write( 96, * ) 0, 0
!
      ix1 = 1
      ix2 = nx
      ix3 = +1
!
      do iz = 1, nz
         do ix = ix1, ix2, ix3
            write( 96, '(1p,2e15.6)' ) x(ix,iz), z(ix,iz)
         end do
         ix = ix1
         ix1 = ix2
         ix2 = ix
         ix3 = -ix3
      end do
!
      iz1 = nz
      iz2 = 1
      iz3 = -1
!
      do ix = ix1, ix2, ix3
         do iz = iz1, iz2, iz3
            write( 96, * ) x(ix,iz), z(ix,iz)
         end do
         iz = iz1
         iz1 = iz2
         iz2 = iz
         iz3 = -iz3
      end do
!
      close( 96 )
!
   end subroutine dumpGrid
!
!=======================================================================
!
end module Target_Mechanics2d_mod

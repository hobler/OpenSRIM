module Geom_Mesh_mod
!
!  Manages a mesh containing the closest distance from the boundary
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, TWO, HALF, QUART

   use Param_Sn_mod,          only: sn
   use Geom_Polygon_Old_mod,  only: mcdis4, mcloc
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
   public :: mcmesh                          ! initializes the mesh
!            ======
   public :: isNearBnd                       ! checks whether position is near boundary
!            =========
!
!  Public Data Type:
!  ================
!
   integer(INT4), parameter :: NIM = 1024    ! maximum number of intervals in x- and
                                             !    z-direction
   integer(INT4), parameter :: NLXZM = 10000 ! maximum number of mesh elements

   type, public :: Mesh_type
      real(REAL8) :: x0                      ! minimum x-value of the region
      real(REAL8) :: z0                      ! minimum z-value of the region
      real(REAL8) :: dx                      ! interval length in x-direction
      real(REAL8) :: dz                      ! interval length in z-direction
      integer(INT4) :: ni                    ! actual number of intervals in x- and
                                             !    z-direction (<=NIM)
      integer(INT4) :: jx(NIM)               ! mesh index as a function of interval
                                             !    index in x-direction
      integer(INT4) :: jz(NIM)               ! mesh index as a function of interval 
                                             !    index in z-direction
      integer(INT4) :: nlx                   ! number of mesh lines in x-direction (=jx(ni))
      integer(INT4) :: nlz                   ! number of mesh lines in z-direction (=jz(ni))
      real(REAL8) :: distmin(NLXZM)          ! minimum over all points inside rectangle of
                                             !    the minimum distance from boundary
      real(REAL8) :: distmax(NLXZM)          ! maximum over all points inside rectangle of
                                             !    the minimum distance from boundary
      logical :: inside(NLXZM)               ! flag indicating whether center of rectangle
                                             !    is inside the polygon
   end type Mesh_type

contains
!
!=======================================================================
!
   logical function isNearBnd( x, z, mesh, dsput )
!                   =========  i  i   i      i
!
!     Determine whether point (x,z) is close to the region boundary
!
!  Arguments:
!
      real(REAL8), intent(in) :: x              ! x-coordinate of point to be checked
      real(REAL8), intent(in) :: z              ! z-coordinate of point to be checked
      type(Mesh_type), intent(in) :: mesh       ! mesh
      real(REAL8), intent(in) :: dsput          ! upper limit to sputtered distance
!
!  Local Variables:
!
      integer(INT4) :: imeshx, imeshz, imesh
      real(REAL8) :: dxact, dzact
!
!  Start of Executable Code:
!  ------------------------
!
!
!  Determine if (x,z) is potentially located in an interface layer. If so, 
!  determine the mesh indices. We know that the point is located in layer ir.

      dxact = x - mesh%x0
      dzact = z - mesh%z0
                                             ! check if point is really inside mesh
      if( dxact < 0 .or. dxact > mesh%ni * mesh%dx .or.  &
          dzact < 0 .or. dzact > mesh%ni * mesh%dz ) then
            error stop 'getRegion2d: the point is not inside the mesh!'
      end if
      imeshx = int( dxact/mesh%dx ) + 1
      imeshz = int( dzact/mesh%dz ) + 1
      imesh = (mesh%jz(imeshz)-1) * mesh%nlx + mesh%jx(imeshx)
                                             ! check if interface layer is potentially
                                             !    entered
      isNearBnd = (mesh%distmin(imesh) - dsput <= sn%pmax) 

   end function isNearBnd
!
!=======================================================================
!
   subroutine mcmesh( px, pz, coseta, sineta, ffpm, mesh )
!             ======  i   i     i       i      i     o
!
!     Generates a mesh (mesh) for each region using the boundaries and a minimum mesh 
!     spacing  (ffpm)
!
!  Arguments:
!
      real(REAL8), intent(in) :: px(:)          ! x-coordinates of polygon (first=last point)
      real(REAL8), intent(in) :: pz(:)          ! z-coordinates of polygon (first=last point)
      real(REAL8), intent(in) :: coseta         ! cosine of direction angle of auxiliary line
      real(REAL8), intent(in) :: sineta         ! sine of direction angle of auxiliary line
      real(REAL8), intent(in) :: ffpm           ! maximum free flight path = minimum grid 
                                                !    spacing
      type(Mesh_type) :: mesh                   ! mesh for region
!
!  Local Variables:
!
      real(REAL8), parameter :: DFAC = 0.98

      logical :: reentr(size( px ))             ! flag indicating reentrant corner
      logical :: far(NIM)                       ! flag indicating distmin is larger than the
                                                !    minimum of the interval lengths in x- 
                                                !    and z-direction
      logical :: farold, insert, inside

      integer(INT4) :: jp, ip, iloc, ii, j, jxx, jzz, np
      integer(INT4) :: ixmax, izmax, ix, iz
      integer(INT4) :: index0, index, index1, index2
      integer(INT4) :: jxold(NIM), jzold(NIM)   ! old values of jx, jz
      integer(INT4) :: ix1(NIM), iz1(NIM)
      integer(INT4) :: ix2(NIM), iz2(NIM)

      real(REAL8) :: xmax, zmax, ffpmsq, deltsq, rsqi, rsqj, distminsq, distmaxsq
      real(REAL8) :: x1, x2, x3, z1, z2, z3, dx, dz
      real(REAL8) :: distminold(NLXZM)          ! old values of distmin
      real(REAL8) :: distmaxold(NLXZM)          ! old values of distmax
      logical :: insideold(NLXZM)               ! old values of inside
      real(REAL8) :: rx(size( px ))
      real(REAL8) :: rz(size( px ))
!
!  Start of Executable Code:
!  ------------------------
!
      np = size( px ) - 1                       ! number of distinct points
!
!  Determine reentrant corners
! 
      do jp = 1, np+1

         do ip = 1, jp-1
            rx(ip) = px(ip)
            rz(ip) = pz(ip)
         end do
         do ip = jp+1, np+1
            rx(ip-1) = px(ip)
            rz(ip-1) = pz(ip)
         end do

         if( jp == 1 ) then
            rx(np) = rx(1)
            rz(np) = rz(1)
         else if( jp == np+1 ) then
            rx(1) = rx(np)
            rz(1) = rz(np)
         end if

         iloc = mcloc( px(jp), pz(jp), rx(1:np), rz(1:np), coseta, sineta )
!               =====   i       i        i         i         i       i

         reentr(jp) = (iloc == 1)

      end do
!
!  Calculate the boundaries of the mesh which will be overlayed to each region
!
      mesh%x0 = px(1)
      xmax = px(1)
      mesh%z0 = pz(1)
      zmax = pz(1)

      do ip = 1, np
         mesh%x0 = min( mesh%x0, px(ip) )
         xmax = max( xmax, px(ip) )
         mesh%z0 = min( mesh%z0, pz(ip) )
         zmax = max( zmax, pz(ip) )
      end do

      mesh%dx = xmax - mesh%x0
      mesh%dz = zmax - mesh%z0

      mesh%x0 = mesh%x0 - 0.00001 * mesh%dx
      mesh%z0 = mesh%z0 - 0.00001 * mesh%dz
      mesh%dx = 1.00002 * mesh%dx
      mesh%dz = 1.00002 * mesh%dz
! 
!  Determine maximum free flight path      
!
      ffpmsq = ffpm**2 
!
!  Initialize mesh
!
      mesh%ni = 1
      mesh%jx(1) = 1
      mesh%jz(1) = 1
      mesh%nlx = 1
      mesh%nlz = 1
      mesh%distmin(1) = ZERO
      mesh%distmax(1) = 0.25 * (mesh%dx**2 + mesh%dz**2)    ! only rough upper limit
      mesh%inside(1) = .true.                               ! only a guess

      ix1(1) = 1
      ix2(1) = 1
      iz1(1) = 1
      iz2(1) = 1
!
!  Loop over mesh refinements
!
      Refine_loop: do 
!
!  Save old data
!
         do ii = 1, mesh%ni
            jxold(ii) = mesh%jx(ii)
            jzold(ii) = mesh%jz(ii)
         end do

         do j=1,mesh%nlx*mesh%nlz
            distminold(j) = mesh%distmin(j)
            distmaxold(j) = mesh%distmax(j)
            insideold(j) = mesh%inside(j)
         end do

         mesh%ni = 2 * mesh%ni
         mesh%dx = HALF * mesh%dx
         mesh%dz = HALF * mesh%dz
         deltsq = ( DFAC * min( mesh%dx, mesh%dz ) )**2
         rsqi = QUART * (mesh%dx**2 + mesh%dz**2)
!
! Expand the mesh
!
         do ii = 1, mesh%ni
            mesh%jx(ii) = jxold((ii+1)/2)
            mesh%jz(ii) = jzold((ii+1)/2)
         end do

         do jxx = 1, mesh%nlx
            ix1(jxx) = 2*ix1(jxx) - 1
            ix2(jxx) = 2*ix2(jxx)
         end do

         do jzz = 1, mesh%nlz
            iz1(jzz) = 2*iz1(jzz) - 1
            iz2(jzz) = 2*iz2(jzz)
         end do
!
!  Insert new x-lines
!  ==================
! 
         izmax = 0
!
!  Loop over intervals
!  -------------------
!
         IntervalX_loop: do ix = 2, mesh%ni
            ! 
            ! if there is already a line between ix and ix-1, no line can be
            ! inserted. it is essential to skip the rest of the loop.
            ! izmax is the maximum index for which far(iz) has been 
            ! calculated for ix-1.
            !
            if( mesh%jx(ix-1) /= mesh%jx(ix) ) then
               izmax = 0
               cycle IntervalX_loop
            end if
!
!  Check if a line is to be inserted
! 
            x1 = mesh%x0 + (ix-1) * mesh%dx
            x2 = x1 + mesh%dx
            x3 = x1 - mesh%dx
            insert = .false.
            ! 
            ! check if far=farold.
            ! far=.true., if element (ix,iz) is inside and "far" from 
            !             boundaries
            ! farold=.true., if element (ix-1,iz) is inside and "far" 
            !                from boundaries.
            ! if any far(iz).neqv.farold(iz), a line will be inserted.
            !
            ! remark 1: index(iz,ix) = index(iz,ix-1).
            ! remark 2: if iz <= izmax, farold may be obtained from the 
            !           old far (far(iz) for ix-1).
            ! remark 3: no check is necessary for iz=1 and iz=ni, because
            !           in these cases always far=.false. and farold=.false.
            !
            index0 = mesh%jx(ix) - mesh%nlx

            do iz = 1, mesh%ni

               index = mesh%jz(iz) * mesh%nlx + index0
               if( mesh%distmin(index) >= deltsq ) then
                  far(iz) = .true.
                  cycle
               end if
               if( mesh%distmin(index) < ZERO ) then
                  far(iz) = .false.
                  cycle
               end if

               z1 = mesh%z0 + (iz-1) * mesh%dz
               z2 = z1 + mesh%dz

               if( iz <= izmax ) then
                  farold = far(iz)
               else
                  call mcdis4( distminsq, distmaxsq, inside, x3, x1, z1, z2, rsqi,  &
!                      ======      o          o        o     i   i   i   i    i
                               px(1:np+1), pz(1:np+1), reentr(1:np+1), coseta, sineta )
!                                   i           i           i            i       i
                  farold = (distminsq >= deltsq)
               end if

               call mcdis4( distminsq, distmaxsq, inside, x1, x2, z1, z2, rsqi,  &
!                   ======      o          o        o     i   i   i   i    i    
                            px(1:np+1), pz(1:np+1), reentr(1:np+1), coseta, sineta )
!                               i           i            i            i       i
               far(iz) = (distminsq >= deltsq)

               if( far(iz) .neqv. farold ) then
                  insert = .true.
                  exit
               end if
            end do

            if( .not.insert ) then
               izmax = mesh%ni
               cycle IntervalX_loop
            end if

            izmax = iz
!
!  Insert a new line between ix and ix-1
!
            ! The insertion is performed by increasing jx for indices >= ix.
            !   ix1, ix2, and distmin must also be adjusted.
            ! Remark: this point is only reached, if a line is to be 
            !         inserted.
            !
            do ii = ix, mesh%ni
               mesh%jx(ii) = mesh%jx(ii) + 1
            end do
            mesh%nlx = mesh%jx(mesh%ni)
            if( mesh%nlx * mesh%nlz > NLXZM)  exit Refine_loop

            do jxx = mesh%nlx, mesh%jx(ix)+1, -1
               ix1(jxx) = ix1(jxx-1)
               ix2(jxx) = ix2(jxx-1)
            end do
            ix2(mesh%jx(ix)) = ix2(mesh%jx(ix)-1)
            ix1(mesh%jx(ix)) = ix
            ix2(mesh%jx(ix)-1) = ix-1

            do jzz = mesh%nlz, 1, -1
               do jxx = mesh%nlx, mesh%jx(ix)+1, -1
                  index2 = (jzz-1) * mesh%nlx + jxx
                  index1 = (jzz-1) * (mesh%nlx-1) + (jxx-1)
                  mesh%distmin(index2) = mesh%distmin(index1)
                  mesh%distmax(index2) = mesh%distmax(index1)
                  mesh%inside(index2) = mesh%inside(index1)
               end do
               z1 = mesh%z0 + (iz1(jzz)-1) * mesh%dz
               z2 = mesh%z0 + iz2(jzz) * mesh%dz
               dz = z2 - z1
               do jxx = mesh%jx(ix), mesh%jx(ix)-1, -1
                  x1 = mesh%x0 + (ix1(jxx)-1) * mesh%dx
                  x2 = mesh%x0 + ix2(jxx) * mesh%dx
                  dx = x2 - x1
                  rsqj = QUART * (dx**2 + dz**2)
                  index = (jzz-1) * mesh%nlx + jxx
                  call mcdis4( mesh%distmin(index), mesh%distmax(index), mesh%inside(index),  &
!                      ======        o                    o                    o     
                               x1, x2, z1, z2, rsqj,  &
!                              i   i   i   i    i
                               px(1:np+1), pz(1:np+1), reentr(1:np+1), coseta, sineta )
!                                   i           i           i            i       i
               end do
               do jxx = mesh%jx(ix)-2, 1, -1
                  index2 = (jzz-1) * mesh%nlx + jxx
                  index1 = (jzz-1) * (mesh%nlx-1) + jxx
                  mesh%distmin(index2) = mesh%distmin(index1)
                  mesh%distmax(index2) = mesh%distmax(index1)
                  mesh%inside(index2) = mesh%inside(index1)
               end do
            end do

         end do IntervalX_loop
!
!  Insert new z-lines
!  ==================
!
         ixmax = 0                           ! the first column may never be far from the
!                                               !    boundary
!  Loop over intervals
!  -------------------
!
         IntervalZ_loop: do iz = 2, mesh%ni
            !
            ! If there is already a line between iz and iz-1, no line can be
            !   inserted. it is essential to skip the rest of the loop.
            ! ixmax is the maximum index for which far(ix) has been 
            !   calculated for iz-1.
            !
            if( mesh%jz(iz-1) /= mesh%jz(iz) ) then
               ixmax = 0
               cycle IntervalZ_loop
            end if
!
!  Check, if a line is to be inserted
! 
            z1 = mesh%z0 + (iz-1) * mesh%dz
            z2 = z1 + mesh%dz
            z3 = z1 - mesh%dz
            insert = .false.
            ! 
            ! Check, if far(ix)=farold(ix)
            !   far=.true., if element (ix,iz) is inside and "far" from 
            !               boundaries
            !   farold=.true., if element (ix,iz-1) is inside and "far" 
            !                  from boundaries.
            ! If any far(ix).neqv.farold(ix), a line will be inserted.
            !
            ! Remark 1: index(ix,iz) = index(ix,iz-1)
            ! Remark 2: if ix <= ixmax, farold may be obtained from the old
            !           far (far(ix) for iz-1)
            ! Remark 3: no check is necessary for ix=1 and ix=ni, because in
            !           these cases always far=.false. and farold=.false.
            !
            index0 = (mesh%jz(iz)-1) * mesh%nlx

            do ix = 1, mesh%ni
               index = index0 + mesh%jx(ix)
               if( mesh%distmin(index) >= deltsq ) then
                  far(ix) = .true.
                  cycle
               end if
               if( mesh%distmin(index) < ZERO ) then
                  far(ix) = .false.
                  cycle
               end if

               x1 = mesh%x0 + (ix-1) * mesh%dx
               x2 = x1 + mesh%dx
               if( ix <= ixmax ) then
                  farold = far(ix)
               else
                  call mcdis4( distminsq, distmaxsq, inside, x1, x2, z3, z1, rsqi,  & 
!                      ======      o          o        o     i   i   i   i    i    
                               px(1:np+1), pz(1:np+1), reentr(1:np+1), coseta, sineta )
!                                  i           i          i              i       i
                  farold = (distminsq >= deltsq)
               end if
               call mcdis4( distminsq, distmaxsq, inside, x1, x2, z1, z2, rsqi,  &
!                   ======      o          o        o     i   i   i   i    i    
                            px(1:np+1), pz(1:np+1), reentr(1:np+1), coseta, sineta)
!                               i           i          i              i       i
               far(ix) = (distminsq >= deltsq)

               if( far(ix) .neqv. farold ) then
                  insert = .true.
                  exit
               end if
            end do

            if( .not.insert ) then
               ixmax = mesh%ni
               cycle IntervalZ_loop
            end if

            ixmax = ix
! 
!  Insert a new line between iz and iz-1
! 
            ! The insertion is performed by increasing jz for indices >= iz.
            !   iz1, iz2, and distmin must also be adjusted.
            ! Remark: this point will only be reached, if a line is to be 
            !         inserted.
            !
            do ii = iz, mesh%ni
               mesh%jz(ii) = mesh%jz(ii) + 1
            end do
            mesh%nlz = mesh%jz(mesh%ni)
            if( mesh%nlx * mesh%nlz > NLXZM)  exit Refine_loop

            do jzz = mesh%nlz, mesh%jz(iz)+1, -1
               iz1(jzz) = iz1(jzz-1)
               iz2(jzz) = iz2(jzz-1)
            end do
            iz2(mesh%jz(iz)) = iz2(mesh%jz(iz)-1)
            iz1(mesh%jz(iz)) = iz
            iz2(mesh%jz(iz)-1) = iz - 1

            do jxx = mesh%nlx, 1, -1
               do jzz = mesh%nlz, mesh%jz(iz)+1, -1
                  index2 = (jzz-1) * mesh%nlx + jxx
                  index1 = (jzz-2) * mesh%nlx + jxx
                  mesh%distmin(index2) = mesh%distmin(index1)
                  mesh%distmax(index2) = mesh%distmax(index1)
                  mesh%inside(index2) = mesh%inside(index1)
               end do
               x1 = mesh%x0 + (ix1(jxx)-1) * mesh%dx
               x2 = mesh%x0 + ix2(jxx) * mesh%dx
               dx = x2 - x1
               do jzz = mesh%jz(iz), mesh%jz(iz)-1, -1
                  z1 = mesh%z0 + (iz1(jzz)-1) * mesh%dz
                  z2 = mesh%z0 + iz2(jzz) * mesh%dz
                  dz = z2 - z1
                  rsqj = QUART * (dx**2 + dz**2)
                  index = (jzz-1) * mesh%nlx + jxx
                  call mcdis4( mesh%distmin(index), mesh%distmax(index), mesh%inside(index),  &
!                      ======        o                    o                    o     
                               x1, x2, z1, z2, rsqj,  &
!                              i   i   i   i    i 
                               px(1:np+1), pz(1:np+1), reentr(1:np+1), coseta, sineta )
!                                  i           i            i            i       i
               end do
            end do

         end do IntervalZ_loop

         if( 2*mesh%ni > NIM .or. deltsq <= ffpmsq )  exit Refine_loop

      end do Refine_loop
!
!  End of mesh refinement loop
!
      if( 2*mesh%ni <= NIM ) then
         mesh%ni = mesh%ni / 2
         do ii = 1, mesh%ni
            mesh%jx(ii) = jxold(ii)
            mesh%jz(ii) = jzold(ii)
         end do
         mesh%nlx = mesh%jx(mesh%ni)
         mesh%nlz = mesh%jz(mesh%ni)
         mesh%dx = TWO * mesh%dx
         mesh%dz = TWO * mesh%dz
         do j = 1, mesh%nlx * mesh%nlz
            mesh%distmin(j) = distminold(j)
            mesh%distmax(j) = distmaxold(j)
            mesh%inside(j) = insideold(j)
         end do
      end if

      do j = 1, mesh%nlx * mesh%nlz
         mesh%distmin(j) = sqrt(mesh%distmin(j))
         mesh%distmax(j) = sqrt(mesh%distmax(j))
      end do
!
!  End of loop over regions
!
   end subroutine mcmesh
!
!=======================================================================
!
end module Geom_Mesh_mod

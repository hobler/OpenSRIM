module Cinterface_Mechanics1d_mod
!
!  Provides functions that relax a 1D target in response to the changes in cell contents
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, ONE
   use Array_mod,             only: alloc

   use Param_Setup_mod,       only: nelem
   use Param_Elem_mod,        only: elem
   use Param_Geom_mod,        only: geom
   use Param_Dam_mod,         only: dam
   use param_Out_mod,         only: out

   use Target_Surface1d_mod,  only: updateSurface1d

   use, intrinsic :: iso_c_binding
!
!  Default Type and Access:
!  =======================
!
   implicit none
   private
!
!  Public Procedures:
!  =================
!
   public :: relaxViscoElasticEuler1d           ! relax densities according to a viscoelastic
!            ========================           !    model on an Eulerian grid
   public :: initSurface1dFromFF                ! set surface from fill fractions
!            ===================
!
contains
!
!=======================================================================
!
   subroutine relaxViscoElasticEuler1d( z, dens, damdens, dens0, ffgeom, dpa,  &
!             ======================== io   io     io       i      io     i
                                        pressure, stress, istep )
!                                          io       io      i
!
!     Relax 1D densities according to a viscoelastic model
!
      real(REAL8), intent(inout) :: z(0:)       ! z-coordinates of vertices
      real(REAL8), intent(inout) :: dens(:,:)   ! densities
                                                !    on entry: before relaxation
                                                !    on exit: after relaxation
      real(REAL8), intent(inout) :: damdens(:,:)! damage densities
      real(REAL8), intent(in) :: dens0(:,:)     ! densities after previous relaxation
      real(REAL8), intent(inout) :: ffgeom(:)   ! geometric fill fraction:
                                                !    on entry: after last relaxation
                                                !    on exit: after this relaxation (?)
      real(REAL8), intent(in) :: dpa(:)         ! displacements per atom since previous
                                                !    relaxation
      real(REAL8), intent(inout) :: pressure(:) ! pressure
      real(REAL8), intent(inout) :: stress(:)   ! zz component of stress tensor
      integer(INT4), intent(in) :: istep        ! step counter
!
!  Local Variables:
!
      integer(INT4) :: ielem

      interface
         subroutine relax1d( nelem_c, nz_c, nother_c,  &
                             z_c, posfix_c, ff0_c, ffgeom_c, ff_c, dpa_c,  &
                             pressure_c, stress_c, other_c,  &
                             modbulk_c, modshear_c, viscosity_c, istep_c )  &
               bind(c) !, name='relax1d')
            use, intrinsic :: iso_c_binding
            integer(c_int), value :: nelem_c
            integer(c_int), value :: nz_c
            integer(c_int), value :: nother_c
            real(c_double) :: z_c(*)
            real(c_double), value :: posfix_c
            real(c_double) :: ff0_c(*)
            real(c_double) :: ffgeom_c(*)
            real(c_double) :: ff_c(*)
            real(c_double) :: dpa_c(*)
            real(c_double) :: pressure_c(*)
            real(c_double) :: stress_c(*)
            real(c_double) :: other_c(*)
            real(c_double), value :: modbulk_c
            real(c_double), value :: modshear_c
            real(c_double), value :: viscosity_c
            integer(c_int), value :: istep_c
         end subroutine relax1d
      end interface

      integer(c_int) :: nelem_c, nz_c, nother_c
      real(c_double), allocatable :: z_c(:)
      real(c_double) :: posfix_c
      real(c_double), allocatable :: ff0_c(:,:)
      real(c_double), allocatable :: ffgeom_c(:)
      real(c_double), allocatable :: ff_c(:,:)
      real(c_double), allocatable :: dpa_c(:)
      real(c_double), allocatable :: pressure_c(:)
      real(c_double), allocatable :: stress_c(:)
      real(c_double), allocatable :: other_c(:,:)
      real(c_double) :: modbulk_c
      real(c_double) :: modshear_c
      real(c_double) :: viscosity_c
      integer(c_int) :: istep_c
!
!  Start of Executable Code:
!  ------------------------
!
      nelem_c = nelem
      nz_c = size( z )
      z_c = z
      posfix_c = geom%posfix

      call alloc( ff0_c, shape( dens0 ), 'ff0_c', 'relaxViscoElasticEuler1d' )
!          =====    o          i            i                i
      call alloc( ff_c, shape( dens ), 'ff_c', 'relaxViscoElasticEuler1d' )
!          =====    o        i           i                i
      do ielem = 1, nelem
         ff0_c(ielem,:) = dens0(ielem,:) / elem(ielem)%dens
         ff_c(ielem,:) = dens(ielem,:) / elem(ielem)%dens
      end do

      ffgeom_c = ffgeom
      dpa_c = dpa
      pressure_c = pressure
      stress_c = stress
      modbulk_c = dam%modbulk
      modshear_c = dam%modshear
      viscosity_c = dam%viscosity
      if( out%lhis ) then
         nother_c = 2*nelem
         other_c = damdens
      else
         nother_c = 0
         call alloc( other_c, [1, 1, 1], 'other_c', 'relaxViscoElasticEuler1d' )
!             =====    io         i          i                 i
      end if
      istep_c = istep

      open( 96, file='relax1d_input.dat')
      if( .true. ) then                         ! long format (+ comma at end of line!)
         write( 96, '(a)' ) 'nelem, nz, nother:'
         write( 96, '(3(i5,:'',''))' )  nelem, nz_c, nother_c
         write( 96, '(a)' )  'z:'
         write( 96, '(100(es24.16e3,:'',''))' )  z_c(:)
         write( 96, '(a)' )  'ff0:'
         write( 96, '(10000(es24.16e3,:'',''))' ) ff0_c
         write( 96, '(a)' )  'ffgeom:'
         write( 96, '(10000(es24.16e3,:'',''))' ) ffgeom_c
         write( 96, '(a)' )  'ff:'
         write( 96, '(10000(es24.16e3,:'',''))' ) ff_c
         write( 96, '(a)' )  'dpa:'
         write( 96, '(10000(es24.16e3,:'',''))' ) dpa_c
         write( 96, '(a)' )  'pressure:'
         write( 96, '(10000(es24.16e3,:'',''))' ) pressure_c
         write( 96, '(a)' )  'stress:  '
         write( 96, '(40000(es24.16e3,:'',''))' ) stress_c
         write( 96, '(a)' )  'other:'
         write( 96, '(10000(es24.16e3,:'',''))' ) other_c
      end if

      call relax1d( nelem_c, nz_c, nother_c, z_c, posfix_c, ff0_c, ffgeom_c,  &
!          =======     i      i        i      i      i        i       io
                    ff_c, dpa_c, pressure_c, stress_c, other_c,  &
!                    io    io       io         io        io
                    modbulk_c, modshear_c, viscosity_c, istep_c )
!                       i           i           i          i

      do ielem = 1, nelem
         dens(ielem,:) = ff_c(ielem,:) * elem(ielem)%dens
      end do
      if( out%lhis )  damdens(:,:) = other_c(:,:)
      ffgeom = ffgeom_c
      pressure = pressure_c
      stress = stress_c

   end subroutine relaxViscoElasticEuler1d
!
!=======================================================================
!
   subroutine initSurface1dFromFF( z, ffgeom, dens, istep )
!             ===================  i    i      i      i
!
!     Set surface contours from fill fractions
!
!  Arguments:
!
      real(REAL8), intent(in) :: z(0:)          ! z-coordinates of vertices
      real(REAL8), intent(in) :: ffgeom(:)      ! geometric fill fractions
      real(REAL8), intent(in) :: dens(:,:)      ! densities
      integer(INT4), intent(in) :: istep        ! step counter
!
!  Local Variables
!
      integer(INT4) :: nz, ielem
      real(REAL8) :: zmin, zmax

      interface
         subroutine contour1d( nelem_c, nz_c, z_c, ffgeom_c, ff_c, zsurf_c, istep_c )  &
               bind(c) !, name='contour1d')
            use, intrinsic :: iso_c_binding
            integer(c_int), value :: nelem_c
            integer(c_int), value :: nz_c
            real(c_double) :: z_c(*)
            real(c_double) :: ffgeom_c(*)
            real(c_double) :: ff_c(*)
            real(c_double) :: zsurf_c(2)
            integer(c_int), value :: istep_c
         end subroutine contour1d
      end interface

      integer(c_int) :: nelem_c
      integer(c_int) :: nz_c
      real(c_double), allocatable :: z_c(:)
      real(c_double), allocatable :: ffgeom_c(:)
      real(c_double), allocatable :: ff_c(:,:)
      real(c_double) :: zsurf_c(2)
      integer(c_int) :: istep_c
!
!  Start of Executable Code:
!  ------------------------
!
      nelem_c = nelem
      nz = ubound( z, 1 )
      nz_c = nz + 1
      z_c = z
      ffgeom_c = ffgeom
      call alloc( ff_c, shape( dens ), 'ff_c', 'initSurface1dFromFF' )
!          =====    o        i           i               i
      do ielem = 1, nelem
         ff_c(ielem,:) = dens(ielem,:) / elem(ielem)%dens
      end do
      istep_c = istep

      open( 96, file='contour1d_input.dat')
      write( 96, '(a)' ) 'nz:'
      write( 96, '(i5,:'','')' )  nz_c
      write( 96, '(a)' )  'z:'
      write( 96, '(100(es24.16e3,:'',''))' )  z_c(:)
      write( 96, '(a)' )  'ffgeom:'
      write( 96, '(10000(es24.16e3,:'',''))' ) ffgeom_c
      write( 96, '(a)' )  'ff:'
      write( 96, '(10000(es24.16e3,:'',''))' ) ff_c
      close( 96 )

      call contour1d( nelem_c, nz_c, z_c, ffgeom_c, ff_c, zsurf_c, istep_c )
!          =========     i      i     i      i       i       o        i

      zmin = zsurf_c(1)
      zmax = zsurf_c(2)
      !print *, 'zmin, zmax=', zmin, zmax

      call updateSurface1d( zmin, zmax )
!          ===============   i     i

   end subroutine initSurface1dFromFF
!
!=======================================================================
!
end module Cinterface_Mechanics1d_mod

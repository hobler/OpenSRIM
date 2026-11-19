module Cinterface_Mechanics2d_mod
!
!  Provides functions that relax a 1D target in response to the changes in cell contents
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, ONE, TWO, EPSMB,  &
                                    ENUM
   use String_mod,            only: str
   use Array_mod,             only: alloc, changeShape

   use Param_Setup_mod,       only: nelem
   use Param_Elem_mod,        only: elem
   use Param_Geom_mod,        only: geom
   use Param_Dam_mod,         only: dam
   use param_Out_mod,         only: out
   use Param_Sn_mod,          only: sn

   use Target_Surface2d_mod,  only: initSurface2dFromSegments

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
   public :: relaxViscoElasticEuler2d           ! relax densities according to a viscoelastic
!            ========================           !    model on an Eulerian grid
   public :: initSurface2dFromFF                ! set surface from fill fractions
!            ===================
!
contains
!
!=======================================================================
!
   subroutine relaxViscoElasticEuler2d( x, z, dens, damdens, dens0, ffgeom, dpa,  &
!             ======================== io io   io     io       i      io     i
                                        pressure, stress, istep )
!                                          io       io      i
!
!     Relax 1D densities according to a viscoelastic model
!
      real(REAL8), intent(inout) :: x(0:)             ! x coordinates of vertices
      real(REAL8), intent(inout) :: z(0:)             ! z coordinates of vertices
      real(REAL8), intent(inout) :: dens(:,:,:)       ! densities
                                                      !    on entry: before relaxation
                                                      !    on exit: after relaxation
      real(REAL8), intent(inout) :: damdens(:,:,:)    ! damage densities
      real(REAL8), intent(in) :: dens0(:,:,:)         ! densities after previous relaxation
      real(REAL8), intent(inout) :: ffgeom(:,:)       ! geometric fill fraction:
                                                      !    on entry: after last relaxation
                                                      !    on exit: after this relaxation (?)
      real(REAL8), intent(in) :: dpa(:,:)             ! displacements per atom since previous
                                                      !    relaxation
      real(REAL8), intent(inout) :: pressure(:,:)     ! pressure
      real(REAL8), intent(inout) :: stress(:,:,:,:)   ! stress tensor
      integer(INT4), intent(in) :: istep              ! step counter
!
!  Local Variables:
!
      integer(INT4) :: iz, ielem, ix
      integer :: rss
      real(REAL8) :: sumff0, sumff

      interface
         subroutine relax2d( nelem_c, nx_c, nz_c, nother_c,  &
                             x_c, z_c, posfix_c, ff0_c, ffgeom_c, ff_c, dpa_c,  &
                             pressure_c, stress_c, other_c,  &
                             modbulk_c, modshear_c, viscosity_c, gammasurf_c, istep_c )  &
               bind(c) !, name='relax2d')
            use, intrinsic :: iso_c_binding
            integer(c_int), value :: nelem_c
            integer(c_int), value :: nx_c
            integer(c_int), value :: nz_c
            integer(c_int), value :: nother_c
            real(c_double) :: x_c(*)
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
            real(c_double), value :: gammasurf_c
            integer(c_int), value :: istep_c
         end subroutine relax2d
      end interface

      integer(c_int) :: nelem_c, nx_c, nz_c, nother_c
      real(c_double), allocatable :: x_c(:)
      real(c_double), allocatable :: z_c(:)
      real(c_double) :: posfix_c
      real(c_double), allocatable :: ff0_c(:,:,:)
      real(c_double), allocatable :: ffgeom_c(:,:)
      real(c_double), allocatable :: ff_c(:,:,:)
      real(c_double), allocatable :: dpa_c(:,:)
      real(c_double), allocatable :: pressure_c(:,:)
      real(c_double), allocatable :: stress_c(:,:,:,:)
      real(c_double), allocatable :: other_c(:,:,:)
      real(c_double) :: modbulk_c
      real(c_double) :: modshear_c
      real(c_double) :: viscosity_c
      real(c_double) :: gammasurf_c
      integer(c_int) :: istep_c
!
!  Start of Executable Code:
!  ------------------------
!
!  Move SURF file into istep-1 folder
!
      if( istep > 1 )  &
            call execute_command_line( 'mv SURF ' // str( istep-1 ) // ' 2>/dev/null' )
!
!  Copy to C tape variables
!
      nelem_c = nelem
      nx_c = size( x )
      nz_c = size( z )
      x_c = x
      z_c = z
      posfix_c = geom%posfix
      call alloc( ff0_c, shape( dens0 ), 'ff0_c', 'relaxViscoElasticEuler2d' )
!          =====    o          i            i                i
      call alloc( ff_c, shape( dens ), 'ff_c', 'relaxViscoElasticEuler2d' )
!          =====    o        i           i                i
      do ielem = 1, nelem
         ff0_c(ielem,:,:) = dens0(ielem,:,:) / elem(ielem)%dens
         ff_c(ielem,:,:) = dens(ielem,:,:) / elem(ielem)%dens
      end do

      ffgeom_c = ffgeom
      dpa_c = dpa
      pressure_c = pressure
      stress_c = stress
      modbulk_c = dam%modbulk
      modshear_c = dam%modshear
      viscosity_c = dam%viscosity
      gammasurf_c = dam%gammasurf
      if( out%lhis ) then
         nother_c = 2*nelem
         other_c = damdens
      else
         nother_c = 0
         call alloc( other_c, [1, 1, 1], 'other_c', 'relaxViscoElasticEuler2d' )
!             =====    io         i          i                 i
      end if
      istep_c = istep

      open( 96, file='relax2d_input.dat')
      if( .true. ) then                         ! long format (+ comma at end of line!)
         write( 96, '(a)' ) 'nelem, nx, nz, nother:'
         write( 96, '(4(i5,:'',''))' )  nelem, nx_c, nz_c, nother_c
         write( 96, '(a)' )  'x:'
         write( 96, '(100(es24.16e3,:'',''))' )  x_c(:)
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
      else                                      ! short format
         write( 96, '(a)' ) 'nelem, nx, nz, nother:'
         write( 96, '(4(i5,:'',''))' )  nelem_c, nx_c, nz_c, nother_c
         write( 96, '(a)' )  'x:'
         write( 96, '(100(f8.2,:'',''))' )  x_c(:)
         write( 96, '(a)' )  'z:'
         write( 96, '(100(f8.2,:'',''))' )  z_c(:)
         write( 96, '(a)' )  'ff0_elem:'
         write( 96, '(10000(f8.5,:'',''))' ) ff0_c
         !do iz = 1, nz_c-1
         !   write( 96, '(100f8.5)' ) ff0_c(:,iz)
         !end do
         write( 96, '(a)' )  'ffgeom:'
         write( 96, '(10000(f8.5,:'',''))' ) ffgeom_c
         !do iz = 1, nz_c-1
         !   write( 96, '(100f8.5)' ) ffgeom_c(:,iz)
         !end do
         write( 96, '(a)' )  'ff:'
         write( 96, '(10000(f8.5,:'',''))' ) ff_c
         !do iz = 1, nz_c-1
         !   write( 96, '(100f8.5)' ) ff_c(:,iz)
         !end do
         write( 96, '(a)' )  'dpa:'
         write( 96, '(10000(f8.5,:'',''))' ) dpa_c
         !do iz = 1, nz_c-1
         !   write( 96, '(100f8.5)' ) dpa_c(:,iz)
         !end do
         write( 96, '(a)' )  'pressure:'
         write( 96, '(10000(es13.6,:'',''))' ) pressure_c
         write( 96, '(a)' )  'stress:  '
         write( 96, '(40000(es13.6,:'',''))' ) stress_c
         write( 96, '(a)' )  'other:'
         write( 96, '(10000(es13.6,:'',''))' ) other_c
      end if
      close( 96 )

      ! debug
      if( .false. ) then
         ix = (25 - x(0)) / (x(1) - x(0)) + 1
         iz = (225 - z(0)) / (z(1) - z(0)) + 1
         write( 10, '(a)' )  'step ' // str( istep ) // ' (before): ff0(Xe,25,225)=' //  &
                             str( ff0_c(1,ix,iz) )
         write( 10, '(a)' )  'step ' // str( istep ) // ' (before): ff(Xe,25,225) =' //  &
                             str( ff_c(1,ix,iz) )
      end if

      !print '(10000f8.5)', ff0_c(:,:)
      !print '(a, 10000f8.5)', 'ffgeom=', ffgeom_c(:,:)
      !print *, 'ff=    ', ff_c(:,:)
      !print *, 'dpa=    ', dpa_c(:,:)

      sumff0 = sum( ff0_c(:,:,:) )
      sumff = sum( ff_c(:,:,:) )
      rss = systMemUsage()
      call relax2d( nelem_c, nx_c, nz_c, nother_c, x_c, z_c, posfix_c, ff0_c, ffgeom_c,  &
!          =======     i      i     i       i       i    i      i        i       io
                    ff_c, dpa_c, pressure_c, stress_c, other_c,  &
!                    io    io       io         io        io
                    modbulk_c, modshear_c, viscosity_c, gammasurf_c, istep_c )
!                       i           i           i            i          i
      ! debug
      if( .false. ) then
         write( 10, '(a)' )  'step ' // str( istep ) // ' (after) : ff0(Xe,25,225)=' //  &
                             str( ff0_c(1,ix,iz) )
         write( 10, '(a)' )  'step ' // str( istep ) // ' (after) : ff(Xe,25,225) =' //  &
                             str( ff_c(1,ix,iz) )
         write( 10, '()' )
      end if

      print *, 'FORTRAN Memory before relax2d:', rss
      rss = systMemUsage()
      print *, 'FORTRAN Memory after relax2d:', rss
      !print '(a, 10f8.5)', 'ff0=   ', ff0_c(:10)
      !print '(a, 10f8.5)', 'ffreom=', ffgeom_c(:10)
      !print '(a, 10f8.5)', 'ff=    ', ff_c(:10)

      do ielem = 1, nelem
         dens(ielem,:,:) = ff_c(ielem,:,:) * elem(ielem)%dens
      end do
      if( out%lhis )  damdens(:,:,:) = other_c(:,:,:)
      ffgeom = ffgeom_c
      pressure = pressure_c
      stress = stress_c

      print *, 'sum(ff0) before relax2d:', sumff0
      print *, 'sum(ff) before relax2d:', sumff
      print *, 'sum(ff) after relax2d:', sum( ff_c(:,:,:) )

      call execute_command_line( 'mv relax2d_input.dat ' // str( istep ) )

   end subroutine relaxViscoElasticEuler2d
!
!=======================================================================
!
   subroutine initSurface2dFromFF( x, z, ffgeom, dens, istep )
!             ===================  i  i    i      i      i
!
!     Set surface contours from fill fractions
!
!  Arguments:
!
      real(REAL8), intent(in) :: x(0:)          ! x coordinates of grid lines
      real(REAL8), intent(in) :: z(0:)          ! z coordinates of grid lines
      real(REAL8), intent(in) :: ffgeom(:,:)    ! geometric fill fractions
      real(REAL8), intent(in) :: dens(:,:,:)    ! densities
      integer(INT4), intent(in) :: istep        ! step counter
!
!  Local Variables
!
      integer(INT4) :: nx, nz, ielem
      real(REAL8), pointer :: segments(:,:,:)

      integer :: rss

      interface
         subroutine contour2d( nelem_c, nx_c, nz_c, x_c, z_c, ffgeom_c, ff_c, ns_c, s_ptr_c,  &
                               istep_c )  &
               bind(c)
            use, intrinsic :: iso_c_binding
            integer(c_int), value :: nelem_c
            integer(c_int), value :: nx_c
            integer(c_int), value :: nz_c
            real(c_double) :: x_c(*)
            real(c_double) :: z_c(*)
            real(c_double) :: ffgeom_c(*)
            real(c_double) :: ff_c(*)
            integer(c_int) :: ns_c
            type(c_ptr) :: s_ptr_c
            integer(c_int), value :: istep_c
         end subroutine contour2d
      end interface

      integer(c_int) :: nelem_c
      integer(c_int) :: nx_c
      integer(c_int) :: nz_c
      real(c_double), allocatable :: x_c(:)
      real(c_double), allocatable :: z_c(:)
      real(c_double), allocatable :: ffgeom_c(:,:)
      real(c_double), allocatable :: ff_c(:,:,:)
      integer(c_int) :: ns_c                       ! 1/4 size of *s_ptr_c
      type(c_ptr) :: s_ptr_c
      integer(c_int) :: istep_c
!
!  Start of Executable Code:
!  ------------------------
!
      nelem_c = nelem
      nx = ubound( x, 1 )
      nz = ubound( z, 1 )
      nx_c = nx + 1
      nz_c = nz + 1
      x_c = x
      z_c = z
      ffgeom_c = ffgeom
      call alloc( ff_c, shape( dens ), 'ff_c', 'initSurface2dFromFF' )
!          =====    o        i           i               i
      do ielem = 1, nelem
         ff_c(ielem,:,:) = dens(ielem,:,:) / elem(ielem)%dens
      end do
      istep_c = istep

      open( 96, file='contour2d_input.dat')
      write( 96, '(a)' ) 'nx, nz:'
      write( 96, '(2(i5,:'',''))' )  nx_c, nz_c
      write( 96, '(a)' )  'x:'
      write( 96, '(100(es24.16e3,:'',''))' )  x_c(:)
      write( 96, '(a)' )  'z:'
      write( 96, '(100(es24.16e3,:'',''))' )  z_c(:)
      write( 96, '(a)' )  'ffgeom:'
      write( 96, '(10000(es24.16e3,:'',''))' ) ffgeom_c
      write( 96, '(a)' )  'ff:'
      write( 96, '(10000(es24.16e3,:'',''))' ) ff_c
      close( 96 )

      rss = systMemUsage()
      call contour2d( nelem_c, nx_c, nz_c, x_c, z_c, ffgeom_c, ff_c, ns_c, s_ptr_c, istep_c )
!          =========     i      i     i      i   i      i       i     o       o        i
      print *, 'FORTRAN Memory before contour2d:', rss
      rss = systMemUsage()
      print *, 'FORTRAN Memory after contour2d:', rss

      call c_f_pointer( s_ptr_c, segments, [2,2,ns_c] )
!
!  Initialize the surface
!
      call initSurface2dFromSegments( x(:), z(:), ffgeom(:,:), segments(:,:,:) )
!          =========================   i     i        i              i
!      call initSurface2dFromSegmentsOld( x(:), z(:), ffgeom(:,:), segments(:,:,:) )
!!          ============================   i     i        i              i

   end subroutine initSurface2dFromFF
!
!=======================================================================
!
    function systMemUsage( ) result( valueRSS )
       ! Ref: https://stackoverflow.com/questions/22028571/track-memory-usage-in-fortran-90
        character(len=80) :: line
        integer ::  ios, fu, valueRSS
        valueRSS = -1

        open( newunit=fu, file='/proc/self/status', action='read' )
        do
            read( fu, '(a)', iostat=ios )  line
            if( ios /= 0 )  exit
            if( line(1:6) == 'VmRSS:' ) then
                read( line(7:), * )  valueRSS
                exit
            endif
        enddo
        close( fu )
    end
!
!=======================================================================
!
end module Cinterface_Mechanics2d_mod

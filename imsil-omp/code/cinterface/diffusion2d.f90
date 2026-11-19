module Cinterface_Diffusion2d_mod
!
!  Provides functions that relax a 1D target in response to the changes in cell contents
!
!  Used Modules:
!  ============
!
   use Const_mod,             only: INT4, REAL8,  &
                                    ZERO, ONE, TWO, EPSMB
   use String_mod,            only: str
   use Array_mod,             only: alloc

   use Param_Setup_mod,       only: nelem
   use Param_Elem_mod,        only: elem
   use Param_Geom_mod,        only: geom
   use Param_Dam_mod,         only: dam
   use param_Out_mod,         only: out
   use Param_Sn_mod,          only: sn

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
   public :: diffuseCell2dReg                   ! diffuse densities
!            ================
!
contains
!
!=======================================================================
!
   subroutine diffuseCell2dReg( x, z, ffgeom, dens, dpa, istep )
!             ================  i  i    i      io    i     i
!
!     Relax 1D densities according to a viscoelastic model
!
      real(REAL8), intent(in) :: x(0:)          ! x coordinates of vertices
      real(REAL8), intent(in) :: z(0:)          ! z coordinates of vertices
      real(REAL8), intent(inout) :: dens(:,:)   ! densities
                                                !    on entry: before diffusion
                                                !    on exit: after diffusion
      real(REAL8), intent(inout) :: ffgeom(:)   ! geometric fill fraction
      real(REAL8), intent(in) :: dpa(:)         ! displacements per atom since previous
                                                !    relaxation
      integer(INT4), intent(in) :: istep        ! step counter
!
!  Local Variables:
!
      interface
         subroutine diffuse2d( nelem_c, nx_c, nz_c,  &
                               x_c, z_c, ffgeom_c, dens_c, dpa_c,  &
                               diffusivity_c, istep_c ) bind(c)
            use, intrinsic :: iso_c_binding
            integer(c_int), value :: nelem_c
            integer(c_int), value :: nx_c
            integer(c_int), value :: nz_c
            real(c_double) :: x_c(*)
            real(c_double) :: z_c(*)
            real(c_double) :: ffgeom_c(*)
            real(c_double) :: dens_c(*)
            real(c_double) :: dpa_c(*)
            real(c_double) :: diffusivity_c(*)
            integer(c_int), value :: istep_c
         end subroutine diffuse2d
      end interface

      integer(c_int) :: nelem_c, nx_c, nz_c
      real(c_double), allocatable :: x_c(:)
      real(c_double), allocatable :: z_c(:)
      real(c_double), allocatable :: dens_c(:,:)
      real(c_double), allocatable :: ffgeom_c(:)
      real(c_double), allocatable :: dpa_c(:)
      real(c_double), allocatable :: diffusivity_c(:)
      integer(c_int) :: istep_c

      integer :: rss
!
!  Start of Executable Code:
!  ------------------------
!
      nelem_c = nelem
      nx_c = size( x )
      nz_c = size( z )
      x_c = x
      z_c = z
      ffgeom_c = ffgeom
      dens_c = dens
      dpa_c = dpa
      diffusivity_c = dam%diffusivity
      istep_c = istep

      open( 96, file='diffuse2d_input.dat')
      if( .true. ) then                         ! long format (+ comma at end of line!)
         write( 96, '(a)' ) 'nelem, nx, nz:'
         write( 96, '(3(i5,:'',''))' )  nelem, nx_c, nz_c
         write( 96, '(a)' )  'x:'
         write( 96, '(100(es24.16e3,:'',''))' )  x_c(:)
         write( 96, '(a)' )  'z:'
         write( 96, '(100(es24.16e3,:'',''))' )  z_c(:)
         write( 96, '(a)' )  'ffgeom:'
         write( 96, '(10000(es24.16e3,:'',''))' ) ffgeom_c
         write( 96, '(a)' )  'dens:'
         write( 96, '(10000(es24.16e3,:'',''))' ) dens_c
         write( 96, '(a)' )  'dpa:'
         write( 96, '(10000(es24.16e3,:'',''))' ) dpa_c
      else                                      ! short format
         write( 96, '(a)' ) 'nelem, nx, nz:'
         write( 96, '(3(i5,:'',''))' )  nelem_c, nx_c, nz_c
         write( 96, '(a)' )  'x:'
         write( 96, '(100(f8.2,:'',''))' )  x_c(:)
         write( 96, '(a)' )  'z:'
         write( 96, '(100(f8.2,:'',''))' )  z_c(:)
         write( 96, '(a)' )  'ffgeom:'
         write( 96, '(10000(f8.5,:'',''))' ) ffgeom_c
         write( 96, '(a)' )  'dens:'
         write( 96, '(10000(f8.5,:'',''))' ) dens_c
         write( 96, '(a)' )  'dpa:'
         write( 96, '(10000(f8.5,:'',''))' ) dpa_c
      end if
      close( 96 )

      rss = systMemUsage()

      call diffuse2d( nelem_c, nx_c, nz_c, x_c, z_c, ffgeom_c, dens_c,  &
!          =========     i      i     i     i    i      i        io
                      dpa_c, diffusivity_c, istep_c )
!                       i          i          i

      print *, 'FORTRAN Memory before diffuse2d:', rss
      rss = systMemUsage()
      print *, 'FORTRAN Memory after diffuse2d:', rss
      !print '(a, 10f8.5)', 'ff0=   ', ff0_c(:10)
      !print '(a, 10f8.5)', 'ffgeom=', ffgeom_c(:10)
      !print '(a, 10f8.5)', 'ff=    ', ff_c(:10)

      dens = dens_c

      call execute_command_line( 'mv diffuse2d_input.dat ' // str( istep ) )

   end subroutine diffuseCell2dReg
!
!=======================================================================
!
    function systMemUsage() result(valueRSS)
       ! Ref: https://stackoverflow.com/questions/22028571/track-memory-usage-in-fortran-90
        character(len=80) :: line
        integer ::  ios, fu, valueRSS
        valueRSS = -1

        open(newunit=fu, file='/proc/self/status', action='read')
        do
            read(fu, '(a)',iostat=ios ) line
            if(ios /=0) exit
            if(line(1:6) == 'VmRSS:') then
                read(line(7:), *) valueRSS
                exit
            endif
        enddo
        close(fu)
    end
!
!=======================================================================
!
end module Cinterface_Diffusion2d_mod

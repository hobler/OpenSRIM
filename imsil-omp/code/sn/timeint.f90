module Sn_TimeInt_mod
!
!  Calculate the time integral using the universal ZBL interatomic potential
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8,  &
                              ZERO, BIG

   use Param_Setup_mod, only: nelem
   use Param_Sn_mod,    only: sn

   use Open_mod,        only: openFile
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
   public :: initTimeInt                        ! initialize interpolation coefficients
!            ===========
   public :: timeint                            ! calculate angle and time integral
!            =======
!
!  Private Data:
!  ============
!
   integer(INT4) :: ne
   integer(INT4) :: np
!
!  Grid
!                                                   
   real(REAL8), allocatable :: et(:)      ! (ne)
   real(REAL8), allocatable :: pt(:)      ! (np)
   real(REAL8) :: pcutoff
!
!  Interpolation coefficients
!
   real(REAL8), allocatable :: b00(:,:)   ! (ne,np)
   real(REAL8), allocatable :: b10(:,:)   ! (ne,np)
   real(REAL8), allocatable :: b20(:,:)   ! (ne-1,np-1)
   real(REAL8), allocatable :: b30(:,:)   ! (ne-1,np-1)
   real(REAL8), allocatable :: b01(:,:)   ! (ne,np)
   real(REAL8), allocatable :: b02(:,:)   ! (ne-1,np-1)
   real(REAL8), allocatable :: b03(:,:)   ! (ne-1,np-1)
   real(REAL8), allocatable :: b11(:,:)   ! (ne-1,np-1)
   real(REAL8), allocatable :: b21(:,:)   ! (ne-1,np-1)
   real(REAL8), allocatable :: b31(:,:)   ! (ne-1,np-1)
   real(REAL8), allocatable :: b12(:,:)   ! (ne-1,np-1)
   real(REAL8), allocatable :: b13(:,:)   ! (ne-1,np-1)
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initTimeInt( )
!             ===========
!
!     Initialize the interpolation coefficients
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      character(80) :: scatfile

      integer(INT4) :: ie, ip

      real(REAL8) :: f00, f10, f01, f11, fx00, fx10, fx01, fx11
      real(REAL8) :: fy00, fy10, fy01, fy11, s
!
!  Start of Executable Code:
!  ------------------------
!
      scatfile = 'scattab/ZBLuniv'
!
!  Determine table dimensions
!
      call openFile( 10, scatfile, 'formatted', 'read' )
!          ========  i     i          i         i

      read( 10, *, err=900 )  ne, np
      close( 10 )
!
!  Allocate arrays
!
      allocate( et( ne ) )
      allocate( pt( np ) )

      allocate( b00( ne, np ) )
      allocate( b10( ne, np ) )
      allocate( b20( ne-1, np-1 ) )
      allocate( b30( ne-1, np-1 ) )
      allocate( b01( ne, np ) )
      allocate( b02( ne-1, np-1 ) )
      allocate( b03( ne-1, np-1 ) )
      allocate( b11( ne-1, np-1 ) )
      allocate( b21( ne-1, np-1 ) )
      allocate( b31( ne-1, np-1 ) )
      allocate( b12( ne-1, np-1 ) )
      allocate( b13( ne-1, np-1 ) )
!
!  Load the tables (b00=tau,  b10=dtau/de,  b01=dtau/dp) skipping the thetas
!
      call openFile( 10, scatfile, 'formatted', 'read' )
!          ========  i     i          i         i
      read( 10,*,err=900) ne, np
      read( 10,*,err=900) (et(ie), ie=1,ne)
      read( 10,*,err=900) (pt(ip), ip=1,np)
      read( 10,*,err=900) pcutoff
      read( 10,*,err=900) ((b00(ie,ip), ie=1,ne), ip=1,np)
      read( 10,*,err=900) ((b10(ie,ip), ie=1,ne), ip=1,np)
      read( 10,*,err=900) ((b01(ie,ip), ie=1,ne), ip=1,np)
      read( 10,*,err=900) ((b00(ie,ip), ie=1,ne), ip=1,np)
      read( 10,*,err=900) ((b10(ie,ip), ie=1,ne), ip=1,np)
      read( 10,*,err=900) ((b01(ie,ip), ie=1,ne), ip=1,np)
      close (10)
!
!  Initialize the parameters for bicubic interpolation (aij, bij)
!
      do ie = 1, ne-1                           ! ie and ip must be increasing in these loops
         do ip = 1, np-1

            f00 = b00(ie,ip)
            f01 = b00(ie,ip+1)
            f10 = b00(ie+1,ip)
            f11 = b00(ie+1,ip+1)
            fx00 = b10(ie,ip) * (et(ie+1)-et(ie))
            fx01 = b10(ie,ip+1) * (et(ie+1)-et(ie))
            fx10 = b10(ie+1,ip) * (et(ie+1)-et(ie))
            fx11 = b10(ie+1,ip+1) * (et(ie+1)-et(ie))
            fy00 = b01(ie,ip) * (pt(ip+1)-pt(ip))
            fy01 = b01(ie,ip+1) * (pt(ip+1)-pt(ip))
            fy10 = b01(ie+1,ip) * (pt(ip+1)-pt(ip))
            fy11 = b01(ie+1,ip+1) * (pt(ip+1)-pt(ip))

            b00(ie,ip) = f00
            b10(ie,ip) = fx00
            b01(ie,ip) = fy00

            b30(ie,ip) = fx10 - f10-f10 + f00+f00 + fx00
            b20(ie,ip) = f10 - f00 - fx00 - b30(ie,ip)
            b03(ie,ip) = fy01 - f01-f01 + f00+f00 + fy00
            b02(ie,ip) = f01 - f00 - fy00 - b03(ie,ip)

            s = f10 + f01 - f11 - f00
            b11(ie,ip) = fx01 + fy10 - fx00 - fy00 + s
            s = s + s
            b13(ie,ip) = fy11 - fy01 + fy10 - fy00 + s
            b12(ie,ip) = fx01 - fx00 - b11(ie,ip) - b13(ie,ip)
            b31(ie,ip) = fx11 - fx10 + fx01 - fx00 + s
            b21(ie,ip) = fy10 - fy00 - b11(ie,ip) - b31(ie,ip)

         end do
      end do

      return

900   error stop 'Error during reading scattering table file'

   end subroutine initTimeInt
!
!=======================================================================
!
   subroutine timeint( tau, e, p )
!             =======   o   i  i
!
!     Calculates the scattering angle in the center-of-mass system and the
!     time integral
!
!     -) uses bicubic interpolation in a table for cot(th/2)
!     -) the table is read from the scattering table file scatfile
!
!  Arguments:
!
      real(REAL8), intent(out) :: tau           ! time integral [rnorm]
      real(REAL8), intent(in) :: e              ! ion energy [enorm]
      real(REAL8), intent(in) :: p              ! impact parameter [rnorm]
!
!  Local Variables:
!
      integer(INT4) :: il, im, iu, ie, ip
      real(REAL8) :: eq, pq
!
!  Start of Executable Code:
!  ------------------------
!
!  no deflection if impact parameter larger than pcutoff
!
      if( p >= pcutoff ) then
         tau = ZERO
         return
      end if
!  
!  binary search in the pt-vector
!  
      il = 1
      iu = np
 30   im = (il+iu) / 2
      if( p <= pt(im) ) then
         iu = im
      else
         il = im
      end if
      if( iu-il > 1 )  goto 30
      ip = il
!  
!  reduced impact parameter
!  
      pq = (p - pt(ip)) / (pt(ip+1) - pt(ip))
!  
!  binary search in the e-vector
!  
      il = 1
      iu = ne
 40   im = (il+iu)/2
      if( e <= et(im) ) then
         iu = im
      else
         il = im
      end if
      if( iu-il > 1 )  goto 40
      ie = il
!
!  reduced energy
!
      eq = (e - et(ie)) / (et(ie+1) - et(ie))
!
!  bicubic interpolation of tau in the table
!
      if( e < et(1) ) then
         tau = ( (b03(1,ip)*pq+b02(1,ip))*pq + b01(1,ip) ) * pq + b00(1,ip)
      else
         tau = (  ( (b13(ie,ip)*eq+b03(ie,ip))*pq + &
                    (b12(ie,ip)*eq+b02(ie,ip)) ) * pq + &
                  ( (b31(ie,ip)*eq+b21(ie,ip))*eq+b11(ie,ip) ) * eq + &
                  b01(ie,ip)  )  *  pq    +  &
               (  ( (b30(ie,ip)*eq+b20(ie,ip))*eq + &
                    b10(ie,ip) ) * eq + b00(ie,ip)  )
      end if

   end subroutine timeint
!
!=======================================================================
!
end module Sn_TimeInt_mod

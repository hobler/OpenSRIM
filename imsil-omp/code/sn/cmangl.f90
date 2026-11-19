module Sn_Cmangl_mod
!
!  Calculate the scattering angle in the center-of-mass system
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
   public :: initCmangl             ! initialize interpolation coefficients
!            ==========
   public :: cmangl                 ! calculate angle and time integral
!            ======
!
!  Private Data:
!  ============
!
!  Table index
!
   integer(INT4), allocatable :: isca(:,:)   ! (nelem,nelem)
!
!  Table dimensions
!
   integer(INT4), allocatable :: ne(:)       ! (nsca)
   integer(INT4), allocatable :: np(:)       ! (nsca)
!
!  Grid
!                                                   
   real(REAL8), allocatable :: etab(:,:)     ! (nemax,nsca)
   real(REAL8), allocatable :: ptab(:,:)       ! (npmax,nsca)
   real(REAL8), allocatable :: eminln(:)     ! (nsca)
   real(REAL8), allocatable :: pcutoff(:)    ! (nsca)
!
!  Interpolation coefficients
!
   real(REAL8), allocatable :: a00(:,:,:)    ! (nemax,npmax,nsca)
   real(REAL8), allocatable :: a10(:,:,:)    ! (nemax,npmax,nsca)
   real(REAL8), allocatable :: a20(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: a30(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: a01(:,:,:)    ! (nemax,npmax,nsca)
   real(REAL8), allocatable :: a02(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: a03(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: a11(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: a21(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: a31(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: a12(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: a13(:,:,:)    ! (nemax-1,npmax-1,nsca)
                                                
   real(REAL8), allocatable :: b00(:,:,:)    ! (nemax,npmax,nsca)
   real(REAL8), allocatable :: b10(:,:,:)    ! (nemax,npmax,nsca)
   real(REAL8), allocatable :: b20(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: b30(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: b01(:,:,:)    ! (nemax,npmax,nsca)
   real(REAL8), allocatable :: b02(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: b03(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: b11(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: b21(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: b31(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: b12(:,:,:)    ! (nemax-1,npmax-1,nsca)
   real(REAL8), allocatable :: b13(:,:,:)    ! (nemax-1,npmax-1,nsca)
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine initCmangl( )
!             ==========
!
!     Initialize the interpolation coefficients
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem1, ielem2, jelem1, jelem2, ie, ip, ne1, nemax, np1, npmax
      integer(INT4) :: is, nsca

      real(REAL8) :: f00, f10, f01, f11, fx00, fx10, fx01, fx11
      real(REAL8) :: fy00, fy10, fy01, fy11, s
!
!  Start of Executable Code:
!  ------------------------
!
!  Determine table dimensions
!
      nsca = 0
      nemax = 0
      npmax = 0

      do_i1a: do ielem1 = 1, nelem
         do_i2a: do ielem2 = 1, nelem
                                                ! new scattering table file ?
            if( sn%lquad(ielem2,ielem1) )  cycle

            do_j1a: do jelem1 = 1, nelem
               do_j2a: do jelem2 = 1, nelem
                  if( ielem1==jelem1 .and. ielem2==jelem2 )  exit do_j1a
                  if( sn%scatfile(ielem2,ielem1) == sn%scatfile(jelem2,jelem1) )  cycle do_i2a
               end do do_j2a
            end do do_j1a

            nsca = nsca + 1

            call openFile( 10, sn%scatfile(ielem2,ielem1), 'formatted', 'read' )
!                ========  i             i                    i         i

            read( 10,*,err=900) ne1, np1
            close (10)
            nemax = max (ne1,nemax)
            npmax = max (np1,npmax)

         end do do_i2a
      end do do_i1a
!
!  Allocate arrays
!
      allocate( isca( nelem, nelem ) )

      allocate( ne( nsca ) )
      allocate( np( nsca ) )
             
      allocate( etab( nemax, nsca ) )
      allocate( ptab( npmax, nsca ) )
      allocate( eminln( nsca ) )
      allocate( pcutoff( nsca ) )

      allocate( a00( nemax, npmax, nsca ) )
      allocate( a10( nemax, npmax, nsca ) )
      allocate( a20( nemax-1, npmax-1, nsca ) )
      allocate( a30( nemax-1, npmax-1, nsca ) )
      allocate( a01( nemax, npmax, nsca ) )
      allocate( a02( nemax-1, npmax-1, nsca ) )
      allocate( a03( nemax-1, npmax-1, nsca ) )
      allocate( a11( nemax-1, npmax-1, nsca ) )
      allocate( a21( nemax-1, npmax-1, nsca ) )
      allocate( a31( nemax-1, npmax-1, nsca ) )
      allocate( a12( nemax-1, npmax-1, nsca ) )
      allocate( a13( nemax-1, npmax-1, nsca ) )
             
      allocate( b00( nemax, npmax, nsca ) )
      allocate( b10( nemax, npmax, nsca ) )
      allocate( b20( nemax-1, npmax-1, nsca ) )
      allocate( b30( nemax-1, npmax-1, nsca ) )
      allocate( b01( nemax, npmax, nsca ) )
      allocate( b02( nemax-1, npmax-1, nsca ) )
      allocate( b03( nemax-1, npmax-1, nsca ) )
      allocate( b11( nemax-1, npmax-1, nsca ) )
      allocate( b21( nemax-1, npmax-1, nsca ) )
      allocate( b31( nemax-1, npmax-1, nsca ) )
      allocate( b12( nemax-1, npmax-1, nsca ) )
      allocate( b13( nemax-1, npmax-1, nsca ) )
!
!  Load the tables (a00=cot(th/2),  a10=dcot(th/2)/de,  a01=dcot(th/2)/dp)
!                  (b00=tau,  b10=dtau/de,  b01=dtau/dp)
!
      is = 0

      do_i1b: do ielem1 = 1, nelem
         do_i2b: do ielem2 = 1, nelem
                                                ! new scattering table file ?
            if( sn%lquad(ielem2,ielem1) ) then
               isca(ielem2,ielem1) = 0
               cycle
            end if

            do_j1b: do jelem1 = 1, nelem
               do_j2b: do jelem2 = 1, nelem
                  if( ielem1==jelem1 .and. ielem2==jelem2 )  exit do_j1b
                  if( sn%scatfile(ielem2,ielem1) == sn%scatfile(jelem2,jelem1) ) then
                     isca(ielem2,ielem1) = isca(jelem2,jelem1)
                     cycle do_i2b
                  end if
               end do do_j2b
            end do do_j1b

            is = is + 1
            isca(ielem2,ielem1) = is

            call openFile( 10, sn%scatfile(ielem2,ielem1), 'formatted', 'read' )
!                ========  i        i                i          i
            read( 10,*,err=900) ne(is), np(is)
            read( 10,*,err=900) (etab(ie,is), ie=1,ne(is))
            read( 10,*,err=900) (ptab(ip,is), ip=1,np(is))
            read( 10,*,err=900) pcutoff(is)
            read( 10,*,err=900) ((a00(ie,ip,is), ie=1,ne(is)), ip=1,np(is))
            read( 10,*,err=900) ((a10(ie,ip,is), ie=1,ne(is)), ip=1,np(is))
            read( 10,*,err=900) ((a01(ie,ip,is), ie=1,ne(is)), ip=1,np(is))
            read( 10,*,err=900) ((b00(ie,ip,is), ie=1,ne(is)), ip=1,np(is))
            read( 10,*,err=900) ((b10(ie,ip,is), ie=1,ne(is)), ip=1,np(is))
            read( 10,*,err=900) ((b01(ie,ip,is), ie=1,ne(is)), ip=1,np(is))
            close (10)

         end do do_i2b
      end do do_i1b

      if( is /= nsca )  error stop 'cmangl: inconsistency in nsca'
!
!  Initialize the parameters for bicubic interpolation (aij, bij)
!
      do is = 1, nsca

         do ie = 1, ne(is)-1                       ! ie and ip must be increasing in these loops
            do ip = 1, np(is)-1

               f00 = a00(ie,ip,is)
               f01 = a00(ie,ip+1,is)
               f10 = a00(ie+1,ip,is)
               f11 = a00(ie+1,ip+1,is)
               fx00 = a10(ie,ip,is) * (etab(ie+1,is)-etab(ie,is))
               fx01 = a10(ie,ip+1,is) * (etab(ie+1,is)-etab(ie,is))
               fx10 = a10(ie+1,ip,is) * (etab(ie+1,is)-etab(ie,is))
               fx11 = a10(ie+1,ip+1,is) * (etab(ie+1,is)-etab(ie,is))
               fy00 = a01(ie,ip,is) * (ptab(ip+1,is)-ptab(ip,is))
               fy01 = a01(ie,ip+1,is) * (ptab(ip+1,is)-ptab(ip,is))
               fy10 = a01(ie+1,ip,is) * (ptab(ip+1,is)-ptab(ip,is))
               fy11 = a01(ie+1,ip+1,is) * (ptab(ip+1,is)-ptab(ip,is))

               a00(ie,ip,is) = f00
               a10(ie,ip,is) = fx00
               a01(ie,ip,is) = fy00

               a30(ie,ip,is) = fx10 - f10-f10 + f00+f00 + fx00 
               a20(ie,ip,is) = f10 - f00 - fx00 - a30(ie,ip,is)
               a03(ie,ip,is) = fy01 - f01-f01 + f00+f00 + fy00 
               a02(ie,ip,is) = f01 - f00 - fy00 - a03(ie,ip,is)

               s = f10 + f01 - f11 - f00
               a11(ie,ip,is) = fx01 + fy10 - fx00 - fy00 + s
               s = s + s
               a13(ie,ip,is) = fy11 - fy01 + fy10 - fy00 + s
               a12(ie,ip,is) = fx01 - fx00 - a11(ie,ip,is) - a13(ie,ip,is)
               a31(ie,ip,is) = fx11 - fx10 + fx01 - fx00 + s
               a21(ie,ip,is) = fy10 - fy00 - a11(ie,ip,is) - a31(ie,ip,is)

               eminln(is) = log(etab(1,is))

               f00 = b00(ie,ip,is)              ! now the same with the b's
               f01 = b00(ie,ip+1,is)
               f10 = b00(ie+1,ip,is)
               f11 = b00(ie+1,ip+1,is)
               fx00 = b10(ie,ip,is) * (etab(ie+1,is)-etab(ie,is))
               fx01 = b10(ie,ip+1,is) * (etab(ie+1,is)-etab(ie,is))
               fx10 = b10(ie+1,ip,is) * (etab(ie+1,is)-etab(ie,is))
               fx11 = b10(ie+1,ip+1,is) * (etab(ie+1,is)-etab(ie,is))
               fy00 = b01(ie,ip,is) * (ptab(ip+1,is)-ptab(ip,is))
               fy01 = b01(ie,ip+1,is) * (ptab(ip+1,is)-ptab(ip,is))
               fy10 = b01(ie+1,ip,is) * (ptab(ip+1,is)-ptab(ip,is))
               fy11 = b01(ie+1,ip+1,is) * (ptab(ip+1,is)-ptab(ip,is))

               b00(ie,ip,is) = f00
               b10(ie,ip,is) = fx00
               b01(ie,ip,is) = fy00

               b30(ie,ip,is) = fx10 - f10-f10 + f00+f00 + fx00 
               b20(ie,ip,is) = f10 - f00 - fx00 - b30(ie,ip,is)
               b03(ie,ip,is) = fy01 - f01-f01 + f00+f00 + fy00 
               b02(ie,ip,is) = f01 - f00 - fy00 - b03(ie,ip,is)

               s = f10 + f01 - f11 - f00
               b11(ie,ip,is) = fx01 + fy10 - fx00 - fy00 + s
               s = s + s
               b13(ie,ip,is) = fy11 - fy01 + fy10 - fy00 + s
               b12(ie,ip,is) = fx01 - fx00 - b11(ie,ip,is) - b13(ie,ip,is)
               b31(ie,ip,is) = fx11 - fx10 + fx01 - fx00 + s
               b21(ie,ip,is) = fy10 - fy00 - b11(ie,ip,is) - b31(ie,ip,is)

            end do
         end do

      end do

      return

 900  error stop 'Error during reading scattering table file'

   end subroutine initCmangl
!
!=======================================================================
!
   subroutine cmangl( cotthh, tau, e, p, ielem1, ielem2)
!             ======    o      o   i  i    i       i
!
!     Calculates the scattering angle in the center-of-mass system and the
!     time integral
!
!     -) uses bicubic interpolation in a table for cot(th/2)
!     -) the table is read from the scattering table file scatfile
!
!  Arguments:
!
      real(REAL8), intent(out) :: cotthh        ! cot(scattering angle/2)
      real(REAL8), intent(out) :: tau           ! time integral [rnorm]
      real(REAL8), intent(in) :: e              ! ion energy [enorm]
      real(REAL8), intent(in) :: p              ! impact parameter [rnorm]
      integer(INT4), intent(in) :: ielem1       ! index of projectile
      integer(INT4), intent(in) :: ielem2       ! index of recoil
!
!  Local Variables:
!
      integer(INT4) :: il, im, iu, ie, ip, is
      real(REAL8) :: eq, pq
!
!  Start of Executable Code:
!  ------------------------
!
!  Select table
!
      is = isca(ielem2,ielem1)
      if( is == 0 )  error stop 'cmangl: Missing table for specified elements'
!
!  no deflection if impact parameter larger than pcutoff
!
      if( p >= pcutoff(is) ) then
         cotthh = BIG
         tau = ZERO
         return
      end if
!  
!  binary search in the ptab-vector
!  
      il = 1
      iu = np(is)
 30   im = (il+iu) / 2
      if( p <= ptab(im,is) ) then
         iu = im
      else
         il = im
      end if
      if( iu-il > 1 )  goto 30
      ip = il
!  
!  reduced impact parameter
!  
      pq = (p - ptab(ip,is)) / (ptab(ip+1,is) - ptab(ip,is))
!  
!  binary search in the e-vector
!  
      il = 1
      iu = ne(is)
 40   im = (il+iu)/2
      if( e <= etab(im,is) ) then
         iu = im
      else
         il = im
      end if
      if( iu-il > 1 )  goto 40
      ie = il
!
!  reduced energy
!
      eq = (e - etab(ie,is)) / (etab(ie+1,is) - etab(ie,is))
!
!  bicubic interpolation of cot(th/2) in the table
!
      if( e < etab(1,is) ) then                 ! special treatment for e<etab(1,is)
         cotthh = ( (a03(1,ip,is)*pq+a02(1,ip,is))*pq + a01(1,ip,is) ) * pq + a00(1,ip,is)
         cotthh = cotthh * eminln(is) / log(e)
      else
         cotthh = (  ( (a13(ie,ip,is)*eq+a03(ie,ip,is))*pq +  &
                       (a12(ie,ip,is)*eq+a02(ie,ip,is)) ) * pq +  &
                     ( (a31(ie,ip,is)*eq+a21(ie,ip,is))*eq+a11(ie,ip,is) ) * eq +  &
                     a01(ie,ip,is)  )  *  pq    +  &
                  (  ( (a30(ie,ip,is)*eq+a20(ie,ip,is))*eq +  &
                       a10(ie,ip,is) ) * eq + a00(ie,ip,is)  )
      end if  
!
!  bicubic interpolation of tau in the table
!
      if( e < etab(1,is) ) then
         tau = ( (b03(1,ip,is)*pq+b02(1,ip,is))*pq + b01(1,ip,is) ) * pq + b00(1,ip,is)
      else
         tau = (  ( (b13(ie,ip,is)*eq+b03(ie,ip,is))*pq + &
                    (b12(ie,ip,is)*eq+b02(ie,ip,is)) ) * pq + &
                  ( (b31(ie,ip,is)*eq+b21(ie,ip,is))*eq+b11(ie,ip,is) ) * eq + &
                  b01(ie,ip,is)  )  *  pq    +  &
               (  ( (b30(ie,ip,is)*eq+b20(ie,ip,is))*eq + &
                    b10(ie,ip,is) ) * eq + b00(ie,ip,is)  )
      end if

   end subroutine cmangl
!
!=======================================================================
!
end module Sn_Cmangl_mod

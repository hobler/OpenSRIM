module Undefine_mod
!
!  Set and check undefined variables or arrays. In case of arrays a logical array is 
!  returned indicating the definition status of each array element.
!
!  Note: works by using an unlikely value to indicate an undefined variable. Therefore 
!        might fail if this value is identical to a desired defined value.
!
!  Used Modules:
!  ============
!
   use Const_mod, only: INT4, REAL8
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
   public :: undefine                        ! set variable undefined
!            ========
   interface undefine
      module procedure iundefine             ! undefine integer variable
      module procedure rundefine             ! undefine real variable
      module procedure cundefine             ! undefine character variable

      module procedure iaundefine            ! undefine 1-d integer array
      module procedure iaaundefine           ! undefine 2-d integer array
      module procedure raundefine            ! undefine 1-d real array
      module procedure raaundefine           ! undefine 2-d real array
      module procedure raaaundefine          ! undefine 3-d real array
      module procedure caundefine            ! undefine 1-d array of character variables
      module procedure caaundefine           ! undefine 2-d array of character variables
   end interface

   public :: defined                         ! check whether variable is undefined
!            =======
   interface defined
      module procedure idefined              ! check definition of integer variable
      module procedure rdefined              ! check definition of real variable
      module procedure cdefined              ! check definition of character variable

      module procedure iadefined             ! check definition of 1-d integer array
      module procedure iaadefined            ! check definition of 2-d integer array
      module procedure radefined             ! check definition of 1-d real array
      module procedure raadefined            ! check definition of 2-d real array
      module procedure raaadefined           ! check definition of 3-d real array
      module procedure cadefined             ! check definition of 1-d character array
      module procedure caadefined            ! check definition of 2-d character array
   end interface

   public :: definedlength                   ! get number of defined contiguous integer array
!            =============                   !    elements starting from index 1
   interface definedlength
      module procedure iadefinedlength       ! get number of defined contiguous integer array
                                             !    elements starting from index 1
      module procedure radefinedlength       ! get number of defined contiguous real array
                                             !    elements starting from index 1
   end interface
!
!  Private Data:
!  ============
!
   integer(INT4), parameter :: IUNDEF = -1367130542
   real(REAL8), parameter ::   RUNDEF = -8.97518E+11
   character(10), parameter :: CUNDEF = '#@/&^%>;~*'
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine iundefine( i )
!             =========  o

      integer(INT4), intent(out) :: i

      i = IUNDEF

   end subroutine iundefine
!
!-----------------------------------------------------------------------
!
   logical function idefined( i )
!                   ========  i

      integer(INT4), intent(in) :: i

      idefined = ( i /= IUNDEF )

   end function idefined
!
!=======================================================================
!
   subroutine rundefine( r )
!             =========  o

      real(REAL8), intent(out) :: r

      r = RUNDEF

   end subroutine rundefine
!
!-----------------------------------------------------------------------
!
   logical function rdefined( r )
!                   ========  i

      real(REAL8), intent(in) :: r

      rdefined = ( r /= RUNDEF )

   end function rdefined
!
!=======================================================================
!
   subroutine cundefine( c )
!             =========  o

      character(*), intent(out) :: c

      c = CUNDEF

   end subroutine cundefine
!
!-----------------------------------------------------------------------
!
   logical function cdefined( c )
!                   ========  i

      character(*), intent(in) :: c

      integer(INT4) :: lenmin

      lenmin = min( len(c), len(CUNDEF) )

      cdefined = ( c(1:lenmin) /= CUNDEF(1:lenmin) )

   end function cdefined
!
!=======================================================================
!
   subroutine iaundefine( i )
!             ==========  o

      integer(INT4), intent(out) :: i(:)
 
      i = IUNDEF

   end subroutine iaundefine
!
!-----------------------------------------------------------------------
!
   function iadefined( i ) result( idef )
!           =========  i            o

      integer(INT4), intent(in) :: i(:)
      logical :: idef(size( i ))

      idef = ( i /= IUNDEF )

   end function iadefined
!
!-----------------------------------------------------------------------
!
   function iadefinedlength( i ) result( ideflength )
!           ===============  i               o

      integer(INT4), intent(in) :: i(:)
      integer(INT4) :: ideflength

      do ideflength = 1, size( i )
         if( i(ideflength) == IUNDEF ) exit
      end do

      ideflength = ideflength - 1

   end function iadefinedlength
!
!=======================================================================
!
   subroutine raundefine( r )
!             ==========  o

      real(REAL8), intent(out) :: r(:)

      r = RUNDEF

   end subroutine raundefine
!
!-----------------------------------------------------------------------
!
   function radefined( r ) result( rdef )
!           =========  i            o

      real(REAL8), intent(in) :: r(:)
      logical :: rdef(size( r ))

      rdef = ( r /= RUNDEF )

   end function radefined
!
!-----------------------------------------------------------------------
!
   function radefinedlength( r ) result( rdeflength )
!           ===============  i               o

      real(REAL8), intent(in) :: r(:)
      integer(INT4) :: rdeflength

      do rdeflength = 1, size( r )
         if( r(rdeflength) == RUNDEF ) exit
      end do

      rdeflength = rdeflength - 1

   end function radefinedlength
!
!=======================================================================
!
   subroutine iaaundefine( i )
!             ===========  o

      integer(INT4), intent(out) :: i(:,:)

      i = IUNDEF

   end subroutine iaaundefine
!
!-----------------------------------------------------------------------
!
   function iaadefined( i ) result( idef )
!           ==========  i            o

      integer(INT4), intent(in) :: i(:,:)
      logical :: idef(size( i, 1 ),size( i, 2 ))

      idef = ( i /= IUNDEF )

   end function iaadefined
!
!=======================================================================
!
   subroutine raaundefine( r )
!             ===========  o

      real(REAL8), intent(out) :: r(:,:)

      r = RUNDEF

   end subroutine raaundefine
!
!-----------------------------------------------------------------------
!
   function raadefined( r ) result( rdef )
!           ==========  i            o

      real(REAL8), intent(in) :: r(:,:)
      logical :: rdef(size( r, 1 ),size( r, 2 ))

      rdef = ( r /= RUNDEF )

   end function raadefined
!
!=======================================================================
!
   subroutine raaaundefine( r )
!             ============  o

      real(REAL8), intent(out) :: r(:,:,:)

      r = RUNDEF

   end subroutine raaaundefine
!
!-----------------------------------------------------------------------
!
   function raaadefined( r ) result( rdef )
!           ===========  i            o

      real(REAL8), intent(in) :: r(:,:,:)
      logical :: rdef(size( r, 1 ),size( r, 2 ),size( r, 3 ))

      rdef = ( r /= RUNDEF )

   end function raaadefined
!
!=======================================================================
!
   subroutine caundefine( c )
!             ==========  o

      character(*), intent(out) :: c(:)

      c = CUNDEF

   end subroutine caundefine
!
!-----------------------------------------------------------------------
!
   function cadefined( c ) result( cdef )
!           =========  i            o

      character(*), intent(in) :: c(:)
      logical :: cdef(size( c ))

      integer(INT4) :: lenmin

      lenmin = min( len(c), len(CUNDEF) )

      cdef = ( c(:)(1:lenmin) /= CUNDEF(1:lenmin) )

   end function cadefined
!
!=======================================================================
!
   subroutine caaundefine( c )
!             ===========  o

      character(*), intent(out) :: c(:,:)

      c = CUNDEF

   end subroutine caaundefine
!
!-----------------------------------------------------------------------
!
   function caadefined( c ) result( cdef )
!           ==========  i            o

      character(*), intent(in) :: c(:,:)
      logical :: cdef(size( c, 1 ),size( c, 2 ))

      integer(INT4) :: lenmin

      lenmin = min( len(c), len(CUNDEF) )

      cdef = ( c(:,:)(1:lenmin) /= CUNDEF(1:lenmin) )

   end function caadefined
!
!=======================================================================
!
end module Undefine_mod


!*******************************************************************************
!*******************************************************************************
!      program test_defined
!
!      use Undefine_mod
!
!      implicit none
!
!      integer(INT4) i
!      real(REAL8) d
!      character*5  c5
!      character*15 c15
!      integer(INT4) ia(5)
!      real(REAL8) d1(4), d2(3,3), d3(2,2,2)
!
!**    start of executable code
!*     ************************
!
!      i=2
!      print *, 'i=',i,'defined=',defined(i)
!
!      call undefine(i)
!      print *, 'i=',i,'defined=',defined(i)
!
!      i=4
!      print *, 'i=',i,'defined=',defined(i)
!
!
!      d=2
!      print *, 'd=',d,'defined=',defined(d)
!
!      call undefine(d)
!      print *, 'd=',d,'defined=',defined(d)
!
!      d=4
!      print *, 'd=',d,'defined=',defined(d)
!
!
!      c5='12345'
!      print *, 'c5=',c5,'defined=',defined(c5)
!
!      call undefine(c5)
!      print *, 'c5=',c5,'defined=',defined(c5)
!
!      c5='54321'
!      print *, 'c5=',c5,'defined=',defined(c5)
!
!
!      c15='1234567890abcde'
!      print *, 'c15=',c15,'defined=',defined(c15)
!
!      call undefine(c15)
!      print *, 'c15=',c15,'defined=',defined(c15)
!
!      c15='edcba0987654321'
!      print *, 'c15=',c15,'defined=',defined(c15)
!
!
!      ia(1:5)=1
!      print *, 'ia=',ia,'defined=',defined(ia)
!
!      call undefine(ia)
!      print *, 'ia=',ia,'defined=',defined(ia)
!
!      ia(1:5)=5
!      print *, 'ia=',ia,'defined=',defined(ia)
!
!
!      d1(1:4)=1
!      print *, 'd1=',d1,'defined=',defined(d1)
!
!      call undefine(d1)
!      print *, 'd1=',d1,'defined=',defined(d1)
!
!      d1(1:4)=4
!      print *, 'd1=',d1,'defined=',defined(d1)
!
!
!      d2(1:3,1:3)=1
!      print *, 'd2=',d2,'defined=',defined(d2)
!
!      call undefine(d2)
!      print *, 'd2=',d2,'defined=',defined(d2)
!
!      d2(1:3,1:3)=3
!      print *, 'd2=',d2,'defined=',defined(d2)
!
!
!      d3(1:2,1:2,1:2)=1
!      print *, 'd3=',d3,'defined=',defined(d3)
!
!      call undefine(d3)
!      print *, 'd3=',d3,'defined=',defined(d3)
!
!      d3(1:2,1:2,1:2)=2
!      print *, 'd3=',d3,'defined=',defined(d3)
!
!      call undefine(d3(1,1,2))
!      print *, 'd3=',d3,'defined=',defined(d3)
!
!      call undefine(d3(:,1,1))
!      print *, 'd3=',d3,'defined=',defined(d3)
!
!      end program test_defined

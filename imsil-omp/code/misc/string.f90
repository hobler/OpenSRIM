module String_mod
!
!  Defines various operations on strings
!
!  Used Modules:
!  ============
!
   use Const_mod,    only: INT4, REAL8
   use Undefine_mod, only: defined
!
!  Default Type and Access
!  =======================
!
   implicit none
   private
   save
!
!  Public functions
!  ================
!
   public :: str                             ! returns a variable of (any) type as a string
!            ===                             !    with leading and trailing blanks removed
   interface str
      module procedure strc
      module procedure stri
      module procedure strl
      module procedure strr
   end interface str

   public :: adjustc                         ! center string
!            =======
   public :: to_lower                        ! convert string to lowercase
!            ========
   public :: to_upper                        ! convert string to uppercase
!            ========
   public :: capitalize                      ! convert first character to uppercase
!            ==========

contains
!
!=======================================================================
!
   function strc( var, fmt_ )
!           ====   i    i
!
!     Write string with given format to other string and remove blanks from both sides
!
!  Arguments:
!
      character(:), allocatable :: strc         ! result string
      character(*), intent(in) :: var           ! variable to be written to string
      character(*), intent(in), optional :: fmt_! format
 !
 !  Local variables:
 !
      character(256) :: string
 !
 !  Start of Executable Code:
 !  ------------------------
 !       
      if( present( fmt_ ) ) then
         write( string, '(' // fmt_ // ')' )  var
      else
         write( string, * )  var
      end if
        
      strc = trim( adjustl( string ) )

   end function strc
!
!=======================================================================
!
   function stri( var, fmt_ )
!           ====   i    i
!
!     Write an integer number with given format to other string and remove blanks from both 
!     sides
!
!  Arguments:
!
      character(:), allocatable :: stri         ! result string
      integer(INT4), intent(in) :: var          ! variable to be written to string
      character(*), intent(in), optional :: fmt_! format
 !
 !  Local variables:
 !
      character(256) :: string
 !
 !  Start of Executable Code:
 !  ------------------------
 !       
      if( present( fmt_ ) ) then
         write( string, '(' // fmt_ // ')' )  var
      else
         write( string, * )  var
      end if
        
      stri = trim( adjustl( string ) )

   end function stri
!
!=======================================================================
!
   function strr( var, fmt_ )
!           ====   i    i
!
!     Write a real number with given format to other string and remove blanks from both sides
!
!  Arguments:
!
      character(:), allocatable :: strr         ! result string
      real(REAL8), intent(in) :: var            ! variable to be written to string
      character(*), intent(in), optional :: fmt_! format
 !
 !  Local variables:
 !
      character(256) :: string
 !
 !  Start of Executable Code:
 !  ------------------------
 !       
      if( present( fmt_ ) ) then
         write( string, '(' // fmt_ // ')' )  var
      else
         write( string, * )  var
      end if
        
      strr = trim( adjustl( string ) )

   end function strr
!
!=======================================================================
!
   function strl( var, fmt_ )
!           ====   i    i
!
!     Write a logical with given format to other string and remove blanks from both sides
!
!  Arguments:
!
      character(:), allocatable :: strl         ! result string
      logical, intent(in) :: var                ! variable to be written to string
      character(*), intent(in), optional :: fmt_! format
 !
 !  Local variables:
 !
      character(256) :: string
 !
 !  Start of Executable Code:
 !  ------------------------
 !       
      if( present( fmt_ ) ) then
         write( string, '(' // fmt_ // ')' )  var
      else
         write( string, * )  var
      end if
        
      strl = trim( adjustl( string ) )
      
      if( strl == 'T' )  strl = 't'
      if( strl == 'F' )  strl = 'f'

   end function strl
!
!=======================================================================
!
   function adjustc( string )
!           =======    i
!
!     Center the string so that number of leading blanks equals or is one more than
!     number of trailing blanks
!
!     In addition, replace "undefined" string with blank
!
!     Note: Fortran provides the intrinsic procedures adjustl and adjustr
!
!  Arguments:
!
      character(:), allocatable :: adjustc   ! centered string
      character(*), intent(in) :: string     ! string to be centered
!
!  Local Variables:
!
      integer(INT4) :: indent
!
!  Start of Executable Code:
!  ------------------------
!
      adjustc = adjustl( string )

      if( defined( string ) ) then
         indent = (len(adjustc) - len_trim(adjustc) + 1) / 2
         adjustc(1+indent:) = adjustc
         adjustc(1:indent) = ' '
      else
         adjustc = ' '
      end if

   end function adjustc
!
!=======================================================================
!
   function to_lower( string )
!           ========    i
!
!     Convert a string to lowercase
!
!     Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
!     Original author: Clive Page
!
      character(:), allocatable :: to_lower
      character(*), intent(in) :: string
      integer :: i,j

      to_lower = string
      do i = 1, len(string)
          j = iachar(string(i:i))
          if (j>= iachar("A") .and. j<=iachar("Z") ) then
               to_lower(i:i) = achar(iachar(string(i:i))+32)
          else
               to_lower(i:i) = string(i:i)
          end if
      end do

   end function to_lower
!
!=======================================================================
!
   function to_upper( string )
!           ========    i
!
!     Convert a string to uppercase
!
!     Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
!     Original author: Clive Page
!
      character(:), allocatable :: to_upper
      character(*), intent(in) :: string
      integer :: i,j

      to_upper = string
      do i = 1, len(string)
          j = iachar(string(i:i))
          if (j>= iachar("a") .and. j<=iachar("z") ) then
               to_upper(i:i) = achar(iachar(string(i:i))-32)
          else
               to_upper(i:i) = string(i:i)
          end if
      end do

   end function to_upper
!
!=======================================================================
!
   function capitalize( string )
!           ==========    i
!
!     Capitalize a string, i.e. convert the first letter to uppercase
!
      character(:), allocatable :: capitalize
      character(*), intent(in) :: string

      capitalize = string
      capitalize(1:1) = to_upper(string(1:1))

   end function capitalize
!
!=======================================================================
!
end module String_mod

module Param_WriteTab_mod
!
!  Write an array as a table to the output file
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: INT4, REAL8
   use String_mod,      only: adjustc
   use Undefine_mod,    only: defined, undefine

   use Param_Setup_mod, only: nelem
   use Param_Elem_mod,  only: elem
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
   public :: writeTab                        ! write an array to the output file
!            ========
   interface writeTab
      module procedure cWriteTab2d           ! write table of character array
      module procedure rWriteTab2d           ! write table of real array
      module procedure rWriteTab1d           ! write table (list) of 1D real array
   end interface
!
!=======================================================================
!
contains
!
!=======================================================================
!
   subroutine cWriteTab2d( carr, title, ntriml_, indent_ )
!             ===========   i      i      i         i
!
!     Write a 2d character array as a table to the output file (unit 21)
!
!     This version assumes index1 = ion atom index, index2 = target atom index
!
!  Arguments:
!
      character(*), intent(in) :: carr(:,:)           ! array to be written to output file
      character(*), intent(in) :: title               ! title to be written above table
      integer(INT4), intent(in), optional :: ntriml_  ! number of characters at the beginning
                                                      !    of carr to be ignores
      logical, intent(in), optional :: indent_        ! flag indicating indentation of first
                                                      !    line
!
!  Local Variables:
!
      logical :: indent
      !character( max( maxval( len_trim(carr(:,:))-ntriml ), 10 ) ) :: bar, text
      character(:), allocatable :: bar, text
      character(:), allocatable :: carr_tmp(:,:)
      integer(INT4) :: i, ielem1, ielem2, nelem1, nelem2, ielemmin, ielemmax
      integer(INT4) :: ntriml, wrow, nrows
!
!  Start of Executable Code:
!  ------------------------
!
      if( present( ntriml_ ) ) then
         ntriml = ntriml_
      else
         ntriml = 0
      end if

      if( present( indent_ ) ) then
         indent = indent_
      else
         indent = .true.
      end if
!
!  No table if all entries are undefined
!
      if( .not.any( defined( carr(:,:) ) ) )  return
!
!  One or two lines only insted of table if all entries the same
!
      if( indent )  write( 21, '(t13)', advance='no' )

      if( all( carr(:,:) == carr(1,1) ) ) then
         if( len_trim(title) + len_trim( carr(1,1)(ntriml+1:) ) > 66 ) then
            write( 21, '(a)' )  trim( title ) // ':'
            write( 21, '(t16,a)' )  trim( carr(1,1)(ntriml+1:) )
         else
            write( 21, '(a)' )  trim( title ) // ': ' // trim( carr(1,1)(ntriml+1:) )
         end if
         write( 21, '()' )
         return
      end if
!
!  Otherwise write table...
!
      wrow = maxval( len_trim( carr(:,:) ) )
      wrow = wrow - ntriml                   ! subtract ntriml here instead of the line above for every entry
      wrow = max( wrow, 10 )                 ! 10 characters needed by header
      nrows = 66 / (wrow+1) - 1
      nrows = max( nrows, 1 )

      bar = repeat( '-', wrow )
      allocate( character(wrow) :: text )

      nelem2 = size( carr, 1 )
      nelem1 = size( carr, 2 )
      allocate( character(wrow) :: carr_tmp( nelem2, nelem1 ) )
      do ielem1 = 1, nelem1
         do ielem2 = 1, nelem2
            if( .not.defined(carr(ielem2,ielem1)) ) then
               carr_tmp(ielem2,ielem1) = '-'
            else
               carr_tmp(ielem2,ielem1) = carr(ielem2,ielem1)(ntriml+1:)
            end if
         end do
      end do
!
!  Title line
!
      write( 21, '(a)' )  trim(title) // ':'
!
!  Loop over subtables
!
      do i = 0, (nelem-1)/nrows
         ielemmin = nrows*i + 1
         ielemmax = min( nrows*i+nrows, nelem )
!
!  Header
!
         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) (bar // '+', ielem1=ielemmin,ielemmax)

         if( nelem1 == 1 ) then
            write( 21, '(t13,a)', advance='no' )  '|Target\   |'
            write( 21, '(t13,a)' )  bar // '|'
         else
            if( nelem2 == 1 ) then
               write( 21, '(t13,a)', advance='no' )  '|      \Ion|'
            else
               write( 21, '(t13,a)', advance='no' )  '|Target\Ion|'
            end if
            do ielem1 = ielemmin, ielemmax
               write( text, '(a,i2,a,a2,a)' )  'Atom', ielem1, '(', elem(ielem1)%name, ')'
               write( 21, '(a,a)', advance='no' ) adjustc( text ), '|'
            end do
            write( 21, '()' )
         end if

         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) (bar // '+', ielem1=ielemmin,ielemmax)
!
!  Table rows
!
         if( nelem2 == 1 ) then
            write( 21, '(t13,a,i2,a,a2,a)', advance='no' )  &
                  '|          |'
            do ielem1 = ielemmin, ielemmax
               write( 21, '(a,a)', advance='no' )  &
                     adjustc( carr_tmp(1,ielem1)(1:wrow) ), '|'
            end do
            write( 21, '()' )
         else
            do ielem2 = 1, nelem
               write( 21, '(t13,a,i2,a,a2,a)', advance='no' )  &
                     '|Atom', ielem2, '(', elem(ielem2)%name, ')|'
               do ielem1 = ielemmin, ielemmax
                  write( 21, '(a,a)', advance='no' )  &
                        adjustc( carr_tmp(ielem2,ielem1)(1:wrow) ), '|'
               end do
               write( 21, '()' )
            end do
         end if
!
!  Footer
!
         write( 21, '(t13,a)', advance='no' ) '+----------+'
         write( 21, '(5a)' ) (bar // '+', ielem1=ielemmin,ielemmax)
         write( 21, '()' )

      end do

      deallocate( carr_tmp )

   end subroutine cWriteTab2d
!
!=======================================================================
!
   subroutine rWriteTab2d( rarr, title, undeftext_, indent_ )
!             ===========   i      i        i          i
!
!     Write a 2d real array as a table to the output file (unit 21)
!
!     This version assumes index1 = ion atom index, index2 = target atom index
!
!  Arguments:
!
      real(REAL8), intent(in) :: rarr(:,:)            ! array to be written to output file
      character(*), intent(in) :: title               ! title to be written above table
      character(*), intent(in), optional :: undeftext_! text for undefined values
      logical, intent(in), optional :: indent_        ! flag indicating indentation of first
                                                      !    line
!
!  Local Variables:
!
      character(10) :: carr(size( rarr, 1 ),size( rarr, 2 ))
      integer(INT4) :: ielem1, ielem2
!
!  Start of Executable Code:
!  ------------------------
!
!  No table if all entries are undefined
!
      if( .not.any( defined( carr(:,:) ) ) )  return
!
!  Convert to character array, replacing undefined values by undeftext_
!
      do ielem1 = 1, size( carr, 2 )
         do ielem2 = 1, size( carr, 1 )
            if( defined( rarr(ielem2,ielem1) ) ) then
               write( carr(ielem2,ielem1), '(es10.3e2)' )  rarr(ielem2,ielem1)
            else
               if( present( undeftext_ ) ) then
                  carr(ielem2,ielem1) = undeftext_
               else
                  call undefine( carr(ielem2,ielem1) )
               end if
            end if
         end do
      end do
!
!  Write table
!
      call cWriteTab2d( carr, title )
!          ===========   i      i

   end subroutine rWriteTab2d
!
!=======================================================================
!
   subroutine rWriteTab1d( rarr, title, undeftext )
!             ===========   i      i        i
!
!     Write a 2d real array as a table to the output file (unit 21)
!
!     This version assumes index1 = ion atom index, index2 = target atom index
!
!  Arguments:
!
      real(REAL8), intent(in) :: rarr(:)     ! array to be written to output file
      character(*), intent(in) :: title      ! title to be written above table
      character(*), intent(in), optional :: undeftext ! text for undefined values
!
!  Local Variables:
!
      real(REAL8), allocatable :: rarr2d(:,:)
!
!  Start of Executable Code:
!  ------------------------
!
   allocate( rarr2d( 1, size( rarr ) ) )
   rarr2d(1,:) = rarr(:)

   call rWriteTab2d( rarr2d, title, undeftext )
!       ===========    i       i        i

   end subroutine rWriteTab1d
!
!=======================================================================
!
end module Param_WriteTab_mod

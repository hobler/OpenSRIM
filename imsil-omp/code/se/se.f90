module Se_mod
!
!  Handle electronic stopping
!
!  Used modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, HALF
   use Undefine_mod
   use Types_mod,          only: Proj_type

   use Param_Se_mod,       only: se

   use Param_Setup_mod,    only: nr, nelemIon, nelem
   use Param_Elem_mod,     only: elem
   use Param_Ion_mod,      only: ion
   use Param_Mat_mod,      only: mat
   use Param_Sn_mod,       only: sn

   use Se_Mixed_mod,       only: initSeMixed, evalSeMixed

   use Target_Surface_mod, only: getDistSurface
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
   public :: initSeModels           ! initialize stopping parameters
!            ============
   public :: evalSe                 ! evaluate electronic stopping
!            ====== 
   public :: tabSe                  ! calls evalSe repeatedly and writes a table
!            =====                  !   
!
!  Private Data:
!  ============
!
   logical :: init = .false.
!
!========================================================================
!
contains
!
!========================================================================
!
   subroutine initSeModels( )
!             ============
!
!     Initialize stopping power parameters
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
!     None
!
!  Start of Executable Code:
!  ------------------------
!
      if( all( se%semodel(:,:) == 'none' ) )  return

      init = .true.
!
!  Call initialization subroutines of specialized electronic stopping modules
!
      call initSeMixed( )                       ! mixed local/nonlocal model
!          ===========

   end subroutine initSeModels
!
!=======================================================================
!
   subroutine evalSe( dee, dedxnl, deeloc, proj, posc, ffp, pmax, ef, p_, ielem2_ )
!             ======   o     o       o      i     i    io    i    i   i    i
!
!     Evaluates the electronic energy loss
!
!  Arguments:
!
      real(REAL8), intent(out) :: dee                 ! electronic energy loss [eV]
      real(REAL8), intent(out) :: dedxnl              ! nonlocal energy loss per path length
                                                      !    [eV/A] (at end of ffp if iterated)
      real(REAL8), intent(out) :: deeloc              ! local energy loss [eV]
      type(Proj_type), intent(in) :: proj             ! projectile
      real(REAL8), intent(in) :: posc(:)              ! collision point (= where p is normal to
                                                      !    ffp)
      real(REAL8), intent(inout) :: ffp               ! free flight path [A] (reduced on exit if
                                                      !    projectile stops)
      real(REAL8), intent(in) :: pmax                 ! maximum impact param [A]
      real(REAL8), intent(in) :: ef                   ! trajectory cutoff energy [eV]
      real(REAL8), optional, intent(in) :: p_         ! impact parameters [A]
      integer(INT4), optional, intent(in) :: ielem2_  ! atom index of recoil
!
!  Local Variables
!
      real(REAL8) :: e                          ! projectile energy [eV]
      integer(INT4) :: ielem1                   ! index of projectile
      integer(INT4) :: id                       ! domain index

      logical :: inside_beg, inside_end
      real(REAL8) :: distsurfc, ffp1, ffp2, dee2
!
!  Start of Executable Code:
!  ------------------------
!
      e = proj%e
      ielem1 = proj%ielem
      id = proj%id

      if( all( se%semodel(:,ielem1) == 'none' ) ) then
         dee = ZERO
         dedxnl = ZERO
         deeloc = ZERO
         return
      end if

      if( proj%e < ef ) then
         dee = ZERO
         dedxnl = ZERO
         deeloc = ZERO
         ffp = ZERO
         return
      end if  

      if( .not.init )  goto 901
!
!  Take care that nonlocal electronic stopping is applied only inside the target. For that,
!  potentially divide the calculation in two parts and distinguish depending on whether the
!  start and endpoint of the free flight path are inside or outside the target (distsurf>/<0)
!
!  - inside -> inside   ... use evalSeMixed like normal
!  - outside -> outside ... use evalSeMixed with FFP=0
!  - inside -> outside  ... use evalSeMixed with FFP part inside the target and no local
!                           contribution; in a second step use evalSeMixed with FFP=0
!  - outside -> inside  ... use evalSeMixed with FFP part inside the target
!
!  Before the final call of evalSeMixed, ffp2 denotes the FFP for the final call of evalSeMixed,
!  and dee and ffp are set to values that yield the correct final values when the values
!  returned by the final call of evalSeMixed are added
!
      if( proj%distsurf > ffp ) then
         inside_beg = .true.
         inside_end = .true.
      else
         inside_beg = ( proj%distsurf > ZERO )
         call getDistSurface( posc(:), distsurfc )
!             ==============    i          o
         inside_end = ( distsurfc > ZERO )
      end if

      if( inside_beg ) then
         if( inside_end ) then                  ! inside -> inside
            ffp2 = ffp
            ffp = ZERO
            dee = ZERO
         else                                   ! inside -> outside
            ffp1 = ffp * proj%distsurf / (proj%distsurf - distsurfc) ! ffp inside target
            ffp2 = ZERO                         ! no nonlocal Se outside target
            call evalSeMixed( dee, dedxnl, deeloc, e, ielem1, id, ffp1, pmax, ef )
!                ===========   o     o       o     i    i     i    io    i    i
            e = e - dee
            if( e < ef ) then
               dedxnl = ZERO
               deeloc = ZERO
               ffp = ffp1
               return
            end if
         end if
      else
         if( inside_end ) then                  ! outside - > inside
            ffp2 = ffp * distsurfc / (distsurfc - proj%distsurf)
            ffp = ffp - ffp2
            dee = ZERO
         else                                   ! outside -> outside
            ffp2 = ZERO
            dee = ZERO
         end if
      end if
!
!  Mixed nonlocal and impact parameter dependent models
!
      call evalSeMixed( dee2, dedxnl, deeloc, e, ielem1, id, ffp2, pmax, ef, p_, ielem2_ )
!          ===========    o     o       o     i    i     i    io    i    i   i     i

      dee = dee + dee2
      if( ffp2 == ZERO )  dedxnl = ZERO         ! ffp2 was used to force nonlocal Se to zero
      ffp = ffp + ffp2

      return

 901  error stop 'evalSe: not initialized'

   end subroutine evalSe
!
!=======================================================================
!
   subroutine tabSe( )
!             =====
!
!     Makes a table of the nonlocal electronic stopping power as a function
!     of energy between e=ef and e=e0 and writes it to file SE
!
!  Arguments:
!
!     None
!
!  Local Variables:
!
      integer(INT4) :: ielem1, ielem2, ir, i
      real(REAL8) :: ffp, p, e, pmax
      real(REAL8) :: dee, dedxnl, deeloc
!
!  Start of Executable Code:
!  ------------------------
!
      if( .not.se%ltab ) return

      if( .not.init )  goto 901

      open( 122, file='SE', status='replace', form='formatted' )

      ffp = 1

      p = se%ptab

      do ir = 1, nr
         write( 122,'(a,a,a,a)' ) &
               '* Ion: ', trim(ion%name),', Target: ', trim(mat(ir)%name)
         do ielem1 = 1, nelemIon
            do ielem2 = 1, nelem
               write( 122,'(a,a,a,a)' ) &
                     '* Ion Atom: ', elem(ielem1)%name,', Target Atom: ', elem(ielem2)%name
               write( 122,*) 3,101
               write( 122,*) 0,0,0
               do i = 0, 100
                  e = sn%ef * (ion%e / sn%ef) ** (i/100d0)
                  if( mat(ir)%lxtal ) then
                     pmax = sn%pmax
                  else
                     pmax = sn%pmaxpsi
                  end if
                  call evalSeMixed( dee, dedxnl, deeloc, e, ielem1, ir, ffp, pmax, ZERO, p,  &
!                      ===========   o     o       o     i    i     i   io    i     i    i
                                    ielem2 )
!                                     i
                  write( 122,*) e, dedxnl, deeloc
               end do
            end do
         end do
      end do

      close (122)

      return

 901  error stop 'tabSe: not initialized'

   end subroutine tabSe
!
!=======================================================================
!
end module Se_mod

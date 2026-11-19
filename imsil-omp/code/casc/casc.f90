module Casc_mod
!
!  Simulate the cascade(s) originating from one impinging ion.
!
!  In case of molecular ions the (overlapping) cascades of all ion atoms are simulated
!
!  Used Modules:
!  ============
!
   use Const_mod,          only: INT4, REAL8,  &
                                 ZERO, STATE
   use Array_mod,          only: check_allocation
   use Types_mod,          only: Proj_type
   use String_mod,         only: str

   use Random_mod,         only: rnnset

   use Param_Setup_mod,    only: nelemIon
   use Param_Ion_mod,      only: ion
   use Param_Mat_mod,      only: irxtal

   use Sn_RestFFP_mod,     only: resetRestFFP
   use Xtal_OldSites_mod,  only: initOldSites, updateOldSites

   use Ion_Proj_mod,       only: getProjectile
   use Casc_Branch_mod,    only: branch
   use Casc_Step_mod,      only: step

   use Eval_mod,           only: eval11, eval12
   !debug
   !   use Eval_Stat_mod,         only: sumdee, sumedam, sumeb, sumet
   !end debug
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
   public :: casc                            ! simulate the cascade(s)
!
contains
!
!=======================================================================
!
   subroutine casc( iion, done )
!             ====   i     o
!
!     Simulate the cascade(s) originating from one impinging ion.
!
!  Arguments:
!
      integer(INT4), intent(in) :: iion         ! index of ion
      logical, intent(out) :: done              ! flag indicating last ion has been processed
!
!  Local Variables:
!
      logical :: newmolecule, found
      character(256) :: errmsg
      integer(INT4) :: ielem, ielem2, iatom, iatom2
      integer(INT4) :: ig, istat
      type(Proj_type), allocatable :: proj(:)   ! recoil generation stack

      integer(INT4) :: icoll                    ! collision counter (for debugging)
!     debug
!        real(REAL8) :: sume
!        real(REAL8) :: sumeion = ZERO
!     end debug
!
!  Start of Executable Code:
!  ------------------------
!
      if( ion%usetra ) then
         ielem2 = 1
      else
         ielem2 = nelemIon
      end if

      Ion_Elements_loop: do ielem = 1, ielem2
!     ---------------------------------------

         if( ion%usetra ) then
            iatom2 = 1
         else
            iatom2 = ion%natom(ielem)
         end if

         Element_Atoms_loop: do iatom = 1, iatom2
!        ----------------------------------------
!
!  Get projectile and do some initializations needed for each atom of the ion
!
            ig = 1
            call rnnset( )
!                ======

            if( allocated( proj ) )  deallocate( proj )
            allocate( proj( 1 ), stat=istat, errmsg=errmsg )
            call check_allocation( 'main', 'proj', istat, errmsg )

            newmolecule = (ielem == 1 .and. iatom == 1)
            call getProjectile( ielem, newmolecule, proj(1), found, done )
!                =============    i         i          o       o     o
            if( done )  return
            if( .not.found )  cycle Element_Atoms_loop

            proj(1)%iion = iion

            if( irxtal /= 0 ) then
               call initOldSites( )
!                   ============
               if( proj(1)%origin%lattice )  call updateOldSites( 1, proj_=proj(1) )
!                                                 ==============  i       i
            end if

            call resetRestFFP( 1 )           ! reset old ffp list
!                ============  i

!              debug
!                  sumeion = sumeion + proj(1)%e
!              end debug
!
! Initialization that are done only for each new molecule
!
            if( newmolecule )  call eval11( )
!                                   ======
!
! Initializations done for each atom of the ion
!
            call eval12( proj(1) )
!                ======    i
!
! Check if ion trajectory is terminated immediately (indicated by ig=0)
!
            call branch( proj, ig, ZERO )  ! results in ig=0 if trajectory terminates
!                ======   io   io   i
!
! Brute force implementation of molecular dynamics option
!
!            if( sn%lmd ) then
!               call moldyn (proj(1)%pos,proj(1)%dir,proj(1)%e,proj(1)%ielem)
!!                   ======   io       io      io    i
!               cycle
!            end if

            icoll = 0                        ! for debugging

            Collision_loop: do
!           ------------------

               !debug
                  icoll = icoll + 1
                  !print *,  'Thread=', omp_get_thread_num(), ' icoll=', icoll
                  !write( 998, *)  'icoll=', icoll
                  !if( proj(1)%emin < 1.03 ) then
                  !   print *, icoll, proj(1)%e, proj(1)%pos(3)
                  !end if
               ! end debug
!
!     Check if end of cascade is reached
!
               if( ig == 0 )  exit Collision_loop
!
!     Perform one step in the cascade development
!
               call step( proj, ig )
!                   ====   io   io

               !debug
!                     sume = sum( proj(1:ig)%e, mask=.not.proj(1:ig)%virtual .or.  &
!                                                    (.not.proj(1:ig)%virtualparent .and.  &
!                                                     proj(1:ig)%state==STATE%SUBTHRESHOLD) ) - &
!                            sum( proj(1:ig)%ededuct ) +  &
!                            sum( sumdee(:)/ion%ddose ) + sum( sumedam(:)/ion%ddose ) +  &
!                            sum( sumeb(:)/ion%ddose ) + sum( sumet(:)/ion%ddose )
                  !if( abs( sume - sumeion ) > EPSMB * max( sume, sumeion ) ) then
                  !   write( 998, * ) '>>>>> ', str(icoll),  &
                  !         ' ', str( sume, 'f99.3' ),  &
                  !         ' ', str( sum( proj(1:ig)%e, mask=.not.proj(1:ig)%virtual .or.  &
                  !                               (.not.proj(1:ig)%virtualparent .and.  &
                  !                                proj(1:ig)%state==STATE%SUBTHRESHOLD) ), &
                  !                   'f99.3' ),  &
                  !         ' ', str( sum( proj(1:ig)%ededuct ), 'f99.3' ),  &
                  !         ' ', str( sum( sumdee(:)/ion%ddose ), 'f99.3' ),  &
                  !         ' ', str( sum( sumedam(:)/ion%ddose ), 'f99.3' ),  &
                  !         ' ', str( sum( sumeb(:)/ion%ddose ), 'f99.3' ),  &
                  !         ' ', str( sum( sumet(:)/ion%ddose ), 'f99.3' )
                  !end if
               !end debug

            end do Collision_loop
!           ---------------------

         end do Element_Atoms_loop
!        -------------------------

      end do Ion_Elements_loop
!     ------------------------

   end subroutine casc
!
!=======================================================================
!
end module Casc_mod

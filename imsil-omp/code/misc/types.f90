module Types_mod
!
!  Define the projectiles data type
!
!  Used Modules:
!  ============
!
   use Const_mod,       only: REAL8, INT4, ZERO
!
!  Default Type and Access:
!  =======================
!
   implicit none
   public
   save
!
!  Public Data Types:
!  =================
!
   type :: Ipos_type                         ! lattice site index type
!          =========
      integer(INT4) :: ipos(3)               ! coordinates in units of quarter lattice
                                             !    constant (only for old diamond lattice)
      integer(INT4) :: iuc(3)                ! coordinates of unit cell in units of lattice
                                             !    constant
      integer(INT4) :: iucs                  ! index of lattice site within the unit cell
   end type

   type Origin_type                          ! Projectile origin
!       ===========
      logical :: internal                    ! indicates whether projectile originates from
                                             !    inside the target
      logical :: amorphous                   ! indicates whether projectile originates from
                                             !    an amorphous region
      logical :: lattice                     ! indicates whether projectile originates from
                                             !    a lattice site
      type(Ipos_type) :: ipos                ! original position in lattice coordinates
                                             !    (only set if lattice=t)
      integer(INT4) :: isa                   ! site index of lattice atom
      real(REAL8) :: e                       ! original energy
      real(REAL8) :: pos(3)                  ! original position
      real(REAL8) :: dir(3)                  ! original direction
   end type

   type :: Proj_type                         ! Projectile data type
!          =========
      real(REAL8) :: e                       ! energy [eV]
      real(REAL8) :: pos(3)                  ! position [A]
      real(REAL8) :: dir(3)                  ! direction vector (|dir|=1)
      integer(INT4) :: ielem                 ! chemical element index of the projectile
      integer(INT4) :: id                    ! index of domain pos(:) is in or closest to
      integer(INT4) :: ig                    ! recoil generation (index in projectile stack)
      integer(INT4) :: iion                  ! index of ion
      integer(INT4) :: nreflect              ! number of surface reflections
      logical :: reentered                   ! flag indicating that projectile or one of its
                                             !    parents has reentered the target

      type(Origin_type) :: origin            ! origin of projectile
      real(REAL8) :: posold(3)               ! previous position of projectile
      real(REAL8) :: ffpmin                  ! minimum free flight path
      real(REAL8) :: distsurf                ! distance from the surface
      real(REAL8) :: emin                    ! minimum kinetic energy during collision [eV]
      real(REAL8) :: ededuct                 ! energy of a virtual projectile to be deducted
                                             !    from its parent (due to sputtering)
      integer(INT4) :: state                 ! one of the STATEs defined in Const_mod
      logical :: virtual = .false.           ! flag indicating projectile finally exists
                                             !    only if it leaves the target
      logical :: virtualparent = .false.     ! virtuality of projectile that has recoiled
                                             !    current projectile
   end type

   integer(INT4), parameter :: NATOMM = 4    ! maximum number of atoms per lattice site

   type :: Site_type                         ! Lattice site type
!          =========
      type(Ipos_type) :: ipos                ! position (crystal coordinate system)
      logical :: latticeAtom = .true.        ! lattice atom flag
      integer(INT4) :: ielem                 ! index of atom species normally on this site
      integer(INT4) :: nsa = 1               ! number of atoms at lattice site
      integer(INT4) :: iatom(NATOMM)         ! (nsa) index of atom in atom(:)
      logical :: flag                        ! auxiliary flag
   end type

   type :: Atom_type                         ! Lattice site atom type
!          =========
      real(REAL8) :: pos(3)                  ! position (target coordinate system)
      integer(INT4) :: ielem                 ! atom species index
      integer(INT4) :: isite                 ! index of lattice site in site(:)
      integer(INT4) :: isa                   ! index among atoms belonging to site(isite)
      logical :: lvib                        ! flag indicating whether pos(:) contains
   end type                                  !   lattice vibrations

   type :: LSDefect_type                     ! Lattice site defect type (to be modified)
!          =============
      integer(INT4) :: ipos(3) = 0           ! lattice site position in the column (crystal
                                             !    coordinate system)
      integer(INT4) :: nsa  = 1              ! number of atoms at lattice site
      integer(INT4) :: ielem(NATOMM) = 0     ! (nsa) atom species index     
      real(REAL8) :: dpos(3,NATOMM) = ZERO   ! (nsa) displacement vector of atom wrt lattice
                                             !    site (crystal coordinate system)
      integer(INT4) :: ior = 0               ! defect orientation
      integer(INT4) :: type = 0              ! type of defect
   end type
!
!=======================================================================
!
end module Types_mod

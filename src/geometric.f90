module GEOMETRIC

  use PARAMETERS

  implicit none

  private

  public TIBPoint, TIBPoint_MomentumSource, TIBPoint_ScalFlSource, TIBPoint_Viscosity, &
         UIBPoints, VIBPoints, WIBPoints, ScalFlIBPoints, &
         InitSolidBodies,GetSolidBodiesBC, &
         obstaclefile
#ifdef CUSTOMSB
  public TLine, TPlane, TPolyhedron, TBall, TCylJacket, TCylinder, TTerrainPoint, TTerrain, TSolidBody, &
         NoneBody, Polyhedron, Ball, Cylinder, Terrain, &
         AddSolidBody
#endif

  type TLine                            !These object could be implemented using Fortran's 2003 inheritance.
    real(KND) xc,yc,zc                  !This approach using pointers is more portable and less safe.
    real(KND) a,b,c
 !  contains
 !    procedure Nearest => TLine_Nearest
  end type TLine


  type TPlane
    real(KND) a,b,c,d      !ax+by+cz+d/=0 for inner half-space
    logical gl             !T > in ineq. above F < in ineq. above
    logical :: rough = .false.!T rough surface, F flat surface
    real(KND) z0           !roughness parameter
 !  contains
 !    procedure Inside => TPlane_Inside
 !    procedure Nearest => TPlane_Nearest
  end type TPlane


  type TPolyhedron
    integer :: nplanes = 0
    type(TPlane),dimension(:),allocatable :: Planes !intersection of half-spaces
 !  contains
 !    procedure Inside => TPolyhedron_Inside
 !    procedure Nearest => TPolyhedron_Nearest
 !    procedure NearestOut => TPolyhedron_NearestOut
  end type TPolyhedron


  type TBall
    real(KND) xc,yc,zc,r
    logical :: rough = .false. !T rough surface, F flat surface
    real(KND) z0            !roughness parameter
 !  contains
 !    procedure Inside => TBall_Inside
 !    procedure Nearest => TBall_Nearest
 !    procedure NearestOut => TBall_NearestOut
  end type TBall


  type TCylJacket
    real(KND) xc,yc,zc
    real(KND) a,b,c
    real(KND) r
    logical :: rough = .false. !T rough surface, F flat surface
    real(KND) z0            !roughness parameter
 !  contains
 !    procedure Inside => TCylJacket_Inside
 !    procedure Nearest => TCylJacket_Nearest
  end type TCylJacket


  type TCylinder
    type(TCylJacket) Jacket
    type(TPlane),pointer :: Plane1 => null() ,Plane2 => null()
 !  contains
 !    procedure Inside => TCylinder_Inside
 !    procedure Nearest => TCylinder_Nearest
 !    procedure NearestOut => TCylinder_NearestOut
  end type TCylinder


  type TTerrainPoint
    real(KND) :: elev = 0
    logical :: rough = .false.
    real(KND) z0
  end type TTerrainPoint


   type TTerrain
     type(TTerrainPoint),dimension(:,:),allocatable :: UPoints,VPoints,PrPoints !allocate with a buffer of width 1 (i.e. 0:Xnx)
 !  contains
 !    procedure Inside => TTerrain_Inside
 !    procedure Nearest => TTerrain_Nearest
 !    procedure NearestOut => TTerrain_NearestOut
   end type TTerrain


  type TSolidBody
    integer numofbody
    integer :: typeofbody = 0                              !0.. none, 1..polyhedron, 2.. ball, 3.. cylinder, 4.. terrain
    type(TPolyhedron),pointer :: Polyhedron => null()    !asociated will be only part writen in typeofbody
    type(TBall),pointer :: Ball => null()
    type(TCylinder),pointer :: Cylinder => null()
    type(TTerrain),pointer :: Terrain => null()
    logical   :: rough = .false.                             !T rough surface, F flat surface
    real(KND) :: z0 = 0                                      !roughness parameter
    real(KND) :: temperatureflux = 0
    type(TSolidBody),pointer :: next =>null()
 !  contains
 !    procedure Inside => TSolidbody_Inside
 !    procedure Nearest => TSolidbody_Nearest
 !    procedure NearestOut => TSolidbody_NearestOut
  end type TSolidbody


  type(TSolidBody),pointer :: FirstSB => null()         !First member of the linked list of solid objects


  integer, parameter :: NoneBody = 0, Polyhedron = 1, Ball = 2, Cylinder = 3, Terrain = 4

  interface Inside
    module procedure TPlane_Inside
    module procedure TPolyhedron_Inside
    module procedure TBall_Inside
    module procedure TCylJacket_Inside
    module procedure TCylinder_Inside
    module procedure TTerrain_Inside
    module procedure TSolidBody_Inside
  end interface Inside

  interface Nearest   !shadows intrinsic nearest(), use intrinsic statement if needed
    module procedure TLine_Nearest
    module procedure TPlane_Nearest
    module procedure TPolyhedron_Nearest
    module procedure TBall_Nearest
    module procedure TCylJacket_Nearest
    module procedure TCylinder_Nearest
    module procedure TTerrain_Nearest
    module procedure TSolidBody_Nearest
  end interface Nearest

  interface NearestOut  !from outside
    module procedure TPolyhedron_NearestOut
    module procedure TBall_NearestOut
    module procedure TCylinder_NearestOut
    module procedure TTerrain_NearestOut
    module procedure TSolidBody_NearestOut
  end interface NearestOut



  type TInterpolationPoint
    integer   :: xi                !coordinates of the interpolation points
    integer   :: yj
    integer   :: zk
    real(KND) :: coef              !interpolation coefficients for the interpolation points
  endtype TInterpolationPoint

  type TVelIBPoint    !NOTE: this implementation is not very efficient, because it computes the int. coefficients again and again. It should be changed similarly to the scalar variant below.
    integer   :: component      !1..U, 2..V, 3..W
    integer   :: xi             !coordinates of the grid point
    integer   :: yj
    integer   :: zk
    real(KND) :: distx          !vector to the nearest boundary point
    real(KND) :: disty
    real(KND) :: distz
    integer   :: dirx           !integer form of the above vector (~sign(distx))
    integer   :: diry
    integer   :: dirz
    integer   :: interp         !kind of interpolation 0.. none (boundarypoint), 1..linear, 2..bilinear, 3..trilinear
    integer   :: interpdir      !direction of interpolation in the linear case, for bilinear it is normal direction to the interpolation plane
    type(TInterpolationPoint),dimension(:),allocatable :: IntPoints !array of interpolation points
    type(TVelIBPoint),pointer :: next => null()  !pointer to the next item in the list, not used when in the array in any way
 !  contains
 !    procedure Create         => TVelIBPoint_Create
 !    procedure AddToList      => TVelIBPoint_AddToList
 !    procedure DeallocateList => TVelIBPoint_DeallocateList
  end type TVelIBPoint



  type TScalFlIBPoint
    integer                :: xi                       !coordinates of the grid point
    integer                :: yj
    integer                :: zk
    real(KND)              :: dist                     !distance to the boundary
    type(TInterpolationPoint),dimension(:),allocatable :: IntPoints !array of interpolation points
    integer                :: interp                   !kind of interpolation 1.. none (1 point outside), 2..linear, 4..bilinear  other values not allowed
    real(KND)              :: temperatureflux = 0      !desired temperature flux
    type(TScalFlIBPoint),pointer :: next => null()     !pointer to the next item in the list, not used when in the array in any way
 !  contains
 !    procedure Create         => TScalFlIBPoint_Create
 !    procedure AddToList      => TScalFlIBPoint_AddToList
 !    procedure DeallocateList => TScalFlIBPoint_DeallocateList
  end type TScalFlIBPoint


  interface Create
    module procedure TVelIBPoint_Create
    module procedure TScalFlIBPoint_Create
  end interface Create

  interface AddToList
    module procedure TVelIBPoint_AddToList
    module procedure TScalFlIBPoint_AddToList
  end interface AddToList

  interface DeallocateList
    module procedure TVelIBPoint_DeallocateList
    module procedure TScalFlIBPoint_DeallocateList
  end interface DeallocateList



  type TIBPoint
    integer   :: xi
    integer   :: yj
    integer   :: zk
    real(KND) :: dist
    real(KND) :: temperatureflux = 0
    integer   :: interp
    type(TInterpolationPoint),dimension(:),allocatable :: IntPoints !array of interpolation points
 !   contains
 !     procedure Interpolate      => TIBPoint_Interpolate
 !     procedure InterpolateTDiff => TIBPoint_Interpolate_TDiff
 !     procedure ScalFlSource     => TIBPoint_ScalFlSource
 !     procedure MomentumSource   => TIBPoint_MomentumSource
 !     procedure Viscosity        => TIBPoint_Viscosity
  end type TIBPoint

  type(TVelIBPoint),pointer,save    :: FirstIBPoint => null(), LastIBPoint => null()
  type(TScalFlIBPoint),pointer,save :: FirstScalFlIBPoint => null(), LastScalFlIBPoint => null()

  type(TIBPoint),dimension(:),allocatable,save :: UIBPoints, VIBPoints, WIBPoints
  type(TIBPoint),dimension(:),allocatable,save :: ScalFlIBPoints

  integer, save :: NUIBPoints = 0 ,NVIBPoints = 0, NWIBPoints = 0, NScalFlIBPoints = 0

  character(80) :: obstaclefile = ''

  interface assignment (=)
    module procedure VelIBPtoIBP
  end interface

  interface assignment (=)
    module procedure ScalFlIBPtoIBP
  end interface











contains












  subroutine VelIBPtoIBP(IBP,VelIBP)
    type(TIBPoint),intent(out)    :: IBP
    type(TVelIBPoint),intent(in)  :: VelIBP

    IBP%xi = VelIBP%xi
    IBP%yj = VelIBP%yj
    IBP%zk = VelIBP%zk
    IBP%interp = VelIBP%interp

    allocate(IBP%IntPoints(size(VelIBP%IntPoints)))

    IBP%IntPoints = VelIBP%IntPoints
  end subroutine VelIBPtoIBP

  subroutine ScalFlIBPtoIBP(IBP,ScalFlIBP)
    type(TIBPoint),intent(out)    :: IBP
    type(TScalFlIBPoint),intent(in)  :: ScalFlIBP

    IBP%xi = ScalFlIBP%xi
    IBP%yj = ScalFlIBP%yj
    IBP%zk = ScalFlIBP%zk
    IBP%dist = ScalFlIBP%dist
    IBP%temperatureflux = ScalFlIBP%temperatureflux
    IBP%interp = ScalFlIBP%interp

    allocate(IBP%IntPoints(size(ScalFlIBP%IntPoints)))

    IBP%IntPoints = ScalFlIBP%IntPoints
  end subroutine ScalFlIBPtoIBP











  pure function TIBPoint_Interpolate(IBP,U,lb) result(Uint)
    real(KND) :: Uint
    type(TIBPoint),intent(in) :: IBP
    integer,intent(in) :: lb
    real(KND),dimension(lb:,lb:,lb:),intent(in) :: U
    integer i

    Uint = 0

    do i=1,IBP%interp
      Uint = Uint + IBP%IntPoints(i)%coef * U(IBP%IntPoints(i)%xi,&
                                              IBP%IntPoints(i)%yj,&
                                              IBP%IntPoints(i)%zk)
    enddo
  end function TIBPoint_Interpolate


  pure function TIBPoint_InterpolateTDiff(IBP,U) result(Uint)
    real(KND) :: Uint
    type(TIBPoint),intent(in) :: IBP
    real(KND),dimension(-1:,-1:,-1:),intent(in) :: U
    integer i,n

    n = 0
    Uint = 0

    do i=1,IBP%interp
      if (abs(IBP%IntPoints(i)%coef-0._KND)>epsilon(1._KND)) then
        Uint = Uint + U(IBP%IntPoints(i)%xi,&
                        IBP%IntPoints(i)%yj,&
                        IBP%IntPoints(i)%zk)
        n = n + 1
      endif
    enddo

    Uint = Uint + U(IBP%xi,IBP%yj,IBP%zk)
    Uint = Uint / (n+1)
  end function TIBPoint_InterpolateTDiff


  pure function TIBPoint_ScalFlSource(IBP,Scalar,sctype)  result(src)  !Virtual scalar source for the Immersed Boundary Method with prescribed scalar flux on the boundary
    real(KND) :: src

    type(TIBPoint),intent(in) :: IBP
    real(KND),intent(in)      :: Scalar(-1:,-1:,-1:)
    integer,intent(in)        :: sctype

    real(KND) intscal,intTDiff

    intscal = TIBPoint_Interpolate(IBP,Scalar,-1)

    if (sctype==1) then
      intTDiff = TIBPoint_InterpolateTDiff(IBP,TDiff)

      if (intTDiff>0)  intscal = intscal + IBP%temperatureflux * IBP%dist / intTDiff
    endif

    src = (intscal - Scalar(IBP%xi,IBP%yj,IBP%zk)) / dt

  end function TIBPoint_ScalFlSource



  pure function TIBPoint_Viscosity(IBP,Viscosity)  result(src)  !Virtual scalar source for the Immersed Boundary Method with prescribed scalar flux on the boundary
    real(KND) :: src

    type(TIBPoint),intent(in) :: IBP
    real(KND),intent(in)      :: Viscosity(-1:,-1:,-1:)

    src = TIBPoint_Interpolate(IBP,Viscosity,-1)

  end function TIBPoint_Viscosity



  pure function TIBPoint_MomentumSource(IBP,U) result(src)
    real(KND) :: src

    type(TIBpoint),intent(in) :: IBP
    real(KND),intent(in)                      :: U(-2:,-2:,-2:)

    src = (TIBPoint_Interpolate(IBP,U,-2) - U(IBP%xi,IBP%yj,IBP%zk)) / dt

  end function TIBPoint_MomentumSource





  subroutine TVelIBPoint_Create(IBP,xi,yj,zk,xU,yU,zU,Utype,component)
    type(TVelIBPoint),intent(out)               :: IBP
    integer,intent(in)                          :: xi,yj,zk
    real(KND),dimension(-2:),intent(in)         :: xU,yU,zU
    integer,dimension(-2:,-2:,-2:),intent(in)   :: Utype
    integer,intent(in)                          :: component

    type(TSolidBody),pointer :: SB
    integer dirx,diry,dirz,n1,n2,nx,ny,nz
    real(KND) x,y,z,xnear,ynear,znear,t
    logical free100,free010,free001

    x = xU(xi)                                !real coordinates of the IB forcing point
    y = yU(yj)
    z = zU(zk)
    call SetCurrentSB(SB,Utype(xi,yj,zk))
    call NearestOut(SB,xnear,ynear,znear,x,y,z)
    IBP%component = component
    IBP%xi = xi                                !integer grid coordinates
    IBP%yj = yj
    IBP%zk = zk
    IBP%distx = xnear-x                       !real distance to the boundary in the x,y,z direction
    IBP%disty = ynear-y
    IBP%distz = znear-z
    IBP%dirx = nint(sign(1.0_KND,IBP%distx))  !integer value denoting direction to the boundary
    IBP%diry = nint(sign(1.0_KND,IBP%disty))
    IBP%dirz = nint(sign(1.0_KND,IBP%distz))

    dirx = abs(IBP%dirx)                      !local temporary variable with abs(dir)
    diry = abs(IBP%diry)
    dirz = abs(IBP%dirz)

    if (.not.Inside(SB,xU(xi),yPr(yj),zPr(zk)))  then  !For now, if actually outside the body, set an artificial boundary
      IBP%interp = 0                                   !point here. In future we can use another interpolation.
      IBP%interpdir = 0

      allocate(IBP%IntPoints(0))

      return

    endif


    if (abs(IBP%distx)<(xU(xi+1)-xU(xi-1))/1000._KND) then      !if too close to the boundary, set the distance to 0
      IBP%distx = 0
      dirx = 0
      IBP%dirx = 0
    endif

    if (abs(IBP%disty)<(yU(yj+1)-yU(yj-1))/1000._KND) then
      IBP%disty = 0
      diry = 0
      IBP%diry = 0
    endif

    if (abs(IBP%distz)<(zU(zk+1)-zU(zk-1))/1000._KND) then
      IBP%distz = 0
      dirz = 0
      IBP%dirz = 0
    endif


    if (dirx==0) then                                               !If there is a cell free of solid bodies in the direction
      free100=.false.
    else
      free100 = (Utype(xi+IBP%dirx,yj,zk)==0)
    endif

    if (diry==0) then
      free010=.false.
    else
      free010 = (Utype(xi,yj+IBP%diry,zk)==0)
    endif

    if (dirz==0) then
      free001=.false.
    else
      free001 = (Utype(xi,yj,zk+IBP%dirz)==0)
    endif


    nx = 0
    ny = 0
    nz = 0
    n1 = 0                                                          ! n1 number of neighbouring cells free of solid bodies
    n2 = 0                                                          ! n2 number of nonzero coordinate directions to boundary

    if (Utype(xi+1,yj,zk)<=0) then
      n1 = n1+1
      nx = nx+1
    endif

    if (Utype(xi-1,yj,zk)<=0) then
      n1 = n1+1
      nx = nx+1
    endif

    if (Utype(xi,yj+1,zk)<=0) then
      n1 = n1+1
      ny = ny+1
    endif

    if (Utype(xi,yj-1,zk)<=0) then
      n1 = n1+1
      ny = ny+1
    endif

    if (Utype(xi,yj,zk+1)<=0) then
      n1 = n1+1
      nz = nz+1
    endif

    if (Utype(xi,yj,zk-1)<=0) then
      n1 = n1+1
      nz = nz+1
    endif

    if (dirx/=0) n2 = n2+1
    if (diry/=0) n2 = n2+1
    if (dirz/=0) n2 = n2+1


    if (n1>n2) then   ! If too many free directions, treat as directly on the boundary,
      IBP%interp = 0  ! because we are probably at some edge.
      IBP%distx = 0
      IBP%disty = 0
      IBP%distz = 0
      dirx = 0
      diry = 0
      dirz = 0
      IBP%dirx = 0
      IBP%diry = 0
      IBP%dirz = 0
    endif

      !At least in one direction free  cells on both sides.
      !Unresolvably small object or a thin wall -> treat as directly on the boundary.
    if (nx>1.or.ny>1.or.nz>1) then
      IBP%interp = 0
      IBP%interpdir = 0
      IBP%distx = 0
      IBP%disty = 0
      IBP%distz = 0
      dirx = 0
      diry = 0
      dirz = 0
      IBP%dirx = 0
      IBP%diry = 0
      IBP%dirz = 0
    endif


    if ((dirx==0.and.diry==0.and.dirz==0).or.((.not.free001).and.(.not.free010).and.(.not.free100))) then
                              !If no free direction, the boundary is here.
      IBP%interp = 0
      IBP%interpdir = 0
      IBP%distx = 0
      IBP%disty = 0
      IBP%distz = 0
      dirx = 0
      diry = 0
      dirz = 0
      IBP%dirx = 0
      IBP%diry = 0
      IBP%dirz = 0

    elseif ((free100.and.dirx==1).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then !Three free directions,

      IBP%interp = 9                                                                           !use trilinear interpolation.
      IBP%interpdir = 0

    elseif ((.not.free100).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then !Two free directions y and z,
                                                                                     !use bilinear interpolation normal to x.
      IBP%interp = 6
      IBP%interpdir = 1

      if (dirx==1) then                                       !If dirx /= 0 then forget the x component of dist vector
        t = TSolidBody_NearestOnLineOut(SB,x,y,z,x,y+IBP%disty,z+IBP%distz) ! and find an intersection of the new vector with the boundary
        dirx = 0
        IBP%dirx = 0
        IBP%distx = 0
        IBP%disty = IBP%disty*t
        IBP%distz = IBP%distz*t
      endif

    elseif ((.not.free010).and.(free100.and.dirx==1).and.(free001.and.dirz==1)) then  !the same normal to y

      IBP%interp = 6
      IBP%interpdir = 2

      if (diry==1) then
        t = TSolidBody_NearestOnLineOut(SB,x,y,z,x+IBP%distx,y,z+IBP%distz)
        diry = 0
        IBP%diry = 0
        IBP%distx = IBP%distx*t
        IBP%disty = 0
        IBP%distz = IBP%distz*t
      endif

    elseif ((.not.free001).and.(free100.and.dirx==1).and.(free010.and.diry==1)) then  !the same normal to z

      IBP%interp = 6
      IBP%interpdir = 3

      if (dirz==1) then
        t = TSolidBody_NearestOnLineOut(SB,x,y,z,x+IBP%distx,y+IBP%disty,z)
        dirz = 0
        IBP%dirz = 0
        IBP%distx = IBP%distx*t
        IBP%disty = IBP%disty*t
        IBP%distz = 0
      endif

    elseif (free100.and.dirx==1) then  !Only one free direction, use linear interpolation in direction x.

      IBP%interp = 3
      IBP%interpdir = 1

      if (diry==1.or.dirz==1) then                   !If other dir components nonzero, delete them and
        t = TSolidBody_NearestOnLineOut(SB,x,y,z,x+IBP%distx,y,z)  !find an intersection of the new vector with the boundary
        diry = 0
        dirz = 0
        IBP%diry = 0
        IBP%dirz = 0
        IBP%distx = t*IBP%distx
        IBP%disty = 0
        IBP%distz = 0
      endif

    elseif (free010.and.diry==1) then               !the same in y

      IBP%interp = 3
      IBP%interpdir = 2

      if (dirx==1.or.dirz==1) then
        t = TSolidBody_NearestOnLineOut(SB,x,y,z,x,y+IBP%disty,z)
        dirx = 0
        dirz = 0
        IBP%dirx = 0
        IBP%dirz = 0
        IBP%distx = 0
        IBP%disty = t*IBP%disty
        IBP%distz = 0
      endif

    elseif (free001.and.dirz==1) then               !the same in z

      IBP%interp = 3
      IBP%interpdir = 3

      if (dirx==1.or.diry==1) then
        t = TSolidBody_NearestOnLineOut(SB,x,y,z,x,y,z+IBP%distz)
        dirx = 0
        diry = 0
        IBP%dirx = 0
        IBP%diry = 0
        IBP%distx = 0
        IBP%disty = 0
        IBP%distz = t*IBP%distz
      endif

    else                                           ! We should have find some interpolation and not come here!
      write(*,*) "Assertion error"
      write(*,*) "free100",free100
      write(*,*) "free010",free010
      write(*,*) "free001",free001
      write(*,*) "dirx",dirx
      write(*,*) "diry",diry
      write(*,*) "dirz",dirz
      write(*,*) "i",IBP%xi
      write(*,*) "j",IBP%yj
      write(*,*) "j",IBP%zk
      write(*,*) "distx",IBP%distx
      write(*,*) "disty",IBP%disty
      write(*,*) "distz",IBP%distz
      write(*,*) "x",x
      write(*,*) "y",y
      write(*,*) "z",z
      write(*,*) "component",IBP%component
      stop
    endif

    allocate(IBP%IntPoints(IBP%interp))

    call TVelIBPoint_InterpolationCoefs(IBP,xU,yU,zU)

  end subroutine TVelIBPoint_Create








!
!   subroutine TVelIBPoint_InterpolationCoefs(IBP,xU,yU,zU)
!     type(TVelIBpoint),intent(inout)     :: IBP
!     real(KND),dimension(-2:),intent(in) :: xU,yU,zU
!     real(KND) p0,p1,p2,p3,p4,p5
!     integer xi,yj,zk,dirx,diry,dirz
!
!     xi = IBP%xi
!     yj = IBP%yj
!     zk = IBP%zk
!
!     dirx = IBP%dirx
!     diry = IBP%diry
!     dirz = IBP%dirz
!
!     if (IBP%interp==2) then
!
!
!         IBP%IntPoints(1)%xi = xi+dirx
!         IBP%IntPoints(1)%yj = yj+diry
!         IBP%IntPoints(1)%zk = zk+dirz
!
!         IBP%IntPoints(2)%xi = xi+2*dirx
!         IBP%IntPoints(2)%yj = yj+2*diry
!         IBP%IntPoints(2)%zk = zk+2*dirz
!
!         if (IBP%interpdir==1) then
!           p0 = abs(IBP%distx)
!           p1 = abs(xU(xi+dirx)-xU(xi))
!           p2 = abs(xU(xi+2*dirx)-xU(xi))
!         elseif (IBP%interpdir==2) then
!           p0 = abs(IBP%disty)
!           p1 = abs(yU(yj+diry)-yU(yj))
!           p2 = abs(yU(yj+2*diry)-yU(yj))
!         else
!           p0 = abs(IBP%distz)
!           p1 = abs(zU(zk+dirz)-zU(zk))
!           p2 = abs(zU(zk+2*dirz)-zU(zk))
!         endif
!
!         call IBLinInterpolationCoefs(IBP%IntPoints%coef,p0,p1,p2)
!
!
!     elseif (IBP%interp==4) then
!
!
!         if (IBP%interpdir==1) then
!
!           IBP%IntPoints(1)%xi = xi
!           IBP%IntPoints(1)%yj = yj + diry
!           IBP%IntPoints(1)%zk = zk
!
!           IBP%IntPoints(2)%xi = xi
!           IBP%IntPoints(2)%yj = yj + 2*diry
!           IBP%IntPoints(2)%zk = zk
!
!           IBP%IntPoints(3)%xi = xi
!           IBP%IntPoints(3)%yj = yj
!           IBP%IntPoints(3)%zk = zk + dirz
!
!           IBP%IntPoints(4)%xi = xi
!           IBP%IntPoints(4)%yj = yj
!           IBP%IntPoints(4)%zk = zk + 2*dirz
!
!           p0 = abs(IBP%disty)
!           p1 = abs(yU(yj+diry)-yU(yj))
!           p2 = abs(yU(yj+2*diry)-yU(yj))
!
!           p3 = abs(IBP%distz)
!           p4 = abs(zU(zk+dirz)-zU(zk))
!           p5 = abs(zU(zk+2*dirz)-zU(zk))
!
!
!         elseif (IBP%interpdir==2) then
!
!           IBP%IntPoints(1)%xi = xi
!           IBP%IntPoints(1)%yj = yj
!           IBP%IntPoints(1)%zk = zk + dirz
!
!           IBP%IntPoints(2)%xi = xi
!           IBP%IntPoints(2)%yj = yj
!           IBP%IntPoints(2)%zk = zk + 2*dirz
!
!           IBP%IntPoints(3)%xi = xi + dirx
!           IBP%IntPoints(3)%yj = yj
!           IBP%IntPoints(3)%zk = zk
!
!           IBP%IntPoints(4)%xi = xi + 2*dirx
!           IBP%IntPoints(4)%yj = yj
!           IBP%IntPoints(4)%zk = zk
!
!           p0 = abs(IBP%distz)
!           p1 = abs(zU(zk+dirz)-zU(zk))
!           p2 = abs(zU(zk+2*dirz)-zU(zk))
!
!           p3 = abs(IBP%distx)
!           p4 = abs(xU(xi+dirx)-xU(xi))
!           p5 = abs(xU(xi+2*dirx)-xU(xi))
!
!         else
!
!           IBP%IntPoints(1)%xi = xi + dirx
!           IBP%IntPoints(1)%yj = yj
!           IBP%IntPoints(1)%zk = zk
!
!           IBP%IntPoints(2)%xi = xi + 2*dirx
!           IBP%IntPoints(2)%yj = yj
!           IBP%IntPoints(2)%zk = zk
!
!           IBP%IntPoints(3)%xi = xi
!           IBP%IntPoints(3)%yj = yj + diry
!           IBP%IntPoints(3)%zk = zk
!
!           IBP%IntPoints(4)%xi = xi
!           IBP%IntPoints(4)%yj = yj + 2*diry
!           IBP%IntPoints(4)%zk = zk
!
!           p0 = abs(IBP%distx)
!           p1 = abs(xU(xi+dirx)-xU(xi))
!           p2 = abs(xU(xi+2*dirx)-xU(xi))
!
!           p3 = abs(IBP%disty)
!           p4 = abs(yU(yj+diry)-yU(yj))
!           p5 = abs(yU(yj+2*diry)-yU(yj))
!
!         endif
!
!         call IBLinInterpolationCoefs(IBP%IntPoints(1:2)%coef,p0,p1,p2)
!
!         call IBLinInterpolationCoefs(IBP%IntPoints(3:4)%coef,p3,p4,p5)
!
!         !beta 1 and  2 in Peller et al., doi:1.1002/fld.1227
!         p1 = p3/p0
!         p2 = p0/p3
!         !sum of betas
!         p4 = p1 + p2
!
!         IBP%IntPoints(1:2)%coef = IBP%IntPoints(1:2)%coef * p1/p4
!
!         IBP%IntPoints(3:4)%coef = IBP%IntPoints(3:4)%coef * p2/p4
!
!     elseif (IBP%interp==7) then
!
!
!         IBP%IntPoints(1)%xi = xi+dirx
!         IBP%IntPoints(1)%xi = yj
!         IBP%IntPoints(1)%xi = zk
!
!         IBP%IntPoints(2)%xi = xi
!         IBP%IntPoints(2)%xi = yj+diry
!         IBP%IntPoints(2)%xi = zk
!
!         IBP%IntPoints(3)%xi = xi+dirx
!         IBP%IntPoints(3)%xi = yj+diry
!         IBP%IntPoints(3)%xi = zk
!
!         IBP%IntPoints(4)%xi = xi
!         IBP%IntPoints(4)%xi = yj
!         IBP%IntPoints(4)%xi = zk+dirz
!
!         IBP%IntPoints(5)%xi = xi+dirx
!         IBP%IntPoints(5)%xi = yj
!         IBP%IntPoints(5)%xi = zk+dirz
!
!         IBP%IntPoints(6)%xi = xi
!         IBP%IntPoints(6)%xi = yj+diry
!         IBP%IntPoints(6)%xi = zk+dirz
!
!         IBP%IntPoints(7)%xi = xi+dirx
!         IBP%IntPoints(7)%xi = yj+diry
!         IBP%IntPoints(7)%xi = zk+dirz
!
!         p0 = abs(IBP%distx)
!         p1 = abs(IBP%disty)
!         p2 = abs(IBP%distz)
!         p3 = abs(xU(xi+dirx)-xU(xi))
!         p4 = abs(yU(yj+diry)-yU(yj))
!         p5 = abs(zU(zk+dirz)-zU(zk))
!
!         call IBTriLinInterpolationCoefs(IBP%IntPoints%coef,p0,p1,p2,p3,p4,p5)
!
!     endif
!
!     IBP%Intpoints%coef = 0
!   end subroutine TVelIBPoint_InterpolationCoefs









  subroutine TVelIBPoint_InterpolationCoefs(IBP,xU,yU,zU)
    type(TVelIBpoint),intent(inout)     :: IBP
    real(KND),dimension(-2:),intent(in) :: xU,yU,zU
    real(KND) xr, yr, zr, x(0:3), y(0:3), z(0:3)
    real(KND) b1, b2, b3, c
    integer xi,yj,zk,dirx,diry,dirz

    xi = IBP%xi
    yj = IBP%yj
    zk = IBP%zk

    dirx = IBP%dirx
    diry = IBP%diry
    dirz = IBP%dirz

    if (IBP%interp==3) then


        IBP%IntPoints(1)%xi = xi+dirx
        IBP%IntPoints(1)%yj = yj+diry
        IBP%IntPoints(1)%zk = zk+dirz

        IBP%IntPoints(2)%xi = xi+2*dirx
        IBP%IntPoints(2)%yj = yj+2*diry
        IBP%IntPoints(2)%zk = zk+2*dirz

        IBP%IntPoints(3)%xi = xi+3*dirx
        IBP%IntPoints(3)%yj = yj+3*diry
        IBP%IntPoints(3)%zk = zk+3*dirz

        if (IBP%interpdir==1) then
          xr = xU(xi) + IBP%distx
          x  = (/ xU(xi), xU(xi+dirx), xU(xi+2*dirx), xU(xi+3*dirx)  /)
        elseif (IBP%interpdir==2) then
          xr = yU(yj) + IBP%disty
          x  = (/ yU(yj), yU(yj+diry), yU(yj+2*diry), yU(yj+3*diry)  /)
        else
          xr = zU(zk) + IBP%distz
          x  = (/ zU(zk), zU(zk+dirz), zU(zk+2*dirz), zU(zk+3*dirz)  /)
        endif

        call IBLeastSquare2InterpolationCoefs(IBP%IntPoints%coef,xr,x)


    elseif (IBP%interp==6) then


        if (IBP%interpdir==1) then

          IBP%IntPoints(1)%xi = xi
          IBP%IntPoints(1)%yj = yj + diry
          IBP%IntPoints(1)%zk = zk

          IBP%IntPoints(2)%xi = xi
          IBP%IntPoints(2)%yj = yj + 2*diry
          IBP%IntPoints(2)%zk = zk

          IBP%IntPoints(3)%xi = xi
          IBP%IntPoints(3)%yj = yj + 3*diry
          IBP%IntPoints(3)%zk = zk

          IBP%IntPoints(4)%xi = xi
          IBP%IntPoints(4)%yj = yj
          IBP%IntPoints(4)%zk = zk + dirz

          IBP%IntPoints(5)%xi = xi
          IBP%IntPoints(5)%yj = yj
          IBP%IntPoints(5)%zk = zk + 2*dirz

          IBP%IntPoints(6)%xi = xi
          IBP%IntPoints(6)%yj = yj
          IBP%IntPoints(6)%zk = zk + 3*dirz

          xr = yU(yj) + IBP%disty
          x  = (/ yU(yj), yU(yj+diry), yU(yj+2*diry), yU(yj+3*diry)  /)

          yr = zU(zk) + IBP%distz
          y  = (/ zU(zk), zU(zk+dirz), zU(zk+2*dirz), zU(zk+3*dirz)  /)


        elseif (IBP%interpdir==2) then

          IBP%IntPoints(1)%xi = xi
          IBP%IntPoints(1)%yj = yj
          IBP%IntPoints(1)%zk = zk + dirz

          IBP%IntPoints(2)%xi = xi
          IBP%IntPoints(2)%yj = yj
          IBP%IntPoints(2)%zk = zk + 2*dirz

          IBP%IntPoints(3)%xi = xi
          IBP%IntPoints(3)%yj = yj
          IBP%IntPoints(3)%zk = zk + 3*dirz

          IBP%IntPoints(4)%xi = xi + dirx
          IBP%IntPoints(4)%yj = yj
          IBP%IntPoints(4)%zk = zk

          IBP%IntPoints(5)%xi = xi + 2*dirx
          IBP%IntPoints(5)%yj = yj
          IBP%IntPoints(5)%zk = zk

          IBP%IntPoints(6)%xi = xi + 3*dirx
          IBP%IntPoints(6)%yj = yj
          IBP%IntPoints(6)%zk = zk

          xr = zU(zk) + IBP%distz
          x  = (/ zU(zk), zU(zk+dirz), zU(zk+2*dirz), zU(zk+3*dirz)  /)

          yr = xU(xi) + IBP%distx
          y  = (/ xU(xi), xU(xi+dirx), xU(xi+2*dirx), xU(xi+3*dirx)  /)

        else

          IBP%IntPoints(1)%xi = xi + dirx
          IBP%IntPoints(1)%yj = yj
          IBP%IntPoints(1)%zk = zk

          IBP%IntPoints(2)%xi = xi + 2*dirx
          IBP%IntPoints(2)%yj = yj
          IBP%IntPoints(2)%zk = zk

          IBP%IntPoints(3)%xi = xi + 3*dirx
          IBP%IntPoints(3)%yj = yj
          IBP%IntPoints(3)%zk = zk

          IBP%IntPoints(4)%xi = xi
          IBP%IntPoints(4)%yj = yj + diry
          IBP%IntPoints(4)%zk = zk

          IBP%IntPoints(5)%xi = xi
          IBP%IntPoints(5)%yj = yj + 2*diry
          IBP%IntPoints(5)%zk = zk

          IBP%IntPoints(6)%xi = xi
          IBP%IntPoints(6)%yj = yj + 3*diry
          IBP%IntPoints(6)%zk = zk

          xr = xU(xi) + IBP%distx
          x  = (/ xU(xi), xU(xi+dirx), xU(xi+2*dirx), xU(xi+3*dirx)  /)

          yr = yU(yj) + IBP%disty
          y  = (/ yU(yj), yU(yj+diry), yU(yj+2*diry), yU(yj+3*diry)  /)

        endif

        call IBLeastSquare2InterpolationCoefs(IBP%IntPoints(1:3)%coef,xr,x)

        call IBLeastSquare2InterpolationCoefs(IBP%IntPoints(4:6)%coef,yr,y)

        !beta 1 and  2 in eq. 17-19 in Peller et al., doi:1.1002/fld.1227
        b1 = abs(y(0)-yr)/abs(x(0)-xr)
        b2 = 1/b1
        !sum of betas
        c = b1 + b2

        IBP%IntPoints(1:3)%coef = IBP%IntPoints(1:3)%coef * b1/c

        IBP%IntPoints(4:6)%coef = IBP%IntPoints(4:6)%coef * b2/c

    elseif (IBP%interp==9) then

        IBP%IntPoints(1)%xi = xi + dirx
        IBP%IntPoints(1)%yj = yj
        IBP%IntPoints(1)%zk = zk

        IBP%IntPoints(2)%xi = xi + 2*dirx
        IBP%IntPoints(2)%yj = yj
        IBP%IntPoints(2)%zk = zk

        IBP%IntPoints(3)%xi = xi + 3*dirx
        IBP%IntPoints(3)%yj = yj
        IBP%IntPoints(3)%zk = zk

        IBP%IntPoints(4)%xi = xi
        IBP%IntPoints(4)%yj = yj + diry
        IBP%IntPoints(4)%zk = zk

        IBP%IntPoints(5)%xi = xi
        IBP%IntPoints(5)%yj = yj + 2*diry
        IBP%IntPoints(5)%zk = zk

        IBP%IntPoints(6)%xi = xi
        IBP%IntPoints(6)%yj = yj + 3*diry
        IBP%IntPoints(6)%zk = zk

        IBP%IntPoints(7)%xi = xi
        IBP%IntPoints(7)%yj = yj
        IBP%IntPoints(7)%zk = zk + dirz

        IBP%IntPoints(8)%xi = xi
        IBP%IntPoints(8)%yj = yj
        IBP%IntPoints(8)%zk = zk + 2*dirz

        IBP%IntPoints(9)%xi = xi
        IBP%IntPoints(9)%yj = yj
        IBP%IntPoints(9)%zk = zk + 3*dirz

        xr = xU(xi) + IBP%distx
        x  = (/ xU(xi), xU(xi+dirx), xU(xi+2*dirx), xU(xi+3*dirx)  /)

        yr = yU(yj) + IBP%disty
        y  = (/ yU(yj), yU(yj+diry), yU(yj+2*diry), yU(yj+3*diry)  /)

        zr = zU(zk) + IBP%distz
        z  = (/ yU(zk), zU(zk+dirz), zU(zk+2*dirz), zU(zk+3*dirz)  /)

        call IBLeastSquare2InterpolationCoefs(IBP%IntPoints(1:3)%coef,xr,x)

        call IBLeastSquare2InterpolationCoefs(IBP%IntPoints(4:6)%coef,yr,y)

        !beta 1 and  2 in eq. 17-19 in Peller et al., doi:1.1002/fld.1227
        b1 = abs(y(0)-yr)*abs(z(0)-zr)/abs(x(0)-xr)
        b2 = abs(z(0)-zr)*abs(x(0)-xr)/abs(y(0)-yr)
        b3 = abs(x(0)-xr)*abs(y(0)-yr)/abs(z(0)-zr)
        !sum of betas
        c = b1 + b2 + b3

        IBP%IntPoints(1:3)%coef = IBP%IntPoints(1:3)%coef * b1/c

        IBP%IntPoints(4:6)%coef = IBP%IntPoints(4:6)%coef * b2/c

        IBP%IntPoints(7:9)%coef = IBP%IntPoints(7:9)%coef * b3/c

    endif

    IBP%Intpoints%coef = 0
  end subroutine TVelIBPoint_InterpolationCoefs












  subroutine TScalFlIBPoint_Create(IBP,xi,yj,zk)
    type(TScalFlIBPoint),intent(out) :: IBP
    integer,intent(in) :: xi,yj,zk               !grid coordinates of the forcing point
    type(TSolidBody),pointer :: SB
    integer dirx,diry,dirz,dirx2,diry2,dirz2,nfreedirs,ndirs,i
    real(KND) x,y,z,xnear,ynear,znear,distx,disty,distz,t,tx,ty,tz
    logical freep00,free0p0,free00p,freem00,free0m0,free00m

    x = xPr(xi)                                   !physical coordinates of the forcing point
    y = yPr(yj)
    z = zPr(zk)
    call SetCurrentSB(SB,Prtype(xi,yj,zk))
    call NearestOut(SB,xnear,ynear,znear,x,y,z)

    IBP%xi = xi
    IBP%yj = yj
    IBP%zk = zk

    distx = xnear-x                               !distance to the nearest point on the boundary
    disty = ynear-y
    distz = znear-z
    dirx = nint(sign(1.0_KND,distx))              !direction to the boundary point
    diry = nint(sign(1.0_KND,disty))
    dirz = nint(sign(1.0_KND,distz))


    if (abs(distx)<(xPr(xi+1)-xPr(xi-1))/100._KND) then  !If very close to the boundary, set the boundary here
      distx = 0
      dirx = 0
    endif

    if (abs(disty)<(yPr(yj+1)-yPr(yj-1))/100._KND) then
      disty = 0
      diry = 0
    endif

    if (abs(distz)<(zW(zk+1)-zW(zk-1))/100._KND) then
      distz = 0
      dirz = 0
    endif


    freep00 = (Prtype(xi+1,yj,zk)<=0)   !logicals denoting if the cell in plus x direction is free of SB
    freem00 = (Prtype(xi-1,yj,zk)<=0)
    free0p0 = (Prtype(xi,yj+1,zk)<=0)
    free0m0 = (Prtype(xi,yj-1,zk)<=0)
    free00p = (Prtype(xi,yj,zk+1)<=0)
    free00m = (Prtype(xi,yj,zk-1)<=0)



    if (.not.(freep00.or.freem00)) then
      dirx = 0
      distx = 0
    end if
    if (.not.(free0p0.or.free0m0)) then
      diry = 0
      disty = 0
    end if
    if (.not.(free00p.or.free00m)) then
      dirz = 0
      distz = 0
    end if

    nfreedirs = 0
    ndirs = 0

    if (freep00) nfreedirs = nfreedirs+1
    if (free0p0) nfreedirs = nfreedirs+1
    if (free00p) nfreedirs = nfreedirs+1
    if (freem00) nfreedirs = nfreedirs+1
    if (free0m0) nfreedirs = nfreedirs+1
    if (free00m) nfreedirs = nfreedirs+1

    if (dirx/=0) ndirs = ndirs+1
    if (diry/=0) ndirs = ndirs+1
    if (dirz/=0) ndirs = ndirs+1

    if (nfreedirs>ndirs.or.(freep00.and.freem00).or.(free0p0.and.free0m0).or.(free00p.and.free00m)) then
      IBP%interp = 0    !If more free spaces than directions, or free space in both oposite directions,
      distx = 0         !Assume is an unresolvably small feature and treat as on the boundary.
      disty = 0
      distz = 0
      dirx = 0
      diry = 0
      dirz = 0
      dirx = 0
      diry = 0
      dirz = 0
    endif


    if (dirx==0.and.diry==0.and.dirz==0) then   !If dir>0 treat as on the boundary.

      dirx2 = 0
      diry2 = 0
      dirz2 = 0
      if (freep00) dirx2 = dirx2+1     !Find a free neighbouring cell and interpolate from there.
      if (freem00) dirx2 = dirx2-1
      if (free0p0) diry2 = diry2+1
      if (free0m0) diry2 = diry2-1
      if (free00p) dirz2 = dirz2+1
      if (free00m) dirz2 = dirz2-1

      if (dirx2==0.and.diry2==0.and.dirz2==0) then  !probably a very thin body, just average free points around

        IBP%interp = nfreedirs

        allocate(IBP%IntPoints(nfreedirs))

        i = 1
        IBP%dist = 0
        if (freep00) then
          IBP%IntPoints(i) = TInterpolationPoint(xi + 1, yj, zk, 1._KND / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + xPr(xi+1)-xPr(xi)
        endif
        if (freem00) then
          IBP%IntPoints(i) = TInterpolationPoint(xi - 1, yj, zk, 1._KND / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + xPr(xi)-xPr(xi-1)
        endif
        if (free0p0) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj + 1, zk, 1._KND / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + yPr(yj+1)-yPr(yj)
        endif
        if (free0m0) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj - 1, zk, 1._KND / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + yPr(yj)-yPr(yj-1)
        endif
        if (free00p) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj, zk + 1, 1._KND / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + zPr(zk+1)-zPr(zk)
        endif
        if (free00m) then
          IBP%IntPoints(i) = TInterpolationPoint(xi, yj, zk - 1, 1._KND / nfreedirs)
          i = i + 1
          IBP%dist = IBP%dist + zPr(zk)-zPr(zk-1)
        endif
        IBP%dist = IBP%dist / nfreedirs

      else  !some edge

        allocate(IBP%IntPoints(1))

        IBP%IntPoints(1)%xi = IBP%xi+dirx2
        IBP%IntPoints(1)%yj = IBP%yj+diry2
        IBP%IntPoints(1)%zk = IBP%zk+dirz2
        IBP%IntPoints%coef = 1._KND
        IBP%interp = 1
        IBP%dist = sqrt((x-xPr(IBP%xi+dirx2))**2+(y-yPr(IBP%yj+diry2))**2+(z-zPr(IBP%zk+dirz2))**2)

      endif

    elseif (nfreedirs==1.or.ndirs==1) then     !Only one free point outside, interpolate from there.

      allocate(IBP%IntPoints(1))

      IBP%IntPoints(1)%xi = IBP%xi+dirx
      IBP%IntPoints(1)%yj = IBP%yj+diry
      IBP%IntPoints(1)%zk = IBP%zk+dirz
      IBP%IntPoints%coef = 1._KND
      IBP%interp = 1
      IBP%dist = sqrt((x-xPr(IBP%xi+dirx))**2+(y-yPr(IBP%yj+diry))**2+(z-zPr(IBP%zk+dirz))**2)

    elseif (nfreedirs==2.or.ndirs==2) then    !Two free directions, use bilinear interpolation in the plane contaning these two neigbours.

      allocate(IBP%IntPoints(2))

      if (dirx==0) then     !plane yz

        if (abs(disty/distz)<abs(yPr(yj+diry)-y)/abs(zPr(zk+dirz)-z)) then  !Which gridline does the line from this point
          t = (zPr(zk+dirz)-z)/distz                                           !to the boundary intersect? Then use the points
          IBP%IntPoints(1)%xi = IBP%xi
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(disty*t)/abs(yPr(IBP%yj+diry)-yPr(IBP%yj))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk+dirz
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((disty*t)**2+(distz*t)**2)
        else
          t = (yPr(yj+diry)-y)/disty
          IBP%IntPoints(1)%xi = IBP%xi
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(distz*t)/abs(zPr(IBP%zk+dirz)-zPr(IBP%zk))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj+diry
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((disty*t)**2+(distz*t)**2)
        endif

      elseif (diry==0) then     !plane xz

        if (abs(distx/distz)<abs(xPr(xi+dirx)-x)/abs(zPr(zk+dirz)-z)) then
          t = (zPr(zk+dirz)-z)/distz
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(distx*t)/abs(xPr(IBP%xi+dirx)-xPr(IBP%xi))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk+dirz
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(distz*t)**2)
        else
          t = (xPr(xi+dirx)-x)/distx
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj
          IBP%IntPoints(1)%zk = IBP%zk+dirz
          IBP%IntPoints(1)%coef = abs(distz*t)/abs(zPr(IBP%zk+dirz)-zPr(IBP%zk))
          IBP%IntPoints(2)%xi = IBP%xi+dirx
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(distz*t)**2)
        endif

      else                 !plane xy

        if (abs(distx/disty)<abs(xPr(xi+dirx)-x)/abs(yPr(yj+diry)-y)) then
          t = (yPr(yj+diry)-y)/disty
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk
          IBP%IntPoints(1)%coef = abs(distx*t)/abs(xPr(IBP%xi+dirx)-xPr(IBP%xi))
          IBP%IntPoints(2)%xi = IBP%xi
          IBP%IntPoints(2)%yj = IBP%yj+diry
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(disty*t)**2)
        else
          t = (xPr(xi+dirx)-x)/distx
          IBP%IntPoints(1)%xi = IBP%xi+dirx
          IBP%IntPoints(1)%yj = IBP%yj+diry
          IBP%IntPoints(1)%zk = IBP%zk
          IBP%IntPoints(1)%coef = abs(disty*t)/abs(yPr(IBP%yj+diry)-yPr(IBP%yj))
          IBP%IntPoints(2)%xi = IBP%xi+dirx
          IBP%IntPoints(2)%yj = IBP%yj
          IBP%IntPoints(2)%zk = IBP%zk
          IBP%IntPoints(2)%coef = 1-IBP%IntPoints(1)%coef
          IBP%interp = 2
          IBP%dist = sqrt((distx*t)**2+(disty*t)**2)
        endif
      endif

    else                         !more than two free directions

      tx = (xPr(xi+dirx)-x)/distx   !a vector pointing from a free point to this forcing point.
      ty = (yPr(yj+diry)-y)/disty
      tz = (zPr(zk+dirz)-z)/distz

      IBP%interp = 4

      allocate(IBP%IntPoints(4))

      if (tx<=ty.and.tx<=tz) then
        !coordinates of interpolation point are therefore
        !xPr(xi+dirx) = x+tx*distx
        !y+tx*disty
        !z+tx*distz
        IBP%dist = sqrt((distx*tx)**2+(disty*tx)**2+(distz*tx)**2)
        IBP%IntPoints(1)%xi = IBP%xi+dirx
        IBP%IntPoints(1)%yj = IBP%yj+diry
        IBP%IntPoints(1)%zk = IBP%zk+dirz
        IBP%IntPoints(1)%coef = abs(disty*tx)/abs(yPr(yj+diry)-y)*&
                                  abs(distz*tx)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(2)%xi = IBP%xi+dirx
        IBP%IntPoints(2)%yj = IBP%yj
        IBP%IntPoints(2)%zk = IBP%zk+dirz
        IBP%IntPoints(2)%coef = (1-abs(disty*tx)/abs(yPr(yj+diry)-y))*&
                                  abs(distz*tx)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(3)%xi = IBP%xi+dirx
        IBP%IntPoints(3)%yj = IBP%yj+diry
        IBP%IntPoints(3)%zk = IBP%zk
        IBP%IntPoints(3)%coef = abs(disty*tx)/abs(yPr(yj+diry)-y)*&
                                  (1-abs(distz*tx)/abs(zPr(zk+dirz)-z))
        IBP%IntPoints(4)%xi = IBP%xi+dirx
        IBP%IntPoints(4)%yj = IBP%yj
        IBP%IntPoints(4)%zk = IBP%zk
        IBP%IntPoints(4)%coef = (1-abs(disty*tx)/abs(yPr(yj+diry)-y))*&
                                  (1-abs(distz*tx)/abs(zPr(zk+dirz)-z))
      elseif (ty<=tx.and.ty<=tz) then
        !the same with ty
        IBP%dist = sqrt((distx*ty)**2+(disty*ty)**2+(distz*ty)**2)
        IBP%IntPoints(1)%xi = IBP%xi+dirx
        IBP%IntPoints(1)%yj = IBP%yj+diry
        IBP%IntPoints(1)%zk = IBP%zk+dirz
        IBP%IntPoints(1)%coef = abs(distx*ty)/abs(xPr(xi+dirx)-x)*&
                                  abs(distz*ty)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(2)%xi = IBP%xi
        IBP%IntPoints(2)%yj = IBP%yj+diry
        IBP%IntPoints(2)%zk = IBP%zk+dirz
        IBP%IntPoints(2)%coef = (1-abs(distx*ty)/abs(xPr(xi+dirx)-x))*&
                                  abs(distz*ty)/abs(zPr(zk+dirz)-z)
        IBP%IntPoints(3)%xi = IBP%xi+dirx
        IBP%IntPoints(3)%yj = IBP%yj+diry
        IBP%IntPoints(3)%zk = IBP%zk
        IBP%IntPoints(3)%coef = abs(distx*ty)/abs(xPr(xi+dirx)-x)*&
                                  (1-abs(distz*ty)/abs(zPr(zk+dirz)-z))
        IBP%IntPoints(4)%xi = IBP%xi
        IBP%IntPoints(4)%yj = IBP%yj+diry
        IBP%IntPoints(4)%zk = IBP%zk
        IBP%IntPoints(4)%coef = (1-abs(distx*ty)/abs(xPr(xi+dirx)-x))*&
                                  (1-abs(distz*ty)/abs(zPr(zk+dirz)-z))
      else
        !the same with tz
        IBP%dist = sqrt((distx*tz)**2+(disty*tz)**2+(distz*tz)**2)
        IBP%IntPoints(1)%xi = IBP%xi+dirx
        IBP%IntPoints(1)%yj = IBP%yj+diry
        IBP%IntPoints(1)%zk = IBP%zk+dirz
        IBP%IntPoints(1)%coef = abs(distx*tz)/abs(xPr(xi+dirx)-x)*&
                                  abs(disty*tz)/abs(yPr(yj+diry)-y)
        IBP%IntPoints(2)%xi = IBP%xi
        IBP%IntPoints(2)%yj = IBP%yj+diry
        IBP%IntPoints(2)%zk = IBP%zk+dirz
        IBP%IntPoints(2)%coef = (1-abs(distx*tz)/abs(xPr(xi+dirx)-x))*&
                                  abs(disty*tz)/abs(yPr(yj+diry)-y)
        IBP%IntPoints(3)%xi = IBP%xi+dirx
        IBP%IntPoints(3)%yj = IBP%yj
        IBP%IntPoints(3)%zk = IBP%zk+dirz
        IBP%IntPoints(3)%coef = abs(distx*tz)/abs(xPr(xi+dirx)-x)*&
                                  (1-abs(disty*tz)/abs(yPr(yj+diry)-y))
        IBP%IntPoints(4)%xi = IBP%xi
        IBP%IntPoints(4)%yj = IBP%yj
        IBP%IntPoints(4)%zk = IBP%zk+dirz
        IBP%IntPoints(4)%coef = (1-abs(distx*tz)/abs(xPr(xi+dirx)-x))*&
                                  (1-abs(disty*tz)/abs(yPr(yj+diry)-y))
      endif

    endif

    IBP%temperatureflux = SB%temperatureflux

  end subroutine TScalFlIBPoint_Create






  subroutine TVelIBPoint_DeallocateList(first)
    type(TVelIBPoint),pointer :: first,node,tmp

    node => first

    if (associated(node)) then
      tmp => node
      node => node%next

      deallocate(node)
    endif

    node => null()

  end subroutine TVelIBPoint_DeallocateList



  subroutine TScalFlIBPoint_DeallocateList(first)
    type(TScalFlIBPoint),pointer ::first,node,tmp

    node => first

    if (associated(node)) then
      tmp => node
      node => node%next

      deallocate(node)
    endif

    node => null()

  end subroutine TScalFlIBPoint_DeallocateList




  subroutine TVelIBPoint_AddToList(IBP)
    type(TVelIBPoint),intent(in) :: IBP

    if (.not.associated(LastIBPoint)) then
     allocate(FirstIBPoint)
     FirstIBPoint = IBP
     LastIBPoint => FirstIBPoint
    else
     allocate(LastIBPoint%next)
     LastIBPoint%next = IBP
     LastIBPoint => LastIBPoint%next
    endif
  end subroutine TVelIBPoint_AddToList



  subroutine TScalFlIBPoint_AddToList(SIBP)
    type(TScalFlIBPoint),intent(in) :: SIBP

    if (.not.associated(LastScalFlIBPoint)) then
     allocate(FirstScalFlIBPoint)
     FirstScalFlIBPoint = SIBP
     LastScalFlIBPoint => FirstScalFlIBPoint
    else
     allocate(LastScalFlIBPoint%next)
     LastScalFlIBPoint%next = SIBP
     LastScalFlIBPoint => LastScalFlIBPoint%next
    endif
  end subroutine TScalFlIBPoint_AddToList







  subroutine MoveIBPointsToArray
    type(TVelIBPoint),pointer :: CurrentIBPoint
    type(TScalFlIBPoint),pointer :: CurrentScalFlIBPoint
    integer :: i,iU,iV,iW

    NUIBPoints = 0
    NVIBPoints = 0
    NWIBPoints = 0
    NScalFlIBPoints = 0

    CurrentIBPoint => FirstIBPoint
    do
     if (associated(CurrentIBPoint)) then
       if (CurrentIBPoint%component==1) then
         NUIBPoints = NUIBPoints + 1
       elseif (CurrentIBPoint%component==2) then
         NVIBPoints = NVIBPoints + 1
       else
         NWIBPoints = NWIBPoints + 1
       endif
     else
       exit
     endif
     CurrentIBPoint => CurrentIBPoint%next
    enddo

    CurrentScalFlIBPoint => FirstScalFlIBPoint
    do
     if (associated(CurrentScalFlIBPoint)) then
       NScalFlIBPoints = NScalFlIBPoints + 1
     else
       exit
     endif
     CurrentScalFlIBPoint => CurrentScalFlIBPoint%next
    enddo

    allocate(UIBPoints(NUIBPoints))
    allocate(VIBPoints(NVIBPoints))
    allocate(WIBPoints(NWIBPoints))
    allocate(ScalFlIBPoints(NScalFlIBPoints))

    CurrentIBPoint => FirstIBPoint
    iU = 0
    iV = 0
    iW = 0
    do
     if (associated(CurrentIBPoint)) then
       if (CurrentIBPoint%component==1) then
         iU = iU + 1
         UIBPoints(iU) = CurrentIBPoint
       elseif (CurrentIBPoint%component==2) then
         iV = iV + 1
         VIBPoints(iV) = CurrentIBPoint
       else
         iW = iW + 1
         WIBPoints(iW) = CurrentIBPoint
       endif
     else
       exit
     endif
     CurrentIBPoint => CurrentIBPoint%next
    enddo

    CurrentScalFlIBPoint => FirstScalFlIBPoint
    i = 0
    do
     if (associated(CurrentScalFlIBPoint)) then
       i = i + 1
       ScalFlIBPoints(i) = CurrentScalFlIBPoint
     else
       exit
     endif
     CurrentScalFlIBPoint => CurrentScalFlIBPoint%next
    enddo

  end subroutine MoveIBPointsToArray

























  subroutine SetCurrentSB(SB,n)
    type(TSolidBody),pointer :: SB
    integer(SINT),intent(in) :: n
    integer i

    SB => FirstSB

    do i = 2,n
      if (associated(SB%next)) then
        SB => SB%next
      else
        SB => null()
        exit
      endif
    enddo
  end subroutine SetCurrentSB







  real(KND) function PointDist(x1,y1,z1,x2,y2,z2)
    real(KND),intent(in) :: x1,y1,z1,x2,y2,z2

    PointDist = SQRT((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
  end function PointDist



  subroutine TLine_Nearest(L,xnear,ynear,znear,x,y,z)
    type(TLine),intent(in) :: L
    real(KND),intent(out)  :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    real(KND) t

    if (L%a/=0 .or. L%b/=0 .or. L%c/=0) then
      t = ( L%a*(x-L%xc) + L%b*(y-L%yc) + L%c*(z-L%zc) ) / (L%a**2 + L%b**2 + L%c**2)
    else
      t = 0
    endif

    xnear = L%xc + L%a * t
    ynear = L%yc + L%b * t
    znear = L%zc + L%c * t
  end subroutine TLine_Nearest












  pure logical function TPlane_Inside(PL,x,y,z,eps)
    type(TPlane),intent(in) :: PL
    real(KND),intent(in) :: x,y,z
    real(KND),optional,intent(in) ::eps
    real(KND) eps2
    logical ins

    if (present(eps)) then
      eps2 = eps
    else
      eps2 = MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND)
    endif

    if (PL%GL) then
      if (PL%a*x+PL%b*y+PL%c*z+PL%d>=-eps2) then
       ins = .true.
      else
       ins = .false.
      endif
    else
      if (PL%a*x+PL%b*y+PL%c*z+PL%d<=eps2) then
       ins = .true.
      else
       ins = .false.
      endif
    endif
    TPlane_Inside = ins
  end function TPlane_Inside



  pure logical function TPolyhedron_Inside(PH,x,y,z,eps)
    type(TPolyhedron),intent(in) :: PH
    real(KND),intent(in) :: x,y,z
    real(KND),optional,intent(in) ::eps
    logical ins
    integer i

    if (PH%nplanes>0) then
      ins = .true.
      do i = 1,PH%nplanes
       if (present(eps)) then
        ins = Inside(PH%Planes(i),x,y,z,eps)
       else
        ins = Inside(PH%Planes(i),x,y,z)
       endif
       if (.not. ins) exit
      enddo
    else
      ins = .false.
    endif

    TPolyhedron_Inside = ins
  end function TPolyhedron_Inside



  pure logical function TBall_Inside(B,x,y,z,eps)
    type(TBall),intent(in) :: B
    real(KND),intent(in) :: x,y,z
    real(KND),intent(in),optional ::eps
    real(KND) :: eps2
    logical ins

    if (present(eps)) then
     eps2 = eps
    else
     eps2 = MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND)
    endif

    if ((B%xc-x)**2+(B%yc-y)**2+(B%zc-z)**2<=(B%r+eps2)**2) then
     ins = .true.
    else
     ins = .false.
    endif

    TBall_Inside = ins
  end function TBall_Inside



  pure real(KND) function LineDist(x,y,z,xl,yl,zl,a,b,c)
    real(KND),intent(in) :: x,y,z,xl,yl,zl,a,b,c
    real(KND) t

    if (((a/=0).or.(b/=0)).or.(c/=0)) then
     t = (a*(x-xl)+b*(y-yl)+c*(z-zl))/(a**2+b**2+c**2)
    else
     t = 0
    endif

    LineDist = SQRT((xl+a*t-x)**2+(yl+b*t-y)**2+(zl+c*t-z)**2)
  end function LineDist



  pure logical function TCylJacket_Inside(J,x,y,z,eps)
    type(TCylJacket),intent(in) :: J
    real(KND),intent(in) :: x,y,z
    real(KND),intent(in),optional ::eps
    real(KND) eps2
    logical ins


    if (present(eps)) then
     eps2 = eps
    else
     eps2 = MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND)
    endif

    if (LineDist(x,y,z,j%xc,j%yc,j%zc,J%a,J%b,J%c)<=j%r+eps2) then
     ins = .true.
    else
     ins = .false.
    endif

    TCylJacket_Inside = ins
  end function TCylJacket_Inside



  pure logical function TCylinder_Inside(C,x,y,z,eps)
   type(TCylinder),intent(in) :: C
   real(KND),intent(in) :: x,y,z
   real(KND),intent(in),optional ::eps
   logical ins

    ins = .true.

    if (.not.Inside(C%Jacket,x,y,z,eps)) ins = .false.

    if (ins.and.associated(C%Plane1)) then
           if (.not.Inside(C%Plane1,x,y,z,eps)) ins = .false.
    endif
    if (ins.and.associated(C%Plane2)) then
          if (.not.Inside(C%Plane2,x,y,z,eps)) ins = .false.
    endif

    TCylinder_Inside = ins
  end function TCylinder_Inside



   subroutine TTerrain_GridCoords(x2,y2,xi,yj,comp)
    real(KND),intent(in) :: x2,y2
    integer,intent(out) :: xi,yj,comp
    real(KND) x,y,distPr,distU,distV
    integer xPri,yPrj,xUi,yVj,i

    x = x2
    y = y2

    if (gridtype==uniformgrid) then

        xPri = min(max(nint( (x - xU(0))/dxmin + 0.5_KND ),1),Prnx+1)

        yPrj = min(max(nint( (y - yV(0))/dymin + 0.5_KND ),1),Prny+1)

        xUi = min(max(nint( (x-xU(0))/dxmin ),0), Unx+1)

        yVj = min(max(nint( (y-yV(0))/dymin ),0), Vny+1)
    else


        xPri = Prnx+1
        do i=0,Prnx+1
         if (xU(i)>=x) then
                       xPri = i
                       exit
                      endif
        enddo

        yPrj = Prny+1
        do i=0,Prny+1
         if (yV(i)>=y) then
                       yPrj = i
                       exit
                      endif
        enddo

        xUi = Prnx+1
        do i = 0,Prnx+1
         if (xPr(i+1)>=x) then
                       xUi = i
                       exit
                      endif
        enddo

        yVj = Prny+1
        do i = 0,Prny+1
         if (yPr(i+1)>=y) then
                       yVj = i
                       exit
                      endif
        enddo

    endif

    distPr = (x-xPr(xPri))**2+(y-yPr(yPrj))**2
    distU = (x-xU(xUi))**2+(y-yPr(yPrj))**2
    distV = (x-xPr(xPri))**2+(y-yV(yVj))**2

    if (distU<distPr.and.distU<distV) then
     xi = xUi
     yj = yPrj
     comp = 1
    elseif (distV<distPr.and.distV<distU) then
     xi = xPri
     yj = yVj
     comp = 2
    else
     xi = xPri
     yj = yPrj
     comp = 3
    endif
  end subroutine TTerrain_GridCoords



   logical function TTerrain_Inside(T,x,y,z,eps)
    type(TTerrain),intent(in) :: T
    real(KND),intent(in) :: x,y,z
    real(KND),intent(in),optional ::eps
    real(KND) eps2
    logical ins
    integer xi,yj,comp

    if (present(eps)) then
     eps2 = eps
    else
     eps2 = MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND)
    endif

    ins = .false.
    call TTerrain_GridCoords(x,y,xi,yj,comp)

    if (comp==1) then
     if (z<=T%UPoints(xi,yj)%elev+eps2) ins = .true.
    elseif (comp==2) then
     if (z<=T%VPoints(xi,yj)%elev+eps2) ins = .true.
    elseif (comp==3) then
     if (z<=T%PrPoints(xi,yj)%elev+eps2) ins = .true.
    endif

    TTerrain_Inside = ins
  end function TTerrain_Inside



   logical function TSolidBody_Inside(SB,x,y,z,eps)
    type(TSolidBody),intent(in) :: SB
    real(KND),intent(in) :: x,y,z
    real(KND),intent(in),optional :: eps
    real(KND) x2,y2,z2
    real(KND) lx,ly,lz

    x2 = x
    y2 = y
    z2 = z

    lx = xU(Prnx) - xU(0)
    ly = yV(Prny) - yV(0)
    lz = zW(Prnz) - zW(0)

    if (Btype(Ea)==PERIODIC.and.x2>xU(Prnx+1)) x2 = x2-lx
    if (Btype(No)==PERIODIC.and.y2>yV(Prny+1)) y2 = y2-ly
    if (Btype(To)==PERIODIC.and.z2>zW(Prnz+1)) z2 = z2-lz

    if (Btype(We)==PERIODIC.and.x2<xU(0)) x2 = x2+lx
    if (Btype(So)==PERIODIC.and.y2<yV(0)) y2 = y2+ly
    if (Btype(Bo)==PERIODIC.and.z2<zW(0)) z2 = z2+lz

    select case (SB%typeofbody)

     case (Polyhedron)

       if (present(eps)) then
         TSolidBody_Inside = Inside(SB%Polyhedron,x2,y2,z2,eps)
       else
         TSolidBody_Inside = Inside(SB%Polyhedron,x2,y2,z2)
       endif

     case (Ball)

       if (present(eps)) then
         TSolidBody_Inside = Inside(SB%Ball,x2,y2,z2,eps)
       else
         TSolidBody_Inside = Inside(SB%Ball,x2,y2,z2)
       endif

     case (Cylinder)

       if (present(eps)) then
         TSolidBody_Inside = Inside(SB%Cylinder,x2,y2,z2,eps)
       else
         TSolidBody_Inside = Inside(SB%Cylinder,x2,y2,z2)
       endif

     case (Terrain)

       if (present(eps)) then
         TSolidBody_Inside = Inside(SB%Terrain,x2,y2,z2,eps)
       else
         TSolidBody_Inside = Inside(SB%Terrain,x2,y2,z2)
       endif

     case default
       TSolidBody_Inside = .false.

    end select

  end function TSolidBody_Inside


  subroutine TPlane_Nearest(PL,xnear,ynear,znear,x,y,z)
    type(TPlane),intent(in) :: PL
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) ::x,y,z
    real(KND) t

    if (((PL%a/=0).or.(PL%b/=0)).or.(PL%c/=0)) then
     t = -(PL%a*x+PL%b*y+PL%c*z+PL%d)/(PL%a**2+PL%b**2+PL%c**2)
    else
     t = 0
    endif
    xnear = x+PL%a*t
    ynear = y+PL%b*t
    znear = z+PL%c*t
  end subroutine TPlane_Nearest



  subroutine TCylJacket_Nearest(J,xnear,ynear,znear,x,y,z)
    type(TCylJacket),intent(in) :: J
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    real(KND) t,xl,yl,zl,a,b,c

    call Nearest(TLine(J%xc,J%yc,J%zc,J%a,J%b,J%c),xl,yl,zl,x,y,z)

    a = x-xl
    b = y-yl
    c = z-zl
    t = J%r/SQRT(a**2+b**2+c**2)

    xnear = a*t+xl
    ynear = b*t+yl
    znear = c*t+zl
  end subroutine TCylJacket_Nearest



  subroutine TCylinder_Nearest(C,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
   type(TCylinder),intent(in) :: C
   real(KND),intent(out) :: xnear,ynear,znear
   real(KND),intent(in) :: x,y,z
   real(KND) xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2

   !!!Only for Planes perpendicular to jacket!!!!

   if (associated(C%Plane1)) then
       if (associated(C%Plane2)) then
          call Nearest(C%Jacket,xJ,yJ,zJ,x,y,z)
          call Nearest(C%Plane1,xP1,yP1,zP1,x,y,z)
          call Nearest(C%Plane2,xP2,yP2,zP2,x,y,z)
          if (Inside(C%Jacket,x,y,z)) then
             if (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
                xnear = xP1
                ynear = yP1
                znear = zP1
             else
                xnear = xP2
                ynear = yP2
                znear = zP2
             endif
          elseif (Inside(C%Plane1,x,y,z).and.Inside(C%Plane2,x,y,z)) then
                xnear = xJ
                ynear = yJ
                znear = zJ
          elseif (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
           call Nearest(C%Jacket,xnear,ynear,znear,xP1,yP1,zP1)
          else
           call Nearest(C%Jacket,xnear,ynear,znear,xP2,yP2,zP2)
          endif
       else
          call Nearest(C%Jacket,xJ,yJ,zJ,x,y,z)
          call Nearest(C%Plane1,xP1,yP1,zP1,x,y,z)
          if (Inside(C%Jacket,x,y,z)) then
                xnear = xP1
                ynear = yP1
                znear = zP1
          elseif (Inside(C%Plane1,x,y,z)) then
                xnear = xJ
                ynear = yJ
                znear = zJ
          else
           call Nearest(C%Jacket,xnear,ynear,znear,xP1,yP1,zP1)
         endif
       endif

    else

      call Nearest(C%Jacket,xnear,ynear,znear,x,y,z)

    endif
  end subroutine TCylinder_Nearest

  subroutine TBall_Nearest(Bl,xnear,ynear,znear,x,y,z)
    type(TBall),intent(in) :: Bl
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    real(KND) t,a,b,c

    a = x-Bl%xc
    b = y-Bl%yc
    c = z-Bl%zc
    t = Bl%r/SQRT(a**2+b**2+c**2)
    xnear = a*t+Bl%xc
    ynear = b*t+Bl%yc
    znear = c*t+Bl%zc
  end subroutine TBall_Nearest


  subroutine TPolyhedron_Nearest(PH,xnear,ynear,znear,x,y,z)
    type(TPolyhedron),intent(in) :: PH
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    real(KND) :: dists(PH%nplanes),xP(PH%nplanes),yP(PH%nplanes),zP(PH%nplanes)
    real(KND) :: minv
    real(KND) :: ailine,biline,ciline
    real(KND) :: x0iline,y0iline,z0iline,xln,yln,zln,p
    integer   :: inearest,inearest2,inearest3
    integer   :: i

    integer,dimension(3)    :: ipivot,work2             !arguments of the LAPACK
    real(KND),dimension(3)  :: xg,bg,R,C,ferr,berr      !  routine xGESVX
    real(KND),dimension(12) :: work
    real(KND),dimension(3,3) :: ag,af                    !
    integer   :: info                                   !
    real(KND) :: rcond                                  !
    character :: equed                                  !

    external :: DGESVX, SGESVX

   !Vzdalenosti od rovin, pokud je nejbl. bod roviny uvnitr jine, nebo puv. bod na vnitrni strane -> vzd. *-1
   !Pokud je nejbl. body +- eps. (norm. vektor!,dxmin/100) uvnitr a vne polyh -> hotovo
   !Jinak 2. nejbl. rovina v abs. hodnote -> prusecnice a nejbl bod na ni
   !Nejbl. bod na prusecnici. Pokud +-eps. uvnitr, (najit vekt. v rovine  kolme na prusecnic)-:hotovo
   !Jinak iterativne najit bod na prusecnici uvnitr

    do i = 1,PH%nplanes
     call Nearest(PH%Planes(i),xP(i),yP(i),zP(i),x,y,z)
     dists(i) = SQRT((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
     if (Inside(PH%Planes(i),x,y,z)) dists(i) = -ABS(dists(i))
    enddo
    !find nearest plane with
    inearest = 0
    minv = huge(minv)
    do i = 1,PH%nplanes
     if (dists(i)>=0.and.dists(i)<minv) then
       inearest = i
       minv = dists(i)
     endif
    enddo

    if (inearest==0) then
      write(*,*) "no nearest"
      return
    endif

    if (Inside(PH,xP(inearest),yP(inearest),zP(inearest),MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND))) then
     xnear = xP(inearest)
     ynear = yP(inearest)
     znear = zP(inearest)
     return
    endif

    dists = abs(dists)

    inearest2 = 0
    minv = huge(minv)
    do i = 1,PH%nplanes
     if (i/=inearest.and.dists(i)<minv) then
       inearest2 = i
       minv = dists(i)
     endif
    enddo

    inearest3 = 0
    minv = huge(minv)
    do i = 1,PH%nplanes
     if (i/=inearest.and.i/=inearest2.and.dists(i)<minv) then
       inearest3 = i
       minv = dists(i)
     endif
    enddo

    ailine = PH%Planes(inearest)%b*PH%Planes(inearest2)%c-PH%Planes(inearest)%c*PH%Planes(inearest2)%b
    biline = PH%Planes(inearest)%c*PH%Planes(inearest2)%a-PH%Planes(inearest)%a*PH%Planes(inearest2)%c
    ciline = PH%Planes(inearest)%a*PH%Planes(inearest2)%b-PH%Planes(inearest)%b*PH%Planes(inearest2)%a

    if (abs(ailine)<=epsilon(ailine).and.abs(biline)<=epsilon(biline).and.abs(ciline)<=epsilon(ciline)) then
     write(*,*) "cross product 0"
     write(*,*) "numbers of planes",inearest,inearest2
     write(*,*) "x,y,z",x,y,z
     write(*,*) "pl1",PH%Planes(inearest)%a,PH%Planes(inearest)%b,PH%Planes(inearest)%c,PH%Planes(inearest)%d
     write(*,*) "pl2",PH%Planes(inearest2)%a,PH%Planes(inearest2)%b,PH%Planes(inearest2)%c,PH%Planes(inearest2)%d

     call Nearest(PH%Planes(inearest),xln,yln,zln,x,y,z)
     call Nearest(PH%Planes(inearest2),xnear,ynear,znear,x,y,z)

     if ((xln-x)**2+(yln-y)**2+(zln-z)**2<(xln-x)**2+(yln-y)**2+(zln-z)**2) then
       xnear = xln
       ynear = yln
       znear = zln
     end if

     return

    endif

    p = SQRT(ailine**2+biline**2+ciline**2)
    ailine = ailine/p
    biline = biline/p
    ciline = ciline/p

    if (abs(ciline)>=0.1_KND) then

       ag = reshape(source = (/ PH%Planes(inearest)%a,PH%Planes(inearest2)%a,0._KND,&
                           PH%Planes(inearest)%b,PH%Planes(inearest2)%b,0._KND,&
                           PH%Planes(inearest)%c,PH%Planes(inearest2)%c,1._KND /),&
                 shape = (/ 3,3 /))

       bg = (/ -PH%Planes(inearest)%d,-PH%Planes(inearest2)%d,(zW(Wnz+1)+zW(0))/2._KND /)

       if (KND==DBL) then
        call DGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       else
        call SGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       endif

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    elseif (abs(biline)>=0.1_KND) then

       ag = reshape(source = (/ PH%Planes(inearest)%a,PH%Planes(inearest2)%a,0._KND,&
                           PH%Planes(inearest)%b,PH%Planes(inearest2)%b,1._KND,&
                           PH%Planes(inearest)%c,PH%Planes(inearest2)%c,0._KND /),&
                 shape = (/ 3,3 /))

       bg = (/ -PH%Planes(inearest)%d,-PH%Planes(inearest2)%d,(yV(Vny+1)+yV(0))/2._KND /)

       if (KND==DBL) then
        call DGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       else
        call SGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       endif

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    else

       ag = reshape(source = (/ PH%Planes(inearest)%a,PH%Planes(inearest2)%a,1._KND,&
                           PH%Planes(inearest)%b,PH%Planes(inearest2)%b,0._KND,&
                           PH%Planes(inearest)%c,PH%Planes(inearest2)%c,0._KND /),&
                  shape = (/ 3,3 /))

       bg = (/ -PH%Planes(inearest)%d,-PH%Planes(inearest2)%d,(xU(Unx+1)+xU(0))/2._KND /)

       if (KND==DBL) then
        call DGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       else
        call SGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       endif

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    endif

    call Nearest(TLine(x0iline,y0iline,z0iline,ailine,biline,ciline),xln,yln,zln,x,y,z)

    if (Inside(PH,xln,yln,zln,min(dxmin/1000._KND,dymin/10000._KND,dzmin/10000._KND))) then
     xnear = xln
     ynear = yln
     znear = zln
     return
    endif



    ag = reshape(source = (/ PH%Planes(inearest)%a,PH%Planes(inearest2)%a,PH%Planes(inearest3)%a,&
                        PH%Planes(inearest)%b,PH%Planes(inearest2)%b,PH%Planes(inearest3)%b,&
                        PH%Planes(inearest)%c,PH%Planes(inearest2)%c,PH%Planes(inearest3)%c /),&
               shape = (/ 3,3 /))
    bg = (/ -PH%Planes(inearest)%d,-PH%Planes(inearest2)%d,-PH%Planes(inearest3)%d /)

    if (KND==DBL) then
     call DGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
    else
     call SGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
    endif

    xnear = xg(1)
    ynear = xg(2)
    znear = xg(3)

  end subroutine TPolyhedron_Nearest




  subroutine TTerrain_Nearest(T,xnear,ynear,znear,x,y,z)
    type(TTerrain),intent(in) :: T
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    type(TPlane) :: Pl
    integer     :: xi,yj,comp
    real(KND)   :: a,b,zloc
    xnear = x
    ynear = y
    znear = z

    call TTerrain_GridCoords(x,y,xi,yj,comp)

    if (comp==1) then  !Construct a tangent plane using first derivatives

     zloc = T%UPoints(xi,yj)%elev
     a = (T%PrPoints(xi+1,yj)%elev - T%PrPoints(xi,yj)%elev) / dxU(xi)
     b = (T%UPoints(xi,yj+1)%elev - T%UPoints(xi,yj-1)%elev) / (yPr(yj+1)-yPr(yj-1))

     Pl%a = a
     Pl%b = b
     Pl%c = -1._KND
     Pl%d = -a*x-b*y+zloc

     call Nearest(Pl,xnear,ynear,znear,x,y,z)

    elseif (comp==2) then

     zloc = T%VPoints(xi,yj)%elev
     a = (T%VPoints(xi+1,yj)%elev - T%VPoints(xi-1,yj)%elev) / (xPr(xi+1)-xPr(xi-1))
     b = (T%PrPoints(xi,yj+1)%elev - T%PrPoints(xi,yj)%elev) / dyV(yj)

     Pl%a = a
     Pl%b = b
     Pl%c = -1._KND
     Pl%d = -a*x-b*y+zloc

     call Nearest(Pl,xnear,ynear,znear,x,y,z)

    elseif (comp==3) then

     zloc = T%PrPoints(xi,yj)%elev
     a = (T%UPoints(xi,yj)%elev - T%UPoints(xi-1,yj)%elev) / dxPr(xi)
     b = (T%VPoints(xi,yj)%elev - T%VPoints(xi,yj-1)%elev) / dyPr(yj)

     Pl%a = a
     Pl%b = b
     Pl%c = -1._KND
     Pl%d = - a*x - b*y + zloc

     call Nearest(Pl,xnear,ynear,znear,x,y,z)

    endif
  end subroutine TTerrain_Nearest


  subroutine TSolidBody_Nearest(SB,xnear,ynear,znear,x,y,z)
    type(TSolidBody),intent(in) :: SB
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    select case (SB%typeofbody)
     case (1)
      call Nearest(SB%Polyhedron,xnear,ynear,znear,x,y,z)
     case (2)
      call Nearest(SB%Ball,xnear,ynear,znear,x,y,z)
     case (3)
      call Nearest(SB%Cylinder,xnear,ynear,znear,x,y,z)
     case (4)
      call Nearest(SB%Terrain,xnear,ynear,znear,x,y,z)
     case default
      xnear = -huge(znear)
      ynear = -huge(znear)
      znear = -huge(znear)
    end select
  end subroutine TSolidBody_Nearest















  subroutine TPolyhedron_NearestOut(PH,xnear,ynear,znear,x,y,z)
    type(TPolyhedron),intent(in) :: PH
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    real(KND) dists(PH%nplanes),xP(PH%nplanes),yP(PH%nplanes),zP(PH%nplanes),minv
    integer inearest,i

    dists = huge(minv)
    do i = 1,PH%nplanes
     call Nearest(PH%Planes(i),xP(i),yP(i),zP(i),x,y,z)
     dists(i) = SQRT((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
    enddo

    inearest = 0
    minv = huge(minv)
    do i = 1,PH%nplanes
     if (dists(i)>=0.and.dists(i)<minv) then
       inearest = i
       minv = dists(i)
     endif
    enddo

    xnear = xP(inearest)
    ynear = yP(inearest)
    znear = zP(inearest)
  end subroutine TPolyhedron_NearestOut



  subroutine TCylinder_NearestOut(C,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
    type(TCylinder),intent(in) :: C
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    real(KND) xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2

    if (associated(C%Plane1)) then
      call Nearest(C%Plane1,xP1,yP1,zP1,x,y,z)
    else
     xP1 = sqrt(huge(xP1))/10
     yP1 = sqrt(huge(yP1))/10
     zP1 = sqrt(huge(zP1))/10
    endif
    if (associated(C%Plane2)) then
      call Nearest(C%Plane2,xP2,yP2,zP2,x,y,z)
    else
     xP2 = sqrt(huge(xP2))/10
     yP2 = sqrt(huge(yP2))/10
     zP2 = sqrt(huge(zP2))/10
    endif
    call Nearest(C%Jacket,xJ,yJ,zJ,x,y,z)

    if (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2).and.&
         PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xJ,yJ,zJ)) then
                 xnear = xP1
                 ynear = yP1
                 znear = zP1
    elseif (PointDist(x,y,z,xP2,yP2,zP2)<=PointDist(x,y,z,xP1,yP1,zP1).and.&
             PointDist(x,y,z,xP2,yP2,zP2)<=PointDist(x,y,z,xJ,yJ,zJ)) then
                 xnear = xP2
                 ynear = yP2
                 znear = zP2
    elseif (PointDist(x,y,z,xJ,yJ,zJ)<=PointDist(x,y,z,xP2,yP2,zP2).and.&
             PointDist(x,y,z,xJ,yJ,zJ)<=PointDist(x,y,z,xP1,yP1,zP1)) then
                 xnear = xJ
                 ynear = yJ
                 znear = zJ
    endif
  end subroutine TCylinder_NearestOut



  subroutine TBall_NearestOut(Bl,xnear,ynear,znear,x,y,z)
    type(TBall),intent(in) :: Bl
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    call Nearest(Bl,xnear,ynear,znear,x,y,z)
  end subroutine TBall_NearestOut



  subroutine TTerrain_NearestOut(T,xnear,ynear,znear,x,y,z)
    type(TTerrain),intent(in) :: T
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    call Nearest(T,xnear,ynear,znear,x,y,z)
  end subroutine TTerrain_NearestOut




  subroutine TSolidBody_NearestOut(SB,xnear,ynear,znear,x,y,z)
    type(TSolidBody),intent(in) :: SB
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    xnear = -1e+9;ynear = -1e+9;znear = -1e+9

    select case (SB%typeofbody)
     case (1)
      call NearestOut(SB%Polyhedron,xnear,ynear,znear,x,y,z)
     case (2)
      call NearestOut(SB%Ball,xnear,ynear,znear,x,y,z)
     case (3)
      call NearestOut(SB%Cylinder,xnear,ynear,znear,x,y,z)
     case (4)
      call NearestOut(SB%Terrain,xnear,ynear,znear,x,y,z)
     case default
      xnear = -huge(znear);ynear = -huge(znear);znear = -huge(znear)
    end select
  end subroutine TSolidBody_NearestOut





  real(KND) function TSolidBody_NearestOnLineOut(SB,x,y,z,x2,y2,z2) !Find t, such that x+(x2-x)*t lies on the boundary of the SB
    type(TSolidBody),intent(in) :: SB
    real(KND),intent(in) :: x,y,z,x2,y2,z2

    real(KND) t,t1,t2
    integer i

    t1 = 0
    t2 = 1
    if (Inside(SB,x2,y2,z2)) then
     do
      t2 = t2*2._KND
      if (.not.Inside(SB,x+(x2-x)*t2,y+(y2-y)*t2,z+(z2-z)*t2)) exit
     enddo
    endif
    t = (t1+t2)/2._KND

    do i = 1,20         !The bisection method with maximum 20 iterations (should be well enough)
     if (Inside(SB,x+(x2-x)*t,y+(y2-y)*t,z+(z2-z)*t)) then
      t1 = t
     else
      t2 = t
     endif
     t = (t1+t2)/2._KND
     if (abs(t1-t2)<MIN(dxmin/1000._KND,dymin/1000._KND,dzmin/1000._KND))   exit
    enddo
    TSolidBody_NearestOnLineOut = t
  end function TSolidBody_NearestOnLineOut




  recursive subroutine TSolidBody_DeallocateList(SB)
    type(TSolidbody),pointer :: SB

    if (associated(SB%next)) call TSolidBody_DeallocateList(SB%next)

    if (associated(SB%Terrain)) then
           if (allocated(SB%Terrain%UPoints)) deallocate(SB%Terrain%UPoints)
           if (allocated(SB%Terrain%VPoints)) deallocate(SB%Terrain%VPoints)
           if (allocated(SB%Terrain%PrPoints)) deallocate(SB%Terrain%PrPoints)
           deallocate(SB%Ball)
    endif

    if (associated(SB%Ball)) deallocate(SB%Ball)

    if (associated(SB%Cylinder)) then
           if (associated(SB%Cylinder%Plane1)) deallocate(SB%Cylinder%Plane1)
           if (associated(SB%Cylinder%Plane2)) deallocate(SB%Cylinder%Plane2)
           deallocate(SB%Cylinder)
    endif

    if (associated(SB%Polyhedron)) then
           if (allocated(SB%Polyhedron%Planes)) deallocate(SB%Polyhedron%Planes)
           deallocate(SB%Polyhedron)
    endif

    deallocate(SB)

  end subroutine TSolidBody_DeallocateList




















  subroutine GetSolidBodiesBC

    call FindInsideCells

    call InitImBoundaries

    call GetSolidBodiesWM

  end subroutine GetSolidBodiesBC




  subroutine FindInsideCells
    integer i,j,k
    type(TSolidBody),pointer :: CurrentSB => null()

    !find if the gridpoints lie ins a solid body and write it's number
    !do not nullify the .type arrays, they could have been made nonzero by other unit
    CurrentSB => FirstSB

    do
     if (associated(CurrentSB)) then
      !$omp parallel do private(i,j,k)
      do k = 0,Prnz+1
       do j = 0,Prny+1
        do i = 0,Prnx+1
           if (Inside(CurrentSB,xPr(i),yPr(j),zPr(k))) Prtype(i,j,k) = CurrentSB%numofbody
        enddo
       enddo
      enddo
      !$omp end parallel do
     else
       exit
     endif
     CurrentSB => CurrentSB%next
    enddo

    CurrentSB => FirstSB
    do
     if (associated(CurrentSB)) then
      !$omp parallel do private(i,j,k)
      do k = 0,Unz+1
       do j = 0,Uny+1
         do i = 0,Unx+1
            if (Inside(CurrentSB,xU(i),yPr(j),zPr(k),(dxmin*dymin*dzmin)**(1._KND/3)/1000))&
                     Utype(i,j,k) = CurrentSB%numofbody
         enddo
        enddo
       enddo
      !$omp end parallel do
     else
       exit
     endif
     CurrentSB => CurrentSB%next
    enddo

    CurrentSB => FirstSB
    do
     if (associated(CurrentSB)) then
      !$omp parallel do private(i,j,k)
      do k = 0,Vnz+1
       do j = 0,Vny+1
        do i = 0,Vnx+1
           if (Inside(CurrentSB,xPr(i),yV(j),zPr(k),(dxmin*dymin*dzmin)**(1._KND/3)/1000))&
                     Vtype(i,j,k) = CurrentSB%numofbody
        enddo
       enddo
      enddo
      !$omp end parallel do
     else
       exit
     endif
     CurrentSB => CurrentSB%next
    enddo

    CurrentSB => FirstSB
    do
     if (associated(CurrentSB)) then
      !$omp parallel do private(i,j,k)
      do k = 0,Wnz+1
       do j = 0,Wny+1
        do i = 0,Wnx+1
           if (Inside(CurrentSB,xPr(i),yPr(j),zW(k),(dxmin*dymin*dzmin)**(1._KND/3)/1000))&
                     Wtype(i,j,k) = CurrentSB%numofbody
        enddo
       enddo
      enddo
      !$omp end parallel do
     else
       exit
     endif
     CurrentSB => CurrentSB%next
    enddo

  end subroutine FindInsideCells




  subroutine InitImBoundaries
    type(TVelIBPoint) IBP
    type(TScalFlIBPoint) SIBP
    integer i,j,k


    do k = 1,Unz
     do j = 1,Uny
      do i = 1,Unx
       if (Utype(i,j,k)>0) then
        if (Utype(i+1,j,k)<=0.or.Utype(i-1,j,k)<=0.or.Utype(i,j+1,k)<=0&
          .or.Utype(i,j-1,k)<=0.or.Utype(i,j,k+1)<=0.or.Utype(i,j,k-1)<=0)  then
            call  Create(IBP,i,j,k,xU(-2:),yPr(-2:),zPr(-2:),Utype,1)
            call  AddToList(IBP)
        endif
       endif
      enddo
     enddo
    enddo

    do k = 1,Vnz
     do j = 1,Vny
      do i = 1,Vnx
       if (Vtype(i,j,k)>0) then
        if (Vtype(i+1,j,k)<=0.or.Vtype(i-1,j,k)<=0.or.Vtype(i,j+1,k)<=0&
          .or.Vtype(i,j-1,k)<=0.or.Vtype(i,j,k+1)<=0.or.Vtype(i,j,k-1)<=0)  then
            call  Create(IBP,i,j,k,xPr(-2:),yV(-2:),zPr(-2:),Vtype,2)
            call  AddToList(IBP)
        endif
       endif
      enddo
     enddo
    enddo

    do k = 1,Wnz
     do j = 1,Wny
      do i = 1,Wnx
       if (Wtype(i,j,k)>0) then
        if (Wtype(i+1,j,k)<=0.or.Wtype(i-1,j,k)<=0.or.Wtype(i,j+1,k)<=0&
          .or.Wtype(i,j-1,k)<=0.or.Wtype(i,j,k+1)<=0.or.Wtype(i,j,k-1)<=0)  then
            call  Create(IBP,i,j,k,xPr(-2:),yPr(-2:),zW(-2:),Wtype,3)
            call  AddToList(IBP)
        endif
       endif
      enddo
     enddo
    enddo

    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       if (Prtype(i,j,k)>0) then
        if (Prtype(i+1,j,k)<=0.or.Prtype(i-1,j,k)<=0.or.Prtype(i,j+1,k)<=0&
          .or.Prtype(i,j-1,k)<=0.or.Prtype(i,j,k+1)<=0.or.Prtype(i,j,k-1)<=0)  then
            call  Create(SIBP,i,j,k)
            call  AddToList(SIBP)
        endif
       endif
      enddo
     enddo
    enddo

    call MoveIBPointsToArray

    call DeallocateList(FirstIBPoint)
    call DeallocateList(FirstScalFlIBPoint)

  end subroutine InitImBoundaries




  subroutine GetSolidBodiesWM
    use WallModels, only: WMPoint, AddWMPoint
    type(WMPOINT)            :: WMP
    type(TSolidBody),pointer :: CurrentSB => null()
    integer                  :: neighbours(3,6)
    real(KND)     :: dist,nearx,neary,nearz
    integer       :: i,j,k,m,n,o,p
    integer(SINT) :: nb

    allocate(WMP%depscalar(computescalars))

    neighbours = 0
    neighbours(1,1) =  1
    neighbours(1,2) = -1
    neighbours(2,3) =  1
    neighbours(2,4) = -1
    neighbours(3,5) =  1
    neighbours(3,6) = -1

    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       if (Prtype(i,j,k)<=0) then
        dist = huge(dist)
        nb = 0

        do p=1,6
           m=neighbours(1,p)
           n=neighbours(2,p)
           o=neighbours(3,p)
           if ((Prtype(i+m,j+n,k+o)>0).and.Prtype(i+m,j+n,k+o)/=nb.and.(sum(abs((/m,n,o/)))==1)) then
             call SetCurrentSB(CurrentSB,Prtype(i+m,j+n,k+o))
             call Nearest(CurrentSB,nearx,neary,nearz,xPr(i),yPr(j),zPr(k))
             if (SQRT((nearx-xPr(i))**2+(neary-yPr(j))**2+(nearz-zPr(k))**2)<dist) then
              dist = SQRT((nearx-xPr(i))**2+(neary-yPr(j))**2+(nearz-zPr(k))**2)
              nb = Prtype(i+m,j+n,k+o)
             endif
           endif
        enddo

        if (nb>0) then
          call SetCurrentSB(CurrentSB,nb)
          WMP%xi = i
          WMP%yj = j
          WMP%zk = k
          WMP%distx = nearx-xPr(i)
          WMP%disty = neary-yPr(j)
          WMP%distz = nearz-zPr(k)
          WMP%ustar = 1
          if (CurrentSB%typeofbody==4) then
           if (CurrentSB%Terrain%PrPoints(i,j)%rough) then
             WMP%z0 = CurrentSB%Terrain%PrPoints(i,j)%z0
           else
             WMP%z0 = 0
           endif
          else
           if (CurrentSB%rough) then
             WMP%z0 = CurrentSB%z0
           else
             WMP%z0 = 0
           endif
          endif
          call AddWMPoint(WMP)
        endif
       endif
      enddo
     enddo
    enddo

  end subroutine GetSolidBodiesWM








  subroutine InitSolidBodies

    if (len_trim(obstaclefile)>0) then
      call ReadSolidBodiesFromFile(obstaclefile)
    end if

#ifdef CUSTOMSB
    !An external subroutine, it should use this module and use AddSolidBody to supply
    ! pointers to the new solid bodies.
    call CustomSolidBodies
#endif
  end subroutine InitSolidBodies




  subroutine ReadSolidBodiesFromFile(filename)
    use Strings
    character(*),intent(in) :: filename
    integer unit,io
    character(180) :: line
    type(TSolidBody),pointer :: SB => null()

    unit=20
    !provisional, Solaris Studio still doesn't support newunit.
    open(unit=unit,file=filename,action='read',status='old',iostat=io)

    if (io==0) then
      do
        read(unit,'(a)',iostat=io) line

        if (io/=0) exit

        line = adjustl(line)

        if (len_trim(line)>0) then

          if (upcase(line(1:10))=='POLYHEDRON') then

            call ReadPolyhedron(unit,line(11:),SB)
          end if

          if (associated(SB)) then

            call AddSolidBody(SB)
            SB => null()
          end if

        end if

      end do

      close(unit)

    else

      write(*,*) "Could not open file",filename
      stop

    end if
  end subroutine ReadSolidBodiesFromFile


  subroutine ReadPolyhedron(unit,restline,SB)
    integer,intent(in)       :: unit
    character(*),intent(in)  :: restline
    type(TSolidBody),pointer :: SB
    integer nPlanes,i,io

    read(restline,*,iostat=io) nPlanes

    if (io/=0) then
      write(*,*) "Expected number of planes in polyhedron, received '",trim(restline),"' instead."
      stop
    end if

    allocate(SB)
    allocate(SB%Polyhedron)
    allocate(SB%Polyhedron%Planes(nPlanes))

    SB%Polyhedron%nplanes = nPlanes

    SB%typeofbody = Polyhedron

    do i=1,nPlanes
      call ReadPlane(unit,SB%Polyhedron%Planes(i))
    end do

  end subroutine ReadPolyhedron


  subroutine ReadPlane(unit,Pl)
    use Strings
    integer,intent(in) :: unit
    type(TPlane),intent(inout) :: Pl
    character(180) :: line
    integer io

    read(unit,'(a)',iostat=io) line
    if (io/=0) then
      write(*,*) "Error reading the line with the plane definition."
      stop
    end if

    if (count_multispaces(line) == 4) then
      read(line,*,iostat=io) Pl%a,Pl%b,Pl%c,Pl%d,Pl%gl
    else if (count_multispaces(line) == 6) then
      read(line,*,iostat=io) Pl%a,Pl%b,Pl%c,Pl%d,Pl%gl,Pl%rough, Pl%z0
    else
      io = 999
    end if
    if (io/=0) then
      write(*,*) "Error parsing the line with the plane definition."
      stop
    end if
  end subroutine ReadPlane



  subroutine AddSolidBody(SB)
    type(TSolidBody),pointer :: SB
    type(TSolidBody),pointer :: p

    if (.not. associated(FirstSB)) then
      FirstSB => SB
      SB%numofbody = 1
    else
      p => FirstSB
      do
        if (associated(p%next)) then

          p => p%next

        else

          p%next => SB
          SB%numofbody = p%numofbody + 1
          exit

        end if

      end do
    end if
  end subroutine AddSolidBody
















  pure subroutine IBLeastSquare2InterpolationCoefs(Coefs,xr,x)
    real(KND),intent(out) :: Coefs(:)
    real(KND),intent(in)  :: xr,x(0:)
    real(KND) :: A1, A2, A4
    integer   :: n

    n = size(Coefs)
    if ( n<3 .or. size(x)<=n ) then
      Coefs = 0
    else
      A1 = sum( x(1:) - xr )**2
      A2 = sum( (x(1:)**2 - xr**2) * (x(1:) - xr) )
      A4 = sum( x(1:)**2 - xr**2 )**2

      Coefs = ( ( A2 * (x(1:)**2 - xr**2) - A4 * (x(1:) - xr) ) * (x(0)    - xr)&
            +   ( A2 * (x(1:) - xr) - A1 * (x(1:)**2 - xr**2) ) * (x(0)**2 - xr**2) )&
            / (A2**2 - A1*A4)
    end if

   end subroutine IBLeastSquare2InterpolationCoefs

endmodule GEOMETRIC

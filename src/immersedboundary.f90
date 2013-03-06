module ImmersedBoundaryWM
  use Parameters
  use SolidBodies
  use GeometricShapes, only: TTerrain
  use WallModels, only: WMPoint, AddWMPoint

  private
  public GetSolidBodiesWM
 
contains
  
  subroutine GetSolidBodiesWM
    type(WMPOINT)            :: WMP
    type(TSolidBody),pointer :: CurrentSB => null()
    integer                  :: neighbours(3,6)
    real(KND)     :: dist,nearx,neary,nearz
    integer       :: i,j,k,m,n,o,p
    integer       :: nb

    allocate(WMP%depscalar(computescalars))

    !six triplets [1,0,0], [-1,0,0], [0,1,0],...
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
             call CurrentSB%Nearest(nearx,neary,nearz,xPr(i),yPr(j),zPr(k))
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

          select type (geomshape => CurrentSB%GeometricShape)
            type is (TTerrain)
              if (geomshape%PrPoints(i,j)%rough) then
                WMP%z0 = geomshape%PrPoints(i,j)%z0
              else
                WMP%z0 = 0
              endif
            class default
              if (CurrentSB%rough) then
                WMP%z0 = CurrentSB%z0
              else
                WMP%z0 = 0
              endif
          end select

          call AddWMPoint(WMP)
        endif
       endif
      enddo
     enddo
    enddo

  end subroutine GetSolidBodiesWM

end module ImmersedBoundaryWM




module ImmersedBoundary

  use Parameters
  use Lists
  use SolidBodies
  use ImmersedBoundaryWM

  implicit none

  private

  public TIBPoint, TIBPoint_MomentumSource, TIBPoint_ScalFlSource, TIBPoint_Viscosity, &
         UIBPoints, VIBPoints, WIBPoints, ScalFlIBPoints, &
         GetSolidBodiesBC
         !InitSolidBodies imported from SolidBodies



  type TInterpolationPoint
    integer   :: xi                !coordinates of the interpolation points
    integer   :: yj
    integer   :: zk
    real(KND) :: coef              !interpolation coefficients for the interpolation points
  endtype TInterpolationPoint

  type,extends(TListable) :: TVelIBPoint
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
 !  contains
 !    procedure Create         => TVelIBPoint_Create
  end type TVelIBPoint



  type,extends(TListable) :: TScalFlIBPoint
    integer                :: xi                       !coordinates of the grid point
    integer                :: yj
    integer                :: zk
    real(KND)              :: dist                     !distance to the boundary
    type(TInterpolationPoint),dimension(:),allocatable :: IntPoints !array of interpolation points
    integer                :: interp                   !kind of interpolation 1.. none (1 point outside), 2..linear, 4..bilinear  other values not allowed
    real(KND)              :: temperatureflux = 0      !desired temperature flux
 !  contains
 !    procedure Create         => TScalFlIBPoint_Create
  end type TScalFlIBPoint


  interface Create
    module procedure TVelIBPoint_Create
    module procedure TScalFlIBPoint_Create
  end interface Create



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

  type(TList)  :: UIBPointsList, VIBPointsList, WIBPointsList, ScalFlIBPointsList

  type(TIBPoint),dimension(:),allocatable,save :: UIBPoints, VIBPoints, WIBPoints
  type(TIBPoint),dimension(:),allocatable,save :: ScalFlIBPoints

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



  pure recursive function TIBPoint_Interpolate(IBP,U,lb) result(Uint)
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


  pure recursive function TIBPoint_InterpolateTDiff(IBP,U) result(Uint)
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


  pure recursive function TIBPoint_ScalFlSource(IBP,Scalar,sctype)  result(src)  !Virtual scalar source for the Immersed Boundary Method with prescribed scalar flux on the boundary
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



  pure recursive function TIBPoint_Viscosity(IBP,Viscosity)  result(src)  !Virtual scalar source for the Immersed Boundary Method with prescribed scalar flux on the boundary
    real(KND) :: src

    type(TIBPoint),intent(in) :: IBP
    real(KND),intent(in)      :: Viscosity(-1:,-1:,-1:)

    src = TIBPoint_Interpolate(IBP,Viscosity,-1)

  end function TIBPoint_Viscosity



  pure recursive function TIBPoint_MomentumSource(IBP,U) result(src)
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
    call SB%NearestOut(xnear,ynear,znear,x,y,z)
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

    if (.not.SB%Inside(xU(xi),yPr(yj),zPr(zk)))  then  !For now, if actually outside the body, set an artificial boundary
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
        t = SB%NearestOnLineOut(x,y,z,x,y+IBP%disty,z+IBP%distz) ! and find an intersection of the new vector with the boundary
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
        t = SB%NearestOnLineOut(x,y,z,x+IBP%distx,y,z+IBP%distz)
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
        t = SB%NearestOnLineOut(x,y,z,x+IBP%distx,y+IBP%disty,z)
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
        t = SB%NearestOnLineOut(x,y,z,x+IBP%distx,y,z)  !find an intersection of the new vector with the boundary
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
        t = SB%NearestOnLineOut(x,y,z,x,y+IBP%disty,z)
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
        t = SB%NearestOnLineOut(x,y,z,x,y,z+IBP%distz)
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
    call SB%NearestOut(xnear,ynear,znear,x,y,z)

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


  subroutine MoveIBPointsToArray
    !It would be posible call a generic procedure with the list and array
    ! as an argument, but we want the final arrays not polymorphic.
    integer :: i

    allocate(UIBPoints(UIBPointsList%Len()))
    allocate(VIBPoints(VIBPointsList%Len()))
    allocate(WIBPoints(WIBPointsList%Len()))
    allocate(ScalFlIBPoints(ScalFlIBPointsList%Len()))

    i = 0
    call UIBPointsList%ForEach(CopyIBPoint)
    i = 0
    call VIBPointsList%ForEach(CopyIBPoint)
    i = 0
    call WIBPointsList%ForEach(CopyIBPoint)
    i = 0
    call ScalFlIBPointsList%ForEach(CopyIBPoint)

    contains

      subroutine CopyIBPoint(CurrentIBPoint)
        class(TListable) :: CurrentIBPoint

        i = i + 1

        select type (CurrentIBPoint)
          type is (TVelIBPoint)
            if (CurrentIBPoint%component==1) then
              UIBPoints(i) = CurrentIBPoint
            elseif (CurrentIBPoint%component==2) then
              VIBPoints(i) = CurrentIBPoint
            else
              WIBPoints(i) = CurrentIBPoint
            endif
          type is (TScalFlIBPoint)
            ScalFlIBPoints(i) = CurrentIBPoint
          class default
            stop "Type error in immersed boundary points list."
        end select
      end subroutine

  end subroutine MoveIBPointsToArray


  subroutine AuxNeighbours(Xtype,nx,ny,nz,s) !helper procedure to FindNeighbouringCells
    integer,intent(in)    :: nx,ny,nz,s
    integer,intent(inout) :: Xtype(s:,s:,s:)
    integer i,j,k
    !$omp parallel do private(i,j,k)
    do k = 1,nz
      do j = 1,ny
        do i = 1,nx
          if (Xtype(i,j,k)==0.and.any(Xtype(i-1:i+1,j-1:j+1,k-1:k+1)>0)) &
              Xtype(i,j,k) = -1
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine

  subroutine FindNeighbouringCells
    !sets type of the cells closest to the solid body to -1
    call AuxNeighbours(Prtype,Prnx,Prny,Prnz,0)
    call AuxNeighbours(Utype,Unx,Uny,Unz,-2)
    call AuxNeighbours(Vtype,Vnx,Vny,Vnz,-2)
    call AuxNeighbours(Wtype,Wnx,Wny,Wnz,-2)
  end subroutine FindNeighbouringCells

  



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
            call  UIBPointsList%Add(IBP)
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
            call  VIBPointsList%Add(IBP)
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
            call  WIBPointsList%Add(IBP)
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
            call  ScalFlIBPointsList%Add(SIBP)
        endif
       endif
      enddo
     enddo
    enddo

    call MoveIBPointsToArray

    call UIBPointsList%Deallocate
    call VIBPointsList%Deallocate
    call WIBPointsList%Deallocate
    call ScalFlIBPointsList%Deallocate

  end subroutine InitImBoundaries


  subroutine GetSolidBodiesBC

    call FindInsideCells

    call FindNeighbouringCells

    call InitImBoundaries

    call GetSolidBodiesWM

  end subroutine GetSolidBodiesBC


end module ImmersedBoundary
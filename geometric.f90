module GEOMETRIC

 use PARAMETERS

 implicit none

 private

 public TIBPoint,TScalFlIBPoint,UIBPoints,VIBPoints,WIBPoints,ScalFlIBPoints,&
        InitSolidBodies,GetSolidBodiesBC,&
        IBLinInterpolation,IBBilinInterpolation,IBTriLinInterpolation


 type TLine                            !These object could be implemented using Fortran's 2003 inheritance, or at least using
   real(KND) xc,yc,zc                  !allocatable components. This approach using pointers is more portable and less safe.
   real(KND) a,b,c
 end type TLine

 type TPlane
   real(KND) a,b,c,d      !ax+by+cz+d/=0 for inner half-space
   logical gl             !T > in ineq. above F < in ineq. above
   logical :: rough=.false.!T rough surface, F flat surface
   real(KND) z0           !roughness parameter
  end type TPlane

  type TPolyhedron
    integer :: nplanes = 0
    type(TPlane),dimension(:),allocatable :: Planes !intersection of half-spaces
  end type TPolyhedron

  type TBall
    real(KND) xc,yc,zc,r
    logical :: rough=.false. !T rough surface, F flat surface
    real(KND) z0            !roughness parameter
  end type TBall

  type TCylJacket
    real(KND) xc,yc,zc
    real(KND) a,b,c
    real(KND) r
    logical :: rough=.false. !T rough surface, F flat surface
    real(KND) z0            !roughness parameter
  end type TCylJacket

  type TCylinder
    type(TCylJacket) Jacket
    type(TPlane),pointer :: Plane1 => null() ,Plane2 => null()
  end type TCylinder


 type TTerrainPoint
   real(KND) :: elev = 0
   logical :: rough=.false.
   real(KND) z0
 end type TTerrainPoint


  type TTerrain
    type(TTerrainPoint),dimension(:,:),allocatable :: UPoints,VPoints,PrPoints !allocate with a buffer of width 1 (i.e. 0:Xnx)
  end type TTerrain

 type TSolidBody
   integer numofbody
   integer :: typeofbody = 0                              !0.. none, 1..polyhedron, 2.. ball, 3.. cylinder, 4.. terrain
   type(TPolyhedron),pointer :: Polyhedron => null()    !asociated will be only part writen in typeofbody
   type(TBall),pointer :: Ball => null()
   type(TCylinder),pointer :: Cylinder => null()
   type(TTerrain),pointer :: Terrain => null()
   logical :: rough=.false.                             !T rough surface, F flat surface
   real(KND) z0                                        !roughness parameter
   type(TSolidBody),pointer :: next =>null()
 end type TSolidbody

 type(TSolidBody),pointer :: FirstSB => null()         !First member of the linked list of solid objects

 integer, parameter :: NoneBody = 0, Polyhedron = 1, Ball = 2, Cylinder = 3, Terrain = 4



 type TIBPoint    !NOTE: this implementation is not very efficient, because it computes the int. coefficients again and again. It should be changed similarly to the scalar variant below.
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
   real(KND) :: momsrc         !virtual momentum source
   type(TIBPoint),pointer :: next => null()  !pointer to the next item in the list, not used when in the array in any way
 end type TIBPoint

 type TScalFlIBPoint
   integer                :: xi                       !coordinates of the grid point
   integer                :: yj
   integer                :: zk
   real(KND)              :: dist                     !distance to the boundary
   integer,dimension(4)   :: intpointi                !coordinates of the interpolation points
   integer,dimension(4)   :: intpointj
   integer,dimension(4)   :: intpointk
   real(KND),dimension(4) :: intcoef                  !interpolation coefficients for the interpolation points
   integer                :: interp                   !kind of interpolation 1.. none (1 point outside), 2..linear, 4..bilinear  other values not allowed
   real(KND)              :: scalsrc                  !virtual scalar source
   real(KND)              :: flux = 0                 !desired scalar flux
   type(TScalFlIBPoint),pointer :: next => null()     !pointer to the next item in the list, not used when in the array in any way
 end type TScalFlIBPoint


 type(TIBPoint),pointer,save       :: FirstIBPoint => null(), LastIBPoint => null()
 type(TScalFlIBPoint),pointer,save :: FirstScalFlIBPoint => null(), LastScalFlIBPoint => null()

 type(TIBPoint),dimension(:),allocatable,save       :: UIBPoints, VIBPoints, WIBPoints
 type(TScalFLIBPoint),dimension(:),allocatable,save :: ScalFLIBPoints

 integer, save :: NUIBPoints = 0 ,NVIBPoints = 0, NWIBPoints = 0, NScalFlIBPoints = 0

contains

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

  pure logical function PlaneInside(PL,x,y,z,eps)
    type(TPlane),intent(in) :: PL
    real(KND),intent(in) :: x,y,z
    real(KND),optional,intent(in) ::eps
    real(KND) eps2
    logical interior

    if (present(eps)) then
      eps2 = eps
    else
      eps2 = MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND)
    endif

    if (PL%GL) then
      if (PL%a*x+PL%b*y+PL%c*z+PL%d>=-eps2) then
       interior=.true.
      else
       interior=.false.
      endif
    else
      if (PL%a*x+PL%b*y+PL%c*z+PL%d<=eps2) then
       interior=.true.
      else
       interior=.false.
      endif
    endif
    PlaneInside = interior
  endfunction PlaneInside

  pure logical function PolyhedronInside(PH,x,y,z,eps)
    type(TPolyhedron),intent(in) :: PH
    real(KND),intent(in) :: x,y,z
    real(KND),optional,intent(in) ::eps
    logical interior
    integer i

    if (PH%nplanes>0) then
      interior=.true.
      do i = 1,PH%nplanes
       if (present(eps)) then
        interior = PlaneInside(PH%Planes(i),x,y,z,eps)
       else
        interior = PlaneInside(PH%Planes(i),x,y,z)
       endif
       if (.not. interior) exit
      enddo
    else
      interior=.false.
    endif

    PolyhedronInside = interior
  endfunction PolyhedronInside

  pure logical function BallInside(B,x,y,z,eps)
    type(TBall),intent(in) :: B
    real(KND),intent(in) :: x,y,z
    real(KND),intent(in),optional ::eps
    real(KND) :: eps2
    logical interior

    if (present(eps)) then
     eps2 = eps
    else
     eps2 = MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND)
    endif

    if ((B%xc-x)**2+(B%yc-y)**2+(B%zc-z)**2<=(B%r+eps2)**2) then
     interior=.true.
    else
     interior=.false.
    endif

    BallInside = interior
  endfunction BallInside


  pure real(KND) function LineDist(x,y,z,xl,yl,zl,a,b,c)
    real(KND),intent(in) :: x,y,z,xl,yl,zl,a,b,c
    real(KND) t

    if (((a/=0).or.(b/=0)).or.(c/=0)) then
     t = (a*(x-xl)+b*(y-yl)+c*(z-zl))/(a**2+b**2+c**2)
    else
     t = 0
    endif

    LineDist = SQRT((xl+a*t-x)**2+(yl+b*t-y)**2+(zl+c*t-z)**2)
  endfunction LineDist


  pure logical function JacketInside(J,x,y,z,eps)
    type(TCylJacket),intent(in) :: J
    real(KND),intent(in) :: x,y,z
    real(KND),intent(in),optional ::eps
    real(KND) eps2
    logical interior


    if (present(eps)) then
     eps2 = eps
    else
     eps2 = MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND)
    endif

    if (LineDist(x,y,z,j%xc,j%yc,j%zc,J%a,J%b,J%c)<=j%r+eps2) then
     interior=.true.
    else
     interior=.false.
    endif

    JacketInside = interior
  endfunction JacketInside


  pure logical function CylinderInside(C,x,y,z,eps)
   type(TCylinder),intent(in) :: C
   real(KND),intent(in) :: x,y,z
   real(KND),intent(in),optional ::eps
   logical interior

    interior=.true.

    if (.not.JacketInside(C%Jacket,x,y,z,eps)) interior=.false.

    if (interior.and.associated(C%Plane1)) then
           if (.not.PlaneInside(C%Plane1,x,y,z,eps)) interior=.false.
    endif
    if (interior.and.associated(C%Plane2)) then
          if (.not.PlaneInside(C%Plane2,x,y,z,eps)) interior=.false.
    endif

    CylinderInside = interior
  endfunction CylinderInside



  pure subroutine TerrGridCoords(x2,y2,xi,yj,comp)
    real(KND),intent(in) :: x2,y2
    integer,intent(out) :: xi,yj,comp
    real(KND) x,y,distPr,distU,distV
    integer xPri,yPrj,xUi,yVj,i

    x = x2
    y = y2

    if (gridtype==uniformgrid) then

        xPri = min(max(nint( (x-xU(0)) / (xU(Prnx+1)-xU(0)) * real(Prnx,KND) + 0.5_KND ),1),Prnx)

        yPrj = min(max(nint( (y-yV(0)) / (yV(Prny+1)-yV(0)) * real(Prny,KND) + 0.5_KND ),1),Prny)

        xUi = min(max(nint( (x-xU(0)) / (xU(Prnx+1)-xU(0)) * real(Prnx,KND) ),0), Prnx+1)

        yVj = min(max(nint( (y-yV(0)) / (yV(Prny+1)-yV(0)) * real(Prny,KND) ),0),Prny+1)
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
  end subroutine TerrGridCoords




  pure logical function TerrainInside(T,x,y,z,eps)
    type(TTerrain),intent(in) :: T
    real(KND),intent(in) :: x,y,z
    real(KND),intent(in),optional ::eps
    real(KND) eps2
    logical interior
    integer xi,yj,comp

    if (present(eps)) then
     eps2 = eps
    else
     eps2 = MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND)
    endif

    interior=.false.
    call TerrGridCoords(x,y,xi,yj,comp)
    if (comp==1) then
     if (z<=T%UPoints(xi,yj)%elev+eps2) interior=.true.
    elseif (comp==2) then
     if (z<=T%VPoints(xi,yj)%elev+eps2) interior=.true.
    elseif (comp==3) then
     if (z<=T%PrPoints(xi,yj)%elev+eps2) interior=.true.
    endif

    TerrainInside = interior
  endfunction TerrainInside


  pure logical function SolidBodyInside(SB,x,y,z,eps)
    type(TSolidBody),intent(in) :: SB
    real(KND),intent(in) :: x,y,z
    real(KND),intent(in),optional :: eps
    real(KND) x2,y2,z2

    x2 = x
    y2 = y
    z2 = z

    if (Btype(Ea)==PERIODIC.and.x2>xU(Prnx+1)) x2 = x2-lx
    if (Btype(No)==PERIODIC.and.y2>yV(Prny+1)) y2 = y2-ly
    if (Btype(To)==PERIODIC.and.z2>zW(Prnz+1)) z2 = z2-lz

    if (Btype(We)==PERIODIC.and.x2<xU(0)) x2 = x2+lx
    if (Btype(So)==PERIODIC.and.y2<yV(0)) y2 = y2+ly
    if (Btype(Bo)==PERIODIC.and.z2<zW(0)) z2 = z2+lz

    select case (SB%typeofbody)

     case (Polyhedron)

       if (present(eps)) then
         SolidBodyInside = PolyhedronInside(SB%Polyhedron,x2,y2,z2,eps)
       else
         SolidBodyInside = PolyhedronInside(SB%Polyhedron,x2,y2,z2)
       endif

     case (Ball)

       if (present(eps)) then
         SolidBodyInside = BallInside(SB%Ball,x2,y2,z2,eps)
       else
         SolidBodyInside = BallInside(SB%Ball,x2,y2,z2)
       endif

     case (Cylinder)

       if (present(eps)) then
         SolidBodyInside = CylinderInside(SB%Cylinder,x2,y2,z2,eps)
       else
         SolidBodyInside = CylinderInside(SB%Cylinder,x2,y2,z2)
       endif

     case (Terrain)

       if (present(eps)) then
         SolidBodyInside = TerrainInside(SB%Terrain,x2,y2,z2,eps)
       else
         SolidBodyInside = TerrainInside(SB%Terrain,x2,y2,z2)
       endif

     case default
       SolidBodyInside=.false.

    end select

  endfunction SolidBodyInside





  real(KND) function PointDist(x1,y1,z1,x2,y2,z2)
    real(KND),intent(in) :: x1,y1,z1,x2,y2,z2

    PointDist = SQRT((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
  endfunction PointDist


  subroutine LineNearest(xnear,ynear,znear,x,y,z,xl,yl,zl,a,b,c)
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z,xl,yl,zl,a,b,c
    real(KND) t

    if (((a/=0).or.(b/=0)).or.(c/=0)) then
     t = (a*(x-xl)+b*(y-yl)+c*(z-zl))/(a**2+b**2+c**2)
    else
     t = 0
    endif
    xnear = xl+a*t
    ynear = yl+b*t
    znear = zl+c*t
  end subroutine LineNearest

  subroutine PlaneNearest(PL,xnear,ynear,znear,x,y,z)
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
  end subroutine PlaneNearest



  subroutine JacketNearest(J,xnear,ynear,znear,x,y,z)
    type(TCylJacket),intent(in) :: J
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    real(KND) t,xl,yl,zl,a,b,c

    call LineNearest(xl,yl,zl,x,y,z,J%xc,J%yc,J%zc,J%a,J%b,J%c)
    a = x-xl
    b = y-yl
    c = z-zl
    t = J%r/SQRT(a**2+b**2+c**2)
    xnear = a*t+xl
    ynear = b*t+yl
    znear = c*t+zl
  end subroutine JacketNearest



  subroutine CylinderNearest(C,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
   type(TCylinder),intent(in) :: C
   real(KND),intent(out) :: xnear,ynear,znear
   real(KND),intent(in) :: x,y,z
   real(KND) xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2

   !!!Only for Planes perpendicular to jacket!!!!

   if (associated(C%Plane1)) then
       if (associated(C%Plane2)) then
          call JacketNearest(C%Jacket,xJ,yJ,zJ,x,y,z)
          call PlaneNearest(C%Plane1,xP1,yP1,zP1,x,y,z)
          call PlaneNearest(C%Plane2,xP2,yP2,zP2,x,y,z)
          if (JacketInside(C%Jacket,x,y,z)) then
             if (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
                xnear = xP1
                ynear = yP1
                znear = zP1
             else
                xnear = xP2
                ynear = yP2
                znear = zP2
             endif
          elseif (PlaneInside(C%Plane1,x,y,z).and.PlaneInside(C%Plane2,x,y,z)) then
                xnear = xJ
                ynear = yJ
                znear = zJ
          elseif (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
           call JacketNearest(C%Jacket,xnear,ynear,znear,xP1,yP1,zP1)
          else
           call JacketNearest(C%Jacket,xnear,ynear,znear,xP2,yP2,zP2)
          endif
       else
          call JacketNearest(C%Jacket,xJ,yJ,zJ,x,y,z)
          call PlaneNearest(C%Plane1,xP1,yP1,zP1,x,y,z)
          if (JacketInside(C%Jacket,x,y,z)) then
                xnear = xP1
                ynear = yP1
                znear = zP1
          elseif (PlaneInside(C%Plane1,x,y,z)) then
                xnear = xJ
                ynear = yJ
                znear = zJ
          else
           call JacketNearest(C%Jacket,xnear,ynear,znear,xP1,yP1,zP1)
         endif
       endif

    else

      call JacketNearest(C%Jacket,xnear,ynear,znear,x,y,z)

    endif
  end subroutine CylinderNearest

  subroutine BallNearest(Bl,xnear,ynear,znear,x,y,z)
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
  end subroutine BallNearest


  subroutine PolyhedronNearest(PH,xnear,ynear,znear,x,y,z)
    type(TPolyhedron),intent(in) :: PH
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    real(KND) :: dists(PH%nplanes),xP(PH%nplanes),yP(PH%nplanes),zP(PH%nplanes)
    real(KND) :: minv
    real(KND) :: ailine,biline,ciline
    real(KND) :: x0iline,y0iline,z0iline,xln,yln,zln,p
    integer   :: nearest,nearest2,nearest3
    integer   :: i,j

    integer,dimension(3)    :: ipivot,work2             !arguments of the LAPACK
    real(KND),dimension(3)  :: xg,bg,R,C,ferr,berr      !  routine xGESVX
    real(KND),dimension(12) :: work
    real(KND),dimension(3,3) :: ag,af                    !
    integer   :: info                                   !
    real(KND) :: rcond                                  !
    character :: equed                                  !

   !Vzdalenosti od rovin, pokud je nejbl. bod roviny uvnitr jine, nebo puv. bod na vnitrni strane -> vzd. *-1
   !Pokud je nejbl. body +- eps. (norm. vektor!,dxmin/100) uvnitr a vne polyh -> hotovo
   !Jinak 2. nejbl. rovina v abs. hodnote -> prusecnice a nejbl bod na ni
   !Nejbl. bod na prusecnici. Pokud +-eps. uvnitr, (najit vekt. v rovine  kolme na prusecnic)-:hotovo
   !Jinak iterativne najit bod na prusecnici uvnitr

    do i = 1,PH%nplanes
     call PlaneNearest(PH%Planes(i),xP(i),yP(i),zP(i),x,y,z)
     dists(i) = SQRT((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
     if (PlaneInside(PH%Planes(i),x,y,z)) dists(i) = -ABS(dists(i))
    enddo
    !find nearest plane with
    nearest = 0
    minv = huge(minv)
    do i = 1,PH%nplanes
     if (dists(i)>=0.and.dists(i)<minv) then
       nearest = i
       minv = dists(i)
     endif
    enddo

    if (nearest==0) then
      write(*,*) "no nearest"
      return
    endif

    if (PolyhedronInside(PH,xP(nearest),yP(nearest),zP(nearest),MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND))) then
     xnear = xP(nearest)
     ynear = yP(nearest)
     znear = zP(nearest)
     return
    endif

    dists = abs(dists)

    nearest2 = 0
    minv = huge(minv)
    do i = 1,PH%nplanes
     if (i/=nearest.and.dists(i)<minv) then
       nearest2 = i
       minv = dists(i)
     endif
    enddo

    nearest3 = 0
    minv = huge(minv)
    do i = 1,PH%nplanes
     if (i/=nearest.and.i/=nearest2.and.dists(i)<minv) then
       nearest3 = i
       minv = dists(i)
     endif
    enddo

    ailine = PH%Planes(nearest)%b*PH%Planes(nearest2)%c-PH%Planes(nearest)%c*PH%Planes(nearest2)%b
    biline = PH%Planes(nearest)%c*PH%Planes(nearest2)%a-PH%Planes(nearest)%a*PH%Planes(nearest2)%c
    ciline = PH%Planes(nearest)%a*PH%Planes(nearest2)%b-PH%Planes(nearest)%b*PH%Planes(nearest2)%a

    if (abs(ailine)<=epsilon(ailine).and.abs(biline)<=epsilon(biline).and.abs(ciline)<=epsilon(ciline)) then
     write(*,*) "cross product 0"
     write(*,*) "numbers of planes",nearest,nearest2
     write(*,*) "x,y,z",x,y,z
     write(*,*) "pl1",PH%Planes(nearest)%a,PH%Planes(nearest)%b,PH%Planes(nearest)%c,PH%Planes(nearest)%d
     write(*,*) "pl2",PH%Planes(nearest2)%a,PH%Planes(nearest2)%b,PH%Planes(nearest2)%c,PH%Planes(nearest2)%d
    endif

    p = SQRT(ailine**2+biline**2+ciline**2)
    ailine = ailine/p
    biline = biline/p
    ciline = ciline/p

    if (abs(ciline)>=0.1_KND) then

       ag = reshape(source = (/ PH%Planes(nearest)%a,PH%Planes(nearest2)%a,0._KND,&
                           PH%Planes(nearest)%b,PH%Planes(nearest2)%b,0._KND,&
                           PH%Planes(nearest)%c,PH%Planes(nearest2)%c,1._KND /),&
                 shape = (/ 3,3 /))
 !       ag = reshape(source = (/ PH%Planes(nearest)%a,PH%Planes(nearest)%b,PH%Planes(nearest)%c,&
 !                           PH%Planes(nearest2)%a,PH%Planes(nearest2)%b,PH%Planes(nearest2)%c,&
 !                           0._KND,0._KND,1._KND /),&
 !                 shape = (/ 3,3 /))
       bg = (/ -PH%Planes(nearest)%d,-PH%Planes(nearest2)%d,(zW(Wnz+1)+zW(0))/2._KND /)

       if (KND==DBL) then
        call DGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       else
        call SGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       endif

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    elseif (abs(biline)>=0.1_KND) then

       ag = reshape(source = (/ PH%Planes(nearest)%a,PH%Planes(nearest2)%a,0._KND,&
                           PH%Planes(nearest)%b,PH%Planes(nearest2)%b,1._KND,&
                           PH%Planes(nearest)%c,PH%Planes(nearest2)%c,0._KND /),&
                 shape = (/ 3,3 /))
       bg = (/ -PH%Planes(nearest)%d,-PH%Planes(nearest2)%d,(yV(Vny+1)+yV(0))/2._KND /)

       if (KND==DBL) then
        call DGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       else
        call SGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       endif

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    else

       ag = reshape(source = (/ PH%Planes(nearest)%a,PH%Planes(nearest2)%a,1._KND,&
                           PH%Planes(nearest)%b,PH%Planes(nearest2)%b,0._KND,&
                           PH%Planes(nearest)%c,PH%Planes(nearest2)%c,0._KND /),&
                  shape = (/ 3,3 /))
       bg = (/ -PH%Planes(nearest)%d,-PH%Planes(nearest2)%d,(xU(Unx+1)+xU(0))/2._KND /)

       if (KND==DBL) then
        call DGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       else
        call SGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
       endif

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    endif

    call LineNearest(xln,yln,zln,x,y,z,x0iline,y0iline,z0iline,ailine,biline,ciline)

    if (PolyhedronInside(PH,xln,yln,zln,min(dxmin/1000._KND,dymin/10000._KND,dzmin/10000._KND))) then
     xnear = xln
     ynear = yln
     znear = zln
     return
    endif



    ag = reshape(source = (/ PH%Planes(nearest)%a,PH%Planes(nearest2)%a,PH%Planes(nearest3)%a,&
                        PH%Planes(nearest)%b,PH%Planes(nearest2)%b,PH%Planes(nearest3)%b,&
                        PH%Planes(nearest)%c,PH%Planes(nearest2)%c,PH%Planes(nearest3)%c /),&
               shape = (/ 3,3 /))
    bg = (/ -PH%Planes(nearest)%d,-PH%Planes(nearest2)%d,-PH%Planes(nearest3)%d /)

    if (KND==DBL) then
     call DGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
    else
     call SGESVX("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3,rcond,ferr,berr,work,work2,info)
    endif

    xnear = xg(1)
    ynear = xg(2)
    znear = xg(3)

  end subroutine PolyhedronNearest




  subroutine TerrainNearest(T,xnear,ynear,znear,x,y,z)
    type(TPlane) :: Pl
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    type(TTerrain),intent(in) :: T
    integer     :: xi,yj,comp
    real(KND)   :: a,b,zloc
    xnear = x
    ynear = y
    znear = z0

    call TerrGridCoords(x,y,xi,yj,comp)

    if (comp==1) then  !Construct a tangent plane using first derivatives

     zloc = T%UPoints(xi,yj)%elev
     a = (T%PrPoints(xi+1,yj)%elev - T%PrPoints(xi,yj)%elev) / dxU(xi)
     b = (T%UPoints(xi,yj+1)%elev - T%UPoints(xi,yj-1)%elev) / (yPr(yj+1)-yPr(yj-1))

     Pl%a = a
     Pl%b = b
     Pl%c = -1._KND
     Pl%d = -a*x-b*y+zloc

     call PlaneNearest(Pl,xnear,ynear,znear,x,y,z)

    elseif (comp==2) then

     zloc = T%VPoints(xi,yj)%elev
     a = (T%VPoints(xi+1,yj)%elev - T%VPoints(xi-1,yj)%elev) / (xPr(xi+1)-xPr(xi-1))
     b = (T%PrPoints(xi,yj+1)%elev - T%PrPoints(xi,yj)%elev) / dyV(yj)

     Pl%a = a
     Pl%b = b
     Pl%c = -1._KND
     Pl%d = -a*x-b*y+zloc

     call PlaneNearest(Pl,xnear,ynear,znear,x,y,z)

    elseif (comp==3) then

     zloc = T%PrPoints(xi,yj)%elev
     a = (T%UPoints(xi,yj)%elev - T%UPoints(xi-1,yj)%elev) / dxPr(xi)
     b = (T%VPoints(xi,yj)%elev - T%VPoints(xi,yj-1)%elev) / dyPr(yj)

     Pl%a = a
     Pl%b = b
     Pl%c = -1._KND
     Pl%d = -a*x-b*y+zloc

     call PlaneNearest(Pl,xnear,ynear,znear,x,y,z)

    endif
  end subroutine TerrainNearest


  subroutine SolidBodyNearest(SB,xnear,ynear,znear,x,y,z)
    type(TSolidBody),intent(in) :: SB
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    select case (SB%typeofbody)
     case (1)
      call PolyhedronNearest(SB%Polyhedron,xnear,ynear,znear,x,y,z)
     case (2)
      call BallNearest(SB%Ball,xnear,ynear,znear,x,y,z)
     case (3)
      call CylinderNearest(SB%Cylinder,xnear,ynear,znear,x,y,z)
     case (4)
      call TerrainNearest(SB%Terrain,xnear,ynear,znear,x,y,z)
     case default
      xnear = -huge(znear)
      ynear = -huge(znear)
      znear = -huge(znear)
    end select
  end subroutine SolidBodyNearest



  subroutine CylinderNearestOut(C,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
    type(TCylinder),intent(in) :: C
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    real(KND) xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2

    if (associated(C%Plane1)) then
      call PlaneNearest(C%Plane1,xP1,yP1,zP1,x,y,z)
    else
     xP1 = sqrt(huge(xP1))/10
     yP1 = sqrt(huge(yP1))/10
     zP1 = sqrt(huge(zP1))/10
    endif
    if (associated(C%Plane2)) then
      call PlaneNearest(C%Plane2,xP2,yP2,zP2,x,y,z)
    else
     xP2 = sqrt(huge(xP2))/10
     yP2 = sqrt(huge(yP2))/10
     zP2 = sqrt(huge(zP2))/10
    endif
    call JacketNearest(C%Jacket,xJ,yJ,zJ,x,y,z)

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
  end subroutine CylinderNearestOut


  subroutine PolyhedronNearestOut(PH,xnear,ynear,znear,x,y,z)
    type(TPolyhedron),intent(in) :: PH
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    real(KND) dists(PH%nplanes),xP(PH%nplanes),yP(PH%nplanes),zP(PH%nplanes),minv
    integer nearest,i

    dists = huge(minv)
    do i = 1,PH%nplanes
     call PlaneNearest(PH%Planes(i),xP(i),yP(i),zP(i),x,y,z)
     dists(i) = SQRT((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
    enddo

    nearest = 0
    minv = huge(minv)
    do i = 1,PH%nplanes
     if (dists(i)>=0.and.dists(i)<minv) then
       nearest = i
       minv = dists(i)
     endif
    enddo

    xnear = xP(nearest)
    ynear = yP(nearest)
    znear = zP(nearest)
  end subroutine PolyhedronNearestOut





  subroutine SolidBodyNearestOut(SB,xnear,ynear,znear,x,y,z)
    type(TSolidBody),intent(in) :: SB
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    xnear = -1e+9;ynear = -1e+9;znear = -1e+9

    select case (SB%typeofbody)
     case (1)
      call PolyhedronNearestOut(SB%Polyhedron,xnear,ynear,znear,x,y,z)
     case (2)
      call BallNearest(SB%Ball,xnear,ynear,znear,x,y,z)
     case (3)
      call CylinderNearestOut(SB%Cylinder,xnear,ynear,znear,x,y,z)
     case (4)
      call TerrainNearest(SB%Terrain,xnear,ynear,znear,x,y,z)
     case default
      xnear = -huge(znear);ynear = -huge(znear);znear = -huge(znear)
    end select
  end subroutine SolidBodyNearestOut






  subroutine GetSolidBodiesBC
    use WallModels, only: WMPoint, AddWMPoint
    integer i,j,k,m,n,o
    integer(SINT) nb
    real(KND) dist,nearx,neary,nearz
    type(TSolidBody),pointer :: CurrentSB => null()
    type(WMPOINT) :: WMP

    !find if the gridpoints lie inside a solid body and write it's number
    !do not nullify the .type arrays, they could have been made nonzero by other unit
    CurrentSB => FirstSB

    do
     if (associated(CurrentSB)) then
      !$omp parallel do
      do k = 0,Prnz+1
       do j = 0,Prny+1
        do i = 0,Prnx+1
           if (SolidBodyInside(CurrentSB,xPr(i),yPr(j),zPr(k))) Prtype(i,j,k) = CurrentSB%numofbody
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
      !$omp parallel do
      do k = 0,Unz+1
       do j = 0,Uny+1
         do i = 0,Unx+1
            if (SolidBodyInside(CurrentSB,xU(i),yPr(j),zPr(k),max(dxU(i),dyPr(j),dzPr(k))/5._KND))&
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
      !$omp parallel do
      do k = 0,Vnz+1
       do j = 0,Vny+1
        do i = 0,Vnx+1
           if (SolidBodyInside(CurrentSB,xPr(i),yV(j),zPr(k),max(dxPr(i),dyV(j),dzPr(k))/5._KND))&
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
      !$omp parallel do
      do k = 0,Wnz+1
       do j = 0,Wny+1
        do i = 0,Wnx+1
           if (SolidBodyInside(CurrentSB,xPr(i),yPr(j),zW(k),max(dxPr(i),dyPr(j),dzW(k))/5._KND))&
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


    call InitImBoundaries

    allocate(WMP%depscalar(computescalars))

    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       if (Prtype(i,j,k)==0) then
        dist = huge(dist)
        nb = 0
        do o = -1,1
         do n = -1,1
          do m = -1,1
           if ((m/=0.or.n/=0.or.o/=0).and.i+m>0.and.j+n>0.and.k+o>0.and.i+m<=Prnx.and.j+n<=Prny.and.k+o<=Prnz.and.&
            (Prtype(i+m,j+n,k+o)>0).and.Prtype(i+m,j+n,k+o)/=nb) then
             call SetCurrentSB(CurrentSB,Prtype(i+m,j+n,k+o))
             call SolidBodyNearest(CurrentSB,nearx,neary,nearz,xPr(i),yPr(j),zPr(k))
             if (SQRT((nearx-xPr(i))**2+(neary-yPr(j))**2+(nearz-zPr(k))**2)<dist) then
              dist = SQRT((nearx-xPr(i))**2+(neary-yPr(j))**2+(nearz-zPr(k))**2)
              nb = Prtype(i+m,j+n,k+o)
             endif
           endif
          enddo
         enddo
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

  end subroutine GetSolidBodiesBC



  ! File customBC.f90 should contain user generated subroutine InitSolidBodies that inicializes
  ! the linked list of solid bodies objects. This is to avoid unnecessary changes in the history
  ! of this file, when only boundary conditions changes, until a runtime BC generation is developed.
  subroutine InitSolidBodies
    include "customBC.f90"
  end subroutine InitSolidBodies




  subroutine AddIBpoint(IBP)
    type(TIBPoint),intent(in) :: IBP

    if (.not.associated(LastIBPoint)) then
     allocate(FirstIBPoint)
     FirstIBPoint = IBP
     LastIBPoint => FirstIBPoint
    else
     allocate(LastIBPoint%next)
     LastIBPoint%next = IBP
     LastIBPoint => LastIBPoint%next
    endif
  end subroutine AddIBPoint

  subroutine AddScalIBpoint(SIBP)
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
  end subroutine AddScalIBPoint


  real(KND) function NearestOnLineOut(SB,x,y,z,x2,y2,z2) !Find t, such that x+(x2-x)*t lies on the boundary of the SB
    type(TSolidBody),intent(in) :: SB
    real(KND),intent(in) :: x,y,z,x2,y2,z2

    real(KND) t,t1,t2
    integer i

    t1 = 0
    t2 = 1
    if (SolidBodyInside(SB,x2,y2,z2)) then
     do
      t2 = t2*2._KND
      if (.not.SolidBodyInside(SB,x+(x2-x)*t2,y+(y2-y)*t2,z+(z2-z)*t2)) exit
     enddo
    endif
    t = (t1+t2)/2._KND

    do i = 1,20         !The bisection method with maximum 20 iterations (should be well enough)
     if (SolidBodyInside(SB,x+(x2-x)*t,y+(y2-y)*t,z+(z2-z)*t)) then
      t1 = t
     else
      t2 = t
     endif
     t = (t1+t2)/2._KND
     if (abs(t1-t2)<MIN(dxmin/1000._KND,dymin/1000._KND,dzmin/1000._KND))   exit
    enddo
    NearestOnLineOut = t
  endfunction NearestOnLineOut



  subroutine GetUIBPoint(IBP,xi,yj,zk)
    type(TIBPoint),intent(out) :: IBP
    integer,intent(in) :: xi,yj,zk

    type(TSolidBody),pointer :: SB
    integer dirx,diry,dirz,n1,n2,nx,ny,nz
    real(KND) x,y,z,xnear,ynear,znear,t
    logical free100,free010,free001

    x = xU(xi)                                !real coordinates of the IB forcing point
    y = yPr(yj)
    z = zPr(zk)
    call SetCurrentSB(SB,Utype(xi,yj,zk))
    call SolidBodyNearestOut(SB,xnear,ynear,znear,x,y,z)
    IBP%component = 1
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

    if (.not.SolidBodyInside(SB,xU(xi),yPr(yj),zPr(zk)))  then  !For now, if actually outside the body, set an artificial boundary
     IBP%interp = 0                                            !point here. In future we can use another interpolation.
     IBP%interpdir = 0
     return
    endif


    if (abs(IBP%distx)<(xU(xi+1)-xU(xi-1))/100._KND) then      !if too close to the boundary, set the distance to 0
     IBP%distx = 0
     dirx = 0
     IBP%dirx = 0
    endif
    if (abs(IBP%disty)<(yPr(yj+1)-yPr(yj-1))/100._KND) then
     IBP%disty = 0
     diry = 0
     IBP%diry = 0
    endif
    if (abs(IBP%distz)<(zPr(zk+1)-zPr(zk-1))/100._KND) then
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
    if (Utype(xi+1,yj,zk)==0) then
     n1 = n1+1
     nx = nx+1
    endif
    if (Utype(xi-1,yj,zk)==0) then
     n1 = n1+1
     nx = nx+1
    endif
    if (Utype(xi,yj+1,zk)==0) then
     n1 = n1+1
     ny = ny+1
    endif
    if (Utype(xi,yj-1,zk)==0) then
     n1 = n1+1
     ny = ny+1
    endif
    if (Utype(xi,yj,zk+1)==0) then
     n1 = n1+1
     nz = nz+1
    endif
    if (Utype(xi,yj,zk-1)==0) then
     n1 = n1+1
     nz = nz+1
    endif

    if (dirx/=0) n2 = n2+1
    if (diry/=0) n2 = n2+1
    if (dirz/=0) n2 = n2+1


    if (n1>n2) then   ! If too many free directions, treat as directly on the boundary,
     IBP%interp = 0     !because we are probably at some edge.
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
     IBP%interp = 0          !If no free direction, the boundary is here.
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
     IBP%interp = 3                                                                           !use trilinear interpolation.
     IBP%interpdir = 0
    elseif ((.not.free100).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then !Two free directions y and z,
     IBP%interp = 2                                                                    !use bilinear interpolation normal to x.
     IBP%interpdir = 1
     if (dirx==1) then                                       !If dirx /= 0 then forget the x component of dist vector
      t = NearestOnLineOut(SB,x,y,z,x,y+IBP%disty,z+IBP%distz) ! and find an intersection of the new vector with the boundary

      dirx = 0
      IBP%dirx = 0
      IBP%distx = 0
      IBP%disty = IBP%disty*t
      IBP%distz = IBP%distz*t
     endif
    elseif ((.not.free010).and.(free100.and.dirx==1).and.(free001.and.dirz==1)) then  !the same normal to y
     IBP%interp = 2
     IBP%interpdir = 2
     if (diry==1) then
      t = NearestOnLineOut(SB,x,y,z,x+IBP%distx,y,z+IBP%distz)

      diry = 0
      IBP%diry = 0
      IBP%distx = IBP%distx*t
      IBP%disty = 0
      IBP%distz = IBP%distz*t
     endif
    elseif ((.not.free001).and.(free100.and.dirx==1).and.(free010.and.diry==1)) then  !the same normal to z
     IBP%interp = 2
     IBP%interpdir = 3
     if (dirz==1) then
      t = NearestOnLineOut(SB,x,y,z,x+IBP%distx,y+IBP%disty,z)

      dirz = 0
      IBP%dirz = 0
      IBP%distx = IBP%distx*t
      IBP%disty = IBP%disty*t
      IBP%distz = 0
     endif
    elseif (free100.and.dirx==1) then  !Only one free direction, use linear interpolation in direction x.
     IBP%interp = 1
     IBP%interpdir = 1
     if (diry==1.or.dirz==1) then                   !If other dir components nonzero, delete them and
      t = NearestOnLineOut(SB,x,y,z,x+IBP%distx,y,z)  !find an intersection of the new vector with the boundary

      diry = 0
      dirz = 0
      IBP%diry = 0
      IBP%dirz = 0
      IBP%distx = t*IBP%distx
      IBP%disty = 0
      IBP%distz = 0
     endif
    elseif (free010.and.diry==1) then               !the same in y
     IBP%interp = 1
     IBP%interpdir = 2
     if (dirx==1.or.dirz==1) then
      t = NearestOnLineOut(SB,x,y,z,x,y+IBP%disty,z)

      dirx = 0
      dirz = 0
      IBP%dirx = 0
      IBP%dirz = 0
      IBP%distx = 0
      IBP%disty = t*IBP%disty
      IBP%distz = 0
     endif
    elseif (free001.and.dirz==1) then               !the same in z
     IBP%interp = 1
     IBP%interpdir = 3
     if (dirx==1.or.diry==1) then
      t = NearestOnLineOut(SB,x,y,z,x,y,z+IBP%distz)

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
  end subroutine GetUIBPoint


  subroutine GetVIBPoint(IBP,xi,yj,zk)
    type(TIBPoint),intent(out) :: IBP
    integer,intent(in) :: xi,yj,zk

    type(TSolidBody),pointer :: SB
    integer dirx,diry,dirz,n1,n2,nx,ny,nz
    real(KND) x,y,z,xnear,ynear,znear,t
    logical free100,free010,free001


    x = xPr(xi)
    y = yV(yj)
    z = zPr(zk)
    call SetCurrentSB(SB,Vtype(xi,yj,zk))
    call SolidBodyNearestOut(SB,xnear,ynear,znear,x,y,z)

    IBP%component = 2
    IBP%xi = xi
    IBP%yj = yj
    IBP%zk = zk
    IBP%distx = xnear-x
    IBP%disty = ynear-y
    IBP%distz = znear-z
    IBP%dirx = nint(sign(1.0_KND,IBP%distx))
    IBP%diry = nint(sign(1.0_KND,IBP%disty))
    IBP%dirz = nint(sign(1.0_KND,IBP%distz))

    dirx = abs(IBP%dirx)
    diry = abs(IBP%diry)
    dirz = abs(IBP%dirz)

    if (.not.SolidBodyInside(SB,xPr(xi),yV(yj),zPr(zk)))  then  !For now, if actually outside the body, set an artificial boundary
     IBP%interp = 0                                            !point here. In future we can use another interpolation.
     IBP%interpdir = 0
     return
    endif


    if (abs(IBP%distx)<(xPr(xi+1)-xPr(xi-1))/100._KND) then
     IBP%distx = 0
     dirx = 0
     IBP%dirx = 0
    endif
    if (abs(IBP%disty)<(yV(yj+1)-yV(yj-1))/100._KND) then
     IBP%disty = 0
     diry = 0
     IBP%diry = 0
    endif
    if (abs(IBP%distz)<(zPr(zk+1)-zPr(zk-1))/100._KND) then
     IBP%distz = 0
     dirz = 0
     IBP%dirz = 0
    endif

    if (dirx==0) then
     free100=.false.
    else
     free100 = (Vtype(xi+IBP%dirx,yj,zk)==0)
    endif
    if (diry==0) then
     free010=.false.
    else
     free010 = (Vtype(xi,yj+IBP%diry,zk)==0)
    endif
    if (dirz==0) then
     free001=.false.
    else
     free001 = (Vtype(xi,yj,zk+IBP%dirz)==0)
    endif


    nx = 0
    ny = 0
    nz = 0
    n1 = 0                                                          ! n1 number of neighbouring cells free of solid bodies
    n2 = 0                                                          ! n2 number of nonzero coordinate directions to boundary
    if (Vtype(xi+1,yj,zk)==0) then
     n1 = n1+1
     nx = nx+1
    endif
    if (Vtype(xi-1,yj,zk)==0) then
     n1 = n1+1
     nx = nx+1
    endif
    if (Vtype(xi,yj+1,zk)==0) then
     n1 = n1+1
     ny = ny+1
    endif
    if (Vtype(xi,yj-1,zk)==0) then
     n1 = n1+1
     ny = ny+1
    endif
    if (Vtype(xi,yj,zk+1)==0) then
     n1 = n1+1
     nz = nz+1
    endif
    if (Vtype(xi,yj,zk-1)==0) then
     n1 = n1+1
     nz = nz+1
    endif

    if (dirx/=0) n2 = n2+1
    if (diry/=0) n2 = n2+1
    if (dirz/=0) n2 = n2+1

    if (n1>n2) then
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

    if (nx>1.or.ny>1.or.nz>1) then
     IBP%interp = 0
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
    elseif ((free100.and.dirx==1).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then
     IBP%interp = 3
     IBP%interpdir = 0
    elseif ((.not.free100).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then
     IBP%interp = 2
     IBP%interpdir = 1
     if (dirx==1) then
      t = NearestOnLineOut(SB,x,y,z,x,y+IBP%disty,z+IBP%distz)

      dirx = 0
      IBP%dirx = 0
      IBP%distx = 0
      IBP%disty = IBP%disty*t
      IBP%distz = IBP%distz*t
     endif
    elseif ((.not.free010).and.(free100.and.dirx==1).and.(free001.and.dirz==1)) then
     IBP%interp = 2
     IBP%interpdir = 2
     if (diry==1) then
      t = NearestOnLineOut(SB,x,y,z,x+IBP%distx,y,z+IBP%distz)

      diry = 0
      IBP%diry = 0
      IBP%distx = IBP%distx*t
      IBP%disty = 0
      IBP%distz = IBP%distz*t
     endif
    elseif ((.not.free001).and.(free100.and.dirx==1).and.(free010.and.diry==1)) then
     IBP%interp = 2
     IBP%interpdir = 3
     if (dirz==1) then
      t = NearestOnLineOut(SB,x,y,z,x+IBP%distx,y+IBP%disty,z)

      dirz = 0
      IBP%dirz = 0
      IBP%distx = IBP%distx*t
      IBP%disty = IBP%disty*t
      IBP%distz = 0
     endif
    elseif (free100.and.dirx==1) then
     IBP%interp = 1
     IBP%interpdir = 1
     if (diry==1.or.dirz==1) then
      t = NearestOnLineOut(SB,x,y,z,x+IBP%distx,y,z)

      diry = 0
      dirz = 0
      IBP%diry = 0
      IBP%dirz = 0
      IBP%distx = t*IBP%distx
      IBP%disty = 0
      IBP%distz = 0
     endif
    elseif (free010.and.diry==1) then
     IBP%interp = 1
     IBP%interpdir = 2
     if (dirx==1.or.dirz==1) then
      t = NearestOnLineOut(SB,x,y,z,x,y+IBP%disty,z)

      dirx = 0
      dirz = 0
      IBP%dirx = 0
      IBP%dirz = 0
      IBP%distx = 0
      IBP%disty = t*IBP%disty
      IBP%distz = 0
     endif
    elseif (free001.and.dirz==1) then
     IBP%interp = 1
     IBP%interpdir = 3
     if (diry==1.or.dirx==1) then
      t = NearestOnLineOut(SB,x,y,z,x,y,z+IBP%distz)

      dirx = 0
      diry = 0
      IBP%dirx = 0
      IBP%diry = 0
      IBP%distx = 0
      IBP%disty = 0
      IBP%distz = t*IBP%distz
     endif
    else
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
     write(*,*) "x",x
     write(*,*) "y",y
     write(*,*) "z",z
     write(*,*) "IBP%component",IBP%component
     stop
    endif

  end subroutine GetVIBPoint


  subroutine GetWIBPoint(IBP,xi,yj,zk)
    type(TIBPoint),intent(out) :: IBP
    integer,intent(in) :: xi,yj,zk

    type(TSolidBody),pointer :: SB
    integer dirx,diry,dirz,n1,n2,nx,ny,nz
    real(KND) x,y,z,xnear,ynear,znear,t
    logical free100,free010,free001


    x = xPr(xi)
    y = yPr(yj)
    z = zW(zk)
    call SetCurrentSB(SB,Wtype(xi,yj,zk))
    call SolidBodyNearestOut(SB,xnear,ynear,znear,x,y,z)

    IBP%component = 3
    IBP%xi = xi
    IBP%yj = yj
    IBP%zk = zk
    IBP%distx = xnear-x
    IBP%disty = ynear-y
    IBP%distz = znear-z
    IBP%dirx = nint(sign(1.0_KND,IBP%distx))
    IBP%diry = nint(sign(1.0_KND,IBP%disty))
    IBP%dirz = nint(sign(1.0_KND,IBP%distz))

    dirx = abs(IBP%dirx)
    diry = abs(IBP%diry)
    dirz = abs(IBP%dirz)

    if (.not.SolidBodyInside(SB,xPr(xi),yPr(yj),zW(zk)))  then  !For now, if actually outside the body, set an artificial boundary
     IBP%interp = 0                                            !point here. In future we can use another interpolation.
     IBP%interpdir = 0
     return
    endif


    if (abs(IBP%distx)<(xPr(xi+1)-xPr(xi-1))/100._KND) then
     IBP%distx = 0
     dirx = 0
     IBP%dirx = 0
    endif
    if (abs(IBP%disty)<(yPr(yj+1)-yPr(yj-1))/100._KND) then
     IBP%disty = 0
     diry = 0
     IBP%diry = 0
    endif
    if (abs(IBP%distz)<(zW(zk+1)-zW(zk-1))/100._KND) then
     IBP%distz = 0
     dirz = 0
     IBP%dirz = 0
    endif


    if (dirx==0) then
     free100=.false.
    else
     free100 = (Wtype(xi+IBP%dirx,yj,zk)==0)
    endif
    if (diry==0) then
     free010=.false.
    else
     free010 = (Wtype(xi,yj+IBP%diry,zk)==0)
    endif
    if (dirz==0) then
     free001=.false.
    else
     free001 = (Wtype(xi,yj,zk+IBP%dirz)==0)
    endif


    nx = 0
    ny = 0
    nz = 0
    n1 = 0                                                          ! n1 number of neighbouring cells free of solid bodies
    n2 = 0                                                          ! n2 number of nonzero coordinate directions to boundary
    if (Wtype(xi+1,yj,zk)==0) then
     n1 = n1+1
     nx = nx+1
    endif
    if (Wtype(xi-1,yj,zk)==0) then
     n1 = n1+1
     nx = nx+1
    endif
    if (Wtype(xi,yj+1,zk)==0) then
     n1 = n1+1
     ny = ny+1
    endif
    if (Wtype(xi,yj-1,zk)==0) then
     n1 = n1+1
     ny = ny+1
    endif
    if (Wtype(xi,yj,zk+1)==0) then
     n1 = n1+1
     nz = nz+1
    endif
    if (Wtype(xi,yj,zk-1)==0) then
     n1 = n1+1
     nz = nz+1
    endif

    if (dirx/=0) n2 = n2+1
    if (diry/=0) n2 = n2+1
    if (dirz/=0) n2 = n2+1

    if (n1>n2) then
     IBP%interp = 0
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
    elseif ((free100.and.dirx==1).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then
     IBP%interp = 3
     IBP%interpdir = 0
    elseif ((.not.free100).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then
     IBP%interp = 2
     IBP%interpdir = 1
     if (dirx==1) then
      t = NearestOnLineOut(SB,x,y,z,x,y+IBP%disty,z+IBP%distz)

      dirx = 0
      IBP%dirx = 0
      IBP%distx = 0
      IBP%disty = IBP%disty*t
      IBP%distz = IBP%distz*t
     endif
    elseif ((.not.free010).and.(free100.and.dirx==1).and.(free001.and.dirz==1)) then
     IBP%interp = 2
     IBP%interpdir = 2
     if (diry==1) then
      t = NearestOnLineOut(SB,x,y,z,x+IBP%distx,y,z+IBP%distz)

      diry = 0
      IBP%diry = 0
      IBP%distx = IBP%distx*t
      IBP%disty = 0
      IBP%distz = IBP%distz*t
     endif
    elseif ((.not.free001).and.(free100.and.dirx==1).and.(free010.and.diry==1)) then
     IBP%interp = 2
     IBP%interpdir = 3
     if (dirz==1) then
      t = NearestOnLineOut(SB,x,y,z,x+IBP%distx,y+IBP%disty,z)

      dirz = 0
      IBP%dirz = 0
      IBP%distx = IBP%distx*t
      IBP%disty = IBP%disty*t
      IBP%distz = 0
     endif
    elseif (free100.and.dirx==1) then
     IBP%interp = 1
     IBP%interpdir = 1
     if (diry==1.or.dirz==1) then
      t = NearestOnLineOut(SB,x,y,z,x+IBP%distx,y,z)

      diry = 0
      dirz = 0
      IBP%diry = 0
      IBP%dirz = 0
      IBP%distx = t*IBP%distx
      IBP%disty = 0
      IBP%distz = 0
     endif
    elseif (free010.and.diry==1) then
     IBP%interp = 1
     IBP%interpdir = 2
     if (dirx==1.or.dirz==1) then
      t = NearestOnLineOut(SB,x,y,z,x,y+IBP%disty,z)

      dirx = 0
      dirz = 0
      IBP%dirx = 0
      IBP%dirz = 0
      IBP%distx = 0
      IBP%disty = t*IBP%disty
      IBP%distz = 0
     endif
    elseif (free001.and.dirz==1) then
     IBP%interp = 1
     IBP%interpdir = 3
     if (dirx==1.or.diry==1) then
      t = NearestOnLineOut(SB,x,y,z,x,y,z+IBP%distz)

      dirx = 0
      diry = 0
      IBP%dirx = 0
      IBP%diry = 0
      IBP%distx = 0
      IBP%disty = 0
      IBP%distz = t*IBP%distz
     endif
    else
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
     write(*,*) "x",x
     write(*,*) "y",y
     write(*,*) "z",z
     write(*,*) "IBP%component",IBP%component
     stop
    endif
  end subroutine GetWIBPoint


  subroutine GeTScalIFlBPoint(IBP,xi,yj,zk)
    type(TScalFlIBPoint),intent(out) :: IBP
    integer,intent(in) :: xi,yj,zk               !grid coordinates of the forcing point

    type(TSolidBody),pointer :: SB
    integer dirx,diry,dirz,dirx2,diry2,dirz2,n1,n2
    real(KND) x,y,z,xnear,ynear,znear,distx,disty,distz,t,tx,ty,tz
    logical freep00,free0p0,free00p,freem00,free0m0,free00m


    x = xPr(xi)                                   !physical coordinates of the forcing point
    y = yPr(yj)
    z = zPr(zk)
    call SetCurrentSB(SB,Prtype(xi,yj,zk))
    call SolidBodyNearestOut(SB,xnear,ynear,znear,x,y,z)

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

    freep00 = (Prtype(xi+1,yj,zk)==0)   !logicals denoting if the cell in plus x direction is free of SB
    freem00 = (Prtype(xi-1,yj,zk)==0)
    free0p0 = (Prtype(xi,yj+1,zk)==0)
    free0m0 = (Prtype(xi,yj-1,zk)==0)
    free00p = (Prtype(xi,yj,zk+1)==0)
    free00m = (Prtype(xi,yj,zk-1)==0)

    n1 = 0
    n2 = 0
    if (freep00) n1 = n1+1
    if (free0p0) n1 = n1+1
    if (free00p) n1 = n1+1
    if (freem00) n1 = n1+1
    if (free0m0) n1 = n1+1
    if (free00m) n1 = n1+1

    if (dirx/=0) n2 = n2+1
    if (diry/=0) n2 = n2+1
    if (dirz/=0) n2 = n2+1

    if (n1>n2.or.(freep00.and.freem00).or.(free0p0.and.free0m0).or.(free00p.and.free00m)) then
     IBP%interp = 0    !If moore free spaces than directions, or free space in both oposite directions,
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
      IBP%intpointi(1) = IBP%xi+dirx2
      IBP%intpointj(1) = IBP%yj+diry2
      IBP%intpointk(1) = IBP%zk+dirz2
      IBP%intcoef = 1._KND
      IBP%interp = 1
      IBP%dist = sqrt((x-xPr(IBP%xi+dirx2))**2+(y-yPr(IBP%yj+diry2))**2+(z-zPr(IBP%zk+dirz2))**2)

    elseif (n2==1) then     !Only one free point outside, interpolate from there.

      IBP%intpointi(1) = IBP%xi+dirx
      IBP%intpointj(1) = IBP%yj+diry
      IBP%intpointk(1) = IBP%zk+dirz
      IBP%intcoef = 1._KND
      IBP%interp = 1
      IBP%dist = sqrt((x-xPr(IBP%xi+dirx))**2+(y-yPr(IBP%yj+diry))**2+(z-zPr(IBP%zk+dirz))**2)

    elseif (n2==2) then    !Two free directions, use bilinear interpolation in the plane contaning these two neigbours.

     if (dirx==0) then     !plane yz

      if (abs(disty/distz)>abs(yPr(yj+diry)-y)/abs(zPr(zk+dirz)-z)) then  !Which gridline dous the line from this point
       t = (zPr(zk+dirz)-z)/distz                                           !to the boundary intersect? Then use the points
       IBP%intpointi(1) = IBP%xi                                             !on both end s of the grid line.
       IBP%intpointj(1) = IBP%yj+diry
       IBP%intpointk(1) = IBP%zk+dirz
       IBP%intcoef(1) = abs(disty*t)/abs(yPr(IBP%yj+diry)-yPr(IBP%yj))
       IBP%intpointi(2) = IBP%xi
       IBP%intpointj(2) = IBP%yj
       IBP%intpointk(2) = IBP%zk+dirz
       IBP%intcoef(2) = 1-IBP%intcoef(1)
       IBP%interp = 2
       IBP%dist = sqrt((disty*t)**2+(distz*t)**2)
      else
       t = (yPr(yj+diry)-y)/disty
       IBP%intpointi(1) = IBP%xi
       IBP%intpointj(1) = IBP%yj+diry
       IBP%intpointk(1) = IBP%zk+dirz
       IBP%intcoef(1) = abs(distz*t)/abs(zPr(IBP%zk+dirz)-zPr(IBP%zk))
       IBP%intpointi(2) = IBP%xi
       IBP%intpointj(2) = IBP%yj+diry
       IBP%intpointk(2) = IBP%zk
       IBP%intcoef(2) = 1-IBP%intcoef(1)
       IBP%interp = 2
       IBP%dist = sqrt((disty*t)**2+(distz*t)**2)
      endif

     elseif (diry==0) then     !plane xz

      if (abs(distx/distz)>abs(xPr(xi+dirx)-x)/abs(zPr(zk+dirz)-z)) then
       t = (zPr(zk+dirz)-z)/distz
       IBP%intpointi(1) = IBP%xi+dirx
       IBP%intpointj(1) = IBP%yj
       IBP%intpointk(1) = IBP%zk+dirz
       IBP%intcoef(1) = abs(distx*t)/abs(xPr(IBP%xi+dirx)-xPr(IBP%xi))
       IBP%intpointi(2) = IBP%xi
       IBP%intpointj(2) = IBP%yj
       IBP%intpointk(2) = IBP%zk+dirz
       IBP%intcoef(2) = 1-IBP%intcoef(1)
       IBP%interp = 2
       IBP%dist = sqrt((distx*t)**2+(distz*t)**2)
      else
       t = (xPr(xi+dirx)-x)/distx
       IBP%intpointi(1) = IBP%xi+dirx
       IBP%intpointj(1) = IBP%yj
       IBP%intpointk(1) = IBP%zk+dirz
       IBP%intcoef(1) = abs(distz*t)/abs(zPr(IBP%zk+dirz)-zPr(IBP%zk))
       IBP%intpointi(2) = IBP%xi+dirx
       IBP%intpointj(2) = IBP%yj
       IBP%intpointk(2) = IBP%zk
       IBP%intcoef(2) = 1-IBP%intcoef(1)
       IBP%interp = 2
       IBP%dist = sqrt((distx*t)**2+(distz*t)**2)
      endif

     else                 !plane xy

      if (abs(distx/disty)>abs(xPr(xi+dirx)-x)/abs(yPr(yj+diry)-y)) then
       t = (yPr(yj+diry)-y)/disty
       IBP%intpointi(1) = IBP%xi+dirx
       IBP%intpointj(1) = IBP%yj+diry
       IBP%intpointk(1) = IBP%zk
       IBP%intcoef(1) = abs(distx*t)/abs(xPr(IBP%xi+dirx)-xPr(IBP%xi))
       IBP%intpointi(2) = IBP%xi
       IBP%intpointj(2) = IBP%yj+diry
       IBP%intpointk(2) = IBP%zk
       IBP%intcoef(2) = 1-IBP%intcoef(1)
       IBP%interp = 2
       IBP%dist = sqrt((distx*t)**2+(disty*t)**2)
      else
       t = (xPr(xi+dirx)-x)/distx
       IBP%intpointi(1) = IBP%xi+dirx
       IBP%intpointj(1) = IBP%yj+diry
       IBP%intpointk(1) = IBP%zk
       IBP%intcoef(1) = abs(disty*t)/abs(yPr(IBP%yj+diry)-yPr(IBP%yj))
       IBP%intpointi(2) = IBP%xi+dirx
       IBP%intpointj(2) = IBP%yj
       IBP%intpointk(2) = IBP%zk
       IBP%intcoef(2) = 1-IBP%intcoef(1)
       IBP%interp = 2
       IBP%dist = sqrt((distx*t)**2+(disty*t)**2)
      endif
     endif

    else                         !more than two free directions

     tx = (xPr(xi+dirx)-x)/distx   !a vector pointing from a free point to this forcing point.
     ty = (yPr(yj+diry)-y)/disty
     tz = (zPr(zk+dirz)-z)/distz

     IBP%interp = 4
     if (tx<=ty.and.tx<=tz) then
      !coordinates of interpolation point are therefore
      !xPr(xi+dirx) = x+tx*distx
      !y+tx*disty
      !z+tx*distz
      IBP%dist = sqrt((distx*tx)**2+(disty*tx)**2+(distz*tx)**2)
      IBP%intpointi(1) = IBP%xi+dirx
      IBP%intpointj(1) = IBP%yj+diry
      IBP%intpointk(1) = IBP%zk+dirz
      IBP%intcoef(1) = abs(disty*tx)/abs(yPr(yj+diry)-y)*&
                                abs(distz*tx)/abs(zPr(zk+dirz)-z)
      IBP%intpointi(2) = IBP%xi+dirx
      IBP%intpointj(2) = IBP%yj
      IBP%intpointk(2) = IBP%zk+dirz
      IBP%intcoef(2) = (1-abs(disty*tx)/abs(yPr(yj+diry)-y))*&
                                abs(distz*tx)/abs(zPr(zk+dirz)-z)
      IBP%intpointi(3) = IBP%xi+dirx
      IBP%intpointj(3) = IBP%yj+diry
      IBP%intpointk(3) = IBP%zk
      IBP%intcoef(3) = abs(disty*tx)/abs(yPr(yj+diry)-y)*&
                                (1-abs(distz*tx)/abs(zPr(zk+dirz)-z))
      IBP%intpointi(4) = IBP%xi+dirx
      IBP%intpointj(4) = IBP%yj
      IBP%intpointk(4) = IBP%zk
      IBP%intcoef(4) = (1-abs(disty*tx)/abs(yPr(yj+diry)-y))*&
                                (1-abs(distz*tx)/abs(zPr(zk+dirz)-z))
     elseif (ty<=tx.and.ty<=tz) then
      !the same with ty
      IBP%dist = sqrt((distx*ty)**2+(disty*ty)**2+(distz*ty)**2)
      IBP%intpointi(1) = IBP%xi+dirx
      IBP%intpointj(1) = IBP%yj+diry
      IBP%intpointk(1) = IBP%zk+dirz
      IBP%intcoef(1) = abs(distx*ty)/abs(xPr(xi+dirx)-x)*&
                                abs(distz*ty)/abs(zPr(zk+dirz)-z)
      IBP%intpointi(2) = IBP%xi
      IBP%intpointj(2) = IBP%yj+diry
      IBP%intpointk(2) = IBP%zk+dirz
      IBP%intcoef(2) = (1-abs(distx*ty)/abs(xPr(xi+dirx)-x))*&
                                abs(distz*ty)/abs(zPr(zk+dirz)-z)
      IBP%intpointi(3) = IBP%xi+dirx
      IBP%intpointj(3) = IBP%yj+diry
      IBP%intpointk(3) = IBP%zk
      IBP%intcoef(3) = abs(distx*ty)/abs(xPr(xi+dirx)-x)*&
                                (1-abs(distz*ty)/abs(zPr(zk+dirz)-z))
      IBP%intpointi(4) = IBP%xi
      IBP%intpointj(4) = IBP%yj+diry
      IBP%intpointk(4) = IBP%zk
      IBP%intcoef(4) = (1-abs(distx*ty)/abs(xPr(xi+dirx)-x))*&
                                (1-abs(distz*ty)/abs(zPr(zk+dirz)-z))
     else
      !the same with tz
      IBP%dist = sqrt((distx*tz)**2+(disty*tz)**2+(distz*tz)**2)
      IBP%intpointi(1) = IBP%xi+dirx
      IBP%intpointj(1) = IBP%yj+diry
      IBP%intpointk(1) = IBP%zk+dirz
      IBP%intcoef(1) = abs(distx*tz)/abs(xPr(xi+dirx)-x)*&
                                abs(disty*tz)/abs(yPr(yj+diry)-y)
      IBP%intpointi(2) = IBP%xi
      IBP%intpointj(2) = IBP%yj+diry
      IBP%intpointk(2) = IBP%zk+dirz
      IBP%intcoef(2) = (1-abs(distx*tz)/abs(xPr(xi+dirx)-x))*&
                                abs(disty*tz)/abs(yPr(yj+diry)-y)
      IBP%intpointi(3) = IBP%xi+dirx
      IBP%intpointj(3) = IBP%yj
      IBP%intpointk(3) = IBP%zk+dirz
      IBP%intcoef(3) = abs(distx*tz)/abs(xPr(xi+dirx)-x)*&
                                (1-abs(disty*tz)/abs(yPr(yj+diry)-y))
      IBP%intpointi(4) = IBP%xi
      IBP%intpointj(4) = IBP%yj
      IBP%intpointk(4) = IBP%zk+dirz
      IBP%intcoef(4) = (1-abs(distx*tz)/abs(xPr(xi+dirx)-x))*&
                                (1-abs(disty*tz)/abs(yPr(yj+diry)-y))
     endif

    endif
  end subroutine GeTScalIFlBPoint



  subroutine InitImBoundaries
    type(TIBPoint) IBP
    type(TScalFlIBPoint) SIBP
    integer i,j,k


    do k = 1,Unz
     do j = 1,Uny
      do i = 1,Unx
       if (Utype(i,j,k)>0) then
        if (Utype(i+1,j,k)==0.or.Utype(i-1,j,k)==0.or.Utype(i,j+1,k)==0&
          .or.Utype(i,j-1,k)==0.or.Utype(i,j,k+1)==0.or.Utype(i,j,k-1)==0)  then
            call  GetUIBPoint(IBP,i,j,k)
            call  AddIBPoint(IBP)
        endif
       endif
      enddo
     enddo
    enddo

    do k = 1,Vnz
     do j = 1,Vny
      do i = 1,Vnx
       if (Vtype(i,j,k)>0) then
        if (Vtype(i+1,j,k)==0.or.Vtype(i-1,j,k)==0.or.Vtype(i,j+1,k)==0&
          .or.Vtype(i,j-1,k)==0.or.Vtype(i,j,k+1)==0.or.Vtype(i,j,k-1)==0)  then
            call  GetVIBPoint(IBP,i,j,k)
            call  AddIBPoint(IBP)
        endif
       endif
      enddo
     enddo
    enddo

    do k = 1,Wnz
     do j = 1,Wny
      do i = 1,Wnx
       if (Wtype(i,j,k)>0) then
        if (Wtype(i+1,j,k)==0.or.Wtype(i-1,j,k)==0.or.Wtype(i,j+1,k)==0&
          .or.Wtype(i,j-1,k)==0.or.Wtype(i,j,k+1)==0.or.Wtype(i,j,k-1)==0)  then
            call  GetWIBPoint(IBP,i,j,k)
            call  AddIBPoint(IBP)
        endif
       endif
      enddo
     enddo
    enddo

    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       if (Prtype(i,j,k)>0) then
        if (Prtype(i+1,j,k)==0.or.Prtype(i-1,j,k)==0.or.Prtype(i,j+1,k)==0&
          .or.Prtype(i,j-1,k)==0.or.Prtype(i,j,k+1)==0.or.Prtype(i,j,k-1)==0)  then
            call  GeTScalIFlBPoint(SIBP,i,j,k)
            call  AddScalIBPoint(SIBP)
        endif
       endif
      enddo
     enddo
    enddo

    call MoveIBPointsToArray

    call DeallIBP(FirstIBPoint)
    call DeallScalFLIBP(FirstScalFlIBPoint)

  end subroutine InitImBoundaries


  subroutine MoveIBPointsToArray
    type(TIBPoint),pointer :: CurrentIBPoint
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


  pure function IBLinInterpolation(h0,h1,h2,vel1,vel2)
    real(KND) :: IBLinInterpolation
    real(KND),intent(in) :: h0,h1,h2,vel1,vel2

    if (h0<=h1-h0) then
     IBLinInterpolation = -(h0/(h1-h0))*vel1
  !   elseif  (h1-h0<h1/10._KND) then
  !    IBLinInt= 0
    else
     IBLinInterpolation = -((h2-2*h0)*vel1+(2*h0-h1)*vel2)/(h2-h1)
    endif
  endfunction IBLinInterpolation


  pure function IBBiLinInterpolation(x0,y0,x1,y1,velx,vely,velxy)
    real(KND) :: IBBiLinInterpolation
    real(KND),intent(in) :: x0,y0,x1,y1
    real(KND),intent(in) :: velx,vely,velxy

    real(KND) :: a,b

    a = (x1-x0)/x1
    b = (y1-y0)/y1

    IBBiLinInterpolation = -(a*(1-b)*vely+(1-a)*(1-b)*velxy+(1-a)*b*velx)/(a*b)

  endfunction IBBiLinInterpolation


  pure function IBTriLinInterpolation(x0,y0,z0,x1,y1,z1,velx,vely,velxy,velz,velxz,velyz,velxyz)
    real(KND) :: IBTriLinInterpolation
    real(KND),intent(in) :: x0,y0,z0,x1,y1,z1
    real(KND),intent(in) :: velx,vely,velxy,velz,velxz,velyz,velxyz
    real(KND) :: a,b,c

    a = x0/x1
    b = y0/y1
    c = z0/z1

    IBTriLinInterpolation = - (a*(1-b)*(1-c)*velx+&
                   (1-a)*b*(1-c)*vely+&
                   (1-a)*(1-b)*c*velz+&
                   a*b*(1-c)*velxy+&
                   a*(1-b)*c*velxz+&
                   (1-a)*b*c*velyz+&
                   a*b*c*velxyz)/((1-a)*(1-b)*(1-c))

  endfunction IBTriLinInterpolation



  recursive subroutine DeallIBP(IBP)
    type(TIBPoint),pointer :: IBP

    if (associated(IBP%next)) call DeallIBP(IBP%next)

    deallocate(IBP)

  end subroutine DeallIBP


  recursive subroutine DeallScalFlIBP(IBP)
    type(TScalFlIBPoint),pointer :: IBP

    if (associated(IBP%next)) call DeallScalFlIBP(IBP%next)

    deallocate(IBP)

  end subroutine DeallScalFlIBP



  recursive subroutine DeallSB(SB)
    type(TSolidbody),pointer :: SB

    if (associated(SB%next)) call DeallSB(SB%next)

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

  end subroutine DeallSB

endmodule GEOMETRIC

module GeometricShapes
  use iso_c_binding, only: c_ptr
  use Kinds
  use Parameters
  use CGAL_Polyhedra
  
  implicit none

  private

  public TGeometricShape, TLine, TPlane, TPolyhedron, TCGALPolyhedron, &
         TSphere, TEllipsoid, TCylJacket, TCylinder, TTerrainPoint, TTerrain

         
  type r3
    real(knd) :: x,y,z
  end type
  
  type TBbox
    real(knd) :: xmin,ymin,zmin,xmax,ymax,zmax
  end type

         
  type,abstract :: TGeometricShape
    private
    type(TBbox) :: bbox
  contains
    procedure,private :: in_bbox
    procedure :: Inside => TGeometricShape_Inside
    procedure(Inside_interface),private,deferred :: InsideEps
    procedure(Closest_interface), deferred :: Closest
    procedure(Closest_interface), deferred :: ClosestOut
  end type

  abstract interface
     pure logical function Inside_interface(self,x,y,z,eps)
      import
      class(TGeometricShape),intent(in) :: self
      real(knd),intent(in) :: x,y,z
      real(knd),intent(in) ::eps
    end function
    subroutine Closest_interface(self,xnear,ynear,znear,x,y,z)
      import
      class(TGeometricShape),intent(in) :: self
      real(knd),intent(out) :: xnear,ynear,znear
      real(knd),intent(in) :: x,y,z
    end subroutine
  end interface

  type,extends(TGeometricShape) :: TLine
    real(knd) xc,yc,zc
    real(knd) a,b,c
  contains
    procedure :: InsideEps => TLine_Inside
    procedure :: Closest => TLine_Closest
    procedure :: ClosestOut => TLine_Closest
  end type TLine


  type,extends(TGeometricShape) :: TPlane
    real(knd) a,b,c,d      !ax+by+cz+d/=0 for inner half-space
    logical gl             !T > in ineq. above F < in ineq. above
    logical :: rough = .false.!T rough surface, F flat surface
    real(knd) z0           !roughness parameter
  contains
    procedure :: InsideEps => TPlane_Inside
    procedure :: Closest => TPlane_Closest
    procedure :: ClosestOut => TPlane_Closest
  end type TPlane


  type,extends(TGeometricShape) :: TPolyhedron
    integer :: nplanes = 0
    type(TPlane),dimension(:),allocatable :: Planes !intersection of half-spaces
  contains
    procedure,private :: InsideEps => TPolyhedron_Inside
    procedure :: Closest => TPolyhedron_Closest
    procedure :: ClosestOut => TPolyhedron_ClosestOut
  end type TPolyhedron


  type,extends(TGeometricShape) :: TCGALPolyhedron
    type(c_ptr) :: cgalptr
    type(r3)    :: ref
  contains
    procedure,private :: InsideEps => TCGALPolyhedron_Inside
    procedure :: Closest => TCGALPolyhedron_Closest
    procedure :: ClosestOut => TCGALPolyhedron_ClosestOut
    procedure :: ReadOff => TCGALPolyhedron_ReadOff
    procedure :: InitBbox => TCGALPolyhedron_InitBbox
  end type TCGALPolyhedron


  type,extends(TGeometricShape) :: TSphere
    real(knd) xc,yc,zc,r
    logical :: rough = .false. !T rough surface, F flat surface
    real(knd) z0            !roughness parameter
  contains
    procedure,private :: InsideEps => TSphere_Inside
    procedure :: Closest => TSphere_Closest
    procedure :: ClosestOut => TSphere_ClosestOut
  end type TSphere


  type,extends(TGeometricShape) :: TEllipsoid
    real(knd) xc,yc,zc,a,b,c
    logical :: rough = .false. !T rough surface, F flat surface
    real(knd) z0            !roughness parameter
  contains
    procedure,private :: InsideEps => TEllipsoid_Inside
    procedure :: Closest => TEllipsoid_Closest
    procedure :: ClosestOut => TEllipsoid_ClosestOut
  end type TEllipsoid


  type,extends(TGeometricShape) :: TCylJacket
    real(knd) xc,yc,zc
    real(knd) a,b,c
    real(knd) r
    logical :: rough = .false. !T rough surface, F flat surface
    real(knd) z0            !roughness parameter
  contains
    procedure :: InsideEps => TCylJacket_Inside
    procedure :: Closest => TCylJacket_Closest
    procedure :: ClosestOut => TCylJacket_Closest
  end type TCylJacket


  type,extends(TGeometricShape) :: TCylinder
    type(TCylJacket) Jacket
    type(TPlane),allocatable :: Plane1 ,Plane2
  contains
    procedure,private :: InsideEps => TCylinder_Inside
    procedure :: Closest => TCylinder_Closest
    procedure :: ClosestOut => TCylinder_ClosestOut
  end type TCylinder


  type TTerrainPoint
    real(knd) :: elev = 0
    logical :: rough = .false.
    real(knd) z0
  end type TTerrainPoint


  type,extends(TGeometricShape) :: TTerrain
    type(TTerrainPoint),dimension(:,:),allocatable :: UPoints,VPoints,PrPoints !allocate with a buffer of width 1 (i.e. 0:Xnx)
  contains
    procedure,private :: InsideEps => TTerrain_Inside
    procedure :: Closest => TTerrain_Closest
    procedure :: ClosestOut => TTerrain_ClosestOut
  end type TTerrain
   
  interface Closest
    module procedure TLine_Closest
  end interface
 
  !initializers
  interface TLine
    module procedure TLine_Init
  end interface

  interface TEllipsoid
    module procedure TEllipsoid_Init
  end interface

contains

  pure logical function in_bbox(self,x,y,z,eps)
    class(TGeometricShape),intent(in) :: self
    real(knd),intent(in) :: x,y,z,eps
    
    associate(b=>self%bbox)
      in_bbox = (x+eps>=b%xmin).and. &
                (x-eps<=b%xmax).and. &
                (y+eps>=b%ymin).and. &
                (y-eps<=b%ymax).and. &
                (z+eps>=b%zmin).and. &
                (z-eps<=b%zmax)
    end associate
  end function
  
 pure logical function TGeometricShape_Inside(self,x,y,z,eps) result(ins)
    class(TGeometricShape),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),optional,intent(in) ::eps
    
    if (present(eps)) then
      ins = self%InsideEps(x,y,z,eps)
    else
      ins = self%InsideEps(x,y,z,epsilon(1._knd))
    end if
    
  end function
  

  real(knd) function PointDist(x1,y1,z1,x2,y2,z2)
    real(knd),intent(in) :: x1,y1,z1,x2,y2,z2

    PointDist = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
  end function PointDist



  subroutine TLine_Closest(self,xnear,ynear,znear,x,y,z)
    class(TLine),intent(in) :: self
    real(knd),intent(out)  :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) t

    if (self%a/=0 .or. self%b/=0 .or. self%c/=0) then
      t = ( self%a*(x-self%xc) + self%b*(y-self%yc) + self%c*(z-self%zc) ) / (self%a**2 + self%b**2 + self%c**2)
    else
      t = 0
    endif

    xnear = self%xc + self%a * t
    ynear = self%yc + self%b * t
    znear = self%zc + self%c * t
  end subroutine TLine_Closest


  
  pure logical function TLine_Inside(self,x,y,z,eps) result(ins)
    class(TLine),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps

    ins = .false.
  end function TLine_Inside
  
  pure logical function TPlane_Inside(self,x,y,z,eps) result(ins)
    class(TPlane),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps

    if (self%GL) then
      if (self%a*x+self%b*y+self%c*z+self%d>=-eps) then
       ins = .true.
      else
       ins = .false.
      endif
    else
      if (self%a*x+self%b*y+self%c*z+self%d<=eps) then
       ins = .true.
      else
       ins = .false.
      endif
    endif
  end function TPlane_Inside



  pure logical function TPolyhedron_Inside(self,x,y,z,eps) result(ins)
    class(TPolyhedron),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps
    integer i

    if (self%nplanes>0) then
      ins = .true.
      do i = 1,self%nplanes

       ins = self%Planes(i)%Inside(x,y,z,eps)

       if (.not. ins) exit
      enddo
    else
      ins = .false.
    endif
  end function TPolyhedron_Inside

  
  pure logical function TCGALPolyhedron_Inside(self,x,y,z,eps) result(ins)
    class(TCGALPolyhedron),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps
    
    !it is probably better to use one fixed reference point, 
    !because test segment aligned with probable walls could cause problems
    if (self%in_bbox(x,y,z,eps)) then
      ins = &
           cgal_polyhedron_inside(self%cgalptr, &
                                  x,y,z, &
                                  self%ref%x,self%ref%y,self%ref%z)
      if (eps>0) then                            
        if (.not.ins) ins = &
             cgal_polyhedron_inside(self%cgalptr, &
                                    x-eps,y,z, &
                                    self%ref%x,self%ref%y,self%ref%z)
        if (.not.ins) ins = &
             cgal_polyhedron_inside(self%cgalptr, &
                                    x+eps,y,z, &
                                    self%ref%x,self%ref%y,self%ref%z)
        if (.not.ins) ins = &
             cgal_polyhedron_inside(self%cgalptr, &
                                    x,y-eps,z, &
                                    self%ref%x,self%ref%y,self%ref%z)
        if (.not.ins) ins = &
             cgal_polyhedron_inside(self%cgalptr, &
                                    x,y-eps,z, &
                                    self%ref%x,self%ref%y,self%ref%z)
        if (.not.ins) ins = &
             cgal_polyhedron_inside(self%cgalptr, &
                                    x,y,z-eps, &
                                    self%ref%x,self%ref%y,self%ref%z)
        if (.not.ins) ins = &
             cgal_polyhedron_inside(self%cgalptr, &
                                    x,y,z+eps, &
                                    self%ref%x,self%ref%y,self%ref%z)
      end if
    else
      ins = .false.
    end if
  
  end function TCGALPolyhedron_Inside



  pure logical function TSphere_Inside(self,x,y,z,eps) result(ins)
    class(TSphere),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps

    if ((self%xc-x)**2+(self%yc-y)**2+(self%zc-z)**2<=(self%r+eps)**2) then
     ins = .true.
    else
     ins = .false.
    endif

  end function TSphere_Inside



  pure logical function TEllipsoid_Inside(self,x,y,z,eps) result(ins)
    class(TEllipsoid),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps
    
    if (self%in_bbox(x,y,z,eps)) then
      ins =  ((self%xc-x)**2/(self%a)**2 + &
            (self%yc-y)**2/(self%b)**2 + &
            (self%zc-z)**2/(self%c)**2 <= 1._knd+sqrt(eps))
    else
        ins = .false.
    end if
  end function TEllipsoid_Inside



  pure real(knd) function LineDist(x,y,z,xl,yl,zl,a,b,c)
    real(knd),intent(in) :: x,y,z,xl,yl,zl,a,b,c
    real(knd) t

    if (((a/=0).or.(b/=0)).or.(c/=0)) then
     t = (a*(x-xl)+b*(y-yl)+c*(z-zl))/(a**2+b**2+c**2)
    else
     t = 0
    endif

    LineDist = sqrt((xl+a*t-x)**2+(yl+b*t-y)**2+(zl+c*t-z)**2)
  end function LineDist



  pure logical function TCylJacket_Inside(self,x,y,z,eps) result(ins)
    class(TCylJacket),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps

    if (LineDist(x,y,z,self%xc,self%yc,self%zc,self%a,self%b,self%c)<=self%r+eps) then
     ins = .true.
    else
     ins = .false.
    endif
  end function TCylJacket_Inside



  pure logical function TCylinder_Inside(self,x,y,z,eps) result(ins)
    class(TCylinder),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps

    ins = .true.

    if (.not.self%Jacket%Inside(x,y,z,eps)) ins = .false.

    if (ins.and.allocated(self%Plane1)) then
           if (.not.self%Plane1%Inside(x,y,z,eps)) ins = .false.
    endif
    if (ins.and.allocated(self%Plane2)) then
          if (.not.self%Plane2%Inside(x,y,z,eps)) ins = .false.
    endif
  end function TCylinder_Inside



  pure subroutine TTerrain_GridCoords(x2,y2,xi,yj,comp)
    real(knd),intent(in) :: x2,y2
    integer,intent(out) :: xi,yj,comp
    real(knd) x,y,distPr,distU,distV
    integer xPri,yPrj,xUi,yVj,i

    x = x2
    y = y2

    if (gridtype==uniformgrid) then

        xPri = min(max(nint( (x - xU(0))/dxmin + 0.5_knd ),1),Prnx+1)

        yPrj = min(max(nint( (y - yV(0))/dymin + 0.5_knd ),1),Prny+1)

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



  pure logical function TTerrain_Inside(self,x,y,z,eps)
    class(TTerrain),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps
    logical ins
    integer xi,yj,comp

    ins = .false.
    call TTerrain_GridCoords(x,y,xi,yj,comp)

    if (comp==1) then
     if (z<=self%UPoints(xi,yj)%elev+eps) ins = .true.
    elseif (comp==2) then
     if (z<=self%VPoints(xi,yj)%elev+eps) ins = .true.
    elseif (comp==3) then
     if (z<=self%PrPoints(xi,yj)%elev+eps) ins = .true.
    endif

    TTerrain_Inside = ins
  end function TTerrain_Inside


  subroutine TPlane_Closest(self,xnear,ynear,znear,x,y,z)
    class(TPlane),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) ::x,y,z
    real(knd) t

    if (((self%a/=0).or.(self%b/=0)).or.(self%c/=0)) then
     t = -(self%a*x+self%b*y+self%c*z+self%d)/(self%a**2+self%b**2+self%c**2)
    else
     t = 0
    endif
    xnear = x+self%a*t
    ynear = y+self%b*t
    znear = z+self%c*t
  end subroutine TPlane_Closest



  subroutine TCylJacket_Closest(self,xnear,ynear,znear,x,y,z)
    class(TCylJacket),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) t,xl,yl,zl,a,b,c

    call Closest(TLine(self%xc,self%yc,self%zc,self%a,self%b,self%c),xl,yl,zl,x,y,z)

    a = x-xl
    b = y-yl
    c = z-zl
    t = self%r/sqrt(a**2+b**2+c**2)

    xnear = a*t+xl
    ynear = b*t+yl
    znear = c*t+zl
  end subroutine TCylJacket_Closest



  subroutine TCylinder_Closest(self,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
   class(TCylinder),intent(in) :: self
   real(knd),intent(out) :: xnear,ynear,znear
   real(knd),intent(in) :: x,y,z
   real(knd) xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2

   !!!Only for Planes perpendicular to jacket!!!!

   if (allocated(self%Plane1)) then
       if (allocated(self%Plane2)) then
          call self%Jacket%Closest(xJ,yJ,zJ,x,y,z)
          call self%Plane1%Closest(xP1,yP1,zP1,x,y,z)
          call self%Plane2%Closest(xP2,yP2,zP2,x,y,z)
          if (self%Jacket%Inside(x,y,z)) then
             if (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
                xnear = xP1
                ynear = yP1
                znear = zP1
             else
                xnear = xP2
                ynear = yP2
                znear = zP2
             endif
          elseif (self%Plane1%Inside(x,y,z).and.self%Plane2%Inside(x,y,z)) then
                xnear = xJ
                ynear = yJ
                znear = zJ
          elseif (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
           call self%Jacket%Closest(xnear,ynear,znear,xP1,yP1,zP1)
          else
           call self%Jacket%Closest(xnear,ynear,znear,xP2,yP2,zP2)
          endif
       else
          call self%Jacket%Closest(xJ,yJ,zJ,x,y,z)
          call self%Plane1%Closest(xP1,yP1,zP1,x,y,z)
          if (self%Jacket%Inside(x,y,z)) then
                xnear = xP1
                ynear = yP1
                znear = zP1
          elseif (self%Plane1%Inside(x,y,z)) then
                xnear = xJ
                ynear = yJ
                znear = zJ
          else
           call self%Jacket%Closest(xnear,ynear,znear,xP1,yP1,zP1)
         endif
       endif

    else

      call self%Jacket%Closest(xnear,ynear,znear,x,y,z)

    endif
  end subroutine TCylinder_Closest

  subroutine TSphere_Closest(self,xnear,ynear,znear,x,y,z)
    class(TSphere),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) t,a,b,c

    a = x-self%xc
    b = y-self%yc
    c = z-self%zc
    t = self%r/sqrt(a**2+b**2+c**2)
    xnear = a*t+self%xc
    ynear = b*t+self%yc
    znear = c*t+self%zc
  end subroutine TSphere_Closest


  subroutine TEllipsoid_Closest(self,xnear,ynear,znear,x,y,z)
    class(TEllipsoid),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) t,a,b,c

    stop "Error, closest point on an ellipsoid not implemented!"
    a = (x-self%xc)/self%a
    b = (y-self%yc)/self%b
    c = (z-self%zc)/self%c
    t = 1._knd/sqrt(a**2+b**2+c**2)
    xnear = self%a*a*t+self%xc
    ynear = self%b*b*t+self%yc
    znear = self%c*c*t+self%zc
  end subroutine TEllipsoid_Closest


  subroutine TPolyhedron_Closest(self,xnear,ynear,znear,x,y,z)
    use Lapack
    class(TPolyhedron),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    real(knd) :: dists(self%nplanes),xP(self%nplanes),yP(self%nplanes),zP(self%nplanes)
    real(knd) :: minv
    real(knd) :: ailine,biline,ciline
    real(knd) :: x0iline,y0iline,z0iline,xln,yln,zln,p
    integer   :: inearest,inearest2,inearest3
    integer   :: i

    integer,dimension(3)    :: ipivot,work2             !arguments of the LAPACK
    real(knd),dimension(3)  :: xg,bg,R,C,ferr,berr      !  routine GESVX
    real(knd),dimension(12) :: work
    real(knd),dimension(3,3) :: ag,af                    !
    integer   :: info                                   !
    real(knd) :: rcond                                  !
    character :: equed                                  !

   !Vzdalenosti od rovin, pokud je nejbl. bod roviny uvnitr jine, nebo puv. bod na vnitrni strane -> vzd. *-1
   !Pokud je nejbl. body +- eps. (norm. vektor!,dxmin/100) uvnitr a vne polyh -> hotovo
   !Jinak 2. nejbl. rovina v abs. hodnote -> prusecnice a nejbl bod na ni
   !Nejbl. bod na prusecnici. Pokud +-eps. uvnitr, (najit vekt. v rovine  kolme na prusecnic)-:hotovo
   !Jinak iterativne najit bod na prusecnici uvnitr

    do i = 1,self%nplanes
     call self%Planes(i)%Closest(xP(i),yP(i),zP(i),x,y,z)
     dists(i) = sqrt((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
     if (self%Planes(i)%Inside(x,y,z)) dists(i) = -ABS(dists(i))
    enddo
    !find nearest plane with
    inearest = 0
    minv = huge(minv)
    do i = 1,self%nplanes
     if (dists(i)>=0.and.dists(i)<minv) then
       inearest = i
       minv = dists(i)
     endif
    enddo

    if (inearest==0) then
      write(*,*) "no nearest"
      return
    endif

    if (self%Inside(xP(inearest),yP(inearest),zP(inearest), &
                    MIN(dxmin/10000._knd,dymin/10000._knd,dzmin/10000._knd) &
                   )) then
     xnear = xP(inearest)
     ynear = yP(inearest)
     znear = zP(inearest)
     return
    endif

    dists = abs(dists)

    inearest2 = 0
    minv = huge(minv)
    do i = 1,self%nplanes
     if (i/=inearest.and.dists(i)<minv) then
       inearest2 = i
       minv = dists(i)
     endif
    enddo

    inearest3 = 0
    minv = huge(minv)
    do i = 1,self%nplanes
     if (i/=inearest.and.i/=inearest2.and.dists(i)<minv) then
       inearest3 = i
       minv = dists(i)
     endif
    enddo

    ailine = self%Planes(inearest)%b*self%Planes(inearest2)%c-self%Planes(inearest)%c*self%Planes(inearest2)%b
    biline = self%Planes(inearest)%c*self%Planes(inearest2)%a-self%Planes(inearest)%a*self%Planes(inearest2)%c
    ciline = self%Planes(inearest)%a*self%Planes(inearest2)%b-self%Planes(inearest)%b*self%Planes(inearest2)%a

    if ( abs(ailine)<=epsilon(ailine).and. &
         abs(biline)<=epsilon(biline).and. &
         abs(ciline)<=epsilon(ciline) )   then
     write(*,*) "cross product 0"
     write(*,*) "numbers of planes",inearest,inearest2
     write(*,*) "x,y,z",x,y,z
     write(*,*) "pl1", self%Planes(inearest)%a, self%Planes(inearest)%b, &
                       self%Planes(inearest)%c, self%Planes(inearest)%d
     write(*,*) "pl2", self%Planes(inearest2)%a, self%Planes(inearest2)%b, &
                       self%Planes(inearest2)%c, self%Planes(inearest2)%d

     call self%Planes(inearest)%Closest(xln,yln,zln,x,y,z)
     call self%Planes(inearest2)%Closest(xnear,ynear,znear,x,y,z)

     if ((xln-x)**2+(yln-y)**2+(zln-z)**2<(xln-x)**2+(yln-y)**2+(zln-z)**2) then
       xnear = xln
       ynear = yln
       znear = zln
     end if

     return

    endif

    p = sqrt(ailine**2+biline**2+ciline**2)
    ailine = ailine/p
    biline = biline/p
    ciline = ciline/p

    if (abs(ciline)>=0.1_knd) then

       ag = reshape(source = (/ self%Planes(inearest)%a,self%Planes(inearest2)%a,0._knd, &
                           self%Planes(inearest)%b,self%Planes(inearest2)%b,0._knd, &
                           self%Planes(inearest)%c,self%Planes(inearest2)%c,1._knd /), &
                 shape = (/ 3,3 /))

       bg = (/ -self%Planes(inearest)%d,-self%Planes(inearest2)%d,(zW(Wnz+1)+zW(0))/2._knd /)

       call gesvx("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3, &
                  rcond,ferr,berr,work,work2,info)

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    elseif (abs(biline)>=0.1_knd) then

       ag = reshape(source = (/ self%Planes(inearest)%a,self%Planes(inearest2)%a,0._knd, &
                           self%Planes(inearest)%b,self%Planes(inearest2)%b,1._knd, &
                           self%Planes(inearest)%c,self%Planes(inearest2)%c,0._knd /), &
                 shape = (/ 3,3 /))

       bg = (/ -self%Planes(inearest)%d,-self%Planes(inearest2)%d,(yV(Vny+1)+yV(0))/2._knd /)

       call gesvx("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3, &
                  rcond,ferr,berr,work,work2,info)

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    else

       ag = reshape(source = (/ self%Planes(inearest)%a,self%Planes(inearest2)%a,1._knd, &
                           self%Planes(inearest)%b,self%Planes(inearest2)%b,0._knd, &
                           self%Planes(inearest)%c,self%Planes(inearest2)%c,0._knd /), &
                  shape = (/ 3,3 /))

       bg = (/ -self%Planes(inearest)%d,-self%Planes(inearest2)%d,(xU(Unx+1)+xU(0))/2._knd /)

       call gesvx("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3, &
                  rcond,ferr,berr,work,work2,info)

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    endif

    call Closest(TLine(x0iline,y0iline,z0iline,ailine,biline,ciline),xln,yln,zln,x,y,z)

    if (self%Inside(xln,yln,zln,min(dxmin/1000._knd,dymin/10000._knd,dzmin/10000._knd))) then
     xnear = xln
     ynear = yln
     znear = zln
     return
    endif



    ag = reshape(source = (/ self%Planes(inearest)%a,self%Planes(inearest2)%a,self%Planes(inearest3)%a, &
                        self%Planes(inearest)%b,self%Planes(inearest2)%b,self%Planes(inearest3)%b, &
                        self%Planes(inearest)%c,self%Planes(inearest2)%c,self%Planes(inearest3)%c /), &
               shape = (/ 3,3 /))
    bg = (/ -self%Planes(inearest)%d,-self%Planes(inearest2)%d,-self%Planes(inearest3)%d /)

    call gesvx("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3, &
               rcond,ferr,berr,work,work2,info)

    xnear = xg(1)
    ynear = xg(2)
    znear = xg(3)

  end subroutine TPolyhedron_Closest

  
  
  subroutine TCGALPolyhedron_Closest(self,xnear,ynear,znear,x,y,z)
    class(TCGALPolyhedron),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call cgal_polyhedron_closest(self%cgalptr, &
                                 x,y,z, &
                                 xnear, ynear, znear)
  
  end subroutine TCGALPolyhedron_Closest
  
  
  

  subroutine TTerrain_Closest(self,xnear,ynear,znear,x,y,z)
    class(TTerrain),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    type(TPlane) :: Pl
    integer     :: xi,yj,comp
    real(knd)   :: a,b,zloc
    xnear = x
    ynear = y
    znear = z

    call TTerrain_GridCoords(x,y,xi,yj,comp)

    if (comp==1) then  !Construct a tangent plane using first derivatives

     zloc = self%UPoints(xi,yj)%elev
     a = (self%PrPoints(xi+1,yj)%elev - self%PrPoints(xi,yj)%elev) / dxU(xi)
     b = (self%UPoints(xi,yj+1)%elev - self%UPoints(xi,yj-1)%elev) / (yPr(yj+1)-yPr(yj-1))

     Pl%a = a
     Pl%b = b
     Pl%c = -1._knd
     Pl%d = -a*x-b*y+zloc

     call Pl%Closest(xnear,ynear,znear,x,y,z)

    elseif (comp==2) then

     zloc = self%VPoints(xi,yj)%elev
     a = (self%VPoints(xi+1,yj)%elev - self%VPoints(xi-1,yj)%elev) / (xPr(xi+1)-xPr(xi-1))
     b = (self%PrPoints(xi,yj+1)%elev - self%PrPoints(xi,yj)%elev) / dyV(yj)

     Pl%a = a
     Pl%b = b
     Pl%c = -1._knd
     Pl%d = -a*x-b*y+zloc

     call Pl%Closest(xnear,ynear,znear,x,y,z)

    elseif (comp==3) then

     zloc = self%PrPoints(xi,yj)%elev
     a = (self%UPoints(xi,yj)%elev - self%UPoints(xi-1,yj)%elev) / dxPr(xi)
     b = (self%VPoints(xi,yj)%elev - self%VPoints(xi,yj-1)%elev) / dyPr(yj)

     Pl%a = a
     Pl%b = b
     Pl%c = -1._knd
     Pl%d = - a*x - b*y + zloc

     call Pl%Closest(xnear,ynear,znear,x,y,z)

    endif
  end subroutine TTerrain_Closest





  subroutine TPolyhedron_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(TPolyhedron),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) dists(self%nplanes),xP(self%nplanes),yP(self%nplanes),zP(self%nplanes),minv
    integer inearest,i

    dists = huge(minv)
    do i = 1,self%nplanes
     call self%Planes(i)%Closest(xP(i),yP(i),zP(i),x,y,z)
     dists(i) = sqrt((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
    enddo

    inearest = 0
    minv = huge(minv)
    do i = 1,self%nplanes
     if (dists(i)>=0.and.dists(i)<minv) then
       inearest = i
       minv = dists(i)
     endif
    enddo

    xnear = xP(inearest)
    ynear = yP(inearest)
    znear = zP(inearest)
  end subroutine TPolyhedron_ClosestOut


  subroutine TCGALPolyhedron_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(TCGALPolyhedron),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%Closest(xnear,ynear,znear,x,y,z)
  
  end subroutine TCGALPolyhedron_ClosestOut

  
  subroutine TCylinder_ClosestOut(self,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
    class(TCylinder),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2

    if (allocated(self%Plane1)) then
      call self%Plane1%Closest(xP1,yP1,zP1,x,y,z)
    else
     xP1 = sqrt(huge(xP1))/10
     yP1 = sqrt(huge(yP1))/10
     zP1 = sqrt(huge(zP1))/10
    endif
    if (allocated(self%Plane2)) then
      call self%Plane2%Closest(xP2,yP2,zP2,x,y,z)
    else
     xP2 = sqrt(huge(xP2))/10
     yP2 = sqrt(huge(yP2))/10
     zP2 = sqrt(huge(zP2))/10
    endif
    call self%Jacket%Closest(xJ,yJ,zJ,x,y,z)

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
  end subroutine TCylinder_ClosestOut



  subroutine TSphere_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(TSphere),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%Closest(xnear,ynear,znear,x,y,z)
  end subroutine TSphere_ClosestOut



  subroutine TEllipsoid_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(TEllipsoid),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%Closest(xnear,ynear,znear,x,y,z)
  end subroutine TEllipsoid_ClosestOut



  subroutine TTerrain_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(TTerrain),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%Closest(xnear,ynear,znear,x,y,z)
  end subroutine TTerrain_ClosestOut


  subroutine TCGALPolyhedron_ReadOff(self,filename)
    !reads geometry from an .off file
    use iso_c_binding, only: c_ptr,c_associated
    class(TCGALPolyhedron),intent(out) :: self
    character(*),intent(in) :: filename

    call cgal_polyhedron_read(self%cgalptr, filename)
    
    if (.not.c_associated(self%cgalptr)) then
      write(*,*) "Error reading polyhedron from ",filename
      stop
    end if
  end subroutine TCGALPolyhedron_ReadOff
  
  subroutine TCGALPolyhedron_InitBbox(self)
    class(TCGALPolyhedron),intent(inout) :: self

    associate(b=>self%bbox)
      call cgal_polyhedron_bbox(self%cgalptr, b%xmin, b%ymin, b%zmin, b%xmax, b%ymax, b%zmax)

      b%xmin=max(b%xmin,xU(-2))
      b%ymin=max(b%ymin,yV(-2))
      b%zmin=max(b%zmin,zW(-2))
      b%xmax=min(b%xmax,xU(Prnx+2))
      b%ymax=min(b%ymax,yV(Prny+2))
      b%zmax=min(b%zmax,zW(Prnz+2))
      
      self%ref = r3((b%xmin+b%xmax)/2, (b%ymin+b%ymax)/2, b%zmax + (b%zmax-b%zmin))
    end associate
  end subroutine TCGALPolyhedron_InitBbox
  
  
  
  
  function TLine_Init(xc,yc,zc,a,b,c) result(res)
    type(TLine) :: res
    real(knd),intent(in) :: xc, yc, zc, a, b, c
    
    res%xc = xc
    res%yc = yc
    res%zc = zc
    res%a  = a
    res%b  = b
    res%c  = c

  end function
  
  function TEllipsoid_Init(xc,yc,zc,a,b,c,rough,z0) result(res)
    type(TEllipsoid) :: res
    real(knd),intent(in) :: xc, yc, zc, a, b, c
    logical,intent(in),optional :: rough
    real(knd),intent(in),optional :: z0
    
    res%xc = xc
    res%yc = yc
    res%zc = zc
    res%a  = a
    res%b  = b
    res%c  = c
    if (present(rough)) res%rough = rough
    if (present(z0)) res%z0 = z0
    
    res%bbox%xmin = xc - a
    res%bbox%xmax = xc + a
    res%bbox%ymin = yc - b
    res%bbox%ymax = yc + b
    res%bbox%zmin = zc - c
    res%bbox%zmax = zc + c    

  end function
  
end module GeometricShapes







module TBody_class
  use Kinds, only: knd
  use Lists, only: TListable
  use GeometricShapes, only: TGeometricShape
  use Parameters
  
  implicit none

  private

  public TBody,Inside


  type, extends(TListable) :: TBody
     integer numofbody
     class(TGeometricShape),allocatable :: GeometricShape
  contains
     procedure :: Inside => CInside  !Hack around yet inidentified problem in GCC.
     procedure :: Closest
     procedure :: ClosestOut
     procedure :: ClosestOnLineOut
  end type
  interface Inside
    module procedure CInside
  end interface

contains


   logical function CInside(self,x,y,z,eps)
    class(TBody),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in),optional :: eps
    real(knd) x2,y2,z2
    real(knd) lx,ly,lz

    if (.not.allocated(self%GeometricShape)) then
      CInside = .false.
    else
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

      CInside = self%GeometricShape%Inside(x2,y2,z2,eps)

    end if

  end function CInside

  

  subroutine Closest(self,xnear,ynear,znear,x,y,z)
    class(TBody),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%GeometricShape%Closest(xnear,ynear,znear,x,y,z)

  end subroutine Closest

  

  subroutine ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(TBody),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%GeometricShape%ClosestOut(xnear,ynear,znear,x,y,z)

  end subroutine ClosestOut



  real(knd) function ClosestOnLineOut(self,x,y,z,x2,y2,z2) !Find t, such that x+(x2-x)*t lies on the boundary of the SB
    class(TBody),intent(in) :: self
    real(knd),intent(in) :: x,y,z,x2,y2,z2

    real(knd) t,t1,t2
    integer i

    t1 = 0
    t2 = 1
    if (self%Inside(x2,y2,z2)) then
     do
      t2 = t2*2._knd
      if (.not. self%Inside(x+(x2-x)*t2,y+(y2-y)*t2,z+(z2-z)*t2)) exit
     enddo
    endif
    t = (t1+t2)/2._knd

    do i = 1,20         !The bisection method with maximum 20 iterations (should be well enough)
     if (self%Inside(x+(x2-x)*t,y+(y2-y)*t,z+(z2-z)*t)) then
      t1 = t
     else
      t2 = t
     endif
     t = (t1+t2)/2._knd
     if (abs(t1-t2)<MIN(dxmin/1000._knd,dymin/1000._knd,dzmin/1000._knd))   exit
    enddo
    ClosestOnLineOut = t
  end function ClosestOnLineOut

  
end module TBody_class







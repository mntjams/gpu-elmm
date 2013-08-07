module r3_type
  use Kinds

  type r3
    real(knd) :: x,y,z
  end type
  
  interface assignment(=)
    module procedure r_to_r3
    module procedure v3_to_r3
    module procedure r3_to_v3
  end interface
  
contains
  
  elemental subroutine r_to_r3(lhs,rhs)
    type(r3),intent(out) :: lhs
    real(knd),intent(in) :: rhs
    lhs = r3(rhs,rhs,rhs)
  end subroutine
  
  pure subroutine v3_to_r3(lhs,rhs)
    type(r3),intent(out) :: lhs
    real(knd),intent(in) :: rhs(3)
    lhs = r3(rhs(1),rhs(2),rhs(3))
  end subroutine

  pure subroutine r3_to_v3(lhs,rhs)
    real(knd),intent(out) :: lhs(3)
    type(r3),intent(in)   :: rhs
    lhs = [rhs%x,rhs%y,rhs%z]
  end subroutine
  
  pure function v3(rhs) result(lhs)
    real(knd) :: lhs(3)
    type(r3),intent(in)   :: rhs
    lhs = [rhs%x,rhs%y,rhs%z]
  end function
  
end module r3_type


module GeometricShapes
  use iso_c_binding, only: c_ptr
  use Kinds
  use r3_type
  use Parameters
  use CGAL_Polyhedra
  
  implicit none

  private

  public GeometricShape, Line, Ray, Plane, ConvexPolyhedron, Polyhedron, &
         Sphere, Ellipsoid, CylJacket, Cylinder, TerrainPoint, Terrain, &
         Translation, Scaling, LinearTransform, Union, ArrayFromObst

  type Bbox
    real(knd) :: xmin = 0, ymin = 0, zmin = 0, xmax = 0, ymax = 0, zmax = 0
  end type

         
  type,abstract :: GeometricShape
    private
    type(Bbox) :: bbox = Bbox()
  contains
    procedure,private :: in_bbox
    procedure :: Inside => GeometricShape_Inside
    procedure(Inside_interface),private,deferred :: InsideEps
    procedure(Closest_interface), deferred :: Closest
    procedure :: ClosestOut => GeometricShape_ClosestOut
    procedure :: IntersectsRay => GeometricShape_IntersectsRay
  end type

  abstract interface
     logical function Inside_interface(self,x,y,z,eps)
      import
      class(GeometricShape),intent(in) :: self
      real(knd),intent(in) :: x,y,z
      real(knd),intent(in) ::eps
    end function
    subroutine Closest_interface(self,xnear,ynear,znear,x,y,z)
      import
      class(GeometricShape),intent(in) :: self
      real(knd),intent(out) :: xnear,ynear,znear
      real(knd),intent(in) :: x,y,z
    end subroutine
  end interface

  type,extends(GeometricShape) :: Line
    real(knd) xc,yc,zc
    real(knd) a,b,c
  contains
    procedure :: InsideEps => Line_Inside
    procedure :: Closest => Line_Closest
  end type

  type,extends(GeometricShape) :: Ray
    real(knd) xc,yc,zc
    real(knd) a,b,c
  contains
    procedure :: InsideEps => Ray_Inside
    procedure :: Closest => Ray_Closest
  end type
  
  type,extends(GeometricShape) :: Plane
    real(knd) a,b,c,d      !ax+by+cz+d/=0 for inner half-space
    logical gl             !T > in ineq. above F < in ineq. above
    logical :: rough = .false.!T rough surface, F flat surface
    real(knd) z0           !roughness parameter
  contains
    procedure :: InsideEps => Plane_Inside
    procedure :: Closest => Plane_Closest
  end type


  type,extends(GeometricShape) :: ConvexPolyhedron
    integer :: nplanes = 0
    type(Plane),dimension(:),allocatable :: Planes !intersection of half-spaces
  contains
    procedure,private :: InsideEps => ConvexPolyhedron_Inside
    procedure :: Closest => ConvexPolyhedron_Closest
    procedure :: ClosestOut => ConvexPolyhedron_ClosestOut
  end type


  type,extends(GeometricShape) :: Polyhedron
    type(c_ptr) :: cgalptr
    type(r3)    :: ref
  contains
    procedure,private :: ReadOff => Polyhedron_ReadOff
    procedure,private :: InitBbox => Polyhedron_InitBbox
    procedure,private :: InsideEps => Polyhedron_Inside
    procedure :: Closest => Polyhedron_Closest
    procedure :: IntersectsRay => Polyhedron_IntersectsRay
  end type


  type,extends(GeometricShape) :: Sphere
    real(knd) xc,yc,zc,r
    logical :: rough = .false. !T rough surface, F flat surface
    real(knd) :: z0 = 0           !roughness parameter
  contains
    procedure,private :: InsideEps => Sphere_Inside
    procedure :: Closest => Sphere_Closest
    procedure :: IntersectsRay => Sphere_IntersectsRay
  end type


  type,extends(GeometricShape) :: Ellipsoid
    real(knd) xc,yc,zc,a,b,c
    logical :: rough = .false. !T rough surface, F flat surface
    real(knd) z0            !roughness parameter
  contains
    procedure,private :: InsideEps => Ellipsoid_Inside
    procedure :: Closest => Ellipsoid_Closest
    procedure :: IntersectsRay => Ellipsoid_IntersectsRay
  end type


  type,extends(GeometricShape) :: CylJacket
    real(knd) xc,yc,zc
    real(knd) a,b,c
    real(knd) r
    logical :: rough = .false. !T rough surface, F flat surface
    real(knd) z0            !roughness parameter
  contains
    procedure :: InsideEps => CylJacket_Inside
    procedure :: Closest => CylJacket_Closest
  end type


  type,extends(GeometricShape) :: Cylinder
    type(CylJacket) Jacket
    type(Plane),allocatable :: Plane1 ,Plane2
  contains
    procedure,private :: InsideEps => Cylinder_Inside
    procedure :: Closest => Cylinder_Closest
    procedure :: ClosestOut => Cylinder_ClosestOut
  end type


  type TerrainPoint
    real(knd) :: elev = 0
    logical :: rough = .false.
    real(knd) z0
  end type


  type,extends(GeometricShape) :: Terrain
    type(TerrainPoint),dimension(:,:),allocatable :: UPoints,VPoints,PrPoints !allocate with a buffer of width 1 (i.e. 0:Xnx)
  contains
    procedure,private :: InsideEps => Terrain_Inside
    procedure,private,nopass,non_overridable :: GridCoords => Terrain_GridCoords
    procedure :: Closest => Terrain_Closest
  end type
  
  type,extends(GeometricShape) :: Translation
    class(GeometricShape),allocatable :: original
    type(r3) :: shift
  contains
    procedure,private :: Translation_in_bbox 
    procedure,private :: InsideEps => Translation_Inside
    procedure :: Closest => Translation_Closest
    procedure :: ClosestOut => Translation_ClosestOut
    procedure :: IntersectsRay => Translation_IntersectsRay
  end type
   
  type,extends(GeometricShape) :: Scaling
    class(GeometricShape),allocatable :: original
    type(r3) :: factor
  contains
    procedure,private :: in_bbox => Scaling_in_bbox 
    procedure,private :: InsideEps => Scaling_Inside
    procedure :: Closest => Scaling_Closest
    procedure :: ClosestOut => Scaling_ClosestOut
    procedure :: IntersectsRay => Scaling_IntersectsRay
  end type
   
  type,extends(GeometricShape) :: LinearTransform
    class(GeometricShape),allocatable,private :: original
    real(knd),private :: matrix(3,3), inv_matrix(3,3)
  contains
    procedure,private :: in_bbox => LinearTransform_in_bbox 
    procedure,private :: InsideEps => LinearTransform_Inside
    procedure :: Closest => LinearTransform_Closest
    procedure :: ClosestOut => LinearTransform_ClosestOut
    procedure :: IntersectsRay => LinearTransform_IntersectsRay
  end type
   
  type,extends(GeometricShape) :: Union
    class(GeometricShape),allocatable,private :: items(:)
    integer,private :: size
  contains
    procedure,private :: in_bbox => Union_in_bbox 
    procedure,private :: InsideEps => Union_Inside
    procedure :: Closest => Union_Closest
    procedure :: ClosestOut => Union_ClosestOut
    procedure :: IntersectsRay => Union_IntersectsRay
  end type
   
  interface Closest
    module procedure Line_Closest
  end interface
 
  !initializers
  interface Line
    module procedure Line_Init
  end interface

  interface Ray
    module procedure Ray_Init
    module procedure Ray_Init_v3
    module procedure Ray_Init_r3
  end interface

  interface Plane
    module procedure Plane_Init_3r_32
    module procedure Plane_Init_3r_64
  end interface

  interface ConvexPolyhedron
    module procedure ConvexPolyhedron_Init
  end interface

  interface Polyhedron
    module procedure Polyhedron_Init
  end interface
  
  interface Ellipsoid
    module procedure Ellipsoid_Init
  end interface
  
  interface Translation
    module procedure Translation_Init_3r
    module procedure Translation_Init_3r3
    module procedure Translation_Init_v3
  end interface

  interface Scaling
    module procedure Scaling_Init_r
    module procedure Scaling_Init_r3
    module procedure Scaling_Init_v3
  end interface
  
  interface LinearTransform
    module procedure LinearTransform_Init_scale_r
    module procedure LinearTransform_Init_scale_r3
    module procedure LinearTransform_Init_scale_v3
    module procedure LinearTransform_Init_rot
  end interface
  
  interface Union
    !TODO change to a polymorphic function that changes result according to data read
    module procedure Union_Init_Obst
  end interface
  
  real(knd),parameter :: unit_matrix_3(3,3) = reshape(source=[0,0,1,0,1,0,0,0,1], &
                                                      shape=[3,3])

contains

  !defaults
  
  logical function GeometricShape_Inside(self,x,y,z,eps) result(ins)
    class(GeometricShape),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),optional,intent(in) ::eps
    
    if (present(eps)) then
      ins = self%InsideEps(x,y,z,eps)
    else
      ins = self%InsideEps(x,y,z,epsilon(1._knd))
    end if
    
  end function
  
  subroutine GeometricShape_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(GeometricShape),intent(in) :: self
    real(knd),intent(out)  :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    
    call self%Closest(xnear,ynear,znear,x,y,z)
  end subroutine  
    
  logical function GeometricShape_IntersectsRay(self,r) result(intersects)
    class(GeometricShape),intent(in) :: self
    class(Ray),intent(in) :: r
    
    !shapes will be transparent for (solar) rays if they do not override this method.
    intersects = .false.
 
  end function
  

  
  
  !helpers
  
  real(knd) function PointDist(x1,y1,z1,x2,y2,z2)
    real(knd),intent(in) :: x1,y1,z1,x2,y2,z2

    PointDist = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
  end function
  
   real(knd) function LineDist(x,y,z,xl,yl,zl,a,b,c)
    real(knd),intent(in) :: x,y,z,xl,yl,zl,a,b,c
    real(knd) t

    if (((a/=0).or.(b/=0)).or.(c/=0)) then
     t = (a*(x-xl)+b*(y-yl)+c*(z-zl))/(a**2+b**2+c**2)
    else
     t = 0
    endif

    LineDist = sqrt((xl+a*t-x)**2+(yl+b*t-y)**2+(zl+c*t-z)**2)
  end function LineDist
 
  
  

  
  logical function in_bbox(self,x,y,z,eps)
    class(GeometricShape),intent(in) :: self
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
  
  
  

  
  



  function Line_Init(xc,yc,zc,a,b,c) result(res)
    type(Line) :: res
    real(knd),intent(in) :: xc, yc, zc, a, b, c
    
    res%xc = xc
    res%yc = yc
    res%zc = zc
    res%a  = a
    res%b  = b
    res%c  = c

  end function
  
  logical function Line_Inside(self,x,y,z,eps) result(ins)
    class(Line),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps

    ins = .false.
  end function
  

  subroutine Line_Closest(self,xnear,ynear,znear,x,y,z)
    class(Line),intent(in) :: self
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
  end subroutine


  
  
  
  
  
  
  
  function Ray_Init(xc,yc,zc,a,b,c) result(res)
    type(Ray) :: res
    real(knd),intent(in) :: xc, yc, zc, a, b, c
    
    res%xc = xc
    res%yc = yc
    res%zc = zc
    res%a  = a
    res%b  = b
    res%c  = c

  end function
  
  function Ray_Init_v3(c,v) result(res)
    type(Ray) :: res
    real(knd),intent(in) :: c(3),v(3)
    
    res%xc = c(1)
    res%yc = c(2)
    res%zc = c(3)
    res%a  = v(1)
    res%b  = v(2)
    res%c  = v(3)

  end function
  
  function Ray_Init_r3(c,v) result(res)
    type(Ray) :: res
    type(r3),intent(in) :: c,v
    
    res%xc = c%x
    res%yc = c%y
    res%zc = c%z
    res%a  = v%x
    res%b  = v%y
    res%c  = v%z

  end function
  
  
  logical function Ray_Inside(self,x,y,z,eps) result(ins)
    class(Ray),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps

    ins = .false.
  end function
  
  
  subroutine Ray_Closest(self,xnear,ynear,znear,x,y,z)
    class(Ray),intent(in) :: self
    real(knd),intent(out)  :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) t

    if (self%a/=0 .or. self%b/=0 .or. self%c/=0) then
      t = ( self%a*(x-self%xc) + self%b*(y-self%yc) + self%c*(z-self%zc) ) / (self%a**2 + self%b**2 + self%c**2)
    else
      t = 0
    endif

    xnear = self%xc + self%a * max(t,0._knd)
    ynear = self%yc + self%b * max(t,0._knd)
    znear = self%zc + self%c * max(t,0._knd)
  end subroutine


  
  
  
  
  
  
  
  
  
  
  function Plane_Init_3r(xc,yc,zc,a,b,c) result(res)
    !Point and an outward normal
    type(Plane) :: res
    real(knd),intent(in) :: xc, yc, zc, a, b, c
    
    res%a  = a
    res%b  = b
    res%c  = c
    res%gl = .false.
    res%d = - (a*xc + b*yc + c*zc)

  end function
  
  function Plane_Init_3r_32(xc,yc,zc,a,b,c) result(res)
    !Point and an outward normal
    type(Plane) :: res
    real(real32),intent(in) :: xc, yc, zc, a, b, c
    
    res = Plane_Init_3r(real(xc,knd), real(yc,knd), real(zc,knd), &
                        real(a,knd),  real(b,knd),  real(c,knd))
  end function
  
  function Plane_Init_3r_64(xc,yc,zc,a,b,c) result(res)
    !Point and an outward normal
    type(Plane) :: res
    real(real64),intent(in) :: xc, yc, zc, a, b, c
    
    res = Plane_Init_3r(real(xc,knd), real(yc,knd), real(zc,knd), &
                        real(a,knd),  real(b,knd),  real(c,knd))
  end function
  
  logical function Plane_Inside(self,x,y,z,eps) result(ins)
    class(Plane),intent(in) :: self
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
  end function

  subroutine Plane_Closest(self,xnear,ynear,znear,x,y,z)
    class(Plane),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) ::x,y,z
    real(knd) t

    if (abs(self%a)>tiny(1._knd).and. &
        abs(self%b)>tiny(1._knd).and. &
        abs(self%c)>tiny(1._knd)) then
      t = -(self%a*x+self%b*y+self%c*z+self%d)/(self%a**2+self%b**2+self%c**2)
    else
      t = 0
    endif
    xnear = x+self%a*t
    ynear = y+self%b*t
    znear = z+self%c*t
  end subroutine Plane_Closest



  
  
  
  
  
  
  
  function ConvexPolyhedron_Init(Planes) result(res)
    !Point and an outward normal
    type(ConvexPolyhedron) :: res
    type(Plane),intent(in) :: Planes(:)
    !limitation of gfortran 4.8 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=44672
    allocate(res%Planes(size(Planes)), source=Planes)
    res%nplanes = size(Planes)

  end function
  
  logical function ConvexPolyhedron_Inside(self,x,y,z,eps) result(ins)
    class(ConvexPolyhedron),intent(in) :: self
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
  end function
  
  subroutine ConvexPolyhedron_Closest(self,xnear,ynear,znear,x,y,z)
    use Lapack
    class(ConvexPolyhedron),intent(in) :: self
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

    call Closest(Line(x0iline,y0iline,z0iline,ailine,biline,ciline),xln,yln,zln,x,y,z)

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

  end subroutine ConvexPolyhedron_Closest

  subroutine ConvexPolyhedron_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(ConvexPolyhedron),intent(in) :: self
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
  end subroutine


  
  
  
  
  
  
  
  function Polyhedron_Init(filename) result (res)
    type(Polyhedron) :: res
    character(*) :: filename
    
    call res%ReadOff(filename)
    call res%InitBbox
  end function
  
  
  subroutine Polyhedron_ReadOff(self,filename)
    !reads geometry from an .off file
    use iso_c_binding, only: c_ptr,c_associated
    class(Polyhedron),intent(out) :: self
    character(*),intent(in) :: filename

    call cgal_polyhedron_read(self%cgalptr, filename)
    
    if (.not.c_associated(self%cgalptr)) then
      write(*,*) "Error reading polyhedron from ",filename
      stop
    end if
  end subroutine
  
  subroutine Polyhedron_InitBbox(self)
    class(Polyhedron),intent(inout) :: self

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
  end subroutine
  
  logical function Polyhedron_Inside(self,x,y,z,eps) result(ins)
    class(Polyhedron),intent(in) :: self
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
  
  end function Polyhedron_Inside

  subroutine Polyhedron_Closest(self,xnear,ynear,znear,x,y,z)
    class(Polyhedron),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call cgal_polyhedron_closest(self%cgalptr, &
                                 x,y,z, &
                                 xnear, ynear, znear)
  
  end subroutine
  
  logical function Polyhedron_IntersectsRay(self,r) result(intersects)
    class(Polyhedron),intent(in) :: self
    class(Ray),intent(in) :: r

      intersects = &
           cgal_polyhedron_intersects_ray(self%cgalptr, &
                                          r%xc, r%yc, r%zc, &
                                          r%a,  r%b,  r%c )
 
  end function



  
  
  
  
  
  
  
  
  
  

  logical function Sphere_Inside(self,x,y,z,eps) result(ins)
    class(Sphere),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps

    if ((self%xc-x)**2+(self%yc-y)**2+(self%zc-z)**2<=(self%r+eps)**2) then
     ins = .true.
    else
     ins = .false.
    endif

  end function

  subroutine Sphere_Closest(self,xnear,ynear,znear,x,y,z)
    class(Sphere),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) t,a,b,c

    a = x - self%xc
    b = y - self%yc
    c = z - self%zc
    t = self%r / sqrt(a**2+b**2+c**2)
    xnear = a*t + self%xc
    ynear = b*t + self%yc
    znear = c*t + self%zc
  end subroutine
  
  logical function Sphere_IntersectsRay(self,r) result(res)
     class(Sphere),intent(in) :: self
     class(Ray),intent(in) :: r
     real(knd) rc(3),rv(3) !transformed ray center and vector
     real(knd) a,b,c,D,t1,t2
     
     rc = [r%xc - self%xc, r%yc - self%yc, r%zc - self%zc]
     rv = [r%a, r%b, r%c]
     
     a = dot_product(rv,rv)
     b = 2 * dot_product(rc,rv)
     c = dot_product(rc,rc) - self%r**2
     
     D = b**2 - 4*a*c
     
     if (D<0) then
       res = .false.
     else
       t1 = (-b - sqrt(D)) / (2*a)
       t2 = (-b + sqrt(D)) / (2*a)
       res = (t1>=0 .or. t2>=0)
     end if
  end function
  
  
  
  
  
  
  
  
  

  function Ellipsoid_Init(xc,yc,zc,a,b,c,rough,z0) result(res)
    type(Ellipsoid) :: res
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
  
  logical function Ellipsoid_Inside(self,x,y,z,eps) result(ins)
    class(Ellipsoid),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps
    
    if (self%in_bbox(x,y,z,eps)) then
      ins =  ((self%xc-x)**2/(self%a)**2 + &
            (self%yc-y)**2/(self%b)**2 + &
            (self%zc-z)**2/(self%c)**2 <= 1._knd+sqrt(eps))
    else
        ins = .false.
    end if
  end function

  subroutine Ellipsoid_Closest(self,xnear,ynear,znear,x,y,z)
    class(Ellipsoid),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) t,a,b,c !auxiliary ray parameters

    ! NOT EXACT!
    a = (x - self%xc)/self%a
    b = (y - self%yc)/self%b
    c = (z - self%zc)/self%c
    t = 1._knd / sqrt(a**2+b**2+c**2)
    xnear = self%a * a*t + self%xc
    ynear = self%b * b*t + self%yc
    znear = self%c * c*t + self%zc
  end subroutine
  
  logical function Ellipsoid_IntersectsRay(self,r) result(res)
     class(Ellipsoid),intent(in) :: self
     class(Ray),intent(in) :: r
     type(Ray) :: r2
     !FIXME constructors contains irelevant properties and no keywords due to problem in ifort 13.1
     type(Sphere),parameter :: unit_sphere = Sphere( Bbox(0,0,0,0,0,0), 0._knd,  0._knd,&
                                                      0._knd,  1._knd,  .false.,  0)
     r2 = Ray((r%xc - self%xc)/self%a, (r%yc - self%yc)/self%b, (r%zc - self%zc)/self%c, &
                  r%a/self%a, r%b/self%b, r%c/self%c )
     
     res = unit_sphere%IntersectsRay(r2)

  end function
  






  logical function CylJacket_Inside(self,x,y,z,eps) result(ins)
    class(CylJacket),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps

    if (LineDist(x,y,z,self%xc,self%yc,self%zc,self%a,self%b,self%c)<=self%r+eps) then
     ins = .true.
    else
     ins = .false.
    endif
  end function

  subroutine CylJacket_Closest(self,xnear,ynear,znear,x,y,z)
    class(CylJacket),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) t,xl,yl,zl,a,b,c

    call Closest(Line(self%xc,self%yc,self%zc,self%a,self%b,self%c),xl,yl,zl,x,y,z)

    a = x-xl
    b = y-yl
    c = z-zl
    t = self%r/sqrt(a**2+b**2+c**2)

    xnear = a*t+xl
    ynear = b*t+yl
    znear = c*t+zl
  end subroutine

  
  
  
  
  logical function Cylinder_Inside(self,x,y,z,eps) result(ins)
    class(Cylinder),intent(in) :: self
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
  end function

  subroutine Cylinder_Closest(self,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
   class(Cylinder),intent(in) :: self
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
  end subroutine Cylinder_Closest

  subroutine Cylinder_ClosestOut(self,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
    class(Cylinder),intent(in) :: self
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
  end subroutine Cylinder_ClosestOut


  
  
  
  
  
  
  
  subroutine Terrain_GridCoords(x2,y2,xi,yj,comp)
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
  end subroutine Terrain_GridCoords

  logical function Terrain_Inside(self,x,y,z,eps)
    class(Terrain),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) :: eps
    logical ins
    integer xi,yj,comp

    ins = .false.
    call self%GridCoords(x,y,xi,yj,comp)

    if (comp==1) then
     if (z<=self%UPoints(xi,yj)%elev+eps) ins = .true.
    elseif (comp==2) then
     if (z<=self%VPoints(xi,yj)%elev+eps) ins = .true.
    elseif (comp==3) then
     if (z<=self%PrPoints(xi,yj)%elev+eps) ins = .true.
    endif

    Terrain_Inside = ins
  end function

  subroutine Terrain_Closest(self,xnear,ynear,znear,x,y,z)
    class(Terrain),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    type(Plane) :: Pl
    integer     :: xi,yj,comp
    real(knd)   :: a,b,zloc
    xnear = x
    ynear = y
    znear = z

    call Terrain_GridCoords(x,y,xi,yj,comp)

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
  end subroutine Terrain_Closest



  
  logical function Translation_in_bbox(self,x,y,z,eps) result(in)
    class(Translation),intent(in) :: self
    real(knd),intent(in) :: x,y,z,eps
    
    in = self%original%in_bbox(x - self%shift%x, & 
                               y - self%shift%y, &
                               z - self%shift%z, &
                               eps)
  end function


  logical function Translation_Inside(self,x,y,z,eps) result(ins)
    class(Translation),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps
    
    ins = self%original%InsideEps(x - self%shift%x, & 
                                  y - self%shift%y, &
                                  z - self%shift%z, &
                                  eps)
    
  end function
  
  
  subroutine Translation_Closest(self,xnear,ynear,znear,x,y,z)
    class(Translation),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%original%Closest(xnear, &
                               ynear, &
                               znear, &
                               x - self%shift%x, & 
                               y - self%shift%y, &
                               z - self%shift%z)
    xnear = xnear + self%shift%x
    ynear = ynear + self%shift%y
    znear = znear + self%shift%z
  end subroutine
  
  subroutine Translation_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(Translation),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%original%ClosestOut(xnear, &
                                  ynear, &
                                  znear, &
                                  x - self%shift%x, & 
                                  y - self%shift%y, &
                                  z - self%shift%z)
    xnear = xnear + self%shift%x
    ynear = ynear + self%shift%y
    znear = znear + self%shift%z
  end subroutine
  
  
  logical function Translation_IntersectsRay(self,r) result(res)
    class(Translation),intent(in) :: self
    class(Ray),intent(in) :: r
    
    res = self%original%IntersectsRay(Ray(r%xc - self%shift%x, &
                               r%yc - self%shift%y, &
                               r%zc - self%shift%z, &
                               r%a, &
                               r%b, &
                               r%c))
  end function
  
  elemental function Translation_Init_3r(original,sx,sy,sz) result(res)
    type(Translation) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: sx,sy,sz
   
    allocate(res%original, source=original)
    res%shift = [sx,sy,sz]
  end function
  
  elemental function Translation_Init_3r3(original,shift) result(res)
    type(Translation) :: res
    class(GeometricShape),intent(in) :: original
    type(r3),intent(in) :: shift
   
    allocate(res%original, source=original)
    res%shift = shift
  end function
  
  function Translation_Init_v3(original,shift) result(res)
    type(Translation) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: shift(3)
   
    allocate(res%original, source=original)
    res%shift = shift
  end function
  
  
  
  logical function Scaling_in_bbox(self,x,y,z,eps) result(in)
    class(Scaling),intent(in) :: self
    real(knd),intent(in) :: x,y,z,eps
    
    in = self%original%in_bbox(x / self%factor%x, & 
                               y / self%factor%y, &
                               z / self%factor%z, &
                               eps / (self%factor%x*self%factor%y*self%factor%z)**(1._knd/3))
  end function


  logical function Scaling_Inside(self,x,y,z,eps) result(ins)
    class(Scaling),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps

    ins = self%original%InsideEps(x / self%factor%x, & 
                                  y / self%factor%y, &
                                  z / self%factor%z, &
                                  eps / (self%factor%x*self%factor%y*self%factor%z)**(1._knd/3))
    
  end function
  
  
  subroutine Scaling_Closest(self,xnear,ynear,znear,x,y,z)
    class(Scaling),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%original%Closest(xnear, &
                               ynear, &
                               znear, &
                               x / self%factor%x, & 
                               y / self%factor%y, &
                               z / self%factor%z)
    xnear = xnear * self%factor%x
    ynear = ynear * self%factor%y
    znear = znear * self%factor%z
  end subroutine
  
  subroutine Scaling_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(Scaling),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%original%ClosestOut(xnear, &
                                  ynear, &
                                  znear, &
                                  x / self%factor%x, & 
                                  y / self%factor%y, &
                                  z / self%factor%z)
    xnear = xnear * self%factor%x
    ynear = ynear * self%factor%y
    znear = znear * self%factor%z
  end subroutine
  
  
  logical function Scaling_IntersectsRay(self,r) result(res)
    class(Scaling),intent(in) :: self
    class(Ray),intent(in) :: r
    
    res = self%original%IntersectsRay(Ray(r%xc / self%factor%x, &
                               r%yc / self%factor%y, &
                               r%zc / self%factor%z, &
                               r%a, &
                               r%b, &
                               r%c))
  end function
  
  elemental function Scaling_Init_r(original,scalar_factor) result(res)
    type(Scaling) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: scalar_factor
   
    allocate(res%original, source=original)
    res%factor = scalar_factor
  end function
  
  elemental function Scaling_Init_r3(original,factor) result(res)
    type(Scaling) :: res
    class(GeometricShape),intent(in) :: original
    type(r3),intent(in) :: factor
   
    allocate(res%original, source=original)
    res%factor = factor
  end function
  
  function Scaling_Init_v3(original,factor) result(res)
    type(Scaling) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: factor(3)
   
    allocate(res%original, source=original)
    res%factor = factor
  end function
  

  
  
  logical function LinearTransform_in_bbox(self,x,y,z,eps) result(in)
    class(LinearTransform),intent(in) :: self
    real(knd),intent(in) :: x,y,z,eps
    real(knd) :: xyz(3)
    
    xyz = matmul(self%inv_matrix, [x,y,z])
    
    in = self%original%in_bbox(xyz(1), & 
                               xyz(2), &
                               xyz(3), &
                               eps)
  end function


  logical function LinearTransform_Inside(self,x,y,z,eps) result(ins)
    class(LinearTransform),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps
    real(knd) :: xyz(3)
    
    xyz = matmul(self%inv_matrix, [x,y,z])
    
    ins = self%original%InsideEps(xyz(1), & 
                                  xyz(2), &
                                  xyz(3), &
                                  eps)
    
  end function
  
  
  subroutine LinearTransform_Closest(self,xnear,ynear,znear,x,y,z)
    class(LinearTransform),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) xyz(3),xyznear(3)

    xyz = matmul(self%inv_matrix, [x,y,z])
    
    call self%original%Closest(xyznear(1), &
                               xyznear(2), &
                               xyznear(3), &
                               xyz(1), & 
                               xyz(2), &
                               xyz(3))
                               
    xyznear = matmul(self%matrix, xyznear)
                        
    xnear = xyznear(1)
    ynear = xyznear(2)
    znear = xyznear(3)
  end subroutine
  
  subroutine LinearTransform_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(LinearTransform),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) xyz(3),xyznear(3)

    xyz = matmul(self%inv_matrix, [x,y,z])
    
    call self%original%Closest(xyznear(1), &
                               xyznear(2), &
                               xyznear(3), &
                               xyz(1), & 
                               xyz(2), &
                               xyz(3))
                               
    xyznear = matmul(self%matrix, xyznear)
                        
    xnear = xyznear(1)
    ynear = xyznear(2)
    znear = xyznear(3)
  end subroutine
  
  
  logical function LinearTransform_IntersectsRay(self,r) result(res)
    class(LinearTransform),intent(in) :: self
    class(Ray),intent(in) :: r
    real(knd) xyz(3),abc(3)
    
    xyz = matmul(self%inv_matrix, [r%xc, r%yc, r%zc])
    abc = matmul(self%inv_matrix, [r%a, r%b, r%c])
    
    res = self%original%IntersectsRay(Ray(xyz, abc))
  end function
  
  elemental function LinearTransform_Init_scale_r(original,scalar_factor) result(res)
    type(LinearTransform) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: scalar_factor
   
    allocate(res%original, source=original)
    
    res%matrix = unit_matrix_3 * scalar_factor
    res%inv_matrix = unit_matrix_3 / scalar_factor
  end function
  
  function LinearTransform_Init_scale_v3(original,factor) result(res)
    type(LinearTransform) :: res
    class(GeometricShape),intent(in) :: original
    real(knd),intent(in) :: factor(3)
    integer :: i
   
    allocate(res%original, source=original)
    
    res%matrix = 0
    res%inv_matrix = 0

    forall(i=1:3)  res%matrix(i,i) = factor(i)
    forall(i=1:3)  res%matrix(i,i) = 1._knd/factor(i)
  end function
  
  elemental function LinearTransform_Init_scale_r3(original,factor) result(res)
    type(LinearTransform) :: res
    class(GeometricShape),intent(in) :: original
    type(r3),intent(in) :: factor
    real(knd) :: v(3)
    integer :: i
   
    allocate(res%original, source=original)
    
    v = factor
    
    res%matrix = 0
    res%inv_matrix = 0

    forall(i=1:3)  res%matrix(i,i) = v(i)
    forall(i=1:3)  res%matrix(i,i) = 1._knd/v(i)
  end function
  
  pure function rotation_matrix_x(phi) result(res)
    real(knd) :: res(3,3)
    real(knd),intent(in) :: phi
    real(knd) :: c,s
    
    c = cos(phi)
    s = sin(phi)
    res(1,1) = 1
    res(2,1) = 0
    res(3,1) = 0
    res(1,2) = 0
    res(2,2) = c
    res(3,2) = s
    res(1,3) = 0
    res(2,3) = -s
    res(3,3) = c
  end function
  
  pure function rotation_matrix_y(phi) result(res)
    real(knd) :: res(3,3)
    real(knd),intent(in) :: phi
    real(knd) :: c,s
    
    c = cos(phi)
    s = sin(phi)
    res(1,1) = c
    res(2,1) = 0
    res(3,1) = -s
    res(1,2) = 0
    res(2,2) = 1
    res(3,2) = 0
    res(1,3) = s
    res(2,3) = 0
    res(3,3) = c
  end function
  
  pure function rotation_matrix_z(phi) result(res)
    real(knd) :: res(3,3)
    real(knd),intent(in) :: phi
    real(knd) :: c,s
    
    c = cos(phi)
    s = sin(phi)
    res(1,1) = c
    res(2,1) = s
    res(3,1) = 0
    res(1,2) = -s
    res(2,2) = c
    res(3,2) = 0
    res(1,3) = 0
    res(2,3) = 0
    res(3,3) = 1
  end function
  
  
  elemental function LinearTransform_Init_rot(original,axis,phi) result(res)
    type(LinearTransform) :: res
    class(GeometricShape),intent(in) :: original
    integer,intent(in) :: axis
    real(knd),intent(in) :: phi
  
    allocate(res%original, source=original)

    if (axis==1) then
      res%matrix = rotation_matrix_x(phi)
      res%inv_matrix = rotation_matrix_x(-phi)
    else if (axis==2) then
      res%matrix = rotation_matrix_y(phi)
      res%inv_matrix = rotation_matrix_y(-phi)
    else if (axis==3) then
      res%matrix = rotation_matrix_z(phi)
      res%inv_matrix = rotation_matrix_z(-phi)
    end if
  end function


  
  
  logical function Union_in_bbox(self,x,y,z,eps) result(in)
    class(Union),intent(in) :: self
    real(knd),intent(in) :: x,y,z,eps
    integer i
    
    in = any([ ( self%items(i)%in_bbox(x,y,z,eps), i=1,self%size ) ])

  end function


  logical function Union_Inside(self,x,y,z,eps) result(ins)
    class(Union),intent(in) :: self
    real(knd),intent(in) :: x,y,z
    real(knd),intent(in) ::eps
    integer i
    
    ins = any([ ( self%items(i)%InsideEps(x,y,z,eps), i=1,self%size ) ])
    
  end function
  
  
  subroutine Union_Closest(self,xnear,ynear,znear,x,y,z)
    class(Union),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: xs(self%size),ys(self%size),zs(self%size)
    integer i

    do i=1,self%size
      call self%items(i)%Closest(xs(i),ys(i),zs(i),x,y,z)
    end do
    
    associate (j => minloc(hypot(xs,hypot(ys,zs))))
      xnear = xs(j(1))
      ynear = ys(j(1))
      znear = zs(j(1))
    end associate
  end subroutine
  
  subroutine Union_ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(Union),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z
    real(knd) :: xs(self%size),ys(self%size),zs(self%size)
    integer i
    !FIXME: this will work assuming we are close enough to the boundary
    do i=1,self%size
      call self%items(i)%Closest(xs(i),ys(i),zs(i),x,y,z)
    end do
    
    associate (j => minloc(hypot(xs,hypot(ys,zs))))
      xnear = xs(j(1))
      ynear = ys(j(1))
      znear = zs(j(1))
    end associate
  end subroutine  
  
  logical function Union_IntersectsRay(self,r) result(res)
    class(Union),intent(in) :: self
    class(Ray),intent(in) :: r
    integer i
    
    res = any([ ( self%items(i)%IntersectsRay(Ray(r%xc,r%yc,r%zc,r%a,r%b,r%c)), &
                  i=1,self%size ) ])
  end function
  
  function Union_Init(items) result(res)
     type(Union) :: res
     class(GeometricShape),intent(in) :: items(:)
    
     allocate(res%items(size(items)), source=items)
     res%size = size(res%items)
  end function
  
  function Union_Init_Obst(filename) result(res)
    use Strings, only: upcase
    type(Union) :: res    
    character(*),intent(in) :: filename
    integer unit,io
    character(180) :: line
    type(ConvexPolyhedron) :: poly
    type(ConvexPolyhedron),allocatable :: items(:)

    open(newunit=unit,file=filename,action='read',status='old',iostat=io)
    
    allocate(items(0))

    if (io==0) then
      do
        read(unit,'(a)',iostat=io) line

        if (io/=0) exit

        line = adjustl(line)

        if (len_trim(line)>0) then

          if (upcase(line(1:10))=='POLYHEDRON') then
          
            call ReadPolyhedron(poly,line(11:))
            
           call add_element(items,poly)
!             items = [items, poly]
            
          end if

        end if

      end do

      close(unit)

    else

      write(*,*) "Could not open file ",filename
      stop

    end if
    
    call move_alloc(items,res%items)
    
    res%size = size(res%items)
    
    contains
    
      subroutine ReadPolyhedron(poly,restline)
        type(ConvexPolyhedron),intent(out) :: poly
        character(*),intent(in)  :: restline
        integer nPlanes,i,io

        read(restline,*,iostat=io) nPlanes

        if (io/=0) then
          write(*,*) "Expected number of planes in polyhedron, received '",trim(restline),"' instead."
          stop
        end if

        allocate(poly%Planes(nPlanes))

        poly%nplanes = nPlanes

        do i=1,nPlanes
          call ReadPlane(poly%Planes(i))
        end do

      end subroutine ReadPolyhedron


      subroutine ReadPlane(Pl)
        use Strings
        type(Plane),intent(out) :: Pl
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
      
      subroutine add_element(a,e)
        type(ConvexPolyhedron),allocatable,intent(inout) :: a(:)
        type(ConvexPolyhedron),intent(in) :: e
        type(ConvexPolyhedron),allocatable :: tmp(:)

        if (.not.allocated(a)) then
          a = [e]
        else
          call move_alloc(a,tmp)
          allocate(a(size(tmp)+1))
          a(1:size(tmp)) = tmp
          a(size(tmp)+1) = e
        end if
      end subroutine
  end function
  
  
  
  subroutine ArrayFromObst(res,filename)
    use Strings, only: upcase
    class(GeometricShape),allocatable,intent(out) :: res(:)
    character(*),intent(in) :: filename
    integer unit,io
    character(180) :: line
    type(ConvexPolyhedron) :: poly
    type(ConvexPolyhedron),allocatable :: items(:)

    open(newunit=unit,file=filename,action='read',status='old',iostat=io)
    
    allocate(items(0))

    if (io==0) then
      do
        read(unit,'(a)',iostat=io) line

        if (io/=0) exit

        line = adjustl(line)

        if (len_trim(line)>0) then

          if (upcase(line(1:10))=='POLYHEDRON') then
          
            call ReadPolyhedron(poly,line(11:))
            
            call add_element(items,poly)

!               items = [items, poly]

           end if
        end if

      end do

      close(unit)

    else

      write(*,*) "Could not open file ",filename
      stop

    end if
    
    call move_alloc(items,res)
    
    contains
    
      subroutine ReadPolyhedron(poly,restline)
        type(ConvexPolyhedron),intent(out) :: poly
        character(*),intent(in)  :: restline
        integer nPlanes,i,io

        read(restline,*,iostat=io) nPlanes

        if (io/=0) then
          write(*,*) "Expected number of planes in polyhedron, received '",trim(restline),"' instead."
          stop
        end if

        allocate(poly%Planes(nPlanes))

        poly%nplanes = nPlanes

        do i=1,nPlanes
          call ReadPlane(poly%Planes(i))
        end do

      end subroutine ReadPolyhedron


      subroutine ReadPlane(Pl)
        use Strings
        type(Plane),intent(out) :: Pl
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
      
      subroutine add_element(a,e)
        type(ConvexPolyhedron),allocatable,intent(inout) :: a(:)
        type(ConvexPolyhedron),intent(in) :: e
        type(ConvexPolyhedron),allocatable :: tmp(:)

        if (.not.allocated(a)) then
          a = [e]
        else
          call move_alloc(a,tmp)
          allocate(a(size(tmp)+1))
          a(1:size(tmp)) = tmp
          a(size(tmp)+1) = e
        end if
      end subroutine
  end subroutine

end module GeometricShapes







module Body_class
  use Kinds, only: knd
  use Lists, only: Listable
  use GeometricShapes, only: GeometricShape, Ray
  use Parameters
  
  implicit none

  private

  public Body,Inside


  type, extends(Listable),abstract :: Body
     integer numofbody
     class(GeometricShape),allocatable :: GeometricShape
  contains
     procedure :: Inside => CInside  !Hack around yet inidentified problem in GCC.
     procedure :: Closest
     procedure :: ClosestOut
     procedure :: ClosestOnLineOut
     procedure :: IntersectsRay
  end type
  interface Inside
    module procedure CInside
  end interface

contains


  logical function CInside(self,x,y,z,eps)
    class(Body),intent(in) :: self
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
    class(Body),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%GeometricShape%Closest(xnear,ynear,znear,x,y,z)

  end subroutine Closest

  

  subroutine ClosestOut(self,xnear,ynear,znear,x,y,z)
    class(Body),intent(in) :: self
    real(knd),intent(out) :: xnear,ynear,znear
    real(knd),intent(in) :: x,y,z

    call self%GeometricShape%ClosestOut(xnear,ynear,znear,x,y,z)

  end subroutine ClosestOut



  real(knd) function ClosestOnLineOut(self,x,y,z,x2,y2,z2) !Find t, such that x+(x2-x)*t lies on the boundary of the SB
    class(Body),intent(in) :: self
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

  logical function IntersectsRay(self,r) result(intersects)
    class(Body),intent(in) :: self
    class(Ray),intent(in) :: r

    intersects = self%GeometricShape%IntersectsRay(r)
 
  end function

end module Body_class







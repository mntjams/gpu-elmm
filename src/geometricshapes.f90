module GeometricShapes
  use Kinds
  use Parameters
  
  implicit none

  private

  public TGeometricShape, TLine, TPlane, TPolyhedron, TBall, TCylJacket, &
         TCylinder, TTerrainPoint, TTerrain

  type,abstract :: TGeometricShape
  contains
    procedure(Inside_interface),  deferred :: Inside
    procedure(Nearest_interface), deferred :: Nearest
    procedure(Nearest_interface), deferred :: NearestOut
  end type

  abstract interface
    pure logical function Inside_interface(self,x,y,z,eps)
      import
      class(TGeometricShape),intent(in) :: self
      real(KND),intent(in) :: x,y,z
      real(KND),optional,intent(in) ::eps
    end function
    subroutine Nearest_interface(self,xnear,ynear,znear,x,y,z)
      import
      class(TGeometricShape),intent(in) :: self
      real(KND),intent(out) :: xnear,ynear,znear
      real(KND),intent(in) :: x,y,z
    end subroutine
  end interface

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
!   contains
!     procedure Inside => TPlane_Inside
 !    procedure Nearest => TPlane_Nearest
  end type TPlane


  type,extends(TGeometricShape) :: TPolyhedron
    integer :: nplanes = 0
    type(TPlane),dimension(:),allocatable :: Planes !intersection of half-spaces
  contains
    procedure :: Inside => TPolyhedron_Inside
    procedure :: Nearest => TPolyhedron_Nearest
    procedure :: NearestOut => TPolyhedron_NearestOut
  end type TPolyhedron


  type,extends(TGeometricShape) :: TBall
    real(KND) xc,yc,zc,r
    logical :: rough = .false. !T rough surface, F flat surface
    real(KND) z0            !roughness parameter
  contains
    procedure :: Inside => TBall_Inside
    procedure :: Nearest => TBall_Nearest
    procedure :: NearestOut => TBall_NearestOut
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


  type,extends(TGeometricShape) :: TCylinder
    type(TCylJacket) Jacket
    type(TPlane),allocatable :: Plane1 ,Plane2
  contains
    procedure :: Inside => TCylinder_Inside
    procedure :: Nearest => TCylinder_Nearest
    procedure :: NearestOut => TCylinder_NearestOut
  end type TCylinder


  type TTerrainPoint
    real(KND) :: elev = 0
    logical :: rough = .false.
    real(KND) z0
  end type TTerrainPoint


   type,extends(TGeometricShape) :: TTerrain
     type(TTerrainPoint),dimension(:,:),allocatable :: UPoints,VPoints,PrPoints !allocate with a buffer of width 1 (i.e. 0:Xnx)
   contains
     procedure :: Inside => TTerrain_Inside
     procedure :: Nearest => TTerrain_Nearest
     procedure :: NearestOut => TTerrain_NearestOut
   end type TTerrain


  interface Inside
    module procedure TPlane_Inside
    module procedure TCylJacket_Inside
  end interface Inside

  interface Nearest   !shadows intrinsic nearest(), use intrinsic statement if needed
    module procedure TLine_Nearest
    module procedure TPlane_Nearest
    module procedure TCylJacket_Nearest
  end interface Nearest



contains

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



  pure logical function TPolyhedron_Inside(self,x,y,z,eps)
    class(TPolyhedron),intent(in) :: self
    real(KND),intent(in) :: x,y,z
    real(KND),optional,intent(in) ::eps
    logical ins
    integer i

    if (self%nplanes>0) then
      ins = .true.
      do i = 1,self%nplanes
       if (present(eps)) then
        ins = Inside(self%Planes(i),x,y,z,eps)
       else
        ins = Inside(self%Planes(i),x,y,z)
       endif
       if (.not. ins) exit
      enddo
    else
      ins = .false.
    endif

    TPolyhedron_Inside = ins
  end function TPolyhedron_Inside



  pure logical function TBall_Inside(self,x,y,z,eps)
    class(TBall),intent(in) :: self
    real(KND),intent(in) :: x,y,z
    real(KND),intent(in),optional ::eps
    real(KND) :: eps2
    logical ins

    if (present(eps)) then
     eps2 = eps
    else
     eps2 = MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND)
    endif

    if ((self%xc-x)**2+(self%yc-y)**2+(self%zc-z)**2<=(self%r+eps2)**2) then
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



  pure logical function TCylinder_Inside(self,x,y,z,eps)
   class(TCylinder),intent(in) :: self
   real(KND),intent(in) :: x,y,z
   real(KND),intent(in),optional ::eps
   logical ins

    ins = .true.

    if (.not.Inside(self%Jacket,x,y,z,eps)) ins = .false.

    if (ins.and.allocated(self%Plane1)) then
           if (.not.Inside(self%Plane1,x,y,z,eps)) ins = .false.
    endif
    if (ins.and.allocated(self%Plane2)) then
          if (.not.Inside(self%Plane2,x,y,z,eps)) ins = .false.
    endif

    TCylinder_Inside = ins
  end function TCylinder_Inside



  pure subroutine TTerrain_GridCoords(x2,y2,xi,yj,comp)
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



  pure logical function TTerrain_Inside(self,x,y,z,eps)
    class(TTerrain),intent(in) :: self
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
     if (z<=self%UPoints(xi,yj)%elev+eps2) ins = .true.
    elseif (comp==2) then
     if (z<=self%VPoints(xi,yj)%elev+eps2) ins = .true.
    elseif (comp==3) then
     if (z<=self%PrPoints(xi,yj)%elev+eps2) ins = .true.
    endif

    TTerrain_Inside = ins
  end function TTerrain_Inside


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



  subroutine TCylinder_Nearest(self,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
   class(TCylinder),intent(in) :: self
   real(KND),intent(out) :: xnear,ynear,znear
   real(KND),intent(in) :: x,y,z
   real(KND) xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2

   !!!Only for Planes perpendicular to jacket!!!!

   if (allocated(self%Plane1)) then
       if (allocated(self%Plane2)) then
          call Nearest(self%Jacket,xJ,yJ,zJ,x,y,z)
          call Nearest(self%Plane1,xP1,yP1,zP1,x,y,z)
          call Nearest(self%Plane2,xP2,yP2,zP2,x,y,z)
          if (Inside(self%Jacket,x,y,z)) then
             if (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
                xnear = xP1
                ynear = yP1
                znear = zP1
             else
                xnear = xP2
                ynear = yP2
                znear = zP2
             endif
          elseif (Inside(self%Plane1,x,y,z).and.Inside(self%Plane2,x,y,z)) then
                xnear = xJ
                ynear = yJ
                znear = zJ
          elseif (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
           call Nearest(self%Jacket,xnear,ynear,znear,xP1,yP1,zP1)
          else
           call Nearest(self%Jacket,xnear,ynear,znear,xP2,yP2,zP2)
          endif
       else
          call Nearest(self%Jacket,xJ,yJ,zJ,x,y,z)
          call Nearest(self%Plane1,xP1,yP1,zP1,x,y,z)
          if (Inside(self%Jacket,x,y,z)) then
                xnear = xP1
                ynear = yP1
                znear = zP1
          elseif (Inside(self%Plane1,x,y,z)) then
                xnear = xJ
                ynear = yJ
                znear = zJ
          else
           call Nearest(self%Jacket,xnear,ynear,znear,xP1,yP1,zP1)
         endif
       endif

    else

      call Nearest(self%Jacket,xnear,ynear,znear,x,y,z)

    endif
  end subroutine TCylinder_Nearest

  subroutine TBall_Nearest(self,xnear,ynear,znear,x,y,z)
    class(TBall),intent(in) :: self
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    real(KND) t,a,b,c

    a = x-self%xc
    b = y-self%yc
    c = z-self%zc
    t = self%r/SQRT(a**2+b**2+c**2)
    xnear = a*t+self%xc
    ynear = b*t+self%yc
    znear = c*t+self%zc
  end subroutine TBall_Nearest


  subroutine TPolyhedron_Nearest(self,xnear,ynear,znear,x,y,z)
    use Lapack
    class(TPolyhedron),intent(in) :: self
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    real(KND) :: dists(self%nplanes),xP(self%nplanes),yP(self%nplanes),zP(self%nplanes)
    real(KND) :: minv
    real(KND) :: ailine,biline,ciline
    real(KND) :: x0iline,y0iline,z0iline,xln,yln,zln,p
    integer   :: inearest,inearest2,inearest3
    integer   :: i

    integer,dimension(3)    :: ipivot,work2             !arguments of the LAPACK
    real(KND),dimension(3)  :: xg,bg,R,C,ferr,berr      !  routine GESVX
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

    do i = 1,self%nplanes
     call Nearest(self%Planes(i),xP(i),yP(i),zP(i),x,y,z)
     dists(i) = SQRT((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
     if (Inside(self%Planes(i),x,y,z)) dists(i) = -ABS(dists(i))
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
                    MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND) &
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

     call Nearest(self%Planes(inearest),xln,yln,zln,x,y,z)
     call Nearest(self%Planes(inearest2),xnear,ynear,znear,x,y,z)

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

       ag = reshape(source = (/ self%Planes(inearest)%a,self%Planes(inearest2)%a,0._KND, &
                           self%Planes(inearest)%b,self%Planes(inearest2)%b,0._KND, &
                           self%Planes(inearest)%c,self%Planes(inearest2)%c,1._KND /), &
                 shape = (/ 3,3 /))

       bg = (/ -self%Planes(inearest)%d,-self%Planes(inearest2)%d,(zW(Wnz+1)+zW(0))/2._KND /)

       call gesvx("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3, &
                  rcond,ferr,berr,work,work2,info)

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    elseif (abs(biline)>=0.1_KND) then

       ag = reshape(source = (/ self%Planes(inearest)%a,self%Planes(inearest2)%a,0._KND, &
                           self%Planes(inearest)%b,self%Planes(inearest2)%b,1._KND, &
                           self%Planes(inearest)%c,self%Planes(inearest2)%c,0._KND /), &
                 shape = (/ 3,3 /))

       bg = (/ -self%Planes(inearest)%d,-self%Planes(inearest2)%d,(yV(Vny+1)+yV(0))/2._KND /)

       call gesvx("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3, &
                  rcond,ferr,berr,work,work2,info)

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    else

       ag = reshape(source = (/ self%Planes(inearest)%a,self%Planes(inearest2)%a,1._KND, &
                           self%Planes(inearest)%b,self%Planes(inearest2)%b,0._KND, &
                           self%Planes(inearest)%c,self%Planes(inearest2)%c,0._KND /), &
                  shape = (/ 3,3 /))

       bg = (/ -self%Planes(inearest)%d,-self%Planes(inearest2)%d,(xU(Unx+1)+xU(0))/2._KND /)

       call gesvx("E","N",3,1,ag,3,af,3,ipivot,EQUED,R,C,bg,3,xg,3, &
                  rcond,ferr,berr,work,work2,info)

       x0iline = xg(1)
       y0iline = xg(2)
       z0iline = xg(3)

    endif

    call Nearest(TLine(x0iline,y0iline,z0iline,ailine,biline,ciline),xln,yln,zln,x,y,z)

    if (self%Inside(xln,yln,zln,min(dxmin/1000._KND,dymin/10000._KND,dzmin/10000._KND))) then
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

  end subroutine TPolyhedron_Nearest




  subroutine TTerrain_Nearest(self,xnear,ynear,znear,x,y,z)
    class(TTerrain),intent(in) :: self
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

     zloc = self%UPoints(xi,yj)%elev
     a = (self%PrPoints(xi+1,yj)%elev - self%PrPoints(xi,yj)%elev) / dxU(xi)
     b = (self%UPoints(xi,yj+1)%elev - self%UPoints(xi,yj-1)%elev) / (yPr(yj+1)-yPr(yj-1))

     Pl%a = a
     Pl%b = b
     Pl%c = -1._KND
     Pl%d = -a*x-b*y+zloc

     call Nearest(Pl,xnear,ynear,znear,x,y,z)

    elseif (comp==2) then

     zloc = self%VPoints(xi,yj)%elev
     a = (self%VPoints(xi+1,yj)%elev - self%VPoints(xi-1,yj)%elev) / (xPr(xi+1)-xPr(xi-1))
     b = (self%PrPoints(xi,yj+1)%elev - self%PrPoints(xi,yj)%elev) / dyV(yj)

     Pl%a = a
     Pl%b = b
     Pl%c = -1._KND
     Pl%d = -a*x-b*y+zloc

     call Nearest(Pl,xnear,ynear,znear,x,y,z)

    elseif (comp==3) then

     zloc = self%PrPoints(xi,yj)%elev
     a = (self%UPoints(xi,yj)%elev - self%UPoints(xi-1,yj)%elev) / dxPr(xi)
     b = (self%VPoints(xi,yj)%elev - self%VPoints(xi,yj-1)%elev) / dyPr(yj)

     Pl%a = a
     Pl%b = b
     Pl%c = -1._KND
     Pl%d = - a*x - b*y + zloc

     call Nearest(Pl,xnear,ynear,znear,x,y,z)

    endif
  end subroutine TTerrain_Nearest





  subroutine TPolyhedron_NearestOut(self,xnear,ynear,znear,x,y,z)
    class(TPolyhedron),intent(in) :: self
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    real(KND) dists(self%nplanes),xP(self%nplanes),yP(self%nplanes),zP(self%nplanes),minv
    integer inearest,i

    dists = huge(minv)
    do i = 1,self%nplanes
     call Nearest(self%Planes(i),xP(i),yP(i),zP(i),x,y,z)
     dists(i) = SQRT((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
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
  end subroutine TPolyhedron_NearestOut



  subroutine TCylinder_NearestOut(self,xnear,ynear,znear,x,y,z) !only for planes perpendicular to the axis
    class(TCylinder),intent(in) :: self
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z
    real(KND) xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2

    if (allocated(self%Plane1)) then
      call Nearest(self%Plane1,xP1,yP1,zP1,x,y,z)
    else
     xP1 = sqrt(huge(xP1))/10
     yP1 = sqrt(huge(yP1))/10
     zP1 = sqrt(huge(zP1))/10
    endif
    if (allocated(self%Plane2)) then
      call Nearest(self%Plane2,xP2,yP2,zP2,x,y,z)
    else
     xP2 = sqrt(huge(xP2))/10
     yP2 = sqrt(huge(yP2))/10
     zP2 = sqrt(huge(zP2))/10
    endif
    call Nearest(self%Jacket,xJ,yJ,zJ,x,y,z)

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



  subroutine TBall_NearestOut(self,xnear,ynear,znear,x,y,z)
    class(TBall),intent(in) :: self
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    call self%Nearest(xnear,ynear,znear,x,y,z)
  end subroutine TBall_NearestOut



  subroutine TTerrain_NearestOut(self,xnear,ynear,znear,x,y,z)
    class(TTerrain),intent(in) :: self
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    call self%Nearest(xnear,ynear,znear,x,y,z)
  end subroutine TTerrain_NearestOut




end module GeometricShapes







module TBody_class
  use Kinds, only: KND
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
     procedure :: Nearest
     procedure :: NearestOut
     procedure :: NearestOnLineOut
  end type
  interface Inside
    module procedure CInside
  end interface

contains


  pure logical function CInside(self,x,y,z,eps)
    class(TBody),intent(in) :: self
    real(KND),intent(in) :: x,y,z
    real(KND),intent(in),optional :: eps
    real(KND) x2,y2,z2
    real(KND) lx,ly,lz

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

      if (present(eps)) then
        CInside = self%GeometricShape%Inside(x2,y2,z2,eps)
      else
        CInside = self%GeometricShape%Inside(x2,y2,z2)
      endif
    end if

  end function CInside

  

  subroutine Nearest(self,xnear,ynear,znear,x,y,z)
    class(TBody),intent(in) :: self
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    call self%GeometricShape%Nearest(xnear,ynear,znear,x,y,z)

  end subroutine Nearest

  

  subroutine NearestOut(self,xnear,ynear,znear,x,y,z)
    class(TBody),intent(in) :: self
    real(KND),intent(out) :: xnear,ynear,znear
    real(KND),intent(in) :: x,y,z

    call self%GeometricShape%NearestOut(xnear,ynear,znear,x,y,z)

  end subroutine NearestOut



  real(KND) function NearestOnLineOut(self,x,y,z,x2,y2,z2) !Find t, such that x+(x2-x)*t lies on the boundary of the SB
    class(TBody),intent(in) :: self
    real(KND),intent(in) :: x,y,z,x2,y2,z2

    real(KND) t,t1,t2
    integer i

    t1 = 0
    t2 = 1
    if (self%Inside(x2,y2,z2)) then
     do
      t2 = t2*2._KND
      if (.not. self%Inside(x+(x2-x)*t2,y+(y2-y)*t2,z+(z2-z)*t2)) exit
     enddo
    endif
    t = (t1+t2)/2._KND

    do i = 1,20         !The bisection method with maximum 20 iterations (should be well enough)
     if (self%Inside(x+(x2-x)*t,y+(y2-y)*t,z+(z2-z)*t)) then
      t1 = t
     else
      t2 = t
     endif
     t = (t1+t2)/2._KND
     if (abs(t1-t2)<MIN(dxmin/1000._KND,dymin/1000._KND,dzmin/1000._KND))   exit
    enddo
    NearestOnLineOut = t
  end function NearestOnLineOut


end module TBody_class







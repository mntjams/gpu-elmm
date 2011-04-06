module GEOMETRIC
use PARAMETERS
use WALLMODELS

implicit none


 TYPE TLine
   real(KND) xc,yc,zc
   real(KND) a,b,c
 ENDTYPE TLine

 TYPE TPlane
   real(KND) a,b,c,d !ax+by+cz+d/=0 for inner half-space
   logical gl !T > in ineq. above F < in ineq. above
   logical rough !T rough surface, F flat surface
   real(KND) z0 !roughness parameter
  ENDTYPE TPlane

  TYPE TPolyhedron
   integer nplanes
   TYPE(TPlane),dimension(:),allocatable:: Planes !intersection of half-spaces
  ENDTYPE TPolyhedron

  TYPE TBall
   real(KND) xc,yc,zc,r
   logical rough !T rough surface, F flat surface
   real(KND) z0 !roughness parameter
  ENDTYPE TBall

  TYPE TCylJacket
   real(KND) xc,yc,zc
   real(KND) a,b,c
   real(KND) r
   logical rough !T rough surface, F flat surface
   real(KND) z0 !roughness parameter
  ENDTYPE TCylJacket

  TYPE TCylinder
   TYPE(TCylJacket) Jacket
   TYPE(TPlane),pointer:: Plane1 => null() ,Plane2 => null()
  ENDTYPE TCylinder


 TYPE TTerrainPoint
  real(KND) elev
  logical:: rough=.false.
  real(KND) z0
 ENDTYPE TTerrainPoint


  TYPE TTerrain
   TYPE(TTerrainPoint),DIMENSION(:,:),ALLOCATABLE:: UTerrPoints,VTerrPoints,PrTerrPoints
  ENDTYPE TTerrain

 TYPE TSolidBody
  integer numofbody
  integer:: typeofbody=0 !0.. none, 1..polyhedron, 2.. ball, 3.. cylinder, 4.. terrain
  TYPE(TPolyhedron),pointer:: Polyhedron => null() !asociated will be only part writen in typeofbody
  TYPE(TBall),pointer:: Ball => null()
  TYPE(TCylinder),pointer:: Cylinder => null()
  TYPE(TTerrain),pointer:: Terrain => null()
  logical:: rough=.false. !T rough surface, F flat surface
  real(KND) z0 !roughness parameter
  TYPE(TSolidBody),pointer:: next =>null()
 ENDTYPE TSolidbody

 TYPE(TSolidBody),pointer:: FirstSB => null() !First member of the linked list of solid objects
  
 


 TYPE TIBPoint
  integer component !1..U, 2..V, 3..W
  integer x
  integer y
  integer z
  real(KND) distx !coords of nearest boundary point (in the mesh units!)
  real(KND) disty
  real(KND) distz
  integer dirx
  integer diry
  integer dirz
  integer interp !kind of interp. 0.. none (boundarypint), 1..linear, 2..bilinear, 3..trilinear
  integer interpdir
  real(KND) MSourc
  TYPE(TIBPoint),pointer:: next => null()
 ENDTYPE TIBPoint

 TYPE TScalFlIBPoint
  integer x
  integer y
  integer z
  real(KND) dist
  integer,dimension(4):: intpointi
  integer,dimension(4):: intpointj
  integer,dimension(4):: intpointk
  real(KND),dimension(4):: intcoef
  integer interp !kind of interp. 1.. none (1 point outside), 2..linear, 4..bilinear  other values not allowed
  real(KND) ScalSourc !virtual scalar source
  real(KND):: Flux=0 !desired scalar flux
  type(TScalFlIBPoint),pointer::next=>null()
 ENDTYPE TScalFlIBPoint


 type(TIBPoint),pointer::FirstIBPoint=>null(), LastIBPoint=>null()
 type(TScalFlIBPoint),pointer::FirstScalFlIBPoint=>null(), LastScalFlIBPoint=>null()


contains

  subroutine SetCurrentSB(SB,n)
  type(TSolidBody),pointer:: SB
  integer(SINT),intent(IN):: n
  integer i

   SB=>FirstSB
   do i=2,n
    if (associated(SB%next)) then
     SB=>SB%next
    else
     SB=>null()
     EXIT
    endif
   enddo
  endsubroutine SetCurrentSB
  
  pure logical function PlaneInt(x,y,z,PL,eps)
   real(KND),intent(IN):: x,y,z
   TYPE(TPlane),intent(IN):: PL
   real(KND),optional,intent(IN)::eps
   real(KND) eps2
   logical interior

   if (PRESENT(eps)) then
    eps2=ABS(eps)
   else
    eps2=0
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
   PlaneInt=interior
  endfunction PlaneInt

  pure logical function PolyhedronInt(x,y,z,PH,eps)
   real(KND),intent(IN):: x,y,z
   TYPE(TPolyhedron),intent(IN):: PH
   real(KND),optional,intent(IN)::eps
   logical interior
   integer i

   if (PH%nplanes>0) then
    interior=.true.
    do i=1,PH%nplanes
     if (PRESENT(eps)) then
      interior=PlaneInt(x,y,z,PH%Planes(i),eps)
     else
      interior=PlaneInt(x,y,z,PH%Planes(i))
     endif
     if (.not. INTERIOR) EXIT
    enddo
   else
    interior=.false.
   endif
   PolyhedronInt=interior
  endfunction PolyhedronInt

  pure logical function BallInt(x,y,z,B)
   real(KND),intent(IN):: x,y,z
   TYPE(TBall),intent(IN):: B
   logical interior

   if ((B%xc-x)**2+(B%yc-y)**2+(B%zc-z)**2<=(B%r)**2) then
    interior=.true.
   else
    interior=.false.
   endif
   BallInt=interior
  endfunction BallInt

  pure real(KND) function LineDist(x,y,z,xl,yl,zl,a,b,c)
   real(KND),intent(IN):: x,y,z,xl,yl,zl,a,b,c
   real(KND) t
   if (((a/=0).or.(b/=0)).or.(c/=0)) then
    t=(a*(x-xl)+b*(y-yl)+c*(z-zl))/(a**2+b**2+c**2)
   else
    t=0
   endif
   LineDist=SQRT((xl+a*t-x)**2+(yl+b*t-y)**2+(zl+c*t-z)**2)
  endfunction LineDist

  pure logical function JacketInt(x,y,z,J)
   real(KND),intent(IN):: x,y,z
   TYPE(TCylJacket),intent(IN):: J
   real(KND) eps
   logical interior
! write (*,*) x
! write (*,*)y
! write (*,*)z
! write (*,*)j%xc
! write (*,*)j%yc
! write (*,*)j%zc
! write (*,*)J%a
! write (*,*)J%b
! write (*,*)J%c
! write (*,*)j%r


   eps=MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND)
! write (*,*)eps
   if (LineDist(x,y,z,j%xc,j%yc,j%zc,J%a,J%b,J%c)<=j%r+eps) then
    interior=.true.
   else
    interior=.false.
   endif
   JacketInt=interior
  endfunction JacketInt

  pure logical function CylinderInt(x,y,z,C)
   real(KND),intent(IN):: x,y,z
   TYPE(TCylinder),intent(IN):: C
   logical interior
    interior=.true.
    if (.not.JacketInt(x,y,z,C%Jacket)) interior=.false.
    if (interior.and.associated(C%Plane1)) then
           if (.not.PlaneInt(x,y,z,C%Plane1)) interior=.false.
    endif
    if (interior.and.associated(C%Plane2)) then
          if (.not.PlaneInt(x,y,z,C%Plane2)) interior=.false.
    endif
    CylinderInt=interior
  endfunction CylinderInt



  pure subroutine TerrGridCoords(x2,y2,xi,yj,comp)
  real(KND),intent(IN):: x2,y2
  integer,intent(OUT):: xi,yj,comp
  real(KND) x,y,distPr,distU,distV
  integer xPri,yPrj,xUi,yVj,i

   x=x2
   y=y2

   xPri=Prnx
   do i=1,Prnx
    if (xU(i-1)>=x) then
                  xPri=i
                  EXIT
                 endif
   enddo

   yPrj=Prny
   do i=1,Prny
    if (yV(i-1)>=y) then
                  yPrj=i
                  EXIT
                 endif
   enddo

   xUi=Prnx+1
   do i=0,Prnx+1
    if (xPr(i)>=x) then
                  xUi=i
                  EXIT
                 endif
   enddo

   yVj=Prny+1
   do i=0,Prny+1
    if (yPr(i)>=y) then
                  yVj=i
                  EXIT
                 endif
   enddo

   distPr=(x-xPr(xPri))**2+(y-yPr(yPrj))**2
   distU=(x-xU(xUi))**2+(y-yPr(yPrj))**2
   distV=(x-xPr(xPri))**2+(y-yV(yVj))**2

   if (distU<distPr.and.distU<distV) then
    xi=xUi
    yj=yPrj
    comp=1
   elseif (distV<distPr.and.distV<distU) then
    xi=xPri
    yj=yVj
    comp=2
   else
    xi=xPri
    yj=yPrj
    comp=3
   endif
   endsubroutine TerrGridCoords




  pure logical function TerrainInt(x,y,z,T)
   real(KND),intent(IN):: x,y,z
   TYPE(TTerrain),intent(IN):: T
   logical interior
   integer xi,yj,comp
   interior=.false.
   call TerrGridCoords(x,y,xi,yj,comp)
   if (comp==1) then
    if (z<=T%UTerrPoints(xi,yj)%elev) interior=.true.
   elseif (comp==2) then
    if (z<=T%VTerrPoints(xi,yj)%elev) interior=.true.
   elseif (comp==3) then
    if (z<=T%PrTerrPoints(xi,yj)%elev) interior=.true.
   endif

   TerrainInt=interior
  endfunction TerrainInt

    
  pure logical function SolidBodyInt(x,y,z,SB)
   real(KND),intent(IN):: x,y,z
   real(KND) x2,y2,z2
   TYPE(TSolidBody),intent(IN):: SB

   x2=x
   y2=y
   z2=z
   if (BtypeE==PERIODIC.and.x2>xU(Prnx+1)) x2=x2-lx
   if (BtypeN==PERIODIC.and.y2>yV(Prny+1)) y2=y2-ly
   if (BtypeT==PERIODIC.and.z2>zW(Prnz+1)) z2=z2-lz

   if (BtypeW==PERIODIC.and.x2<xU(0)) x2=x2+lx
   if (BtypeS==PERIODIC.and.y2<yV(0)) y2=y2+ly
   if (BtypeB==PERIODIC.and.z2<zW(0)) z2=z2+lz

   select case (SB%typeofbody)
    case (1)
     SolidBodyInt=PolyhedronInt(x2,y2,z2,SB%Polyhedron)
    case (2)
     SolidBodyInt=BallInt(x2,y2,z2,SB%Ball)
    case (3)
     SolidBodyInt=CylinderInt(x2,y2,z2,SB%Cylinder)
    case (4)
     SolidBodyInt=TerrainInt(x2,y2,z2,SB%Terrain)
    case default
     SolidBodyInt=.false.
    end select
  endfunction SolidBodyInt





  real(KND) function PointDist(x1,y1,z1,x2,y2,z2)
  real(KND) x1,y1,z1,x2,y2,z2
   PointDist=SQRT((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
  endfunction PointDist


  subroutine LineNearest(xnear,ynear,znear,x,y,z,xl,yl,zl,a,b,c)
  real(KND) xnear,ynear,znear,x,y,z,xl,yl,zl,a,b,c
  real(KND) t
   if (((a/=0).or.(b/=0)).or.(c/=0)) then
    t=(a*(x-xl)+b*(y-yl)+c*(z-zl))/(a**2+b**2+c**2)
   else
    t=0
   endif
   xnear=xl+a*t
   ynear=yl+b*t
   znear=zl+c*t
  endsubroutine LineNearest

  subroutine PlaneNearest(xnear,ynear,znear,x,y,z,PL)
   real(KND) xnear,ynear,znear,x,y,z
   TYPE(TPlane):: PL
   real(KND) t
   
   if (((PL%a/=0).or.(PL%b/=0)).or.(PL%c/=0)) then
    t=-(PL%a*x+PL%b*y+PL%c*z+PL%d)/(PL%a**2+PL%b**2+PL%c**2)
   else
    t=0
   endif
   xnear=x+PL%a*t
   ynear=y+PL%b*t
   znear=z+PL%c*t
  endsubroutine PlaneNearest



  subroutine JacketNearest(xnear,ynear,znear,x,y,z,J)
   real(KND) xnear,ynear,znear,x,y,z
   TYPE(TCylJacket):: J
   real(KND) t,xl,yl,zl,a,b,c

   call LineNearest(xl,yl,zl,x,y,z,J%xc,J%yc,J%zc,J%a,J%b,J%c)
   a=x-xl
   b=y-yl
   c=z-zl
   t=J%r/SQRT(a**2+b**2+c**2)
   xnear=a*t+xl
   ynear=b*t+yl
   znear=c*t+zl
  endsubroutine JacketNearest



  subroutine CylinderNearest(xnear,ynear,znear,x,y,z,C) !only for planes perpendicular to the axis
   real(KND) xnear,ynear,znear,x,y,z
   TYPE(TCylinder):: C
   real(KND) xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2
  
   !!!Only for Planes perpendicular to jacket!!!!

   if (associated(C%Plane1)) then
       if (associated(C%Plane2)) then
          call JacketNearest(xJ,yJ,zJ,x,y,z,C%Jacket)
          call PlaneNearest(xP1,yP1,zP1,x,y,z,C%Plane1)
          call PlaneNearest(xP2,yP2,zP2,x,y,z,C%Plane2)
          if (JacketInt(x,y,z,C%Jacket)) then
             if (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
                xnear=xP1
                ynear=yP1
                znear=zP1
             else
                xnear=xP2
                ynear=yP2
                znear=zP2
             endif
          elseif (PlaneInt(x,y,z,C%Plane1).and.PlaneInt(x,y,z,C%Plane2)) then
                xnear=xJ
                ynear=yJ
                znear=zJ
          elseif (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2)) then
           call JacketNearest(xnear,ynear,znear,xP1,yP1,zP1,C%Jacket)
          else
           call JacketNearest(xnear,ynear,znear,xP2,yP2,zP2,C%Jacket)
          endif
       else
          call JacketNearest(xJ,yJ,zJ,x,y,z,C%Jacket)
          call PlaneNearest(xP1,yP1,zP1,x,y,z,C%Plane1)
          if (JacketInt(x,y,z,C%Jacket)) then
                xnear=xP1
                ynear=yP1
                znear=zP1
          elseif (PlaneInt(x,y,z,C%Plane1)) then
                xnear=xJ
                ynear=yJ
                znear=zJ
          else
           call JacketNearest(xnear,ynear,znear,xP1,yP1,zP1,C%Jacket)
         endif
       endif

   else
    call JacketNearest(xnear,ynear,znear,x,y,z,C%Jacket)
   endif
  endsubroutine CylinderNearest

  subroutine BallNearest(xnear,ynear,znear,x,y,z,Bl)
   real(KND) xnear,ynear,znear,x,y,z
   TYPE(TBall):: Bl
   real(KND) t,a,b,c

   a=x-Bl%xc
   b=y-Bl%yc
   c=z-Bl%zc
   t=Bl%r/SQRT(a**2+b**2+c**2)
   xnear=a*t+Bl%xc
   ynear=b*t+Bl%yc
   znear=c*t+Bl%zc
  endsubroutine BallNearest


  subroutine PolyhedronNearest(xnear,ynear,znear,x,y,z,PH)
  real(KND) xnear,ynear,znear,x,y,z
  type(TPolyhedron) PH
  real(KND) dists(PH%nplanes),xP(PH%nplanes),yP(PH%nplanes),zP(PH%nplanes),minv
  real(KND) ailine,biline,ciline,x0iline,y0iline,z0iline,ag(3,3),bg(3),xg(3),xln,yln,zln,c
  integer nearest,nearest2,nearest3,i,j,INDX(3)
  !Vzdalenesti od rovin, pokud je nejbl. bod roviny uvnitr jine, nebo puv. bod na vnitrni strane -> vzd. *-1
  !Pokud je nejbl. body +- eps. (norm. vektor!,dxmin/100) uvnitr a vne polyh -> hotovo
  !Jinak 2. nejbl. rovina v abs. hodnote -> prusecnice a nejbl bod na ni
  !Nejbl. bod na prusecnici. Pokud +-eps. uvnitr, (najit vekt. v rovine  kolme na prusecnic)-:hotovo
  !Jinak iterativne najit bod na prusecnici uvnitr



   do i=1,PH%nplanes
    call PlaneNearest(xP(i),yP(i),zP(i),x,y,z,PH%Planes(i))
    dists(i)=SQRT((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
    if (PlaneInt(x,y,z,PH%Planes(i))) dists(i)=-ABS(dists(i))
   enddo
   !find nearest plane with
   write (*,*) "dists",dists(:)
   nearest=0
   minv=huge(minv)
   do i=1,PH%nplanes
    if (dists(i)>=0.and.dists(i)<minv) then
      nearest=i
      minv=dists(i)
    endif
   enddo

   if (nearest==0) then
     write(*,*) "no nearest"
     goto 10
   endif

   if (PolyhedronInt(xP(nearest),yP(nearest),zP(nearest),PH,MIN(dxmin/10000._KND,dymin/10000._KND,dzmin/10000._KND))) then
    xnear=xP(nearest)
    ynear=yP(nearest)
    znear=zP(nearest)
    goto 10
   endif

   dists=abs(dists)

   nearest2=0
   minv=huge(minv)
   do i=1,PH%nplanes
    if (i/=nearest.and.dists(i)<minv) then
      nearest2=i
      minv=dists(i)
    endif
   enddo

   nearest3=0
   minv=huge(minv)
   do i=1,PH%nplanes
    if (i/=nearest.and.i/=nearest2.and.dists(i)<minv) then
      nearest3=i
      minv=dists(i)
    endif
   enddo

   ailine=PH%Planes(nearest)%b*PH%Planes(nearest2)%c-PH%Planes(nearest)%c*PH%Planes(nearest2)%b
   biline=PH%Planes(nearest)%c*PH%Planes(nearest2)%a-PH%Planes(nearest)%a*PH%Planes(nearest2)%c
   ciline=PH%Planes(nearest)%a*PH%Planes(nearest2)%b-PH%Planes(nearest)%b*PH%Planes(nearest2)%a
 
   if (abs(ailine)<=epsilon(ailine).and.abs(biline)<=epsilon(biline).and.abs(ciline)<=epsilon(ciline)) then
    write(*,*) PH%Planes(nearest)%c*PH%Planes(nearest2)%a,PH%Planes(nearest)%a*PH%Planes(nearest2)%c
    write(*,*) "cross product 0"
    write(*,*) "pl1",PH%Planes(nearest)%a,PH%Planes(nearest)%b,PH%Planes(nearest)%c,PH%Planes(nearest)%d
    write(*,*) "pl2",PH%Planes(nearest2)%a,PH%Planes(nearest2)%b,PH%Planes(nearest2)%c,PH%Planes(nearest2)%d
   endif

   c=SQRT(ailine**2+biline**2+ciline**2)
   ailine=ailine/c
   biline=biline/c
   ciline=ciline/c

   if (abs(ciline)>=0.1_KND) then
      ag=reshape(source=(/ PH%Planes(nearest)%a,PH%Planes(nearest2)%a,0._KND,&
                          PH%Planes(nearest)%b,PH%Planes(nearest2)%b,0._KND,&
                          PH%Planes(nearest)%c,PH%Planes(nearest2)%c,1._KND /),&
                shape=(/ 3,3 /))
      bg=(/ -PH%Planes(nearest)%d,-PH%Planes(nearest2)%d,(zW(Wnz+1)+zW(0))/2._KND /)

      call LEGS(ag,3,bg,xg,INDX)

      x0iline=xg(1)
      y0iline=xg(2)
      z0iline=xg(3)
   elseif (abs(biline)>=0.1_KND) then
      ag=reshape(source=(/ PH%Planes(nearest)%a,PH%Planes(nearest2)%a,0._KND,&
                          PH%Planes(nearest)%b,PH%Planes(nearest2)%b,1._KND,&
                          PH%Planes(nearest)%c,PH%Planes(nearest2)%c,0._KND /),&
                shape=(/ 3,3 /))
      bg=(/ -PH%Planes(nearest)%d,-PH%Planes(nearest2)%d,(yV(Vny+1)+yV(0))/2._KND /)

      call LEGS(ag,3,bg,xg,INDX)

      x0iline=xg(1)
      y0iline=xg(2)
      z0iline=xg(3)

   else
      ag=reshape(source=(/ PH%Planes(nearest)%a,PH%Planes(nearest2)%a,1._KND,&
                          PH%Planes(nearest)%b,PH%Planes(nearest2)%b,0._KND,&
                          PH%Planes(nearest)%c,PH%Planes(nearest2)%c,0._KND /),&
                 shape=(/ 3,3 /))
      bg=(/ -PH%Planes(nearest)%d,-PH%Planes(nearest2)%d,(xU(Unx+1)+xU(0))/2._KND /)

      call LEGS(ag,3,bg,xg,INDX)

       x0iline=xg(1)
       y0iline=xg(2)
       z0iline=xg(3)
   endif
   call LineNearest(xln,yln,zln,x,y,z,x0iline,y0iline,z0iline,ailine,biline,ciline)
   if (PolyhedronInt(xln,yln,zln,PH,min(dxmin/1000._KND,dymin/10000._KND,dzmin/10000._KND))) then
    xnear=xln
    ynear=yln
    znear=zln
    goto 10
   endif
  


     ag=reshape(source=(/ PH%Planes(nearest)%a,PH%Planes(nearest2)%a,PH%Planes(nearest3)%a,&
                         PH%Planes(nearest)%b,PH%Planes(nearest2)%b,PH%Planes(nearest3)%b,&
                         PH%Planes(nearest)%c,PH%Planes(nearest2)%c,PH%Planes(nearest3)%c /),&
                shape=(/ 3,3 /))
     bg=(/ -PH%Planes(nearest)%d,-PH%Planes(nearest2)%d,-PH%Planes(nearest3)%d /)

     call LEGS(ag,3,bg,xg,INDX)
       xnear=xg(1)
       ynear=xg(2)
       znear=xg(3)


   10 i=1
  endsubroutine PolyhedronNearest


   

  subroutine TerrainNearest(xnear,ynear,znear,x,y,z,T)
   real(KND) xnear,ynear,znear
   real(KND) x,y,z
   TYPE(TTerrain):: T
   integer xi,yj,comp

   xnear=x
   ynear=y
   znear=z0
   call TerrGridCoords(x,y,xi,yj,comp)
   if (comp==1) znear=T%UTerrPoints(xi,yj)%elev
   if (comp==2) znear=T%VTerrPoints(xi,yj)%elev
   if (comp==3) znear=T%PrTerrPoints(xi,yj)%elev
  endsubroutine TerrainNearest


  subroutine SolidBodyNearest(xnear,ynear,znear,x,y,z,SB)
   real(KND),INTENT(OUT):: xnear,ynear,znear
   real(KND) x,y,z
   TYPE(TSolidBody):: SB
   
   select case (SB%typeofbody)
    case (1)
     call PolyhedronNearest(xnear,ynear,znear,x,y,z,SB%Polyhedron)
    case (2)
     call BallNearest(xnear,ynear,znear,x,y,z,SB%Ball)
    case (3)
     call CylinderNearest(xnear,ynear,znear,x,y,z,SB%Cylinder)
    case (4)
     call TerrainNearest(xnear,ynear,znear,x,y,z,SB%Terrain)
    case default
     xnear=-huge(znear);ynear=-huge(znear);znear=-huge(znear)
   end select
  endsubroutine SolidBodyNearest



  subroutine CylinderNearestOut(xnear,ynear,znear,x,y,z,C) !only for planes perpendicular to the axis
   real(KND) xnear,ynear,znear,x,y,z
   TYPE(TCylinder):: C
   real(KND) xJ,yJ,zJ,xP1,yP1,zP1,xP2,yP2,zP2
  
   if (associated(C%Plane1)) then
     call PlaneNearest(xP1,yP1,zP1,x,y,z,C%Plane1)
   else
    xP1=sqrt(huge(xP1))/10
    yP1=sqrt(huge(yP1))/10
    zP1=sqrt(huge(zP1))/10
   endif
   if (associated(C%Plane2)) then
     call PlaneNearest(xP2,yP2,zP2,x,y,z,C%Plane2)
   else
    xP2=sqrt(huge(xP2))/10
    yP2=sqrt(huge(yP2))/10
    zP2=sqrt(huge(zP2))/10
   endif
   call JacketNearest(xJ,yJ,zJ,x,y,z,C%Jacket)

   if (PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xP2,yP2,zP2).and.&
        PointDist(x,y,z,xP1,yP1,zP1)<=PointDist(x,y,z,xJ,yJ,zJ)) then
                xnear=xP1
                ynear=yP1
                znear=zP1
   elseif (PointDist(x,y,z,xP2,yP2,zP2)<=PointDist(x,y,z,xP1,yP1,zP1).and.&
            PointDist(x,y,z,xP2,yP2,zP2)<=PointDist(x,y,z,xJ,yJ,zJ)) then
                xnear=xP2
                ynear=yP2
                znear=zP2
   elseif (PointDist(x,y,z,xJ,yJ,zJ)<=PointDist(x,y,z,xP2,yP2,zP2).and.&
            PointDist(x,y,z,xJ,yJ,zJ)<=PointDist(x,y,z,xP1,yP1,zP1)) then
                xnear=xJ
                ynear=yJ
                znear=zJ
   endif
  endsubroutine CylinderNearestOut


  subroutine PolyhedronNearestOut(xnear,ynear,znear,x,y,z,PH)
  real(KND) xnear,ynear,znear,x,y,z
  type(TPolyhedron) PH
  real(KND) dists(PH%nplanes),xP(PH%nplanes),yP(PH%nplanes),zP(PH%nplanes),minv
  integer nearest,i

   dists=huge(minv)
   do i=1,PH%nplanes
    call PlaneNearest(xP(i),yP(i),zP(i),x,y,z,PH%Planes(i))
    dists(i)=SQRT((x-xP(i))**2+(y-yP(i))**2+(z-zp(i))**2)
   enddo

   nearest=0
   minv=huge(minv)
   do i=1,PH%nplanes
    if (dists(i)>=0.and.dists(i)<minv) then
      nearest=i
      minv=dists(i)
    endif
   enddo

   xnear=xP(nearest)
   ynear=yP(nearest)
   znear=zP(nearest)
  endsubroutine PolyhedronNearestOut





  subroutine SolidBodyNearestOut(xnear,ynear,znear,x,y,z,SB)
   real(KND),INTENT(OUT):: xnear,ynear,znear
   real(KND) x,y,z
   TYPE(TSolidBody):: SB

  xnear=-1E+9;ynear=-1E+9;znear=-1E+9

  select case (SB%typeofbody)
    case (1)
     call PolyhedronNearestOut(xnear,ynear,znear,x,y,z,SB%Polyhedron)
    case (2)
     call BallNearest(xnear,ynear,znear,x,y,z,SB%Ball)
    case (3)
     call CylinderNearestOut(xnear,ynear,znear,x,y,z,SB%Cylinder)
    case (4)
     call TerrainNearest(xnear,ynear,znear,x,y,z,SB%Terrain)
    case default
     xnear=-huge(znear);ynear=-huge(znear);znear=-huge(znear)
    end select
  endsubroutine SolidBodyNearestOut






  subroutine GetSolidBodiesBC
  integer i,j,k,m,n,o
  integer(SINT) nb
  real(KND) dist,nearx,neary,nearz
  TYPE(TSolidBody),pointer:: CurrentSB => null()
  TYPE(WMPOINT):: WMP

   !find if the gridpoints lie inside a solid body and write it's number
   !do not nullifie the .type arrays, they could have been made nonzero by other unit
   CurrentSB=>FirstSB

   do
    if (associated(CurrentSB)) then
     do k=0,Prnz+1
      do j=0,Prny+1
       do i=0,Prnx+1
          if (SolidBodyInt(xPr(i),yPr(j),zPr(k),CurrentSB)) Prtype(i,j,k)=CurrentSB%numofbody
       enddo
      enddo
     enddo
    else
      EXIT
    endif
    CurrentSB=>CurrentSB%next
   enddo
      
   CurrentSB=>FirstSB
   do
    if (associated(CurrentSB)) then
    do k=-1,Unz+2
     do j=-1,Uny+2
       do i=-1,Unx+2
          if (SolidBodyInt(xU(i),yPr(j),zPr(k),CurrentSB)) Utype(i,j,k)=CurrentSB%numofbody
       enddo
      enddo
     enddo
    else
      EXIT
    endif
    CurrentSB=>CurrentSB%next
   enddo

   CurrentSB=>FirstSB
   do
    if (associated(CurrentSB)) then
     do k=-1,Vnz+2
      do j=-1,Vny+2
       do i=-1,Vnx+2
          if (SolidBodyInt(xPr(i),yV(j),zPr(k),CurrentSB)) Vtype(i,j,k)=CurrentSB%numofbody
       enddo
      enddo
     enddo
    else
      EXIT
    endif
    CurrentSB=>CurrentSB%next
   enddo

   CurrentSB=>FirstSB
   do
    if (associated(CurrentSB)) then
     do k=-1,Wnz+2
      do j=-1,Wny+2
       do i=-1,Wnx+2
          if (SolidBodyInt(xPr(i),yPr(j),zW(k),CurrentSB)) Wtype(i,j,k)=CurrentSB%numofbody
       enddo
      enddo
     enddo
    else
      EXIT
    endif
    CurrentSB=>CurrentSB%next
   enddo


   call InitImBoundaries
   
   allocate(WMP%depscalar(computescalars))

   do k=1,Prnz 
    do j=1,Prny
     do i=1,Prnx
      if (Prtype(i,j,k)==0) then
       dist=huge(dist)
       nb=0
       do o=-1,1
        do n=-1,1
         do m=-1,1
          if ((m/=0.or.n/=0.or.o/=0).and.(Prtype(i+m,j+n,k+o)>0)) then
            call SetCurrentSB(CurrentSB,Prtype(i+m,j+n,k+o))
            call SolidBodyNearest(nearx,neary,nearz,xPr(i),yPr(j),zPr(k),CurrentSB)
            if (SQRT((nearx-xPr(i))**2+(neary-yPr(j))**2+(nearz-zPr(k))**2)<dist) then
             dist=SQRT((nearx-xPr(i))**2+(neary-yPr(j))**2+(nearz-zPr(k))**2)
             nb=Prtype(i+m,j+n,k+o)            
            endif
          endif
         enddo
        enddo
       enddo

       if (nb>0) then
         call SetCurrentSB(CurrentSB,nb)
         WMP%x=i
         WMP%y=j
         WMP%z=k
         WMP%distx=nearx-xPr(i)
         WMP%disty=neary-yPr(j)
         WMP%distz=nearz-zPr(k)
         WMP%ustar=1
         if (CurrentSB%typeofbody==4) then
          if (CurrentSB%Terrain%PrTerrPoints(i,j)%rough) then
            WMP%z0=CurrentSB%Terrain%PrTerrPoints(i,j)%z0
          else
            WMP%z0=0
          endif
         else
          if (CurrentSB%rough) then
            WMP%z0=CurrentSB%z0
          else
            WMP%z0=0
          endif
         endif
         call AddWMPoint(WMP)
       endif
      endif
     enddo
    enddo
   enddo

 
  endsubroutine GetSolidBodiesBC








  subroutine InitSolidBodies
  real(KND) xnear,ynear,znear
  integer i,j,k
  character(70):: str
  character(8):: fname
  TYPE(TSolidBody),pointer:: CurrentSB => null()
  real(KND) a,b,a2,b2,phi

     
! !     CEDVAL A1-3
! !     allocate(FirstSB%Polyhedron)
! !     allocate(FirstSB%Polyhedron%Planes(3))
! !     
! !     FirstSB%Polyhedron%nplanes=3
! !     FirstSB%Polyhedron%Planes(1)%a=0
! !     FirstSB%Polyhedron%Planes(1)%b=0
! !     FirstSB%Polyhedron%Planes(1)%c=1
! !     FirstSB%Polyhedron%Planes(1)%d=-1
! !     FirstSB%Polyhedron%Planes(1)%GL=.false.
! !     FirstSB%Polyhedron%Planes(2)%a=1
! !     FirstSB%Polyhedron%Planes(2)%b=0
! !     FirstSB%Polyhedron%Planes(2)%c=0
! !     FirstSB%Polyhedron%Planes(2)%d=-0.5
! !     FirstSB%Polyhedron%Planes(2)%GL=.false.
! ! 
! !     FirstSB%Polyhedron%Planes(3)%a=-1
! !     FirstSB%Polyhedron%Planes(3)%b=0
! !     FirstSB%Polyhedron%Planes(3)%c=0
! !     FirstSB%Polyhedron%Planes(3)%d=-0.5
! !     FirstSB%Polyhedron%Planes(3)%GL=.false.




!!BASIC LAYOUT OF SURO EXPERIMENTS
! ! !      allocate(FirstSB)
! ! !   
! ! !      FirstSB%numofbody=1
! ! !      FirstSB%typeofbody=1
! ! ! 
! ! !     allocate(FirstSB%Polyhedron)
! ! !     allocate(FirstSB%Polyhedron%Planes(5))
! ! !     
! ! !     FirstSB%Polyhedron%nplanes=5
! ! !     FirstSB%Polyhedron%Planes(1)%a=0
! ! !     FirstSB%Polyhedron%Planes(1)%b=0
! ! !     FirstSB%Polyhedron%Planes(1)%c=1
! ! !     FirstSB%Polyhedron%Planes(1)%d=-3.2
! ! !     FirstSB%Polyhedron%Planes(1)%GL=.false.
! ! ! 
! ! !     FirstSB%Polyhedron%Planes(2)%a=1
! ! !     FirstSB%Polyhedron%Planes(2)%b=0
! ! !     FirstSB%Polyhedron%Planes(2)%c=0
! ! !     FirstSB%Polyhedron%Planes(2)%d=-14.5
! ! !     FirstSB%Polyhedron%Planes(2)%GL=.false.
! ! ! 
! ! !     FirstSB%Polyhedron%Planes(3)%a=-1
! ! !     FirstSB%Polyhedron%Planes(3)%b=0
! ! !     FirstSB%Polyhedron%Planes(3)%c=0
! ! !     FirstSB%Polyhedron%Planes(3)%d=12
! ! !     FirstSB%Polyhedron%Planes(3)%GL=.false.
! ! ! 
! ! !     FirstSB%Polyhedron%Planes(4)%a=0
! ! !     FirstSB%Polyhedron%Planes(4)%b=1
! ! !     FirstSB%Polyhedron%Planes(4)%c=0
! ! !     FirstSB%Polyhedron%Planes(4)%d=-5.2
! ! !     FirstSB%Polyhedron%Planes(4)%GL=.false.
! ! ! 
! ! !     FirstSB%Polyhedron%Planes(5)%a=0
! ! !     FirstSB%Polyhedron%Planes(5)%b=-1
! ! !     FirstSB%Polyhedron%Planes(5)%c=0
! ! !     FirstSB%Polyhedron%Planes(5)%d=-5.2
! ! !     FirstSB%Polyhedron%Planes(5)%GL=.false.
! ! ! 
! ! !     FirstSB%rough=.true.
! ! !     FirstSB%z0=0.005
! ! ! 
! ! !     allocate(FirstSB%next)
! ! !     CurrentSB=>FirstSB%next
! ! ! 
! ! !     allocate(CurrentSB%Polyhedron)
! ! !     allocate(CurrentSB%Polyhedron%Planes(5))
! ! !     
! ! !     CurrentSB%numofbody=2
! ! !     CurrentSB%typeofbody=1
! ! ! 
! ! !     CurrentSB%Polyhedron%nplanes=5
! ! !     CurrentSB%Polyhedron%Planes(1)%a=0
! ! !     CurrentSB%Polyhedron%Planes(1)%b=0
! ! !     CurrentSB%Polyhedron%Planes(1)%c=1
! ! !     CurrentSB%Polyhedron%Planes(1)%d=-3
! ! !     CurrentSB%Polyhedron%Planes(1)%GL=.false.
! ! ! 
! ! !     CurrentSB%Polyhedron%Planes(2)%a=1
! ! !     CurrentSB%Polyhedron%Planes(2)%b=-1
! ! !     CurrentSB%Polyhedron%Planes(2)%c=0
! ! !     CurrentSB%Polyhedron%Planes(2)%d=-22.1
! ! !     CurrentSB%Polyhedron%Planes(2)%GL=.false.
! ! ! 
! ! !     CurrentSB%Polyhedron%Planes(3)%a=1
! ! !     CurrentSB%Polyhedron%Planes(3)%b=1
! ! !     CurrentSB%Polyhedron%Planes(3)%c=0
! ! !     CurrentSB%Polyhedron%Planes(3)%d=0
! ! !     CurrentSB%Polyhedron%Planes(3)%GL=.false.
! ! ! 
! ! !     CurrentSB%Polyhedron%Planes(4)%a=-1
! ! !     CurrentSB%Polyhedron%Planes(4)%b=-1
! ! !     CurrentSB%Polyhedron%Planes(4)%c=0
! ! !     CurrentSB%Polyhedron%Planes(4)%d=-4.24
! ! !     CurrentSB%Polyhedron%Planes(4)%GL=.false.
! ! ! 
! ! !     CurrentSB%Polyhedron%Planes(5)%a=-1
! ! !     CurrentSB%Polyhedron%Planes(5)%b=1
! ! !     CurrentSB%Polyhedron%Planes(5)%c=0
! ! !     CurrentSB%Polyhedron%Planes(5)%d=+20
! ! !     CurrentSB%Polyhedron%Planes(5)%GL=.false.
! ! ! 
! ! !     CurrentSB%rough=.true.
! ! !     CurrentSB%z0=0.005





!SURO EXPERIMENTS IN FRAME OF REFERENCE FOR WIND DIRECTION FROM PHI
! !     xheading=288
! !     phi=108! 4.5.2010....108 (po staru 252),     22.6.2010.....191
! !     phi=MOD((phi-250._KND)/180._KND*pi,2._KND*pi)
! ! 
! !      allocate(FirstSB)
! !   
! !      FirstSB%numofbody=1
! !      FirstSB%typeofbody=1
! ! 
! !     allocate(FirstSB%Polyhedron)
! !     allocate(FirstSB%Polyhedron%Planes(5))
! !     
! !     FirstSB%Polyhedron%nplanes=5
! !     FirstSB%Polyhedron%Planes(1)%a=0
! !     FirstSB%Polyhedron%Planes(1)%b=0
! !     FirstSB%Polyhedron%Planes(1)%c=1
! !     FirstSB%Polyhedron%Planes(1)%d=-6
! !     FirstSB%Polyhedron%Planes(1)%GL=.false.
! ! 
! !           a=1
! !           b=0
! !           a2=a*cos(phi)-b*sin(phi)
! !           b2=a*sin(phi)+b*cos(phi)
! !     FirstSB%Polyhedron%Planes(2)%a=a2
! !     FirstSB%Polyhedron%Planes(2)%b=b2
! !     FirstSB%Polyhedron%Planes(2)%c=0
! !     FirstSB%Polyhedron%Planes(2)%d=-14.5
! !     FirstSB%Polyhedron%Planes(2)%GL=.false.
! ! 
! !           a=-1
! !           b=0
! !           a2=a*cos(phi)-b*sin(phi)
! !           b2=a*sin(phi)+b*cos(phi)
! !     FirstSB%Polyhedron%Planes(3)%a=a2
! !     FirstSB%Polyhedron%Planes(3)%b=b2
! !     FirstSB%Polyhedron%Planes(3)%c=0
! !     FirstSB%Polyhedron%Planes(3)%d=12
! !     FirstSB%Polyhedron%Planes(3)%GL=.false.
! ! 
! !           a=0
! !           b=1
! !           a2=a*cos(phi)-b*sin(phi)
! !           b2=a*sin(phi)+b*cos(phi)
! !     FirstSB%Polyhedron%Planes(4)%a=a2
! !     FirstSB%Polyhedron%Planes(4)%b=b2
! !     FirstSB%Polyhedron%Planes(4)%c=0
! !     FirstSB%Polyhedron%Planes(4)%d=-5.2
! !     FirstSB%Polyhedron%Planes(4)%GL=.false.
! ! 
! !           a=0
! !           b=-1
! !           a2=a*cos(phi)-b*sin(phi)
! !           b2=a*sin(phi)+b*cos(phi)
! !     FirstSB%Polyhedron%Planes(5)%a=a2
! !     FirstSB%Polyhedron%Planes(5)%b=b2
! !     FirstSB%Polyhedron%Planes(5)%c=0
! !     FirstSB%Polyhedron%Planes(5)%d=-5.2
! !     FirstSB%Polyhedron%Planes(5)%GL=.false.
! ! 
! !     FirstSB%rough=.true.
! !     FirstSB%z0=0.005
! 
!     allocate(FirstSB%next)
!     CurrentSB=>FirstSB%next
! 
!     allocate(CurrentSB%Polyhedron)
!     allocate(CurrentSB%Polyhedron%Planes(5))
!     
!     CurrentSB%numofbody=2
!     CurrentSB%typeofbody=1
! 
!     CurrentSB%Polyhedron%nplanes=5
!     CurrentSB%Polyhedron%Planes(1)%a=0
!     CurrentSB%Polyhedron%Planes(1)%b=0
!     CurrentSB%Polyhedron%Planes(1)%c=1
!     CurrentSB%Polyhedron%Planes(1)%d=-3
!     CurrentSB%Polyhedron%Planes(1)%GL=.false.
! 
!           a=1
!           b=-1
!           a2=a*cos(phi)-b*sin(phi)
!           b2=a*sin(phi)+b*cos(phi)
!     CurrentSB%Polyhedron%Planes(2)%a=a2
!     CurrentSB%Polyhedron%Planes(2)%b=b2
!     CurrentSB%Polyhedron%Planes(2)%c=0
!     CurrentSB%Polyhedron%Planes(2)%d=-22.1
!     CurrentSB%Polyhedron%Planes(2)%GL=.false.
! 
!           a=1
!           b=1
!           a2=a*cos(phi)-b*sin(phi)
!           b2=a*sin(phi)+b*cos(phi)
!     CurrentSB%Polyhedron%Planes(3)%a=a2
!     CurrentSB%Polyhedron%Planes(3)%b=b2
!     CurrentSB%Polyhedron%Planes(3)%c=0
!     CurrentSB%Polyhedron%Planes(3)%d=0
!     CurrentSB%Polyhedron%Planes(3)%GL=.false.
! 
!           a=-1
!           b=-1
!           a2=a*cos(phi)-b*sin(phi)
!           b2=a*sin(phi)+b*cos(phi)
!     CurrentSB%Polyhedron%Planes(4)%a=a2
!     CurrentSB%Polyhedron%Planes(4)%b=b2
!     CurrentSB%Polyhedron%Planes(4)%c=0
!     CurrentSB%Polyhedron%Planes(4)%d=-4.24
!     CurrentSB%Polyhedron%Planes(4)%GL=.false.
! 
!           a=-1
!           b=1
!           a2=a*cos(phi)-b*sin(phi)
!           b2=a*sin(phi)+b*cos(phi)
!     CurrentSB%Polyhedron%Planes(5)%a=a2
!     CurrentSB%Polyhedron%Planes(5)%b=b2
!     CurrentSB%Polyhedron%Planes(5)%c=0
!     CurrentSB%Polyhedron%Planes(5)%d=+20
!     CurrentSB%Polyhedron%Planes(5)%GL=.false.
! 
!     CurrentSB%rough=.true.
!     CurrentSB%z0=0.005












!!Thompson's data wide4 building
!      allocate(FirstSB)
!   
!      FirstSB%numofbody=1
!      FirstSB%typeofbody=1
! 
!     allocate(FirstSB%Polyhedron)
!     allocate(FirstSB%Polyhedron%Planes(5))
!     
!     FirstSB%Polyhedron%nplanes=5
!     FirstSB%Polyhedron%Planes(1)%a=0
!     FirstSB%Polyhedron%Planes(1)%b=0
!     FirstSB%Polyhedron%Planes(1)%c=1
!     FirstSB%Polyhedron%Planes(1)%d=-150
!     FirstSB%Polyhedron%Planes(1)%GL=.false.
! 
!     FirstSB%Polyhedron%Planes(2)%a=1
!     FirstSB%Polyhedron%Planes(2)%b=0
!     FirstSB%Polyhedron%Planes(2)%c=0
!     FirstSB%Polyhedron%Planes(2)%d=-150
!     FirstSB%Polyhedron%Planes(2)%GL=.false.
! 
!     FirstSB%Polyhedron%Planes(3)%a=-1
!     FirstSB%Polyhedron%Planes(3)%b=0
!     FirstSB%Polyhedron%Planes(3)%c=0
!     FirstSB%Polyhedron%Planes(3)%d=0
!     FirstSB%Polyhedron%Planes(3)%GL=.false.
! 
!     FirstSB%Polyhedron%Planes(4)%a=0
!     FirstSB%Polyhedron%Planes(4)%b=1
!     FirstSB%Polyhedron%Planes(4)%c=0
!     FirstSB%Polyhedron%Planes(4)%d=-300
!     FirstSB%Polyhedron%Planes(4)%GL=.false.
! 
!     FirstSB%Polyhedron%Planes(5)%a=0
!     FirstSB%Polyhedron%Planes(5)%b=-1
!     FirstSB%Polyhedron%Planes(5)%c=0
!     FirstSB%Polyhedron%Planes(5)%d=-300
!     FirstSB%Polyhedron%Planes(5)%GL=.false.
!      FirstSB%z0=0.015













!     allocate(FirstSB)
!   
!     FirstSB%numofbody=1
!     FirstSB%typeofbody=3
!     allocate(FirstSB%Cylinder)
!     FirstSB%Cylinder%Jacket%xc=0
!     FirstSB%Cylinder%Jacket%yc=0
!     FirstSB%Cylinder%Jacket%zc=0
!     FirstSB%Cylinder%Jacket%a=0
!     FirstSB%Cylinder%Jacket%b=1
!     FirstSB%Cylinder%Jacket%c=0
!     FirstSB%Cylinder%Jacket%r=.5
!     FirstSB%rough=.false.
!     FirstSB%z0=0

! ! !     allocate(FirstSB)
! ! !     CurrentSB=>FirstSB
! ! !     
! ! !     CurrentSB%numofbody=1
! ! !     CurrentSB%typeofbody=4
! ! !     allocate(CurrentSB%Terrain)
! ! !     allocate(CurrentSB%Terrain%UTerrPoints(0:Prnx+1,1:Prny))
! ! !     allocate(CurrentSB%Terrain%VTerrPoints(1:Prnx,0:Prny+1))
! ! !     allocate(CurrentSB%Terrain%PrTerrPoints(1:Prnx,1:Prny))
! ! !     do j=1,Prny
! ! !      do i=1,Prnx
! ! !       CurrentSB%Terrain%PrTerrPoints(i,j)%elev=1._KND/(1._KND+(sqrt(xPr(i)**2)/0.6_KND))!+yPr(j)**2
! ! !       CurrentSB%Terrain%PrTerrPoints(i,j)%rough=.false.
! ! !      enddo
! ! !     enddo
! ! !     do j=1,Prny
! ! !      do i=0,Prnx+1
! ! !       CurrentSB%Terrain%UTerrPoints(i,j)%elev=1._KND/(1+(sqrt(xU(i)**2)/0.6_KND))!+yPr(j)**2
! ! !       CurrentSB%Terrain%UTerrPoints(i,j)%rough=.false.
! ! !      enddo
! ! !     enddo
! ! !     do j=0,Prny+1
! ! !      do i=1,Prnx
! ! !       CurrentSB%Terrain%VTerrPoints(i,j)%elev=1._KND/(1+(sqrt(xPr(i)**2)/0.6_KND))!+yV(j)**2
! ! !       CurrentSB%Terrain%VTerrPoints(i,j)%rough=.false.
! ! !      enddo
! ! !     enddo
! ! !    


!     allocate(FirstSB)
!   
!     FirstSB%numofbody=1
!     FirstSB%typeofbody=2
!     allocate(FirstSB%Ball)
!     FirstSB%Ball%xc=0
!     FirstSB%Ball%yc=0
!     FirstSB%Ball%zc=0
!     FirstSB%Ball%r=1
!     FirstSB%rough=.true.
!     FirstSB%z0=z0inlet





! !   OPEN(11,file="near.vtk")
! !   write (11,"(A)") "# vtk DataFile Version 2.0"
! !   write (11,"(A)") "diplomka output file"
! !   write (11,"(A)") "ASCII"
! !   write (11,"(A)") "DATASET RECTILINEAR_GRID"
! !   str="DIMENSIONS"
! !   write (str(12:),*) Prnx,Prny,Prnz
! !   write (11,"(A)") str
! !   str="X_COORDINATES"
! !   write (str(15:),*) Prnx,"float"
! !   write (11,"(A)") str
! !   write (11,*) xPr(1:Prnx)
! !   str="Y_COORDINATES"
! !   write (str(15:),*) Prny,"float"
! !   write (11,"(A)") str
! !   write (11,*) yPr(1:Prny)
! !   str="Z_COORDINATES"
! !   write (str(15:),*) Prnz,"float"
! !   write (11,"(A)") str
! !   write (11,*) zPr(1:Prnz)
! !   str="POINT_DATA"
! !   write (str(12:),*) Prnx*Prny*Prnz
! !   write (11,"(A)") str
! ! 
! ! 
! !   write (11,"(A)") "VECTORS near float"
! !   do k=1,Prnz
! !    do j=1,Prny
! !     do i=1,Prnx
! !       write(*,*) "-------------"
! !       write(*,*) xPr(i),yPr(j),zPr(k)
! !       if (.not.SolidBodyInt(xPr(i),yPr(j),zPr(k),FirstSB)) then
! !         call SolidBodyNearest(xnear,ynear,znear,xPr(i),yPr(j),zPr(k),FirstSB)
! !       else
! !        xnear=xPr(i);ynear=yPr(j);znear=zPr(k)
! !       endif
! !       Write (11,*) xnear-xPr(i),ynear-yPr(j),znear-zPr(k)
! !     enddo
! !    enddo
! !   enddo
! !  CLOSE(11)



! !      allocate(FirstSB)
! !   
! !      FirstSB%numofbody=1
! !      FirstSB%typeofbody=1
! ! 
! !     allocate(FirstSB%Polyhedron)
! !     allocate(FirstSB%Polyhedron%Planes(4))
! !     
! !     FirstSB%Polyhedron%nplanes=4
! !     FirstSB%Polyhedron%Planes(1)%a=1
! !     FirstSB%Polyhedron%Planes(1)%b=1
! !     FirstSB%Polyhedron%Planes(1)%c=0
! !     FirstSB%Polyhedron%Planes(1)%d=-1
! !     FirstSB%Polyhedron%Planes(1)%GL=.false.
! ! 
! !     FirstSB%Polyhedron%Planes(2)%a=-1
! !     FirstSB%Polyhedron%Planes(2)%b=1
! !     FirstSB%Polyhedron%Planes(2)%c=0
! !     FirstSB%Polyhedron%Planes(2)%d=-1
! !     FirstSB%Polyhedron%Planes(2)%GL=.false.
! ! 
! !     FirstSB%Polyhedron%Planes(3)%a=1
! !     FirstSB%Polyhedron%Planes(3)%b=-1
! !     FirstSB%Polyhedron%Planes(3)%c=0
! !     FirstSB%Polyhedron%Planes(3)%d=-1
! !     FirstSB%Polyhedron%Planes(3)%GL=.false.
! ! 
! !     FirstSB%Polyhedron%Planes(4)%a=-1
! !     FirstSB%Polyhedron%Planes(4)%b=-1
! !     FirstSB%Polyhedron%Planes(4)%c=0
! !     FirstSB%Polyhedron%Planes(4)%d=-1
! !     FirstSB%Polyhedron%Planes(4)%GL=.false.
! 
!     FirstSB%Polyhedron%Planes(5)%a=0
!     FirstSB%Polyhedron%Planes(5)%b=0
!     FirstSB%Polyhedron%Planes(5)%c=1
!     FirstSB%Polyhedron%Planes(5)%d=-1
!     FirstSB%Polyhedron%Planes(5)%GL=.false.
! 
!     FirstSB%Polyhedron%Planes(6)%a=0
!     FirstSB%Polyhedron%Planes(6)%b=0
!     FirstSB%Polyhedron%Planes(6)%c=-1
!     FirstSB%Polyhedron%Planes(6)%d=-1
!     FirstSB%Polyhedron%Planes(6)%GL=.false.

!     FirstSB%Polyhedron%Planes(7)%a=-1
!     FirstSB%Polyhedron%Planes(7)%b=1
!     FirstSB%Polyhedron%Planes(7)%c=-1
!     FirstSB%Polyhedron%Planes(7)%d=-1
!     FirstSB%Polyhedron%Planes(7)%GL=.false.
! 
!     FirstSB%Polyhedron%Planes(8)%a=-1
!     FirstSB%Polyhedron%Planes(8)%b=-1
!     FirstSB%Polyhedron%Planes(8)%c=-1
!     FirstSB%Polyhedron%Planes(8)%d=-1
!     FirstSB%Polyhedron%Planes(8)%GL=.false.



!      allocate(FirstSB)
!   
!      FirstSB%numofbody=1
!      FirstSB%typeofbody=1
! 
!     allocate(FirstSB%Polyhedron)
!     allocate(FirstSB%Polyhedron%Planes(3))
!     
!     FirstSB%Polyhedron%nplanes=3
!     FirstSB%Polyhedron%Planes(1)%a=0
!     FirstSB%Polyhedron%Planes(1)%b=0
!     FirstSB%Polyhedron%Planes(1)%c=1
!     FirstSB%Polyhedron%Planes(1)%d=-0.125
!     FirstSB%Polyhedron%Planes(1)%GL=.false.
! 
!     FirstSB%Polyhedron%Planes(2)%a=1
!     FirstSB%Polyhedron%Planes(2)%b=0
!     FirstSB%Polyhedron%Planes(2)%c=0
!     FirstSB%Polyhedron%Planes(2)%d=-0.075
!     FirstSB%Polyhedron%Planes(2)%GL=.false.
! 
!     FirstSB%Polyhedron%Planes(3)%a=-1
!     FirstSB%Polyhedron%Planes(3)%b=0
!     FirstSB%Polyhedron%Planes(3)%c=0
!     FirstSB%Polyhedron%Planes(3)%d=-0.075
!     FirstSB%Polyhedron%Planes(3)%GL=.false.
! 
!     FirstSB%rough=.false.
!     FirstSB%z0=0

!     allocate(FirstSB)
!     allocate(FirstSB%Polyhedron)
!     allocate(FirstSB%Polyhedron%Planes(4))
!     
!      FirstSB%numofbody=1
!      FirstSB%typeofbody=1
! 
!     FirstSB%Polyhedron%nplanes=4
!     FirstSB%Polyhedron%Planes(1)%a=-1
!     FirstSB%Polyhedron%Planes(1)%b=0
!     FirstSB%Polyhedron%Planes(1)%c=0
!     FirstSB%Polyhedron%Planes(1)%d=-1
!     FirstSB%Polyhedron%Planes(1)%GL=.false.
!     FirstSB%Polyhedron%Planes(2)%a=1
!     FirstSB%Polyhedron%Planes(2)%b=0
!     FirstSB%Polyhedron%Planes(2)%c=0
!     FirstSB%Polyhedron%Planes(2)%d=-1
!     FirstSB%Polyhedron%Planes(2)%GL=.false.
! 
!     FirstSB%Polyhedron%Planes(3)%a=0
!     FirstSB%Polyhedron%Planes(3)%b=0
!     FirstSB%Polyhedron%Planes(3)%c=-1
!     FirstSB%Polyhedron%Planes(3)%d=-1
!     FirstSB%Polyhedron%Planes(3)%GL=.false.
!     FirstSB%Polyhedron%Planes(4)%a=0
!     FirstSB%Polyhedron%Planes(4)%b=0
!     FirstSB%Polyhedron%Planes(4)%c=1
!     FirstSB%Polyhedron%Planes(4)%d=-1
!     FirstSB%Polyhedron%Planes(4)%GL=.false.

   call InitChannelIT2

  endsubroutine InitSolidBodies





  subroutine InitChannelIT1
  TYPE(TSolidBody),pointer:: CurrentSB => null()
  integer i

    allocate(CurrentSB)
    FirstSB=>CurrentSB

   do i=1,20
 
    CurrentSB%numofbody=i
    CurrentSB%typeofbody=1

    allocate(CurrentSB%Polyhedron)
    allocate(CurrentSB%Polyhedron%Planes(3))
    
    CurrentSB%Polyhedron%nplanes=3
    CurrentSB%Polyhedron%Planes(1)%a=0
    CurrentSB%Polyhedron%Planes(1)%b=0
    CurrentSB%Polyhedron%Planes(1)%c=1
    CurrentSB%Polyhedron%Planes(1)%d=-0.05
    CurrentSB%Polyhedron%Planes(1)%GL=.false.

    CurrentSB%Polyhedron%Planes(2)%a=1
    CurrentSB%Polyhedron%Planes(2)%b=0
    CurrentSB%Polyhedron%Planes(2)%c=0
    CurrentSB%Polyhedron%Planes(2)%d=0.025+(10-i)*0.1
    CurrentSB%Polyhedron%Planes(2)%GL=.false.

    CurrentSB%Polyhedron%Planes(3)%a=1
    CurrentSB%Polyhedron%Planes(3)%b=0
    CurrentSB%Polyhedron%Planes(3)%c=0
    CurrentSB%Polyhedron%Planes(3)%d=0.075+(10-i)*0.1
    CurrentSB%Polyhedron%Planes(3)%GL=.true.

    CurrentSB%rough=.false.
    CurrentSB%z0=0

    if (i<20) then
     allocate(CurrentSB%next)
     CurrentSB=>CurrentSB%next
    endif
   enddo

   endsubroutine InitChannelIT1


  subroutine InitChannelIT2
  TYPE(TSolidBody),pointer:: CurrentSB => null()
  integer i
  real(KND) x, d1, d2
    allocate(CurrentSB)
    FirstSB=>CurrentSB

   do i=1,20
 
    CurrentSB%numofbody=i
    CurrentSB%typeofbody=1

    allocate(CurrentSB%Polyhedron)
    allocate(CurrentSB%Polyhedron%Planes(3))
    
    CurrentSB%Polyhedron%nplanes=3
    CurrentSB%Polyhedron%Planes(1)%a=0
    CurrentSB%Polyhedron%Planes(1)%b=0
    CurrentSB%Polyhedron%Planes(1)%c=1
    CurrentSB%Polyhedron%Planes(1)%d=-0.03
    CurrentSB%Polyhedron%Planes(1)%GL=.false.

    CurrentSB%Polyhedron%Planes(2)%a=1
    CurrentSB%Polyhedron%Planes(2)%b=0
    CurrentSB%Polyhedron%Planes(2)%c=0
    CurrentSB%Polyhedron%Planes(2)%d=0.025+(10-i)*0.1
    CurrentSB%Polyhedron%Planes(2)%GL=.false.

    CurrentSB%Polyhedron%Planes(3)%a=1
    CurrentSB%Polyhedron%Planes(3)%b=0
    CurrentSB%Polyhedron%Planes(3)%c=0
    CurrentSB%Polyhedron%Planes(3)%d=0.075+(10-i)*0.1
    CurrentSB%Polyhedron%Planes(3)%GL=.true.

    CurrentSB%rough=.false.
    CurrentSB%z0=0

     allocate(CurrentSB%next)
     CurrentSB=>CurrentSB%next
   enddo

   do i=21,40
 
    CurrentSB%numofbody=i
    CurrentSB%typeofbody=1

    allocate(CurrentSB%Polyhedron)
    allocate(CurrentSB%Polyhedron%Planes(3))
    
    x=0.05+(30-i)*0.1
    d1=-2*x-0.125
    d2=2*x-0.125

    CurrentSB%Polyhedron%nplanes=3
    CurrentSB%Polyhedron%Planes(1)%a=0
    CurrentSB%Polyhedron%Planes(1)%b=0
    CurrentSB%Polyhedron%Planes(1)%c=1
    CurrentSB%Polyhedron%Planes(1)%d=-0.03
    CurrentSB%Polyhedron%Planes(1)%GL=.true.

    CurrentSB%Polyhedron%Planes(2)%a=2
    CurrentSB%Polyhedron%Planes(2)%b=0
    CurrentSB%Polyhedron%Planes(2)%c=2.5
    CurrentSB%Polyhedron%Planes(2)%d=d1
    CurrentSB%Polyhedron%Planes(2)%GL=.false.

    CurrentSB%Polyhedron%Planes(3)%a=-2
    CurrentSB%Polyhedron%Planes(3)%b=0
    CurrentSB%Polyhedron%Planes(3)%c=2.5
    CurrentSB%Polyhedron%Planes(3)%d=d2
    CurrentSB%Polyhedron%Planes(3)%GL=.false.

    CurrentSB%rough=.false.
    CurrentSB%z0=0

    if (i<40) then
     allocate(CurrentSB%next)
     CurrentSB=>CurrentSB%next
    endif
   enddo

   endsubroutine InitChannelIT2



  subroutine AddIBpoint(IBP)
  type(TIBPoint):: IBP
  
  if (.not.associated(LastIBPoint)) then
   allocate(FirstIBPoint)
   FirstIBPoint=IBP
   LastIBPoint=>FirstIBPoint
  else
   allocate(LastIBPoint%next)
   LastIBPoint%next=IBP
   LastIBPoint=>LastIBPoint%next
  endif
  endsubroutine AddIBPoint

  subroutine AddScalIBpoint(SIBP)
  type(TScalFlIBPoint):: SIBP
  
  if (.not.associated(LastScalFlIBPoint)) then
   allocate(FirstScalFlIBPoint)
   FirstScalFlIBPoint=SIBP
   LastScalFlIBPoint=>FirstScalFlIBPoint
  else
   allocate(LastScalFlIBPoint%next)
   LastScalFlIBPoint%next=SIBP
   LastScalFlIBPoint=>LastScalFlIBPoint%next
  endif
  endsubroutine AddScalIBPoint


  subroutine NearestOnLineOut(x,y,z,x2,y2,z2,t,SB)
  TYPE(TSolidBody):: SB
  real(KND):: x,y,z,x2,y2,z2,t,t1,t2
  integer i

  t1=0
  t2=1
  if (SolidBodyInt(x2,y2,z2,SB)) then
   do
    t2=t2*2
    if (.not.SolidBodyInt(x+(x2-x)*t2,y+(y2-y)*t2,z+(z2-z)*t2,SB)) EXIT
   enddo
  endif

  do i=1,20
   t=(t1+t2)/2._KND
   if (SolidBodyInt(x+(x2-x)*t,y+(y2-y)*t,z+(z2-z)*t,SB)) then
    t1=t
   else
    t2=t
   endif
   if (abs(t1-t2)<MIN(dxmin/1000._KND,dymin/1000._KND,dzmin/1000._KND)) EXIT
  enddo
  t=1
  endsubroutine NearestOnLineOut



  subroutine GetUIBPoint(IBP,xi,yj,zk)
  TYPE(TIBPoint) IBP
  TYPE(TSolidBody),pointer:: SB
  integer xi,yj,zk,dirx,diry,dirz,n1,n2
  real(KND) x,y,z,xnear,ynear,znear,t
  logical free100,free010,free001

  x=xU(xi)
  y=yPr(yj)
  z=zPr(zk)
  call SetCurrentSB(SB,Utype(xi,yj,zk))
  call SolidBodyNearestOut(xnear,ynear,znear,x,y,z,SB)
  IBP%component=1
  IBP%x=xi
  IBP%y=yj
  IBP%z=zk
  IBP%distx=xnear-x
  IBP%disty=ynear-y
  IBP%distz=znear-z
  IBP%dirx=nint(sign(1.0_KND,IBP%distx))
  IBP%diry=nint(sign(1.0_KND,IBP%disty))
  IBP%dirz=nint(sign(1.0_KND,IBP%distz))
  dirx=abs(IBP%dirx)
  diry=abs(IBP%diry)
  dirz=abs(IBP%dirz)


  if (abs(IBP%distx)<(xU(xi+1)-xU(xi-1))/20._KND) then
   IBP%distx=0
   dirx=0
   IBP%dirx=0
  endif
  if (abs(IBP%disty)<(yPr(yj+1)-yPr(yj-1))/20._KND) then
   IBP%disty=0
   diry=0
   IBP%diry=0
  endif
  if (abs(IBP%distz)<(zPr(zk+1)-zPr(zk-1))/20._KND) then
   IBP%distz=0
   dirz=0
   IBP%dirz=0
  endif

  if (dirx==0) then
   free100=.false.
  else
   free100=.not.SolidBodyInt(xU(xi+IBP%dirx),yPr(yj),zPr(zk),SB)
  endif
  if (diry==0) then
   free010=.false.
  else
   free010=.not.SolidBodyInt(xU(xi),yPr(yj+IBP%diry),zPr(zk),SB)
  endif
  if (dirz==0) then
   free001=.false.
  else
   free001=.not.SolidBodyInt(xU(xi),yPr(yj),zPr(zk+IBP%dirz),SB)
  endif


  n1=0
  n2=0
  if (.not.SolidBodyInt(xU(xi+1),yPr(yj),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xU(xi-1),yPr(yj),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xU(xi),yPr(yj+1),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xU(xi),yPr(yj-1),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xU(xi),yPr(yj),zPr(zk+1),SB)) n1=n1+1
  if (.not.SolidBodyInt(xU(xi),yPr(yj),zPr(zk-1),SB)) n1=n1+1
  if (dirx/=0) n2=n2+1
  if (diry/=0) n2=n2+1
  if (dirz/=0) n2=n2+1
  if (n1>n2) then
   IBP%interp=0
   IBP%distx=0
   IBP%disty=0
   IBP%distz=0
   dirx=0
   diry=0
   dirz=0
   IBP%dirx=0
   IBP%diry=0
   IBP%dirz=0
  endif


  
  if ((dirx==0.and.diry==0.and.dirz==0).or.((.not.free001).and.(.not.free010).and.(.not.free100))) then
   IBP%interp=0
   IBP%interpdir=0
   IBP%distx=0
   IBP%disty=0
   IBP%distz=0
   dirx=0
   diry=0
   dirz=0
   IBP%dirx=0
   IBP%diry=0
   IBP%dirz=0
   
  elseif ((free100.and.dirx==1).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then
   IBP%interp=3
   IBP%interpdir=0
  elseif ((.not.free100).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then
   IBP%interp=2
   IBP%interpdir=1
   if (dirx==1) then
    call NearestOnLineOut(x,y,z,x,y+IBP%disty,z+IBP%distz,t,SB)

    dirx=0
    IBP%dirx=0
    IBP%distx=0
    IBP%disty=IBP%disty*t
    IBP%distz=IBP%distz*t
   endif
  elseif ((.not.free010).and.(free100.and.dirx==1).and.(free001.and.dirz==1)) then
   IBP%interp=2
   IBP%interpdir=2
   if (diry==1) then
    call NearestOnLineOut(x,y,z,x+IBP%distx,y,z+IBP%distz,t,SB)

    diry=0
    IBP%diry=0
    IBP%distx=IBP%distx*t
    IBP%disty=0
    IBP%distz=IBP%distz*t
   endif
  elseif ((.not.free001).and.(free100.and.dirx==1).and.(free010.and.diry==1)) then
   IBP%interp=2
   IBP%interpdir=3
   if (dirz==1) then
    call NearestOnLineOut(x,y,z,x+IBP%distx,y+IBP%disty,z,t,SB)

    dirz=0
    IBP%dirz=0
    IBP%distx=IBP%distx*t
    IBP%disty=IBP%disty*t
    IBP%distz=0
   endif
  elseif (free100.and.dirx==1) then
   IBP%interp=1
   IBP%interpdir=1
   if (diry==1.or.dirz==1) then
    call NearestOnLineOut(x,y,z,x+IBP%distx,y,z,t,SB)

    diry=0
    dirz=0
    IBP%diry=0
    IBP%dirz=0
    IBP%distx=t*IBP%distx
    IBP%disty=0
    IBP%distz=0
   endif
  elseif (free010.and.diry==1) then
   IBP%interp=1
   IBP%interpdir=2
   if (dirx==1.or.dirz==1) then
    call NearestOnLineOut(x,y,z,x,y+IBP%disty,z,t,SB)

    dirx=0
    dirz=0
    IBP%dirx=0
    IBP%dirz=0
    IBP%distx=0
    IBP%disty=t*IBP%disty
    IBP%distz=0
   endif
  elseif (free001.and.dirz==1) then
   IBP%interp=1
   IBP%interpdir=3
   if (dirx==1.or.diry==1) then
    call NearestOnLineOut(x,y,z,x,y,z+IBP%distz,t,SB)

    dirx=0
    diry=0
    IBP%dirx=0
    IBP%diry=0
    IBP%distx=0
    IBP%disty=0
    IBP%distz=t*IBP%distz
   endif
  else
   write(*,*) "Assertion error"
   write(*,*) "free100",free100
   write(*,*) "free010",free010
   write(*,*) "free001",free001
   write(*,*) "dirx",dirx
   write(*,*) "diry",diry
   write(*,*) "dirz",dirz
   write(*,*) "i",IBP%x
   write(*,*) "j",IBP%y
   write(*,*) "j",IBP%z
   write(*,*) "distx",IBP%distx
   write(*,*) "disty",IBP%disty
   write(*,*) "distz",IBP%distz
   write(*,*) "x",x
   write(*,*) "y",y
   write(*,*) "z",z
   write(*,*) "component",IBP%component
   stop 
  endif
  endsubroutine GetUIBPoint


  subroutine GetVIBPoint(IBP,xi,yj,zk)
  TYPE(TIBPoint) IBP
  TYPE(TSolidBody),pointer:: SB
  integer xi,yj,zk,dirx,diry,dirz,n1,n2
  real(KND) x,y,z,xnear,ynear,znear,t
  logical free100,free010,free001


  x=xPr(xi)
  y=yV(yj)
  z=zPr(zk)
  call SetCurrentSB(SB,Vtype(xi,yj,zk))
  call SolidBodyNearestOut(xnear,ynear,znear,x,y,z,SB)

  IBP%component=2
  IBP%x=xi
  IBP%y=yj
  IBP%z=zk
  IBP%distx=xnear-x
  IBP%disty=ynear-y
  IBP%distz=znear-z
  IBP%dirx=nint(sign(1.0_KND,IBP%distx))
  IBP%diry=nint(sign(1.0_KND,IBP%disty))
  IBP%dirz=nint(sign(1.0_KND,IBP%distz))
  dirx=abs(IBP%dirx)
  diry=abs(IBP%diry)
  dirz=abs(IBP%dirz)
    
  if (abs(IBP%distx)<(xPr(xi+1)-xPr(xi-1))/20._KND) then
   IBP%distx=0
   dirx=0
   IBP%dirx=0
  endif
  if (abs(IBP%disty)<(yV(yj+1)-yV(yj-1))/20._KND) then
   IBP%disty=0
   diry=0
   IBP%diry=0
  endif
  if (abs(IBP%distz)<(zPr(zk+1)-zPr(zk-1))/20._KND) then
   IBP%distz=0
   dirz=0
   IBP%dirz=0
  endif

  if (dirx==0) then
   free100=.false.
  else
   free100=.not.SolidBodyInt(xPr(xi+IBP%dirx),yV(yj),zPr(zk),SB)
  endif
  if (diry==0) then
   free010=.false.
  else
   free010=.not.SolidBodyInt(xPr(xi),yV(yj+IBP%diry),zPr(zk),SB)
  endif
  if (dirz==0) then
   free001=.false.
  else
   free001=.not.SolidBodyInt(xPr(xi),yV(yj),zPr(zk+IBP%dirz),SB)
  endif


  n1=0
  n2=0
  if (.not.SolidBodyInt(xPr(xi+1),yV(yj),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi-1),yV(yj),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yV(yj+1),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yV(yj-1),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yV(yj),zPr(zk+1),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yV(yj),zPr(zk-1),SB)) n1=n1+1
  if (dirx/=0) n2=n2+1
  if (diry/=0) n2=n2+1
  if (dirz/=0) n2=n2+1
  if (n1>n2) then
   IBP%interp=0
   IBP%distx=0
   IBP%disty=0
   IBP%distz=0
   dirx=0
   diry=0
   dirz=0
   IBP%dirx=0
   IBP%diry=0
   IBP%dirz=0
  endif

  if ((dirx==0.and.diry==0.and.dirz==0).or.((.not.free001).and.(.not.free010).and.(.not.free100))) then
   IBP%interp=0
   IBP%distx=0
   IBP%disty=0
   IBP%distz=0
   dirx=0
   diry=0
   dirz=0
   IBP%dirx=0
   IBP%diry=0
   IBP%dirz=0
  elseif ((free100.and.dirx==1).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then
   IBP%interp=3
  elseif ((.not.free100).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then
   IBP%interp=2
   IBP%interpdir=1
   if (dirx==1) then
    call NearestOnLineOut(x,y,z,x,y+IBP%disty,z+IBP%distz,t,SB)

    dirx=0
    IBP%dirx=0
    IBP%distx=0
    IBP%disty=IBP%disty*t
    IBP%distz=IBP%distz*t
   endif
  elseif ((.not.free010).and.(free100.and.dirx==1).and.(free001.and.dirz==1)) then
   IBP%interp=2
   IBP%interpdir=2
   if (diry==1) then
    call NearestOnLineOut(x,y,z,x+IBP%distx,y,z+IBP%distz,t,SB)

    diry=0
    IBP%diry=0
    IBP%distx=IBP%distx*t
    IBP%disty=0
    IBP%distz=IBP%distz*t
   endif
  elseif ((.not.free001).and.(free100.and.dirx==1).and.(free010.and.diry==1)) then
   IBP%interp=2
   IBP%interpdir=3
   if (dirz==1) then
    call NearestOnLineOut(x,y,z,x+IBP%distx,y+IBP%disty,z,t,SB)

    dirz=0
    IBP%dirz=0
    IBP%distx=IBP%distx*t
    IBP%disty=IBP%disty*t
    IBP%distz=0
   endif
  elseif (free100.and.dirx==1) then
   IBP%interp=1
   IBP%interpdir=1
   if (diry==1.or.dirz==1) then
    call NearestOnLineOut(x,y,z,x+IBP%distx,y,z,t,SB)

    diry=0
    dirz=0
    IBP%diry=0
    IBP%dirz=0
    IBP%distx=t*IBP%distx
    IBP%disty=0
    IBP%distz=0
   endif
  elseif (free010.and.diry==1) then
   IBP%interp=1
   IBP%interpdir=2
   if (dirx==1.or.dirz==1) then
    call NearestOnLineOut(x,y,z,x,y+IBP%disty,z,t,SB)

    dirx=0
    dirz=0
    IBP%dirx=0
    IBP%dirz=0
    IBP%distx=0
    IBP%disty=t*IBP%disty
    IBP%distz=0
   endif
  elseif (free001.and.dirz==1) then
   IBP%interp=1
   IBP%interpdir=3
   if (diry==1.or.dirx==1) then
    call NearestOnLineOut(x,y,z,x,y,z+IBP%distz,t,SB)

    dirx=0
    diry=0
    IBP%dirx=0
    IBP%diry=0
    IBP%distx=0
    IBP%disty=0
    IBP%distz=t*IBP%distz
   endif
  else
   write(*,*) "Assertion error"
   write(*,*) "free100",free100
   write(*,*) "free010",free010
   write(*,*) "free001",free001
   write(*,*) "dirx",dirx
   write(*,*) "diry",diry
   write(*,*) "dirz",dirz
   write(*,*) "i",IBP%x
   write(*,*) "j",IBP%y
   write(*,*) "j",IBP%z
   write(*,*) "x",x
   write(*,*) "y",y
   write(*,*) "z",z
   write(*,*) "IBP%component",IBP%component
   stop
  endif

  endsubroutine GetVIBPoint


  subroutine GetWIBPoint(IBP,xi,yj,zk)
  TYPE(TIBPoint) IBP
  TYPE(TSolidBody),pointer:: SB
  integer xi,yj,zk,dirx,diry,dirz,n1,n2
  real(KND) x,y,z,xnear,ynear,znear,t
  logical free100,free010,free001


  x=xPr(xi)
  y=yPr(yj)
  z=zW(zk)
  call SetCurrentSB(SB,Wtype(xi,yj,zk))
  call SolidBodyNearestOut(xnear,ynear,znear,x,y,z,SB)

  IBP%component=3
  IBP%x=xi
  IBP%y=yj
  IBP%z=zk
  IBP%distx=xnear-x
  IBP%disty=ynear-y
  IBP%distz=znear-z
  IBP%dirx=nint(sign(1.0_KND,IBP%distx))
  IBP%diry=nint(sign(1.0_KND,IBP%disty))
  IBP%dirz=nint(sign(1.0_KND,IBP%distz))
  dirx=abs(IBP%dirx)
  diry=abs(IBP%diry)
  dirz=abs(IBP%dirz)
  
  if (abs(IBP%distx)<(xPr(xi+1)-xPr(xi-1))/20._KND) then
   IBP%distx=0
   dirx=0
   IBP%dirx=0
  endif
  if (abs(IBP%disty)<(yPr(yj+1)-yPr(yj-1))/20._KND) then
   IBP%disty=0
   diry=0
   IBP%diry=0
  endif
  if (abs(IBP%distz)<(zW(zk+1)-zW(zk-1))/20._KND) then
   IBP%distz=0
   dirz=0
   IBP%dirz=0
  endif


  if (dirx==0) then
   free100=.false.
  else
   free100=.not.SolidBodyInt(xPr(xi+IBP%dirx),yPr(yj),zW(zk),SB)
  endif
  if (diry==0) then
   free010=.false.
  else
   free010=.not.SolidBodyInt(xPr(xi),yPr(yj+IBP%diry),zW(zk),SB)
  endif
  if (dirz==0) then
   free001=.false.
  else
   free001=.not.SolidBodyInt(xPr(xi),yPr(yj),zW(zk+IBP%dirz),SB)
  endif


  n1=0
  n2=0
  if (.not.SolidBodyInt(xPr(xi+1),yPr(yj),zW(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi-1),yPr(yj),zW(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yPr(yj+1),zW(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yPr(yj-1),zW(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yPr(yj),zW(zk+1),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yPr(yj),zW(zk-1),SB)) n1=n1+1
  if (dirx/=0) n2=n2+1
  if (diry/=0) n2=n2+1
  if (dirz/=0) n2=n2+1
  if (n1>n2) then
   IBP%interp=0
   IBP%distx=0
   IBP%disty=0
   IBP%distz=0
   dirx=0
   diry=0
   dirz=0
   IBP%dirx=0
   IBP%diry=0
   IBP%dirz=0
  endif


  
  if ((dirx==0.and.diry==0.and.dirz==0).or.((.not.free001).and.(.not.free010).and.(.not.free100))) then
   IBP%interp=0
   IBP%distx=0
   IBP%disty=0
   IBP%distz=0
   dirx=0
   diry=0
   dirz=0
   IBP%dirx=0
   IBP%diry=0
   IBP%dirz=0
  elseif ((free100.and.dirx==1).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then
   IBP%interp=3
  elseif ((.not.free100).and.(free010.and.diry==1).and.(free001.and.dirz==1)) then
   IBP%interp=2
   IBP%interpdir=1
   if (dirx==1) then
    call NearestOnLineOut(x,y,z,x,y+IBP%disty,z+IBP%distz,t,SB)

    dirx=0
    IBP%dirx=0
    IBP%distx=0
    IBP%disty=IBP%disty*t
    IBP%distz=IBP%distz*t
   endif
  elseif ((.not.free010).and.(free100.and.dirx==1).and.(free001.and.dirz==1)) then
   IBP%interp=2
   IBP%interpdir=2
   if (diry==1) then
    call NearestOnLineOut(x,y,z,x+IBP%distx,y,z+IBP%distz,t,SB)

    diry=0
    IBP%diry=0
    IBP%distx=IBP%distx*t
    IBP%disty=0
    IBP%distz=IBP%distz*t
   endif
  elseif ((.not.free001).and.(free100.and.dirx==1).and.(free010.and.diry==1)) then
   IBP%interp=2
   IBP%interpdir=3
   if (dirz==1) then
    call NearestOnLineOut(x,y,z,x+IBP%distx,y+IBP%disty,z,t,SB)

    dirz=0
    IBP%dirz=0
    IBP%distx=IBP%distx*t
    IBP%disty=IBP%disty*t
    IBP%distz=0
   endif
  elseif (free100.and.dirx==1) then
   IBP%interp=1
   IBP%interpdir=1
   if (diry==1.or.dirz==1) then
    call NearestOnLineOut(x,y,z,x+IBP%distx,y,z,t,SB)

    diry=0
    dirz=0
    IBP%diry=0
    IBP%dirz=0
    IBP%distx=t*IBP%distx
    IBP%disty=0
    IBP%distz=0
   endif
  elseif (free010.and.diry==1) then
   IBP%interp=1
   IBP%interpdir=2
   if (dirx==1.or.dirz==1) then
    call NearestOnLineOut(x,y,z,x,y+IBP%disty,z,t,SB)

    dirx=0
    dirz=0
    IBP%dirx=0
    IBP%dirz=0
    IBP%distx=0
    IBP%disty=t*IBP%disty
    IBP%distz=0
   endif
  elseif (free001.and.dirz==1) then
   IBP%interp=1
   IBP%interpdir=3
   if (dirx==1.or.diry==1) then
    call NearestOnLineOut(x,y,z,x,y,z+IBP%distz,t,SB)

    dirx=0
    diry=0
    IBP%dirx=0
    IBP%diry=0
    IBP%distx=0
    IBP%disty=0
    IBP%distz=t*IBP%distz
   endif
  else
   write(*,*) "Assertion error"
   write(*,*) "free100",free100
   write(*,*) "free010",free010
   write(*,*) "free001",free001
   write(*,*) "dirx",dirx
   write(*,*) "diry",diry
   write(*,*) "dirz",dirz
   write(*,*) "i",IBP%x
   write(*,*) "j",IBP%y
   write(*,*) "j",IBP%z
   write(*,*) "x",x
   write(*,*) "y",y
   write(*,*) "z",z
   write(*,*) "IBP%component",IBP%component
   stop
  endif
  endsubroutine GetWIBPoint


  subroutine GeTScalIFlBPoint(IBP,xi,yj,zk)
  TYPE(TScalFlIBPoint) IBP
  TYPE(TSolidBody),pointer:: SB
  integer xi,yj,zk,dirx,diry,dirz,dirx2,diry2,dirz2,n1,n2
  real(KND) x,y,z,xnear,ynear,znear,distx,disty,distz,t,tx,ty,tz
  logical freep00,free0p0,free00p,freem00,free0m0,free00m


  x=xPr(xi)
  y=yPr(yj)
  z=zPr(zk)
  call SetCurrentSB(SB,Prtype(xi,yj,zk))
  call SolidBodyNearestOut(xnear,ynear,znear,x,y,z,SB)

  IBP%x=xi
  IBP%y=yj
  IBP%z=zk
  distx=xnear-x
  disty=ynear-y
  distz=znear-z
  dirx=nint(sign(1.0_KND,distx))
  diry=nint(sign(1.0_KND,disty))
  dirz=nint(sign(1.0_KND,distz))

  if (abs(distx)<(xPr(xi+1)-xPr(xi-1))/20._KND) then
   distx=0
   dirx=0
  endif
  if (abs(disty)<(yPr(yj+1)-yPr(yj-1))/20._KND) then
   disty=0
   diry=0
  endif
  if (abs(distz)<(zW(zk+1)-zW(zk-1))/20._KND) then
   distz=0
   dirz=0
  endif

  freep00=.not.SolidBodyInt(xPr(xi+1),yPr(yj),zPr(zk),SB)
  free0p0=.not.SolidBodyInt(xPr(xi),yPr(yj+1),zPr(zk),SB)
  free00p=.not.SolidBodyInt(xPr(xi),yPr(yj),zPr(zk+1),SB)
  freem00=.not.SolidBodyInt(xPr(xi-1),yPr(yj),zPr(zk),SB)
  free0m0=.not.SolidBodyInt(xPr(xi),yPr(yj-1),zPr(zk),SB)
  free00m=.not.SolidBodyInt(xPr(xi),yPr(yj),zPr(zk-1),SB)

  n1=0
  n2=0
  if (.not.SolidBodyInt(xPr(xi+1),yPr(yj),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi-1),yPr(yj),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yPr(yj+1),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yPr(yj-1),zPr(zk),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yPr(yj),zPr(zk+1),SB)) n1=n1+1
  if (.not.SolidBodyInt(xPr(xi),yPr(yj),zPr(zk-1),SB)) n1=n1+1
  if (dirx/=0) n2=n2+1
  if (diry/=0) n2=n2+1
  if (dirz/=0) n2=n2+1
  if (n1>n2) then
   IBP%interp=0
   distx=0
   disty=0
   distz=0
   dirx=0
   diry=0
   dirz=0
   dirx=0
   diry=0
   dirz=0
  endif


  if (dirx==0.and.diry==0.and.dirz==0) then
    dirx2=0
    diry2=0
    dirz2=0
    if (freep00) dirx2=1
    if (freem00) dirx2=-1
    if (free0p0) diry2=1
    if (free0m0) diry2=-1
    if (free00p) dirz2=1
    if (free00m) dirz2=-1
    IBP%intpointi(1)=IBP%x+dirx2
    IBP%intpointj(1)=IBP%y+diry2
    IBP%intpointk(1)=IBP%z+dirz2
    IBP%intcoef=1._KND
    IBP%interp=1
    IBP%dist=sqrt((x-xPr(IBP%x+dirx2))**2+(y-yPr(IBP%y+diry2))**2+(z-zPr(IBP%z+dirz2))**2)
  elseif (n2==1) then
    IBP%intpointi(1)=IBP%x+dirx
    IBP%intpointj(1)=IBP%y+diry
    IBP%intpointk(1)=IBP%z+dirz
    IBP%intcoef=1._KND
    IBP%interp=1
    IBP%dist=sqrt((x-xPr(IBP%x+dirx))**2+(y-yPr(IBP%y+diry))**2+(z-zPr(IBP%z+dirz))**2)
  elseif (n2==2) then
   if (dirx==0) then
    if (abs(disty/distz)>abs(yPr(yj+diry)-y)/abs(zPr(zk+dirz)-z)) then
     t=(zPr(zk+dirz)-z)/distz
     IBP%intpointi(1)=IBP%x
     IBP%intpointj(1)=IBP%y+diry
     IBP%intpointk(1)=IBP%z+dirz
     IBP%intcoef(1)=abs(disty*t)/abs(yPr(IBP%y+diry)-yPr(IBP%y))
     IBP%intpointi(2)=IBP%x
     IBP%intpointj(2)=IBP%y
     IBP%intpointk(2)=IBP%z+dirz
     IBP%intcoef(2)=1-IBP%intcoef(1)
     IBP%interp=2
     IBP%dist=sqrt((disty*t)**2+(distz*t)**2)
    else
     t=(yPr(yj+diry)-y)/disty
     IBP%intpointi(1)=IBP%x
     IBP%intpointj(1)=IBP%y+diry
     IBP%intpointk(1)=IBP%z+dirz
     IBP%intcoef(1)=abs(distz*t)/abs(zPr(IBP%z+dirz)-zPr(IBP%z))
     IBP%intpointi(2)=IBP%x
     IBP%intpointj(2)=IBP%y+diry
     IBP%intpointk(2)=IBP%z
     IBP%intcoef(2)=1-IBP%intcoef(1)
     IBP%interp=2
     IBP%dist=sqrt((disty*t)**2+(distz*t)**2)
    endif
   elseif (diry==0) then
    if (abs(distx/distz)>abs(xPr(xi+dirx)-x)/abs(zPr(zk+dirz)-z)) then
     t=(zPr(zk+dirz)-z)/distz
     IBP%intpointi(1)=IBP%x+dirx
     IBP%intpointj(1)=IBP%y
     IBP%intpointk(1)=IBP%z+dirz
     IBP%intcoef(1)=abs(distx*t)/abs(xPr(IBP%x+dirx)-xPr(IBP%x))
     IBP%intpointi(2)=IBP%x
     IBP%intpointj(2)=IBP%y
     IBP%intpointk(2)=IBP%z+dirz
     IBP%intcoef(2)=1-IBP%intcoef(1)
     IBP%interp=2
     IBP%dist=sqrt((distx*t)**2+(distz*t)**2)
    else
     t=(xPr(xi+dirx)-x)/distx
     IBP%intpointi(1)=IBP%x+dirx
     IBP%intpointj(1)=IBP%y
     IBP%intpointk(1)=IBP%z+dirz
     IBP%intcoef(1)=abs(distz*t)/abs(zPr(IBP%z+dirz)-zPr(IBP%z))
     IBP%intpointi(2)=IBP%x+dirx
     IBP%intpointj(2)=IBP%y
     IBP%intpointk(2)=IBP%z
     IBP%intcoef(2)=1-IBP%intcoef(1)
     IBP%interp=2
     IBP%dist=sqrt((distx*t)**2+(distz*t)**2)
    endif
   else
    if (abs(distx/disty)>abs(xPr(xi+dirx)-x)/abs(yPr(yj+diry)-y)) then
     t=(yPr(yj+diry)-y)/disty
     IBP%intpointi(1)=IBP%x+dirx
     IBP%intpointj(1)=IBP%y+diry
     IBP%intpointk(1)=IBP%z
     IBP%intcoef(1)=abs(distx*t)/abs(xPr(IBP%x+dirx)-xPr(IBP%x))
     IBP%intpointi(2)=IBP%x
     IBP%intpointj(2)=IBP%y+diry
     IBP%intpointk(2)=IBP%z
     IBP%intcoef(2)=1-IBP%intcoef(1)
     IBP%interp=2
     IBP%dist=sqrt((distx*t)**2+(disty*t)**2)
    else
     t=(xPr(xi+dirx)-x)/distx
     IBP%intpointi(1)=IBP%x+dirx
     IBP%intpointj(1)=IBP%y+diry
     IBP%intpointk(1)=IBP%z
     IBP%intcoef(1)=abs(disty*t)/abs(yPr(IBP%y+diry)-yPr(IBP%y))
     IBP%intpointi(2)=IBP%x+dirx
     IBP%intpointj(2)=IBP%y
     IBP%intpointk(2)=IBP%z
     IBP%intcoef(2)=1-IBP%intcoef(1)
     IBP%interp=2
     IBP%dist=sqrt((distx*t)**2+(disty*t)**2)
    endif
   endif
  else
   tx=(xPr(xi+dirx)-x)/distx
   ty=(yPr(yj+diry)-y)/disty
   tz=(zPr(zk+dirz)-z)/distz

   IBP%interp=4
   if (tx<=ty.and.tx<=tz) then
    !coordinatess if interpolation point are therefore
    !xPr(xi+dirx)=x+tx*distx
    !y+tx*disty
    !z+tz*distz
    IBP%dist=sqrt((distx*tx)**2+(disty*tx)**2+(distz*tx)**2)
    IBP%intpointi(1)=IBP%x+dirx
    IBP%intpointj(1)=IBP%y+diry
    IBP%intpointk(1)=IBP%z+dirz
    IBP%intcoef(1)=abs(disty*tx)/abs(yPr(yj+diry)-y)*&
                              abs(distz*tx)/abs(zPr(zk+dirz)-z)
    IBP%intpointi(2)=IBP%x+dirx
    IBP%intpointj(2)=IBP%y
    IBP%intpointk(2)=IBP%z+dirz
    IBP%intcoef(2)=(1-abs(disty*tx)/abs(yPr(yj+diry)-y))*&
                              abs(distz*tx)/abs(zPr(zk+dirz)-z)
    IBP%intpointi(3)=IBP%x+dirx
    IBP%intpointj(3)=IBP%y+diry
    IBP%intpointk(3)=IBP%z
    IBP%intcoef(3)=abs(disty*tx)/abs(yPr(yj+diry)-y)*&
                              (1-abs(distz*tx)/abs(zPr(zk+dirz)-z))
    IBP%intpointi(4)=IBP%x+dirx
    IBP%intpointj(4)=IBP%y
    IBP%intpointk(4)=IBP%z
    IBP%intcoef(4)=(1-abs(disty*tx)/abs(yPr(yj+diry)-y))*&
                              (1-abs(distz*tx)/abs(zPr(zk+dirz)-z))
   elseif (ty<=tx.and.ty<=tz) then
    !the same with ty
    IBP%dist=sqrt((distx*ty)**2+(disty*ty)**2+(distz*ty)**2)
    IBP%intpointi(1)=IBP%x+dirx
    IBP%intpointj(1)=IBP%y+diry
    IBP%intpointk(1)=IBP%z+dirz
    IBP%intcoef(1)=abs(distx*ty)/abs(xPr(xi+dirx)-x)*&
                              abs(distz*ty)/abs(zPr(zk+dirz)-z)
    IBP%intpointi(2)=IBP%x
    IBP%intpointj(2)=IBP%y+diry
    IBP%intpointk(2)=IBP%z+dirz
    IBP%intcoef(2)=(1-abs(distx*ty)/abs(xPr(xi+dirx)-x))*&
                              abs(distz*ty)/abs(zPr(zk+dirz)-z)
    IBP%intpointi(3)=IBP%x+dirx
    IBP%intpointj(3)=IBP%y+diry
    IBP%intpointk(3)=IBP%z
    IBP%intcoef(3)=abs(distx*ty)/abs(xPr(xi+dirx)-x)*&
                              (1-abs(distz*ty)/abs(zPr(zk+dirz)-z))
    IBP%intpointi(4)=IBP%x
    IBP%intpointj(4)=IBP%y+diry
    IBP%intpointk(4)=IBP%z
    IBP%intcoef(4)=(1-abs(distx*ty)/abs(xPr(xi+dirx)-x))*&
                              (1-abs(distz*ty)/abs(zPr(zk+dirz)-z))
   else
    !the same with tz
    IBP%dist=sqrt((distx*tz)**2+(disty*tz)**2+(distz*tz)**2)
    IBP%intpointi(1)=IBP%x+dirx
    IBP%intpointj(1)=IBP%y+diry
    IBP%intpointk(1)=IBP%z+dirz
    IBP%intcoef(1)=abs(distx*tz)/abs(xPr(xi+dirx)-x)*&
                              abs(disty*tz)/abs(yPr(yj+diry)-y)
    IBP%intpointi(2)=IBP%x
    IBP%intpointj(2)=IBP%y+diry
    IBP%intpointk(2)=IBP%z+dirz
    IBP%intcoef(2)=(1-abs(distx*tz)/abs(xPr(xi+dirx)-x))*&
                              abs(disty*tz)/abs(yPr(yj+diry)-y)
    IBP%intpointi(3)=IBP%x+dirx
    IBP%intpointj(3)=IBP%y
    IBP%intpointk(3)=IBP%z+dirz
    IBP%intcoef(3)=abs(distx*tz)/abs(xPr(xi+dirx)-x)*&
                              (1-abs(disty*tz)/abs(yPr(yj+diry)-y))
    IBP%intpointi(4)=IBP%x
    IBP%intpointj(4)=IBP%y
    IBP%intpointk(4)=IBP%z+dirz
    IBP%intcoef(4)=(1-abs(distx*tz)/abs(xPr(xi+dirx)-x))*&
                              (1-abs(disty*tz)/abs(yPr(yj+diry)-y))
   endif
  endif
  endsubroutine GeTScalIFlBPoint



  subroutine InitImBoundaries
  TYPE(TIBPoint) IBP
  TYPE(TScalFlIBPoint) SIBP
  integer i,j,k


  do k=1,Unz
   do j=1,Uny
    do i=1,Unx
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

  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
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

  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
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

  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
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

  endsubroutine InitImBoundaries



  real(KND) function IBLinInt(h0,h1,h2,vel1,vel2)
  real(KND) h0,h1,h2,vel1,vel2
  if (h0<=h1-h0) then
   IBLinInt=-(h0/(h1-h0))*vel1
  elseif  (h1-h0<h1/10._KND) then
   IBLinInt= 0
  else   
   IBLinInt=-((h2-2*h0)*vel1+(2*h0-h1)*vel2)/(h2-h1)
  endif
  endfunction IBLinInt

  real(KND) function IBBiLinInt(x0,y0,x1,y1,velx,vely,velxy)
  real(KND) x0,y0,x1,y1,velx,vely,velxy,a,b

  a=(x1-x0)/x1
  b=(y1-y0)/y1
   IBBiLinInt=-(a*(1-b)*vely+(1-a)*(1-b)*velxy+(1-a)*b*velx)/(a*b)
  endfunction IBBiLinInt

  real(KND) function IBTriLinInt(x0,y0,z0,x1,y1,z1,velx,vely,velxy,velz,velxz,velyz,velxyz)
  real(KND) x0,y0,z0,x1,y1,z1,velx,vely,velxy,velz,velxz,velyz,velxyz,a,b,c

  a=x0/x1
  b=y0/y1
  c=z0/z1
    IBTriLinInt=- (a*(1-b)*(1-c)*velx+&
                 (1-a)*b*(1-c)*vely+&
                 (1-a)*(1-b)*c*velz+&
                 a*b*(1-c)*velxy+&
                 a*(1-b)*c*velxz+&
                 (1-a)*b*c*velyz+&
                 a*b*c*velxyz)/((1-a)*(1-b)*(1-c))

  endfunction IBTriLinInt



  recursive subroutine DeallIBP(IBP)
  type(TIBPoint),pointer:: IBP

  if (associated(IBP%next)) call DeallIBP(IBP%next)
  deallocate(IBP)
  endsubroutine DeallIBP



  recursive subroutine DeallSB(SB)
  type(TSolidbody),pointer:: SB

  if (associated(SB%next)) call DeallSB(SB%next)
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

  endsubroutine DeallSB



!Nahradit vlastni implementaci!!!!!
SUBROUTINE LEGS (A,N,B,X,INDX)
!
! Subroutine to solve the equation A(N,N)*X(N) = B(N) with the
! partial-pivoting Gaussian elimination scheme.
! Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I,J
  INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
  real(KND), INTENT (INOUT), DIMENSION (N,N) :: A
  real(KND), INTENT (INOUT), DIMENSION (N) :: B
  real(KND), INTENT (OUT), DIMENSION (N) :: X
!
  CALL ELGS (A,N,INDX)
!
  DO I = 1, N-1
    DO J = I+1, N
      B(INDX(J)) = B(INDX(J))-A(INDX(J),I)*B(INDX(I))
    END DO
  END DO
!
  X(N) = B(INDX(N))/A(INDX(N),N)
  DO I = N-1, 1, -1
    X(I) = B(INDX(I))
    DO J = I+1, N
      X(I) = X(I)-A(INDX(I),J)*X(J)
    END DO
    X(I) =  X(I)/A(INDX(I),I)
  END DO
!
END SUBROUTINE LEGS
!
SUBROUTINE ELGS (A,N,INDX)
!
! Subroutine to perform the partial-pivoting Gaussian elimination.
! A(N,N) is the original matrix in the input and transformed matrix
! plus the pivoting element ratios below the diagonal in the output.
! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I,J,K,ITMP
  INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
  real(KND) :: C1,PI,PI1,PJ
  real(KND), INTENT (INOUT), DIMENSION (N,N) :: A
  real(KND), DIMENSION (N) :: C
!
! Initialize the index
!
  DO I = 1, N
    INDX(I) = I
  END DO
!
! Find the rescaling factors, one from each row
!
  DO I = 1, N
    C1= 0.0
    DO J = 1, N
      C1 = MAX(C1,ABS(A(I,J)))
    END DO
    C(I) = C1
  END DO
!
! Search the pivoting (largest) element from each column
!
  DO J = 1, N-1
    PI1 = 0.0
    DO I = J, N
      PI = ABS(A(INDX(I),J))/C(INDX(I))
      IF (PI.GT.PI1) THEN
        PI1 = PI
        K   = I
      ENDIF
    END DO
!
! Interchange the rows via INDX(N) to record pivoting order
!
    ITMP    = INDX(J)
    INDX(J) = INDX(K)
    INDX(K) = ITMP
    DO I = J+1, N
      PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! Record pivoting ratios below the diagonal
!
      A(INDX(I),J) = PJ
!
! Modify other elements accordingly
!
      DO K = J+1, N
        A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      END DO
    END DO
  END DO
!
END SUBROUTINE ELGS

endmodule GEOMETRIC
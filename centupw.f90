!     
! File:   centupw2d.f90
! Author: lada
!
! Created on 18. květen 2006, 14:26
!

MODULE centupw
 use PARAMETERS
 use BOUNDARIES
 implicit none

 interface MINMOD
         module procedure MINMOD2, MINMOD3
  end interface
  real(KND),PARAMETER:: theta=2.0

 contains
  !znaceni: UL(i)~U-(i-1/2),UR(i)~U+(i-1/2) jsou interpolovane rychlosti na stene
  !s indexem jako ma blizky bod Pr(i)  
  subroutine CWENOU(U2,U,V,coef)
  real(KND) U2(-3:,-3:),U(-3:,-3:),V(-3:,-3:)
  real(KND),DIMENSION(LBOUND(U,1):UBOUND(U,1),LBOUND(U,2):UBOUND(U,2)):: UL,UR,VL,VR,VtoU
  real(KND) coef


  call Bound_CondU(U)

  call Bound_CondV(V)
  VtoU=0
  call TRANSVtoU(VtoU,V)

  call INTERPUUx(U,UL,UR)

  call FLUXESUx(U2,UL,UR,coef)

  call INTERPUUy(U,UL,UR)
  
  call INTERPUVy(VtoU,VL,VR)

  call FlUXESUy(U2,UL,UR,VL,VR,coef)
  endsubroutine CWENOU
  
  subroutine CWENOV(V2,U,V,coef)
  real(KND) V2(-3:,-3:),U(-3:,-3:),V(-3:,-3:)
  real(KND),DIMENSION(LBOUND(V,1):UBOUND(V,1),LBOUND(V,2):UBOUND(V,2)):: UL,UR,VL,VR,UtoV
  real(KND) coef


  call Bound_CondU(U)

  call Bound_CondV(V)
  UtoV=0
  call TRANSUtoV(UtoV,U)

  call INTERPVUx(UtoV,UL,UR)

  call INTERPVVx(V,VL,VR)

  call FLUXESVx(V2,UL,UR,VL,VR,coef)

  call INTERPVVy(V,VL,VR)

  call FlUXESVy(V2,VL,VR,coef)
  endsubroutine CWENOV
  
  subroutine FLUXESUX(U2,UL,UR,coef)
  real(KND),DIMENSION(-3:,-3:):: U2,UL,UR
  real(KND),DIMENSION(LBOUND(U2,1):UBOUND(U2,1),LBOUND(U2,2):UBOUND(U2,2)):: H 
  real(KND) coef,p
  integer i,j
  !hrany maji stejny index jako odp. P
  H=0
  
  
   do j=0,Uny
    H(0,j)=-(((UL(0,j)*UL(0,j)+UR(0,j)*UR(0,j))&
        -MAX(ABS(UL(0,j)),ABS(UR(0,j)))*(UR(0,j)-UL(0,j))))/dxU(0)
    do i=1,Unx
      p=((UL(i,j)*UL(i,j)+UR(i,j)*UR(i,j))&
         -MAX(ABS(UL(i,j)),ABS(UR(i,j)))*(UR(i,j)-UL(i,j)))
      H(i-1,j)=H(i-1,j)+p/dxU(i-1)
      H(i,j)=H(i,j)-p/dxU(i)
    enddo
    H(Unx,j)=H(Unx,j)+((UL(Unx+1,j)*UL(Unx+1,j)+UR(Unx+1,j)*UR(Unx+1,j))&
      -MAX(ABS(UL(Unx+1,j)),ABS(UR(Unx+1,j)))*(UR(Unx+1,j)-UL(Unx+1,j)))/dxU(Unx)
   enddo
  
  H=dt*coef*H/2.
  U2=U2-H
  endsubroutine FLUXESUX
  
  subroutine FLUXESUY(U2,UL,UR,VL,VR,coef)
  real(KND),DIMENSION(-3:,-3:):: U2,UL,UR,VL,VR
  real(KND),DIMENSION(LBOUND(U2,1):UBOUND(U2,1),LBOUND(U2,2):UBOUND(U2,2)):: H 
  real(KND) coef,p
  integer i,j
  !hrany maji stejny index jako odp. P
  H=0
  
  
   do i=0,Unx
    H(i,0)=-(((UL(i,0)*VL(i,0)+UR(i,0)*VR(i,0))&
        -MAX(ABS(VL(i,0)),ABS(VR(i,0)))*(UR(i,0)-UL(i,0))))/dyPr(0)
    do j=1,Uny
      p=((UL(i,j)*VL(i,j)+UR(i,j)*VR(i,j))&
         -MAX(ABS(VL(i,j)),ABS(VR(i,j)))*(UR(i,j)-UL(i,j)))
      H(i,j-1)=H(i,j-1)+p/dyPr(j-1)
      H(i,j)=H(i,j)-p/dyPr(j)
    enddo
    H(i,Uny)=H(i,Uny)+((UL(i,Uny+1)*VL(i,Uny+1)+UR(i,Uny+1)*VR(i,Uny+1))&
      -MAX(ABS(VL(i,Uny+1)),ABS(VR(i,Uny+1)))*(UR(i,Uny+1)-UL(i,Uny+1)))/dyPr(Uny)
   enddo
 
  H=dt*coef*H/2.
  U2=U2-H
  endsubroutine FLUXESUY
  
  subroutine FLUXESVX(V2,UL,UR,VL,VR,coef)
  real(KND),DIMENSION(-3:,-3:):: V2,UL,UR,VL,VR
  real(KND),DIMENSION(LBOUND(V2,1):UBOUND(V2,1),LBOUND(V2,2):UBOUND(V2,2)):: H 
  real(KND) coef,p
  integer i,j
  !hrany maji stejny index jako odp. P
  H=0
  
  
   do j=0,Vny
    H(0,j)=-(((UL(0,j)*VL(0,j)+UR(0,j)*VR(0,j))&
        -MAX(ABS(UL(0,j)),ABS(UR(0,j)))*(VR(0,j)-VL(0,j))))/dxPr(0)
    do i=1,Vnx
      p=((UL(i,j)*VL(i,j)+UR(i,j)*VR(i,j))&
         -MAX(ABS(UL(i,j)),ABS(UR(i,j)))*(VR(i,j)-VL(i,j)))
      H(i-1,j)=H(i-1,j)+p/dxPr(i-1)
      H(i,j)=H(i,j)-p/dxPr(i)
    enddo
    H(Vnx,j)=H(Vnx,j)+((UL(Vnx+1,j)*VL(Vnx+1,j)+UR(Vnx+1,j)*VR(Vnx+1,j))&
      -MAX(ABS(UL(Vnx+1,j)),ABS(UR(Vnx+1,j)))*(VR(Vnx+1,j)-VL(Vnx+1,j)))/dxPr(Vnx)
   enddo
  
  H=dt*coef*H/2.
  V2=V2-H
  endsubroutine FLUXESVX
  
  subroutine FLUXESVY(V2,VL,VR,coef)
  real(KND),DIMENSION(-3:,-3:):: V2,VL,VR
  real(KND),DIMENSION(LBOUND(V2,1):UBOUND(V2,1),LBOUND(V2,2):UBOUND(V2,2)):: H 
  real(KND) coef,p
  integer i,j
  !hrany maji stejny index jako odp. P
  H=0
 
   do i=0,Vnx
    H(i,0)=-(((VL(i,0)*VL(i,0)+VR(i,0)*VR(i,0))&
        -MAX(ABS(VL(i,0)),ABS(VR(i,0)))*(VR(i,0)-VL(i,0))))/dyV(0)
    do j=1,Vny
      p=((VL(i,j)*VL(i,j)+VR(i,j)*VR(i,j))&
         -MAX(ABS(VL(i,j)),ABS(VR(i,j)))*(VR(i,j)-VL(i,j)))
      H(i,j-1)=H(i,j-1)+p/dyV(j-1)
      H(i,j)=H(i,j)-p/dyV(j)
    enddo
    H(i,Vny)=H(i,Vny)+((VL(i,Vny+1)*VL(i,Vny+1)+VR(i,Vny+1)*VR(i,Vny+1))&
      -MAX(ABS(VL(i,Vny+1)),ABS(VR(i,Vny+1)))*(VR(i,Vny+1)-VL(i,Vny+1)))/dyV(Vny)
   enddo
  
  H=dt*coef*H/2.
  V2=V2-H
  endsubroutine FLUXESVY
  
  subroutine INTERPUUx(U,UL,UR)
  real(KND),DIMENSION(-3:,-3:):: U,UL,UR
  real(KND) sL,sC,sR
  integer i,j

  do j=0,Uny
   do i=-1,Unx+2
    sL=(U(i,j)-U(i-1,j))/dxPr(i)
    !sC=(U(i+1,j)-U(i-1,j))/(dxPr(i)+dxPr(i+1))
    sR=(U(i+1,j)-U(i,j))/dxPr(i+1)
    sL=sL!*theta
    sR=sR!*theta
    sC=SUPERBEE(sL,sR)!MINMOD(sL,sC,sR)
    UR(i,j)=U(i,j)+sC*(xPr(i)-xU(i))
    UL(i+1,j)=U(i,j)+sC*(xPr(i+1)-xU(i))
   enddo
  enddo
 
  endsubroutine INTERPUUx
  
  subroutine INTERPUUy(U,UL,UR)
  real(KND),DIMENSION(-3:,-3:):: U,UL,UR
  real(KND) sL,sC,sR
  integer i,j

  do i=0,Unx
   do j=-1,Uny+2
    sL=(U(i,j)-U(i,j-1))/dyV(j-1)
    !sC=(U(i,j+1)-U(i,j-1))/(dyV(j-1)+dyV(j))
    sR=(U(i,j+1)-U(i,j))/dyV(j)
    sL=sL!*theta
    sR=sR!*theta
    sC=SUPERBEE(sL,sR)!MINMOD(sL,sC,sR)
    UR(i,j)=U(i,j)+sC*(yV(j-1)-yPr(j))
    UL(i,j+1)=U(i,j)+sC*(yV(j)-yPr(j))
   enddo
  enddo
  endsubroutine INTERPUUy

  subroutine INTERPUVy(V,VL,VR)
  real(KND),DIMENSION(-3:,-3:):: V,VL,VR
  real(KND) sL,sC,sR
  integer i,j
  
  do i=0,Unx
   do j=-1,Uny+2
    sL=(V(i,j)-V(i,j-1))/dyPr(j)
    !sC=(V(i,j+1)-V(i,j-1))/(dyPr(j)+dyPr(j+1))
    sR=(V(i,j+1)-V(i,j))/dyPr(j+1)
    sL=sL!*theta
    sR=sR!*theta
    sC=SUPERBEE(sL,sR)!MINMOD(sL,sC,sR)
    VR(i,j)=V(i,j)+sC*(yV(j-1)-yPr(j))
    VL(i,j+1)=V(i,j)+sC*(yV(j)-yPr(j))
   enddo
  enddo
  endsubroutine INTERPUVy
  
  subroutine INTERPVVx(V,VL,VR)
  real(KND),DIMENSION(-3:,-3:):: V,VL,VR
  real(KND) sL,sC,sR
  integer i,j

  do j=0,Vny
   do i=-1,Vnx+2
    sL=(V(i,j)-V(i-1,j))/dxU(i-1)
    !sC=(V(i+1,j)-V(i-1,j))/(dxU(i-1)+dxU(i))
    sR=(V(i+1,j)-V(i,j))/dxU(i)
    sL=sL!*theta
    sR=sR!*theta
    sC=SUPERBEE(sL,sR)!MINMOD(sL,sC,sR)
    VR(i,j)=V(i,j)+sC*(xU(i-1)-xPr(i))
    VL(i+1,j)=V(i,j)+sC*(xU(i)-xPr(i))
   enddo
  enddo
  endsubroutine INTERPVVx

  subroutine INTERPVUx(U,UL,UR)
  real(KND),DIMENSION(-3:,-3:):: U,UL,UR
  real(KND) sL,sC,sR
  integer i,j

  do j=0,Vny
   do i=-1,Vnx+2
    sL=(U(i,j)-U(i-1,j))/dxU(i-1)
    sC=(U(i+1,j)-U(i-1,j))/(dxU(i-1)+dxU(i))
    sR=(U(i+1,j)-U(i,j))/dxU(i)
    sL=sL!*theta
    sR=sR!*theta
    sC=SUPERBEE(sL,sR)!MINMOD(sL,sC,sR)
    UR(i,j)=U(i,j)+sC*(xU(i-1)-xPr(i))
    UL(i+1,j)=U(i,j)+sC*(xU(i)-xPr(i))
    enddo
  enddo
  endsubroutine INTERPVUx
  
  subroutine INTERPVVy(V,VL,VR)
  real(KND),DIMENSION(-3:,-3:):: V,VL,VR
  real(KND) sL,sC,sR
  integer i,j

  do i=0,Vnx
   do j=-1,Vny+2
    sL=(V(i,j)-V(i,j-1))/dyPr(j)
    sC=(V(i,j+1)-V(i,j-1))/(dyPr(j)+dyPr(j+1))
    sR=(V(i,j+1)-V(i,j))/dyPr(j+1)
    sL=sL!*theta
    sR=sR!*theta
    sC=SUPERBEE(sL,sR)!MINMOD(sL,sC,sR)
    VR(i,j)=V(i,j)+sC*(yPr(j)-yV(j))
    VL(i,j+1)=V(i,j)+sC*(yPr(j+1)-yV(j))
   enddo
  enddo
  endsubroutine INTERPVVy
  
  real(KND) function SUPERBEE(a,b)
  real(KND) a,b,ab,b2
  ab=a*b
  b2=b*b
  if ((a*b)<=0._KND) then
   SUPERBEE=0._KND
  else if (theta*ab<b2) then
   SUPERBEE=a*theta
  else if (ab<b2) then
   SUPERBEE=b
  else if (ab<theta*b2) then
   SUPERBEE=a
  else
   SUPERBEE=b*theta
  endif
  endfunction SUPERBEE  

  
  real(KND) function MINMOD2(a,b)
  real(KND) a,b
      MINMOD2=(SIGN(1.,a)+SIGN(1.,b))*MIN(ABS(a),ABS(b))/2.
  endfunction MINMOD2

  real(KND) function MINMOD3(a,b,c)
  real(KND) a,b,c
    if ((a>0).and.(b>0).and.(c>0)) then  
        MINMOD3=MIN(a,b,c)
    elseif ((a<0).and.(b<0).and.(c<0)) then
        MINMOD3=MAX(a,b,c)
    else
        MINMOD3=0
    endif
  endfunction MINMOD3

  subroutine TRANSVtoU(VtoU,V)    
  real(KND) VtoU(-3:,-3:),V(-3:,-3:)    
  integer i,j    
  do j=-2,Uny+2
   do i=-2,Unx+2
    if (.not. (dxU(i)*dyPr(j)==0)) VtoU(i,j)=&
                ((xU(i)-xPr(i))*(yV(j)-yPr(j))*V(i,j)&
               +(xU(i)-xPr(i))*(yPr(j)-yV(j-1))*V(i,j-1)&
               +(xPr(i+1)-xU(i))*(yV(j)-yPr(j))*V(i+1,j)&
               +(xPr(i+1)-xU(i))*(yPr(j)-yV(j-1))*V(i+1,j-1)&
               )/(dxU(i)*dyPr(j))
   enddo
  enddo
  endsubroutine TRANSVtoU
  
  subroutine TRANSUtoV(UtoV,U)
  real(KND) UtoV(-3:,-3:),U(-3:,-3:)
  integer i,j
  do j=-2,Vny+2
   do i=-2,Vnx+2
    if (.not. (dxPr(i)*dyV(j)==0)) UtoV(i,j)=&
               ((xPr(i)-xU(i-1))*(yV(j)-yPr(j))*U(i-1,j)&
              +(xPr(i)-xU(i-1))*(yPr(j+1)-yV(j))*U(i-1,j+1)&
              +(xU(i)-xPr(i))*(yV(j)-yPr(j))*U(i,j)&
              +(xU(i)-xPr(i))*(yPr(j+1)-yV(j))*U(i,j+1)&
              )/(dxPr(i)*dyV(j))
   enddo
  enddo
  endsubroutine TRANSUtoV

END MODULE centupw

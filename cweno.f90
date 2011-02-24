module CWENO
    
 use PARAMETERS
 use BOUNDARIES
 implicit none

 contains
  !znaceni: UL(i)~U-(i-1/2),UR(i)~U+(i-1/2) jsou interpolovane rychlosti na stene
  !s indexem jako ma blizky bod Pr(i)  
  subroutine CWENOU(U2,U,V,coef)
  real(KND) U2(-3:,-3:),U(-3:,-3:),V(-3:,-3:)
  real(KND),DIMENSION(LBOUND(U,1):UBOUND(U,1),LBOUND(U,2):UBOUND(U,2)):: UL,UR,VL,VR,VtoU
  real(KND) coef

  call Bound_CondU(U)
  call Bound_CondV(V)
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
  
  
   do j=1,Uny
    H(0,j)=-((UL(0,j)*UL(0,j)+UR(0,j)*UR(0,j))&
        -MAX(ABS(UL(0,j)),ABS(UR(0,j)))*(UR(0,j)-UL(0,j)))
    do i=1,Unx
      p=((UL(i,j)*UL(i,j)+UR(i,j)*UR(i,j))&
         -MAX(ABS(UL(i,j)),ABS(UR(i,j)))*(UR(i,j)-UL(i,j)))
      H(i-1,j)=H(i-1,j)+p
      H(i,j)=H(i,j)-p
    enddo
    H(Unx,j)=H(Unx,j)+((UL(Unx+1,j)*UL(Unx+1,j)+UR(Unx+1,j)*UR(Unx+1,j))&
      -MAX(ABS(UL(Unx+1,j)),ABS(UR(Unx+1,j)))*(UR(Unx+1,j)-UL(Unx+1,j)))
   enddo
  
  H=dt*coef*H/(2*dxmin)
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
    H(i,1)=-((UL(i,1)*VL(i,1)+UR(i,1)*VR(i,1))&
        -MAX(ABS(VL(i,1)),ABS(VR(i,1)))*(UR(i,1)-UL(i,1)))
    do j=2,Uny
      p=((UL(i,j)*VL(i,j)+UR(i,j)*VR(i,j))&
         -MAX(ABS(VL(i,j)),ABS(VR(i,j)))*(UR(i,j)-UL(i,j)))
      H(i,j-1)=H(i,j-1)+p
      H(i,j)=H(i,j)-p
    enddo
    H(i,Uny)=H(i,Uny)+((UL(i,Uny+1)*VL(i,Uny+1)+UR(i,Uny+1)*VR(i,Uny+1))&
      -MAX(ABS(VL(i,Uny+1)),ABS(VR(i,Uny+1)))*(UR(i,Uny+1)-UL(i,Uny+1)))
   enddo
 
  H=dt*coef*H/(2.*dymin)
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
    H(1,j)=-((UL(1,j)*VL(1,j)+UR(1,j)*VR(1,j))&
        -MAX(ABS(UL(1,j)),ABS(UR(1,j)))*(VR(1,j)-VL(1,j)))
    do i=2,Vnx
      p=((UL(i,j)*VL(i,j)+UR(i,j)*VR(i,j))&
         -MAX(ABS(UL(i,j)),ABS(UR(i,j)))*(VR(i,j)-VL(i,j)))
      H(i-1,j)=H(i-1,j)+p
      H(i,j)=H(i,j)-p
    enddo
    H(Vnx,j)=H(Vnx,j)+((UL(Vnx+1,j)*VL(Vnx+1,j)+UR(Vnx+1,j)*VR(Vnx+1,j))&
      -MAX(ABS(UL(Vnx+1,j)),ABS(UR(Vnx+1,j)))*(VR(Vnx+1,j)-VL(Vnx+1,j)))
   enddo
  
  H=dt*coef*H/(2.*dxmin)
  V2=V2-H
  endsubroutine FLUXESVX
  
  subroutine FLUXESVY(V2,VL,VR,coef)
  real(KND),DIMENSION(-3:,-3:):: V2,VL,VR
  real(KND),DIMENSION(LBOUND(V2,1):UBOUND(V2,1),LBOUND(V2,2):UBOUND(V2,2)):: H 
  real(KND) coef,p
  integer i,j
  !hrany maji stejny index jako odp. P
  H=0
 
   do i=1,Vnx
    H(i,0)=-((VL(i,0)*VL(i,0)+VR(i,0)*VR(i,0))&
        -MAX(ABS(VL(i,0)),ABS(VR(i,0)))*(VR(i,0)-VL(i,0)))
    do j=1,Vny
      p=((VL(i,j)*VL(i,j)+VR(i,j)*VR(i,j))&
         -MAX(ABS(VL(i,j)),ABS(VR(i,j)))*(VR(i,j)-VL(i,j)))
      H(i,j-1)=H(i,j-1)+p
      H(i,j)=H(i,j)-p
    enddo
    H(i,Vny)=H(i,Vny)+((VL(i,Vny+1)*VL(i,Vny+1)+VR(i,Vny+1)*VR(i,Vny+1))&
      -MAX(ABS(VL(i,Vny+1)),ABS(VR(i,Vny+1)))*(VR(i,Vny+1)-VL(i,Vny+1)))
   enddo
  
  H=dt*coef*H/(2*dymin)
  V2=V2-H
  endsubroutine FLUXESVY
  
  subroutine INTERPUUx(U,UL,UR)
  real(KND),DIMENSION(-3:,-3:):: U,UL,UR
  real(KND),DIMENSION(1:3,LBOUND(UL,1):UBOUND(UL,1),LBOUND(UL,2):UBOUND(UL,2)):: a,w
  real(KND),PARAMETER:: e=1E-6,f=2.
  real(KND) p
  integer i,j,m


   do j=1,Uny
    a(1,-1,j)=0.25/((U(-1,j)-U(-2,j))*(U(-1,j)-U(-2,j))+e)**f
    a(2,-1,j)=0.5/((13/3.)*(U(-2,j)-2*U(-1,j)+U(0,j))**2+&
          0.25*(U(0,j)-U(-2,j)*U(0,j)-U(-2,j))+e)**f
    p=0.25/((U(0,j)-U(1,j))*(U(0,j)-U(-1,j))+e)**f
    a(3,-1,j)=p
    do i=0,Unx+1
     a(1,i,j)=p
     a(2,i,j)=1./((13/3.)*(U(i-1,j)-2*U(i,j)+U(i+1,j))**2+&
          0.25*(U(i+1,j)-U(i-1,j)*U(i+1,j)-U(i-1,j))+e)**f
          
     p=1./(4.*(U(i+1,j)-U(i,j))*(U(i+1,j)-U(i,j))+e)**f
     a(3,i,j)=p
    enddo
   enddo

    
   do j=1,Uny
    do i=-1,Unx+1
     p=a(1,i,j)+a(2,i,j)+a(3,i,j)
     do m=1,3     
         w(m,i,j)=a(m,i,j)/p
     enddo
    enddo
   enddo
  UL=0
  UR=0
  do j=1,Uny
   do i=-1,Unx
    do m=1,3
        UL(i+1,j)=UL(i+1,j)+P(U,(i+0.5)*dxmin,i,j,m,1)*w(m,i,j)
        UR(i+1,j)=UR(i+1,j)+P(U,(i+0.5)*dxmin,i+1,j,m,1)*w(m,i+1,j)
    enddo
   enddo
  enddo
  endsubroutine INTERPUUx
  
  subroutine INTERPUUy(U,UL,UR)
  real(KND),DIMENSION(-3:,-3:):: U,UL,UR
  real(KND),DIMENSION(1:3,LBOUND(UL,1):UBOUND(UL,1),LBOUND(UL,2):UBOUND(UL,2)):: a,w
  real(KND),PARAMETER:: e=1E-6,f=2.
  real(KND) p
  integer i,j,m


   do i=0,Unx
    a(1,i,0)=0.25/((U(i,0)-U(i,-1))*(U(i,0)-U(i,-1))+e)**f
    a(2,i,0)=0.5/((13/3.)*(U(i,-1)-2*U(i,0)+U(i,1))**2+&
          0.25*(U(i,1)-U(i,-1)*U(i,1)-U(i,-1))+e)**f
    p=0.25/((U(i,1)-U(i,0))*(U(i,1)-U(i,0))+e)**f
    a(3,i,0)=p
    do j=1,Uny+1
     a(1,i,j)=p
     a(2,i,j)=1./((13/3.)*(U(i,j-1)-2*U(i,j)+U(i,j+1))**2+&
          0.25*(U(i,j+1)-U(i,j-1)*U(i,j+1)-U(i,j-1))+e)**f
          
     p=1./(4.*(U(i,j+1)-U(i,j))*(U(i,j+1)-U(i,j))+e)**f
     a(3,i,j)=p
    enddo
   enddo

    
   do j=0,Uny+1
    do i=0,Unx
     p=a(1,i,j)+a(2,i,j)+a(3,i,j)
     do m=1,3     
         w(m,i,j)=a(m,i,j)/p
     enddo
    enddo
   enddo
  UL=0
  UR=0
  do j=0,Uny
   do i=0,Unx
    do m=1,3
        UL(i,j+1)=UL(i,j+1)+P(U,(j+0.5)*dymin,i,j,m,2)*w(m,i,j)
        UR(i,j+1)=UR(i,j+1)+P(U,(j+0.5)*dymin,i,j+1,m,2)*w(m,i,j+1)
    enddo
   enddo
  enddo
  endsubroutine INTERPUUy

  subroutine INTERPUVy(V,VL,VR)
  real(KND),DIMENSION(-3:,-3:):: V,VL,VR
  real(KND),DIMENSION(1:3,LBOUND(VL,1):UBOUND(VL,1),LBOUND(VL,2):UBOUND(VL,2)):: a,w
  real(KND),PARAMETER:: e=1E-6,f=2.
  real(KND) p
  integer i,j,m
  
   do i=0,Unx
    a(1,i,0)=0.25/((V(i,0)-V(i,-1))*(V(i,0)-V(i,-1))+e)**f
    a(2,i,0)=0.5/((13/3.)*(V(i,-1)-2*V(i,0)+V(i,1))**2+&
          0.25*(V(i,1)-V(i,-1)*V(i,1)-V(i,-1))+e)**f
    p=0.25/((V(i,1)-V(i,0))*(V(i,1)-V(i,0))+e)**f
    a(3,i,0)=p
    do j=1,Vny+1
     a(1,i,j)=p
     a(2,i,j)=1./((13/3.)*(V(i,j-1)-2*V(i,j)+V(i,j+1))**2+&
          0.25*(V(i,j+1)-V(i,j-1)*V(i,j+1)-V(i,j-1))+e)**f
          
     p=1./(4.*(V(i,j+1)-V(i,j))*(V(i,j+1)-V(i,j))+e)**f
     a(3,i,j)=p
    enddo
   enddo

    
   do j=0,Uny+1
    do i=0,Unx
     p=a(1,i,j)+a(2,i,j)+a(3,i,j)
     do m=1,3     
         w(m,i,j)=a(m,i,j)/p
     enddo
    enddo
   enddo
  VL=0
  VR=0
  do j=0,Uny
   do i=0,Unx
    do m=1,3
        VL(i,j+1)=VL(i,j+1)+P(V,(j+0.5)*dymin,i,j,m,2)*w(m,i,j)
        VR(i,j+1)=VR(i,j+1)+P(V,(j+0.5)*dymin,i,j+1,m,2)*w(m,i,j+1)
    enddo
   enddo
  enddo
  endsubroutine INTERPUVy
  
  subroutine INTERPVVx(V,VL,VR)
  real(KND),DIMENSION(-3:,-3:):: V,VL,VR
  real(KND),DIMENSION(1:3,LBOUND(VL,1):UBOUND(VL,1),LBOUND(VL,2):UBOUND(VL,2)):: a,w
  real(KND),PARAMETER:: e=1E-6,f=2.
  real(KND) p
  integer i,j,m

   do j=0,Vny
    a(1,0,j)=0.25/((V(0,j)-V(-1,j))*(V(0,j)-V(-1,j))+e)**f
    a(2,0,j)=0.5/((13/3.)*(V(-1,j)-2*V(0,j)+V(1,j))**2+&
          0.25*(V(1,j)-V(-1,j)*V(1,j)-V(-1,j))+e)**f
    p=0.25/((V(1,j)-V(0,j))*(V(1,j)-V(0,j))+e)**f
    a(3,0,j)=p
    do i=1,Vnx+1
     a(1,i,j)=p
     a(2,i,j)=1./((13/3.)*(V(i-1,j)-2*V(i,j)+V(i+1,j))**2+&
          0.25*(V(i+1,j)-V(i-1,j)*V(i+1,j)-V(i-1,j))+e)**f
          
     p=1./(4.*(V(i+1,j)-V(i,j))*(V(i+1,j)-V(i,j))+e)**f
     a(3,i,j)=p
    enddo
   enddo

    
   do j=0,Vny
    do i=0,Vnx+1
     p=a(1,i,j)+a(2,i,j)+a(3,i,j)
     do m=1,3     
         w(m,i,j)=a(m,i,j)/p
     enddo
    enddo
   enddo
  VL=0
  VR=0
  do j=0,Vny
   do i=0,Vnx
    do m=1,3
        VL(i+1,j)=VL(i+1,j)+P(V,(i+0.5)*dxmin,i,j,m,1)*w(m,i,j)
        VR(i+1,j)=VR(i+1,j)+P(V,(i+0.5)*dxmin,i+1,j,m,1)*w(m,i+1,j)
    enddo
   enddo
  enddo
  endsubroutine INTERPVVx

  subroutine INTERPVUx(U,UL,UR)
  real(KND),DIMENSION(-3:,-3:):: U,UL,UR
  real(KND),DIMENSION(1:3,LBOUND(UL,1):UBOUND(UL,1),LBOUND(UL,2):UBOUND(UL,2)):: a,w
  real(KND),PARAMETER:: e=1E-6,f=2.
  real(KND) p
  integer i,j,m

   do j=0,Vny
    a(1,0,j)=0.25/((U(0,j)-U(-1,j))*(U(0,j)-U(-1,j))+e)**f
    a(2,0,j)=0.5/((13/3.)*(U(-1,j)-2*U(0,j)+U(1,j))**2+&
          0.25*(U(1,j)-U(-1,j)*U(1,j)-U(-1,j))+e)**f
    p=0.25/((U(1,j)-U(0,j))*(U(1,j)-U(0,j))+e)**f
    a(3,0,j)=p
    do i=1,Vnx+1
     a(1,i,j)=p
     a(2,i,j)=1./((13/3.)*(U(i-1,j)-2*U(i,j)+U(i+1,j))**2+&
          0.25*(U(i+1,j)-U(i-1,j)*U(i+1,j)-U(i-1,j))+e)**f
          
     p=1./(4.*(U(i+1,j)-U(i,j))*(U(i+1,j)-U(i,j))+e)**f
     a(3,i,j)=p
    enddo
   enddo

    
   do j=1,Vny
    do i=0,Vnx+1
     p=a(1,i,j)+a(2,i,j)+a(3,i,j)
     do m=1,3     
         w(m,i,j)=a(m,i,j)/p
     enddo
    enddo
   enddo
  UL=0
  UR=0
  do j=1,Vny
   do i=0,Vnx
    do m=1,3
        UL(i+1,j)=UL(i+1,j)+P(U,(i+0.5)*dxmin,i,j,m,1)*w(m,i,j)
        UR(i+1,j)=UR(i+1,j)+P(U,(i+0.5)*dxmin,i+1,j,m,1)*w(m,i+1,j)
    enddo
   enddo
  enddo
  endsubroutine INTERPVUx
  
  subroutine INTERPVVy(V,VL,VR)
  real(KND),DIMENSION(-3:,-3:):: V,VL,VR
  real(KND),DIMENSION(1:3,LBOUND(VL,1):UBOUND(VL,1),LBOUND(VL,2):UBOUND(VL,2)):: a,w
  real(KND),PARAMETER:: e=1E-6,f=2.
  real(KND) p
  integer i,j,m

   do i=1,Vnx
    a(1,i,-1)=0.25/((V(i,-1)-V(i,-2))*(V(i,-1)-V(i,-2))+e)**f
    a(2,i,-1)=0.5/((13/3.)*(V(i,-2)-2*V(i,-1)+V(i,0))**2+&
          0.25*(V(i,0)-V(i,-2)*V(i,0)-V(i,-3))+e)**f
    p=0.25/((V(i,0)-V(i,-1))*(V(i,0)-V(i,-1))+e)**f
    a(3,i,-1)=p
    do j=0,Vny+1
     a(1,i,j)=p
     a(2,i,j)=1./((13/3.)*(V(i,j-1)-2*V(i,j)+V(i,j+1))**2+&
          0.25*(V(i,j+1)-V(i,j-1)*V(i,j+1)-V(i,j-1))+e)**f
          
     p=1./(4.*(V(i,j+1)-V(i,j))*(V(i,j+1)-V(i,j))+e)**f
     a(3,i,j)=p
    enddo
   enddo

    
   do j=0,Vny+1
    do i=0,Vnx
     p=a(1,i,j)+a(2,i,j)+a(3,i,j)
     do m=1,3     
         w(m,i,j)=a(m,i,j)/p
     enddo
    enddo
   enddo
  VL=0
  VR=0
  do j=0,Vny
   do i=0,Vnx
    do m=1,3
        VL(i,j+1)=VL(i,j+1)+P(V,(j+0.5)*dymin,i,j,m,2)*w(m,i,j)
        VR(i,j+1)=VR(i,j+1)+P(V,(j+0.5)*dymin,i,j+1,m,2)*w(m,i,j+1)
    enddo
   enddo
  enddo
  endsubroutine INTERPVVy
  


  real(KND) function P(U,x,i,j,m,dir)
  real(KND) U(-3:,-3:)
  real(KND) x
  integer i,j,m,dir
  
  if (m==1) then
   if (dir==1) then
    P=U(i,j)+(U(i,j)-U(i-1,j))*(x-dxmin*i)/dxmin 
   else
    P=U(i,j)+(U(i,j)-U(i,j-1))*(x-dymin*j)/dymin
   endif
  elseif (m==3) then
   if (dir==1) then
    P=U(i,j)+(U(i+1,j)-U(i,j))*(x-dxmin*i)/dxmin 
   else
    P=U(i,j)+(U(i,j+1)-U(i,j))*(x-dymin*j)/dymin
   endif
  else
   if (dir==1) then
    P=U(i,j)-(1./12)*(U(i+1,j)-2*U(i,j)+U(i-1,j))-(1./12)*(U(i,j-1)-2*U(i,j)+U(i,j-1))&
      +(U(i+1,j)-U(i-1,j))*(x-dxmin*i)/(2*dxmin)+(U(i+1,j)-2*U(i,j)+U(i-1,j))*(x-dxmin*i)*(x-dxmin*i)/(dxmin*dxmin) 
   else
    P=U(i,j)-(1./12)*(U(i+1,j)-2*U(i,j)+U(i-1,j))-(1./12)*(U(i,j-1)-2*U(i,j)+U(i,j-1))&
      +(U(i,j+1)-U(i,j-1))*(x-dymin*j)/(2*dymin)+(U(i,j+1)-2*U(i,j)+U(i,j-1))*(x-dymin*j)*(x-dymin*j)/(dymin*dymin)
   endif
  endif
  endfunction P



  subroutine TRANSVtoU(VtoU,V)    
  real(KND) VtoU(-3:,-3:),V(-3:,-3:)    
  integer i,j    
  do j=-1,Uny+2
   do i=-2,Unx+2
    VtoU(i,j)=((xU(i)-xPr(i))*(yV(j)-yPr(j))*V(i,j)&
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
   do i=-1,Vnx+2
    UtoV(i,j)=((xPr(i)-xU(i-1))*(yV(j)-yPr(j))*U(i-1,j)&
              +(xPr(i)-xU(i-1))*(yPr(j+1)-yV(j))*U(i-1,j+1)&
              +(xU(i)-xPr(i))*(yV(j)-yPr(j))*U(i,j)&
              +(xU(i)-xPr(i))*(yPr(j+1)-yV(j))*U(i,j+1)&
              )/(dxPr(i)*dyV(j))
   enddo
  enddo
  endsubroutine TRANSUtoV
endmodule CWENO2d
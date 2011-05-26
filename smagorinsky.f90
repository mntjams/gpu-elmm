module SMAGORINSKY

 use PARAMETERS
 use BOUNDARIES, only : Bound_CondU, Bound_CondV, Bound_CondW

 implicit none

 private
 public :: Smag, Smag2, StabSmag, DynSmag, Vreman, Filter
 
 real(KND),parameter:: CSmag=0.1_KND
 
 contains
 
  subroutine Smag(U,V,W)  !Standard Smagorinsky model with implicit filtering
   integer i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(inout):: U,V,W

   call Bound_CondU(U)
   call Bound_CondV(V)
   call Bound_CondW(W)

   do k=-1,Prnz+2
    do j=-1,Prny+2
     do i=-1,Prnx+2
       Visc(i,j,k)=NuSmag(i,j,k,U,V,W)
     enddo
    enddo
   enddo
   if (Re>0) then
     Visc=Visc+1._KND/Re
   endif
  endsubroutine Smag

  
  subroutine Smag2(U,V,W,width)  !Standard Smagorinsky model with explicit filtering
   integer,intent(in):: width
   real(KND),dimension(-2:,-2:,-2:):: U,V,W
   integer i,j,k

   call Bound_CondU(U)
   call Bound_CondV(V)
   call Bound_CondW(W)

   do k=-1,Prnz+2
    do j=-1,Prny+2
     do i=-1,Prnx+2
      if (Re>0) then
       Visc(i,j,k)=1._KND/Re+NuSmag2(i,j,k,U,V,W,width*2+1)
      else
       Visc(i,j,k)=NuSmag2(i,j,k,U,V,W,width*2+1)
      endif
     enddo
    enddo
   enddo
  endsubroutine Smag2




  real(KND) function NuSmag(i,j,k,U,V,W)    !subgrid viscosity for Smagorinsky with implicit filtering
   integer,intent(in):: i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND) S(1:3,1:3)
   real(KND) width,Sbar

   width=(dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND)
   call StrainIJ(i,j,k,U,V,W,S)
   Sbar=Strainu(S)
   NuSmag=Sbar*(width*CSmag)**2
  endfunction NuSmag


  real(KND) function NuSmag2(i,j,k,U,V,W,width2)    !subgrid viscosity for Smagorinsky with implicit filtering
   integer,intent(in):: i,j,k,width2
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND) S(1:3,1:3)
   real(KND) width,Sbar

   width=(dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND)
   call StrainIJ(i,j,k,U,V,W,S)
   Sbar=Strainu(S)
   NuSmag2=Sbar*((width2)*width*CSmag)**2
  endfunction NuSmag2
 



  subroutine StabSmag(U,V,W)  !Smagorinsky with a stability correction Brown et al. (1994)
   integer i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(inout):: U,V,W
   real(KND) Ri,l,l0
   real(KND) width,Sbar
   real(KND),parameter:: CS=0.17_KND
   real(KND) S(1:3,1:3)

   call Bound_CondU(U)
   call Bound_CondV(V)
   call Bound_CondW(W)

   do k=-1,Prnz+2
    do j=-1,Prny+2
     do i=-1,Prnx+2
       width=(dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND)

       call StrainIJ(i,j,k,U,V,W,S)
       Sbar=Strainu(S)

       Ri=Rig(i,j,k,U,V,temperature)


       l0=CS*width
       l=WallDamp(l0,z0B,zPr(k))

       Visc(i,j,k)=Sbar*Fm(Ri)*l
       TDiff(i,j,k)=Sbar*Fh(Ri)*l
     enddo
    enddo
   enddo
   if (Re>0) then
     Visc=Visc+1._KND/Re
   endif
  endsubroutine StabSmag


  pure real(KND) function Fm(Ri)       !Adjustment function for stability
   real(KND),intent(in):: Ri           !Pointwise Richarson number
   real(KND),parameter :: Ric=0.25_KND

   if (Ri>0) then
    Fm  =  (1._KND-Ri/Ric)**4
   else
    Fm  =  sqrt(1._KND - 16._KND*Ri)
   endif
  end function Fm

  pure real(KND) function Fh(Ri)       !Adjustment function for stability
   real(KND),intent(in):: Ri           !Pointwise Richarson number
   real(KND),parameter :: Ric=0.25_KND

   if (Ri>0) then
    Fh  =  (1./0.7_KND)  *  (1._KND - Ri/Ric)**4  *  (1._KND - 1.2_KND*Ri)
   else
    Fh  =  sqrt(1._KND - 40._KND*Ri) / 0.7_KND
   endif
  end function Fh

  pure real(KND) function WallDamp(l0,z0,z) !Wall damping for Smagorinsky model
   real(KND),intent(in):: l0,z0,z

   WallDamp= 1._KND / (&
                 1._KND/l0**2  +  1._KND/(0.4_KND*(z+z0))**2 &
                 )
  end function WallDamp
  

  pure real(KND) function Rig(i,j,k,U,V,temperature)
   integer,intent(in):: i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V
   real(KND),dimension(-1:,-1:,-1:),intent(in):: temperature
   real(KND) num,denom

   num=(grav_acc/temperature_ref)*(temperature(i,j,k+1)-temperature(i,j,k-1))/(zPr(k+1)-zPr(k-1))
   denom=((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._KND*(zPr(k+1)-zPr(k-1))))**2
   denom=denom+((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(zPr(k+1)-zPr(k-1))))**2
   if (abs(denom)>tiny(1._KND)*100) then
    Rig=num/denom
   else
    Rig=0
   endif
  endfunction Rig






  pure real(KND) function Strainu(S)      !Magnitude of the strain rate tensor.
   real(KND),intent(in):: S(1:3,1:3)
   integer::ii,jj

   Strainu=0
   do jj=1,3
    do ii=1,3
     Strainu=Strainu+S(ii,jj)*S(ii,jj)
    enddo
   enddo
   Strainu=SQRT(2._KND*Strainu)
  endfunction Strainu

  
  
  pure subroutine StrainIJ(i,j,k,U,V,W,S)     !Computes components of the strain rate tensor.
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND),intent(out):: S(1:3,1:3)
   integer,intent(in):: i,j,k
   real(KND) D(1:3,1:3)
   integer::ii,jj

   D=0

   if (i>-2.and.i<Prnx+3) D(1,1)=(U(i,j,k)-U(i-1,j,k))/dxPr(i)
   if (j>-2.and.j<Prny+3) D(2,2)=(V(i,j,k)-V(i,j-1,k))/dyPr(j)
   if (k>-2.and.j<Prnz+3) D(3,3)=(W(i,j,k)-W(i,j,k-1))/dzPr(k)
   if (j>-1.and.i>-1.and.j<Uny+3.and.i<Unx+4)&
      D(1,2)=(U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(2._KND*(yPr(j+1)-yPr(j-1)))
   if (k>-1.and.i>-1.and.k<Uny+3.and.i<Unx+4)&
      D(1,3)=(U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
   if (j>-1.and.i>-1.and.i<Vnx+3.and.j<Vny+4)&
      D(2,1)=(V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(2._KND*(xPr(i+1)-xPr(i-1)))
   if (j>-1.and.k>-1.and.k<Vnz+3.and.j<Vny+4)&
      D(2,3)=(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
   if (k>-1.and.i>-1.and.i<Wnx+3.and.k<Wnz+4)&
      D(3,1)=(W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(2._KND*(xPr(i+1)-xPr(i-1)))
   if (j>-1.and.k>-1.and.j<Wny+3.and.k<Wnz+4)&
      D(3,2)=(W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(2._KND*(yPr(j+1)-yPr(j-1)))

   do jj=1,3
    do ii=1,3
     S(ii,jj)=(D(ii,jj)+D(jj,ii))
    enddo
   enddo
   S=S/2._KND
  endsubroutine StrainIJ



  

  subroutine DynSmag(U,V,W)    ! Dynamic Smagorinsky
   integer ii,jj,i,j,k
   real(KND),dimension(-2:,-2:,-2:):: U,V,W
   real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3):: FiltU
   real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3):: FiltV
   real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3):: FiltW
   real(KND),dimension(-2:Prnx+3,-2:Prny+3,-2:Prnz+3):: P,CSDyn
   real(KND),dimension(1:3,1:3,-2:Prnx+3,-2:Prny+3,-2:Prnz+3):: M,L
   real(KND),dimension(1:3,1:3):: S
   real(KND) Strain,A,B,sums
   character(70):: str

   call Bound_CondU(U)
   call Bound_CondV(V)
   call Bound_CondW(W)

   FiltU=0
   FiltV=0
   FiltW=0

   call Filter(U,FILTU,Unx,Uny,Unz,1)
   call Filter(V,FILTV,Vnx,Vny,Vnz,1)
   call Filter(W,FILTW,Wnx,Wny,Wnz,1)




   do k=-2,Prnz+3           !Aij -> Lij
    do j=-2,Prny+3
     do i=-2,Prnx+3
      call StrainIJ(i,j,k,U,V,W,S)
      Strain=Strainu(S)
      L(:,:,i,j,k)=Abs(Strain)*S*((dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND))**2
     enddo
    enddo
   enddo
   M=0
   call Filter(L(1,1,:,:,:),M(1,1,:,:,:),Prnx,Prny,Prnz,1)         !F(Aij) -> Mij
   call Filter(L(2,2,:,:,:),M(2,2,:,:,:),Prnx,Prny,Prnz,1)
   call Filter(L(3,3,:,:,:),M(3,3,:,:,:),Prnx,Prny,Prnz,1)
   call Filter(L(1,2,:,:,:),M(1,2,:,:,:),Prnx,Prny,Prnz,1)
   call Filter(L(1,3,:,:,:),M(1,3,:,:,:),Prnx,Prny,Prnz,1)
   call Filter(L(2,3,:,:,:),M(2,3,:,:,:),Prnx,Prny,Prnz,1)
   M(2,1,:,:,:)=M(1,2,:,:,:)
   M(3,1,:,:,:)=M(1,3,:,:,:)
   M(3,2,:,:,:)=M(2,3,:,:,:)
   P=0
   do k=-1,Prnz+2 !FiltS used here for UiUj               F(UiUj) -> L
    do j=-1,Prny+2
     do i=-1,Prnx+2
      P(i,j,k)=((U(i,j,k)+U(i-1,j,k))/2._KND)*((U(i,j,k)+U(i-1,j,k))/2._KND)
     enddo
    enddo
   enddo

   call Filter(P,L(1,1,:,:,:),Prnx,Prny,Prnz,1)

   do k=-1,Prnz+2 !FiltS used here for UiUj
    do j=-1,Prny+2
     do i=-1,Prnx+2
      P(i,j,k)=((V(i,j,k)+V(i,j-1,k))/2._KND)*((V(i,j,k)+V(i,j-1,k))/2._KND)
     enddo
    enddo
   enddo

   call Filter(P,L(2,2,:,:,:),Prnx,Prny,Prnz,1)

   do k=-1,Prnz+2 !FiltS used here for UiUj
    do j=-1,Prny+2
     do i=-1,Prnx+2
      P(i,j,k)=((W(i,j,k)+W(i,j,k-1))/2._KND)*((W(i,j,k)+W(i,j,k-1))/2._KND)
    enddo
    enddo
   enddo

   call Filter(P,L(3,3,:,:,:),Prnx,Prny,Prnz,1)

   do k=-1,Prnz+2 !FiltS used here for UiUj
    do j=-1,Prny+2
     do i=-1,Prnx+2
      P(i,j,k)=((U(i,j,k)+U(i-1,j,k))/2._KND)*((V(i,j,k)+V(i,j-1,k))/2._KND)
     enddo
    enddo
   enddo

   call Filter(P,L(1,2,:,:,:),Prnx,Prny,Prnz,1)

   do k=-1,Prnz+2 !FiltS used here for UiUj
    do j=-1,Prny+2
     do i=-1,Prnx+2
      P(i,j,k)=((U(i,j,k)+U(i-1,j,k))/2._KND)*((W(i,j,k)+W(i,j,k-1))/2._KND)
     enddo
    enddo
   enddo

   call Filter(P,L(1,3,:,:,:),Prnx,Prny,Prnz,1)

   do k=-1,Prnz+2 !FiltS used here for UiUj
    do j=-1,Prny+2
     do i=-1,Prnx+2
      P(i,j,k)=((V(i,j,k)+V(i,j-1,k))/2._KND)*((W(i,j,k)+W(i,j,k-1))/2._KND)
     enddo
    enddo
   enddo

   call Filter(P,L(2,3,:,:,:),Prnx,Prny,Prnz,1)

   do k=-1,Prnz+2 !Modified Leonard stresses         F(Ui)F(Uj)-F(UiUj) -> Lij
    do j=-1,Prny+2
     do i=-1,Prnx+2
      L(1,1,i,j,k)=-L(1,1,i,j,k)+&
                (FiltU(i,j,k)+FiltU(i-1,j,k)/2._KND)*(FiltU(i,j,k)+FiltU(i-1,j,k)/2._KND)
      L(2,2,i,j,k)=-L(2,2,i,j,k)+&
                (FiltV(i,j,k)+FiltV(i,j-1,k)/2._KND)*(FiltV(i,j,k)+FiltV(i,j-1,k)/2._KND)
      L(3,3,i,j,k)=-L(3,3,i,j,k)+&
                (FiltW(i,j,k)+FiltW(i,j,k-1)/2._KND)*(FiltW(i,j,k)+FiltW(i,j,k-1)/2._KND)
      L(1,2,i,j,k)=-L(1,2,i,j,k)+&
                (FiltU(i,j,k)+FiltU(i-1,j,k)/2._KND)*(FiltV(i,j,k)+FiltV(i,j-1,k)/2._KND)
      L(1,3,i,j,k)=-L(1,3,i,j,k)+&
                (FiltU(i,j,k)+FiltU(i-1,j,k)/2._KND)*(FiltW(i,j,k)+FiltW(i,j,k-1)/2._KND)
      L(2,3,i,j,k)=-L(2,3,i,j,k)+&
                (FiltV(i,j,k)+FiltV(i,j-1,k)/2._KND)*(FiltW(i,j,k)+FiltW(i,j,k-1)/2._KND)
     enddo
    enddo
   enddo

   L(2,1,:,:,:)=L(1,2,:,:,:)
   L(3,1,:,:,:)=L(1,3,:,:,:)
   L(3,2,:,:,:)=L(2,3,:,:,:)

   !write(*,*) "i",SUM(L)/(9*Prnx*Prny*Prnz)
   do k=-1,Prnz+2                                ! Bij)-F(Aij) ->Mij
    do j=-1,Prny+2
     do i=-1,Prnx+2
      call StrainIJ(i,j,k,FiltU,FiltV,FiltW,S)
      Strain=Strainu(S)
      M(:,:,i,j,k)=((1._KND/debugparam)**2)*(((dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND))**2)*ABS(Strain)*S(:,:)-M(:,:,i,j,k)
     enddo
    enddo
   enddo

   do k=-1,Prnz+2
    do j=-1,Prny+2
     do i=-1,Prnx+2
      A=0
      B=0
      do jj=1,3
       do ii=1,3
        A=A+L(ii,jj,i,j,k)*M(ii,jj,i,j,k)
        B=B+M(ii,jj,i,j,k)*M(ii,jj,i,j,k)
       enddo
      enddo
      CSDyn(i,j,k)=A
      P(i,j,k)=B
     enddo
    enddo
   enddo



    L(1,1,:,:,:)=CSDyn
    M(1,1,:,:,:)=P

 !     do i=1,10
 !      call  Filter(CSDyn,L(1,1,:,:,:),Prnx,Prny,Prnz,1)
 !      call  Filter(P,M(1,1,:,:,:),Prnx,Prny,Prnz,1)
 !      CSDyn=L(1,1,:,:,:)
 !      P=M(1,1,:,:,:)
 !     enddo

 !   L(1,1,:,:,:)=SUM(L(1:Prnx,1:Prny,1:Prnz,1,1))/(Prnx*Prny*Prnz)
 !   M(1,1,:,:,:)=SUM(M(1:Prnx,1:Prny,1:Prnz,1,1))/(Prnx*Prny*Prnz)





   Sums=0
   do k=-1,Prnz+2
    do j=-1,Prny+2
     do i=-1,Prnx+2
      if (abs(M(1,1,i,j,k))<1.e-4_KND*abs(L(1,1,i,j,k))) then
                              CSDyn(i,j,k)=0
                              if (i>0.and.i<Prnx+1.and.j>0.and.j<Prny+1.and.k>0.and.k<Prnz+1.) Sums=Sums+1
                             else
                              CSDyn(i,j,k)=0.5_KND*L(1,1,i,j,k)/M(1,1,i,j,k)
                              !if (i>0.and.i<Prnx+1.and.j>0.and.j<Prny+1.and.k>0.and.k<Prnz+1.) write (*,*) A,"/",B,"=",CSDyn(i,j,k)
                             endif
     enddo
    enddo
   enddo



   write(*,*)  100*Sums/(Prnx*Prny*Prnz),"% of zero nominator"

   !do j=-1,Prny+2
   ! CSDyn(1:Prnx,j,1:Prnz)=SUM(CSDyn(1:Prnx,j,1:Prnz))/(Prnx*Prnz)
   ! write(*,*) j,CSDyn(1,j,1)
   !enddo

   Strain=0
   do k=-1,Prnz+2
    do j=-1,Prny+2
     do i=-1,Prnx+2
      CSDyn(i,j,k)=MIN(CSDyn(i,j,k),10000._KND)
      CSDyn(i,j,k)=MAX(-0.001_KND,CSDyn(i,j,k))
      if (i>0.and.i<Prnx+1.and.j>0.and.j<Prny+1.and.k>0.and.k<Prnz+1.)  Strain=Strain+CSDyn(i,j,k)
     enddo
    enddo
   enddo
   write(*,*) "Average dynamic constant:",Strain/(Prnx*Prny*Prnz)
   write(*,*) "Maximal dynamic constant:",MAXVAL(CSDyn(1:Prnx,1:Prny,1:Prnz))



   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
      call StrainIJ(i,j,k,U,V,W,S)
      Strain=Strainu(S)
      if (Re>0) then
       Visc(i,j,k)=1._KND/Re+CSDyn(i,j,k)*(((dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND))**2)*Strain
      else
       Visc(i,j,k)=CSDyn(i,j,k)*(((dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND))**2)*Strain
      endif
     enddo
    enddo
   enddo
  endsubroutine DynSmag





  subroutine Vreman(U,V,W)   !Vreman subgrid model (Physics of Fluids, 2004)
   real(KND),dimension(-2:,-2:,-2:),intent(inout):: U,V,W
   integer i,j,k,ii,jj
   real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)::bb
   real(KND),dimension(1:3,1:3,-1:Prnx+2,-1:Prny+2,-1:Prnz+2)::a,b
   real(KND),parameter::c=0.05

   call Bound_CondU(U)
   call Bound_CondV(V)
   call Bound_CondW(W)

   do k=-1,Prnz+2
    do j=-1,Prny+2
     do i=-1,Prnx+2
      a(1,1,i,j,k)=(U(i,j,k)-U(i-1,j,k))/dxPr(i)
      a(2,1,i,j,k)=(U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(2._KND*(dyV(j)+dyV(j-1)))
      a(3,1,i,j,k)=(U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._KND*(dzW(k)+dzW(k-1)))
      a(2,2,i,j,k)=(V(i,j,k)-V(i,j-1,k))/dyPr(j)
      a(1,2,i,j,k)=(V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(2._KND*(dxU(i)+dxU(i-1)))
      a(3,2,i,j,k)=(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(dzW(k)+dzW(k-1)))
      a(3,3,i,j,k)=(W(i,j,k)-W(i,j,k-1))/dzPr(k)
      a(1,3,i,j,k)=(W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(2._KND*(dxU(i)+dxU(i-1)))
      a(2,3,i,j,k)=(W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(2._KND*(dyV(j)+dyV(j-1)))
     enddo
    enddo
   enddo


   forall(k=-1:Prnz+2,j=-1:Prny+2,i=-1:Prnx+2,jj=1:3,ii=1:3)
    b(ii,jj,i,j,k)=(dxPr(i))**2*a(1,ii,i,j,k)*a(1,jj,i,j,k)+&
                   (dyPr(j))**2*a(2,ii,i,j,k)*a(2,jj,i,j,k)+&
                   (dzPr(k))**2*a(3,ii,i,j,k)*a(3,jj,i,j,k)
   endforall

   bb(:,:,:)=          b(1,1,:,:,:)*b(2,2,:,:,:)-b(1,2,:,:,:)**2
   bb(:,:,:)=bb(:,:,:)+b(1,1,:,:,:)*b(3,3,:,:,:)-b(1,3,:,:,:)**2
   bb(:,:,:)=bb(:,:,:)+b(2,2,:,:,:)*b(3,3,:,:,:)-b(2,3,:,:,:)**2

   Visc=0._KND
   forall(k=-1:Prnz+2,j=-1:Prny+2,i=-1:Prnx+2)
        a(1,1,i,j,k)=sum(a(:,:,i,j,k)**2)
   endforall
   forall(k=-1:Prnz+2,j=-1:Prny+2,i=-1:Prnx+2,abs(a(1,1,i,j,k))>1e-5.and.bb(i,j,k)>0)
        Visc(i,j,k)=c*sqrt(bb(i,j,k)/a(1,1,i,j,k))
   endforall

   if (Re>0) then
     Visc=Visc+1._KND/Re
   endif
  endsubroutine Vreman













 subroutine TrapesField(U1,U2,nx,ny,nz)   !Trapesoidal filter
  real(KND),dimension(-2:,-2:,-2:),intent(in):: U1
  real(KND),dimension(-2:,-2:,-2:),intent(out):: U2
  integer,intent(in):: nx,ny,nz
  integer i,j,k

  U2=0
  do k=-1,nz+2
   do j=-1,ny+2
    do i=-1,nx+2
     U2(i,j,k)=(U1(i-1,j,k)+U1(i+1,j,k)+U1(i,j-1,k)+U1(i,j-1,k)+U1(i,j,k-1)+U1(i,j,k+1)+6*U1(i,j,k))/12._KND
    enddo
   enddo
  enddo

 endsubroutine TrapesField



 subroutine TrapesField2(U1,U2,nx,ny,nz) !Filter with weighting
  real(KND),dimension(-2:,-2:,-2:),intent(in):: U1
  real(KND),dimension(-2:,-2:,-2:),intent(out):: U2
  integer,intent(in):: nx,ny,nz
  real(KND),dimension(-1:1):: F
  real(KND) S,S1
  integer ii,jj,i,j,k,kk

  F=(/ 1,2,1 /)
  U2=0
  S1=0
  do kk=-1,1
   do jj=-1,1
    do ii=-1,1
     S1=S1+F(ii)*F(jj)*F(kk)
    enddo
   enddo
  enddo
  do k=-1,nz+2
   do j=-1,ny+2
    do i=-1,nx+2
     S=0
     do kk=-1,1
      do jj=-1,1
       do ii=-1,1
        S=S+F(ii)*F(jj)*F(kk)*U1(i+ii,j+jj,k+kk)
       enddo
      enddo
     enddo
     U2(i,j,k)=S/S1
    enddo
   enddo
  enddo
 endsubroutine TrapesField2



 subroutine Filter(U1,U2,nx,ny,nz,width)    !Calls a selected filter
  real(KND),dimension(-2:,-2:,-2:),intent(in):: U1
  real(KND),dimension(-2:,-2:,-2:),intent(out):: U2
  integer,intent(in):: nx,ny,nz
  integer width

  call TrapesField(U1,U2,nx,ny,nz)

 endsubroutine Filter
 
endmodule SMAGORINSKY

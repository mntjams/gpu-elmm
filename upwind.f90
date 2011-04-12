!     
! File:   upwind2d.f90
! Author: lada
!
! Created on 11. rijen 2006, 21:56
!


module UPWIND
 use PARAMETERS
 use BOUNDARIES
 implicit none

 private minmod,minmod2,minmod3,superbee,limiter,heav,gamma,vanalbada,vanleer

 interface MINMOD
         module procedure MINMOD2, MINMOD3
  end interface
  real(KND),dimension(:,:,:),allocatable:: Uco,Vco,Wco,Ux,Uy,Uz,Vx,Vy,Vz,Wx,Wy,Wz,Ut,Vt,Wt
 contains

  

 subroutine Tau(U2,V2,W2,U,V,W,Pr)
  real(KND),dimension(-2:,-2:,-2:):: U2,V2,W2,U,V,W
  real(KND),dimension(1:,1:,1:):: Pr

  if (step==1) then
  allocate(Ux(-1:Unx+2,-1:Uny+2,-1:Unz+2))
  allocate(Uy(-1:Unx+2,-1:Uny+2,-1:Unz+2))
  allocate(Uz(-1:Unx+2,-1:Uny+2,-1:Unz+2))
  allocate(Vx(-1:Vnx+2,-1:Vny+2,-1:Vnz+2))
  allocate(Vy(-1:Vnx+2,-1:Vny+2,-1:Vnz+2))
  allocate(Vz(-1:Vnx+2,-1:Vny+2,-1:Vnz+2))
  allocate(Wx(-1:Wnx+2,-1:Wny+2,-1:Wnz+2))
  allocate(Wy(-1:Wnx+2,-1:Wny+2,-1:Wnz+2))
  allocate(Wz(-1:Wnx+2,-1:Wny+2,-1:Wnz+2))
  endif
  Ux=10000000
  Uy=10000000
  Uz=10000000
  Vx=10000000
  Vy=10000000
  Vz=10000000
  Wx=10000000
  Wy=10000000
  Wz=10000000
  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)
  call Bound_Pr(Pr)
  if (gridtype==uniformgrid) then
   call SLOPESUNIF(U,V,W)
  else
   call SLOPES(U,V,W)
  endif

  !Ux=0000000
  !Uy=0000000
  !Uz=0000000
  !Vx=0000000
  !Vy=0000000
  !Vz=0000000
  !Wx=0000000
  !Wy=0000000
  !Wz=0000000


  if (step==1) allocate(Ut(-2:Unx+3,-2:Uny+3,-2:Unz+3),Vt(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),Wt(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
  Ut=99999999
  Vt=99999999
  Wt=99999999
  call Utim(U,V,W,Pr)
  call Vtim(U,V,W,Pr)
  call Wtim(U,V,W,Pr)
  Ut=Ut*0.99
  Vt=Vt*0.99
  Wt=Wt*0.99

  if (step==1) allocate(Uco(-2:Unx+3,-2:Uny+3,-2:Unz+3),Vco(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),Wco(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
  Uco=0!-1000000
  Vco=0!1000000
  Wco=0!~1000000
  call Corners(U,V,W)
  !Uco=0
  !Vco=0
  !Wco=0
  !deallocate(Uy,Uz,Vx,Vz,Wx,Wy)

  call TauU(U2,U,V,W)
  !deallocate(Ux,Ut)
  call TauV(V2,U,V,W)
  !deallocate(Vy,Vt)
  call TauW(W2,U,V,W)
  !deallocate(Wz,Wt,Uco,Vco,Wco)
 endsubroutine Tau

  

 subroutine SLOPES(U,V,W)
 real(KND) U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
 integer i,j,k
 real(KND) Ul,Ur,Uc,Vl,Vr,Vc,Wl,Wr,Wc

 !$OMP PARALLEL PRIVATE(Ul,Ur,Uc,Vl,Vr,Vc,Wl,Wr,Wc)
 !$OMP DO
  do k=-1,Unz+2
   do j=-1,Uny+2
    do i=-1,Unx+2
     Uc=(U(i+1,j,k)-U(i-1,j,k))/(dxPr(i+1)+dxPr(i))
     Ur=(U(i+1,j,k)-U(i,j,k))/(dxPr(i+1))
     Ul=(U(i,j,k)-U(i-1,j,k))/(dxPr(i))
     Ux(i,j,k)=LIMITER(Ul,Ur,Uc)
    enddo
   enddo
  enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
  do k=-1,Vnz+2
   do j=-1,Vny+2
    do i=-1,Vnx+2
     Vc=(V(i,j+1,k)-V(i,j-1,k))/(dyPr(j+1)+dyPr(j))
     Vr=(V(i,j+1,k)-V(i,j,k))/(dyPr(j+1))
     Vl=(V(i,j,k)-V(i,j-1,k))/(dyPr(j))
     Vy(i,j,k)=LIMITER(Vl,Vr,Vc)
    enddo
   enddo
  enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
  do k=-1,Wnz+2
   do j=-1,Wny+2
    do i=-1,Wnx+2
     Wc=(W(i,j,k+1)-W(i,j,k-1))/(dzPr(k+1)+dzPr(k))
     Wr=(W(i,j,k+1)-W(i,j,k))/(dzPr(k+1))
     Wl=(W(i,j,k)-W(i,j,k-1))/(dzPr(k))
     Wz(i,j,k)=LIMITER(Wl,Wr,Wc)
    enddo
   enddo
  enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Unz+2
  do j=-1,Uny+2
   do i=-1,Unx+2
    Uc=(U(i,j+1,k)-U(i,j-1,k))/(dyV(j)+dyV(j-1))
    Ur=(U(i,j+1,k)-U(i,j,k))/(dyV(j))
    Ul=(U(i,j,k)-U(i,j-1,k))/(dyV(j-1))
    Uy(i,j,k)=LIMITER(Ul,Ur,Uc)
   enddo
  enddo
 enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Unz+2
  do j=-1,Uny+2
   do i=-1,Unx+2
    Uc=(U(i,j,k+1)-U(i,j,k-1))/(dzW(k)+dzW(k-1))
    Ur=(U(i,j,k+1)-U(i,j,k))/(dzW(k))
    Ul=(U(i,j,k)-U(i,j,k-1))/(dzW(k-1))
    Uz(i,j,k)=LIMITER(Ul,Ur,Uc)
   enddo
  enddo
 enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Vnz+2
  do j=-1,Vny+2
   do i=-1,Vnx+2
    Vc=(V(i+1,j,k)-V(i-1,j,k))/(dxU(i)+dxU(i-1))
    Vr=(V(i+1,j,k)-V(i,j,k))/(dxU(i))
    Vl=(V(i,j,k)-V(i-1,j,k))/(dxU(i-1))
    Vx(i,j,k)=LIMITER(Vl,Vr,Vc)
   enddo
  enddo
 enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Vnz+2
  do j=-1,Vny+2
   do i=-1,Vnx+2
    Vc=(V(i,j,k+1)-V(i,j,k-1))/(dzW(k)+dzW(k-1))
    Vr=(V(i,j,k+1)-V(i,j,k))/(dzW(k))
    Vl=(V(i,j,k)-V(i,j,k-1))/(dzW(k-1))
    Vz(i,j,k)=LIMITER(Vl,Vr,Vc)
   enddo
  enddo
 enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Wnz+2
  do j=-1,Wny+2
   do i=-1,Wnx+2
    Wc=(W(i+1,j,k)-W(i-1,j,k))/(dxU(i)+dxU(i-1))
    Wr=(W(i+1,j,k)-W(i,j,k))/(dxU(i))
    Wl=(W(i,j,k)-W(i-1,j,k))/(dxU(i-1))
    Wx(i,j,k)=LIMITER(Wl,Wr,Wc)
   enddo
  enddo
 enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Wnz+2
  do j=-1,Wny+2
   do i=-1,Wnx+2
    Wc=(W(i,j+1,k)-W(i,j-1,k))/(dyV(j)+dyV(j-1))
    Wr=(W(i,j+1,k)-W(i,j,k))/(dyV(j))
    Wl=(W(i,j,k)-W(i,j-1,k))/(dyV(j-1))
    Wy(i,j,k)=LIMITER(Wl,Wr,Wc)
   enddo
  enddo
 enddo
 !$OMP ENDDO
 !$OMP END PARALLEL
 end subroutine SLOPES

  
 subroutine Utim(U,V,W,Pr)
  real(KND) U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND) Vloc,Wloc,Uyc,Uzc,dx2,dy2,dz2
  integer i,j,k,x,y,z,nx,ny,nz
  type(TIBPoint),pointer:: IBP
     !!!!!!!!PORESIT IBM!!!!!!!!!!
  nx=Unx
  ny=Vny+1
  nz=Wnz+1

  if (gridtype==uniformgrid) then
   dx2=dxmin*dxmin
   dy2=dymin*dymin
   dz2=dzmin*dzmin
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Ut(i,j,k)=((U(i+1,j,k)-2*U(i,j,k)+U(i-1,j,k))/(dx2)+&     
       (U(i,j+1,k)-2*U(i,j,k)+U(i,j-1,k))/(dy2)+&
       (U(i,j,k+1)-2*U(i,j,k)+U(i,j,k-1))/(dz2))/Re
     enddo
    enddo
   enddo
  else
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Ut(i,j,k)=(Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))/dxPr(i+1)-&
       Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k))/dxPr(i))/dxU(i)+&
       (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))/dyV(j)-&
       0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k))/dyV(j-1))/dyPr(j)+&
       (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)-&
       0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1))/dzW(k-1))/dzPr(k)
     enddo
    enddo
   enddo
  endif

  do k=1,nz
   do j=1,ny
    do i=1,nx
     Ut(i,j,k)=Ut(i,j,k)-(Pr(i+1,j,k)-Pr(i,j,k))/dxU(i)
    enddo
   enddo
  enddo
  do k=1,nz
   do j=1,ny
    do i=0,nx+1
     Ut(i,j,k)=Ut(i,j,k)-Ux(i,j,k)*U(i,j,k)
    enddo
   enddo
  enddo

  if (gridtype==uniformgrid) then
   do k=1,nz
    do j=1,ny
     do i=1,nx 
      Vloc=(V(i,j-1,k)+V(i,j,k)+V(i+1,j-1,k)+V(i+1,j,k))/4._KND
      if (Vloc>=0) then
                   Uyc=(U(i,j,k)-U(i,j-1,k))/dymin+(0.5_KND)*(1-Vloc*dt/dymin)&
                        *(Uy(i,j,k)-Uy(i,j-1,k))
      else
                   Uyc=(U(i,j+1,k)-U(i,j,k))/dyV(j)-(0.5_KND)*(1+Vloc*dt/dymin)&
                        *(Uy(i,j+1,k)-Uy(i,j,k))
      endif
      Ut(i,j,k)=Ut(i,j,k)-Uyc*Vloc
     enddo
    enddo
   enddo
  else
   do k=1,nz
    do j=1,ny
     do i=1,nx 
      Vloc=(V(i,j-1,k)+V(i,j,k)+V(i+1,j-1,k)+V(i+1,j,k))/4._KND
      !Vloc=(V(i,j-1,k)+V(i,j,k)+V(i+1,j-1,k)+V(i+1,j,k)+V(i,j-1,k)+V(i,j,k-1)+V(i+1,j-1,k-1)+V(i+1,j,k-1))/8._KND
      if (Vloc>=0) then
                   Uyc=(U(i,j,k)-U(i,j-1,k))/dyV(j-1)+(0.5_KND)*(1-Vloc*dt/dyV(j-1))&
                        *(Uy(i,j,k)*dyPr(j)-Uy(i,j-1,k)*dyPr(j-1))/dyV(j-1)
      else
                   Uyc=(U(i,j+1,k)-U(i,j,k))/dyV(j)-(0.5_KND)*(1+Vloc*dt/dyV(j))&
                        *(Uy(i,j+1,k)*dyPr(j+1)-Uy(i,j,k)*dyPr(j))/dyV(j)
      endif
      Ut(i,j,k)=Ut(i,j,k)-Uyc*Vloc
     enddo
    enddo
   enddo
  endif

  if (gridtype==uniformgrid) then
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Wloc=(W(i,j,k-1)+W(i,j,k)+W(i+1,j,k-1)+W(i+1,j,k))/4._KND
      if (Wloc>=0) then
                   Uzc=(U(i,j,k)-U(i,j,k-1))/dzmin+(0.5_KND)*(1-Wloc*dt/dzmin)&
                        *(Uz(i,j,k)*dzPr(k)-Uz(i,j,k-1)*dzPr(k-1))/dzW(k-1)
      else
                   Uzc=(U(i,j,k+1)-U(i,j,k))/dzmin-(0.5_KND)*(1+Wloc*dt/dzmin)&
                        *(Uz(i,j,k+1)-Uz(i,j,k))
      endif
      Ut(i,j,k)=Ut(i,j,k)-Uzc*Wloc
     enddo
    enddo
   enddo
  else
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Wloc=(W(i,j,k-1)+W(i,j,k)+W(i+1,j,k-1)+W(i+1,j,k))/4._KND
      !Wloc=(W(i,j-1,k)+W(i,j,k)+W(i+1,j-1,k)+W(i+1,j,k)+W(i,j-1,k)+W(i,j,k-1)+W(i+1,j-1,k-1)+W(i+1,j,k-1))/8._KND
      if (Wloc>=0) then
                   Uzc=(U(i,j,k)-U(i,j,k-1))/dzW(k-1)+(0.5_KND)*(1-Wloc*dt/dzW(k-1))&
                        *(Uz(i,j,k)*dzPr(k)-Uz(i,j,k-1)*dzPr(k-1))/dzW(k-1)
      else
                   Uzc=(U(i,j,k+1)-U(i,j,k))/dzW(k)-(0.5_KND)*(1+Wloc*dt/dzW(k))&
                        *(Uz(i,j,k+1)*dzPr(k+1)-Uz(i,j,k)*dzPr(k))/dzW(k)
      endif
      Ut(i,j,k)=Ut(i,j,k)-Uzc*Wloc
     enddo
    enddo
   enddo
  endif

  if (associated(FirstIBPoint)) then
   IBP => FirstIBPoint
   do
    if (IBP%component==1) then
     x=IBP%x
     y=IBP%y
     z=IBP%z
     Ut(x,y,z)=Ut(x,y,z)+IBP%MSourc
    endif
    if (associated(IBP%next)) then
     IBP=>IBP%next
    else
     exit
    endif
   enddo
  endif

  do k=1,nz
   do j=1,ny
    do i=1,nx
     Ut(i,j,k)=Ut(i,j,k)*(dt/2._KND)
    enddo
   enddo
  enddo

  do k=-2,nz+2
   do j=-2,ny+2
    do i=-2,nx+2
     if (i<1.or.i>nx.or.j<1.or.j>ny.or.k<1.or.k>nz) Ut(i,j,k)=-100000
    enddo
   enddo
  enddo

  if (BtypeE==PERIODIC) then
   do k=1,nz
    do j=1,ny
     Ut(nx+1,j,k)=Ut(1,j,k)
    enddo
   enddo
  else
   do k=1,nz
    do j=1,ny
     Ut(nx+1,j,k)=0
    enddo
   enddo
 endif
 endsubroutine Utim


 subroutine Vtim(U,V,W,Pr)
  real(KND) U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND) Vxc,Vzc,Uloc,Wloc,dx2,dy2,dz2
  integer i,j,k,x,y,z,nx,ny,nz
  type(TIBPoint),pointer:: IBP
  
  nx=Unx+1
  ny=Vny
  nz=Wnz+1


  if (gridtype==uniformgrid) then
   dx2=dxmin*dxmin
   dy2=dymin*dymin
   dz2=dzmin*dzmin
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Vt(i,j,k)=((V(i+1,j,k)-2*V(i,j,k)+V(i-1,j,k))/dx2+&     
       (V(i,j+1,k)-2*V(i,j,k)+V(i,j-1,k))/dy2+&
       (V(i,j,k+1)-2*V(i,j,k)+V(i,j,k-1))/dz2)/Re
     enddo
    enddo
   enddo
  else
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Vt(i,j,k)=(0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V(i+1,j,k)-V(i,j,k))/dxU(i)-&
       0.25_KND*(Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(V(i,j,k)-V(i-1,j,k))/dxU(i-1))/dxPr(i)+&     
       (Visc(i,j+1,k)*(V(i,j+1,k)-V(i,j,k))/dyPr(j+1)-&
       Visc(i,j,k)*(V(i,j,k)-V(i,j-1,k))/dyPr(j))/dyV(j)+&
       (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))/dzW(k)-&
       0.25_KND*(Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(V(i,j,k)-V(i,j,k-1))/dzW(k-1))/dzPr(k)
     enddo
    enddo
   enddo
 endif

  do k=1,nz
   do j=1,ny
    do i=1,nx
     Vt(i,j,k)=Vt(i,j,k)-(Pr(i,j+1,k)-Pr(i,j,k))/dyV(j)
    enddo
   enddo
  enddo

  do k=1,nz
   do j=1,ny
    do i=1,nx
     Vt(i,j,k)=Vt(i,j,k)-Vy(i,j,k)*V(i,j,k)
    enddo
   enddo
  enddo


  if (gridtype==uniformgrid) then
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Uloc=(U(i-1,j,k)+U(i,j,k)+U(i-1,j+1,k)+U(i,j+1,k))/4._KND
      !Uloc=(U(i-1,j,k)+U(i,j,k)+U(i-1,j+1,k)+U(i,j+1,k)+U(i-1,j,k-1)+U(i,j,k-1)+U(i-1,j+1,k-1)+U(i,j+1,k-1))/8._KND
      if (Uloc>=0) then
                   Vxc=(V(i,j,k)-V(i-1,j,k))/dxmin+(0.5_KND)*(1-Uloc*dt/dxmin)&
                        *(Vx(i,j,k)-Vx(i-1,j,k))
      else
                   Vxc=(V(i+1,j,k)-V(i,j,k))/dxmin-(0.5_KND)*(1+Uloc*dt/dxmin)&
                        *(Vx(i+1,j,k)-Vx(i,j,k))
      endif
                  
      Vt(i,j,k)=Vt(i,j,k)-Vxc*Uloc
     enddo
    enddo
   enddo
  else 
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Uloc=(U(i-1,j,k)+U(i,j,k)+U(i-1,j+1,k)+U(i,j+1,k))/4._KND
      !Uloc=(U(i-1,j,k)+U(i,j,k)+U(i-1,j+1,k)+U(i,j+1,k)+U(i-1,j,k-1)+U(i,j,k-1)+U(i-1,j+1,k-1)+U(i,j+1,k-1))/8._KND
      if (Uloc>=0) then
                   Vxc=(V(i,j,k)-V(i-1,j,k))/dxU(i-1)+(0.5_KND)*(1-Uloc*dt/dxU(i-1))&
                        *(Vx(i,j,k)*dxPr(i)-Vx(i-1,j,k)*dxPr(i-1))/dxU(i-1)
      else
                   Vxc=(V(i+1,j,k)-V(i,j,k))/dxU(i)-(0.5_KND)*(1+Uloc*dt/dxU(i))&
                        *(Vx(i+1,j,k)*dxPr(i+1)-Vx(i,j,k)*dxPr(i))/dxU(i)
      endif
                  
      Vt(i,j,k)=Vt(i,j,k)-Vxc*Uloc
     enddo
    enddo
   enddo
  endif

  if (gridtype==uniformgrid) then
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Wloc=(W(i,j,k-1)+W(i,j,k)+W(i,j+1,k-1)+W(i,j+1,k))/4._KND
      if (Wloc>=0) then
                   Vzc=(V(i,j,k)-V(i,j,k-1))/dzmin+(0.5_KND)*(1-Wloc*dt/dzmin)&
                        *(Vz(i,j,k)-Vz(i,j,k-1))
      else
                   Vzc=(V(i,j,k+1)-V(i,j,k))/dzmin-(0.5_KND)*(1+Wloc*dt/dzmin)&
                        *(Vz(i,j,k+1)-Vz(i,j,k))
      endif
                  
      Vt(i,j,k)=Vt(i,j,k)-Vzc*Wloc
     enddo
    enddo
   enddo
  else
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Wloc=(W(i,j,k-1)+W(i,j,k)+W(i,j+1,k-1)+W(i,j+1,k))/4._KND
      !Wloc=(W(i-1,j,k)+W(i,j,k)+W(i-1,j+1,k)+W(i,j+1,k)+W(i-1,j,k-1)+W(i,j,k-1)+W(i-1,j+1,k-1)+W(i,j+1,k-1))/8._KND
      if (Wloc>=0) then
                   Vzc=(V(i,j,k)-V(i,j,k-1))/dzW(k-1)+(0.5_KND)*(1-Wloc*dt/dzW(k-1))&
                        *(Vz(i,j,k)*dzPr(k)-Vz(i,j,k-1)*dzPr(k-1))/dzW(k-1)
      else
                   Vzc=(V(i,j,k+1)-V(i,j,k))/dzW(k)-(0.5_KND)*(1+Wloc*dt/dzW(k))&
                        *(Vz(i,j,k+1)*dzPr(k+1)-Vz(i,j,k)*dzPr(k))/dzW(k)
      endif
                  
      Vt(i,j,k)=Vt(i,j,k)-Vzc*Wloc
     enddo
    enddo
   enddo
  endif
 
  if (associated(FirstIBPoint)) then
   IBP => FirstIBPoint
   do
    if (IBP%component==2) then
     x=IBP%x
     y=IBP%y
     z=IBP%z
     Vt(x,y,z)=Vt(x,y,z)+IBP%MSourc
    endif
    if (associated(IBP%next)) then
     IBP=>IBP%next
    else
     exit
    endif
   enddo
  endif

  do k=1,nz
   do j=1,ny
    do i=1,nx
     Vt(i,j,k)=Vt(i,j,k)*(dt/2._KND)
    enddo
   enddo
  enddo

  do k=-2,nz+2
   do j=-2,ny+2
    do i=-2,nx+2
     if (i<1.or.i>nx.or.j<1.or.j>ny.or.k<1.or.k>nz) Vt(i,j,k)=-100000
    enddo
   enddo
  enddo
  
  if (BtypeN==PERIODIC) then
   do k=1,nz
    do i=1,nx
     Vt(i,ny+1,k)=Vt(i,1,k)
    enddo
   enddo
  else
   do k=1,nz
    do i=1,nx
     Vt(i,ny+1,k)=0
    enddo
   enddo
  endif
 endsubroutine Vtim


 subroutine Wtim(U,V,W,Pr)
  real(KND) U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND) Wxc,Wyc,Uloc,Vloc,dx2,dy2,dz2
  integer i,j,k,x,y,z,nx,ny,nz
  type(TIBPoint),pointer:: IBP
  
  nx=Unx+1
  ny=Vny+1
  nz=Wnz

  if (gridtype==uniformgrid) then
   dx2=dxmin*dxmin
   dy2=dymin*dymin
   dz2=dzmin*dzmin
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Wt(i,j,k)=((W(i+1,j,k)-2*W(i,j,k)+W(i-1,j,k))/dx2+&     
        (W(i,j+1,k)-2*W(i,j,k)+W(i,j-1,k))/dy2+&
        (W(i,j,k+1)-2*W(i,j,k)+W(i,j,k-1))/dz2)/Re
     enddo
    enddo
   enddo
  else
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Wt(i,j,k)=(0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W(i+1,j,k)-W(i,j,k))/dxU(i)-&
        0.25_KND*(Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(W(i,j,k)-W(i-1,j,k))/dxU(i-1))/dxPr(i)+&     
        (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W(i,j+1,k)-W(i,j,k))/dyV(j)-&
        0.25_KND*(Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(W(i,j,k)-W(i,j-1,k))/dyV(j-1))/dyPr(j)+&
        (Visc(i,j,k+1)*(W(i,j,k+1)-W(i,j,k))/dzPr(k+1)-&
        Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1))/dzPr(k))/dzW(k)
     enddo
    enddo
   enddo
  endif


  do k=1,nz
   do j=1,ny
    do i=1,nx
     Wt(i,j,k)=Wt(i,j,k)-(Pr(i,j,k+1)-Pr(i,j,k))/dzW(k)
    enddo
   enddo
  enddo

  do k=1,nz
   do j=1,ny
    do i=1,nx
      Wt(i,j,k)=Wt(i,j,k)-Wz(i,j,k)*W(i,j,k)
    enddo
   enddo
  enddo

  if (gridtype==uniformgrid) then
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Uloc=(U(i-1,j,k)+U(i,j,k)+U(i-1,j,k+1)+U(i,j,k+1))/4._KND
      if (Uloc>=0) then
                   Wxc=(W(i,j,k)-W(i-1,j,k))/dxmin+(0.5_KND)*(1-Uloc*dt/dxmin)&
                        *(Wx(i,j,k)-Wx(i-1,j,k))
      else
                   Wxc=(W(i+1,j,k)-W(i,j,k))/dxmin-(0.5_KND)*(1+Uloc*dt/dxmin)&
                        *(Wx(i+1,j,k)-Wx(i,j,k))
      endif
                  
      Wt(i,j,k)=Wt(i,j,k)-Wxc*Uloc
     enddo
    enddo
   enddo
  else
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Uloc=(U(i-1,j,k)+U(i,j,k)+U(i-1,j,k+1)+U(i,j,k+1))/4._KND
      !Uloc=(U(i-1,j,k)+U(i,j,k)+U(i-1,j,k+1)+U(i,j,k+1)+U(i-1,j-1,k)+U(i,j-1,k)+U(i-1,j-1,k+1)+U(i,j-1,k+1))/8._KND
      if (Uloc>=0) then
                   Wxc=(W(i,j,k)-W(i-1,j,k))/dxU(i-1)+(0.5_KND)*(1-Uloc*dt/dxU(i-1))&
                        *(Wx(i,j,k)*dxPr(i)-Wx(i-1,j,k)*dxPr(i-1))/dxU(i-1)
      else
                   Wxc=(W(i+1,j,k)-W(i,j,k))/dxU(i)-(0.5_KND)*(1+Uloc*dt/dxU(i))&
                        *(Wx(i+1,j,k)*dxPr(i+1)-Wx(i,j,k)*dxPr(i))/dxU(i)
      endif
                  
      Wt(i,j,k)=Wt(i,j,k)-Wxc*Uloc
     enddo
    enddo
   enddo
  endif

  if (gridtype==uniformgrid) then
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Vloc=(V(i,j-1,k)+V(i,j,k)+V(i,j-1,k+1)+V(i,j,k+1))/4._KND
      !Vloc=(V(i,j-1,k)+V(i,j,k)+V(i,j-1,k+1)+V(i,j,k+1)+V(i,j-1,k)+V(i-1,j,k)+V(i-1,j-1,k+1)+V(i-1,j,k+1))/8._KND
      if (Vloc>=0) then
                   Wyc=(W(i,j,k)-W(i,j-1,k))/dyV(j-1)+(0.5_KND)*(1-Vloc*dt/dyV(j-1))&
                        *(Wy(i,j,k)*dyPr(j)-Wy(i,j-1,k)*dyPr(j-1))/dyV(j-1)
      else
                   Wyc=(W(i,j+1,k)-W(i,j,k))/dyV(j)-(0.5_KND)*(1+Vloc*dt/dyV(j))&
                        *(Wy(i,j+1,k)*dyPr(j+1)-Wy(i,j,k)*dyPr(j))/dyV(j)
      endif
     
      Wt(i,j,k)=Wt(i,j,k)-Wyc*Vloc
     enddo
    enddo
   enddo
  else
   do k=1,nz
    do j=1,ny
     do i=1,nx
      Vloc=(V(i,j-1,k)+V(i,j,k)+V(i,j-1,k+1)+V(i,j,k+1))/4._KND
      !Vloc=(V(i,j-1,k)+V(i,j,k)+V(i,j-1,k+1)+V(i,j,k+1)+V(i,j-1,k)+V(i-1,j,k)+V(i-1,j-1,k+1)+V(i-1,j,k+1))/8._KND
      if (Vloc>=0) then
                   Wyc=(W(i,j,k)-W(i,j-1,k))/dymin+(0.5_KND)*(1-Vloc*dt/dymin)&
                        *(Wy(i,j,k)-Wy(i,j-1,k))
      else
                   Wyc=(W(i,j+1,k)-W(i,j,k))/dymin-(0.5_KND)*(1+Vloc*dt/dymin)&
                        *(Wy(i,j+1,k)-Wy(i,j,k))
      endif
     
      Wt(i,j,k)=Wt(i,j,k)-Wyc*Vloc
     enddo
    enddo
   enddo
  endif

  if (associated(FirstIBPoint)) then
   IBP => FirstIBPoint
   do
    if (IBP%component==3) then
     x=IBP%x
     y=IBP%y
     z=IBP%z
     Wt(x,y,z)=Wt(x,y,z)+IBP%MSourc
    endif
    if (associated(IBP%next)) then
     IBP=>IBP%next
    else
     exit
    endif
   enddo
  endif

  do k=1,nz
   do j=1,ny
    do i=1,nx
     Wt(i,j,k)=Wt(i,j,k)*(dt/2._KND)
    enddo
   enddo
  enddo

  do k=-2,nz+2
   do j=-2,ny+2
    do i=-2,nx+2
     if (i<1.or.i>nx.or.j<1.or.j>ny.or.k<1.or.k>nz) Wt(i,j,k)=-100000
    enddo
   enddo
  enddo

  if (BtypeT==PERIODIC) then
   do j=1,ny
    do i=1,nx
     Wt(i,j,nz+1)=Wt(i,j,1)
    enddo
   enddo
  else
   do j=1,ny
    do i=1,nx
     Wt(i,j,nz+1)=0
    enddo
   enddo
 endif
 endsubroutine Wtim


 subroutine Corners(U,V,W)
 real(KND),dimension(-2:,-2:,-2:):: U,V,W
 real(KND) USB,UST,UNB,UNT
 real(KND) VWB,VWT,VEB,VET
 real(KND) WWS,WWN,WES,WEN
 integer i,j,k
 integer normalupwind,xupw,yupw,zupw


 
 normalupwind=1
 if (normalupwind==0) then
  do k=1,Wnz
   do j=1,Vny
    do i=1,Unx
     USB=U(i,j,k)+(yV(j)-yPr(j))*Uy(i,j,k)+(zW(k)-zPr(k))*Uz(i,j,k)+Ut(i,j,k)
     UST=U(i,j,k+1)+(yV(j)-yPr(j))*Uy(i,j,k+1)-(zPr(k+1)-zW(k))*Uz(i,j,k+1)+Ut(i,j,k+1)
     UNB=U(i,j+1,k)-(yPr(j+1)-yV(j))*Uy(i,j+1,k)+(zW(k)-zPr(k))*Uz(i,j+1,k)+Ut(i,j+1,k)
     UNT=U(i,j+1,k+1)-(yPr(j+1)-yV(j))*Uy(i,j+1,k+1)-(zPr(k+1)-zW(k))*Uz(i,j+1,k+1)+Ut(i,j+1,k+1)

     VWB=V(i,j,k)+(xU(i)-xPr(i))*Vx(i,j,k)+(zW(k)-zPr(k))*Vz(i,j,k)+Vt(i,j,k)
     VWT=V(i,j,k+1)+(xU(i)-xPr(i))*Vx(i,j,k+1)-(zPr(k+1)-zW(k))*Vz(i,j,k+1)+Vt(i,j,k+1)
     VEB=V(i+1,j,k)-(xPr(i+1)-xU(i))*Vx(i+1,j,k)+(zW(k)-zPr(k))*Vz(i+1,j,k)+Vt(i+1,j,k)
     VET=V(i+1,j,k+1)-(xPr(i+1)-xU(i))*Vx(i+1,j,k+1)-(zPr(k+1)-zW(k))*Vz(i+1,j,k+1)+Vt(i+1,j,k+1)

     WWS=W(i,j,k)+(xU(i)-xPr(i))*Wx(i,j,k)+(yV(j)-yPr(j))*Wy(i,j,k)+Wt(i,j,k)
     WWN=W(i,j+1,k)+(xU(i)-xPr(i))*Wx(i,j+1,k)-(yPr(j+1)-yV(j))*Wy(i,j+1,k)+Wt(i,j+1,k)
     WES=W(i+1,j,k)-(xPr(i+1)-xU(i))*Wx(i+1,j,k)+(yV(j)-yPr(j))*Wy(i+1,j,k)+Wt(i+1,j,k)
     WEN=W(i+1,j+1,k)-(xPr(i+1)-xU(i))*Wx(i+1,j+1,k)-(yPr(j+1)-yV(j))*Wy(i+1,j+1,k)+Wt(i+1,j+1,k)
   
     Uco(i,j,k)=(USB+UST+UNB+UNT)/4._KND
     Vco(i,j,k)=(VWB+VWT+VEB+VET)/4._KND
     Wco(i,j,k)=(WWS+WWN+WES+WEN)/4._KND
    enddo
   enddo
  enddo
 else
  do k=1,Wnz
   do j=1,Vny
    do i=1,Unx
     USB=U(i,j,k)+(yV(j)-yPr(j))*Uy(i,j,k)+(zW(k)-zPr(k))*Uz(i,j,k)+Ut(i,j,k)
     UST=U(i,j,k+1)+(yV(j)-yPr(j))*Uy(i,j,k+1)-(zPr(k+1)-zW(k))*Uz(i,j,k+1)+Ut(i,j,k+1)
     UNB=U(i,j+1,k)-(yPr(j+1)-yV(j))*Uy(i,j+1,k)+(zW(k)-zPr(k))*Uz(i,j+1,k)+Ut(i,j+1,k)
     UNT=U(i,j+1,k+1)-(yPr(j+1)-yV(j))*Uy(i,j+1,k+1)-(zPr(k+1)-zW(k))*Uz(i,j+1,k+1)+Ut(i,j+1,k+1)

     VWB=V(i,j,k)+(xU(i)-xPr(i))*Vx(i,j,k)+(zW(k)-zPr(k))*Vz(i,j,k)+Vt(i,j,k)
     VWT=V(i,j,k+1)+(xU(i)-xPr(i))*Vx(i,j,k+1)-(zPr(k+1)-zW(k))*Vz(i,j,k+1)+Vt(i,j,k+1)
     VEB=V(i+1,j,k)-(xPr(i+1)-xU(i))*Vx(i+1,j,k)+(zW(k)-zPr(k))*Vz(i+1,j,k)+Vt(i+1,j,k)
     VET=V(i+1,j,k+1)-(xPr(i+1)-xU(i))*Vx(i+1,j,k+1)-(zPr(k+1)-zW(k))*Vz(i+1,j,k+1)+Vt(i+1,j,k+1)

     WWS=W(i,j,k)+(xU(i)-xPr(i))*Wx(i,j,k)+(yV(j)-yPr(j))*Wy(i,j,k)+Wt(i,j,k)
     WWN=W(i,j+1,k)+(xU(i)-xPr(i))*Wx(i,j+1,k)-(yPr(j+1)-yV(j))*Wy(i,j+1,k)+Wt(i,j+1,k)
     WES=W(i+1,j,k)-(xPr(i+1)-xU(i))*Wx(i+1,j,k)+(yV(j)-yPr(j))*Wy(i+1,j,k)+Wt(i+1,j,k)
     WEN=W(i+1,j+1,k)-(xPr(i+1)-xU(i))*Wx(i+1,j+1,k)-(yPr(j+1)-yV(j))*Wy(i+1,j+1,k)+Wt(i+1,j+1,k)
   
     if (USB>0.and.UST>0.and.UNB>0.and.UNT>0) then
      xupw=-1
     elseif (USB<0.and.UST<0.and.UNB<0.and.UNT<0) then
      xupw=+1
     else
      xupw=0
     endif

     if (VWB>0.and.VWT>0.and.VEB>0.and.VET>0) then
      yupw=-1
     elseif (VWB<0.and.VWT<0.and.VEB<0.and.VET<0) then
      yupw=+1
     else
      yupw=0
     endif

     if (WWS>0.and.WWN>0.and.WES>0.and.WEN>0) then
      zupw=-1
     elseif (WWS<0.and.WWN<0.and.WES<0.and.WEN<0) then
      zupw=+1
     else
      zupw=0
     endif

     if (yupw==-1.and.zupw==-1) then
      Uco(i,j,k)=USB
     elseif (yupw==-1.and.zupw==+1) then
      Uco(i,j,k)=UST
     elseif (yupw==+1.and.zupw==-1) then
      Uco(i,j,k)=UNB
     elseif (yupw==+1.and.zupw==+1) then
      Uco(i,j,k)=UNT
     elseif (yupw==-1) then
      Uco(i,j,k)=(USB+UST)/2._KND
     elseif (yupw==+1) then
      Uco(i,j,k)=(UNB+UNT)/2._KND
     elseif (zupw==-1) then
      Uco(i,j,k)=(USB+UNB)/2._KND
     elseif (zupw==+1) then
      Uco(i,j,k)=(UST+UNT)/2._KND
     else
      Uco(i,j,k)=(USB+UST+UNB+UNT)/4._KND
     endif

     if (xupw==-1.and.zupw==-1) then
      Vco(i,j,k)=VWB
     elseif (xupw==-1.and.zupw==+1) then
      Vco(i,j,k)=VWT
     elseif (xupw==+1.and.zupw==-1) then
      Vco(i,j,k)=VEB
     elseif (xupw==+1.and.zupw==+1) then
      Vco(i,j,k)=VET
     elseif (xupw==-1) then
      Vco(i,j,k)=(VWB+VWT)/2._KND
     elseif (xupw==+1) then
      Vco(i,j,k)=(VEB+VET)/2._KND
     elseif (zupw==-1) then
      Vco(i,j,k)=(VWB+VEB)/2._KND
     elseif (zupw==+1) then
      Vco(i,j,k)=(VWT+VET)/2._KND
     else
      Vco(i,j,k)=(VWB+VWT+VEB+VET)/4._KND
     endif

     if (xupw==-1.and.yupw==-1) then
      Wco(i,j,k)=WWS
     elseif (xupw==-1.and.yupw==+1) then
      Wco(i,j,k)=WWN
     elseif (xupw==+1.and.yupw==-1) then
      Wco(i,j,k)=WES
     elseif (xupw==+1.and.yupw==+1) then
      Wco(i,j,k)=WEN
     elseif (xupw==-1) then
      Wco(i,j,k)=(WWS+WWN)/2._KND
     elseif (xupw==+1) then
      Wco(i,j,k)=(WES+WEN)/2._KND
     elseif (yupw==-1) then
      Wco(i,j,k)=(WWS+WES)/2._KND
     elseif (yupw==+1) then
      Wco(i,j,k)=(WWN+WEN)/2._KND
     else
      Wco(i,j,k)=(WWS+WWN+WES+WEN)/4._KND
     endif
    enddo
   enddo
  enddo
 endif
 
 call Bound_CondUco(Uco)
 call Bound_CondVco(Vco)
 call Bound_CondWco(Wco)
 endsubroutine Corners

  
 subroutine TauU(U2,U,V,W)
 real(KND),dimension(-2:,-2:,-2:):: U2,U,V,W
 real(KND),dimension(LBOUND(U,1):UBOUND(U,1),LBOUND(U,2):UBOUND(U,2),LBOUND(U,3):UBOUND(U,3))::Uce
 real(KND) UL,UR
 integer i,j,k
 
 Uce=10000
 
 do k=1,Unz
  i=1            
  do j=1,Uny
    if (BtypeW==PERIODIC) then
       UL=U(i-1,j,k)+(xPr(i)-xU(i-1))*Ux(i-1,j,k)+Ut(Unx,j,k)
    else
       UL=U(i-1,j,k)+(xPr(i)-xU(i-1))*Ux(i-1,j,k)
    endif   
    UR=U(i,j,k)-(xU(i)-xPr(i))*Ux(i,j,k)+Ut(i,j,k) 
    if ((UL>=0).and.(UL+UR>=0)) then
                                  Uce(i,j,k)=UL
    elseif ((UL<0).and.(UR>0)) then
                                  Uce(i,j,k)=(UL+UR)/2._KND
    else
                                  Uce(i,j,k)=UR
    endif
  enddo
 enddo

 do k=1,Unz
  do j=1,Uny
   do i=2,Unx
    UL=U(i-1,j,k)+(xPr(i)-xU(i-1))*Ux(i-1,j,k)+Ut(i-1,j,k)
    UR=U(i,j,k)-(xU(i)-xPr(i))*Ux(i,j,k)+Ut(i,j,k)
    if ((UL>=0).and.(UL+UR>=0)) then
                                  Uce(i,j,k)=UL
    elseif ((UL<0).and.(UR>0)) then
                                  Uce(i,j,k)=(UL+UR)/2._KND
    else
                                  Uce(i,j,k)=UR
    endif
   enddo
  enddo
 enddo

 do k=1,Unz
  i=Unx+1 
  do j=1,Uny
    UL=U(i-1,j,k)+(xPr(i)-xU(i-1))*Ux(i-1,j,k)+Ut(i-1,j,k)
    if (BtypeE==PERIODIC) then
       UR=U(i,j,k)-(xU(i)-xPr(i))*Ux(i,j,k)+Ut(1,j,k)
    else
       UR=U(i,j,k)-(xU(i)-xPr(i))*Ux(i,j,k)
    endif   
    if ((UL>=0).and.(UL+UR>=0)) then
                                  Uce(i,j,k)=UL
    elseif ((UL<0).and.(UR>0)) then
                                  Uce(i,j,k)=(UL+UR)/2._KND
    else
                                  Uce(i,j,k)=UR
    endif
  enddo
 enddo

 U2(1:Unx,1:Uny,1:Unz)=0
 do k=1,Unz
  do j=1,Uny
   do i=1,Unx
    U2(i,j,k)=U2(i,j,k)-dt*((Uce(i+1,j,k)+Uce(i,j,k))/2._KND)*(Uce(i+1,j,k)-Uce(i,j,k))/dxU(i)
   enddo
  enddo  
 enddo
 do k=1,Unz
  do j=1,Uny
   do i=1,Unx
    U2(i,j,k)=U2(i,j,k)-dt*((Vco(i,j,k)+Vco(i,j-1,k)+Vco(i,j,k-1)+Vco(i,j-1,k-1))/4._KND)*&
                           ((Uco(i,j,k)+Uco(i,j,k-1)-Uco(i,j-1,k)-Uco(i,j-1,k-1))/2._KND)/dyPr(j)
   enddo
  enddo
 enddo
 do k=1,Unz
  do j=1,Uny
   do i=1,Unx
    U2(i,j,k)=U2(i,j,k)-dt*((Wco(i,j,k)+Wco(i,j-1,k)+Wco(i,j,k-1)+Wco(i,j-1,k-1))/4._KND)*&
                           ((Uco(i,j,k)+Uco(i,j-1,k)-Uco(i,j,k-1)-Uco(i,j-1,k-1))/2._KND)/dzPr(k)
   enddo
  enddo
 enddo
 endsubroutine TauU


 subroutine TauV(V2,U,V,W)
 real(KND),dimension(-2:,-2:,-2:):: V2,U,V,W
 real(KND),dimension(LBOUND(V,1):UBOUND(V,1),LBOUND(V,2):UBOUND(V,2),LBOUND(V,3):UBOUND(V,3))::Vce
 real(KND) VL,VR
 integer i,j,k

 Vce=10000


 do k=1,Vnz
  j=1
  do i=1,Vnx
    if (BtypeS==PERIODIC) then
       VL=V(i,j-1,k)+(yPr(j)-yV(j-1))*Vy(i,j-1,k)+Vt(i,Vny,k)
    else
       VL=V(i,j-1,k)+(yPr(j)-yV(j-1))*Vy(i,j-1,k) 
    endif
    VR=V(i,j,k)-(yV(j)-yPr(j))*Vy(i,j,k)+Vt(i,j,k)
    if ((VL>=0).and.(VL+VR>=0)) then
                                  Vce(i,j,k)=VL
    elseif ((VL<0).and.(VR>0)) then
                                  Vce(i,j,k)=(VL+VR)/2._KND
    else
                                  Vce(i,j,k)=VR
    endif
  enddo
 enddo
 do k=1,Vnz
  do j=2,Vny
   do i=1,Vnx
    VL=V(i,j-1,k)+(yPr(j)-yV(j-1))*Vy(i,j-1,k)+Vt(i,j-1,k)
    VR=V(i,j,k)-(yV(j)-yPr(j))*Vy(i,j,k)+Vt(i,j,k)
    if ((VL>=0).and.(VL+VR>=0)) then
                                  Vce(i,j,k)=VL
    elseif ((VL<0).and.(VR>0)) then
                                  Vce(i,j,k)=(VL+VR)/2._KND
    else
                                  Vce(i,j,k)=VR
    endif
   enddo
  enddo
 enddo
 do k=1,Vnz
  j=Vny+1
  do i=1,Vnx
    VL=V(i,j-1,k)+(yPr(j)-yV(j-1))*Vy(i,j-1,k)+Vt(i,j-1,k) 
    if (BtypeN==PERIODIC) then
       VR=V(i,j,k)-(yV(j)-yPr(j))*Vy(i,j,k)+Vt(i,1,k)
    else
       VR=V(i,j,k)-(yV(j)-yPr(j))*Vy(i,j,k)
    endif   
    if ((VL>=0).and.(VL+VR>=0)) then
                                  Vce(i,j,k)=VL
    elseif ((VL<0).and.(VR>0)) then
                                  Vce(i,j,k)=(VL+VR)/2._KND
    else
                                  Vce(i,j,k)=VR
    endif
  enddo
 enddo
 V2(1:Vnx,1:Vny,1:Vnz)=0
 do k=1,Vnz
  do j=1,Vny
   do i=1,Vnx
    V2(i,j,k)=V2(i,j,k)-dt*((Vce(i,j+1,k)+Vce(i,j,k))/2._KND)*(Vce(i,j+1,k)-Vce(i,j,k))/dyV(j)
   enddo
  enddo
 enddo 
 do k=1,Vnz
  do j=1,Vny
   do i=1,Vnx
    V2(i,j,k)=V2(i,j,k)-dt*((Uco(i-1,j,k)+Uco(i,j,k-1)+Uco(i-1,j,k-1)+Uco(i,j,k))/4._KND)*&
                           ((Vco(i,j,k)+Vco(i,j,k-1)-Vco(i-1,j,k)-Vco(i-1,j,k-1))/2._KND)/dxPr(i)
   enddo
  enddo
 enddo
 do k=1,Vnz
  do j=1,Vny
   do i=1,Vnx
    V2(i,j,k)=V2(i,j,k)-dt*((Wco(i-1,j,k)+Wco(i,j,k-1)+Wco(i-1,j,k-1)+Wco(i,j,k))/4._KND)*&
                           ((Vco(i,j,k)+Vco(i-1,j,k)-Vco(i,j,k-1)-Vco(i-1,j,k-1))/2._KND)/dzPr(k)
   enddo
  enddo
 enddo
 endsubroutine TauV



 subroutine TauW(W2,U,V,W)
 real(KND),dimension(-2:,-2:,-2:):: W2,U,V,W
 real(KND),dimension(LBOUND(W,1):UBOUND(W,1),LBOUND(W,2):UBOUND(W,2),LBOUND(W,3):UBOUND(W,3))::Wce
 real(KND) WL,WR
 integer i,j,k

 Wce=10000


 do j=1,Wny
  k=1
  do i=1,Wnx
    if (BtypeB==PERIODIC) then    
       WL=W(i,j,k-1)+(zPr(k)-zW(k-1))*Wz(i,j,k-1)+Wt(Wnx,j,k)
    else
       WL=W(i,j,k-1)+(zPr(k)-zW(k-1))*Wz(i,j,k-1)
    endif
    WR=W(i,j,k)-(zW(k)-zPr(k))*Wz(i,j,k)+Wt(i,j,k)
    if ((WL>=0).and.(WL+WR>=0)) then
                                  Wce(i,j,k)=WL
    elseif ((WL<0).and.(WR>0)) then
                                  Wce(i,j,k)=(WL+WR)/2._KND
    else
                                  Wce(i,j,k)=WR
    endif
  enddo
 enddo
 do k=2,Wnz
  do j=1,Wny
   do i=1,Wnx
    WL=W(i,j,k-1)+(zPr(k)-zW(k-1))*Wz(i,j,k-1)+Wt(i,j,k-1)
    WR=W(i,j,k)-(zW(k)-zPr(k))*Wz(i,j,k)+Wt(i,j,k)
    if ((WL>=0).and.(WL+WR>=0)) then
                                  Wce(i,j,k)=WL
    elseif ((WL<0).and.(WR>0)) then
                                  Wce(i,j,k)=(WL+WR)/2._KND
    else
                                  Wce(i,j,k)=WR
    endif
   enddo
  enddo
 enddo
 do j=1,Wny
  k=Wnz+1
  do i=1,Wnx
    WL=W(i,j,k-1)+(zPr(k)-zW(k-1))*Wz(i,j,k-1)+Wt(i,j,k-1)
    if (BtypeT==PERIODIC) then
       WR=W(i,j,k)-(zW(k)-zPr(k))*Wz(i,j,k)+Wt(i,j,1)
    else
       WR=W(i,j,k)-(zW(k)-zPr(k))*Wz(i,j,k)
    endif   
    if ((WL>=0).and.(WL+WR>=0)) then
                                  Wce(i,j,k)=WL
    elseif ((WL<0).and.(WR>0)) then
                                  Wce(i,j,k)=(WL+WR)/2._KND
    else
                                  Wce(i,j,k)=WR
    endif
  enddo
 enddo
 W2(1:Wnx,1:Wny,1:Wnz)=0
 do k=1,Wnz
  do j=1,Wny
   do i=1,Wnx
    W2(i,j,k)=W2(i,j,k)-dt*((Wce(i,j,k+1)+Wce(i,j,k))/2._KND)*(Wce(i,j,k+1)-Wce(i,j,k))/dzW(k)
   enddo
  enddo
 enddo 
 do k=1,Wnz
  do j=1,Wny
   do i=1,Wnx
    W2(i,j,k)=W2(i,j,k)-dt*((Uco(i-1,j,k)+Uco(i,j-1,k)+Uco(i-1,j-1,k)+Uco(i,j,k))/4._KND)*&
                           ((Wco(i,j,k)+Wco(i,j-1,k)-Wco(i-1,j,k)-Wco(i-1,j-1,k))/2._KND)/dxPr(i)
   enddo
  enddo
 enddo
 do k=1,Wnz
  do j=1,Wny
   do i=1,Wnx
    W2(i,j,k)=W2(i,j,k)-dt*((Vco(i-1,j,k)+Vco(i,j-1,k)+Vco(i-1,j-1,k)+Vco(i,j,k))/4._KND)*&
                           ((Wco(i,j,k)+Wco(i-1,j,k)-Wco(i,j-1,k)-Wco(i-1,j-1,k))/2._KND)/dyPr(j)
   enddo
  enddo
 enddo
 endsubroutine TauW






  subroutine Bound_CondUco(Uco)
  real(KND),dimension(-2:,-2:,-2:):: Uco
  integer i,j,k,nx,ny,nz

  nx=Unx
  ny=Vny
  nz=Wnz   
 
  if (BtypeS==DIRICHLET) then
   do k=1,nz
    do i=1,nx                       !Dirichlet inlet
     Uco(i,0,k)=SsideU
     Uco(i,-1,k)=SsideU+(SsideU-Uco(i,1,k))
     Uco(i,-2,k)=SsideU+(SsideU-Uco(i,2,k))
    enddo       
   enddo   
  elseif (BtypeS==NOSLIP) then
   do k=1,nz
    do i=1,nx                       !Solid wall
     Uco(i,0,k)=0
     Uco(i,-1,k)=-Uco(i,1,k)
     Uco(i,-2,k)=-Uco(i,2,k)
    enddo
   enddo
  elseif (BtypeS==NEUMANN) then
   do k=1,nz
    do i=1,nx                       !Neumann inlet
     Uco(i,0,k)=Uco(i,1,k)
     Uco(i,-1,k)=Uco(i,1,k)
     Uco(i,-2,k)=Uco(i,1,k)
    enddo
   enddo
  elseif (BtypeS==FREESLIP) then  !FREESLIP
   do k=1,nz
    do i=1,nx
     Uco(i,0,k)=Uco(i,1,k)
     Uco(i,-1,k)=Uco(i,1,k)
     Uco(i,-2,k)=Uco(i,1,k)
    enddo
   enddo
  elseif (BtypeS==PERIODIC) then  !Periodic BC
   do k=1,nz
    do i=1,nx
     Uco(i,0,k)=Uco(i,ny,k)
     Uco(i,-1,k)=Uco(i,ny-1,k)
     Uco(i,-2,k)=Uco(i,ny-2,k)
    enddo
   enddo      
  endif   

  if (BtypeN==DIRICHLET) then
   do k=1,nz
    do i=1,nx                       !Dirichlet inlet
     Uco(i,ny+1,k)=NsideU
     Uco(i,ny+2,k)=NsideU+(NsideU-Uco(i,ny,k))
     Uco(i,ny+3,k)=NsideU+(NsideU-Uco(i,ny-1,k))
    enddo       
   enddo   
  elseif (BtypeN==NOSLIP) then
   do k=1,nz
    do i=1,nx                       !Solid wall
     Uco(i,ny+1,k)=0
     Uco(i,ny+2,k)=-Uco(i,ny,k)
     Uco(i,ny+3,k)=-Uco(i,ny-1,k)
    enddo
   enddo
  elseif (BtypeN==NEUMANN) then
   do k=1,nz
    do i=1,nx                       !Neumann inlet
     Uco(i,ny+1,k)=Uco(i,ny,k)
     Uco(i,ny+2,k)=Uco(i,ny,k)
     Uco(i,ny+3,k)=Uco(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==FREESLIP) then  !FREESLIP
   do k=1,nz
    do i=1,nx
     Uco(i,ny+1,k)=Uco(i,ny,k)
     Uco(i,ny+2,k)=Uco(i,ny,k)
     Uco(i,ny+3,k)=Uco(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==PERIODIC) then  !Periodic BC
   do k=1,nz
    do i=1,nx
     Uco(i,ny+1,k)=Uco(i,1,k)
     Uco(i,ny+2,k)=Uco(i,2,k)
     Uco(i,ny+3,k)=Uco(i,3,k)
    enddo
   enddo      
  endif

  if (BtypeB==DIRICHLET) then
   do j=1,ny
    do i=1,nx                       !Dirichlet inlet
     Uco(i,j,0)=BsideU
     Uco(i,j,-1)=BsideU+(BsideU-Uco(i,j,1))
     Uco(i,j,-2)=BsideU+(BsideU-Uco(i,j,2))
    enddo       
   enddo   
  elseif (BtypeB==NOSLIP) then
   do j=1,ny
    do i=1,nx                       !Solid wall
     Uco(i,j,0)=0
     Uco(i,j,-1)=-Uco(i,j,1)
     Uco(i,j,-2)=-Uco(i,j,2)
    enddo
   enddo
  elseif (BtypeB==NEUMANN) then
   do j=1,ny
    do i=1,nx                       !Neumann inlet
     Uco(i,j,0)=Uco(i,j,1)
     Uco(i,j,-1)=Uco(i,j,1)
     Uco(i,j,-2)=Uco(i,j,1)
    enddo
   enddo
  elseif (BtypeB==FREESLIP) then  !FREESLIP
   do j=1,ny
    do i=1,nx
     Uco(i,j,0)=Uco(i,j,1)
     Uco(i,j,-1)=Uco(i,j,1)
     Uco(i,j,-2)=Uco(i,j,1)
    enddo
   enddo
  elseif (BtypeB==PERIODIC) then  !Periodic BC
   do j=1,ny
    do i=1,nx
     Uco(i,j,0)=Uco(i,j,nz)
     Uco(i,j,-1)=Uco(i,j,nz-1)
     Uco(i,j,-2)=Uco(i,j,nz-2)
    enddo
   enddo      
  endif   
  
  if (BtypeT==DIRICHLET) then
   do j=1,ny
    do i=1,nx                       !Dirichlet inlet
     Uco(i,j,nz+1)=TsideU
     Uco(i,j,nz+2)=TsideU+(TsideU-Uco(i,j,nz))
     Uco(i,j,nz+3)=TsideU+(TsideU-Uco(i,j,nz-1))
    enddo       
   enddo   
  elseif (BtypeT==NOSLIP) then
   do j=1,ny
    do i=1,nx                       !Solid wall
     Uco(i,j,nz+1)=0
     Uco(i,j,nz+2)=-Uco(i,j,nz)
     Uco(i,j,nz+3)=-Uco(i,j,nz-1)
    enddo
   enddo
  elseif (BtypeT==NEUMANN) then
   do j=1,ny
    do i=1,nx                       !Neumann inlet
     Uco(i,j,nz+1)=Uco(i,j,nz)
     Uco(i,j,nz+2)=Uco(i,j,nz)
     Uco(i,j,nz+3)=Uco(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==FREESLIP) then  !FREESLIP
   do j=1,ny
    do i=1,nx
     Uco(i,j,nz+1)=Uco(i,j,nz)
     Uco(i,j,nz+2)=Uco(i,j,nz)
     Uco(i,j,nz+3)=Uco(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==PERIODIC) then  !Periodic BC
   do j=1,ny
    do i=1,nx
     Uco(i,j,nz+1)=Uco(i,j,1)
     Uco(i,j,nz+2)=Uco(i,j,2)
     Uco(i,j,nz+3)=Uco(i,j,3)
    enddo
   enddo      
  endif
   
  if (BtypeW==DIRICHLET) then
   do k=1,nz
    do j=1,ny                       !Dirichlet inlet
     Uco(0,j,k)=(Uin(j,k)+Uin(j+1,k)+Uin(j,k+1)+Uin(j+1,k+1))/4._KND
     Uco(-1,j,k)=(Uin(j,k)+Uin(j+1,k)+Uin(j,k+1)+Uin(j+1,k+1))/4._KND
     Uco(-2,j,k)=(Uin(j,k)+Uin(j+1,k)+Uin(j,k+1)+Uin(j+1,k+1))/4._KND
    enddo       
   enddo
  elseif (BtypeW==NOSLIP) then
   do k=1,nz
    do j=1,ny                       !Solid wall
     Uco(0,j,k)=0
     Uco(-1,j,k)=-Uco(1,j,k)
     Uco(-2,j,k)=-Uco(2,j,k)
    enddo
   enddo
  elseif (BtypeW==NEUMANN) then
   do k=1,nz
    do j=1,ny                       !Neumann inlet
     Uco(0,j,k)=Uco(1,j,k)
     Uco(-1,j,k)=Uco(1,j,k)
     Uco(-2,j,k)=Uco(1,j,k)
    enddo
   enddo
  elseif (BtypeW==FREESLIP) then  !FREESLIP
   do k=1,nz
    do j=1,ny                       
     Uco(0,j,k)=0
     Uco(-1,j,k)=-Uco(1,j,k)
     Uco(-2,j,k)=-Uco(2,j,k)
    enddo
   enddo
  elseif (BtypeW==PERIODIC) then  !Periodic BC
   do k=1,nz
    do j=1,ny                       
     Uco(0,j,k)=Uco(nx,j,k)
     Uco(-1,j,k)=Uco(nx-1,j,k)
     Uco(-2,j,k)=Uco(nx-2,j,k)
    enddo
   enddo      
  endif      
     
  if (BtypeE==DIRICHLET) then
   do k=1,nz
    do j=1,ny                       !Dirichlet inlet
     Uco(nx+1,j,k)=(Uin(j,k)+Uin(j+1,k)+Uin(j,k+1)+Uin(j+1,k+1))/4._KND
     Uco(nx+2,j,k)=(Uin(j,k)+Uin(j+1,k)+Uin(j,k+1)+Uin(j+1,k+1))/4._KND
     Uco(nx+3,j,k)=(Uin(j,k)+Uin(j+1,k)+Uin(j,k+1)+Uin(j+1,k+1))/4._KND
    enddo       
   enddo
  elseif (BtypeE==NOSLIP) then
   do k=1,nz
    do j=1,ny                       !Solid wall
     Uco(nx+1,j,k)=0
     Uco(nx+2,j,k)=-Uco(nx,j,k)
     Uco(nx+3,j,k)=-Uco(nx-1,j,k)
    enddo
   enddo
  elseif (BtypeE==NEUMANN) then   !Neumann outlet
   do k=-2,nz+3
    do j=-2,ny+3
     Uco(nx+1,j,k)=Uco(nx,j,k)
     Uco(nx+2,j,k)=Uco(nx,j,k)
     Uco(nx+3,j,k)=Uco(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==FREESLIP) then  !FREESLIP
   do k=1,nz
    do j=1,ny                       
     Uco(nx+1,j,k)=0
     Uco(nx+2,j,k)=Uco(nx,j,k)
     Uco(nx+3,j,k)=Uco(nx-1,j,k)
    enddo
   enddo
  elseif (BtypeE==PERIODIC) then  !Periodic BC
   do k=1,nz
    do j=1,ny                       
     Uco(nx+1,j,k)=Uco(1,j,k)
     Uco(nx+2,j,k)=Uco(2,j,k)
     Uco(nx+3,j,k)=Uco(3,j,k)
    enddo
   enddo      
  endif



  !!!!!!ROHY

  if (BtypeS==NOSLIP.or.BtypeB==NOSLIP) then
   do i=1,nx
    Uco(i,0,0)=0
   enddo
  elseif (BtypeS==PERIODIC.and.BtypeB==PERIODIC) then
   do i=1,nx
    Uco(i,0,0)=(Uco(i,1,1)+Uco(i,1,nz)+Uco(i,ny,1)+Uco(i,ny,nz))/4._KND
   enddo
  elseif (BtypeS==PERIODIC) then
   do i=1,nx
    Uco(i,0,0)=(Uco(i,1,0)+Uco(i,ny,0))/2._KND
   enddo
  elseif (BtypeB==PERIODIC) then
   do i=1,nx
    Uco(i,0,0)=(Uco(i,0,1)+Uco(i,0,nz))/2._KND
   enddo
  else
   do i=1,nx
    Uco(i,0,0)=(Uco(i,1,0)+Uco(i,0,1))/2._KND
   enddo
  endif

  if (BtypeN==NOSLIP.or.BtypeB==NOSLIP) then
   do i=1,nx
    Uco(i,ny+1,0)=0
   enddo
  elseif (BtypeN==PERIODIC.and.BtypeB==PERIODIC) then
   do i=1,nx
    Uco(i,ny+1,0)=(Uco(i,1,1)+Uco(i,1,nz)+Uco(i,ny,1)+Uco(i,ny,nz))/4._KND
   enddo
  elseif (BtypeN==PERIODIC) then
   do i=1,nx
    Uco(i,ny+1,0)=(Uco(i,1,0)+Uco(i,ny,0))/2._KND
   enddo
  elseif (BtypeB==PERIODIC) then
   do i=1,nx
    Uco(i,ny+1,0)=(Uco(i,ny+1,1)+Uco(i,ny+1,nz))/2._KND
   enddo
  else
   do i=1,nx
    Uco(i,ny+1,0)=(Uco(i,ny,0)+Uco(i,ny+1,1))/2._KND
   enddo
  endif

  if (BtypeS==NOSLIP.or.BtypeT==NOSLIP) then
   do i=1,nx
    Uco(i,0,nz+1)=0
   enddo
  elseif (BtypeS==PERIODIC.and.BtypeT==PERIODIC) then
   do i=1,nx
    Uco(i,0,nz+1)=(Uco(i,1,1)+Uco(i,1,nz)+Uco(i,ny,1)+Uco(i,ny,nz))/4._KND
   enddo
  elseif (BtypeS==PERIODIC) then
   do i=1,nx
    Uco(i,0,nz+1)=(Uco(i,1,nz+1)+Uco(i,ny,nz+1))/2._KND
   enddo
  elseif (BtypeT==PERIODIC) then
   do i=1,nx
    Uco(i,0,nz+1)=(Uco(i,0,1)+Uco(i,0,nz))/2._KND
   enddo
  else
   do i=1,nx
    Uco(i,0,nz+1)=(Uco(i,1,nz+1)+Uco(i,0,nz))/2._KND
   enddo
  endif

  if (BtypeN==NOSLIP.or.BtypeT==NOSLIP) then
   do i=1,nx
    Uco(i,ny+1,nz+1)=0
   enddo
  elseif (BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do i=1,nx
    Uco(i,ny+1,nz+1)=(Uco(i,1,1)+Uco(i,1,nz)+Uco(i,ny,1)+Uco(i,ny,nz))/4._KND
   enddo
  elseif (BtypeN==PERIODIC) then
   do i=1,nx
    Uco(i,ny+1,nz+1)=(Uco(i,1,nz+1)+Uco(i,ny,nz+1))/2._KND
   enddo
  elseif (BtypeT==PERIODIC) then
   do i=1,nx
    Uco(i,ny+1,nz+1)=(Uco(i,ny+1,1)+Uco(i,ny+1,nz))/2._KND
   enddo
  else
   do i=1,nx
    Uco(i,ny+1,nz+1)=(Uco(i,ny,nz+1)+Uco(i,ny+1,nz))/2._KND
   enddo
  endif


  if (BtypeE==NOSLIP.or.BtypeE==FREESLIP) then
   do k=1,nz
    Uco(nx+1,0,k)=0
    Uco(nx+1,ny+1,k)=0
   enddo
   do j=1,ny
    Uco(nx+1,j,0)=0
    Uco(nx+1,j,nz+1)=0
   enddo
  elseif (BtypeE==DIRICHLET) then
   do k=1,nz
    Uco(nx+1,0,k)=(Uin(0,k)+Uin(1,k)+Uin(0,k+1)+Uin(1,k+1))/4._KND
    Uco(nx+1,ny+1,k)=(Uin(ny+1,k)+Uin(ny+2,k)+Uin(ny+1,k+1)+Uin(ny+2,k+1))/4._KND
   enddo
   do j=1,ny
    Uco(nx+1,j,0)=(Uin(j,0)+Uin(j+1,0)+Uin(j,1)+Uin(j+1,1))/4._KND
    Uco(nx+1,j,nz+1)=(Uin(j,nz+1)+Uin(j+1,nz+1)+Uin(j,nz+2)+Uin(j+1,nz+2))/4._KND
   enddo
  elseif (BtypeE==PERIODIC) then
   do k=1,nz
    Uco(nx+1,0,k)=(Uco(1,0,k)+Uco(nx,0,k))/2._KND
    Uco(nx+1,ny+1,k)=(Uco(1,ny+1,k)+Uco(nx,ny+1,k))/2._KND
   enddo
   do j=1,ny
    Uco(nx+1,j,0)=(Uco(1,j,0)+Uco(nx,j,1))/2._KND
    Uco(nx+1,j,nz+1)=(Uco(1,j,nz+1)+Uco(nx,j,nz+1))/2._KND
   enddo
  else
   do k=1,nz
    Uco(nx+1,0,k)=(Uco(nx+1,1,k)+Uco(nx,0,k))/2._KND
    Uco(nx+1,ny+1,k)=(Uco(nx+1,ny,k)+Uco(nx,ny+1,k))/2._KND
   enddo
   do j=1,ny
    Uco(nx+1,j,0)=(Uco(nx+1,j,1)+Uco(nx,j,0))/2._KND
    Uco(nx+1,j,nz+1)=(Uco(nx+1,j,nz)+Uco(nx,j,nz+1))/2._KND
   enddo
  endif


  if (BtypeW==NOSLIP.or.BtypeW==FREESLIP) then
   do k=1,nz
    Uco(0,0,k)=0
    Uco(0,ny+1,k)=0
   enddo
   do j=1,ny
    Uco(0,j,0)=0
    Uco(0,j,nz+1)=0
   enddo
  elseif (BtypeW==DIRICHLET) then
   do k=1,nz
    Uco(0,0,k)=(Uin(0,k)+Uin(1,k)+Uin(0,k+1)+Uin(1,k+1))/4._KND
    Uco(0,ny+1,k)=(Uin(ny+1,k)+Uin(ny+2,k)+Uin(ny+1,k+1)+Uin(ny+2,k+1))/4._KND
   enddo
   do j=1,ny
    Uco(0,j,0)=(Uin(j,0)+Uin(j+1,0)+Uin(j,1)+Uin(j+1,1))/4._KND
    Uco(0,j,nz+1)=(Uin(j,nz+1)+Uin(j+1,nz+1)+Uin(j,nz+2)+Uin(j+1,nz+2))/4._KND
   enddo
  elseif (BtypeW==PERIODIC) then
   do k=1,nz
    Uco(0,0,k)=(Uco(1,0,k)+Uco(nx,0,k))/2._KND
    Uco(0,ny+1,k)=(Uco(1,ny+1,k)+Uco(nx,ny+1,k))/2._KND
   enddo
   do j=1,ny
    Uco(0,j,0)=(Uco(1,j,0)+Uco(nx,j,0))/2._KND
    Uco(0,j,nz+1)=(Uco(1,j,nz+1)+Uco(nx,j,nz+1))/2._KND
   enddo
  else
   do k=1,nz
    Uco(0,0,k)=(Uco(0,1,k)+Uco(1,0,k))/2._KND
    Uco(0,ny+1,k)=(Uco(0,ny,k)+Uco(1,ny+1,k))/2._KND
   enddo
   do j=1,ny
    Uco(0,j,0)=(Uco(0,j,1)+Uco(1,j,0))/2._KND
    Uco(0,j,nz+1)=(Uco(0,j,nz)+Uco(1,j,nz+1))/2._KND
   enddo
  endif
  
  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   Uco(0,0,0)=(Uco(1,0,0)+Uco(0,1,0)+Uco(0,0,1)+Uco(nx,0,0)+Uco(0,ny,0)+Uco(0,0,nz))/6._KND
   Uco(nx,0,0)=(Uco(1,0,0)+Uco(0,1,0)+Uco(0,0,1)+Uco(nx,0,0)+Uco(0,ny,0)+Uco(0,0,nz))/6._KND
   Uco(0,ny,0)=(Uco(1,0,0)+Uco(0,1,0)+Uco(0,0,1)+Uco(nx,0,0)+Uco(0,ny,0)+Uco(0,0,nz))/6._KND
   Uco(0,0,nz)=(Uco(1,0,0)+Uco(0,1,0)+Uco(0,0,1)+Uco(nx,0,0)+Uco(0,ny,0)+Uco(0,0,nz))/6._KND
   Uco(nx,ny,0)=(Uco(1,0,0)+Uco(0,1,0)+Uco(0,0,1)+Uco(nx,0,0)+Uco(0,ny,0)+Uco(0,0,nz))/6._KND
   Uco(0,ny,nz)=(Uco(1,0,0)+Uco(0,1,0)+Uco(0,0,1)+Uco(nx,0,0)+Uco(0,ny,0)+Uco(0,0,nz))/6._KND
   Uco(nx,0,nz)=(Uco(1,0,0)+Uco(0,1,0)+Uco(0,0,1)+Uco(nx,0,0)+Uco(0,ny,0)+Uco(0,0,nz))/6._KND
   Uco(nx,ny,nz)=(Uco(1,0,0)+Uco(0,1,0)+Uco(0,0,1)+Uco(nx,0,0)+Uco(0,ny,0)+Uco(0,0,nz))/6._KND
  endif


  endsubroutine Bound_CondUco

  
  subroutine Bound_CondVco(Vco)
  real(KND),dimension(-2:,-2:,-2:):: Vco
  integer i,j,k,nx,ny,nz

  nx=Unx
  ny=Vny
  nz=Wnz   

  if (BtypeW==DIRICHLET) then
   do k=1,nz
    do j=1,ny                       !Dirichlet inlet
     Vco(0,j,k)=0
     Vco(-1,j,k)=0
     Vco(-2,j,k)=0
    enddo       
   enddo
  elseif (BtypeW==NOSLIP) then
   do k=1,nz
    do j=1,ny                       !Solid wall
     Vco(0,j,k)=0
     Vco(-1,j,k)=-Vco(1,j,k)
     Vco(-2,j,k)=-Vco(2,j,k)
    enddo
   enddo
  elseif (BtypeW==NEUMANN) then
   do k=1,nz
    do j=1,ny                       !Neumann inlet
     Vco(0,j,k)=Vco(1,j,k)
     Vco(-1,j,k)=Vco(1,j,k)
     Vco(-2,j,k)=Vco(1,j,k)
    enddo
   enddo
  elseif (BtypeW==FREESLIP) then  !FREESLIP
   do k=1,nz
    do j=1,ny                       
     Vco(0,j,k)=Vco(1,j,k)
     Vco(-1,j,k)=Vco(1,j,k)
     Vco(-2,j,k)=Vco(1,j,k)
    enddo
   enddo
  elseif (BtypeW==PERIODIC) then  !Periodic BC
   do k=1,nz
    do j=1,ny                       
     Vco(0,j,k)=Vco(nx,j,k)
     Vco(-1,j,k)=Vco(nx-1,j,k)
     Vco(-2,j,k)=Vco(nx-2,j,k)
    enddo
   enddo      
  endif      
     
  if (BtypeE==DIRICHLET) then
   do k=1,nz
    do j=1,ny                       !Dirichlet inlet
     Vco(nx+1,j,k)=0
     Vco(nx+2,j,k)=0
     Vco(nx+3,j,k)=0
    enddo       
   enddo
  elseif (BtypeE==NOSLIP) then
   do k=1,nz
    do j=1,ny                       !Solid wall
     Vco(nx+1,j,k)=0
     Vco(nx+2,j,k)=-Vco(nx,j,k)
     Vco(nx+3,j,k)=-Vco(nx-1,j,k)
    enddo
   enddo
  elseif (BtypeE==NEUMANN) then
   do k=1,nz
    do j=1,ny                       !Neumann inlet
     Vco(nx+1,j,k)=Vco(nx,j,k)
     Vco(nx+2,j,k)=Vco(nx,j,k)
     Vco(nx+3,j,k)=Vco(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==FREESLIP) then  !FREESLIP
   do k=1,nz
    do j=1,ny                       
     Vco(nx+1,j,k)=Vco(nx,j,k)
     Vco(nx+2,j,k)=Vco(nx,j,k)
     Vco(nx+3,j,k)=Vco(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==PERIODIC) then  !Periodic BC
   do k=1,nz
    do j=1,ny                       
     Vco(nx+1,j,k)=Vco(1,j,k)
     Vco(nx+2,j,k)=Vco(2,j,k)
     Vco(nx+3,j,k)=Vco(3,j,k)
    enddo
   enddo      
  endif
      
  if (BtypeS==DIRICHLET) then
   do k=1,nz
    do i=1,nx                       !Dirichlet inlet
     Vco(i,0,k)=SsideV
     Vco(i,-1,k)=SsideV+(SsideV-Vco(i,1,k))
     Vco(i,-2,k)=SsideV+(SsideV-Vco(i,2,k))
    enddo       
   enddo   
  elseif (BtypeS==NOSLIP) then
   do k=1,nz
    do i=1,nx                       !Solid wall
     Vco(i,0,k)=0
     Vco(i,-1,k)=-Vco(i,1,k)
     Vco(i,-2,k)=-Vco(i,2,k)
    enddo
   enddo
  elseif (BtypeS==NEUMANN) then
   do k=1,nz
    do i=1,nx                       !Neumann inlet
     Vco(i,0,k)=Vco(i,1,k)
     Vco(i,-1,k)=Vco(i,1,k)
     Vco(i,-2,k)=Vco(i,1,k)
    enddo
   enddo
  elseif (BtypeS==FREESLIP) then  !FREESLIP
   do k=1,nz
    do i=1,nx
     Vco(i,0,k)=0
     Vco(i,-1,k)=-Vco(i,1,k)
     Vco(i,-2,k)=-Vco(i,2,k)
    enddo
   enddo
  elseif (BtypeS==PERIODIC) then  !Periodic BC
   do k=1,nz
    do i=1,nx
     Vco(i,0,k)=Vco(i,ny,k)
     Vco(i,-1,k)=Vco(i,ny-1,k)
     Vco(i,-2,k)=Vco(i,ny-2,k)
    enddo
   enddo      
  endif   

  if (BtypeN==DIRICHLET) then
   do k=1,nz
    do i=1,nx                       !Dirichlet inlet
     Vco(i,ny+1,k)=NsideV
     Vco(i,ny+2,k)=NsideV+(NsideV-Vco(i,ny,k))
     Vco(i,ny+3,k)=NsideV+(NsideV-Vco(i,ny-1,k))
    enddo       
   enddo   
  elseif (BtypeN==NOSLIP) then
   do k=1,nz
    do i=1,nx                       !Solid wall
     Vco(i,ny+1,k)=0
     Vco(i,ny+2,k)=-Vco(i,ny,k)
     Vco(i,ny+3,k)=-Vco(i,ny-1,k)
    enddo
   enddo
  elseif (BtypeN==NEUMANN) then
   do k=1,nz
    do i=1,nx                       !Neumann inlet
     Vco(i,ny+1,k)=Vco(i,ny,k)
     Vco(i,ny+2,k)=Vco(i,ny,k)
     Vco(i,ny+3,k)=Vco(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==FREESLIP) then  !FREESLIP
   do k=1,nz
    do i=1,nx
     Vco(i,ny+1,k)=0
     Vco(i,ny+2,k)=-Vco(i,ny,k)
     Vco(i,ny+3,k)=-Vco(i,ny-1,k)
    enddo
   enddo
  elseif (BtypeN==PERIODIC) then  !Periodic BC
   do k=1,nz
    do i=1,nx
     Vco(i,ny+1,k)=Vco(i,1,k)
     Vco(i,ny+2,k)=Vco(i,2,k)
     Vco(i,ny+3,k)=Vco(i,3,k)
    enddo
   enddo      
  endif

  if (BtypeB==DIRICHLET) then
   do j=1,ny
    do i=1,nx                       !Dirichlet inlet
     Vco(i,j,0)=BsideV
     Vco(i,j,-1)=BsideV+(BsideV-Vco(i,j,1))
     Vco(i,j,-2)=BsideV+(BsideV-Vco(i,j,2))
    enddo       
   enddo   
  elseif (BtypeB==NOSLIP) then
   do j=1,ny
    do i=1,nx                       !Solid wall
     Vco(i,j,0)=0
     Vco(i,j,-1)=-Vco(i,j,1)
     Vco(i,j,-2)=-Vco(i,j,2)
    enddo
   enddo
  elseif (BtypeB==NEUMANN) then
   do j=1,ny
    do i=1,nx                       !Neumann inlet
     Vco(i,j,0)=Vco(i,j,1)
     Vco(i,j,-1)=Vco(i,j,1)
     Vco(i,j,-2)=Vco(i,j,1)
    enddo
   enddo
  elseif (BtypeB==FREESLIP) then  !FREESLIP
   do j=1,ny
    do i=1,nx
     Vco(i,j,0)=Vco(i,j,1)
     Vco(i,j,-1)=Vco(i,j,1)
     Vco(i,j,-2)=Vco(i,j,1)
    enddo
   enddo
  elseif (BtypeB==PERIODIC) then  !Periodic BC
   do j=1,ny
    do i=1,nx
     Vco(i,j,0)=Vco(i,j,nz)
     Vco(i,j,-1)=Vco(i,j,nz-1)
     Vco(i,j,-2)=Vco(i,j,nz-2)
    enddo
   enddo      
  endif   
  
  if (BtypeT==DIRICHLET) then
   do j=1,ny
    do i=1,nx                       !Dirichlet inlet
     Vco(i,j,nz+1)=TsideV
     Vco(i,j,nz+2)=TsideV+(TsideV-Vco(i,j,nz))
     Vco(i,j,nz+3)=TsideV+(TsideV-Vco(i,j,nz-1))
    enddo       
   enddo   
  elseif (BtypeT==NOSLIP) then
   do j=1,ny
    do i=1,nx                       !Solid wall
     Vco(i,j,nz+1)=0
     Vco(i,j,nz+2)=-Vco(i,j,nz)
     Vco(i,j,nz+3)=-Vco(i,j,nz-1)
    enddo
   enddo
  elseif (BtypeT==NEUMANN) then
   do j=1,ny
    do i=1,nx                       !Neumann inlet
     Vco(i,j,nz+1)=Vco(i,j,nz)
     Vco(i,j,nz+2)=Vco(i,j,nz)
     Vco(i,j,nz+3)=Vco(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==FREESLIP) then  !FREESLIP
   do j=1,ny
    do i=1,nx
     Vco(i,j,nz+1)=Vco(i,j,nz)
     Vco(i,j,nz+2)=Vco(i,j,nz)
     Vco(i,j,nz+3)=Vco(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==PERIODIC) then  !Periodic BC
   do j=1,ny
    do i=1,nx
     Vco(i,j,nz+1)=Vco(i,j,1)
     Vco(i,j,nz+2)=Vco(i,j,2)
     Vco(i,j,nz+3)=Vco(i,j,3)
    enddo
   enddo      
  endif

  !!!!!!ROHY

  if (BtypeW==NOSLIP.or.BtypeB==NOSLIP) then
   do j=1,ny
    Vco(0,j,0)=0
   enddo
  elseif (BtypeW==PERIODIC.and.BtypeB==PERIODIC) then
   do j=1,ny
    Vco(0,j,0)=(Vco(1,j,1)+Vco(1,j,nz)+Vco(nx,j,1)+Vco(nx,j,nz))/4._KND
   enddo
  elseif (BtypeW==PERIODIC) then
   do j=1,ny
    Vco(0,j,0)=(Vco(1,j,0)+Vco(nx,j,0))/2._KND
   enddo
  elseif (BtypeB==PERIODIC) then
   do j=1,ny
    Vco(0,j,0)=(Vco(0,j,1)+Vco(0,j,nz))/2._KND
   enddo
  else
   do j=1,ny
    Vco(0,j,0)=(Vco(1,j,0)+Vco(0,j,1))/2._KND
   enddo
  endif

  if (BtypeE==NOSLIP.or.BtypeB==NOSLIP) then
   do j=1,ny
    Vco(nx+1,j,0)=0
   enddo
  elseif (BtypeE==PERIODIC.and.BtypeB==PERIODIC) then
   do j=1,ny
    Vco(nx+1,j,0)=(Vco(1,j,1)+Vco(1,j,nz)+Vco(nx,j,1)+Vco(nx,j,nz))/4._KND
   enddo
  elseif (BtypeE==PERIODIC) then
   do j=1,ny
    Vco(nx+1,j,0)=(Vco(1,j,0)+Vco(nx,j,0))/2._KND
   enddo
  elseif (BtypeB==PERIODIC) then
   do j=1,ny
    Vco(nx+1,j,0)=(Vco(nx+1,j,1)+Vco(nx+1,j,nz))/2._KND
   enddo
  else
   do j=1,ny
    Vco(nx+1,j,0)=(Vco(nx,j,0)+Vco(nx+1,j,1))/2._KND
   enddo
  endif

  if (BtypeW==NOSLIP.or.BtypeT==NOSLIP) then
   do j=1,ny
    Vco(0,j,nz+1)=0
   enddo
  elseif (BtypeW==PERIODIC.and.BtypeT==PERIODIC) then
   do j=1,ny
    Vco(0,j,nz+1)=(Vco(1,j,1)+Vco(1,j,nz)+Vco(nx,j,1)+Vco(nx,j,nz))/4._KND
   enddo
  elseif (BtypeW==PERIODIC) then
   do j=1,ny
    Vco(0,j,nz+1)=(Vco(1,j,nz+1)+Vco(nx,j,nz+1))/2._KND
   enddo
  elseif (BtypeT==PERIODIC) then
   do j=1,ny
    Vco(0,j,nz+1)=(Vco(0,j,1)+Vco(0,j,nz))/2._KND
   enddo
  else
   do j=1,ny
    Vco(0,j,nz+1)=(Vco(1,j,nz+1)+Vco(0,j,nz))/2._KND
   enddo
  endif

  if (BtypeE==NOSLIP.or.BtypeT==NOSLIP) then
   do j=1,ny
    Vco(nx+1,j,nz+1)=0
   enddo
  elseif (BtypeE==PERIODIC.and.BtypeT==PERIODIC) then
   do j=1,ny
    Vco(nx+1,j,nz+1)=(Vco(1,j,1)+Vco(1,j,nz)+Vco(nx,j,1)+Vco(nx,j,nz))/4._KND
   enddo
  elseif (BtypeE==PERIODIC) then
   do j=1,ny
    Vco(nx+1,j,nz+1)=(Vco(1,j,nz+1)+Vco(nx,j,nz+1))/2._KND
   enddo
  elseif (BtypeT==PERIODIC) then
   do j=1,ny
    Vco(nx+1,j,nz+1)=(Vco(nx+1,j,1)+Vco(nx+1,j,nz))/2._KND
   enddo
  else
   do j=1,ny
    Vco(nx+1,j,nz+1)=(Vco(nx,j,nz+1)+Vco(nx+1,j,nz))/2._KND
   enddo
  endif


  if (BtypeN==NOSLIP.or.BtypeN==FREESLIP) then
   do k=1,nz
    Vco(0,ny+1,k)=0
    Vco(nx+1,ny+1,k)=0
   enddo
   do i=1,nx
    Vco(i,ny+1,0)=0
    Vco(i,ny+1,nz+1)=0
   enddo
  elseif (BtypeN==DIRICHLET) then
   do k=1,nz
    Vco(0,ny+1,k)=NsideV
    Vco(nx+1,ny+1,k)=NsideV
   enddo
   do i=1,nx
    Vco(i,ny+1,0)=NsideV
    Vco(i,ny+1,nz+1)=NsideV
   enddo
  elseif (BtypeN==PERIODIC) then
   do k=1,nz
    Vco(0,ny+1,k)=(Vco(0,1,k)+Vco(0,ny,k))/2._KND
    Vco(nx+1,ny+1,k)=(Vco(nx+1,1,k)+Vco(nx+1,ny,k))/2._KND
   enddo
   do i=1,nx
    Vco(i,ny+1,0)=(Vco(i,1,0)+Vco(i,ny,1))/2._KND
    Vco(i,ny+1,nz+1)=(Vco(i,1,nz+1)+Vco(i,ny,nz+1))/2._KND
   enddo
  else
   do k=1,nz
    Vco(0,ny+1,k)=(Vco(1,ny+1,k)+Vco(0,ny,k))/2._KND
    Vco(nx+1,ny+1,k)=(Vco(nx,ny+1,k)+Vco(nx+1,ny,k))/2._KND
   enddo
   do i=1,nx
    Vco(i,ny+1,0)=(Vco(i,ny+1,1)+Vco(i,ny,0))/2._KND
    Vco(i,ny+1,nz+1)=(Vco(i,ny+1,nz)+Vco(i,ny,nz+1))/2._KND
   enddo
  endif


  if (BtypeS==NOSLIP.or.BtypeS==FREESLIP) then
   do k=1,nz
    Vco(0,0,k)=0
    Vco(nx+1,0,k)=0
   enddo
   do i=1,nx
    Vco(i,0,0)=0
    Vco(i,0,nz+1)=0
   enddo
  elseif (BtypeS==DIRICHLET) then
   do k=1,nz
    Vco(0,0,k)=SsideV
    Vco(nx+1,0,k)=SsideV
   enddo
   do i=1,nx
    Vco(i,0,0)=SsideV
    Vco(i,0,nz+1)=SsideV
   enddo
  elseif (BtypeS==PERIODIC) then
   do k=1,nz
    Vco(0,0,k)=(Vco(0,1,k)+Vco(0,ny,k))/2._KND
    Vco(nx+1,0,k)=(Vco(nx+1,1,k)+Vco(nx+1,ny,k))/2._KND
   enddo
   do i=1,nx
    Vco(i,0,0)=(Vco(i,1,0)+Vco(i,ny,0))/2._KND
    Vco(i,0,nz+1)=(Vco(i,1,nz+1)+Vco(i,ny,nz+1))/2._KND
   enddo
  else
   do k=1,nz
    Vco(0,0,k)=(Vco(1,0,k)+Vco(0,1,k))/2._KND
    Vco(nx+1,0,k)=(Vco(nx,0,k)+Vco(nx+1,1,k))/2._KND
   enddo
   do i=1,nx
    Vco(i,0,0)=(Vco(i,0,1)+Vco(i,1,0))/2._KND
    Vco(i,0,nz+1)=(Vco(i,0,nz)+Vco(i,1,nz+1))/2._KND
   enddo
  endif

  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   Vco(0,0,0)=(Vco(1,0,0)+Vco(0,1,0)+Vco(0,0,1)+Vco(nx,0,0)+Vco(0,ny,0)+Vco(0,0,nz))/6._KND
   Vco(nx,0,0)=(Vco(1,0,0)+Vco(0,1,0)+Vco(0,0,1)+Vco(nx,0,0)+Vco(0,ny,0)+Vco(0,0,nz))/6._KND
   Vco(0,ny,0)=(Vco(1,0,0)+Vco(0,1,0)+Vco(0,0,1)+Vco(nx,0,0)+Vco(0,ny,0)+Vco(0,0,nz))/6._KND
   Vco(0,0,nz)=(Vco(1,0,0)+Vco(0,1,0)+Vco(0,0,1)+Vco(nx,0,0)+Vco(0,ny,0)+Vco(0,0,nz))/6._KND
   Vco(nx,ny,0)=(Vco(1,0,0)+Vco(0,1,0)+Vco(0,0,1)+Vco(nx,0,0)+Vco(0,ny,0)+Vco(0,0,nz))/6._KND
   Vco(0,ny,nz)=(Vco(1,0,0)+Vco(0,1,0)+Vco(0,0,1)+Vco(nx,0,0)+Vco(0,ny,0)+Vco(0,0,nz))/6._KND
   Vco(nx,0,nz)=(Vco(1,0,0)+Vco(0,1,0)+Vco(0,0,1)+Vco(nx,0,0)+Vco(0,ny,0)+Vco(0,0,nz))/6._KND
   Vco(nx,ny,nz)=(Vco(1,0,0)+Vco(0,1,0)+Vco(0,0,1)+Vco(nx,0,0)+Vco(0,ny,0)+Vco(0,0,nz))/6._KND
  endif

  endsubroutine Bound_CondVco


  subroutine Bound_CondWco(Wco)
  real(KND),dimension(-2:,-2:,-2:):: Wco
  integer i,j,k,nx,ny,nz

  nx=Unx
  ny=Vny
  nz=Wnz   

  if (BtypeW==DIRICHLET) then
   do k=1,nz
    do j=1,ny                       !Dirichlet inlet
     Wco(0,j,k)=0
     Wco(-1,j,k)=0
     Wco(-2,j,k)=0
    enddo       
   enddo
  elseif (BtypeW==NOSLIP) then
   do k=1,nz
    do j=1,ny                       !Solid wall
     Wco(0,j,k)=0
     Wco(-1,j,k)=-Wco(1,j,k)
     Wco(-2,j,k)=-Wco(2,j,k)
    enddo
   enddo
  elseif (BtypeW==NEUMANN) then
   do k=1,nz
    do j=1,ny                       !Neumann inlet
     Wco(0,j,k)=Wco(1,j,k)
     Wco(-1,j,k)=Wco(1,j,k)
     Wco(-2,j,k)=Wco(1,j,k)
    enddo
   enddo
  elseif (BtypeW==FREESLIP) then  !FREESLIP
   do k=1,nz
    do j=1,ny                       
     Wco(0,j,k)=Wco(1,j,k)
     Wco(-1,j,k)=Wco(1,j,k)
     Wco(-2,j,k)=Wco(1,j,k)
    enddo
   enddo
  elseif (BtypeW==PERIODIC) then  !Periodic BC
   do k=1,nz
    do j=1,ny                       
     Wco(0,j,k)=Wco(nx,j,k)
     Wco(-1,j,k)=Wco(nx-1,j,k)
     Wco(-2,j,k)=Wco(nx-2,j,k)
    enddo
   enddo      
  endif      
     
  if (BtypeE==DIRICHLET) then
   do k=1,nz
    do j=1,ny                       !Dirichlet inlet
     Wco(nx+1,j,k)=0
     Wco(nx+2,j,k)=0
     Wco(nx+3,j,k)=0
    enddo       
   enddo
  elseif (BtypeE==NOSLIP) then
   do k=1,nz
    do j=1,ny                       !Solid wall
     Wco(nx+1,j,k)=0
     Wco(nx+2,j,k)=-Wco(nx,j,k)
     Wco(nx+3,j,k)=-Wco(nx-1,j,k)
    enddo
   enddo
  elseif (BtypeE==NEUMANN) then
   do k=1,nz
    do j=1,ny                       !Neumann inlet
     Wco(nx+1,j,k)=Wco(nx,j,k)
     Wco(nx+2,j,k)=Wco(nx,j,k)
     Wco(nx+3,j,k)=Wco(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==FREESLIP) then  !FREESLIP
   do k=1,nz
    do j=1,ny                       
     Wco(nx+1,j,k)=Wco(nx,j,k)
     Wco(nx+2,j,k)=Wco(nx,j,k)
     Wco(nx+3,j,k)=Wco(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==PERIODIC) then  !Periodic BC
   do k=1,nz
    do j=1,ny                       
     Wco(nx+1,j,k)=Wco(1,j,k)
     Wco(nx+2,j,k)=Wco(2,j,k)
     Wco(nx+3,j,k)=Wco(3,j,k)
    enddo
   enddo      
  endif
      
  if (BtypeS==DIRICHLET) then
   do k=1,nz
    do i=1,nx                       !Dirichlet inlet
     Wco(i,0,k)=SsideW
     Wco(i,-1,k)=SsideW+(SsideW-Wco(i,1,k))
     Wco(i,-2,k)=SsideW+(SsideW-Wco(i,2,k))
    enddo       
   enddo   
  elseif (BtypeS==NOSLIP) then
   do k=1,nz
    do i=1,nx                       !Solid wall
     Wco(i,0,k)=0
     Wco(i,-1,k)=-Wco(i,1,k)
     Wco(i,-2,k)=-Wco(i,2,k)
    enddo
   enddo
  elseif (BtypeS==NEUMANN) then
   do k=1,nz
    do i=1,nx                       !Neumann inlet
     Wco(i,0,k)=Wco(i,1,k)
     Wco(i,-1,k)=Wco(i,1,k)
     Wco(i,-2,k)=Wco(i,1,k)
    enddo
   enddo
  elseif (BtypeS==FREESLIP) then  !FREESLIP
   do k=1,nz
    do i=1,nx
     Wco(i,0,k)=Wco(i,1,k)
     Wco(i,-1,k)=Wco(i,1,k)
     Wco(i,-2,k)=Wco(i,1,k)
    enddo
   enddo
  elseif (BtypeS==PERIODIC) then  !Periodic BC
   do k=1,nz
    do i=1,nx
     Wco(i,0,k)=Wco(i,ny,k)
     Wco(i,-1,k)=Wco(i,ny-1,k)
     Wco(i,-2,k)=Wco(i,ny-2,k)
    enddo
   enddo      
  endif   

  if (BtypeN==DIRICHLET) then
   do k=1,nz
    do i=1,nx                       !Dirichlet inlet
     Wco(i,ny+1,k)=NsideW
     Wco(i,ny+2,k)=NsideW+(NsideW-Wco(i,ny,k))
     Wco(i,ny+3,k)=NsideW+(NsideW-Wco(i,ny-1,k))
    enddo       
   enddo   
  elseif (BtypeN==NOSLIP) then
   do k=1,nz
    do i=1,nx                       !Solid wall
     Wco(i,ny+1,k)=0
     Wco(i,ny+2,k)=-Wco(i,ny,k)
     Wco(i,ny+3,k)=-Wco(i,ny-1,k)
    enddo
   enddo
  elseif (BtypeN==NEUMANN) then
   do k=1,nz
    do i=1,nx                       !Neumann inlet
     Wco(i,ny+1,k)=Wco(i,ny,k)
     Wco(i,ny+2,k)=Wco(i,ny,k)
     Wco(i,ny+3,k)=Wco(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==FREESLIP) then  !FREESLIP
   do k=1,nz
    do i=1,nx
     Wco(i,ny+1,k)=Wco(i,ny,k)
     Wco(i,ny+2,k)=Wco(i,ny,k)
     Wco(i,ny+3,k)=Wco(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==PERIODIC) then  !Periodic BC
   do k=1,nz
    do i=1,nx
     Wco(i,ny+1,k)=Wco(i,1,k)
     Wco(i,ny+2,k)=Wco(i,2,k)
     Wco(i,ny+3,k)=Wco(i,3,k)
    enddo
   enddo      
  endif

  if (BtypeB==DIRICHLET) then
   do j=1,ny
    do i=1,nx                       !Dirichlet inlet
     Wco(i,j,0)=BsideW
     Wco(i,j,-1)=BsideW+(BsideW-Wco(i,j,1))
     Wco(i,j,-2)=BsideW+(BsideW-Wco(i,j,2))
    enddo       
   enddo   
  elseif (BtypeB==NOSLIP) then
   do j=1,ny
    do i=1,nx                       !Solid wall
     Wco(i,j,0)=0
     Wco(i,j,-1)=-Wco(i,j,1)
     Wco(i,j,-2)=-Wco(i,j,2)
    enddo
   enddo
  elseif (BtypeB==NEUMANN) then
   do j=1,ny
    do i=1,nx                       !Neumann inlet
     Wco(i,j,0)=Wco(i,j,1)
     Wco(i,j,-1)=Wco(i,j,1)
     Wco(i,j,-2)=Wco(i,j,1)
    enddo
   enddo
  elseif (BtypeB==FREESLIP) then  !FREESLIP
   do j=1,ny
    do i=1,nx
     Wco(i,j,0)=0
     Wco(i,j,-1)=-Wco(i,j,1)
     Wco(i,j,-2)=-Wco(i,j,2)
    enddo
   enddo
  elseif (BtypeB==PERIODIC) then  !Periodic BC
   do j=1,ny
    do i=1,nx
     Wco(i,j,0)=Wco(i,j,nz)
     Wco(i,j,-1)=Wco(i,j,nz-1)
     Wco(i,j,-2)=Wco(i,j,nz-2)
    enddo
   enddo      
  endif   
  
  if (BtypeT==DIRICHLET) then
   do j=1,ny
    do i=1,nx                       !Dirichlet inlet
     Wco(i,j,nz+1)=TsideW
     Wco(i,j,nz+2)=TsideW+(TsideW-Wco(i,j,nz))
     Wco(i,j,nz+3)=TsideW+(TsideW-Wco(i,j,nz-1))
    enddo       
   enddo   
  elseif (BtypeT==NOSLIP) then
   do j=1,ny
    do i=1,nx                       !Solid wall
     Wco(i,j,nz+1)=0
     Wco(i,j,nz+2)=-Wco(i,j,nz)
     Wco(i,j,nz+3)=-Wco(i,j,nz-1)
    enddo
   enddo
  elseif (BtypeT==NEUMANN) then
   do j=1,ny
    do i=1,nx                       !Neumann inlet
     Wco(i,j,nz+1)=Wco(i,j,nz)
     Wco(i,j,nz+2)=Wco(i,j,nz)
     Wco(i,j,nz+3)=Wco(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==FREESLIP) then  !FREESLIP
   do j=1,ny
    do i=1,nx
     Wco(i,j,nz+1)=0
     Wco(i,j,nz+2)=-Wco(i,j,nz)
     Wco(i,j,nz+3)=-Wco(i,j,nz-1)
    enddo
   enddo
  elseif (BtypeT==PERIODIC) then  !Periodic BC
   do j=1,ny
    do i=1,nx
     Wco(i,j,nz+1)=Wco(i,j,1)
     Wco(i,j,nz+2)=Wco(i,j,2)
     Wco(i,j,nz+3)=Wco(i,j,3)
    enddo
   enddo      
  endif

  

  if (BtypeW==NOSLIP.or.BtypeS==NOSLIP) then
   do k=1,nz
    Wco(0,0,k)=0
   enddo
  elseif (BtypeW==PERIODIC.and.BtypeS==PERIODIC) then
   do k=1,nz
    Wco(0,0,k)=(Wco(1,1,k)+Wco(1,ny,k)+Wco(nx,1,k)+Wco(nx,ny,k))/4._KND
   enddo
  elseif (BtypeW==PERIODIC) then
   do k=1,nz
    Wco(0,0,k)=(Wco(1,0,k)+Wco(nx,0,k))/2._KND
   enddo
  elseif (BtypeS==PERIODIC) then
   do k=1,nz
    Wco(0,0,k)=(Wco(0,1,k)+Wco(0,ny,k))/2._KND
   enddo
  else
   do k=1,nz
    Wco(0,0,k)=(Wco(1,0,k)+Wco(0,1,k))/2._KND
   enddo
  endif

  if (BtypeE==NOSLIP.or.BtypeS==NOSLIP) then
   do k=1,nz
    Wco(nx+1,0,k)=0
   enddo
  elseif (BtypeE==PERIODIC.and.BtypeS==PERIODIC) then
   do k=1,nz
    Wco(nx+1,0,k)=(Wco(1,1,k)+Wco(1,ny,k)+Wco(nx,1,k)+Wco(nx,ny,k))/4._KND
   enddo
  elseif (BtypeE==PERIODIC) then
   do k=1,nz
    Wco(nx+1,0,k)=(Wco(1,0,k)+Wco(nx,0,k))/2._KND
   enddo
  elseif (BtypeS==PERIODIC) then
   do k=1,nz
    Wco(nx+1,0,k)=(Wco(nx+1,1,k)+Wco(nx+1,ny,k))/2._KND
   enddo
  else
   do k=1,nz
    Wco(nx+1,0,k)=(Wco(nx,0,k)+Wco(nx+1,1,k))/2._KND
   enddo
  endif

  if (BtypeW==NOSLIP.or.BtypeN==NOSLIP) then
   do k=1,nz
    Wco(0,ny+1,k)=0
   enddo
  elseif (BtypeW==PERIODIC.and.BtypeN==PERIODIC) then
   do k=1,nz
    Wco(0,ny+1,k)=(Wco(1,1,k)+Wco(1,ny,k)+Wco(nx,1,k)+Wco(nx,ny,k))/4._KND
   enddo
  elseif (BtypeW==PERIODIC) then
   do k=1,nz
    Wco(0,ny+1,k)=(Wco(1,ny+1,k)+Wco(nx,ny+1,k))/2._KND
   enddo
  elseif (BtypeN==PERIODIC) then
   do k=1,nz
    Wco(0,ny+1,k)=(Wco(0,1,k)+Wco(0,ny,k))/2._KND
   enddo
  else
   do k=1,nz
    Wco(0,ny+1,k)=(Wco(1,ny+1,k)+Wco(0,ny,k))/2._KND
   enddo
  endif

  if (BtypeE==NOSLIP.or.BtypeN==NOSLIP) then
   do k=1,nz
    Wco(nx+1,ny+1,k)=0
   enddo
  elseif (BtypeE==PERIODIC.and.BtypeN==PERIODIC) then
   do k=1,nz
    Wco(nx+1,ny+1,k)=(Wco(1,1,k)+Wco(1,ny,k)+Wco(nx,1,k)+Wco(nx,ny,k))/4._KND
   enddo
  elseif (BtypeE==PERIODIC) then
   do k=1,nz
    Wco(nx+1,ny+1,k)=(Wco(1,ny+1,k)+Wco(nx,ny+1,k))/2._KND
   enddo
  elseif (BtypeN==PERIODIC) then
   do k=1,nz
    Wco(nx+1,ny+1,k)=(Wco(nx+1,1,k)+Wco(nx+1,ny,k))/2._KND
   enddo
  else
   do k=1,nz
    Wco(nx+1,ny+1,k)=(Wco(nx,ny+1,k)+Wco(nx+1,ny,k))/2._KND
   enddo
  endif


  if (BtypeT==NOSLIP.or.BtypeT==FREESLIP) then
   do j=1,ny
    Wco(0,j,nz+1)=0
    Wco(nx+1,j,nz+1)=0
   enddo
   do i=1,nx
    Wco(i,0,nz+1)=0
    Wco(i,ny+1,nz+1)=0
   enddo
  elseif (BtypeT==DIRICHLET) then
   do j=1,ny
    Wco(0,j,nz+1)=TsideW
    Wco(nx+1,j,nz+1)=TsideW
   enddo
   do i=1,nx
    Wco(i,0,nz+1)=TsideW
    Wco(i,ny+1,nz+1)=TsideW
   enddo
  elseif (BtypeT==PERIODIC) then
   do j=1,ny
    Wco(0,j,nz+1)=(Wco(0,j,1)+Wco(0,j,nz))/2._KND
    Wco(nx+1,j,nz+1)=(Wco(nx+1,j,1)+Wco(nx+1,j,nz))/2._KND
   enddo
   do i=1,nx
    Wco(i,0,nz+1)=(Wco(i,0,1)+Wco(i,0,nz))/2._KND
    Wco(i,ny+1,nz+1)=(Wco(i,ny+1,1)+Wco(i,ny+1,nz))/2._KND
   enddo
  else
   do j=1,ny
    Wco(0,j,nz+1)=(Wco(1,j,nz+1)+Wco(0,j,nz))/2._KND
    Wco(nx+1,j,nz+1)=(Wco(nx,j,nz+1)+Wco(nx+1,j,nz))/2._KND
   enddo
   do i=1,nx
    Wco(i,0,nz+1)=(Wco(i,1,nz+1)+Wco(i,0,nz))/2._KND
    Wco(i,ny+1,nz+1)=(Wco(i,ny,nz+1)+Wco(i,ny+1,nz))/2._KND
   enddo
  endif


  if (BtypeB==NOSLIP.or.BtypeB==FREESLIP) then
   do j=1,ny
    Wco(0,j,0)=0
    Wco(nx+1,j,0)=0
   enddo
   do i=1,nx
    Wco(i,0,0)=0
    Wco(i,ny+1,0)=0
   enddo
  elseif (BtypeB==DIRICHLET) then
   do j=1,ny
    Wco(0,j,0)=BsideW
    Wco(nx+1,j,0)=BsideW
   enddo
   do i=1,nx
    Wco(i,0,0)=BsideW
    Wco(i,ny+1,0)=BsideW
   enddo
  elseif (BtypeB==PERIODIC) then
   do j=1,ny
    Wco(0,j,0)=(Wco(0,j,1)+Wco(0,j,nz))/2._KND
    Wco(nx+1,j,0)=(Wco(nx+1,j,1)+Wco(nx+1,j,nz))/2._KND
   enddo
   do i=1,nx
    Wco(i,0,0)=(Wco(i,0,1)+Wco(i,0,nz))/2._KND
    Wco(i,ny+1,0)=(Wco(i,ny+1,1)+Wco(i,ny+1,nz))/2._KND
   enddo
  else
   do j=1,ny
    Wco(0,j,0)=(Wco(1,j,0)+Wco(0,j,1))/2._KND
    Wco(nx+1,j,0)=(Wco(nx,j,0)+Wco(nx+1,j,1))/2._KND
   enddo
   do i=1,nx
    Wco(i,0,0)=(Wco(i,1,0)+Wco(i,0,1))/2._KND
    Wco(i,ny+1,0)=(Wco(i,ny,0)+Wco(i,ny+1,1))/2._KND
   enddo
  endif

  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   Wco(0,0,0)=(Wco(1,0,0)+Wco(0,1,0)+Wco(0,0,1)+Wco(nx,0,0)+Wco(0,ny,0)+Wco(0,0,nz))/6._KND
   Wco(nx,0,0)=(Wco(1,0,0)+Wco(0,1,0)+Wco(0,0,1)+Wco(nx,0,0)+Wco(0,ny,0)+Wco(0,0,nz))/6._KND
   Wco(0,ny,0)=(Wco(1,0,0)+Wco(0,1,0)+Wco(0,0,1)+Wco(nx,0,0)+Wco(0,ny,0)+Wco(0,0,nz))/6._KND
   Wco(0,0,nz)=(Wco(1,0,0)+Wco(0,1,0)+Wco(0,0,1)+Wco(nx,0,0)+Wco(0,ny,0)+Wco(0,0,nz))/6._KND
   Wco(nx,ny,0)=(Wco(1,0,0)+Wco(0,1,0)+Wco(0,0,1)+Wco(nx,0,0)+Wco(0,ny,0)+Wco(0,0,nz))/6._KND
   Wco(0,ny,nz)=(Wco(1,0,0)+Wco(0,1,0)+Wco(0,0,1)+Wco(nx,0,0)+Wco(0,ny,0)+Wco(0,0,nz))/6._KND
   Wco(nx,0,nz)=(Wco(1,0,0)+Wco(0,1,0)+Wco(0,0,1)+Wco(nx,0,0)+Wco(0,ny,0)+Wco(0,0,nz))/6._KND
   Wco(nx,ny,nz)=(Wco(1,0,0)+Wco(0,1,0)+Wco(0,0,1)+Wco(nx,0,0)+Wco(0,ny,0)+Wco(0,0,nz))/6._KND
  endif

  endsubroutine Bound_CondWco


  elemental real(KND) function LIMITER(a,b,c)
  real(KND),intent(in):: a,b,c
  real(KND) R
  real(KND),parameter:: epsil=0.00001_KND
  
  if (limitertype==minmodlim) then
   LIMITER=MINMOD(a,b)
  elseif (limitertype==extminmodlim) then
   LIMITER=MINMOD(limparam*a,limparam*b,c)
  elseif (limitertype==gammalim) then
   if (abs(b)>epsil.and.abs(a)>epsil) then
    R=a/b
    LIMITER=b*GAMMA(R)
   else
    LIMITER=MINMOD(a,b)
   endif
  elseif (limitertype==vanalbadalim) then
   if (abs(b)>epsil) then
    R=a/b
    LIMITER=b*VANALBADA(R)
   else
    LIMITER=MINMOD(a,b)
   endif
  elseif (limitertype==vanleerlim) then
   if (abs(b)>epsil) then
    R=a/b
    LIMITER=b*VANLEER(R)
   else
    LIMITER=MINMOD(a,b)
   endif
  elseif (limitertype==superbeelim) then
   if (abs(b)>epsil) then
    R=a/b
    LIMITER=b*SUPERBEE(R)
   else
    LIMITER=MINMOD(a,b)
   endif
  elseif (limitertype==0) then
   LIMITER=0
  else
   LIMITER=(a+b)/2._KND
  endif
  endfunction LIMITER

  elemental real(KND) function HEAV(x)
  real(KND),intent(in):: x
  HEAV=(1._KND+SIGN(1._KND,x))/2._KND
  endfunction HEAV

  elemental real(KND) function GAMMA(R)
  real(KND),intent(in):: R
  GAMMA=((1-limparam)/limparam)*R*(HEAV(R) - HEAV(r-(limparam/(1-limparam)))) + HEAV(r-(limparam/(1-limparam)))
  endfunction GAMMA

  elemental real(KND) function VANALBADA(R)
  real(KND),intent(in):: R
  VANALBADA=(R+R*R)/(1+R*R)
  endfunction VANALBADA

  elemental real(KND) function VANLEER(R)
  real(KND),intent(in):: R
  VANLEER=(R+abs(R))/(1+abs(R))
  endfunction VANLEER

  elemental real(KND) function SUPERBEE(R)
  real(KND),intent(in):: R
  SUPERBEE=MAX(0._KND,MAX(MIN(2*R,1._KND),MIN(R,2._KND)))
  endfunction SUPERBEE

  elemental real(KND) function MINMOD2(a,b)
  real(KND),intent(in):: a,b
      MINMOD2=(SIGN(1._KND,a)+SIGN(1._KND,b))*MIN(ABS(a),ABS(b))/2._KND
  endfunction MINMOD2

  elemental real(KND) function MINMOD3(a,b,c)
  real(KND),intent(in):: a,b,c
    if ((a>0).and.(b>0).and.(c>0)) then  
        MINMOD3=MIN(a,b,c)
    elseif ((a<0).and.(b<0).and.(c<0)) then
        MINMOD3=MAX(a,b,c)
    else
        MINMOD3=0
    endif
  endfunction MINMOD3




 subroutine SLOPESUNIF(U,V,W)
 real(KND) U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
 integer i,j,k
 real(KND) Ul,Ur,Uc,Vl,Vr,Vc,Wl,Wr,Wc

 !$OMP PARALLEL PRIVATE(Ul,Ur,Uc,Vl,Vr,Vc,Wl,Wr,Wc)
 !$OMP DO
  do k=-1,Unz+2
   do j=-1,Uny+2
    do i=-1,Unx+2
     Ur=(U(i+1,j,k)-U(i,j,k))
     Ul=(U(i,j,k)-U(i-1,j,k))
     Uc=(Ul+Ur)/2._KND
     Ux(i,j,k)=LIMITER(Ul,Ur,Uc)
    enddo
   enddo
  enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
  do k=-1,Vnz+2
   do j=-1,Vny+2
    do i=-1,Vnx+2
     Vr=(V(i,j+1,k)-V(i,j,k))
     Vl=(V(i,j,k)-V(i,j-1,k))
     Vc=(Vl+Vr)/2._KND
     Vy(i,j,k)=LIMITER(Vl,Vr,Vc)
    enddo
   enddo
  enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
  do k=-1,Wnz+2
   do j=-1,Wny+2
    do i=-1,Wnx+2
     Wr=(W(i,j,k+1)-W(i,j,k))
     Wl=(W(i,j,k)-W(i,j,k-1))
     Wc=(Wl+Wr)/2._KND
     Wz(i,j,k)=LIMITER(Wl,Wr,Wc)
    enddo
   enddo
  enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Unz+2
  do j=-1,Uny+2
   do i=-1,Unx+2
    Ur=(U(i,j+1,k)-U(i,j,k))
    Ul=(U(i,j,k)-U(i,j-1,k))
    Uc=(Ul+Ur)/2._KND
    Uy(i,j,k)=LIMITER(Ul,Ur,Uc)
   enddo
  enddo
 enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Unz+2
  do j=-1,Uny+2
   do i=-1,Unx+2
    Ur=(U(i,j,k+1)-U(i,j,k))
    Ul=(U(i,j,k)-U(i,j,k-1))
    Uc=(Ul+Ur)/2._KND
    Uz(i,j,k)=LIMITER(Ul,Ur,Uc)
   enddo
  enddo
 enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Vnz+2
  do j=-1,Vny+2
   do i=-1,Vnx+2
    Vr=(V(i+1,j,k)-V(i,j,k))
    Vl=(V(i,j,k)-V(i-1,j,k))
    Vc=(Vl+Vr)/2._KND
    Vx(i,j,k)=LIMITER(Vl,Vr,Vc)
   enddo
  enddo
 enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Vnz+2
  do j=-1,Vny+2
   do i=-1,Vnx+2
    Vr=(V(i,j,k+1)-V(i,j,k))
    Vl=(V(i,j,k)-V(i,j,k-1))
    Vc=(Vl+Vr)/2._KND
    Vz(i,j,k)=LIMITER(Vl,Vr,Vc)
   enddo
  enddo
 enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Wnz+2
  do j=-1,Wny+2
   do i=-1,Wnx+2
    Wr=(W(i+1,j,k)-W(i,j,k))
    Wl=(W(i,j,k)-W(i-1,j,k))
    Wc=(Wl+Wr)/2._KND
    Wx(i,j,k)=LIMITER(Wl,Wr,Wc)
   enddo
  enddo
 enddo
 !$OMP ENDDO NOWAIT
 !$OMP DO
 do k=-1,Wnz+2
  do j=-1,Wny+2
   do i=-1,Wnx+2
    Wr=(W(i,j+1,k)-W(i,j,k))
    Wl=(W(i,j,k)-W(i,j-1,k))
    Wc=(Wl+Wr)/2._KND
    Wy(i,j,k)=LIMITER(Wl,Wr,Wc)
   enddo
  enddo
 enddo
 !$OMP ENDDO
 !$OMP WORKSHARE
 Ux=Ux/dxmin
 Vx=Vx/dxmin
 Wx=Wx/dxmin
 Uy=Uy/dymin
 Vy=Vy/dymin
 Wy=Wy/dymin
 Uz=Uz/dzmin
 Vz=Vz/dzmin
 Wz=Wz/dzmin
 !$OMP ENDWORKSHARE
 !$OMP END PARALLEL
 end subroutine SLOPESUNIF


 


end module upwind

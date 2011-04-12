module LAXFRIED
 use PARAMETERS
 use BOUNDARIES
 implicit none

 contains
   

 subroutine LF(U2,V2,W2,U,V,W)
  real(KND),dimension(-2:,-2:,-2:):: U2,V2,W2,U,V,W

  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)

  call FLUXESU(U2,U,V,W)
  call FLUXESV(V2,U,V,W)
  call FLUXESW(W2,U,V,W)
 endsubroutine LF
  
 subroutine FLUXESU(U2,U,V,W)
 real(KND),dimension(-2:,-2:,-2:):: U2,U,V,W
 real(KND) Vlocn,Vlocs,Wloct,Wlocb
 integer i,j,k
 
  do k=1,Unz
   do j=1,Uny
    do i=1,Unx
     Vlocn=(V(i,j,k)+V(i,j+1,k)+V(i+1,j,k)+V(i+1,j+1,k))/4._KND
     Vlocs=(V(i,j-2,k)+V(i,j-1,k)+V(i+1,j-2,k)+V(i+1,j-1,k))/4._KND
     Wloct=(W(i,j,k)+W(i,j,k+1)+W(i+1,j,k)+W(i+1,j,k+1))/4._KND
     Wlocb=(W(i,j,k-2)+W(i,j,k-1)+W(i+1,j,k-2)+W(i+1,j,k-1))/4._KND
     U2(i,j,k)=(U(i+1,j,k)+U(i-1,j,k)+U(i,j+1,k)+U(i,j-1,k)+U(i,j,k-1)+U(i,j,k+1))/6._KND&
               -(dt/dxmin)*((U(i+1,j,k)*U(i+1,j,k))-(U(i-1,j,k)*U(i-1,j,k)))/2._KND&
               -(dt/dymin)*((U(i,j+1,k)*Vlocn)-(U(i,j-1,k)*Vlocs))/2._KND&
               -(dt/dzmin)*((U(i,j,k+1)*Wloct)-(U(i,j,k-1)*Wlocb))/2._KND
    enddo
   enddo
  enddo             
  U2=U2-U
 endsubroutine FLUXESU

 subroutine FLUXESV(V2,U,V,W)
 real(KND),dimension(-2:,-2:,-2:):: V2,U,V,W
 real(KND) Uloce,Ulocw,Wloct,Wlocb
 integer i,j,k

  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
     Uloce=(U(i,j,k)+U(i+1,j,k)+U(i,j+1,k)+U(i+1,j+1,k))/4._KND
     Ulocw=(U(i-2,j,k)+U(i-1,j,k)+U(i-2,j+1,k)+U(i-1,j+1,k))/4._KND
     Wloct=(W(i,j,k)+W(i,j,k+1)+W(i,j+1,k)+W(i,j+1,k+1))/4._KND
     Wlocb=(W(i,j,k-2)+W(i,j,k-1)+W(i,j+1,k-2)+W(i,j+1,k-1))/4._KND
     V2(i,j,k)=(V(i+1,j,k)+V(i-1,j,k)+V(i,j+1,k)+V(i,j-1,k)+V(i,j,k-1)+V(i,j,k+1))/6._KND&
               -(dt/dxmin)*((V(i+1,j,k)*Uloce)-(V(i-1,j,k)*Ulocw))/2._KND&
               -(dt/dymin)*((V(i,j+1,k)*V(i,j+1,k))-(V(i,j-1,k)*V(i,j-1,k)))/2._KND&
               -(dt/dzmin)*((V(i,j,k+1)*Wloct)-(V(i,j,k-1)*Wlocb))/2._KND
    enddo
   enddo
  enddo             
  V2=V2-V
 endsubroutine FLUXESV
 
 subroutine FLUXESW(W2,U,V,W)
 real(KND),dimension(-2:,-2:,-2:):: W2,U,V,W
 real(KND) Uloce,Ulocw,Vlocn,Vlocs
 integer i,j,k

  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
     Uloce=(U(i,j,k)+U(i+1,j,k)+U(i,j,k+1)+U(i+1,j,k+1))/4._KND
     Ulocw=(U(i-2,j,k)+U(i-1,j,k)+U(i-2,j,k+1)+U(i-1,j,k+1))/4._KND
     Vlocn=(V(i,j,k)+V(i,j+1,k)+V(i,j,k+1)+V(i,j+1,k+1))/4._KND
     Vlocs=(V(i,j-2,k)+V(i,j-1,k)+V(i,j-2,k+1)+V(i,j-1,k+1))/4._KND
     W2(i,j,k)=(W(i+1,j,k)+W(i-1,j,k)+W(i,j+1,k)+W(i,j-1,k)+W(i,j,k-1)+W(i,j,k+1))/6._KND&
               -(dt/dxmin)*((W(i+1,j,k)*Uloce)-(W(i-1,j,k)*Ulocw))/2._KND&
               -(dt/dymin)*((W(i,j+1,k)*Vlocn)-(W(i,j-1,k)*Vlocs))/2._KND&
               -(dt/dzmin)*((W(i,j,k+1)*W(i,j,k+1))-(W(i,j,k-1)*W(i,j,k-1)))/2._KND
    enddo
   enddo
  enddo             
  W2=W2-W
 endsubroutine FLUXESW


  !!!!Jsou ted spravne indexy hran a objemu????
 

end module LAXFRIED

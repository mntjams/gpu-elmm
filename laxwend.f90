module LAXWEND
 use PARAMETERS
 use BOUNDARIES
 implicit none

 contains


 subroutine MC1(U2,V2,W2,U,V,W)
 real(KND),dimension(-2:,-2:,-2:),intent(inout):: U2,V2,W2,U,V,W  !U,V,W changed only by BC
 real(KND) U1(LBOUND(U,1):UBOUND(U,1),LBOUND(U,2):UBOUND(U,2),LBOUND(U,3):UBOUND(U,3))
 real(KND) V1(LBOUND(V,1):UBOUND(V,1),LBOUND(V,2):UBOUND(V,2),LBOUND(V,3):UBOUND(V,3))
 real(KND) W1(LBOUND(W,1):UBOUND(W,1),LBOUND(W,2):UBOUND(W,2),LBOUND(W,3):UBOUND(W,3))
 integer i,j,k
 real(KND) Uloc,Uloce,Ulocw,Vloc,Vlocn,Vlocs,Wloc,Wloct,Wlocb

  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)

  do k=1,Unz
   do j=1,Uny
    do i=1,Unx
     Vlocn=(V(i,j,k)+V(i,j+1,k)+V(i+1,j,k)+V(i+1,j+1,k))/4._KND
     Vloc=(V(i,j-1,k)+V(i,j,k)+V(i+1,j-1,k)+V(i+1,j,k))/4._KND
     Wloct=(W(i,j,k)+W(i,j,k+1)+W(i+1,j,k)+W(i+1,j,k+1))/4._KND
     Wloc=(W(i,j,k-1)+W(i,j,k)+W(i+1,j,k-1)+W(i+1,j,k))/4._KND
     U1(i,j,k)=U(i,j,k)&
                       -(dt/dxmin)*((U(i+1,j,k)*U(i+1,j,k)-U(i,j,k)*U(i,j,k)))&
                       -(dt/dymin)*((U(i,j+1,k)*Vlocn)-(U(i,j,k)*Vloc))&
                       -(dt/dzmin)*((U(i,j,k+1)*Wloct)-(U(i,j,k)*Wloc))
    enddo
   enddo
  enddo

  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
     Uloce=(U(i,j,k)+U(i+1,j,k)+U(i,j+1,k)+U(i+1,j+1,k))/4._KND
     Uloc=(U(i-1,j,k)+U(i,j,k)+U(i-1,j+1,k)+U(i,j+1,k))/4._KND
     Wloct=(W(i,j,k)+W(i,j,k+1)+W(i,j+1,k)+W(i,j+1,k+1))/4._KND
     Wloc=(W(i,j,k-1)+W(i,j,k)+W(i,j+1,k-1)+W(i,j+1,k))/4._KND
     V1(i,j,k)=V(i,j,k)&
               -(dt/dxmin)*((V(i+1,j,k)*Uloce)-(V(i,j,k)*Uloc))&
               -(dt/dymin)*((V(i,j+1,k)*V(i,j+1,k))-(V(i,j,k)*V(i,j,k)))&
               -(dt/dzmin)*((V(i,j,k+1)*Wloct)-(V(i,j,k)*Wloc))
    enddo
   enddo
  enddo             

  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
     Uloce=(U(i,j,k)+U(i+1,j,k)+U(i,j,k+1)+U(i+1,j,k+1))/4._KND
     Uloc=(U(i-1,j,k)+U(i,j,k)+U(i-1,j,k+1)+U(i,j,k+1))/4._KND
     Vlocn=(V(i,j,k)+V(i,j+1,k)+V(i,j,k+1)+V(i,j+1,k+1))/4._KND
     Vloc=(V(i,j-1,k)+V(i,j,k)+V(i,j-1,k+1)+V(i,j,k+1))/4._KND
     W1(i,j,k)=W(i,j,k)&
               -(dt/dxmin)*((W(i+1,j,k)*Uloce)-(W(i,j,k)*Uloc))&
               -(dt/dymin)*((W(i,j+1,k)*Vlocn)-(W(i,j,k)*Vloc))&
               -(dt/dzmin)*((W(i,j,k+1)*W(i,j,k+1))-(W(i,j,k)*W(i,j,k)))
    enddo
   enddo
  enddo             
 
  call Bound_CondU(U1)
  call Bound_CondV(V1)
  call Bound_CondW(W1)

  do k=1,Unz
   do j=1,Uny
    do i=1,Unx
     Vlocs=(V1(i,j-2,k)+V1(i,j-1,k)+V1(i+1,j-2,k)+V1(i+1,j-1,k))/4._KND
     Vloc=(V1(i,j-1,k)+V1(i,j,k)+V1(i+1,j-1,k)+V1(i+1,j,k))/4._KND
     Wlocb=(W1(i,j,k-2)+W1(i,j,k-1)+W1(i+1,j,k-2)+W1(i+1,j,k-1))/4._KND
     Wloc=(W1(i,j,k-1)+W1(i,j,k)+W1(i+1,j,k-1)+W1(i+1,j,k))/4._KND
     U2(i,j,k)=(U(i,j,k)+U1(i,j,k)&
                       -(dt/dxmin)*((U1(i,j,k)*U1(i,j,k)-U1(i-1,j,k)*U1(i-1,j,k)))&
                       -(dt/dymin)*((U1(i,j,k)*Vloc)-(U1(i,j-1,k)*Vlocs))&
                       -(dt/dzmin)*((U1(i,j,k)*Wloc)-(U1(i,j,k-1)*Wlocb)))/2._KND
    enddo
   enddo
  enddo

  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
     Ulocw=(U1(i-2,j,k)+U1(i-1,j,k)+U1(i-2,j+1,k)+U1(i-1,j+1,k))/4._KND
     Uloc=(U1(i-1,j,k)+U1(i,j,k)+U1(i-1,j+1,k)+U1(i,j+1,k))/4._KND
     Wlocb=(W1(i,j,k-2)+W1(i,j,k-1)+W1(i+1,j,k-2)+W1(i+1,j,k-1))/4._KND
     Wloc=(W1(i,j,k-1)+W1(i,j,k)+W1(i,j+1,k-1)+W1(i,j+1,k))/4._KND
     V2(i,j,k)=(V(i,j,k)+V1(i,j,k)&
               -(dt/dxmin)*((V1(i,j,k)*Uloc)-(V1(i-2,j,k)*Ulocw))&
               -(dt/dymin)*((V1(i,j,k)*V1(i,j,k))-(V1(i,j-1,k)*V1(i,j-1,k)))&
               -(dt/dzmin)*((V1(i,j,k)*Wloc)-(V1(i,j,k-1)*Wlocb)))/2._KND
    enddo
   enddo
  enddo             

  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
     Ulocw=(U1(i-2,j,k)+U1(i-1,j,k)+U1(i-2,j+1,k)+U1(i-1,j+1,k))/4._KND
     Uloc=(U1(i-1,j,k)+U1(i,j,k)+U1(i-1,j,k+1)+U1(i,j,k+1))/4._KND
     Vlocs=(V1(i,j-2,k)+V1(i,j-1,k)+V1(i+1,j-2,k)+V1(i+1,j-1,k))/4._KND
     Vloc=(V1(i,j-1,k)+V1(i,j,k)+V1(i,j-1,k+1)+V1(i,j,k+1))/4._KND
     W2(i,j,k)=(W(i,j,k)+W1(i,j,k)&
               -(dt/dxmin)*((W1(i,j,k)*Uloc)-(W1(i-1,j,k)*Ulocw))&
               -(dt/dymin)*((W1(i,j,k)*Vloc)-(W1(i,j-1,k)*Vlocs))&
               -(dt/dzmin)*((W1(i,j,k)*W1(i,j,k))-(W1(i,j,k-1)*W1(i,j,k-1))))/2._KND
    enddo
   enddo
  enddo             


  U2=U2-U
  V2=V2-V
  W2=W2-W
 endsubroutine MC1



 subroutine MC2(U2,V2,W2,U,V,W)
 real(KND),dimension(-2:,-2:,-2:),intent(inout):: U2,V2,W2,U,V,W !U,V,W changed only by BC
 real(KND) U1(LBOUND(U,1):UBOUND(U,1),LBOUND(U,2):UBOUND(U,2),LBOUND(U,3):UBOUND(U,3))
 real(KND) V1(LBOUND(V,1):UBOUND(V,1),LBOUND(V,2):UBOUND(V,2),LBOUND(V,3):UBOUND(V,3))
 real(KND) W1(LBOUND(W,1):UBOUND(W,1),LBOUND(W,2):UBOUND(W,2),LBOUND(W,3):UBOUND(W,3))
 integer i,j,k
 real(KND) Uloc,Uloce,Ulocw,Vloc,Vlocn,Vlocs,Wloc,Wloct,Wlocb

  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)

  do k=1,Unz
   do j=1,Uny
    do i=1,Unx
     Vlocs=(V1(i,j-2,k)+V1(i,j-1,k)+V1(i+1,j-2,k)+V1(i+1,j-1,k))/4._KND
     Vloc=(V1(i,j-1,k)+V1(i,j,k)+V1(i+1,j-1,k)+V1(i+1,j,k))/4._KND
     Wlocb=(W1(i,j,k-2)+W1(i,j,k-1)+W1(i+1,j,k-2)+W1(i+1,j,k-1))/4._KND
     Wloc=(W1(i,j,k-1)+W1(i,j,k)+W1(i+1,j,k-1)+W1(i+1,j,k))/4._KND
     U1(i,j,k)=U(i,j,k)&
                       -(dt/dxmin)*((U1(i,j,k)*U1(i,j,k)-U1(i-1,j,k)*U1(i-1,j,k)))&
                       -(dt/dymin)*((U1(i,j,k)*Vloc)-(U1(i,j-1,k)*Vlocs))&
                       -(dt/dzmin)*((U1(i,j,k)*Wloc)-(U1(i,j,k-1)*Wlocb))
    enddo
   enddo
  enddo

  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
     Ulocw=(U1(i-2,j,k)+U1(i-1,j,k)+U1(i-2,j+1,k)+U1(i-1,j+1,k))/4._KND
     Uloc=(U1(i-1,j,k)+U1(i,j,k)+U1(i-1,j+1,k)+U1(i,j+1,k))/4._KND
     Wlocb=(W1(i,j,k-2)+W1(i,j,k-1)+W1(i+1,j,k-2)+W1(i+1,j,k-1))/4._KND
     Wloc=(W1(i,j,k-1)+W1(i,j,k)+W1(i,j+1,k-1)+W1(i,j+1,k))/4._KND
     V1(i,j,k)=V(i,j,k)&
               -(dt/dxmin)*((V1(i,j,k)*Uloc)-(V1(i-2,j,k)*Ulocw))&
               -(dt/dymin)*((V1(i,j,k)*V1(i,j,k))-(V1(i,j-1,k)*V1(i,j-1,k)))&
               -(dt/dzmin)*((V1(i,j,k)*Wloc)-(V1(i,j,k-1)*Wlocb))
    enddo
   enddo
  enddo             

  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
     Ulocw=(U1(i-2,j,k)+U1(i-1,j,k)+U1(i-2,j+1,k)+U1(i-1,j+1,k))/4._KND
     Uloc=(U1(i-1,j,k)+U1(i,j,k)+U1(i-1,j,k+1)+U1(i,j,k+1))/4._KND
     Vlocs=(V1(i,j-2,k)+V1(i,j-1,k)+V1(i+1,j-2,k)+V1(i+1,j-1,k))/4._KND
     Vloc=(V1(i,j-1,k)+V1(i,j,k)+V1(i,j-1,k+1)+V1(i,j,k+1))/4._KND
     W1(i,j,k)=W(i,j,k)&
               -(dt/dxmin)*((W1(i,j,k)*Uloc)-(W1(i-1,j,k)*Ulocw))&
               -(dt/dymin)*((W1(i,j,k)*Vloc)-(W1(i,j-1,k)*Vlocs))&
               -(dt/dzmin)*((W1(i,j,k)*W1(i,j,k))-(W1(i,j,k-1)*W1(i,j,k-1)))
    enddo
   enddo
  enddo             
 
  call Bound_CondU(U1)
  call Bound_CondV(V1)
  call Bound_CondW(W1)

  do k=1,Unz
   do j=1,Uny
    do i=1,Unx
     Vlocn=(V(i,j,k)+V(i,j+1,k)+V(i+1,j,k)+V(i+1,j+1,k))/4._KND
     Vloc=(V(i,j-1,k)+V(i,j,k)+V(i+1,j-1,k)+V(i+1,j,k))/4._KND
     Wloct=(W(i,j,k)+W(i,j,k+1)+W(i+1,j,k)+W(i+1,j,k+1))/4._KND
     Wloc=(W(i,j,k-1)+W(i,j,k)+W(i+1,j,k-1)+W(i+1,j,k))/4._KND
     U2(i,j,k)=(U(i,j,k)+U1(i,j,k)&
                       -(dt/dxmin)*((U(i+1,j,k)*U(i+1,j,k)-U(i,j,k)*U(i,j,k)))&
                       -(dt/dymin)*((U(i,j+1,k)*Vlocn)-(U(i,j,k)*Vloc))&
                       -(dt/dzmin)*((U(i,j,k+1)*Wloct)-(U(i,j,k)*Wloc)))/2._KND
    enddo
   enddo
  enddo

  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
     Uloce=(U(i,j,k)+U(i+1,j,k)+U(i,j+1,k)+U(i+1,j+1,k))/4._KND
     Ulocw=(U(i-1,j,k)+U(i,j,k)+U(i-1,j+1,k)+U(i,j+1,k))/4._KND
     Wloct=(W(i,j,k)+W(i,j,k+1)+W(i,j+1,k)+W(i,j+1,k+1))/4._KND
     Wloc=(W(i,j,k-1)+W(i,j,k)+W(i,j+1,k-1)+W(i,j+1,k))/4._KND
     V2(i,j,k)=(V(i,j,k)+V1(i,j,k)&
               -(dt/dxmin)*((V(i+1,j,k)*Uloce)-(V(i,j,k)*Uloc))&
               -(dt/dymin)*((V(i,j+1,k)*V(i,j+1,k))-(V(i,j,k)*V(i,j,k)))&
               -(dt/dzmin)*((V(i,j,k+1)*Wloct)-(V(i,j,k)*Wloc)))/2._KND
    enddo
   enddo
  enddo             

  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
     Uloce=(U(i,j,k)+U(i+1,j,k)+U(i,j,k+1)+U(i+1,j,k+1))/4._KND
     Uloc=(U(i-1,j,k)+U(i,j,k)+U(i-1,j,k+1)+U(i,j,k+1))/4._KND
     Vlocn=(V(i,j,k)+V(i,j+1,k)+V(i,j,k+1)+V(i,j+1,k+1))/4._KND
     Vloc=(V(i,j-1,k)+V(i,j,k)+V(i,j-1,k+1)+V(i,j,k+1))/4._KND
     W2(i,j,k)=(W(i,j,k)+W1(i,j,k)&
               -(dt/dxmin)*((W(i+1,j,k)*Uloce)-(W(i,j,k)*Uloc))&
               -(dt/dymin)*((W(i,j+1,k)*Vlocn)-(W(i,j,k)*Vloc))&
               -(dt/dzmin)*((W(i,j,k+1)*W(i,j,k+1))-(W(i,j,k)*W(i,j,k))))/2._KND
    enddo
   enddo
  enddo             


  U2=U2-U
  V2=V2-V
  W2=W2-W
 endsubroutine MC2

endmodule LAXWEND








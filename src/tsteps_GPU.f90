  !$hmpp <tsteps> UnifRedBlack codelet
  subroutine UNIFREDBLACK_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,&
                              Btype,sideU,&
                              dt,dxmin,dymin,dzmin,&
                              Uin,Vin,Win,U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                              coef,maxCNiter,epsCN,iters,residuum)
  implicit none

#ifdef __HMPP
   integer, parameter:: KND=4,TIM=4
#endif


   integer,intent(in):: Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz
   integer,intent(in):: Btype(6)
   real(KND),intent(in):: sideU(3,6)
   real(KND),dimension(-2:Prny+3,-2:Prnz+3),intent(in)          :: Uin,Vin,Win
   real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(in):: U
   real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(in):: V
   real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(in):: W
   real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(inout):: U2,U3
   real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(inout):: V2,V3
   real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(inout):: W2,W3
   real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(in):: Visc
   real(TIM),intent(in):: dt
   real(KND),intent(in):: dxmin,dymin,dzmin,coef,epsCN
   integer,intent(in):: maxCNiter
   integer,intent(out):: iters
   real(KND),intent(out):: residuum

   real(KND),dimension(1:Unx,1:Uny,1:Unz):: Apu
   real(KND),dimension(1:Vnx,1:Vny,1:Vnz):: ApV
   real(KND),dimension(1:Wnx,1:Wny,1:Wnz):: ApW
   real(KND) recdxmin2,recdymin2,recdzmin2                                                               !reciprocal values of dx**2
   real(KND) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
   integer i,j,k,l

   intrinsic mod,abs,max

       Ap=coef*dt/(2._KND)
       S=0
       l=0

       recdxmin2=1./dxmin**2
       recdymin2=1./dymin**2
       recdzmin2=1./dzmin**2

       !$hmppcg grid blocksize 512x1
       !$hmppcg permute (k,i,j)
       do k=1,Unz    !The explicit part, which doesn't have to be changed inside the loop
        do j=1,Uny
         do i=1,Unx
          U2(i,j,k)=U2(i,j,k)+Ap*(&
          ((Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))-&
          Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k)))*recdxmin2+0.25_KND*(&
           ((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))-&
           (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k)))*recdymin2+&
           ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))-&
           (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1)))*recdzmin2)))
         enddo
        enddo
       enddo
       !$hmppcg grid blocksize 512x1
       !$hmppcg permute (k,i,j)
       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          V2(i,j,k)=V2(i,j,k)+Ap*(&
          (0.25_KND*((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V(i+1,j,k)-V(i,j,k))-&
          (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(V(i,j,k)-V(i-1,j,k)))*recdxmin2+&
           (Visc(i,j+1,k)*(V(i,j+1,k)-V(i,j,k))-&
           Visc(i,j,k)*(V(i,j,k)-V(i,j-1,k)))*recdymin2+&
           0.25_KND*((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))-&
           (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(V(i,j,k)-V(i,j,k-1)))*recdzmin2))
         enddo
        enddo
       enddo
       !$hmppcg grid blocksize 512x1
       !$hmppcg permute (k,i,j)
       do k=1,Wnz
       do j=1,Wny
        do i=1,Wnx
         W2(i,j,k)=W2(i,j,k)+Ap*(&
         (0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W(i+1,j,k)-W(i,j,k))-&
         (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(W(i,j,k)-W(i-1,j,k)))*recdxmin2+&
          ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W(i,j+1,k)-W(i,j,k))-&
          (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(W(i,j,k)-W(i,j-1,k)))*recdymin2)+&
          (Visc(i,j,k+1)*(W(i,j,k+1)-W(i,j,k))-&
          Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1)))*recdzmin2))
        enddo
       enddo
      enddo

       !$hmppcg grid blocksize 512x1
       !$hmppcg permute (k,i,j)
       do k=1,Unz         !Auxiliary coefficients to better efficiency in loops
        do j=1,Uny
         do i=1,Unx
          ApU(i,j,k)=((Visc(i+1,j,k)+&
                      Visc(i,j,k))*recdxmin2+&
                      0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))+&
                      (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k)))*recdymin2+&
                      ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                      (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1)))*recdzmin2))
          ApU(i,j,k)=1._KND/(1._KND+Ap*ApU(i,j,k))
         enddo
        enddo
       enddo


       !$hmppcg grid blocksize 512x1
       !$hmppcg permute (k,i,j)
       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          ApV(i,j,k)=((Visc(i,j+1,k)+&
                     Visc(i,j,k))*recdymin2+&
                     0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))+&
                      (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k)))*recdxmin2+&
                     ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                     (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1)))*recdzmin2))
          ApV(i,j,k)=1._KND/(1._KND+Ap*ApV(i,j,k))
         enddo
        enddo
       enddo


       !$hmppcg grid blocksize 512x1
       !$hmppcg permute (k,i,j)
       do k=1,Wnz
        do j=1,Wny
         do i=1,Wnx
          ApW(i,j,k)=(0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                      (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k)))*recdxmin2+&
                     ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))+&
                     (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k)))*recdymin2)+&
                     (Visc(i,j,k+1)+&
                     Visc(i,j,k))*recdzmin2)
          ApW(i,j,k)=1._KND/(1._KND+Ap*ApW(i,j,k))
         enddo
        enddo
       enddo


       Suavg=0    !maximum values of velocities to norm the residues.
       !$hmppcg grid blocksize 512x1
       !$hmppcg gridify (k,i), reduce (max:Suavg)
       do k=1,Unz
        do i=1,Unx
         do j=1,Uny
          Suavg=max(Suavg,abs(U3(i,j,k)))
         enddo
        enddo
       enddo
       Svavg=0
       !$hmppcg grid blocksize 512x1
       !$hmppcg gridify (k,i), reduce (max:Svavg)
       do k=1,Vnz
        do i=1,Vnx
         do j=1,Vny
          Svavg=max(Svavg,abs(V3(i,j,k)))
         enddo
        enddo
       enddo
       Swavg=0
       !$hmppcg grid blocksize 512x1
       !$hmppcg gridify (k,i), reduce (max:Swavg)
       do k=1,Wnz
        do i=1,Wnx
         do j=1,Wny
          Swavg=max(Swavg,abs(W3(i,j,k)))
         enddo
        enddo
       enddo
       if (Suavg<=1e-3_KND) Suavg=1
       if (Svavg<=1e-3_KND) Svavg=1
       if (Swavg<=1e-3_KND) Swavg=1


       l=1
       S=epsCN+1.
       do while (S>epsCN.and.l<=maxCNiter)          !Gauss-Seidel iteration for Crank-Nicolson result
        call BoundU_GPU(1,Unx,Uny,Unz,Prny,Prnz,&
                             Btype,sideU,&
                             Uin,U3,0)
        call BoundU_GPU(2,Vnx,Vny,Vnz,Prny,Prnz,&
                             Btype,sideU,&
                             Vin,V3,0)
        call BoundU_GPU(3,Wnx,Wny,Wnz,Prny,Prnz,&
                             Btype,sideU,&
                             Win,W3,0)
        S=0
        Su=0
        Sv=0
        Sw=0
        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Su)
        do k=1,Unz
         do i=1,Unx
          do j=1+mod(i+k,2),Uny,2
            p=((Visc(i+1,j,k)*(U3(i+1,j,k))-&
             Visc(i,j,k)*(-U3(i-1,j,k)))*recdxmin2+&
             0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U3(i,j+1,k))-&
             (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(-U3(i,j-1,k)))*recdymin2+&
             ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U3(i,j,k+1))-&
             (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-U3(i,j,k-1)))*recdzmin2))
            p=Ap*p+U2(i,j,k)+U(i,j,k)
            p=p*ApU(i,j,k)


            Su=max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Sv)
        do k=1,Vnz
         do i=1,Vnx
          do j=1+mod(i+k,2),Vny,2
            p=((Visc(i,j+1,k)*(V3(i,j+1,k))-&
             Visc(i,j,k)*(-V3(i,j-1,k)))*recdymin2+&
             0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V3(i+1,j,k))-&
             (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(-V3(i-1,j,k)))*recdxmin2+&
             ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V3(i,j,k+1))-&
             (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-V3(i,j,k-1)))*recdzmin2))
            p=Ap*p+V2(i,j,k)+V(i,j,k)
            p=p*ApV(i,j,k)
            Sv=max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Sw)
        do k=1,Wnz
         do i=1,Wnx
          do j=1+mod(i+k,2),Wny,2
            p=(0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W3(i+1,j,k))-&
             (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(-W3(i-1,j,k)))*recdxmin2+&
             ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W3(i,j+1,k))-&
             (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(-W3(i,j-1,k)))*recdymin2)+&
             (Visc(i,j,k+1)*(W3(i,j,k+1))-&
             Visc(i,j,k)*(-W3(i,j,k-1)))*recdzmin2)
            p=Ap*p+W2(i,j,k)+W(i,j,k)
            p=p*ApW(i,j,k)
            Sw=max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo


        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Su)
        do k=1,Unz
         do i=1,Unx
          do j=1+mod(i+k+1,2),Uny,2
            p=((Visc(i+1,j,k)*(U3(i+1,j,k))-&
             Visc(i,j,k)*(-U3(i-1,j,k)))*recdxmin2+&
             0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U3(i,j+1,k))-&
             (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(-U3(i,j-1,k)))*recdymin2+&
             ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U3(i,j,k+1))-&
             (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-U3(i,j,k-1)))*recdzmin2))
            p=Ap*p+U2(i,j,k)+U(i,j,k)
            p=p*ApU(i,j,k)


            Su=max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Sv)
        do k=1,Vnz
         do i=1,Vnx
          do j=1+mod(i+k+1,2),Vny,2
            p=((Visc(i,j+1,k)*(V3(i,j+1,k))-&
             Visc(i,j,k)*(-V3(i,j-1,k)))*recdymin2+&
             0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V3(i+1,j,k))-&
             (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(-V3(i-1,j,k)))*recdxmin2+&
             ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V3(i,j,k+1))-&
             (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-V3(i,j,k-1)))*recdzmin2))
            p=Ap*p+V2(i,j,k)+V(i,j,k)
            p=p*ApV(i,j,k)
            Sv=max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Sw)
        do k=1,Wnz
         do i=1,Wnx
          do j=1+mod(i+k+1,2),Wny,2
            p=(0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W3(i+1,j,k))-&
             (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(-W3(i-1,j,k)))*recdxmin2+&
             ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W3(i,j+1,k))-&
             (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(-W3(i,j-1,k)))*recdymin2)+&
             (Visc(i,j,k+1)*(W3(i,j,k+1))-&
             Visc(i,j,k)*(-W3(i,j,k-1)))*recdzmin2)
            p=Ap*p+W2(i,j,k)+W(i,j,k)
            p=p*ApW(i,j,k)
            Sw=max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        S=max(Su/Suavg,Sv/Svavg,Sw/Swavg)
        l=l+1
       enddo
    iters=l-1
    residuum=S


    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
       U2(i,j,k)=U3(i,j,k)
      enddo
     enddo
    enddo
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
       V2(i,j,k)=V3(i,j,k)
      enddo
     enddo
    enddo
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
       W2(i,j,k)=W3(i,j,k)
      enddo
     enddo
    enddo

    call BoundU_GPU(1,Unx,Uny,Unz,Prny,Prnz,&
                         Btype,sideU,&
                         Uin,U2,0)
    call BoundU_GPU(2,Vnx,Vny,Vnz,Prny,Prnz,&
                         Btype,sideU,&
                         Vin,V2,0)
    call BoundU_GPU(3,Wnx,Wny,Wnz,Prny,Prnz,&
                         Btype,sideU,&
                         Win,W2,0)
  endsubroutine UNIFREDBLACK_GPU


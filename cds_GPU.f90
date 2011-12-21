
    !$hmpp CDS_GPU codelet, target=CUDA
    subroutine CDS_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dx,dy,dz,dt,U2,V2,W2,U,V,W)
    implicit none
    integer, parameter:: KND=4

    integer,intent(in)    :: Unx, Uny, Unz, Vnx, Vny, Vnz, Wnx, Wny, Wnz
    real(KND),intent(in)  :: dx, dy, dz, dt
    real(KND),intent(out) :: U2(-2:Unx+3,-2:Uny+3,-2:Unz+3)
    real(KND),intent(in)  :: U( -2:Unx+3,-2:Uny+3,-2:Unz+3)
    real(KND),intent(out) :: V2(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
    real(KND),intent(in)  :: V( -2:Vnx+3,-2:Vny+3,-2:Vnz+3)
    real(KND),intent(out) :: W2(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
    real(KND),intent(in)  :: W( -2:Wnx+3,-2:Wny+3,-2:Wnz+3)
    integer i,j,k
    real(KND) Ax,Ay,Az


        Ax=0.25_KND*dt/dx
        Ay=0.25_KND*dt/dy
        Az=0.25_KND*dt/dz
       !$hmppcg grid blocksize 512x1
       !$hmppcg permute (k,i,j)
        do k=1,Unz
            do j=1,Uny
                do i=1,Unx
                    U2(i,j,k)= - ((Ax*(U(i+1,j,k)+U(i,j,k))*(U(i+1,j,k)+U(i,j,k))&
                    -Ax*(U(i,j,k)+U(i-1,j,k))*(U(i,j,k)+U(i-1,j,k)))&
                    +(Ay*(U(i,j+1,k)+U(i,j,k))*(V(i+1,j,k)+V(i,j,k))&
                    -Ay*(U(i,j,k)+U(i,j-1,k))*(V(i+1,j-1,k)+V(i,j-1,k)))&
                    +(Az*(U(i,j,k+1)+U(i,j,k))*(W(i+1,j,k)+W(i,j,k))&
                    -Az*(U(i,j,k)+U(i,j,k-1))*(W(i+1,j,k-1)+W(i,j,k-1))))
                enddo
            enddo
        enddo

       !$hmppcg grid blocksize 512x1
       !$hmppcg permute (k,i,j)
        do k=1,Vnz
            do j=1,Vny
                do i=1,Vnx
                    V2(i,j,k)= - ((Ay*(V(i,j+1,k)+V(i,j,k))*(V(i,j+1,k)+V(i,j,k))&
                    -Ay*(V(i,j,k)+V(i,j-1,k))*(V(i,j,k)+V(i,j-1,k)))&
                    +(Ax*(V(i+1,j,k)+V(i,j,k))*(U(i,j+1,k)+U(i,j,k))&
                    -Ax*(V(i,j,k)+V(i-1,j,k))*(U(i-1,j+1,k)+U(i-1,j,k)))&
                    +(Az*(V(i,j,k+1)+V(i,j,k))*(W(i,j+1,k)+W(i,j,k))&
                    -Az*(V(i,j,k)+V(i,j,k-1))*(W(i,j+1,k-1)+W(i,j,k-1))))
                enddo
            enddo
        enddo

       !$hmppcg grid blocksize 512x1
       !$hmppcg permute (k,i,j)
        do k=1,Wnz
            do j=1,Wny
                do i=1,Wnx
                    W2(i,j,k)= - ((Az*(W(i,j,k+1)+W(i,j,k))*(W(i,j,k+1)+W(i,j,k))&
                    -Az*(W(i,j,k)+W(i,j,k-1))*(W(i,j,k)+W(i,j,k-1)))&
                    +(Ay*(W(i,j+1,k)+W(i,j,k))*(V(i,j,k+1)+V(i,j,k))&
                    -Ay*(W(i,j,k)+W(i,j-1,k))*(V(i,j-1,k)+V(i,j-1,k+1)))&
                    +(Ax*(W(i+1,j,k)+W(i,j,k))*(U(i,j,k+1)+U(i,j,k))&
                    -Ax*(W(i,j,k)+W(i-1,j,k))*(U(i-1,j,k+1)+U(i-1,j,k))))
                enddo
            enddo
        enddo

    end subroutine CDS_GPU


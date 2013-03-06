module CDS
    use PARAMETERS
    use BOUNDARIES
    use LIMITERS, only: FluxLimiter
    use ArrayUtilities
    use Tiling

    implicit none

!     private
!     public CDU, CDV, CDW, CDS4U, CDS4V, CDS4W

    contains





    subroutine CDUdiv(U2,U,V,W)
    real(knd) :: U2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer nx,ny,nz,i,j,k
    real(knd) Ax,Ay,Az,Utmp
    real(knd),allocatable,dimension(:),save:: fe,fp,gn,ht
    integer,save:: called=0

     nx=Unx
     ny=Uny
     nz=Unz

     if (gridtype==UNIFORMGRID) then
        Ax=0.25*dt/dxmin
        Ay=0.25*dt/dymin
        Az=0.25*dt/dzmin

        !$omp parallel do private(i,j,k)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    U2(i,j,k)= - ((Ax*(U(i+1,j,k)+U(i,j,k))*(U(i+1,j,k)+U(i,j,k)) &
                    -Ax*(U(i,j,k)+U(i-1,j,k))*(U(i,j,k)+U(i-1,j,k))) &
                    +(Ay*(U(i,j+1,k)+U(i,j,k))*(V(i+1,j,k)+V(i,j,k)) &
                    -Ay*(U(i,j,k)+U(i,j-1,k))*(V(i+1,j-1,k)+V(i,j-1,k))) &
                    +(Az*(U(i,j,k+1)+U(i,j,k))*(W(i+1,j,k)+W(i,j,k)) &
                    -Az*(U(i,j,k)+U(i,j,k-1))*(W(i+1,j,k-1)+W(i,j,k-1))))
                end do
            end do
        end do
        !$omp end parallel do
     else
        if (called==0) then
         allocate(fe(0:nx),fp(1:nx),gn(0:ny),ht(0:nz))
         called=1
         forall (i=0:nx)      fe(i)=(xPr(i+1)-xU(i))/(xU(i+1)-xU(i))
         forall (i=1:nx)      fp(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
         forall (j=0:ny)      gn(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
         forall (k=0:nz)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
        end if

        !$omp parallel do private(i,j,k)
        do k=1,nz
         do j=1,ny
          do i=1,nx
           Utmp=( ((fe(i)*U(i+1,j,k)+(1-fe(i))*U(i,j,k))*(fe(i)*U(i+1,j,k)+(1-fe(i))*U(i,j,k))) &
                      -((fe(i-1)*U(i,j,k)+(1-fe(i-1))*U(i-1,j,k))*(fe(i-1)*U(i,j,k)+(1-fe(i-1))*U(i-1,j,k))))/dxU(i)
           Utmp=Utmp+( (gn(j)*U(i,j+1,k)+(1-gn(j))*U(i,j,k))*(fp(i)*V(i+1,j,k)+(1-fp(i))*V(i,j,k)) &
                      -(gn(j-1)*U(i,j,k)+(1-gn(j-1))*U(i,j-1,k))*(fp(i)*V(i+1,j-1,k)+(1-fp(i))*V(i,j-1,k)))/dyPr(j)
           Utmp=Utmp+( (ht(k)*U(i,j,k+1)+(1-ht(k))*U(i,j,k))*(fp(i)*W(i+1,j,k)+(1-fp(i))*W(i,j,k)) &
                      -(ht(k-1)*U(i,j,k)+(1-ht(k-1))*U(i,j,k-1))*(fp(i)*W(i+1,j,k-1)+(1-fp(i))*W(i,j,k-1)))/dzPr(k)
           U2(i,j,k)=-dt*Utmp
          end do
         end do
        end do
        !$omp end parallel do
     end if
    end subroutine CDUdiv






    subroutine CDVdiv(V2,U,V,W)
    real(knd) :: V2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer nx,ny,nz,i,j,k
    real(knd) Ax,Ay,Az,Vtmp
    real(knd),allocatable,dimension(:),save:: gn,gp,fe,ht
    integer,save:: called=0


     nx=Vnx
     ny=Vny
     nz=Vnz

     if (gridtype==UNIFORMGRID) then
        Ax=0.25*dt/dxmin
        Ay=0.25*dt/dymin
        Az=0.25*dt/dzmin

        !$omp parallel do private(i,j,k)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    V2(i,j,k)= - ((Ay*(V(i,j+1,k)+V(i,j,k))*(V(i,j+1,k)+V(i,j,k)) &
                    -Ay*(V(i,j,k)+V(i,j-1,k))*(V(i,j,k)+V(i,j-1,k))) &
                    +(Ax*(V(i+1,j,k)+V(i,j,k))*(U(i,j+1,k)+U(i,j,k)) &
                    -Ax*(V(i,j,k)+V(i-1,j,k))*(U(i-1,j+1,k)+U(i-1,j,k))) &
                    +(Az*(V(i,j,k+1)+V(i,j,k))*(W(i,j+1,k)+W(i,j,k)) &
                    -Az*(V(i,j,k)+V(i,j,k-1))*(W(i,j+1,k-1)+W(i,j,k-1))))
                end do
            end do
        end do
        !$omp end parallel do
     else
        if (called==0) then
         allocate(gn(0:ny),gp(1:ny),fe(0:nx),ht(0:nz))
         called=1
         forall (j=0:ny)      gn(j)=(yPr(j+1)-yV(j))/(yV(j+1)-yV(j))
         forall (j=1:ny)      gp(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
         forall (i=0:nx)      fe(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
         forall (k=0:nz)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
        end if

        !$omp parallel do private(i,j,k)
        do k=1,nz
         do j=1,ny
          do i=1,nx
           Vtmp=( ((gn(j)*V(i,j+1,k)+(1-gn(j))*V(i,j,k))*(gn(j)*V(i,j+1,k)+(1-gn(j))*V(i,j,k))) &
                      -((gn(j-1)*V(i,j,k)+(1-gn(j-1))*V(i,j-1,k))*(gn(j-1)*V(i,j,k)+(1-gn(j-1))*V(i,j-1,k))))/dyV(j)
           Vtmp=Vtmp+( (fe(i)*V(i+1,j,k)+(1-fe(i))*V(i,j,k))*(gp(j)*U(i,j+1,k)+(1-gp(j))*U(i,j,k)) &
                      -(fe(i-1)*V(i,j,k)+(1-fe(i-1))*V(i-1,j,k))*(gp(j)*U(i-1,j+1,k)+(1-gp(j))*U(i-1,j,k)))/dxPr(i)
           Vtmp=Vtmp+( (ht(k)*V(i,j,k+1)+(1-ht(k))*V(i,j,k))*(gp(j)*W(i,j+1,k)+(1-gp(j))*W(i,j,k)) &
                      -(ht(k-1)*V(i,j,k)+(1-ht(k-1))*V(i,j,k-1))*(gp(j)*W(i,j+1,k-1)+(1-gp(j))*W(i,j,k-1)))/dzPr(k)
           V2(i,j,k)=-dt*Vtmp
          end do
         end do
        end do
        !$omp end parallel do
     end if
    end subroutine CDVdiv





    subroutine CDWdiv(W2,U,V,W)
    real(knd) :: W2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer nx,ny,nz,i,j,k
    real(knd) Ax,Ay,Az,Wtmp
    real(knd),allocatable,dimension(:),save:: ht,hp,fe,gn
    integer,save:: called=0

     nx=Wnx
     ny=Wny
     nz=Wnz

     if (gridtype==UNIFORMGRID) then
        Ax=0.25*dt/dxmin
        Ay=0.25*dt/dymin
        Az=0.25*dt/dzmin

        !$omp parallel do private(i,j,k)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    W2(i,j,k)= - ((Az*(W(i,j,k+1)+W(i,j,k))*(W(i,j,k+1)+W(i,j,k)) &
                    -Az*(W(i,j,k)+W(i,j,k-1))*(W(i,j,k)+W(i,j,k-1))) &
                    +(Ay*(W(i,j+1,k)+W(i,j,k))*(V(i,j,k+1)+V(i,j,k)) &
                    -Ay*(W(i,j,k)+W(i,j-1,k))*(V(i,j-1,k)+V(i,j-1,k+1))) &
                    +(Ax*(W(i+1,j,k)+W(i,j,k))*(U(i,j,k+1)+U(i,j,k)) &
                    -Ax*(W(i,j,k)+W(i-1,j,k))*(U(i-1,j,k+1)+U(i-1,j,k))))
                end do
            end do
        end do
        !$omp end parallel do
     else
        if (called==0) then
         allocate(ht(0:nz),hp(1:nz),fe(0:nx),gn(0:ny))
         called=1
         forall (k=0:nz)      ht(k)=(zPr(k+1)-zW(k))/(zW(k+1)-zW(k))
         forall (k=1:nz)      hp(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
         forall (i=0:nx)      fe(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
         forall (j=0:ny)      gn(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
        end if

        !$omp parallel do private(i,j,k)
        do k=1,nz
         do j=1,ny
          do i=1,nx
           Wtmp=( ((ht(k)*W(i,j,k+1)+(1-ht(k))*W(i,j,k))*(ht(k)*W(i,j,k+1)+(1-ht(k))*W(i,j,k))) &
                      -((ht(k-1)*W(i,j,k)+(1-ht(k-1))*W(i,j,k-1))*(ht(k-1)*W(i,j,k)+(1-ht(k-1))*W(i,j,k-1))))/dzW(k)
           Wtmp=Wtmp+( (fe(i)*W(i+1,j,k)+(1-fe(i))*W(i,j,k))*(hp(k)*U(i,j,k+1)+(1-hp(k))*U(i,j,k)) &
                      -(fe(i-1)*W(i,j,k)+(1-fe(i-1))*W(i-1,j,k))*(hp(k)*U(i-1,j,k+1)+(1-hp(k))*U(i-1,j,k)))/dxPr(i)
           Wtmp=Wtmp+( (gn(j)*W(i,j+1,k)+(1-gn(j))*W(i,j,k))*(hp(k)*V(i,j,k+1)+(1-hp(k))*V(i,j,k)) &
                      -(gn(j-1)*W(i,j,k)+(1-gn(j-1))*W(i,j-1,k))*(hp(k)*V(i,j-1,k+1)+(1-hp(k))*V(i,j-1,k)))/dzPr(k)
           W2(i,j,k)=-dt*Wtmp
          end do
         end do
        end do
        !$omp end parallel do
     end if
    end subroutine CDWdiv











    subroutine CDUadv(U2,U,V,W)
    real(knd) :: U2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer i,j,k
    real(knd) Ax,Ay,Az,Vadv,Wadv


     if (gridtype==UNIFORMGRID) then
        Ax=0.5*dt/dxmin
        Ay=0.125*dt/dymin
        Az=0.125*dt/dzmin

        !$omp parallel do private(i,j,k,Vadv,Wadv)
        do k=1,Unz
            do j=1,Uny
                do i=1,Unx
                    Vadv = ( V(i,j,k) + V(i+1,j,k) + V(i,j-1,k) + V(i+1,j-1,k) )
                    Wadv = ( W(i,j,k) + W(i+1,j,k) + W(i,j,k-1) + W(i+1,j,k-1) )
                    U2(i,j,k)= U2(i,j,k) &
                               - (Ax*(U(i+1,j,k)-U(i-1,j,k)) * U(i,j,k) &
                               +  Ay*(U(i,j+1,k)-U(i,j-1,k)) * Vadv&
                               +  Az*(U(i,j,k+1)-U(i,j,k-1)) * Wadv )
                end do
            end do
        end do

     end if
    end subroutine CDUadv






    subroutine CDVadv(V2,U,V,W)
    real(knd) :: V2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer i,j,k
    real(knd) Ax,Ay,Az,Uadv,Wadv


     if (gridtype==UNIFORMGRID) then
        Ax=0.125*dt/dxmin
        Ay=0.5*dt/dymin
        Az=0.125*dt/dzmin

        !$omp parallel do private(i,j,k,Uadv,Wadv)
        do k=1,Vnz
            do j=1,Vny
                do i=1,Vnx
                    Uadv = ( U(i,j,k) + U(i,j+1,k) + U(i-1,j,k) + U(i-1,j+1,k) )
                    Wadv = ( W(i,j,k) + W(i,j+1,k) + W(i,j,k-1) + W(i,j+1,k-1) )
                    V2(i,j,k)= V2(i,j,k) &
                               - (Ax*(V(i+1,j,k)-V(i-1,j,k)) * Uadv&
                               +  Ay*(V(i,j+1,k)-V(i,j-1,k)) * V(i,j,k) &
                               +  Az*(V(i,j,k+1)-V(i,j,k-1)) * Wadv )
                end do
            end do
        end do
        !$omp end parallel do

     end if
    end subroutine CDVadv





    subroutine CDWadv(W2,U,V,W)
    real(knd) :: W2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer i,j,k
    real(knd) Ax,Ay,Az,Uadv,Vadv


     if (gridtype==UNIFORMGRID) then
        Ax=0.125*dt/dxmin
        Ay=0.125*dt/dymin
        Az=0.5*dt/dzmin

        !$omp parallel do private(i,j,k,Uadv,Vadv)
        do k=1,Wnz
            do j=1,Wny
                do i=1,Wnx
                    Uadv = ( U(i,j,k) + U(i,j,k+1) + U(i-1,j,k) + U(i-1,j,k+1) )
                    Vadv = ( V(i,j,k) + V(i,j,k+1) + V(i,j-1,k) + V(i,j-1,k+1) )
                    W2(i,j,k)= W2(i,j,k) &
                               - (Ax*(W(i+1,j,k)-W(i-1,j,k)) * Uadv&
                               +  Ay*(W(i,j+1,k)-W(i,j-1,k)) * Vadv&
                               +  Az*(W(i,j,k+1)-W(i,j,k-1)) * W(i,j,k) )
                end do
            end do
        end do
        !$omp end parallel do

     end if
    end subroutine CDWadv















    subroutine CDU(U2,U,V,W)
    real(knd) :: U2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)

     U2 = 0
     call CDUdiv(U2,U,V,W)
     call CDUadv(U2,U,V,W)
     U2 = U2 / 2
    end subroutine CDU






    subroutine CDV(V2,U,V,W)
    real(knd) :: V2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)

     V2 = 0
     call CDVdiv(V2,U,V,W)
     call CDVadv(V2,U,V,W)
     V2 = V2 / 2

    end subroutine CDV





    subroutine CDW(W2,U,V,W)
    real(knd) :: W2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)

     W2 = 0
     call CDWdiv(W2,U,V,W)
     call CDWadv(W2,U,V,W)
     W2 = W2 / 2

    end subroutine CDW
















    subroutine CDS4U(U2,U,V,W)
      real(knd),dimension(-2:,-2:,-2:),intent(out) :: U2
      real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: mask4(4) = (/ -1, 9, 9, -1 /)
      integer,parameter :: divcoef = 256
      integer,parameter :: mask2(2) = (/ 8, 8 /)
      integer   :: i,j,k,l
      real(knd) :: flux

      !$omp parallel private(i,j,k,l,flux)
      !$omp workshare
      U2 = 0
      !$omp end workshare

      !$omp do
      do k=1,Unz
        do j=1,Uny
          do i=1,Unx+1
            if (Utype(i-1,j,k)<=0 .and. Utype(i,j,k)<=0) then

              flux = ( ( sum(mask4 * U(i-2:i+1,j,k)) )**2 ) / dxmin
              U2(i-1,j,k) = U2(i-1,j,k) - flux
              U2(i,j,k)   = U2(i,j,k)   + flux

            else if (Utype(i-1,j,k)<=0 .or. Utype(i,j,k)<=0) then

              flux = ( ( sum(mask2 * U(i-1:i,j,k)) )**2 ) / dxmin
              U2(i-1,j,k) = U2(i-1,j,k) - flux
              U2(i,j,k)   = U2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do

      !$omp do
      do k=1,Unz
        do j=1,Uny+1
          do i=1,Unx
            if (Utype(i,j-1,k)<=0 .and. Utype(i,j,k)<=0) then

              flux = ( sum(mask4 * U(i,j-2:j+1,k)) * sum(mask4 * V(i-1:i+2,j-1,k)) ) / dymin
              U2(i,j-1,k) = U2(i,j-1,k) - flux
              U2(i,j,k)   = U2(i,j,k)   + flux

            else if (Utype(i,j-1,k)<=0 .or. Utype(i,j,k)<=0) then

              flux = ( sum(mask2 * U(i,j-1:j,k)) * sum(mask2 * V(i:i+1,j-1,k)) ) / dymin
              U2(i,j-1,k) = U2(i,j-1,k) - flux
              U2(i,j,k)   = U2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do

      do l=1,2
        !$omp do
        do k=l,Unz+1,2
          do j=1,Uny
            do i=1,Unx
              if (Utype(i,j,k-1)<=0 .and. Utype(i,j,k)<=0) then

                  flux = ( sum(mask4 * U(i,j,k-2:k+1))) * sum(mask4 * W(i-1:i+2,j,k-1)) / dzmin
                  U2(i,j,k-1) = U2(i,j,k-1) - flux
                  U2(i,j,k)   = U2(i,j,k)   + flux

              else if (Utype(i,j,k-1)<=0 .or. Utype(i,j,k)<=0) then

                  flux = ( sum(mask2 * U(i,j,k-1:k))) * sum(mask2 * W(i:i+1,j,k-1)) / dzmin
                  U2(i,j,k-1) = U2(i,j,k-1) - flux
                  U2(i,j,k)   = U2(i,j,k)   + flux

              end if
            end do
          end do
        end do
        !$omp end do
      end do

      !$omp workshare
      U2 = U2 * dt / divcoef
      !$omp end workshare
      !$omp end parallel

    end subroutine CDS4U


    subroutine CDS4V(V2,U,V,W)
      real(knd),dimension(-2:,-2:,-2:),intent(out) :: V2
      real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: mask4(4) = (/ -1, 9, 9, -1 /)
      integer,parameter :: divcoef = 256
      integer,parameter :: mask2(2) = (/ 8, 8 /)
      integer   :: i,j,k,l
      real(knd) :: flux

      !$omp parallel private(i,j,k,l,flux)
      !$omp workshare
      V2 = 0
      !$omp end workshare

      !$omp do
      do k=1,Vnz
        do j=1,Vny
          do i=1,Vnx+1
            if (Vtype(i-1,j,k)<=0 .and. Vtype(i,j,k)<=0) then

              flux = ( sum(mask4 * V(i-2:i+1,j,k)) * sum(mask4 * U(i-1,j-1:j+2,k)) ) / dxmin
              V2(i-1,j,k) = V2(i-1,j,k) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            else if (Vtype(i-1,j,k)<=0 .or. Vtype(i,j,k)<=0) then

              flux = ( sum(mask2 * V(i-1:i,j,k)) * sum(mask2 * U(i-1,j:j+1,k)) ) / dxmin
              V2(i-1,j,k) = V2(i-1,j,k) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do


      !$omp do
      do k=1,Vnz
        do j=1,Vny+1
          do i=1,Vnx
            if (Vtype(i,j-1,k)<=0 .and. Vtype(i,j,k)<=0) then

              flux = ( ( sum(mask4 * V(i,j-2:j+1,k)) )**2 ) / dymin
              V2(i,j-1,k) = V2(i,j-1,k) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            else if (Vtype(i,j-1,k)<=0 .or. Vtype(i,j,k)<=0) then

              flux = ( ( sum(mask2 * V(i,j-1:j,k)) )**2 ) / dymin
              V2(i,j-1,k) = V2(i,j-1,k) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do

      do l=1,2
        !$omp do
        do k=l,Vnz+1,2
          do j=1,Vny
            do i=1,Vnx
              if (Vtype(i,j,k-1)<=0 .and. Vtype(i,j,k)<=0) then

                flux = ( sum(mask4 * V(i,j,k-2:k+1)) * sum(mask4 * W(i,j-1:j+2,k-1)) ) / dzmin
                V2(i,j,k-1) = V2(i,j,k-1) - flux
                V2(i,j,k)   = V2(i,j,k)   + flux

              else if (Vtype(i,j,k-1)<=0 .or. Vtype(i,j,k)<=0) then

                flux = ( sum(mask2 * V(i,j,k-1:k)) * sum(mask2 * W(i,j:j+1,k-1)) ) / dzmin
                V2(i,j,k-1) = V2(i,j,k-1) - flux
                V2(i,j,k)   = V2(i,j,k)   + flux

              end if
            end do
          end do
        end do
        !$omp end do
      end do

      !$omp workshare
      V2 = V2 * dt / divcoef
      !$omp end workshare
      !$omp end parallel

    end subroutine CDS4V


    subroutine CDS4W(W2,U,V,W)
      real(knd),dimension(-2:,-2:,-2:),intent(out) :: W2
      real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: mask4(4) = (/ -1, 9, 9, -1 /)
      integer,parameter :: divcoef = 256
      integer,parameter :: mask2(2) = (/  8, 8 /)
      integer   :: i,j,k,l
      real(knd) :: flux

      !$omp parallel private(i,j,k,l,flux)
      !$omp workshare
      W2 = 0
      !$omp end workshare

      !$omp do
      do k=1,Wnz
        do j=1,Wny
          do i=1,Wnx+1
            if (Wtype(i-1,j,k)<=0 .and. Wtype(i,j,k)<=0) then

              flux = ( sum(mask4 * W(i-2:i+1,j,k)) * sum(mask4 * U(i-1,j,k-1:k+2)) ) / dxmin
              W2(i-1,j,k) = W2(i-1,j,k) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            else if (Wtype(i-1,j,k)<=0 .or. Wtype(i,j,k)<=0) then

              flux = ( sum(mask2 * W(i-1:i,j,k)) * sum(mask2 * U(i-1,j,k:k+1)) ) / dxmin
              W2(i-1,j,k) = W2(i-1,j,k) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do

      !$omp do
      do k=1,Wnz
        do j=1,Wny+1
          do i=1,Wnx
            if (Wtype(i,j-1,k)<=0 .and. Wtype(i,j,k)<=0) then

              flux = ( sum(mask4 * W(i,j-2:j+1,k)) * sum(mask4 * V(i,j-1,k-1:k+2)) ) / dymin
              W2(i,j-1,k) = W2(i,j-1,k) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            else if (Wtype(i,j-1,k)<=0 .or. Wtype(i,j,k)<=0) then

              flux = ( sum(mask2 * W(i,j-1:j,k)) * sum(mask2 * V(i,j-1,k:k+1)) ) / dymin
              W2(i,j-1,k) = W2(i,j-1,k) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do

      do l=1,2
        !$omp do
        do k=l,Wnz+1,2
          do j=1,Wny
            do i=1,Wnx
              if (Wtype(i,j,k-1)<=0 .and. Wtype(i,j,k)<=0) then

                flux = ( ( sum(mask4 * W(i,j,k-2:k+1)) )**2 ) / dzmin
                W2(i,j,k-1) = W2(i,j,k-1) - flux
                W2(i,j,k)   = W2(i,j,k)   + flux

              else if (Wtype(i,j,k-1)<=0 .or. Wtype(i,j,k)<=0) then

                flux = ( ( sum(mask2 * W(i,j,k-1:k)) )**2 ) / dzmin
                W2(i,j,k-1) = W2(i,j,k-1) - flux
                W2(i,j,k)   = W2(i,j,k)   + flux

              end if
            end do
          end do
        end do
        !$omp end do
      end do

      !$omp workshare
      W2 = W2 * dt / divcoef
      !$omp end workshare
      !$omp end parallel

    end subroutine CDS4W

















    subroutine CDS4_2U(U2,U,V,W)
      real(knd),dimension(-2:,-2:,-2:),intent(out) :: U2
      real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: mask4(4) = (/ -1, 9, 9, -1 /)
      integer,parameter :: divcoef = 256
      integer,parameter :: mask2(2) = (/ 8, 8 /)
      integer   :: i,j,k,l
      real(knd) :: flux

      !$omp parallel private(i,j,k,l,flux)
      !$omp workshare
      U2 = 0
      !$omp end workshare

      !$omp do
      do k=1,Unz
        do j=1,Uny
          do i=1,Unx+1
            if (Utype(i-1,j,k)<=0 .and. Utype(i,j,k)<=0) then

              flux = ( ( sum(mask4 * U(i-2:i+1,j,k)) )**2 ) / dxmin
              U2(i-1,j,k) = U2(i-1,j,k) - flux
              U2(i,j,k)   = U2(i,j,k)   + flux

            else if (Utype(i-1,j,k)<=0 .or. Utype(i,j,k)<=0) then

              flux = ( ( sum(mask2 * U(i-1:i,j,k)) )**2 ) / dxmin
              U2(i-1,j,k) = U2(i-1,j,k) - flux
              U2(i,j,k)   = U2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do

      !$omp do
      do k=1,Unz
        do j=1,Uny+1
          do i=1,Unx
            if (Utype(i,j-1,k)<=0 .and. Utype(i,j,k)<=0) then

              flux = ( sum(mask4 * U(i,j-2:j+1,k)) * sum(mask2 * V(i:i+1,j-1,k)) ) / dymin
              U2(i,j-1,k) = U2(i,j-1,k) - flux
              U2(i,j,k)   = U2(i,j,k)   + flux

            else if (Utype(i,j-1,k)<=0 .or. Utype(i,j,k)<=0) then

              flux = ( sum(mask2 * U(i,j-1:j,k)) * sum(mask2 * V(i:i+1,j-1,k)) ) / dymin
              U2(i,j-1,k) = U2(i,j-1,k) - flux
              U2(i,j,k)   = U2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do

      do l=1,2
        !$omp do
        do k=l,Unz+1,2
          do j=1,Uny
            do i=1,Unx
              if (Utype(i,j,k-1)<=0 .and. Utype(i,j,k)<=0) then

                  flux = ( sum(mask4 * U(i,j,k-2:k+1))) * sum(mask2 * W(i:i+1,j,k-1)) / dzmin
                  U2(i,j,k-1) = U2(i,j,k-1) - flux
                  U2(i,j,k)   = U2(i,j,k)   + flux

              else if (Utype(i,j,k-1)<=0 .or. Utype(i,j,k)<=0) then

                  flux = ( sum(mask2 * U(i,j,k-1:k))) * sum(mask2 * W(i:i+1,j,k-1)) / dzmin
                  U2(i,j,k-1) = U2(i,j,k-1) - flux
                  U2(i,j,k)   = U2(i,j,k)   + flux

              end if
            end do
          end do
        end do
        !$omp end do
      end do

      !$omp workshare
      U2 = U2 * dt / divcoef
      !$omp end workshare
      !$omp end parallel

    end subroutine CDS4_2U


    subroutine CDS4_2V(V2,U,V,W)
      real(knd),dimension(-2:,-2:,-2:),intent(out) :: V2
      real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: mask4(4) = (/ -1, 9, 9, -1 /)
      integer,parameter :: divcoef = 256
      integer,parameter :: mask2(2) = (/ 8, 8 /)
      integer   :: i,j,k,l
      real(knd) :: flux

      !$omp parallel private(i,j,k,l,flux)
      !$omp workshare
      V2 = 0
      !$omp end workshare

      !$omp do
      do k=1,Vnz
        do j=1,Vny
          do i=1,Vnx+1
            if (Vtype(i-1,j,k)<=0 .and. Vtype(i,j,k)<=0) then

              flux = ( sum(mask4 * V(i-2:i+1,j,k)) * sum(mask2 * U(i-1,j:j+1,k)) ) / dxmin
              V2(i-1,j,k) = V2(i-1,j,k) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            else if (Vtype(i-1,j,k)<=0 .or. Vtype(i,j,k)<=0) then

              flux = ( sum(mask2 * V(i-1:i,j,k)) * sum(mask2 * U(i-1,j:j+1,k)) ) / dxmin
              V2(i-1,j,k) = V2(i-1,j,k) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do


      !$omp do
      do k=1,Vnz
        do j=1,Vny+1
          do i=1,Vnx
            if (Vtype(i,j-1,k)<=0 .and. Vtype(i,j,k)<=0) then

              flux = ( ( sum(mask4 * V(i,j-2:j+1,k)) )**2 ) / dymin
              V2(i,j-1,k) = V2(i,j-1,k) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            else if (Vtype(i,j-1,k)<=0 .or. Vtype(i,j,k)<=0) then

              flux = ( ( sum(mask2 * V(i,j-1:j,k)) )**2 ) / dymin
              V2(i,j-1,k) = V2(i,j-1,k) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do

      do l=1,2
        !$omp do
        do k=l,Vnz+1,2
          do j=1,Vny
            do i=1,Vnx
              if (Vtype(i,j,k-1)<=0 .and. Vtype(i,j,k)<=0) then

                flux = ( sum(mask4 * V(i,j,k-2:k+1)) * sum(mask2 * W(i,j:j+1,k-1)) ) / dzmin
                V2(i,j,k-1) = V2(i,j,k-1) - flux
                V2(i,j,k)   = V2(i,j,k)   + flux

              else if (Vtype(i,j,k-1)<=0 .or. Vtype(i,j,k)<=0) then

                flux = ( sum(mask2 * V(i,j,k-1:k)) * sum(mask2 * W(i,j:j+1,k-1)) ) / dzmin
                V2(i,j,k-1) = V2(i,j,k-1) - flux
                V2(i,j,k)   = V2(i,j,k)   + flux

              end if
            end do
          end do
        end do
        !$omp end do
      end do

      !$omp workshare
      V2 = V2 * dt / divcoef
      !$omp end workshare
      !$omp end parallel

    end subroutine CDS4_2V


    subroutine CDS4_2W(W2,U,V,W)
      real(knd),dimension(-2:,-2:,-2:),intent(out) :: W2
      real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: mask4(4) = (/ -1, 9, 9, -1 /)
      integer,parameter :: divcoef = 256
      integer,parameter :: mask2(2) = (/  8, 8 /)
      integer   :: i,j,k,l
      real(knd) :: flux

      !$omp parallel private(i,j,k,l,flux)
      !$omp workshare
      W2 = 0
      !$omp end workshare

      !$omp do
      do k=1,Wnz
        do j=1,Wny
          do i=1,Wnx+1
            if (Wtype(i-1,j,k)<=0 .and. Wtype(i,j,k)<=0) then

              flux = ( sum(mask4 * W(i-2:i+1,j,k)) * sum(mask2 * U(i-1,j,k:k+1)) ) / dxmin
              W2(i-1,j,k) = W2(i-1,j,k) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            else if (Wtype(i-1,j,k)<=0 .or. Wtype(i,j,k)<=0) then

              flux = ( sum(mask2 * W(i-1:i,j,k)) * sum(mask2 * U(i-1,j,k:k+1)) ) / dxmin
              W2(i-1,j,k) = W2(i-1,j,k) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do

      !$omp do
      do k=1,Wnz
        do j=1,Wny+1
          do i=1,Wnx
            if (Wtype(i,j-1,k)<=0 .and. Wtype(i,j,k)<=0) then

              flux = ( sum(mask4 * W(i,j-2:j+1,k)) * sum(mask2 * V(i,j-1,k:k+1)) ) / dymin
              W2(i,j-1,k) = W2(i,j-1,k) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            else if (Wtype(i,j-1,k)<=0 .or. Wtype(i,j,k)<=0) then

              flux = ( sum(mask2 * W(i,j-1:j,k)) * sum(mask2 * V(i,j-1,k:k+1)) ) / dymin
              W2(i,j-1,k) = W2(i,j-1,k) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do

      do l=1,2
        !$omp do
        do k=l,Wnz+1,2
          do j=1,Wny
            do i=1,Wnx
              if (Wtype(i,j,k-1)<=0 .and. Wtype(i,j,k)<=0) then

                flux = ( ( sum(mask4 * W(i,j,k-2:k+1)) )**2 ) / dzmin
                W2(i,j,k-1) = W2(i,j,k-1) - flux
                W2(i,j,k)   = W2(i,j,k)   + flux

              else if (Wtype(i,j,k-1)<=0 .or. Wtype(i,j,k)<=0) then

                flux = ( ( sum(mask2 * W(i,j,k-1:k)) )**2 ) / dzmin
                W2(i,j,k-1) = W2(i,j,k-1) - flux
                W2(i,j,k)   = W2(i,j,k)   + flux

              end if
            end do
          end do
        end do
        !$omp end do
      end do

      !$omp workshare
      W2 = W2 * dt / divcoef
      !$omp end workshare
      !$omp end parallel

    end subroutine CDS4_2W













    subroutine CD4divU(U2,U,V,W)
      !Morinishi et al., JCP 143, http://dx.doi.org/10.1006/jcph.1998.5962
      real(knd),dimension(-2:,-2:,-2:),intent(out) :: U2
      real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: i_mask4(4) = [ -1, 9, 9, -1 ]
      integer,parameter :: divcoef = 256
      integer,parameter :: coef2ord = divcoef / 4
      integer,parameter :: narr = 4
      !UV1 and UV3 are sized to fit to the L1 cache, so no stack overflow should be possible
      real(knd) :: UV1(-1:tilenx(narr)+2,-1:tileny(narr)+2,-1:tilenz(narr)+2)
      real(knd) :: UV3(-1:tilenx(narr)+2,-1:tileny(narr)+2,-1:tilenz(narr)+2)
      integer   :: bi,bj,bk,i,j,k,li,lj,lk
      real(knd) :: Uint,Vint,Wint,dU

      call set(U2,0._knd)

      !$omp parallel private(bi,bj,bk,i,j,k,li,lj,lk,Uint,Vint,Wint,dU,UV1,UV3)
      !$omp do
      do bk = 1,Unz,tilenz(narr)
       do bj = 1,Uny,tileny(narr)
        do bi = 1,Unx,tilenx(narr)

          do k = bk,min(bk+tilenz(narr)-1,Unz)
           do j = bj,min(bj+tileny(narr)-1,Uny)
            do i = bi-1,min(bi+tilenx(narr)-1,Unx)+2
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 Uj_1 - 1/8 Uj_3
              Uint = sum( i_mask4 * U(i-2:i+1,j,k) )
              UV1(li,lj,lk) = Uint * (U(i-1,j,k)+U(i  ,j,k))
              UV3(li,lj,lk) = Uint * (U(i-2,j,k)+U(i+1,j,k))
            end do
           end do
          end do

          do k = bk,min(bk+tilenz(narr)-1,Unz)
           do j = bj,min(bj+tileny(narr)-1,Uny)
            do i = bi,min(bi+tilenx(narr)-1,Unx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 d1 - 1/8 d3
              if (Utype(i,j,k)==0) then
                dU =     9 * (UV1(li+1,lj,lk) - UV1(li  ,lj,lk))
                dU = dU -    (UV3(li+2,lj,lk) - UV3(li-1,lj,lk))/3
              else if (Utype(i,j,k)<0) then  !near a boundary - 2nd order
                dU = coef2ord * ((U(i+1,j,k)+U(i  ,j,k)) * (U(i+1,j,k)+U(i  ,j,k)) &
                                -(U(i  ,j,k)+U(i-1,j,k)) * (U(i  ,j,k)+U(i-1,j,k)))
              end if
              U2(i,j,k) = U2(i,j,k) + dU / dxmin
            end do
           end do
          end do

        end do
       end do
      end do
      !$omp end do

      !$omp do
      do bk = 1,Unz,tilenz(narr)
       do bj = 1,Uny,tileny(narr)
        do bi = 1,Unx,tilenx(narr)

          do k = bk,min(bk+tilenz(narr)-1,Unz)
           do j = bj-2,min(bj+tileny(narr)-1,Uny)+1
            do i = bi,min(bi+tilenx(narr)-1,Unx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 Uj_1 - 1/8 Uj_3
              Vint = sum( i_mask4 * V(i-1:i+2,j,k) )
              UV1(li,lj,lk) = Vint * (U(i,j  ,k)+U(i,j+1,k))
              UV3(li,lj,lk) = Vint * (U(i,j-1,k)+U(i,j+2,k))
            end do
           end do
          end do

          do k = bk,min(bk+tilenz(narr)-1,Unz)
           do j = bj,min(bj+tileny(narr)-1,Uny)
            do i = bi,min(bi+tilenx(narr)-1,Unx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 d1 - 1/8 d3
              if (Utype(i,j,k)==0) then
                dU =     9 * (UV1(li,lj  ,lk) - UV1(li,lj-1,lk))
                dU = dU -    (UV3(li,lj+1,lk) - UV3(li,lj-2,lk))/3
              else if (Utype(i,j,k)<0) then  !near a boundary - 2nd order
                dU = coef2ord * ((U(i,j+1,k)+U(i,j  ,k)) * (V(i+1,j  ,k)+V(i,j  ,k)) &
                                -(U(i,j  ,k)+U(i,j-1,k)) * (V(i+1,j-1,k)+V(i,j-1,k)))
              end if
              U2(i,j,k) = U2(i,j,k) + dU / dymin
            end do
           end do
          end do

        end do
       end do
      end do
      !$omp end do

      !$omp do
      do bk = 1,Unz,tilenz(narr)
       do bj = 1,Uny,tileny(narr)
        do bi = 1,Unx,tilenx(narr)

          do k = bk-2,min(bk+tilenz(narr),Unz)+1
           do j = bj,min(bj+tileny(narr)-1,Uny)
            do i = bi,min(bi+tilenx(narr)-1,Unx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 Uj_1 - 1/8 Uj_3
              Wint = sum( i_mask4 * W(i-1:i+2,j,k) )
              UV1(li,lj,lk) = Wint * (U(i,j,k  )+U(i,j,k+1))
              UV3(li,lj,lk) = Wint * (U(i,j,k-1)+U(i,j,k+2))
            end do
           end do
          end do

          do k = bk,min(bk+tilenz(narr)-1,Unz)
           do j = bj,min(bj+tileny(narr)-1,Uny)
            do i = bi,min(bi+tilenx(narr)-1,Unx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 d1 - 1/8 d3
              if (Utype(i,j,k)==0) then
                dU =     9 * (UV1(li,lj,lk  ) - UV1(li,lj,lk-1))
                dU = dU -    (UV3(li,lj,lk+1) - UV3(li,lj,lk-2))/3
              else if (Utype(i,j,k)<0) then  !near a boundary - 2nd order
                dU = coef2ord * ((U(i,j,k+1)+U(i,j,k  )) * (W(i+1,j,k  )+W(i,j,k  )) &
                                -(U(i,j,k  )+U(i,j,k-1)) * (W(i+1,j,k-1)+W(i,j,k-1)))
              end if
              U2(i,j,k) = U2(i,j,k) + dU / dzmin
            end do
           end do
          end do

        end do
       end do
      end do
      !$omp end do
      !$omp end parallel

      call multiply(U2,-dt/divcoef)

    end subroutine CD4divU











    subroutine CD4divV(V2,U,V,W)
      !Morinishi et al., JCP 143, http://dx.doi.org/10.1006/jcph.1998.5962
      real(knd),dimension(-2:,-2:,-2:),intent(out) :: V2
      real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: i_mask4(4) = [ -1, 9, 9, -1 ]
      integer,parameter :: divcoef = 256
      integer,parameter :: coef2ord = divcoef / 4
      integer,parameter :: narr = 4
      !UV1 and UV3 are sized to fit to the L1 cache, so no stack overflow should be possible
      real(knd) :: UV1(-1:tilenx(narr)+2,-1:tileny(narr)+2,-1:tilenz(narr)+2)
      real(knd) :: UV3(-1:tilenx(narr)+2,-1:tileny(narr)+2,-1:tilenz(narr)+2)
      integer   :: bi,bj,bk,i,j,k,li,lj,lk
      real(knd) :: Uint,Vint,Wint,dV

      call set(V2,0._knd)

      !$omp parallel private(bi,bj,bk,i,j,k,li,lj,lk,Uint,Vint,Wint,dV,UV1,UV3)
      !$omp do
      do bk = 1,Vnz,tilenz(narr)
       do bj = 1,Vny,tileny(narr)
        do bi = 1,Vnx,tilenx(narr)

          do k = bk,min(bk+tilenz(narr)-1,Vnz)
           do j = bj,min(bj+tileny(narr)-1,Vny)
            do i = bi-2,min(bi+tilenx(narr)-1,Vnx)+1
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 Uj_1 - 1/8 Uj_3
              Uint = sum( i_mask4 * U(i,j-1:j+2,k) )
              UV1(li,lj,lk) = Uint * (V(i  ,j,k)+V(i+1,j,k))
              UV3(li,lj,lk) = Uint * (V(i-1,j,k)+V(i+2,j,k))
            end do
           end do
          end do

          do k = bk,min(bk+tilenz(narr)-1,Vnz)
           do j = bj,min(bj+tileny(narr)-1,Vny)
            do i = bi,min(bi+tilenx(narr)-1,Vnx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 d1 - 1/8 d3
              if (Vtype(i,j,k)==0) then
                dV =     9 * (UV1(li  ,lj,lk) - UV1(li-1,lj,lk))
                dV = dV -    (UV3(li+1,lj,lk) - UV3(li-2,lj,lk))/3
              else if (Vtype(i,j,k)<0) then  !near a boundary - 2nd order
                dV = coef2ord * ((V(i+1,j,k)+V(i  ,j,k)) * (U(i  ,j+1,k)+U(i  ,j,k)) &
                                -(V(i  ,j,k)+V(i-1,j,k))  *(U(i-1,j+1,k)+U(i-1,j,k)))
              end if
              V2(i,j,k) = V2(i,j,k) + dV / dxmin
            end do
           end do
          end do

        end do
       end do
      end do
      !$omp end do

      !$omp do
      do bk = 1,Vnz,tilenz(narr)
       do bj = 1,Vny,tileny(narr)
        do bi = 1,Vnx,tilenx(narr)

          do k = bk,min(bk+tilenz(narr)-1,Vnz)
           do j = bj-1,min(bj+tileny(narr)-1,Vny)+2
            do i = bi,min(bi+tilenx(narr)-1,Vnx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 Uj_1 - 1/8 Uj_3
              Vint = sum( i_mask4 * V(i,j-2:j+1,k) )
              UV1(li,lj,lk) = Vint * (V(i,j-1,k)+V(i,j  ,k))
              UV3(li,lj,lk) = Vint * (V(i,j-2,k)+V(i,j+1,k))
            end do
           end do
          end do

          do k = bk,min(bk+tilenz(narr)-1,Vnz)
           do j = bj,min(bj+tileny(narr)-1,Vny)
            do i = bi,min(bi+tilenx(narr)-1,Vnx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 d1 - 1/8 d3
              if (Vtype(i,j,k)==0) then
                dV =     9 * (UV1(li,lj+1,lk) - UV1(li,lj  ,lk))
                dV = dV -    (UV3(li,lj+2,lk) - UV3(li,lj-1,lk))/3
              else if (Vtype(i,j,k)<0) then  !near a boundary - 2nd order
                dV = coef2ord * ((V(i,j+1,k)+V(i,j  ,k)) * (V(i,j+1,k)+V(i,j  ,k)) &
                                -(V(i,j  ,k)+V(i,j-1,k)) * (V(i,j  ,k)+V(i,j-1,k)))
              end if
              V2(i,j,k) = V2(i,j,k) + dV / dymin
            end do
           end do
          end do

        end do
       end do
      end do
      !$omp end do

      !$omp do
      do bk = 1,Vnz,tilenz(narr)
       do bj = 1,Vny,tileny(narr)
        do bi = 1,Vnx,tilenx(narr)

          do k = bk-2,min(bk+tilenz(narr),Vnz)+1
           do j = bj,min(bj+tileny(narr)-1,Vny)
            do i = bi,min(bi+tilenx(narr)-1,Vnx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 Uj_1 - 1/8 Uj_3
              Wint = sum( i_mask4 * W(i,j-1:j+2,k) )
              UV1(li,lj,lk) = Wint * (V(i,j,k  )+V(i,j,k+1))
              UV3(li,lj,lk) = Wint * (V(i,j,k-1)+V(i,j,k+2))
            end do
           end do
          end do

          do k = bk,min(bk+tilenz(narr)-1,Vnz)
           do j = bj,min(bj+tileny(narr)-1,Vny)
            do i = bi,min(bi+tilenx(narr)-1,Vnx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 d1 - 1/8 d3
              if (Vtype(i,j,k)==0) then
                dV =     9 * (UV1(li,lj,lk  ) - UV1(li,lj,lk-1))
                dV = dV -    (UV3(li,lj,lk+1) - UV3(li,lj,lk-2))/3
              else if (Vtype(i,j,k)<0) then  !near a boundary - 2nd order
                dV = coef2ord * ((V(i,j,k+1)+V(i,j,k  )) * (W(i,j+1,k  )+W(i,j,k  )) &
                                -(V(i,j,k  )+V(i,j,k-1)) * (W(i,j+1,k-1)+W(i,j,k-1)))
              end if
              V2(i,j,k) = V2(i,j,k) + dV / dzmin
            end do
           end do
          end do

        end do
       end do
      end do
      !$omp end do
      !$omp end parallel

      call multiply(V2,-dt/divcoef)

    end subroutine CD4divV
















    subroutine CD4divW(W2,U,V,W)
      !Morinishi et al., JCP 143, http://dx.doi.org/10.1006/jcph.1998.5962
      real(knd),dimension(-2:,-2:,-2:),intent(out) :: W2
      real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: i_mask4(4) = [ -1, 9, 9, -1 ]
      integer,parameter :: divcoef = 256
      integer,parameter :: coef2ord = divcoef / 4
      integer,parameter :: narr = 4
      !UV1 and UV3 are sized to fit to the L1 cache, so no stack overflow should be possible
      real(knd) :: UV1(-1:tilenx(narr)+2,-1:tileny(narr)+2,-1:tilenz(narr)+2)
      real(knd) :: UV3(-1:tilenx(narr)+2,-1:tileny(narr)+2,-1:tilenz(narr)+2)
      integer   :: bi,bj,bk,i,j,k,li,lj,lk
      real(knd) :: Uint,Vint,Wint,dW

      call set(W2,0._knd)

      !$omp parallel private(bi,bj,bk,i,j,k,li,lj,lk,Uint,Vint,Wint,dW,UV1,UV3)
      !$omp do
      do bk = 1,Wnz,tilenz(narr)
       do bj = 1,Wny,tileny(narr)
        do bi = 1,Wnx,tilenx(narr)

          do k = bk,min(bk+tilenz(narr)-1,Wnz)
           do j = bj,min(bj+tileny(narr)-1,Wny)
            do i = bi-2,min(bi+tilenx(narr)-1,Wnx)+1
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 Uj_1 - 1/8 Uj_3
              Uint = sum( i_mask4 * U(i,j,k-1:k+2) )
              UV1(li,lj,lk) = Uint * (W(i  ,j,k)+W(i+1,j,k))
              UV3(li,lj,lk) = Uint * (W(i-1,j,k)+W(i+2,j,k))
            end do
           end do
          end do

          do k = bk,min(bk+tilenz(narr)-1,Wnz)
           do j = bj,min(bj+tileny(narr)-1,Wny)
            do i = bi,min(bi+tilenx(narr)-1,Wnx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 d1 - 1/8 d3
              if (Wtype(i,j,k)==0) then
                dW =     9 * (UV1(li  ,lj,lk) - UV1(li-1,lj,lk))
                dW = dW -    (UV3(li+1,lj,lk) - UV3(li-2,lj,lk))/3
              else if (Wtype(i,j,k)<0) then  !near a boundary - 2nd order
                dW = coef2ord * ((W(i+1,j,k)+W(i  ,j,k)) * (U(i  ,j,k+1)+U(i  ,j,k)) &
                                -(W(i  ,j,k)+W(i-1,j,k)) * (U(i-1,j,k+1)+U(i-1,j,k)))
              end if
              W2(i,j,k) = W2(i,j,k) + dW / dxmin
            end do
           end do
          end do

        end do
       end do
      end do
      !$omp end do

      !$omp do
      do bk = 1,Wnz,tilenz(narr)
       do bj = 1,Wny,tileny(narr)
        do bi = 1,Wnx,tilenx(narr)

          do k = bk,min(bk+tilenz(narr)-1,Wnz)
           do j = bj-2,min(bj+tileny(narr)-1,Wny)+1
            do i = bi,min(bi+tilenx(narr)-1,Wnx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 Uj_1 - 1/8 Uj_3
              Vint = sum( i_mask4 * V(i,j,k-1:k+2) )
              UV1(li,lj,lk) = Vint * (W(i,j  ,k)+W(i,j+1,k))
              UV3(li,lj,lk) = Vint * (W(i,j-1,k)+W(i,j+2,k))
            end do
           end do
          end do

          do k = bk,min(bk+tilenz(narr)-1,Wnz)
           do j = bj,min(bj+tileny(narr)-1,Wny)
            do i = bi,min(bi+tilenx(narr)-1,Wnx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 d1 - 1/8 d3
              if (Wtype(i,j,k)==0) then
                dW =     9 * (UV1(li,lj  ,lk) - UV1(li,lj-1,lk))
                dW = dW -    (UV3(li,lj+1,lk) - UV3(li,lj-2,lk))/3
              else if (Wtype(i,j,k)<0) then  !near a boundary - 2nd order
                dW = coef2ord * ((W(i,j+1,k)+W(i,j  ,k)) * (V(i,j,k+1)+V(i,j  ,k  )) &
                                -(W(i,j  ,k)+W(i,j-1,k)) * (V(i,j-1,k)+V(i,j-1,k+1)))
              end if
              W2(i,j,k) = W2(i,j,k) + dW / dymin
            end do
           end do
          end do

        end do
       end do
      end do
      !$omp end do

      !$omp do
      do bk = 1,Wnz,tilenz(narr)
       do bj = 1,Wny,tileny(narr)
        do bi = 1,Wnx,tilenx(narr)

          do k = bk-1,min(bk+tilenz(narr)-1,Wnz)+2
           do j = bj,min(bj+tileny(narr)-1,Wny)
            do i = bi,min(bi+tilenx(narr)-1,Wnx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 Uj_1 - 1/8 Uj_3
              Wint = sum( i_mask4 * W(i,j,k-2:k+1) )
              UV1(li,lj,lk) = Wint * (W(i,j,k-1)+W(i,j,k  ))
              UV3(li,lj,lk) = Wint * (W(i,j,k-2)+W(i,j,k+1))
            end do
           end do
          end do

          do k = bk,min(bk+tilenz(narr)-1,Wnz)
           do j = bj,min(bj+tileny(narr)-1,Wny)
            do i = bi,min(bi+tilenx(narr)-1,Wnx)
              li = i-bi+1
              lj = j-bj+1
              lk = k-bk+1
              !9/8 d1 - 1/8 d3
              if (Wtype(i,j,k)==0) then
                dW =     9 * (UV1(li,lj,lk+1) - UV1(li,lj,lk  ))
                dW = dW -    (UV3(li,lj,lk+2) - UV3(li,lj,lk-1))/3
              else if (Wtype(i,j,k)<0) then  !near a boundary - 2nd order
                dW = coef2ord * ((W(i,j,k+1)+W(i,j,k  )) * (W(i,j,k+1)+W(i,j,k  )) &
                                -(W(i,j,k  )+W(i,j,k-1)) * (W(i,j,k  )+W(i,j,k-1)))
              end if
              W2(i,j,k) = W2(i,j,k) + dW / dzmin
            end do
           end do
          end do

        end do
       end do
      end do
      !$omp end do
      !$omp end parallel

      call multiply(W2,-dt/divcoef)

    end subroutine CD4divW




























!     subroutine JST4U(U2,U,V,W)
!       real(knd),dimension(-2:,-2:,-2:),intent(out) :: U2
!       real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
!       integer,parameter :: mask4(4) = (/ -1, 9, 9, -1 /)
!       integer,parameter :: divcoef = 16
!       integer,parameter :: mask2(4) = (/ 0, 8, 8, 0 /)
!       integer,parameter :: diffmask(4) = (/ -1, +3, -3, +1 /)
!       real(knd),parameter :: k4 = 0.01_knd
!       integer   :: i,j,k,l,mask(4)
!       real(knd) :: flux,nubar,eps4,Uadv,Vadv,Wadv,d
!       real(knd) :: nu(-2:ubound(U,1),-2:ubound(U,2),-2:ubound(U,3))
!
!       !$omp parallel private(i,j,k,l,flux,nubar,eps4,Uadv,Vadv,Wadv,mask)
!       !$omp workshare
!       U2 = 0
!       !$omp end workshare
!
!       !$omp do
!       do k=1,Unz
!         do j=1,Uny
!           do i=-1,Unx+2
!             d = sum(U(i-1:i+1,j,k)*[1,2,1])  !divisor
!             if (abs(d) > 100*tiny(1._knd)) then
!               nu(i,j,k) = abs( sum(U(i-1:i+1,j,k)*[1,-2,1]) / d )
!             else
!               nu(i,j,k) = 0
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       !$omp do
!       do k=1,Unz
!         do j=1,Uny
!           do i=1,Unx+1
!             if (Utype(i-1,j,k)<=0 .or. Utype(i,j,k)<=0) then
!
!               nubar = maxval(nu(i-2:i+1,j,k))
!               eps4 = 0.1!max(0._knd, k4 * max(nubar,0.1_knd))
!
!               if (Utype(i-1,j,k)<=0 .and. Utype(i,j,k)<=0) then
!                 mask = mask4
!               else
!                 mask = mask2
!               end if
!
!               Uadv = sum(mask * U(i-2:i+1,j,k))
!               flux = eps4 * Uadv * sum(U(i-2:i+1,j,k) * diffmask)
!
!               flux = flux + ( ( Uadv )**2 ) / divcoef
!               flux = flux / dxmin
!               U2(i-1,j,k) = U2(i-1,j,k) - flux
!               U2(i,j,k)   = U2(i,j,k)   + flux
!
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       !$omp do
!       do k=1,Unz
!         do j=-1,Uny+2
!           do i=1,Unx
!             d = sum(U(i,j-1:j+1,k)*[1,2,1])  !divisor
!             if (abs(d) > 100*tiny(1._knd)) then
!               nu(i,j,k) = abs( sum(U(i,j-1:j+1,k)*[1,-2,1]) / d )
!             else
!               nu(i,j,k) = 0
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       !$omp do
!       do k=1,Unz
!         do j=1,Uny+1
!           do i=1,Unx
!             if (Utype(i,j-1,k)<=0 .or. Utype(i,j,k)<=0) then
!
!               nubar = maxval(nu(i,j-2:j+1,k))
!               eps4 = 0.1!max(0._knd, k4 * max(nubar,.1_knd))
!
!               if (Utype(i,j-1,k)<=0 .and. Utype(i,j,k)<=0) then
!                 mask = mask4
!               else
!                 mask = mask2
!               end if
!
!               Vadv = sum(mask * V(i-1:i+2,j-1,k))
!               flux = eps4 * Vadv * sum(U(i-2:i+1,j,k) * diffmask)
!
!               flux = flux + ( sum(mask * U(i,j-2:j+1,k)) * Vadv ) / divcoef
!               flux = flux / dymin
!               U2(i,j-1,k) = U2(i,j-1,k) - flux
!               U2(i,j,k)   = U2(i,j,k)   + flux
!
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       !$omp do
!       do k=-1,Unz+2
!         do j=1,Uny
!           do i=1,Unx
!             d = sum(U(i,j,k-1:k+1)*[1,2,1])  !divisor
!             if (abs(d) > 100*tiny(1._knd)) then
!               nu(i,j,k) = abs( sum(U(i,j,k-1:k+1)*[1,-2,1]) / d )
!             else
!               nu(i,j,k) = 0
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       do l=1,2
!         !$omp do
!         do k=l,Unz+1,2
!           do j=1,Uny
!             do i=1,Unx
!               if (Utype(i,j,k-1)<=0 .or. Utype(i,j,k)<=0) then
!
!                 nubar = maxval(nu(i,j,k-2:k+1))
!                 eps4 = 0.1!max(0._knd, k4 * max(nubar,.1_knd))
!
!                 if (Utype(i,j,k-1)<=0 .and. Utype(i,j,k)<=0) then
!                   mask = mask4
!                 else
!                   mask = mask2
!                 end if
!
!                 Wadv = sum(mask * W(i-1:i+2,j,k-1))
!                 flux = eps4 * Wadv * sum(U(i-2:i+1,j,k) * diffmask)
!
!                 flux = flux + ( sum(mask * U(i,j,k-2:k+1)) * Wadv ) / divcoef
!                 flux = flux / dzmin
!                 U2(i,j,k-1) = U2(i,j,k-1) - flux
!                 U2(i,j,k)   = U2(i,j,k)   + flux
!
!               end if
!             end do
!           end do
!         end do
!         !$omp end do
!       end do
!
!       !$omp workshare
!       U2 = U2 * dt / divcoef
!       !$omp end workshare
!       !$omp end parallel
!
!     end subroutine JST4U
!
!
!     subroutine JST4V(V2,U,V,W)
!       real(knd),dimension(-2:,-2:,-2:),intent(out) :: V2
!       real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
!       integer,parameter :: mask4(4) = (/ -1, 9, 9, -1 /)
!       integer,parameter :: divcoef = 16
!       integer,parameter :: mask2(4) = (/ 0, 8, 8, 0 /)
!       integer,parameter :: diffmask(4) = (/ -1, +3, -3, +1 /)
!       real(knd),parameter :: k4 = 0.01_knd
!       integer   :: i,j,k,l,mask(4)
!       real(knd) :: flux,nubar,eps4,Uadv,Vadv,Wadv,d
!       real(knd) :: nu(-2:ubound(V,1),-2:ubound(V,2),-2:ubound(V,3))
!
!       !$omp parallel private(i,j,k,l,flux,nubar,eps4,Uadv,Vadv,Wadv,mask)
!       !$omp workshare
!       V2 = 0
!       !$omp end workshare
!
!
!       !$omp do
!       do k=1,Vnz
!         do j=1,Vny
!           do i=-1,Vnx+2
!             d = sum(V(i-1:i+1,j,k)*[1,2,1])  !divisor
!             if (abs(d) > 100*tiny(1._knd)) then
!               nu(i,j,k) = abs( sum(V(i-1:i+1,j,k)*[1,-2,1]) / d )
!             else
!               nu(i,j,k) = 0
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!      !$omp do
!       do k=1,Vnz
!         do j=1,Vny
!           do i=1,Vnx+1
!             if (Vtype(i-1,j,k)<=0 .or. Vtype(i,j,k)<=0) then
!
!               nubar = maxval(nu(i-2:i+1,j,k))
!               eps4 = 0.1!max(0._knd, k4 * max(nubar,.1_knd))
!
!               if (Vtype(i-1,j,k)<=0 .and. Vtype(i,j,k)<=0) then
!                 mask = mask4
!               else
!                 mask = mask2
!               end if
!
!               Uadv = sum(mask * U(i-1,j-1:j+2,k))
!               flux = eps4 * Uadv * sum(V(i-2:i+1,j,k) * diffmask)
!
!               flux = flux + ( sum(mask * V(i-2:i+1,j,k)) * Uadv ) / divcoef
!               flux = flux / dxmin
!               V2(i-1,j,k) = V2(i-1,j,k) - flux
!               V2(i,j,k)   = V2(i,j,k)   + flux
!
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       !$omp do
!       do k=1,Vnz
!         do j=-1,Vny+2
!           do i=1,Vnx
!             d = sum(V(i,j-1:j+1,k)*[1,2,1])  !divisor
!             if (abs(d) > 100*tiny(1._knd)) then
!               nu(i,j,k) = abs( sum(V(i,j-1:j+1,k)*[1,-2,1]) / d )
!             else
!               nu(i,j,k) = 0
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       !$omp do
!       do k=1,Vnz
!         do j=1,Vny+1
!           do i=1,Vnx
!             if (Vtype(i,j-1,k)<=0 .or. Vtype(i,j,k)<=0) then
!
!               nubar = maxval(nu(i,j-2:j+1,k))
!               eps4 = 0.1!max(0._knd, k4 * max(nubar,.1_knd))
!
!               if (Vtype(i,j-1,k)<=0 .and. Vtype(i,j,k)<=0) then
!                 mask = mask4
!               else
!                 mask = mask2
!               end if
!
!               Vadv = sum(mask * V(i,j-2:j+1,k))
!               flux = eps4 * Vadv * sum(V(i,j-2:j+1,k) * diffmask)
!
!               flux = flux + ( ( Vadv )**2 ) / divcoef
!               flux = flux / dymin
!               V2(i,j-1,k) = V2(i,j-1,k) - flux
!               V2(i,j,k)   = V2(i,j,k)   + flux
!
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       !$omp do
!       do k=-1,Vnz+2
!         do j=1,Vny
!           do i=1,Vnx
!             d = sum(V(i,j,k-1:k+1)*[1,2,1])  !divisor
!             if (abs(d) > 100*tiny(1._knd)) then
!               nu(i,j,k) = abs( sum(V(i,j,k-1:k+1)*[1,-2,1]) / d )
!             else
!               nu(i,j,k) = 0
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       do l=1,2
!         !$omp do
!         do k=l,Vnz+1,2
!           do j=1,Vny
!             do i=1,Vnx
!               if (Vtype(i,j,k-1)<=0 .or. Vtype(i,j,k)<=0) then
!
!                 nubar = maxval(nu(i,j,k-2:k+1))
!                 eps4 = 0.1!max(0._knd, k4 * max(nubar,.1_knd))
!
!                 if (Vtype(i,j,k-1)<=0 .and. Vtype(i,j,k)<=0) then
!                   mask = mask4
!                 else
!                   mask = mask2
!                 end if
!
!                 Wadv = sum(mask * W(i,j-1:j+2,k-1))
!                 flux = eps4 * Wadv * sum(V(i,j,k-2:k+1) * diffmask)
!
!                 flux = flux + ( sum(mask * V(i,j,k-2:k+1)) * Wadv ) / divcoef
!                 flux = flux / dzmin
!                 V2(i,j,k-1) = V2(i,j,k-1) - flux
!                 V2(i,j,k)   = V2(i,j,k)   + flux
!
!               end if
!             end do
!           end do
!         end do
!         !$omp end do
!       end do
!
!       !$omp workshare
!       V2 = V2 * dt / divcoef
!       !$omp end workshare
!       !$omp end parallel
!
!     end subroutine JST4V
!
!
!     subroutine JST4W(W2,U,V,W)
!       real(knd),dimension(-2:,-2:,-2:),intent(out) :: W2
!       real(knd),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
!       integer,parameter :: mask4(4) = (/ -1, 9, 9, -1 /)
!       integer,parameter :: divcoef = 16
!       integer,parameter :: mask2(4) = (/ 0, 8, 8, 0 /)
!       integer,parameter :: diffmask(4) = (/ -1, +3, -3, +1 /)
!       real(knd),parameter :: k4 = 0.01_knd
!       integer   :: i,j,k,l,mask(4)
!       real(knd) :: flux,nubar,eps4,Uadv,Vadv,Wadv,d
!       real(knd) :: nu(-2:ubound(U,1),-2:ubound(U,2),-2:ubound(U,3))
!
!       !$omp parallel private(i,j,k,l,flux,nubar,eps4,Uadv,Vadv,Wadv,mask)
!       !$omp workshare
!       W2 = 0
!       !$omp end workshare
!
!       !$omp do
!       do k=1,Wnz
!         do j=1,Wny
!           do i=-1,Wnx+2
!             d = sum(W(i-1:i+1,j,k)*[1,2,1])  !divisor
!             if (abs(d) > 100*tiny(1._knd)) then
!               nu(i,j,k) = abs( sum(W(i-1:i+1,j,k)*[1,-2,1]) / d )
!             else
!               nu(i,j,k) = 0
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       !$omp do
!       do k=1,Wnz
!         do j=1,Wny
!           do i=1,Wnx+1
!             if (Wtype(i-1,j,k)<=0 .or. Wtype(i,j,k)<=0) then
!
!               nubar = maxval(nu(i-2:i+1,j,k))
!               eps4 = 0.1!max(0._knd, k4 * max(nubar,.1_knd))
!
!               if (Wtype(i-1,j,k)<=0 .and. Wtype(i,j,k)<=0) then
!                 mask = mask4
!               else
!                 mask = mask2
!               end if
!
!               Uadv = sum(mask * U(i-1,j,k-1:k+2))
!               flux = eps4 * Uadv * sum(W(i-2:i+1,j,k) * diffmask)
!
!               flux = flux + ( sum(mask * W(i-2:i+1,j,k)) * Uadv ) / divcoef
!               flux = flux /dxmin
!               W2(i-1,j,k) = W2(i-1,j,k) - flux
!               W2(i,j,k)   = W2(i,j,k)   + flux
!
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       !$omp do
!       do k=1,Wnz
!         do j=-1,Wny+2
!           do i=1,Wnx
!             d = sum(W(i,j-1:j+1,k)*[1,2,1])  !divisor
!             if (abs(d) > 100*tiny(1._knd)) then
!               nu(i,j,k) = abs( sum(W(i,j-1:j+1,k)*[1,-2,1]) / d )
!             else
!               nu(i,j,k) = 0
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       !$omp do
!       do k=1,Wnz
!         do j=1,Wny+1
!           do i=1,Wnx
!             if (Wtype(i,j-1,k)<=0 .or. Wtype(i,j,k)<=0) then
!
!               nubar = maxval(nu(i,j-2:j+1,k))
!               eps4 = 0.1!max(0._knd, k4 * max(nubar,.1_knd))
!
!               if (Wtype(i,j-1,k)<=0 .and. Wtype(i,j,k)<=0) then
!                 mask = mask4
!               else
!                 mask = mask2
!               end if
!
!               Vadv = sum(mask * V(i,j-1,k-1:k+2))
!               flux = eps4 * Vadv * sum(W(i,j-2:j+1,k) * diffmask)
!
!               flux = flux + ( sum(mask * W(i,j-2:j+1,k)) * Vadv ) / divcoef
!               flux = flux / dymin
!               W2(i,j-1,k) = W2(i,j-1,k) - flux
!               W2(i,j,k)   = W2(i,j,k)   + flux
!
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       !$omp do
!       do k=-1,Wnz+2
!         do j=1,Wny
!           do i=1,Wnx
!             d = sum(W(i,j,k-1:k+1)*[1,2,1])  !divisor
!             if (abs(d) > 100*tiny(1._knd)) then
!               nu(i,j,k) = abs( sum(W(i,j,k-1:k+1)*[1,-2,1]) / d )
!             else
!               nu(i,j,k) = 0
!             end if
!           end do
!         end do
!       end do
!       !$omp end do
!
!       do l=1,2
!         !$omp do
!         do k=l,Wnz+1,2
!           do j=1,Wny
!             do i=1,Wnx
!               if (Wtype(i,j,k-1)<=0 .or. Wtype(i,j,k)<=0) then
!
!                 nubar = maxval(nu(i,j,k-2:k+1))
!                 eps4 = 0.1!max(0._knd, k4 * max(nubar,.1_knd))
!
!                 if (Wtype(i,j,k-1)<=0 .and. Wtype(i,j,k)<=0) then
!                   mask = mask4
!                 else
!                   mask = mask2
!                 end if
!
!                 Wadv = sum(mask * W(i,j,k-2:k+1))
!                 flux = eps4 * Wadv * sum(W(i,j,k-2:k+1) * diffmask)
!
!
!                 flux = flux + ( ( Wadv )**2 ) / divcoef
!                 flux = flux / dzmin
!                 W2(i,j,k-1) = W2(i,j,k-1) - flux
!                 W2(i,j,k)   = W2(i,j,k)   + flux
!
!               end if
!             end do
!           end do
!         end do
!         !$omp end do
!       end do
!
!       !$omp workshare
!       W2 = W2 * dt / divcoef
!       !$omp end workshare
!       !$omp end parallel
!
!     end subroutine JST4W



end module CDS

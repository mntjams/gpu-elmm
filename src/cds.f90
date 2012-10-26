module CDS
    use PARAMETERS
    use BOUNDARIES
    use LIMITERS, only: FluxLimiter

    implicit none

    private
    public CDU, CDV, CDW, CDS4U, CDS4V, CDS4W, KAPPAU, KAPPAV, KAPPAW

    contains





    subroutine CDU(U2,U,V,W)
    real(KND) :: U2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer nx,ny,nz,i,j,k
    real(KND) Ax,Ay,Az,Utmp
    real(KND),allocatable,dimension(:),save:: fe,fp,gn,ht
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
                    U2(i,j,k)= - ((Ax*(U(i+1,j,k)+U(i,j,k))*(U(i+1,j,k)+U(i,j,k))&
                    -Ax*(U(i,j,k)+U(i-1,j,k))*(U(i,j,k)+U(i-1,j,k)))&
                    +(Ay*(U(i,j+1,k)+U(i,j,k))*(V(i+1,j,k)+V(i,j,k))&
                    -Ay*(U(i,j,k)+U(i,j-1,k))*(V(i+1,j-1,k)+V(i,j-1,k)))&
                    +(Az*(U(i,j,k+1)+U(i,j,k))*(W(i+1,j,k)+W(i,j,k))&
                    -Az*(U(i,j,k)+U(i,j,k-1))*(W(i+1,j,k-1)+W(i,j,k-1))))
                enddo
            enddo
        enddo
        !$omp end parallel do
     else
        if (called==0) then
         allocate(fe(0:nx),fp(1:nx),gn(0:ny),ht(0:nz))
         called=1
         forall (i=0:nx)      fe(i)=(xPr(i+1)-xU(i))/(xU(i+1)-xU(i))
         forall (i=1:nx)      fp(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
         forall (j=0:ny)      gn(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
         forall (k=0:nz)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
        endif

        !$omp parallel do private(i,j,k)
        do k=1,nz
         do j=1,ny
          do i=1,nx
           Utmp=( ((fe(i)*U(i+1,j,k)+(1-fe(i))*U(i,j,k))*(fe(i)*U(i+1,j,k)+(1-fe(i))*U(i,j,k)))&
                      -((fe(i-1)*U(i,j,k)+(1-fe(i-1))*U(i-1,j,k))*(fe(i-1)*U(i,j,k)+(1-fe(i-1))*U(i-1,j,k))))/dxU(i)
           Utmp=Utmp+( (gn(j)*U(i,j+1,k)+(1-gn(j))*U(i,j,k))*(fp(i)*V(i+1,j,k)+(1-fp(i))*V(i,j,k))&
                      -(gn(j-1)*U(i,j,k)+(1-gn(j-1))*U(i,j-1,k))*(fp(i)*V(i+1,j-1,k)+(1-fp(i))*V(i,j-1,k)))/dyPr(j)
           Utmp=Utmp+( (ht(k)*U(i,j,k+1)+(1-ht(k))*U(i,j,k))*(fp(i)*W(i+1,j,k)+(1-fp(i))*W(i,j,k))&
                      -(ht(k-1)*U(i,j,k)+(1-ht(k-1))*U(i,j,k-1))*(fp(i)*W(i+1,j,k-1)+(1-fp(i))*W(i,j,k-1)))/dzPr(k)
           U2(i,j,k)=-dt*Utmp
          enddo
         enddo
        enddo
        !$omp end parallel do
     endif
    end subroutine CDU






    subroutine CDV(V2,U,V,W)
    real(KND) :: V2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer nx,ny,nz,i,j,k
    real(KND) Ax,Ay,Az,Vtmp
    real(KND),allocatable,dimension(:),save:: gn,gp,fe,ht
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
                    V2(i,j,k)= - ((Ay*(V(i,j+1,k)+V(i,j,k))*(V(i,j+1,k)+V(i,j,k))&
                    -Ay*(V(i,j,k)+V(i,j-1,k))*(V(i,j,k)+V(i,j-1,k)))&
                    +(Ax*(V(i+1,j,k)+V(i,j,k))*(U(i,j+1,k)+U(i,j,k))&
                    -Ax*(V(i,j,k)+V(i-1,j,k))*(U(i-1,j+1,k)+U(i-1,j,k)))&
                    +(Az*(V(i,j,k+1)+V(i,j,k))*(W(i,j+1,k)+W(i,j,k))&
                    -Az*(V(i,j,k)+V(i,j,k-1))*(W(i,j+1,k-1)+W(i,j,k-1))))
                enddo
            enddo
        enddo
        !$omp end parallel do
     else
        if (called==0) then
         allocate(gn(0:ny),gp(1:ny),fe(0:nx),ht(0:nz))
         called=1
         forall (j=0:ny)      gn(j)=(yPr(j+1)-yV(j))/(yV(j+1)-yV(j))
         forall (j=1:ny)      gp(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
         forall (i=0:nx)      fe(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
         forall (k=0:nz)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
        endif

        !$omp parallel do private(i,j,k)
        do k=1,nz
         do j=1,ny
          do i=1,nx
           Vtmp=( ((gn(j)*V(i,j+1,k)+(1-gn(j))*V(i,j,k))*(gn(j)*V(i,j+1,k)+(1-gn(j))*V(i,j,k)))&
                      -((gn(j-1)*V(i,j,k)+(1-gn(j-1))*V(i,j-1,k))*(gn(j-1)*V(i,j,k)+(1-gn(j-1))*V(i,j-1,k))))/dyV(j)
           Vtmp=Vtmp+( (fe(i)*V(i+1,j,k)+(1-fe(i))*V(i,j,k))*(gp(j)*U(i,j+1,k)+(1-gp(j))*U(i,j,k))&
                      -(fe(i-1)*V(i,j,k)+(1-fe(i-1))*V(i-1,j,k))*(gp(j)*U(i-1,j+1,k)+(1-gp(j))*U(i-1,j,k)))/dxPr(i)
           Vtmp=Vtmp+( (ht(k)*V(i,j,k+1)+(1-ht(k))*V(i,j,k))*(gp(j)*W(i,j+1,k)+(1-gp(j))*W(i,j,k))&
                      -(ht(k-1)*V(i,j,k)+(1-ht(k-1))*V(i,j,k-1))*(gp(j)*W(i,j+1,k-1)+(1-gp(j))*W(i,j,k-1)))/dzPr(k)
           V2(i,j,k)=-dt*Vtmp
          enddo
         enddo
        enddo
        !$omp end parallel do
     endif
    end subroutine CDV





    subroutine CDW(W2,U,V,W)
    real(KND) :: W2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer nx,ny,nz,i,j,k
    real(KND) Ax,Ay,Az,Wtmp
    real(KND),allocatable,dimension(:),save:: ht,hp,fe,gn
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
                    W2(i,j,k)= - ((Az*(W(i,j,k+1)+W(i,j,k))*(W(i,j,k+1)+W(i,j,k))&
                    -Az*(W(i,j,k)+W(i,j,k-1))*(W(i,j,k)+W(i,j,k-1)))&
                    +(Ay*(W(i,j+1,k)+W(i,j,k))*(V(i,j,k+1)+V(i,j,k))&
                    -Ay*(W(i,j,k)+W(i,j-1,k))*(V(i,j-1,k)+V(i,j-1,k+1)))&
                    +(Ax*(W(i+1,j,k)+W(i,j,k))*(U(i,j,k+1)+U(i,j,k))&
                    -Ax*(W(i,j,k)+W(i-1,j,k))*(U(i-1,j,k+1)+U(i-1,j,k))))
                enddo
            enddo
        enddo
        !$omp end parallel do
     else
        if (called==0) then
         allocate(ht(0:nz),hp(1:nz),fe(0:nx),gn(0:ny))
         called=1
         forall (k=0:nz)      ht(k)=(zPr(k+1)-zW(k))/(zW(k+1)-zW(k))
         forall (k=1:nz)      hp(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
         forall (i=0:nx)      fe(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
         forall (j=0:ny)      gn(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
        endif

        !$omp parallel do private(i,j,k)
        do k=1,nz
         do j=1,ny
          do i=1,nx
           Wtmp=( ((ht(k)*W(i,j,k+1)+(1-ht(k))*W(i,j,k))*(ht(k)*W(i,j,k+1)+(1-ht(k))*W(i,j,k)))&
                      -((ht(k-1)*W(i,j,k)+(1-ht(k-1))*W(i,j,k-1))*(ht(k-1)*W(i,j,k)+(1-ht(k-1))*W(i,j,k-1))))/dzW(k)
           Wtmp=Wtmp+( (fe(i)*W(i+1,j,k)+(1-fe(i))*W(i,j,k))*(hp(k)*U(i,j,k+1)+(1-hp(k))*U(i,j,k))&
                      -(fe(i-1)*W(i,j,k)+(1-fe(i-1))*W(i-1,j,k))*(hp(k)*U(i-1,j,k+1)+(1-hp(k))*U(i-1,j,k)))/dxPr(i)
           Wtmp=Wtmp+( (gn(j)*W(i,j+1,k)+(1-gn(j))*W(i,j,k))*(hp(k)*V(i,j,k+1)+(1-hp(k))*V(i,j,k))&
                      -(gn(j-1)*W(i,j,k)+(1-gn(j-1))*W(i,j-1,k))*(hp(k)*V(i,j-1,k+1)+(1-hp(k))*V(i,j-1,k)))/dzPr(k)
           W2(i,j,k)=-dt*Wtmp
          enddo
         enddo
        enddo
        !$omp end parallel do
     endif
    end subroutine CDW



    subroutine CDS4U(U2,U,V,W)
      real(KND),dimension(-2:,-2:,-2:),intent(out) :: U2
      real(KND),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: mask(4) = (/ -1, 9, 9, -1 /)
      integer,parameter :: divcoef = 256
      real(KND),parameter :: mask2(2) = (/ 0.5_KND, 0.5_KND /)
      integer   :: i,j,k,l
      real(KND) :: flux

      !$omp parallel private(i,j,k,l,flux)
      !$omp workshare
      U2 = 0
      !$omp end workshare

      !$omp do
      do k=1,Unz
        do j=1,Uny
          do i=1,Unx+1
            flux = ( ( sum(mask * U(i-2:i+1,j,k)) )**2 ) / dxmin
            U2(i-1,j,k) = U2(i-1,j,k) - flux
            U2(i,j,k)   = U2(i,j,k)   + flux
          enddo
        enddo
      enddo
      !$omp end do

      !$omp do
      do k=1,Unz
        do j=1,Uny+1
          do i=1,Unx
            flux = ( sum(mask * U(i,j-2:j+1,k)) * sum(mask * V(i-1:i+2,j-1,k)) ) / dymin
            U2(i,j-1,k) = U2(i,j-1,k) - flux
            U2(i,j,k)   = U2(i,j,k)   + flux
          enddo
        enddo
      enddo
      !$omp end do

      do l=1,2
        !$omp do
        do k=l,Unz+1,2
          do j=1,Uny
            do i=1,Unx
              flux = ( sum(mask * U(i,j,k-2:k+1))) * sum(mask * W(i-1:i+2,j,k-1)) / dzmin
              U2(i,j,k-1) = U2(i,j,k-1) - flux
              U2(i,j,k)   = U2(i,j,k)   + flux
            enddo
          enddo
        enddo
        !$omp end do
      end do

      !$omp workshare
      U2 = U2 * dt / divcoef
      !$omp end workshare
      !$omp end parallel

    end subroutine CDS4U


    subroutine CDS4V(V2,U,V,W)
      real(KND),dimension(-2:,-2:,-2:),intent(out) :: V2
      real(KND),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: mask(4) = (/ -1, 9, 9, -1 /)
      integer,parameter :: divcoef = 256
      real(KND),parameter :: mask2(2) = (/ 0.5_KND, 0.5_KND /)
      integer   :: i,j,k,l
      real(KND) :: flux

      !$omp parallel private(i,j,k,l,flux)
      !$omp workshare
      V2 = 0
      !$omp end workshare

      !$omp do
      do k=1,Vnz
        do j=1,Vny
          do i=1,Vnx+1
            flux = ( sum(mask * V(i-2:i+1,j,k)) * sum(mask * U(i-1,j-1:j+2,k)) ) / dxmin
            V2(i-1,j,k) = V2(i-1,j,k) - flux
            V2(i,j,k)   = V2(i,j,k)   + flux
          enddo
        enddo
      enddo
      !$omp end do

      !$omp do
      do k=1,Vnz
        do j=1,Vny+1
          do i=1,Vnx
            flux = ( ( sum(mask * V(i,j-2:j+1,k)) )**2 ) / dymin
            V2(i,j-1,k) = V2(i,j-1,k) - flux
            V2(i,j,k)   = V2(i,j,k)   + flux
          enddo
        enddo
      enddo
      !$omp end do

      do l=1,2
        !$omp do
        do k=l,Vnz+1,2
          do j=1,Vny
            do i=1,Vnx
              flux = ( sum(mask * V(i,j,k-2:k+1)) * sum(mask * W(i,j-1:j+2,k-1)) ) / dzmin
              V2(i,j,k-1) = V2(i,j,k-1) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux
            enddo
          enddo
        enddo
        !$omp end do
      end do

      !$omp workshare
      V2 = V2 * dt / divcoef
      !$omp end workshare
      !$omp end parallel

    end subroutine CDS4V


    subroutine CDS4W(W2,U,V,W)
      real(KND),dimension(-2:,-2:,-2:),intent(out) :: W2
      real(KND),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
      integer,parameter :: mask(4) = (/ -1, 9, 9, -1 /)
      integer,parameter :: divcoef = 256
      real(KND),parameter :: mask2(2) = (/  0.5_KND, 0.5_KND/)
      integer   :: i,j,k,l
      real(KND) :: flux

      !$omp parallel private(i,j,k,l,flux)
      !$omp workshare
      W2 = 0
      !$omp end workshare

      !$omp do
      do k=1,Wnz
        do j=1,Wny
          do i=1,Wnx+1
            flux = ( sum(mask * W(i-2:i+1,j,k)) * sum(mask * U(i-1,j,k-1:k+2)) ) / dxmin
            W2(i-1,j,k) = W2(i-1,j,k) - flux
            W2(i,j,k)   = W2(i,j,k)   + flux
          enddo
        enddo
      enddo
      !$omp end do

      !$omp do
      do k=1,Wnz
        do j=1,Wny+1
          do i=1,Wnx
            flux = ( sum(mask * W(i,j-2:j+1,k)) * sum(mask * V(i,j-1,k-1:k+2)) ) / dymin
            W2(i,j-1,k) = W2(i,j-1,k) - flux
            W2(i,j,k)   = W2(i,j,k)   + flux
          enddo
        enddo
      enddo
      !$omp end do

      do l=1,2
        !$omp do
        do k=l,Wnz+1,2
          do j=1,Wny
            do i=1,Wnx
              flux = ( ( sum(mask * W(i,j,k-2:k+1)) )**2 ) / dzmin
              W2(i,j,k-1) = W2(i,j,k-1) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux
            enddo
          enddo
        enddo
        !$omp end do
      end do

      !$omp workshare
      W2 = W2 * dt / divcoef
      !$omp end workshare
      !$omp end parallel

    end subroutine CDS4W


  subroutine KAPPAU(U2,U,V,W)
  real(KND),intent(out)::U2(-2:,-2:,-2:) !Hunsdorfer et al. 1995, JCP
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer i,j,k
  real(KND) A,SL,SR,FLUX
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3):: SLOPE
  real(KND),parameter::eps=1e-8


  A=dt/2._KND

  SLOPE=0
  do i=0,Unx
   do j=1,Uny
    do k=1,Unz
     if (U(i,j,k)+U(i+1,j,k)>0) then
      SR=(U(i+1,j,k)-U(i,j,k))!/dxU(i)
      SL=(U(i,j,k)-U(i-1,j,k))!/dxU(i-1)
     else
      SR=(U(i,j,k)-U(i+1,j,k))!/dxU(i)
      SL=(U(i+1,j,k)-U(i+2,j,k))!/dxU(i-1)
     endif
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps)/(SL+eps))
    enddo
   enddo
  enddo

  do i=1,Unx
   do j=1,Uny
    do k=0,Unz
     if ((U(i,j,k)+U(i+1,j,k))>0) then
      FLUX=(U(i,j,k)+U(i+1,j,k))*(U(i,j,k)+(U(i,j,k)-U(i-1,j,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=(U(i,j,k)+U(i+1,j,k))*(U(i+1,j,k)+(U(i+1,j,k)-U(i+2,j,k))*SLOPE(i,j,k)/2._KND)
     endif

     U2(i,j,k)=-A*FLUX/dxPr(i)
     U2(i+1,j,k)=U2(i+1,j,k)+A*FLUX/dxPr(i+1)
    enddo
   enddo
  enddo

  SLOPE=0
  do i=1,Unx
   do j=0,Uny
    do k=1,Unz
     if (V(i,j,k)+V(i+1,j,k)>0) then
      SR=(U(i,j+1,k)-U(i,j,k))!/dxU(i)
      SL=(U(i,j,k)-U(i,j-1,k))!/dxU(i-1)
     else
      SR=(U(i,j,k)-U(i,j+1,k))!/dxU(i)
      SL=(U(i,j+1,k)-U(i,j+2,k))!/dxU(i-1)
     endif
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps)/(SL+eps))
    enddo
   enddo
  enddo


  do i=1,Unx
   do j=0,Uny
    do k=1,Unz
     if (V(i,j,k)+V(i+1,j,k)>0) then
      FLUX=(V(i,j,k)+V(i+1,j,k))*(U(i,j,k)+(U(i,j,k)-U(i,j-1,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=(V(i,j,k)+V(i+1,j,k))*(U(i,j+1,k)+(U(i,j+1,k)-U(i,j+2,k))*SLOPE(i,j,k)/2._KND)
     endif

     U2(i,j,k)=U2(i,j,k)-A*FLUX/dyPr(j)
     U2(i,j+1,k)=U2(i,j+1,k)+A*FLUX/dyPr(j+1)
    enddo
   enddo
  enddo

  SLOPE=0
  do i=1,Unx
   do j=1,Uny
    do k=0,Unz
    if (W(i,j,k)+W(i+1,j,k)>0) then
      SR=(U(i,j,k+1)-U(i,j,k))!/dxU(i)
      SL=(U(i,j,k)-U(i,j,k-1))!/dxU(i-1)
     else
      SR=(U(i,j,k)-U(i,j,k+1))!/dxU(i)
      SL=(U(i,j,k+1)-U(i,j,k+2))!/dxU(i-1)
     endif
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps)/(SL+eps))
    enddo
   enddo
  enddo

  do i=1,Unx
   do j=1,Uny
    do k=0,Unz
     if (W(i,j,k)+W(i+1,j,k)>0) then
      FLUX=(W(i,j,k)+W(i+1,j,k))*(U(i,j,k)+(U(i,j,k)-U(i,j,k-1))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=(W(i,j,k)+W(i+1,j,k))*(U(i,j,k+1)+(U(i,j,k+1)-U(i,j,k+2))*SLOPE(i,j,k)/2._KND)
     endif

     U2(i,j,k)=U2(i,j,k)-A*FLUX/dzPr(k)
     U2(i,j,k+1)=U2(i,j,k+1)+A*FLUX/dzPr(k+1)
    enddo
   enddo
  enddo
  endsubroutine KAPPAU



  pure subroutine KAPPAV(V2,U,V,W) !Kappa scheme with flux limiter
  real(KND),intent(out)::V2(-2:,-2:,-2:) !Hunsdorfer et al. 1995, JCP
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer i,j,k
  real(KND) A,SL,SR,FLUX
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3):: SLOPE
  real(KND),parameter::eps=1e-8


  A=dt/2._KND

  SLOPE=0
  do i=0,Vnx
   do j=1,Vny
    do k=1,Vnz
     if (U(i,j,k)+U(i,j+1,k)>0) then
      SR=(V(i+1,j,k)-V(i,j,k))!/dxU(i)
      SL=(V(i,j,k)-V(i-1,j,k))!/dxU(i-1)
     else
      SR=(V(i,j,k)-V(i+1,j,k))!/dxU(i)
      SL=(V(i+1,j,k)-V(i+2,j,k))!/dxU(i-1)
     endif
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps)/(SL+eps))
    enddo
   enddo
  enddo


  do i=0,Vnx
   do j=1,Vny
    do k=1,Vnz
     if (U(i,j,k)+U(i,j+1,k)>0) then
      FLUX=(U(i,j,k)+U(i,j+1,k))*(V(i,j,k)+(V(i,j,k)-V(i-1,j,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=(U(i,j,k)+U(i,j+1,k))*(V(i+1,j,k)+(V(i+1,j,k)-V(i+2,j,k))*SLOPE(i,j,k)/2._KND)
     endif
     V2(i,j,k)=-A*FLUX/dxPr(i)
     V2(i+1,j,k)=V2(i+1,j,k)+A*FLUX/dxPr(i+1)
    enddo
   enddo
  enddo

  SLOPE=0
  do i=1,Vnx
   do j=0,Vny
    do k=1,Vnz
     if (V(i,j,k)+V(i,j+1,k)>0) then
      SR=(V(i,j+1,k)-V(i,j,k))!/dxU(i)
      SL=(V(i,j,k)-V(i,j-1,k))!/dxU(i-1)
     else
      SR=(V(i,j,k)-V(i,j+1,k))!/dxU(i)
      SL=(V(i,j+1,k)-V(i,j+2,k))!/dxU(i-1)
     endif
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps)/(SL+eps))
    enddo
   enddo
  enddo


  do i=1,Vnx
   do j=0,Vny
    do k=1,Vnz
     if (V(i,j,k)+V(i,j+1,k)>0) then
      FLUX=(V(i,j,k)+V(i,j+1,k))*(V(i,j,k)+(V(i,j,k)-V(i,j-1,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=(V(i,j,k)+V(i,j+1,k))*(V(i,j+1,k)+(V(i,j+1,k)-V(i,j+2,k))*SLOPE(i,j,k)/2._KND)
     endif

     V2(i,j,k)=V2(i,j,k)-A*FLUX/dyPr(j)
     V2(i,j+1,k)=V2(i,j+1,k)+A*FLUX/dyPr(j+1)
    enddo
   enddo
  enddo

  SLOPE=0
  do i=1,Vnx
   do j=1,Vny
    do k=0,Vnz
    if (W(i,j,k)+W(i,j+1,k)>0) then
      SR=(V(i,j,k+1)-V(i,j,k))!/dxU(i)
      SL=(V(i,j,k)-V(i,j,k-1))!/dxU(i-1)
     else
      SR=(V(i,j,k)-V(i,j,k+1))!/dxU(i)
      SL=(V(i,j,k+1)-V(i,j,k+2))!/dxU(i-1)
     endif
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps)/(SL+eps))
    enddo
   enddo
  enddo

  do i=1,Vnx
   do j=1,Vny
    do k=0,Vnz
     if (W(i,j,k)+W(i,j+1,k)>0) then
      FLUX=(W(i,j,k)+W(i,j+1,k))*(V(i,j,k)+(V(i,j,k)-V(i,j,k-1))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=(W(i,j,k)+W(i,j+1,k))*(V(i,j,k+1)+(V(i,j,k+1)-V(i,j,k+2))*SLOPE(i,j,k)/2._KND)
     endif

     V2(i,j,k)=V2(i,j,k)-A*FLUX/dzPr(k)
     V2(i,j,k+1)=V2(i,j,k+1)+A*FLUX/dzPr(k+1)
    enddo
   enddo
  enddo
  endsubroutine KAPPAV

  pure subroutine KAPPAW(W2,U,V,W) !Kappa scheme with flux limiter
  real(KND),intent(out)::W2(-2:,-2:,-2:) !Hunsdorfer et al. 1995, JCP
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer i,j,k
  real(KND) A,SL,SR,FLUX
  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3):: SLOPE
  real(KND),parameter::eps=1e-8


  A=dt/2._KND

  SLOPE=0
  do i=0,Wnx
   do j=1,Wny
    do k=1,Wnz
     if (U(i,j,k)+U(i,j,k+1)>0) then
      SR=(W(i+1,j,k)-W(i,j,k))!/dxU(i)
      SL=(W(i,j,k)-W(i-1,j,k))!/dxU(i-1)
     else
      SR=(W(i,j,k)-W(i+1,j,k))!/dxU(i)
      SL=(W(i+1,j,k)-W(i+2,j,k))!/dxU(i-1)
     endif
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps)/(SL+eps))
    enddo
   enddo
  enddo


  do i=0,Wnx
   do j=1,Wny
    do k=1,Wnz
     if (U(i,j,k)+U(i,j,k+1)>0) then
      FLUX=(U(i,j,k)+U(i,j,k+1))*(W(i,j,k)+(W(i,j,k)-W(i-1,j,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=(U(i,j,k)+U(i,j,k+1))*(W(i+1,j,k)+(W(i+1,j,k)-W(i+2,j,k))*SLOPE(i,j,k)/2._KND)
     endif
     W2(i,j,k)=-A*FLUX/dxPr(i)
     W2(i+1,j,k)=W2(i+1,j,k)+A*FLUX/dxPr(i+1)
    enddo
   enddo
  enddo

  SLOPE=0
  do i=1,Wnx
   do j=0,Wny
    do k=1,Wnz
     if (V(i,j,k)+V(i,j,k+1)>0) then
      SR=(W(i,j+1,k)-W(i,j,k))!/dxU(i)
      SL=(W(i,j,k)-W(i,j-1,k))!/dxU(i-1)
     else
      SR=(W(i,j,k)-W(i,j+1,k))!/dxU(i)
      SL=(W(i,j+1,k)-W(i,j+2,k))!/dxU(i-1)
     endif
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps)/(SL+eps))
    enddo
   enddo
  enddo


  do i=1,Wnx
   do j=0,Wny
    do k=1,Wnz
     if (V(i,j,k)+V(i,j,k+1)>0) then
      FLUX=(V(i,j,k)+V(i,j,k+1))*(W(i,j,k)+(W(i,j,k)-W(i,j-1,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=(V(i,j,k)+V(i,j,k+1))*(W(i,j+1,k)+(W(i,j+1,k)-W(i,j+2,k))*SLOPE(i,j,k)/2._KND)
     endif

     W2(i,j,k)=W2(i,j,k)-A*FLUX/dyPr(j)
     W2(i,j+1,k)=W2(i,j+1,k)+A*FLUX/dyPr(j+1)
    enddo
   enddo
  enddo

  SLOPE=0
  do i=1,Wnx
   do j=1,Wny
    do k=0,Wnz
    if (W(i,j,k)+W(i,j,k+1)>0) then
      SR=(W(i,j,k+1)-W(i,j,k))!/dxU(i)
      SL=(W(i,j,k)-W(i,j,k-1))!/dxU(i-1)
     else
      SR=(W(i,j,k)-W(i,j,k+1))!/dxU(i)
      SL=(W(i,j,k+1)-W(i,j,k+2))!/dxU(i-1)
     endif
     SLOPE(i,j,k)=FLUXLIMITER((SR+eps)/(SL+eps))
    enddo
   enddo
  enddo

  do i=1,Wnx
   do j=1,Wny
    do k=0,Wnz
     if (W(i,j,k)+W(i,j,k+1)>0) then
      FLUX=(W(i,j,k)+W(i,j,k+1))*(W(i,j,k)+(W(i,j,k)-W(i,j,k-1))*SLOPE(i,j,k)/2._KND)
     else
      FLUX=(W(i,j,k)+W(i,j,k+1))*(W(i,j,k+1)+(W(i,j,k+1)-W(i,j,k+2))*SLOPE(i,j,k)/2._KND)
     endif

     W2(i,j,k)=W2(i,j,k)-A*FLUX/dzPr(k)
     W2(i,j,k+1)=W2(i,j,k+1)+A*FLUX/dzPr(k+1)
    enddo
   enddo
  enddo
  endsubroutine KAPPAW



end module CDS

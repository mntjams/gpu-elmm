!
! File:   CDS.f90
! Author: lada
!
! Created on 15. leden 2006, 12:03
!

module CDS
    use PARAMETERS
    use BOUNDARIES
    implicit none
    
    contains
    
    
    subroutine CDU(U2,U,V,W,coef)        
    real(KND) :: U2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(KND) coef
    integer nx,ny,nz,i,j,k
    real(KND) Ax,Ay,Az,Utmp
    real(KND),allocatable,dimension(:),save:: fe,fp,gn,ht
    integer,save:: called=0

     nx=Unx
     ny=Uny
     nz=Unz
     call BOUND_CONDU(U)
     call BOUND_CONDV(V)
     call BOUND_CONDW(W)

     if (gridtype==UNIFORMGRID) then
        Ax=0.25*coef*dt/dxmin
        Ay=0.25*coef*dt/dymin
        Az=0.25*coef*dt/dzmin
        
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    U2(i,j,k)=U2(i,j,k) - coef*((Ax*((U(i+1,j,k)*U(i+1,j,k))+(U(i,j,k)*U(i,j,k)))&
                    -Ax*((U(i,j,k)*U(i,j,k))+(U(i-1,j,k)*U(i-1,j,k))))&
                    +(Ay*(U(i,j+1,k)+U(i,j,k))*(V(i+1,j,k)+V(i,j,k))&
                    -Ay*(U(i,j,k)+U(i,j-1,k))*(V(i+1,j-1,k)+V(i,j-1,k)))&
                    +(Az*(U(i,j,k+1)+U(i,j,k))*(W(i+1,j,k)+W(i,j,k))&
                    -Az*(U(i,j,k)+U(i,j,k-1))*(W(i+1,j,k-1)+W(i,j,k-1))))
                enddo
            enddo
        enddo
     else
        if (called==0) then
         allocate(fe(0:nx),fp(1:nx),gn(0:ny),ht(0:nz))
         called=1
         forall (i=0:nx)      fe(i)=(xPr(i+1)-xU(i-1))/(xU(i)-xU(i-1))
         forall (i=1:nx)      fp(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))        
         forall (j=0:ny)      gn(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
         forall (k=0:nz)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
        endif

        do k=1,nz
         do j=1,ny
          do i=1,nx
           Utmp=Utmp+( (fe(i)*U(i+1,j,k)*U(i+1,j,k) + (1-fe(i))*U(i,j,k)*U(i,j,k))&
                      -(fe(i-1)*U(i,j,k)*U(i,j,k) + (1-fe(i-1))*U(i-1,j,k)*U(i-1,j,k)))/dxU(i)
           Utmp=Utmp+( (gn(j)*U(i,j+1,k)+(1-gn(j))*U(i,j,k))*(fp(i)*V(i+1,j,k)+(1-fp(i))*V(i,j,k))&
                      -(gn(j-1)*U(i,j,k)+(1-gn(j-1))*U(i,j-1,k))*(fp(i)*V(i+1,j-1,k)+(1-fp(i))*V(i,j-1,k)))/dyPr(j)
           Utmp=Utmp+( (ht(k)*U(i,j,k+1)+(1-ht(k))*U(i,j,k))*(fp(i)*W(i+1,j,k)+(1-fp(i))*W(i,j,k))&
                      -(ht(k-1)*U(i,j,k)+(1-ht(k-1))*U(i,j,k-1))*(fp(i)*W(i+1,j,k-1)+(1-fp(i))*W(i,j,k-1)))/dzPr(k)
           U2(i,j,k)=U2(i,j,k)-dt*coef*Utmp
          enddo
         enddo
        enddo
     endif
    end subroutine CDU
        
    subroutine CDV(V2,U,V,W,coef)
    real(KND) :: V2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(KND) coef
    integer nx,ny,nz,i,j,k
    real(KND) Ax,Ay,Az,Vtmp
    real(KND),allocatable,dimension(:),save:: gn,gp,fe,ht
    integer,save:: called=0


     nx=Vnx
     ny=Vny
     nz=Vnz
     call BOUND_CONDU(U)
     call BOUND_CONDV(V)
     call BOUND_CONDW(W)

     if (gridtype==UNIFORMGRID) then
        Ax=0.25*coef*dt/dxmin
        Ay=0.25*coef*dt/dymin
        Az=0.25*coef*dt/dzmin
        
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    V2(i,j,k)=V2(i,j,k) - coef*((Ay*((V(i,j+1,k)*V(i,j+1,k))+(V(i,j,k)*V(i,j,k)))&
                    -Ay*((V(i,j,k)*V(i,j,k))+(V(i,j-1,k)*V(i,j-1,k))))&
                    +(Ax*(V(i+1,j,k)+V(i,j,k))*(U(i,j+1,k)+U(i,j,k))&
                    -Ax*(V(i,j,k)+V(i-1,j,k))*(U(i-1,j+1,k)+U(i-1,j,k)))&
                    +(Az*(V(i,j,k+1)+V(i,j,k))*(W(i,j+1,k)+W(i,j,k))&
                    -Az*(V(i,j,k)+V(i,j,k-1))*(W(i,j+1,k-1)+W(i,j,k-1))))
                enddo
            enddo
        enddo
     else
        if (called==0) then
         allocate(gn(0:ny),gp(1:ny),fe(0:nx),ht(0:nz))
         called=1
         forall (j=0:ny)      gn(j)=(yPr(j+1)-yV(j-1))/(yV(j)-yV(j-1))
         forall (j=1:ny)      gp(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))        
         forall (i=0:nx)      fe(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
         forall (k=0:nz)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
        endif

        do k=1,nz
         do j=1,ny
          do i=1,nx
           Vtmp=Vtmp+( (gn(j)*V(i,j+1,k)*V(i,j+1,k) + (1-gn(j))*V(i,j,k)*V(i,j,k))&
                      -(gn(j-1)*V(i,j,k)*V(i,j,k) + (1-gn(j-1))*V(i,j-1,k)*V(i,j-1,k)))/dyV(j)
           Vtmp=Vtmp+( (fe(i)*V(i+1,j,k)+(1-fe(i))*V(i,j,k))*(gp(j)*U(i,j+1,k)+(1-gp(j))*U(i,j,k))&
                      -(fe(i-1)*V(i,j,k)+(1-fe(i-1))*V(i-1,j,k))*(gp(j)*U(i-1,j+1,k)+(1-gp(j))*U(i-1,j,k)))/dxPr(i)
           Vtmp=Vtmp+( (ht(k)*V(i,j,k+1)+(1-ht(k))*V(i,j,k))*(gp(j)*W(i,j+1,k)+(1-gp(j))*W(i,j,k))&
                      -(ht(k-1)*V(i,j,k)+(1-ht(k-1))*V(i,j,k-1))*(gp(j)*W(i,j+1,k-1)+(1-gp(j))*W(i,j,k-1)))/dzPr(k)
           V2(i,j,k)=V2(i,j,k)-dt*coef*Vtmp
          enddo
         enddo
        enddo
     endif
    end subroutine CDV            
    
    subroutine CDW(W2,U,V,W,coef)
    real(KND) :: W2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(KND) coef
    integer nx,ny,nz,i,j,k
    real(KND) Ax,Ay,Az,Wtmp
    real(KND),allocatable,dimension(:),save:: ht,hp,fe,gn
    integer,save:: called=0

     nx=Wnx
     ny=Wny
     nz=Wnz
     call BOUND_CONDU(U)
     call BOUND_CONDV(V)
     call BOUND_CONDW(W)

     if (gridtype==UNIFORMGRID) then
        Ax=0.25*coef*dt/dxmin
        Ay=0.25*coef*dt/dymin
        Az=0.25*coef*dt/dzmin
        
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    W2(i,j,k)=W2(i,j,k) - coef*((Az*((W(i,j,k+1)*W(i,j,k+1))+(W(i,j,k)*W(i,j,k)))&
                    -Az*((W(i,j,k)*W(i,j,k))+(W(i,j,k-1)*W(i,j,k-1))))&
                    +(Ay*(W(i,j+1,k)+W(i,j,k))*(V(i,j,k+1)+V(i,j,k))&
                    -Ay*(W(i,j,k)+W(i,j-1,k))*(V(i,j-1,k)+V(i,j-1,k+1)))&
                    +(Ax*(W(i+1,j,k)+W(i,j,k))*(U(i,j,k+1)+U(i,j,k))&
                    -Ax*(W(i,j,k)+W(i-1,j,k))*(U(i-1,j,k+1)+U(i-1,j,k))))
                enddo
            enddo
        enddo
     else
        if (called==0) then
         allocate(ht(0:nz),hp(1:nz),fe(0:nx),gn(0:ny))
         called=1
         forall (k=0:nz)      ht(k)=(zPr(k+1)-zW(k-1))/(zW(k)-zW(k-1))
         forall (k=1:nz)      hp(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))        
         forall (i=0:nx)      fe(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
         forall (j=0:ny)      gn(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
        endif

        do k=1,nz
         do j=1,ny
          do i=1,nx
           Wtmp=Wtmp+( (ht(k)*W(i,j,k+1)*W(i,j,k+1) + (1-ht(k))*W(i,j,k)*W(i,j,k))&
                      -(ht(k-1)*W(i,j,k)*W(i,j,k) + (1-ht(k-1))*W(i,j,k-1)*W(i,j,k-1)))/dzW(k)
           Wtmp=Wtmp+( (fe(i)*W(i+1,j,k)+(1-fe(i))*W(i,j,k))*(hp(k)*U(i,j,k+1)+(1-hp(k))*U(i,j,k))&
                      -(fe(i-1)*W(i,j,k)+(1-fe(i-1))*W(i-1,j,k))*(hp(k)*U(i-1,j,k+1)+(1-hp(k))*U(i-1,j,k)))/dxPr(i)
           Wtmp=Wtmp+( (gn(k)*V(i,j+1,k)+(1-gn(k))*V(i,j,k))*(hp(k)*V(i,j,k+1)+(1-hp(k))*V(i,j,k))&
                      -(gn(k-1)*V(i,j,k)+(1-gn(k-1))*V(i,j-1,k))*(hp(j)*V(i,j-1,k+1)+(1-hp(k))*V(i,j-1,k)))/dzPr(k)
           W2(i,j,k)=W2(i,j,k)-dt*coef*Wtmp
          enddo
         enddo
        enddo
     endif
    end subroutine CDW            
            
          












  
    subroutine CDU2(U2,U,V,W,coef)
    real(KND) :: U2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(KND) coef
    integer nx,ny,nz,i,j,k
    real(KND) Ax,Ay,Az,Utmp
    real(KND),allocatable,dimension(:),save:: fe,fp,gn,ht
    integer,save:: called=0

     nx=Unx
     ny=Uny
     nz=Unz
     call BOUND_CONDU(U)
     call BOUND_CONDV(V)
     call BOUND_CONDW(W)

     if (gridtype==UNIFORMGRID) then
        Ax=0.25*coef*dt/dxmin
        Ay=0.25*coef*dt/dymin
        Az=0.25*coef*dt/dzmin
        
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    U2(i,j,k)=U2(i,j,k) - coef*((Ax*(U(i+1,j,k)+U(i,j,k))*(U(i+1,j,k)+U(i,j,k))&
                    -Ax*(U(i,j,k)+U(i-1,j,k))*(U(i,j,k)+U(i-1,j,k)))&
                    +(Ay*(U(i,j+1,k)+U(i,j,k))*(V(i+1,j,k)+V(i,j,k))&
                    -Ay*(U(i,j,k)+U(i,j-1,k))*(V(i+1,j-1,k)+V(i,j-1,k)))&
                    +(Az*(U(i,j,k+1)+U(i,j,k))*(W(i+1,j,k)+W(i,j,k))&
                    -Az*(U(i,j,k)+U(i,j,k-1))*(W(i+1,j,k-1)+W(i,j,k-1))))
                enddo
            enddo
        enddo
     else
        if (called==0) then
         allocate(fe(0:nx),fp(1:nx),gn(0:ny),ht(0:nz))
         called=1
         forall (i=0:nx)      fe(i)=(xPr(i+1)-xU(i))/(xU(i+1)-xU(i))
         forall (i=1:nx)      fp(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
         forall (j=0:ny)      gn(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
         forall (k=0:nz)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
        endif

        do k=1,nz
         do j=1,ny
          do i=1,nx
           Utmp=( ((fe(i)*U(i+1,j,k)+(1-fe(i))*U(i,j,k))*(fe(i)*U(i+1,j,k)+(1-fe(i))*U(i,j,k)))&
                      -((fe(i-1)*U(i,j,k)+(1-fe(i-1))*U(i-1,j,k))*(fe(i-1)*U(i,j,k)+(1-fe(i-1))*U(i-1,j,k))))/dxU(i)
           Utmp=Utmp+( (gn(j)*U(i,j+1,k)+(1-gn(j))*U(i,j,k))*(fp(i)*V(i+1,j,k)+(1-fp(i))*V(i,j,k))&
                      -(gn(j-1)*U(i,j,k)+(1-gn(j-1))*U(i,j-1,k))*(fp(i)*V(i+1,j-1,k)+(1-fp(i))*V(i,j-1,k)))/dyPr(j)
           Utmp=Utmp+( (ht(k)*U(i,j,k+1)+(1-ht(k))*U(i,j,k))*(fp(i)*W(i+1,j,k)+(1-fp(i))*W(i,j,k))&
                      -(ht(k-1)*U(i,j,k)+(1-ht(k-1))*U(i,j,k-1))*(fp(i)*W(i+1,j,k-1)+(1-fp(i))*W(i,j,k-1)))/dzPr(k)
           U2(i,j,k)=U2(i,j,k)-dt*coef*Utmp
          enddo
         enddo
        enddo
     endif
    end subroutine CDU2
        





    subroutine CDV2(V2,U,V,W,coef)
    real(KND) :: V2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(KND) coef
    integer nx,ny,nz,i,j,k
    real(KND) Ax,Ay,Az,Vtmp
    real(KND),allocatable,dimension(:),save:: gn,gp,fe,ht
    integer,save:: called=0


     nx=Vnx
     ny=Vny
     nz=Vnz
     call BOUND_CONDU(U)
     call BOUND_CONDV(V)
     call BOUND_CONDW(W)

     if (gridtype==UNIFORMGRID) then
        Ax=0.25*coef*dt/dxmin
        Ay=0.25*coef*dt/dymin
        Az=0.25*coef*dt/dzmin
        
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    V2(i,j,k)=V2(i,j,k) - coef*((Ay*(V(i,j+1,k)+V(i,j,k))*(V(i,j+1,k)+V(i,j,k))&
                    -Ay*(V(i,j,k)+V(i,j-1,k))*(V(i,j,k)+V(i,j-1,k)))&
                    +(Ax*(V(i+1,j,k)+V(i,j,k))*(U(i,j+1,k)+U(i,j,k))&
                    -Ax*(V(i,j,k)+V(i-1,j,k))*(U(i-1,j+1,k)+U(i-1,j,k)))&
                    +(Az*(V(i,j,k+1)+V(i,j,k))*(W(i,j+1,k)+W(i,j,k))&
                    -Az*(V(i,j,k)+V(i,j,k-1))*(W(i,j+1,k-1)+W(i,j,k-1))))
                enddo
            enddo
        enddo
     else
        if (called==0) then
         allocate(gn(0:ny),gp(1:ny),fe(0:nx),ht(0:nz))
         called=1
         forall (j=0:ny)      gn(j)=(yPr(j+1)-yV(j))/(yV(j+1)-yV(j))
         forall (j=1:ny)      gp(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
         forall (i=0:nx)      fe(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
         forall (k=0:nz)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
        endif

        do k=1,nz
         do j=1,ny
          do i=1,nx
           Vtmp=( ((gn(j)*V(i,j+1,k)+(1-gn(j))*V(i,j,k))*(gn(j)*V(i,j+1,k)+(1-gn(j))*V(i,j,k)))&
                      -((gn(j-1)*V(i,j,k)+(1-gn(j-1))*V(i,j-1,k))*(gn(j-1)*V(i,j,k)+(1-gn(j-1))*V(i,j-1,k))))/dyV(j)
           Vtmp=Vtmp+( (fe(i)*V(i+1,j,k)+(1-fe(i))*V(i,j,k))*(gp(j)*U(i,j+1,k)+(1-gp(j))*U(i,j,k))&
                      -(fe(i-1)*V(i,j,k)+(1-fe(i-1))*V(i-1,j,k))*(gp(j)*U(i-1,j+1,k)+(1-gp(j))*U(i-1,j,k)))/dxPr(i)
           Vtmp=Vtmp+( (ht(k)*V(i,j,k+1)+(1-ht(k))*V(i,j,k))*(gp(j)*W(i,j+1,k)+(1-gp(j))*W(i,j,k))&
                      -(ht(k-1)*V(i,j,k)+(1-ht(k-1))*V(i,j,k-1))*(gp(j)*W(i,j+1,k-1)+(1-gp(j))*W(i,j,k-1)))/dzPr(k)
           V2(i,j,k)=V2(i,j,k)-dt*coef*Vtmp
          enddo
         enddo
        enddo
     endif
    end subroutine CDV2     
    




    subroutine CDW2(W2,U,V,W,coef)
    real(KND) :: W2(-2:,-2:,-2:),U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(KND) coef
    integer nx,ny,nz,i,j,k
    real(KND) Ax,Ay,Az,Wtmp
    real(KND),allocatable,dimension(:),save:: ht,hp,fe,gn
    integer,save:: called=0

     nx=Wnx
     ny=Wny
     nz=Wnz
     call BOUND_CONDU(U)
     call BOUND_CONDV(V)
     call BOUND_CONDW(W)

     if (gridtype==UNIFORMGRID) then
        Ax=0.25*coef*dt/dxmin
        Ay=0.25*coef*dt/dymin
        Az=0.25*coef*dt/dzmin
        
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    W2(i,j,k)=W2(i,j,k) - coef*((Az*(W(i,j,k+1)+W(i,j,k))*(W(i,j,k+1)+W(i,j,k))&
                    -Az*(W(i,j,k)+W(i,j,k-1))*(W(i,j,k)+W(i,j,k-1)))&
                    +(Ay*(W(i,j+1,k)+W(i,j,k))*(V(i,j,k+1)+V(i,j,k))&
                    -Ay*(W(i,j,k)+W(i,j-1,k))*(V(i,j-1,k)+V(i,j-1,k+1)))&
                    +(Ax*(W(i+1,j,k)+W(i,j,k))*(U(i,j,k+1)+U(i,j,k))&
                    -Ax*(W(i,j,k)+W(i-1,j,k))*(U(i-1,j,k+1)+U(i-1,j,k))))
                enddo
            enddo
        enddo
     else
        if (called==0) then
         allocate(ht(0:nz),hp(1:nz),fe(0:nx),gn(0:ny))
         called=1
         forall (k=0:nz)      ht(k)=(zPr(k+1)-zW(k))/(zW(k+1)-zW(k))
         forall (k=1:nz)      hp(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))        
         forall (i=0:nx)      fe(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
         forall (j=0:ny)      gn(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
        endif

        do k=1,nz
         do j=1,ny
          do i=1,nx
           Wtmp=( ((ht(k)*W(i,j,k+1)+(1-ht(k))*W(i,j,k))*(ht(k)*W(i,j,k+1)+(1-ht(k))*W(i,j,k)))&
                      -((ht(k-1)*W(i,j,k)+(1-ht(k-1))*W(i,j,k-1))*(ht(k-1)*W(i,j,k)+(1-ht(k-1))*W(i,j,k-1))))/dzW(k)
           Wtmp=Wtmp+( (fe(i)*W(i+1,j,k)+(1-fe(i))*W(i,j,k))*(hp(k)*U(i,j,k+1)+(1-hp(k))*U(i,j,k))&
                      -(fe(i-1)*W(i,j,k)+(1-fe(i-1))*W(i-1,j,k))*(hp(k)*U(i-1,j,k+1)+(1-hp(k))*U(i-1,j,k)))/dxPr(i)
           Wtmp=Wtmp+( (gn(j)*W(i,j+1,k)+(1-gn(j))*W(i,j,k))*(hp(k)*V(i,j,k+1)+(1-hp(k))*V(i,j,k))&
                      -(gn(j-1)*W(i,j,k)+(1-gn(j-1))*W(i,j-1,k))*(hp(k)*V(i,j-1,k+1)+(1-hp(k))*V(i,j-1,k)))/dzPr(k)
           W2(i,j,k)=W2(i,j,k)-dt*coef*Wtmp
          enddo
         enddo
        enddo
     endif
    end subroutine CDW2            
            
            
  subroutine KAPPAU(U2,U,V,W,coef)
  real(KND),intent(INOUT)::U2(-2:,-2:,-2:),U(-2:,-2:,-2:) !Hunsdorfer et al. 1995, JCP
  real(KND),intent(IN):: V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer i,j,k
  real(KND) A,SL,SR,FLUX
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3):: SLOPE
  real(KND),parameter::eps=1E-8

  call BOUND_CONDU(U)

  A=coef*dt/2._KND

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

     U2(i,j,k)=U2(i,j,k)-A*FLUX/dxPr(i)
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



  pure subroutine KAPPAV(V2,U,V,W,coef) !Kappa scheme with flux limiter
  real(KND),intent(INOUT)::V2(-2:,-2:,-2:),V(-2:,-2:,-2:) !Hunsdorfer et al. 1995, JCP
  real(KND),intent(IN):: U(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer i,j,k
  real(KND) A,SL,SR,FLUX
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3):: SLOPE
  real(KND),parameter::eps=1E-8

  call BOUND_CONDV(V)
  
  A=coef*dt/2._KND

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
     V2(i,j,k)=V2(i,j,k)-A*FLUX/dxPr(i)
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

  pure subroutine KAPPAW(W2,U,V,W,coef) !Kappa scheme with flux limiter
  real(KND),intent(INOUT)::W2(-2:,-2:,-2:),W(-2:,-2:,-2:) !Hunsdorfer et al. 1995, JCP
  real(KND),intent(IN):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),coef
  integer i,j,k
  real(KND) A,SL,SR,FLUX
  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3):: SLOPE
  real(KND),parameter::eps=1E-8

  call BOUND_CONDW(W)
  
  A=coef*dt/2._KND

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
     W2(i,j,k)=W2(i,j,k)-A*FLUX/dxPr(i)
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

module POISSON

 use PARAMETERS
 use BOUNDARIES
 use MULTIGRID,   only: POISSMG
 use MULTIGRID2d, only: POISSMG2d
 use FISHPOISSON

implicit none

 integer MUDbnx,MUDbny,MUDbnz,MUDlmg

contains


subroutine PR_CORRECT(U,V,W,Pr,coef,Q)                   !Pressure correction
  real(KND),dimension(-2:,-2:,-2:),intent(inout):: U,V,W !Phi is computed in Poisson eq. with div of U in RHS
  real(KND),dimension(1:,1:,1:),intent(inout):: Pr       !Depending on active projection method Phi becomes new pressure
  real(KND),dimension(1:,1:,1:),optional,intent(in)::Q   !or is added to last pressure
  real(KND),intent(in):: coef
  real(KND),allocatable,dimension(:,:,:):: RHS,RHS2      !U,V,W velocity field for correction U(-2:Unx,-2:Uny,-2:Unz)
  real(KND),save,allocatable:: Phi(:,:,:)                !Pr(1:Prnx+1,1:Prny+1,Prnz+1) pressure
  real(KND) Phiref                                       !coef cofficient from Runge Kutta, Q mass sources from immersed boundary
  real(TIM) dt2,dt3                                      !RHS right hand side of eq. with divergence of U
  integer i,j,k                                          !Phi computed pseudopressure, saved as first guess for next time
  integer nx,ny,nz
  real(KND) S,S2,Apr
  integer maxI,maxJ,maxK
  character(70):: str
  integer,save:: called=0
  real(KND),save:: lastcoef
  integer maxiterbackup

  if (called==0) then
    allocate(Phi(0:Prnx+1,0:Prny+1,0:Prnz+1))
    Phi=0
    called=1
    lastcoef=coef
  else
   called=1
  endif

  Phi=lastcoef*Phi/coef
  lastcoef=coef

  projectiontype=1

  write (*,*) "Preparing pressure correction"

  if (Btype(Ea)==PERIODIC.and.poissmet==2) then
   nx=Prnx+1
  else
   nx=Prnx
  endif
  if (Btype(No)==PERIODIC.and.poissmet==2) then
   ny=Prny+1
  else
   ny=Prny
  endif
  if (Btype(To)==PERIODIC.and.poissmet==2) then
   nz=Prnz+1
  else
   nz=Prnz
  endif

 allocate(RHS(nx,ny,nz))

 if (poissmet==3)  allocate(RHS2(1:nx,1:ny,1:nz))


  if (projectiontype==2) then
   Apr=coef*dt
   !Subtract pressure gradient terms in case of pressure correction method no. 2
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
          U(i,j,k)=U(i,j,k)+Apr*(Pr(i+1,j,k)-Pr(i,j,k))/dxU(i)
     enddo
    enddo
   enddo
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
          V(i,j,k)=V(i,j,k)+Apr*(Pr(i,j+1,k)-Pr(i,j,k))/dyV(j)
     enddo
    enddo
   enddo
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
          W(i,j,k)=W(i,j,k)+Apr*(Pr(i,j,k+1)-Pr(i,j,k))/dzW(k)
     enddo
    enddo
   enddo
  endif


   dt2=coef*dt
   if (Re>0)  then
    dt3=coef*dt/(2._TIM*Re)
   else
    dt3=0
   endif

   call BoundU(1,U)
   call BoundU(2,V)
   call BoundU(3,W)

   S=0
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
      S=S-((U(i,j,k)-U(i-1,j,k))/(dxPr(i))+(V(i,j,k)-V(i,j-1,k))/(dyPr(j))+(W(i,j,k)-W(i,j,k-1))/(dzPr(k)))*dxPr(i)*dyPr(j)*dzPr(k)
     enddo
    enddo
   enddo

   S=S/(ly*lz)
   write(*,*) "Compatibility correction",S
   U(Unx+1,:,:)=U(Unx+1,:,:)+S

   S=0
   S2=0
   maxI=-5
   maxJ=-5
   maxK=-5


   do k=1,Prnz            !divergence of U -> RHS
    do j=1,Prny
     do i=1,Prnx
        if (present(Q)) then
          RHS(i,j,k)=(U(i,j,k)-U(i-1,j,k))/(dxPr(i))+(V(i,j,k)-V(i,j-1,k))/(dyPr(j))+(W(i,j,k)-W(i,j,k-1))/(dzPr(k))-Q(i,j,k)
          S2=S2+abs(RHS(i,j,k)*dxPr(i)*dyPr(j)*dzPr(k))
          RHS(i,j,k)=RHS(i,j,k)/(dt2)
       else
          RHS(i,j,k)=(U(i,j,k)-U(i-1,j,k))/(dxPr(i))+(V(i,j,k)-V(i,j-1,k))/(dyPr(j))+(W(i,j,k)-W(i,j,k-1))/(dzPr(k))
          S2=S2+abs(RHS(i,j,k)*dxPr(i)*dyPr(j)*dzPr(k))
          RHS(i,j,k)=RHS(i,j,k)/(dt2)
       endif
     enddo
    enddo
   enddo
   write (*,*) "avgRHS",S2/(lx*ly*lz)


   if (poissmet==1) then

       call POISSSOR(Phi,RHS)

   elseif (poissmet==2) then

       call POISSFISH(Phi,RHS)

   elseif (poissmet==3) then

       RHS2=RHS
       call POISSFISH(Phi,RHS)

       Phiref=Phi(Prnx/2,Prny/2,Prnz/2)
       do k=1,Prnz
        do j=1,Prny
         do i=1,Prnx
          Phi(i,j,k)=Phi(i,j,k)-Phiref
         enddo
        enddo
       enddo

       RHS=RHS2
       if (Prny==1) then
        call POISSMG2d(Phi,RHS)
       else
        call POISSMG(Phi,RHS)
       endif

   elseif (poissmet==4) then

       if (Prny==1) then
        call POISSMG2d(Phi,RHS)
       else
        call POISSMG(Phi,RHS)
       endif

   elseif (poissmet==5) then

       write(*,*) "MUDPACK support disabled."
   endif



   Phiref=Phi(Prnx/2,Prny/2,Prnz/2)
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
      Phi(i,j,k)=Phi(i,j,k)-Phiref
     enddo
    enddo
   enddo


   call Bound_Phi(Phi)




 do k=1,Unz
  do j=1,Uny
    do i=1,Unx
      U(i,j,k)=U(i,j,k)-dt2*(Phi(i+1,j,k)-Phi(i,j,k))/dxU(i)
    enddo
   enddo
  enddo
  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
      V(i,j,k)=V(i,j,k)-dt2*(Phi(i,j+1,k)-Phi(i,j,k))/dyV(j)
    enddo
   enddo
  enddo
  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
      W(i,j,k)=W(i,j,k)-dt2*(Phi(i,j,k+1)-Phi(i,j,k))/dzW(k)
    enddo
   enddo
  enddo
  if (projectiontype==1) then
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
       Pr(i,j,k)=Pr(i,j,k)+Phi(i,j,k)-dt3*(((Phi(i+1,j,k)-Phi(i,j,k))/dxU(i)-(Phi(i,j,k)-Phi(i-1,j,k))/dxU(i-1))/dxPr(i)+&
                                            ((Phi(i,j+1,k)-Phi(i,j,k))/dyV(j)-(Phi(i,j,k)-Phi(i,j-1,k))/dyV(j-1))/dyPr(j)+&
                                            ((Phi(i,j,k+1)-Phi(i,j,k))/dzW(k)-(Phi(i,j,k)-Phi(i,j,k-1))/dzW(k-1))/dzPr(k))
     enddo
    enddo
   enddo
  elseif (projectiontype==2) then
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
       Pr(i,j,k)=Phi(i,j,k)
     enddo
    enddo
   enddo
  endif


  S=0
  S2=0
  maxI=-5
  maxJ=-5
  maxK=-5
  call BoundU(1,U)
  call BoundU(2,V)
  call BoundU(3,W)
  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
        if (present(Q)) then
          RHS(i,j,k)=(U(i,j,k)-U(i-1,j,k))/(dxPr(i))+(V(i,j,k)-V(i,j-1,k))/(dyPr(j))+(W(i,j,k)-W(i,j,k-1))/(dzPr(k))-Q(i,j,k)
          S2=S2+abs(RHS(i,j,k)*dxPr(i)*dyPr(j)*dzPr(k))
          if (abs(RHS(i,j,k))>abs(S).and.i>1) then
                                      S=RHS(i,j,k)
                                      maxI=i
                                      maxJ=j
                                      maxK=k
                                     endif
          RHS(i,j,k)=RHS(i,j,k)/(dt2)
        else
          RHS(i,j,k)=(U(i,j,k)-U(i-1,j,k))/(dxPr(i))+(V(i,j,k)-V(i,j-1,k))/(dyPr(j))+(W(i,j,k)-W(i,j,k-1))/(dzPr(k))
          S2=S2+abs(RHS(i,j,k)*dxPr(i)*dyPr(j)*dzPr(k))
          if (abs(RHS(i,j,k))>abs(S).and.i>1) then
                                      S=RHS(i,j,k)
                                      maxI=i
                                      maxJ=j
                                      maxK=k
                                     endif
          RHS(i,j,k)=RHS(i,j,k)/(dt2)
        endif
    enddo
   enddo
  enddo

  if ((maxI>-1).and.(maxJ>-1)) write (*,*) "maxRHS",maxI,maxJ,maxK,S
  write (*,*) "avgRHS",S2/(lx*ly*lz)


!     GOTO 30
!   OPEN(11,file="Phi.vtk")
!   write (11,"(A)") "# vtk DataFile Version 2.0"
!   write (11,"(A)") "diplomka output file"
!   write (11,"(A)") "ASCII"
!   write (11,"(A)") "DATASET RECTILINEAR_GRID"
!   str="DIMENSIONS"
!   write (str(12:),*) Prnx,Prny,Prnz
!   write (11,"(A)") str
!   str="X_COORDINATES"
!   write (str(15:),*) Prnx,"float"
!   write (11,"(A)") str
!   write (11,*) xPr(1:Prnx)
!   str="Y_COORDINATES"
!   write (str(15:),*) Prny,"float"
!   write (11,"(A)") str
!   write (11,*) yPr(1:Prny)
!   str="Z_COORDINATES"
!   write (str(15:),*) Prnz,"float"
!   write (11,"(A)") str
!   write (11,*) zPr(1:Prnz)
!   str="POINT_DATA"
!   write (str(12:),*) Prnx*Prny*Prnz
!   write (11,"(A)") str
!
!
!   write (11,"(A)") "SCALARS RHS float"
!   write (11,"(A)") "LOOKUP_TABLE default"
!   do k=1,Prnz
!    do j=1,Prny
!     do i=1,Prnx
!       Write (11,*) Phi(i,j,k)
!     enddo
!    enddo
!   enddo
!   write (11,*)
!   CLOSE(11)
!   20 continue
!   OPEN(11,file="RHS2.vtk")
!   write (11,"(A)") "# vtk DataFile Version 2.0"
!   write (11,"(A)") "diplomka output file"
!   write (11,"(A)") "ASCII"
!   write (11,"(A)") "DATASET RECTILINEAR_GRID"
!   str="DIMENSIONS"
!   write (str(12:),*) Prnx,Prny,Prnz
!   write (11,"(A)") str
!   str="X_COORDINATES"
!   write (str(15:),*) Prnx,"float"
!   write (11,"(A)") str
!   write (11,*) xPr(1:Prnx)
!   str="Y_COORDINATES"
!   write (str(15:),*) Prny,"float"
!   write (11,"(A)") str
!   write (11,*) yPr(1:Prny)
!   str="Z_COORDINATES"
!   write (str(15:),*) Prnz,"float"
!   write (11,"(A)") str
!   write (11,*) zPr(1:Prnz)
!   str="POINT_DATA"
!   write (str(12:),*) Prnx*Prny*Prnz
!   write (11,"(A)") str
!
!
!   write (11,"(A)") "SCALARS RHS float"
!   write (11,"(A)") "LOOKUP_TABLE default"
!   do k=1,Prnz
!    do j=1,Prny
!     do i=1,Prnx
!       Write (11,*) coef*RHS(i,j,k)*dt*3./(1./dxmin+1./dymin+1./dzmin)
!     enddo
!    enddo
!   enddo
!   write (11,*)
!   CLOSE(11)
!   30 continue
  endsubroutine PR_CORRECT



  subroutine POISSSOR(Phi,RHS)        !Solves Poisson equation using Successive over-relaxation

  real(KND),dimension(0:,0:,0:),intent(inout):: Phi
  real(KND),dimension(1:,1:,1:),intent(in)::RHS
  integer,save:: called=0
  integer nx,ny,nz,i,j,k,l,blockit
  real(KND) S,P,Ap
  real(KND),dimension(:),allocatable,save::Aw,Ae,As,An,Ab,At


   write (*,*) "Computing poisson equation"
   S=0
   nx=Prnx
   ny=Prny
   nz=Prnz
   if (called==0) then                             !coefficients based on grid spacing computed only once
       allocate(Aw(1:nx),Ae(1:nx))
       allocate(As(1:ny),An(1:ny))
       allocate(Ab(1:nz),At(1:nz))
       forall(i=1:nx)
        Ae(i)=1._KND/(dxU(i)*dxPr(i))
        Aw(i)=1._KND/(dxU(i-1)*dxPr(i))
       endforall
       forall(j=1:ny)
        An(j)=1._KND/(dyV(j)*dyPr(j))
        As(j)=1._KND/(dyV(j-1)*dyPr(j))
       endforall
       forall(k=1:nz)
        At(k)=1._KND/(dzW(k)*dzPr(k))
        Ab(k)=1._KND/(dzW(k-1)*dzPr(k))
       endforall
       called=1
   endif

   if (GPU>0) then
     l=1
     blockit=200
     S=huge(1.0)
     if (blockit>maxPoissoniter) blockit=maxPoissoniter
     !$hmpp <SOR> allocate
     !$hmpp <SOR> advancedload, args[GS_GPU::Phi,GS_GPU::RHS,GS_GPU::nx,GS_GPU::ny,GS_GPU::nz,GS_GPU::nit,&
     !$hmpp <SOR>    GS_GPU::Aw,GS_GPU::Ae,GS_GPU::As,GS_GPU::An,GS_GPU::Ab,GS_GPU::At,&
     !$hmpp <SOR>    GS_GPU::Btype]
     do while (l<=maxPoissoniter.and.S>epsPoisson)
      !$hmpp  <SOR> GS_GPU callsite
      call GS_GPU(nx,ny,nz,blockit,&
                      Phi,RHS,&
                      Aw,Ae,As,An,Ab,At,&
                      Btype)
      !$hmpp  <SOR> Res_GPU callsite
      call Res_GPU(nx,ny,nz,&
                      Phi,RHS,&
                      Aw,Ae,As,An,Ab,At,&
                      Btype,S)
      !$hmpp <SOR> delegatedStore,Args[Res_GPU::R]
      l=l+blockit
      write (*,*) "GPU Poisson iter: ",l,S
     enddo
      !$hmpp <SOR> delegatedStore,Args[Res_GPU::Phi]
      !$hmpp <SOR> release
   else
    l=1
    S=huge(1.0_KND)
        do while (l<=maxPoissoniter.and.S>maxPoissoniter)
          S=0
          !$OMP PARALLEL PRIVATE(i,j,k,p) REDUCTION(max:S)
          !$OMP DO
          do k=1,nz
             do j=1,ny
                do i=1+mod(j+k,2),nx,2
                  p=0
                  Ap=0
                  if (i>1) then
                            p=p+Phi(i-1,j,k)*Aw(i)
                            Ap=Ap+Aw(i)
                  elseif (Btype(We)==PERIODIC) then
                            p=p+Phi(nx,j,k)*Aw(i)
                            Ap=Ap+Aw(i)
                  endif
                  if (i<nx) then
                            p=p+Phi(i+1,j,k)*Ae(i)
                            Ap=Ap+Ae(i)
                  elseif (Btype(We)==PERIODIC) then
                            p=p+Phi(1,j,k)*Ae(i)
                            Ap=Ap+Ae(i)
                  endif
                  if (j>1) then
                            p=p+Phi(i,j-1,k)*As(j)
                            Ap=Ap+As(j)
                  elseif (Btype(No)==PERIODIC) then
                            p=p+Phi(i,ny,k)*As(j)
                            Ap=Ap+As(j)
                  endif
                  if (j<ny) then
                            p=p+Phi(i,j+1,k)*An(j)
                            Ap=Ap+An(j)
                  elseif (Btype(No)==PERIODIC) then
                            p=p+Phi(i,1,k)*An(j)
                            Ap=Ap+An(j)
                  endif
                  if (k>1) then
                            p=p+Phi(i,j,k-1)*Ab(k)
                            Ap=Ap+Ab(k)
                  elseif (Btype(To)==PERIODIC) then
                            p=p+Phi(i,j,nz)*Ab(k)
                            Ap=Ap+Ab(k)
                  endif
                  if (k<nz) then
                            p=p+Phi(i,j,k+1)*At(k)
                            Ap=Ap+At(k)
                  elseif (Btype(To)==PERIODIC) then
                            p=p+Phi(i,j,1)*At(k)
                            Ap=Ap+At(k)
                  endif
                  p=p-RHS(i,j,k)
                  p=p/Ap
                  S=max(abs(p-Phi(i,j,k)),S)
                  Phi(i,j,k)=p
                enddo
             enddo
          enddo
          !$OMP ENDDO
          !$OMP DO
          do k=1,nz
             do j=1,ny
                do i=1+mod(j+k+1,2),nx,2
                  p=0
                  Ap=0
                  if (i>1) then
                            p=p+Phi(i-1,j,k)*Aw(i)
                            Ap=Ap+Aw(i)
                  elseif (Btype(We)==PERIODIC) then
                            p=p+Phi(nx,j,k)*Aw(i)
                            Ap=Ap+Aw(i)
                  endif
                  if (i<nx) then
                            p=p+Phi(i+1,j,k)*Ae(i)
                            Ap=Ap+Ae(i)
                  elseif (Btype(We)==PERIODIC) then
                            p=p+Phi(1,j,k)*Ae(i)
                            Ap=Ap+Ae(i)
                  endif
                  if (j>1) then
                            p=p+Phi(i,j-1,k)*As(j)
                            Ap=Ap+As(j)
                  elseif (Btype(No)==PERIODIC) then
                            p=p+Phi(i,ny,k)*As(j)
                            Ap=Ap+As(j)
                  endif
                  if (j<ny) then
                            p=p+Phi(i,j+1,k)*An(j)
                            Ap=Ap+An(j)
                  elseif (Btype(No)==PERIODIC) then
                            p=p+Phi(i,1,k)*An(j)
                            Ap=Ap+An(j)
                  endif
                  if (k>1) then
                            p=p+Phi(i,j,k-1)*Ab(k)
                            Ap=Ap+Ab(k)
                  elseif (Btype(To)==PERIODIC) then
                            p=p+Phi(i,j,nz)*Ab(k)
                            Ap=Ap+Ab(k)
                  endif
                  if (k<nz) then
                            p=p+Phi(i,j,k+1)*At(k)
                            Ap=Ap+At(k)
                  elseif (Btype(To)==PERIODIC) then
                            p=p+Phi(i,j,1)*At(k)
                            Ap=Ap+At(k)
                  endif
                  p=p-RHS(i,j,k)
                  p=p/Ap
                  S=max(abs(p-Phi(i,j,k)),S)
                  Phi(i,j,k)=p
                enddo
             enddo
          enddo
          !$OMP ENDDO
        !$OMP ENDPARALLEL
          p=abs(maxval(Phi(1:nx,1:ny,1:nz)))
          if (p>0) S=S/p
          if (MOD(l,10)==0)  write (*,*) "   Poisson iter: ",l,S
        enddo
   endif
 end subroutine POISSSOR










  subroutine  POISSVEL(Phi,RHS)      !Calls Bicgstab solver from LAPACK, does not work

  real(KND),dimension(0:,0:,0:),intent(inout):: Phi
  real(KND),dimension(1:,1:,1:),intent(in)::RHS
  real(KND),dimension(1:Prnx+2,1:Prny+2,1:Prnz+2)::RHS2,Phi2
  real(KND),dimension(1:Prnx+2)::x
  real(KND),dimension(1:Prny+2)::y
  real(KND),dimension(1:Prnz+2)::z
  integer nx,ny,nz,i0,j0,k0,itmax,itol,np
  real(KND) tol,err
  logical,dimension(2):: dbcx,dbcy,dbcz

  nx=Prnx+2
  ny=Prny+2
  nz=Prnz+2
  np=1
  rhs2=0
  rhs2(2:Prnx+1,2:Prny+1,2:Prnz+1)=RHS(1:Prnx,1:Prny,1:Prnz)
  x=xPR(0:Prnx+1)
  y=yPR(0:Prny+1)
  z=zPR(0:Prnz+1)
  dbcx=.false.
  dbcy=.false.
  dbcz=.false.
  i0=1
  j0=1
  k0=1
  Phi2=Phi(0:Prnx+1,0:Prny+1,0:Prnz+1)
  itmax=1000000
  itol=1
  tol=5.e-5


    !call  poisson2(nx,ny,nz,np,x,y,z,rhs2,dbcx,dbcy,dbcz,i0,j0,k0,phi2,itmax,itol,tol,err)

  Phi(0:Prnx+1,0:Prny+1,0:Prnz+1)=Phi2
  endsubroutine POISSVEL




  subroutine POISSADI(Phi,RHS)                      !Alternating direction implicit with tridiag solver from LINPACK
  real(KND),dimension(0:,0:,0:),intent(inout):: Phi !does not converge
  real(KND),dimension(1:,1:,1:),intent(in)::RHS
  real(KND) Phi2(0:Prnx+2,0:Prny+2,0:Prnz+2)
  real(KND),parameter:: adirelaxpar=10._KND
  real(KND),allocatable,dimension(:),save:: bx,cx,dx,ex
  real(KND),allocatable,dimension(:),save:: by,cy,dy,ey
  real(KND),allocatable,dimension(:),save:: bz,cz,dz,ez
  integer i,j,k,l,info,nx,ny,nz
  real(KND),dimension(:),allocatable,save::Apx,Apy,Apz,Aw,Ae,As,An,Ab,At
  integer,save:: called=0
  real(KND) S


  nx=Prnx
  ny=Prny
  nz=Prnz
  if (called==0) then
       allocate(Aw(1:nx),Ae(1:nx),Apx(1:nx))
       allocate(As(1:ny),An(1:ny),Apy(1:ny))
       allocate(Ab(1:nz),At(1:nz),Apz(1:nz))
       allocate(bx(1:nx+2),cx(1:nx+2),dx(1:nx+2),ex(1:nx+2))
       allocate(by(1:ny+2),cy(1:ny+2),dy(1:ny+2),ey(1:ny+2))
       allocate(bz(1:nz+2),cz(1:nz+2),dz(1:nz+2),ez(1:nz+2))
       forall(i=1:nx)
        Ae(i)=1._KND/(dxU(i)*dxPr(i))
        Aw(i)=1._KND/(dxU(i-1)*dxPr(i))
        Apx(i)=(1._KND/dxU(i)+1._KND/dxU(i-1))/dxPr(i)
       endforall
       forall(j=1:ny)
        An(j)=1._KND/(dyV(j)*dyPr(j))
        As(j)=1._KND/(dyV(j-1)*dyPr(j))
        Apy(j)=(1._KND/dyV(j)+1._KND/dyV(j-1))/dyPr(j)
       endforall
       forall(k=1:nz)
        At(k)=1._KND/(dzW(k)*dzPr(k))
        Ab(k)=1._KND/(dzW(k-1)*dzPr(k))
        Apz(k)=(1._KND/dzW(k)+1._KND/dzW(k-1))/dzPr(k)
       endforall
       called=1

  endif
!  call POISSSOR(Phi,RHS)
  do l=1,maxpoissoniter
  S=0
  Phi2=Phi
   call Bound_Phi(Phi)
   do k=1,nz
    do j=1,ny
       forall(i=2:nx+1)
        dx(i)=-Apx(i-1)-adirelaxpar
       endforall
       forall(i=2:nx)
        cx(i)=Aw(i)
       endforall
       forall(i=3:nx+1)
        ex(i)=Ae(i-2)
       endforall
       if (Btype(We)==PERIODIC) then
        dx(1)=1
        dx(nx+2)=1
        ex(2)=0
        cx(nx+1)=0
       else
        dx(1)=1
        dx(nx+2)=1
        ex(2)=-1
        cx(nx+1)=-1
       endif
     if (Btype(We)==PERIODIC) then
      bx(1)=Phi(nx,j,k)
      bx(nx+2)=Phi(1,j,k)
     else
      bx(1)=0
      bx(nx+2)=0
     endif
     do i=1,nx
      bx(i+1)=RHS(i,j,k)+(Apy(j)+Apz(k)-adirelaxpar)*Phi(i,j,k)-(Phi(i,j+1,k)*An(j)+Phi(i,j-1,k)*As(j)+&
                      Phi(i,j,k+1)*At(k)+Phi(i,j,k-1)*Ab(k))
     enddo
!     call SGTSL(nx+2,cx,dx,ex,bx,info)
     if (info/=0)  write (*,*) info

     Phi(1:nx,j,k)=bx(2:nx+1)
    enddo
   enddo

   call Bound_Phi(Phi)
   do k=1,nz
    do i=1,nx
       forall(j=2:ny+1)
        dy(j)=-Apy(j-1)-adirelaxpar
       endforall
       forall(j=2:ny)
        cy(j)=As(j)
       endforall
       forall(j=3:ny+1)
        ey(j)=An(j-2)
       endforall
       if (Btype(So)==PERIODIC) then
        dy(1)=1
        dy(ny+2)=1
        ey(2)=0
        cy(ny+1)=0
       else
        dy(1)=1
        dy(ny+2)=1
        ey(2)=-1
        cy(ny+1)=-1
       endif
     if (Btype(So)==PERIODIC) then
      by(1)=Phi(i,ny,k)
      by(ny+2)=Phi(i,1,k)
     else
      by(1)=0
      by(ny+2)=0
     endif
     do j=1,ny
      by(j+1)=(Apy(j)-adirelaxpar)*Phi(i,j,k)+(Phi(i,j+1,k)*An(j)+Phi(i,j-1,k)*As(j))
     enddo
     ! call SGTSL(ny+2,cy,dy,ey,by,info)

     if (info/=0)  write (*,*) "info",info
    Phi(i,1:ny,k)=by(2:ny+1)
    enddo
   enddo


   call Bound_Phi(Phi)
   do j=1,ny
    do i=1,nx
       forall(k=2:nz+1)
        dz(k)=-Apz(k-1)-adirelaxpar
       endforall
       forall(k=2:nz)
        cz(k)=Ab(k)
       endforall
       forall(k=3:nz+1)
        ez(k)=At(k-2)
       endforall
       if (Btype(Bo)==PERIODIC) then
        dz(1)=1
        dz(nz+2)=1
        ez(2)=0
        cz(nz+1)=0
       else
        dz(1)=1
        dz(nz+2)=1
        ez(2)=-1
        cz(nz+1)=-1
       endif
     if (Btype(Bo)==PERIODIC) then
      bz(1)=Phi(i,j,nz)
      bz(nz+2)=Phi(i,j,1)
     else
      bz(1)=0
      bz(nz+2)=0
     endif
     do k=1,nz
      bz(k+1)=(Apz(k)-adirelaxpar)*Phi(i,j,k)+(Phi(i,j,k+1)*At(k)+Phi(i,j,k-1)*Ab(k))
     enddo
!     call SGTSL(nz+2,cz,dz,ez,bz,info)

     if (info/=0) write (*,*) "info",info
     Phi(i,j,1:nz)=bz(2:nz+1)
    enddo
   enddo
   S=SUM(Abs(Phi(1:Prnx,1:Prny,1:Prnz)-Phi2(1:Prnx,1:Prny,1:Prnz)))/(Prnx*Prny*Prnz)
   write (*,*) l,S
   if (Abs(S)<epspoisson) exit
  enddo

  endsubroutine POISSADI







  subroutine POISSADIGS(Phi,RHS)                    !Alternating direction implicit with Gauss-Seidel
  real(KND),dimension(0:,0:,0:),intent(inout):: Phi !does not converge
  real(KND),dimension(1:,1:,1:),intent(in)::RHS
  real(KND) Phi2(0:Prnx+2,0:Prny+2,0:Prnz+2)
  real(KND):: adirelaxpar=0.78_KND
  real(KND),allocatable,dimension(:),save:: bx,cx,dx,ex
  real(KND),allocatable,dimension(:),save:: by,cy,dy,ey
  real(KND),allocatable,dimension(:),save:: bz,cz,dz,ez
  integer i,j,k,l,nx,ny,nz
  real(KND),dimension(:),allocatable,save::Apx,Apy,Apz,Aw,Ae,As,An,Ab,At
  integer,save:: called=0
  real(KND) S


  nx=Prnx
  ny=Prny
  nz=Prnz
  if (called==0) then
       allocate(Aw(1:nx),Ae(1:nx),Apx(1:nx))
       allocate(As(1:ny),An(1:ny),Apy(1:ny))
       allocate(Ab(1:nz),At(1:nz),Apz(1:nz))
       allocate(bx(1:nx+2),cx(1:nx+2),dx(1:nx+2),ex(1:nx+2))
       allocate(by(1:ny+2),cy(1:ny+2),dy(1:ny+2),ey(1:ny+2))
       allocate(bz(1:nz+2),cz(1:nz+2),dz(1:nz+2),ez(1:nz+2))
       forall(i=1:nx)
        Ae(i)=1._KND/(dxU(i)*dxPr(i))
        Aw(i)=1._KND/(dxU(i-1)*dxPr(i))
        Apx(i)=-(1._KND/dxU(i)+1._KND/dxU(i-1))/dxPr(i)
       endforall
       forall(j=1:ny)
        An(j)=1._KND/(dyV(j)*dyPr(j))
        As(j)=1._KND/(dyV(j-1)*dyPr(j))
        Apy(j)=-(1._KND/dyV(j)+1._KND/dyV(j-1))/dyPr(j)
       endforall
       forall(k=1:nz)
        At(k)=1._KND/(dzW(k)*dzPr(k))
        Ab(k)=1._KND/(dzW(k-1)*dzPr(k))
        Apz(k)=-(1._KND/dzW(k)+1._KND/dzW(k-1))/dzPr(k)
       endforall
       called=1

  endif
 ! call POISSSOR(Phi,RHS)
  write (*,*) "---------"
  do l=1,maxpoissoniter
  write (*,*) "l=",l
  S=0
   Phi2=Phi
   call Bound_Phi(Phi)
   do k=1,nz
    do j=1,ny
     do i=1,nx
      cx(i)=2*RHS(i,j,k)-(Apx(i)+2*Apy(j)+2*Apz(k)+adirelaxpar)*Phi(i,j,k)-(Phi(i+1,j,k)*Ae(i)+Phi(i-1,j,k)*Aw(i)+&
                                                                     2*Phi(i,j+1,k)*An(j)+2*Phi(i,j-1,k)*As(j)+&
                                                                     2*Phi(i,j,k+1)*At(k)+2*Phi(i,j,k-1)*Ab(k))
     enddo
     bx(2:nx+1)=cx(1:nx)
     dx(1:nx)=Apx(1:nx)-adirelaxpar
!     write (*,*) "-",k,j
     call TRIDAGGS(nx,Btype(We),Btype(Ea),Aw(1:nx),Ae(1:nx),dx(1:nx),bx(1:nx+2),cx(1:nx))
      !write(*,*) "***"
     !write(*,*) bx
     S=S+sum(abs(Phi(1:nx,j,k)-bx(2:nx+1)))/(Prnx*Prny*Prnz)
     Phi2(1:nx,j,k)=bx(2:nx+1)
    enddo
   enddo
   write(*,*)"x", S

   S=0
   Phi=Phi2
   call Bound_Phi(Phi)
   do k=1,nz
    do i=1,nx
     do j=1,ny
      cy(j)=(Apy(j)-adirelaxpar)*Phi(i,j,k)+(Phi(i,j+1,k)*An(j)+Phi(i,j-1,k)*As(j))
     enddo
     by(2:ny+1)=cy(1:ny)
     dy(1:ny)=Apy(1:ny)-adirelaxpar

     call TRIDAGGS(ny,Btype(So),Btype(No),As(1:ny),An(1:ny),dy(1:ny),by(1:ny+2),cy(1:ny))

     S=S+sum(abs(Phi(i,1:ny,k)-by(2:ny+1)))/(Prnx*Prny*Prnz)
     Phi2(i,1:ny,k)=by(2:ny+1)
    enddo
   enddo
   write(*,*)"y", S


   Phi=Phi2
   call Bound_Phi(Phi)
   do j=1,ny
    do i=1,nx
     do k=1,nz
      cz(k)=(Apz(k)-adirelaxpar)*Phi(i,j,k)+(Phi(i,j,k+1)*At(k)+Phi(i,j,k-1)*Ab(k))
     enddo
     bz(2:nz+1)=cz(1:nz)
     dz(1:nz)=Apz(1:nz)-adirelaxpar

     call TRIDAGGS(nz,Btype(Bo),Btype(To),Ab(1:nz),At(1:nz),dz(1:nz),bz(1:nz+2),cz(1:nz))

     S=S+sum(abs(Phi(i,j,1:nz)-bz(2:nz+1)))/(Prnx*Prny*Prnz)
     Phi2(i,j,1:nz)=bz(2:nz+1)
    enddo
   enddo
   Phi=Phi2
   write (*,*) l,S

  enddo

  endsubroutine POISSADIGS


  subroutine TRIDAGGS(n,btw,bte,aw,ae,ap,Phi,RHS)  !Solves tridiagonal system using Gauss Seidel
  !tridiagonal system
  integer,intent(in):: n,btw,bte
  real(KND),dimension(1:),intent(in)::aw,ap,ae,RHS
  real(KND),dimension(0:),intent(inout):: Phi
  integer,parameter:: tridagiter=1000
  real(KND),parameter::eps=1e-5
  integer i,j
  real(KND) S,p
  do j=1,tridagiter
   S=0
   if (btw==PERIODIC) then
               Phi(0)=Phi(n)
   else
               Phi(0)=Phi(1)
   endif
   if (bte==PERIODIC) then
               Phi(n+1)=Phi(1)
   else
               Phi(n+1)=Phi(n)
   endif
   do i=1,n
    p=(RHS(i)-aw(i)*Phi(i-1)-ae(i)*Phi(i+1))
    p=p/ap(i)
    S=S+(p-Phi(i))**2
    p=Phi(i)
   enddo
  ! write (*,*) "tridag",j,S
   if (S<eps) exit
  enddo
  endsubroutine TRIDAGGS















!GPU Codelets. More or less like externals. Inside the module only because of KND.

  !$hmpp <SOR> group, target=CUDA
  !$hmpp <SOR> mapbyname, nx,ny,nz,Phi,RHS,Aw,Ae,As,An,Ab,At,Btype

!$hmpp <SOR> GS_GPU codelet
subroutine GS_GPU(nx,ny,nz,nit,Phi,RHS,&
                     Aw,Ae,As,An,Ab,At,&
                     Btype)
   implicit none

   integer,parameter:: KND=4,PERIODIC=3
   integer, parameter :: Ea=1,We=2,So=3,No=4,Bo=5,To=6

   integer,intent(in)   :: nx,ny,nz,nit
   real(KND),dimension(0:nx+1,0:ny+1,0:nz+1),intent(inout)::Phi
   real(KND),dimension(1:nx,1:ny,1:nz),intent(in)::RHS
   real(KND),dimension(1:nx),intent(in) :: Aw,Ae
   real(KND),dimension(1:ny),intent(in) :: As,An
   real(KND),dimension(1:nz),intent(in) :: Ab,At
   integer(KND),intent(in) :: Btype(6)
   integer i,j,k,l
   real(KND) :: p,Ap
   intrinsic mod

  do l=1,nit
   !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
     do k=1,nz
        do i=1,nx
            do j=1+mod(i+k,2),ny,2
             p=0
             Ap=0
             if (i>1) then
                       p=p+Phi(i-1,j,k)*Aw(i)
                       Ap=Ap+Aw(i)
             elseif (Btype(We)==PERIODIC) then
                       p=p+Phi(nx-1,j,k)*Aw(i)
                       Ap=Ap+Aw(i)
             endif
             if (i<nx) then
                       p=p+Phi(i+1,j,k)*Ae(i)
                       Ap=Ap+Ae(i)
             elseif (Btype(Ea)==PERIODIC) then
                       p=p+Phi(1,j,k)*Ae(i)
                       Ap=Ap+Ae(i)
             endif
             if (j>1) then
                       p=p+Phi(i,j-1,k)*As(j)
                       Ap=Ap+As(j)
             elseif (Btype(So)==PERIODIC) then
                       p=p+Phi(i,ny-1,k)*As(j)
                       Ap=Ap+As(j)
             endif
             if (j<ny) then
                       p=p+Phi(i,j+1,k)*An(j)
                       Ap=Ap+An(j)
             elseif (Btype(No)==PERIODIC) then
                       p=p+Phi(i,1,k)*An(j)
                       Ap=Ap+An(j)
             endif
             if (k>1) then
                       p=p+Phi(i,j,k-1)*Ab(k)
                       Ap=Ap+Ab(k)
             elseif (Btype(Bo)==PERIODIC) then
                       p=p+Phi(i,j,nz-1)*Ab(k)
                       Ap=Ap+Ab(k)
             endif
             if (k<nz) then
                       p=p+Phi(i,j,k+1)*At(k)
                       Ap=Ap+At(k)
             elseif (Btype(To)==PERIODIC) then
                       p=p+Phi(i,j,1)*At(k)
                       Ap=Ap+At(k)
             endif
             p=p-RHS(i,j,k)

             p=p/Ap
             Phi(i,j,k)=p
            enddo
        enddo
    enddo
  !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
     do k=1,nz
        do i=1,nx
            do j=1+mod(i+k+1,2),ny,2
             p=0
             Ap=0
             if (i>1) then
                       p=p+Phi(i-1,j,k)*Aw(i)
                       Ap=Ap+Aw(i)
             elseif (Btype(We)==PERIODIC) then
                       p=p+Phi(nx-1,j,k)*Aw(i)
                       Ap=Ap+Aw(i)
             endif
             if (i<nx) then
                       p=p+Phi(i+1,j,k)*Ae(i)
                       Ap=Ap+Ae(i)
             elseif (Btype(Ea)==PERIODIC) then
                       p=p+Phi(1,j,k)*Ae(i)
                       Ap=Ap+Ae(i)
             endif
             if (j>1) then
                       p=p+Phi(i,j-1,k)*As(j)
                       Ap=Ap+As(j)
             elseif (Btype(So)==PERIODIC) then
                       p=p+Phi(i,ny-1,k)*As(j)
                       Ap=Ap+As(j)
             endif
             if (j<ny) then
                       p=p+Phi(i,j+1,k)*An(j)
                       Ap=Ap+An(j)
             elseif (Btype(No)==PERIODIC) then
                       p=p+Phi(i,1,k)*An(j)
                       Ap=Ap+An(j)
             endif
             if (k>1) then
                       p=p+Phi(i,j,k-1)*Ab(k)
                       Ap=Ap+Ab(k)
             elseif (Btype(Bo)==PERIODIC) then
                       p=p+Phi(i,j,nz-1)*Ab(k)
                       Ap=Ap+Ab(k)
             endif
             if (k<nz) then
                       p=p+Phi(i,j,k+1)*At(k)
                       Ap=Ap+At(k)
             elseif (Btype(To)==PERIODIC) then
                       p=p+Phi(i,j,1)*At(k)
                       Ap=Ap+At(k)
             endif
             p=p-RHS(i,j,k)

             p=p/Ap
             Phi(i,j,k)=p
            enddo
        enddo
    enddo
  enddo
endsubroutine GS_GPU


!$hmpp <SOR> Res_GPU codelet
subroutine Res_GPU(nx,ny,nz,Phi,RHS,&
                     Aw,Ae,As,An,Ab,At,&
                     Btype,R)

   implicit none
   intrinsic mod, abs, max


   integer,parameter:: KND=4,PERIODIC=3
   integer, parameter :: Ea=1,We=2,So=3,No=4,Bo=5,To=6

   integer,intent(in)   :: nx,ny,nz
   real(KND),dimension(0:nx+1,0:ny+1,0:nz+1),intent(inout)::Phi
   real(KND),dimension(1:nx,1:ny,1:nz),intent(in)::RHS
   real(KND),dimension(1:nx),intent(in) :: Aw,Ae
   real(KND),dimension(1:ny),intent(in) :: As,An
   real(KND),dimension(1:nz),intent(in) :: Ab,At
   integer(KND),intent(in) :: Btype(6)
   real(KND),intent(out) :: R
   integer i,j,k,l
   real(KND) :: p,Ap

   R=0
   !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i), reduce(max:R)
     do k=1,nz
        do i=1,nx
            do j=1,ny
             p=0
             Ap=0
             if (i>1) then
                       p=p+Phi(i-1,j,k)*Aw(i)
                       Ap=Ap+Aw(i)
             elseif (Btype(We)==PERIODIC) then
                       p=p+Phi(nx-1,j,k)*Aw(i)
                       Ap=Ap+Aw(i)
             endif
             if (i<nx) then
                       p=p+Phi(i+1,j,k)*Ae(i)
                       Ap=Ap+Ae(i)
             elseif (Btype(Ea)==PERIODIC) then
                       p=p+Phi(1,j,k)*Ae(i)
                       Ap=Ap+Ae(i)
             endif
             if (j>1) then
                       p=p+Phi(i,j-1,k)*As(j)
                       Ap=Ap+As(j)
             elseif (Btype(So)==PERIODIC) then
                       p=p+Phi(i,ny-1,k)*As(j)
                       Ap=Ap+As(j)
             endif
             if (j<ny) then
                       p=p+Phi(i,j+1,k)*An(j)
                       Ap=Ap+An(j)
             elseif (Btype(No)==PERIODIC) then
                       p=p+Phi(i,1,k)*An(j)
                       Ap=Ap+An(j)
             endif
             if (k>1) then
                       p=p+Phi(i,j,k-1)*Ab(k)
                       Ap=Ap+Ab(k)
             elseif (Btype(Bo)==PERIODIC) then
                       p=p+Phi(i,j,nz-1)*Ab(k)
                       Ap=Ap+Ab(k)
             endif
             if (k<nz) then
                       p=p+Phi(i,j,k+1)*At(k)
                       Ap=Ap+At(k)
             elseif (Btype(To)==PERIODIC) then
                       p=p+Phi(i,j,1)*At(k)
                       Ap=Ap+At(k)
             endif
             p=p-RHS(i,j,k)
             p=abs(-p +Ap*Phi(i,j,k))
             R=max(R,abs(p))
            enddo
        enddo
     enddo
endsubroutine Res_GPU




end module POISSON

module POISSON

 use PARAMETERS
 use BOUNDARIES
 use MULTIGRID,   only: POISSMG
 use MULTIGRID2d, only: POISSMG2d

implicit none

 integer MUDbnx,MUDbny,MUDbnz,MUDlmg

contains


subroutine PR_CORRECT(U,V,W,Pr,coef,Q)                   !Pressure correction
  real(KND),dimension(-2:,-2:,-2:),intent(inout):: U,V,W !Phi is computed in Poisson eq. with div of U in RHS
  real(KND),dimension(1:,1:,1:),intent(inout):: Pr       !Depending on active projection method Phi becomes new pressure
  real(KND),dimension(0:,0:,0:),optional,intent(in)::Q   !or is added to last pressure
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

   !$omp parallel do private(i,j,k) reduction(+:s2)
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
   !$omp end parallel do
   write (*,*) "avgRHS",S2/(lx*ly*lz)


   if (poissmet==1) then

       call POISSSOR(Phi,RHS)

   elseif (poissmet==2) then

       call POISS_POISFFT(Phi,RHS)

   elseif (poissmet==3) then

       if (Prny==1) then
        call POISSMG2d(Phi,RHS)
       else
        call POISSMG(Phi,RHS)
       endif

   endif



   Phiref=Phi(Prnx/2,Prny/2,Prnz/2)
   Phi = Phi-Phiref



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
  call BoundU(1,U)
  call BoundU(2,V)
  call BoundU(3,W)

!   !$omp parallel do private(i,j,k) reduction(max:S) reduction(+:S2)
!   do k=1,Prnz
!    do j=1,Prny
!     do i=1,Prnx
!         if (present(Q)) then
!           RHS(i,j,k)=(U(i,j,k)-U(i-1,j,k))/(dxPr(i))+(V(i,j,k)-V(i,j-1,k))/(dyPr(j))+(W(i,j,k)-W(i,j,k-1))/(dzPr(k))-Q(i,j,k)
!           S2=S2+abs(RHS(i,j,k)*dxPr(i)*dyPr(j)*dzPr(k))
!
!           if (i>2) then
!             S=max(abs(RHS(i,j,k)),S)
!           endif
!
!         else
!           RHS(i,j,k)=(U(i,j,k)-U(i-1,j,k))/(dxPr(i))+(V(i,j,k)-V(i,j-1,k))/(dyPr(j))+(W(i,j,k)-W(i,j,k-1))/(dzPr(k))
!           S2=S2+abs(RHS(i,j,k)*dxPr(i)*dyPr(j)*dzPr(k))
!
!           if (i>2) then
!             S=max(abs(RHS(i,j,k)),S)
!           endif
!
!         endif
!     enddo
!    enddo
!   enddo
!   !$omp end parallel do
!   write (*,*) "maxRHS",S
!   write (*,*) "avgRHS",S2/(lx*ly*lz)
!   if (present(Q)) write (*,*) "maxQ",maxval(Q)
  endsubroutine PR_CORRECT





  subroutine POISS_POISFFT(Phi,RHS)
    use poisfft
    type(PoisFFT_Solver3D),save :: Solver
    real(KND),dimension(0:,0:,0:),intent(inout):: Phi
    real(KND),dimension(1:,1:,1:),intent(in)::RHS
    integer i
    logical, save :: called = .false.
    integer bounds(6)
    integer(DBL), save :: trate
    integer(DBL)       :: t1, t2

    do i=1,6
      if (Btype(i)==PERIODIC) then
        bounds(i) = PoisFFT_PERIODIC
      else
        bounds(i) = PoisFFT_NeumannStag
      endif
    enddo

    if (.not.called) then
      Solver = PoisFFT_Solver3D_New(Prnx,Prny,Prnz,dxmin,dymin,dzmin,bounds)
      called = .true.
      call system_clock(count_rate=trate)
    endif

    call system_clock(count=t1)

    call PoisFFT_Solver3D_Execute(Solver,Phi,RHS)

    call system_clock(count=t2)
    write(*,*) "solver cpu time", real(t2-t1)/real(trate)

  endsubroutine POISS_POISFFT





  subroutine POISSSOR(Phi,RHS)        !Solves Poisson equation using Successive over-relaxation

  real(KND),dimension(0:,0:,0:),intent(inout):: Phi
  real(KND),dimension(1:,1:,1:),intent(in)::RHS
  integer,save:: called=0
  integer nx,ny,nz,i,j,k,l,blockit
  real(KND) S,P,Ap
  real(KND),dimension(:),allocatable,save::Aw,Ae,As,An,Ab,At


   write (*,*) "Computing Poisson equation"
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
    l=0
    S=huge(1.0_KND)
        do while (l<=maxPoissoniter.and.S>epsPoisson)
          l=l+1
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






!GPU Codelets. More or less like externals. Inside the module only because of KND.

  !$hmpp <SOR> group, target=CUDA
  !$hmpp <SOR> mapbyname, nx,ny,nz,Phi,RHS,Aw,Ae,As,An,Ab,At,Btype

!$hmpp <SOR> GS_GPU codelet
subroutine GS_GPU(nx,ny,nz,nit,Phi,RHS,&
                     Aw,Ae,As,An,Ab,At,&
                     Btype)
   implicit none
#ifdef __HMPP
   integer,parameter:: KND=4,PERIODIC=3
   integer, parameter :: Ea=1,We=2,So=3,No=4,Bo=5,To=6
#endif

   integer,intent(in)   :: nx,ny,nz,nit
   real(KND),dimension(0:nx+1,0:ny+1,0:nz+1),intent(inout)::Phi
   real(KND),dimension(1:nx,1:ny,1:nz),intent(in)::RHS
   real(KND),dimension(1:nx),intent(in) :: Aw,Ae
   real(KND),dimension(1:ny),intent(in) :: As,An
   real(KND),dimension(1:nz),intent(in) :: Ab,At
   integer,intent(in) :: Btype(6)
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

#ifdef __HMPP
   integer,parameter:: KND=4,PERIODIC=3
   integer, parameter :: Ea=1,We=2,So=3,No=4,Bo=5,To=6
#endif

   integer,intent(in)   :: nx,ny,nz
   real(KND),dimension(0:nx+1,0:ny+1,0:nz+1),intent(inout)::Phi
   real(KND),dimension(1:nx,1:ny,1:nz),intent(in)::RHS
   real(KND),dimension(1:nx),intent(in) :: Aw,Ae
   real(KND),dimension(1:ny),intent(in) :: As,An
   real(KND),dimension(1:nz),intent(in) :: Ab,At
   integer,intent(in) :: Btype(6)
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

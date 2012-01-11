module MULTIGRID

 use PARAMETERS

 implicit none

 private

 public POISSMG, SetMGParams

  type TMGArr
    real(KND),allocatable,dimension(:,:,:):: Arr
  endtype

  type TCoefs
   real(KND) dx,dy,dz,Aw,Ae,As,An,Ab,At
   integer nx,ny,nz
  endtype


  type(TMGArr),allocatable,dimension(:):: PhiMG,RHSMG,ResMG
  type(TCoefs),allocatable,dimension(:):: CoefMg

  integer bnx,bny,bnz !Base cube dimensions for Multigrid, Prnx=bnx*2**level
  integer LMG !depth of multigrid
  integer minmglevel !innermost MG level
  integer,save :: minGPUlevel=5
  integer mgncgc !type of cycling 1..V cycle, 2..W cycle
  integer mgnpre
  integer mgnpost
  integer mgmaxinnerGSiter
  real(KND) mgepsinnerGS

  integer,dimension(0:8) :: nxa,nya,nza
  real(KND),dimension(0:8) :: Aw,Ae,As,An,Ab,At

contains


 subroutine SetMGParams(llmg, lminmglevel, lmingpulevel, lbnx, lbny, lbnz,&
                          lmgncgc, lmgnpre, lmgnpost, lmgmaxinnerGSiter, lmgepsinnerGS)

  integer,intent(in)::   llmg, lminmglevel, lmingpulevel, lbnx, lbny, lbnz, lmgncgc, lmgnpre, lmgnpost, lmgmaxinnerGSiter
  real(KND),intent(in):: lmgepsinnerGS

  lmg=llmg
  minmglevel=lminmglevel
  minGPUlevel=lminGPUlevel
  bnx=lbnx
  bny=lbny
  bnz=lbnz
  mgncgc=lmgncgc
  mgnpre=lmgnpre
  mgnpost=lmgnpost
  mgmaxinnerGSiter=lmgmaxinnerGSiter
  mgepsinnerGS=lmgepsinnerGS
 endsubroutine SetMGParams

 subroutine Prolongate(AFine,ACoarse,level)
  integer,intent(in):: level
  real(KND),dimension(-1:,-1:,-1:),intent(in):: ACoarse
  real(KND),dimension(-1:,-1:,-1:),intent(inout):: AFine
  integer:: i,j,k,nx,ny,nz
  real(KND):: A1,A2,A4,A8


   nx=bnx*2**level !level means from which grid we interpolate
   ny=bny*2**level
   nz=bnz*2**level

   if (GPU>0.and.level>=minGPUlevel-1) then!.and.level>3


      !$hmpp <GSKernels> advancedload, args[Pr::level]
      !$hmpp <GSKernels> Pr callsite,args[*].noupdate=true
      call Pr_GPU(nxa,nya,nza,PhiMG(0)%Arr,PhiMG(1)%Arr,PhiMG(2)%Arr,PhiMG(3)%Arr,PhiMG(4)%Arr,&
                             PhiMG(5)%Arr,PhiMG(6)%Arr,PhiMG(7)%Arr,PhiMG(8)%Arr,&
                             Btype,level)

   else
    do k=0,nz
     do j=0,ny
      do i=0,nx
                      A1=ACoarse(i,j,k)
                      A2=0.5_KND*A1
                      A4=0.25_KND*A1
                      A8=0.125*A1
                      AFine(2*i,2*j,2*k)=AFine(2*i,2*j,2*k)+A1
                      AFine(2*i+1,2*j,2*k)=AFine(2*i+1,2*j,2*k)+A2
                      AFine(2*i-1,2*j,2*k)=AFine(2*i-1,2*j,2*k)+A2
                      AFine(2*i,2*j+1,2*k)=AFine(2*i,2*j+1,2*k)+A2
                      AFine(2*i,2*j-1,2*k)=AFine(2*i,2*j-1,2*k)+A2
                      AFine(2*i,2*j,2*k+1)=AFine(2*i,2*j,2*k+1)+A2
                      AFine(2*i,2*j,2*k-1)=AFine(2*i,2*j,2*k-1)+A2
                      AFine(2*i+1,2*j+1,2*k)=AFine(2*i+1,2*j+1,2*k)+A4
                      AFine(2*i+1,2*j,2*k+1)=AFine(2*i+1,2*j,2*k+1)+A4
                      AFine(2*i,2*j+1,2*k+1)=AFine(2*i,2*j+1,2*k+1)+A4
                      AFine(2*i-1,2*j-1,2*k)=AFine(2*i-1,2*j-1,2*k)+A4
                      AFine(2*i-1,2*j,2*k-1)=AFine(2*i-1,2*j,2*k-1)+A4
                      AFine(2*i,2*j-1,2*k-1)=AFine(2*i,2*j-1,2*k-1)+A4
                      AFine(2*i-1,2*j+1,2*k)=AFine(2*i-1,2*j+1,2*k)+A4
                      AFine(2*i-1,2*j,2*k+1)=AFine(2*i-1,2*j,2*k+1)+A4
                      AFine(2*i+1,2*j-1,2*k)=AFine(2*i+1,2*j-1,2*k)+A4
                      AFine(2*i,2*j-1,2*k+1)=AFine(2*i,2*j-1,2*k+1)+A4
                      AFine(2*i+1,2*j,2*k-1)=AFine(2*i+1,2*j,2*k-1)+A4
                      AFine(2*i,2*j+1,2*k-1)=AFine(2*i,2*j+1,2*k-1)+A4
                      AFine(2*i+1,2*j+1,2*k+1)=AFine(2*i+1,2*j+1,2*k+1)+A8
                      AFine(2*i-1,2*j+1,2*k+1)=AFine(2*i-1,2*j+1,2*k+1)+A8
                      AFine(2*i+1,2*j-1,2*k+1)=AFine(2*i+1,2*j-1,2*k+1)+A8
                      AFine(2*i+1,2*j+1,2*k-1)=AFine(2*i+1,2*j+1,2*k-1)+A8
                      AFine(2*i-1,2*j-1,2*k+1)=AFine(2*i-1,2*j-1,2*k+1)+A8
                      AFine(2*i-1,2*j+1,2*k-1)=AFine(2*i-1,2*j+1,2*k-1)+A8
                      AFine(2*i+1,2*j-1,2*k-1)=AFine(2*i+1,2*j-1,2*k-1)+A8
                      AFine(2*i-1,2*j-1,2*k-1)=AFine(2*i-1,2*j-1,2*k-1)+A8
      enddo
      enddo
    enddo
   endif
   call Bound_Phi_MG(Afine,2*nx,2*ny,2*nz)
 endsubroutine Prolongate


 subroutine Restrict(ACoarse,AFine,level)
 integer,intent(in):: level
 real(KND),dimension(-1:,-1:,-1:),intent(out):: ACoarse
 real(KND),dimension(-1:,-1:,-1:),intent(inout):: AFine
 real(KND) q
 integer:: i,j,k,nx,ny,nz

   nx=bnx*2**level !level means on which grid we restrict
   ny=bny*2**level
   nz=bnz*2**level

   if (GPU>0.and.level>=minGPUlevel-1) then
      !$hmpp <GSKernels> advancedload, args[Re::level]
      !$hmpp <GSKernels> Re callsite,args[*].noupdate=true
          call Re_GPU(nxa,nya,nza,RHSMG(0)%Arr,RHSMG(1)%Arr,RHSMG(2)%Arr,RHSMG(3)%Arr,RHSMG(4)%Arr,&
                             RHSMG(5)%Arr,RHSMG(6)%Arr,RHSMG(7)%Arr,RHSMG(8)%Arr,&
                             ResMG(0)%Arr,ResMG(1)%Arr,ResMG(2)%Arr,ResMG(3)%Arr,ResMG(4)%Arr,&
                             ResMG(5)%Arr,ResMG(6)%Arr,ResMG(7)%Arr,ResMG(8)%Arr,&
                             Btype,level)
   else


      call Bound_Phi_MG(Afine,2*nx,2*ny,2*nz)

      do k=0,nz
       do j=0,ny
        do i=0,nx
            ACoarse(i,j,k)=AFine(2*i,2*j,2*k); q=1.0

            if (i<nx) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i+1,2*j,2*k); q=q+0.5
            endif

            if (i>0) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i-1,2*j,2*k); q=q+0.5
            endif

            if (j<ny) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i,2*j+1,2*k); q=q+0.5
            endif

            if (j>0) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i,2*j-1,2*k); q=q+0.5
            endif

            if (k<nz) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i,2*j,2*k+1); q=q+0.5
            endif

            if (k>0) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i,2*j,2*k-1); q=q+0.5
            endif

            if ((i<nx).and.(j<ny)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i+1,2*j+1,2*k); q=q+0.25
            endif

            if ((i<nx).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i+1,2*j,2*k+1); q=q+0.25
            endif

            if ((j<ny).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i,2*j+1,2*k+1); q=q+0.25
            endif

            if ((i>0).and.(j>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i-1,2*j-1,2*k); q=q+0.25
            endif

            if ((i>0).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i-1,2*j,2*k-1); q=q+0.25
            endif

            if ((j>0).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i,2*j-1,2*k-1); q=q+0.25
            endif

            if ((i>0).and.(j<ny)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i-1,2*j+1,2*k); q=q+0.25
            endif

            if ((i>0).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i-1,2*j,2*k+1); q=q+0.25
            endif

            if ((i<nx).and.(j>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i+1,2*j-1,2*k); q=q+0.25
            endif

            if ((j>0).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i,2*j-1,2*k+1); q=q+0.25
            endif

            if ((i<nx).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i+1,2*j,2*k-1); q=q+0.25
            endif

            if ((j<ny).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i,2*j+1,2*k-1); q=q+0.25
            endif

            if ((i<nx).and.(j<ny).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i+1,2*j+1,2*k+1); q=q+0.125
            endif

            if ((i>0).and.(j<ny).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i-1,2*j+1,2*k+1); q=q+0.125
            endif

            if ((i<nx).and.(j>0).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i+1,2*j-1,2*k+1); q=q+0.125
            endif

            if ((i<nx).and.(j<ny).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i+1,2*j+1,2*k-1); q=q+0.125
            endif

            if ((i>0).and.(j>0).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i-1,2*j-1,2*k+1); q=q+0.125
            endif

            if ((i>0).and.(j<ny).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i-1,2*j+1,2*k-1); q=q+0.125
            endif

            if ((i<nx).and.(j>0).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i+1,2*j-1,2*k-1); q=q+0.125
            endif

            if ((i>0).and.(j>0).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i-1,2*j-1,2*k-1); q=q+0.125
            endif

            ACoarse(i,j,k)=ACoarse(i,j,k)/q
        enddo
        enddo
      enddo
        if (Btype(Ea)==PERIODIC) ACoarse(0,:,:)=ACoarse(nx,:,:)
        if (Btype(No)==PERIODIC) ACoarse(:,0,:)=ACoarse(:,ny,:)
        if (Btype(To)==PERIODIC) ACoarse(:,:,0)=ACoarse(:,:,nz)
   endif
endsubroutine Restrict







pure subroutine BOUND_Phi_MG(Phi,nx,ny,nz)
  real(KND),intent(inout):: Phi(-1:,-1:,-1:)
  integer,intent(in):: nx,ny,nz
  integer i,j,k


   if (Btype(We)==PERIODIC) then
    do k=0,nz
     do j=0,ny                      !Periodic BC
      Phi(-1,j,k)=Phi(nx-1,j,k)
     enddo
    enddo
  else
    do k=0,nz
     do j=0,ny                      !Other BCs
      Phi(-1,j,k)=Phi(0,j,k)
     enddo
    enddo
   endif

   if (Btype(Ea)==PERIODIC) then
    do k=0,nz
     do j=0,ny                      !Periodic BC
      Phi(nx+1,j,k)=Phi(1,j,k)
     enddo
    enddo
   else
    do k=0,nz
     do j=0,ny                      !Other BCs
      Phi(nx+1,j,k)=Phi(nx,j,k)
     enddo
    enddo
   endif

   if (Btype(So)==PERIODIC) then
   do k=0,nz
     do i=-1,nx+1                      !Periodic BC
      Phi(i,-1,k)=Phi(i,ny-1,k)
     enddo
    enddo
   else
    do k=0,nz
     do i=-1,nx+1                      !Other BCs
      Phi(i,-1,k)=Phi(i,0,k)
     enddo
    enddo
   endif

   if (Btype(No)==PERIODIC) then
   do k=0,nz
     do i=-1,nx+1                      !Periodic BC
      Phi(i,ny+1,k)=Phi(i,1,k)
     enddo
    enddo
   else
    do k=0,nz
     do i=-1,nx+1                      !Other BCs
      Phi(i,ny+1,k)=Phi(i,ny,k)
     enddo
    enddo
   endif

   if (Btype(Bo)==PERIODIC) then
    do j=-1,ny+1
     do i=-1,nx+1                      !Periodic BC
      Phi(i,j,-1)=Phi(i,j,nz-1)
     enddo
    enddo
   else
    do j=-1,ny+1
     do i=-1,nx+1                      !Other BCs
      Phi(i,j,-1)=Phi(i,j,0)
     enddo
    enddo
   endif

   if (Btype(To)==PERIODIC) then
    do j=-1,ny+1
     do i=-1,nx+1                      !Periodic BC
      Phi(i,j,nz+1)=Phi(i,j,1)
     enddo
    enddo
   else
    do j=-1,ny+1
     do i=-1,nx+1                      !Other BCs
      Phi(i,j,nz+1)=Phi(i,j,nz)
     enddo
    enddo
   endif
endsubroutine BOUND_Phi_MG



pure function ind(level,nulx,nuly,nulz,i,j,k)
integer ind
integer,intent(in):: level,i,j,k,nulx,nuly,nulz
 ind=(1-nulx)+i+(j-nuly)*(CoefMG(level)%nx+(1-nulx))+(k-nulz)*(CoefMG(level)%nx+1-nulx)*(CoefMG(level)%ny+1-nuly)
endfunction ind




subroutine MG_GE(level)
integer,intent(in)::level
integer i,j,k,l,info
integer,allocatable,dimension(:),save:: ige,ipivot,work2
real(KND),allocatable,dimension(:),save:: xge,bge,work,R,C,ferr,berr
real(KND),allocatable,dimension(:,:),save:: age,af
real(KND) Ap,rcond
integer,save:: nx,ny,nz,nulx,nuly,nulz,nxyz,called=0
character(1),save:: equed

 if (called==0) then
  nx=CoefMG(level)%nx
  ny=CoefMG(level)%ny
  nz=CoefMG(level)%nz

  if (Btype(Ea)==PERIODIC) then
   nulx=1
  else
   nulx=0
  endif
  if (Btype(No)==PERIODIC) then
   nuly=1
  else
   nuly=0
  endif
  if (Btype(To)==PERIODIC) then
   nulz=1
  else
   nulz=0
  endif
  nxyz=(nx+(1-nulx))*(ny+(1-nuly))*(nz+(1-nulz))
  allocate(xge(nxyz))
  allocate(bge(nxyz))
  allocate(ige(nxyz))
  allocate(ipivot(nxyz))
  allocate(age(nxyz,nxyz))
  allocate(af(nxyz,nxyz))
  allocate(work(4*nxyz))
  allocate(work2(nxyz))
  allocate(R(nxyz))
  allocate(C(nxyz))
  allocate(ferr(nxyz))
  allocate(berr(nxyz))

  called=1


  age=0
  bge=0
  xge=0
  ige=0
     do k=nulz,nz
        do j=nuly,ny
           do i=nulx,nx
             l=ind(level,nulx,nuly,nulz,i,j,k)
             Ap=0
             if (i>nulx) then
                       age(l,ind(level,nulx,nuly,nulz,i-1,j,k))=CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             elseif (Btype(We)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nx,j,k))=CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<nx) then
                       age(l,ind(level,nulx,nuly,nulz,i+1,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (Btype(Ea)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nulx,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>nuly) then
                       age(l,ind(level,nulx,nuly,nulz,i,j-1,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (Btype(So)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,ny,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<ny) then
                       age(l,ind(level,nulx,nuly,nulz,i,j+1,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (Btype(No)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,nuly,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>nulz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k-1))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (Btype(Bo)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,nz))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<nz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k+1))=CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (Btype(To)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,nulz))=CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             endif

             age(l,l)=-Ap
             bge(l)=RHSMG(level)%Arr(i,j,k)
            enddo
        enddo
    enddo

    l=ind(level,nulx,nuly,nulz,nx,ny,nz)
    age(1,:)=0
    age(1,1)=age(2,2)
    bge(1)=0


    if (KND==DBL) then
     call DGESVX("E","N",nxyz,1,age,nxyz,af,nxyz,ipivot,EQUED,R,C,bge,nxyz,xge,nxyz,rcond,ferr,berr,work,work2,info)
    else
     call SGESVX("E","N",nxyz,1,age,nxyz,af,nxyz,ipivot,EQUED,R,C,bge,nxyz,xge,nxyz,rcond,ferr,berr,work,work2,info)
    endif


    write (*,*) "RCOND",rcond
    if (info/=0) then
     stop
     write (*,*) "info",info
    endif
  else


    do k=nulz,nz
       do j=nuly,ny
          do i=nulx,nx
            l=ind(level,nulx,nuly,nulz,i,j,k)
            bge(l)=RHSMG(level)%Arr(i,j,k)
          enddo
       enddo
    enddo
    bge(1)=0


    if (KND==DBL) then
     call DGESVX("F","N",nxyz,1,age,nxyz,af,nxyz,ipivot,EQUED,R,C,bge,nxyz,xge,nxyz,rcond,ferr,berr,work,work2,info)
    else
     call SGESVX("F","N",nxyz,1,age,nxyz,af,nxyz,ipivot,EQUED,R,C,bge,nxyz,xge,nxyz,rcond,ferr,berr,work,work2,info)
    endif


    if (info/=0) then
     stop
     write (*,*) "info",info
    endif
  endif
!  do i=1,l
!   write (*,*) "diag",i,l,age(i,i)
!  enddo


!  call DGESV(nxyz,1,age,nxyz,ipivot,bge,nxyz,info)

!  call DGETRI(nxyz,age,nxyz,ipivot,work,nxyz,info)

!  do i=1,l
!   do j=1,l
!    write (*,*) age(i,j)
!   enddo
!  enddo
!  xge=bge
!  call LEGS(age,nxyz,bge,xge,ige)

     do k=nulz,nz
        do j=nuly,ny
           do i=nulx,nx
             l=ind(level,nulx,nuly,nulz,i,j,k)
             PhiMG(level)%Arr(i,j,k)=xge(l)
           enddo
        enddo
     enddo

 if (Btype(Ea)==PERIODIC) then
  PhiMG(level)%Arr(0,:,:)=PhiMG(level)%Arr(nx,:,:)
 endif
 if (Btype(No)==PERIODIC) then
  PhiMG(level)%Arr(:,0,:)=PhiMG(level)%Arr(:,ny,:)
 endif
 if (Btype(To)==PERIODIC) then
  PhiMG(level)%Arr(:,:,0)=PhiMG(level)%Arr(:,:,nz)
 endif

endsubroutine MG_GE


subroutine MG_LU(level)
integer,intent(in)::level
integer i,j,k,l,info
integer,allocatable,dimension(:),save:: ige,ipivot,work2
real(KND),allocatable,dimension(:),save:: xge,bge,work,R,C,ferr,berr
real(KND),allocatable,dimension(:,:),save:: age,af
real(KND) Ap,rcond
integer,save:: nx,ny,nz,nulx,nuly,nulz,nxyz,called=0
character(1),save:: equed

 if (called==0) then
  nx=CoefMG(level)%nx
  ny=CoefMG(level)%ny
  nz=CoefMG(level)%nz

  if (Btype(Ea)==PERIODIC) then
   nulx=1
  else
   nulx=0
  endif
  if (Btype(No)==PERIODIC) then
   nuly=1
  else
   nuly=0
  endif
  if (Btype(To)==PERIODIC) then
   nulz=1
  else
   nulz=0
  endif
  nxyz=(nx+(1-nulx))*(ny+(1-nuly))*(nz+(1-nulz))
  allocate(xge(nxyz))
  allocate(bge(nxyz))
  allocate(ige(nxyz))
  allocate(ipivot(nxyz))
  allocate(age(nxyz,nxyz))
  allocate(af(nxyz,nxyz))
  allocate(work(4*nxyz))
  allocate(work2(nxyz))
  allocate(R(nxyz))
  allocate(C(nxyz))
  allocate(ferr(nxyz))
  allocate(berr(nxyz))

  called=1


  age=0
!   bge=0
  xge=0
  ige=0
     do k=nulz,nz
        do j=nuly,ny
           do i=nulx,nx
             l=ind(level,nulx,nuly,nulz,i,j,k)
             Ap=0
             if (i>nulx) then
                       age(l,ind(level,nulx,nuly,nulz,i-1,j,k))=CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             elseif (Btype(We)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nx,j,k))=CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<nx) then
                       age(l,ind(level,nulx,nuly,nulz,i+1,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (Btype(Ea)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nulx,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>nuly) then
                       age(l,ind(level,nulx,nuly,nulz,i,j-1,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (Btype(So)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,ny,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<ny) then
                       age(l,ind(level,nulx,nuly,nulz,i,j+1,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (Btype(No)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,nuly,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>nulz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k-1))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (Btype(Bo)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,nz))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<nz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k+1))=CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (Btype(To)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,nulz))=CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             endif

             age(l,l)=-Ap
!              bge(l)=RHSMG(level)%Arr(i,j,k)
            enddo
        enddo
    enddo

    l=ind(level,nulx,nuly,nulz,nx,ny,nz)
    age(1,:)=0
    age(1,1)=age(2,2)
!     bge(1)=0


    if (KND==DBL) then
     call  DGETRF(nxyz,nxyz,age,nxyz,ipivot,info)
    else
     call  SGETRF(nxyz,nxyz,age,nxyz,ipivot,info)
    endif

    if (info/=0) then
     stop
     write (*,*) "info",info
    endif
 endif


    do k=nulz,nz
       do j=nuly,ny
          do i=nulx,nx
            l=ind(level,nulx,nuly,nulz,i,j,k)
            bge(l)=RHSMG(level)%Arr(i,j,k)
          enddo
       enddo
    enddo
    bge(1)=0

    if (KND==DBL) then
     call DGETRS("N",nxyz,1,age,nxyz,ipivot,bge,nxyz,info)
    else
     call SGETRS("N",nxyz,1,age,nxyz,ipivot,bge,nxyz,info)
    endif

    if (info/=0) then
     stop
     write (*,*) "info",info
    endif

!  do i=1,l
!   write (*,*) "diag",i,l,age(i,i)
!  enddo


!  call DGESV(nxyz,1,age,nxyz,ipivot,bge,nxyz,info)

!  call DGETRI(nxyz,age,nxyz,ipivot,work,nxyz,info)

!  do i=1,l
!   do j=1,l
!    write (*,*) age(i,j)
!   enddo
!  enddo
!  xge=bge
!  call LEGS(age,nxyz,bge,xge,ige)

     do k=nulz,nz
        do j=nuly,ny
           do i=nulx,nx
             l=ind(level,nulx,nuly,nulz,i,j,k)
             PhiMG(level)%Arr(i,j,k)=bge(l)
           enddo
        enddo
     enddo

 if (Btype(Ea)==PERIODIC) then
  PhiMG(level)%Arr(0,:,:)=PhiMG(level)%Arr(nx,:,:)
 endif
 if (Btype(No)==PERIODIC) then
  PhiMG(level)%Arr(:,0,:)=PhiMG(level)%Arr(:,ny,:)
 endif
 if (Btype(To)==PERIODIC) then
  PhiMG(level)%Arr(:,:,0)=PhiMG(level)%Arr(:,:,nz)
 endif

endsubroutine MG_LU





subroutine MG_INV(level)
integer,intent(in)::level
integer i,j,k,l,info,ldwork
integer,allocatable,dimension(:):: ipivot
real(KND),allocatable,dimension(:),save:: xge,bge,work
real(KND),allocatable,dimension(:,:),save:: age
real(KND) Ap
integer,save:: nx,ny,nz,nulx,nuly,nulz,nxyz,called=0

 if (called==0) then
  nx=CoefMG(level)%nx
  ny=CoefMG(level)%ny
  nz=CoefMG(level)%nz

  if (Btype(Ea)==PERIODIC) then
   nulx=1
  else
   nulx=0
  endif
  if (Btype(No)==PERIODIC) then
   nuly=1
  else
   nuly=0
  endif
  if (Btype(To)==PERIODIC) then
   nulz=1
  else
   nulz=0
  endif
  nxyz=(nx+(1-nulx))*(ny+(1-nuly))*(nz+(1-nulz))
  allocate(xge(nxyz))
  allocate(bge(nxyz))
  allocate(ipivot(nxyz))
  allocate(age(nxyz,nxyz))


  called=1


  age=0
!   bge=0
  xge=0
     do k=nulz,nz
        do j=nuly,ny
           do i=nulx,nx
             l=ind(level,nulx,nuly,nulz,i,j,k)
             Ap=0
             if (i>nulx) then
                       age(l,ind(level,nulx,nuly,nulz,i-1,j,k))=CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             elseif (Btype(We)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nx,j,k))=CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<nx) then
                       age(l,ind(level,nulx,nuly,nulz,i+1,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (Btype(Ea)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nulx,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>nuly) then
                       age(l,ind(level,nulx,nuly,nulz,i,j-1,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (Btype(So)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,ny,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<ny) then
                       age(l,ind(level,nulx,nuly,nulz,i,j+1,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (Btype(No)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,nuly,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>nulz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k-1))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (Btype(Bo)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,nz))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<nz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k+1))=CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (Btype(To)==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,nulz))=CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             endif

             age(l,l)=-Ap
!              bge(l)=RHSMG(level)%Arr(i,j,k)
            enddo
        enddo
    enddo

    l=ind(level,nulx,nuly,nulz,nx,ny,nz)
    age(1,:)=0
    age(1,1)=age(2,2)
!     bge(1)=0


    if (KND==DBL) then
     call  DGETRF(nxyz,nxyz,age,nxyz,ipivot,info)
    else
     call  SGETRF(nxyz,nxyz,age,nxyz,ipivot,info)
    endif

    if (info/=0) then
     stop
     write (*,*) "info",info
    endif

    allocate(work(1))

    if (KND==DBL) then
     call DGETRI(nxyz,age,nxyz,ipivot,work,-1,info)
    else
     call SGETRI(nxyz,age,nxyz,ipivot,work,-1,info)
    endif

    if (info/=0) then
     stop
     write (*,*) "info",info
    endif

    ldwork=nint(work(1))

    deallocate(work)
    allocate(work(ldwork))

    if (KND==DBL) then
     call DGETRI(nxyz,age,nxyz,ipivot,work,ldwork,info)
    else
     call SGETRI(nxyz,age,nxyz,ipivot,work,ldwork,info)
    endif

    if (info/=0) then
     stop
     write (*,*) "info",info
    endif

 endif


    do k=nulz,nz
       do j=nuly,ny
          do i=nulx,nx
            l=ind(level,nulx,nuly,nulz,i,j,k)
            bge(l)=RHSMG(level)%Arr(i,j,k)
          enddo
       enddo
    enddo
    bge(1)=0


!     do j=1,nxyz !rows of X{j}
!      xge(j)=0
!      do i=1,nxyz !columns of A(i,j)
!       xge(j)=xge(i,j)+age(j,i)*b(i)
!      enddo
!     enddo

    xge=0
    do j=1,nxyz !rows of X{j}
     do i=1,nxyz !columns of A(i,j) if (age(i,j)>100*tiny(1._KND))
      xge(i)=xge(i)+age(i,j)*bge(j)
     enddo
    enddo



     do k=nulz,nz
        do j=nuly,ny
           do i=nulx,nx
             l=ind(level,nulx,nuly,nulz,i,j,k)
             PhiMG(level)%Arr(i,j,k)=xge(l)
           enddo
        enddo
     enddo

 if (Btype(Ea)==PERIODIC) then
  PhiMG(level)%Arr(0,:,:)=PhiMG(level)%Arr(nx,:,:)
 endif
 if (Btype(No)==PERIODIC) then
  PhiMG(level)%Arr(:,0,:)=PhiMG(level)%Arr(:,ny,:)
 endif
 if (Btype(To)==PERIODIC) then
  PhiMG(level)%Arr(:,:,0)=PhiMG(level)%Arr(:,:,nz)
 endif

endsubroutine MG_INV








subroutine MG_GS(level,niter)
integer,intent(in)::level,niter
integer i,j,k,l
real(KND) p,Ap

  if (GPU>0.and.level>=minGPUlevel) then!.and.level>3
!    write (*,*) "Gauss-Seidel GPU call level", level
   !$hmpp  <GSKernels> advancedload, args[GS::nit,GS::level]
   !$hmpp  <GSKernels> Gs callsite,args[*].noupdate=true
   call Gs_GPU(nxa,nya,nza,niter,PhiMG(0)%Arr,PhiMG(1)%Arr,PhiMG(2)%Arr,PhiMG(3)%Arr,PhiMG(4)%Arr,&
                             PhiMG(5)%Arr,PhiMG(6)%Arr,PhiMG(7)%Arr,PhiMG(8)%Arr,&
                             RHSMG(0)%Arr,RHSMG(1)%Arr,RHSMG(2)%Arr,RHSMG(3)%Arr,RHSMG(4)%Arr,&
                             RHSMG(5)%Arr,RHSMG(6)%Arr,RHSMG(7)%Arr,RHSMG(8)%Arr,&
                     Aw,Ae,As,An,Ab,At,&
                     Btype,level)

  else

   do l=1,niter
    !$OMP PARALLEL PRIVATE(i,j,k,p,Ap)
    !$OMP DO
    do k=0,CoefMG(level)%nz
        do j=0,CoefMG(level)%ny
            do i=0+mod(j+k,2),CoefMG(level)%nx,2
             p=0
             Ap=0
             if (i>0) then
                       p=p+PhiMG(level)%Arr(i-1,j,k)*CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             elseif (Btype(We)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(CoefMG(level)%nx-1,j,k)*CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<CoefMG(level)%nx) then
                       p=p+PhiMG(level)%Arr(i+1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (Btype(Ea)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>0) then
                       p=p+PhiMG(level)%Arr(i,j-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (Btype(So)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,CoefMG(level)%ny-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<CoefMG(level)%ny) then
                       p=p+PhiMG(level)%Arr(i,j+1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (Btype(No)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>0) then
                       p=p+PhiMG(level)%Arr(i,j,k-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (Btype(Bo)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,j,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<CoefMG(level)%nz) then
                       p=p+PhiMG(level)%Arr(i,j,k+1)*CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (Btype(To)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,j,1)*CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             endif
             p=p-RHSMG(level)%Arr(i,j,k)

             p=p/Ap
             PhiMG(level)%Arr(i,j,k)=p
            enddo
        enddo
    enddo
    !$OMP ENDDO
    !$OMP DO
    do k=0,CoefMG(level)%nz
        do j=0,CoefMG(level)%ny
            do i=0+mod(j+k+1,2),CoefMG(level)%nx,2
             p=0
             Ap=0
             if (i>0) then
                       p=p+PhiMG(level)%Arr(i-1,j,k)*CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             elseif (Btype(We)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(CoefMG(level)%nx-1,j,k)*CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<CoefMG(level)%nx) then
                       p=p+PhiMG(level)%Arr(i+1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (Btype(Ea)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>0) then
                       p=p+PhiMG(level)%Arr(i,j-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (Btype(So)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,CoefMG(level)%ny-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<CoefMG(level)%ny) then
                       p=p+PhiMG(level)%Arr(i,j+1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (Btype(No)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>0) then
                       p=p+PhiMG(level)%Arr(i,j,k-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (Btype(Bo)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,j,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<CoefMG(level)%nz) then
                       p=p+PhiMG(level)%Arr(i,j,k+1)*CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (Btype(To)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,j,1)*CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             endif
             p=p-RHSMG(level)%Arr(i,j,k)
             p=p/Ap

             PhiMG(level)%Arr(i,j,k)=p
            enddo
        enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL

   enddo
  endif
endsubroutine MG_GS

subroutine MG_res(level,R)
integer,intent(in)::level
real(KND),intent(out)::R
integer i,j,k
real(KND),save:: p,Ap

  if (GPU>0.and.level>=minGPUlevel) then!.and.level>3
    !$hmpp <GSKernels> advancedload, args[Res::level]
   !$hmpp  <GSKernels> Res callsite,args[*].noupdate=true
   call Res_GPU(nxa,nya,nza,PhiMG(0)%Arr,PhiMG(1)%Arr,PhiMG(2)%Arr,PhiMG(3)%Arr,PhiMG(4)%Arr,&
                             PhiMG(5)%Arr,PhiMG(6)%Arr,PhiMG(7)%Arr,PhiMG(8)%Arr,&
                             RHSMG(0)%Arr,RHSMG(1)%Arr,RHSMG(2)%Arr,RHSMG(3)%Arr,RHSMG(4)%Arr,&
                             RHSMG(5)%Arr,RHSMG(6)%Arr,RHSMG(7)%Arr,RHSMG(8)%Arr,&
                             ResMG(0)%Arr,ResMG(1)%Arr,ResMG(2)%Arr,ResMG(3)%Arr,ResMG(4)%Arr,&
                             ResMG(5)%Arr,ResMG(6)%Arr,ResMG(7)%Arr,ResMG(8)%Arr,&
                     Aw,Ae,As,An,Ab,At,&
                     Btype,R,level)
  else

     do k=0,CoefMG(level)%nz
        do j=0,CoefMG(level)%ny
            do i=0,CoefMG(level)%nx
             p=0
             Ap=0
             if (i>0) then
                       p=p+PhiMG(level)%Arr(i-1,j,k)*CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             elseif (Btype(We)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(CoefMG(level)%nx-1,j,k)*CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<CoefMG(level)%nx) then
                       p=p+PhiMG(level)%Arr(i+1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (Btype(Ea)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>0) then
                       p=p+PhiMG(level)%Arr(i,j-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (Btype(So)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,CoefMG(level)%ny-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<CoefMG(level)%ny) then
                       p=p+PhiMG(level)%Arr(i,j+1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (Btype(No)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>0) then
                       p=p+PhiMG(level)%Arr(i,j,k-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (Btype(Bo)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,j,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<CoefMG(level)%nz) then
                       p=p+PhiMG(level)%Arr(i,j,k+1)*CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (Btype(To)==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,j,1)*CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             endif
             p=p-RHSMG(level)%Arr(i,j,k)
             p=-p +Ap*PhiMG(level)%Arr(i,j,k)
             ResMG(level)%Arr(i,j,k)=p
            enddo
        enddo
    enddo
    R=MAXVAL(ResMG(level)%Arr(0:CoefMG(level)%nx,0:CoefMG(level)%ny,0:CoefMG(level)%nz))
   endif
endsubroutine MG_res


subroutine MG_clear(level)
integer,intent(in):: level
integer l

 if (GPU>0.and.level>=minGPUlevel) then!.and.level>2
  !$hmpp <GSKernels> advancedload, args[Cl::level]
  !$hmpp <GSKernels> Cl callsite,args[*].noupdate=true
  call Cl_GPU(nxa,nya,nza,PhiMG(0)%Arr,PhiMG(1)%Arr,PhiMG(2)%Arr,PhiMG(3)%Arr,PhiMG(4)%Arr,&
                             PhiMG(5)%Arr,PhiMG(6)%Arr,PhiMG(7)%Arr,&
                             RHSMG(0)%Arr,RHSMG(1)%Arr,RHSMG(2)%Arr,RHSMG(3)%Arr,RHSMG(4)%Arr,&
                             RHSMG(5)%Arr,RHSMG(6)%Arr,RHSMG(7)%Arr,&
                             ResMG(0)%Arr,ResMG(1)%Arr,ResMG(2)%Arr,ResMG(3)%Arr,ResMG(4)%Arr,&
                             ResMG(5)%Arr,ResMG(6)%Arr,ResMG(7)%Arr,level,minmglevel)
 else
  do l=level,minmglevel,-1
   PhiMG(l)%Arr=0
   RHSMG(l)%Arr=0
   ResMG(l)%Arr=0
  enddo
 endif
endsubroutine



recursive subroutine MG_CGC(level,eps,ncgc,npre,npost,R)
integer,intent(in):: level,ncgc,npre,npost
real(KND),intent(in):: eps
real(KND),intent(out):: R
real(KND) R1
integer k

  if (level==minGPUlevel-1) then
    if (level==0) then
    !$hmpp <GSKernels> delegatedstore, args[GS::Phi0,GS::RHS0]
    elseif (level==1) then
    !$hmpp <GSKernels> delegatedstore, args[GS::Phi1,GS::RHS1]
    elseif (level==2) then
    !$hmpp <GSKernels> delegatedstore, args[GS::Phi2,GS::RHS2]
    elseif (level==3) then
    !$hmpp <GSKernels> delegatedstore, args[GS::Phi3,GS::RHS3]
    elseif (level==4) then
    !$hmpp <GSKernels> delegatedstore, args[GS::Phi4,GS::RHS4]
    elseif (level==5) then
    !$hmpp <GSKernels> delegatedstore, args[GS::Phi5,GS::RHS5]
    elseif (level==6) then
    !$hmpp <GSKernels> delegatedstore, args[GS::Phi6,GS::RHS6]
    elseif (level==7) then
    !$hmpp <GSKernels> delegatedstore, args[GS::Phi7,GS::RHS7]
    elseif (level==8) then
    !$hmpp <GSKernels> delegatedstore, args[GS::Phi8,GS::RHS8]
    endif
  endif

  if (level == minmglevel) then !1
     R=huge(mgepsinnerGS)/2.

     if (.not.(GPU>0.and.level>=minGPUlevel)) then
          call MG_GE(level)
          call MG_res(level,R)
!        write (*,*) "GEres",sqrt(R)
     endif

     if (R>mgepsinnerGS**2) then
      R1=R
      do k=1,mgmaxinnerGSiter

        call MG_GS(level, 10)
        call MG_Res(level, R)

         !$hmpp <GSKernels> delegatedStore,Args[Res::R]


 !        write (*,*) "GSres",R
        if (R < mgepsinnerGS) exit
        if (R1<R*1.05_KND) exit
        R1=R
      enddo
     endif
  else
   do k = 1, ncgc !number of recurrent calls

    call MG_GS(level,npre)
    call MG_res(level,R)

    call MG_clear(level-1)

    call Restrict(RHSMG(level-1)%Arr,ResMG(level)%Arr,level-1)


    call MG_CGC(level-1,eps,ncgc,npre,npost,R)

    call Prolongate(PhiMG(level)%Arr,PhiMG(level-1)%Arr,level-1)


    call MG_GS(level,npost)
    call MG_res(level,R)

   enddo
  endif

  if (level==minGPUlevel-1) then
    if (level==0) then
    !$hmpp <GSKernels> advancedload, args[GS::Phi0]
    elseif (level==1) then
    !$hmpp <GSKernels> advancedload, args[GS::Phi1]
    elseif (level==2) then
    !$hmpp <GSKernels> advancedload, args[GS::Phi2]
    elseif (level==3) then
    !$hmpp <GSKernels> advancedload, args[GS::Phi3]
    elseif (level==4) then
    !$hmpp <GSKernels> advancedload, args[GS::Phi4]
    elseif (level==5) then
    !$hmpp <GSKernels> advancedload, args[GS::Phi5]
    elseif (level==6) then
    !$hmpp <GSKernels> advancedload, args[GS::Phi6]
    elseif (level==7) then
    !$hmpp <GSKernels> advancedload, args[GS::Phi7]
    elseif (level==8) then
    !$hmpp <GSKernels> advancedload, args[GS::Phi8]
    endif
  endif
endsubroutine MG_CGC


subroutine POISSMG(Phi,RHS)        !Solves Poisson equation using Successive over-relaxation
real(KND),dimension(0:,0:,0:),intent(inout):: Phi
real(KND),dimension(1:,1:,1:),intent(in)::RHS
integer i,j,k,l,nx,ny,nz,sx,sy,sz
real(KND) mgeps,R,Phiref
real(KND),save:: called=0

 mgeps=epsPoisson
 Phi=0

 if (Btype(Ea)==PERIODIC) then
  sx=1
 else
  sx=0
 endif
 if (Btype(No)==PERIODIC) then
  sy=1
 else
  sy=0
 endif
 if (Btype(To)==PERIODIC) then
  sz=1
 else
  sz=0
 endif


 if (called==0) then
  allocate(PhiMG(0:8),RHSMG(0:8),ResMG(0:8),CoefMG(0:8))
  do l=8,0,-1
   if (l>LMG.or.l<minmglevel) then
    CoefMG(l)%nx=0
    CoefMG(l)%ny=0
    CoefMG(l)%nz=0
    allocate(PhiMG(l)%Arr(1,1,1))
    allocate(RHSMG(l)%Arr(1,1,1))
    allocate(ResMG(l)%Arr(1,1,1))
   else
    CoefMG(l)%nx=bnx*2**l
    CoefMG(l)%ny=bny*2**l
    CoefMG(l)%nz=bnz*2**l
    CoefMG(l)%dx=dxmin*2**(LMG-l)
    CoefMG(l)%dy=dymin*2**(LMG-l)
    CoefMG(l)%dz=dzmin*2**(LMG-l)
    allocate(PhiMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%ny+1,-1:CoefMG(l)%nz+1))
    allocate(RHSMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%ny+1,-1:CoefMG(l)%nz+1))
    allocate(ResMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%ny+1,-1:CoefMG(l)%nz+1))
    CoefMG(l)%Ae=1._KND/(CoefMG(l)%dx*CoefMG(l)%dx)
    CoefMG(l)%Aw=1._KND/(CoefMG(l)%dx*CoefMG(l)%dx)
    CoefMG(l)%An=1._KND/(CoefMG(l)%dy*CoefMG(l)%dy)
    CoefMG(l)%As=1._KND/(CoefMG(l)%dy*CoefMG(l)%dy)
    CoefMG(l)%At=1._KND/(CoefMG(l)%dz*CoefMG(l)%dz)
    CoefMG(l)%Ab=1._KND/(CoefMG(l)%dz*CoefMG(l)%dz)
   endif
  enddo
  called=1
 endif

 nx=bnx*2**LMG
 ny=bny*2**LMG
 nz=bnz*2**LMG

 if (nx-sx/=Prnx-1.or.ny-sy/=Prny-1.or.nz-sz/=Prnz-1) then
    write (*,*) "Incorrect dimensions, multigrid, vs. grid defined in grid.conf:"
    write (*,*) 0+sx,":",nx,"--",1,":",Prnx
    write (*,*) 0+sy,":",ny,"--",1,":",Prny
    write (*,*) 0+sz,":",nz,"--",1,":",Prnz
    stop
 endif


 PhiMG(LMG)%Arr=0
 PhiMG(LMG)%Arr(0+sx:nx,0+sy:ny,0+sz:nz)=Phi(1:Prnx,1:Prny,1:Prnz)
 RHSMG(LMG)%Arr(0+sx:nx,0+sy:ny,0+sz:nz)=RHS(1:Prnx,1:Prny,1:Prnz)

 if (Btype(Ea)==PERIODIC) then
  RHSMG(LMG)%Arr(0,:,:)=RHSMG(LMG)%Arr(Prnx,:,:)
  PhiMG(LMG)%Arr(0,:,:)=PhiMG(LMG)%Arr(Prnx,:,:)
 endif
 if (Btype(No)==PERIODIC) then
  RHSMG(LMG)%Arr(:,0,:)=RHSMG(LMG)%Arr(:,Prny,:)
  PhiMG(LMG)%Arr(:,0,:)=PhiMG(LMG)%Arr(:,Prny,:)
 endif
 if (Btype(To)==PERIODIC) then
  RHSMG(LMG)%Arr(:,:,0)=RHSMG(LMG)%Arr(:,:,Prnz)
  PhiMG(LMG)%Arr(:,:,0)=PhiMG(LMG)%Arr(:,:,Prnz)
 endif





  nxa=(/ CoefMG(0)%nx, CoefMG(1)%nx, CoefMG(2)%nx, CoefMG(3)%nx, CoefMG(4)%nx,&
         CoefMG(5)%nx, CoefMG(6)%nx, CoefMG(7)%nx, CoefMG(8)%nx /)
  nya=(/ CoefMG(0)%ny, CoefMG(1)%ny, CoefMG(2)%ny, CoefMG(3)%ny, CoefMG(4)%ny,&
         CoefMG(5)%ny, CoefMG(6)%ny, CoefMG(7)%ny, CoefMG(8)%ny /)
  nza=(/ CoefMG(0)%nz, CoefMG(1)%nz, CoefMG(2)%nz, CoefMG(3)%nz, CoefMG(4)%nz,&
         CoefMG(5)%nz, CoefMG(6)%nz, CoefMG(7)%nz, CoefMG(8)%nz /)

  Aw=(/ CoefMG(0)%Aw, CoefMG(1)%Aw, CoefMG(2)%Aw, CoefMG(3)%Aw, CoefMG(4)%Aw,&
         CoefMG(5)%Aw, CoefMG(6)%Aw, CoefMG(7)%Aw, CoefMG(8)%Aw /)
  Ae=(/ CoefMG(0)%Ae, CoefMG(1)%Ae, CoefMG(2)%Ae, CoefMG(3)%Ae, CoefMG(4)%Ae,&
         CoefMG(5)%Ae, CoefMG(6)%Ae, CoefMG(7)%Ae, CoefMG(8)%Ae /)
  As=(/ CoefMG(0)%As, CoefMG(1)%As, CoefMG(2)%As, CoefMG(3)%As, CoefMG(4)%As,&
         CoefMG(5)%As, CoefMG(6)%As, CoefMG(7)%As, CoefMG(8)%As /)
  An=(/ CoefMG(0)%An, CoefMG(1)%An, CoefMG(2)%An, CoefMG(3)%An, CoefMG(4)%An,&
         CoefMG(5)%An, CoefMG(6)%An, CoefMG(7)%An, CoefMG(8)%An /)
  Ab=(/ CoefMG(0)%Ab, CoefMG(1)%Ab, CoefMG(2)%Ab, CoefMG(3)%Ab, CoefMG(4)%Ab,&
         CoefMG(5)%Ab, CoefMG(6)%Ab, CoefMG(7)%Ab, CoefMG(8)%Ab /)
  At=(/ CoefMG(0)%At, CoefMG(1)%At, CoefMG(2)%At, CoefMG(3)%At, CoefMG(4)%At,&
         CoefMG(5)%At, CoefMG(6)%At, CoefMG(7)%At, CoefMG(8)%At /)

  !$hmpp <GSKernels> allocate
  !$hmpp <GSkernels> advancedload, args[GS::nx,GS::ny,GS::nz,GS::Aw,GS::Ae,GS::As,GS::An,GS::Ab,GS::At,&
  !$hmpp <GSkernels>  GS::Btype]




  if (LMG==0) then
   !$hmpp <GSKernels> advancedload, args[GS::Phi0,GS::RHS0]
  elseif (LMG==1) then
   !$hmpp <GSKernels> advancedload, args[GS::Phi1,GS::RHS1]
  elseif (LMG==2) then
   !$hmpp <GSKernels> advancedload, args[GS::Phi2,GS::RHS2]
  elseif (LMG==3) then
   !$hmpp <GSKernels> advancedload, args[GS::Phi3,GS::RHS3]
  elseif (LMG==4) then
   !$hmpp <GSKernels> advancedload, args[GS::Phi4,GS::RHS4]
  elseif (LMG==5) then
   !$hmpp <GSKernels> advancedload, args[GS::Phi5,GS::RHS5]
  elseif (LMG==6) then
   !$hmpp <GSKernels> advancedload, args[GS::Phi6,GS::RHS6]
  elseif (LMG==7) then
   !$hmpp <GSKernels> advancedload, args[GS::Phi7,GS::RHS7]
  elseif (LMG==8) then
   !$hmpp <GSKernels> advancedload, args[GS::Phi8,GS::RHS8]
  endif

 do l=1,maxPoissoniter
  write (*,*) "MG iteration:",l
!    Phiref=SUM(PhiMG(LMG)%Arr(0:CoefMG(LMG)%nx,0:CoefMG(LMG)%ny,0:CoefMG(LMG)%nz))/&
!           ((CoefMG(LMG)%nx+1)*(CoefMG(LMG)%ny+1)*(CoefMG(LMG)%nz+1))
!    do k=0,CoefMG(LMG)%nz
!     do j=0,CoefMG(LMG)%ny
!      do i=0,CoefMG(LMG)%nx
!       PhiMG(LMG)%Arr(i,j,k)=PhiMG(LMG)%Arr(i,j,k)-Phiref
!      enddo
!     enddo
!    enddo
!    call Bound_Phi_MG(PhiMG(LMG)%Arr,CoefMG(LMG)%nx,CoefMG(LMG)%ny,CoefMG(LMG)%nz)




  call MG_CGC(LMG,mgeps,mgncgc,mgnpre,mgnpost,R)


  if (LMG==0) then
   !$hmpp <GSKernels> delegatedstore, args[Res::R]
  elseif (LMG==1) then
   !$hmpp <GSKernels> delegatedstore, args[Res::R]
  elseif (LMG==2) then
   !$hmpp <GSKernels> delegatedstore, args[Res::R]
  elseif (LMG==3) then
   !$hmpp <GSKernels> delegatedstore, args[Res::R]
  elseif (LMG==4) then
   !$hmpp <GSKernels> delegatedstore, args[Res::R]
  elseif (LMG==5) then
   !$hmpp <GSKernels> delegatedstore, args[Res::R]
  elseif (LMG==6) then
   !$hmpp <GSKernels> delegatedstore, args[Res::R]
  elseif (LMG==7) then
   !$hmpp <GSKernels> delegatedstore, args[Res::R]
  elseif (LMG==8) then
   !$hmpp <GSKernels> delegatedstore, args[Res::R]
  endif


  write (*,*) "MG residuum",R
  if (R<mgeps)   exit
 enddo

  if (LMG==0) then
   !$hmpp <GSKernels> delegatedstore, args[GS::Phi0,Res::R]
  elseif (LMG==1) then
   !$hmpp <GSKernels> delegatedstore, args[GS::Phi1,Res::R]
  elseif (LMG==2) then
   !$hmpp <GSKernels> delegatedstore, args[GS::Phi2,Res::R]
  elseif (LMG==3) then
   !$hmpp <GSKernels> delegatedstore, args[GS::Phi3,Res::R]
  elseif (LMG==4) then
   !$hmpp <GSKernels> delegatedstore, args[GS::Phi4,Res::R]
  elseif (LMG==5) then
   !$hmpp <GSKernels> delegatedstore, args[GS::Phi5,Res::R]
  elseif (LMG==6) then
   !$hmpp <GSKernels> delegatedstore, args[GS::Phi6,Res::R]
  elseif (LMG==7) then
   !$hmpp <GSKernels> delegatedstore, args[GS::Phi7,Res::R]
  elseif (LMG==8) then
   !$hmpp <GSKernels> delegatedstore, args[GS::Phi8,Res::R]
  endif

  !$hmpp <GSKernels> release

 Phi(1:Prnx,1:Prny,1:Prnz)=PhiMG(LMG)%Arr(0+sx:nx,0+sy:ny,0+sz:nz)

endsubroutine POISSMG
















!GPU Codelets.

  !$hmpp <GSKernels> group, target=CUDA


  !$hmpp <GSKernels> mapbyname, nx,ny,nz,Btype
  !$hmpp <GSKernels> mapbyname, Phi0,Phi1,Phi2,Phi3,Phi4,Phi5,Phi6,Phi7,Phi8
  !$hmpp <GSKernels> mapbyname, RHS0,RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,RHS7,RHS8
  !$hmpp <GSKernels> mapbyname, Res0,Res1,Res2,Res3,Res4,Res5,Res6,Res7,Res8
  !$hmpp <GSKernels> mapbyname, Aw,Ae,As,An,Ab,At


 !$hmpp <GSKernels> Pr codelet
 subroutine Pr_GPU(nx,ny,nz,Phi0,Phi1,Phi2,Phi3,Phi4,Phi5,Phi6,Phi7,Phi8,&
                                 Btype,level)
 implicit none
 integer,parameter:: KND=4,PERIODIC=3
  integer,intent(in),dimension(0:8) :: nx,ny,nz
  real(KND),dimension(-1:nx(0)+1,-1:ny(0)+1,-1:nz(0)+1),intent(inout)::Phi0
  real(KND),dimension(-1:nx(1)+1,-1:ny(1)+1,-1:nz(1)+1),intent(inout)::Phi1
  real(KND),dimension(-1:nx(2)+1,-1:ny(2)+1,-1:nz(2)+1),intent(inout)::Phi2
  real(KND),dimension(-1:nx(3)+1,-1:ny(3)+1,-1:nz(3)+1),intent(inout)::Phi3
  real(KND),dimension(-1:nx(4)+1,-1:ny(4)+1,-1:nz(4)+1),intent(inout)::Phi4
  real(KND),dimension(-1:nx(5)+1,-1:ny(5)+1,-1:nz(5)+1),intent(inout)::Phi5
  real(KND),dimension(-1:nx(6)+1,-1:ny(6)+1,-1:nz(6)+1),intent(inout)::Phi6
  real(KND),dimension(-1:nx(7)+1,-1:ny(7)+1,-1:nz(7)+1),intent(inout)::Phi7
  real(KND),dimension(-1:nx(8)+1,-1:ny(8)+1,-1:nz(8)+1),intent(inout)::Phi8
  integer,intent(in) :: Btype(6),level
  integer :: nxl, nyl, nzl

  nxl=nx(level)
  nyl=ny(level)
  nzl=nz(level)

    if (level==0) then
     call Prolongate_GPU(nxl,nyl,nzl,Phi1,Phi0,Btype)
    elseif (level==1) then
     call Prolongate_GPU(nxl,nyl,nzl,Phi2,Phi1,Btype)
    elseif (level==2) then
     call Prolongate_GPU(nxl,nyl,nzl,Phi3,Phi2,Btype)
    elseif (level==3) then
     call Prolongate_GPU(nxl,nyl,nzl,Phi4,Phi3,Btype)
    elseif (level==4) then
     call Prolongate_GPU(nxl,nyl,nzl,Phi5,Phi4,Btype)
    elseif (level==5) then
     call Prolongate_GPU(nxl,nyl,nzl,Phi6,Phi5,Btype)
    elseif (level==6) then
     call Prolongate_GPU(nxl,nyl,nzl,Phi7,Phi6,Btype)
    elseif (level==7) then
     call Prolongate_GPU(nxl,nyl,nzl,Phi8,Phi7,Btype)
    endif
 endsubroutine Pr_GPU


 subroutine Prolongate_GPU(nx,ny,nz,AFine,ACoarse,Btype)

 implicit none
 integer,parameter:: KND=4,PERIODIC=3
 integer,intent(in):: nx,ny,nz,Btype(6)
 real(KND),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(in):: ACoarse
 real(KND),dimension(-1:2*nx+1,-1:2*ny+1,-1:2*nz+1),intent(inout):: AFine
 real(KND),dimension(-1:1,-1:1,-1:1):: Cf
 integer:: i,j,k,ii,jj,kk
 intrinsic abs

 do kk=-1,1
  do jj=-1,1
   do ii=-1,1
    if (ii==0.and.jj==0.and.kk==0) then
        Cf(ii,jj,kk)=1._KND
    else if (abs(ii)+abs(jj)+abs(kk)==1) then
        Cf(ii,jj,kk)=1/2._KND
    else if (abs(ii)+abs(jj)+abs(kk)==2) then
        Cf(ii,jj,kk)=1/4._KND
    else
        Cf(ii,jj,kk)=1/8._KND
    end if
   enddo
  enddo
 enddo

 do kk=-1,1
  do jj=-1,1
   do ii=-1,1
   !$hmppcg grid blocksize 512x1
   !$hmppcg permute (k,i,j)
   !$hmppcg gridify(k,i)
    do k=0,nz
     do j=0,ny
      do i=0,nx
                      AFine(2*i+ii,2*j+jj,2*k+kk)=AFine(2*i+ii,2*j+jj,2*k+kk)+Cf(ii,jj,kk)*ACoarse(i,j,k)
      enddo
     enddo
    enddo

   enddo
  enddo
 enddo

   if (Btype(We)==PERIODIC) then
    do k=0,2*nz
     do j=0,2*ny                      !Periodic BC
      Afine(-1,j,k)=Afine(2*nx-1,j,k)
     enddo
    enddo
  else
    do k=0,2*nz
     do j=0,2*ny                      !Other BCs
      Afine(-1,j,k)=Afine(0,j,k)
     enddo
    enddo
   endif

   if (Btype(Ea)==PERIODIC) then
    do k=0,2*nz
     do j=0,2*ny                      !Periodic BC
      Afine(2*nx+1,j,k)=Afine(1,j,k)
     enddo
    enddo
   else
    do k=0,2*nz
     do j=0,2*ny                      !Other BCs
      Afine(2*nx+1,j,k)=Afine(2*nx,j,k)
     enddo
    enddo
   endif

   if (Btype(So)==PERIODIC) then
   do k=0,2*nz
     do i=-1,2*nx+1                      !Periodic BC
      Afine(i,-1,k)=Afine(i,2*ny-1,k)
     enddo
    enddo
   else
    do k=0,2*nz
     do i=-1,2*nx+1                      !Other BCs
      Afine(i,-1,k)=Afine(i,0,k)
     enddo
    enddo
   endif

   if (Btype(No)==PERIODIC) then
   do k=0,2*nz
     do i=-1,2*nx+1                      !Periodic BC
      Afine(i,2*ny+1,k)=Afine(i,1,k)
     enddo
    enddo
   else
    do k=0,2*nz
     do i=-1,2*nx+1                      !Other BCs
      Afine(i,2*ny+1,k)=Afine(i,2*ny,k)
     enddo
    enddo
   endif

   if (Btype(Bo)==PERIODIC) then
    do j=-1,2*ny+1
     do i=-1,2*nx+1                      !Periodic BC
      Afine(i,j,-1)=Afine(i,j,2*nz-1)
     enddo
    enddo
   else
    do j=-1,2*ny+1
     do i=-1,2*nx+1                      !Other BCs
      Afine(i,j,-1)=Afine(i,j,0)
     enddo
    enddo
   endif

   if (Btype(To)==PERIODIC) then
    do j=-1,2*ny+1
     do i=-1,2*nx+1                      !Periodic BC
      Afine(i,j,2*nz+1)=Afine(i,j,1)
     enddo
    enddo
   else
    do j=-1,2*ny+1
     do i=-1,2*nx+1                      !Other BCs
      Afine(i,j,2*nz+1)=Afine(i,j,2*nz)
     enddo
    enddo
   endif
 endsubroutine Prolongate_GPU





 !$hmpp <GSKernels> Re codelet
 subroutine Re_GPU(nx,ny,nz,RHS0,RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,RHS7,RHS8,&
                            Res0,Res1,Res2,Res3,Res4,Res5,Res6,Res7,Res8,&
                                 Btype,level)
 implicit none
 integer,parameter:: KND=4,PERIODIC=3
  integer,intent(in),dimension(0:8) :: nx,ny,nz
  real(KND),dimension(-1:nx(0)+1,-1:ny(0)+1,-1:nz(0)+1),intent(inout)::Res0,RHS0
  real(KND),dimension(-1:nx(1)+1,-1:ny(1)+1,-1:nz(1)+1),intent(inout)::Res1,RHS1
  real(KND),dimension(-1:nx(2)+1,-1:ny(2)+1,-1:nz(2)+1),intent(inout)::Res2,RHS2
  real(KND),dimension(-1:nx(3)+1,-1:ny(3)+1,-1:nz(3)+1),intent(inout)::Res3,RHS3
  real(KND),dimension(-1:nx(4)+1,-1:ny(4)+1,-1:nz(4)+1),intent(inout)::Res4,RHS4
  real(KND),dimension(-1:nx(5)+1,-1:ny(5)+1,-1:nz(5)+1),intent(inout)::Res5,RHS5
  real(KND),dimension(-1:nx(6)+1,-1:ny(6)+1,-1:nz(6)+1),intent(inout)::Res6,RHS6
  real(KND),dimension(-1:nx(7)+1,-1:ny(7)+1,-1:nz(7)+1),intent(inout)::Res7,RHS7
  real(KND),dimension(-1:nx(8)+1,-1:ny(8)+1,-1:nz(8)+1),intent(inout)::Res8,RHS8
  integer,intent(in) :: Btype(6),level
  integer :: nxl, nyl, nzl

  nxl=nx(level)
  nyl=ny(level)
  nzl=nz(level)

    if (level==0) then
     call Restrict_GPU(nxl,nyl,nzl,RHS0,Res1,Btype)
    elseif (level==1) then
     call Restrict_GPU(nxl,nyl,nzl,RHS1,Res2,Btype)
    elseif (level==2) then
     call Restrict_GPU(nxl,nyl,nzl,RHS2,Res3,Btype)
    elseif (level==3) then
     call Restrict_GPU(nxl,nyl,nzl,RHS3,Res4,Btype)
    elseif (level==4) then
     call Restrict_GPU(nxl,nyl,nzl,RHS4,Res5,Btype)
    elseif (level==5) then
     call Restrict_GPU(nxl,nyl,nzl,RHS5,Res6,Btype)
    elseif (level==6) then
     call Restrict_GPU(nxl,nyl,nzl,RHS6,Res7,Btype)
    elseif (level==7) then
     call Restrict_GPU(nxl,nyl,nzl,RHS6,Res7,Btype)
    endif
 endsubroutine Re_GPU




 subroutine Restrict_GPU(nx,ny,nz,ACoarse,AFine,Btype)
 implicit none
 integer,parameter:: KND=4,PERIODIC=3
 integer,intent(in):: nx,ny,nz,Btype(6)
 real(KND),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(out):: ACoarse
 real(KND),dimension(-1:2*nx+1,-1:2*ny+1,-1:2*nz+1),intent(inout):: AFine
 real(KND) q
 integer:: i,j,k


   if (Btype(We)==PERIODIC) then
    do k=0,2*nz
     do j=0,2*ny                      !Periodic BC
      Afine(-1,j,k)=Afine(2*nx-1,j,k)
     enddo
    enddo
  else
    do k=0,2*nz
     do j=0,2*ny                      !Other BCs
      Afine(-1,j,k)=Afine(0,j,k)
     enddo
    enddo
   endif

   if (Btype(Ea)==PERIODIC) then
    do k=0,2*nz
     do j=0,2*ny                      !Periodic BC
      Afine(2*nx+1,j,k)=Afine(1,j,k)
     enddo
    enddo
   else
    do k=0,2*nz
     do j=0,2*ny                      !Other BCs
      Afine(2*nx+1,j,k)=Afine(2*nx,j,k)
     enddo
    enddo
   endif

   if (Btype(So)==PERIODIC) then
   do k=0,2*nz
     do i=-1,2*nx+1                      !Periodic BC
      Afine(i,-1,k)=Afine(i,2*ny-1,k)
     enddo
    enddo
   else
    do k=0,2*nz
     do i=-1,2*nx+1                      !Other BCs
      Afine(i,-1,k)=Afine(i,0,k)
     enddo
    enddo
   endif

   if (Btype(No)==PERIODIC) then
   do k=0,2*nz
     do i=-1,2*nx+1                      !Periodic BC
      Afine(i,2*ny+1,k)=Afine(i,1,k)
     enddo
    enddo
   else
    do k=0,2*nz
     do i=-1,2*nx+1                      !Other BCs
      Afine(i,2*ny+1,k)=Afine(i,2*ny,k)
     enddo
    enddo
   endif

   if (Btype(Bo)==PERIODIC) then
    do j=-1,2*ny+1
     do i=-1,2*nx+1                      !Periodic BC
      Afine(i,j,-1)=Afine(i,j,2*nz-1)
     enddo
    enddo
   else
    do j=-1,2*ny+1
     do i=-1,2*nx+1                      !Other BCs
      Afine(i,j,-1)=Afine(i,j,0)
     enddo
    enddo
   endif

   if (Btype(To)==PERIODIC) then
    do j=-1,2*ny+1
     do i=-1,2*nx+1                      !Periodic BC
      Afine(i,j,2*nz+1)=Afine(i,j,1)
     enddo
    enddo
   else
    do j=-1,2*ny+1
     do i=-1,2*nx+1                      !Other BCs
      Afine(i,j,2*nz+1)=Afine(i,j,2*nz)
     enddo
    enddo
   endif


   !$hmppcg permute (k,i,j)
   !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
      do k=0,nz
       do j=0,ny
        do i=0,nx
            ACoarse(i,j,k)=AFine(2*i,2*j,2*k); q=1.0

            if (i<nx) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i+1,2*j,2*k); q=q+0.5
            endif

            if (i>0) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i-1,2*j,2*k); q=q+0.5
            endif

            if (j<ny) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i,2*j+1,2*k); q=q+0.5
            endif

            if (j>0) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i,2*j-1,2*k); q=q+0.5
            endif

            if (k<nz) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i,2*j,2*k+1); q=q+0.5
            endif

            if (k>0) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.5*AFine(2*i,2*j,2*k-1); q=q+0.5
            endif

            if ((i<nx).and.(j<ny)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i+1,2*j+1,2*k); q=q+0.25
            endif

            if ((i<nx).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i+1,2*j,2*k+1); q=q+0.25
            endif

            if ((j<ny).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i,2*j+1,2*k+1); q=q+0.25
            endif

            if ((i>0).and.(j>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i-1,2*j-1,2*k); q=q+0.25
            endif

            if ((i>0).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i-1,2*j,2*k-1); q=q+0.25
            endif

            if ((j>0).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i,2*j-1,2*k-1); q=q+0.25
            endif

            if ((i>0).and.(j<ny)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i-1,2*j+1,2*k); q=q+0.25
            endif

            if ((i>0).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i-1,2*j,2*k+1); q=q+0.25
            endif

            if ((i<nx).and.(j>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i+1,2*j-1,2*k); q=q+0.25
            endif

            if ((j>0).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i,2*j-1,2*k+1); q=q+0.25
            endif

            if ((i<nx).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i+1,2*j,2*k-1); q=q+0.25
            endif

            if ((j<ny).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.25*AFine(2*i,2*j+1,2*k-1); q=q+0.25
            endif

            if ((i<nx).and.(j<ny).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i+1,2*j+1,2*k+1); q=q+0.125
            endif

            if ((i>0).and.(j<ny).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i-1,2*j+1,2*k+1); q=q+0.125
            endif

            if ((i<nx).and.(j>0).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i+1,2*j-1,2*k+1); q=q+0.125
            endif

            if ((i<nx).and.(j<ny).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i+1,2*j+1,2*k-1); q=q+0.125
            endif

            if ((i>0).and.(j>0).and.(k<nz)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i-1,2*j-1,2*k+1); q=q+0.125
            endif

            if ((i>0).and.(j<ny).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i-1,2*j+1,2*k-1); q=q+0.125
            endif

            if ((i<nx).and.(j>0).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i+1,2*j-1,2*k-1); q=q+0.125
            endif

            if ((i>0).and.(j>0).and.(k>0)) then
              ACoarse(i,j,k)=ACoarse(i,j,k)+0.125*AFine(2*i-1,2*j-1,2*k-1); q=q+0.125
            endif

            ACoarse(i,j,k)=ACoarse(i,j,k)/q
        enddo
       enddo
      enddo

      do k=-1,nz+1
       do j=-1,ny+1
        if (Btype(Ea)==PERIODIC) ACoarse(0,j,k)=ACoarse(nx,j,k)
       enddo
      enddo
      do k=-1,nz+1
       do i=-1,nx+1
        if (Btype(No)==PERIODIC) ACoarse(i,0,k)=ACoarse(i,ny,k)
       enddo
      enddo
      do j=-1,ny+1
       do i=-1,nx+1
        if (Btype(To)==PERIODIC) ACoarse(i,j,0)=ACoarse(i,j,nz)
       enddo
      enddo
endsubroutine Restrict_GPU





  !$hmpp <GSKernels> GS codelet
  subroutine GS_GPU(nx,ny,nz,nit,Phi0,Phi1,Phi2,Phi3,Phi4,Phi5,Phi6,Phi7,Phi8,&
                              RHS0,RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,RHS7,RHS8,&
                     Aw,Ae,As,An,Ab,At,&
                     Btype,level)
  implicit none
  integer,parameter:: KND=4
  integer,intent(in),dimension(0:8) :: nx,ny,nz
  integer,intent(in) :: nit
  real(KND),dimension(-1:nx(0)+1,-1:ny(0)+1,-1:nz(0)+1),intent(inout)::Phi0,RHS0
  real(KND),dimension(-1:nx(1)+1,-1:ny(1)+1,-1:nz(1)+1),intent(inout)::Phi1,RHS1
  real(KND),dimension(-1:nx(2)+1,-1:ny(2)+1,-1:nz(2)+1),intent(inout)::Phi2,RHS2
  real(KND),dimension(-1:nx(3)+1,-1:ny(3)+1,-1:nz(3)+1),intent(inout)::Phi3,RHS3
  real(KND),dimension(-1:nx(4)+1,-1:ny(4)+1,-1:nz(4)+1),intent(inout)::Phi4,RHS4
  real(KND),dimension(-1:nx(5)+1,-1:ny(5)+1,-1:nz(5)+1),intent(inout)::Phi5,RHS5
  real(KND),dimension(-1:nx(6)+1,-1:ny(6)+1,-1:nz(6)+1),intent(inout)::Phi6,RHS6
  real(KND),dimension(-1:nx(7)+1,-1:ny(7)+1,-1:nz(7)+1),intent(inout)::Phi7,RHS7
  real(KND),dimension(-1:nx(8)+1,-1:ny(8)+1,-1:nz(8)+1),intent(inout)::Phi8,RHS8
  integer,intent(in) :: level
  real(KND),dimension(0:8),intent(in) :: Aw,Ae,As,An,Ab,At
  integer,intent(in)   :: Btype(6)
  integer :: nxl, nyl, nzl

  nxl=nx(level)
  nyl=ny(level)
  nzl=nz(level)

   if (level==0) then
      call MG_Gs_GPU(nxl,nyl,nzl,nit,&
                      Phi0,RHS0,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype)
   elseif (level==1) then
      call MG_Gs_GPU(nxl,nyl,nzl,nit,&
                      Phi1,RHS1,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype)
   elseif (level==2) then
      call MG_Gs_GPU(nxl,nyl,nzl,nit,&
                      Phi2,RHS2,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype)
   elseif (level==3) then
      call MG_Gs_GPU(nxl,nyl,nzl,nit,&
                      Phi3,RHS3,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype)
   elseif (level==4) then
      call MG_Gs_GPU(nxl,nyl,nzl,nit,&
                      Phi4,RHS4,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype)
   elseif (level==5) then
      call MG_Gs_GPU(nxl,nyl,nzl,nit,&
                      Phi5,RHS5,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype)
   elseif (level==6) then
      call MG_Gs_GPU(nxl,nyl,nzl,nit,&
                      Phi6,RHS6,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype)
   elseif (level==7) then
      call MG_Gs_GPU(nxl,nyl,nzl,nit,&
                      Phi7,RHS7,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype)
   elseif (level==8) then
      call MG_Gs_GPU(nxl,nyl,nzl,nit,&
                      Phi8,RHS8,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype)
   endif

  endsubroutine GS_GPU



  subroutine MG_GS_GPU(nx,ny,nz,nit,Phi,RHS,&
                     Aw,Ae,As,An,Ab,At,&
                     Btype)
  implicit none

  integer,parameter:: KND=4,PERIODIC=3

  integer,intent(in)   :: nx,ny,nz,nit,Btype(6)
  real(KND),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(inout)::Phi
  real(KND),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(inout)::RHS
  real(KND),intent(in) :: Aw,Ae,As,An,Ab,At
  integer i,j,k,l
  real(KND) :: p,Ap
  intrinsic mod

  do l=1,nit
   !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
     do k=0,nz
        do i=0,nx
            do j=0+mod(i+k,2),ny,2
             p=0
             Ap=0
             if (i>0) then
                       p=p+Phi(i-1,j,k)*Aw
                       Ap=Ap+Aw
             elseif (Btype(We)==PERIODIC) then
                       p=p+Phi(nx-1,j,k)*Aw
                       Ap=Ap+Aw
             endif
             if (i<nx) then
                       p=p+Phi(i+1,j,k)*Ae
                       Ap=Ap+Ae
             elseif (Btype(Ea)==PERIODIC) then
                       p=p+Phi(1,j,k)*Ae
                       Ap=Ap+Ae
             endif
             if (j>0) then
                       p=p+Phi(i,j-1,k)*As
                       Ap=Ap+As
             elseif (Btype(So)==PERIODIC) then
                       p=p+Phi(i,ny-1,k)*As
                       Ap=Ap+As
             endif
             if (j<ny) then
                       p=p+Phi(i,j+1,k)*An
                       Ap=Ap+An
             elseif (Btype(No)==PERIODIC) then
                       p=p+Phi(i,1,k)*An
                       Ap=Ap+An
             endif
             if (k>0) then
                       p=p+Phi(i,j,k-1)*Ab
                       Ap=Ap+Ab
             elseif (Btype(Bo)==PERIODIC) then
                       p=p+Phi(i,j,nz-1)*Ab
                       Ap=Ap+Ab
             endif
             if (k<nz) then
                       p=p+Phi(i,j,k+1)*At
                       Ap=Ap+At
             elseif (Btype(To)==PERIODIC) then
                       p=p+Phi(i,j,1)*At
                       Ap=Ap+At
             endif
             p=p-RHS(i,j,k)

             p=p/Ap
             Phi(i,j,k)=p
            enddo
        enddo
    enddo
  !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
     do k=0,nz
        do i=0,nx
            do j=0+mod(i+k+1,2),ny,2
             p=0
             Ap=0
             if (i>0) then
                       p=p+Phi(i-1,j,k)*Aw
                       Ap=Ap+Aw
             elseif (Btype(We)==PERIODIC) then
                       p=p+Phi(nx-1,j,k)*Aw
                       Ap=Ap+Aw
             endif
             if (i<nx) then
                       p=p+Phi(i+1,j,k)*Ae
                       Ap=Ap+Ae
             elseif (Btype(Ea)==PERIODIC) then
                       p=p+Phi(1,j,k)*Ae
                       Ap=Ap+Ae
             endif
             if (j>0) then
                       p=p+Phi(i,j-1,k)*As
                       Ap=Ap+As
             elseif (Btype(So)==PERIODIC) then
                       p=p+Phi(i,ny-1,k)*As
                       Ap=Ap+As
             endif
             if (j<ny) then
                       p=p+Phi(i,j+1,k)*An
                       Ap=Ap+An
             elseif (Btype(No)==PERIODIC) then
                       p=p+Phi(i,1,k)*An
                       Ap=Ap+An
             endif
             if (k>0) then
                       p=p+Phi(i,j,k-1)*Ab
                       Ap=Ap+Ab
             elseif (Btype(Bo)==PERIODIC) then
                       p=p+Phi(i,j,nz-1)*Ab
                       Ap=Ap+Ab
             endif
             if (k<nz) then
                       p=p+Phi(i,j,k+1)*At
                       Ap=Ap+At
             elseif (Btype(To)==PERIODIC) then
                       p=p+Phi(i,j,1)*At
                       Ap=Ap+At
             endif
             p=p-RHS(i,j,k)

             p=p/Ap
             Phi(i,j,k)=p
            enddo
        enddo
    enddo
  enddo
  endsubroutine MG_GS_GPU


!$hmpp <GSKernels> Res codelet
  subroutine Res_GPU(nx,ny,nz,Phi0,Phi1,Phi2,Phi3,Phi4,Phi5,Phi6,Phi7,Phi8,&
                              RHS0,RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,RHS7,RHS8,&
                              Res0,Res1,Res2,Res3,Res4,Res5,Res6,Res7,Res8,&
                     Aw,Ae,As,An,Ab,At,&
                     Btype,R,level)
  implicit none
  integer,parameter:: KND=4
  integer,intent(in),dimension(0:8) :: nx,ny,nz
  real(KND),dimension(-1:nx(0)+1,-1:ny(0)+1,-1:nz(0)+1),intent(inout)::Phi0,Res0,RHS0
  real(KND),dimension(-1:nx(1)+1,-1:ny(1)+1,-1:nz(1)+1),intent(inout)::Phi1,Res1,RHS1
  real(KND),dimension(-1:nx(2)+1,-1:ny(2)+1,-1:nz(2)+1),intent(inout)::Phi2,Res2,RHS2
  real(KND),dimension(-1:nx(3)+1,-1:ny(3)+1,-1:nz(3)+1),intent(inout)::Phi3,Res3,RHS3
  real(KND),dimension(-1:nx(4)+1,-1:ny(4)+1,-1:nz(4)+1),intent(inout)::Phi4,Res4,RHS4
  real(KND),dimension(-1:nx(5)+1,-1:ny(5)+1,-1:nz(5)+1),intent(inout)::Phi5,Res5,RHS5
  real(KND),dimension(-1:nx(6)+1,-1:ny(6)+1,-1:nz(6)+1),intent(inout)::Phi6,Res6,RHS6
  real(KND),dimension(-1:nx(7)+1,-1:ny(7)+1,-1:nz(7)+1),intent(inout)::Phi7,Res7,RHS7
  real(KND),dimension(-1:nx(8)+1,-1:ny(8)+1,-1:nz(8)+1),intent(inout)::Phi8,Res8,RHS8
  integer,intent(in) :: level
  real(KND),dimension(0:8),intent(in) :: Aw,Ae,As,An,Ab,At
  integer,intent(in)   :: Btype(6)
  real(KND),intent(out) :: R
  integer :: nxl, nyl, nzl

  nxl=nx(level)
  nyl=ny(level)
  nzl=nz(level)


   if (level==0) then
      call MG_Res_GPU(nxl,nyl,nzl,&
                      Phi0,Res0,RHS0,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype,R)
   elseif (level==1) then
      call MG_Res_GPU(nxl,nyl,nzl,&
                      Phi1,Res1,RHS1,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype,R)
   elseif (level==2) then
      call MG_Res_GPU(nxl,nyl,nzl,&
                      Phi2,Res2,RHS2,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype,R)
   elseif (level==3) then
      call MG_Res_GPU(nxl,nyl,nzl,&
                      Phi3,Res3,RHS3,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype,R)
   elseif (level==4) then
      call MG_Res_GPU(nxl,nyl,nzl,&
                      Phi4,Res4,RHS4,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype,R)
   elseif (level==5) then
      call MG_Res_GPU(nxl,nyl,nzl,&
                      Phi5,Res5,RHS5,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype,R)
   elseif (level==6) then
      call MG_Res_GPU(nxl,nyl,nzl,&
                      Phi6,Res6,RHS6,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype,R)
   elseif (level==7) then
      call MG_Res_GPU(nxl,nyl,nzl,&
                      Phi7,Res7,RHS7,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype,R)
   elseif (level==8) then
      call MG_Res_GPU(nxl,nyl,nzl,&
                      Phi8,Res8,RHS8,&
                      Aw(level),Ae(level),As(level),An(level),Ab(level),At(level),&
                      Btype,R)
   endif

  end subroutine Res_GPU





  subroutine MG_Res_GPU(nx,ny,nz,Phi,Res,RHS,&
                             Aw,Ae,As,An,Ab,At,&
                             Btype,R)

   implicit none


   integer,parameter:: KND=4,PERIODIC=3

   integer,intent(in)   :: nx,ny,nz,Btype(6)
   real(KND),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(in)::Phi
   real(KND),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(out)::Res
   real(KND),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(in)::RHS
   real(KND),intent(in) :: Aw,Ae,As,An,Ab,At
   real(KND),intent(out) :: R
   integer i,j,k,l
   real(KND) :: p,Ap
   intrinsic max, abs

    !$hmppcg grid blocksize 512x1
    !$hmppcg gridify(k,i)
     do k=0,nz
        do i=0,nx
            do j=0,ny
             p=0
             Ap=0
             if (i>0) then
                       p=p+Phi(i-1,j,k)*Aw
                       Ap=Ap+Aw
             elseif (Btype(We)==PERIODIC) then
                       p=p+Phi(nx-1,j,k)*Aw
                       Ap=Ap+Aw
             endif
             if (i<nx) then
                       p=p+Phi(i+1,j,k)*Ae
                       Ap=Ap+Ae
             elseif (Btype(Ea)==PERIODIC) then
                       p=p+Phi(1,j,k)*Ae
                       Ap=Ap+Ae
             endif
             if (j>0) then
                       p=p+Phi(i,j-1,k)*As
                       Ap=Ap+As
             elseif (Btype(So)==PERIODIC) then
                       p=p+Phi(i,ny-1,k)*As
                       Ap=Ap+As
             endif
             if (j<ny) then
                       p=p+Phi(i,j+1,k)*An
                       Ap=Ap+An
             elseif (Btype(No)==PERIODIC) then
                       p=p+Phi(i,1,k)*An
                       Ap=Ap+An
             endif
             if (k>0) then
                       p=p+Phi(i,j,k-1)*Ab
                       Ap=Ap+Ab
             elseif (Btype(Bo)==PERIODIC) then
                       p=p+Phi(i,j,nz-1)*Ab
                       Ap=Ap+Ab
             endif
             if (k<nz) then
                       p=p+Phi(i,j,k+1)*At
                       Ap=Ap+At
             elseif (Btype(To)==PERIODIC) then
                       p=p+Phi(i,j,1)*At
                       Ap=Ap+At
             endif
             p=p-RHS(i,j,k)
             p=-p +Ap*Phi(i,j,k)
             Res(i,j,k)=p
            enddo
        enddo
     enddo

     R=0
     !$hmppcg grid blocksize 512x1
     !$hmppcg gridify(k,i), reduce(max:R)
     do k=0,nz
        do i=0,nx
            do j=0,ny
             R=max(R,abs(Res(i,j,k)))
            enddo
        enddo
    enddo
  endsubroutine MG_Res_GPU






  !$hmpp <GSKernels> Cl codelet
  subroutine Cl_GPU(nx,ny,nz,Phi0,Phi1,Phi2,Phi3,Phi4,Phi5,Phi6,Phi7,&
                             RHS0,RHS1,RHS2,RHS3,RHS4,RHS5,RHS6,RHS7,&
                             Res0,Res1,Res2,Res3,Res4,Res5,Res6,Res7,level,minmglevel)
  implicit none
  integer,parameter:: KND=4
  integer,intent(in),dimension(0:8) :: nx,ny,nz
  real(KND),dimension(-1:nx(0)+1,-1:ny(0)+1,-1:nz(0)+1),intent(inout)::Phi0,Res0,RHS0
  real(KND),dimension(-1:nx(1)+1,-1:ny(1)+1,-1:nz(1)+1),intent(inout)::Phi1,Res1,RHS1
  real(KND),dimension(-1:nx(2)+1,-1:ny(2)+1,-1:nz(2)+1),intent(inout)::Phi2,Res2,RHS2
  real(KND),dimension(-1:nx(3)+1,-1:ny(3)+1,-1:nz(3)+1),intent(inout)::Phi3,Res3,RHS3
  real(KND),dimension(-1:nx(4)+1,-1:ny(4)+1,-1:nz(4)+1),intent(inout)::Phi4,Res4,RHS4
  real(KND),dimension(-1:nx(5)+1,-1:ny(5)+1,-1:nz(5)+1),intent(inout)::Phi5,Res5,RHS5
  real(KND),dimension(-1:nx(6)+1,-1:ny(6)+1,-1:nz(6)+1),intent(inout)::Phi6,Res6,RHS6
  real(KND),dimension(-1:nx(7)+1,-1:ny(7)+1,-1:nz(7)+1),intent(inout)::Phi7,Res7,RHS7
  integer,intent(in) :: level,minmglevel
  integer :: l
  integer :: nxl, nyl, nzl

  nxl=nx(level)
  nyl=ny(level)
  nzl=nz(level)

  do l=level,minmglevel,-1
   if (l==0) then
    call Clear_GPU(nx(l),ny(l),nz(l),Phi0,RHS0,Res0)
   elseif (l==1) then
    call Clear_GPU(nx(l),ny(l),nz(l),Phi1,RHS1,Res1)
   elseif (l==2) then
    call Clear_GPU(nx(l),ny(l),nz(l),Phi2,RHS2,Res2)
   elseif (l==3) then
    call Clear_GPU(nx(l),ny(l),nz(l),Phi3,RHS3,Res3)
   elseif (l==4) then
    call Clear_GPU(nx(l),ny(l),nz(l),Phi4,RHS4,Res4)
   elseif (l==5) then
    call Clear_GPU(nx(l),ny(l),nz(l),Phi5,RHS5,Res5)
   elseif (l==6) then
    call Clear_GPU(nx(l),ny(l),nz(l),Phi6,RHS6,Res6)
   elseif (l==7) then
    call Clear_GPU(nx(l),ny(l),nz(l),Phi7,RHS7,Res7)
   endif
  enddo
  end subroutine Cl_GPU








  subroutine Clear_GPU(nx,ny,nz,Phi,RHS,Res)
  implicit none
  integer,parameter:: KND=4
  integer,intent(in) :: nx,ny,nz
  real(KND),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(out)::Phi
  real(KND),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(out)::Res
  real(KND),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(out)::RHS
  integer :: i,j,k
   !$hmppcg permute (k,i,j)
   !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
  do k=-1,nz+1
   do j=-1,ny+1
    do i=-1,nx+1
     Phi(i,j,k)=0
    enddo
   enddo
  enddo

   !$hmppcg permute (k,i,j)
   !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
  do k=-1,nz+1
   do j=-1,ny+1
    do i=-1,nx+1
     Res(i,j,k)=0
    enddo
   enddo
  enddo

   !$hmppcg permute (k,i,j)
   !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
  do k=-1,nz+1
   do j=-1,ny+1
    do i=-1,nx+1
     RHS(i,j,k)=0
    enddo
   enddo
  enddo

  endsubroutine Clear_GPU


endmodule MULTIGRID

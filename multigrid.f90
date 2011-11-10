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
  integer mgncgc !type of cycling 1..V cycle, 2..W cycle
  integer mgnpre
  integer mgnpost
  integer mgmaxinnerGSiter
  real(KND) mgepsinnerGS


contains 


 subroutine SetMGParams(llmg, lminmglevel, lbnx, lbny, lbnz,&
                          lmgncgc, lmgnpre, lmgnpost, lmgmaxinnerGSiter, lmgepsinnerGS)

  integer,intent(in)::   llmg, lminmglevel, lbnx, lbny, lbnz, lmgncgc, lmgnpre, lmgnpost, lmgmaxinnerGSiter
  real(KND),intent(in):: lmgepsinnerGS

  lmg=llmg
  minmglevel=lminmglevel
  bnx=lbnx
  bny=lbny
  bnz=lbnz
  mgncgc=lmgncgc
  mgnpre=lmgnpre
  mgnpost=lmgnpost
  mgmaxinnerGSiter=lmgmaxinnerGSiter
  mgepsinnerGS=lmgepsinnerGS
 endsubroutine SetMGParams

 pure subroutine Prolongate(AFine,ACoarse,level)
  implicit none
  integer,intent(in):: level
  real(KND),dimension(-1:,-1:,-1:),intent(in):: ACoarse
  real(KND),dimension(-1:,-1:,-1:),intent(inout):: AFine
  integer:: i,j,k,nx,ny,nz
 
   nx=bnx*2**level !level means from which grid we interpolate
   ny=bny*2**level
   nz=bnz*2**level

   do k=0,nz
    do j=0,ny
     do i=0,nx
                     AFine(2*i,2*j,2*k)=AFine(2*i,2*j,2*k)+ACoarse(i,j,k)
                     AFine(2*i+1,2*j,2*k)=AFine(2*i+1,2*j,2*k)+0.5*ACoarse(i,j,k)
                     AFine(2*i-1,2*j,2*k)=AFine(2*i-1,2*j,2*k)+0.5*ACoarse(i,j,k)
                     AFine(2*i,2*j+1,2*k)=AFine(2*i,2*j+1,2*k)+0.5*ACoarse(i,j,k)
                     AFine(2*i,2*j-1,2*k)=AFine(2*i,2*j-1,2*k)+0.5*ACoarse(i,j,k)
                     AFine(2*i,2*j,2*k+1)=AFine(2*i,2*j,2*k+1)+0.5*ACoarse(i,j,k)
                     AFine(2*i,2*j,2*k-1)=AFine(2*i,2*j,2*k-1)+0.5*ACoarse(i,j,k)
                     AFine(2*i+1,2*j+1,2*k)=AFine(2*i+1,2*j+1,2*k)+0.25*ACoarse(i,j,k)
                     AFine(2*i+1,2*j,2*k+1)=AFine(2*i+1,2*j,2*k+1)+0.25*ACoarse(i,j,k)
                     AFine(2*i,2*j+1,2*k+1)=AFine(2*i,2*j+1,2*k+1)+0.25*ACoarse(i,j,k)
                     AFine(2*i-1,2*j-1,2*k)=AFine(2*i-1,2*j-1,2*k)+0.25*ACoarse(i,j,k)
                     AFine(2*i-1,2*j,2*k-1)=AFine(2*i-1,2*j,2*k-1)+0.25*ACoarse(i,j,k)
                     AFine(2*i,2*j-1,2*k-1)=AFine(2*i,2*j-1,2*k-1)+0.25*ACoarse(i,j,k)
                     AFine(2*i-1,2*j+1,2*k)=AFine(2*i-1,2*j+1,2*k)+0.25*ACoarse(i,j,k)
                     AFine(2*i-1,2*j,2*k+1)=AFine(2*i-1,2*j,2*k+1)+0.25*ACoarse(i,j,k)
                     AFine(2*i+1,2*j-1,2*k)=AFine(2*i+1,2*j-1,2*k)+0.25*ACoarse(i,j,k)
                     AFine(2*i,2*j-1,2*k+1)=AFine(2*i,2*j-1,2*k+1)+0.25*ACoarse(i,j,k)
                     AFine(2*i+1,2*j,2*k-1)=AFine(2*i+1,2*j,2*k-1)+0.25*ACoarse(i,j,k)
                     AFine(2*i,2*j+1,2*k-1)=AFine(2*i,2*j+1,2*k-1)+0.25*ACoarse(i,j,k)
                     AFine(2*i+1,2*j+1,2*k+1)=AFine(2*i+1,2*j+1,2*k+1)+0.125*ACoarse(i,j,k)
                     AFine(2*i-1,2*j+1,2*k+1)=AFine(2*i-1,2*j+1,2*k+1)+0.125*ACoarse(i,j,k)
                     AFine(2*i+1,2*j-1,2*k+1)=AFine(2*i+1,2*j-1,2*k+1)+0.125*ACoarse(i,j,k)
                     AFine(2*i+1,2*j+1,2*k-1)=AFine(2*i+1,2*j+1,2*k-1)+0.125*ACoarse(i,j,k)
                     AFine(2*i-1,2*j-1,2*k+1)=AFine(2*i-1,2*j-1,2*k+1)+0.125*ACoarse(i,j,k)
                     AFine(2*i-1,2*j+1,2*k-1)=AFine(2*i-1,2*j+1,2*k-1)+0.125*ACoarse(i,j,k)
                     AFine(2*i+1,2*j-1,2*k-1)=AFine(2*i+1,2*j-1,2*k-1)+0.125*ACoarse(i,j,k)
                     AFine(2*i-1,2*j-1,2*k-1)=AFine(2*i-1,2*j-1,2*k-1)+0.125*ACoarse(i,j,k)
     enddo
    enddo
   enddo

   call Bound_Phi_MG(Afine,2*nx,2*ny,2*nz)
 endsubroutine Prolongate


 subroutine Restrict(ACoarse,AFine,level)
 implicit none
 integer,intent(in):: level
 real(KND),dimension(-1:,-1:,-1:),intent(out):: ACoarse
 real(KND),dimension(-1:,-1:,-1:),intent(inout):: AFine
 real(KND) q
 integer:: i,j,k,nx,ny,nz
 
   nx=bnx*2**level !level means on which grid we restrict
   ny=bny*2**level
   nz=bnz*2**level
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
    if (BtypeE==PERIODIC) ACoarse(0,:,:)=ACoarse(nx,:,:)
    if (BtypeN==PERIODIC) ACoarse(:,0,:)=ACoarse(:,ny,:)
    if (BtypeT==PERIODIC) ACoarse(:,:,0)=ACoarse(:,:,nz)
endsubroutine Restrict







pure subroutine BOUND_Phi_MG(Phi,nx,ny,nz)
  real(KND),intent(inout):: Phi(-1:,-1:,-1:)
  integer,intent(in):: nx,ny,nz
  integer i,j,k


   if (BtypeW==PERIODIC) then
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

   if (BtypeE==PERIODIC) then
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

   if (BtypeS==PERIODIC) then
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

   if (BtypeN==PERIODIC) then
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

   if (BtypeB==PERIODIC) then
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

   if (BtypeT==PERIODIC) then
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

  if (BtypeE==PERIODIC) then
   nulx=1
  else
   nulx=0
  endif
  if (BtypeN==PERIODIC) then
   nuly=1
  else
   nuly=0
  endif
  if (BtypeT==PERIODIC) then
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
             elseif (BtypeW==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nx,j,k))=CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<nx) then
                       age(l,ind(level,nulx,nuly,nulz,i+1,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (BtypeE==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nulx,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>nuly) then
                       age(l,ind(level,nulx,nuly,nulz,i,j-1,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (BtypeS==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,ny,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<ny) then
                       age(l,ind(level,nulx,nuly,nulz,i,j+1,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (BtypeN==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,nuly,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>nulz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k-1))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (BtypeB==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,nz))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<nz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k+1))=CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (BtypeT==PERIODIC) then
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

 if (BtypeE==PERIODIC) then
  PhiMG(level)%Arr(0,:,:)=PhiMG(level)%Arr(nx,:,:)
 endif
 if (BtypeN==PERIODIC) then
  PhiMG(level)%Arr(:,0,:)=PhiMG(level)%Arr(:,ny,:)
 endif
 if (BtypeT==PERIODIC) then
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

  if (BtypeE==PERIODIC) then
   nulx=1
  else
   nulx=0
  endif
  if (BtypeN==PERIODIC) then
   nuly=1
  else
   nuly=0
  endif
  if (BtypeT==PERIODIC) then
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
             elseif (BtypeW==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nx,j,k))=CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<nx) then
                       age(l,ind(level,nulx,nuly,nulz,i+1,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (BtypeE==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nulx,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>nuly) then
                       age(l,ind(level,nulx,nuly,nulz,i,j-1,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (BtypeS==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,ny,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<ny) then
                       age(l,ind(level,nulx,nuly,nulz,i,j+1,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (BtypeN==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,nuly,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>nulz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k-1))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (BtypeB==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,nz))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<nz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k+1))=CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (BtypeT==PERIODIC) then
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

 if (BtypeE==PERIODIC) then
  PhiMG(level)%Arr(0,:,:)=PhiMG(level)%Arr(nx,:,:)
 endif
 if (BtypeN==PERIODIC) then
  PhiMG(level)%Arr(:,0,:)=PhiMG(level)%Arr(:,ny,:)
 endif
 if (BtypeT==PERIODIC) then
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

  if (BtypeE==PERIODIC) then
   nulx=1
  else
   nulx=0
  endif
  if (BtypeN==PERIODIC) then
   nuly=1
  else
   nuly=0
  endif
  if (BtypeT==PERIODIC) then
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
             elseif (BtypeW==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nx,j,k))=CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<nx) then
                       age(l,ind(level,nulx,nuly,nulz,i+1,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (BtypeE==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,nulx,j,k))=CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>nuly) then
                       age(l,ind(level,nulx,nuly,nulz,i,j-1,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (BtypeS==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,ny,k))=CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<ny) then
                       age(l,ind(level,nulx,nuly,nulz,i,j+1,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (BtypeN==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,nuly,k))=CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>nulz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k-1))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (BtypeB==PERIODIC) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,nz))=CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<nz) then
                       age(l,ind(level,nulx,nuly,nulz,i,j,k+1))=CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (BtypeT==PERIODIC) then
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

 if (BtypeE==PERIODIC) then
  PhiMG(level)%Arr(0,:,:)=PhiMG(level)%Arr(nx,:,:)
 endif
 if (BtypeN==PERIODIC) then
  PhiMG(level)%Arr(:,0,:)=PhiMG(level)%Arr(:,ny,:)
 endif
 if (BtypeT==PERIODIC) then
  PhiMG(level)%Arr(:,:,0)=PhiMG(level)%Arr(:,:,nz)
 endif

endsubroutine MG_INV








subroutine MG_GS(level,niter)
integer,intent(in)::level,niter
integer i,j,k,l
real(KND) p,Ap

  if (GPU>0.and.level>4) then
   write (*,*) "Gauss-Seidel GPU call level", level

   if (level==0) then
      !$hmpp <GsKernels> advancedload, args[GS0::Phi,GS0::RHS]
      !$hmpp  <GSKernels> GS0 callsite
      call MG_GS_GPU(CoefMG(0)%nx,CoefMG(0)%ny,CoefMG(0)%nz,niter,&
                    PhiMG(0)%Arr,RHSMG(0)%Arr,&
                    CoefMG(0)%Aw,CoefMG(0)%Ae,CoefMG(0)%As,CoefMG(0)%An,CoefMG(0)%Ab,CoefMG(0)%At,&
                    BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT)
      !$hmpp <GSKernels> delegatedStore,Args[GS0::Phi]

   elseif (level==1) then
      !$hmpp <GsKernels> advancedload, args[GS1::Phi,GS1::RHS]
      !$hmpp  <GSKernels> GS1 callsite
      call MG_GS_GPU(CoefMG(1)%nx,CoefMG(1)%ny,CoefMG(1)%nz,niter,&
                    PhiMG(1)%Arr,RHSMG(1)%Arr,&
                    CoefMG(1)%Aw,CoefMG(1)%Ae,CoefMG(1)%As,CoefMG(1)%An,CoefMG(1)%Ab,CoefMG(1)%At,&
                    BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT)
      !$hmpp <GSKernels> delegatedStore,Args[GS1::Phi]

   elseif (level==2) then
      !$hmpp <GsKernels> advancedload, args[GS2::Phi,GS2::RHS]
      !$hmpp  <GSKernels> GS2 callsite
      call MG_GS_GPU(CoefMG(2)%nx,CoefMG(2)%ny,CoefMG(2)%nz,niter,&
                    PhiMG(2)%Arr,RHSMG(2)%Arr,&
                    CoefMG(2)%Aw,CoefMG(2)%Ae,CoefMG(2)%As,CoefMG(2)%An,CoefMG(2)%Ab,CoefMG(2)%At,&
                    BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT)
      !$hmpp <GSKernels> delegatedStore,Args[GS5::Phi]

   elseif (level==3) then
      !$hmpp <GsKernels> advancedload, args[GS3::Phi,GS3::RHS]
      !$hmpp  <GSKernels> GS3 callsite
      call MG_GS_GPU(CoefMG(3)%nx,CoefMG(3)%ny,CoefMG(3)%nz,niter,&
                    PhiMG(3)%Arr,RHSMG(3)%Arr,&
                    CoefMG(3)%Aw,CoefMG(3)%Ae,CoefMG(3)%As,CoefMG(3)%An,CoefMG(3)%Ab,CoefMG(3)%At,&
                    BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT)
      !$hmpp <GSKernels> delegatedStore,Args[GS3::Phi]

   elseif (level==4) then
      !$hmpp <GsKernels> advancedload, args[GS4::Phi,GS4::RHS]
      !$hmpp  <GSKernels> GS4 callsite
      call MG_GS_GPU(CoefMG(4)%nx,CoefMG(4)%ny,CoefMG(4)%nz,niter,&
                    PhiMG(4)%Arr,RHSMG(4)%Arr,&
                    CoefMG(4)%Aw,CoefMG(4)%Ae,CoefMG(4)%As,CoefMG(4)%An,CoefMG(4)%Ab,CoefMG(4)%At,&
                    BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT)
      !$hmpp <GSKernels> delegatedStore,Args[GS4::Phi]

   elseif (level==5) then
      !$hmpp <GsKernels> advancedload, args[GS5::Phi,GS5::RHS]
      !$hmpp  <GSKernels> GS5 callsite
      call MG_GS_GPU(CoefMG(5)%nx,CoefMG(5)%ny,CoefMG(5)%nz,niter,&
                    PhiMG(5)%Arr,RHSMG(5)%Arr,&
                    CoefMG(5)%Aw,CoefMG(5)%Ae,CoefMG(5)%As,CoefMG(5)%An,CoefMG(5)%Ab,CoefMG(5)%At,&
                    BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT)
      !$hmpp <GSKernels> delegatedStore,Args[GS5::Phi]

   elseif (level==6) then
      !$hmpp <GsKernels> advancedload, args[GS6::Phi,GS6::RHS]
      !$hmpp  <GSKernels> GS6 callsite
      call MG_GS_GPU(CoefMG(6)%nx,CoefMG(6)%ny,CoefMG(6)%nz,niter,&
                    PhiMG(6)%Arr,RHSMG(6)%Arr,&
                    CoefMG(6)%Aw,CoefMG(6)%Ae,CoefMG(6)%As,CoefMG(6)%An,CoefMG(6)%Ab,CoefMG(6)%At,&
                    BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT)
      !$hmpp <GSKernels> delegatedStore,Args[GS6::Phi]


   elseif (level==7) then
      !$hmpp <GsKernels> advancedload, args[GS7::Phi,GS7::RHS]
      !$hmpp  <GSKernels> GS7 callsite
      call MG_GS_GPU(CoefMG(7)%nx,CoefMG(7)%ny,CoefMG(7)%nz,niter,&
                    PhiMG(7)%Arr,RHSMG(7)%Arr,&
                    CoefMG(7)%Aw,CoefMG(7)%Ae,CoefMG(7)%As,CoefMG(7)%An,CoefMG(7)%Ab,CoefMG(7)%At,&
                    BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT)
      !$hmpp <GSKernels> delegatedStore,Args[GS7::Phi]

   elseif (level==8) then
      !$hmpp <GsKernels> advancedload, args[GS8::Phi,GS8::RHS]
      !$hmpp  <GSKernels> GS8 callsite
      call MG_GS_GPU(CoefMG(8)%nx,CoefMG(8)%ny,CoefMG(8)%nz,niter,&
                    PhiMG(8)%Arr,RHSMG(8)%Arr,&
                    CoefMG(8)%Aw,CoefMG(8)%Ae,CoefMG(8)%As,CoefMG(8)%An,CoefMG(8)%Ab,CoefMG(8)%At,&
                    BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT)
      !$hmpp <GSKernels> delegatedStore,Args[GS8::Phi]

   endif


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
             elseif (BtypeW==PERIODIC) then
                       p=p+PhiMG(level)%Arr(CoefMG(level)%nx-1,j,k)*CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<CoefMG(level)%nx) then
                       p=p+PhiMG(level)%Arr(i+1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (BtypeE==PERIODIC) then
                       p=p+PhiMG(level)%Arr(1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>0) then
                       p=p+PhiMG(level)%Arr(i,j-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (BtypeS==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,CoefMG(level)%ny-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<CoefMG(level)%ny) then
                       p=p+PhiMG(level)%Arr(i,j+1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (BtypeN==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>0) then
                       p=p+PhiMG(level)%Arr(i,j,k-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (BtypeB==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,j,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<CoefMG(level)%nz) then
                       p=p+PhiMG(level)%Arr(i,j,k+1)*CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (BtypeT==PERIODIC) then
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
             elseif (BtypeW==PERIODIC) then
                       p=p+PhiMG(level)%Arr(CoefMG(level)%nx-1,j,k)*CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<CoefMG(level)%nx) then
                       p=p+PhiMG(level)%Arr(i+1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (BtypeE==PERIODIC) then
                       p=p+PhiMG(level)%Arr(1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>0) then
                       p=p+PhiMG(level)%Arr(i,j-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (BtypeS==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,CoefMG(level)%ny-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<CoefMG(level)%ny) then
                       p=p+PhiMG(level)%Arr(i,j+1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (BtypeN==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>0) then
                       p=p+PhiMG(level)%Arr(i,j,k-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (BtypeB==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,j,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<CoefMG(level)%nz) then
                       p=p+PhiMG(level)%Arr(i,j,k+1)*CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (BtypeT==PERIODIC) then
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

  if (GPU>0.and.level>4) then
   if (level==0) then
      !$hmpp  <GSKernels> Res0 callsite
      call MG_Res_GPU(CoefMG(0)%nx,CoefMG(0)%ny,CoefMG(0)%nz,&
                      PhiMG(0)%Arr,ResMG(0)%Arr,RHSMG(0)%Arr,&
                      CoefMG(0)%Aw,CoefMG(0)%Ae,CoefMG(0)%As,CoefMG(0)%An,CoefMG(0)%Ab,CoefMG(0)%At,&
                      BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT,R)
      !$hmpp <GSKernels> delegatedStore,Args[Res0::Res]
   elseif (level==1) then
      !$hmpp  <GSKernels> Res1 callsite
      call MG_Res_GPU(CoefMG(1)%nx,CoefMG(1)%ny,CoefMG(1)%nz,&
                      PhiMG(1)%Arr,ResMG(1)%Arr,RHSMG(1)%Arr,&
                      CoefMG(1)%Aw,CoefMG(1)%Ae,CoefMG(1)%As,CoefMG(1)%An,CoefMG(1)%Ab,CoefMG(1)%At,&
                      BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT,R)
      !$hmpp <GSKernels> delegatedStore,Args[Res1::Res]
   elseif (level==2) then
      !$hmpp  <GSKernels> Res2 callsite
      call MG_Res_GPU(CoefMG(2)%nx,CoefMG(2)%ny,CoefMG(2)%nz,&
                      PhiMG(2)%Arr,ResMG(2)%Arr,RHSMG(2)%Arr,&
                      CoefMG(2)%Aw,CoefMG(2)%Ae,CoefMG(2)%As,CoefMG(2)%An,CoefMG(2)%Ab,CoefMG(2)%At,&
                      BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT,R)
      !$hmpp <GSKernels> delegatedStore,Args[Res2::Res]
   elseif (level==3) then
      !$hmpp  <GSKernels> Res3 callsite
      call MG_Res_GPU(CoefMG(3)%nx,CoefMG(3)%ny,CoefMG(3)%nz,&
                      PhiMG(3)%Arr,ResMG(3)%Arr,RHSMG(3)%Arr,&
                      CoefMG(3)%Aw,CoefMG(3)%Ae,CoefMG(3)%As,CoefMG(3)%An,CoefMG(3)%Ab,CoefMG(3)%At,&
                      BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT,R)
      !$hmpp <GSKernels> delegatedStore,Args[Res3::Res]
   elseif (level==4) then
      !$hmpp  <GSKernels> Res4 callsite
      call MG_Res_GPU(CoefMG(4)%nx,CoefMG(4)%ny,CoefMG(4)%nz,&
                      PhiMG(4)%Arr,ResMG(4)%Arr,RHSMG(4)%Arr,&
                      CoefMG(4)%Aw,CoefMG(4)%Ae,CoefMG(4)%As,CoefMG(4)%An,CoefMG(4)%Ab,CoefMG(4)%At,&
                      BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT,R)
      !$hmpp <GSKernels> delegatedStore,Args[Res4::Res]
   elseif (level==5) then
      !$hmpp  <GSKernels> Res5 callsite
      call MG_Res_GPU(CoefMG(5)%nx,CoefMG(5)%ny,CoefMG(5)%nz,&
                      PhiMG(5)%Arr,ResMG(5)%Arr,RHSMG(5)%Arr,&
                      CoefMG(5)%Aw,CoefMG(5)%Ae,CoefMG(5)%As,CoefMG(5)%An,CoefMG(5)%Ab,CoefMG(5)%At,&
                      BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT,R)
      !$hmpp <GSKernels> delegatedStore,Args[Res5::Res]
   elseif (level==6) then
      !$hmpp  <GSKernels> Res6 callsite
      call MG_Res_GPU(CoefMG(6)%nx,CoefMG(6)%ny,CoefMG(6)%nz,&
                      PhiMG(6)%Arr,ResMG(6)%Arr,RHSMG(6)%Arr,&
                      CoefMG(6)%Aw,CoefMG(6)%Ae,CoefMG(6)%As,CoefMG(6)%An,CoefMG(6)%Ab,CoefMG(6)%At,&
                      BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT,R)
      !$hmpp <GSKernels> delegatedStore,Args[Res6::Res]
   elseif (level==7) then
      !$hmpp  <GSKernels> Res7 callsite
      call MG_Res_GPU(CoefMG(7)%nx,CoefMG(7)%ny,CoefMG(7)%nz,&
                      PhiMG(7)%Arr,ResMG(7)%Arr,RHSMG(7)%Arr,&
                      CoefMG(7)%Aw,CoefMG(7)%Ae,CoefMG(7)%As,CoefMG(7)%An,CoefMG(7)%Ab,CoefMG(7)%At,&
                      BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT,R)
      !$hmpp <GSKernels> delegatedStore,Args[Res7::Res]
   elseif (level==8) then
      !$hmpp  <GSKernels> Res8 callsite
      call MG_Res_GPU(CoefMG(8)%nx,CoefMG(8)%ny,CoefMG(8)%nz,&
                      PhiMG(8)%Arr,ResMG(8)%Arr,RHSMG(8)%Arr,&
                      CoefMG(8)%Aw,CoefMG(8)%Ae,CoefMG(8)%As,CoefMG(8)%An,CoefMG(8)%Ab,CoefMG(8)%At,&
                      BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT,R)
      !$hmpp <GSKernels> delegatedStore,Args[Res8::Res]
   endif
  else

     do k=0,CoefMG(level)%nz
        do j=0,CoefMG(level)%ny
            do i=0,CoefMG(level)%nx
             p=0
             Ap=0
             if (i>0) then
                       p=p+PhiMG(level)%Arr(i-1,j,k)*CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             elseif (BtypeW==PERIODIC) then
                       p=p+PhiMG(level)%Arr(CoefMG(level)%nx-1,j,k)*CoefMG(level)%Aw
                       Ap=Ap+CoefMG(level)%Aw
             endif
             if (i<CoefMG(level)%nx) then
                       p=p+PhiMG(level)%Arr(i+1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             elseif (BtypeE==PERIODIC) then
                       p=p+PhiMG(level)%Arr(1,j,k)*CoefMG(level)%Ae
                       Ap=Ap+CoefMG(level)%Ae
             endif
             if (j>0) then
                       p=p+PhiMG(level)%Arr(i,j-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             elseif (BtypeS==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,CoefMG(level)%ny-1,k)*CoefMG(level)%As
                       Ap=Ap+CoefMG(level)%As
             endif
             if (j<CoefMG(level)%ny) then
                       p=p+PhiMG(level)%Arr(i,j+1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             elseif (BtypeN==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,1,k)*CoefMG(level)%An
                       Ap=Ap+CoefMG(level)%An
             endif
             if (k>0) then
                       p=p+PhiMG(level)%Arr(i,j,k-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             elseif (BtypeB==PERIODIC) then
                       p=p+PhiMG(level)%Arr(i,j,CoefMG(level)%nz-1)*CoefMG(level)%Ab
                       Ap=Ap+CoefMG(level)%Ab
             endif
             if (k<CoefMG(level)%nz) then
                       p=p+PhiMG(level)%Arr(i,j,k+1)*CoefMG(level)%At
                       Ap=Ap+CoefMG(level)%At
             elseif (BtypeT==PERIODIC) then
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

  do l=level,minmglevel,-1
   PhiMG(l)%Arr=0
   RHSMG(l)%Arr=0
   ResMG(l)%Arr=0
  enddo
endsubroutine



recursive subroutine MG_CGC(level,eps,ncgc,npre,npost,R)
integer,intent(in):: level,ncgc,npre,npost
real(KND),intent(in):: eps
real(KND),intent(out):: R
real(KND) R1
integer k

  if (level == minmglevel) then !1
        call MG_GE(level)
        call MG_res(level,R)
!        write (*,*) "GEres",sqrt(R)

     if (R>mgepsinnerGS**2) then
      R1=R
      do k=1,mgmaxinnerGSiter
        call MG_GS(level, 10)
        call MG_res(level,R)
!         write (*,*) "GSres",sqrt(R)
        if (R < mgepsinnerGS**2) exit
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

!   write (*,*) "cgc res",level,R
endsubroutine MG_CGC


subroutine POISSMG(Phi,RHS)        !Solves Poisson equation using Successive over-relaxation
real(KND),dimension(0:,0:,0:),intent(inout):: Phi
real(KND),dimension(1:,1:,1:),intent(in)::RHS
integer i,j,k,l,nx,ny,nz,sx,sy,sz
real(KND) mgeps,R,Phiref
real(KND),save:: called=0

 mgeps=epsPoisson
 Phi=0

 if (BtypeE==PERIODIC) then
  sx=1
 else
  sx=0
 endif
 if (BtypeN==PERIODIC) then
  sy=1
 else
  sy=0
 endif
 if (BtypeT==PERIODIC) then
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
    allocate(PhiMG(l)%Arr(0,0,0))
    allocate(RHSMG(l)%Arr(0,0,0))
    allocate(ResMG(l)%Arr(0,0,0))
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

 if (BtypeE==PERIODIC) then
  RHSMG(LMG)%Arr(0,:,:)=RHSMG(LMG)%Arr(Prnx,:,:)
  PhiMG(LMG)%Arr(0,:,:)=PhiMG(LMG)%Arr(Prnx,:,:)
 endif
 if (BtypeN==PERIODIC) then
  RHSMG(LMG)%Arr(:,0,:)=RHSMG(LMG)%Arr(:,Prny,:)
  PhiMG(LMG)%Arr(:,0,:)=PhiMG(LMG)%Arr(:,Prny,:)
 endif
 if (BtypeT==PERIODIC) then
  RHSMG(LMG)%Arr(:,:,0)=RHSMG(LMG)%Arr(:,:,Prnz)
  PhiMG(LMG)%Arr(:,:,0)=PhiMG(LMG)%Arr(:,:,Prnz)
 endif

 do l=1,maxPoissoniter
  write (*,*) "MG iteration:",l
   Phiref=SUM(PhiMG(LMG)%Arr(0:CoefMG(LMG)%nx,0:CoefMG(LMG)%ny,0:CoefMG(LMG)%nz))/&
          ((CoefMG(LMG)%nx+1)*(CoefMG(LMG)%ny+1)*(CoefMG(LMG)%nz+1))
   do k=0,CoefMG(LMG)%nz
    do j=0,CoefMG(LMG)%ny
     do i=0,CoefMG(LMG)%nx
      PhiMG(LMG)%Arr(i,j,k)=PhiMG(LMG)%Arr(i,j,k)-Phiref
     enddo
    enddo
   enddo
   call Bound_Phi_MG(PhiMG(LMG)%Arr,CoefMG(LMG)%nx,CoefMG(LMG)%ny,CoefMG(LMG)%nz)

  !$hmpp <GSKernels> allocate
  call MG_CGC(LMG,mgeps,mgncgc,mgnpre,mgnpost,R)
  !$hmpp <GSKernels> release

  if (R<mgeps)   exit
 write (*,*) "MG residuum",R
 enddo
 write (*,*) "MG residuum",R

 Phi(1:Prnx,1:Prny,1:Prnz)=PhiMG(LMG)%Arr(0+sx:nx,0+sy:ny,0+sz:nz)

endsubroutine POISSMG














endmodule MULTIGRID


!GPU Codelets. More or less like externals. Inside the module only because of KND.

  !$hmpp <GSKernels> group, target=CUDA
  !$hmpp <GSKernels> map args[GS0::Phi;Res0::Phi]
  !$hmpp <GSKernels> map args[GS0::RHS;Res0::RHS]

  !$hmpp <GSKernels> map args[GS1::Phi;Res1::Phi]
  !$hmpp <GSKernels> map args[GS1::RHS;Res1::RHS]
 
  !$hmpp <GSKernels> map args[GS2::Phi;Res2::Phi]
  !$hmpp <GSKernels> map args[GS2::RHS;Res2::RHS]

  !$hmpp <GSKernels> map args[GS3::Phi;Res3::Phi]
  !$hmpp <GSKernels> map args[GS3::RHS;Res3::RHS]

  !$hmpp <GSKernels> map args[GS4::Phi;Res4::Phi]
  !$hmpp <GSKernels> map args[GS4::RHS;Res4::RHS]

  !$hmpp <GSKernels> map args[GS5::Phi;Res5::Phi]
  !$hmpp <GSKernels> map args[GS5::RHS;Res5::RHS]

  !$hmpp <GSKernels> map args[GS6::Phi;Res6::Phi]
  !$hmpp <GSKernels> map args[GS6::RHS;Res6::RHS]

  !$hmpp <GSKernels> map args[GS7::Phi;Res7::Phi]
  !$hmpp <GSKernels> map args[GS7::RHS;Res7::RHS]

  !$hmpp <GSKernels> map args[GS8::Phi;Res8::Phi]
  !$hmpp <GSKernels> map args[GS8::RHS;Res8::RHS]

  !$hmpp <GSKernels> mapbyname, BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT

!$hmpp <GSKernels> GS0 codelet
!$hmpp <GSKernels> GS1 codelet
!$hmpp <GSKernels> GS2 codelet
!$hmpp <GSKernels> GS3 codelet
!$hmpp <GSKernels> GS4 codelet
!$hmpp <GSKernels> GS5 codelet
!$hmpp <GSKernels> GS6 codelet
!$hmpp <GSKernels> GS7 codelet
!$hmpp <GSKernels> GS8 codelet
subroutine MG_GS_GPU(nx,ny,nz,nit,Phi,RHS,&
                     Aw,Ae,As,An,Ab,At,&
                     BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT)
   implicit none

   integer,parameter:: PERIODIC=3

   integer,intent(in)   :: nx,ny,nz,nit
   real(4),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(inout)::Phi
   real(4),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(in)::RHS
   real(4),intent(in) :: Aw,Ae,As,An,Ab,At
   integer(4),intent(in) :: BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT
   integer i,j,k,l
   real(4) :: p,Ap

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
             elseif (BtypeW==PERIODIC) then
                       p=p+Phi(nx-1,j,k)*Aw
                       Ap=Ap+Aw
             endif
             if (i<nx) then
                       p=p+Phi(i+1,j,k)*Ae
                       Ap=Ap+Ae
             elseif (BtypeE==PERIODIC) then
                       p=p+Phi(1,j,k)*Ae
                       Ap=Ap+Ae
             endif
             if (j>0) then
                       p=p+Phi(i,j-1,k)*As
                       Ap=Ap+As
             elseif (BtypeS==PERIODIC) then
                       p=p+Phi(i,ny-1,k)*As
                       Ap=Ap+As
             endif
             if (j<ny) then
                       p=p+Phi(i,j+1,k)*An
                       Ap=Ap+An
             elseif (BtypeN==PERIODIC) then
                       p=p+Phi(i,1,k)*An
                       Ap=Ap+An
             endif
             if (k>0) then
                       p=p+Phi(i,j,k-1)*Ab
                       Ap=Ap+Ab
             elseif (BtypeB==PERIODIC) then
                       p=p+Phi(i,j,nz-1)*Ab
                       Ap=Ap+Ab
             endif
             if (k<nz) then
                       p=p+Phi(i,j,k+1)*At
                       Ap=Ap+At
             elseif (BtypeT==PERIODIC) then
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
             elseif (BtypeW==PERIODIC) then
                       p=p+Phi(nx-1,j,k)*Aw
                       Ap=Ap+Aw
             endif
             if (i<nx) then
                       p=p+Phi(i+1,j,k)*Ae
                       Ap=Ap+Ae
             elseif (BtypeE==PERIODIC) then
                       p=p+Phi(1,j,k)*Ae
                       Ap=Ap+Ae
             endif
             if (j>0) then
                       p=p+Phi(i,j-1,k)*As
                       Ap=Ap+As
             elseif (BtypeS==PERIODIC) then
                       p=p+Phi(i,ny-1,k)*As
                       Ap=Ap+As
             endif
             if (j<ny) then
                       p=p+Phi(i,j+1,k)*An
                       Ap=Ap+An
             elseif (BtypeN==PERIODIC) then
                       p=p+Phi(i,1,k)*An
                       Ap=Ap+An
             endif
             if (k>0) then
                       p=p+Phi(i,j,k-1)*Ab
                       Ap=Ap+Ab
             elseif (BtypeB==PERIODIC) then
                       p=p+Phi(i,j,nz-1)*Ab
                       Ap=Ap+Ab
             endif
             if (k<nz) then
                       p=p+Phi(i,j,k+1)*At
                       Ap=Ap+At
             elseif (BtypeT==PERIODIC) then
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


!$hmpp <GSKernels> Res0 codelet
!$hmpp <GSKernels> Res1 codelet
!$hmpp <GSKernels> Res2 codelet
!$hmpp <GSKernels> Res3 codelet
!$hmpp <GSKernels> Res4 codelet
!$hmpp <GSKernels> Res5 codelet
!$hmpp <GSKernels> Res6 codelet
!$hmpp <GSKernels> Res7 codelet
!$hmpp <GSKernels> Res8 codelet
subroutine MG_Res_GPU(nx,ny,nz,Phi,Res,RHS,&
                     Aw,Ae,As,An,Ab,At,&
                     BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT,R)

   implicit none


   integer,parameter:: PERIODIC=3

   integer,intent(in)   :: nx,ny,nz
   real(4),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(in)::Phi
   real(4),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(out)::Res
   real(4),dimension(-1:nx+1,-1:ny+1,-1:nz+1),intent(in)::RHS
   real(4),intent(in) :: Aw,Ae,As,An,Ab,At
   integer(4),intent(in) :: BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT
   real(4),intent(out) :: R
   integer i,j,k,l
   real(4) :: p,Ap

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
             elseif (BtypeW==PERIODIC) then
                       p=p+Phi(nx-1,j,k)*Aw
                       Ap=Ap+Aw
             endif
             if (i<nx) then
                       p=p+Phi(i+1,j,k)*Ae
                       Ap=Ap+Ae
             elseif (BtypeE==PERIODIC) then
                       p=p+Phi(1,j,k)*Ae
                       Ap=Ap+Ae
             endif
             if (j>0) then
                       p=p+Phi(i,j-1,k)*As
                       Ap=Ap+As
             elseif (BtypeS==PERIODIC) then
                       p=p+Phi(i,ny-1,k)*As
                       Ap=Ap+As
             endif
             if (j<ny) then
                       p=p+Phi(i,j+1,k)*An
                       Ap=Ap+An
             elseif (BtypeN==PERIODIC) then
                       p=p+Phi(i,1,k)*An
                       Ap=Ap+An
             endif
             if (k>0) then
                       p=p+Phi(i,j,k-1)*Ab
                       Ap=Ap+Ab
             elseif (BtypeB==PERIODIC) then
                       p=p+Phi(i,j,nz-1)*Ab
                       Ap=Ap+Ab
             endif
             if (k<nz) then
                       p=p+Phi(i,j,k+1)*At
                       Ap=Ap+At
             elseif (BtypeT==PERIODIC) then
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





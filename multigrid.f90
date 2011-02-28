module MULTIGRID

use PARAMETERS
implicit none

  type TMGArr
    real(KND),allocatable,dimension(:,:,:):: Arr
  endtype

  type TCoefs
   real(KND) dx,dy,dz,Aw,Ae,As,An,Ab,At
   integer nx,ny,nz
  endtype

  type TRedBlack
   logical,allocatable,dimension(:,:,:):: Arr
  endtype

  type(TMGArr),allocatable,dimension(:):: PhiMG,RHSMG,ResMG
  type(TCoefs),allocatable,dimension(:):: CoefMg
  type(TRedBlack),allocatable,dimension(:):: RedBlackMG

  integer bnx,bny,bnz !Base cube dimensions for Multigrid, Prnx=bnx*2**level
  integer LMG !depth of multigrid
  integer minmglevel !innermost MG level
  integer mgncgc !type of cycling 1..V cycle, 2..W cycle
  integer mgnpre
  integer mgnpost
  integer mgmaxinnnerGSiter
  real(KND) mgepsinnerGS


!type (fishworkspace) :: WFISH
!real(KND),allocatable::WFISH(:)
contains 


pure subroutine Prolongate(AFine,ACoarse,level)
 integer,intent(IN):: level
 real(KND),dimension(-1:,-1:,-1:),intent(IN):: ACoarse
 real(KND),dimension(-1:,-1:,-1:),intent(INOUT):: AFine
 integer:: i,j,k,nx,ny,nz
 
   nx=bnx*2**level !level means from which grid we interpolate
   ny=bny*2**level
   nz=bnz*2**level

   do k=0,nz
    do j=0,ny
     do i=0,nx
                     AFine(2*i,2*j,2*k)=AFine(2*i,2*j,2*k)+ACoarse(i,j,k)
!                      if (i<nx) AFine(2*i+1,2*j,2*k)=AFine(2*i+1,2*j,2*k)+0.5*ACoarse(i,j,k)
!                      if (i>0) AFine(2*i-1,2*j,2*k)=AFine(2*i-1,2*j,2*k)+0.5*ACoarse(i,j,k)
!                      if (j<ny) AFine(2*i,2*j+1,2*k)=AFine(2*i,2*j+1,2*k)+0.5*ACoarse(i,j,k)
!                      if (j>0) AFine(2*i,2*j-1,2*k)=AFine(2*i,2*j-1,2*k)+0.5*ACoarse(i,j,k)
!                      if (k<nz) AFine(2*i,2*j,2*k+1)=AFine(2*i,2*j,2*k+1)+0.5*ACoarse(i,j,k)
!                      if (k>0) AFine(2*i,2*j,2*k-1)=AFine(2*i,2*j,2*k-1)+0.5*ACoarse(i,j,k)
!                      if((i<nx).and.(j<ny)) AFine(2*i+1,2*j+1,2*k)=AFine(2*i+1,2*j+1,2*k)+0.25*ACoarse(i,j,k)
!                      if((i<nx).and.(k<nz)) AFine(2*i+1,2*j,2*k+1)=AFine(2*i+1,2*j,2*k+1)+0.25*ACoarse(i,j,k)
!                      if((j<ny).and.(k<nz)) AFine(2*i,2*j+1,2*k+1)=AFine(2*i,2*j+1,2*k+1)+0.25*ACoarse(i,j,k)
!                      if((i>0).and.(j>0)) AFine(2*i,2*j+1,2*k+1)=AFine(2*i,2*j+1,2*k+1)+0.25*ACoarse(i,j,k)
!                      if((i>0).and.(k>0)) AFine(2*i,2*j+1,2*k+1)=AFine(2*i,2*j+1,2*k+1)+0.25*ACoarse(i,j,k)
!                      if((j>0).and.(k>0)) AFine(2*i,2*j-1,2*k-1)=AFine(2*i,2*j-1,2*k-1)+0.25*ACoarse(i,j,k)
!                      if((i>0).and.(j<ny)) AFine(2*i-1,2*j+1,2*k)=AFine(2*i-1,2*j+1,2*k)+0.25*ACoarse(i,j,k)
!                      if((i>0).and.(k<nz)) AFine(2*i-1,2*j,2*k+1)=AFine(2*i-1,2*j,2*k+1)+0.25*ACoarse(i,j,k)
!                      if((i<nx).and.(j>0)) AFine(2*i+1,2*j-1,2*k)=AFine(2*i+1,2*j-1,2*k)+0.25*ACoarse(i,j,k)
!                      if((j>0).and.(k<nz)) AFine(2*i,2*j-1,2*k+1)=AFine(2*i,2*j-1,2*k+1)+0.25*ACoarse(i,j,k)
!                      if((i<nx).and.(k>0)) AFine(2*i+1,2*j,2*k-1)=AFine(2*i+1,2*j,2*k-1)+0.25*ACoarse(i,j,k)
!                      if((j<ny).and.(k>0)) AFine(2*i,2*j+1,2*k-1)=AFine(2*i,2*j+1,2*k-1)+0.25*ACoarse(i,j,k)
!                      if((i<nx).and.(j<ny).and.(k<nz)) AFine(2*i+1,2*j+1,2*k+1)=AFine(2*i+1,2*j+1,2*k+1)+0.125*ACoarse(i,j,k)
!                      if((i>0).and.(j<ny).and.(k<nz)) AFine(2*i-1,2*j+1,2*k+1)=AFine(2*i-1,2*j+1,2*k+1)+0.125*ACoarse(i,j,k)
!                      if((i<nx).and.(j>0).and.(k<nz)) AFine(2*i+1,2*j-1,2*k+1)=AFine(2*i+1,2*j-1,2*k+1)+0.125*ACoarse(i,j,k)
!                      if((i<nx).and.(j<ny).and.(k>0)) AFine(2*i+1,2*j+1,2*k-1)=AFine(2*i+1,2*j+1,2*k-1)+0.125*ACoarse(i,j,k)
!                      if((i>0).and.(j>0).and.(k<nz)) AFine(2*i-1,2*j-1,2*k+1)=AFine(2*i-1,2*j-1,2*k+1)+0.125*ACoarse(i,j,k)
!                      if((i>0).and.(j<ny).and.(k>0)) AFine(2*i-1,2*j+1,2*k-1)=AFine(2*i-1,2*j+1,2*k-1)+0.125*ACoarse(i,j,k)
!                      if((i<nx).and.(j>0).and.(k>0)) AFine(2*i+1,2*j-1,2*k-1)=AFine(2*i+1,2*j-1,2*k-1)+0.125*ACoarse(i,j,k)
!                      if((i>0).and.(j>0).and.(k>0)) AFine(2*i-1,2*j-1,2*k-1)=AFine(2*i-1,2*j-1,2*k-1)+0.125*ACoarse(i,j,k)
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
!    if (BtypeE==PERIODIC) Afine(0,:,:)=Afine(2*nx,:,:)
!    if (BtypeN==PERIODIC) Afine(:,0,:)=Afine(:,2*ny,:)
!    if (BtypeT==PERIODIC) Afine(:,:,0)=Afine(:,:,2*nz)
   call Bound_Phi_MG(Afine,2*nx,2*ny,2*nz)
endsubroutine Prolongate



 subroutine Restrict(ACoarse,AFine,level)
 integer,intent(IN):: level
 real(KND),dimension(-1:,-1:,-1:),intent(OUT):: ACoarse
 real(KND),dimension(-1:,-1:,-1:),intent(INOUT):: AFine
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





pure subroutine Prolongate2(AFine,ACoarse,level)
 integer,intent(IN):: level
 real(KND),dimension(-1:,-1:,-1:),intent(IN):: ACoarse
 real(KND),dimension(-1:,-1:,-1:),intent(INOUT):: AFine
 real(KND),dimension(LBOUND(AFine,1):UBOUND(AFine,1),LBOUND(AFine,1):UBOUND(AFine,1),LBOUND(AFine,1):UBOUND(AFine,1)):: Mask
 integer:: i,j,k,nx,ny,nz
 
   nx=bnx*2**level !level means from which grid we interpolate
   ny=bny*2**level
   nz=bnz*2**level

   do k=0,nz
    do j=0,ny
     do i=0,nx
                     AFine(2*i,2*j,2*k)=AFine(2*i,2*j,2*k)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i+1,2*j,2*k)=AFine(2*i+1,2*j,2*k)+ACoarse(i,j,k)
                     Mask(2*i+1,2*j,2*k)=Mask(2*i+1,2*j,2*k)+1
                     AFine(2*i-1,2*j,2*k)=AFine(2*i-1,2*j,2*k)+ACoarse(i,j,k)
                     Mask(2*i-1,2*j,2*k)=Mask(2*i-1,2*j,2*k)+1
                     AFine(2*i,2*j+1,2*k)=AFine(2*i,2*j+1,2*k)+ACoarse(i,j,k)
                     Mask(2*i,2*j+1,2*k)=Mask(2*i,2*j+1,2*k)+1
                     AFine(2*i,2*j-1,2*k)=AFine(2*i,2*j-1,2*k)+ACoarse(i,j,k)
                     Mask(2*i,2*j-1,2*k)=Mask(2*i,2*j-1,2*k)+1
                     AFine(2*i,2*j,2*k+1)=AFine(2*i,2*j,2*k+1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k+1)=Mask(2*i,2*j,2*k+1)+1
                     AFine(2*i,2*j,2*k-1)=AFine(2*i,2*j,2*k-1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k-1)=Mask(2*i,2*j,2*k-1)+1
                     AFine(2*i+1,2*j+1,2*k)=AFine(2*i+1,2*j+1,2*k)+ACoarse(i,j,k)
                     Mask(2*i+1,2*j+1,2*k)=Mask(2*i+1,2*j+1,2*k)+1
                     AFine(2*i+1,2*j,2*k+1)=AFine(2*i+1,2*j,2*k+1)+ACoarse(i,j,k)
                     Mask(2*i+1,2*j,2*k+1)=Mask(2*i+1,2*j,2*k+1)+1
                     AFine(2*i,2*j+1,2*k+1)=AFine(2*i,2*j+1,2*k+1)+ACoarse(i,j,k)
                     Mask(2*i,2*j+1,2*k+1)=Mask(2*i,2*j+1,2*k+1)+1
                     AFine(2*i,2*j+1,2*k+1)=AFine(2*i,2*j+1,2*k+1)+ACoarse(i,j,k)
                     Mask(2*i,2*j+1,2*k+1)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i,2*j+1,2*k+1)=AFine(2*i,2*j+1,2*k+1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i,2*j-1,2*k-1)=AFine(2*i,2*j-1,2*k-1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i-1,2*j+1,2*k)=AFine(2*i-1,2*j+1,2*k)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i-1,2*j,2*k+1)=AFine(2*i-1,2*j,2*k+1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i+1,2*j-1,2*k)=AFine(2*i+1,2*j-1,2*k)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i,2*j-1,2*k+1)=AFine(2*i,2*j-1,2*k+1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i+1,2*j,2*k-1)=AFine(2*i+1,2*j,2*k-1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i,2*j+1,2*k-1)=AFine(2*i,2*j+1,2*k-1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i+1,2*j+1,2*k+1)=AFine(2*i+1,2*j+1,2*k+1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i-1,2*j+1,2*k+1)=AFine(2*i-1,2*j+1,2*k+1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i+1,2*j-1,2*k+1)=AFine(2*i+1,2*j-1,2*k+1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i+1,2*j+1,2*k-1)=AFine(2*i+1,2*j+1,2*k-1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i-1,2*j-1,2*k+1)=AFine(2*i-1,2*j-1,2*k+1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i-1,2*j+1,2*k-1)=AFine(2*i-1,2*j+1,2*k-1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i+1,2*j-1,2*k-1)=AFine(2*i+1,2*j-1,2*k-1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
                     AFine(2*i-1,2*j-1,2*k-1)=AFine(2*i-1,2*j-1,2*k-1)+ACoarse(i,j,k)
                     Mask(2*i,2*j,2*k)=Mask(2*i,2*j,2*k)+1
     enddo
    enddo
   enddo
!    if (BtypeE==PERIODIC) Afine(0,:,:)=Afine(2*nx,:,:)
!    if (BtypeN==PERIODIC) Afine(:,0,:)=Afine(:,2*ny,:)
!    if (BtypeT==PERIODIC) Afine(:,:,0)=Afine(:,:,2*nz)
   call Bound_Phi_MG(Afine,2*nx,2*ny,2*nz)
endsubroutine Prolongate2



pure subroutine Restrict2(ACoarse,AFine,level)
 integer,intent(IN):: level
 real(KND),dimension(-1:,-1:,-1:),intent(OUT):: ACoarse
 real(KND),dimension(-1:,-1:,-1:),intent(INOUT):: AFine
 real(KND) q
 integer:: i,j,k,nx,ny,nz
 
   nx=bnx*2**level !level means on which grid we restrict
   ny=bny*2**level
   nz=bnz*2**level

   call Bound_Phi_MG(Afine,2*nx,2*ny,2*nz)

   do k=0,nz
    do j=0,ny
     do i=0,nx
        ACoarse(i,j,k)=AFine(2*i,2*j,2*k); q=1
!         if (i<nx) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i+1,2*j,2*k); q=q+1
! endif
!         if (i>0) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i-1,2*j,2*k); q=q+1
! endif
!         if (j<ny) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i,2*j+1,2*k); q=q+1
! endif
!         if (j>0) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i,2*j-1,2*k); q=q+1
! endif
!         if (k<nz) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i,2*j,2*k+1); q=q+1
! endif
!         if (k>0) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i,2*j,2*k-1); q=q+1
! endif
!         if ((i<nx).and.(j<ny)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i+1,2*j+1,2*k); q=q+1
! endif
!         if ((i<nx).and.(k<nz)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i+1,2*j,2*k+1); q=q+1
! endif
!         if ((j<ny).and.(k<nz)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i,2*j+1,2*k+1); q=q+1
! endif
!         if ((i>0).and.(j>0)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i-1,2*j-1,2*k); q=q+1
! endif
!         if ((i>0).and.(k>0)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i-1,2*j,2*k-1); q=q+1
! endif
!         if ((j>0).and.(k>0)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i,2*j-1,2*k-1); q=q+1
! endif
!         if ((i>0).and.(j<ny)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i-1,2*j+1,2*k); q=q+1
! endif
!         if ((i>0).and.(k<nz)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i-1,2*j,2*k+1); q=q+1
! endif
!         if ((i<nx).and.(j>0)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i+1,2*j-1,2*k); q=q+1
! endif
!         if ((j>0).and.(k<nz)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i,2*j-1,2*k+1); q=q+1
! endif
!         if ((i<nx).and.(k>0)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i+1,2*j,2*k-1); q=q+1
! endif
!         if ((j<ny).and.(k>0)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i,2*j+1,2*k-1); q=q+1
! endif
!         if ((i<nx).and.(j<ny).and.(k<nz)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i+1,2*j+1,2*k+1); q=q+1
! endif
!         if ((i>0).and.(j<ny).and.(k<nz)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i-1,2*j+1,2*k+1); q=q+1
! endif
!         if ((i<nx).and.(j>0).and.(k<nz)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i+1,2*j-1,2*k+1); q=q+1
! endif
!         if ((i<nx).and.(j<ny).and.(k>0)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i+1,2*j+1,2*k-1); q=q+1
! endif
!         if ((i>0).and.(j>0).and.(k<nz)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i-1,2*j-1,2*k+1); q=q+1
! endif
!         if ((i>0).and.(j<ny).and.(k>0)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i-1,2*j+1,2*k-1); q=q+1
! endif
!         if ((i<nx).and.(j>0).and.(k>0)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i+1,2*j-1,2*k-1); q=q+1
! endif
!         if ((i>0).and.(j>0).and.(k>0)) then
 ACoarse(i,j,k)=ACoarse(i,j,k)+AFine(2*i-1,2*j-1,2*k-1); q=q+1
! endif
        ACoarse(i,j,k)=ACoarse(i,j,k)/q             
     enddo
    enddo
   enddo
!    if (BtypeE==PERIODIC) ACoarse(0,:,:)=ACoarse(nx,:,:)
!    if (BtypeN==PERIODIC) ACoarse(:,0,:)=ACoarse(:,ny,:)
!    if (BtypeT==PERIODIC) ACoarse(:,:,0)=ACoarse(:,:,nz)
endsubroutine Restrict2




pure subroutine BOUND_Phi_MG(Phi,nx,ny,nz)
  real(KND),intent(INOUT):: Phi(-1:,-1:,-1:)
  integer,intent(IN):: nx,ny,nz
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
integer,intent(IN):: level,i,j,k,nulx,nuly,nulz
 ind=(1-nulx)+i+(j-nuly)*(CoefMG(level)%nx+(1-nulx))+(k-nulz)*(CoefMG(level)%nx+1-nulx)*(CoefMG(level)%ny+1-nuly)
endfunction ind




subroutine MG_GE(level)
integer,intent(IN)::level
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
     STOP
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
     STOP
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
integer,intent(IN)::level
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
     STOP
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
     STOP
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
integer,intent(IN)::level
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
     STOP
     write (*,*) "info",info
    endif

    allocate(work(1))

    if (KND==DBL) then
     call DGETRI(nxyz,age,nxyz,ipivot,work,-1,info)
    else
     call SGETRI(nxyz,age,nxyz,ipivot,work,-1,info)
    endif
     
    if (info/=0) then
     STOP
     write (*,*) "info",info
    endif

    ldwork=work(1)
    deallocate(work)
    allocate(work(ldwork))

    if (KND==DBL) then
     call DGETRI(nxyz,age,nxyz,ipivot,work,ldwork,info)
    else
     call SGETRI(nxyz,age,nxyz,ipivot,work,ldwork,info)
    endif
     
    if (info/=0) then
     STOP
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
integer,intent(IN)::level,niter
! real(KND),intent(OUT):: R
integer i,j,k,l
real(KND) p,Ap

   do l=1,niter
!     S=0
    !$OMP PARALLEL PRIVATE(i,j,k,p,Ap)
    !$OMP DO
    do k=0,CoefMG(level)%nz
        do j=0,CoefMG(level)%ny
            do i=0,CoefMG(level)%nx
             if (RedBlackMG(level)%Arr(i,j,k)) then
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
!              S=S+(p-Ap*PhiMG(level)%Arr(i,j,k))**2
!              if (l==niter) ResMG(level)%Arr(i,j,k)=-p +Ap*PhiMG(level)%Arr(i,j,k)
             p=p/Ap
             PhiMG(level)%Arr(i,j,k)=p
             endif
            enddo
        enddo
    enddo
    !$OMP ENDDO
    !$OMP DO
    do k=0,CoefMG(level)%nz
        do j=0,CoefMG(level)%ny
            do i=0,CoefMG(level)%nx
             if (.not.RedBlackMG(level)%Arr(i,j,k)) then
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
!              S=S+(p-Ap*PhiMG(level)%Arr(i,j,k))**2
!              if (l==niter) ResMG(level)%Arr(i,j,k)=-p +Ap*PhiMG(level)%Arr(i,j,k)
             p=p/Ap
             PhiMG(level)%Arr(i,j,k)=p
             endif
            enddo
        enddo
    enddo
    !$OMP ENDDO
    !$OMP END PARALLEL
!     write(*,*) "GS residuum",S
!     call BOUND_Phi_MG(PhiMG(level)%Arr,nx-1,ny-1,nz-1)
   enddo  
!    R=MAXVAL(ResMG(level)%Arr(0:CoefMG(level)%nx,0:CoefMG(level)%ny,0:CoefMG(level)%nz))

!     write(*,*) "GS residuum",sqrt(S)
!      call BOUND_Phi_MG(PhiMG(level)%Arr,nx-1,ny-1,nz-1)
endsubroutine MG_GS

subroutine MG_res(level,R)
integer,intent(IN)::level
real(KND),intent(OUT)::R
integer i,j,k
real(KND),save:: p,Ap
character(70):: str

     call BOUND_Phi_MG(PhiMG(level)%Arr,CoefMG(level)%nx,CoefMG(level)%ny,CoefMG(level)%nz)
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
!      write (*,*) "residuum",Sqrt(R)
!  if (level<LMG) goto 50
!   OPEN(11,file="Res.vtk")
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
!       Write (11,*) ResMG(level)%Arr(i,j,k)
!     enddo
!    enddo
!   enddo 
!   write (11,*)
!   CLOSE(11)
!   50 i=1
endsubroutine MG_res


subroutine MG_clear(level)
integer,intent(IN):: level
integer l

  do l=level,minmglevel,-1
   PhiMG(l)%Arr=0
   RHSMG(l)%Arr=0
   ResMG(l)%Arr=0
  enddo
endsubroutine


subroutine Filter(A,level)
 real(KND),dimension(-1:,-1:,-1:),intent(INOUT):: A
 real(KND),dimension(LBOUND(A,1):UBOUND(A,1),LBOUND(A,1):UBOUND(A,1),LBOUND(A,1):UBOUND(A,1)):: B
 integer,intent(IN):: level
 integer nx,ny,nz,i,j,k,ii,jj,kk
 real(KND) S 
 nx=bnx*2**level
 ny=bny*2**level
 nz=bnz*2**level

 do k=0,nz
  do j=0,ny
   do i=0,nx
    S=0
    do kk=-1,1
     do jj=-1,1
      do ii=-1,1
       if ((i/=0).or.(j/=0).or.(k/=0))  S=S+A(i+ii,j+jj,k+kk)
      enddo
     enddo
    enddo
    S=S/26._KND
    B(i,j,k)=A(i,j,k)/.8_KND+.2_KND*S
   enddo
  enddo
 enddo
 A=B
 call Bound_Phi_MG(A,nx,ny,nz)
endsubroutine


recursive subroutine MG_CGC(level,eps,ncgc,npre,npost,R)
integer,intent(IN):: level,ncgc,npre,npost
real(KND),intent(IN):: eps
real(KND),intent(OUT):: R
real(KND) R1
integer k
!   write (*,*) "cgc",level
  if (level == minmglevel) then !1
        call MG_INV(level)
        call MG_res(level,R)
!        write (*,*) "GEres",sqrt(R)

     if (R>mgepsinnerGS**2) then
      R1=R
      do k=1,mgmaxinnnerGSiter
        call MG_GS(level, 10)
        call MG_res(level,R)
!         write (*,*) "GSres",sqrt(R)
        if (R < mgepsinnerGS**2) EXIT
        if (R1<R*1.05_KND) EXIT
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
real(KND),dimension(0:,0:,0:),intent(INOUT):: Phi
real(KND),dimension(1:,1:,1:),intent(IN)::RHS
integer i,j,k,l,nx,ny,nz,sx,sy,sz
real(KND) mgeps,R,Phiref
real(KND),save:: called=0

!  if (KND==SNG) then
!     mgeps=1E-7
!  else
!     mgeps=1E-14
!  endif
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
  allocate(PhiMG(minmglevel:LMG),RHSMG(minmglevel:LMG),ResMG(minmglevel:LMG),CoefMG(minmglevel:LMG),RedBlackMG(minmglevel:LMG))
  do l=LMG,minmglevel,-1
   CoefMG(l)%nx=bnx*2**l
   CoefMG(l)%ny=bny*2**l
   CoefMG(l)%nz=bnz*2**l
   CoefMG(l)%dx=dxmin*2**(LMG-l)
   CoefMG(l)%dy=dymin*2**(LMG-l)
   CoefMG(l)%dz=dzmin*2**(LMG-l)
   allocate(PhiMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%ny+1,-1:CoefMG(l)%nz+1))
   allocate(RHSMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%ny+1,-1:CoefMG(l)%nz+1))
   allocate(ResMG(l)%Arr(-1:CoefMG(l)%nx+1,-1:CoefMG(l)%ny+1,-1:CoefMG(l)%nz+1))
   allocate(RedBlackMG(l)%Arr(-1:CoefMG(l)%nx,-1:CoefMG(l)%ny,-1:CoefMG(l)%nz))
   CoefMG(l)%Ae=1._KND/(CoefMG(l)%dx*CoefMG(l)%dx)
   CoefMG(l)%Aw=1._KND/(CoefMG(l)%dx*CoefMG(l)%dx)
   CoefMG(l)%An=1._KND/(CoefMG(l)%dy*CoefMG(l)%dy)
   CoefMG(l)%As=1._KND/(CoefMG(l)%dy*CoefMG(l)%dy)
   CoefMG(l)%At=1._KND/(CoefMG(l)%dz*CoefMG(l)%dz)
   CoefMG(l)%Ab=1._KND/(CoefMG(l)%dz*CoefMG(l)%dz)

   RedBlackMG(l)%Arr(-1,-1,-1)=.true.
   do k=0,CoefMG(l)%nz
    RedBlackMG(l)%Arr(-1,-1,k)=.not.RedBlackMG(l)%Arr(-1,-1,k-1)
    do j=0,CoefMG(l)%ny
     RedBlackMG(l)%Arr(-1,j,k)=.not.RedBlackMG(l)%Arr(-1,j-1,k)
     do i=0,CoefMG(l)%nx
      RedBlackMG(l)%Arr(i,j,k)=.not.RedBlackMG(l)%Arr(i-1,j,k)
     enddo
    enddo
   enddo

  enddo
  called=1
 endif

 nx=bnx*2**LMG
 ny=bny*2**LMG
 nz=bnz*2**LMG

 if (nx-sx/=Prnx-1.or.ny-sy/=Prny-1.or.nz-sz/=Prnz-1) then
    write (*,*) 0+sx,":",nx,"--",1,":",Prnx
    write (*,*) 0+sy,":",ny,"--",1,":",Prny
    write (*,*) 0+sz,":",nz,"--",1,":",Prnz
    STOP
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
          ((CoefMG(LMG)%nx+1)*(CoefMG(LMG)%ny+1)*(CoefMG(LMG)%ny+1))
   do k=0,CoefMG(LMG)%nz
    do j=0,CoefMG(LMG)%ny
     do i=0,CoefMG(LMG)%nx
      PhiMG(LMG)%Arr(i,j,k)=PhiMG(LMG)%Arr(i,j,k)-Phiref
     enddo
    enddo
   enddo
   call Bound_Phi_MG(PhiMG(LMG)%Arr,CoefMG(LMG)%nx,CoefMG(LMG)%ny,CoefMG(LMG)%nz)
  call MG_CGC(LMG,mgeps,mgncgc,mgnpre,mgnpost,R)

  if (R<mgeps)   EXIT
 write (*,*) "MG residuum",R
 enddo
 write (*,*) "MG residuum",R
 Phi(1:Prnx,1:Prny,1:Prnz)=PhiMG(LMG)%Arr(0+sx:nx,0+sy:ny,0+sz:nz)

endsubroutine POISSMG





endmodule MULTIGRID
   nx=bnx*2**level !level means from which grid we interpolate
   ny=bny*2**level
   nz=bnz*2**level

   do k=0,nz
    do j=0,ny
     do i=0,nx
                     AFine(2*i,2*j,2*k)+=ACoarse(i,j,k)
                     if (i<nx) AFine(2*i+1,2*j,2*k)+=0.5*ACoarse(i,j,k)
                     if (i>0) AFine(2*i-1,2*j,2*k)+=0.5*ACoarse(i,j,k)
                     if (j<ny) AFine(2*i,2*j+1,2*k)+=0.5*ACoarse(i,j,k)
                     if (j>0) AFine(2*i,2*j-1,2*k)+=0.5*ACoarse(i,j,k)
                     if (k<nz) AFine(2*i,2*j,2*k+1)+=0.5*ACoarse(i,j,k)
                     if (k>0) AFine(2*i,2*j,2*k-1)+=0.5*ACoarse(i,j,k)
                     if ((i<nx).and.(j<ny)) AFine(2*i+1,2*j+1,2*k)+=0.25*ACoarse(i,j,k)
                     if ((i<nx).and.(k<nz)) AFine(2*i+1,2*j,2*k+1)+=0.25*ACoarse(i,j,k)
                     if ((j<ny).and.(k<nz)) AFine(2*i,2*j+1,2*k+1)+=0.25*ACoarse(i,j,k)
                     if ((i>0).and.(j>0)) AFine(2*i-1,2*j-1,2*k)+=0.25*ACoarse(i,j,k)
                     if ((i>0).and.(k>0)) AFine(2*i-1,2*j,2*k-1)+=0.25*ACoarse(i,j,k)
                     if ((j>0).and.(k>0)) AFine(2*i,2*j-1,2*k-1)+=0.25*ACoarse(i,j,k)
                     if ((i>0).and.(j<ny)) AFine(2*i-1,2*j+1,2*k)+=0.25*ACoarse(i,j,k)
                     if ((i>0).and.(k<nz)) AFine(2*i-1,2*j,2*k+1)+=0.25*ACoarse(i,j,k)
                     if ((i<nx).and.(j>0)) AFine(2*i+1,2*j-1,2*k)+=0.25*ACoarse(i,j,k)
                     if ((j>0).and.(k<nz)) AFine(2*i,2*j-1,2*k+1)+=0.25*ACoarse(i,j,k)
                     if ((i<nx).and.(k>0)) AFine(2*i+1,2*j,2*k-1)+=0.25*ACoarse(i,j,k)
                     if ((j<ny).and.(k>0)) AFine(2*i,2*j+1,2*k-1)+=0.25*ACoarse(i,j,k)
                     if ((i<nx).and.(j<ny).and.(k<nz)) AFine(2*i+1,2*j+1,2*k+1)+=0.125*ACoarse(i,j,k)
                     if ((i>0).and.(j<ny).and.(k<nz)) AFine(2*i-1,2*j+1,2*k+1)+=0.125*ACoarse(i,j,k)
                     if ((i<nx).and.(j>0).and.(k<nz)) AFine(2*i+1,2*j-1,2*k+1)+=0.125*ACoarse(i,j,k)
                     if ((i<nx).and.(j<ny).and.(k>0)) AFine(2*i+1,2*j+1,2*k-1)+=0.125*ACoarse(i,j,k)
                     if ((i>0).and.(j>0).and.(k<nz)) AFine(2*i-1,2*j-1,2*k+1)+=0.125*ACoarse(i,j,k)
                     if ((i>0).and.(j<ny).and.(k>0)) AFine(2*i-1,2*j+1,2*k-1)+=0.125*ACoarse(i,j,k)
                     if ((i<nx).and.(j>0).and.(k>0)) AFine(2*i+1,2*j-1,2*k-1)+=0.125*ACoarse(i,j,k)
                     if ((i>0).and.(j>0).and.(k>0)) AFine(2*i-1,2*j-1,2*k-1)+=0.125*ACoarse(i,j,k)
     enddo
    enddo
   enddo
end


pure subroutine Restrict(level,ACoarse,AFine)
 integer,intent(IN):: level
 real(KND),dimension(0:,0:,0:),intent(OUT):: ACoarse
 real(KND),dimension(0:,0:,0:),intent(IN):: AFine
 real(KND) q
 integer:: i,j,k,nx,ny,nz
 
   nx=bnx*2**level !level means on which grid we restrict
   ny=bny*2**level
   nz=bnz*2**level

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
endsubroutine Restrict





subroutine MG_CGC(level,eps)
  integer k
  real(KND) R
  real(KND),INTENT(IN):: eps
  if (level == 0) then
    k = 0
    R=1 !R greater then eps**2
    do while ( ( k < 20 ) && ( R > eps*eps ) )    
      call MG_GS(level, 10)
      call get_rez_dif(level,R)
      k=k+1
    enddo
  else
   do k = 1, ncgc !number of recurrent calls
    call MG_GS(level,npre)
    call get_rez_dif(level,R)
    call grid_clear_ad(level-1)
    call Restrict(MGRHS(level-1),MGRes(level),level-1)
    call MG_CGC(level-1,eps)
    call Prolongate(MGPhi(level),MGPhi(level-1),level-1)
    call MG_GS(level,npost)
    call get_rez_dif(level,R)
   enddo
  endif
endsubroutine MG_CGC







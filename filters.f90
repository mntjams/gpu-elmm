module FILTERS

use PARAMETERS

 subroutine TOPHATFIELD(U1,U2,nx,ny,nz,width)
 real(KND),dimension(-2:,-2:,-2:) U1,U2
 integer nx,ny,nz,width,i,j,k,ii,jj,kk

 do k=1,nz
  do j=1,ny
   do i=1,nx
    S=0
    do kk=k-width,k+width
     do jj=j-width,j+width
      do kk=k-width,k+width
       S=S+U1(ii,jj,kk)
      enddo
     enddo
    enddo
    U2(i,j,k)=S/((2*width+1)**3)
   enddo
  enddo
 enddo
 endsubroutine TOPHATFIELD

endmodule FILTERS

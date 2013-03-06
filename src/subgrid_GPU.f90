

  !$hmpp <tsteps> Vreman codelet
  subroutine Vreman_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dx,dy,dz,dt,Re,U,V,W,Visc)
   !Vreman subgrid model (Physics of Fluids, 2004)

   implicit none

#include "hmpp-include.f90"

   integer,intent(in)    :: Prnx, Prny, Prnz, Unx, Uny, Unz, Vnx, Vny, Vnz, Wnx, Wny, Wnz
   real(knd),intent(in)  :: dx, dy, dz, dt, Re
   real(knd),dimension(-2:Prnx+3,-2:Prny+3,-2:Prnz+3),intent(in):: U,V,W
   real(knd),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(out):: Visc

   integer i,j,k,ii,jj
   real(knd)::bb,aa
   real(knd),dimension(1:3,1:3)::a,b
   real(knd) :: dx2,dy2,dz2
   real(knd),parameter::c=0.05

   intrinsic abs,sqrt



   dx2=dx**2
   dy2=dy**2
   dz2=dz**2
   !$hmppcg grid blocksize myblocksize2
   !$hmppcg gridify (k,i)  private(aa,bb,a,b,i,j,k,ii,jj)
   do k=0,Prnz+1
    do i=0,Prnx+1
     do j=0,Prny+1
      a(1,1)=(U(i,j,k)-U(i-1,j,k))/dx
      a(2,1)=(U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(4._knd*dy)
      a(3,1)=(U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(4._knd*dz)

      a(2,2)=(V(i,j,k)-V(i,j-1,k))/dy
      a(1,2)=(V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(4._knd*dx)
      a(3,2)=(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(4._knd*dz)

      a(3,3)=(W(i,j,k)-W(i,j,k-1))/dz
      a(1,3)=(W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(4._knd*dx)
      a(2,3)=(W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(4._knd*dy)

      !$hmppcg unroll jj:3, ii:3, noremainder
      do jj=1,3
       do ii=1,3
        b(ii,jj)=dx2*a(1,ii)*a(1,jj)+&
                       dy2*a(2,ii)*a(2,jj)+&
                       dz2*a(3,ii)*a(3,jj)
       enddo
      enddo

      bb=          b(1,1)*b(2,2)-b(1,2)**2
      bb=bb+b(1,1)*b(3,3)-b(1,3)**2
      bb=bb+b(2,2)*b(3,3)-b(2,3)**2

      aa=0



      Visc(i,j,k)=0._knd

      !$hmppcg unroll jj:3, ii:3, noremainder
      do jj=1,3
       do ii=1,3
        aa=aa+(a(ii,jj)**2)
       enddo
      enddo



      if (abs(aa)>1e-5.and.bb>0)  Visc(i,j,k)=c*sqrt(bb/aa)


      if (Re>0)  Visc(i,j,k)=Visc(i,j,k)+1._knd/Re

     enddo
    enddo
   enddo
  endsubroutine Vreman_GPU


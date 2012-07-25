  subroutine BOUND_PASSSCALAR_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,ScalBtype,sideScal,TDiff,SCAL)
  implicit none
#include "hmpp-include.f90"

  integer, intent(in) :: Prnx,Prny,Prnz,ScalBtype(6)
  real(KND),intent(in) :: dxmin,dymin,dzmin,Re,sideScal(6),TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  real(KND),intent(inout) :: SCAL(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  integer i,j,k,nx,ny,nz

   nx = Prnx
   ny = Prny
   nz = Prnz
   if (ScalBtype(We)==DIRICHLET) then
      do k = 1,nz
       do j = 1,ny
        SCAL(0,j,k) = sideScal(We)-(SCAL(1,j,k)-sideScal(We))
        SCAL(-1,j,k) = sideScal(We)-(SCAL(2,j,k)-sideScal(We))
       enddo
      enddo
   else if (ScalBtype(We)==PERIODIC) then
    do k = 1,nz
     do j = 1,ny
      SCAL(0,j,k) = SCAL(nx,j,k)
      SCAL(-1,j,k) = SCAL(nx-1,j,k)
     enddo
    enddo
   else if (ScalBtype(We)==NEUMANN) then
    do k = 1,nz
     do j = 1,ny
      SCAL(0,j,k) = SCAL(1,j,k)-sideScal(We)*dxmin
      SCAL(-1,j,k) = SCAL(1,j,k)-sideScal(We)*(dxmin+dxmin)
     enddo
    enddo
   else if (ScalBtype(We)==CONSTFLUX) then
    do k = 1,nz
     do j = 1,ny
      SCAL(0,j,k) = (SCAL(1,j,k)*((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxmin))+sideScal(We))/&
                                  ((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxmin))
      SCAL(-1,j,k) = SCAL(0,j,k)-(SCAL(1,j,k)-SCAL(0,j,k))
     enddo
    enddo
   else
    do k = 1,nz
     do j = 1,ny
      SCAL(0,j,k) = SCAL(1,j,k)-sideScal(We)*dxmin
      SCAL(-1,j,k) = SCAL(1,j,k)-sideScal(We)*(dxmin+dxmin)
     enddo
    enddo
   endif

   if (ScalBtype(Ea)==DIRICHLET) then
    do k = 1,nz
     do j = 1,ny
      SCAL(nx+1,j,k) = sideScal(Ea)-(SCAL(nx,j,k)-sideScal(Ea))
      SCAL(nx+2,j,k) = sideScal(Ea)-(SCAL(nx-1,j,k)-sideScal(Ea))
     enddo
    enddo
   else if (ScalBtype(Ea)==PERIODIC) then
    do k = 1,nz
     do j = 1,ny
      SCAL(nx+1,j,k) = SCAL(1,j,k)
      SCAL(nx+2,j,k) = SCAL(2,j,k)
     enddo
    enddo
   else if (ScalBtype(Ea)==NEUMANN) then
    do k = 1,nz
     do j = 1,ny
      SCAL(nx+1,j,k) = SCAL(nx,j,k)+sideScal(Ea)*dxmin
      SCAL(nx+2,j,k) = SCAL(nx,j,k)+sideScal(Ea)*(dxmin+dxmin)
     enddo
    enddo
   else if (ScalBtype(Ea)==CONSTFLUX) then
    do k = 1,nz
     do j = 1,ny
      SCAL(nx+1,j,k) = (SCAL(nx,j,k)*((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxmin))+sideScal(Ea))/&
       ((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxmin))
      SCAL(nx+2,j,k) = SCAL(nx+1,j,k)-(SCAL(nx,j,k)-SCAL(nx+1,j,k))
     enddo
    enddo
   else
    do k = 1,nz
     do j = 1,ny
      SCAL(nx+1,j,k) = SCAL(nx,j,k)+sideScal(Ea)*dxmin
      SCAL(nx+2,j,k) = SCAL(nx,j,k)+sideScal(Ea)*(dxmin+dxmin)
     enddo
    enddo
   endif


   if (ScalBtype(So)==DIRICHLET) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,0,k) = sideScal(So)-(SCAL(i,1,k)-sideScal(So))
      SCAL(i,-1,k) = sideScal(So)-(SCAL(i,2,k)-sideScal(So))
     enddo
    enddo
   else if (ScalBtype(So)==PERIODIC) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,0,k) = SCAL(i,ny,k)
      SCAL(i,-1,k) = SCAL(i,ny-1,k)
     enddo
    enddo
   else if (ScalBtype(So)==NEUMANN) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,0,k) = SCAL(i,1,k)-sideScal(So)*dymin
      SCAL(i,-1,k) = SCAL(i,1,k)-sideScal(So)*(dymin+dymin)
     enddo
    enddo
   else if (ScalBtype(So)==CONSTFLUX) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,0,k) = (SCAL(i,1,k)*((TDiff(i,1,k)+TDiff(i,0,k))/(2*dymin))+sideScal(So))/&
                                  ((TDiff(i,1,k)+TDiff(i,0,k))/(2*dymin))
      SCAL(i,-1,k) = SCAL(i,0,k)-(SCAL(i,1,k)-SCAL(i,0,k))
     enddo
    enddo
   else
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,0,k) = SCAL(i,1,k)-sideScal(So)*dymin
      SCAL(i,-1,k) = SCAL(i,1,k)-sideScal(So)*(dymin+dymin)
     enddo
    enddo
   endif

   if (ScalBtype(No)==DIRICHLET) then
    do k = 1,nz
     do i=-1,nx+2
       SCAL(i,ny+1,k) = sideScal(No)-(SCAL(i,ny,k)-sideScal(No))
       SCAL(i,ny+2,k) = sideScal(No)-(SCAL(i,ny-1,k)-sideScal(No))
      enddo
     enddo
   else if (ScalBtype(No)==PERIODIC) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,ny+1,k) = SCAL(i,1,k)
      SCAL(i,ny+2,k) = SCAL(i,2,k)
     enddo
    enddo
   else if (ScalBtype(No)==NEUMANN) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,ny+1,k) = SCAL(i,ny,k)+sideScal(No)*dymin
      SCAL(i,ny+2,k) = SCAL(i,ny,k)+sideScal(No)*(dymin+dymin)
     enddo
    enddo
   else if (ScalBtype(No)==CONSTFLUX) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,ny+1,k) = (SCAL(i,ny,k)*((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dymin))+sideScal(No))/&
       ((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dymin))
      SCAL(i,ny+2,k) = SCAL(i,ny+1,k)-(SCAL(i,ny,k)-SCAL(i,ny+1,k))
     enddo
    enddo
   else
    do k = 1,nz
     do i=-1,ny+2
      SCAL(i,ny+1,k) = SCAL(i,ny,k)+sideScal(No)*dymin
      SCAL(i,ny+2,k) = SCAL(i,ny,k)+sideScal(No)*(dymin+dymin)
     enddo
    enddo
   endif



    if (ScalBtype(Bo)==DIRICHLET)  then
     do j=-1,ny+2
      do i=-1,nx+2
       SCAL(i,j,0) = sideScal(Bo)-(SCAL(i,j,1)-sideScal(Bo))
       SCAL(i,j,-1) = sideScal(Bo)-(SCAL(i,j,2)-sideScal(Bo))
      enddo
     enddo
   else if (ScalBtype(Bo)==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,0) = SCAL(i,j,nz)
      SCAL(i,j,-1) = SCAL(i,j,nz-1)
     enddo
    enddo
   else if (ScalBtype(Bo)==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,0) = SCAL(i,j,1)-sideScal(Bo)*dzmin
      SCAL(i,j,-1) = SCAL(i,j,1)-sideScal(Bo)*(dzmin+dzmin)
     enddo
    enddo
   else if (ScalBtype(Bo)==CONSTFLUX.or.ScalBtype(Bo)==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,0) = SCAL(i,j,1)+sideScal(Bo)*dzmin/((TDiff(i,j,1)+TDiff(i,j,0))/(2._KND))
      SCAL(i,j,-1) = SCAL(i,j,0)-(SCAL(i,j,1)-SCAL(i,j,0))
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,0) = SCAL(i,j,1)-sideScal(Bo)*dzmin
      SCAL(i,j,-1) = SCAL(i,j,1)-sideScal(Bo)*(dzmin+dzmin)
     enddo
    enddo
   endif

   if (ScalBtype(To)==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
       SCAL(i,j,nz+1) = sideScal(To)-(SCAL(i,j,nz)-sideScal(To))
       SCAL(i,j,nz+2) = sideScal(To)-(SCAL(i,j,nz-1)-sideScal(To))
      enddo
     enddo
   else if (ScalBtype(To)==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,nz+1) = SCAL(i,j,1)
      SCAL(i,j,nz+2) = SCAL(i,j,2)
     enddo
    enddo
   else if (ScalBtype(To)==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,nz+1) = SCAL(i,j,nz)+sideScal(To)*dzmin
      SCAL(i,j,nz+2) = SCAL(i,j,nz)+sideScal(To)*(dzmin+dzmin)
     enddo
    enddo
   else if (ScalBtype(To)==CONSTFLUX) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,nz+1) = SCAL(i,j,nz)-sideScal(To)*dzmin/((TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2._KND))
      SCAL(i,j,nz+2) = SCAL(i,j,nz+1)-(SCAL(i,j,nz)-SCAL(i,j,nz+1))
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,nz+1) = SCAL(i,j,nz)+sideScal(To)*dzmin
      SCAL(i,j,nz+2) = SCAL(i,j,nz)+sideScal(To)*(dzmin+dzmin)
     enddo
    enddo
   endif
  end subroutine BOUND_PASSSCALAR_GPU






  subroutine BOUND_TEMP_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,TBtype,sideTemp,BsideTArr,BsideTFLArr,TDiff,TempIn,Theta)
  implicit none
#include "hmpp-include.f90"

  integer, intent(in) :: Prnx,Prny,Prnz,TBtype(6)
  real(KND),intent(in) :: Re,dxmin,dymin,dzmin,sideTemp(6),TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),TempIn(-1:Prny+2,-1:Prnz+2)
  real(KND),intent(in),dimension(-1:Prnx+2,-1:Prny+2) :: BsideTArr,BsideTFLArr
  real(KND),intent(inout) :: theta(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  integer i,j,k,nx,ny,nz

  intrinsic abs

   nx = Prnx
   ny = Prny
   nz = Prnz
   if (TBtype(We)==DIRICHLET) then
     do k = 1,nz
      do j = 1,ny
       theta(0,j,k) = Tempin(j,k)
       theta(-1,j,k) = Tempin(j,k)
      enddo
     enddo
   else if (TBtype(We)==PERIODIC) then
    do k = 1,nz
     do j = 1,ny
      theta(0,j,k) = theta(nx,j,k)
      theta(-1,j,k) = theta(nx-1,j,k)
     enddo
    enddo
   else if (TBtype(We)==NEUMANN) then
    do k = 1,nz
     do j = 1,ny
      theta(0,j,k) = theta(1,j,k)-sideTemp(We)*dxmin
      theta(-1,j,k) = theta(1,j,k)-sideTemp(We)*(dxmin+dxmin)
     enddo
    enddo
   else if (TBtype(We)==CONSTFLUX) then
    do k = 1,nz
     do j = 1,ny
      theta(0,j,k) = (theta(1,j,k)*((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxmin))+sideTemp(We))/&
                                  ((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxmin))
      theta(-1,j,k) = theta(0,j,k)-(theta(1,j,k)-theta(0,j,k))
     enddo
    enddo
   else
    do k = 1,nz
     do j = 1,ny
      theta(0,j,k) = theta(1,j,k)-sideTemp(We)*dxmin
      theta(-1,j,k) = theta(1,j,k)-sideTemp(We)*(dxmin+dxmin)
     enddo
    enddo
   endif

   if (TBtype(Ea)==DIRICHLET) then
    do k = 1,nz
     do j = 1,ny
      theta(nx+1,j,k) = sideTemp(Ea)
      theta(nx+2,j,k) = sideTemp(Ea)
     enddo
    enddo
   else if (TBtype(Ea)==PERIODIC) then
    do k = 1,nz
     do j = 1,ny
      theta(nx+1,j,k) = theta(1,j,k)
      theta(nx+2,j,k) = theta(2,j,k)
     enddo
    enddo
   else if (TBtype(Ea)==NEUMANN) then
    do k = 1,nz
     do j = 1,ny
      theta(nx+1,j,k) = theta(nx,j,k)+sideTemp(Ea)*dxmin
      theta(nx+2,j,k) = theta(nx,j,k)+sideTemp(Ea)*(dxmin+dxmin)
     enddo
    enddo
   else if (TBtype(Ea)==CONSTFLUX) then
    do k = 1,nz
     do j = 1,ny
      theta(nx+1,j,k) = (theta(nx,j,k)*((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxmin))+sideTemp(Ea))/&
       ((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxmin))
      theta(nx+2,j,k) = theta(nx+1,j,k)-(theta(nx,j,k)-theta(nx+1,j,k))
     enddo
    enddo
   else
    do k = 1,nz
     do j = 1,ny
      theta(nx+1,j,k) = theta(nx,j,k)+sideTemp(Ea)*dxmin
      theta(nx+2,j,k) = theta(nx,j,k)+sideTemp(Ea)*(dxmin+dxmin)
     enddo
    enddo
   endif


   if (TBtype(So)==DIRICHLET) then
    do k = 1,nz
     do i=-1,nx+2
      theta(i,0,k) = sideTemp(So)
      theta(i,-1,k) = sideTemp(So)
     enddo
    enddo
   else if (TBtype(So)==PERIODIC) then
    do k = 1,nz
     do i=-1,nx+2
      theta(i,0,k) = theta(i,ny,k)
      theta(i,-1,k) = theta(i,ny-1,k)
     enddo
    enddo
   else if (TBtype(So)==NEUMANN) then
    do k = 1,nz
     do i=-1,nx+2
      theta(i,0,k) = theta(i,1,k)-sideTemp(So)*dymin
      theta(i,-1,k) = theta(i,1,k)-sideTemp(So)*(dymin+dymin)
     enddo
    enddo
   else if (TBtype(So)==CONSTFLUX) then
    do k = 1,nz
     do i=-1,nx+2
     theta(i,0,k) = (theta(i,1,k)*((TDiff(i,1,k)+TDiff(i,0,k))/(2*dymin))+sideTemp(So))/&
                                  ((TDiff(i,1,k)+TDiff(i,0,k))/(2*dymin))
      theta(i,-1,k) = theta(i,0,k)-(theta(i,1,k)-theta(i,0,k))
     enddo
    enddo
   else
    do k = 1,nz
     do i=-1,nx+2
      theta(i,0,k) = theta(i,1,k)-sideTemp(So)*dymin
      theta(i,-1,k) = theta(i,1,k)-sideTemp(So)*(dymin+dymin)
     enddo
    enddo
   endif

   if (TBtype(No)==DIRICHLET) then
    do k = 1,nz
     do i=-1,nx+2
       theta(i,ny+1,k) = sideTemp(No)
       theta(i,ny+2,k) = sideTemp(No)
      enddo
     enddo
   else if (TBtype(No)==PERIODIC) then
    do k = 1,nz
     do i=-1,nx+2
      theta(i,ny+1,k) = theta(i,1,k)
      theta(i,ny+2,k) = theta(i,2,k)
     enddo
    enddo
   else if (TBtype(No)==NEUMANN) then
    do k = 1,nz
     do i=-1,nx+2
      theta(i,ny+1,k) = theta(i,ny,k)+sideTemp(No)*dymin
      theta(i,ny+2,k) = theta(i,ny,k)+sideTemp(No)*(dymin+dymin)
     enddo
    enddo
   else if (TBtype(No)==CONSTFLUX) then
    do k = 1,nz
     do i=-1,nx+2
      theta(i,ny+1,k) = (theta(i,ny,k)*((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dymin))+sideTemp(No))/&
       ((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dymin))
      theta(i,ny+2,k) = theta(i,ny+1,k)-(theta(i,ny,k)-theta(i,ny+1,k))
     enddo
    enddo
   else
    do k = 1,nz
     do i=-1,ny+2
      theta(i,ny+1,k) = theta(i,ny,k)+sideTemp(No)*dymin
      theta(i,ny+2,k) = theta(i,ny,k)+sideTemp(No)*(dymin+dymin)
     enddo
    enddo
   endif



   if (TBtype(Bo)==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,0) = theta(i,j,nz)
      theta(i,j,-1) = theta(i,j,nz-1)
     enddo
    enddo
   else if (TBtype(Bo)==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,0) = theta(i,j,1)-sideTemp(Bo)*dzmin
      theta(i,j,-1) = theta(i,j,1)-sideTemp(Bo)*(dzmin+dzmin)
     enddo
    enddo
   else if (TBtype(Bo)==CONSTFLUX.or.TBtype(Bo)==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
      if (abs(BsideTFLArr(i,j))<=0.and.TDiff(i,j,1)<1.5_KND/(Re).and.TBtype(Bo)==DIRICHLET) then
       theta(i,j,0) = BsideTArr(i,j)
      else
       theta(i,j,0) = theta(i,j,1)+BsideTFLArr(i,j)*dzmin/((TDiff(i,j,1)+TDiff(i,j,0))/(2._KND))
      endif
      theta(i,j,-1) = theta(i,j,0)-(theta(i,j,1)-theta(i,j,0))
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,0) = theta(i,j,1)-sideTemp(Bo)*dzmin
      theta(i,j,-1) = theta(i,j,1)-sideTemp(Bo)*(dzmin+dzmin)
     enddo
    enddo
   endif

   if (TBtype(To)==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
       theta(i,j,nz+1) = sideTemp(To)
       theta(i,j,nz+2) = sideTemp(To)
      enddo
     enddo
   else if (TBtype(To)==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,nz+1) = theta(i,j,1)
      theta(i,j,nz+2) = theta(i,j,2)
     enddo
    enddo
   else if (TBtype(To)==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,nz+1) = theta(i,j,nz)+sideTemp(To)*dzmin
      theta(i,j,nz+2) = theta(i,j,nz)+sideTemp(To)*(dzmin+dzmin)
     enddo
    enddo
   else if (TBtype(To)==CONSTFLUX) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,nz+1) = theta(i,j,nz)-sideTemp(To)*dzmin/((TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2._KND))
      theta(i,j,nz+2) = theta(i,j,nz+1)-(theta(i,j,nz)-theta(i,j,nz+1))
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,nz+1) = theta(i,j,nz)+sideTemp(To)*dzmin
      theta(i,j,nz+2) = theta(i,j,nz)+sideTemp(To)*(dzmin+dzmin)
     enddo
    enddo
   endif
  endsubroutine BOUND_TEMP_GPU


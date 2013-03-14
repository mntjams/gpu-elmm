  subroutine BoundScalar_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,ScalBtype,sideScal,TDiff,SCAL)
  implicit none
#include "hmpp-include.f90"

  integer, intent(in) :: Prnx,Prny,Prnz,ScalBtype(6)
  real(knd),intent(in) :: dxmin,dymin,dzmin,Re,sideScal(6),TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  real(knd),intent(inout) :: SCAL(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
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
      SCAL(i,j,0) = SCAL(i,j,1)+sideScal(Bo)*dzmin/((TDiff(i,j,1)+TDiff(i,j,0))/(2._knd))
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
      SCAL(i,j,nz+1) = SCAL(i,j,nz)-sideScal(To)*dzmin/((TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2._knd))
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
  end subroutine BoundScalar_GPU





!   !$hmpp <tsteps> Bound_Temperature codelet
  subroutine BoundTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,&
                                  Re,TempBtype,sideTemp,BsideTArr,BsideTFLArr,&
                                  TDiff,TempIn,Temperature)
  implicit none
#include "hmpp-include.f90"

  integer, intent(in) :: Prnx,Prny,Prnz,TempBtype(6)
  real(knd),intent(in) :: Re,dxmin,dymin,dzmin,sideTemp(6),TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),TempIn(-1:Prny+2,-1:Prnz+2)
  real(knd),intent(in),dimension(-1:Prnx+2,-1:Prny+2) :: BsideTArr,BsideTFLArr
  real(knd),intent(inout) :: Temperature(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  integer i,j,k,nx,ny,nz

  intrinsic abs

   nx = Prnx
   ny = Prny
   nz = Prnz
   if (TempBtype(We)==DIRICHLET) then
     do k = 1,nz
      do j = 1,ny
       Temperature(0,j,k) = TempIn(j,k)
       Temperature(-1,j,k) = TempIn(j,k)
      enddo
     enddo
   else if (TempBtype(We)==PERIODIC) then
    do k = 1,nz
     do j = 1,ny
      Temperature(0,j,k) = Temperature(nx,j,k)
      Temperature(-1,j,k) = Temperature(nx-1,j,k)
     enddo
    enddo
   else if (TempBtype(We)==NEUMANN) then
    do k = 1,nz
     do j = 1,ny
      Temperature(0,j,k) = Temperature(1,j,k)-sideTemp(We)*dxmin
      Temperature(-1,j,k) = Temperature(1,j,k)-sideTemp(We)*(dxmin+dxmin)
     enddo
    enddo
   else if (TempBtype(We)==CONSTFLUX) then
    do k = 1,nz
     do j = 1,ny
      Temperature(0,j,k) = (Temperature(1,j,k)*((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxmin))+sideTemp(We))/&
                                  ((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxmin))
      Temperature(-1,j,k) = Temperature(0,j,k)-(Temperature(1,j,k)-Temperature(0,j,k))
     enddo
    enddo
   else
    do k = 1,nz
     do j = 1,ny
      Temperature(0,j,k) = Temperature(1,j,k)-sideTemp(We)*dxmin
      Temperature(-1,j,k) = Temperature(1,j,k)-sideTemp(We)*(dxmin+dxmin)
     enddo
    enddo
   endif

   if (TempBtype(Ea)==DIRICHLET) then
    do k = 1,nz
     do j = 1,ny
      Temperature(nx+1,j,k) = sideTemp(Ea)
      Temperature(nx+2,j,k) = sideTemp(Ea)
     enddo
    enddo
   else if (TempBtype(Ea)==PERIODIC) then
    do k = 1,nz
     do j = 1,ny
      Temperature(nx+1,j,k) = Temperature(1,j,k)
      Temperature(nx+2,j,k) = Temperature(2,j,k)
     enddo
    enddo
   else if (TempBtype(Ea)==NEUMANN) then
    do k = 1,nz
     do j = 1,ny
      Temperature(nx+1,j,k) = Temperature(nx,j,k)+sideTemp(Ea)*dxmin
      Temperature(nx+2,j,k) = Temperature(nx,j,k)+sideTemp(Ea)*(dxmin+dxmin)
     enddo
    enddo
   else if (TempBtype(Ea)==CONSTFLUX) then
    do k = 1,nz
     do j = 1,ny
      Temperature(nx+1,j,k) = (Temperature(nx,j,k)*((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxmin))+sideTemp(Ea))/&
       ((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxmin))
      Temperature(nx+2,j,k) = Temperature(nx+1,j,k)-(Temperature(nx,j,k)-Temperature(nx+1,j,k))
     enddo
    enddo
   else
    do k = 1,nz
     do j = 1,ny
      Temperature(nx+1,j,k) = Temperature(nx,j,k)+sideTemp(Ea)*dxmin
      Temperature(nx+2,j,k) = Temperature(nx,j,k)+sideTemp(Ea)*(dxmin+dxmin)
     enddo
    enddo
   endif


   if (TempBtype(So)==DIRICHLET) then
    do k = 1,nz
     do i=-1,nx+2
      Temperature(i,0,k) = sideTemp(So)
      Temperature(i,-1,k) = sideTemp(So)
     enddo
    enddo
   else if (TempBtype(So)==PERIODIC) then
    do k = 1,nz
     do i=-1,nx+2
      Temperature(i,0,k) = Temperature(i,ny,k)
      Temperature(i,-1,k) = Temperature(i,ny-1,k)
     enddo
    enddo
   else if (TempBtype(So)==NEUMANN) then
    do k = 1,nz
     do i=-1,nx+2
      Temperature(i,0,k) = Temperature(i,1,k)-sideTemp(So)*dymin
      Temperature(i,-1,k) = Temperature(i,1,k)-sideTemp(So)*(dymin+dymin)
     enddo
    enddo
   else if (TempBtype(So)==CONSTFLUX) then
    do k = 1,nz
     do i=-1,nx+2
     Temperature(i,0,k) = (Temperature(i,1,k)*((TDiff(i,1,k)+TDiff(i,0,k))/(2*dymin))+sideTemp(So))/&
                                  ((TDiff(i,1,k)+TDiff(i,0,k))/(2*dymin))
      Temperature(i,-1,k) = Temperature(i,0,k)-(Temperature(i,1,k)-Temperature(i,0,k))
     enddo
    enddo
   else
    do k = 1,nz
     do i=-1,nx+2
      Temperature(i,0,k) = Temperature(i,1,k)-sideTemp(So)*dymin
      Temperature(i,-1,k) = Temperature(i,1,k)-sideTemp(So)*(dymin+dymin)
     enddo
    enddo
   endif

   if (TempBtype(No)==DIRICHLET) then
    do k = 1,nz
     do i=-1,nx+2
       Temperature(i,ny+1,k) = sideTemp(No)
       Temperature(i,ny+2,k) = sideTemp(No)
      enddo
     enddo
   else if (TempBtype(No)==PERIODIC) then
    do k = 1,nz
     do i=-1,nx+2
      Temperature(i,ny+1,k) = Temperature(i,1,k)
      Temperature(i,ny+2,k) = Temperature(i,2,k)
     enddo
    enddo
   else if (TempBtype(No)==NEUMANN) then
    do k = 1,nz
     do i=-1,nx+2
      Temperature(i,ny+1,k) = Temperature(i,ny,k)+sideTemp(No)*dymin
      Temperature(i,ny+2,k) = Temperature(i,ny,k)+sideTemp(No)*(dymin+dymin)
     enddo
    enddo
   else if (TempBtype(No)==CONSTFLUX) then
    do k = 1,nz
     do i=-1,nx+2
      Temperature(i,ny+1,k) = (Temperature(i,ny,k)*((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dymin))+sideTemp(No))/&
       ((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dymin))
      Temperature(i,ny+2,k) = Temperature(i,ny+1,k)-(Temperature(i,ny,k)-Temperature(i,ny+1,k))
     enddo
    enddo
   else
    do k = 1,nz
     do i=-1,ny+2
      Temperature(i,ny+1,k) = Temperature(i,ny,k)+sideTemp(No)*dymin
      Temperature(i,ny+2,k) = Temperature(i,ny,k)+sideTemp(No)*(dymin+dymin)
     enddo
    enddo
   endif



   if (TempBtype(Bo)==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      Temperature(i,j,0) = Temperature(i,j,nz)
      Temperature(i,j,-1) = Temperature(i,j,nz-1)
     enddo
    enddo
   else if (TempBtype(Bo)==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      Temperature(i,j,0) = Temperature(i,j,1)-sideTemp(Bo)*dzmin
      Temperature(i,j,-1) = Temperature(i,j,1)-sideTemp(Bo)*(dzmin+dzmin)
     enddo
    enddo
   else if (TempBtype(Bo)==CONSTFLUX.or.TempBtype(Bo)==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
      if (abs(BsideTFLArr(i,j))<=0.and.TDiff(i,j,1)<1.5_knd/(Re).and.TempBtype(Bo)==DIRICHLET) then
       Temperature(i,j,0) = BsideTArr(i,j)
      else
       Temperature(i,j,0) = Temperature(i,j,1)+BsideTFLArr(i,j)*dzmin/((TDiff(i,j,1)+TDiff(i,j,0))/(2._knd))
      endif
      Temperature(i,j,-1) = Temperature(i,j,0)-(Temperature(i,j,1)-Temperature(i,j,0))
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2
      Temperature(i,j,0) = Temperature(i,j,1)-sideTemp(Bo)*dzmin
      Temperature(i,j,-1) = Temperature(i,j,1)-sideTemp(Bo)*(dzmin+dzmin)
     enddo
    enddo
   endif

   if (TempBtype(To)==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
       Temperature(i,j,nz+1) = sideTemp(To)
       Temperature(i,j,nz+2) = sideTemp(To)
      enddo
     enddo
   else if (TempBtype(To)==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      Temperature(i,j,nz+1) = Temperature(i,j,1)
      Temperature(i,j,nz+2) = Temperature(i,j,2)
     enddo
    enddo
   else if (TempBtype(To)==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      Temperature(i,j,nz+1) = Temperature(i,j,nz)+sideTemp(To)*dzmin
      Temperature(i,j,nz+2) = Temperature(i,j,nz)+sideTemp(To)*(dzmin+dzmin)
     enddo
    enddo
   else if (TempBtype(To)==CONSTFLUX) then
    do j=-1,ny+2
     do i=-1,nx+2
      Temperature(i,j,nz+1) = Temperature(i,j,nz)-sideTemp(To)*dzmin/((TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2._knd))
      Temperature(i,j,nz+2) = Temperature(i,j,nz+1)-(Temperature(i,j,nz)-Temperature(i,j,nz+1))
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2
      Temperature(i,j,nz+1) = Temperature(i,j,nz)+sideTemp(To)*dzmin
      Temperature(i,j,nz+2) = Temperature(i,j,nz)+sideTemp(To)*(dzmin+dzmin)
     enddo
    enddo
   endif
  endsubroutine BoundTemperature_GPU


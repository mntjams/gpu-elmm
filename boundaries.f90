module BOUNDARIES
use PARAMETERS
use WALLMODELS

implicit none


  private
  public GridCoords, GridCoordsU, GridCoordsV, GridCoordsW,&
         Bound_CondU, Bound_CondV, Bound_CondW,&
         Bound_CondU2, Bound_CondV2, Bound_CondW2,&
         Bound_Phi, Bound_Pr, Bound_Q,&
         ShearInlet, ParInlet, ConstInlet


 contains


  elemental subroutine GridCoords(xi,yj,zk,x,y,z)
  integer,intent(out):: xi,yj,zk
  real(KND),intent(in):: x,y,z
  integer i

  xi=Prnx
  do i=1,Prnx
   if (xU(i)>=x) then
                 xi=i
                 exit
                endif
  enddo

  yj=Prny
  do i=1,Prny
   if (yV(i)>=y) then
                 yj=i
                 exit
                endif
  enddo
  zk=Prnz
  do i=1,Prnz
   if (zW(i)>=z) then
                 zk=i
                 exit
                endif
  enddo
  endsubroutine GridCoords

  elemental subroutine GridCoordsU(xi,yj,zk,x,y,z)
  integer,intent(out) :: xi,yj,zk
  real(KND),intent(in) :: x,y,z
  integer i

  xi=Unx+1
  do i=1,Unx+1
   if (xPr(i+1)>=x) then
                 xi=i-1
                 exit
                endif
  enddo

  yj=Vny
  do i=1,Vny
   if (yV(i)>=y) then
                 yj=i
                 exit
                endif
  enddo
  zk=Wnz
  do i=1,Wnz
   if (zW(i)>=z) then
                 zk=i
                 exit
                endif
  enddo
  end subroutine GridCoordsU


  elemental subroutine GridCoordsV(xi,yj,zk,x,y,z)
  integer,intent(out) :: xi,yj,zk
  real(KND),intent(in) :: x,y,z
  integer i

  xi=Unx
  do i=1,Unx
   if (xU(i)>=x) then
                 xi=i
                 exit
                endif
  enddo

  yj=Vny+1
  do i=1,Vny+1
   if (yPr(i)>=y) then
                 yj=i-1
                 exit
                endif
  enddo
  zk=Wnz
  do i=1,Wnz
   if (zW(i)>=z) then
                 zk=i
                 exit
                endif
  enddo
  end subroutine GridCoordsV


  elemental subroutine GridCoordsW(xi,yj,zk,x,y,z)
  integer,intent(out) :: xi,yj,zk
  real(KND),intent(in) :: x,y,z
  integer i

  xi=Unx
  do i=1,Unx
   if (xU(i)>=x) then
                 xi=i
                 exit
                endif
  enddo

  yj=Vny
  do i=1,Vny
   if (yV(i)>=y) then
                 yj=i
                 exit
                endif
  enddo

  zk=Wnz+1
  do i=1,Wnz+1
   if (zPr(i)>=z) then
                 zk=i-1
                 exit
                endif
  enddo
  end subroutine GridCoordsW

  subroutine BOUND_CONDU(U)
  real(KND),intent(inout):: U(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Unx
  ny=Uny
  nz=Unz

  !!! corners and edges for periodic conditions
  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=-2,0
      U(i,j,k)=U(i+nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j+ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=-2,0
      U(i,j,k)=U(i+nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=-2,0
      U(i,j,k)=U(i+nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=-2,0
      U(i,j,k)=U(i+nx,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif

  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC) then
   do k=1,nz
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j+ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=ny+1,ny+3
     do i=-2,0
      U(i,j,k)=U(i+nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=-2,0
      U(i,j,k)=U(i+nx,j+ny,k)
     enddo
    enddo
   enddo
   endif


  if (BtypeE==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=1,ny
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=1,ny
     do i=-2,0
      U(i,j,k)=U(i+nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=-2,0
      U(i,j,k)=U(i+nx,j,k+nz)
     enddo
    enddo
   enddo
  endif


  if (BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=1,nx
      U(i,j,k)=U(i,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=1,nx
      U(i,j,k)=U(i,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=1,nx
      U(i,j,k)=U(i,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=1,nx
      U(i,j,k)=U(i,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif


  if (BtypeW==DIRICHLET) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     U(0,j,k)=Uin(j,k)
     U(-1,j,k)=Uin(j,k)
     U(-2,j,k)=Uin(j,k)
    enddo
   enddo
  elseif (BtypeW==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     U(0,j,k)=0
     U(-1,j,k)=-U(1,j,k)
     U(-2,j,k)=-U(2,j,k)
    enddo
   enddo
  elseif (BtypeW==NEUMANN) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Neumann inlet
     U(0,j,k)=U(1,j,k)
     U(-1,j,k)=U(1,j,k)
     U(-2,j,k)=U(1,j,k)
    enddo
   enddo
  elseif (BtypeW==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     U(0,j,k)=0
     U(-1,j,k)=-U(1,j,k)
     U(-2,j,k)=-U(2,j,k)
    enddo
   enddo
  elseif (BtypeW==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     U(0,j,k)=U(nx,j,k)
     U(-1,j,k)=U(nx-1,j,k)
     U(-2,j,k)=U(nx-2,j,k)
    enddo
   enddo
  elseif (BtypeW==TURBULENTINLET) then
   do k=1,Prnz
    do j=1,Prny                       !Dirichlet inlet
     U(0,j,k)=Uin(j,k)
     U(-1,j,k)=Uin(j,k)
     U(-2,j,k)=Uin(j,k)
    enddo
   enddo
  elseif (BtypeW==INLETFROMFILE) then
   do k=1,Prnz
    do j=1,Prny                       !Dirichlet inlet
     U(0,j,k)=Uin(j,k)
     U(-1,j,k)=Uin(j,k)
     U(-2,j,k)=Uin(j,k)
    enddo
   enddo
  endif


  if (BtypeE==DIRICHLET) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     U(nx+1,j,k)=Uin(j,k)
     U(nx+2,j,k)=Uin(j,k)
     U(nx+3,j,k)=Uin(j,k)
    enddo
   enddo
  elseif (BtypeE==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     U(nx+1,j,k)=0
     U(nx+2,j,k)=-U(nx,j,k)
     U(nx+3,j,k)=-U(nx-1,j,k)
    enddo
   enddo
  elseif (BtypeE==NEUMANN.or.BtypeE==OUTLETBUFF) then   !Neumann outlet
   do k=-2,nz+3
    do j=-2,ny+3
     U(nx+1,j,k)=U(nx,j,k)
     U(nx+2,j,k)=U(nx,j,k)
     U(nx+3,j,k)=U(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     U(nx+1,j,k)=0
     U(nx+2,j,k)=-U(nx,j,k)
     U(nx+3,j,k)=-U(nx-1,j,k)
    enddo
   enddo
  elseif (BtypeE==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     U(nx+1,j,k)=U(1,j,k)
     U(nx+2,j,k)=U(2,j,k)
     U(nx+3,j,k)=U(3,j,k)
    enddo
   enddo
  elseif (BtypeE==TURBULENTINLET) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     U(nx+1,j,k)=Uin(j,k)
     U(nx+2,j,k)=Uin(j,k)
     U(nx+3,j,k)=Uin(j,k)
    enddo
   enddo
  endif

  if (BtypeS==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     U(i,0,k)=SsideU+(SsideU-U(i,1,k))
     U(i,-1,k)=SsideU+(SsideU-U(i,2,k))
     U(i,-2,k)=SsideU+(SsideU-U(i,3,k))
    enddo
   enddo
  elseif (BtypeS==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     U(i,0,k)=-U(i,1,k)
     U(i,-1,k)=-U(i,2,k)
     U(i,-2,k)=-U(i,3,k)
    enddo
   enddo
  elseif (BtypeS==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     U(i,0,k)=U(i,1,k)
     U(i,-1,k)=U(i,1,k)
     U(i,-2,k)=U(i,1,k)
    enddo
   enddo
  elseif (BtypeS==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     U(i,0,k)=U(i,1,k)
     U(i,-1,k)=U(i,1,k)
     U(i,-2,k)=U(i,1,k)
    enddo
   enddo
  elseif (BtypeS==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     U(i,0,k)=U(i,ny,k)
     U(i,-1,k)=U(i,ny-1,k)
     U(i,-2,k)=U(i,ny-2,k)
    enddo
   enddo
  endif


  if (BtypeN==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     U(i,ny+1,k)=NsideU+(NsideU-U(i,ny,k))
     U(i,ny+2,k)=NsideU+(NsideU-U(i,ny-1,k))
     U(i,ny+3,k)=NsideU+(NsideU-U(i,ny-2,k))
    enddo
   enddo
  elseif (BtypeN==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     U(i,ny+1,k)=-U(i,ny,k)
     U(i,ny+2,k)=-U(i,ny-1,k)
     U(i,ny+3,k)=-U(i,ny-2,k)
    enddo
   enddo
  elseif (BtypeN==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     U(i,ny+1,k)=U(i,ny,k)
     U(i,ny+2,k)=U(i,ny,k)
     U(i,ny+3,k)=U(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     U(i,ny+1,k)=U(i,ny,k)
     U(i,ny+2,k)=U(i,ny,k)
     U(i,ny+3,k)=U(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     U(i,ny+1,k)=U(i,1,k)
     U(i,ny+2,k)=U(i,2,k)
     U(i,ny+3,k)=U(i,3,k)
    enddo
   enddo
  endif


  if (BtypeB==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     U(i,j,0)=BsideU+(BsideU-U(i,j,1))
     U(i,j,-1)=BsideU+(BsideU-U(i,j,2))
     U(i,j,-2)=BsideU+(BsideU-U(i,j,3))
    enddo
   enddo
  elseif (BtypeB==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     U(i,j,0)=-U(i,j,1)
     U(i,j,-1)=-U(i,j,2)
     U(i,j,-2)=-U(i,j,3)
    enddo
   enddo
  elseif (BtypeB==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     U(i,j,0)=U(i,j,1)
     U(i,j,-1)=U(i,j,1)
     U(i,j,-2)=U(i,j,1)
    enddo
   enddo
  elseif (BtypeB==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     U(i,j,0)=U(i,j,1)
     U(i,j,-1)=U(i,j,1)
     U(i,j,-2)=U(i,j,1)
    enddo
   enddo
  elseif (BtypeB==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     U(i,j,0)=U(i,j,nz)
     U(i,j,-1)=U(i,j,nz-1)
     U(i,j,-2)=U(i,j,nz-2)
    enddo
   enddo
  endif

  if (BtypeT==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     U(i,j,nz+1)=TsideU+(TsideU-U(i,j,nz))
     U(i,j,nz+2)=TsideU+(TsideU-U(i,j,nz-1))
     U(i,j,nz+3)=TsideU+(TsideU-U(i,j,nz-2))
    enddo
   enddo
  elseif (BtypeT==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     U(i,j,nz+1)=-U(i,j,nz)
     U(i,j,nz+2)=-U(i,j,nz-1)
     U(i,j,nz+3)=-U(i,j,nz-2)
    enddo
   enddo
  elseif (BtypeT==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     U(i,j,nz+1)=U(i,j,nz)
     U(i,j,nz+2)=U(i,j,nz)
     U(i,j,nz+3)=U(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==FREESLIPBUFF.or.BtypeT==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     U(i,j,nz+1)=U(i,j,nz)
     U(i,j,nz+2)=U(i,j,nz)
     U(i,j,nz+3)=U(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     U(i,j,nz+1)=U(i,j,1)
     U(i,j,nz+2)=U(i,j,2)
     U(i,j,nz+3)=U(i,j,3)
    enddo
   enddo
  endif

  end subroutine BOUND_CONDU






  pure subroutine BOUND_CONDV(V)
  real(KND),intent(inout):: V(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Vnx
  ny=Vny
  nz=Vnz

  !!! corners and edges for periodic conditions
  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=-2,0
      V(i,j,k)=V(i+nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j+ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=-2,0
      V(i,j,k)=V(i+nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=-2,0
      V(i,j,k)=V(i+nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=-2,0
      V(i,j,k)=V(i+nx,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif

  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC) then
   do k=1,nz
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j+ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=ny+1,ny+3
     do i=-2,0
      V(i,j,k)=V(i+nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=-2,0
      V(i,j,k)=V(i+nx,j+ny,k)
     enddo
    enddo
   enddo
   endif

  if (BtypeE==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=1,ny
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=1,ny
     do i=-2,0
      V(i,j,k)=V(i+nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=-2,0
      V(i,j,k)=V(i+nx,j,k+nz)
     enddo
    enddo
   enddo
  endif

  if (BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=1,nx
      V(i,j,k)=V(i,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=1,nx
      V(i,j,k)=V(i,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=1,nx
      V(i,j,k)=V(i,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=1,nx
      V(i,j,k)=V(i,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif


  if (BtypeW==DIRICHLET) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     V(0,j,k)=0
     V(-1,j,k)=0
     V(-2,j,k)=0
    enddo
   enddo
  elseif (BtypeW==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     V(0,j,k)=-V(1,j,k)
     V(-1,j,k)=-V(2,j,k)
     V(-2,j,k)=-V(3,j,k)
    enddo
   enddo
  elseif (BtypeW==NEUMANN) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Neumann inlet
     V(0,j,k)=V(1,j,k)
     V(-1,j,k)=V(1,j,k)
     V(-2,j,k)=V(1,j,k)
    enddo
   enddo
  elseif (BtypeW==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     V(0,j,k)=V(1,j,k)
     V(-1,j,k)=V(1,j,k)
     V(-2,j,k)=V(1,j,k)
    enddo
   enddo
  elseif (BtypeW==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     V(0,j,k)=V(nx,j,k)
     V(-1,j,k)=V(nx-1,j,k)
     V(-2,j,k)=V(nx-2,j,k)
    enddo
   enddo
  elseif (BtypeW==TURBULENTINLET) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     V(0,j,k)=Vin(j,k)
     V(-1,j,k)=Vin(j,k)
     V(-2,j,k)=Vin(j,k)
    enddo
   enddo
  elseif (BtypeW==INLETFROMFILE) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     V(0,j,k)=Vin(j,k)
     V(-1,j,k)=Vin(j,k)
     V(-2,j,k)=Vin(j,k)
    enddo
   enddo
  endif

  if (BtypeE==DIRICHLET) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     V(nx+1,j,k)=0
     V(nx+2,j,k)=0
     V(nx+3,j,k)=0
    enddo
   enddo
  elseif (BtypeE==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     V(nx+1,j,k)=-V(nx,j,k)
     V(nx+2,j,k)=-V(nx-1,j,k)
     V(nx+3,j,k)=-V(nx-2,j,k)
    enddo
   enddo
  elseif (BtypeE==NEUMANN.or.BtypeE==OUTLETBUFF) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Neumann inlet
     V(nx+1,j,k)=V(nx,j,k)
     V(nx+2,j,k)=V(nx,j,k)
     V(nx+3,j,k)=V(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     V(nx+1,j,k)=V(nx,j,k)
     V(nx+2,j,k)=V(nx,j,k)
     V(nx+3,j,k)=V(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     V(nx+1,j,k)=V(1,j,k)
     V(nx+2,j,k)=V(2,j,k)
     V(nx+3,j,k)=V(3,j,k)
    enddo
   enddo
  elseif (BtypeE==TURBULENTINLET) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     V(nx+1,j,k)=Vin(j,k)
     V(nx+2,j,k)=Vin(j,k)
     V(nx+3,j,k)=Vin(j,k)
    enddo
   enddo
  endif

  if (BtypeS==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     V(i,0,k)=SsideV
     V(i,-1,k)=SsideV+(SsideV-V(i,1,k))
     V(i,-2,k)=SsideV+(SsideV-V(i,2,k))
    enddo
   enddo
  elseif (BtypeS==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     V(i,0,k)=0
     V(i,-1,k)=-V(i,1,k)
     V(i,-2,k)=-V(i,2,k)
    enddo
   enddo
  elseif (BtypeS==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     V(i,0,k)=V(i,1,k)
     V(i,-1,k)=V(i,1,k)
     V(i,-2,k)=V(i,1,k)
    enddo
   enddo
  elseif (BtypeS==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     V(i,0,k)=0
     V(i,-1,k)=-V(i,1,k)
     V(i,-2,k)=-V(i,2,k)
    enddo
   enddo
  elseif (BtypeS==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     V(i,0,k)=V(i,ny,k)
     V(i,-1,k)=V(i,ny-1,k)
     V(i,-2,k)=V(i,ny-2,k)
    enddo
   enddo
  endif

  if (BtypeN==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     V(i,ny+1,k)=NsideV
     V(i,ny+2,k)=NsideV+(NsideV-V(i,ny,k))
     V(i,ny+3,k)=NsideV+(NsideV-V(i,ny-1,k))
    enddo
   enddo
  elseif (BtypeN==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     V(i,ny+1,k)=0
     V(i,ny+2,k)=-V(i,ny,k)
     V(i,ny+3,k)=-V(i,ny-1,k)
    enddo
   enddo
  elseif (BtypeN==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     V(i,ny+1,k)=V(i,ny,k)
     V(i,ny+2,k)=V(i,ny,k)
     V(i,ny+3,k)=V(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     V(i,ny+1,k)=0
     V(i,ny+2,k)=-V(i,ny,k)
     V(i,ny+3,k)=-V(i,ny-1,k)
    enddo
   enddo
  elseif (BtypeN==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     V(i,ny+1,k)=V(i,1,k)
     V(i,ny+2,k)=V(i,2,k)
     V(i,ny+3,k)=V(i,3,k)
    enddo
   enddo
  endif

  if (BtypeB==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     V(i,j,0)=BsideV+(BsideV-V(i,j,1))
     V(i,j,-1)=BsideV+(BsideV-V(i,j,2))
     V(i,j,-2)=BsideV+(BsideV-V(i,j,3))
    enddo
   enddo
  elseif (BtypeB==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     V(i,j,0)=-V(i,j,1)
     V(i,j,-1)=-V(i,j,2)
     V(i,j,-2)=-V(i,j,3)
    enddo
   enddo
  elseif (BtypeB==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     V(i,j,0)=V(i,j,1)
     V(i,j,-1)=V(i,j,1)
     V(i,j,-2)=V(i,j,1)
    enddo
   enddo
  elseif (BtypeB==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     V(i,j,0)=V(i,j,1)
     V(i,j,-1)=V(i,j,1)
     V(i,j,-2)=V(i,j,1)
    enddo
   enddo
  elseif (BtypeB==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     V(i,j,0)=V(i,j,nz)
     V(i,j,-1)=V(i,j,nz-1)
     V(i,j,-2)=V(i,j,nz-2)
    enddo
   enddo
  endif

  if (BtypeT==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     V(i,j,nz+1)=TsideV+(TsideV-V(i,j,nz))
     V(i,j,nz+2)=TsideV+(TsideV-V(i,j,nz-1))
     V(i,j,nz+3)=TsideV+(TsideV-V(i,j,nz-2))
    enddo
   enddo
  elseif (BtypeT==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     V(i,j,nz+1)=-V(i,j,nz)
     V(i,j,nz+2)=-V(i,j,nz-1)
     V(i,j,nz+3)=-V(i,j,nz-2)
    enddo
   enddo
  elseif (BtypeT==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     V(i,j,nz+1)=V(i,j,nz)
     V(i,j,nz+2)=V(i,j,nz)
     V(i,j,nz+3)=V(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==FREESLIPBUFF.or.BtypeT==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     V(i,j,nz+1)=V(i,j,nz)
     V(i,j,nz+2)=V(i,j,nz)
     V(i,j,nz+3)=V(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     V(i,j,nz+1)=V(i,j,1)
     V(i,j,nz+2)=V(i,j,2)
     V(i,j,nz+3)=V(i,j,3)
    enddo
   enddo
  endif
  end subroutine BOUND_CONDV







  pure subroutine BOUND_CONDW(W)
  real(KND),intent(inout):: W(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Wnx
  ny=Wny
  nz=Wnz

  !!! corners and edges for periodic conditions
  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=-2,0
      W(i,j,k)=W(i+nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j+ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=-2,0
      W(i,j,k)=W(i+nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=-2,0
      W(i,j,k)=W(i+nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=-2,0
      W(i,j,k)=W(i+nx,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif

  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC) then
   do k=1,nz
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j+ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=ny+1,ny+3
     do i=-2,0
      W(i,j,k)=W(i+nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=-2,0
      W(i,j,k)=W(i+nx,j+ny,k)
     enddo
    enddo
   enddo
   endif

  if (BtypeE==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=1,ny
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=1,ny
     do i=-2,0
      W(i,j,k)=W(i+nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=-2,0
      W(i,j,k)=W(i+nx,j,k+nz)
     enddo
    enddo
   enddo
  endif

  if (BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=1,nx
      W(i,j,k)=W(i,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=1,nx
      W(i,j,k)=W(i,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=1,nx
      W(i,j,k)=W(i,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=1,nx
      W(i,j,k)=W(i,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif


  if (BtypeW==DIRICHLET) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     W(0,j,k)=0
     W(-1,j,k)=0
     W(-2,j,k)=0
    enddo
   enddo
  elseif (BtypeW==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     W(0,j,k)=-W(1,j,k)
     W(-1,j,k)=-W(2,j,k)
     W(-2,j,k)=-W(3,j,k)
    enddo
   enddo
  elseif (BtypeW==NEUMANN) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Neumann inlet
     W(0,j,k)=W(1,j,k)
     W(-1,j,k)=W(1,j,k)
     W(-2,j,k)=W(1,j,k)
    enddo
   enddo
  elseif (BtypeW==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     W(0,j,k)=W(1,j,k)
     W(-1,j,k)=W(1,j,k)
     W(-2,j,k)=W(1,j,k)
    enddo
   enddo
  elseif (BtypeW==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     W(0,j,k)=W(nx,j,k)
     W(-1,j,k)=W(nx-1,j,k)
     W(-2,j,k)=W(nx-2,j,k)
    enddo
   enddo
  elseif (BtypeW==TURBULENTINLET) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     W(0,j,k)=Win(j,k)
     W(-1,j,k)=Win(j,k)
     W(-2,j,k)=Win(j,k)
    enddo
   enddo
  elseif (BtypeW==INLETFROMFILE) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     W(0,j,k)=Win(j,k)
     W(-1,j,k)=Win(j,k)
     W(-2,j,k)=Win(j,k)
    enddo
   enddo
  endif

  if (BtypeE==DIRICHLET) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     W(nx+1,j,k)=0
     W(nx+2,j,k)=0
     W(nx+3,j,k)=0
    enddo
   enddo
  elseif (BtypeE==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     W(nx+1,j,k)=-W(nx,j,k)
     W(nx+2,j,k)=-W(nx-1,j,k)
     W(nx+3,j,k)=-W(nx-2,j,k)
    enddo
   enddo
  elseif (BtypeE==NEUMANN.or.BtypeE==OUTLETBUFF) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Neumann inlet
     W(nx+1,j,k)=W(nx,j,k)
     W(nx+2,j,k)=W(nx,j,k)
     W(nx+3,j,k)=W(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     W(nx+1,j,k)=W(nx,j,k)
     W(nx+2,j,k)=W(nx,j,k)
     W(nx+3,j,k)=W(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     W(nx+1,j,k)=W(1,j,k)
     W(nx+2,j,k)=W(2,j,k)
     W(nx+3,j,k)=W(3,j,k)
    enddo
   enddo
  elseif (BtypeE==TURBULENTINLET) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     W(nx+1,j,k)=Win(j,k)
     W(nx+2,j,k)=Win(j,k)
     W(nx+3,j,k)=Win(j,k)
    enddo
   enddo
  endif

  if (BtypeS==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     W(i,0,k)=SsideW+(SsideW-W(i,1,k))
     W(i,-1,k)=SsideW+(SsideW-W(i,2,k))
     W(i,-2,k)=SsideW+(SsideW-W(i,3,k))
    enddo
   enddo
  elseif (BtypeS==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     W(i,0,k)=-W(i,1,k)
     W(i,-1,k)=-W(i,2,k)
     W(i,-2,k)=-W(i,3,k)
    enddo
   enddo
  elseif (BtypeS==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     W(i,0,k)=W(i,1,k)
     W(i,-1,k)=W(i,1,k)
     W(i,-2,k)=W(i,1,k)
    enddo
   enddo
  elseif (BtypeS==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     W(i,0,k)=W(i,1,k)
     W(i,-1,k)=W(i,1,k)
     W(i,-2,k)=W(i,1,k)
    enddo
   enddo
  elseif (BtypeS==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     W(i,0,k)=W(i,ny,k)
     W(i,-1,k)=W(i,ny-1,k)
     W(i,-2,k)=W(i,ny-2,k)
    enddo
   enddo
  endif

  if (BtypeN==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     W(i,ny+1,k)=NsideW+(NsideW-W(i,ny,k))
     W(i,ny+2,k)=NsideW+(NsideW-W(i,ny-1,k))
     W(i,ny+3,k)=NsideW+(NsideW-W(i,ny-2,k))
    enddo
   enddo
  elseif (BtypeN==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     W(i,ny+1,k)=-W(i,ny,k)
     W(i,ny+2,k)=-W(i,ny-1,k)
     W(i,ny+3,k)=-W(i,ny-2,k)
    enddo
   enddo
  elseif (BtypeN==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     W(i,ny+1,k)=W(i,ny,k)
     W(i,ny+2,k)=W(i,ny,k)
     W(i,ny+3,k)=W(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     W(i,ny+1,k)=W(i,ny,k)
     W(i,ny+2,k)=W(i,ny,k)
     W(i,ny+3,k)=W(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     W(i,ny+1,k)=W(i,1,k)
     W(i,ny+2,k)=W(i,2,k)
     W(i,ny+3,k)=W(i,3,k)
    enddo
   enddo
  endif

  if (BtypeB==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     W(i,j,0)=BsideW
     W(i,j,-1)=BsideW+(BsideW-W(i,j,1))
     W(i,j,-2)=BsideW+(BsideW-W(i,j,2))
    enddo
   enddo
  elseif (BtypeB==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     W(i,j,0)=0
     W(i,j,-1)=-W(i,j,1)
     W(i,j,-2)=-W(i,j,2)
    enddo
   enddo
  elseif (BtypeB==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     W(i,j,0)=W(i,j,1)
     W(i,j,-1)=W(i,j,1)
     W(i,j,-2)=W(i,j,1)
    enddo
   enddo
  elseif (BtypeB==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     W(i,j,0)=0
     W(i,j,-1)=-W(i,j,1)
     W(i,j,-2)=-W(i,j,2)
    enddo
   enddo
  elseif (BtypeB==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     W(i,j,0)=W(i,j,nz)
     W(i,j,-1)=W(i,j,nz-1)
     W(i,j,-2)=W(i,j,nz-2)
    enddo
   enddo
  endif

  if (BtypeT==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     W(i,j,nz+1)=TsideW
     W(i,j,nz+2)=TsideW+(TsideW-W(i,j,nz))
     W(i,j,nz+3)=TsideW+(TsideW-W(i,j,nz-1))
    enddo
   enddo
  elseif (BtypeT==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     W(i,j,nz+1)=0
     W(i,j,nz+2)=-W(i,j,nz)
     W(i,j,nz+3)=-W(i,j,nz-1)
    enddo
   enddo
  elseif (BtypeT==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     W(i,j,nz+1)=W(i,j,nz)
     W(i,j,nz+2)=W(i,j,nz)
     W(i,j,nz+3)=W(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==FREESLIPBUFF.or.BtypeT==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     W(i,j,nz+1)=0
     W(i,j,nz+2)=-W(i,j,nz)
     W(i,j,nz+3)=-W(i,j,nz-1)
    enddo
   enddo
  elseif (BtypeT==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     W(i,j,nz+1)=W(i,j,1)
     W(i,j,nz+2)=W(i,j,2)
     W(i,j,nz+3)=W(i,j,3)
    enddo
   enddo
  endif
  end subroutine BOUND_CONDW





  pure subroutine BOUND_CONDU2(U)
  real(KND),intent(inout):: U(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Unx
  ny=Uny
  nz=Unz

  !!! corners and edges for periodic conditions
  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=-2,0
      U(i,j,k)=U(i+nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j+ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=-2,0
      U(i,j,k)=U(i+nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=-2,0
      U(i,j,k)=U(i+nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=-2,0
      U(i,j,k)=U(i+nx,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif

  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC) then
   do k=1,nz
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j+ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=ny+1,ny+3
     do i=-2,0
      U(i,j,k)=U(i+nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=-2,0
      U(i,j,k)=U(i+nx,j+ny,k)
     enddo
    enddo
   enddo
   endif

  if (BtypeE==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=1,ny
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=1,ny
     do i=-2,0
      U(i,j,k)=U(i+nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=-2,0
      U(i,j,k)=U(i+nx,j,k+nz)
     enddo
    enddo
   enddo
  endif

  if (BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=1,nx
      U(i,j,k)=U(i,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=1,nx
      U(i,j,k)=U(i,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=1,nx
      U(i,j,k)=U(i,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=1,nx
      U(i,j,k)=U(i,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif


  if (BtypeS==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     U(i,0,k)=-U(i,1,k)
     U(i,-1,k)=-U(i,2,k)
     U(i,-2,k)=-U(i,3,k)
    enddo
   enddo
  elseif (BtypeS==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     U(i,0,k)=-U(i,1,k)
     U(i,-1,k)=-U(i,2,k)
     U(i,-2,k)=-U(i,3,k)
    enddo
   enddo
  elseif (BtypeS==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     U(i,0,k)=U(i,1,k)
     U(i,-1,k)=U(i,1,k)
     U(i,-2,k)=U(i,1,k)
    enddo
   enddo
  elseif (BtypeS==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     U(i,0,k)=U(i,1,k)
     U(i,-1,k)=U(i,1,k)
     U(i,-2,k)=U(i,1,k)
    enddo
   enddo
  elseif (BtypeS==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     U(i,0,k)=U(i,ny,k)
     U(i,-1,k)=U(i,ny-1,k)
     U(i,-2,k)=U(i,ny-2,k)
    enddo
   enddo
  endif

  if (BtypeN==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     U(i,ny+1,k)=-U(i,ny,k)
     U(i,ny+2,k)=-U(i,ny-1,k)
     U(i,ny+3,k)=-U(i,ny-2,k)
    enddo
   enddo
  elseif (BtypeN==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     U(i,ny+1,k)=-U(i,ny,k)
     U(i,ny+2,k)=-U(i,ny-1,k)
     U(i,ny+3,k)=-U(i,ny-2,k)
    enddo
   enddo
  elseif (BtypeN==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     U(i,ny+1,k)=U(i,ny,k)
     U(i,ny+2,k)=U(i,ny,k)
     U(i,ny+3,k)=U(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     U(i,ny+1,k)=U(i,ny,k)
     U(i,ny+2,k)=U(i,ny,k)
     U(i,ny+3,k)=U(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     U(i,ny+1,k)=U(i,1,k)
     U(i,ny+2,k)=U(i,2,k)
     U(i,ny+3,k)=U(i,3,k)
    enddo
   enddo
  endif

  if (BtypeB==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     U(i,j,0)=-U(i,j,1)
     U(i,j,-1)=-U(i,j,2)
     U(i,j,-2)=-U(i,j,3)
    enddo
   enddo
  elseif (BtypeB==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     U(i,j,0)=-U(i,j,1)
     U(i,j,-1)=-U(i,j,2)
     U(i,j,-2)=-U(i,j,3)
    enddo
   enddo
  elseif (BtypeB==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     U(i,j,0)=U(i,j,1)
     U(i,j,-1)=U(i,j,1)
     U(i,j,-2)=U(i,j,1)
    enddo
   enddo
  elseif (BtypeB==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     U(i,j,0)=U(i,j,1)
     U(i,j,-1)=U(i,j,1)
     U(i,j,-2)=U(i,j,1)
    enddo
   enddo
  elseif (BtypeB==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     U(i,j,0)=U(i,j,nz)
     U(i,j,-1)=U(i,j,nz-1)
     U(i,j,-2)=U(i,j,nz-2)
    enddo
   enddo
  endif

  if (BtypeT==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     U(i,j,nz+1)=-U(i,j,nz)
     U(i,j,nz+2)=-U(i,j,nz-1)
     U(i,j,nz+3)=-U(i,j,nz-2)
    enddo
   enddo
  elseif (BtypeT==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     U(i,j,nz+1)=-U(i,j,nz)
     U(i,j,nz+2)=-U(i,j,nz-1)
     U(i,j,nz+3)=-U(i,j,nz-2)
    enddo
   enddo
  elseif (BtypeT==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     U(i,j,nz+1)=U(i,j,nz)
     U(i,j,nz+2)=U(i,j,nz)
     U(i,j,nz+3)=U(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==FREESLIPBUFF.or.BtypeT==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     U(i,j,nz+1)=U(i,j,nz)
     U(i,j,nz+2)=U(i,j,nz)
     U(i,j,nz+3)=U(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     U(i,j,nz+1)=U(i,j,1)
     U(i,j,nz+2)=U(i,j,2)
     U(i,j,nz+3)=U(i,j,3)
    enddo
   enddo
  endif

  if ((BtypeW==DIRICHLET).or.(BtypeW==TURBULENTINLET).or.(BtypeW==INLETFROMFILE)) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     U(0,j,k)=0
     U(-1,j,k)=0
     U(-2,j,k)=0
    enddo
   enddo
  elseif (BtypeW==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     U(0,j,k)=0
     U(-1,j,k)=-U(1,j,k)
     U(-2,j,k)=-U(2,j,k)
    enddo
   enddo
  elseif (BtypeW==NEUMANN) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Neumann inlet
     U(0,j,k)=U(1,j,k)
     U(-1,j,k)=U(1,j,k)
     U(-2,j,k)=U(1,j,k)
    enddo
   enddo
  elseif (BtypeW==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     U(0,j,k)=0
     U(-1,j,k)=-U(1,j,k)
     U(-2,j,k)=-U(2,j,k)
    enddo
   enddo
  elseif (BtypeW==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     U(0,j,k)=U(nx,j,k)
     U(-1,j,k)=U(nx-1,j,k)
     U(-2,j,k)=U(nx-2,j,k)
    enddo
   enddo
  endif

  if ((BtypeE==DIRICHLET).or.(BtypeE==TURBULENTINLET)) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     U(nx+1,j,k)=0
     U(nx+2,j,k)=0
     U(nx+3,j,k)=0
    enddo
   enddo
  elseif (BtypeE==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     U(nx+1,j,k)=0
     U(nx+2,j,k)=-U(nx,j,k)
     U(nx+3,j,k)=-U(nx-1,j,k)
    enddo
   enddo
  elseif (BtypeE==NEUMANN.or.BtypeE==OUTLETBUFF) then   !Neumann outlet

   do k=-2,nz+3
    do j=-2,ny+3
     U(nx+1,j,k)=U(nx,j,k)
     U(nx+2,j,k)=U(nx,j,k)
     U(nx+3,j,k)=U(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     U(nx+1,j,k)=0
     U(nx+2,j,k)=-U(nx,j,k)
     U(nx+3,j,k)=-U(nx-1,j,k)
    enddo
   enddo
  elseif (BtypeE==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     U(nx+1,j,k)=U(1,j,k)
     U(nx+2,j,k)=U(2,j,k)
     U(nx+3,j,k)=U(3,j,k)
    enddo
   enddo
  endif

  end subroutine BOUND_CONDU2


  pure subroutine BOUND_CONDV2(V)
  real(KND),intent(inout):: V(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Vnx
  ny=Vny
  nz=Vnz

  !!! corners and edges for periodic conditions
  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=-2,0
      V(i,j,k)=V(i+nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j+ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=-2,0
      V(i,j,k)=V(i+nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=-2,0
      V(i,j,k)=V(i+nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=-2,0
      V(i,j,k)=V(i+nx,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif

  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC) then
   do k=1,nz
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j+ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=ny+1,ny+3
     do i=-2,0
      V(i,j,k)=V(i+nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=-2,0
      V(i,j,k)=V(i+nx,j+ny,k)
     enddo
    enddo
   enddo
   endif

  if (BtypeE==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=1,ny
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=nx+1,nx+3
      V(i,j,k)=V(i-nx,j,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=1,ny
     do i=-2,0
      V(i,j,k)=V(i+nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=-2,0
      V(i,j,k)=V(i+nx,j,k+nz)
     enddo
    enddo
   enddo
  endif

  if (BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=1,nx
      V(i,j,k)=V(i,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=1,nx
      V(i,j,k)=V(i,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=1,nx
      V(i,j,k)=V(i,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=1,nx
      V(i,j,k)=V(i,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif


  if ((BtypeW==DIRICHLET).or.(BtypeW==TURBULENTINLET).or.(BtypeW==INLETFROMFILE)) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     V(0,j,k)=0
     V(-1,j,k)=0
     V(-2,j,k)=0
    enddo
   enddo
  elseif (BtypeW==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     V(0,j,k)=-V(1,j,k)
     V(-1,j,k)=-V(2,j,k)
     V(-2,j,k)=-V(3,j,k)
    enddo
   enddo
  elseif (BtypeW==NEUMANN) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Neumann inlet
     V(0,j,k)=V(1,j,k)
     V(-1,j,k)=V(1,j,k)
     V(-2,j,k)=V(1,j,k)
    enddo
   enddo
  elseif (BtypeW==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     V(0,j,k)=V(1,j,k)
     V(-1,j,k)=V(1,j,k)
     V(-2,j,k)=V(1,j,k)
    enddo
   enddo
  elseif (BtypeW==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     V(0,j,k)=V(nx,j,k)
     V(-1,j,k)=V(nx-1,j,k)
     V(-2,j,k)=V(nx-2,j,k)
    enddo
   enddo
  endif

  if ((BtypeE==DIRICHLET).or.(BtypeE==TURBULENTINLET)) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     V(nx+1,j,k)=0
     V(nx+2,j,k)=0
     V(nx+3,j,k)=0
    enddo
   enddo
  elseif (BtypeE==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     V(nx+1,j,k)=-V(nx,j,k)
     V(nx+2,j,k)=-V(nx-1,j,k)
     V(nx+3,j,k)=-V(nx-2,j,k)
    enddo
   enddo
  elseif (BtypeE==NEUMANN.or.BtypeE==OUTLETBUFF) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Neumann inlet
     V(nx+1,j,k)=V(nx,j,k)
     V(nx+2,j,k)=V(nx,j,k)
     V(nx+3,j,k)=V(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     V(nx+1,j,k)=V(nx,j,k)
     V(nx+2,j,k)=V(nx,j,k)
     V(nx+3,j,k)=V(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     V(nx+1,j,k)=V(1,j,k)
     V(nx+2,j,k)=V(2,j,k)
     V(nx+3,j,k)=V(3,j,k)
    enddo
   enddo
  endif

  if (BtypeS==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     V(i,0,k)=0
     V(i,-1,k)=-V(i,1,k)
     V(i,-2,k)=-V(i,2,k)
    enddo
   enddo
  elseif (BtypeS==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     V(i,0,k)=0
     V(i,-1,k)=-V(i,1,k)
     V(i,-2,k)=-V(i,2,k)
    enddo
   enddo
  elseif (BtypeS==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     V(i,0,k)=V(i,1,k)
     V(i,-1,k)=V(i,1,k)
     V(i,-2,k)=V(i,1,k)
    enddo
   enddo
  elseif (BtypeS==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     V(i,0,k)=0
     V(i,-1,k)=-V(i,1,k)
     V(i,-2,k)=-V(i,2,k)
    enddo
   enddo
  elseif (BtypeS==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     V(i,0,k)=V(i,ny,k)
     V(i,-1,k)=V(i,ny-1,k)
     V(i,-2,k)=V(i,ny-2,k)
    enddo
   enddo
  endif

  if (BtypeN==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     V(i,ny+1,k)=0
     V(i,ny+2,k)=-V(i,ny,k)
     V(i,ny+3,k)=-V(i,ny-1,k)
    enddo
   enddo
  elseif (BtypeN==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     V(i,ny+1,k)=0
     V(i,ny+2,k)=-V(i,ny,k)
     V(i,ny+3,k)=-V(i,ny-1,k)
    enddo
   enddo
  elseif (BtypeN==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     V(i,ny+1,k)=V(i,ny,k)
     V(i,ny+2,k)=V(i,ny,k)
     V(i,ny+3,k)=V(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     V(i,ny+1,k)=0
     V(i,ny+2,k)=-V(i,ny,k)
     V(i,ny+3,k)=-V(i,ny-1,k)
    enddo
   enddo
  elseif (BtypeN==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     V(i,ny+1,k)=V(i,1,k)
     V(i,ny+2,k)=V(i,2,k)
     V(i,ny+3,k)=V(i,3,k)
    enddo
   enddo
  endif

  if (BtypeB==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     V(i,j,0)=-V(i,j,1)
     V(i,j,-1)=-V(i,j,2)
     V(i,j,-2)=-V(i,j,3)
    enddo
   enddo
  elseif (BtypeB==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     V(i,j,0)=-V(i,j,1)
     V(i,j,-1)=-V(i,j,2)
     V(i,j,-2)=-V(i,j,3)
    enddo
   enddo
  elseif (BtypeB==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     V(i,j,0)=V(i,j,1)
     V(i,j,-1)=V(i,j,1)
     V(i,j,-2)=V(i,j,1)
    enddo
   enddo
  elseif (BtypeB==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     V(i,j,0)=V(i,j,1)
     V(i,j,-1)=V(i,j,1)
     V(i,j,-2)=V(i,j,1)
    enddo
   enddo
  elseif (BtypeB==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     V(i,j,0)=V(i,j,nz)
     V(i,j,-1)=V(i,j,nz-1)
     V(i,j,-2)=V(i,j,nz-2)
    enddo
   enddo
  endif

  if (BtypeT==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     V(i,j,nz+1)=-V(i,j,nz)
     V(i,j,nz+2)=-V(i,j,nz-1)
     V(i,j,nz+3)=-V(i,j,nz-2)
    enddo
   enddo
  elseif (BtypeT==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     V(i,j,nz+1)=-V(i,j,nz)
     V(i,j,nz+2)=-V(i,j,nz-1)
     V(i,j,nz+3)=-V(i,j,nz-2)
    enddo
   enddo
  elseif (BtypeT==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     V(i,j,nz+1)=V(i,j,nz)
     V(i,j,nz+2)=V(i,j,nz)
     V(i,j,nz+3)=V(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==FREESLIPBUFF.or.BtypeT==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     V(i,j,nz+1)=V(i,j,nz)
     V(i,j,nz+2)=V(i,j,nz)
     V(i,j,nz+3)=V(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     V(i,j,nz+1)=V(i,j,1)
     V(i,j,nz+2)=V(i,j,2)
     V(i,j,nz+3)=V(i,j,3)
    enddo
   enddo
  endif

  end subroutine BOUND_CONDV2


  pure subroutine BOUND_CONDW2(W)
  real(KND),intent(inout):: W(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Wnx
  ny=Wny
  nz=Wnz

  !!! corners and edges for periodic conditions
  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=-2,0
      W(i,j,k)=W(i+nx,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j+ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=-2,0
      W(i,j,k)=W(i+nx,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=-2,0
      W(i,j,k)=W(i+nx,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=-2,0
      W(i,j,k)=W(i+nx,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif

  if (BtypeE==PERIODIC.and.BtypeN==PERIODIC) then
   do k=1,nz
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j+ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=ny+1,ny+3
     do i=-2,0
      W(i,j,k)=W(i+nx,j-ny,k)
     enddo
    enddo
   enddo
   do k=1,nz
    do j=-2,0
     do i=-2,0
      W(i,j,k)=W(i+nx,j+ny,k)
     enddo
    enddo
   enddo
   endif

  if (BtypeE==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=1,ny
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=nx+1,nx+3
      W(i,j,k)=W(i-nx,j,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=1,ny
     do i=-2,0
      W(i,j,k)=W(i+nx,j,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=1,ny
     do i=-2,0
      W(i,j,k)=W(i+nx,j,k+nz)
     enddo
    enddo
   enddo
  endif

  if (BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=1,nx
      W(i,j,k)=W(i,j-ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=ny+1,ny+3
     do i=1,nx
      W(i,j,k)=W(i,j-ny,k+nz)
     enddo
    enddo
   enddo
   do k=nz+1,nz+3
    do j=-2,0
     do i=1,nx
      W(i,j,k)=W(i,j+ny,k-nz)
     enddo
    enddo
   enddo
   do k=-2,0
    do j=-2,0
     do i=1,nx
      W(i,j,k)=W(i,j+ny,k+nz)
     enddo
    enddo
   enddo
  endif


  if ((BtypeW==DIRICHLET).or.(BtypeW==TURBULENTINLET).or.(BtypeW==INLETFROMFILE)) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     W(0,j,k)=0
     W(-1,j,k)=0
     W(-2,j,k)=0
    enddo
   enddo
  elseif (BtypeW==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     W(0,j,k)=-W(1,j,k)
     W(-1,j,k)=-W(2,j,k)
     W(-2,j,k)=-W(3,j,k)
    enddo
   enddo
  elseif (BtypeW==NEUMANN) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Neumann inlet
     W(0,j,k)=W(1,j,k)
     W(-1,j,k)=W(1,j,k)
     W(-2,j,k)=W(1,j,k)
    enddo
   enddo
  elseif (BtypeW==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     W(0,j,k)=W(1,j,k)
     W(-1,j,k)=W(1,j,k)
     W(-2,j,k)=W(1,j,k)
    enddo
   enddo
  elseif (BtypeW==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     W(0,j,k)=W(nx,j,k)
     W(-1,j,k)=W(nx-1,j,k)
     W(-2,j,k)=W(nx-2,j,k)
    enddo
   enddo
  endif

  if ((BtypeE==DIRICHLET).or.(BtypeE==TURBULENTINLET)) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     W(nx+1,j,k)=0
     W(nx+2,j,k)=0
     W(nx+3,j,k)=0
    enddo
   enddo
  elseif (BtypeE==NOSLIP) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Solid wall
     W(nx+1,j,k)=-W(nx,j,k)
     W(nx+2,j,k)=-W(nx-1,j,k)
     W(nx+3,j,k)=-W(nx-2,j,k)
    enddo
   enddo
  elseif (BtypeE==NEUMANN.or.BtypeE==OUTLETBUFF) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Neumann inlet
     W(nx+1,j,k)=W(nx,j,k)
     W(nx+2,j,k)=W(nx,j,k)
     W(nx+3,j,k)=W(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do j=-2,ny+3
     W(nx+1,j,k)=W(nx,j,k)
     W(nx+2,j,k)=W(nx,j,k)
     W(nx+3,j,k)=W(nx,j,k)
    enddo
   enddo
  elseif (BtypeE==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do j=-2,ny+3
     W(nx+1,j,k)=W(1,j,k)
     W(nx+2,j,k)=W(2,j,k)
     W(nx+3,j,k)=W(3,j,k)
    enddo
   enddo
  endif

  if (BtypeS==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     W(i,0,k)=-W(i,1,k)
     W(i,-1,k)=-W(i,2,k)
     W(i,-2,k)=-W(i,3,k)
    enddo
   enddo
  elseif (BtypeS==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     W(i,0,k)=-W(i,1,k)
     W(i,-1,k)=-W(i,2,k)
     W(i,-2,k)=-W(i,3,k)
    enddo
   enddo
  elseif (BtypeS==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     W(i,0,k)=W(i,1,k)
     W(i,-1,k)=W(i,1,k)
     W(i,-2,k)=W(i,1,k)
    enddo
   enddo
  elseif (BtypeS==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     W(i,0,k)=W(i,1,k)
     W(i,-1,k)=W(i,1,k)
     W(i,-2,k)=W(i,1,k)
    enddo
   enddo
  elseif (BtypeS==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     W(i,0,k)=W(i,ny,k)
     W(i,-1,k)=W(i,ny-1,k)
     W(i,-2,k)=W(i,ny-2,k)
    enddo
   enddo
  endif

  if (BtypeN==DIRICHLET) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Dirichlet inlet
     W(i,ny+1,k)=-W(i,ny,k)
     W(i,ny+2,k)=-W(i,ny-1,k)
     W(i,ny+3,k)=-W(i,ny-2,k)
    enddo
   enddo
  elseif (BtypeN==NOSLIP) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Solid wall
     W(i,ny+1,k)=-W(i,ny,k)
     W(i,ny+2,k)=-W(i,ny-1,k)
     W(i,ny+3,k)=-W(i,ny-2,k)
    enddo
   enddo
  elseif (BtypeN==NEUMANN) then
   do k=-2,nz+3
    do i=-2,nx+3                       !Neumann inlet
     W(i,ny+1,k)=W(i,ny,k)
     W(i,ny+2,k)=W(i,ny,k)
     W(i,ny+3,k)=W(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==FREESLIP) then  !FREESLIP
   do k=-2,nz+3
    do i=-2,nx+3
     W(i,ny+1,k)=W(i,ny,k)
     W(i,ny+2,k)=W(i,ny,k)
     W(i,ny+3,k)=W(i,ny,k)
    enddo
   enddo
  elseif (BtypeN==PERIODIC) then  !Periodic BC
   do k=-2,nz+3
    do i=-2,nx+3
     W(i,ny+1,k)=W(i,1,k)
     W(i,ny+2,k)=W(i,2,k)
     W(i,ny+3,k)=W(i,3,k)
    enddo
   enddo
  endif

  if (BtypeB==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     W(i,j,0)=0
     W(i,j,-1)=-W(i,j,1)
     W(i,j,-2)=-W(i,j,2)
    enddo
   enddo
  elseif (BtypeB==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     W(i,j,0)=0
     W(i,j,-1)=-W(i,j,1)
     W(i,j,-2)=-W(i,j,2)
    enddo
   enddo
  elseif (BtypeB==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     W(i,j,0)=W(i,j,1)
     W(i,j,-1)=W(i,j,1)
     W(i,j,-2)=W(i,j,1)
    enddo
   enddo
  elseif (BtypeB==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     W(i,j,0)=0
     W(i,j,-1)=-W(i,j,1)
     W(i,j,-2)=-W(i,j,2)
    enddo
   enddo
  elseif (BtypeB==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     W(i,j,0)=W(i,j,nz)
     W(i,j,-1)=W(i,j,nz-1)
     W(i,j,-2)=W(i,j,nz-2)
    enddo
   enddo
  endif

  if (BtypeT==DIRICHLET) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Dirichlet inlet
     W(i,j,nz+1)=0
     W(i,j,nz+2)=-W(i,j,nz)
     W(i,j,nz+3)=-W(i,j,nz-1)
    enddo
   enddo
  elseif (BtypeT==NOSLIP) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Solid wall
     W(i,j,nz+1)=0
     W(i,j,nz+2)=-W(i,j,nz)
     W(i,j,nz+3)=-W(i,j,nz-1)
    enddo
   enddo
  elseif (BtypeT==NEUMANN) then
   do j=-2,ny+3
    do i=-2,nx+3                       !Neumann inlet
     W(i,j,nz+1)=W(i,j,nz)
     W(i,j,nz+2)=W(i,j,nz)
     W(i,j,nz+3)=W(i,j,nz)
    enddo
   enddo
  elseif (BtypeT==FREESLIPBUFF.or.BtypeT==FREESLIP) then  !FREESLIP
   do j=-2,ny+3
    do i=-2,nx+3
     W(i,j,nz+1)=0
     W(i,j,nz+2)=-W(i,j,nz)
     W(i,j,nz+3)=-W(i,j,nz-1)
    enddo
   enddo
  elseif (BtypeT==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     W(i,j,nz+1)=W(i,j,1)
     W(i,j,nz+2)=W(i,j,2)
     W(i,j,nz+3)=W(i,j,3)
    enddo
   enddo
  endif

  end subroutine BOUND_CONDW2




  pure subroutine BOUND_Phi(Phi)
  real(KND),intent(inout):: Phi(0:,0:,0:)
  integer i,j,k,nx,ny,nz
   nx=Prnx
   ny=Prny
   nz=Prnz


   if (BtypeW==PERIODIC) then
    do k=1,nz
     do j=1,ny                      !Periodic BC
      Phi(0,j,k)=Phi(nx,j,k)
     enddo
    enddo
  else
    do k=1,nz
     do j=1,ny                      !Other BCs
      Phi(0,j,k)=Phi(1,j,k)
     enddo
    enddo
   endif

   if (BtypeE==PERIODIC) then
    do k=1,nz
     do j=1,ny                      !Periodic BC
      Phi(nx+1,j,k)=Phi(1,j,k)
     enddo
    enddo
   else
    do k=1,nz
     do j=1,ny                      !Other BCs
      Phi(nx+1,j,k)=Phi(nx,j,k)
     enddo
    enddo
   endif

   if (BtypeS==PERIODIC) then
   do k=1,nz
     do i=1,nx                      !Periodic BC
      Phi(i,0,k)=Phi(i,ny,k)
     enddo
    enddo
   else
    do k=1,nz
     do i=1,nx                      !Other BCs
      Phi(i,0,k)=Phi(i,1,k)
     enddo
    enddo
   endif

   if (BtypeN==PERIODIC) then
   do k=1,nz
     do i=1,nx                      !Periodic BC
      Phi(i,ny+1,k)=Phi(i,1,k)
     enddo
    enddo
   else
    do k=1,nz
     do i=1,nx                      !Other BCs
      Phi(i,ny+1,k)=Phi(i,ny,k)
     enddo
    enddo
   endif

   if (BtypeB==PERIODIC) then
    do j=1,ny
     do i=1,nx                      !Periodic BC
      Phi(i,j,0)=Phi(i,j,nz)
     enddo
    enddo
   else
    do j=1,ny
     do i=1,nx                      !Other BCs
      Phi(i,j,0)=Phi(i,j,1)
     enddo
    enddo
   endif

   if (BtypeT==PERIODIC) then
    do j=1,ny
     do i=1,nx                      !Periodic BC
      Phi(i,j,nz+1)=Phi(i,j,1)
     enddo
    enddo
   else
    do j=1,ny
     do i=1,nx                      !Other BCs
      Phi(i,j,nz+1)=Phi(i,j,nz)
     enddo
    enddo
   endif
  end subroutine BOUND_Phi

  pure subroutine BOUND_Pr(Pr)
  real(KND),intent(inout):: Pr(1:,1:,1:)
  integer i,j,k,nx,ny,nz

  nx=Prnx
  ny=Prny
  nz=Prnz

  if (BtypeW==PERIODIC) then
   do k=1,nz
    do j=1,ny
     Pr(nx+1,j,k)=Pr(1,j,k)
    enddo
   enddo
  endif
  if (BtypeN==PERIODIC) then
   do k=1,nz
    do i=1,nx
     Pr(i,ny+1,k)=Pr(i,1,k)
    enddo
   enddo
  endif
  if (BtypeT==PERIODIC) then
   do j=1,ny
    do i=1,nx
     Pr(i,j,nz+1)=Pr(i,j,1)
    enddo
   enddo
  endif

  if (BtypeW==PERIODIC.and.BtypeN==PERIODIC) then
   do k=1,nz
      Pr(nx+1,ny+1,k)=Pr(1,1,k)
   enddo
  endif

   if (BtypeW==PERIODIC.and.BtypeT==PERIODIC) then
   do j=1,ny
      Pr(nx+1,j,nz+1)=Pr(1,j,1)
   enddo
  endif

  if (BtypeN==PERIODIC.and.BtypeT==PERIODIC) then
   do i=1,nx
      Pr(i,ny+1,nz+1)=Pr(i,1,1)
   enddo
  endif

  if (BtypeW==PERIODIC.and.BtypeN==PERIODIC.and.BtypeT==PERIODIC)  Pr(nx+1,ny+1,nz+1)=Pr(1,1,1)

  end subroutine BOUND_Pr



  pure subroutine BOUND_Q(Phi)
  real(KND),intent(inout):: Phi(0:,0:,0:)
  integer i,j,k,nx,ny,nz
   nx=Prnx
   ny=Prny
   nz=Prnz


   if (BtypeW==PERIODIC) then
    do k=1,nz
     do j=1,ny                      !Periodic BC
      Phi(nx,j,k)=Phi(0,j,k)+Phi(nx,j,k)
      Phi(0,j,k)=0
     enddo
    enddo
  else
    do k=1,nz
     do j=1,ny                      !Other BCs
      Phi(1,j,k)=Phi(1,j,k)+Phi(0,j,k)
      Phi(0,j,k)=0
     enddo
    enddo
   endif

   if (BtypeE==PERIODIC) then
    do k=1,nz
     do j=1,ny                      !Periodic BC
      Phi(1,j,k)=Phi(1,j,k)+Phi(nx+1,j,k)
      Phi(nx+1,j,k)=0
     enddo
    enddo
   else
    do k=1,nz
     do j=1,ny                      !Other BCs
      Phi(nx,j,k)=Phi(nx,j,k)+Phi(nx+1,j,k)
      Phi(nx+1,j,k)=0
     enddo
    enddo
   endif

   if (BtypeS==PERIODIC) then
   do k=1,nz
     do i=1,nx                      !Periodic BC
      Phi(i,ny,k)=Phi(i,ny,k)+Phi(i,0,k)
      Phi(i,0,k)=0
     enddo
    enddo
   else
    do k=1,nz
     do i=1,nx                      !Other BCs
      Phi(i,1,k)=Phi(i,1,k)+Phi(i,0,k)
      Phi(i,0,k)=0
     enddo
    enddo
   endif

   if (BtypeN==PERIODIC) then
   do k=1,nz
     do i=1,nx                      !Periodic BC
      Phi(i,1,k)=Phi(i,1,k)+Phi(i,ny+1,k)
      Phi(i,ny+1,k)=0
     enddo
    enddo
   else
    do k=1,nz
     do i=1,nx                      !Other BCs
      Phi(i,ny,k)=Phi(i,ny,k)+Phi(i,ny+1,k)
      Phi(i,ny+1,k)=0
     enddo
    enddo
   endif

   if (BtypeB==PERIODIC) then
    do j=1,ny
     do i=1,nx                      !Periodic BC
      Phi(i,j,nz)=Phi(i,j,nz)+Phi(i,j,0)
      Phi(i,j,0)=0
     enddo
    enddo
   else
    do j=1,ny
     do i=1,nx                      !Other BCs
      Phi(i,j,1)=Phi(i,j,1)+Phi(i,j,0)
      Phi(i,j,0)=0
     enddo
    enddo
   endif

   if (BtypeT==PERIODIC) then
    do j=1,ny
     do i=1,nx                      !Periodic BC
      Phi(i,j,1)=Phi(i,j,1)+Phi(i,j,nz+1)
      Phi(i,j,nz+1)=0
     enddo
    enddo
   else
    do j=1,ny
     do i=1,nx                      !Other BCs
      Phi(i,j,nz)=Phi(i,j,nz)+Phi(i,j,nz+1)
      Phi(i,j,nz+1)=0
     enddo
    enddo
   endif
  end subroutine BOUND_Q















  subroutine CONSTINLET
  !Dirichlet inlet condition
  !Uin must be allocated beforehand
  integer j,k
  do k=LBOUND(Uin,2),UBOUND(Uin,2)
   do j=LBOUND(Uin,1),UBOUND(Uin,1)
     Uin(j,k)=Uinlet
   enddo
  enddo
  do k=LBOUND(Vin,2),UBOUND(Vin,2)
   do j=LBOUND(Vin,1),UBOUND(Vin,1)
     Vin(j,k)=0
   enddo
  enddo
  do k=LBOUND(Win,2),UBOUND(Win,2)
   do j=LBOUND(Win,1),UBOUND(Win,1)
     Win(j,k)=0
   enddo
  enddo
  endsubroutine CONSTINLET

  subroutine SHEARINLET(G)
  real(KND) G
  integer j,k
  do k=LBOUND(Uin,2),UBOUND(Uin,2)
   do j=LBOUND(Uin,1),UBOUND(Uin,1)
     Uin(j,k)=G*(zPr(k)-((zW(Wnz+1)+zW(0))/2._KND))
   enddo
  enddo
  Vin=0
  Win=0
  endsubroutine SHEARINLET

  subroutine PARINLET
  integer j,k
  do k=LBOUND(Uin,2),UBOUND(Uin,2)
   do j=LBOUND(Uin,1),UBOUND(Uin,1)
     Uin(j,k)=1.5*Uinlet*(1-((lz/2._KND-(zPr(k)-zW(0)))/(lz/2._KND))**2)
   enddo
  enddo
  Vin=0
  Win=0
  endsubroutine PARINLET


end module BOUNDARIES

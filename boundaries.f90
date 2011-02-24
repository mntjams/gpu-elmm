module BOUNDARIES
use PARAMETERS
!use BLASIUS
use WALLMODELS
use SCALARS
use TURBINLET
use GEOMETRIC

implicit none

  integer nBoundUx,nBoundUy,nBoundUz,nBoundVx,nBoundVy,nBoundVz,nBoundWx,nBoundWy,nBoundWz
  integer iXupa,iXup,jYha,jYh


 contains 
  

   elemental subroutine GridCoords(xi,yj,zk,x,y,z)
   integer,intent(out):: xi,yj,zk
   real(KND),intent(in):: x,y,z
   integer i

   xi=Prnx
   do i=1,Prnx
    if (xU(i)>=x) then
                  xi=i
                  EXIT
                 endif
   enddo

   yj=Prny
   do i=1,Prny
    if (yV(i)>=y) then
                  yj=i
                  EXIT
                 endif
   enddo
   zk=Prnz
   do i=1,Prnz
    if (zW(i)>=z) then
                  zk=i
                  EXIT
                 endif
   enddo
   endsubroutine GridCoords


  subroutine BOUND_CONDU(U)
  real(KND),intent(INOUT):: U(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Unx
  ny=Uny
  nz=Unz   
 
  !!! oblasti nepokryte predchozim pri periodickych podminkach
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
  elseif (BtypeW==PERIODIC) then
   do k=-2,nz+3
    do j=-2,ny+3                       !Dirichlet inlet
     U(0,j,k)=Uin(j,k)
     U(-1,j,k)=Uin(j,k)
     U(-2,j,k)=Uin(j,k)
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
  elseif (BtypeE==NEUMANN) then   !Neumann outlet
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
  elseif (BtypeT==FREESLIP) then  !FREESLIP
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
  real(KND),intent(INOUT):: V(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Vnx
  ny=Vny
  nz=Vnz

  !!! oblasti nepokryte predchozim pri periodickych podminkach
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
  elseif (BtypeE==NEUMANN) then
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
  elseif (BtypeT==FREESLIP) then  !FREESLIP
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
  real(KND),intent(INOUT):: W(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Wnx
  ny=Wny
  nz=Wnz
    
  !!! oblasti nepokryte predchozim pri periodickych podminkach
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
  elseif (BtypeE==NEUMANN) then
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
  elseif (BtypeT==FREESLIP) then  !FREESLIP
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
  real(KND),intent(INOUT):: U(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Unx
  ny=Uny
  nz=Unz

  !!! oblasti nepokryte predchozim pri periodickych podminkach
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
  elseif (BtypeT==FREESLIP) then  !FREESLIP
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
  elseif (BtypeE==NEUMANN) then   !Neumann outlet

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
  real(KND),intent(INOUT):: V(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Vnx
  ny=Vny
  nz=Vnz

  !!! oblasti nepokryte predchozim pri periodickych podminkach
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
  elseif (BtypeE==NEUMANN) then
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
  elseif (BtypeT==FREESLIP) then  !FREESLIP
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
  real(KND),intent(INOUT):: W(-2:,-2:,-2:)
  integer i,j,k,nx,ny,nz

  nx=Wnx
  ny=Wny
  nz=Wnz

  !!! oblasti nepokryte predchozim pri periodickych podminkach
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
  elseif (BtypeE==NEUMANN) then
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
  elseif (BtypeT==FREESLIP) then  !FREESLIP
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
  real(KND),intent(INOUT):: Phi(0:,0:,0:)
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
  real(KND),intent(INOUT):: Pr(1:,1:,1:)
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
  real(KND),intent(INOUT):: Phi(0:,0:,0:)
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


  

  subroutine READBOUNDS
  real(KND),allocatable:: xU2(:),yV2(:),zW2(:)
  integer i,j,k,nx,ny,nz,nxup,nxdown,nyup,nydown,nzup,nzdown,io
  real(KND) P
  type(WMPoint):: WMP
  
  nx=Prnx-1
  ny=Prny-1
  nz=Prnz-1
  
  if (xgridfromfile) then
   OPEN(11,file="xgrid.txt")
   j=-1
   do
    read (11,*,iostat=io) P
    if (io==0) then
      j=j+1
    else
      EXIT
    endif
   enddo
   nx=j
   Prnx=nx
   Vnx=Prnx
   Wnx=Prnx
   if (BtypeE==PERIODIC) then
                         Unx=Prnx
                        else
                         Unx=Prnx-1
   endif
   CLOSE(11)
  endif

  if (ygridfromfile) then
   OPEN(11,file="ygrid.txt")
   j=-1
   do
    read (11,*,iostat=io) P
    if (io==0) then
      j=j+1
    else
      EXIT
    endif
   enddo
   ny=j
   Prny=ny
   Uny=Prny
   Wny=Prny
   if (BtypeN==PERIODIC) then
                         Vny=Prny
                        else
                         Vny=Prny-1
   endif
   CLOSE(11)
  endif

  if (zgridfromfile) then
   OPEN(11,file="zgrid.txt")
   j=-1
   do
    read (11,*,iostat=io) P
    if (io==0) then
      j=j+1
      write(*,*) j
    else
      EXIT
    endif
   enddo
   nz=j
   Prnz=nz
   Unz=Prnz
   Vnz=Prnz
   if (BtypeT==PERIODIC) then
                         Wnz=Prnz
                        else
                         Wnz=Prnz-1
   endif
   CLOSE(11)
  endif 





  allocate(xU2(-3:nx+4))
  allocate(yV2(-3:ny+4))
  allocate(zW2(-3:nz+4))



  if (xgridfromfile) then
   OPEN(11,file="xgrid.txt")
   do j=0,nx
    write(*,*) j
    read(11,*) xU2(j)
   enddo
   CLOSE(11)

   if (BtypeW==PERIODIC) then
    do j=-1,-3,-1
     xU2(j)=xU2(0)-(xU2(nx)-xU2(nx+j))
    enddo
   else
    do j=-1,-3,-1
     xU2(j)=xU2(0)-(xU2(0-j)-xU2(0))
    enddo
   endif

   if (BtypeE==PERIODIC) then   
    do j=nx+1,nx+4
     xU2(j)=xU2(nx)+(xU2(j-nx)-xU2(0))
    enddo
   else
    do j=nx+1,nx+4
     xU2(j)=xU2(nx)+(xU2(nx)-xU2(nx-(j-nx)))
    enddo
   endif

   x0=xU2(0)
  else
   forall (i=-3:nx+4)
    xU2(i)=(i)*dxmin+x0
   endforall
  endif


  if (ygridfromfile) then
   OPEN(11,file="ygrid.txt")
   do j=0,ny
    read(11,*) yV2(j)
   enddo
   CLOSE(11)

   if (BtypeS==PERIODIC) then
    do j=-1,-3,-1
     yV2(j)=yV2(0)-(yV2(ny)-yV2(ny+j))
    enddo
   else
    do j=-1,-3,-1
     yV2(j)=yV2(0)-(yV2(0-j)-yV2(0))
    enddo
   endif

   if (BtypeN==PERIODIC) then   
    do j=ny+1,ny+4
     yV2(j)=yV2(ny)+(yV2(j-ny)-yV2(0))
    enddo
   else
    do j=ny+1,ny+4
     yV2(j)=yV2(ny)+(yV2(ny)-yV2(ny-(j-ny)))
    enddo
   endif

   y0=yV2(0)
  else
   forall (j=-3:ny+4)
     yV2(j)=(j)*dymin+y0
   endforall
  endif


  if (zgridfromfile) then
   OPEN(11,file="zgrid.txt")
   do j=0,nz
    read(11,*) zW2(j)
   enddo
   CLOSE(11)

   if (BtypeB==PERIODIC) then
    do j=-1,-3,-1
     zW2(j)=zW2(0)-(zW2(nz)-zW2(nz+j))
    enddo
   else
    do j=-1,-3,-1
     zW2(j)=zW2(0)-(zW2(0-j)-zW2(0))
    enddo
   endif

   if (BtypeT==PERIODIC) then   
    do j=nz+1,nz+4
     zW2(j)=zW2(nz)+(zW2(j-nz)-zW2(0))
    enddo
   else
    do j=nz+1,nz+4
     zW2(j)=zW2(nz)+(zW2(nz)-zW2(nz-(j-nz)))
    enddo
   endif

   z0=zW2(0)
  else
   forall (k=-3:nz+4)
    zW2(k)=(k)*dzmin+z0
   endforall
  endif


  nxup=nx+1
  nxdown=0
  nyup=ny+1
  nydown=0
  nzup=nz+1
  nzdown=0

  
  allocate(xU(-3:nx+4))
  allocate(yV(-3:ny+4))
  allocate(zW(-3:nz+4))
  allocate(dxU(-2:nx+3))
  allocate(dyV(-2:ny+3))
  allocate(dzW(-2:nz+3))
  allocate(xPr(-2:nx+4),dxPr(-2:nx+4))
  allocate(yPr(-2:ny+4),dyPr(-2:ny+4))
  allocate(zPr(-2:nz+4),dzPr(-2:nz+4))

  xU=xU2(nxdown-3:nxup+3)
  yV=yV2(nydown-3:nyup+3)
  zW=zW2(nzdown-3:nzup+3)
  
  forall (i=-2:nx+4)
   xPr(i)=(xU(i-1)+xU(i))/2._KND
   dxPr(i)=xU(i)-xU(i-1)
  endforall
  forall (j=-2:ny+4)
   yPr(j)=(yV(j-1)+yV(j))/2._KND
   dyPr(j)=yV(j)-yV(j-1)
  endforall
  forall (k=-2:nz+4)
   zPr(k)=(zW(k-1)+zW(k))/2._KND
   dzPr(k)=zW(k)-zW(k-1)
  endforall
  forall (i=-2:nx+3)
   dxU(i)=xPr(i+1)-xPr(i)
  endforall 
  forall (j=-2:ny+3)
   dyV(j)=yPr(j+1)-yPr(j)
  endforall 
  forall (k=-2:nz+3)
   dzW(k)=zPr(k+1)-zPr(k)
  endforall 
    
  deallocate(xU2)
  deallocate(yV2)
  deallocate(zW2)

  
  allocate(Utype(-2:Unx+3,-2:Uny+3,-2:Unz+3))
  allocate(Vtype(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
  allocate(Wtype(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
  allocate(Prtype(0:Prnx+1,0:Prny+1,0:Prnz+1))
  Utype=0
  Vtype=0
  Wtype=0
  Prtype=0
  
!   write (*,*) Prnx
!   write (*,*) Prny
!   write (*,*) Prnz
!   write (*,*) Unx
!   write (*,*) Uny
!   write (*,*) Unz


  allocate(Uin(-2:Uny+3,-2:Unz+3),Vin(-2:Vny+3,-2:Vnz+3),Win(-2:Wny+3,-2:Wnz+3))
  if (buoyancy>0) allocate(Tempin(-2:Prny+3,-2:Prnz+3))
 
  select case (inlettype)
   case (NOINLET)
    Uin=0
    Vin=0
    Win=0
   case (SHEAR)
    call SHEARINLET(SHEARG)
   case (PARABOLIC)
    call PARINLET
   case (TURBULENTINLET)
    call GETTURBINLET
   case (INLETFROMFILE)
    call GETINLETFROMFILE(starttime)
   case default
    call CONSTINLET
  endselect

  nBoundUx=0
  nBoundVx=0
  nBoundWx=0
  nBoundUy=0
  nBoundVy=0
  nBoundWy=0
  nBoundUz=0
  nBoundVz=0
  nBoundWz=0
  
  


   allocate(REDBLACKU(0:Unx,0:Uny,0:Unz),REDBLACKV(0:Vnx,0:Vny,0:Vnz),REDBLACKW(0:Wnx,0:Wny,0:Wnz),REDBLACKPR(0:Prnx,0:Prny,0:Prnz))

   REDBLACKU(0,0,0)=.true.
   do k=1,Unz
    REDBLACKU(0,0,k)=.not.REDBLACKU(0,0,k-1)
    do j=1,Uny
     REDBLACKU(0,j,k)=.not.REDBLACKU(0,j-1,k)
     do i=1,Unx
      REDBLACKU(i,j,k)=.not.REDBLACKU(i-1,j,k)
     enddo
    enddo
   enddo

   REDBLACKV(0,0,0)=.true.
   do k=1,Vnz
    REDBLACKV(0,0,k)=.not.REDBLACKV(0,0,k-1)
    do j=1,Vny
     REDBLACKV(0,j,k)=.not.REDBLACKV(0,j-1,k)
     do i=1,Vnx
      REDBLACKV(i,j,k)=.not.REDBLACKV(i-1,j,k)
     enddo
    enddo
   enddo

   REDBLACKW(0,0,0)=.true.
   do k=1,Wnz
    REDBLACKW(0,0,k)=.not.REDBLACKW(0,0,k-1)
    do j=1,Wny
     REDBLACKW(0,j,k)=.not.REDBLACKW(0,j-1,k)
     do i=1,Wnx
      REDBLACKW(i,j,k)=.not.REDBLACKW(i-1,j,k)
     enddo
    enddo
   enddo
   RedBlackPr=.false.
   REDBLACKPr(0,0,0)=.true.
   do k=1,Prnz
    REDBLACKPr(0,0,k)=.not.REDBLACKPr(0,0,k-1)
    do j=1,Prny
     REDBLACKPr(0,j,k)=.not.REDBLACKPr(0,j-1,k)
     do i=1,Prnx
      REDBLACKPr(i,j,k)=.not.REDBLACKPr(i-1,j,k)
     enddo
    enddo
   enddo

   
!    if (tasktype==1)   then
!      Zu=Prnz/2
!      Yu=5*Prny/9+1!3*Prny/5+1
!      Yd=4*Prny/9+1!2*Prny/5
!      Xup=8*Prnx/20!5*Prnx/10
!      Xd=7*Prnx/20!4*Prnx/10
!      write(*,*) 1,Zu
!      write(*,*) Yd,Yu
!      write(*,*) Xd,Xup
!    elseif (tasktype==6) then
!      Zu=boxnz
!      Yu=(5*boxny)/4
!      Yd=boxny/4+1
!      Xup=(boxnx*3)/2
!      Xd=boxnx/2+1
!    endif





   if (buoyancy==1) then
    allocate(temperature(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    temperature=huge(1.0_KND)
      if (TBtypeB==CONSTFLUX.or.TBtypeB==DIRICHLET) then
       allocate(BsideTFLArr(-1:Prnx+2,-1:Prny+2))
       if (TBtypeB==CONSTFLUX) then
        BsideTFLArr=BsideTemp
       else
        BsideTFLArr=0
       endif
       if (TBtypeB==DIRICHLET) then
        allocate(BsideTArr(-1:Prnx+2,-1:Prny+2))
        BsideTArr=BsideTemp
       endif
      endif
   else
    allocate(temperature(0,0,0))
   endif
   
   if (computescalars>0) then
    allocate(Scalar(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
    Scalar=huge(1.0_KND)
    if (averaging==1) then
     allocate(Scalaravg(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
     Scalaravg=0
    endif
   else
    allocate(Scalar(0,0,0,0))
    if (averaging==1) then
     allocate(Scalaravg(0,0,0,0))
    endif
   endif





    allocate(WMP%depscalar(computescalars))
    WMP%depscalar=0

    if (BtypeW==NOSLIP) then    
     do k=1,Prnz
      do j=1,Prny
       WMP%x=1
       WMP%y=j
       WMP%z=k
       WMP%distx=(xPr(1)-xU(0))
       WMP%disty=0
       WMP%distz=0
       WMP%ustar=1
       WMP%z0=z0W
       call AddWMPoint(WMP)
      enddo
     enddo     
    endif
    
    if (BtypeE==NOSLIP) then
     do k=1,Prnz
      do j=1,Prny
       WMP%x=Prnx
       WMP%y=j
       WMP%z=k
       WMP%distx=(xPr(Prnx)-xU(Unx+1))
       WMP%disty=0
       WMP%distz=0
       WMP%ustar=1
       WMP%z0=z0E
       call AddWMPoint(WMP)
      enddo
     enddo     
    endif
   
    if (BtypeS==NOSLIP) then
     do k=1,Prnz
      do i=1,Prnx
       WMP%x=i
       WMP%y=1
       WMP%z=k
       WMP%distx=0
       WMP%disty=(yPr(1)-yV(0))
       WMP%distz=0
       WMP%ustar=1
       WMP%z0=z0S
       call AddWMPoint(WMP)
      enddo
     enddo
    elseif (BtypeS==DIRICHLET) then 
     do k=1,Prnz
      do i=1,Prnx
       WMP%x=i
       WMP%y=1
       WMP%z=k
       WMP%distx=0
       WMP%disty=(yPr(1)-yV(0))
       WMP%distz=0
       WMP%ustar=1
       WMP%wallu=SsideU
       WMP%wallv=SsideV
       WMP%wallw=SsideW
       WMP%z0=z0S
       call AddWMPoint(WMP)
      enddo
     enddo     
    endif
    
    if (BtypeN==NOSLIP) then
     do k=1,Prnz
      do i=1,Prnx
       WMP%x=i
       WMP%y=Prny
       WMP%z=k
       WMP%distx=0
       WMP%disty=(yPr(Prny)-yV(Vny+1))
       WMP%distz=0
       WMP%ustar=1
       WMP%z0=z0N
       call AddWMPoint(WMP)
      enddo
     enddo
    elseif (BtypeN==DIRICHLET) then
     do k=1,Prnz
      do i=1,Prnx
       WMP%x=i
       WMP%y=Prny
       WMP%z=k
       WMP%distx=0
       WMP%disty=(yPr(Prny)-yV(Vny+1))
       WMP%distz=0
       WMP%ustar=1
       WMP%wallu=NsideU
       WMP%wallv=NsideV
       WMP%wallw=NsideW
       WMP%z0=z0N
       call AddWMPoint(WMP)
      enddo
     enddo
    endif

    if (BtypeB==NOSLIP) then
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,0)==0) then
        WMP%x=i
        WMP%y=j
        WMP%z=1
        WMP%distx=0
        WMP%disty=0
        WMP%distz=(zPr(1)-zW(0))
        WMP%ustar=1
        WMP%z0=z0B
       if (TBtypeB==CONSTFLUX) then
        WMP%tempfl=BsideTemp
       else
        WMP%temp=0
       endif
       if (TBtypeB==DIRICHLET) then
        WMP%temp=BsideTemp
       endif
        call AddWMPoint(WMP)
       endif
      enddo
     enddo     
    elseif (BtypeB==DIRICHLET) then
     do j=1,Prny
      do i=1,Prnx
       WMP%x=i
       WMP%y=j
       WMP%z=1
       WMP%distx=0
       WMP%disty=0
       WMP%distz=(zPr(1)-zW(0))
       WMP%ustar=1
       WMP%wallu=BsideU
       WMP%wallv=BsideV
       WMP%wallw=BsideW
       WMP%z0=z0B
       if (TBtypeB==CONSTFLUX) then
        WMP%tempfl=BsideTemp
       else
        WMP%temp=0
       endif
       if (TBtypeB==DIRICHLET) then
        WMP%temp=BsideTemp
       endif

       call AddWMPoint(WMP)
      enddo
     enddo     
    endif
    
    if (BtypeT==NOSLIP) then
     do j=1,Prny
      do i=1,Prnx
       WMP%x=i
       WMP%y=j
       WMP%z=Prnz
       WMP%distx=0
       WMP%disty=0
       WMP%distz=(zPr(Prnz)-zW(Wnz+1))
       WMP%ustar=1
       WMP%z0=z0T
       call AddWMPoint(WMP)
      enddo
     enddo     
    elseif (BtypeT==DIRICHLET) then
     do j=1,Prny
      do i=1,Prnx
       WMP%x=i
       WMP%y=j
       WMP%z=Prnz
       WMP%distx=0
       WMP%disty=0
       WMP%distz=(zPr(Prnz)-zW(Wnz+1))
       WMP%ustar=1
       WMP%wallu=TsideU
       WMP%wallv=TsideV
       WMP%wallw=TsideW
       WMP%z0=z0T
       call AddWMPoint(WMP)
      enddo
     enddo     
    endif
        

   if (computescalars>0.and.pointscalsource==1) then
        call Gridcoords(scalsrci(:),scalsrcj(:),scalsrck(:),scalsrcx(:),scalsrcy(:),scalsrcz(:))
   endif

   call InitSolidBodies
   call GetSolidBodiesBC

  write (*,*) "set"
 end subroutine READBOUNDS
 

end module BOUNDARIES

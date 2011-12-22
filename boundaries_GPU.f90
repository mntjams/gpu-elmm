 pure subroutine BoundU_GPU(component,nx,ny,nz,Prny,Prnz,&
                             Btype,sideU,&
                             Uin,U,regime)
  implicit none
#ifdef __HMPP
  integer,parameter :: KND=4
  integer,parameter       :: NOSLIP=1, FREESLIP=2, PERIODIC=3, DIRICHLET=4, NEUMANN=5, CONSTFLUX=6,&  !boundary condition types
                               TURBULENTINLET=7, FREESLIPBUFF=8, OUTLETBUFF=9, INLETFROMFILE=10
  integer, parameter :: Ea=1,We=2,So=3,No=4,Bo=5,To=6
#endif
  integer,parameter       :: interm=2

  integer,intent(in)      :: regime     ! if == 2 then BC for U2, V2 and W2 in module tsteps
  integer,intent(in)      :: component  ! 1,2,3 for U,V,W
  integer,intent(in)      :: nx,ny,nz,Prny,Prnz
  integer,intent(in)      :: Btype(6)
  real(KND),intent(in)    :: sideU(3,6)
  real(KND),intent(in)    :: Uin(-2:Prny+3,-2:Prnz+3)
  real(KND),intent(inout) :: U(-2:nx+3,-2:ny+3,-2:nz+3)
  integer i,j,k

  !!! corners and edges for periodic conditions
  if (Btype(Ea)==PERIODIC.and.Btype(No)==PERIODIC.and.Btype(To)==PERIODIC) then
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

  if (Btype(Ea)==PERIODIC.and.Btype(No)==PERIODIC) then
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


  if (Btype(Ea)==PERIODIC.and.Btype(To)==PERIODIC) then
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


  if (Btype(No)==PERIODIC.and.Btype(To)==PERIODIC) then
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

  if (Btype(We)==DIRICHLET) then
    if (component==1.and.regime/=interm) then
       do k=-2,nz+3
        do j=-2,ny+3                       !Dirichlet inlet
         U(0,j,k)=Uin(j,k)
         U(-1,j,k)=Uin(j,k)
         U(-2,j,k)=Uin(j,k)
        enddo
       enddo
    else
       do k=-2,nz+3
        do j=-2,ny+3                       !Dirichlet inlet
         U(0,j,k)=0
         U(-1,j,k)=0
         U(-2,j,k)=0
        enddo
       enddo
    endif
  elseif (Btype(We)==NOSLIP.or.(component==1.and.Btype(We)==FREESLIP)) then
    if (component==1) then
       do k=-2,nz+3
        do j=-2,ny+3                       !Solid wall
         U(0,j,k)=0
         U(-1,j,k)=-U(1,j,k)
         U(-2,j,k)=-U(2,j,k)
        enddo
       enddo
    else
       do k=-2,nz+3
        do j=-2,ny+3                       !Solid wall
         U(0,j,k)=-U(1,j,k)
         U(-1,j,k)=-U(2,j,k)
         U(-2,j,k)=-U(3,j,k)
        enddo
       enddo
    endif
  elseif (Btype(We)==NEUMANN.or.(component/=1.and.Btype(We)==FREESLIP)) then
    if (component==1) then
       do k=-2,nz+3
        do j=-2,ny+3                       !Neumann inlet
         U(0,j,k)=U(1,j,k)
         U(-1,j,k)=U(1,j,k)
         U(-2,j,k)=U(1,j,k)
        enddo
       enddo
    else
       do k=-2,nz+3
        do j=-2,ny+3                       !Neumann inlet
         U(0,j,k)=U(1,j,k)
         U(-1,j,k)=U(1,j,k)
         U(-2,j,k)=U(1,j,k)
        enddo
       enddo
    endif
  elseif (Btype(We)==PERIODIC) then  !Periodic BC
    do k=-2,nz+3
     do j=-2,ny+3
      U(0,j,k)=U(nx,j,k)
      U(-1,j,k)=U(nx-1,j,k)
      U(-2,j,k)=U(nx-2,j,k)
     enddo
    enddo
  elseif (Btype(We)==TURBULENTINLET.or.Btype(We)==INLETFROMFILE) then
    if (regime/=interm) then
      do k=1,Prnz
       do j=1,Prny                       !Dirichlet inlet
        U(0,j,k)=Uin(j,k)
        U(-1,j,k)=Uin(j,k)
        U(-2,j,k)=Uin(j,k)
       enddo
      enddo
    else
      do k=1,Prnz
       do j=1,Prny                       !Dirichlet inlet
        U(0,j,k)=0
        U(-1,j,k)=0
        U(-2,j,k)=0
       enddo
      enddo
    endif
  endif


  if (Btype(Ea)==DIRICHLET) then
    if (component==1.and.regime/=interm) then
      do k=-2,nz+3
       do j=-2,ny+3                       !Dirichlet inlet
        U(nx+1,j,k)=Uin(j,k)
        U(nx+2,j,k)=Uin(j,k)
        U(nx+3,j,k)=Uin(j,k)
       enddo
      enddo
    else
      do k=-2,nz+3
       do j=-2,ny+3                       !Dirichlet inlet
        U(nx+1,j,k)=0
        U(nx+2,j,k)=0
        U(nx+3,j,k)=0
       enddo
      enddo
    endif
  elseif (Btype(Ea)==NOSLIP.or.(component==1.and.Btype(Ea)==FREESLIP)) then
    if (component==1) then
      do k=-2,nz+3
       do j=-2,ny+3                       !Solid wall
        U(nx+1,j,k)=0
        U(nx+2,j,k)=-U(nx,j,k)
        U(nx+3,j,k)=-U(nx-1,j,k)
       enddo
      enddo
    else
      do k=-2,nz+3
       do j=-2,ny+3                       !Solid wall
        U(nx+1,j,k)=-U(nx,j,k)
        U(nx+2,j,k)=-U(nx-1,j,k)
        U(nx+3,j,k)=-U(nx-2,j,k)
       enddo
      enddo
    endif
  elseif (Btype(Ea)==NEUMANN.or.Btype(Ea)==OUTLETBUFF.or.(component/=1.and.Btype(Ea)==FREESLIP)) then   !Neumann outlet
    do k=-2,nz+3
     do j=-2,ny+3
      U(nx+1,j,k)=U(nx,j,k)
      U(nx+2,j,k)=U(nx,j,k)
      U(nx+3,j,k)=U(nx,j,k)
     enddo
    enddo
  elseif (Btype(Ea)==PERIODIC) then  !Periodic BC
    do k=-2,nz+3
     do j=-2,ny+3
      U(nx+1,j,k)=U(1,j,k)
      U(nx+2,j,k)=U(2,j,k)
      U(nx+3,j,k)=U(3,j,k)
     enddo
    enddo
  elseif (Btype(Ea)==TURBULENTINLET) then
    if (regime/=interm) then
      do k=-2,nz+3
       do j=-2,ny+3                       !Dirichlet inlet
        U(nx+1,j,k)=Uin(j,k)
        U(nx+2,j,k)=Uin(j,k)
        U(nx+3,j,k)=Uin(j,k)
       enddo
      enddo
    else
     do k=-2,nz+3
       do j=-2,ny+3                       !Dirichlet inlet
        U(nx+1,j,k)=0
        U(nx+2,j,k)=0
        U(nx+3,j,k)=0
       enddo
      enddo
    endif
  endif

  if (Btype(So)==DIRICHLET) then
    if (component==2.and.regime/=interm) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,0,k)=sideU(component,So)
        U(i,-1,k)=sideU(component,So)+(sideU(component,So)-U(i,1,k))
        U(i,-2,k)=sideU(component,So)+(sideU(component,So)-U(i,2,k))
       enddo
      enddo
    else if (regime/=interm) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,0,k)=sideU(component,So)+(sideU(component,So)-U(i,1,k))
        U(i,-1,k)=sideU(component,So)+(sideU(component,So)-U(i,2,k))
        U(i,-2,k)=sideU(component,So)+(sideU(component,So)-U(i,3,k))
       enddo
      enddo
    else
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,0,k)=0
        U(i,-1,k)=0
        U(i,-2,k)=0
       enddo
      enddo
    endif
  elseif (Btype(So)==NOSLIP.or.(component==2.and.Btype(So)==FREESLIP)) then
    if (component==2) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Solid wall
        U(i,0,k)=0
        U(i,-1,k)=-U(i,1,k)
        U(i,-2,k)=-U(i,2,k)
       enddo
      enddo
    else
      do k=-2,nz+3
       do i=-2,nx+3                       !Solid wall
        U(i,0,k)=-U(i,1,k)
        U(i,-1,k)=-U(i,2,k)
        U(i,-2,k)=-U(i,3,k)
       enddo
      enddo
    endif
  elseif (Btype(So)==NEUMANN.or.(component/=2.and.Btype(So)==FREESLIP)) then
    do k=-2,nz+3
     do i=-2,nx+3                       !Neumann inlet
      U(i,0,k)=U(i,1,k)
      U(i,-1,k)=U(i,1,k)
      U(i,-2,k)=U(i,1,k)
     enddo
    enddo
  elseif (Btype(So)==PERIODIC) then  !Periodic BC
    do k=-2,nz+3
     do i=-2,nx+3
      U(i,0,k)=U(i,ny,k)
      U(i,-1,k)=U(i,ny-1,k)
      U(i,-2,k)=U(i,ny-2,k)
     enddo
    enddo
  endif


  if (Btype(No)==DIRICHLET) then
    if (component==2.and.regime/=interm) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,ny+1,k)=sideU(component,No)
        U(i,ny+2,k)=sideU(component,No)+(sideU(component,No)-U(i,ny,k))
        U(i,ny+3,k)=sideU(component,No)+(sideU(component,No)-U(i,ny-1,k))
       enddo
      enddo
    elseif (regime/=interm) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,ny+1,k)=sideU(component,No)+(sideU(component,No)-U(i,ny,k))
        U(i,ny+2,k)=sideU(component,No)+(sideU(component,No)-U(i,ny-1,k))
        U(i,ny+3,k)=sideU(component,No)+(sideU(component,No)-U(i,ny-2,k))
       enddo
      enddo
    else
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,ny+1,k)=0
        U(i,ny+2,k)=0
        U(i,ny+3,k)=0
       enddo
      enddo
    endif
  elseif (Btype(No)==NOSLIP.or.(component==2.and.Btype(So)==FREESLIP)) then
    if (component==2) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Solid wall
        U(i,ny+1,k)=0
        U(i,ny+2,k)=-U(i,ny,k)
        U(i,ny+3,k)=-U(i,ny-1,k)
       enddo
      enddo
    else
      do k=-2,nz+3
       do i=-2,nx+3                       !Solid wall
        U(i,ny+1,k)=-U(i,ny,k)
        U(i,ny+2,k)=-U(i,ny-1,k)
        U(i,ny+3,k)=-U(i,ny-2,k)
       enddo
      enddo
    endif
  elseif (Btype(No)==NEUMANN.or.(component/=2.and.Btype(No)==FREESLIP)) then
    do k=-2,nz+3
     do i=-2,nx+3                       !Neumann inlet
      U(i,ny+1,k)=U(i,ny,k)
      U(i,ny+2,k)=U(i,ny,k)
      U(i,ny+3,k)=U(i,ny,k)
     enddo
    enddo
  elseif (Btype(No)==PERIODIC) then  !Periodic BC
    do k=-2,nz+3
     do i=-2,nx+3
      U(i,ny+1,k)=U(i,1,k)
      U(i,ny+2,k)=U(i,2,k)
      U(i,ny+3,k)=U(i,3,k)
     enddo
    enddo
  endif


  if (Btype(Bo)==DIRICHLET) then
    if (component==3.and.regime/=interm) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,0)=sideU(component,Bo)
        U(i,j,-1)=sideU(component,Bo)+(sideU(component,Bo)-U(i,j,1))
        U(i,j,-2)=sideU(component,Bo)+(sideU(component,Bo)-U(i,j,2))
       enddo
      enddo
    elseif (regime/=interm) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,0)=sideU(component,Bo)+(sideU(component,Bo)-U(i,j,1))
        U(i,j,-1)=sideU(component,Bo)+(sideU(component,Bo)-U(i,j,2))
        U(i,j,-2)=sideU(component,Bo)+(sideU(component,Bo)-U(i,j,3))
       enddo
      enddo
    else
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,0)=0
        U(i,j,-1)=0
        U(i,j,-2)=0
       enddo
      enddo
    endif
  elseif (Btype(Bo)==NOSLIP.or.(component==3.and.Btype(Bo)==FREESLIP)) then
    if (component==3) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Solid wall
        U(i,j,0)=0
        U(i,j,-1)=-U(i,j,1)
        U(i,j,-2)=-U(i,j,2)
       enddo
      enddo
    else
      do j=-2,ny+3
       do i=-2,nx+3                       !Solid wall
        U(i,j,0)=-U(i,j,1)
        U(i,j,-1)=-U(i,j,2)
        U(i,j,-2)=-U(i,j,3)
       enddo
      enddo
    endif
  elseif (Btype(Bo)==NEUMANN.or.(component/=3.and.Btype(Bo)==FREESLIP)) then
    do j=-2,ny+3
     do i=-2,nx+3                       !Neumann inlet
      U(i,j,0)=U(i,j,1)
      U(i,j,-1)=U(i,j,1)
      U(i,j,-2)=U(i,j,1)
     enddo
    enddo
  elseif (Btype(Bo)==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     U(i,j,0)=U(i,j,nz)
     U(i,j,-1)=U(i,j,nz-1)
     U(i,j,-2)=U(i,j,nz-2)
    enddo
   enddo
  endif

  if (Btype(To)==DIRICHLET) then
    if (component==3.and.regime/=interm) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,nz+1)=sideU(component,To)
        U(i,j,nz+2)=sideU(component,To)+(sideU(component,To)-U(i,j,nz))
        U(i,j,nz+3)=sideU(component,To)+(sideU(component,To)-U(i,j,nz-1))
       enddo
      enddo
    elseif (regime/=interm) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,nz+1)=sideU(component,To)+(sideU(component,To)-U(i,j,nz))
        U(i,j,nz+2)=sideU(component,To)+(sideU(component,To)-U(i,j,nz-1))
        U(i,j,nz+3)=sideU(component,To)+(sideU(component,To)-U(i,j,nz-2))
       enddo
      enddo
    else
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,nz+1)=0
        U(i,j,nz+2)=0
        U(i,j,nz+3)=0
       enddo
      enddo
    endif
  elseif (Btype(To)==NOSLIP.or.(component==3.and.(Btype(To)==FREESLIP.or.Btype(To)==FREESLIPBUFF))) then
    if (component==3) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Solid wall
        U(i,j,nz+1)=0
        U(i,j,nz+2)=-U(i,j,nz)
        U(i,j,nz+3)=-U(i,j,nz-1)
       enddo
      enddo
    else
      do j=-2,ny+3
       do i=-2,nx+3                       !Solid wall
        U(i,j,nz+1)=-U(i,j,nz)
        U(i,j,nz+2)=-U(i,j,nz-1)
        U(i,j,nz+3)=-U(i,j,nz-2)
       enddo
      enddo
    endif
  elseif (Btype(To)==NEUMANN.or.(component/=3.and.(Btype(To)==FREESLIP.or.Btype(To)==FREESLIPBUFF))) then
    do j=-2,ny+3
     do i=-2,nx+3                       !Neumann inlet
      U(i,j,nz+1)=U(i,j,nz)
      U(i,j,nz+2)=U(i,j,nz)
      U(i,j,nz+3)=U(i,j,nz)
     enddo
    enddo
  elseif (Btype(To)==PERIODIC) then  !Periodic BC
   do j=-2,ny+3
    do i=-2,nx+3
     U(i,j,nz+1)=U(i,j,1)
     U(i,j,nz+2)=U(i,j,2)
     U(i,j,nz+3)=U(i,j,3)
    enddo
   enddo
  endif

  end subroutine BoundU_GPU



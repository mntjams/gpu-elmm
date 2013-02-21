module BOUNDARIES
use PARAMETERS
use WALLMODELS

implicit none


  private
  public GridCoords, GridCoordsU, GridCoordsV, GridCoordsW,&
         BoundU, Bound_Phi, Bound_Pr, Bound_Q,&
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



  recursive subroutine BoundU(component,U,reg)
  integer,intent(in)           :: component
  real(KND),intent(inout)      :: U(-2:,-2:,-2:)
  integer,optional,intent(in) :: reg
  integer regime,i,j,k,nx,ny,nz
  integer,parameter :: interm = 2

  if (present(reg)) then
    regime=reg
  else
    regime=0
  endif

  if (component==1) then
    nx=Unx
    ny=Uny
    nz=Unz
  elseif (component==2) then
    nx=Vnx
    ny=Vny
    nz=Vnz
  else
    nx=Wnx
    ny=Wny
    nz=Wnz
  endif

  !$omp parallel private(i,j,k)

  !!! corners and edges for periodic conditions
  if (Btype(Ea)==PERIODIC.and.Btype(No)==PERIODIC.and.Btype(To)==PERIODIC) then
  !$omp sections
  !$omp section
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j-ny,k-nz)
     enddo
    enddo
   enddo
  !$omp section
   do k=-2,0
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j-ny,k+nz)
     enddo
    enddo
   enddo
  !$omp section
   do k=nz+1,nz+3
    do j=-2,0
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j+ny,k-nz)
     enddo
    enddo
   enddo
  !$omp section
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=-2,0
      U(i,j,k)=U(i+nx,j-ny,k-nz)
     enddo
    enddo
   enddo
  !$omp section
   do k=-2,0
    do j=-2,0
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j+ny,k+nz)
     enddo
    enddo
   enddo
  !$omp section
   do k=nz+1,nz+3
    do j=-2,0
     do i=-2,0
      U(i,j,k)=U(i+nx,j+ny,k-nz)
     enddo
    enddo
   enddo
  !$omp section
   do k=-2,0
    do j=ny+1,ny+3
     do i=-2,0
      U(i,j,k)=U(i+nx,j-ny,k+nz)
     enddo
    enddo
   enddo
  !$omp section
   do k=-2,0
    do j=-2,0
     do i=-2,0
      U(i,j,k)=U(i+nx,j+ny,k+nz)
     enddo
    enddo
   enddo
   !$omp end sections
  endif

  if (Btype(Ea)==PERIODIC.and.Btype(No)==PERIODIC) then
   !$omp sections
   !$omp section
   do k=1,nz
    do j=ny+1,ny+3
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j-ny,k)
     enddo
    enddo
   enddo
   !$omp section
   do k=1,nz
    do j=-2,0
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j+ny,k)
     enddo
    enddo
   enddo
   !$omp section
   do k=1,nz
    do j=ny+1,ny+3
     do i=-2,0
      U(i,j,k)=U(i+nx,j-ny,k)
     enddo
    enddo
   enddo
   !$omp section
   do k=1,nz
    do j=-2,0
     do i=-2,0
      U(i,j,k)=U(i+nx,j+ny,k)
     enddo
    enddo
   enddo
   !$omp end sections
  endif


  if (Btype(Ea)==PERIODIC.and.Btype(To)==PERIODIC) then
   !$omp sections
   !$omp section
   do k=nz+1,nz+3
    do j=1,ny
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j,k-nz)
     enddo
    enddo
   enddo
   !$omp section
   do k=-2,0
    do j=1,ny
     do i=nx+1,nx+3
      U(i,j,k)=U(i-nx,j,k+nz)
     enddo
    enddo
   enddo
   !$omp section
   do k=nz+1,nz+3
    do j=1,ny
     do i=-2,0
      U(i,j,k)=U(i+nx,j,k-nz)
     enddo
    enddo
   enddo
   !$omp section
   do k=-2,0
    do j=1,ny
     do i=-2,0
      U(i,j,k)=U(i+nx,j,k+nz)
     enddo
    enddo
   enddo
   !$omp end sections
  endif


  if (Btype(No)==PERIODIC.and.Btype(To)==PERIODIC) then
   !$omp sections
   !$omp section
   do k=nz+1,nz+3
    do j=ny+1,ny+3
     do i=1,nx
      U(i,j,k)=U(i,j-ny,k-nz)
     enddo
    enddo
   enddo
   !$omp section
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
   !$omp section
   do k=-2,0
    do j=-2,0
     do i=1,nx
      U(i,j,k)=U(i,j+ny,k+nz)
     enddo
    enddo
   enddo
   !$omp end sections
  endif

  !$omp sections
  !$omp section

  if (Btype(Bo)==DIRICHLET .and. regime/=interm) then
    if (component==3) then
      do j=1,ny
       do i=1,nx                       !Dirichlet inlet
        U(i,j,0)=sideU(component,Bo)
        U(i,j,-1)=sideU(component,Bo)+(sideU(component,Bo)-U(i,j,1))
        U(i,j,-2)=sideU(component,Bo)+(sideU(component,Bo)-U(i,j,2))
       enddo
      enddo
    else
      do j=1,ny
       do i=1,nx                       !Dirichlet inlet
        U(i,j,0)=sideU(component,Bo)+(sideU(component,Bo)-U(i,j,1))
        U(i,j,-1)=sideU(component,Bo)+(sideU(component,Bo)-U(i,j,2))
        U(i,j,-2)=sideU(component,Bo)+(sideU(component,Bo)-U(i,j,3))
       enddo
      enddo
    endif
  elseif (Btype(Bo)==NOSLIP .or. (component==3.and.Btype(Bo)==FREESLIP) .or. &
             (Btype(Bo)==DIRICHLET.and.regime==interm)) then
    if (component==3) then
      do j=1,ny
       do i=1,nx                       !Solid wall
        U(i,j,0)=0
        U(i,j,-1)=-U(i,j,1)
        U(i,j,-2)=-U(i,j,2)
       enddo
      enddo
    else
      do j=1,ny
       do i=1,nx                       !Solid wall
        U(i,j,0)=-U(i,j,1)
        U(i,j,-1)=-U(i,j,2)
        U(i,j,-2)=-U(i,j,3)
       enddo
      enddo
    endif
  elseif (Btype(Bo)==NEUMANN.or.(component/=3.and.Btype(Bo)==FREESLIP)) then
    do j=1,ny
     do i=1,nx                       !Neumann inlet
      U(i,j,0)=U(i,j,1)
      U(i,j,-1)=U(i,j,1)
      U(i,j,-2)=U(i,j,1)
     enddo
    enddo
  elseif (Btype(Bo)==PERIODIC) then  !Periodic BC
   do j=1,ny
    do i=1,nx
     U(i,j,0)=U(i,j,nz)
     U(i,j,-1)=U(i,j,nz-1)
     U(i,j,-2)=U(i,j,nz-2)
    enddo
   enddo
  endif

  !$omp section
  if (Btype(To)==DIRICHLET.and.regime/=interm) then
    if (component==3) then
      do j=1,ny
       do i=1,nx                       !Dirichlet inlet
        U(i,j,nz+1)=sideU(component,To)
        U(i,j,nz+2)=sideU(component,To)+(sideU(component,To)-U(i,j,nz))
        U(i,j,nz+3)=sideU(component,To)+(sideU(component,To)-U(i,j,nz-1))
       enddo
      enddo
    else
      do j=1,ny
       do i=1,nx                       !Dirichlet inlet
        U(i,j,nz+1)=sideU(component,To)+(sideU(component,To)-U(i,j,nz))
        U(i,j,nz+2)=sideU(component,To)+(sideU(component,To)-U(i,j,nz-1))
        U(i,j,nz+3)=sideU(component,To)+(sideU(component,To)-U(i,j,nz-2))
       enddo
      enddo
    endif
  elseif (Btype(To)==NOSLIP .or. (component==3.and.(Btype(To)==FREESLIP.or.Btype(To)==FREESLIPBUFF)) .or. &
           (Btype(To)==DIRICHLET.and.regime==interm) ) then
    if (component==3) then
      do j=1,ny
       do i=1,nx                       !Solid wall
        U(i,j,nz+1)=0
        U(i,j,nz+2)=-U(i,j,nz)
        U(i,j,nz+3)=-U(i,j,nz-1)
       enddo
      enddo
    else
      do j=1,ny
       do i=1,nx                       !Solid wall
        U(i,j,nz+1)=-U(i,j,nz)
        U(i,j,nz+2)=-U(i,j,nz-1)
        U(i,j,nz+3)=-U(i,j,nz-2)
       enddo
      enddo
    endif
  elseif (Btype(To)==NEUMANN.or.(component/=3.and.(Btype(To)==FREESLIP.or.Btype(To)==FREESLIPBUFF))) then
    do j=1,ny
     do i=1,nx                       !Neumann inlet
      U(i,j,nz+1)=U(i,j,nz)
      U(i,j,nz+2)=U(i,j,nz)
      U(i,j,nz+3)=U(i,j,nz)
     enddo
    enddo
  elseif (Btype(To)==PERIODIC) then  !Periodic BC
   do j=1,ny
    do i=1,nx
     U(i,j,nz+1)=U(i,j,1)
     U(i,j,nz+2)=U(i,j,2)
     U(i,j,nz+3)=U(i,j,3)
    enddo
   enddo
  endif
  !$omp end sections

  !$omp sections
  !$omp section
  if (Btype(So)==DIRICHLET.and.regime/=interm) then
    if (component==2) then
      do k=-2,nz+3
       do i=1,nx                       !Dirichlet inlet
        U(i,0,k)=sideU(component,So)
        U(i,-1,k)=sideU(component,So)+(sideU(component,So)-U(i,1,k))
        U(i,-2,k)=sideU(component,So)+(sideU(component,So)-U(i,2,k))
       enddo
      enddo
    else
      do k=-2,nz+3
       do i=1,nx                       !Dirichlet inlet
        U(i,0,k)=sideU(component,So)+(sideU(component,So)-U(i,1,k))
        U(i,-1,k)=sideU(component,So)+(sideU(component,So)-U(i,2,k))
        U(i,-2,k)=sideU(component,So)+(sideU(component,So)-U(i,3,k))
       enddo
      enddo
    endif
  elseif (Btype(So)==NOSLIP .or. (component==2.and.Btype(So)==FREESLIP) .or. &
             (Btype(So)==DIRICHLET.and.regime==interm)) then
    if (component==2) then
      do k=-2,nz+3
       do i=1,nx                       !Solid wall
        U(i,0,k)=0
        U(i,-1,k)=-U(i,1,k)
        U(i,-2,k)=-U(i,2,k)
       enddo
      enddo
    else
      do k=-2,nz+3
       do i=1,nx                       !Solid wall
        U(i,0,k)=-U(i,1,k)
        U(i,-1,k)=-U(i,2,k)
        U(i,-2,k)=-U(i,3,k)
       enddo
      enddo
    endif
  elseif (Btype(So)==NEUMANN.or.(component/=2.and.Btype(So)==FREESLIP)) then
    do k=-2,nz+3
     do i=1,nx                       !Neumann inlet
      U(i,0,k)=U(i,1,k)
      U(i,-1,k)=U(i,1,k)
      U(i,-2,k)=U(i,1,k)
     enddo
    enddo
  elseif (Btype(So)==PERIODIC) then  !Periodic BC
    do k=-2,nz+3
     do i=1,nx
      U(i,0,k)=U(i,ny,k)
      U(i,-1,k)=U(i,ny-1,k)
      U(i,-2,k)=U(i,ny-2,k)
     enddo
    enddo
  endif


  !$omp section
  if (Btype(No)==DIRICHLET.and.regime/=interm) then
    if (component==2) then
      do k=-2,nz+3
       do i=1,nx                       !Dirichlet inlet
        U(i,ny+1,k)=sideU(component,No)
        U(i,ny+2,k)=sideU(component,No)+(sideU(component,No)-U(i,ny,k))
        U(i,ny+3,k)=sideU(component,No)+(sideU(component,No)-U(i,ny-1,k))
       enddo
      enddo
    else
      do k=-2,nz+3
       do i=1,nx                       !Dirichlet inlet
        U(i,ny+1,k)=sideU(component,No)+(sideU(component,No)-U(i,ny,k))
        U(i,ny+2,k)=sideU(component,No)+(sideU(component,No)-U(i,ny-1,k))
        U(i,ny+3,k)=sideU(component,No)+(sideU(component,No)-U(i,ny-2,k))
       enddo
      enddo
    endif
  elseif (Btype(No)==NOSLIP .or. (component==2.and.Btype(So)==FREESLIP) .or. &
           (Btype(No)==DIRICHLET.and.regime==interm)) then
    if (component==2) then
      do k=-2,nz+3
       do i=1,nx                       !Solid wall
        U(i,ny+1,k)=0
        U(i,ny+2,k)=-U(i,ny,k)
        U(i,ny+3,k)=-U(i,ny-1,k)
       enddo
      enddo
    else
      do k=-2,nz+3
       do i=1,nx                       !Solid wall
        U(i,ny+1,k)=-U(i,ny,k)
        U(i,ny+2,k)=-U(i,ny-1,k)
        U(i,ny+3,k)=-U(i,ny-2,k)
       enddo
      enddo
    endif
  elseif (Btype(No)==NEUMANN .or. (component/=2.and.Btype(No)==FREESLIP)) then
    do k=-2,nz+3
     do i=1,nx                       !Neumann inlet
      U(i,ny+1,k)=U(i,ny,k)
      U(i,ny+2,k)=U(i,ny,k)
      U(i,ny+3,k)=U(i,ny,k)
     enddo
    enddo
  elseif (Btype(No)==PERIODIC) then  !Periodic BC
    do k=-2,nz+3
     do i=1,nx
      U(i,ny+1,k)=U(i,1,k)
      U(i,ny+2,k)=U(i,2,k)
      U(i,ny+3,k)=U(i,3,k)
     enddo
    enddo
  endif
  !$omp end sections

  !$omp sections
  !$omp section
  if (Btype(We)==DIRICHLET.and.regime/=interm) then
    if (component==1) then
       do k=-2,nz+3
        do j=-2,ny+3                       !Dirichlet inlet
         U(0,j,k)=Uin(j,k)
         U(-1,j,k)=Uin(j,k)+(Uin(j,k)-U(1,j,k))
         U(-2,j,k)=Uin(j,k)+(Uin(j,k)-U(2,j,k))
        enddo
       enddo
    else
       do k=-2,nz+3
        do j=-2,ny+3                       !Dirichlet inlet
         U(0,j,k)=-U(1,j,k)
         U(-1,j,k)=-U(2,j,k)
         U(-2,j,k)=-U(3,j,k)
        enddo
       enddo
    endif
  elseif (Btype(We)==NOSLIP .or. (component==1.and.Btype(We)==FREESLIP) .or. &
           (Btype(We)==DIRICHLET.and.regime==interm)) then
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
      if (component==1) then
        do k=-2,nz+3
         do j=-2,ny+3
          U(0,j,k)=Uin(j,k)
          U(-1,j,k)=Uin(j,k)+(Uin(j,k)-U(1,j,k))
          U(-2,j,k)=Uin(j,k)+(Uin(j,k)-U(2,j,k))
         enddo
        enddo
      else
        do k=-2,nz+3
         do j=-2,ny+3
          U(0,j,k)=Uin(j,k)+(Uin(j,k)-U(1,j,k))
          U(-1,j,k)=Uin(j,k)+(Uin(j,k)-U(2,j,k))
          U(-2,j,k)=Uin(j,k)+(Uin(j,k)-U(3,j,k))
         enddo
        enddo
      endif
    else
      if (component==1) then
         do k=-2,nz+3
          do j=-2,ny+3
           U(0,j,k)=0
           U(-1,j,k)=-U(1,j,k)
           U(-2,j,k)=-U(2,j,k)
          enddo
         enddo
      else
         do k=-2,nz+3
          do j=-2,ny+3
           U(0,j,k)=-U(1,j,k)
           U(-1,j,k)=-U(2,j,k)
           U(-2,j,k)=-U(3,j,k)
          enddo
         enddo
      endif
    endif
  endif


  !$omp section
  if (Btype(Ea)==DIRICHLET.and.regime/=interm) then
    if (component==1) then
      do k=-2,nz+3
       do j=-2,ny+3                       !Dirichlet inlet
        U(nx+1,j,k)=Uin(j,k)
        U(nx+2,j,k)=Uin(j,k)+(Uin(j,k)-U(nx,j,k))
        U(nx+3,j,k)=Uin(j,k)+(Uin(j,k)-U(nx-1,j,k))
       enddo
      enddo
    else
      do k=-2,nz+3
       do j=-2,ny+3                       !Dirichlet inlet
        U(nx+1,j,k)=-U(nx,j,k)
        U(nx+2,j,k)=-U(nx-1,j,k)
        U(nx+3,j,k)=-U(nx-2,j,k)
       enddo
      enddo
    endif
  elseif (Btype(Ea)==NOSLIP .or. (component==1.and.Btype(Ea)==FREESLIP) .or. &
            (Btype(Ea)==DIRICHLET.and.regime==interm)) then
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
      if (component==1) then
        do k=-2,nz+3
         do j=-2,ny+3
          U(nx+1,j,k)=Uin(j,k)
          U(nx+2,j,k)=Uin(j,k)+(Uin(j,k)-U(nx,j,k))
          U(nx+3,j,k)=Uin(j,k)+(Uin(j,k)-U(nx-1,j,k))
         enddo
        enddo
      else
        do k=-2,nz+3
         do j=-2,ny+3
          U(nx+1,j,k)=Uin(j,k)+(Uin(j,k)-U(nx,j,k))
          U(nx+2,j,k)=Uin(j,k)+(Uin(j,k)-U(nx-1,j,k))
          U(nx+3,j,k)=Uin(j,k)+(Uin(j,k)-U(nx-2,j,k))
         enddo
        enddo
      endif
    else
      if (component==1) then
         do k=-2,nz+3
          do j=-2,ny+3
           U(nx+1,j,k)=0
           U(nx+2,j,k)=-U(nx,j,k)
           U(nx+3,j,k)=-U(nx-1,j,k)
          enddo
         enddo
      else
         do k=-2,nz+3
          do j=-2,ny+3
           U(nx+1,j,k)=-U(nx,j,k)
           U(nx+2,j,k)=-U(nx-1,j,k)
           U(nx+3,j,k)=-U(nx-2,j,k)
          enddo
         enddo
      endif
    endif
  endif
  !$omp end sections
  !$omp end parallel

  end subroutine BoundU



  pure subroutine BOUND_Phi(Phi)
  real(KND),intent(inout):: Phi(0:,0:,0:)
  integer i,j,k,nx,ny,nz
   nx=Prnx
   ny=Prny
   nz=Prnz


   if (Btype(We)==PERIODIC) then
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

   if (Btype(Ea)==PERIODIC) then
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

   if (Btype(So)==PERIODIC) then
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

   if (Btype(No)==PERIODIC) then
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

   if (Btype(Bo)==PERIODIC) then
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

   if (Btype(To)==PERIODIC) then
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

  if (Btype(We)==PERIODIC) then
   do k=1,nz
    do j=1,ny
     Pr(nx+1,j,k)=Pr(1,j,k)
    enddo
   enddo
  endif
  if (Btype(No)==PERIODIC) then
   do k=1,nz
    do i=1,nx
     Pr(i,ny+1,k)=Pr(i,1,k)
    enddo
   enddo
  endif
  if (Btype(To)==PERIODIC) then
   do j=1,ny
    do i=1,nx
     Pr(i,j,nz+1)=Pr(i,j,1)
    enddo
   enddo
  endif

  if (Btype(We)==PERIODIC.and.Btype(No)==PERIODIC) then
   do k=1,nz
      Pr(nx+1,ny+1,k)=Pr(1,1,k)
   enddo
  endif

   if (Btype(We)==PERIODIC.and.Btype(To)==PERIODIC) then
   do j=1,ny
      Pr(nx+1,j,nz+1)=Pr(1,j,1)
   enddo
  endif

  if (Btype(No)==PERIODIC.and.Btype(To)==PERIODIC) then
   do i=1,nx
      Pr(i,ny+1,nz+1)=Pr(i,1,1)
   enddo
  endif

  if (Btype(We)==PERIODIC.and.Btype(No)==PERIODIC.and.Btype(To)==PERIODIC)  Pr(nx+1,ny+1,nz+1)=Pr(1,1,1)

  end subroutine BOUND_Pr



  pure subroutine BOUND_Q(Phi)
  real(KND),intent(inout):: Phi(0:,0:,0:)
  integer i,j,k,nx,ny,nz
   nx=Prnx
   ny=Prny
   nz=Prnz


   if (Btype(We)==PERIODIC) then
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

   if (Btype(Ea)==PERIODIC) then
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

   if (Btype(So)==PERIODIC) then
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

   if (Btype(No)==PERIODIC) then
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

   if (Btype(Bo)==PERIODIC) then
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

   if (Btype(To)==PERIODIC) then
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
    real(KND) lz

    lz = zW(Prnz) - zW(0)

    do k=LBOUND(Uin,2),UBOUND(Uin,2)
     do j=LBOUND(Uin,1),UBOUND(Uin,1)
       Uin(j,k)=1.5*Uinlet*(1-((lz/2._KND-(zPr(k)-zW(0)))/(lz/2._KND))**2)
     enddo
    enddo

    Vin=0
    Win=0

  endsubroutine PARINLET


end module BOUNDARIES

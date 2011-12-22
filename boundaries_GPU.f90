 pure subroutine BoundU_GPU(component,nx,ny,nz,Prny,Prnz,&
                             BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT,&
                             SsideU,NsideU,BsideU,TsideU,&
                             SsideV,NsideV,BsideV,TsideV,&
                             SsideW,NsideW,BsideW,TsideW,&
                             Uin,U,regime)
  implicit none
#ifdef __HMPP
  integer,parameter :: KND=4
#endif
  integer,parameter       :: NOSLIP=1, FREESLIP=2, PERIODIC=3, DIRICHLET=4, NEUMANN=5, CONSTFLUX=6,&  !boundary condition types
                               TURBULENTINLET=7, FREESLIPBUFF=8, OUTLETBUFF=9, INLETFROMFILE=10
  integer,parameter       :: interm=2

  integer,intent(in)      :: regime     ! if == 2 then BC for U2, V2 and W2 in module tsteps
  integer,intent(in)      :: component  ! 1,2,3 for U,V,W
  integer,intent(in)      :: nx,ny,nz,Prny,Prnz
  integer,intent(in)      :: BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT
  real(KND),intent(in)    :: SsideU,NsideU,BsideU,TsideU
  real(KND),intent(in)    :: SsideV,NsideV,BsideV,TsideV
  real(KND),intent(in)    :: SsideW,NsideW,BsideW,TsideW
  real(KND),intent(in)    :: Uin(-2:Prny+3,-2:Prnz+3)
  real(KND),intent(inout) :: U(-2:nx+3,-2:ny+3,-2:nz+3)
  integer i,j,k

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
  elseif (BtypeW==NOSLIP.or.(component==1.and.BtypeW==FREESLIP)) then
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
  elseif (BtypeW==NEUMANN.or.(component/=1.and.BtypeW==FREESLIP)) then
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
  elseif (BtypeW==PERIODIC) then  !Periodic BC
    do k=-2,nz+3
     do j=-2,ny+3
      U(0,j,k)=U(nx,j,k)
      U(-1,j,k)=U(nx-1,j,k)
      U(-2,j,k)=U(nx-2,j,k)
     enddo
    enddo
  elseif (BtypeW==TURBULENTINLET.or.BtypeW==INLETFROMFILE) then
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


  if (BtypeE==DIRICHLET) then
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
  elseif (BtypeE==NOSLIP.or.(component==1.and.BtypeE==FREESLIP)) then
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
  elseif (BtypeE==NEUMANN.or.BtypeE==OUTLETBUFF.or.(component/=1.and.BtypeE==FREESLIP)) then   !Neumann outlet
    do k=-2,nz+3
     do j=-2,ny+3
      U(nx+1,j,k)=U(nx,j,k)
      U(nx+2,j,k)=U(nx,j,k)
      U(nx+3,j,k)=U(nx,j,k)
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

  if (BtypeS==DIRICHLET) then
    if (component==1.and.regime/=interm) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,0,k)=SsideU+(SsideU-U(i,1,k))
        U(i,-1,k)=SsideU+(SsideU-U(i,2,k))
        U(i,-2,k)=SsideU+(SsideU-U(i,3,k))
       enddo
      enddo
    elseif (component==2.and.regime/=interm) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,0,k)=SsideV
        U(i,-1,k)=SsideV+(SsideV-U(i,1,k))
        U(i,-2,k)=SsideV+(SsideV-U(i,2,k))
       enddo
      enddo
    else if (component==3.and.regime/=interm) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,0,k)=SsideW+(SsideW-U(i,1,k))
        U(i,-1,k)=SsideW+(SsideW-U(i,2,k))
        U(i,-2,k)=SsideW+(SsideW-U(i,3,k))
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
  elseif (BtypeS==NOSLIP.or.(component==2.and.BtypeS==FREESLIP)) then
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
  elseif (BtypeS==NEUMANN.or.(component/=2.and.BtypeS==FREESLIP)) then
    do k=-2,nz+3
     do i=-2,nx+3                       !Neumann inlet
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
    if (component==1.and.regime/=interm) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,ny+1,k)=NsideU+(NsideU-U(i,ny,k))
        U(i,ny+2,k)=NsideU+(NsideU-U(i,ny-1,k))
        U(i,ny+3,k)=NsideU+(NsideU-U(i,ny-2,k))
       enddo
      enddo
    elseif (component==2.and.regime/=interm) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,ny+1,k)=NsideV
        U(i,ny+2,k)=NsideV+(NsideV-U(i,ny,k))
        U(i,ny+3,k)=NsideV+(NsideV-U(i,ny-1,k))
       enddo
      enddo
    elseif (component==3.and.regime/=interm) then
      do k=-2,nz+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,ny+1,k)=NsideW+(NsideW-U(i,ny,k))
        U(i,ny+2,k)=NsideW+(NsideW-U(i,ny-1,k))
        U(i,ny+3,k)=NsideW+(NsideW-U(i,ny-2,k))
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
  elseif (BtypeN==NOSLIP.or.(component==2.and.BtypeS==FREESLIP)) then
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
  elseif (BtypeN==NEUMANN.or.(component/=2.and.BtypeN==FREESLIP)) then
    do k=-2,nz+3
     do i=-2,nx+3                       !Neumann inlet
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
    if (component==1.and.regime/=interm) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,0)=BsideU+(BsideU-U(i,j,1))
        U(i,j,-1)=BsideU+(BsideU-U(i,j,2))
        U(i,j,-2)=BsideU+(BsideU-U(i,j,3))
       enddo
      enddo
    elseif (component==2.and.regime/=interm) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,0)=BsideV+(BsideV-U(i,j,1))
        U(i,j,-1)=BsideV+(BsideV-U(i,j,2))
        U(i,j,-2)=BsideV+(BsideV-U(i,j,3))
       enddo
      enddo
    elseif (component==3.and.regime/=interm) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,0)=BsideW
        U(i,j,-1)=BsideW+(BsideW-U(i,j,1))
        U(i,j,-2)=BsideW+(BsideW-U(i,j,2))
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
  elseif (BtypeB==NOSLIP.or.(component==3.and.BtypeB==FREESLIP)) then
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
  elseif (BtypeB==NEUMANN.or.(component/=3.and.BtypeB==FREESLIP)) then
    do j=-2,ny+3
     do i=-2,nx+3                       !Neumann inlet
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
    if (component==1.and.regime/=interm) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,nz+1)=TsideU+(TsideU-U(i,j,nz))
        U(i,j,nz+2)=TsideU+(TsideU-U(i,j,nz-1))
        U(i,j,nz+3)=TsideU+(TsideU-U(i,j,nz-2))
       enddo
      enddo
    elseif (component==2.and.regime/=interm) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,nz+1)=TsideV+(TsideV-U(i,j,nz))
        U(i,j,nz+2)=TsideV+(TsideV-U(i,j,nz-1))
        U(i,j,nz+3)=TsideV+(TsideV-U(i,j,nz-2))
       enddo
      enddo
    elseif (component==3.and.regime/=interm) then
      do j=-2,ny+3
       do i=-2,nx+3                       !Dirichlet inlet
        U(i,j,nz+1)=TsideW
        U(i,j,nz+2)=TsideW+(TsideW-U(i,j,nz))
        U(i,j,nz+3)=TsideW+(TsideW-U(i,j,nz-1))
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
  elseif (BtypeT==NOSLIP.or.(component==3.and.(BtypeT==FREESLIP.or.BtypeT==FREESLIPBUFF))) then
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
  elseif (BtypeT==NEUMANN.or.(component/=3.and.(BtypeT==FREESLIP.or.BtypeT==FREESLIPBUFF))) then
    do j=-2,ny+3
     do i=-2,nx+3                       !Neumann inlet
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

  end subroutine BoundU_GPU



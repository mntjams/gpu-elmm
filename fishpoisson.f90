module FISHPOISSON

 use PARAMETERS

 contains

  subroutine POISSFISH(Phi,RHS)         !Calls direct solver from FISHPACK, only for uniform grid
                                        !for some boundary conditions it fails and so does fishpack90
  real(KND),dimension(0:,0:,0:),intent(out):: Phi
  real(KND),dimension(1:,1:,1:),intent(in)::RHS

  real(KND),dimension(1:Prny,1:Prnz):: XBC1,XBC2
  real(KND),dimension(1:Prnx,1:Prnz):: YBC1,YBC2
  real(KND),dimension(1:Prnx,1:Prny):: ZBC1,ZBC2

  integer nx,ny,nz,xb,yb,zb
  real(KND) PERTRB
  integer IERR
  real(KND) W(100+Prnx+Prny+6*Prnz+7*(INT((Prnx+2)/2) + INT((Prny+2)/2)))


   write (*,*) "Computing Poisson using hw3crt"
  if (BtypeE==PERIODIC) then
   nx=Prnx+1
  else
   nx=Prnx
  endif
  if (BtypeN==PERIODIC) then
   ny=Prny+1
  else
   ny=Prny
  endif
  if (BtypeT==PERIODIC) then
   nz=Prnz+1
  else
   nz=Prnz
  endif

   XBC1=0
   XBC2=0
   YBC1=0
   YBC2=0
   ZBC1=0
   ZBC2=0

   if (BtypeW==PERIODIC) then
                           xb=0
   else
    xb=3
   endif
   if (BtypeN==PERIODIC) then
                           yb=0
   else
    yb=3
   endif
   if (BtypeT==PERIODIC) then
                           zb=0
   else
    zb=3
   endif

   call hw3crt(dxmin,nx*dxmin,nx-1,xb,XBC1,XBC2,dymin,(ny)*dymin,&
                 ny-1,yb,YBC1,YBC2,dzmin,nz*dzmin,nz-1,zb,ZBC1,ZBC2,0,nx,ny,RHS,PERTRB,IERR,W)

    Phi(1:Prnx,1:Prny,1:Prnz)=RHS(1:Prnx,1:Prny,1:Prnz)

   if (ierr/=0) then
     write (*,*) "error ",IERR
     stop
   endif
   if (pertrb>1e-3) write (*,*) "pertrb",PERTRB

  end subroutine POISSFISH
end module FISHPOISSON


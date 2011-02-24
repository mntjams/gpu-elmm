module BLASIUS
use PARAMETERS
implicit none

  real(KND),dimension(:),allocatable:: Far,Gar
  integer nXI
  real(KND) dXI
contains
  
  real(KND) function Udvort(x,y,t)
  real(KND) x,y,t
    
   Udvort=-cos(3.1415926*x)*sin(3.1415926*y)*exp(-2.*3.1415926**2*t/Re)
  endfunction

  real(KND) function Vdvort(x,y,t)
  real(KND) x,y,t
    
   Vdvort=sin(3.1415926*x)*cos(3.1415926*y)*exp(-2.*3.1415926**2*t/Re)
  endfunction

  real(KND) function Prdvort(x,y,t)
  real(KND) x,y,t
    
   Prdvort=-(1./4.)*(cos(2*3.1415926*x)+cos(2*3.1415926*y))*exp(-4.*3.1415926**2*t/Re)
  endfunction

  

  real(KND) function Ublas(x,y)
  real(KND) x,y
   
   Ublas=Gblas(XIxy(x,y))*Uinlet
  endfunction Ublas


  real(KND) function Vblas(x,y)
  real(KND) x,y,XI
   XI=XIxy(x,y)
   Vblas=0.5*SQRT(Uinlet/(Re*x))*(XI*Gblas(XI)-FBlas(XI))
  endfunction Vblas

  
  real(KND) function XIxy(x,y)
  real(KND) x,y
   XIxy=y*Sqrt(Uinlet*Re/x)
  endfunction XIxy
  
  
  real(KND) function Fblas(XI)
  real(KND) XI
  integer i
  
  if (XI<nXi*dXI) then
   i=INT(XI/dXI)+1
   Fblas=(XI-dXI*(i-1))*(Far(I)-Far(i-1))/(dXI*(i)-dXI*(i-1))+Far(i-1)
  else
   Fblas=Far(nXI)+1.*(XI-nXI*dXI)
  endif
  
  endfunction Fblas

  real(KND) function Gblas(XI)
  real(KND) XI
  integer i
  if (XI<nXi*dXI) then
   i=INT(XI/dXI)+1
   Gblas=(XI-dXI*(i-1))*(Gar(I)-Gar(i-1))/(dXI*(i)-dXI*(i-1))+Gar(i-1)
  else
   Gblas=1.
  endif
  
  endfunction Gblas
  

  subroutine SOLVEF
  real(KND) F,G,H,F2,G2,H2,Flat,Glat,Hlat
  integer i
  
  F=0
  G=0
  H=0.332
  Far(0)=F
  Gar(0)=G
  F2=F + G*dXI
  G2=G + H*dXI
  H2=H - (F*H/2.)*dXI
  Far(1)=F2
  Gar(1)=G2

  do i=2,nXI
   Flat=F
   Glat=G
   Hlat=H
   F=F2
   G=G2
   H=H2
   F2=F + (3./2.)*(G*dXI) - (1./2.)*(Glat*dXI)
   G2=G + (3./2.)*(H*dXI) - (1./2.)*(Hlat*dXI)
   H2=H - (3./2.)*((F*H/2.)*dXI) + (1./2.)*((Flat*Hlat/2.)*dXI)
   Far(i)=F2
   Gar(i)=G2
  enddo
  
  if ((Gar(nXi) - 1)>0.001) write (*,*) "f'(inf)<>1"
  endsubroutine SOLVEF


  subroutine BLASINIT
  integer i
  nXI=1000
  dXI=12./nXI
  allocate(Far(0:nXI),Gar(0:nXI))

  call SOLVEF
  OPEN(11,file="Uxi.txt")
  do i=0,INT(1./dXI)
   write(11,*) i*dXI,Ublas(1._KND,i*dXI)
  enddo
  CLOSE(11)
  OPEN(11,file="Vxi.txt")
  do i=0,INT(1./dXI)
   write(11,*) i*dXI,Vblas(1._KND,i*dXI)!0.5*(i*dXI*Gblas(i*dXI)-Fblas(i*dXI))
  enddo
  CLOSE(11)
  write(*,*) "Solved analytical Blasius solution"
  endsubroutine BLASINIT
endmodule BLASIUS
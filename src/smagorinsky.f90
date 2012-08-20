module SMAGORINSKY

 use PARAMETERS
 use BOUNDARIES, only : BoundU

 implicit none

 private
 public :: Smag, Smag2, StabSmag, Vreman, Filter

 real(KND),parameter:: CSmag=0.1_KND

 contains

  subroutine Smag(U,V,W)  !Standard Smagorinsky model with implicit filtering
   integer i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W

   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
       Visc(i,j,k)=NuSmag(i,j,k,U,V,W)
     end do
    end do
   end do
   if (Re>0) then
     Visc=Visc+1._KND/Re
   end if
  endsubroutine Smag


  subroutine Smag2(U,V,W,width)  !Standard Smagorinsky model with explicit filtering
   integer,intent(in):: width
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   integer i,j,k

   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
      if (Re>0) then
       Visc(i,j,k)=1._KND/Re+NuSmag2(i,j,k,U,V,W,width*2+1)
      else
       Visc(i,j,k)=NuSmag2(i,j,k,U,V,W,width*2+1)
      end if
     end do
    end do
   end do
  endsubroutine Smag2




  real(KND) function NuSmag(i,j,k,U,V,W)    !subgrid viscosity for Smagorinsky with implicit filtering
   integer,intent(in):: i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND) S(1:3,1:3)
   real(KND) width,Sbar

   width=(dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND)
   call StrainIJ(i,j,k,U,V,W,S)
   Sbar=Strainu(S)
   NuSmag=Sbar*(width*CSmag)**2
  endfunction NuSmag


  real(KND) function NuSmag2(i,j,k,U,V,W,width2)    !subgrid viscosity for Smagorinsky with implicit filtering
   integer,intent(in):: i,j,k,width2
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND) S(1:3,1:3)
   real(KND) width,Sbar

   width=(dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND)
   call StrainIJ(i,j,k,U,V,W,S)
   Sbar=Strainu(S)
   NuSmag2=Sbar*((width2)*width*CSmag)**2
  endfunction NuSmag2




  subroutine StabSmag(U,V,W)  !Smagorinsky with a stability correction Brown et al. (1994)
   integer i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND) Ri,l,l0
   real(KND) width,Sbar
   real(KND),parameter:: CS=0.17_KND
   real(KND) S(1:3,1:3)

   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
       width=(dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND)

       call StrainIJ(i,j,k,U,V,W,S)
       Sbar=Strainu(S)

       Ri=Rig(i,j,k,U,V,temperature)


       l0=CS*width
       l=WallDamp(l0,z0B,zPr(k))

       Visc(i,j,k)=Sbar*Fm(Ri)*l
       TDiff(i,j,k)=Sbar*Fh(Ri)*l
     end do
    end do
   end do
   if (Re>0) then
     Visc=Visc+1._KND/Re
     TDiff=TDiff+1._KND/(Re*Prandtl)
   end if
  endsubroutine StabSmag


  pure real(KND) function Fm(Ri)       !Adjustment function for stability
   real(KND),intent(in):: Ri           !Pointwise Richarson number
   real(KND),parameter :: Ric=0.25_KND

   if (Ri>=Ric) then
    Fm  =  0
   elseif (Ri>0) then
    Fm  =  (1._KND-Ri/Ric)**4
   else
    Fm  =  sqrt(1._KND - 16._KND*Ri)
   end if
  end function Fm

  pure real(KND) function Fh(Ri)       !Adjustment function for stability
   real(KND),intent(in):: Ri           !Pointwise Richarson number
   real(KND),parameter :: Ric=0.25_KND

   if (Ri>=Ric) then
    Fh  =  0
   elseif (Ri>0) then
    Fh  =  (1./0.7_KND)  *  (1._KND - Ri/Ric)**4  *  (1._KND - 1.2_KND*Ri)
   else
    Fh  =  sqrt(1._KND - 40._KND*Ri) / 0.7_KND
   end if
  end function Fh

  pure real(KND) function WallDamp(l0,z0,z) !Wall damping for Smagorinsky model
   real(KND),intent(in):: l0,z0,z

   WallDamp= 1._KND / (&
                 1._KND/l0**2  +  1._KND/(0.4_KND*(z+z0))**2 &
                 )
  end function WallDamp


  pure real(KND) function Rig(i,j,k,U,V,temperature)
   integer,intent(in):: i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V
   real(KND),dimension(-1:,-1:,-1:),intent(in):: temperature
   real(KND) num,denom

   num=(grav_acc/temperature_ref)*(temperature(i,j,k+1)-temperature(i,j,k-1))/(zPr(k+1)-zPr(k-1))
   denom=((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._KND*(zPr(k+1)-zPr(k-1))))**2
   denom=denom+((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(zPr(k+1)-zPr(k-1))))**2
   if (abs(denom)>tiny(1._KND)*100) then
    Rig=num/denom
   else
    Rig=0
   end if
  endfunction Rig






  pure real(KND) function Strainu(S)      !Magnitude of the strain rate tensor.
   real(KND),intent(in):: S(1:3,1:3)
   integer::ii,jj

   Strainu=0
   do jj=1,3
    do ii=1,3
     Strainu=Strainu+S(ii,jj)*S(ii,jj)
    end do
   end do
   Strainu=SQRT(2._KND*Strainu)
  endfunction Strainu



  pure subroutine StrainIJ(i,j,k,U,V,W,S)     !Computes components of the strain rate tensor.
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND),intent(out):: S(1:3,1:3)
   integer,intent(in):: i,j,k
   real(KND) D(1:3,1:3)
   integer::ii,jj

   D=0

   if (i>-2.and.i<Prnx+3) D(1,1)=(U(i,j,k)-U(i-1,j,k))/dxPr(i)
   if (j>-2.and.j<Prny+3) D(2,2)=(V(i,j,k)-V(i,j-1,k))/dyPr(j)
   if (k>-2.and.j<Prnz+3) D(3,3)=(W(i,j,k)-W(i,j,k-1))/dzPr(k)
   if (j>-1.and.i>-1.and.j<Uny+3.and.i<Unx+4)&
      D(1,2)=(U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(2._KND*(yPr(j+1)-yPr(j-1)))
   if (k>-1.and.i>-1.and.k<Uny+3.and.i<Unx+4)&
      D(1,3)=(U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
   if (j>-1.and.i>-1.and.i<Vnx+3.and.j<Vny+4)&
      D(2,1)=(V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(2._KND*(xPr(i+1)-xPr(i-1)))
   if (j>-1.and.k>-1.and.k<Vnz+3.and.j<Vny+4)&
      D(2,3)=(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
   if (k>-1.and.i>-1.and.i<Wnx+3.and.k<Wnz+4)&
      D(3,1)=(W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(2._KND*(xPr(i+1)-xPr(i-1)))
   if (j>-1.and.k>-1.and.j<Wny+3.and.k<Wnz+4)&
      D(3,2)=(W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(2._KND*(yPr(j+1)-yPr(j-1)))

   do jj=1,3
    do ii=1,3
     S(ii,jj)=(D(ii,jj)+D(jj,ii))
    end do
   end do
   S=S/2._KND
  endsubroutine StrainIJ




  subroutine Vreman(U,V,W)   !Vreman subgrid model (Physics of Fluids, 2004)
   use Tiling, only: tilenx, tileny, tilenz
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   integer i,j,k,bi,bj,bk,ii,jj
   real(KND) :: aa,bb
   real(KND),dimension(1:3,1:3)::a,b
   real(KND) :: dx2,dy2,dz2
   real(KND),parameter::c=0.05
   integer,parameter :: narr = 4

   if (gridtype==uniformgrid) then


    dx2=dxmin**2
    dy2=dymin**2
    dz2=dzmin**2

    !$omp parallel do private(aa,bb,a,b,i,j,k,bi,bj,bk,ii,jj)
    do bk = 1,Prnz,tilenz(narr)
     do bj = 1,Prny,tileny(narr)
      do bi = 1,Prnx,tilenx(narr)
       do k = bk,min(bk+tilenz(narr)-1,Prnz)
        do j = bj,min(bj+tileny(narr)-1,Prny)
         do i = bi,min(bi+tilenx(narr)-1,Prnx)
          a(1,1)=(U(i,j,k)-U(i-1,j,k))/dxmin
          a(2,1)=(U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(4._KND*dymin)
          a(3,1)=(U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(4._KND*dzmin)

          a(2,2)=(V(i,j,k)-V(i,j-1,k))/dymin
          a(1,2)=(V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(4._KND*dxmin)
          a(3,2)=(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(4._KND*dzmin)

          a(3,3)=(W(i,j,k)-W(i,j,k-1))/dzmin
          a(1,3)=(W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(4._KND*dxmin)
          a(2,3)=(W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(4._KND*dymin)

          do jj=1,3
           do ii=1,3
            b(ii,jj)=dx2*a(1,ii)*a(1,jj)+&
                     dy2*a(2,ii)*a(2,jj)+&
                     dz2*a(3,ii)*a(3,jj)
           end do
          end do

          bb=          b(1,1)*b(2,2)-b(1,2)**2
          bb=bb+b(1,1)*b(3,3)-b(1,3)**2
          bb=bb+b(2,2)*b(3,3)-b(2,3)**2

          aa=0

          do jj=1,3
           do ii=1,3
            aa=aa+(a(ii,jj)**2)
           end do
          end do


          if (abs(aa)>1e-5.and.bb>0)  then
            Visc(i,j,k)=c*sqrt(bb/aa)
          else
            Visc(i,j,k) = 0
          end if

          if (Re>0)  Visc(i,j,k)=Visc(i,j,k)+1._KND/Re

         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end parallel do

   else !general grid

    !$omp parallel do private(aa,bb,a,b,i,j,k,bi,bj,bk,ii,jj)
    do bk = 1,Prnz,tilenz(narr)
     do bj = 1,Prny,tileny(narr)
      do bi = 1,Prnx,tilenx(narr)
       do k = bk,min(bk+tilenz(narr)-1,Prnz)
        do j = bj,min(bj+tileny(narr)-1,Prny)
         do i = bi,min(bi+tilenx(narr)-1,Prnx)
          a(1,1)=(U(i,j,k)-U(i-1,j,k))/dxPr(i)
          a(2,1)=(U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(2._KND*(dyV(j)+dyV(j-1)))
          a(3,1)=(U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._KND*(dzW(k)+dzW(k-1)))
          a(2,2)=(V(i,j,k)-V(i,j-1,k))/dyPr(j)
          a(1,2)=(V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(2._KND*(dxU(i)+dxU(i-1)))
          a(3,2)=(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(dzW(k)+dzW(k-1)))
          a(3,3)=(W(i,j,k)-W(i,j,k-1))/dzPr(k)
          a(1,3)=(W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(2._KND*(dxU(i)+dxU(i-1)))
          a(2,3)=(W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(2._KND*(dyV(j)+dyV(j-1)))


          forall(jj=1:3,ii=1:3)
           b(ii,jj)=(dxPr(i))**2*a(1,ii)*a(1,jj)+&
                          (dyPr(j))**2*a(2,ii)*a(2,jj)+&
                          (dzPr(k))**2*a(3,ii)*a(3,jj)
          endforall

          bb=   b(1,1)*b(2,2)-b(1,2)**2
          bb=bb+b(1,1)*b(3,3)-b(1,3)**2
          bb=bb+b(2,2)*b(3,3)-b(2,3)**2

          Visc(i,j,k)=0._KND


          aa=sum(a(:,:)**2)

          if (abs(aa)>1e-5.and.bb>0) Visc(i,j,k)=c*sqrt(bb/aa)

          if (Re>0)  Visc(i,j,k)=Visc(i,j,k)+1._KND/Re

         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end parallel do

   end if !general grid
  endsubroutine Vreman













 subroutine TrapesField(U1,U2,nx,ny,nz)   !Trapesoidal filter
  real(KND),dimension(-2:,-2:,-2:),intent(in):: U1
  real(KND),dimension(-2:,-2:,-2:),intent(out):: U2
  integer,intent(in):: nx,ny,nz
  integer i,j,k

  U2=0
  do k=-1,nz+2
   do j=-1,ny+2
    do i=-1,nx+2
     U2(i,j,k)=(U1(i-1,j,k)+U1(i+1,j,k)+U1(i,j-1,k)+U1(i,j-1,k)+U1(i,j,k-1)+U1(i,j,k+1)+6*U1(i,j,k))/12._KND
    end do
   end do
  end do

 endsubroutine TrapesField



 subroutine TrapesField2(U1,U2,nx,ny,nz) !Filter with weighting
  real(KND),dimension(-2:,-2:,-2:),intent(in):: U1
  real(KND),dimension(-2:,-2:,-2:),intent(out):: U2
  integer,intent(in):: nx,ny,nz
  real(KND),dimension(-1:1):: F
  real(KND) S,S1
  integer ii,jj,i,j,k,kk

  F=(/ 1,2,1 /)
  U2=0
  S1=0
  do kk=-1,1
   do jj=-1,1
    do ii=-1,1
     S1=S1+F(ii)*F(jj)*F(kk)
    end do
   end do
  end do
  do k=-1,nz+2
   do j=-1,ny+2
    do i=-1,nx+2
     S=0
     do kk=-1,1
      do jj=-1,1
       do ii=-1,1
        S=S+F(ii)*F(jj)*F(kk)*U1(i+ii,j+jj,k+kk)
       end do
      end do
     end do
     U2(i,j,k)=S/S1
    end do
   end do
  end do
 endsubroutine TrapesField2



 subroutine Filter(U1,U2,nx,ny,nz,width)    !Calls a selected filter
  real(KND),dimension(-2:,-2:,-2:),intent(in):: U1
  real(KND),dimension(-2:,-2:,-2:),intent(out):: U2
  integer,intent(in):: nx,ny,nz
  integer width

  call TrapesField(U1,U2,nx,ny,nz)

 endsubroutine Filter

endmodule SMAGORINSKY

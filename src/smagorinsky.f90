module SMAGORINSKY

 use PARAMETERS
 use BOUNDARIES, only : BoundU

 implicit none

 private
 public :: SGS_Smag, Smag2, SGS_StabSmag, SGS_Vreman, SGS_Sigma, SGS_Filter

 real(KND),parameter :: CSmag = 0.1_KND

 interface Filter3ord
   module procedure Filter3ord_2D
   module procedure Filter3ord_3D
 end interface

 contains

  subroutine SGS_Smag(U,V,W,filter_ratio)  !Standard Smagorinsky model with implicit filtering
     real(KND),intent(in) :: filter_ratio
  integer i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W

   do k = 1,Prnz
    do j = 1,Prny
     do i = 1,Prnx
       Visc(i,j,k) = NuSmag(i,j,k,U,V,W,filter_ratio)
     end do
    end do
   end do
   if (Re>0) then
     Visc = Visc+1._KND/Re
   end if
  endsubroutine SGS_Smag


  subroutine Smag2(U,V,W,filter_ratio)  !Standard Smagorinsky model with explicit filtering
   real(KND),intent(in) :: filter_ratio
   real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
   integer i,j,k

   do k = 1,Prnz
    do j = 1,Prny
     do i = 1,Prnx
      if (Re>0) then
       Visc(i,j,k) = 1._KND/Re+NuSmag2(i,j,k,U,V,W,filter_ratio)
      else
       Visc(i,j,k) = NuSmag2(i,j,k,U,V,W,filter_ratio)
      end if
     end do
    end do
   end do
  endsubroutine Smag2




  real(KND) function NuSmag(i,j,k,U,V,W,filter_ratio)    !subgrid viscosity for Smagorinsky with implicit filtering
   integer,intent(in) :: i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
   real(KND),intent(in) :: filter_ratio
   real(KND) S(1:3,1:3)
   real(KND) width,Sbar

   width = (dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND)
   call StrainIJ(i,j,k,U,V,W,S)
   Sbar = Strainu(S)
   NuSmag = Sbar*(width*CSmag)**2
  endfunction NuSmag


  real(KND) function NuSmag2(i,j,k,U,V,W,filter_ratio)    !subgrid viscosity for Smagorinsky with implicit filtering
   integer,intent(in) :: i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
   real(KND),intent(in) :: filter_ratio
   real(KND) S(1:3,1:3)
   real(KND) width,Sbar

   width = filter_ratio * (dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND)
   call StrainIJ(i,j,k,U,V,W,S)
   Sbar = Strainu(S)
   NuSmag2 = Sbar*(width*CSmag)**2
  endfunction NuSmag2




  subroutine SGS_StabSmag(U,V,W,filter_ratio)  !Smagorinsky with a stability correction Brown et al. (1994)
    real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    real(KND),intent(in) :: filter_ratio
    real(KND) Ri,l,l0
    real(KND) width,Sbar
    real(KND),parameter :: CS = 0.17_KND
    real(KND) S(1:3,1:3)
    integer i,j,k

    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
        width = filter_ratio * (dxPr(i)*dyPr(j)*dzPr(k))**(1._KND/3._KND)

        call StrainIJ(i,j,k,U,V,W,S)
        Sbar = Strainu(S)

        Ri = Rig(i,j,k,U,V,temperature)


        l0 = CS*width
        l = WallDamp(l0,z0B,zPr(k))

        Visc(i,j,k) = Sbar*Fm(Ri)*l
        TDiff(i,j,k) = Sbar*Fh(Ri)*l
      end do
     end do
    end do
    if (Re>0) then
      Visc = Visc+1._KND/Re
      TDiff = TDiff+1._KND/(Re*Prandtl)
    end if
  endsubroutine SGS_StabSmag


  pure real(KND) function Fm(Ri)       !Adjustment function for stability
   real(KND),intent(in) :: Ri           !Pointwise Richarson number
   real(KND),parameter :: Ric = 0.25_KND

   if (Ri>=Ric) then
    Fm  =  0
   elseif (Ri>0) then
    Fm  =  (1._KND-Ri/Ric)**4
   else
    Fm  =  sqrt(1._KND - 16._KND*Ri)
   end if
  end function Fm

  pure real(KND) function Fh(Ri)       !Adjustment function for stability
   real(KND),intent(in) :: Ri           !Pointwise Richarson number
   real(KND),parameter :: Ric = 0.25_KND

   if (Ri>=Ric) then
    Fh  =  0
   elseif (Ri>0) then
    Fh  =  (1./0.7_KND)  *  (1._KND - Ri/Ric)**4  *  (1._KND - 1.2_KND*Ri)
   else
    Fh  =  sqrt(1._KND - 40._KND*Ri) / 0.7_KND
   end if
  end function Fh

  pure real(KND) function WallDamp(l0,z0,z) !Wall damping for Smagorinsky model
   real(KND),intent(in) :: l0,z0,z

   WallDamp= 1._KND / (&
                 1._KND/l0**2  +  1._KND/(0.4_KND*(z+z0))**2 &
                 )
  end function WallDamp


  pure real(KND) function Rig(i,j,k,U,V,temperature)
   integer,intent(in) :: i,j,k
   real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V
   real(KND),dimension(-1:,-1:,-1:),intent(in) :: temperature
   real(KND) num,denom

   num = (grav_acc/temperature_ref)*(temperature(i,j,k+1)-temperature(i,j,k-1))/(zPr(k+1)-zPr(k-1))
   denom = ((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._KND*(zPr(k+1)-zPr(k-1))))**2
   denom = denom+((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(zPr(k+1)-zPr(k-1))))**2
   if (abs(denom)>tiny(1._KND)*100) then
    Rig = num/denom
   else
    Rig = 0
   end if
  endfunction Rig






  pure real(KND) function Strainu(S)      !Magnitude of the strain rate tensor.
   real(KND),intent(in) :: S(1:3,1:3)
   integer ::ii,jj

   Strainu = 0
   do jj = 1,3
    do ii = 1,3
     Strainu = Strainu+S(ii,jj)*S(ii,jj)
    end do
   end do
   Strainu = SQRT(2._KND*Strainu)
  endfunction Strainu



  pure subroutine StrainIJ(i,j,k,U,V,W,S)     !Computes components of the strain rate tensor.
   real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
   real(KND),intent(out) :: S(1:3,1:3)
   integer,intent(in) :: i,j,k
   real(KND) D(1:3,1:3)
   integer ::ii,jj

   D = 0

   if (i>-2.and.i<Prnx+3) D(1,1) = (U(i,j,k)-U(i-1,j,k))/dxPr(i)
   if (j>-2.and.j<Prny+3) D(2,2) = (V(i,j,k)-V(i,j-1,k))/dyPr(j)
   if (k>-2.and.j<Prnz+3) D(3,3) = (W(i,j,k)-W(i,j,k-1))/dzPr(k)
   if (j>-1.and.i>-1.and.j<Uny+3.and.i<Unx+4)&
      D(1,2) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(2._KND*(yPr(j+1)-yPr(j-1)))
   if (k>-1.and.i>-1.and.k<Uny+3.and.i<Unx+4)&
      D(1,3) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
   if (j>-1.and.i>-1.and.i<Vnx+3.and.j<Vny+4)&
      D(2,1) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(2._KND*(xPr(i+1)-xPr(i-1)))
   if (j>-1.and.k>-1.and.k<Vnz+3.and.j<Vny+4)&
      D(2,3) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
   if (k>-1.and.i>-1.and.i<Wnx+3.and.k<Wnz+4)&
      D(3,1) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(2._KND*(xPr(i+1)-xPr(i-1)))
   if (j>-1.and.k>-1.and.j<Wny+3.and.k<Wnz+4)&
      D(3,2) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(2._KND*(yPr(j+1)-yPr(j-1)))

   do jj = 1,3
    do ii = 1,3
     S(ii,jj) = (D(ii,jj)+D(jj,ii))
    end do
   end do
   S = S/2._KND
  endsubroutine StrainIJ




  subroutine SGS_Vreman(U,V,W,filter_ratio)   !Vreman subgrid model (Physics of Fluids, 2004)
    use Tiling, only: tilenx, tileny, tilenz
    real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    real(KND),intent(in) :: filter_ratio
    integer i,j,k,bi,bj,bk,ii,jj
    real(KND) :: aa,bb
    real(KND),dimension(1:3,1:3) ::a,b
    real(KND) :: dx2,dy2,dz2
    real(KND),parameter ::c = 0.05
    integer,parameter :: narr = 4
    real(KND) :: c2

    c2 = c * filter_ratio**2

    if (gridtype==uniformgrid) then


     dx2 = dxmin**2
     dy2 = dymin**2
     dz2 = dzmin**2

     !$omp parallel do private(aa,bb,a,b,i,j,k,bi,bj,bk,ii,jj)
     do bk = 1,Prnz,tilenz(narr)
      do bj = 1,Prny,tileny(narr)
       do bi = 1,Prnx,tilenx(narr)
        do k = bk,min(bk+tilenz(narr)-1,Prnz)
         do j = bj,min(bj+tileny(narr)-1,Prny)
          do i = bi,min(bi+tilenx(narr)-1,Prnx)
           a(1,1) = (U(i,j,k)-U(i-1,j,k))/dxmin
           a(2,1) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(4._KND*dymin)
           a(3,1) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(4._KND*dzmin)

           a(2,2) = (V(i,j,k)-V(i,j-1,k))/dymin
           a(1,2) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(4._KND*dxmin)
           a(3,2) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(4._KND*dzmin)

           a(3,3) = (W(i,j,k)-W(i,j,k-1))/dzmin
           a(1,3) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(4._KND*dxmin)
           a(2,3) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(4._KND*dymin)

           do jj = 1,3
            do ii = 1,3
             b(ii,jj) = dx2*a(1,ii)*a(1,jj)+&
                      dy2*a(2,ii)*a(2,jj)+&
                      dz2*a(3,ii)*a(3,jj)
            end do
           end do

           bb=          b(1,1)*b(2,2)-b(1,2)**2
           bb = bb+b(1,1)*b(3,3)-b(1,3)**2
           bb = bb+b(2,2)*b(3,3)-b(2,3)**2

           aa = 0

           do jj = 1,3
            do ii = 1,3
             aa = aa+(a(ii,jj)**2)
            end do
           end do


           if (abs(aa)>1e-5.and.bb>0)  then
             Visc(i,j,k) = c2 * sqrt(bb/aa)
           else
             Visc(i,j,k) = 0
           end if

           if (Re>0)  Visc(i,j,k) = Visc(i,j,k)+1._KND/Re

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
           a(1,1) = (U(i,j,k)-U(i-1,j,k))/dxPr(i)
           a(2,1) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(2._KND*(dyV(j)+dyV(j-1)))
           a(3,1) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._KND*(dzW(k)+dzW(k-1)))
           a(2,2) = (V(i,j,k)-V(i,j-1,k))/dyPr(j)
           a(1,2) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(2._KND*(dxU(i)+dxU(i-1)))
           a(3,2) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(dzW(k)+dzW(k-1)))
           a(3,3) = (W(i,j,k)-W(i,j,k-1))/dzPr(k)
           a(1,3) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(2._KND*(dxU(i)+dxU(i-1)))
           a(2,3) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(2._KND*(dyV(j)+dyV(j-1)))


           forall(jj = 1:3,ii = 1:3)
            b(ii,jj) = (dxPr(i))**2*a(1,ii)*a(1,jj)+&
                           (dyPr(j))**2*a(2,ii)*a(2,jj)+&
                           (dzPr(k))**2*a(3,ii)*a(3,jj)
           endforall

           bb=   b(1,1)*b(2,2)-b(1,2)**2
           bb = bb+b(1,1)*b(3,3)-b(1,3)**2
           bb = bb+b(2,2)*b(3,3)-b(2,3)**2

           Visc(i,j,k) = 0._KND


           aa = sum(a(:,:)**2)

           if (abs(aa)>1e-5.and.bb>0) Visc(i,j,k) = c2 * sqrt(bb/aa)

           if (Re>0)  Visc(i,j,k) = Visc(i,j,k)+1._KND/Re

          end do
         end do
        end do
       end do
      end do
     end do
     !$omp end parallel do

    end if !general grid
  endsubroutine SGS_Vreman





  subroutine SGS_Sigma(U,V,W,filter_ratio)
    !from Nicoud, Toda, Cabrit, Bose, Lee, http://dx.doi.org/10.1063/1.3623274
    use Tiling, only: tilenx, tileny, tilenz
    real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    real(KND),intent(in) :: filter_ratio
    real(KND),parameter :: Csig = 1.35_KND
    integer,parameter   :: narr = 4
    integer   :: i,j,k,bi,bj,bk
    real(KND) :: width, C, D, g(3,3), s1, s2, s3

    width = filter_ratio * (dxmin*dymin*dzmin)**(1._KND/3._KND)
    C = (Csig*width)**2

    !$omp parallel do private(g,s1,s2,s3,D,i,j,k,bi,bj,bk)
    do bk = 1,Prnz,tilenz(narr)
     do bj = 1,Prny,tileny(narr)
      do bi = 1,Prnx,tilenx(narr)
       do k = bk,min(bk+tilenz(narr)-1,Prnz)
        do j = bj,min(bj+tileny(narr)-1,Prny)
         do i = bi,min(bi+tilenx(narr)-1,Prnx)

          call GradientTensorUG(g,i,j,k,U,V,W)

          call Sigmas(s1,s2,s3,g)

          if (s1>0) then
            D = (s3 * (s1 - s2) * (s2 - s3)) / s1**2
          else
            D = 0
          end if

          Visc(i,j,k) = C * D

          if (Re>0)  Visc(i,j,k) = Visc(i,j,k)+1._KND/Re

         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end parallel do

  end subroutine SGS_Sigma


  pure subroutine GradientTensorUG(g,i,j,k,U,V,W)
    real(KND),intent(out) :: g(3,3)
    integer,intent(in) :: i,j,k
    real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W

    g(1,1) = (U(i,j,k)-U(i-1,j,k))/dxmin
    g(2,1) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(4._KND*dymin)
    g(3,1) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(4._KND*dzmin)

    g(2,2) = (V(i,j,k)-V(i,j-1,k))/dymin
    g(1,2) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(4._KND*dxmin)
    g(3,2) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(4._KND*dzmin)

    g(3,3) = (W(i,j,k)-W(i,j,k-1))/dzmin
    g(1,3) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(4._KND*dxmin)
    g(2,3) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(4._KND*dymin)
  end subroutine GradientTensorUG


  pure subroutine Sigmas(s1,s2,s3,grads)
    !from Hasan, Basser, Parker, Alexander, http://dx.doi.org/10.1006/jmre.2001.2400
    !via Nicoud, Toda, Cabrit, Bose, Lee, http://dx.doi.org/10.1063/1.3623274
    real(KND),intent(out) :: s1,s2,s3
    real(KND),intent(in)  :: grads(3,3)

    real(DBL) :: trG2, i1, i2, i3, a1, a2, a3, c, G(3,3)

    G = matmul(transpose(grads),grads)

    trG2 = dot_product(G(:,1),G(:,1)) +&
           dot_product(G(:,2),G(:,2)) +&
           dot_product(G(:,3),G(:,3))

    i1 = G(1,1) + G(2,2) + G(3,3)

    i2 = (i1**2 - trG2)
    i2 = i2 / 2

    i3 = det3x3(G)

    a1 = max((i1**2)/9 - i2/3,0._DBL)

    a2 = (i1**3)/27 - i1*i2/6 +i3/2

    if (a1>10*epsilon(1._DBL)) then
       c =  a2 / sqrt(a1**3)
    else
       c = 0
    end if
    c = max(-1._DBL,min(1._DBL,c))
    a3 = acos(c) / 3

    c = 2*sqrt(a1)

    s1 = real( sqrt( i1/3 + c*cos(a3) ) , KND )

    s2 = real( sqrt( max( i1/3 - c*cos(pi/3 + a3) , 0._DBL ) ) , KND)

    s3 = real( sqrt( max( i1/3 - c*cos(pi/3 - a3) , 0._DBL ) ) , KND)

  end subroutine Sigmas


  pure function det3x3(A) result (res)
    real(DBL),intent(in) :: A(3,3)
    real(DBL) :: res

    res =   A(1,1)*A(2,2)*A(3,3)  &
          - A(1,1)*A(2,3)*A(3,2)  &
          - A(1,2)*A(2,1)*A(3,3)  &
          + A(1,2)*A(2,3)*A(3,1)  &
          + A(1,3)*A(2,1)*A(3,2)  &
          - A(1,3)*A(2,2)*A(3,1)
  end function det3x3


 subroutine TrapesField(U2,U1)   !Trapesoidal filter
  real(KND),dimension(-2:,-2:,-2:),intent(in) :: U1
  real(KND),dimension(-2:,-2:,-2:),intent(out) :: U2
  integer i,j,k

  U2 = 0
  do k = lbound(U1,3)+1,ubound(U1,3)-1
   do j = lbound(U1,2)+1,ubound(U1,2)-1
    do i = lbound(U1,1)+1,ubound(U1,1)-1
     U2(i,j,k) = (U1(i-1,j,k)+U1(i+1,j,k)+U1(i,j-1,k)+U1(i,j-1,k)+U1(i,j,k-1)+U1(i,j,k+1)+6*U1(i,j,k))/12._KND
    end do
   end do
  end do

 endsubroutine TrapesField



 subroutine TrapesField2(U2,U1) !Filter with weighting
  real(KND),dimension(:,:,:),intent(in) :: U1
  real(KND),dimension(:,:,:),intent(out) :: U2
  real(KND),dimension(-1:1) :: F
  real(KND) S,S1
  integer ii,jj,i,j,k,kk

  F = (/ 1,2,1 /)
  U2 = 0
  S1 = 0
  do kk=-1,1
!    do jj=-1,1
    do ii=-1,1
     S1 = S1+F(ii)*F(kk)!+jj
    end do
!    end do
  end do
  do k = lbound(U1,3)+1,ubound(U1,3)-1
   do j = lbound(U1,2),ubound(U1,2)!lbound(U1,2)+1,ubound(U1,2)-1
    do i = lbound(U1,1)+1,ubound(U1,1)-1
     S = 0
     do kk=-1,1
!       do jj=-1,1
       do ii=-1,1
        S = S+F(ii)*F(kk)*U1(i+ii,j,k+kk) !*F(jj) +jj
       end do
!       end do
     end do
     U2(i,j,k) = S/S1
    end do
   end do
  end do
 endsubroutine TrapesField2


 subroutine Filter3ord_3D(U2,U1)
  !3D filter with 3 vanishing moments width 2 delta x
  !Vasilyev, Lund, Moin, 1998, http://dx.doi.org/10.1006/jcph.1998.6060
  real(KND),dimension(:,:,:),intent(in)  :: U1
  real(KND),dimension(:,:,:),intent(out) :: U2
  real(KND),parameter :: w(-3:3) = (/ -1._KND/32, 0._KND, 9._KND/32, 1._KND/2, 9._KND/32, 0._KND, -1._KND/32 /)
  real(KND),save :: w3(-3:3,-3:3,-3:3)
  real(KND) S
  integer ii,jj,i,j,k,kk
  integer,save :: called = 0

  if (called ==0) then
    do kk=-3,3
     do jj=-3,3
      do ii=-3,3
       w3(ii,jj,kk) = w(ii)*w(jj)*w(kk)
      end do
     end do
    end do
    called = 1
  end if

  !$omp parallel do private(i,j,k,ii,jj,kk,S)
  do k = lbound(U1,3)+3,ubound(U1,3)-3
   do j = lbound(U1,2)+3,ubound(U1,2)-3
    do i = lbound(U1,1)+3,ubound(U1,1)-3
     S = 0
     do kk = -3,3
      do jj = -3,3
       do ii = -3,3
        S = S + W3(ii,jj,kk) * U1(i+ii,j+jj,k+kk)
       end do
      end do
     end do
     U2(i,j,k) = S
    end do
   end do
  end do
  !$omp end parallel do
 endsubroutine Filter3ord_3D


 subroutine Filter3ord_2D(U2,U1,dir)
  !3D filter with 3 vanishing moments width 2 delta x
  !Vasilyev, Lund, Moin, 1998, http://dx.doi.org/10.1006/jcph.1998.6060
  real(KND),dimension(-2:,-2:,-2:),intent(in)  :: U1
  real(KND),dimension(-2:,-2:,-2:),intent(out) :: U2
  integer,intent(in) :: dir
  real(KND),parameter :: w(-3:3) = (/ -1._KND/32, 0._KND, 9._KND/32, 1._KND/2, 9._KND/32, 0._KND, -1._KND/32 /)
  real(KND) :: w2(-3:3,-3:3,-3:3)
  real(KND) S
  integer ii,jj,i,j,k,kk
  integer mini,minii,minj,minjj,mink,minkk
  integer maxi,maxii,maxj,maxjj,maxk,maxkk

  minii = -3
  maxii = 3
  minjj = -3
  maxjj = 3
  minkk = -3
  maxkk = 3

  mini = lbound(U1,1)+3
  maxi = ubound(U1,1)-3
  minj = lbound(U1,2)+3
  maxj = ubound(U1,2)-3
  mink = lbound(U1,3)+3
  maxk = ubound(U1,3)-3

  if (dir==1) then
   minii=0
   maxii=0
   mini=1
   maxi=1
  else if (dir==2) then
   minjj=0
   maxjj=0
   minj=1
   maxj=1
  else
   minkk=0
   maxkk=0
   mink=1
   maxk=1
  end if

  do kk = minkk,maxkk
   do jj = minjj,maxjj
    do ii = minii,maxii
     w2(ii,jj,kk) = 2 * w(ii)*w(jj)*w(kk)
    end do
   end do
  end do

  !$omp parallel do  private(i,j,k,ii,jj,kk,S)
  do k = mink,maxk
   do j = minj,maxj
    do i = mini,maxi
     S = 0
     do kk = minkk,maxkk
      do jj = minjj,maxjj
       do ii = minii,maxii
        S = S + W2(ii,jj,kk) * U1(i+ii, j+jj, k+kk)
!         print *,U1(i+ii, j+jj, k+kk),W2(ii,jj,kk) * U1(i+ii, j+jj, k+kk)
       end do
      end do
     end do
     U2(i,j,k) = S
!      print *,i,j,k,'-----',S,'-----'
    end do
   end do
  end do
  !$omp end parallel do
 endsubroutine Filter3ord_2D


  subroutine SpectralFilter(U2,U1,kfilt)
    use PoisFFT
    type(PoisFFT_Solver3D) :: Solver
    real(KND),dimension(:,:,:),intent(in) :: U1
    real(KND),dimension(:,:,:),intent(out) :: U2
    integer,intent(in) :: kfilt
    integer,parameter :: bounds(6) = PoisFFT_Periodic
    logical,save :: called = .false.

    Solver = PoisFFT_Solver3D_New(size(U1,1),size(U1,2),size(U1,3),dxmin,dymin,dzmin,bounds)
    called = .true.

    print *, U1(1,1,1),U1(10,10,10),U1(20,20,20)
    Solver%cwork = cmplx(U1,0._KND,KND)

    call Execute(Solver%forward, Solver%cwork)

!     Solver%cwork(ubound(Solver%cwork,1)-kfilt:ubound(Solver%cwork,1),:,:) = 0
!     Solver%cwork(:,ubound(Solver%cwork,2)-kfilt:ubound(Solver%cwork,2),:) = 0
!     Solver%cwork(:,:,ubound(Solver%cwork,3)-kfilt:ubound(Solver%cwork,3)) = 0

!     Solver%cwork(ubound(Solver%cwork,1)-kfilt:ubound(Solver%cwork,1),&
!                  ubound(Solver%cwork,2)-kfilt:ubound(Solver%cwork,2),&
!                  ubound(Solver%cwork,3)-kfilt:ubound(Solver%cwork,3)) = 0

    call Execute(Solver%backward, Solver%cwork)

    U2 = real(Solver%cwork,KND) / (size(U1,1)*size(U1,2)*size(U1,3))
print *, U2(1,1,1), U2(10,10,10), U2(20,20,20)
    call DeallocateData(Solver)
  end subroutine SpectralFilter

 subroutine SGS_Filter(U2,U1)    !Calls a selected filter  U2 <- filt(U1)
  real(KND),dimension(:,:,:),intent(in) :: U1
  real(KND),dimension(:,:,:),intent(out) :: U2

  if (Prnx==1) then
    call Filter3ord(U2,U1,1)
  else if (Prny==1) then
    call Filter3ord(U2,U1,2)
  else if (Prnz==1) then
    call Filter3ord(U2,U1,3)
  else
    call Filter3ord(U2,U1)
  end if

 endsubroutine SGS_Filter

endmodule SMAGORINSKY

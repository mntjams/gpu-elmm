module Subgrid

  use Parameters
  use Boundaries, only: BoundU

  implicit none

  private
  public :: SGS_Smag, SGS_StabSmag, SGS_Vreman, SGS_Sigma, sgstype

  real(knd),parameter :: CSmag = 0.1_knd

  integer :: sgstype

  contains

    subroutine SGS_Smag(U,V,W,filter_ratio)  !Standard Smagorinsky model with implicit filtering
     use ArrayUtilities,only: add
     real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
     real(knd),intent(in) :: filter_ratio
     integer i,j,k

      !$omp parallel do private(i,j,k)
      do k = 1,Prnz
       do j = 1,Prny
        do i = 1,Prnx
          Viscosity(i,j,k) = NuSmag(i,j,k,U,V,W,filter_ratio)
        end do
       end do
      end do
      !$omp end parallel do
      if (Re>0) then
        call add(Viscosity,1._knd/Re)
      end if

    endsubroutine SGS_Smag




    real(knd) function NuSmag(i,j,k,U,V,W,filter_ratio)    !subgrid viscosity for Smagorinsky with implicit filtering
      integer,intent(in) :: i,j,k
      real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
      real(knd),intent(in) :: filter_ratio
      real(knd) S(1:3,1:3)
      real(knd) width,Sbar

      width = filter_ratio * (dxPr(i)*dyPr(j)*dzPr(k))**(1._knd/3._knd)
      call StrainIJ(i,j,k,U,V,W,S)
      Sbar = Strainu(S)
      NuSmag = Sbar*(width*CSmag)**2
    endfunction NuSmag





    subroutine SGS_StabSmag(U,V,W,Temperature,filter_ratio)  !Smagorinsky with a stability correction Brown et al. (1994)
      real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
      real(knd),dimension(-1:,-1:,-1:),intent(in) :: Temperature
      real(knd),intent(in) :: filter_ratio
      real(knd) Ri,l,l0
      real(knd) width,Sbar
      real(knd),parameter :: CS = 0.17_knd
      real(knd) S(1:3,1:3)
      integer i,j,k

      do k = 1,Prnz
       do j = 1,Prny
        do i = 1,Prnx
          width = filter_ratio * (dxPr(i)*dyPr(j)*dzPr(k))**(1._knd/3._knd)

          call StrainIJ(i,j,k,U,V,W,S)
          Sbar = Strainu(S)

          Ri = Rig(i,j,k,U,V,temperature)


          l0 = CS*width
          l = WallDamp(l0,z0B,zPr(k))

          Viscosity(i,j,k)  = Sbar*Fm(Ri)*l
          TDiff(i,j,k) = Sbar*Fh(Ri)*l
        end do
       end do
      end do
      if (Re>0) then
        Viscosity  = Viscosity  + 1._knd/Re
        TDiff = TDiff + 1._knd/(Re*Prandtl)
      end if
    endsubroutine SGS_StabSmag


    pure real(knd) function Fm(Ri)       !Adjustment function for stability
      real(knd),intent(in) :: Ri           !Pointwise Richarson number
      real(knd),parameter :: Ric = 0.25_knd

      if (Ri>=Ric) then
       Fm  =  0
      elseif (Ri>0) then
       Fm  =  (1._knd-Ri/Ric)**4
      else
       Fm  =  sqrt(1._knd - 16._knd*Ri)
      end if
    end function Fm

    pure real(knd) function Fh(Ri)       !Adjustment function for stability
      real(knd),intent(in) :: Ri           !Pointwise Richarson number
      real(knd),parameter :: Ric = 0.25_knd

      if (Ri>=Ric) then
       Fh  =  0
      elseif (Ri>0) then
       Fh  =  (1./0.7_knd)  *  (1._knd - Ri/Ric)**4  *  (1._knd - 1.2_knd*Ri)
      else
       Fh  =  sqrt(1._knd - 40._knd*Ri) / 0.7_knd
      end if
    end function Fh

    pure real(knd) function WallDamp(l0,z0,z) !Wall damping for Smagorinsky model
      real(knd),intent(in) :: l0,z0,z

      WallDamp= 1._knd / (&
                    1._knd/l0**2  +  1._knd/(0.4_knd*(z+z0))**2 &
                    )
    end function WallDamp


    pure real(knd) function Rig(i,j,k,U,V,temperature)
      integer,intent(in) :: i,j,k
      real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V
      real(knd),dimension(-1:,-1:,-1:),intent(in) :: temperature
      real(knd) num,denom

      num = (grav_acc/temperature_ref) * &
            (temperature(i,j,k+1)-temperature(i,j,k-1)) / (zPr(k+1)-zPr(k-1))
      
      denom = ((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1)) / (2._knd*(zPr(k+1)-zPr(k-1))))**2
      
      denom = denom + &
              ((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1)) / (2._knd*(zPr(k+1)-zPr(k-1))))**2
      
      if (abs(denom)>tiny(1._knd)*100) then
       Rig = num/denom
      else
       Rig = 0
      end if
    endfunction Rig






    pure real(knd) function Strainu(S)      !Magnitude of the strain rate tensor.
      real(knd),intent(in) :: S(1:3,1:3)
      integer ::ii,jj

      Strainu = 0
      do jj = 1,3
       do ii = 1,3
        Strainu = Strainu+S(ii,jj)*S(ii,jj)
       end do
      end do
      Strainu = SQRT(2._knd*Strainu)
    endfunction Strainu



    pure subroutine StrainIJ(i,j,k,U,V,W,S)     !Computes components of the strain rate tensor.
      real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
      real(knd),intent(out) :: S(1:3,1:3)
      integer,intent(in) :: i,j,k
      real(knd) D(1:3,1:3)
      integer ::ii,jj

      D = 0

      if (i>-2.and.i<Prnx+3) D(1,1) = (U(i,j,k)-U(i-1,j,k))/dxPr(i)
      if (j>-2.and.j<Prny+3) D(2,2) = (V(i,j,k)-V(i,j-1,k))/dyPr(j)
      if (k>-2.and.j<Prnz+3) D(3,3) = (W(i,j,k)-W(i,j,k-1))/dzPr(k)
      if (j>-1.and.i>-1.and.j<Uny+3.and.i<Unx+4)&
         D(1,2) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(2._knd*(yPr(j+1)-yPr(j-1)))
      if (k>-1.and.i>-1.and.k<Uny+3.and.i<Unx+4)&
         D(1,3) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._knd*(zPr(k+1)-zPr(k-1)))
      if (j>-1.and.i>-1.and.i<Vnx+3.and.j<Vny+4)&
         D(2,1) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(2._knd*(xPr(i+1)-xPr(i-1)))
      if (j>-1.and.k>-1.and.k<Vnz+3.and.j<Vny+4)&
         D(2,3) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._knd*(zPr(k+1)-zPr(k-1)))
      if (k>-1.and.i>-1.and.i<Wnx+3.and.k<Wnz+4)&
         D(3,1) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(2._knd*(xPr(i+1)-xPr(i-1)))
      if (j>-1.and.k>-1.and.j<Wny+3.and.k<Wnz+4)&
         D(3,2) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(2._knd*(yPr(j+1)-yPr(j-1)))

      do jj = 1,3
       do ii = 1,3
        S(ii,jj) = (D(ii,jj)+D(jj,ii))
       end do
      end do
      S = S/2._knd
    endsubroutine StrainIJ




    subroutine SGS_Vreman(U,V,W,filter_ratio)   !Vreman subgrid model (Physics of Fluids, 2004)
      use Tiling, only: tilenx, tileny, tilenz
      real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
      real(knd),intent(in) :: filter_ratio
      integer i,j,k,bi,bj,bk,ii,jj
      real(knd) :: aa,bb
      real(knd),dimension(1:3,1:3) ::a,b
      real(knd) :: dx2,dy2,dz2
      real(knd),parameter ::c = 0.05
      integer,parameter :: narr = 4
      real(knd) :: c2

      c2 = c * filter_ratio**2

      if (gridtype==uniformgrid) then


       dx2 = dxmin**2
       dy2 = dymin**2
       dz2 = dzmin**2

       !$omp parallel do private(aa,bb,a,b,i,j,k,bi,bj,bk,ii,jj) schedule(runtime) !collapse(3)
       do bk = 1,Prnz,tilenz(narr)
        do bj = 1,Prny,tileny(narr)
         do bi = 1,Prnx,tilenx(narr)
          do k = bk,min(bk+tilenz(narr)-1,Prnz)
           do j = bj,min(bj+tileny(narr)-1,Prny)
            do i = bi,min(bi+tilenx(narr)-1,Prnx)
             a(1,1) = (U(i,j,k)-U(i-1,j,k))/dxmin
             a(2,1) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(4._knd*dymin)
             a(3,1) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(4._knd*dzmin)

             a(2,2) = (V(i,j,k)-V(i,j-1,k))/dymin
             a(1,2) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(4._knd*dxmin)
             a(3,2) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(4._knd*dzmin)

             a(3,3) = (W(i,j,k)-W(i,j,k-1))/dzmin
             a(1,3) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(4._knd*dxmin)
             a(2,3) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(4._knd*dymin)

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
               Viscosity(i,j,k) = c2 * sqrt(bb/aa)
             else
               Viscosity(i,j,k) = 0
             end if

             if (Re>0)  Viscosity(i,j,k) = Viscosity(i,j,k)+1._knd/Re

            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end parallel do

      else !general grid

       !$omp parallel do private(aa,bb,a,b,i,j,k,bi,bj,bk,ii,jj) schedule(runtime) !collapse(3)
       do bk = 1,Prnz,tilenz(narr)
        do bj = 1,Prny,tileny(narr)
         do bi = 1,Prnx,tilenx(narr)
          do k = bk,min(bk+tilenz(narr)-1,Prnz)
           do j = bj,min(bj+tileny(narr)-1,Prny)
            do i = bi,min(bi+tilenx(narr)-1,Prnx)
             a(1,1) = (U(i,j,k)-U(i-1,j,k))/dxPr(i)
             a(2,1) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(2._knd*(dyV(j)+dyV(j-1)))
             a(3,1) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._knd*(dzW(k)+dzW(k-1)))
             a(2,2) = (V(i,j,k)-V(i,j-1,k))/dyPr(j)
             a(1,2) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(2._knd*(dxU(i)+dxU(i-1)))
             a(3,2) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._knd*(dzW(k)+dzW(k-1)))
             a(3,3) = (W(i,j,k)-W(i,j,k-1))/dzPr(k)
             a(1,3) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(2._knd*(dxU(i)+dxU(i-1)))
             a(2,3) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(2._knd*(dyV(j)+dyV(j-1)))


             forall(jj = 1:3,ii = 1:3)
              b(ii,jj) = (dxPr(i))**2*a(1,ii)*a(1,jj)+&
                             (dyPr(j))**2*a(2,ii)*a(2,jj)+&
                             (dzPr(k))**2*a(3,ii)*a(3,jj)
             endforall

             bb=   b(1,1)*b(2,2)-b(1,2)**2
             bb = bb+b(1,1)*b(3,3)-b(1,3)**2
             bb = bb+b(2,2)*b(3,3)-b(2,3)**2

             Viscosity(i,j,k) = 0._knd


             aa = sum(a(:,:)**2)

             if (abs(aa)>1e-5.and.bb>0) Viscosity(i,j,k) = c2 * sqrt(bb/aa)

             if (Re>0)  Viscosity(i,j,k) = Viscosity(i,j,k)+1._knd/Re

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
      real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
      real(knd),intent(in) :: filter_ratio
      real(knd),parameter :: Csig = 1.35_knd
      integer,parameter   :: narr = 4
      integer   :: i,j,k,bi,bj,bk
      real(knd) :: width, C, D, g(3,3), s1, s2, s3
      integer,parameter :: sigma_knd = knd

      width = filter_ratio * (dxmin*dymin*dzmin)**(1._knd/3._knd)
      C = (Csig*width)**2

      !$omp parallel do private(g,s1,s2,s3,D,i,j,k,bi,bj,bk) schedule(runtime)
      ! !collapse(3)
      do bk = 1,Prnz,tilenz(narr)
       do bj = 1,Prny,tileny(narr)
        do bi = 1,Prnx,tilenx(narr)
         do k = bk,min(bk+tilenz(narr)-1,Prnz)
          do j = bj,min(bj+tileny(narr)-1,Prny)
           do i = bi,min(bi+tilenx(narr)-1,Prnx)

            call GradientTensorUG(g,i,j,k)

            call Sigmas(s1,s2,s3,g)

            if (s1>0) then
              D = (s3 * (s1 - s2) * (s2 - s3)) / s1**2
            else
              D = 0
            end if

            Viscosity(i,j,k) = C * D

            if (Re>0)  Viscosity(i,j,k) = Viscosity(i,j,k)+1._knd/Re

           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end parallel do

    contains

      pure subroutine GradientTensorUG(g,i,j,k)
        real(knd),intent(out) :: g(3,3)
        integer,intent(in) :: i,j,k

        g(1,1) = (U(i,j,k)-U(i-1,j,k))/dxmin
        g(2,1) = (U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(4._knd*dymin)
        g(3,1) = (U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(4._knd*dzmin)

        g(2,2) = (V(i,j,k)-V(i,j-1,k))/dymin
        g(1,2) = (V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(4._knd*dxmin)
        g(3,2) = (V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(4._knd*dzmin)

        g(3,3) = (W(i,j,k)-W(i,j,k-1))/dzmin
        g(1,3) = (W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(4._knd*dxmin)
        g(2,3) = (W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(4._knd*dymin)
      end subroutine GradientTensorUG


      pure subroutine Sigmas(s1,s2,s3,grads)
        !from Hasan, Basser, Parker, Alexander, http://dx.doi.org/10.1006/jmre.2001.2400
        !via Nicoud, Toda, Cabrit, Bose, Lee, http://dx.doi.org/10.1063/1.3623274
        use ieee_arithmetic

        real(knd),intent(out) :: s1,s2,s3
        real(knd),intent(in)  :: grads(3,3)

        real(sigma_knd) :: trG2, i1, i2, i3, a1, a2, a3, c, G(3,3)

        G = real( matmul(transpose(grads),grads) ,sigma_knd)

        trG2 = dot_product(G(:,1),G(:,1)) +&
               dot_product(G(:,2),G(:,2)) +&
               dot_product(G(:,3),G(:,3))

        i1 = G(1,1) + G(2,2) + G(3,3)

        i2 = (i1**2 - trG2)
        i2 = i2 / 2

        i3 = det3x3(G)

        a1 = max((i1**2)/9 - i2/3,0._sigma_knd)

        a2 = (i1**3)/27 - i1*i2/6 +i3/2

        !This requires no FPE trapping is in progress!
        c =  a2 / sqrt(a1**3)

        !If c is NaN let it be 1.
        c = max(-1._sigma_knd, min(1._sigma_knd,c))
        a3 = acos(c) / 3

        c = 2*sqrt(a1)

        s1 = real( sqrt( i1/3 + c*cos(a3) ) , knd )

        s2 = real( sqrt( max( i1/3 - c*cos(pi/3 + a3) , 0._sigma_knd ) ) , knd)

        s3 = real( sqrt( max( i1/3 - c*cos(pi/3 - a3) , 0._sigma_knd ) ) , knd)

      end subroutine Sigmas


      pure function det3x3(A) result (res)
        real(sigma_knd),intent(in) :: A(3,3)
        real(sigma_knd) :: res

        res =   A(1,1)*A(2,2)*A(3,3)  &
              - A(1,1)*A(2,3)*A(3,2)  &
              - A(1,2)*A(2,1)*A(3,3)  &
              + A(1,2)*A(2,3)*A(3,1)  &
              + A(1,3)*A(2,1)*A(3,2)  &
              - A(1,3)*A(2,2)*A(3,1)
      end function det3x3

    end subroutine SGS_Sigma

end module Subgrid

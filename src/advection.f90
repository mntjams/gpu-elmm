module MomentumAdvection
  use Parameters
  use Boundaries
  use ArrayUtilities
  use Tiling

  implicit none

contains





  subroutine CDUdiv(U2, U, V, W)
    real(knd), contiguous, intent(out) :: U2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd) :: Ax, Ay, Az
    integer :: i, j, k

    Ax = 0.25_knd / dxmin
    Ay = 0.25_knd / dymin
    Az = 0.25_knd / dzmin

    !$omp parallel do private(i, j, k)
    do k=1, Unz
        do j=1, Uny
            do i=1, Unx
                U2(i,j,k)= - ((Ax*(U(i+1,j,k) + U(i,j,k)) * (U(i+1,j,k) + U(i,j,k)) &
                             - Ax*(U(i,j,k) + U(i-1,j,k)) * (U(i,j,k) + U(i-1,j,k))) &
                            + (Ay*(U(i,j+1,k) + U(i,j,k)) * (V(i+1,j,k) + V(i,j,k)) &
                             - Ay*(U(i,j,k) + U(i,j-1,k)) * (V(i+1,j-1,k) + V(i,j-1,k))) &
                            + (Az*(U(i,j,k+1) + U(i,j,k)) * (W(i+1,j,k) + W(i,j,k)) &
                             - Az*(U(i,j,k) + U(i,j,k-1)) * (W(i+1,j,k-1) + W(i,j,k-1))))
            end do
        end do
    end do
    !$omp end parallel do
  end subroutine CDUdiv






  subroutine CDVdiv(V2, U, V, W)
    real(knd), contiguous, intent(out) :: V2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd) :: Ax, Ay, Az
    integer :: i, j, k

    Ax = 0.25_knd / dxmin
    Ay = 0.25_knd / dymin
    Az = 0.25_knd / dzmin

    !$omp parallel do private(i, j, k)
    do k=1, Vnz
        do j=1, Vny
            do i=1, Vnx
                V2(i,j,k)= - ((Ay*(V(i,j+1,k) + V(i,j,k)) * (V(i,j+1,k) + V(i,j,k)) &
                              -Ay*(V(i,j,k) + V(i,j-1,k)) * (V(i,j,k) + V(i,j-1,k))) &
                             +(Ax*(V(i+1,j,k) + V(i,j,k)) * (U(i,j+1,k) + U(i,j,k)) &
                              -Ax*(V(i,j,k) + V(i-1,j,k)) * (U(i-1,j+1,k) + U(i-1,j,k))) &
                             +(Az*(V(i,j,k+1) + V(i,j,k)) * (W(i,j+1,k) + W(i,j,k)) &
                              -Az*(V(i,j,k) + V(i,j,k-1)) * (W(i,j+1,k-1) + W(i,j,k-1))))
            end do
        end do
    end do
    !$omp end parallel do
  end subroutine CDVdiv





  subroutine CDWdiv(W2, U, V, W)
    real(knd), contiguous, intent(out) :: W2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd) :: Ax, Ay, Az
    integer  :: i, j, k

    Ax = 0.25_knd / dxmin
    Ay = 0.25_knd / dymin
    Az = 0.25_knd / dzmin

    !$omp parallel do private(i, j, k)
    do k=1, Wnz
        do j=1, Wny
            do i=1, Wnx
                W2(i,j,k)= - ((Az*(W(i,j,k+1) + W(i,j,k)) * (W(i,j,k+1) + W(i,j,k)) &
                             - Az*(W(i,j,k) + W(i,j,k-1)) * (W(i,j,k) + W(i,j,k-1))) &
                            + (Ay*(W(i,j+1,k) + W(i,j,k)) * (V(i,j,k+1) + V(i,j,k)) &
                             - Ay*(W(i,j,k) + W(i,j-1,k)) * (V(i,j-1,k) + V(i,j-1,k+1))) &
                            + (Ax*(W(i+1,j,k) + W(i,j,k)) * (U(i,j,k+1) + U(i,j,k)) &
                             - Ax*(W(i,j,k) + W(i-1,j,k)) * (U(i-1,j,k+1) + U(i-1,j,k))))
            end do
        end do
    end do
    !$omp end parallel do
  end subroutine CDWdiv











  subroutine CDUadv(U2, U, V, W)
  real(knd), contiguous, intent(out) :: U2(-2:,-2:,-2:)
  real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
  integer i, j, k
  real(knd) Ax, Ay, Az, Vadv, Wadv


   if (gridtype==UNIFORMGRID) then
      Ax=0.5/dxmin
      Ay=0.125/dymin
      Az=0.125/dzmin

      !$omp parallel do private(i, j, k, Vadv, Wadv)
      do k=1, Unz
          do j=1, Uny
              do i=1, Unx
                  Vadv = ( V(i,j,k) + V(i+1,j,k) + V(i,j-1,k) + V(i+1,j-1,k) )
                  Wadv = ( W(i,j,k) + W(i+1,j,k) + W(i,j,k-1) + W(i+1,j,k-1) )
                  U2(i,j,k)= U2(i,j,k) &
                             - (Ax*(U(i+1,j,k)-U(i-1,j,k)) * U(i,j,k) &
                             +  Ay*(U(i,j+1,k)-U(i,j-1,k)) * Vadv&
                             +  Az*(U(i,j,k+1)-U(i,j,k-1)) * Wadv )
              end do
          end do
      end do

   end if
  end subroutine CDUadv






  subroutine CDVadv(V2, U, V, W)
  real(knd), contiguous, intent(out) :: V2(-2:,-2:,-2:)
  real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
  integer i, j, k
  real(knd) Ax, Ay, Az, Uadv, Wadv


   if (gridtype==UNIFORMGRID) then
      Ax=0.125/dxmin
      Ay=0.5/dymin
      Az=0.125/dzmin

      !$omp parallel do private(i, j, k, Uadv, Wadv)
      do k=1, Vnz
          do j=1, Vny
              do i=1, Vnx
                  Uadv = ( U(i,j,k) + U(i,j+1,k) + U(i-1,j,k) + U(i-1,j+1,k) )
                  Wadv = ( W(i,j,k) + W(i,j+1,k) + W(i,j,k-1) + W(i,j+1,k-1) )
                  V2(i,j,k)= V2(i,j,k) &
                             - (Ax*(V(i+1,j,k)-V(i-1,j,k)) * Uadv&
                             +  Ay*(V(i,j+1,k)-V(i,j-1,k)) * V(i,j,k) &
                             +  Az*(V(i,j,k+1)-V(i,j,k-1)) * Wadv )
              end do
          end do
      end do
      !$omp end parallel do

   end if
  end subroutine CDVadv





  subroutine CDWadv(W2, U, V, W)
  real(knd), contiguous, intent(out) :: W2(-2:,-2:,-2:)
  real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
  integer i, j, k
  real(knd) Ax, Ay, Az, Uadv, Vadv


   if (gridtype==UNIFORMGRID) then
      Ax=0.125/dxmin
      Ay=0.125/dymin
      Az=0.5/dzmin

      !$omp parallel do private(i, j, k, Uadv, Vadv)
      do k=1, Wnz
          do j=1, Wny
              do i=1, Wnx
                  Uadv = ( U(i,j,k) + U(i,j,k+1) + U(i-1,j,k) + U(i-1,j,k+1) )
                  Vadv = ( V(i,j,k) + V(i,j,k+1) + V(i,j-1,k) + V(i,j-1,k+1) )
                  W2(i,j,k)= W2(i,j,k) &
                             - (Ax*(W(i+1,j,k)-W(i-1,j,k)) * Uadv&
                             +  Ay*(W(i,j+1,k)-W(i,j-1,k)) * Vadv&
                             +  Az*(W(i,j,k+1)-W(i,j,k-1)) * W(i,j,k) )
              end do
          end do
      end do
      !$omp end parallel do

   end if
  end subroutine CDWadv















  subroutine CDU(U2, U, V, W)
  real(knd), contiguous, intent(out) :: U2(-2:,-2:,-2:)
  real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)

   U2 = 0
   call CDUdiv(U2, U, V, W)
   call CDUadv(U2, U, V, W)
   U2 = U2 / 2
  end subroutine CDU






  subroutine CDV(V2, U, V, W)
  real(knd), contiguous, intent(out) :: V2(-2:,-2:,-2:)
  real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)

   V2 = 0
   call CDVdiv(V2, U, V, W)
   call CDVadv(V2, U, V, W)
   V2 = V2 / 2

  end subroutine CDV





  subroutine CDW(W2, U, V, W)
  real(knd), contiguous, intent(out) :: W2(-2:,-2:,-2:)
  real(knd), contiguous, intent(in)  :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)

   W2 = 0
   call CDWdiv(W2, U, V, W)
   call CDWadv(W2, U, V, W)
   W2 = W2 / 2

  end subroutine CDW
















  subroutine CDS4U(U2, U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(out) :: U2
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)  :: U, V, W
    integer, parameter :: mask4(4) = [ -1, 9, 9, -1 ]
    integer, parameter :: divcoef = 256
    integer, parameter :: mask2(2) = [ 8, 8 ]
    integer   :: i, j, k, l
    real(knd) :: flux

    !$omp parallel private(i, j, k, l, flux)
    !$omp workshare
    U2 = 0
    !$omp end workshare

    !$omp do
    do k=1, Unz
      do j=1, Uny
        do i=1, Unx+1
          if (Utype(i-1,j,k)<=0 .and. Utype(i,j,k)<=0) then

            flux = ( ( sum(mask4 * U(i-2:i+1,j,k)) )**2 ) / dxmin
            U2(i-1,j,k) = U2(i-1,j,k) - flux
            U2(i,j,k)   = U2(i,j,k)   + flux

          else if (Utype(i-1,j,k)<=0 .or. Utype(i,j,k)<=0) then

            flux = ( ( sum(mask2 * U(i-1:i,j,k)) )**2 ) / dxmin
            U2(i-1,j,k) = U2(i-1,j,k) - flux
            U2(i,j,k)   = U2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do k=1, Unz
      do j=1, Uny+1
        do i=1, Unx
          if (Utype(i,j-1,k)<=0 .and. Utype(i,j,k)<=0) then

            flux = ( sum(mask4 * U(i,j-2:j+1,k)) * sum(mask4 * V(i-1:i+2,j-1,k)) ) / dymin
            U2(i,j-1,k) = U2(i,j-1,k) - flux
            U2(i,j,k)   = U2(i,j,k)   + flux

          else if (Utype(i,j-1,k)<=0 .or. Utype(i,j,k)<=0) then

            flux = ( sum(mask2 * U(i,j-1:j,k)) * sum(mask2 * V(i:i+1,j-1,k)) ) / dymin
            U2(i,j-1,k) = U2(i,j-1,k) - flux
            U2(i,j,k)   = U2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do

    do l=1, 2
      !$omp do
      do k=l, Unz+1, 2
        do j=1, Uny
          do i=1, Unx
            if (Utype(i,j,k-1)<=0 .and. Utype(i,j,k)<=0) then

                flux = ( sum(mask4 * U(i,j,k-2:k+1))) * sum(mask4 * W(i-1:i+2,j,k-1)) / dzmin
                U2(i,j,k-1) = U2(i,j,k-1) - flux
                U2(i,j,k)   = U2(i,j,k)   + flux

            else if (Utype(i,j,k-1)<=0 .or. Utype(i,j,k)<=0) then

                flux = ( sum(mask2 * U(i,j,k-1:k))) * sum(mask2 * W(i:i+1,j,k-1)) / dzmin
                U2(i,j,k-1) = U2(i,j,k-1) - flux
                U2(i,j,k)   = U2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do
    end do

    !$omp workshare
    U2 = U2 / divcoef
    !$omp end workshare
    !$omp end parallel

  end subroutine CDS4U


  subroutine CDS4V(V2, U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(out) :: V2
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)  :: U, V, W
    integer, parameter :: mask4(4) = [ -1, 9, 9, -1 ]
    integer, parameter :: divcoef = 256
    integer, parameter :: mask2(2) = [ 8, 8 ]
    integer   :: i, j, k, l
    real(knd) :: flux

    !$omp parallel private(i, j, k, l, flux)
    !$omp workshare
    V2 = 0
    !$omp end workshare

    !$omp do
    do k=1, Vnz
      do j=1, Vny
        do i=1, Vnx+1
          if (Vtype(i-1,j,k)<=0 .and. Vtype(i,j,k)<=0) then

            flux = ( sum(mask4 * V(i-2:i+1,j,k)) * sum(mask4 * U(i-1,j-1:j+2,k)) ) / dxmin
            V2(i-1,j,k) = V2(i-1,j,k) - flux
            V2(i,j,k)   = V2(i,j,k)   + flux

          else if (Vtype(i-1,j,k)<=0 .or. Vtype(i,j,k)<=0) then

            flux = ( sum(mask2 * V(i-1:i,j,k)) * sum(mask2 * U(i-1,j:j+1,k)) ) / dxmin
            V2(i-1,j,k) = V2(i-1,j,k) - flux
            V2(i,j,k)   = V2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do


    !$omp do
    do k=1, Vnz
      do j=1, Vny+1
        do i=1, Vnx
          if (Vtype(i,j-1,k)<=0 .and. Vtype(i,j,k)<=0) then

            flux = ( ( sum(mask4 * V(i,j-2:j+1,k)) )**2 ) / dymin
            V2(i,j-1,k) = V2(i,j-1,k) - flux
            V2(i,j,k)   = V2(i,j,k)   + flux

          else if (Vtype(i,j-1,k)<=0 .or. Vtype(i,j,k)<=0) then

            flux = ( ( sum(mask2 * V(i,j-1:j,k)) )**2 ) / dymin
            V2(i,j-1,k) = V2(i,j-1,k) - flux
            V2(i,j,k)   = V2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do

    do l=1, 2
      !$omp do
      do k=l, Vnz+1, 2
        do j=1, Vny
          do i=1, Vnx
            if (Vtype(i,j,k-1)<=0 .and. Vtype(i,j,k)<=0) then

              flux = ( sum(mask4 * V(i,j,k-2:k+1)) * sum(mask4 * W(i,j-1:j+2,k-1)) ) / dzmin
              V2(i,j,k-1) = V2(i,j,k-1) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            else if (Vtype(i,j,k-1)<=0 .or. Vtype(i,j,k)<=0) then

              flux = ( sum(mask2 * V(i,j,k-1:k)) * sum(mask2 * W(i,j:j+1,k-1)) ) / dzmin
              V2(i,j,k-1) = V2(i,j,k-1) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do
    end do

    !$omp workshare
    V2 = V2 / divcoef
    !$omp end workshare
    !$omp end parallel

  end subroutine CDS4V


  subroutine CDS4W(W2, U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(out) :: W2
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)  :: U, V, W
    integer, parameter :: mask4(4) = [ -1, 9, 9, -1 ]
    integer, parameter :: divcoef = 256
    integer, parameter :: mask2(2) = [  8, 8 ]
    integer   :: i, j, k, l
    real(knd) :: flux

    !$omp parallel private(i, j, k, l, flux)
    !$omp workshare
    W2 = 0
    !$omp end workshare

    !$omp do
    do k=1, Wnz
      do j=1, Wny
        do i=1, Wnx+1
          if (Wtype(i-1,j,k)<=0 .and. Wtype(i,j,k)<=0) then

            flux = ( sum(mask4 * W(i-2:i+1,j,k)) * sum(mask4 * U(i-1,j,k-1:k+2)) ) / dxmin
            W2(i-1,j,k) = W2(i-1,j,k) - flux
            W2(i,j,k)   = W2(i,j,k)   + flux

          else if (Wtype(i-1,j,k)<=0 .or. Wtype(i,j,k)<=0) then

            flux = ( sum(mask2 * W(i-1:i,j,k)) * sum(mask2 * U(i-1,j,k:k+1)) ) / dxmin
            W2(i-1,j,k) = W2(i-1,j,k) - flux
            W2(i,j,k)   = W2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do k=1, Wnz
      do j=1, Wny+1
        do i=1, Wnx
          if (Wtype(i,j-1,k)<=0 .and. Wtype(i,j,k)<=0) then

            flux = ( sum(mask4 * W(i,j-2:j+1,k)) * sum(mask4 * V(i,j-1,k-1:k+2)) ) / dymin
            W2(i,j-1,k) = W2(i,j-1,k) - flux
            W2(i,j,k)   = W2(i,j,k)   + flux

          else if (Wtype(i,j-1,k)<=0 .or. Wtype(i,j,k)<=0) then

            flux = ( sum(mask2 * W(i,j-1:j,k)) * sum(mask2 * V(i,j-1,k:k+1)) ) / dymin
            W2(i,j-1,k) = W2(i,j-1,k) - flux
            W2(i,j,k)   = W2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do

    do l=1, 2
      !$omp do
      do k=l, Wnz+1, 2
        do j=1, Wny
          do i=1, Wnx
            if (Wtype(i,j,k-1)<=0 .and. Wtype(i,j,k)<=0) then

              flux = ( ( sum(mask4 * W(i,j,k-2:k+1)) )**2 ) / dzmin
              W2(i,j,k-1) = W2(i,j,k-1) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            else if (Wtype(i,j,k-1)<=0 .or. Wtype(i,j,k)<=0) then

              flux = ( ( sum(mask2 * W(i,j,k-1:k)) )**2 ) / dzmin
              W2(i,j,k-1) = W2(i,j,k-1) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do
    end do

    !$omp workshare
    W2 = W2 / divcoef
    !$omp end workshare
    !$omp end parallel

  end subroutine CDS4W

















  subroutine CDS4_2U(U2, U,V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(out) :: U2
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)  :: U, V, W
    integer, parameter :: mask4(4) = [ -1, 9, 9, -1 ]
    integer, parameter :: divcoef = 256
    integer, parameter :: mask2(2) = [ 8, 8 ]
    integer   :: i, j, k, l
    real(knd) :: flux

    !$omp parallel private(i, j, k, l, flux)
    !$omp workshare
    U2 = 0
    !$omp end workshare

    !$omp do
    do k=1, Unz
      do j=1, Uny
        do i=1, Unx+1
          if (Utype(i-1,j,k)<=0 .and. Utype(i,j,k)<=0) then

            flux = ( ( sum(mask4 * U(i-2:i+1,j,k)) )**2 ) / dxmin
            U2(i-1,j,k) = U2(i-1,j,k) - flux
            U2(i,j,k)   = U2(i,j,k)   + flux

          else if (Utype(i-1,j,k)<=0 .or. Utype(i,j,k)<=0) then

            flux = ( ( sum(mask2 * U(i-1:i,j,k)) )**2 ) / dxmin
            U2(i-1,j,k) = U2(i-1,j,k) - flux
            U2(i,j,k)   = U2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do k=1, Unz
      do j=1, Uny+1
        do i=1, Unx
          if (Utype(i,j-1,k)<=0 .and. Utype(i,j,k)<=0) then

            flux = ( sum(mask4 * U(i,j-2:j+1,k)) * sum(mask2 * V(i:i+1,j-1,k)) ) / dymin
            U2(i,j-1,k) = U2(i,j-1,k) - flux
            U2(i,j,k)   = U2(i,j,k)   + flux

          else if (Utype(i,j-1,k)<=0 .or. Utype(i,j,k)<=0) then

            flux = ( sum(mask2 * U(i,j-1:j,k)) * sum(mask2 * V(i:i+1,j-1,k)) ) / dymin
            U2(i,j-1,k) = U2(i,j-1,k) - flux
            U2(i,j,k)   = U2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do

    do l=1, 2
      !$omp do
      do k=l, Unz+1, 2
        do j=1, Uny
          do i=1, Unx
            if (Utype(i,j,k-1)<=0 .and. Utype(i,j,k)<=0) then

                flux = ( sum(mask4 * U(i,j,k-2:k+1))) * sum(mask2 * W(i:i+1,j,k-1)) / dzmin
                U2(i,j,k-1) = U2(i,j,k-1) - flux
                U2(i,j,k)   = U2(i,j,k)   + flux

            else if (Utype(i,j,k-1)<=0 .or. Utype(i,j,k)<=0) then

                flux = ( sum(mask2 * U(i,j,k-1:k))) * sum(mask2 * W(i:i+1,j,k-1)) / dzmin
                U2(i,j,k-1) = U2(i,j,k-1) - flux
                U2(i,j,k)   = U2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do
    end do

    !$omp workshare
    U2 = U2 / divcoef
    !$omp end workshare
    !$omp end parallel

  end subroutine CDS4_2U


  subroutine CDS4_2V(V2, U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(out) :: V2
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)  :: U, V, W
    integer, parameter :: mask4(4) = [ -1, 9, 9, -1 ]
    integer, parameter :: divcoef = 256
    integer, parameter :: mask2(2) = [ 8, 8 ]
    integer   :: i, j, k, l
    real(knd) :: flux

    !$omp parallel private(i, j, k, l, flux)
    !$omp workshare
    V2 = 0
    !$omp end workshare

    !$omp do
    do k=1, Vnz
      do j=1, Vny
        do i=1, Vnx+1
          if (Vtype(i-1,j,k)<=0 .and. Vtype(i,j,k)<=0) then

            flux = ( sum(mask4 * V(i-2:i+1,j,k)) * sum(mask2 * U(i-1,j:j+1,k)) ) / dxmin
            V2(i-1,j,k) = V2(i-1,j,k) - flux
            V2(i,j,k)   = V2(i,j,k)   + flux

          else if (Vtype(i-1,j,k)<=0 .or. Vtype(i,j,k)<=0) then

            flux = ( sum(mask2 * V(i-1:i,j,k)) * sum(mask2 * U(i-1,j:j+1,k)) ) / dxmin
            V2(i-1,j,k) = V2(i-1,j,k) - flux
            V2(i,j,k)   = V2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do


    !$omp do
    do k=1, Vnz
      do j=1, Vny+1
        do i=1, Vnx
          if (Vtype(i,j-1,k)<=0 .and. Vtype(i,j,k)<=0) then

            flux = ( ( sum(mask4 * V(i,j-2:j+1,k)) )**2 ) / dymin
            V2(i,j-1,k) = V2(i,j-1,k) - flux
            V2(i,j,k)   = V2(i,j,k)   + flux

          else if (Vtype(i,j-1,k)<=0 .or. Vtype(i,j,k)<=0) then

            flux = ( ( sum(mask2 * V(i,j-1:j,k)) )**2 ) / dymin
            V2(i,j-1,k) = V2(i,j-1,k) - flux
            V2(i,j,k)   = V2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do

    do l=1, 2
      !$omp do
      do k=l, Vnz+1, 2
        do j=1, Vny
          do i=1, Vnx
            if (Vtype(i,j,k-1)<=0 .and. Vtype(i,j,k)<=0) then

              flux = ( sum(mask4 * V(i,j,k-2:k+1)) * sum(mask2 * W(i,j:j+1,k-1)) ) / dzmin
              V2(i,j,k-1) = V2(i,j,k-1) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            else if (Vtype(i,j,k-1)<=0 .or. Vtype(i,j,k)<=0) then

              flux = ( sum(mask2 * V(i,j,k-1:k)) * sum(mask2 * W(i,j:j+1,k-1)) ) / dzmin
              V2(i,j,k-1) = V2(i,j,k-1) - flux
              V2(i,j,k)   = V2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do
    end do

    !$omp workshare
    V2 = V2 / divcoef
    !$omp end workshare
    !$omp end parallel

  end subroutine CDS4_2V


  subroutine CDS4_2W(W2, U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(out) :: W2
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)  :: U, V, W
    integer, parameter :: mask4(4) = [ -1, 9, 9, -1 ]
    integer, parameter :: divcoef = 256
    integer, parameter :: mask2(2) = [  8, 8 ]
    integer   :: i, j, k, l
    real(knd) :: flux

    !$omp parallel private(i, j, k, l, flux)
    !$omp workshare
    W2 = 0
    !$omp end workshare

    !$omp do
    do k=1, Wnz
      do j=1, Wny
        do i=1, Wnx+1
          if (Wtype(i-1,j,k)<=0 .and. Wtype(i,j,k)<=0) then

            flux = ( sum(mask4 * W(i-2:i+1,j,k)) * sum(mask2 * U(i-1,j,k:k+1)) ) / dxmin
            W2(i-1,j,k) = W2(i-1,j,k) - flux
            W2(i,j,k)   = W2(i,j,k)   + flux

          else if (Wtype(i-1,j,k)<=0 .or. Wtype(i,j,k)<=0) then

            flux = ( sum(mask2 * W(i-1:i,j,k)) * sum(mask2 * U(i-1,j,k:k+1)) ) / dxmin
            W2(i-1,j,k) = W2(i-1,j,k) - flux
            W2(i,j,k)   = W2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do k=1, Wnz
      do j=1, Wny+1
        do i=1, Wnx
          if (Wtype(i,j-1,k)<=0 .and. Wtype(i,j,k)<=0) then

            flux = ( sum(mask4 * W(i,j-2:j+1,k)) * sum(mask2 * V(i,j-1,k:k+1)) ) / dymin
            W2(i,j-1,k) = W2(i,j-1,k) - flux
            W2(i,j,k)   = W2(i,j,k)   + flux

          else if (Wtype(i,j-1,k)<=0 .or. Wtype(i,j,k)<=0) then

            flux = ( sum(mask2 * W(i,j-1:j,k)) * sum(mask2 * V(i,j-1,k:k+1)) ) / dymin
            W2(i,j-1,k) = W2(i,j-1,k) - flux
            W2(i,j,k)   = W2(i,j,k)   + flux

          end if
        end do
      end do
    end do
    !$omp end do

    do l=1, 2
      !$omp do
      do k=l, Wnz+1, 2
        do j=1, Wny
          do i=1, Wnx
            if (Wtype(i,j,k-1)<=0 .and. Wtype(i,j,k)<=0) then

              flux = ( ( sum(mask4 * W(i,j,k-2:k+1)) )**2 ) / dzmin
              W2(i,j,k-1) = W2(i,j,k-1) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            else if (Wtype(i,j,k-1)<=0 .or. Wtype(i,j,k)<=0) then

              flux = ( ( sum(mask2 * W(i,j,k-1:k)) )**2 ) / dzmin
              W2(i,j,k-1) = W2(i,j,k-1) - flux
              W2(i,j,k)   = W2(i,j,k)   + flux

            end if
          end do
        end do
      end do
      !$omp end do
    end do

    !$omp workshare
    W2 = W2 / divcoef
    !$omp end workshare
    !$omp end parallel

  end subroutine CDS4_2W













  subroutine CD4divU(U2, U, V, W)
    !Morinishi et al., JCP 143, http://dx.doi.org/10.1006/jcph.1998.5962
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(out) :: U2
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)  :: U, V, W
    integer, parameter :: i_mask4(4) = [ -1, 9, 9, -1 ]
    integer, parameter :: divcoef = 256
    integer, parameter :: coef2ord = divcoef / 4
    integer, parameter :: narr = 4
    !UV1 and UV3 are sized to fit to the L1 cache, so no stack overflow should be possible
    real(knd) :: UV1(-1:tilenx(narr)+2, -1:tileny(narr)+2, -1:tilenz(narr)+2)
    real(knd) :: UV3(-1:tilenx(narr)+2, -1:tileny(narr)+2, -1:tilenz(narr)+2)
    integer   :: bi, bj, bk, i, j, k, li, lj, lk
    real(knd) :: Uint, Vint, Wint, dU

    call set(U2, 0._knd)

    !$omp parallel private(bi, bj, bk, i, j, k, li, lj, lk, Uint, Vint, Wint, dU, UV1, UV3)
    !$omp do
    do bk = 1, Unz, tilenz(narr)
     do bj = 1, Uny, tileny(narr)
      do bi = 1, Unx, tilenx(narr)

        do k = bk, min(bk+tilenz(narr)-1, Unz)
         do j = bj, min(bj+tileny(narr)-1, Uny)
          do i = bi-1, min(bi+tilenx(narr)-1, Unx)+2
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 Uj_1 - 1/8 Uj_3
            Uint = sum( i_mask4 * U(i-2:i+1,j,k) )
            UV1(li,lj,lk) = Uint * (U(i-1,j,k) + U(i  ,j,k))
            UV3(li,lj,lk) = Uint * (U(i-2,j,k) + U(i+1,j,k))
          end do
         end do
        end do

        do k = bk, min(bk+tilenz(narr)-1, Unz)
         do j = bj, min(bj+tileny(narr)-1, Uny)
          do i = bi, min(bi+tilenx(narr)-1, Unx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 d1 - 1/8 d3
            if (Utype(i,j,k)==0) then
              dU =     9 * (UV1(li+1,lj,lk) - UV1(li  ,lj,lk))
              dU = dU -    (UV3(li+2,lj,lk) - UV3(li-1,lj,lk))/3
              U2(i,j,k) = U2(i,j,k) + dU / dxmin
            else if (Utype(i,j,k)<0) then  !near a boundary - 2nd order
              dU = coef2ord * ((U(i+1,j,k) + U(i  ,j,k)) * (U(i+1,j,k) + U(i  ,j,k)) &
                              -(U(i  ,j,k) + U(i-1,j,k)) * (U(i  ,j,k) + U(i-1,j,k)))
              U2(i,j,k) = U2(i,j,k) + dU / dxmin
            end if
          end do
         end do
        end do

      end do
     end do
    end do
    !$omp end do

    !$omp do
    do bk = 1, Unz, tilenz(narr)
     do bj = 1, Uny, tileny(narr)
      do bi = 1, Unx, tilenx(narr)

        do k = bk, min(bk+tilenz(narr)-1, Unz)
         do j = bj-2, min(bj+tileny(narr)-1, Uny)+1
          do i = bi, min(bi+tilenx(narr)-1, Unx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 Uj_1 - 1/8 Uj_3
            Vint = sum( i_mask4 * V(i-1:i+2,j,k) )
            UV1(li,lj,lk) = Vint * (U(i,j  ,k) + U(i,j+1,k))
            UV3(li,lj,lk) = Vint * (U(i,j-1,k) + U(i,j+2,k))
          end do
         end do
        end do

        do k = bk, min(bk+tilenz(narr)-1, Unz)
         do j = bj, min(bj+tileny(narr)-1, Uny)
          do i = bi, min(bi+tilenx(narr)-1, Unx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 d1 - 1/8 d3
            if (Utype(i,j,k)==0) then
              dU =     9 * (UV1(li,lj  ,lk) - UV1(li,lj-1,lk))
              dU = dU -    (UV3(li,lj+1,lk) - UV3(li,lj-2,lk))/3
              U2(i,j,k) = U2(i,j,k) + dU / dymin
            else if (Utype(i,j,k)<0) then  !near a boundary - 2nd order
              dU = coef2ord * ((U(i,j+1,k) + U(i,j  ,k)) * (V(i+1,j  ,k) + V(i,j  ,k)) &
                              -(U(i,j  ,k) + U(i,j-1,k)) * (V(i+1,j-1,k) + V(i,j-1,k)))
              U2(i,j,k) = U2(i,j,k) + dU / dymin
            end if
          end do
         end do
        end do

      end do
     end do
    end do
    !$omp end do

    !$omp do
    do bk = 1, Unz, tilenz(narr)
     do bj = 1, Uny, tileny(narr)
      do bi = 1, Unx, tilenx(narr)

        do k = bk-2, min(bk+tilenz(narr), Unz)+1
         do j = bj, min(bj+tileny(narr)-1, Uny)
          do i = bi, min(bi+tilenx(narr)-1, Unx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 Uj_1 - 1/8 Uj_3
            Wint = sum( i_mask4 * W(i-1:i+2,j,k) )
            UV1(li,lj,lk) = Wint * (U(i,j,k  ) + U(i,j,k+1))
            UV3(li,lj,lk) = Wint * (U(i,j,k-1) + U(i,j,k+2))
          end do
         end do
        end do

        do k = bk, min(bk+tilenz(narr)-1, Unz)
         do j = bj, min(bj+tileny(narr)-1, Uny)
          do i = bi, min(bi+tilenx(narr)-1, Unx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 d1 - 1/8 d3
            if (Utype(i,j,k)==0) then
              dU =     9 * (UV1(li,lj,lk  ) - UV1(li,lj,lk-1))
              dU = dU -    (UV3(li,lj,lk+1) - UV3(li,lj,lk-2))/3
              U2(i,j,k) = U2(i,j,k) + dU / dzmin
            else if (Utype(i,j,k)<0) then  !near a boundary - 2nd order
              dU = coef2ord * ((U(i,j,k+1) + U(i,j,k  )) * (W(i+1,j,k  ) + W(i,j,k  )) &
                              -(U(i,j,k  ) + U(i,j,k-1)) * (W(i+1,j,k-1) + W(i,j,k-1)))
              U2(i,j,k) = U2(i,j,k) + dU / dzmin
            end if
          end do
         end do
        end do

      end do
     end do
    end do
    !$omp end do
    !$omp end parallel

    call multiply(U2, -1._knd/divcoef)

  end subroutine CD4divU











  subroutine CD4divV(V2, U, V, W)
    !Morinishi et al., JCP 143, http://dx.doi.org/10.1006/jcph.1998.5962
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(out) :: V2
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)  :: U, V, W
    integer, parameter :: i_mask4(4) = [ -1, 9, 9, -1 ]
    integer, parameter :: divcoef = 256
    integer, parameter :: coef2ord = divcoef / 4
    integer, parameter :: narr = 4
    !UV1 and UV3 are sized to fit to the L1 cache, so no stack overflow should be possible
    real(knd) :: UV1(-1:tilenx(narr)+2, -1:tileny(narr)+2, -1:tilenz(narr)+2)
    real(knd) :: UV3(-1:tilenx(narr)+2, -1:tileny(narr)+2, -1:tilenz(narr)+2)
    integer   :: bi, bj, bk, i, j, k, li, lj, lk
    real(knd) :: Uint, Vint, Wint, dV

    call set(V2, 0._knd)

    !$omp parallel private(bi, bj, bk, i, j, k, li, lj, lk, Uint, Vint, Wint, dV, UV1, UV3)
    !$omp do
    do bk = 1, Vnz, tilenz(narr)
     do bj = 1, Vny, tileny(narr)
      do bi = 1, Vnx, tilenx(narr)

        do k = bk, min(bk+tilenz(narr)-1, Vnz)
         do j = bj, min(bj+tileny(narr)-1, Vny)
          do i = bi-2, min(bi+tilenx(narr)-1, Vnx)+1
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 Uj_1 - 1/8 Uj_3
            Uint = sum( i_mask4 * U(i,j-1:j+2,k) )
            UV1(li,lj,lk) = Uint * (V(i  ,j,k) + V(i+1,j,k))
            UV3(li,lj,lk) = Uint * (V(i-1,j,k) + V(i+2,j,k))
          end do
         end do
        end do

        do k = bk, min(bk+tilenz(narr)-1, Vnz)
         do j = bj, min(bj+tileny(narr)-1, Vny)
          do i = bi, min(bi+tilenx(narr)-1, Vnx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 d1 - 1/8 d3
            if (Vtype(i,j,k)==0) then
              dV =     9 * (UV1(li  ,lj,lk) - UV1(li-1,lj,lk))
              dV = dV -    (UV3(li+1,lj,lk) - UV3(li-2,lj,lk))/3
              V2(i,j,k) = V2(i,j,k) + dV / dxmin
            else if (Vtype(i,j,k)<0) then  !near a boundary - 2nd order
              dV = coef2ord * ((V(i+1,j,k) + V(i  ,j,k)) * (U(i  ,j+1,k) + U(i  ,j,k)) &
                              -(V(i  ,j,k) + V(i-1,j,k))  *(U(i-1,j+1,k) + U(i-1,j,k)))
              V2(i,j,k) = V2(i,j,k) + dV / dxmin
            end if
          end do
         end do
        end do

      end do
     end do
    end do
    !$omp end do

    !$omp do
    do bk = 1, Vnz, tilenz(narr)
     do bj = 1, Vny, tileny(narr)
      do bi = 1, Vnx, tilenx(narr)

        do k = bk, min(bk+tilenz(narr)-1, Vnz)
         do j = bj-1, min(bj+tileny(narr)-1, Vny)+2
          do i = bi, min(bi+tilenx(narr)-1, Vnx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 Uj_1 - 1/8 Uj_3
            Vint = sum( i_mask4 * V(i,j-2:j+1,k) )
            UV1(li,lj,lk) = Vint * (V(i,j-1,k) + V(i,j  ,k))
            UV3(li,lj,lk) = Vint * (V(i,j-2,k) + V(i,j+1,k))
          end do
         end do
        end do

        do k = bk, min(bk+tilenz(narr)-1, Vnz)
         do j = bj, min(bj+tileny(narr)-1, Vny)
          do i = bi, min(bi+tilenx(narr)-1, Vnx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 d1 - 1/8 d3
            if (Vtype(i,j,k)==0) then
              dV =     9 * (UV1(li,lj+1,lk) - UV1(li,lj  ,lk))
              dV = dV -    (UV3(li,lj+2,lk) - UV3(li,lj-1,lk))/3
              V2(i,j,k) = V2(i,j,k) + dV / dymin
            else if (Vtype(i,j,k)<0) then  !near a boundary - 2nd order
              dV = coef2ord * ((V(i,j+1,k) + V(i,j  ,k)) * (V(i,j+1,k) + V(i,j  ,k)) &
                              -(V(i,j  ,k) + V(i,j-1,k)) * (V(i,j  ,k) + V(i,j-1,k)))
              V2(i,j,k) = V2(i,j,k) + dV / dymin
            end if
          end do
         end do
        end do

      end do
     end do
    end do
    !$omp end do

    !$omp do
    do bk = 1, Vnz, tilenz(narr)
     do bj = 1, Vny, tileny(narr)
      do bi = 1, Vnx, tilenx(narr)

        do k = bk-2, min(bk+tilenz(narr), Vnz)+1
         do j = bj, min(bj+tileny(narr)-1, Vny)
          do i = bi, min(bi+tilenx(narr)-1, Vnx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 Uj_1 - 1/8 Uj_3
            Wint = sum( i_mask4 * W(i,j-1:j+2,k) )
            UV1(li,lj,lk) = Wint * (V(i,j,k  ) + V(i,j,k+1))
            UV3(li,lj,lk) = Wint * (V(i,j,k-1) + V(i,j,k+2))
          end do
         end do
        end do

        do k = bk, min(bk+tilenz(narr)-1, Vnz)
         do j = bj, min(bj+tileny(narr)-1, Vny)
          do i = bi, min(bi+tilenx(narr)-1, Vnx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 d1 - 1/8 d3
            if (Vtype(i,j,k)==0) then
              dV =     9 * (UV1(li,lj,lk  ) - UV1(li,lj,lk-1))
              dV = dV -    (UV3(li,lj,lk+1) - UV3(li,lj,lk-2))/3
              V2(i,j,k) = V2(i,j,k) + dV / dzmin
            else if (Vtype(i,j,k)<0) then  !near a boundary - 2nd order
              dV = coef2ord * ((V(i,j,k+1) + V(i,j,k  )) * (W(i,j+1,k  ) + W(i,j,k  )) &
                              -(V(i,j,k  ) + V(i,j,k-1)) * (W(i,j+1,k-1) + W(i,j,k-1)))
              V2(i,j,k) = V2(i,j,k) + dV / dzmin
            end if
          end do
         end do
        end do

      end do
     end do
    end do
    !$omp end do
    !$omp end parallel

    call multiply(V2, -1._knd/divcoef)

  end subroutine CD4divV
















  subroutine CD4divW(W2, U, V, W)
    !Morinishi et al., JCP 143, http://dx.doi.org/10.1006/jcph.1998.5962
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(out) :: W2
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)  :: U, V, W
    integer, parameter :: i_mask4(4) = [ -1, 9, 9, -1 ]
    integer, parameter :: divcoef = 256
    integer, parameter :: coef2ord = divcoef / 4
    integer, parameter :: narr = 4
    !UV1 and UV3 are sized to fit to the L1 cache, so no stack overflow should be possible
    real(knd) :: UV1(-1:tilenx(narr)+2, -1:tileny(narr)+2, -1:tilenz(narr)+2)
    real(knd) :: UV3(-1:tilenx(narr)+2, -1:tileny(narr)+2, -1:tilenz(narr)+2)
    integer   :: bi, bj, bk, i, j, k, li, lj, lk
    real(knd) :: Uint, Vint, Wint, dW

    call set(W2, 0._knd)

    !$omp parallel private(bi, bj, bk, i, j, k, li, lj, lk, Uint, Vint, Wint, dW, UV1, UV3)
    !$omp do
    do bk = 1, Wnz, tilenz(narr)
     do bj = 1, Wny, tileny(narr)
      do bi = 1, Wnx, tilenx(narr)

        do k = bk, min(bk+tilenz(narr)-1, Wnz)
         do j = bj, min(bj+tileny(narr)-1, Wny)
          do i = bi-2, min(bi+tilenx(narr)-1, Wnx)+1
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 Uj_1 - 1/8 Uj_3
            Uint = sum( i_mask4 * U(i,j,k-1:k+2) )
            UV1(li,lj,lk) = Uint * (W(i  ,j,k) + W(i+1,j,k))
            UV3(li,lj,lk) = Uint * (W(i-1,j,k) + W(i+2,j,k))
          end do
         end do
        end do

        do k = bk, min(bk+tilenz(narr)-1, Wnz)
         do j = bj, min(bj+tileny(narr)-1, Wny)
          do i = bi, min(bi+tilenx(narr)-1, Wnx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 d1 - 1/8 d3
            if (Wtype(i,j,k)==0) then
              dW =     9 * (UV1(li  ,lj,lk) - UV1(li-1,lj,lk))
              dW = dW -    (UV3(li+1,lj,lk) - UV3(li-2,lj,lk))/3
              W2(i,j,k) = W2(i,j,k) + dW / dxmin
            else if (Wtype(i,j,k)<0) then  !near a boundary - 2nd order
              dW = coef2ord * ((W(i+1,j,k) + W(i  ,j,k)) * (U(i  ,j,k+1) + U(i  ,j,k)) &
                              -(W(i  ,j,k) + W(i-1,j,k)) * (U(i-1,j,k+1) + U(i-1,j,k)))
              W2(i,j,k) = W2(i,j,k) + dW / dxmin
            end if
          end do
         end do
        end do

      end do
     end do
    end do
    !$omp end do

    !$omp do
    do bk = 1, Wnz, tilenz(narr)
     do bj = 1, Wny, tileny(narr)
      do bi = 1, Wnx, tilenx(narr)

        do k = bk, min(bk+tilenz(narr)-1, Wnz)
         do j = bj-2, min(bj+tileny(narr)-1, Wny)+1
          do i = bi, min(bi+tilenx(narr)-1, Wnx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 Uj_1 - 1/8 Uj_3
            Vint = sum( i_mask4 * V(i,j,k-1:k+2) )
            UV1(li,lj,lk) = Vint * (W(i,j  ,k) + W(i,j+1,k))
            UV3(li,lj,lk) = Vint * (W(i,j-1,k) + W(i,j+2,k))
          end do
         end do
        end do

        do k = bk, min(bk+tilenz(narr)-1, Wnz)
         do j = bj, min(bj+tileny(narr)-1, Wny)
          do i = bi, min(bi+tilenx(narr)-1, Wnx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 d1 - 1/8 d3
            if (Wtype(i,j,k)==0) then
              dW =     9 * (UV1(li,lj  ,lk) - UV1(li,lj-1,lk))
              dW = dW -    (UV3(li,lj+1,lk) - UV3(li,lj-2,lk))/3
              W2(i,j,k) = W2(i,j,k) + dW / dymin
            else if (Wtype(i,j,k)<0) then  !near a boundary - 2nd order
              dW = coef2ord * ((W(i,j+1,k) + W(i,j  ,k)) * (V(i,j,k+1) + V(i,j  ,k  )) &
                              -(W(i,j  ,k) + W(i,j-1,k)) * (V(i,j-1,k) + V(i,j-1,k+1)))
              W2(i,j,k) = W2(i,j,k) + dW / dymin
            end if
          end do
         end do
        end do

      end do
     end do
    end do
    !$omp end do

    !$omp do
    do bk = 1, Wnz, tilenz(narr)
     do bj = 1, Wny, tileny(narr)
      do bi = 1, Wnx, tilenx(narr)

        do k = bk-1, min(bk+tilenz(narr)-1, Wnz)+2
         do j = bj, min(bj+tileny(narr)-1, Wny)
          do i = bi, min(bi+tilenx(narr)-1, Wnx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 Uj_1 - 1/8 Uj_3
            Wint = sum( i_mask4 * W(i,j,k-2:k+1) )
            UV1(li,lj,lk) = Wint * (W(i,j,k-1) + W(i,j,k  ))
            UV3(li,lj,lk) = Wint * (W(i,j,k-2) + W(i,j,k+1))
          end do
         end do
        end do

        do k = bk, min(bk+tilenz(narr)-1, Wnz)
         do j = bj, min(bj+tileny(narr)-1, Wny)
          do i = bi, min(bi+tilenx(narr)-1, Wnx)
            li = i-bi+1
            lj = j-bj+1
            lk = k-bk+1
            !9/8 d1 - 1/8 d3
            if (Wtype(i,j,k)==0) then
              dW =     9 * (UV1(li,lj,lk+1) - UV1(li,lj,lk  ))
              dW = dW -    (UV3(li,lj,lk+2) - UV3(li,lj,lk-1))/3
              W2(i,j,k) = W2(i,j,k) + dW / dzmin
            else if (Wtype(i,j,k)<0) then  !near a boundary - 2nd order
              dW = coef2ord * ((W(i,j,k+1) + W(i,j,k  )) * (W(i,j,k+1) + W(i,j,k  )) &
                              -(W(i,j,k  ) + W(i,j,k-1)) * (W(i,j,k  ) + W(i,j,k-1)))
              W2(i,j,k) = W2(i,j,k) + dW / dzmin
            end if
          end do
         end do
        end do

      end do
     end do
    end do
    !$omp end do
    !$omp end parallel

    call multiply(W2, -1._knd/divcoef)

  end subroutine CD4divW


end module MomentumAdvection

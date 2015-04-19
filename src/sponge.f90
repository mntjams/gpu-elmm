module Sponge

  use Parameters
#ifdef MPI
    use custom_mpi
#endif

  implicit none
  
  private
  
  public :: enable_top_sponge, enable_out_sponge, &
            SpongeTop, SpongeOut, SpongeTopScalar, &
            top_sponge_bottom

  logical :: enable_top_sponge = .false.
  logical :: enable_out_sponge = .false.

  real(knd) :: top_sponge_bottom = huge(1._knd)

  real(knd), dimension(:), allocatable :: DF, avg
  
contains

  pure function DampF(x) result(res)
    real(knd) :: res
    real(knd), intent(in) :: x

    if (x<=0) then
      res = 1
    else if (x>=1) then
      res = 0
    else
      res = (1 - 0.04_knd*x**2) * &
              ( 1 - (1 - exp(10._knd*x**2)) / (1 - exp(10._knd)) )
    end if
  end function
  
  pure function ScalarDampF(x) result(res)
    real(knd) :: res
    real(knd), intent(in)::x

    if (x <= 0) then
      res = 1
    else if (x >= 1) then
      res = 0
    else
      res = (1 - 0.04_knd*x**2) * &
              ( 1 - (1 - exp(10._knd*x**2)) / (1 - exp(10._knd)) )
    end if
  end function
  

  subroutine SpongeTop(U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V, W
    integer   :: i, j, k, bufn
    real(knd) :: ze, zs, zb, p
    
    if (top_sponge_bottom<zW(Prnz)) then
      bufn = min(Prnz, int((zW(Prnz) - top_sponge_bottom) / dzmin) )
      zs = top_sponge_bottom
      ze = gzmax
   
      if (.not.allocated(DF)) allocate(DF(  min(Unz, Vnz, Wnz) - bufn  :  max(Unz, Vnz, Wnz)))
      if (.not.allocated(avg)) allocate(avg( min(Unz, Vnz, Wnz) - bufn  :  max(Unz, Vnz, Wnz)))

      !$omp parallel private(i, j, k, p, zb)
      
      !$omp do
      do k = Unz - bufn, Unz
        avg(k) = 0
      end do
      
      !$omp do
      do k = Unz - bufn, Unz
        p = 0

        do j = 1, Uny
          do i = 1, Unx
            p = p + U(i,j,k)
          end do
        end do
        avg(k) = p
      end do
      
#ifdef MPI
      avg = mpi_co_sum(avg, comm_plane_xy)
#endif

      !$omp do
      do k = Unz - bufn, Unz
        avg(k) = avg(k) / (gUnx*gUny)
      end do

      !$omp do
      do k = Unz - bufn, Unz
        zb = (zPr(k)-zs) / (ze-zs)
        DF(k) = DampF(zb)
      end do

      !$omp do
      do k = Unz - bufn, Unz
        do j = -1, Uny + 1
          do i = -1, Unx + 1
            U(i,j,k) = avg(k) + DF(k) * (U(i,j,k) - avg(k))
          end do
        end do
      end do


      !$omp do
      do k = Vnz - bufn, Vnz
        avg(k) = 0
      end do

      !$omp do
      do k = Vnz - bufn, Vnz
        p = 0

        do j = 1, Vny
          do i = 1, Vnx
            p = p + V(i,j,k)
          end do
        end do
        avg(k) = p
      end do

#ifdef MPI
      avg = mpi_co_sum(avg, comm_plane_xy)
#endif

      !$omp do
      do k = Vnz - bufn, Vnz
        avg(k) = avg(k) / (gVnx*gVny)
      end do

      !$omp do
      do k = Vnz - bufn, Vnz
        zb = (zPr(k)-zs) / (ze-zs)
        DF(k) = DampF(zb)
      end do

      !$omp do
      do k = Vnz - bufn, Vnz
        do j = -1, Vny + 1
          do i = -1, Vnx + 1
            V(i,j,k) = avg(k) + DF(k) * (V(i,j,k) - avg(k))
          end do
        end do
      end do


      !$omp do
      do k = Wnz - bufn, Wnz
        avg(k) = 0
      end do

      !$omp do
      do k = Wnz - bufn, Wnz
        p = 0

        do j = 1, Wny
          do i = 1, Wnx
            p = p + W(i,j,k)
          end do
        end do
        avg(k) = p
      end do

#ifdef MPI
      avg = mpi_co_sum(avg, comm_plane_xy)
#endif

      !$omp do
      do k = Wnz - bufn, Wnz
        avg(k) = avg(k) / (gWnx*gWny)
      end do

      !$omp do
      do k = Wnz - bufn, Wnz
        zb = (zW(k)-zs) / (ze-zs)
        DF(k) = DampF(zb)
      end do

      !$omp do
      do k = Wnz - bufn, Wnz
        do j = -1, Wny + 1
          do i = -1, Wnx + 1
            W(i,j,k) = avg(k) + DF(k) * (W(i,j,k) - avg(k))
          end do
        end do
      end do
      
      !$omp end parallel
    
    end if

  end subroutine SpongeTop
  
  
  
  subroutine SpongeTopScalar(Phi)
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(inout) :: Phi
    integer :: i, j, k, bufn
    real(knd) :: ze, zs, zb, p

    if (top_sponge_bottom<zW(Prnz)) then

      bufn = min(Prnz, int((zW(Prnz) - top_sponge_bottom) / dzmin) )
      zs = top_sponge_bottom
      ze = gzmax

      if (.not.allocated(DF)) allocate(DF(  min(Unz, Vnz, Wnz) - bufn  :  max(Unz, Vnz, Wnz)))
      if (.not.allocated(avg)) allocate(avg( min(Unz, Vnz, Wnz) - bufn  :  max(Unz, Vnz, Wnz)))

      !$omp parallel private(i,j,k,p,zb)
      
      !$omp do
      do k = Prnz-bufn, Prnz
        avg(k) = 0
      end do
      !$omp end do

      !$omp do
      do k = Prnz-bufn, Prnz
        p = 0
        do j = 1, Prny
          do i = 1, Prnx
            p = p + Phi(i,j,k)
          end do
        end do
        avg(k) = p
      end do
      !$omp end do

#ifdef MPI
        avg = mpi_co_sum(avg, comm_plane_xy)
#endif

      !$omp do
      do k = Prnz-bufn, Prnz
        avg(k) = avg(k) / (gPrnx*gPrny)
      end do
      !$omp end do

      !$omp do
      do k = Prnz-bufn, Prnz
        zb = (zPr(k)-zs) / (ze-zs)
        DF(k) = ScalarDampF(zb)
      end do
      !$omp end do

      !$omp do
      do k = Prnz-bufn, Prnz
        do j = -1, Prny+1
          do i = -1, Prnx+1
            Phi(i,j,k) = avg(k) + DF(k) * (Phi(i,j,k)-avg(k))
          end do
        end do
      end do
      !$omp end do

      !$omp end parallel
    end if

  endsubroutine SpongeTopScalar
  
  
  
  

  subroutine SpongeOut(U, V, W, temperature)
    real(knd), contiguous, intent(inout), dimension(-2:,-2:,-2:) :: U, V, W
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(inout) :: temperature
    integer   :: i, j, k, bufn
    real(knd) :: p, xe, xs, xb, DF

    bufn = min(max(10, Prnx/8), Prnx/2)
    xs = xU(Prnx - bufn)
    xe = xU(Prnx)

    !$omp parallel private(i, j, k, p, xb, DF)

    !$omp do
    do k = 1, Unz
      do j = 1, Uny
        p = 0
        do i = 2*Unx/3, Unx - 4
          p = p + U(i,j,k)
        end do
        p = p / (Unx - 4 - 2*Unx/3 + 1)
        do i = Unx - bufn, Unx + 1
          xb = (xU(i)-xs) / (xe-xs)
          DF = DampF(xb)
          U(i,j,k) = p + DF * (U(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do k = 1, Vnz
      do j = 1, Vny
        p = 0
        do i = 2*Vnx/3, Vnx - 4
          p = p + V(i,j,k)
        end do
        p = p / (Vnx - 4 - 2*Vnx/3 + 1)
        do i = Vnx-bufn, Vnx + 1
          xb = (xPr(i)-xs) / (xe-xs)
          DF = DampF(xb)
          V(i,j,k) = p + DF * (V(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    !$omp do
    do k = 1, Wnz
      do j = 1, Wny
        p = 0
        do i = 2*Wnx/3, Wnx - 4
          p = p + W(i,j,k)
        end do
        p = p / (Wnx - 4 - 2*Wnx/3 + 1)
        do i = Wnx - bufn, Wnx + 1
          xb = (xPr(i)-xs) / (xe-xs)
          DF = DampF(xb)
          W(i,j,k) = p + DF * (W(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    if (enable_buoyancy) then
      !$omp do
      do k = 1, Prnz
        do j = 1, Prny
          p = 0
          do i = 2*Prnx/3, Prnx-4
            p = p + temperature(i,j,k)
          end do
          p = p / (Prnx - 4 - 2*Prnx/3 + 1)
          do i = Prnx - bufn, Prnx + 1
            xb = (xPr(i)-xs) / (xe-xs)
            DF = DampF(xb)
            temperature(i,j,k) = p + DF * (temperature(i,j,k) - p)
          end do
        end do
      end do
      !$omp end do
    end if
    !$omp end parallel

  end subroutine SpongeOut



end module
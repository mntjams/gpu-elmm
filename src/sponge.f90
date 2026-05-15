module Sponge

  use Parameters
#ifdef PAR
    use custom_par
#endif

  implicit none
  
  private
  
  public :: enable_top_sponge, enable_top_sponge_scalar, &
            enable_in_sponge_x, enable_out_sponge_x, enable_out_sponge_y, &
            SpongeTop, SpongeOut, SpongeTopScalar, &
            top_sponge_bottom, sponge_to_profiles, &
            U_sponge_avg, V_sponge_avg, W_sponge_avg, Temperature_sponge_avg, Moisture_sponge_avg, &
            rayleigh_damping_uvw, rayleigh_damping_temperature, rayleigh_damping_moisture, &
            init_rayleigh_damping_top, rayleigh_damping_top_properties

  logical :: enable_top_sponge = .false.
  logical :: enable_top_sponge_scalar = .false.
  logical :: enable_in_sponge_x = .false.
  logical :: enable_out_sponge_x = .false.
  logical :: enable_out_sponge_y = .false.
  logical :: enable_top_rayleigh_damping = .false.

  logical :: sponge_to_profiles = .false.

  real(knd), dimension(:), allocatable :: U_sponge_avg, V_sponge_avg, W_sponge_avg, Temperature_sponge_avg, Moisture_sponge_avg
  
  type damping_properties
    real(knd) :: bottom, top !global constants for all images
    real(knd), dimension(:), allocatable :: coefficient(:)
    real(knd), dimension(:), allocatable :: coefficientW(:)
    integer :: kbottom = 1, ktop = 0
    integer :: kWbottom = 1, kWtop = 0
    real(knd), dimension(:), allocatable :: u_prof, v_prof, w_prof, temperature_prof, moisture_prof
  end type
  
  type(damping_properties) :: rayleigh_damping_top_properties

  real(knd) :: top_sponge_bottom = huge(1._knd)

  real(knd), dimension(:), allocatable :: DF, avg
  
contains

  elemental function DampF(x) result(res)
    real(knd) :: res
    real(knd), intent(in) :: x

    if (x<=0) then
      res = 1
    else if (x>=1) then
      res = 0
    else
      res = (1 - 0.1_knd*x**2) * &
              ( 1 - (1 - exp(10._knd*x**2)) / (1 - exp(10._knd)) )
    end if
  end function
  
  elemental function ScalarDampF(x) result(res)
    real(knd) :: res
    real(knd), intent(in)::x

    if (x <= 0) then
      res = 1
    else if (x >= 1) then
      res = 0
    else
      res = (1 - 0.1_knd*x**2) * &
              ( 1 - (1 - exp(10._knd*x**2)) / (1 - exp(10._knd)) )
    end if
  end function
  

  subroutine SpongeTop(U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V, W
    integer   :: i, j, k, lo_U, lo_V, lo_W
    real(knd) :: ze, zs, zb, p
    
    if (top_sponge_bottom<zW(Prnz)) then
      lo_U = Unz - min(Unz, int((zW(Prnz) - top_sponge_bottom) / dzmin) ) + 1
      lo_V = Vnz - min(Vnz, int((zW(Prnz) - top_sponge_bottom) / dzmin) ) + 1
      lo_W = Wnz - min(Wnz, int((zW(Prnz) - top_sponge_bottom) / dzmin) ) + 1
      zs = top_sponge_bottom
      ze = gzmax
   
      if (.not.allocated(DF)) allocate(DF(  min(lo_U, lo_V, lo_W)  :  max(Unz, Vnz, Wnz)))
      if (.not.allocated(avg)) allocate(avg( min(lo_U, lo_V, lo_W)  :  max(Unz, Vnz, Wnz)))

      !$omp parallel private(i, j, k, p, zb)
      
      if (sponge_to_profiles) then
        avg(lo_U:Unz) = U_sponge_avg(lo_U:Unz)
      else
        !$omp do
        do k = lo_U, Unz
          avg(k) = 0
        end do
        
        !$omp do
        do k = lo_U, Unz
          p = 0

          do j = 1, Uny
            do i = 1, Unx
              p = p + U(i,j,k)
            end do
          end do
          avg(k) = p
        end do
        
#ifdef PAR
        avg = par_co_sum(avg, comm_plane_xy)
#endif

        !$omp do
        do k = lo_U, Unz
          avg(k) = avg(k) / (gUnx*gUny)
        end do
      end if

      !$omp do
      do k = lo_U, Unz
        zb = (zPr(k)-zs) / (ze-zs)
        DF(k) = DampF(zb)
      end do

      !$omp do
      do k = lo_U, Unz
        do j = -1, Uny + 1
          do i = -1, Unx + 1
            U(i,j,k) = avg(k) + DF(k) * (U(i,j,k) - avg(k))
          end do
        end do
      end do


      if (sponge_to_profiles) then
        avg(lo_V:Vnz) = V_sponge_avg(lo_V:Vnz)
      else
        !$omp do
        do k = lo_V, Vnz
          avg(k) = 0
        end do

        !$omp do
        do k = lo_V, Vnz
          p = 0

          do j = 1, Vny
            do i = 1, Vnx
              p = p + V(i,j,k)
            end do
          end do
          avg(k) = p
        end do

#ifdef PAR
        avg = par_co_sum(avg, comm_plane_xy)
#endif

        !$omp do
        do k = lo_V, Vnz
          avg(k) = avg(k) / (gVnx*gVny)
        end do
      end if

      !$omp do
      do k = lo_V, Vnz
        zb = (zPr(k)-zs) / (ze-zs)
        DF(k) = DampF(zb)
      end do

      !$omp do
      do k = lo_V, Vnz
        do j = -1, Vny + 1
          do i = -1, Vnx + 1
            V(i,j,k) = avg(k) + DF(k) * (V(i,j,k) - avg(k))
          end do
        end do
      end do


      if (sponge_to_profiles) then
        avg(lo_W:Wnz) = W_sponge_avg(lo_W:Wnz)
      else
        !$omp do
        do k = lo_W, Wnz
          avg(k) = 0
        end do

        !$omp do
        do k = lo_W, Wnz
          p = 0

          do j = 1, Wny
            do i = 1, Wnx
              p = p + W(i,j,k)
            end do
          end do
          avg(k) = p
        end do

#ifdef PAR
        avg = par_co_sum(avg, comm_plane_xy)
#endif

        !$omp do
        do k = lo_W, Wnz
          avg(k) = avg(k) / (gWnx*gWny)
        end do
      end if

      !$omp do
      do k = lo_W, Wnz
        zb = (zW(k)-zs) / (ze-zs)
        DF(k) = DampF(zb)
      end do

      !$omp do
      do k = lo_W, Wnz
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
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: Phi
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

#ifdef PAR
        avg = par_co_sum(avg, comm_plane_xy)
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
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: temperature
    
    !NOTE: currently having both of them enabled will likely lead to strange results
    if (enable_in_sponge_x) call SpongeIn_X(U, V, W, temperature)
    if (enable_out_sponge_x) call SpongeOut_X(U, V, W, temperature)
    if (enable_out_sponge_y) call SpongeOut_Y(U, V, W, temperature)
  end subroutine

  subroutine SpongeIn_X(U, V, W, temperature)
    real(knd), contiguous, intent(inout), dimension(-2:,-2:,-2:) :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: temperature
    integer   :: i, j, k, bufn
    integer   :: hi, lo, loU, hiU
    real(knd) :: p, xe, xs, xb, DF

    !size end extent of the buffer region
    bufn = min(max(5,Prnx/50), Prnx/4)
    xs = xU(0)
    xe = xU(bufn + 2)

    !extent of the probe region where the local average is taken from
    loU = bufn+1
    hiU = max(Unx/3, 50)

    lo = bufn+1
    hi = max(Prnx/3, 50)

    !$omp parallel private(i, j, k, p, xb, DF)

    !$omp do collapse(2)
    do k = 1, Unz
      do j = 1, Uny
        p = 0
        do i = loU, hiU
          p = p + U(i,j,k)
        end do
        p = p / (hiU - loU + 1)
        do i = Unx - bufn, Unx + 1
          xb = (xe-xU(i)) / (xe-xs)
          DF = DampF(xb)
          U(i,j,k) = p + DF * (U(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    !$omp do collapse(2)
    do k = 1, Vnz
      do j = 1, Vny
        p = 0
        do i = lo, hi
          p = p + V(i,j,k)
        end do
        p = p / (hi - lo + 1)
        do i = Vnx-bufn, Vnx + 1
          xb = (xe-xPr(i)) / (xe-xs)
          DF = DampF(xb)
          V(i,j,k) = p + DF * (V(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    !$omp do collapse(2)
    do k = 1, Wnz
      do j = 1, Wny
        p = 0
        do i = lo, hi
          p = p + W(i,j,k)
        end do
        p = p / (hi - lo + 1)
        do i = Wnx - bufn, Wnx + 1
          xb = (xe-xPr(i)) / (xe-xs)
          DF = DampF(xb)
          W(i,j,k) = p + DF * (W(i,j,k) - p)
        end do
      end do
    end do
    !$omp end do

    if (enable_buoyancy) then
      !$omp do collapse(2)
      do k = 1, Prnz
        do j = 1, Prny
          p = 0
          do i = lo, hi
            p = p + temperature(i,j,k)
          end do
          p = p / (hi - lo + 1)
          do i = Prnx - bufn, Prnx + 1
            xb = (xe-xPr(i)) / (xe-xs)
            DF = DampF(xb)
            temperature(i,j,k) = p + DF * (temperature(i,j,k) - p)
          end do
        end do
      end do
      !$omp end do
    end if
    !$omp end parallel

  end subroutine SpongeIn_X


  subroutine SpongeOut_X(U, V, W, temperature)
    real(knd), contiguous, intent(inout), dimension(-2:,-2:,-2:) :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: temperature
    integer   :: i, j, k, bufn
    integer   :: hi, lo, loU, hiU
    real(knd) :: p, xe, xs, xb, DF

    !size end extent of the buffer region
    bufn = min(max(5,Prnx/50), Prnx/4)
    xs = xU(Prnx - bufn-2)
    xe = xU(Prnx)

    !extent of the probe region where the local average is taken from
    loU = max(2*Unx/3, Unx-50)
    hiU = Unx-bufn-1

    lo = max(2*Prnx/3, Prnx-50)
    hi = Prnx-bufn-1

    !$omp parallel private(i, j, k, p, xb, DF)

    !$omp do collapse(2)
    do k = 1, Unz
Uj:   do j = 1, Uny
        p = 0
        do i = loU, hiU
          if (Utype(i,j,k)>0) cycle Uj
          p = p + U(i,j,k)
        end do
        p = p / (hiU - loU + 1)
        do i = Unx - bufn, Unx + 1
          xb = (xU(i)-xs) / (xe-xs)
          DF = DampF(xb)
          U(i,j,k) = p + DF * (U(i,j,k) - p)
        end do
      end do Uj
    end do
    !$omp end do

    !$omp do collapse(2)
    do k = 1, Vnz
Vj: do j = 1, Vny
        p = 0
        do i = lo, hi
          if (Vtype(i,j,k)>0) cycle Vj
          p = p + V(i,j,k)
        end do
        p = p / (hi - lo + 1)
        do i = Vnx-bufn, Vnx + 1
          xb = (xPr(i)-xs) / (xe-xs)
          DF = DampF(xb)
          V(i,j,k) = p + DF * (V(i,j,k) - p)
        end do
      end do Vj
    end do
    !$omp end do

    !$omp do collapse(2)
    do k = 1, Wnz
Wj: do j = 1, Wny
        p = 0
        do i = lo, hi
          if (Wtype(i,j,k)>0) cycle Wj
          p = p + W(i,j,k)
        end do
        p = p / (hi - lo + 1)
        do i = Wnx - bufn, Wnx + 1
          xb = (xPr(i)-xs) / (xe-xs)
          DF = DampF(xb)
          W(i,j,k) = p + DF * (W(i,j,k) - p)
        end do
      end do Wj
    end do
    !$omp end do

    if (enable_buoyancy) then
      !$omp do collapse(2)
      do k = 1, Prnz
Prj:    do j = 1, Prny
          p = 0
          do i = lo, hi
            if (Prtype(i,j,k)>0) cycle Prj
            p = p + temperature(i,j,k)
          end do
          p = p / (hi - lo + 1)
          do i = Prnx - bufn, Prnx + 1
            xb = (xPr(i)-xs) / (xe-xs)
            DF = DampF(xb)
            temperature(i,j,k) = p + DF * (temperature(i,j,k) - p)
          end do
        end do Prj
      end do
      !$omp end do
    end if
    !$omp end parallel

  end subroutine SpongeOut_X


  subroutine SpongeOut_Y(U, V, W, temperature)
    use custom_par
    real(knd), contiguous, intent(inout), dimension(-2:,-2:,-2:) :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: temperature
    integer   :: i, j, k, bufn
    integer   :: hi, lo, loV, hiV
    real(knd) :: p, ye, ys, yb, DF

    !TODO: Properly parallelize for a sponge zone across several images
    if (jim==nyims) then
      !size end extent of the buffer region
      bufn = min(max(5,Prnx/50), gPrny/4, Prny/4)
      ys = yV(Prny - bufn-2)
      ye = yV(Prny)

      !extent of the probe region where the local average is taken from
      loV = max(2*gPrny/3-offset_to_global_y, gPrny-50-offset_to_global_y, 1)
      hiV = Vny-bufn-1

      lo = max(2*gPrny/3-offset_to_global_y, gPrny-50-offset_to_global_y, 1)
      hi = Prny-bufn-1

      !$omp parallel private(i, j, k, p, yb, DF)

      !$omp do collapse(2)
      do k = 1, Unz
        do i = 1, Unx
          p = 0
          do j = lo, hi
            p = p + U(i,j,k)
          end do
          p = p / (hi - lo + 1)
          do j = Uny - bufn, Uny + 1
            yb = (yPr(j)-ys) / (ye-ys)       
            DF = DampF(yb)
            U(i,j,k) = p + DF * (U(i,j,k) - p)
          end do
        end do
      end do
      !$omp end do

      !$omp do collapse(2)
      do k = 1, Vnz
        do i = 1, Vnx
          p = 0
          do j = loV, hiV
            p = p + V(i,j,k)
          end do
          p = p / (hiV - loV + 1)
          do j = Vny-bufn, Vny + 1
            yb = (yV(j)-ys) / (ye-ys)
            DF = DampF(yb)
            V(i,j,k) = p + DF * (V(i,j,k) - p)
          end do           
        end do
      end do
      !$omp end do

      !$omp do collapse(2)
      do k = 1, Wnz
        do i = 1, Wnx
          p = 0
          do j = lo, hi
            p = p + W(i,j,k)
          end do
          p = p / (hi - lo + 1)
          do j = Wny - bufn, Wny + 1
            yb = (yPr(j)-ys) / (ye-ys)
            DF = DampF(yb)
            W(i,j,k) = p + DF * (W(i,j,k) - p)
          end do
        end do
      end do
      !$omp end do

      if (enable_buoyancy) then
        !$omp do collapse(2)
        do k = 1, Prnz
          do i = 1, Prny
            p = 0
            do j = lo, hi
              p = p + temperature(i,j,k)
            end do
            p = p / (hi - lo + 1)
            do j = Prny - bufn, Prny + 1
              yb = (yPr(j)-ys) / (ye-ys)
              DF = DampF(yb)
              temperature(i,j,k) = p + DF * (temperature(i,j,k) - p)
            end do
          end do
        end do
        !$omp end do
      end if
      !$omp end parallel
    end if

  end subroutine SpongeOut_Y
  
  
  
  
  
  
  subroutine rayleigh_damping_uvw(U2, V2, W2, U, V, W)
    real(knd), contiguous, intent(inout) :: U2(-2:,-2:,-2:), V2(-2:,-2:,-2:), W2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
  
    if (enable_top_rayleigh_damping) call rayleigh_damping_top_uvw(U2, V2, W2, U, V, W)
  end subroutine
  
  subroutine rayleigh_damping_temperature(Temperature2, Temperature)
    real(knd), contiguous, intent(inout) :: Temperature2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: Temperature(-2:,-2:,-2:)
  
    if (enable_top_rayleigh_damping) call rayleigh_damping_top_temperature(Temperature2, Temperature)
  end subroutine
  
  subroutine rayleigh_damping_moisture(Moisture2, Moisture)
    real(knd), contiguous, intent(inout) :: Moisture2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: Moisture(-2:,-2:,-2:)
  
    if (enable_top_rayleigh_damping) call rayleigh_damping_top_moisture(Moisture2, Moisture)
  end subroutine
  
  
  subroutine rayleigh_damping_top_uvw(U2, V2, W2, U, V, W)
    !an alternative formulation of the top sponge as a tendency
    real(knd), contiguous, intent(inout) :: U2(-2:,-2:,-2:), V2(-2:,-2:,-2:), W2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    integer :: i, j, k
  
    associate(prop => rayleigh_damping_top_properties)
    
      do k = prop%kbottom, prop%ktop
        do j = 1, Uny
          do i = 1, Unx
            U2(i,j,k) = U2(i,j,k) - prop%coefficient(k) * (U(i,j,k) - prop%u_prof(k))
          end do
        end do      
      end do
    
      do k = prop%kbottom, prop%ktop
        do j = 1, Vny
          do i = 1, Vnx
            V2(i,j,k) = V2(i,j,k) - prop%coefficient(k) * (V(i,j,k) - prop%v_prof(k))
          end do
        end do      
      end do
    
      do k = prop%kWbottom, prop%kWtop
        do j = 1, Wny
          do i = 1, Wnx
            W2(i,j,k) = W2(i,j,k) - prop%coefficientW(k) * (W(i,j,k) - prop%w_prof(k))
          end do
        end do      
      end do
    
    end associate
    
  end subroutine

  subroutine rayleigh_damping_top_Temperature(Temperature2, Temperature)
    !an alternative formulation of the top sponge as a tendency
    real(knd), contiguous, intent(inout) :: Temperature2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: Temperature(-2:,-2:,-2:)
    integer :: i, j, k
  
    associate(prop => rayleigh_damping_top_properties)
    
      do k = prop%kbottom, prop%ktop
        do j = 1, Uny
          do i = 1, Unx
            Temperature2(i,j,k) = Temperature2(i,j,k) - prop%coefficient(k) * (Temperature(i,j,k) - prop%temperature_prof(k))
          end do
        end do      
      end do
    
    end associate
    
  end subroutine

  subroutine rayleigh_damping_top_Moisture(Moisture2, Moisture)
    !an alternative formulation of the top sponge as a tendency
    real(knd), contiguous, intent(inout) :: Moisture2(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: Moisture(-2:,-2:,-2:)
    integer :: i, j, k
  
    associate(prop => rayleigh_damping_top_properties)
    
      do k = prop%kbottom, prop%ktop
        do j = 1, Uny
          do i = 1, Unx
            Moisture2(i,j,k) = Moisture2(i,j,k) - prop%coefficient(k) * (Moisture(i,j,k) - prop%moisture_prof(k))
          end do
        end do      
      end do
    
    end associate
    
  end subroutine
  
  
  
  subroutine init_rayleigh_damping_top(bottom, top, coefficient_fun)
    real(knd), intent(in) :: bottom, top    
    interface
      function coefficient_fun(z) result(res)
        use Kinds
        real(knd) :: res
        real(knd), intent(in) :: z
      end function
    end interface
    
    integer :: k
    
    if (top>bottom .and. top>zW(0) .and. bottom<zW(Prnz)) then
      associate(prop => rayleigh_damping_top_properties)
      
        enable_top_rayleigh_damping = .true.
      
        prop%bottom = bottom
        prop%top = top
        
        do k = 1, Prnz
          if (zPr(k)>=bottom) then
            if (prop%kbottom==0) prop%kbottom = k
            if (zPr(k) <= top) prop%ktop = k
          end if
        end do
        
        if (prop%ktop>=prop%kbottom) allocate(prop%coefficient(prop%kbottom:prop%ktop))
        
        do k = prop%kbottom, prop%ktop
          prop%coefficient(k) = coefficient_fun(zPr(k))
        end do
        
        do k = 1, Wnz
          if (zW(k)>=bottom) then
            if (prop%kWbottom==0) prop%kWbottom = k
            if (zW(k) <= top) prop%kWtop = k
          end if
        end do
        
        if (prop%kWtop>=prop%kWbottom) allocate(prop%coefficientW(prop%kWbottom:prop%kWtop))
        
        do k = prop%kWbottom, prop%kWtop
          prop%coefficientW(k) = coefficient_fun(zW(k))
        end do
        
        ! currently to be set by external code that calls this subroutine
        allocate(prop%u_prof(prop%kbottom:prop%ktop))
        allocate(prop%v_prof(prop%kbottom:prop%ktop))
        allocate(prop%w_prof(prop%kWbottom:prop%kWtop))
        prop%u_prof = 0
        prop%v_prof = 0
        prop%w_prof = 0
        if (enable_buoyancy) then
          allocate(prop%temperature_prof(prop%kbottom:prop%ktop))
          prop%temperature_prof = temperature_ref
        else
          allocate(prop%temperature_prof(1:0))
        end if
        if (enable_moisture) then
          allocate(prop%moisture_prof(prop%kbottom:prop%ktop))
          prop%moisture_prof = moisture_ref
        else
          allocate(prop%moisture_prof(1:0))        
        end if
      end associate
    else
      associate(prop => rayleigh_damping_top_properties)
        allocate(prop%u_prof(1:0))
        allocate(prop%v_prof(1:0))
        allocate(prop%w_prof(1:0))
        allocate(prop%temperature_prof(1:0))
        allocate(prop%moisture_prof(1:0))  
      end associate
    end if
  end subroutine


end module

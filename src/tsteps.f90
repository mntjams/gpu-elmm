module TimeSteps

  use Parameters
  use Dynamics
  use Boundaries, only: BoundU,Bound_Q
  use Pressure, only: PressureCorrection
  use Outputs, only: store, display, current_profiles
  use Scalars, only: ScalarRK3
  use Turbinlet, only: GetTurbulentInlet, GetInletFromFile

  implicit none


  private
  public TMarchRK3


contains


  subroutine TMarchRK3(U, V, W, Pr, Temperature, Moisture, Scalar, dt, delta)
    use RK3
    real(knd), allocatable, intent(inout) :: U(:,:,:), V(:,:,:) ,W(:,:,:), Pr(:,:,:)
    real(knd), allocatable, intent(inout) :: Temperature(:,:,:), Moisture(:,:,:), Scalar(:,:,:,:)
    real(knd), intent(out) :: dt, delta

    integer :: RK_stage
    integer, save :: called = 0
    integer(int64), save :: trate
    integer(int64), save :: time1, time2


    if (called==0) then
      called = 1

      !just to allocate it and make it defined in all points
      U2 = U
      V2 = V
      W2 = W

      Ustar = U
      Vstar = V
      Wstar = W


      call BoundU(1, U, Uin)

      call BoundU(2, V, Vin)

      call BoundU(3, W, Win)

      if (enable_buoyancy) call BoundTemperature(temperature)

      call IBMomentum(U, V, W)


      if (masssourc==1) allocate(Q(0:Prnx+1, 0:Prny+1, 0:Prnz+1))

      if (debugparam>1) call system_clock(count_rate=trate)
    end if


    if ((Btype(We)==TurbulentInlet) .or. (Btype(Ea)==TurbulentInlet)) then
      call GetTurbulentInlet(dt)
    else if (Btype(We)==InletFromFile) then
      call GetInletFromFile(time)
    end if


    call TimeStepLength(U, V, W, dt)



    if (master) write (*,'(a,f12.6,a,es12.4)') " time: ", time," dt: ", dt


    do RK_stage = 1, RK_stages

      if (debugparam>1.and.master) call system_clock(count=time1)


      call SubgridStresses(U, V, W, Pr, Temperature)


      call Convection(U, V, W, &
                      U2, V2, W2, &
                      Ustar, Vstar, Wstar, &
                      Temperature, Moisture, &
                      RK_beta, RK_rho, RK_stage, dt)

      call ScalarRK3(U, V, W, &
                     Temperature, Moisture, Scalar, &
                     RK_stage, dt, &
                     current_profiles%tempfl, current_profiles%moistfl)

      call OtherTerms(U, V, W, &
                      U2, V2, W2, &
                      Pr, &
                      2*RK_alpha(RK_stage)*dt)


      if ((Btype(To) ==FreeSlipBuff) .and. (Prnz>15))  then

          call AttenuateTop(U2, V2, W2)
      end if

      if ((Btype(Ea) ==OutletBuff) .and. (Prnx>15)) then

          call AttenuateOut(U2, V2, W2, temperature)
      end if



      call BoundU(1, U2, Uin)

      call BoundU(2, V2, Vin)

      call BoundU(3, W2, Win)

      call IBMomentum(U2, V2, W2)

      if (masssourc==1) then
          call IBMassSources(Q, U2, V2, W2)
      end if


      if (poissmet>0) then
        call PressureCorrection(U2, V2, W2, Pr, Q, 2*RK_alpha(RK_stage)*dt)
      end if



      if (RK_stage==1) delta = 0

#ifdef DEBUG
      if ( debuglevel>0 .or. steady==1 ) then

        if (Unx*Uny*Unz > 0) &
          delta = delta + sum(abs(U(1:Unx,1:Uny,1:Unz) - U2(1:Unx,1:Uny,1:Unz))) / (Unx*Uny*Unz)

        if (Vnx*Vny*Vnz > 0) &
          delta = delta + sum(abs(V(1:Vnx,1:Vny,1:Vnz) - V2(1:Vnx,1:Vny,1:Vnz))) / (Vnx*Vny*Vnz)

        if (Wnx*Wny*Wnz > 0) &
          delta = delta + sum(abs(W(1:Wnx,1:Wny,1:Wnz) - W2(1:Wnx,1:Wny,1:Wnz))) / (Wnx*Wny*Wnz)

      end if
#endif

      call exchange_alloc(U, U2)
      call exchange_alloc(V, V2)
      call exchange_alloc(W, W2)



      if ((Btype(Ea) ==OutletBuff) .and. (Prnx>15)) then

        call AttenuateOut(U, V, W, temperature)

      end if


      call NullInterior(U, V, W)


      call BoundU(1, U, Uin)

      call BoundU(2, V, Vin)

      call BoundU(3, W, Win)

      call IBMomentum(U, V, W)


      if (debugparam>1.and.master) then
        call system_clock(count=time2)
        write (*,*) "ET of part 1", (time2-time1) / real(trate, int64)
        time1 = time2
      end if

    end do


  end subroutine TMarchRK3



















  subroutine OtherTerms(U, V, W, U2, V2, W2, Pr, coef)
    real(knd), contiguous,  intent(inout) :: U(-2:,-2:,-2:), V(-2:,-2:,-2:), W(-2:,-2:,-2:)
    real(knd), contiguous,  intent(inout) :: Pr(1:,1:,1:)
    real(knd), allocatable, intent(inout) :: U2(:,:,:), V2(:,:,:), W2(:,:,:)
    real(knd), intent(in) :: coef

    real(knd) :: S

    integer :: i, j, k
    integer, save :: called=0

    if (called==0) then
      allocate(U3(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),lbound(U,3):ubound(U,3)))
      allocate(V3(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),lbound(V,3):ubound(V,3)))
      allocate(W3(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),lbound(W,3):ubound(W,3)))
      called=1
    end if


    call PressureGrad(Pr, U2, V2, W2, coef)

    if (.not.explicit_diffusion) then

        !semi-implicit diffusion

        Re_gt_0: if (Re>0) then

          !Diffusion using Crank Nicolson
          !first approximation using forward Euler
          !iteration SOR or Gauss-Seidel

          call ImplicitDiffusion_ForwEul(U, V, W, U2, V2, W2, U3, V3, W3, coef)


          !$omp parallel sections
          !$omp section
          call BoundU(1, U3, Uin)
          !$omp section
          call BoundU(2, V3, Vin)
          !$omp section
          call BoundU(3, W3, Win)
          !$omp end parallel sections
          call IBMomentum(U3, V3, W3)

          !Performs the diffusion terms

          if (gridtype==UNIFORMGRID) then
            call ImplicitDiffusion_Iterations(U, V, W, U2, V2, W2, U3, V3, W3, coef)
          else
            call error_stop("Non-uniform grid support was dropped.")
          end if

          call exchange_alloc(U2, U3)
          call exchange_alloc(V2, V3)
          call exchange_alloc(W2, W3)


        else  Re_gt_0  !Re<=0

          U2 = U + U2
          V2 = V + V2
          W2 = W + W2

        end if   Re_gt_0

        !$omp parallel sections
        !$omp section
        call BoundU(1, U2, Uin)
        !$omp section
        call BoundU(2, V2, Vin)
        !$omp section
        call BoundU(3, W2, Win)
        !$omp end parallel sections
        call IBMomentum(U2, V2, W2)

        if (debuglevel>=2) then  !Compute and output the mean friction in the domain.
          S = 0
          do k = 1, Unz
           do j = 1, Uny
            do i = 1, Unx
              S = S-((Viscosity(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) / dxPr(i+1) - &
              Viscosity(i,j,k) * (U(i,j,k)-U(i-1,j,k)) / dxPr(i)) / dxU(i) + &
                (0.25_knd * (Viscosity(i+1,j+1,k)+Viscosity(i+1,j,k)+Viscosity(i,j+1,k)+Viscosity(i,j,k))* &
                      (U(i,j+1,k)-U(i,j,k)) / dyV(j) - &
                0.25_knd * (Viscosity(i+1,j,k)+Viscosity(i+1,j-1,k)+Viscosity(i,j,k)+Viscosity(i,j-1,k))* &
                      (U(i,j,k)-U(i,j-1,k)) / dyV(j-1)) / dyPr(j) + &
                 (0.25_knd * (Viscosity(i+1,j,k+1)+Viscosity(i+1,j,k)+Viscosity(i,j,k+1)+Viscosity(i,j,k))* &
                      (U(i,j,k+1)-U(i,j,k)) / dzW(k) - &
                0.25_knd * (Viscosity(i+1,j,k)+Viscosity(i+1,j,k-1)+Viscosity(i,j,k)+Viscosity(i,j,k-1))* &
                      (U(i,j,k)-U(i,j,k-1)) / dzW(k-1)) / dzPr(k))
            end do
           end do
          end do

          S = S / (Unx*Uny*Unz)
          write(*,*) "Mean friction:", S
        end if

    else !explicit diffusion

        call add(U2, U)
        call add(V2, V)
        call add(W2, W)

    end if

  end subroutine OtherTerms








  
  subroutine IBMomentum(U, V, W)
    use ImmersedBoundary, only: Up => UIBPoints, &
                                Vp => VIBPoints, &
                                Wp => WIBPoints, &
                                Interpolate => TIBPoint_Interpolate

    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V, W
    integer :: i

    if (size(Up) + size(Vp) + size(Wp)>0) then
      if (.not. allocated(Uwm)) allocate(Uwm(size(Up)))
      if (.not. allocated(Vwm)) allocate(Vwm(size(Vp)))
      if (.not. allocated(Wwm)) allocate(Wwm(size(Wp)))

      !$omp parallel
      !$omp do
      do i = 1, ubound(Up, 1)
          Uwm(i) = Interpolate(Up(i), U, -2)
      end do
      !$omp end do
      !$omp do
      do i = 1, ubound(Up, 1)
          U(Up(i)%xi, Up(i)%yj, Up(i)%zk) = Uwm(i)
      end do
      !$omp end do nowait

      !$omp do
      do i = 1, ubound(Vp, 1)
          Vwm(i) = Interpolate(Vp(i), V, -2)
      end do
      !$omp end do
      !$omp do
      do i = 1, ubound(Vp, 1)
          V(Vp(i)%xi, Vp(i)%yj, Vp(i)%zk) = Vwm(i)
      end do
      !$omp end do nowait

      !$omp do
      do i = 1, ubound(Wp, 1)
          Wwm(i) = Interpolate(Wp(i), W, -2)
      end do
      !$omp end do
      !$omp do
      do i = 1, ubound(Wp, 1)
          W(Wp(i)%xi, Wp(i)%yj, Wp(i)%zk) = Wwm(i)
      end do
      !$omp end do
      !$omp end parallel

    end if
  end subroutine IBMomentum






  subroutine IBMassSources(Q, U, V, W)
    use ImmersedBoundary, only: Up => UIBPoints, &
                                Vp => VIBPoints, &
                                Wp => WIBPoints
    use vtkarray

    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
    real(knd), intent(out) :: Q(0:,0:,0:)
    integer :: i, xi, yj, zk

    Q = 0

    !$omp parallel
    !$omp do private(xi, yj, zk)
    do i = 1, ubound(Up, 1)
      xi = Up(i)%xi
      yj = Up(i)%yj
      zk = Up(i)%zk
      !$omp atomic
      Q(xi,yj,zk)   = Q(xi,yj,zk)   + U(xi,yj,zk) / dxPr(xi)
      !$omp atomic
      Q(xi+1,yj,zk) = Q(xi+1,yj,zk) - U(xi,yj,zk) / dxPr(xi+1)
    end do
    !$omp end do

    !$omp do private(xi, yj, zk)
    do i = 1, ubound(Vp, 1)
      xi = Vp(i)%xi
      yj = Vp(i)%yj
      zk = Vp(i)%zk

      !$omp atomic
      Q(xi,yj,zk)   = Q(xi,yj,zk)   + V(xi,yj,zk) / dyPr(yj)
      !$omp atomic
      Q(xi,yj+1,zk) = Q(xi,yj+1,zk) - V(xi,yj,zk) / dyPr(yj+1)
    end do
    !$omp end do


    !$omp do private(xi, yj, zk)
    do i = 1, ubound(Wp, 1)
      xi = Wp(i)%xi
      yj = Wp(i)%yj
      zk = Wp(i)%zk

      !$omp atomic
      Q(xi,yj,zk)   = Q(xi,yj,zk)   + W(xi,yj,zk) / dzPr(zk)
      !$omp atomic
      Q(xi,yj,zk+1) = Q(xi,yj,zk+1) - W(xi,yj,zk) / dzPr(zk+1)
    end do
    !$omp end do
    !$omp end parallel

    call Bound_Q(Q)


  end subroutine IBMassSources


  subroutine AttenuateTop(U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V, W
    integer   :: i, j, k, bufn, mini, maxi, maxUi
    real(knd) :: ze, zs, zb, p
    real(knd), dimension(:), allocatable :: DF, avg

    if (Btype(We)==DIRICHLET.or.Btype(We)==TURBULENTINLET.or.Btype(We)==INLETFROMFILE) then
      mini = min(5, Unx)
    else
      mini = 1
    end if

    if (Btype(Ea)==DIRICHLET.or.Btype(We)==TURBULENTINLET.or.Btype(We)==OUTLETBUFF) then
      maxi = max(1, Prnx-5)
      maxUi = max(1, Unx-5)
    else
      maxi = Prnx
      maxUi = Unx
    end if

    bufn = max(5, Prnz/4)
    zs = zW(Prnz-bufn)
    ze = zW(Prnz)

    allocate(DF(  min(Unz, Vnz, Wnz) - bufn  :  max(Unz, Vnz, Wnz)))
    allocate(avg( min(Unz, Vnz, Wnz) - bufn  :  max(Unz, Vnz, Wnz)))



    do k = Unz - bufn, Unz
      avg(k) = 0
    end do

    !$omp parallel private(i, j, k, p, zb)
    
    !$omp do
    do k = Unz - bufn, Unz
      p = 0

      do j = 1, Uny
        do i = mini, maxUi
          p = p + U(i,j,k)
        end do
      end do
      avg(k) = p
    end do
    !$omp end do

    !$omp do
    do k = Unz - bufn, Unz
      avg(k) = avg(k) / ((maxUi-mini+1)*Uny)
    end do
    !$omp end do

    !$omp do
    do k = Unz - bufn, Unz
      zb=(zPr(k)-zs) / (ze-zs)
      DF(k) = DampF(zb)
    end do
    !$omp end do

    !$omp do
    do k = Unz - bufn, Unz
      do j = -1, Uny + 1
        do i = -1, Unx + 1
          U(i,j,k) = avg(k) + DF(k) * (U(i,j,k) - avg(k))
        end do
      end do
    end do
    !$omp end do


    !$omp do
    do k = Vnz - bufn, Vnz
      avg(k) = 0
    end do
    !$omp end do

    !$omp do
    do k = Vnz - bufn, Vnz
      p = 0

      do j = 1, Vny
        do i = mini, maxi
          p = p + V(i,j,k)
        end do
      end do
      avg(k) = p
    end do
    !$omp end do

    !$omp do
    do k = Vnz - bufn, Vnz
      avg(k) = avg(k) / ((maxi-mini+1)*Vny)
    end do
    !$omp end do

    !$omp do
    do k = Vnz - bufn, Vnz
      zb = (zPr(k)-zs) / (ze-zs)
      DF(k) = DampF(zb)
    end do
    !$omp end do

    !$omp do
    do k = Vnz - bufn, Vnz
      do j = -1, Vny + 1
        do i = -1, Vnx + 1
          V(i,j,k) = avg(k) + DF(k) * (V(i,j,k) - avg(k))
        end do
      end do
    end do
    !$omp end do


    !$omp do
    do k = Wnz - bufn, Wnz
      avg(k) = 0
    end do
    !$omp end do

    !$omp do
    do k = Wnz - bufn, Wnz
      p = 0

      do j = 1, Wny
        do i = mini, maxi
          p = p + W(i,j,k)
        end do
      end do
      avg(k) = p
    end do
    !$omp end do

    !$omp do
    do k = Wnz - bufn, Wnz
      avg(k) = avg(k) / ((maxi-mini+1)*Wny)
    end do
    !$omp end do

    !$omp do
    do k = Wnz - bufn, Wnz
      zb = (zW(k)-zs) / (ze-zs)
      DF(k) = DampF(zb)
    end do
    !$omp end do

    !$omp do
    do k = Wnz - bufn, Wnz
      do j = -1, Wny + 1
        do i = -1, Wnx + 1
          W(i,j,k) = avg(k) + DF(k) * (W(i,j,k) - avg(k))
        end do
      end do
    end do
    !$omp end do
    
    !$omp end parallel

  end subroutine AttenuateTop



  subroutine AttenuateOut(U, V, W, temperature)
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

  end subroutine AttenuateOut



  pure function DampF(x)
    real(knd) :: DampF
    real(knd), intent(in) :: x

    if (x<=0) then
      DampF = 1
    else if (x>=1) then
      DampF = 0
    else
      DampF = (1 - 0.04_knd*x**2) * &
              ( 1 - (1 - exp(10._knd*x**2)) / (1 - exp(10._knd)) )
    end if
  end function Dampf


  subroutine NullInterior(U, V, W)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V, W
    integer :: i

    !$omp parallel private(i)
    !$omp do
    do i = 1, nUnull
      U(Unull(1,i),Unull(2,i),Unull(3,i)) = 0
    end do
    !$omp end do nowait

    !$omp do
    do i = 1, nVnull
      V(Vnull(1,i),Vnull(2,i),Vnull(3,i)) = 0
    end do
    !$omp end do nowait

    !$omp do
    do i = 1, nWnull
      W(Wnull(1,i),Wnull(2,i),Wnull(3,i)) = 0
    end do
    !$omp end do
    !$omp end parallel

  end subroutine NullInterior


end module TimeSteps











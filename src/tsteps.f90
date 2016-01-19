module TimeSteps

  use Parameters
  use Dynamics
  use Boundaries, only: BoundU, Bound_Q
  use Pressure, only: PressureCorrection
  use Outputs, only: current_profiles
  use Scalars, only: ScalarRK3
  use Turbinlet, only: GetTurbulentInlet, GetBC_INLET_FROM_FILE
  use Sponge, only: enable_top_sponge, enable_out_sponge, SpongeTop, SpongeOut

  implicit none


  private
  public TMarchRK3


contains


  subroutine TMarchRK3(U, V, W, Pr, Temperature, Moisture, Scalar, dt, delta)
    use RK3
#ifdef PAR
    use custom_par
    use domains_bc_par
#endif
    real(knd), allocatable, intent(inout) :: U(:,:,:), V(:,:,:) ,W(:,:,:), Pr(:,:,:)
    real(knd), allocatable, intent(inout) :: Temperature(:,:,:), Moisture(:,:,:), Scalar(:,:,:,:)
    real(knd), intent(out) :: dt, delta

    integer :: RK_stage
    integer, save :: called = 0
    integer(int64), save :: trate
    integer(int64), save :: time1, time2

#ifdef  CUSTOM_TIMESTEP_PROCEDURE
    interface
      subroutine CustomTimeStepProcedure
      end subroutine
    end interface
#endif

    effective_time = time

    !uses previous dt (it should be fixed anyway, but...)
    call par_exchange_domain_bounds(U, V, W, Temperature, Moisture, Scalar, time, dt)

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


      if (enable_ibm_mass_sources) allocate(Q(0:Prnx+1, 0:Prny+1, 0:Prnz+1))

      if (debugparam>1) call system_clock(count_rate=trate)
    end if


    if ((Btype(We)==BC_TURBULENT_INLET) .or. (Btype(Ea)==BC_TURBULENT_INLET)) then
      call GetTurbulentInlet(dt)
    else if (Btype(We)==BC_INLET_FROM_FILE) then
      call GetBC_INLET_FROM_FILE(time)
    end if

    call TimeStepLength(U, V, W, dt)

    if (master) write (*,'(a,f12.6,a,es12.4)') " time: ", time," dt: ", dt
    
#ifdef  CUSTOM_TIMESTEP_PROCEDURE
    call CustomTimeStepProcedure
#endif


    do RK_stage = 1, RK_stages

      if (debugparam>1.and.master) call system_clock(count=time1)

      effective_time = time + 2 * sum(RK_alpha(1:RK_stage-1)) * dt      

      call SubgridStresses(U, V, W, Pr, Temperature)


      call Convection(U, V, W, &
                      U2, V2, W2, &
                      Ustar, Vstar, Wstar, &
                      Temperature, Moisture, &
                      RK_beta, RK_rho, RK_stage, dt)
)
      call ScalarRK3(U, V, W, &
                     Temperature, Moisture, Scalar, &
                     RK_stage, dt, &
                     current_profiles%tempfl, current_profiles%moistfl)

      call OtherTerms(U, V, W, &
                      U2, V2, W2, &
                      Pr, &
                      2*RK_alpha(RK_stage)*dt)

      if (enable_top_sponge)  then

          call SpongeTop(U2, V2, W2)
      end if

      if (enable_out_sponge) then

          call SpongeOut(U2, V2, W2, temperature)
      end if



      call BoundU(1, U2, Uin)

      call BoundU(2, V2, Vin)

      call BoundU(3, W2, Win)

      call IBMomentum(U2, V2, W2)

      if (enable_ibm_mass_sources) then
          call IBMassSources(Q, U2, V2, W2)
      end if


      if (poisson_solver>0) then
        call PressureCorrection(U2, V2, W2, Pr, Q, 2*RK_alpha(RK_stage)*dt)
      end if



      if (RK_stage==1) delta = 0

      if ( debuglevel>0 .or. steady==1 ) then

        if (Unx*Uny*Unz > 0) &
          delta = delta + sum(abs(U(1:Unx,1:Uny,1:Unz) - U2(1:Unx,1:Uny,1:Unz))) / (Unx*Uny*Unz)

        if (Vnx*Vny*Vnz > 0) &
          delta = delta + sum(abs(V(1:Vnx,1:Vny,1:Vnz) - V2(1:Vnx,1:Vny,1:Vnz))) / (Vnx*Vny*Vnz)

        if (Wnx*Wny*Wnz > 0) &
          delta = delta + sum(abs(W(1:Wnx,1:Wny,1:Wnz) - W2(1:Wnx,1:Wny,1:Wnz))) / (Wnx*Wny*Wnz)

        if (RK_stage==RK_stages) then
#ifdef PAR
          delta = par_co_sum(delta)
#endif
          if (master) write(*,*) "delta",delta
        end if
      end if

      call exchange_alloc(U, U2)
      call exchange_alloc(V, V2)
      call exchange_alloc(W, W2)



      if (enable_out_sponge) then

        call SpongeOut(U, V, W, temperature)

      end if


      call NullInterior(U, V, W)


      call BoundU(1, U, Uin)

      call BoundU(2, V, Vin)

      call BoundU(3, W, Win)

      call IBMomentum(U, V, W)


      if (debugparam>1.and.master) then
        call system_clock(count=time2)
        write (*,*) "ET of part 1", (time2-time1) / real(trate, dbl)
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

        Re_gt_0: if (molecular_viscosity > 0) then

          !Diffusion using Crank Nicolson
          !first approximation using forward Euler
          !iteration SOR or Gauss-Seidel

          call ImplicitDiffusion_ForwEul(U, V, W, U2, V2, W2, U3, V3, W3, coef)


          call BoundU(1, U3, Uin)

          call BoundU(2, V3, Vin)

          call BoundU(3, W3, Win)

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


        call BoundU(1, U2, Uin)

        call BoundU(2, V2, Vin)

        call BoundU(3, W2, Win)

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

    call set(Q, 0)

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











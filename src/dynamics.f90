module Dynamics
  use Parameters
  use ArrayUtilities
  use Tiling, only: tilenx, tileny, tilenz
  use Boundaries, only: BoundU
  use ScalarBoundaries, only: BoundTemperature, BoundViscosity

  implicit none

  !module variables to allow their deallocation before the program end
  real(knd), dimension(:,:,:), allocatable:: U3,V3,W3

  real(knd), dimension(:,:,:), allocatable :: Q
  real(knd), dimension(:,:,:), allocatable :: U2,Ustar
  real(knd), dimension(:,:,:), allocatable :: V2,Vstar
  real(knd), dimension(:,:,:), allocatable :: W2,Wstar

  real(knd), dimension(:,:,:), allocatable :: Apu, ApV, ApW, wrk

  real(knd), dimension(:), allocatable :: Uwm, Vwm, Wwm

contains


  subroutine Dynamics_Deallocate
    !Deallocates the working arrays
    use Wallmodels
    if (allocated(U3)) deallocate(U3)
    if (allocated(V3)) deallocate(V3)
    if (allocated(W3)) deallocate(W3)
    if (allocated(Q)) deallocate(Q)
    if (allocated(U2)) deallocate(U2)
    if (allocated(V2)) deallocate(V2)
    if (allocated(W2)) deallocate(W2)
    if (allocated(Ustar)) deallocate(Ustar)
    if (allocated(Vstar)) deallocate(Vstar)
    if (allocated(Wstar)) deallocate(Wstar)
    if (allocated(ApU)) deallocate(ApU)
    if (allocated(ApV)) deallocate(ApV)
    if (allocated(ApW)) deallocate(ApW)
    if (allocated(wrk)) deallocate(wrk)

    if (allocated(Uwm)) deallocate(Uwm)
    if (allocated(Vwm)) deallocate(Vwm)
    if (allocated(Wwm)) deallocate(Wwm)

    !imported from Wallmodels
    if (allocated(Uflx_mask)) deallocate(Uflx_mask)
    if (allocated(Ufly_mask)) deallocate(Ufly_mask)
    if (allocated(Uflz_mask)) deallocate(Uflz_mask)
    if (allocated(Vflx_mask)) deallocate(Vflx_mask)
    if (allocated(Vfly_mask)) deallocate(Vfly_mask)
    if (allocated(Vflz_mask)) deallocate(Vflz_mask)
    if (allocated(Wflx_mask)) deallocate(Wflx_mask)
    if (allocated(Wfly_mask)) deallocate(Wfly_mask)
    if (allocated(Wflz_mask)) deallocate(Wflz_mask)
  end subroutine


  subroutine PressureGrad(Pr, U, V, W, coef)
    real(knd), intent(inout), contiguous :: Pr(1:,1:,1:)
    real(knd), intent(inout), contiguous, dimension(-2:,-2:,-2:) :: U,V,W
    real(knd), intent(in)    :: coef
    real(knd) :: A, Ax, Ay, Az
    integer :: i,j,k

    A = -coef
    Ax = - coef / dxmin
    Ay = - coef / dymin
    Az = - coef / dzmin

    !$omp parallel
    if (enable_pr_gradient_x_profile) then
       !$omp do
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            U(i,j,k) = U(i,j,k) + A * pr_gradient_profile_x(k)
          end do
        end do
      end do
      !$omp end do
    else if (enable_pr_gradient_x_uniform) then
      !$omp do
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            U(i,j,k) = U(i,j,k) + A * pr_gradient_x
          end do
        end do
      end do
      !$omp end do
    end if

    if (enable_pr_gradient_y_profile) then
      !$omp do
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            V(i,j,k) = V(i,j,k) + A * pr_gradient_profile_y(k)
          end do
        end do
      end do
      !$omp end do
    else if (enable_pr_gradient_y_uniform) then
      !$omp do
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            V(i,j,k) = V(i,j,k) + A * pr_gradient_y
          end do
        end do
      end do
      !$omp end do
    end if

    !$omp do
    do k = 1, Unz
      do j = 1, Uny
        do i = 1, Unx
          U(i,j,k) = U(i,j,k) + Ax * (Pr(i+1,j,k)-Pr(i,j,k))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k = 1, Vnz
      do j = 1, Vny
        do i = 1, Vnx
          V(i,j,k) = V(i,j,k) + Ay * (Pr(i,j+1,k)-Pr(i,j,k))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k = 1, Wnz
      do j = 1, Wny
        do i = 1, Wnx
          W(i,j,k) = W(i,j,k) + Az * (Pr(i,j,k+1)-Pr(i,j,k))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp end parallel
  end subroutine PressureGrad


  
  
  subroutine StressBoundaryFlux(U2, V2, dt)
    use ArrayUtilities,only: add
    use Outputs,only: current_profiles
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2, V2
    real(knd), intent(in) :: dt
    real(knd) :: flux
    integer :: first,last
    
    if (Btype(To)==BC_AUTOMATIC_FLUX) then
      first = min(Prnz*5/6,Prnz-5)
      last = Prnz-5
      
      flux = sum(current_profiles%uw(first:last)) + sum(current_profiles%uwsgs(first:last))
      flux = flux / (last-first+1)
      call add(U2(:,:,Unz), -dt*flux/dzPr(Unz))
      
      flux = sum(current_profiles%vw(first:last)) + sum(current_profiles%vwsgs(first:last))
      flux = flux / (last-first+1)
      call add(V2(:,:,Vnz), -dt*flux/dzPr(Vnz))
    end if
  end subroutine





  subroutine Convection(U, V, W, U2, V2, W2, &
                        Ustar, Vstar, Wstar, &
                        Temperature, Moisture, &
                        beta, rho, RK_stage, dt)
    use MomentumAdvection
    use VolumeSources, only: ResistanceForce
    use VTKArray
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(out)   :: U2, V2, W2
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: Ustar, Vstar ,Wstar
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(in)    :: Temperature
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(in)    :: Moisture
    real(knd), dimension(1:3), intent(in) :: beta,rho
    integer,   intent(in) :: RK_stage
    real(knd), intent(in) :: dt
    integer :: i,j,k

    if (RK_stage>1) then
      !$omp parallel private(i,j,k)
      !$omp do
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            U2(i,j,k) = Ustar(i,j,k) * rho(RK_stage) * dt
          end do
        end do
      end do
      !$omp end do nowait
      !$omp do
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            V2(i,j,k) = Vstar(i,j,k) * rho(RK_stage) * dt
          end do
        end do
      end do
      !$omp end do nowait
      !$omp do
      do k = 1, Wnz
        do j = 1, Wny
          do i = 1, Wnx
            W2(i,j,k) = Wstar(i,j,k) * rho(RK_stage) * dt
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
    else
      !$omp parallel private(i,j,k)
      !$omp do
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            Ustar(i,j,k) = 0
          end do
        end do
      end do
      !$omp end do nowait
      !$omp do
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            Vstar(i,j,k) = 0
          end do
        end do
      end do
      !$omp end do nowait
      !$omp do
      do k = 1, Wnz
        do j = 1, Wny
          do i = 1, Wnx
            Wstar(i,j,k) = 0
          end do
        end do
      end do
      !$omp end do nowait

      !$omp do
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            U2(i,j,k) = 0
          end do
        end do
      end do
      !$omp end do nowait
      !$omp do
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            V2(i,j,k) = 0
          end do
        end do
      end do
      !$omp end do nowait
      !$omp do
      do k = 1, Wnz
        do j = 1, Wny
          do i = 1, Wnx
            W2(i,j,k) = 0
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
    end if

    if (advection_method>0) then

      if (advection_method==2) then
        call CDU(Ustar,U,V,W)
        call CDV(Vstar,U,V,W)
        call CDW(Wstar,U,V,W)
      else if (advection_method==4) then
        call CD4divU(Ustar,U,V,W)
        call CD4divV(Vstar,U,V,W)
        call CD4divW(Wstar,U,V,W)
      else if (advection_method==5) then
        call CDUdiv(Ustar,U,V,W)
        call CDVdiv(Vstar,U,V,W)
        call CDWdiv(Wstar,U,V,W)
      else if (advection_method==6) then
        call set(Ustar, 0)
        call set(Vstar, 0)
        call set(Wstar, 0)
        call CDUadv(Ustar,U,V,W)
        call CDVadv(Vstar,U,V,W)
        call CDWadv(Wstar,U,V,W)
      end if

    else

      !$omp parallel private(i,j,k)
      !$omp do
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            Ustar(i,j,k) = 0
          end do
        end do
      end do
      !$omp end do nowait
      !$omp do
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            Vstar(i,j,k) = 0
          end do
        end do
      end do
      !$omp end do nowait
      !$omp do
      do k = 1, Wnz
        do j = 1, Wny
          do i = 1, Wnx
            Wstar(i,j,k) = 0
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel

    end if

    call FilterUstar    

    if (abs(Coriolis_parameter)>tiny(1._knd)) call CoriolisForce(Ustar, Vstar, U, V)

    if (enable_buoyancy) call BuoyancyForce(Wstar, Temperature, Moisture)

    call ResistanceForce(Ustar, Vstar, Wstar, U, V, W)

    call StressBoundaryFlux(Ustar, Vstar, dt)

    if (explicit_diffusion) call MomentumDiffusion_nobranch(Ustar, Vstar, Wstar, U, V, W)


    !$omp parallel private(i,j,k)
    !$omp do
    do k = 1, Unz
      do j = 1, Uny
        do i = 1, Unx
          U2(i,j,k) = U2(i,j,k)  +  Ustar(i,j,k) * beta(RK_stage) * dt
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k = 1, Vnz
      do j = 1, Vny
        do i = 1, Vnx
          V2(i,j,k) = V2(i,j,k)  +  Vstar(i,j,k) * beta(RK_stage) * dt
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k = 1, Wnz
      do j = 1, Wny
        do i = 1, Wnx
          W2(i,j,k) = W2(i,j,k)  +  Wstar(i,j,k) * beta(RK_stage) * dt
        end do
      end do
    end do
    !$omp end do
    !$omp end parallel

  contains

    subroutine FilterUstar

      use Filters, only: filtertype, Filter

      if (filtertype/=0) then
        call BoundU(1,Ustar,Uin,2)
        call BoundU(2,Vstar,Vin,2)
        call BoundU(3,Wstar,Win,2)


        call Filter(Ustar,Utype)

        call Filter(Vstar,Vtype)

        call Filter(Wstar,Wtype)

        call BoundU(1,Ustar,Uin,2)
        call BoundU(2,Vstar,Vin,2)
        call BoundU(3,Wstar,Win,2)
      end if

    end subroutine
  end subroutine Convection











  subroutine TimeStepLength(U, V, W, dt)
    use ieee_arithmetic
#ifdef PAR
    use custom_par
#endif
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)  :: U,V,W
    real(tim), intent(out) :: dt
    integer :: i,j,k
    real(knd) :: m, p
    logical :: nan
    
    nan = .false.
    m = 0
    !$omp parallel do private(i,j,k,p) reduction(max:m) reduction(.or.:nan)
    do k = 1, Prnz
      do j = 1, Prny
        do i = 1, Prnx
          !For scalar advection the sum proved to be necessary when the flow is not aligned to grid.
          p =     max( abs(U(i,j,k)), abs(U(i-1,j,k)) ) / dxmin
          p = p + max( abs(V(i,j,k)), abs(V(i,j-1,k)) ) / dymin
          p = p + max( abs(W(i,j,k)), abs(W(i,j,k-1)) ) / dzmin
          
          m = max(m,p)
          if (ieee_is_nan(p)) nan = .true.
        end do
      end do
    end do
    !$omp end parallel do

    if (nan) then
      dt = tiny(dt)
    else if (m>0) then
      dt = min(CFL/m, min(dxmin,dymin,dzmin)/Uref)
    else
      dt = min(dxmin,dymin,dzmin) / Uref
    end if

    if (steady/=1 .and. dt+time>end_time)  dt = end_time-time
    
#ifdef PAR
    dt = par_co_min(dt)
#endif

  endsubroutine TimeStepLength







  subroutine BuoyancyForce(W, Temperature, Moisture)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: W
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(in) :: Temperature, Moisture
    real(knd) :: A, A2
    integer :: i, j, k


    if (enable_moisture) then
      if (enable_liquid) then
        call error_stop("Liquid water not implemented.")
      else
        A = grav_acc / temperature_ref
        A2 = A / 2._KND

        call apply_moist(1)
        call apply_moist(2)
      end if
    else
      A = grav_acc / temperature_ref
      A2 = A / 2._KND
      !$omp parallel do private(i,j,k)
      do k = 1, Wnz
        do j = 1, Wny
          do i = 1, Wnx
            W(i,j,k) = W(i,j,k) + A2 * (Temperature(i,j,k+1)+Temperature(i,j,k)) - A * temperature_ref
          end do
        end do
      end do
      !$omp end parallel do
    end if

    contains

      subroutine apply_moist(start)
        integer, intent(in) :: start
        integer :: i, j, k
        real(knd) :: temperature_virt

        !$omp parallel do private(i,j,k,temperature_virt)
        do k = start, Wnz+1,2
          do j = 1, Wny
            do i = 1, Wnx
              temperature_virt = theta_v(i,j,k)
              W(i,j,k)   = W(i,j,k)   + A2 * temperature_virt - A * temperature_ref
              W(i,j,k-1) = W(i,j,k-1) + A2 * temperature_virt
            end do
          end do
        end do
        !$omp end parallel do
      end subroutine

      pure real(knd) function theta_v(i, j, k)
        integer, intent(in) :: i, j, k

        theta_v = Temperature(i, j, k) * (1._knd + 0.61_knd * Moisture(i, j, k))
      end function
  end subroutine BuoyancyForce





  subroutine CoriolisForce(U2, V2, U, V)
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U, V
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2, V2
    integer :: i,j,k

    if (abs(Coriolis_parameter)>0) then
    !$omp parallel private(i,j,k)
    !$omp do
      do k = 1, Unz
        do j = 1, Uny
          do i = 1, Unx
            U2(i,j,k) = U2(i,j,k) + &
                  Coriolis_parameter * (V(i,j-1,k)+V(i+1,j-1,k)+V(i,j,k)+V(i+1,j,k))/4._knd
          end do
        end do
      end do
      !$omp end do nowait

      !$omp do
      do k = 1, Vnz
        do j = 1, Vny
          do i = 1, Vnx
            V2(i,j,k) = V2(i,j,k) - &
                  Coriolis_parameter * (U(i-1,j,k)+U(i-1,j+1,k)+U(i,j,k)+U(i,j+1,k))/4._knd
          end do
        end do
      end do
      !$omp end do
      !$omp end parallel
    end if
  end subroutine CoriolisForce





  subroutine SubgridStresses(U,V,W,Pr,Temperature)
    use Subgrid, only: SubgridModel, sgstype, StabSubgridModel
    use ImmersedBoundary, only: ScalFlIBPoints, TIBPoint_Viscosity
    use Wallmodels
    use Scalars, only: ComputeTDiff

    real(knd), contiguous, intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: Pr(1:,1:,1:)
    real(knd), contiguous, intent(in) :: Temperature(-1:,-1:,-1:)
    integer :: i


    if (wallmodeltype>0) then
                      !resulting Viscosity intentionally overwritten
                      call ComputeViscsWM(U,V,W,Pr,Temperature)
    end if
    

    call SubgridModel(U,V,W)


    if (wallmodeltype>0) then
                    call ComputeUVWFluxesWM(U,V,W,Pr,Temperature)
    end if

    call BoundViscosity(Viscosity)

    do i = 1, size(ScalFlIBPoints)
      Viscosity(ScalFlIBPoints(i)%xi,ScalFlIBPoints(i)%yj,ScalFlIBPoints(i)%zk) =  &
                                              TIBPoint_Viscosity(ScalFlIBPoints(i),Viscosity)
    end do

    if (sgstype/=StabSubgridModel.and.enable_buoyancy)  call ComputeTDiff(U,V,W)

    if (size(TDiff)>0) call BoundViscosity(TDiff)

  end subroutine SubgridStresses









  subroutine MomentumDiffusion(U2,V2,W2,U,V,W)
    use Parameters, nu => Viscosity
    use Wallmodels
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U,V,W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2,V2,W2
    real(knd) :: recdxmin2, recdymin2, recdzmin2
    integer :: i,j,k,bi,bj,bk
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 6
       
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    recdxmin2 = 1._knd / dxmin**2
    recdymin2 = 1._knd / dymin**2
    recdzmin2 = 1._knd / dzmin**2

    !$omp parallel private(i,j,k,bi,bj,bk)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
           if (Uflx_mask(i+1,j,k)) &
             U2(i,j,k) = U2(i,j,k) + &
              nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) *recdxmin2
           if (Uflx_mask(i,j,k)) &
             U2(i,j,k) = U2(i,j,k) - &
               nu(i,j,k) * (U(i,j,k)-U(i-1,j,k)) * recdxmin2 
           if (Ufly_mask(i,j+1,k)) &
             U2(i,j,k) = U2(i,j,k) + &
               0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (U(i,j+1,k)-U(i,j,k)) * recdymin2
           if (Ufly_mask(i,j,k)) &
             U2(i,j,k) = U2(i,j,k) - &
               0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * (U(i,j,k)-U(i,j-1,k)) * recdymin2
           if (Uflz_mask(i,j,k+1)) &
             U2(i,j,k) = U2(i,j,k) + &
               0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (U(i,j,k+1)-U(i,j,k)) * recdzmin2
           if (Uflz_mask(i,j,k)) &
             U2(i,j,k) = U2(i,j,k) - &
               0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (U(i,j,k)-U(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 1
#define wrk U2
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           if (Vflx_mask(i+1,j,k)) &
             V2(i,j,k) = V2(i,j,k) + &
               0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (V(i+1,j,k)-V(i,j,k)) * recdxmin2
           if (Vflx_mask(i,j,k)) &
             V2(i,j,k) = V2(i,j,k) - &
               0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * (V(i,j,k)-V(i-1,j,k)) * recdxmin2
           if (Vfly_mask(i,j+1,k)) &
             V2(i,j,k) = V2(i,j,k) + &
               nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k)) * recdymin2
           if (Vfly_mask(i,j,k)) &
             V2(i,j,k) = V2(i,j,k) - &
               nu(i,j,k) * (V(i,j,k)-V(i,j-1,k)) * recdymin2
           if (Vflz_mask(i,j,k+1)) &
             V2(i,j,k) = V2(i,j,k) + &
               0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * (V(i,j,k+1)-V(i,j,k)) * recdzmin2
           if (Vflz_mask(i,j,k)) &
             V2(i,j,k) = V2(i,j,k) - &
               0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (V(i,j,k)-V(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 2
#define wrk V2
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           if (Wflx_mask(i+1,j,k)) &
             W2(i,j,k) = W2(i,j,k) + &
               0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (W(i+1,j,k)-W(i,j,k)) * recdxmin2
           if (Wflx_mask(i,j,k)) &
             W2(i,j,k) = W2(i,j,k) - &
               0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * (W(i,j,k)-W(i-1,j,k)) * recdxmin2
           if (Wfly_mask(i,j+1,k)) &
             W2(i,j,k) = W2(i,j,k) + &
               0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * (W(i,j+1,k)-W(i,j,k)) * recdymin2
           if (Wfly_mask(i,j,k)) &
             W2(i,j,k) = W2(i,j,k) - &
               0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * (W(i,j,k)-W(i,j-1,k)) * recdymin2
           if (Wflz_mask(i,j,k+1)) &
             W2(i,j,k) = W2(i,j,k) + &
               nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k)) * recdzmin2
           if (Wflz_mask(i,j,k)) &
             W2(i,j,k) = W2(i,j,k) - &
               nu(i,j,k) * (W(i,j,k)-W(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 3
#define wrk W2
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
    !$omp end parallel

  end subroutine MomentumDiffusion














  subroutine MomentumDiffusion_nobranch(U2,V2,W2,U,V,W)
    use Parameters, nu => Viscosity
    use Wallmodels
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in)    :: U,V,W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2,V2,W2
    real(knd) :: recdxmin2, recdymin2, recdzmin2
    integer :: i,j,k,bi,bj,bk
    integer :: tnx, tny, tnz
    
    integer, parameter :: narr = 3
       
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)
       
    recdxmin2 = 1._knd / dxmin**2
    recdymin2 = 1._knd / dymin**2
    recdzmin2 = 1._knd / dzmin**2

    !$omp parallel private(i,j,k,bi,bj,bk)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
             U2(i,j,k) = U2(i,j,k) + &
              nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) * recdxmin2
             U2(i,j,k) = U2(i,j,k) - &
               nu(i,j,k) * (U(i,j,k)-U(i-1,j,k)) * recdxmin2 
             U2(i,j,k) = U2(i,j,k) + &
               0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (U(i,j+1,k)-U(i,j,k)) * recdymin2
             U2(i,j,k) = U2(i,j,k) - &
               0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * (U(i,j,k)-U(i,j-1,k)) * recdymin2
             U2(i,j,k) = U2(i,j,k) + &
               0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (U(i,j,k+1)-U(i,j,k)) * recdzmin2
             U2(i,j,k) = U2(i,j,k) - &
               0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (U(i,j,k)-U(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 1
#define wrk U2
#include "wmfluxes-nobranch-U-inc.f90"
#undef wrk
#undef comp


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
             V2(i,j,k) = V2(i,j,k) + &
               0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (V(i+1,j,k)-V(i,j,k)) * recdxmin2
             V2(i,j,k) = V2(i,j,k) - &
               0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * (V(i,j,k)-V(i-1,j,k)) * recdxmin2
             V2(i,j,k) = V2(i,j,k) + &
               nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k)) * recdymin2
             V2(i,j,k) = V2(i,j,k) - &
               nu(i,j,k) * (V(i,j,k)-V(i,j-1,k)) * recdymin2
             V2(i,j,k) = V2(i,j,k) + &
               0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * (V(i,j,k+1)-V(i,j,k)) * recdzmin2
             V2(i,j,k) = V2(i,j,k) - &
               0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (V(i,j,k)-V(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 2
#define wrk V2
#include "wmfluxes-nobranch-V-inc.f90"
#undef wrk
#undef comp


    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
             W2(i,j,k) = W2(i,j,k) + &
               0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (W(i+1,j,k)-W(i,j,k)) * recdxmin2
             W2(i,j,k) = W2(i,j,k) - &
               0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * (W(i,j,k)-W(i-1,j,k)) * recdxmin2
             W2(i,j,k) = W2(i,j,k) + &
               0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * (W(i,j+1,k)-W(i,j,k)) * recdymin2
             W2(i,j,k) = W2(i,j,k) - &
               0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * (W(i,j,k)-W(i,j-1,k)) * recdymin2
             W2(i,j,k) = W2(i,j,k) + &
               nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k)) * recdzmin2
             W2(i,j,k) = W2(i,j,k) - &
               nu(i,j,k) * (W(i,j,k)-W(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 3
#define wrk W2
#include "wmfluxes-nobranch-W-inc.f90"
#undef wrk
#undef comp
    !$omp end parallel

  end subroutine MomentumDiffusion_nobranch



















! Remaining parts of implicit diffusion




  subroutine ImplicitDiffusion_ForwEul(U, V, W, U2, V2, W2, U3, V3, W3, coef)
    use Parameters, nu => Viscosity
    use Wallmodels
    real(knd), intent(in),  dimension(-2:,-2:,-2:), contiguous :: U,V,W
    real(knd), intent(in),  dimension(-2:,-2:,-2:), contiguous :: U2,V2,W2
    real(knd), intent(out), dimension(-2:,-2:,-2:), contiguous :: U3,V3,W3
    real(knd), intent(in) :: coef

    real(knd) :: Ap, recdxmin2, recdymin2, recdzmin2
    integer   :: i,j,k


    Ap = coef

    recdxmin2 = 1._knd / dxmin**2
    recdymin2 = 1._knd / dymin**2
    recdzmin2 = 1._knd / dzmin**2

    !$omp parallel private(i,j,k)

    !$omp do
    do k = 1, Unz
      do j = 1, Uny
        do i = 1, Unx
            U3(i,j,k) = 0
            if (Uflx_mask(i+1,j,k)) &
              U3(i,j,k) = U3(i,j,k) + &
                nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) *recdxmin2
            if (Uflx_mask(i,j,k)) &
              U3(i,j,k) = U3(i,j,k) - &
                nu(i,j,k) * (U(i,j,k)-U(i-1,j,k)) * recdxmin2 
            if (Ufly_mask(i,j+1,k)) &
              U3(i,j,k) = U3(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (U(i,j+1,k)-U(i,j,k)) * recdymin2
            if (Ufly_mask(i,j,k)) &
              U3(i,j,k) = U3(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * (U(i,j,k)-U(i,j-1,k)) * recdymin2
            if (Uflz_mask(i,j,k+1)) &
              U3(i,j,k) = U3(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (U(i,j,k+1)-U(i,j,k)) * recdzmin2
            if (Uflz_mask(i,j,k)) &
              U3(i,j,k) = U3(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (U(i,j,k)-U(i,j,k-1)) * recdzmin2
        end do
      end do
    end do
    !$omp end do
#define comp 1
#define wrk U3
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
    !$omp do
    do k = 1, Unz
      do j = 1, Uny
        do i = 1, Unx
          U3(i,j,k) = U3(i,j,k) * Ap
          U3(i,j,k) = U3(i,j,k) + U(i,j,k) + U2(i,j,k)
        end do
      end do
    end do
    !$omp end do nowait


    !$omp do
    do k = 1, Vnz
      do j = 1, Vny
        do i = 1, Vnx
            V3(i,j,k) = 0
            if (Vflx_mask(i+1,j,k)) &
              V3(i,j,k) = V3(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (V(i+1,j,k)-V(i,j,k)) * recdxmin2
            if (Vflx_mask(i,j,k)) &
              V3(i,j,k) = V3(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * (V(i,j,k)-V(i-1,j,k)) * recdxmin2
            if (Vfly_mask(i,j+1,k)) &
              V3(i,j,k) = V3(i,j,k) + &
                nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k)) * recdymin2
            if (Vfly_mask(i,j,k)) &
              V3(i,j,k) = V3(i,j,k) - &
                nu(i,j,k) * (V(i,j,k)-V(i,j-1,k)) * recdymin2
            if (Vflz_mask(i,j,k+1)) &
              V3(i,j,k) = V3(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * (V(i,j,k+1)-V(i,j,k)) * recdzmin2
            if (Vflz_mask(i,j,k)) &
              V3(i,j,k) = V3(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (V(i,j,k)-V(i,j,k-1)) * recdzmin2
        end do
      end do
    end do
    !$omp end do
#define comp 2
#define wrk V3
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
    !$omp do
    do k = 1, Vnz
      do j = 1, Vny
        do i = 1, Vnx
          V3(i,j,k) = V3(i,j,k) * Ap
          V3(i,j,k) = V3(i,j,k) + V(i,j,k) + V2(i,j,k)
        end do
      end do
    end do
    !$omp end do nowait


    !$omp do
    do k = 1, Wnz
      do j = 1, Wny
        do i = 1, Wnx
            W3(i,j,k) = 0
            if (Wflx_mask(i+1,j,k)) &
              W3(i,j,k) = W3(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (W(i+1,j,k)-W(i,j,k)) * recdxmin2
            if (Wflx_mask(i,j,k)) &
              W3(i,j,k) = W3(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * (W(i,j,k)-W(i-1,j,k)) * recdxmin2
            if (Wfly_mask(i,j+1,k)) &
              W3(i,j,k) = W3(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * (W(i,j+1,k)-W(i,j,k)) * recdymin2
            if (Wfly_mask(i,j,k)) &
              W3(i,j,k) = W3(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * (W(i,j,k)-W(i,j-1,k)) * recdymin2
            if (Wflz_mask(i,j,k+1)) &
              W3(i,j,k) = W3(i,j,k) + &
                nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k)) * recdzmin2
            if (Wflz_mask(i,j,k)) &
              W3(i,j,k) = W3(i,j,k) - &
                nu(i,j,k) * (W(i,j,k)-W(i,j,k-1)) * recdzmin2
        end do
      end do
    end do
    !$omp end do
#define comp 3
#define wrk W3
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
    !$omp do
    do k = 1, Wnz
      do j = 1, Wny
        do i = 1, Wnx
          W3(i,j,k) = W3(i,j,k) * Ap
          W3(i,j,k) = W3(i,j,k) + W(i,j,k) + W2(i,j,k)
        end do
      end do
    end do
    !$omp end do

    !$omp end parallel

  end subroutine ImplicitDiffusion_ForwEul






  subroutine ImplicitDiffusion_Iterations(U, V, W, U2, V2, W2, U3, V3, W3, coef)
    use Parameters, nu => Viscosity
    use Wallmodels
    !$ use omp_lib
#ifdef PAR
    use custom_par, only: par_co_max
#endif
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V, W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2, V2, W2, U3, V3, W3
    real(knd), intent(in) :: coef
    real(knd) recdxmin2,recdymin2,recdzmin2                                                               !reciprocal values of dx**2
    real(knd) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
    integer :: i,j,k,bi,bj,bk,l
    integer, save :: called = 0
    integer :: tnx, tny, tnz, tnx2, tny2, tnz2
    
    integer, parameter :: narr = 3, narr2 = 5
    
    tnx = tilenx(narr)
    tny = tileny(narr)
    tnz = tilenz(narr)

    tnx2 = tilenx(narr2)
    tny2 = tileny(narr2)
    tnz2 = tilenz(narr2)


    if (called==0) then
      allocate(Apu(1:Unx,1:Uny,1:Unz))
      allocate(ApV(1:Vnx,1:Vny,1:Vnz))
      allocate(ApW(1:Wnx,1:Wny,1:Wnz))
      allocate(wrk(0:max(Unx+1,Vnx+1,Wnx+1), &
                  0:max(Uny+1,Vny+1,Wny+1), &
                  0:max(Unz+1,Uny+1,Wnz+1)))
      called = 1
    end if


    Ap = coef / 2
    S = 0
    l = 0

    recdxmin2 = 1._knd / dxmin**2
    recdymin2 = 1._knd / dymin**2
    recdzmin2 = 1._knd / dzmin**2

    !$omp parallel private(i,j,k,bi,bj,bk)

    !The explicit part, which doesn't have to be changed inside the loop
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
            wrk(i,j,k) = 0
            if (Uflx_mask(i+1,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) + &
                nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) *recdxmin2
            if (Uflx_mask(i,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) - &
                nu(i,j,k) * (U(i,j,k)-U(i-1,j,k)) * recdxmin2 
            if (Ufly_mask(i,j+1,k)) &
              wrk(i,j,k) = wrk(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (U(i,j+1,k)-U(i,j,k)) * recdymin2
            if (Ufly_mask(i,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * (U(i,j,k)-U(i,j-1,k)) * recdymin2
            if (Uflz_mask(i,j,k+1)) &
              wrk(i,j,k) = wrk(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (U(i,j,k+1)-U(i,j,k)) * recdzmin2
            if (Uflz_mask(i,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (U(i,j,k)-U(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 1
#include "wmfluxes-inc.f90"
#undef comp
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
              U2(i,j,k) = U2(i,j,k) + Ap * wrk(i,j,k)
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do nowait
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
            wrk(i,j,k) = 0
            if (Vflx_mask(i+1,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (V(i+1,j,k)-V(i,j,k)) * recdxmin2
            if (Vflx_mask(i,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * (V(i,j,k)-V(i-1,j,k)) * recdxmin2
            if (Vfly_mask(i,j+1,k)) &
              wrk(i,j,k) = wrk(i,j,k) + &
                nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k)) * recdymin2
            if (Vfly_mask(i,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) - &
                nu(i,j,k) * (V(i,j,k)-V(i,j-1,k)) * recdymin2
            if (Vflz_mask(i,j,k+1)) &
              wrk(i,j,k) = wrk(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * (V(i,j,k+1)-V(i,j,k)) * recdzmin2
            if (Vflz_mask(i,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (V(i,j,k)-V(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 2
#include "wmfluxes-inc.f90"
#undef comp
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
           V2(i,j,k) = V2(i,j,k) + Ap * wrk(i,j,k)
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
            wrk(i,j,k) = 0
            if (Wflx_mask(i+1,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (W(i+1,j,k)-W(i,j,k)) * recdxmin2
            if (Wflx_mask(i,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * (W(i,j,k)-W(i-1,j,k)) * recdxmin2
            if (Wfly_mask(i,j+1,k)) &
              wrk(i,j,k) = wrk(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * (W(i,j+1,k)-W(i,j,k)) * recdymin2
            if (Wfly_mask(i,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * (W(i,j,k)-W(i,j-1,k)) * recdymin2
            if (Wflz_mask(i,j,k+1)) &
              wrk(i,j,k) = wrk(i,j,k) + &
                nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k)) * recdzmin2
            if (Wflz_mask(i,j,k)) &
              wrk(i,j,k) = wrk(i,j,k) - &
                nu(i,j,k) * (W(i,j,k)-W(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
#define comp 3
#include "wmfluxes-inc.f90"
#undef comp
    !$omp do schedule(runtime)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
           W2(i,j,k) = W2(i,j,k) + Ap * wrk(i,j,k)
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do

    !Auxiliary coefficients to better efficiency in loops
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Unz, tnz
     do bj = 1, Uny, tny
      do bi = 1, Unx, tnx
       do k = bk, min(bk+tnz-1,Unz)
        do j = bj, min(bj+tny-1,Uny)
         do i = bi, min(bi+tnx-1,Unx)
            ApU(i,j,k) = 0
            if (Uflx_mask(i+1,j,k)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                nu(i+1,j,k) * recdxmin2
            if (Uflx_mask(i,j,k)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                nu(i,j,k) * recdxmin2
            if (Ufly_mask(i,j+1,k)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * recdymin2
            if (Ufly_mask(i,j,k)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * recdymin2
            if (Uflz_mask(i,j,k+1)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * recdzmin2
            if (Uflz_mask(i,j,k)) &
              ApU(i,j,k) = ApU(i,j,k) + &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp end parallel

!        ApU = 1._knd/(1._knd+Ap*ApU)
    call multiply_and_add_scalar(ApU, Ap, 1._knd)
    call reciprocal(ApU, 1._knd)

    !$omp parallel private(i,j,k,bi,bj,bk)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Vnz, tnz
     do bj = 1, Vny, tny
      do bi = 1, Vnx, tnx
       do k = bk, min(bk+tnz-1,Vnz)
        do j = bj, min(bj+tny-1,Vny)
         do i = bi, min(bi+tnx-1,Vnx)
            ApV(i,j,k) =0
            if (Vflx_mask(i+1,j,k)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * recdxmin2
            if (Vflx_mask(i,j,k)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * recdxmin2
            if (Vfly_mask(i,j+1,k)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                nu(i,j+1,k) * recdymin2
            if (Vfly_mask(i,j,k)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                nu(i,j,k) * recdymin2
            if (Vflz_mask(i,j,k+1)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * recdzmin2
            if (Vflz_mask(i,j,k)) &
              ApV(i,j,k) = ApV(i,j,k) + &
                0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp end parallel

!        ApV = 1._knd/(1._knd+Ap*ApV)
    call multiply_and_add_scalar(ApV, Ap, 1._knd)
    call reciprocal(ApV, 1._knd)

    !$omp parallel private(i,j,k,bi,bj,bk,Suavg,Svavg,Swavg)
    !$omp do schedule(runtime) collapse(3)
    do bk = 1, Wnz, tnz
     do bj = 1, Wny, tny
      do bi = 1, Wnx, tnx
       do k = bk, min(bk+tnz-1,Wnz)
        do j = bj, min(bj+tny-1,Wny)
         do i = bi, min(bi+tnx-1,Wnx)
            ApW(i,j,k) = 0
            if (Wflx_mask(i+1,j,k)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * recdxmin2
            if (Wflx_mask(i,j,k)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * recdxmin2
            if (Wfly_mask(i,j+1,k)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * recdymin2
            if (Wfly_mask(i,j,k)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * recdymin2
            if (Wflz_mask(i,j,k+1)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                nu(i,j,k+1) * recdzmin2
            if (Wflz_mask(i,j,k)) &
              ApW(i,j,k) = ApW(i,j,k) + &
                nu(i,j,k) * recdzmin2
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
    !$omp end parallel


!        ApW = 1._knd/(1._knd+Ap*ApW)
    call multiply_and_add_scalar(ApW, Ap, 1._knd)
    call reciprocal(ApW, 1._knd)


    do l = 1, maxCNiter               !Gauss-Seidel iteration for Crank-Nicolson result
      call BoundU(1,U3,Uin)
      call BoundU(2,V3,Vin)
      call BoundU(3,W3,Win)

      S = 0
      Su = 0
      Sv = 0
      Sw = 0
      !$omp parallel private(i,j,k,bi,bj,bk,p) reduction(max:Su,Sv,Sw)
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Unz, tnz2
       do bj = 1, Uny, tny2
        do bi = 1, Unx, tnx2
         do k = bk, min(bk+tnz2-1,Unz)
          do j = bj, min(bj+tny2-1,Uny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Unx), 2
            if (Utype(i,j,k)<=0) then
              wrk(i,j,k) = 0
              if (Uflx_mask(i+1,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i+1,j,k) * U3(i+1,j,k) * recdxmin2
              if (Uflx_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i,j,k) * U3(i-1,j,k) * recdxmin2
              if (Ufly_mask(i,j+1,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * U3(i,j+1,k) * recdymin2
              if (Ufly_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * U3(i,j-1,k) * recdymin2
              if (Uflz_mask(i,j,k+1)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * U3(i,j,k+1) * recdzmin2
              if (Uflz_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * U3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 1
#include "wmfluxes-inc.f90"
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Unz, tnz2
       do bj = 1, Uny, tny2
        do bi = 1, Unx, tnx2
         do k = bk, min(bk+tnz2-1,Unz)
          do j = bj, min(bj+tny2-1,Uny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Unx), 2
            if (Utype(i,j,k)<=0) then
              p = Ap * wrk(i,j,k) + U2(i,j,k) + U(i,j,k)
              p = p * ApU(i,j,k)
              Su = max(Su,abs(p-U3(i,j,k)))
              U3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do

      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Vnz, tnz2
       do bj = 1, Vny, tny2
        do bi = 1, Vnx, tnx2
         do k = bk, min(bk+tnz2-1,Vnz)
          do j = bj, min(bj+tny2-1,Vny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Vnx), 2
            if (Vtype(i,j,k)<=0) then
              wrk(i,j,k) = 0
              if (Vflx_mask(i+1,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * V3(i+1,j,k) * recdxmin2
              if (Vflx_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * V3(i-1,j,k) * recdxmin2
              if (Vfly_mask(i,j+1,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i,j+1,k) * V3(i,j+1,k) * recdymin2
              if (Vfly_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i,j,k) * V3(i,j-1,k) * recdymin2
              if (Vflz_mask(i,j,k+1)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * V3(i,j,k+1) * recdzmin2
              if (Vflz_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * V3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 2
#include "wmfluxes-inc.f90"
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Vnz, tnz2
       do bj = 1, Vny, tny2
        do bi = 1, Vnx, tnx2
         do k = bk, min(bk+tnz2-1,Vnz)
          do j = bj, min(bj+tny2-1,Vny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Vnx), 2
            if (Vtype(i,j,k)<=0) then
              p = Ap * wrk(i,j,k) + V2(i,j,k) + V(i,j,k)
              p = p * ApV(i,j,k)
              Sv = max(Sv,abs(p-V3(i,j,k)))
              V3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
         
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Wnz, tnz2
       do bj = 1, Wny, tny2
        do bi = 1, Wnx, tnx2
         do k = bk, min(bk+tnz2-1,Wnz)
          do j = bj, min(bj+tny2-1,Wny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Wnx), 2
            if (Wtype(i,j,k)<=0) then
              wrk(i,j,k) = 0
              if (Wflx_mask(i+1,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * W3(i+1,j,k) * recdxmin2
              if (Wflx_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * W3(i-1,j,k) * recdxmin2
              if (Wfly_mask(i,j+1,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * W3(i,j+1,k) * recdymin2
              if (Wfly_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * W3(i,j-1,k) * recdymin2
              if (Wflz_mask(i,j,k+1)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i,j,k+1) * W3(i,j,k+1) * recdzmin2
              if (Wflz_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i,j,k) * W3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 3
#include "wmfluxes-inc.f90"
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Wnz, tnz2
       do bj = 1, Wny, tny2
        do bi = 1, Wnx, tnx2
         do k = bk, min(bk+tnz2-1,Wnz)
          do j = bj, min(bj+tny2-1,Wny)
           do i = bi+mod(bi+j+k-1,2), min(bi+tnx2-1,Wnx), 2
            if (Wtype(i,j,k)<=0) then
              p = Ap * wrk(i,j,k) + W2(i,j,k) + W(i,j,k)
              p = p * ApW(i,j,k)
              Sw = max(Sw,abs(p-W3(i,j,k)))
              W3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do


      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Unz, tnz2
       do bj = 1, Uny, tny2
        do bi = 1, Unx, tnx2
         do k = bk, min(bk+tnz2-1,Unz)
          do j = bj, min(bj+tny2-1,Uny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Unx), 2
            if (Utype(i,j,k)<=0) then
              wrk(i,j,k) = 0
              if (Uflx_mask(i+1,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i+1,j,k) * U3(i+1,j,k) * recdxmin2
              if (Uflx_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i,j,k) * U3(i-1,j,k) * recdxmin2
              if (Ufly_mask(i,j+1,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * U3(i,j+1,k) * recdymin2
              if (Ufly_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * U3(i,j-1,k) * recdymin2
              if (Uflz_mask(i,j,k+1)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * U3(i,j,k+1) * recdzmin2
              if (Uflz_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * U3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 1
#include "wmfluxes-inc.f90"
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Unz, tnz2
       do bj = 1, Uny, tny2
        do bi = 1, Unx, tnx2
         do k = bk, min(bk+tnz2-1,Unz)
          do j = bj, min(bj+tny2-1,Uny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Unx), 2
            if (Utype(i,j,k)<=0) then
              p = Ap * wrk(i,j,k) + U2(i,j,k) + U(i,j,k)
              p = p * ApU(i,j,k)
              Su = max(Su,abs(p-U3(i,j,k)))
              U3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
         
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Vnz, tnz2
       do bj = 1, Vny, tny2
        do bi = 1, Vnx, tnx2
         do k = bk, min(bk+tnz2-1,Vnz)
          do j = bj, min(bj+tny2-1,Vny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Vnx), 2
            if (Vtype(i,j,k)<=0) then
              wrk(i,j,k) = 0
              if (Vflx_mask(i+1,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * V3(i+1,j,k) * recdxmin2
              if (Vflx_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * V3(i-1,j,k) * recdxmin2
              if (Vfly_mask(i,j+1,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i,j+1,k) * V3(i,j+1,k) * recdymin2
              if (Vfly_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i,j,k) * V3(i,j-1,k) * recdymin2
              if (Vflz_mask(i,j,k+1)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * V3(i,j,k+1) * recdzmin2
              if (Vflz_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * V3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 2
#include "wmfluxes-inc.f90"
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Vnz, tnz2
       do bj = 1, Vny, tny2
        do bi = 1, Vnx, tnx2
         do k = bk, min(bk+tnz2-1,Vnz)
          do j = bj, min(bj+tny2-1,Vny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Vnx), 2
            if (Vtype(i,j,k)<=0) then
              p = Ap * wrk(i,j,k) + V2(i,j,k) + V(i,j,k)
              p = p * ApV(i,j,k)
              Sv = max(Sv,abs(p-V3(i,j,k)))
              V3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
         
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Wnz, tnz2
       do bj = 1, Wny, tny2
        do bi = 1, Wnx, tnx2
         do k = bk, min(bk+tnz2-1,Wnz)
          do j = bj, min(bj+tny2-1,Wny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Wnx), 2
            if (Wtype(i,j,k)<=0) then
              wrk(i,j,k) = 0
              if (Wflx_mask(i+1,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * W3(i+1,j,k) * recdxmin2
              if (Wflx_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * W3(i-1,j,k) * recdxmin2
              if (Wfly_mask(i,j+1,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * W3(i,j+1,k) * recdymin2
              if (Wfly_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * W3(i,j-1,k) * recdymin2
              if (Wflz_mask(i,j,k+1)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i,j,k+1) * W3(i,j,k+1) * recdzmin2
              if (Wflz_mask(i,j,k)) &
                wrk(i,j,k) = wrk(i,j,k) + &
                  nu(i,j,k) * W3(i,j,k-1) * recdzmin2
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
#define comp 3
#include "wmfluxes-inc.f90"
#undef comp
      !$omp do schedule(runtime) collapse(3)
      do bk = 1, Wnz, tnz2
       do bj = 1, Wny, tny2
        do bi = 1, Wnx, tnx2
         do k = bk, min(bk+tnz2-1,Wnz)
          do j = bj, min(bj+tny2-1,Wny)
           do i = bi+mod(bi+j+k,2), min(bi+tnx2-1,Wnx), 2
            if (Wtype(i,j,k)<=0) then
              p = Ap * wrk(i,j,k) + W2(i,j,k) + W(i,j,k)
              p = p * ApW(i,j,k)
              Sw = max(Sw,abs(p-W3(i,j,k)))
              W3(i,j,k) = p
            end if
           end do
          end do
         end do
        end do
       end do
      end do
      !$omp end do
      !$omp end parallel

      S = max(Su,Sv,Sw)
#ifdef PAR
      S = par_co_max(S)
#endif
      if (master) write(*,*) "CN ", l, S

      if (S<=epsCN) exit
    end do


  end subroutine ImplicitDiffusion_Iterations







end module Dynamics
module TSTEPS

  use PARAMETERS
  use ArrayUtilities
  use CDS!, only: CDU, CDV, CDW, CDS4U, CDS4V, CDS4W, CDS4_2U, CDS4_2V, CDS4_2W
  use BOUNDARIES, only: BoundU,Bound_Q
  use ScalarBoundaries, only: BoundTemperature, BoundViscosity
  use Pressure !it exports PressureCorrection and GetPrFromGPU
  use OUTPUTS, only: store,display,proftempfl,profmoistfl
  use SCALARS, only: ScalarRK3, ComputeTDiff
  use Subgrid, only: sgstype, SGS_Smag, SGS_StabSmag, SGS_Vreman, SGS_Sigma
  use TURBINLET, only: GetTurbulentInlet, GetInletFromFile
  use Wallmodels
  use Tiling, only: tilenx, tileny, tilenz
#ifdef __HMPP
  use HMPP_CODELETS
#endif

  implicit none


  private
  public TMarchRK3, TSteps_Deallocate

#ifdef __HMPP
  public GetDataFromGPU
  type(hmppWMPoint),allocatable :: hmppWMPoints(:)
  integer :: nWMPoints = 0
#endif

  !module variables to allow their deallocation before programend
  real(knd), dimension(:,:,:), allocatable:: U3,V3,W3

  real(knd), dimension(:,:,:), allocatable :: Q
  real(knd), dimension(:,:,:), allocatable :: U2,Ustar
  real(knd), dimension(:,:,:), allocatable :: V2,Vstar
  real(knd), dimension(:,:,:), allocatable :: W2,Wstar

  real(knd), dimension(:,:,:), allocatable :: Apu, ApV, ApW, wrk

  real(knd), dimension(:), allocatable :: Uwm, Vwm, Wwm

contains

#ifndef __HMPP
#include "tsteps-shared.f90"
#endif

 subroutine TSteps_Deallocate
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


 subroutine TMarchRK3(U,V,W,Pr,Temperature,Moisture,Scalar,dt,delta)
  use RK3
  real(knd),allocatable,intent(inout) :: U(:,:,:),V(:,:,:),W(:,:,:),Pr(:,:,:)
  real(knd),allocatable,intent(inout) :: Temperature(:,:,:),Moisture(:,:,:),Scalar(:,:,:,:)
  real(knd),intent(out) :: dt, delta

  integer RK_stage
  integer,save:: called = 0
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


!    !$omp parallel sections
!    !$omp section
   call BoundU(1,U,Uin)
!    !$omp section
   call BoundU(2,V,Vin)
!    !$omp section
   call BoundU(3,W,Win)
!    !$omp section
   if (enable_buoyancy) call BoundTemperature(temperature)
!    !$omp end parallel sections
    call IBMomentum(U,V,W)


   if (masssourc==1) allocate(Q(0:Prnx+1,0:Prny+1,0:Prnz+1))

   if (debugparam>1) call system_clock(count_rate=trate)

#ifdef __HMPP
   if (allocated(WMPoints)) then
     nWMPoints = size(WMPoints)
   else
     nWMPoints = 0
   end if
write(*,*) "nhWMPointsHMPP:",nWMPoints
   allocate(hmppWMPoints(nWMPoints))

   if (allocated(WMPoints)) call WMPtoHMPP(hmppWMPoints,WMPoints)


    !$hmpp <tsteps> allocate
    !$hmpp <tsteps> advancedload, args[Vreman::Prnx,BoundU::Prny,BoundU::Prnz]
    !$hmpp <tsteps> advancedload, args[BoundU::nx,BoundU::ny,BoundU::nz, &
    !$hmpp  &     BoundV::nx,BoundV::ny,BoundV::nz,BoundW::nx,BoundW::ny,BoundW::nz]

    !$hmpp <tsteps> advancedload, args[BoundU::sideU,BoundU::Uin,BoundV::Uin,BoundW::Uin]

    !$hmpp <tsteps> advancedload, args[BoundU::component,BoundV::component,BoundW::component]
    !$hmpp <tsteps> advancedload, args[BoundU::regime,BoundV::regime,BoundW::regime]

    !$hmpp <tsteps> advancedload, args[Vreman::dx,Vreman::dy,Vreman::dz,Vreman::Re]
    !$hmpp <tsteps> advancedload, args[Convection::enable_buoyancy,Convection::convmet,Convection::coriolisparam]
    !$hmpp <tsteps> advancedload, args[Convection::grav_acc,Convection::temperature_ref,Convection::RK_beta,Convection::RK_rho]
    !$hmpp <tsteps> advancedload, args[PressureGrad::prgradientx,PressureGrad::prgradienty]

    !$hmpp <tsteps> advancedload, args[BoundU::Btype]

    !$hmpp <tsteps> advancedload, args[ForwEul::dxPr,ForwEul::dyPr,ForwEul::dzPr]
!     !$hmpp <tsteps> advancedload, args[AttenuateOut::xPr,AttenuateTop::zPr]

    !$hmpp <tsteps> advancedload, args[PressureGrad::dxU,PressureGrad::dyV,PressureGrad::dzW]
!     !$hmpp <tsteps> advancedload, args[AttenuateOut::xU,AttenuateTop::zW]

    !$hmpp <tsteps> advancedload, args[PressureGrad::Pr]

    !$hmpp <tsteps> advancedload, args[NullInterior::Unull,NullInterior::Vnull,NullInterior::Wnull]
    !$hmpp <tsteps> advancedload, args[NullInterior::nUnull,NullInterior::nVnull,NullInterior::nWnull]

    !$hmpp <tsteps> advancedload, args[TimeStepEul::CFL,TimeStepEul::Uref,TimeStepEul::steady,TimeStepEul::dxmin,TimeStepEul::endtime]

    !$hmpp <tsteps>  advancedload, args[UnifRedBlack::maxCNiter,UnifRedBlack::epsCN]

    !$hmpp <tsteps> advancedload, args[TimeStepEul::U,TimeStepEul::V,TimeStepEul::W]

    if (enable_buoyancy) then
     !$hmpp <tsteps> advancedload, args[Convection::temperature]
    end if

    !$hmpp <tsteps> advancedload, args[ComputeViscWM::nWMPoints,ComputeViscWM::WMPoints]
    !$hmpp <tsteps> advancedload, args[BoundViscosity_TDiff::Prandtl]

    GPU = 1
#endif
  end if

  if ((Btype(We) ==TurbulentInlet).or.(Btype(Ea) ==TurbulentInlet)) then
    call GetTurbulentInlet(dt)
    !$hmpp <tsteps> advancedload, args[BoundU::sideU,BoundU::Uin,BoundV::Uin,BoundW::Uin]
  else if (Btype(We) ==InletFromFile) then
    call GetInletFromFile(time)
    !$hmpp <tsteps> advancedload, args[BoundU::sideU,BoundU::Uin,BoundV::Uin,BoundW::Uin]
  end if

  !$hmpp <tsteps> advancedload, args[TimeStepEul::time]
  !$hmpp <tsteps> TimeStepEul callsite, args[*].noupdate = true
  call TimeStepEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz, &
               dxmin,dxU,dyV,dzW,CFL,Uref,steady,time,end_time, &
               U,V,W,dt)
  !$hmpp <tsteps> delegatedstore, args[TimeStepEul::dt]


  if (master) write (*,*) "time:",time,"dt: ",dt

  do RK_stage = 1,RK_stages

    if (master) write(*,*) "stage:",RK_stage


    if (debugparam>1.and.master) call system_clock(count=time1)


    call SubgridStresses(U,V,W,Pr,Temperature)

    !$hmpp <tsteps> advancedload, args[Convection::RK_stage]


    !$hmpp <tsteps> Convection callsite, args[*].noupdate = true
    call Convection(U,V,W,U2,V2,W2,Ustar,Vstar,Wstar,Temperature,Moisture,RK_beta,RK_rho,RK_stage,dt)

    call ScalarRK3(U,V,W,Temperature,Moisture,Scalar,RK_stage,dt,proftempfl,profmoistfl)

    call OtherTerms(U,V,W,U2,V2,W2,Pr,2._knd*RK_alpha(RK_stage)*dt)

    !download U2 V2 and W2 to main memory Attenuates below are to slow on GPU, waiting for HMPP update
!     !$hmpp <tsteps> delegatedstore, args[UnifRedblack::U2,UnifRedblack::V2,UnifRedblack::W2]

    !this will havr no effect when doing HMPP
    if ((Btype(To) ==FreeSlipBuff) .and. (Prnz>15))  then
        !$hmpp <tsteps> AttenuateTop callsite, args[*].noupdate = true
        call AttenuateTop(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Btype, &
                          zPr,zW,U2,V2,W2)
    end if

    if ((Btype(Ea) ==OutletBuff) .and. (Prnx>15)) then
        !$hmpp <tsteps> AttenuateOut callsite, args[*].noupdate = true
        call AttenuateOut(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz, &
                          xPr,xU,U2,V2,W2,temperature,enable_buoyancy)
    end if


!     !$omp parallel sections
!     !$omp section
    call BoundU(1,U2,Uin)
!     !$omp section
    call BoundU(2,V2,Vin)
!     !$omp section
    call BoundU(3,W2,Win)
!     !$omp end parallel sections
    call IBMomentum(U2,V2,W2)

    if (masssourc==1) then
        call IBMassSources(Q,U2,V2,W2)
    end if


    if (poissmet>0) then
       call PressureCorrection(U2,V2,W2,Pr,Q,2._knd*RK_alpha(RK_stage)*dt)
    end if


#ifdef __HMPP
    !$hmpp <tsteps> UpdateU callsite, args[*].noupdate=true
    call UpdateU_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,U,V,W,U2,V2,W2)
#else
    if (RK_stage==1) delta = 0
#ifdef DEBUG
    if (debuglevel>0.or.steady==1) then
      if (Unx*Uny*Unz>0) &
        delta = delta+sum(abs(U(1:Unx,1:Uny,1:Unz)-U2(1:Unx,1:Uny,1:Unz)))/(Unx*Uny*Unz)
      if (Vnx*Vny*Vnz>0) &
        delta = delta+sum(abs(V(1:Vnx,1:Vny,1:Vnz)-V2(1:Vnx,1:Vny,1:Vnz)))/(Vnx*Vny*Vnz)
      if (Wnx*Wny*Wnz>0) &
        delta = delta+sum(abs(W(1:Wnx,1:Wny,1:Wnz)-W2(1:Wnx,1:Wny,1:Wnz)))/(Wnx*Wny*Wnz)
    end if
#endif
    call exchange_alloc(U,U2)
    call exchange_alloc(V,V2)
    call exchange_alloc(W,W2)
#endif


    if ((Btype(Ea) ==OutletBuff) .and. (Prnx>15)) then
     !$hmpp <tsteps> AttenuateOut2 callsite, args[*].noupdate = true
      call AttenuateOut(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz, &
                        xPr,xU,U,V,W,temperature,enable_buoyancy)
!       if (enable_buoyancy) then
! #ifdef __HMPP
! !              !$hmpp <tsteps> BoundTemperature callsite
!              call BoundTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,TempBtype,sideTemp,BsideTArr,BsideTFLArr,TDiff,TempIn,temperature)
! #else
!              call BoundTemperature(temperature)
! #endif
!       end if
    end if


    !$hmpp <tsteps> NullInterior callsite, args[*].noupdate = true
    call NullInterior(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz, &
                           nUnull,nVnull,nWnull,Unull,Vnull,Wnull,U,V,W)


#ifdef __HMPP
    !$hmpp <tsteps> BoundU callsite, args[*].noupdate = true
    call BoundU_GPU(1,Unx,Uny,Unz,Prny,Prnz, &
                         Btype,sideU, &
                         Uin,U,0)

    !$hmpp <tsteps> BoundV callsite, args[*].noupdate = true
    call BoundU_GPU(2,Vnx,Vny,Vnz,Prny,Prnz, &
                         Btype,sideU, &
                         Vin,V,0)

    !$hmpp <tsteps> BoundW callsite, args[*].noupdate = true
    call BoundU_GPU(3,Wnx,Wny,Wnz,Prny,Prnz, &
                         Btype,sideU, &
                         Win,W,0)
#else
! !$omp parallel
! !$omp sections
! !$omp section
        call BoundU(1,U,Uin)
! !$omp section
        call BoundU(2,V,Vin)
! !$omp section
        call BoundU(3,W,Win)
! !$omp end sections
! !$omp end parallel
        call IBMomentum(U,V,W)
#endif


    if (debugparam>1.and.master) then
     call system_clock(count=time2)
     write (*,*) "ET of part 1", (time2-time1)/real(trate,int64)
     time1 = time2
    end if

   end do
   
   if (.false.) then
   !$hmpp <tsteps> release
   end if

end subroutine TMarchRK3



















  subroutine OtherTerms(U,V,W,U2,V2,W2,Pr,coef)
   real(knd), contiguous, intent(inout) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
   real(knd), contiguous, intent(inout) :: Pr(1:,1:,1:)
   real(knd), allocatable, intent(inout) :: U2(:,:,:),V2(:,:,:),W2(:,:,:)
   real(knd), intent(in) :: coef

   real(knd) S

   integer i,j,k
   integer,save:: called=0

   if (called==0) then
     allocate(U3(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),lbound(U,3):ubound(U,3)))
     allocate(V3(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),lbound(V,3):ubound(V,3)))
     allocate(W3(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),lbound(W,3):ubound(W,3)))
     called=1
   end if

   !$hmpp <tsteps> advancedload, args[PressureGrad::coef]

   !Pressure gradient terms
   !$hmpp <tsteps> PressureGrad callsite, args[*].noupdate = true
   call PressureGrad(Prnx,Prny,Prnz, &
                     Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz, &
                     dxU,dyV,dzW, &
                     Btype,prgradientx,prgradienty, &
                     Pr,U2,V2,W2, &
                     coef)

    if (explicit_diffusion<=0) then

       !semi-implicit diffusion

       Re_gt_0: if (Re>0) then

         !Diffusion using Crank Nicolson
         !first approximation using forward Euler
         !iteration SOR or Gauss-Seidel


    !      !$hmpp <tsteps> advancedload, args[ForwEul::Visc]

         !$hmpp <tsteps> ForwEul callsite, args[*].noupdate = true
         call ForwEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz, &
                      dxPr,dyPr,dzPr,dxU,dyV,dzW, &
                      U,V,W,U2,V2,W2,U3,V3,W3,Viscosity, &
                      coef)


         !$omp parallel sections
         !$omp section
         call BoundU(1,U3,Uin)
         !$omp section
         call BoundU(2,V3,Vin)
         !$omp section
         call BoundU(3,W3,Win)
         !$omp end parallel sections
         call IBMomentum(U3,V3,W3)

         !Performs the diffusion terms
#ifdef __HMPP
    !       !$hmpp <tsteps>  advancedload, args[UnifRedBlack::maxCNiter,UnifRedBlack::epsCN]

          !$hmpp <tsteps> UNIFREDBLACK callsite, args[*].noupdate=true
          call UNIFREDBLACK_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz, &
                                Btype,sideU, &
                                dt,dxmin,dymin,dzmin, &
                                Uin,Vin,Win, &
                                U,V,W,U2,V2,W2,U3,V3,W3,Visc, &
                                coef,maxCNiter,epsCN,it,S)

         !$hmpp <tsteps> delegatedstore, args[UnifRedBlack::iters,UnifRedBlack::residuum]

          write(*,*) "back from GPU CN", it,S

#else
         if (gridtype==UNIFORMGRID) then
          call UNIFREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
         else
          call error_stop
         end if
         call exchange_alloc(U2,U3)
         call exchange_alloc(V2,V3)
         call exchange_alloc(W2,W3)
#endif

       else  Re_gt_0  !Re<=0

        U2 = U+U2
        V2 = V+V2
        W2 = W+W2

       end if   Re_gt_0

       !$omp parallel sections
       !$omp section
       call BoundU(1,U2,Uin)
       !$omp section
       call BoundU(2,V2,Vin)
       !$omp section
       call BoundU(3,W2,Win)
       !$omp end parallel sections
       call IBMomentum(U2,V2,W2)

       if (debuglevel>=2) then  !Compute and output the mean friction in the domain.
        S = 0
        do k=1,Unz
         do j=1,Uny
          do i=1,Unx
           S = S-((Viscosity(i+1,j,k) * (U(i+1,j,k)-U(i,j,k))/dxPr(i+1) - &
           Viscosity(i,j,k) * (U(i,j,k)-U(i-1,j,k))/dxPr(i))/dxU(i) + &
             (0.25_knd * (Viscosity(i+1,j+1,k)+Viscosity(i+1,j,k)+Viscosity(i,j+1,k)+Viscosity(i,j,k))* &
                   (U(i,j+1,k)-U(i,j,k))/dyV(j) - &
             0.25_knd * (Viscosity(i+1,j,k)+Viscosity(i+1,j-1,k)+Viscosity(i,j,k)+Viscosity(i,j-1,k))* &
                   (U(i,j,k)-U(i,j-1,k))/dyV(j-1))/dyPr(j) + &
              (0.25_knd * (Viscosity(i+1,j,k+1)+Viscosity(i+1,j,k)+Viscosity(i,j,k+1)+Viscosity(i,j,k))* &
                   (U(i,j,k+1)-U(i,j,k))/dzW(k) - &
             0.25_knd * (Viscosity(i+1,j,k)+Viscosity(i+1,j,k-1)+Viscosity(i,j,k)+Viscosity(i,j,k-1))* &
                   (U(i,j,k)-U(i,j,k-1))/dzW(k-1))/dzPr(k))
          end do
         end do
        end do

        S = S/(Unx*Uny*Unz)
        write(*,*) "Mean friction:", S
       end if

    else !explicit diffusion

        call add(U2, U)
        call add(V2, V)
        call add(W2, W)

    end if

  end subroutine OtherTerms








  
  subroutine UnifRedBlack(U,V,W,U2,V2,W2,U3,V3,W3,coef)
    use Parameters, nu => Viscosity
    !$ use omp_lib
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U2,V2,W2,U3,V3,W3
    real(knd), intent(in) :: coef
    real(knd) recdxmin2,recdymin2,recdzmin2                                                               !reciprocal values of dx**2
    real(knd) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
    integer i,j,k,bi,bj,bk,l
    integer, parameter :: narr = 3, narr2 = 5 !number of arrays in the loop
    integer, save :: called = 0

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

!        !$ call omp_set_nested(.true.)
       !$omp parallel private(i,j,k,bi,bj,bk)

       !The explicit part, which doesn't have to be changed inside the loop
       !$omp do schedule(runtime) !collapse(3)
       do bk = 1, Unz, tilenz(narr)
        do bj = 1, Uny, tileny(narr)
         do bi = 1, Unx, tilenx(narr)
          do k = bk, min(bk+tilenz(narr)-1,Unz)
           do j = bj, min(bj+tileny(narr)-1,Uny)
            do i = bi, min(bi+tilenx(narr)-1,Unx)
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
       !$omp do schedule(runtime) !collapse(3)
       do bk = 1, Unz, tilenz(narr)
        do bj = 1, Uny, tileny(narr)
         do bi = 1, Unx, tilenx(narr)
          do k = bk, min(bk+tilenz(narr)-1,Unz)
           do j = bj, min(bj+tileny(narr)-1,Uny)
            do i = bi, min(bi+tilenx(narr)-1,Unx)
              U2(i,j,k) = U2(i,j,k) + Ap * wrk(i,j,k)
            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end do nowait
       !$omp do schedule(runtime) !collapse(3)
       do bk = 1, Vnz, tilenz(narr)
        do bj = 1, Vny, tileny(narr)
         do bi = 1, Vnx, tilenx(narr)
          do k = bk, min(bk+tilenz(narr)-1,Vnz)
           do j = bj, min(bj+tileny(narr)-1,Vny)
            do i = bi, min(bi+tilenx(narr)-1,Vnx)
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
       !$omp do schedule(runtime) !collapse(3)
       do bk = 1, Vnz, tilenz(narr)
        do bj = 1, Vny, tileny(narr)
         do bi = 1, Vnx, tilenx(narr)
          do k = bk, min(bk+tilenz(narr)-1,Vnz)
           do j = bj, min(bj+tileny(narr)-1,Vny)
            do i = bi, min(bi+tilenx(narr)-1,Vnx)
              V2(i,j,k) = V2(i,j,k) + Ap * wrk(i,j,k)
            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end do
       !$omp do schedule(runtime) !collapse(3)
       do bk = 1, Wnz, tilenz(narr)
        do bj = 1, Wny, tileny(narr)
         do bi = 1, Wnx, tilenx(narr)
          do k = bk, min(bk+tilenz(narr)-1,Wnz)
           do j = bj, min(bj+tileny(narr)-1,Wny)
            do i = bi, min(bi+tilenx(narr)-1,Wnx)
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
       do bk = 1, Wnz, tilenz(narr)
        do bj = 1, Wny, tileny(narr)
         do bi = 1, Wnx, tilenx(narr)
          do k = bk, min(bk+tilenz(narr)-1,Wnz)
           do j = bj, min(bj+tileny(narr)-1,Wny)
            do i = bi, min(bi+tilenx(narr)-1,Wnx)
              W2(i,j,k) = W2(i,j,k) + Ap * wrk(i,j,k)
            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end do

       !Auxiliary coefficients to better efficiency in loops
       !$omp do schedule(runtime) !collapse(3)
       do bk = 1, Unz, tilenz(narr)
        do bj = 1, Uny, tileny(narr)
         do bi = 1, Unx, tilenx(narr)
          do k = bk, min(bk+tilenz(narr)-1,Unz)
           do j = bj, min(bj+tileny(narr)-1,Uny)
            do i = bi, min(bi+tilenx(narr)-1,Unx)
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
       !$omp do schedule(runtime) !collapse(3)
       do bk = 1, Vnz, tilenz(narr)
        do bj = 1, Vny, tileny(narr)
         do bi = 1, Vnx, tilenx(narr)
          do k = bk, min(bk+tilenz(narr)-1,Vnz)
           do j = bj, min(bj+tileny(narr)-1,Vny)
            do i = bi, min(bi+tilenx(narr)-1,Vnx)
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
       !$omp do schedule(runtime) !collapse(3)
       do bk = 1, Wnz, tilenz(narr)
        do bj = 1, Wny, tileny(narr)
         do bi = 1, Wnx, tilenx(narr)
          do k = bk, min(bk+tilenz(narr)-1,Wnz)
           do j = bj, min(bj+tileny(narr)-1,Wny)
            do i = bi, min(bi+tilenx(narr)-1,Wnx)
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


       do l=1,maxCNiter               !Gauss-Seidel iteration for Crank-Nicolson result
         call BoundU(1,U3,Uin)
         call BoundU(2,V3,Vin)
         call BoundU(3,W3,Win)

         S = 0
         Su = 0
         Sv = 0
         Sw = 0
         !$omp parallel private(i,j,k,bi,bj,bk,p) reduction(max:Su,Sv,Sw)
         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Unz, tilenz(narr2)
          do bj = 1, Uny, tileny(narr2)
           do bi = 1, Unx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Unz)
             do j = bj, min(bj+tileny(narr2)-1,Uny)
              do i = bi+mod(bi+j+k-1,2), min(bi+tilenx(narr2)-1,Unx), 2
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
         !$omp enddo
#define comp 1
#include "wmfluxes-inc.f90"
#undef comp
         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Unz, tilenz(narr2)
          do bj = 1, Uny, tileny(narr2)
           do bi = 1, Unx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Unz)
             do j = bj, min(bj+tileny(narr2)-1,Uny)
              do i = bi+mod(bi+j+k-1,2), min(bi+tilenx(narr2)-1,Unx), 2
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
         !$omp enddo

         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Vnz, tilenz(narr2)
          do bj = 1, Vny, tileny(narr2)
           do bi = 1, Vnx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Vnz)
             do j = bj, min(bj+tileny(narr2)-1,Vny)
              do i = bi+mod(bi+j+k-1,2), min(bi+tilenx(narr2)-1,Vnx), 2
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
         !$omp enddo
#define comp 2
#include "wmfluxes-inc.f90"
#undef comp
         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Vnz, tilenz(narr2)
          do bj = 1, Vny, tileny(narr2)
           do bi = 1, Vnx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Vnz)
             do j = bj, min(bj+tileny(narr2)-1,Vny)
              do i = bi+mod(bi+j+k-1,2), min(bi+tilenx(narr2)-1,Vnx), 2
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
         !$omp enddo
         
         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Wnz, tilenz(narr2)
          do bj = 1, Wny, tileny(narr2)
           do bi = 1, Wnx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Wnz)
             do j = bj, min(bj+tileny(narr2)-1,Wny)
              do i = bi+mod(bi+j+k-1,2), min(bi+tilenx(narr2)-1,Wnx), 2
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
         !$omp enddo
#define comp 3
#include "wmfluxes-inc.f90"
#undef comp
         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Wnz, tilenz(narr2)
          do bj = 1, Wny, tileny(narr2)
           do bi = 1, Wnx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Wnz)
             do j = bj, min(bj+tileny(narr2)-1,Wny)
              do i = bi+mod(bi+j+k-1,2), min(bi+tilenx(narr2)-1,Wnx), 2
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
         !$omp enddo


         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Unz, tilenz(narr2)
          do bj = 1, Uny, tileny(narr2)
           do bi = 1, Unx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Unz)
             do j = bj, min(bj+tileny(narr2)-1,Uny)
              do i = bi+mod(bi+j+k,2), min(bi+tilenx(narr2)-1,Unx), 2
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
         !$omp enddo
#define comp 1
#include "wmfluxes-inc.f90"
#undef comp
         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Unz, tilenz(narr2)
          do bj = 1, Uny, tileny(narr2)
           do bi = 1, Unx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Unz)
             do j = bj, min(bj+tileny(narr2)-1,Uny)
              do i = bi+mod(bi+j+k,2), min(bi+tilenx(narr2)-1,Unx), 2
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
         !$omp enddo
         
         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Vnz, tilenz(narr2)
          do bj = 1, Vny, tileny(narr2)
           do bi = 1, Vnx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Vnz)
             do j = bj, min(bj+tileny(narr2)-1,Vny)
              do i = bi+mod(bi+j+k,2), min(bi+tilenx(narr2)-1,Vnx), 2
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
         !$omp enddo
#define comp 2
#include "wmfluxes-inc.f90"
#undef comp
         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Vnz, tilenz(narr2)
          do bj = 1, Vny, tileny(narr2)
           do bi = 1, Vnx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Vnz)
             do j = bj, min(bj+tileny(narr2)-1,Vny)
              do i = bi+mod(bi+j+k,2), min(bi+tilenx(narr2)-1,Vnx), 2
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
         !$omp enddo
         
         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Wnz, tilenz(narr2)
          do bj = 1, Wny, tileny(narr2)
           do bi = 1, Wnx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Wnz)
             do j = bj, min(bj+tileny(narr2)-1,Wny)
              do i = bi+mod(bi+j+k,2), min(bi+tilenx(narr2)-1,Wnx), 2
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
         !$omp enddo
#define comp 3
#include "wmfluxes-inc.f90"
#undef comp
         !$omp do schedule(runtime) !collapse(3)
         do bk = 1, Wnz, tilenz(narr2)
          do bj = 1, Wny, tileny(narr2)
           do bi = 1, Wnx, tilenx(narr2)
            do k = bk, min(bk+tilenz(narr2)-1,Wnz)
             do j = bj, min(bj+tileny(narr2)-1,Wny)
              do i = bi+mod(bi+j+k,2), min(bi+tilenx(narr2)-1,Wnx), 2
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
         !$omp enddo
         !$omp end parallel
         S = max(Su,Sv,Sw)
         write (*,*) "CN ",l,S

         if (S<=epsCN) exit
       end do


  end subroutine UnifRedBlack



  subroutine GENREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
   use Parameters, nu => Viscosity
   real(knd), contiguous, dimension(-2:,-2:,-2:), intent(in) :: U,V,W
   real(knd), contiguous, dimension(-2:,-2:,-2:), intent(inout) :: U2,V2,W2,U3,V3,W3
   real(knd), intent(in) :: coef
   real(knd) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
   integer i,j,k,l
   integer, save :: called = 0

       if (called==0) then
         allocate(Apu(1:Unx,1:Uny,1:Unz))
         allocate(ApV(1:Vnx,1:Vny,1:Vnz))
         allocate(ApW(1:Wnx,1:Wny,1:Wnz))
         called = 1
       end if



       Ap = coef/(2._knd)
       S = 0
       l = 0


       do k=1,Unz    !The explicit part, which doesn't have to be changed inside the loop
        do j=1,Uny
         do i=1,Unx
          U2(i,j,k) = U2(i,j,k)+Ap*( &
          ((nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k))/dxPr(i+1) - &
          nu(i,j,k) * (U(i,j,k)-U(i-1,j,k))/dxPr(i))/dxU(i) + &
           (0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k))* &
                     (U(i,j+1,k)-U(i,j,k))/dyV(j) - &
           0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k))* &
                    (U(i,j,k)-U(i,j-1,k))/dyV(j-1))/dyPr(j) + &
           (0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k))* &
                     (U(i,j,k+1)-U(i,j,k))/dzW(k) - &
           0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1))* &
                    (U(i,j,k)-U(i,j,k-1))/dzW(k-1))/dzPr(k)))
         end do
        end do
       end do
       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          V2(i,j,k) = V2(i,j,k)+Ap*( &
          ((0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k))* &
                     (V(i+1,j,k)-V(i,j,k))/dxU(i) - &
          0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k))* &
                   (V(i,j,k)-V(i-1,j,k))/dxU(i-1))/dxPr(i) + &
           (nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k))/dyPr(j+1) - &
           nu(i,j,k) * (V(i,j,k)-V(i,j-1,k))/dyPr(j))/dyV(j) + &
           (0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k))* &
                     (V(i,j,k+1)-V(i,j,k))/dzW(k) - &
           0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1))* &
                    (V(i,j,k)-V(i,j,k-1))/dzW(k-1))/dzPr(k)))
         end do
        end do
       end do
       do k=1,Wnz
       do j=1,Wny
        do i=1,Wnx
         W2(i,j,k) = W2(i,j,k)+Ap*( &
         ((0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k))* &
                    (W(i+1,j,k)-W(i,j,k))/dxU(i) - &
         0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k))* &
                  (W(i,j,k)-W(i-1,j,k))/dxU(i-1))/dxPr(i) + &
          (0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k))* &
                    (W(i,j+1,k)-W(i,j,k))/dyV(j) - &
          0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k))* &
                   (W(i,j,k)-W(i,j-1,k))/dyV(j-1))/dyPr(j) + &
          (nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k))/dzPr(k+1) - &
          nu(i,j,k) * (W(i,j,k)-W(i,j,k-1))/dzPr(k))/dzW(k)))
        end do
       end do
      end do

       do k=1,Unz         !Auxiliary coefficients to better efficiency in loops
        do j=1,Uny
         do i=1,Unx
          ApU(i,j,k) = ((nu(i+1,j,k)/dxPr(i+1) + &
                      nu(i,j,k)/dxPr(i))/dxU(i) + &
                      (0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) &
                        /dyV(j) + &
                      0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) &
                        /dyV(j-1))/dyPr(j) + &
                      (0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) &
                        /dzW(k) + &
                      0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1))/ &
                        dzW(k-1))/dzPr(k))
         end do
        end do
       end do

      ApU = 1._knd/(1._knd+Ap*ApU)


       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          ApV(i,j,k) =  ((0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) &
                            /dxU(i) + &
                      0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) &
                        /dxU(i-1))/dxPr(i) + &
                     (nu(i,j+1,k)/dyPr(j+1) + &
                     nu(i,j,k)/dyPr(j))/dyV(j) + &
                     (0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) &
                        /dzW(k) + &
                     0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) &
                       /dzW(k-1))/dzPr(k))
         end do
        end do
       end do

       ApV = 1._knd/(1._knd+Ap*ApV)


       do k=1,Wnz
        do j=1,Wny
         do i=1,Wnx
          ApW(i,j,k) =  ((0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) &
                            /dxU(i) + &
                      0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) &
                            /dxU(i-1))/dxPr(i) + &
                     (0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) &
                            /dyV(j) + &
                     0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) &
                            /dyV(j-1))/dyPr(j) + &
                     (nu(i,j,k+1)/dzPr(k+1) + &
                     nu(i,j,k)/dzPr(k))/dzW(k))
         end do
        end do
       end do

       ApW = 1._knd/(1._knd+Ap*ApW)

       Suavg = abs(MAXVAL(U3(1:Unx,1:Uny,1:Unz)))  !maximum values of velocities to norm the residues.
       Svavg = abs(MAXVAL(V3(1:Vnx,1:Vny,1:Vnz)))
       Swavg = abs(MAXVAL(W3(1:Wnx,1:Wny,1:Wnz)))
       if (Suavg<=1e-3_knd) Suavg = 1
       if (Svavg<=1e-3_knd) Svavg = 1
       if (Swavg<=1e-3_knd) Swavg = 1



       do l=1,maxCNiter               !Gauss-Seidel iteration for Crank-Nicolson result
        call BoundU(1,U3,Uin)
        call BoundU(2,V3,Vin)
        call BoundU(3,W3,Win)
        S = 0
        Su = 0
        Sv = 0
        Sw = 0
        !$omp parallel private(i,j,k,p) reduction(max:Su,Sv,Sw)
        !$omp do
        do k=1,Unz
         do j=1,Uny
          do i=1+mod(j+k,2),Unx, 2
            p=((nu(i+1,j,k) * (U3(i+1,j,k))/dxPr(i+1) - &
             nu(i,j,k) * (-U3(i-1,j,k))/dxPr(i))/dxU(i) + &
             (0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k))* &
                       (U3(i,j+1,k))/dyV(j) - &
             0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k))* &
                      (-U3(i,j-1,k))/dyV(j-1))/dyPr(j) + &
             (0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k))* &
                      (U3(i,j,k+1))/dzW(k) - &
             0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1))* &
                      (-U3(i,j,k-1))/dzW(k-1))/dzPr(k))
            p = Ap*p+U2(i,j,k)+U(i,j,k)
            p = p*ApU(i,j,k)


            Su = max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k) = U3(i,j,k)+(p-U3(i,j,k))!*1.72_knd
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do
        do k=1,Vnz
         do j=1,Vny
          do i=1+mod(j+k,2),Vnx, 2
            p=((0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k))* &
                         (V3(i+1,j,k))/dxU(i) - &
             0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k))* &
                      (-V3(i-1,j,k))/dxU(i-1))/dxPr(i) + &
             (nu(i,j+1,k) * (V3(i,j+1,k))/dyPr(j+1) - &
             nu(i,j,k) * (-V3(i,j-1,k))/dyPr(j))/dyV(j) + &
             (0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k))* &
                       (V3(i,j,k+1))/dzW(k) - &
             0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1))* &
                      (-V3(i,j,k-1))/dzW(k-1))/dzPr(k))
            p = Ap*p+V2(i,j,k)+V(i,j,k)
            p = p*ApV(i,j,k)
            Sv = max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k) = V3(i,j,k)+(p-V3(i,j,k))!*1.72_knd
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do
        do k=1,Wnz
         do j=1,Wny
          do i=1+mod(j+k,2),Wnx, 2
            p=((0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k))* &
                         (W3(i+1,j,k))/dxU(i) - &
             0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k))* &
                      (-W3(i-1,j,k))/dxU(i-1))/dxPr(i) + &
             (0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k))* &
                       (W3(i,j+1,k))/dyV(j) - &
             0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k))* &
                      (-W3(i,j-1,k))/dyV(j-1))/dyPr(j) + &
             (nu(i,j,k+1) * (W3(i,j,k+1))/dzPr(k+1) - &
             nu(i,j,k) * (-W3(i,j,k-1))/dzPr(k))/dzW(k))
            p = Ap*p+W2(i,j,k)+W(i,j,k)
            p = p*ApW(i,j,k)
            Sw = max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k) = W3(i,j,k)+(p-W3(i,j,k))!*1.72_knd
          end do
         end do
        end do
        !$omp enddo

        !$omp do
        do k=1,Unz
         do j=1,Uny
          do i=1+mod(j+k+1,2),Unx, 2
            p=((nu(i+1,j,k) * (U3(i+1,j,k))/dxPr(i+1) - &
             nu(i,j,k) * (-U3(i-1,j,k))/dxPr(i))/dxU(i) + &
             (0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k))* &
                       (U3(i,j+1,k))/dyV(j) - &
             0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k))* &
                      (-U3(i,j-1,k))/dyV(j-1))/dyPr(j) + &
             (0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k))* &
                       (U3(i,j,k+1))/dzW(k) - &
             0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1))* &
                      (-U3(i,j,k-1))/dzW(k-1))/dzPr(k))
            p = Ap*p+U2(i,j,k)+U(i,j,k)
            p = p*ApU(i,j,k)


            Su = max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k) = U3(i,j,k)+(p-U3(i,j,k))!*1.72_knd
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do
        do k=1,Vnz
         do j=1,Vny
          do i=1+mod(j+k+1,2),Vnx, 2
            p=((0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k))* &
                         (V3(i+1,j,k))/dxU(i) - &
             0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k))* &
                      (-V3(i-1,j,k))/dxU(i-1))/dxPr(i) + &
             (nu(i,j+1,k) * (V3(i,j+1,k))/dyPr(j+1) - &
             nu(i,j,k) * (-V3(i,j-1,k))/dyPr(j))/dyV(j) + &
             (0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k))* &
                       (V3(i,j,k+1))/dzW(k) - &
             0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1))* &
                      (-V3(i,j,k-1))/dzW(k-1))/dzPr(k))
            p = Ap*p+V2(i,j,k)+V(i,j,k)
            p = p*ApV(i,j,k)
            Sv = max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k) = V3(i,j,k)+(p-V3(i,j,k))!*1.72_knd
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do
        do k=1,Wnz
         do j=1,Wny
          do i=1+mod(j+k+1,2),Wnx, 2
            p=((0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k))* &
                         (W3(i+1,j,k))/dxU(i) - &
             0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k))* &
                      (-W3(i-1,j,k))/dxU(i-1))/dxPr(i) + &
             (0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k))* &
                       (W3(i,j+1,k))/dyV(j) - &
             0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k))* &
                      (-W3(i,j-1,k))/dyV(j-1))/dyPr(j) + &
             (nu(i,j,k+1) * (W3(i,j,k+1))/dzPr(k+1) - &
             nu(i,j,k) * (-W3(i,j,k-1))/dzPr(k))/dzW(k))
            p = Ap*p+W2(i,j,k)+W(i,j,k)
            p = p*ApW(i,j,k)
            Sw = max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k) = W3(i,j,k)+(p-W3(i,j,k))!*1.72_knd
          end do
         end do
        end do
        !$omp enddo
        !$omp end parallel
        S = max(Su/Suavg,Sv/Svavg,Sw/Swavg)
        write (*,*) "CN ",l,S
        if (S<=epsCN) exit
       end do

   end subroutine GENREDBLACK





  subroutine SubgridStresses(U,V,W,Pr,Temperature)
    use ImmersedBoundary, only: ScalFlIBPoints, TIBPoint_Viscosity
    use Filters, only: filtertype, filter_ratios

    real(knd), contiguous, intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(knd), contiguous, intent(in) :: Pr(1:,1:,1:)
    real(knd), contiguous, intent(in) :: Temperature(-1:,-1:,-1:)
    integer i

#ifdef __HMPP

     !$hmpp <tsteps> Vreman callsite, args[*].noupdate = true
     call Vreman_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxmin,dymin,dzmin,dt,Re,U,V,W,Viscosity)
!      !$hmpp <tsteps> delegatedstore, args[Vreman::Visc]

2._knd
     if (wallmodeltype>0) then

!        !$hmpp <tsteps> delegatedstore, args[AttenuateTop::Temperature], &
!        !$hmpp  & args[AttenuateTop::Temperature].section={-1:Prnx+2,-1:Prny+2,1:1} of {-1:Prnx+2,-1:Prny+2,-1:Prnz+2}

       !$hmpp <tsteps> ComputeViscWM callsite, args[*].noupdate = true
       call ComputeViscsWM_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz, &
                                nWMPoints,hmppWMPoints, &
                                TempBtype,Re,temperature_ref,grav_acc, &
                                U,V,W,Viscosity)
     end if


     !boundaries to Visc, compute TDiff and boundaries to TDiff

      !$hmpp <tsteps> BoundViscosity_TDiff callsite, args[*].noupdate = true
      call BoundViscosity_TDiff(Prnx,Prny,Prnz,enable_buoyancy,Btype,Re,Prandtl,Visc,TDiff)



#else
     if (wallmodeltype>0) then
                       !resulting Viscosity intentionally overwritten
                       call ComputeViscsWM(U,V,W,Pr,Temperature)
     end if

     if (sgstype==SubgridModel) then
                       call SGS_Smag(U,V,W,filter_ratios(filtertype))
     else if (sgstype==SigmaModel) then
                       call SGS_Sigma(U,V,W,filter_ratios(filtertype))
     else if (sgstype==VremanModel) then
                       call SGS_Vreman(U,V,W,filter_ratios(filtertype))
     else if (sgstype==StabSubgridModel) then
                       call SGS_StabSmag(U,V,W,Temperature,filter_ratios(filtertype))
     else
         if (Re>0) then
           Viscosity = 1._knd/Re
         else
           Viscosity = 0
         end if
     end if


     if (debuglevel>0) then
      write(*,*) "NUt", sum(Viscosity(1:Prnx,1:Prny,1:Prnz))/(Prnx*Prny*Prnz)
      write(*,*) "maxNUt", MAXVAL(Viscosity(1:Prnx,1:Prny,1:Prnz))
      write(*,*) "minNUt", MINVAL(Viscosity(1:Prnx,1:Prny,1:Prnz))
     end if

     if (wallmodeltype>0) then
                     call ComputeUVWFluxesWM(U,V,W,Pr,Temperature)
     end if

     call BoundViscosity(Viscosity)

     do i=1,size(ScalFlIBPoints)
       Viscosity(ScalFlIBPoints(i)%xi,ScalFlIBPoints(i)%yj,ScalFlIBPoints(i)%zk) =  &
                                               TIBPoint_Viscosity(ScalFlIBPoints(i),Viscosity)
     end do

     if (sgstype/=StabSubgridModel.and.enable_buoyancy)  call ComputeTDiff(U,V,W)

     if (size(TDiff)>0) call BoundViscosity(TDiff)
#endif

  end subroutine SubgridStresses

























  subroutine IBMomentum(U,V,W)
    use ImmersedBoundary, only: Up => UIBPoints, &
                                Vp => VIBPoints, &
                                Wp => WIBPoints, &
                                Interpolate => TIBPoint_Interpolate

    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U,V,W
    integer   :: i

    if (size(Up)+size(Vp)+size(Wp)>0) then
      if (.not. allocated(Uwm)) allocate(Uwm(size(Up)))
      if (.not. allocated(Vwm)) allocate(Vwm(size(Vp)))
      if (.not. allocated(Wwm)) allocate(Wwm(size(Wp)))

      !$omp parallel
      !$omp do
      do i=1,ubound(Up,1)
          Uwm(i) = Interpolate(Up(i), U, -2)
      end do
      !$omp end do
      !$omp do
      do i=1,ubound(Up,1)
          U(Up(i)%xi, Up(i)%yj, Up(i)%zk) = Uwm(i)
      end do
      !$omp end do nowait

      !$omp do
      do i=1,ubound(Vp,1)
          Vwm(i) = Interpolate(Vp(i), V, -2)
      end do
      !$omp end do
      !$omp do
      do i=1,ubound(Vp,1)
          V(Vp(i)%xi, Vp(i)%yj, Vp(i)%zk) = Vwm(i)
      end do
      !$omp end do nowait

      !$omp do
      do i=1,ubound(Wp,1)
          Wwm(i) = Interpolate(Wp(i), W, -2)
      end do
      !$omp end do
      !$omp do
      do i=1,ubound(Wp,1)
          W(Wp(i)%xi, Wp(i)%yj, Wp(i)%zk) = Wwm(i)
      end do
      !$omp end do
      !$omp end parallel

    end if
  end subroutine IBMomentum






  subroutine IBMassSources(Q,U,V,W)
    use ImmersedBoundary, only: Up => UIBPoints, &
                                Vp => VIBPoints, &
                                Wp => WIBPoints
    use vtkarray

    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U,V,W
    real(knd), intent(out) :: Q(0:,0:,0:)
    integer i,xi,yj,zk

    Q = 0

    do i=1,ubound(Up,1)
      xi = Up(i)%xi
      yj = Up(i)%yj
      zk = Up(i)%zk

      Q(xi,yj,zk)   = Q(xi,yj,zk)   + U(xi,yj,zk)/dxPr(xi)
      Q(xi+1,yj,zk) = Q(xi+1,yj,zk) - U(xi,yj,zk)/dxPr(xi+1)
    end do


    do i=1,ubound(Vp,1)
      xi = Vp(i)%xi
      yj = Vp(i)%yj
      zk = Vp(i)%zk

      Q(xi,yj,zk)   = Q(xi,yj,zk)   + V(xi,yj,zk)/dyPr(yj)
      Q(xi,yj+1,zk) = Q(xi,yj+1,zk) - V(xi,yj,zk)/dyPr(yj+1)
    end do


    do i=1,ubound(Wp,1)
      xi = Wp(i)%xi
      yj = Wp(i)%yj
      zk = Wp(i)%zk

      Q(xi,yj,zk)   = Q(xi,yj,zk)   + W(xi,yj,zk)/dzPr(zk)
      Q(xi,yj,zk+1) = Q(xi,yj,zk+1) - W(xi,yj,zk)/dzPr(zk+1)
    end do

!     call VtkArrayBin(output_dir//"Q1.vtk",Q,offsets_to_global)!(1:Prnx,1:Prny,1:Prnz)

    call Bound_Q(Q)

!     call VtkArrayBin(output_dir//"Q2.vtk",Q,offsets_to_global)!(1:Prnx,1:Prny,1:Prnz)

  end subroutine IBMassSources




#ifdef __HMPP
  subroutine GetDataFromGPU(getU,getV,getW,getPr,getTemperature,getVisc,getTDiff,U,V,W,Pr,Temperature)
    use HMPP_SCALARS, only:GetTemperatureFromGPU
    logical, intent(in) :: getU,getV,getW,getPr,getTemperature
    real(knd), intent(inout) :: U(-2:,-2:,-2:)
    real(knd), intent(inout) :: V(-2:,-2:,-2:)
    real(knd), intent(inout) :: W(-2:,-2:,-2:)
    real(knd), intent(inout) :: Pr(1:,1:,1:)
    real(knd), intent(inout) :: Temperature(-1:,-1:,-1:)


    if (getU) call GetUFromGPU(U)

    if (getV) call GetVFromGPU(V)

    if (getW) call GetWFromGPU(W)

    if (getPr) call GetPrFromGPU(Pr)

    if (getTemperature) call GetTemperatureFromGPU(Temperature)

    if (getViscosity) call GetViscFromGPU

    if (getTDiff) call GetTDiffFromGPU

  end subroutine GetDataFromGPU


  subroutine GetUFromGPU(U)
    real(knd), intent(inout) :: U(-2:,-2:,-2:)
      !$hmpp <tsteps> delegatedstore, args[UpdateU::U]
  end subroutine GetUFromGPU

  subroutine GetVFromGPU(V)
    real(knd), intent(inout) :: V(-2:,-2:,-2:)
      !$hmpp <tsteps> delegatedstore, args[UpdateU::V]
  end subroutine GetVFromGPU

  subroutine GetWFromGPU(W)
    real(knd), intent(inout) :: W(-2:,-2:,-2:)
      !$hmpp <tsteps> delegatedstore, args[UpdateU::W]
  end subroutine GetWFromGPU

  subroutine GetViscFromGPU
      !$hmpp <tsteps> delegatedstore, args[BoundViscosity_TDiff::Visc]
  end subroutine GetViscFromGPU


  subroutine GetTDiffFromGPU
      !$hmpp <tsteps> delegatedstore, args[BoundViscosity_TDiff::TDiff]
  end subroutine GetTDiffFromGPU
#endif


end module TSTEPS











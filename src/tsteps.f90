module TSTEPS

  use PARAMETERS
  use ArrayUtilities
  use CDS!, only: CDU, CDV, CDW, CDS4U, CDS4V, CDS4W, CDS4_2U, CDS4_2V, CDS4_2W
  use LAXFRIED
  use LAXWend
  use BOUNDARIES, only: BoundU,Bound_Q
  use POISSON !it exports Pr_Correct and GetPrFromGPU
  use OUTPUTS, only: store,display,proftempfl
  use SCALARS, only: ScalarRK3, Bound_Visc, ComputeTDiff, BoundTemperature
  use Subgrid, only: sgstype, SGS_Smag, SGS_StabSmag, SGS_Vreman, SGS_Sigma
  use TURBINLET, only: GetTurbInlet, GetInletFromFile
  use Wallmodels
  use Tiling, only: tilenx, tileny, tilenz
#ifdef __HMPP
  use HMPP_CODELETS
#endif

  implicit none


  private
  public TMarchRK3

#ifdef __HMPP
  public GetDataFromGPU
  type(hmppWMPoint),allocatable,save :: hmppWMPoints(:)
  integer,save :: nWMPoints = 0
#endif

contains

#ifndef __HMPP
#include "tsteps-shared.f90"
#endif


 subroutine TMarchRK3(U,V,W,Pr,Temperature,Scalar,delta)
  use RK3
  real(KND),allocatable,intent(inout):: U(:,:,:),V(:,:,:),W(:,:,:),Pr(:,:,:)  !allocatable to anable move_alloc
  real(KND),allocatable,intent(inout) :: Temperature(:,:,:),Scalar(:,:,:,:)
  real(KND),intent(out):: delta

  real(KND),dimension(:,:,:),allocatable,save   :: Q
  real(KND),dimension(:,:,:),allocatable,save   :: U2,Ustar
  real(KND),dimension(:,:,:),allocatable,save   :: V2,Vstar
  real(KND),dimension(:,:,:),allocatable,save   :: W2,Wstar

  integer l
  integer,save:: called = 0
  integer(int64), save :: trate
  integer(int64), save :: time1, time2


  if (called==0) then
   called = 1

   allocate(U2(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),lbound(U,3):ubound(U,3)))
   allocate(V2(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),lbound(V,3):ubound(V,3)))
   allocate(W2(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),lbound(W,3):ubound(W,3)))

   allocate(Ustar(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),lbound(U,3):ubound(U,3)))
   allocate(Vstar(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),lbound(V,3):ubound(V,3)))
   allocate(Wstar(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),lbound(W,3):ubound(W,3)))


   !$omp parallel
   !$omp sections
   !$omp section
   call BoundU(1,U)
   !$omp section
   call BoundU(2,V)
   !$omp section
   call BoundU(3,W)
   !$omp section
   if (buoyancy==1) call BoundTemperature(temperature)
   !$omp end sections
   !$omp end parallel


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
    !$hmpp <tsteps> advancedload, args[BoundU::nx,BoundU::ny,BoundU::nz,&
    !$hmpp &     BoundV::nx,BoundV::ny,BoundV::nz,BoundW::nx,BoundW::ny,BoundW::nz]

    !$hmpp <tsteps> advancedload, args[BoundU::sideU,BoundU::Uin,BoundV::Uin,BoundW::Uin]

    !$hmpp <tsteps> advancedload, args[BoundU::component,BoundV::component,BoundW::component]
    !$hmpp <tsteps> advancedload, args[BoundU::regime,BoundV::regime,BoundW::regime]

    !$hmpp <tsteps> advancedload, args[Vreman::dx,Vreman::dy,Vreman::dz,Vreman::Re]
    !$hmpp <tsteps> advancedload, args[Convection::buoyancy,Convection::convmet,Convection::coriolisparam]
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

    if (buoyancy==1) then
     !$hmpp <tsteps> advancedload, args[Convection::temperature]
    end if

    !$hmpp <tsteps> advancedload, args[ComputeViscWM::nWMPoints,ComputeViscWM::WMPoints]
    !$hmpp <tsteps> advancedload, args[Bound_Visc_TDiff::Prandtl]

    GPU = 1
#endif
  end if

  if ((Btype(We) ==TurbulentInlet).or.(Btype(Ea) ==TurbulentInlet)) then
    call GetTurbInlet
    !$hmpp <tsteps> advancedload, args[BoundU::sideU,BoundU::Uin,BoundV::Uin,BoundW::Uin]
  else if (Btype(We) ==InletFromFile) then
    call GetInletFromFile(time)
    !$hmpp <tsteps> advancedload, args[BoundU::sideU,BoundU::Uin,BoundV::Uin,BoundW::Uin]
  end if

  !$hmpp <tsteps> advancedload, args[TimeStepEul::time]
  !$hmpp <tsteps> TimeStepEul callsite, args[*].noupdate = true
  call TimeStepEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
               dxmin,dxU,dyV,dzW,CFL,Uref,steady,time,endtime,&
               U,V,W,dt)
  !$hmpp <tsteps> delegatedstore, args[TimeStepEul::dt]


  write (*,*) "time:",time,"dt: ",dt

  do l=1,RKstages

    write(*,*) "stage:",l


    if (debugparam>1) call system_clock(count=time1)



    call SubgridStresses(U,V,W,Pr,Temperature)


    !$hmpp <tsteps> advancedload, args[Convection::lev]


    !$hmpp <tsteps> Convection callsite, args[*].noupdate = true
    call Convection(U,V,W,U2,V2,W2,Ustar,Vstar,Wstar,temperature,RK_beta,RK_rho,l)


    call OtherTerms(U,V,W,U2,V2,W2,Pr,2._KND*RK_alpha(l))


    !download U2 V2 and W2 to main memory Attenuates below are to slow on GPU, waiting for HMPP update
!     !$hmpp <tsteps> delegatedstore, args[UnifRedblack::U2,UnifRedblack::V2,UnifRedblack::W2]

!     !this will havr no effect when doing HMPP
!     if ((Btype(To) ==FreeSlipBuff) .and. (Prnz>15))  then
!         !$hmpp <tsteps> AttenuateTop callsite, args[*].noupdate = true
!         call AttenuateTop(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Btype,&
!                           zPr,zW,U2,V2,W2,temperature,buoyancy)
!     end if
!
!     if ((Btype(Ea) ==OutletBuff) .and. (Prnx>15)) then
!         !$hmpp <tsteps> AttenuateOut callsite, args[*].noupdate = true
!         call AttenuateOut(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
!                           xPr,xU,U2,V2,W2,temperature,buoyancy)
!     end if


!     if ( ((Btype(To) ==FreeSlipBuff) .and.&
!          (Prnz>15)) .or. ((Btype(Ea) ==OutletBuff) .and. (Prnx>15)) ) then
!
!              call BoundU(1,U2)
!              call BoundU(2,V2)
!              call BoundU(3,W2)
!              !Not for temperature, it will be taken care of in ScalarRK3
!     end if

    if (masssourc==1) then
        call IBMassSources(Q,U2,V2,W2)
    end if


    if (poissmet>0) then
       call Pr_Correct(U2,V2,W2,Pr,Q,2._KND*RK_alpha(l))
    end if


#ifdef __HMPP
    !$hmpp <tsteps> UpdateU callsite, args[*].noupdate=true
    call UpdateU_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,U,V,W,U2,V2,W2)
#else
    if (l==1) delta = 0
#ifdef DEBUG
    if (debuglevel>0.or.steady==1) then
      if (Unx*Uny*Unz>0) then &
        delta = delta+sum(abs(U(1:Unx,1:Uny,1:Unz)-U2(1:Unx,1:Uny,1:Unz)))/(Unx*Uny*Unz)
      if (Vnx*Vny*Vnz>0) then &
        delta = delta+sum(abs(V(1:Vnx,1:Vny,1:Vnz)-V2(1:Vnx,1:Vny,1:Vnz)))/(Vnx*Vny*Vnz)
      if (Wnx*Wny*Wnz>0) then &
        delta = delta+sum(abs(W(1:Wnx,1:Wny,1:Wnz)-W2(1:Wnx,1:Wny,1:Wnz)))/(Wnx*Wny*Wnz)
    end if
#endif
    call exchange_alloc(U,U2)
    call exchange_alloc(V,V2)
    call exchange_alloc(W,W2)
#endif



    !Upload the new values from Pr_Correct to GPU


    call ScalarRK3(U,V,W,Temperature,Scalar,l,proftempfl)




    if ((Btype(Ea) ==OutletBuff) .and. (Prnx>15)) then
     !$hmpp <tsteps> AttenuateOut2 callsite, args[*].noupdate = true
      call AttenuateOut(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                        xPr,xU,U,V,W,temperature,buoyancy)
!       if (buoyancy==1) then
! #ifdef __HMPP
! !              !$hmpp <tsteps> BoundTemperature callsite
!              call BoundTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,TBtype,sideTemp,BsideTArr,BsideTFLArr,TDiff,TempIn,temperature)
! #else
!              call BoundTemperature(temperature)
! #endif
!       end if
    end if


     !$hmpp <tsteps> NullInterior callsite, args[*].noupdate = true
     call NullInterior(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                           nUnull,nVnull,nWnull,Unull,Vnull,Wnull,U,V,W)


#ifdef __HMPP
    !$hmpp <tsteps> BoundU callsite, args[*].noupdate = true
    call BoundU_GPU(1,Unx,Uny,Unz,Prny,Prnz,&
                         Btype,sideU,&
                         Uin,U,0)

    !$hmpp <tsteps> BoundV callsite, args[*].noupdate = true
    call BoundU_GPU(2,Vnx,Vny,Vnz,Prny,Prnz,&
                         Btype,sideU,&
                         Vin,V,0)

    !$hmpp <tsteps> BoundW callsite, args[*].noupdate = true
    call BoundU_GPU(3,Wnx,Wny,Wnz,Prny,Prnz,&
                         Btype,sideU,&
                         Win,W,0)
#else
!$omp parallel
!$omp sections
!$omp section
        call BoundU(1,U)
!$omp section
        call BoundU(2,V)
!$omp section
        call BoundU(3,W)
!$omp end sections
!$omp end parallel
#endif


    if (debugparam>1) then
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
   real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
   real(KND),intent(inout):: Pr(1:,1:,1:)
   real(KND),allocatable,intent(inout):: U2(:,:,:),V2(:,:,:),W2(:,:,:)
   real(KND),intent(in):: coef

   real(KND),dimension(:,:,:),allocatable,save:: U3,V3,W3
   real(KND) S

   integer i,j,k
   integer,save:: called=0

   if (called==0) then
     allocate(U3(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),lbound(U,3):ubound(U,3)))
     allocate(V3(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),lbound(V,3):ubound(V,3)))
     allocate(W3(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),lbound(W,3):ubound(W,3)))
     called=1
   end if

   write(*,*) "Otherterms:"

   !$hmpp <tsteps> advancedload, args[PressureGrad::coef]

   !Pressure gradient terms
   !$hmpp <tsteps> PressureGrad callsite, args[*].noupdate = true
   call PressureGrad(Prnx,Prny,Prnz,&
                     Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                     dxU,dyV,dzW,&
                     Btype,prgradientx,prgradienty,&
                     Pr,U2,V2,W2,&
                     dt,coef)


   Re_gt_0: if (Re>0) then

     !Diffusion using Crank Nicolson
     !first approximation using forward Euler
     !iteration SOR or Gauss-Seidel


!      !$hmpp <tsteps> advancedload, args[ForwEul::Visc]

     !$hmpp <tsteps> ForwEul callsite, args[*].noupdate = true
     call ForwEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                  dxPr,dyPr,dzPr,dxU,dyV,dzW,&
                  U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                  dt,coef)


     call IBMomentum(U2,V2,W2,U3,V3,W3)

     !Performs the diffusion terms
#ifdef __HMPP
!       !$hmpp <tsteps>  advancedload, args[UnifRedBlack::maxCNiter,UnifRedBlack::epsCN]

      !$hmpp <tsteps> UNIFREDBLACK callsite, args[*].noupdate=true
      call UNIFREDBLACK_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,&
                            Btype,sideU,&
                            dt,dxmin,dymin,dzmin,&
                            Uin,Vin,Win,&
                            U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                            coef,maxCNiter,epsCN,it,S)

     !$hmpp <tsteps> delegatedstore, args[UnifRedBlack::iters,UnifRedBlack::residuum]

      write(*,*) "back from GPU CN", it,S

#else
     if (gridtype==UNIFORMGRID) then
      call UNIFREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
     else
      call GENREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
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

   if (debuglevel>=2) then  !Compute and output the mean friction in the domain.
    S = 0
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
       S = S-((Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))/dxPr(i+1)-&
       Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k))/dxPr(i))/dxU(i)+&
         (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))/dyV(j)-&
         0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k))/dyV(j-1))/dyPr(j)+&
          (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)-&
         0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1))/dzW(k-1))/dzPr(k))
      end do
     end do
    end do

    S = S/(Unx*Uny*Unz)
    write(*,*) "Mean friction:", S
   end if


  end subroutine OtherTerms




  subroutine UNIFREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND),dimension(-2:,-2:,-2:),intent(inout):: U2,V2,W2,U3,V3,W3
   real(KND),intent(in):: coef
   real(KND) recdxmin2,recdymin2,recdzmin2                                                               !reciprocal values of dx**2
   real(KND) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
   integer i,j,k,bi,bj,bk,l
   integer,parameter :: narr = 3, narr2 = 5 !number of arrays in the loop
   real(KND),allocatable,save:: Apu(:,:,:),ApV(:,:,:),ApW(:,:,:)
   integer,save :: called = 0

       if (called==0) then
         allocate(Apu(1:Unx,1:Uny,1:Unz))
         allocate(ApV(1:Vnx,1:Vny,1:Vnz))
         allocate(ApW(1:Wnx,1:Wny,1:Wnz))
         called = 1
       end if


       Ap = coef*dt/(2._KND)
       S = 0
       l = 0

       recdxmin2=1./dxmin**2
       recdymin2=1./dymin**2
       recdzmin2=1./dzmin**2

       !$omp parallel private(i,j,k,bi,bj,bk)

       !The explicit part, which doesn't have to be changed inside the loop
       !$omp do schedule(runtime)
       do bk = 1,Unz,tilenz(narr)
        do bj = 1,Uny,tileny(narr)
         do bi = 1,Unx,tilenx(narr)
          do k = bk,min(bk+tilenz(narr)-1,Unz)
           do j = bj,min(bj+tileny(narr)-1,Uny)
            do i = bi,min(bi+tilenx(narr)-1,Unx)
              U2(i,j,k) = U2(i,j,k)+Ap*(&
              ((Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))-&
              Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k)))*recdxmin2+0.25_KND*(&
               ((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))-&
               (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k)))*recdymin2+&
               ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))-&
               (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1)))*recdzmin2)))
            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end do nowait
       !$omp do schedule(runtime)
       do bk = 1,Vnz,tilenz(narr)
        do bj = 1,Vny,tileny(narr)
         do bi = 1,Vnx,tilenx(narr)
          do k = bk,min(bk+tilenz(narr)-1,Vnz)
           do j = bj,min(bj+tileny(narr)-1,Vny)
            do i = bi,min(bi+tilenx(narr)-1,Vnx)
              V2(i,j,k) = V2(i,j,k)+Ap*(&
            (0.25_KND*((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V(i+1,j,k)-V(i,j,k))-&
            (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(V(i,j,k)-V(i-1,j,k)))*recdxmin2+&
             (Visc(i,j+1,k)*(V(i,j+1,k)-V(i,j,k))-&
             Visc(i,j,k)*(V(i,j,k)-V(i,j-1,k)))*recdymin2+&
             0.25_KND*((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))-&
             (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(V(i,j,k)-V(i,j,k-1)))*recdzmin2))
            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end do nowait
       !$omp do schedule(runtime)
       do bk = 1,Wnz,tilenz(narr)
        do bj = 1,Wny,tileny(narr)
         do bi = 1,Wnx,tilenx(narr)
          do k = bk,min(bk+tilenz(narr)-1,Wnz)
           do j = bj,min(bj+tileny(narr)-1,Wny)
            do i = bi,min(bi+tilenx(narr)-1,Wnx)
              W2(i,j,k) = W2(i,j,k)+Ap*(&
            (0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W(i+1,j,k)-W(i,j,k))-&
            (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(W(i,j,k)-W(i-1,j,k)))*recdxmin2+&
             ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W(i,j+1,k)-W(i,j,k))-&
             (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(W(i,j,k)-W(i,j-1,k)))*recdymin2)+&
             (Visc(i,j,k+1)*(W(i,j,k+1)-W(i,j,k))-&
             Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1)))*recdzmin2))
            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end do nowait

       !Auxiliary coefficients to better efficiency in loops
       !$omp do schedule(runtime)
       do bk = 1,Unz,tilenz(narr)
        do bj = 1,Uny,tileny(narr)
         do bi = 1,Unx,tilenx(narr)
          do k = bk,min(bk+tilenz(narr)-1,Unz)
           do j = bj,min(bj+tileny(narr)-1,Uny)
            do i = bi,min(bi+tilenx(narr)-1,Unx)
              ApU(i,j,k) = ((Visc(i+1,j,k)+&
                          Visc(i,j,k))*recdxmin2+&
                          0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))+&
                          (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k)))*recdymin2+&
                          ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                          (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1)))*recdzmin2))
            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end do
       !$omp end parallel

!        ApU = 1._KND/(1._KND+Ap*ApU)
       call multiply_and_add_scalar(ApU, Ap, 1._KND)
       call reciprocal(ApU, 1._KND)

       !$omp parallel private(i,j,k,bi,bj,bk)
       !$omp do schedule(runtime)
       do bk = 1,Vnz,tilenz(narr)
        do bj = 1,Vny,tileny(narr)
         do bi = 1,Vnx,tilenx(narr)
          do k = bk,min(bk+tilenz(narr)-1,Vnz)
           do j = bj,min(bj+tileny(narr)-1,Vny)
            do i = bi,min(bi+tilenx(narr)-1,Vnx)
              ApV(i,j,k) = ((Visc(i,j+1,k)+&
                         Visc(i,j,k))*recdymin2+&
                         0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))+&
                          (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k)))*recdxmin2+&
                         ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                         (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1)))*recdzmin2))
            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end do
       !$omp end parallel

!        ApV = 1._KND/(1._KND+Ap*ApV)
       call multiply_and_add_scalar(ApV, Ap, 1._KND)
       call reciprocal(ApV, 1._KND)

       !$omp parallel private(i,j,k,bi,bj,bk,Suavg,Svavg,Swavg)
       !$omp do schedule(runtime)
       do bk = 1,Wnz,tilenz(narr)
        do bj = 1,Wny,tileny(narr)
         do bi = 1,Wnx,tilenx(narr)
          do k = bk,min(bk+tilenz(narr)-1,Wnz)
           do j = bj,min(bj+tileny(narr)-1,Wny)
            do i = bi,min(bi+tilenx(narr)-1,Wnx)
              ApW(i,j,k) = (0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                          (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k)))*recdxmin2+&
                         ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))+&
                         (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k)))*recdymin2)+&
                         (Visc(i,j,k+1)+&
                         Visc(i,j,k))*recdzmin2)
            end do
           end do
          end do
         end do
        end do
       end do
       !$omp end do
       !$omp end parallel


!        ApW = 1._KND/(1._KND+Ap*ApW)
       call multiply_and_add_scalar(ApW, Ap, 1._KND)
       call reciprocal(ApW, 1._KND)



!        !$omp workshare   !problem with gfortran-4.7 - valgrind sees uninitialized values
!        Suavg = abs(MAXVAL(U3(1:Unx,1:Uny,1:Unz)))  !maximum values of velocities to norm the residues.
!        Svavg = abs(MAXVAL(V3(1:Vnx,1:Vny,1:Vnz)))
!        Swavg = abs(MAXVAL(W3(1:Wnx,1:Wny,1:Wnz)))
!        !$omp end workshare

       Suavg=0
       Svavg=0
       Swavg=0

       if (Suavg<=1e-3_KND) Suavg = 1
       if (Svavg<=1e-3_KND) Svavg = 1
       if (Swavg<=1e-3_KND) Swavg = 1


       do l=1,maxCNiter               !Gauss-Seidel iteration for Crank-Nicolson result
        call BoundU(1,U3)
        call BoundU(2,V3)
        call BoundU(3,W3)

        S = 0
        Su = 0
        Sv = 0
        Sw = 0
        !$omp parallel private(i,j,k,bi,bj,bk,p) reduction(max:Su,Sv,Sw)
        !$omp do schedule(runtime)
        do bk = 1,Unz,tilenz(narr2)
         do bj = 1,Uny,tileny(narr2)
          do bi = 1,Unx,tilenx(narr2)
           do k = bk,min(bk+tilenz(narr2)-1,Unz)
            do j = bj,min(bj+tileny(narr2)-1,Uny)
             do i = bi+mod(bi+j+k-1,2),min(bi+tilenx(narr2)-1,Unx),2
               p=((Visc(i+1,j,k)*(U3(i+1,j,k))-&
                Visc(i,j,k)*(-U3(i-1,j,k)))*recdxmin2+&
                0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U3(i,j+1,k))-&
                (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(-U3(i,j-1,k)))*recdymin2+&
                ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U3(i,j,k+1))-&
                (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-U3(i,j,k-1)))*recdzmin2))
               p = Ap*p+U2(i,j,k)+U(i,j,k)
               p = p*ApU(i,j,k)
               Su = max(Su,abs(p-U3(i,j,k)))
               U3(i,j,k) = p
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do schedule(runtime)
        do bk = 1,Vnz,tilenz(narr2)
         do bj = 1,Vny,tileny(narr2)
          do bi = 1,Vnx,tilenx(narr2)
           do k = bk,min(bk+tilenz(narr2)-1,Vnz)
            do j = bj,min(bj+tileny(narr2)-1,Vny)
             do i = bi+mod(bi+j+k-1,2),min(bi+tilenx(narr2)-1,Vnx),2
               p=((Visc(i,j+1,k)*(V3(i,j+1,k))-&
                Visc(i,j,k)*(-V3(i,j-1,k)))*recdymin2+&
                0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V3(i+1,j,k))-&
                (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(-V3(i-1,j,k)))*recdxmin2+&
                ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V3(i,j,k+1))-&
                (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-V3(i,j,k-1)))*recdzmin2))
               p = Ap*p+V2(i,j,k)+V(i,j,k)
               p = p*ApV(i,j,k)
               Sv = max(Sv,abs(p-V3(i,j,k)))
               V3(i,j,k) = p
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do schedule(runtime)
        do bk = 1,Wnz,tilenz(narr2)
         do bj = 1,Wny,tileny(narr2)
          do bi = 1,Wnx,tilenx(narr2)
           do k = bk,min(bk+tilenz(narr2)-1,Wnz)
            do j = bj,min(bj+tileny(narr2)-1,Wny)
             do i = bi+mod(bi+j+k-1,2),min(bi+tilenx(narr2)-1,Wnx),2
               p=(0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W3(i+1,j,k))-&
                (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(-W3(i-1,j,k)))*recdxmin2+&
                ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W3(i,j+1,k))-&
                (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(-W3(i,j-1,k)))*recdymin2)+&
                (Visc(i,j,k+1)*(W3(i,j,k+1))-&
                Visc(i,j,k)*(-W3(i,j,k-1)))*recdzmin2)
               p = Ap*p+W2(i,j,k)+W(i,j,k)
               p = p*ApW(i,j,k)
               Sw = max(Sw,abs(p-W3(i,j,k)))
               W3(i,j,k) = p
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp enddo

        !$omp do schedule(runtime)
        do bk = 1,Unz,tilenz(narr2)
         do bj = 1,Uny,tileny(narr2)
          do bi = 1,Unx,tilenx(narr2)
           do k = bk,min(bk+tilenz(narr2)-1,Unz)
            do j = bj,min(bj+tileny(narr2)-1,Uny)
             do i = bi+mod(bi+j+k,2),min(bi+tilenx(narr2)-1,Unx),2
               p=((Visc(i+1,j,k)*(U3(i+1,j,k))-&
                 Visc(i,j,k)*(-U3(i-1,j,k)))*recdxmin2+&
                 0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U3(i,j+1,k))-&
                 (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(-U3(i,j-1,k)))*recdymin2+&
                 ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U3(i,j,k+1))-&
                 (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-U3(i,j,k-1)))*recdzmin2))
                p = Ap*p+U2(i,j,k)+U(i,j,k)
                p = p*ApU(i,j,k)


                Su = max(Su,abs(p-U3(i,j,k)))
                U3(i,j,k) = p
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do schedule(runtime)
        do bk = 1,Vnz,tilenz(narr2)
         do bj = 1,Vny,tileny(narr2)
          do bi = 1,Vnx,tilenx(narr2)
           do k = bk,min(bk+tilenz(narr2)-1,Vnz)
            do j = bj,min(bj+tileny(narr2)-1,Vny)
             do i = bi+mod(bi+j+k-1,2),min(bi+tilenx(narr2)-1,Vnx),2
               p=((Visc(i,j+1,k)*(V3(i,j+1,k))-&
                Visc(i,j,k)*(-V3(i,j-1,k)))*recdymin2+&
                0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V3(i+1,j,k))-&
                (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(-V3(i-1,j,k)))*recdxmin2+&
                ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V3(i,j,k+1))-&
                (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-V3(i,j,k-1)))*recdzmin2))
               p = Ap*p+V2(i,j,k)+V(i,j,k)
               p = p*ApV(i,j,k)
               Sv = max(Sv,abs(p-V3(i,j,k)))
               V3(i,j,k) = p
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do schedule(runtime)
        do bk = 1,Wnz,tilenz(narr2)
         do bj = 1,Wny,tileny(narr2)
          do bi = 1,Wnx,tilenx(narr2)
           do k = bk,min(bk+tilenz(narr2)-1,Wnz)
            do j = bj,min(bj+tileny(narr2)-1,Wny)
             do i = bi+mod(bi+j+k-1,2),min(bi+tilenx(narr2)-1,Wnx),2
               p=(0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W3(i+1,j,k))-&
                (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(-W3(i-1,j,k)))*recdxmin2+&
                ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W3(i,j+1,k))-&
                (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(-W3(i,j-1,k)))*recdymin2)+&
                (Visc(i,j,k+1)*(W3(i,j,k+1))-&
                Visc(i,j,k)*(-W3(i,j,k-1)))*recdzmin2)
               p = Ap*p+W2(i,j,k)+W(i,j,k)
               p = p*ApW(i,j,k)
               Sw = max(Sw,abs(p-W3(i,j,k)))
               W3(i,j,k) = p
             end do
            end do
           end do
          end do
         end do
        end do
        !$omp enddo
        !$omp end parallel
        S = max(Su/Suavg,Sv/Svavg,Sw/Swavg)
        write (*,*) "CN ",l,S

        if (S<=epsCN) exit
       end do

  end subroutine UNIFREDBLACK



  subroutine GENREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND),dimension(-2:,-2:,-2:),intent(inout):: U2,V2,W2,U3,V3,W3
   real(KND),intent(in):: coef
   real(KND) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
   integer i,j,k,l
   real(KND),allocatable,save:: Apu(:,:,:),ApV(:,:,:),ApW(:,:,:)
   integer,save :: called = 0

       if (called==0) then
         allocate(Apu(1:Unx,1:Uny,1:Unz))
         allocate(ApV(1:Vnx,1:Vny,1:Vnz))
         allocate(ApW(1:Wnx,1:Wny,1:Wnz))
         called = 1
       end if



       Ap = coef*dt/(2._KND)
       S = 0
       l = 0


       do k=1,Unz    !The explicit part, which doesn't have to be changed inside the loop
        do j=1,Uny
         do i=1,Unx
          U2(i,j,k) = U2(i,j,k)+Ap*(&
          ((Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))/dxPr(i+1)-&
          Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k))/dxPr(i))/dxU(i)+&
           (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))/dyV(j)-&
           0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k))/dyV(j-1))/dyPr(j)+&
           (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)-&
           0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1))/dzW(k-1))/dzPr(k)))
         end do
        end do
       end do
       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          V2(i,j,k) = V2(i,j,k)+Ap*(&
          ((0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V(i+1,j,k)-V(i,j,k))/dxU(i)-&
          0.25_KND*(Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(V(i,j,k)-V(i-1,j,k))/dxU(i-1))/dxPr(i)+&
           (Visc(i,j+1,k)*(V(i,j+1,k)-V(i,j,k))/dyPr(j+1)-&
           Visc(i,j,k)*(V(i,j,k)-V(i,j-1,k))/dyPr(j))/dyV(j)+&
           (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))/dzW(k)-&
           0.25_KND*(Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(V(i,j,k)-V(i,j,k-1))/dzW(k-1))/dzPr(k)))
         end do
        end do
       end do
       do k=1,Wnz
       do j=1,Wny
        do i=1,Wnx
         W2(i,j,k) = W2(i,j,k)+Ap*(&
         ((0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W(i+1,j,k)-W(i,j,k))/dxU(i)-&
         0.25_KND*(Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(W(i,j,k)-W(i-1,j,k))/dxU(i-1))/dxPr(i)+&
          (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W(i,j+1,k)-W(i,j,k))/dyV(j)-&
          0.25_KND*(Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(W(i,j,k)-W(i,j-1,k))/dyV(j-1))/dyPr(j)+&
          (Visc(i,j,k+1)*(W(i,j,k+1)-W(i,j,k))/dzPr(k+1)-&
          Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1))/dzPr(k))/dzW(k)))
        end do
       end do
      end do

       do k=1,Unz         !Auxiliary coefficients to better efficiency in loops
        do j=1,Uny
         do i=1,Unx
          ApU(i,j,k) = ((Visc(i+1,j,k)/dxPr(i+1)+&
                      Visc(i,j,k)/dxPr(i))/dxU(i)+&
                      (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))/dyV(j)+&
                      0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))/dyV(j-1))/dyPr(j)+&
                      (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))/dzW(k)+&
                      0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))/dzW(k-1))/dzPr(k))
         end do
        end do
       end do

      ApU = 1._KND/(1._KND+Ap*ApU)


       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          ApV(i,j,k) =  ((0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))/dxU(i)+&
                      0.25_KND*(Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))/dxU(i-1))/dxPr(i)+&
                     (Visc(i,j+1,k)/dyPr(j+1)+&
                     Visc(i,j,k)/dyPr(j))/dyV(j)+&
                     (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))/dzW(k)+&
                     0.25_KND*(Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))/dzW(k-1))/dzPr(k))
         end do
        end do
       end do

       ApV = 1._KND/(1._KND+Ap*ApV)


       do k=1,Wnz
        do j=1,Wny
         do i=1,Wnx
          ApW(i,j,k) =  ((0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))/dxU(i)+&
                      0.25_KND*(Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))/dxU(i-1))/dxPr(i)+&
                     (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))/dyV(j)+&
                     0.25_KND*(Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))/dyV(j-1))/dyPr(j)+&
                     (Visc(i,j,k+1)/dzPr(k+1)+&
                     Visc(i,j,k)/dzPr(k))/dzW(k))
         end do
        end do
       end do

       ApW = 1._KND/(1._KND+Ap*ApW)

       Suavg = abs(MAXVAL(U3(1:Unx,1:Uny,1:Unz)))  !maximum values of velocities to norm the residues.
       Svavg = abs(MAXVAL(V3(1:Vnx,1:Vny,1:Vnz)))
       Swavg = abs(MAXVAL(W3(1:Wnx,1:Wny,1:Wnz)))
       if (Suavg<=1e-3_KND) Suavg = 1
       if (Svavg<=1e-3_KND) Svavg = 1
       if (Swavg<=1e-3_KND) Swavg = 1



       do l=1,maxCNiter               !Gauss-Seidel iteration for Crank-Nicolson result
        call BoundU(1,U3)
        call BoundU(2,V3)
        call BoundU(3,W3)
        S = 0
        Su = 0
        Sv = 0
        Sw = 0
        !$omp parallel private(i,j,k,p) reduction(max:Su,Sv,Sw)
        !$omp do
        do k=1,Unz
         do j=1,Uny
          do i=1+mod(j+k,2),Unx,2
            p=((Visc(i+1,j,k)*(U3(i+1,j,k))/dxPr(i+1)-&
             Visc(i,j,k)*(-U3(i-1,j,k))/dxPr(i))/dxU(i)+&
             (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U3(i,j+1,k))/dyV(j)-&
             0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(-U3(i,j-1,k))/dyV(j-1))/dyPr(j)+&
             (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U3(i,j,k+1))/dzW(k)-&
             0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-U3(i,j,k-1))/dzW(k-1))/dzPr(k))
            p = Ap*p+U2(i,j,k)+U(i,j,k)
            p = p*ApU(i,j,k)


            Su = max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k) = U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do
        do k=1,Vnz
         do j=1,Vny
          do i=1+mod(j+k,2),Vnx,2
            p=((0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V3(i+1,j,k))/dxU(i)-&
             0.25_KND*(Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(-V3(i-1,j,k))/dxU(i-1))/dxPr(i)+&
             (Visc(i,j+1,k)*(V3(i,j+1,k))/dyPr(j+1)-&
             Visc(i,j,k)*(-V3(i,j-1,k))/dyPr(j))/dyV(j)+&
             (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V3(i,j,k+1))/dzW(k)-&
             0.25_KND*(Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-V3(i,j,k-1))/dzW(k-1))/dzPr(k))
            p = Ap*p+V2(i,j,k)+V(i,j,k)
            p = p*ApV(i,j,k)
            Sv = max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k) = V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do
        do k=1,Wnz
         do j=1,Wny
          do i=1+mod(j+k,2),Wnx,2
            p=((0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W3(i+1,j,k))/dxU(i)-&
             0.25_KND*(Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(-W3(i-1,j,k))/dxU(i-1))/dxPr(i)+&
             (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W3(i,j+1,k))/dyV(j)-&
             0.25_KND*(Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(-W3(i,j-1,k))/dyV(j-1))/dyPr(j)+&
             (Visc(i,j,k+1)*(W3(i,j,k+1))/dzPr(k+1)-&
             Visc(i,j,k)*(-W3(i,j,k-1))/dzPr(k))/dzW(k))
            p = Ap*p+W2(i,j,k)+W(i,j,k)
            p = p*ApW(i,j,k)
            Sw = max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k) = W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
          end do
         end do
        end do
        !$omp enddo

        !$omp do
        do k=1,Unz
         do j=1,Uny
          do i=1+mod(j+k+1,2),Unx,2
            p=((Visc(i+1,j,k)*(U3(i+1,j,k))/dxPr(i+1)-&
             Visc(i,j,k)*(-U3(i-1,j,k))/dxPr(i))/dxU(i)+&
             (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U3(i,j+1,k))/dyV(j)-&
             0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(-U3(i,j-1,k))/dyV(j-1))/dyPr(j)+&
             (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U3(i,j,k+1))/dzW(k)-&
             0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-U3(i,j,k-1))/dzW(k-1))/dzPr(k))
            p = Ap*p+U2(i,j,k)+U(i,j,k)
            p = p*ApU(i,j,k)


            Su = max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k) = U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do
        do k=1,Vnz
         do j=1,Vny
          do i=1+mod(j+k+1,2),Vnx,2
            p=((0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V3(i+1,j,k))/dxU(i)-&
             0.25_KND*(Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(-V3(i-1,j,k))/dxU(i-1))/dxPr(i)+&
             (Visc(i,j+1,k)*(V3(i,j+1,k))/dyPr(j+1)-&
             Visc(i,j,k)*(-V3(i,j-1,k))/dyPr(j))/dyV(j)+&
             (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V3(i,j,k+1))/dzW(k)-&
             0.25_KND*(Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-V3(i,j,k-1))/dzW(k-1))/dzPr(k))
            p = Ap*p+V2(i,j,k)+V(i,j,k)
            p = p*ApV(i,j,k)
            Sv = max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k) = V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
          end do
         end do
        end do
        !$omp enddo nowait
        !$omp do
        do k=1,Wnz
         do j=1,Wny
          do i=1+mod(j+k+1,2),Wnx,2
            p=((0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W3(i+1,j,k))/dxU(i)-&
             0.25_KND*(Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(-W3(i-1,j,k))/dxU(i-1))/dxPr(i)+&
             (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W3(i,j+1,k))/dyV(j)-&
             0.25_KND*(Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(-W3(i,j-1,k))/dyV(j-1))/dyPr(j)+&
             (Visc(i,j,k+1)*(W3(i,j,k+1))/dzPr(k+1)-&
             Visc(i,j,k)*(-W3(i,j,k-1))/dzPr(k))/dzW(k))
            p = Ap*p+W2(i,j,k)+W(i,j,k)
            p = p*ApW(i,j,k)
            Sw = max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k) = W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
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

    real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(KND),intent(in):: Pr(1:,1:,1:)
    real(KND),intent(in):: Temperature(-1:,-1:,-1:)
    integer i

#ifdef __HMPP

     !$hmpp <tsteps> Vreman callsite, args[*].noupdate = true
     call Vreman_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxmin,dymin,dzmin,dt,Re,U,V,W,Visc)
!      !$hmpp <tsteps> delegatedstore, args[Vreman::Visc]

2._KND
     if (wallmodeltype>0) then

!        !$hmpp <tsteps> delegatedstore, args[AttenuateTop::Temperature],&
!        !$hmpp & args[AttenuateTop::Temperature].section={-1:Prnx+2,-1:Prny+2,1:1} of {-1:Prnx+2,-1:Prny+2,-1:Prnz+2}

       !$hmpp <tsteps> ComputeViscWM callsite, args[*].noupdate = true
       call ComputeViscsWM_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                                nWMPoints,hmppWMPoints,&
                                TBtype,Re,temperature_ref,grav_acc,&
                                U,V,W,Visc)
     end if


     !boundaries to Visc, compute TDiff and boundaries to TDiff

      !$hmpp <tsteps> Bound_Visc_TDiff callsite, args[*].noupdate = true
      call Bound_Visc_TDiff(Prnx,Prny,Prnz,buoyancy,Btype,Re,Prandtl,Visc,TDiff)



#else


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
           Visc=1._KND/Re
         else
           Visc=0
         end if
     end if


     if (debuglevel>0) then
      write(*,*) "NUt", sum(Visc(1:Prnx,1:Prny,1:Prnz))/(Prnx*Prny*Prnz)
      write(*,*) "maxNUt", MAXVAL(Visc(1:Prnx,1:Prny,1:Prnz))
      write(*,*) "minNUt", MINVAL(Visc(1:Prnx,1:Prny,1:Prnz))
     end if

     if (wallmodeltype>0) then
                     call ComputeViscsWM(U,V,W,Pr,Temperature)
     end if

     do i=1,size(ScalFlIBPoints)
       Visc(ScalFlIBPoints(i)%xi,ScalFlIBPoints(i)%yj,ScalFlIBPoints(i)%zk) = &
                                               TIBPoint_Viscosity(ScalFlIBPoints(i),Visc)
     end do

     call Bound_Visc(Visc)

     if (sgstype/=StabSubgridModel.and.buoyancy==1)  call ComputeTDiff(U,V,W)

     if (buoyancy==1) call Bound_Visc(TDiff)
#endif

  end subroutine SubgridStresses

























  subroutine IBMomentum(U2,V2,W2,U3,V3,W3)
    use ImmersedBoundary, only: UIBPoints, VIBPoints, WIBPoints, TIBPoint_MomentumSource

    real(KND),dimension(-2:,-2:,-2:),intent(inout):: U2,V2,W2,U3,V3,W3
    integer   :: i,xi,yj,zk
    real(KND) :: src

    if (size(UIBPoints)+size(UIBPoints)+size(UIBPoints)>0) then
      !$hmpp <tsteps> delegatedstore, args[ForwEul::U3,ForwEul::V3,ForwEul::W3]
      call BoundU(1,U3)
      call BoundU(2,V3)
      call BoundU(3,W3)

     !$omp parallel
      !$omp do private(xi,yj,zk,src)
      do i=1,ubound(UIBPoints,1)
        xi = UIBPoints(i)%xi
        yj = UIBPoints(i)%yj
        zk = UIBPoints(i)%zk

        src = TIBPoint_MomentumSource(UIBPoints(i),U3)

        U3(xi,yj,zk) = U3(xi,yj,zk) + src * dt
        U2(xi,yj,zk) = U2(xi,yj,zk) + src * dt
      end do
      !$omp end do nowait

      !$omp do private(xi,yj,zk,src)
      do i=1,ubound(VIBPoints,1)
        xi = VIBPoints(i)%xi
        yj = VIBPoints(i)%yj
        zk = VIBPoints(i)%zk

        src = TIBPoint_MomentumSource(VIBPoints(i),V3)

        V3(xi,yj,zk) = V3(xi,yj,zk) + src * dt
        V2(xi,yj,zk) = V2(xi,yj,zk) + src * dt
      end do
      !$omp end do nowait

      !$omp do private(xi,yj,zk,src)
      do i=1,ubound(WIBPoints,1)
        xi = WIBPoints(i)%xi
        yj = WIBPoints(i)%yj
        zk = WIBPoints(i)%zk

        src = TIBPoint_MomentumSource(WIBPoints(i),W3)

        W3(xi,yj,zk) = W3(xi,yj,zk) + src * dt
        W2(xi,yj,zk) = W2(xi,yj,zk) + src * dt
      end do
      !$omp end do
      !$omp end parallel

      !$hmpp <tsteps> advancedload, args[UnifRedBlack::U3,UnifRedBlack::V3,UnifRedBlack::W3]
    end if
  end subroutine IBMomentum






  subroutine IBMassSources(Q,U,V,W)
    use ImmersedBoundary, only: UIBPoints, VIBPoints, WIBPoints

    real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
    real(KND),intent(out):: Q(0:,0:,0:)
    integer i,xi,yj,zk

    Q = 0

    do i=1,ubound(UIBPoints,1)
      xi = UIBPoints(i)%xi
      yj = UIBPoints(i)%yj
      zk = UIBPoints(i)%zk

      Q(xi,yj,zk)   = Q(xi,yj,zk)   + U(xi,yj,zk)/dxPr(xi)
      Q(xi+1,yj,zk) = Q(xi+1,yj,zk) - U(xi,yj,zk)/dxPr(xi+1)
    end do


    do i=1,ubound(VIBPoints,1)
      xi = VIBPoints(i)%xi
      yj = VIBPoints(i)%yj
      zk = VIBPoints(i)%zk

      Q(xi,yj,zk)   = Q(xi,yj,zk)   + V(xi,yj,zk)/dyPr(yj)
      Q(xi,yj+1,zk) = Q(xi,yj+1,zk) - V(xi,yj,zk)/dyPr(yj+1)
    end do


    do i=1,ubound(WIBPoints,1)
      xi = WIBPoints(i)%xi
      yj = WIBPoints(i)%yj
      zk = WIBPoints(i)%zk

      Q(xi,yj,zk)   = Q(xi,yj,zk)   + W(xi,yj,zk)/dzPr(zk)
      Q(xi,yj,zk+1) = Q(xi,yj,zk+1) - W(xi,yj,zk)/dzPr(zk+1)
    end do

    call Bound_Q(Q)

  end subroutine IBMassSources




#ifdef __HMPP
  subroutine GetDataFromGPU(getU,getV,getW,getPr,getTemperature,getVisc,getTDiff,U,V,W,Pr,Temperature)
    use HMPP_SCALARS, only:GetTemperatureFromGPU
    logical,intent(in) :: getU,getV,getW,getPr,getTemperature
    real(KND),intent(inout) :: U(-2:,-2:,-2:)
    real(KND),intent(inout) :: V(-2:,-2:,-2:)
    real(KND),intent(inout) :: W(-2:,-2:,-2:)
    real(KND),intent(inout) :: Pr(1:,1:,1:)
    real(KND),intent(inout) :: Temperature(-1:,-1:,-1:)


    if (getU) call GetUFromGPU(U)

    if (getV) call GetVFromGPU(V)

    if (getW) call GetWFromGPU(W)

    if (getPr) call GetPrFromGPU(Pr)

    if (getTemperature) call GetTemperatureFromGPU(Temperature)

    if (getVisc) call GetViscFromGPU

    if (getTDiff) call GetTDiffFromGPU

  end subroutine GetDataFromGPU


  subroutine GetUFromGPU(U)
    real(KND),intent(inout) :: U(-2:,-2:,-2:)
      !$hmpp <tsteps> delegatedstore, args[UpdateU::U]
  end subroutine GetUFromGPU

  subroutine GetVFromGPU(V)
    real(KND),intent(inout) :: V(-2:,-2:,-2:)
      !$hmpp <tsteps> delegatedstore, args[UpdateU::V]
  end subroutine GetVFromGPU

  subroutine GetWFromGPU(W)
    real(KND),intent(inout) :: W(-2:,-2:,-2:)
      !$hmpp <tsteps> delegatedstore, args[UpdateU::W]
  end subroutine GetWFromGPU

  subroutine GetViscFromGPU
      !$hmpp <tsteps> delegatedstore, args[Bound_Visc_TDiff::Visc]
  end subroutine GetViscFromGPU


  subroutine GetTDiffFromGPU
      !$hmpp <tsteps> delegatedstore, args[Bound_Visc_TDiff::TDiff]
  end subroutine GetTDiffFromGPU
#endif


end module TSTEPS











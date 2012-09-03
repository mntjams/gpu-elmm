module TSTEPS

  use CDS, only: CDU, CDV, CDW, KappaU, KappaV, KappaW
  use LAXFRIED
  use LAXWEND
  use PARAMETERS
  use BOUNDARIES, only: BoundU,Bound_Q
  use POISSON, only: Pr_Correct
  use OUTPUTS, only: store,display,proftempfl
  use SCALARS, only: ScalarRK3, Bound_Visc
  use SMAGORINSKY, only: Smag, StabSmag, Vreman
  use TURBINLET, only: GetTurbInlet, GetInletFromFile
  use Wallmodels, only: ComputeViscsWM
  use Tiling, only: tilenx, tileny, tilenz
#ifdef __HMPP
  use HMPP_CODELETS
#endif

  implicit none


  private
  public TMarchRK3

contains

#ifndef __HMPP
#include "tsteps-shared.f90"
#endif


 subroutine TMarchRK3(U,V,W,Pr,Temperature,Scalar,delta)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),intent(inout) :: Temperature(-1:,-1:,-1:),Scalar(-1:,-1:,-1:,-1:)
  real(KND),intent(out):: delta

  real(KND),dimension(:,:,:),allocatable,save   :: Q
  real(KND),dimension(:,:,:),allocatable,save   ::U2,Ustar
  real(KND),dimension(:,:,:),allocatable,save   ::V2,Vstar
  real(KND),dimension(:,:,:),allocatable,save   ::W2,Wstar

  real(KND),dimension(1:3),parameter:: alpha = (/ 4._KND/15._KND, 1._KND/15._KND, 1._KND/6._KND /)
  real(KND),dimension(1:3),parameter:: beta  = (/ 8._KND/15._KND, 5._KND/12._KND, 3._KND/4._KND /)
  real(KND),dimension(1:3),parameter:: rho   = (/       0._KND, -17._KND/60._KND,-5._KND/12._KND/)

  integer i,l
  integer,save:: called = 0
  integer(DBL), save :: trate
  integer(DBL), save :: time1, time2


  if (called==0) then
   called = 1

   allocate(U2(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),lbound(U,3):ubound(U,3)))
   allocate(V2(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),lbound(V,3):ubound(V,3)))
   allocate(W2(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),lbound(W,3):ubound(W,3)))

   allocate(Ustar(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),lbound(U,3):ubound(U,3)))
   allocate(Vstar(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),lbound(V,3):ubound(V,3)))
   allocate(Wstar(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),lbound(W,3):ubound(W,3)))


   if (masssourc==1) allocate(Q(0:Prnx+1,0:Prny+1,0:Prnz+1))

   if (debugparam>1) call system_clock(count_rate=trate)

    !$hmpp <tsteps> allocate
    !$hmpp <tsteps> advancedload, args[Vreman::Prnx,BoundU::Prny,BoundU::Prnz]
    !$hmpp <tsteps> advancedload, args[BoundU::nx,BoundU::ny,BoundU::nz,&
    !$hmpp &     BoundV::nx,BoundV::ny,BoundV::nz,BoundW::nx,BoundW::ny,BoundW::nz]

    !$hmpp <tsteps> advancedload, args[BoundU::sideU,BoundU::Uin,BoundV::Uin,BoundW::Uin]

    !$hmpp <tsteps> advancedload, args[BoundU::component,BoundV::component,BoundW::component]
    !$hmpp <tsteps> advancedload, args[BoundU::regime,BoundV::regime,BoundW::regime]

    !$hmpp <tsteps> advancedload, args[Vreman::dx,Vreman::dy,Vreman::dz,Vreman::Re]
    !$hmpp <tsteps> advancedload, args[Convection::buoyancy,Convection::convmet,Convection::coriolisparam]
    !$hmpp <tsteps> advancedload, args[Convection::grav_acc,Convection::temperature_ref,Convection::beta,Convection::rho]
    !$hmpp <tsteps> advancedload, args[PressureGrad::prgradientx,PressureGrad::prgradienty]

    !$hmpp <tsteps> advancedload, args[BoundU::Btype]

    !$hmpp <tsteps> advancedload, args[ForwEul::dxPr,ForwEul::dyPr,ForwEul::dzPr]
    !$hmpp <tsteps> advancedload, args[AttenuateOut::xPr,AttenuateTop::zPr]

    !$hmpp <tsteps> advancedload, args[PressureGrad::dxU,PressureGrad::dyV,PressureGrad::dzW]
    !$hmpp <tsteps> advancedload, args[AttenuateOut::xU,AttenuateTop::zW]

    !$hmpp <tsteps> advancedload, args[PressureGrad::Pr]

    !$hmpp <tsteps> advancedload, args[NullInterior::Unull,NullInterior::Vnull,NullInterior::Wnull]
    !$hmpp <tsteps> advancedload, args[NullInterior::nUnull,NullInterior::nVnull,NullInterior::nWnull]

    !$hmpp <tsteps> advancedload, args[TimeStepEul::CFL,TimeStepEul::Uref,TimeStepEul::steady,TimeStepEul::dxmin,TimeStepEul::endtime]


    !$hmpp <tsteps> advancedload, args[TimeStepEul::U,TimeStepEul::V,TimeStepEul::W]

    if (buoyancy==1) then
     !$hmpp <tsteps> advancedload, args[Convection::temperature]
    endif
  endif



  if ((Btype(We) ==TurbulentInlet).or.(Btype(Ea) ==TurbulentInlet)) then
    call GetTurbInlet
    !$hmpp <tsteps> advancedload, args[BoundU::sideU,BoundU::Uin,BoundV::Uin,BoundW::Uin]
  else if (Btype(We) ==InletFromFile) then
    call GetInletFromFile(time)
    !$hmpp <tsteps> advancedload, args[BoundU::sideU,BoundU::Uin,BoundV::Uin,BoundW::Uin]
  endif

  !$hmpp <tsteps> advancedload, args[TimeStepEul::time]
  !$hmpp <tsteps> TimeStepEul callsite, args[*].noupdate = true
  call TimeStepEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
               dxmin,dxU,dyV,dzW,CFL,Uref,steady,time,endtime,&
               U,V,W,dt)
  !$hmpp <tsteps> delegatedstore, args[TimeStepEul::dt]

  write (*,*) "time:",time,"dt: ",dt

  do l=1,3

write(*,*) "stage:",l


    if (debugparam>1) call system_clock(count=time1)

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

    call SubgridStresses(U,V,W,Pr)


    !$hmpp <tsteps> advancedload, args[Convection::lev]


    !$hmpp <tsteps> Convection callsite, args[*].noupdate = true
    call Convection(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,buoyancy,convmet,&
                       dxmin,dymin,dzmin,coriolisparam,grav_acc,temperature_ref,&
                       U,V,W,U2,V2,W2,Ustar,Vstar,Wstar,temperature,beta,rho,l,dt)


    call OtherTerms(U,V,W,U2,V2,W2,Pr,2._KND*alpha(l))



    if ((Btype(To) ==FreeSlipBuff) .and. (Prnz>15))  then
        !$hmpp <tsteps> AttenuateTop callsite, args[*].noupdate = true
        call AttenuateTop(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Btype,&
                          zPr,zW,U2,V2,W2,temperature,buoyancy)
    endif

    if ((Btype(Ea) ==OutletBuff) .and. (Prnx>15)) then
        !$hmpp <tsteps> AttenuateOut callsite, args[*].noupdate = true
        call AttenuateOut(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                          xPr,xU,U2,V2,W2,temperature,buoyancy)
    endif

    !download U2 V2 and W2 to main memory
    !$hmpp <tsteps> delegatedstore, args[AttenuateOut::U,AttenuateOut::V,AttenuateOut::W]
    if (buoyancy==1) then
      !$hmpp <tsteps> delegatedstore, args[AttenuateOut::temperature]
    endif



    if (debugparam>1) then
     call system_clock(count=time2)
     write (*,*) "ET of part 1", (time2-time1)/real(trate)
     time1 = time2
    endif


    if (masssourc==1) then
        call BoundU(1,U2)
        call BoundU(2,V2)
        call BoundU(3,W2)
        call IBMassSources(Q,U2,V2,W2)
    endif




    if (poissmet>0) then
     if (masssourc==1) then
       call Pr_Correct(U2,V2,W2,Pr,2._KND*alpha(l),Q)
     else
       call Pr_Correct(U2,V2,W2,Pr,2._KND*alpha(l))
     endif
    endif


    if (l==1) delta = 0

    !$omp parallel
    if (debuglevel>0) then
      !$omp workshare
      delta = delta+sum(abs(U(1:Unx,1:Uny,1:Unz)-U2(1:Unx,1:Uny,1:Unz)))/(Unx*Uny*Unz)
      delta = delta+sum(abs(V(1:Vnx,1:Vny,1:Vnz)-V2(1:Vnx,1:Vny,1:Vnz)))/(Vnx*Vny*Vnz)
      delta = delta+sum(abs(W(1:Wnx,1:Wny,1:Wnz)-W2(1:Wnx,1:Wny,1:Wnz)))/(Wnx*Wny*Wnz)
      !$omp end workshare
    endif

    !$omp workshare
    U = U2
    V = V2
    W = W2
    !$omp end workshare
    !$omp end parallel



    !Upload the new values from Pr_Correct to GPU

    !$hmpp <tsteps> advancedload, args[PressureGrad::Pr]
    !$hmpp <tsteps> advancedload, args[AttenuateOut2::U,AttenuateOut2::V,AttenuateOut2::W]




   ! Visc should be in memory, as it is computed by CPU for now.

    if (store%BLprofiles>0.and.averaging==1) then
      call ScalarRK3(U,V,W,Temperature,Scalar,l,proftempfl)
    else
      call ScalarRK3(U,V,W,Temperature,Scalar,l)
    end if




    if (debugparam>1) then
     call system_clock(count=time2)
     write (*,*) "ET of part 2", (time2-time1)/real(trate)
     time1 = time2
    endif

    if (buoyancy==1) then
     !$hmpp <tsteps> advancedload, args[AttenuateOut2::temperature]
    endif

    if ((Btype(Ea) ==OutletBuff) .and. (Prnx>15)) then
     !$hmpp <tsteps> AttenuateOut2 callsite, args[*].noupdate = true
      call AttenuateOut(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                        xPr,xU,U,V,W,temperature,buoyancy)
    endif


    !$hmpp <tsteps> NullInterior callsite, args[*].noupdate = true
    call NullInterior(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                          nUnull,nVnull,nWnull,Unull,Vnull,Wnull,U,V,W)




    if (debugparam>1) then
     call system_clock(count=time2)
     write (*,*) "ET of part 3", (time2-time1)/real(trate)
     time1 = time2
    endif

   enddo
   if (.false.) then
   !$hmpp <tsteps> release
   endif
 end subroutine TMarchRK3








!  !$hmpp <tsteps> UpdateU codelet
!  subroutine UpdateU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,U,V,W,U2,V2,W2)
!  implicit none
! #ifdef __HMPP
! #include "hmpp-include.f90"
! #endif
!  integer,intent(in)   :: Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz
!  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(out)  :: U
!  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(in)   :: U2
!  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(out)  :: V
!  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(in)   :: V2
!  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(out)  :: W
!  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(in)   :: W2
!  integer i,j,k
!
!    !$hmppcg grid blocksize 512x1
!    !$hmppcg permute (k,i,j)
!    do k=-2,Unz+3
!     do j=-2,Uny+3
!      do i=-2,Unx+3
!       U(i,j,k) = U2(i,j,k)
!      enddo
!     enddo
!    enddo
!
!    !$hmppcg grid blocksize 512x1
!    !$hmppcg permute (k,i,j)
!    do k=-2,Vnz+3
!     do j=-2,Vny+3
!      do i=-2,Vnx+3
!       V(i,j,k) = V2(i,j,k)
!      enddo
!     enddo
!    enddo
!
!    !$hmppcg grid blocksize 512x1
!    !$hmppcg permute (k,i,j)
!    do k=-2,Wnz+3
!     do j=-2,Wny+3
!      do i=-2,Wnx+3
!       W(i,j,k) = W2(i,j,k)
!      enddo
!     enddo
!    enddo
!  end subroutine UpdateU










  subroutine OtherTerms(U,V,W,U2,V2,W2,Pr,coef)
   real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
   real(KND),intent(inout):: Pr(1:,1:,1:)
   real(KND),intent(inout):: U2(-2:,-2:,-2:),V2(-2:,-2:,-2:),W2(-2:,-2:,-2:)
   real(KND),intent(in):: coef

   real(KND),dimension(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),lbound(U,3):ubound(U,3)):: U3
   real(KND),dimension(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),lbound(V,3):ubound(V,3)):: V3
   real(KND),dimension(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),lbound(W,3):ubound(W,3)):: W3
   real(KND) Ap,Apre,Aprn,Aprt,S

   integer i,j,k,it

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


     !$hmpp <tsteps> advancedload, args[ForwEul::Visc]

     !$hmpp <tsteps> ForwEul callsite, args[*].noupdate = true
     call ForwEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                  dxPr,dyPr,dzPr,dxU,dyV,dzW,&
                  U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                  dt,coef)


     call IBMomentum(U2,V2,W2,U3,V3,W3)


#ifdef __HMPP
     if (gridtype==UNIFORMGRID.and.GPU>0) then                  !Performs the diffusion terms
      write (*,*) "GPU CN call"

      !$hmpp <tsteps>  advancedload, args[UnifRedBlack::maxCNiter,UnifRedBlack::epsCN]

      !$hmpp <tsteps> UNIFREDBLACK callsite, args[*].noupdate=true
      call UNIFREDBLACK_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,&
                            Btype,sideU,&
                            dt,dxmin,dymin,dzmin,&
                            Uin,Vin,Win,&
                            U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                            coef,maxCNiter,epsCN,it,S)

     !$hmpp <tsteps> delegatedstore, args[UnifRedBlack::iters,UnifRedBlack::residuum]

      write(*,*) "back from GPU CN", it,S


     else&
#endif
     if (gridtype==UNIFORMGRID) then
      call UNIFREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
     else
      call GENREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
     endif


   else  Re_gt_0  !Re<=0

    U2 = U+U2
    V2 = V+V2
    W2 = W+W2

   endif   Re_gt_0

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
      enddo
     enddo
    enddo

    S = S/(Unx*Uny*Unz)
    write(*,*) "Mean friction:", S
   endif
  endsubroutine OtherTerms




  subroutine UNIFREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND),dimension(-2:,-2:,-2:),intent(inout):: U2,V2,W2,U3,V3,W3
   real(KND),intent(in):: coef
   real(KND),dimension(1:Unx,1:Uny,1:Unz):: Apu
   real(KND),dimension(1:Vnx,1:Vny,1:Vnz):: ApV
   real(KND),dimension(1:Wnx,1:Wny,1:Wnz):: ApW
   real(KND) recdxmin2,recdymin2,recdzmin2                                                               !reciprocal values of dx**2
   real(KND) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
   integer i,j,k,bi,bj,bk,l
   integer,parameter :: narr = 3, narr2 = 5 !number of arrays in the loop

       Ap = coef*dt/(2._KND)
       S = 0
       l = 0

       recdxmin2=1./dxmin**2
       recdymin2=1./dymin**2
       recdzmin2=1./dzmin**2

       !$omp parallel private(i,j,k,bi,bj,bk,Suavg,Svavg,Swavg)

       !The explicit part, which doesn't have to be changed inside the loop
       !$omp do
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
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
       !$omp end do nowait
       !$omp do
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
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
       !$omp end do nowait
       !$omp do
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
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
       !$omp end do nowait

       !Auxiliary coefficients to better efficiency in loops
       !$omp do
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
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
       !$omp end do

       !$omp workshare
       ApU = 1._KND/(1._KND+Ap*ApU)
       !$omp end workshare nowait

       !$omp do
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
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
       !$omp end do

       !$omp workshare
       ApV = 1._KND/(1._KND+Ap*ApV)
       !$omp end workshare nowait

       !$omp do
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
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
       !$omp end do

       !$omp workshare
       ApW = 1._KND/(1._KND+Ap*ApW)
       !$omp end workshare nowait
       !$omp end parallel


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
        !$OMP PARALLEL PRIVATE(i,j,k,p) REDUCTION(max:Su,Sv,Sw)
        !$OMP DO
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
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
        !$OMP ENDDO NOWAIT
        !$OMP DO
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
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
        !$OMP ENDDO NOWAIT
        !$OMP DO
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
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
        !$OMP ENDDO

        !$OMP DO
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
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
        !$OMP ENDDO NOWAIT
        !$OMP DO
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
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
        !$OMP ENDDO NOWAIT
        !$OMP DO
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
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
        !$OMP ENDDO
        !$OMP END PARALLEL
        S = max(Su/Suavg,Sv/Svavg,Sw/Swavg)
        write (*,*) "CN ",l,S

        if (S<=epsCN) exit
       enddo

       !$omp parallel
       !$omp workshare
       U2 = U3
       V2 = V3
       W2 = W3
       !$omp end workshare
       !$omp end parallel

  endsubroutine UNIFREDBLACK



  subroutine GENREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND),dimension(-2:,-2:,-2:),intent(inout):: U2,V2,W2,U3,V3,W3
   real(KND),intent(in):: coef
   real(KND),dimension(1:Unx,1:Uny,1:Unz):: Apu
   real(KND),dimension(1:Vnx,1:Vny,1:Vnz):: ApV
   real(KND),dimension(1:Wnx,1:Wny,1:Wnz):: ApW
   real(KND) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
   integer i,j,k,l
   integer ind(3)


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
         enddo
        enddo
       enddo
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
         enddo
        enddo
       enddo
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
        enddo
       enddo
      enddo

       do k=1,Unz         !Auxiliary coefficients to better efficiency in loops
        do j=1,Uny
         do i=1,Unx
          ApU(i,j,k) = ((Visc(i+1,j,k)/dxPr(i+1)+&
                      Visc(i,j,k)/dxPr(i))/dxU(i)+&
                      (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))/dyV(j)+&
                      0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))/dyV(j-1))/dyPr(j)+&
                      (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))/dzW(k)+&
                      0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))/dzW(k-1))/dzPr(k))
         enddo
        enddo
       enddo

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
         enddo
        enddo
       enddo

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
         enddo
        enddo
       enddo

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
        !$OMP PARALLEL PRIVATE(i,j,k,p) REDUCTION(max:Su,Sv,Sw)
        !$OMP DO
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
          enddo
         enddo
        enddo
        !$OMP ENDDO NOWAIT
        !$OMP DO
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
          enddo
         enddo
        enddo
        !$OMP ENDDO NOWAIT
        !$OMP DO
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
          enddo
         enddo
        enddo
        !$OMP ENDDO

        !$OMP DO
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
          enddo
         enddo
        enddo
        !$OMP ENDDO NOWAIT
        !$OMP DO
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
          enddo
         enddo
        enddo
        !$OMP ENDDO NOWAIT
        !$OMP DO
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
          enddo
         enddo
        enddo
        !$OMP ENDDO
        !$OMP END PARALLEL
        S = max(Su/Suavg,Sv/Svavg,Sw/Swavg)
        write (*,*) "CN ",l,S
        if (S<=epsCN) exit
       enddo

       U2 = U3
       V2 = V3
       W2 = W3

   endsubroutine GENREDBLACK





  subroutine SubgridStresses(U,V,W,Pr)
  use Geometric, only: ScalFlIBPoints, TIBPoint_Viscosity

  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  integer i

     if (sgstype==SmagorinskyModel) then
                       call Smag(U,V,W)
     elseif (sgstype==VremanModel) then
#ifdef __HMPP
                       if (GPU>0.and.gridtype==uniformgrid.and. Prnx*Prny*Prnz > 50) then
                           !$hmpp <tsteps> Vreman callsite, args[*].noupdate = true
                           call Vreman_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxmin,dymin,dzmin,dt,Re,U,V,W,Visc)
                           !$hmpp <tsteps> delegatedstore, args[Vreman::Visc]
                       else
                           call Vreman(U,V,W)
                       endif
#else
                       call Vreman(U,V,W)
#endif
     elseif (sgstype==StabSmagorinskyModel) then
                       call StabSmag(U,V,W)
     else
                       Visc = 1._KND/Re
     endif


     if (debuglevel>0) then
      write(*,*) "NUt", sum(Visc(1:Prnx,1:Prny,1:Prnz))/(Prnx*Prny*Prnz)
      write(*,*) "maxNUt", MAXVAL(Visc(1:Prnx,1:Prny,1:Prnz))
      write(*,*) "minNUt", MINVAL(Visc(1:Prnx,1:Prny,1:Prnz))
     endif

     if (wallmodeltype>0) then
                     call ComputeViscsWM(U,V,W,Pr)
     endif

     do i=1,size(ScalFlIBPoints)
       Visc(ScalFlIBPoints(i)%xi,ScalFlIBPoints(i)%yj,ScalFlIBPoints(i)%zk) = &
                                               TIBPoint_Viscosity(ScalFlIBPoints(i),Visc)
     enddo

     call Bound_Visc(Visc)

  end subroutine SubgridStresses

























  subroutine IBMomentum(U2,V2,W2,U3,V3,W3)
    use GEOMETRIC, only: UIBPoints, VIBPoints, WIBPoints, TIBPoint_MomentumSource

    real(KND),dimension(-2:,-2:,-2:),intent(inout):: U2,V2,W2,U3,V3,W3
    integer   :: i,xi,yj,zk
    real(KND) :: src


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
     enddo
     !$omp enddo nowait

     !$omp do private(xi,yj,zk,src)
     do i=1,ubound(VIBPoints,1)
       xi = VIBPoints(i)%xi
       yj = VIBPoints(i)%yj
       zk = VIBPoints(i)%zk

       src = TIBPoint_MomentumSource(VIBPoints(i),V3)

       V3(xi,yj,zk) = V3(xi,yj,zk) + src * dt
       V2(xi,yj,zk) = V2(xi,yj,zk) + src * dt
     enddo
     !$omp enddo nowait

     !$omp do private(xi,yj,zk,src)
     do i=1,ubound(WIBPoints,1)
       xi = WIBPoints(i)%xi
       yj = WIBPoints(i)%yj
       zk = WIBPoints(i)%zk

       src = TIBPoint_MomentumSource(WIBPoints(i),W3)

       W3(xi,yj,zk) = W3(xi,yj,zk) + src * dt
       W2(xi,yj,zk) = W2(xi,yj,zk) + src * dt
     enddo
     !$omp enddo
     !$omp end parallel

     !$hmpp <tsteps> advancedload, args[UnifRedBlack::U3,UnifRedBlack::V3,UnifRedBlack::W3]
  end subroutine IBMomentum






  subroutine IBMassSources(Q,U,V,W)
    use GEOMETRIC, only: UIBPoints, VIBPoints, WIBPoints

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
    enddo


    do i=1,ubound(VIBPoints,1)
      xi = VIBPoints(i)%xi
      yj = VIBPoints(i)%yj
      zk = VIBPoints(i)%zk

      Q(xi,yj,zk)   = Q(xi,yj,zk)   + V(xi,yj,zk)/dyPr(yj)
      Q(xi,yj+1,zk) = Q(xi,yj+1,zk) - V(xi,yj,zk)/dyPr(yj+1)
    enddo


    do i=1,ubound(WIBPoints,1)
      xi = WIBPoints(i)%xi
      yj = WIBPoints(i)%yj
      zk = WIBPoints(i)%zk

      Q(xi,yj,zk)   = Q(xi,yj,zk)   + W(xi,yj,zk)/dzPr(zk)
      Q(xi,yj,zk+1) = Q(xi,yj,zk+1) - W(xi,yj,zk)/dzPr(zk+1)
    enddo

    call Bound_Q(Q)

  end subroutine IBMassSources








end module TSTEPS











module TSTEPS

  use CDS, only: CDU, CDV, CDW, KappaU, KappaV, KappaW
  use LAXFRIED
  use LAXWEND
  use PARAMETERS
  use BOUNDARIES, only: BoundU,Bound_Q
  use POISSON, only: Pr_Correct
  use OUTPUTS, only: store,display,proftempfl
  use SCALARS, only: ScalarRK3
  use SMAGORINSKY, only: Smag, StabSmag, Vreman
  use TURBINLET, only: GetTurbInlet, GetInletFromFile
  use Wallmodels, only: ComputeViscsWM

  implicit none


  private
  public TMarchRK3

contains



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
    !$hmpp <tsteps>     BoundV::nx,BoundV::ny,BoundV::nz,BoundW::nx,BoundW::ny,BoundW::nz]

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


    if (debugparam>1) call system_clock(count=time1)

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

    call SubgridStresses(U,V,W,Pr)


    !$hmpp <tsteps> advancedload, args[Convection::lev]


    !$hmpp <tsteps> Convection callsite, args[*].noupdate = true
    call Convection(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,buoyancy,convmet,&
                       dxmin,dymin,dzmin,coriolisparam,grav_acc,temperature_ref,&
                       U,V,W,U2,V2,W2,Ustar,Vstar,Wstar,temperature,beta,rho,l,dt)


    call OtherTerms(U,V,W,U2,V2,W2,Pr,2._KND*alpha(l))



    if (Btype(To) ==FreeSlipBuff)  then
        !$hmpp <tsteps> AttenuateTop callsite, args[*].noupdate = true
        call AttenuateTop(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Btype,&
                          zPr,zW,U2,V2,W2,temperature,buoyancy)
    endif

    if (Btype(Ea) ==OutletBuff) then
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
!     if (debuglevel>0) then
      delta = delta+sum(abs(U(1:Unx,1:Uny,1:Unz)-U2(1:Unx,1:Uny,1:Unz)))/(Unx*Uny*Unz)
      delta = delta+sum(abs(V(1:Vnx,1:Vny,1:Vnz)-V2(1:Vnx,1:Vny,1:Vnz)))/(Vnx*Vny*Vnz)
      delta = delta+sum(abs(W(1:Wnx,1:Wny,1:Wnz)-W2(1:Wnx,1:Wny,1:Wnz)))/(Wnx*Wny*Wnz)
!     endif

    U = U2
    V = V2
    W = W2



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

    if (Btype(Ea) ==OutletBuff) then
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
!  integer,parameter :: KND = 4
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






 !$hmpp <tsteps> Convection codelet
 subroutine Convection(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,buoyancy,convmet,&
                       dxmin,dymin,dzmin,coriolisparam,grav_acc,temperature_ref,&
                       U,V,W,U2,V2,W2,Ustar,Vstar,Wstar,temperature,beta,rho,lev,dt)
  implicit none
#ifdef __HMPP
  integer,parameter :: KND = 4
  intrinsic abs
#endif

  integer,intent(in)   :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,buoyancy,convmet,lev
  real(KND),intent(in) :: dxmin,dymin,dzmin,coriolisparam,grav_acc,temperature_ref
  real(KND),intent(in) :: dt
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(in)    :: U
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(out)   :: U2
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(inout) :: Ustar
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(in)    :: V
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(out)   :: V2
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(inout) :: Vstar
  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(in)    :: W
  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(out)   :: W2
  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(inout) :: Wstar
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(in) :: temperature
  real(KND),dimension(1:3),intent(in) :: beta,rho
  integer i,j,k

      if (lev>1) then
        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Unz
         do j=1,Uny
          do i=1,Unx
           U2(i,j,k) = Ustar(i,j,k)*rho(lev)
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Vnz
         do j=1,Vny
          do i=1,Vnx
           V2(i,j,k) = Vstar(i,j,k)*rho(lev)
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Wnz
         do j=1,Wny
          do i=1,Wnx
           W2(i,j,k) = Wstar(i,j,k)*rho(lev)
          enddo
         enddo
        enddo
      else
        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Unz
         do j=1,Uny
          do i=1,Unx
           Ustar(i,j,k) = 0
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Vnz
         do j=1,Vny
          do i=1,Vnx
           Vstar(i,j,k) = 0
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Wnz
         do j=1,Wny
          do i=1,Wnx
           Wstar(i,j,k) = 0
          enddo
         enddo
        enddo

        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Unz
         do j=1,Uny
          do i=1,Unx
           U2(i,j,k) = 0
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Vnz
         do j=1,Vny
          do i=1,Vnx
           V2(i,j,k) = 0
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Wnz
         do j=1,Wny
          do i=1,Wnx
           W2(i,j,k) = 0
          enddo
         enddo
        enddo
      endif

      if (convmet>0) then

         call CDS_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxmin,dymin,dzmin,dt,Ustar,Vstar,Wstar,U,V,W)

      else

        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Unz
         do j=1,Uny
          do i=1,Unx
           Ustar(i,j,k) = 0
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Vnz
         do j=1,Vny
          do i=1,Vnx
           Vstar(i,j,k) = 0
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg permute (k,i,j)
        do k=1,Wnz
         do j=1,Wny
          do i=1,Wnx
           Wstar(i,j,k) = 0
          enddo
         enddo
        enddo

      endif


      if (abs(coriolisparam)>1E-10_KND) call CoriolisForce(Unx,Uny,Unz,Vnx,Vny,Vnz,&
                                                    coriolisparam,&
                                                    Ustar,Vstar,U,V,dt)

      if (buoyancy==1) call BuoyancyForce(Prnx,Prny,Prnz,Wnx,Wny,Wnz,&
                                grav_acc,temperature_ref,Wstar,temperature,dt)

      !$hmppcg grid blocksize 512x1
      !$hmppcg permute (k,i,j)
      do k=1,Unz
       do j=1,Uny
        do i=1,Unx
         U2(i,j,k) = U2(i,j,k)+Ustar(i,j,k)*beta(lev)
        enddo
       enddo
      enddo
      !$hmppcg grid blocksize 512x1
      !$hmppcg permute (k,i,j)
      do k=1,Vnz
       do j=1,Vny
        do i=1,Vnx
         V2(i,j,k) = V2(i,j,k)+Vstar(i,j,k)*beta(lev)
        enddo
       enddo
      enddo
      !$hmppcg grid blocksize 512x1
      !$hmppcg permute (k,i,j)
      do k=1,Wnz
       do j=1,Wny
        do i=1,Wnx
         W2(i,j,k) = W2(i,j,k)+Wstar(i,j,k)*beta(lev)
        enddo
       enddo
      enddo

  end subroutine Convection











  !$hmpp <tsteps> TimeStepEul codelet
  subroutine TimeStepEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
               dxmin,dxU,dyV,dzW,CFL,Uref,steady,time,endtime,&
               U,V,W,dt)
  implicit none
#ifdef __HMPP
  integer,parameter :: KND = 4,TIM = 4
  intrinsic min,max,abs
#endif
  integer,intent(in)   :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,steady
  real(KND),intent(in) :: dxmin,dxU(-2:Prnx+2),dyV(-2:Prny+2),dzW(-2:Prnz+2)
  real(KND),intent(in) :: CFL,Uref
  real(TIM),intent(in) :: time,endtime
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(in)  :: U
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(in)  :: V
  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(in)  :: W
  real(TIM),intent(out) :: dt
  integer i,j,k
  real(KND) m,p

    m = 0

   !$hmppcg grid blocksize 512x1
   !$hmppcg permute (k,i,j)
   !$hmppcg gridify (k,i), private(p), reduce (max:m)
   do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       p = MAX(abs(U(i,j,k)/dxU(i)),abs(U(i-1,j,k)/dxU(i-1)))
       p = MAX(p,abs(V(i,j,k)/dyV(j)),abs(V(i,j-1,k)/dyV(j-1)))
       p = MAX(p,abs(W(i,j,k)/dzW(k)),abs(W(i,j,k-1)/dzW(k-1)))
       m = max(m,p)
      enddo
     enddo
    enddo


    if (m>0) then
     dt = MIN(CFL/m,dxmin/Uref)
    else
     dt = dxmin/Uref
    endif

    if (steady/=1.and.dt+time>endtime)  dt = endtime-time

  endsubroutine TimeStepEul







  pure subroutine BuoyancyForce(Prnx,Prny,Prnz,Wnx,Wny,Wnz,grav_acc,temperature_ref,W2,theta,dt)
  implicit none
#ifdef __HMPP
  integer,parameter :: KND = 4
#endif
  integer,intent(in) :: Prnx,Prny,Prnz,Wnx,Wny,Wnz
  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(inout) :: W2
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(in) :: theta
  real(KND),intent(in) :: grav_acc,temperature_ref,dt
  real(KND) A
  integer i,j,k

    A = grav_acc*dt/temperature_ref
    !$hmppcg grid blocksize 512x1
    !$hmppcg permute (k,i,j)
    do k=1,Wnz
     do j=1,Wny
      do i=1,Wnx
            W2(i,j,k) = W2(i,j,k)+A*((theta(i,j,k+1)+theta(i,j,k))/2._KND-temperature_ref)
      enddo
     enddo
    enddo
  endsubroutine BuoyancyForce





  pure subroutine CoriolisForce(Unx,Uny,Unz,Vnx,Vny,Vnz,&
                                 coriolisparam,&
                                 U2,V2,U,V,dt)
  implicit none
#ifdef __HMPP
  integer,parameter :: KND = 4
#endif
  integer,intent(in) :: Unx,Uny,Unz,Vnx,Vny,Vnz
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(in)    :: U
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(in)    :: V
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(inout) :: U2
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(inout) :: V2
  real(KND),intent(in):: coriolisparam,dt
  real(KND) A
  integer i,j,k

   A=-dt
   if (coriolisparam>0) then
   !$hmppcg grid blocksize 512x1
   !$hmppcg permute (k,i,j)
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
           U2(i,j,k) = U2(i,j,k)-A*coriolisparam*(V(i,j-1,k)+V(i+1,j-1,k)+V(i,j,k)+V(i+1,j,k))/4._KND
      enddo
     enddo
    enddo
    !$hmppcg grid blocksize 512x1
    !$hmppcg permute (k,i,j)
    do k=1,Vnz
     do j=1,Vny
      do i=1,Vnx
           V2(i,j,k) = V2(i,j,k)+A*coriolisparam*(U(i-1,j,k)+U(i-1,j+1,k)+U(i,j,k)+U(i,j+1,k))/4._KND
      enddo
     enddo
    enddo
   endif
  endsubroutine CoriolisForce







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


     if (gridtype==UNIFORMGRID.and.GPU>0) then                  !Performs the diffusion terms
      write (*,*) "GPU CN call"

      !$hmpp <tsteps>  advancedload, args[UnifRedBlack::maxCNiter,UnifRedBlack::epsCN]


      !$hmpp <tsteps> UNIFREDBLACK callsite, args[*].noupdate = true
      call UNIFREDBLACK_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,&
                            Btype,sideU,&
                            dt,dxmin,dymin,dzmin,&
                            Uin,Vin,Win,&
                            U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                            coef,maxCNiter,epsCN,it,S)

     !$hmpp <tsteps> delegatedstore, args[UnifRedBlack::iters,UnifRedBlack::residuum]

      write(*,*) "back from GPU CN", it,S


     else if (gridtype==UNIFORMGRID) then
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


  !$hmpp <tsteps> PressureGrad codelet
  subroutine PressureGrad(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxU,dyV,dzW,&
                          Btype,prgradientx,prgradienty,&
                          Pr,U,V,W,&
                          dt,coef)
  implicit none
#ifdef __HMPP
  integer,parameter :: KND = 4
#endif
  integer,intent(in)      :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz
  real(KND),intent(in)    :: dxU(-2:Prnx+2),dyV(-2:Prny+2),dzW(-2:Prnz+2),prgradientx,prgradienty,dt,coef
  integer,intent(in)      :: Btype(6)
  real(KND),intent(inout) :: Pr(1:Unx+1,1:Vny+1,1:Wnz+1)
  real(KND),intent(inout) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
  real(KND),intent(inout) :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
  real(KND),intent(inout) :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(KND) :: A
  integer i,j,k

   call BoundPr_GPU(Unx,Vny,Wnz,Prnx,Prny,Prnz,Btype,Pr)

   A=-coef*dt
   A=-coef*dt
   A=-coef*dt

   !$omp parallel
   !$omp do
   !$hmppcg grid blocksize 512x1
   !$hmppcg permute (k,i,j)
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
          U(i,j,k) = U(i,j,k)+A*(Pr(i+1,j,k)-Pr(i,j,k))/dxU(i)+A*prgradientx
     enddo
    enddo
   enddo
   !$omp enddo nowait
   !$omp do
   !$hmppcg grid blocksize 512x1
   !$hmppcg permute (k,i,j)
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
          V(i,j,k) = V(i,j,k)+A*(Pr(i,j+1,k)-Pr(i,j,k))/dyV(j)+A*prgradienty
     enddo
    enddo
   enddo
   !$omp enddo nowait
   !$omp do
   !$hmppcg grid blocksize 512x1
   !$hmppcg permute (k,i,j)
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
          W(i,j,k) = W(i,j,k)+A*(Pr(i,j,k+1)-Pr(i,j,k))/dzW(k)
     enddo
    enddo
   enddo
   !$omp enddo nowait
   !$omp end parallel
  end subroutine PressureGrad



  !$hmpp <tsteps> ForwEul codelet
  subroutine ForwEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                    dxPr,dyPr,dzPr,dxU,dyV,dzW,&
                    U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                    dt,coef)
  implicit none
#ifdef __HMPP
  integer,parameter :: KND = 4
#endif
  integer,intent(in) :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz

  real(KND),intent(in):: U(-2:Unx+3,-2:Uny+3,-2:Unz+3),V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(KND),intent(in):: U2(-2:Unx+3,-2:Uny+3,-2:Unz+3),V2(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),W2(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(KND),intent(out):: U3(-2:Unx+3,-2:Uny+3,-2:Unz+3),V3(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),W3(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(KND),intent(in):: Visc(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  real(KND),intent(in) :: dxU(-2:Prnx+2),dyV(-2:Prny+2),dzW(-2:Prnz+2)
  real(KND),intent(in) :: dxPr(-2:Prnx+3),dyPr(-2:Prny+3),dzPr(-2:Prnz+3),dt,coef
  real(KND) :: Ap
  integer i,j,k

     Ap = coef*dt

     !$hmppcg grid blocksize 512x1
     !$hmppcg permute (k,i,j)
     do k=1,Unz    !Forward Euler for the first approximation
      do j=1,Uny
       do i=1,Unx
        U3(i,j,k) = U(i,j,k)+U2(i,j,k)+Ap*(&
        ((Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))/dxPr(i+1)-&
        Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k))/dxPr(i))/dxU(i)+&
         0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))/dyV(j)-&
         (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k))/dyV(j-1))/dyPr(j)+&
         ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)-&
         (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1))/dzW(k-1))/dzPr(k))))
       enddo
      enddo
     enddo
     !$hmppcg grid blocksize 512x1
     !$hmppcg permute (k,i,j)
     do k=1,Vnz
      do j=1,Vny
       do i=1,Vnx
        V3(i,j,k) = V(i,j,k)+V2(i,j,k)+Ap*(&
        ((Visc(i,j+1,k)*(V(i,j+1,k)-V(i,j,k))/dyPr(j+1)-&
         Visc(i,j,k)*(V(i,j,k)-V(i,j-1,k))/dyPr(j))/dyV(j)+&
         0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V(i+1,j,k)-V(i,j,k))/dxU(i)-&
        (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(V(i,j,k)-V(i-1,j,k))/dxU(i-1))/dxPr(i)+&
         ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))/dzW(k)-&
         (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(V(i,j,k)-V(i,j,k-1))/dzW(k-1))/dzPr(k))))
       enddo
      enddo
     enddo
     !$hmppcg grid blocksize 512x1
     !$hmppcg permute (k,i,j)
     do k=1,Wnz
      do j=1,Wny
       do i=1,Wnx
        W3(i,j,k) = W(i,j,k)+W2(i,j,k)+Ap*(&
        ((0.25_KND*((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W(i+1,j,k)-W(i,j,k))/dxU(i)-&
        (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(W(i,j,k)-W(i-1,j,k))/dxU(i-1))/dxPr(i)+&
         ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W(i,j+1,k)-W(i,j,k))/dyV(j)-&
         (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(W(i,j,k)-W(i,j-1,k))/dyV(j-1))/dyPr(j))+&
         (Visc(i,j,k+1)*(W(i,j,k+1)-W(i,j,k))/dzPr(k+1)-&
         Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1))/dzPr(k))/dzW(k)))
       enddo
      enddo
     enddo
  end subroutine ForwEul
















  subroutine UNIFREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND),dimension(-2:,-2:,-2:),intent(inout):: U2,V2,W2,U3,V3,W3
   real(KND),intent(in):: coef
   real(KND),dimension(1:Unx,1:Uny,1:Unz):: Apu
   real(KND),dimension(1:Vnx,1:Vny,1:Vnz):: ApV
   real(KND),dimension(1:Wnx,1:Wny,1:Wnz):: ApW
   real(KND) recdxmin2,recdymin2,recdzmin2                                                               !reciprocal values of dx**2
   real(KND) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
   integer i,j,k,l

       Ap = coef*dt/(2._KND)
       S = 0
       l = 0

       recdxmin2=1./dxmin**2
       recdymin2=1./dymin**2
       recdzmin2=1./dzmin**2

       do k=1,Unz    !The explicit part, which doesn't have to be changed inside the loop
        do j=1,Uny
         do i=1,Unx
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
       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
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
       do k=1,Wnz
       do j=1,Wny
        do i=1,Wnx
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

       do k=1,Unz         !Auxiliary coefficients to better efficiency in loops
        do j=1,Uny
         do i=1,Unx
          ApU(i,j,k) = ((Visc(i+1,j,k)+&
                      Visc(i,j,k))*recdxmin2+&
                      0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))+&
                      (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k)))*recdymin2+&
                      ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                      (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1)))*recdzmin2))
         enddo
        enddo
       enddo

       ApU = 1._KND/(1._KND+Ap*ApU)

       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          ApV(i,j,k) = ((Visc(i,j+1,k)+&
                     Visc(i,j,k))*recdymin2+&
                     0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))+&
                      (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k)))*recdxmin2+&
                     ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                     (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1)))*recdzmin2))
         enddo
        enddo
       enddo

       ApV = 1._KND/(1._KND+Ap*ApV)

       do k=1,Wnz
        do j=1,Wny
         do i=1,Wnx
          ApW(i,j,k) = (0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                      (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k)))*recdxmin2+&
                     ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))+&
                     (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k)))*recdymin2)+&
                     (Visc(i,j,k+1)+&
                     Visc(i,j,k))*recdzmin2)
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
        !$OMP ENDDO NOWAIT
        !$OMP DO
        do k=1,Vnz
         do j=1,Vny
          do i=1+mod(j+k,2),Vnx,2
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
        !$OMP ENDDO NOWAIT
        !$OMP DO
        do k=1,Wnz
         do j=1,Wny
          do i=1+mod(j+k,2),Wnx,2
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
        !$OMP ENDDO

        !$OMP DO
        do k=1,Unz
         do j=1,Uny
          do i=1+mod(j+k+1,2),Unx,2
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
        !$OMP ENDDO NOWAIT
        !$OMP DO
        do k=1,Vnz
         do j=1,Vny
          do i=1+mod(j+k+1,2),Vnx,2
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
        !$OMP ENDDO NOWAIT
        !$OMP DO
        do k=1,Wnz
         do j=1,Wny
          do i=1+mod(j+k+1,2),Wnx,2
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
        !$OMP ENDDO
        !$OMP END PARALLEL
        S = max(Su/Suavg,Sv/Svavg,Sw/Swavg)
        write (*,*) "CN ",l,S
        if (S<=epsCN) exit
       enddo

       U2 = U3
       V2 = V3
       W2 = W3

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






  !Slower using HMPP than CPU, but if we avoid memory transfer, it is still profitable.

  !$hmpp <tsteps> AttenuateTop codelet
  subroutine AttenuateTop(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Btype,zPr,zW,U,V,W,temperature,buoyancy)
  implicit none
#ifdef __HMPP
   integer, parameter:: KND = 4
   integer,parameter       :: NOSLIP=1, FREESLIP=2, PERIODIC=3, DIRICHLET=4, NEUMANN=5, CONSTFLUX=6,&  !boundary condition types
                                TURBULENTINLET=7, FREESLIPBUFF=8, OUTLETBUFF=9, INLETFROMFILE=10
   integer, parameter :: Ea=1,We=2,So=3,No=4,Bo=5,To=6
#endif
  integer,intent(in)      :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,buoyancy
  integer,intent(in)      :: Btype(6)
  real(KND),intent(in)    :: zPr(-2:Prnz+3)
  real(KND),intent(in)    :: zW(-3:Prnz+3)
  real(KND),intent(inout) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
  real(KND),intent(inout) :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
  real(KND),intent(inout) :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(inout) :: temperature
  integer i,j,k,bufn,mini,maxi,maxUi
  real(KND) ze,zs,zb,p
  real(KND),dimension(:),allocatable :: DF,avg
  intrinsic max,min

    if (Btype(We)==DIRICHLET.or.Btype(We)==TURBULENTINLET.or.Btype(We)==INLETFROMFILE) then
      mini = min(5,Unx)
    else
      mini = 1
    endif

    if (Btype(Ea)==DIRICHLET.or.Btype(We)==TURBULENTINLET.or.Btype(We)==OUTLETBUFF) then
      maxi = max(1,Prnx-5)
      maxUi = max(1,Unx-5)
    else
      maxi = Prnx
      maxUi = Unx
    endif

    bufn = max(5,Prnz/4)
    zs = zW(Prnz-bufn)
    ze = zW(Prnz)

    allocate(DF(min(Unz,Vnz,Wnz)-bufn:max(Unz,Vnz,Wnz)))
    allocate(avg(min(Unz,Vnz,Wnz)-bufn:max(Unz,Vnz,Wnz)))



    do k=Unz-bufn,Unz
      avg(k) = 0
    enddo

    do k=Unz-bufn,Unz
      p = 0
      !$hmppcg grid blocksize 512x1
      !$hmppcg gridify (j,i) global(p), reduce(+:p)
      do j=1,Uny
        do i=mini,maxUi
          p = p+U(i,j,k)
        enddo
      enddo
      avg(k) = p
    enddo

    do k=Unz-bufn,Unz
      avg(k) = avg(k)/((maxUi-mini+1)*Uny)
    enddo

    do k=Unz-bufn,Unz
      zb=(zPr(k)-zs)/(ze-zs)
      DF(k) = DampF(zb)
    enddo

    !$hmppcg grid blocksize 512x1
    !$hmppcg permute(k,i,j)
    !$hmppcg gridify (k,i)
    do k=Unz-bufn,Unz
      do j=-1,Uny+1
        do i=-1,Unx+1
          U(i,j,k) = avg(k)+DF(k)*(U(i,j,k)-avg(k))
        enddo
      enddo
    enddo



    do k=Vnz-bufn,Vnz
      avg(k) = 0
    enddo
    do k=Vnz-bufn,Vnz
      p = 0
      !$hmppcg grid blocksize 512x1
      !$hmppcg gridify (j,i) global(p), reduce(+:p)
      do j=1,Vny
        do i=mini,maxi
          p = p+V(i,j,k)
        enddo
      enddo
      avg(k) = p
    enddo
    do k=Vnz-bufn,Vnz
      avg(k) = avg(k)/((maxi-mini+1)*Vny)
    enddo
    do k=Vnz-bufn,Vnz
      zb=(zPr(k)-zs)/(ze-zs)
      DF(k) = DampF(zb)
    enddo
    !$hmppcg grid blocksize 512x1
    !$hmppcg permute(k,i,j)
    !$hmppcg gridify (k,i)
    do k=Vnz-bufn,Vnz
      do j=-1,Vny+1
        do i=-1,Vnx+1
          V(i,j,k) = avg(k)+DF(k)*(V(i,j,k)-avg(k))
        enddo
      enddo
    enddo



    do k=Wnz-bufn,Wnz
      avg(k) = 0
    enddo

    do k=Wnz-bufn,Wnz
      p = 0
      !$hmppcg grid blocksize 512x1
      !$hmppcg gridify (j,i) global(p), reduce(+:p)
      do j=1,Wny
        do i=mini,maxi
          p = p+W(i,j,k)
        enddo
      enddo
      avg(k) = p
    enddo

    do k=Wnz-bufn,Wnz
      avg(k) = avg(k)/((maxi-mini+1)*Wny)
    enddo

    do k=Wnz-bufn,Wnz
      zb=(zW(k)-zs)/(ze-zs)
      DF(k) = DampF(zb)
    enddo

    !$hmppcg grid blocksize 512x1
    !$hmppcg permute(k,i,j)
    !$hmppcg gridify (k,i)
    do k=Wnz-bufn,Wnz
      do j=-1,Wny+1
        do i=-1,Wnx+1
          W(i,j,k) = avg(k)+DF(k)*(W(i,j,k)-avg(k))
        enddo
      enddo
    enddo



    if (buoyancy==1) then

      do k=Prnz-bufn,Prnz
        avg(k) = 0
      enddo

      do k=Prnz-bufn,Prnz
        p = 0
        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify (j,i) global(p), reduce(+:p)
        do j=1,Prny
          do i=mini,maxi
            p = p+temperature(i,j,k)
          enddo
        enddo
        avg(k) = p
      enddo

      do k=Prnz-bufn,Prnz
        avg(k) = avg(k)/((maxi-mini+1)*Prny)
      enddo

      do k=Prnz-bufn,Prnz
        zb=(zPr(k)-zs)/(ze-zs)
        DF(k) = DampF(zb)
      enddo

      !$hmppcg grid blocksize 512x1
      !$hmppcg permute(k,i,j)
      !$hmppcg gridify (k,i)
      do k=Prnz-bufn,Prnz
        do j=-1,Prny+1
          do i=-1,Prnx+1
            temperature(i,j,k) = avg(k)+DF(k)*(temperature(i,j,k)-avg(k))
          enddo
        enddo
      enddo

    endif

  endsubroutine AttenuateTop



  !$hmpp <tsteps> AttenuateOut codelet
  !$hmpp <tsteps> AttenuateOut2 codelet
  subroutine AttenuateOut(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,xPr,xU,U,V,W,temperature,buoyancy)
  implicit none
#ifdef __HMPP
   integer, parameter:: KND = 4
#endif
  integer,intent(in)      :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,buoyancy
  real(KND),intent(in)    :: xPr(-2:Prnx+3)
  real(KND),intent(in)    :: xU(-3:Prnx+3)
  real(KND),intent(inout) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
  real(KND),intent(inout) :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
  real(KND),intent(inout) :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(inout) :: temperature
  integer i,j,k,bufn
  real(KND) p,xe,xs,xb,DF
  intrinsic max

    bufn = max(10,Prnx/8)
    xs = xU(Prnx-bufn)
    xe = xU(Prnx)

    !$hmppcg grid blocksize 512x1
    !$hmppcg gridify (k,j) private(p,xb,DF)
    do k=1,Unz
      do j=1,Uny
        p = 0
        do i=2*Unx/3,Unx-4
          p = p+U(i,j,k)
        enddo
        p = p/(Unx-4-2*Unx/3+1)
        do i=Unx-bufn,Unx+1
          xb=(xU(i)-xs)/(xe-xs)
          DF = DampF(xb)
          U(i,j,k) = p+DF*(U(i,j,k)-p)
        enddo
      enddo
    enddo

    !$hmppcg grid blocksize 512x1
    !$hmppcg gridify (k,j) private(p,xb,DF)
    do k=1,Vnz
      do j=1,Vny
        p = 0
        do i=2*Vnx/3,Vnx-4
          p = p+V(i,j,k)
        enddo
        p = p/(Vnx-4-2*Vnx/3+1)
        do i=Vnx-bufn,Vnx+1
          xb=(xPr(i)-xs)/(xe-xs)
          DF = DampF(xb)
          V(i,j,k) = p+DF*(V(i,j,k)-p)
        enddo
      enddo
    enddo

    !$hmppcg grid blocksize 512x1
    !$hmppcg gridify (k,j) private(p,xb,DF)
    do k=1,Wnz
      do j=1,Wny
        p = 0
        do i=2*Wnx/3,Wnx-4
          p = p+W(i,j,k)
        enddo
        p = p/((Wnx-4-2*Wnx/3+1))
        do i=Wnx-bufn,Wnx+1
          xb=(xPr(i)-xs)/(xe-xs)
          DF = DampF(xb)
          W(i,j,k) = p+DF*(W(i,j,k)-p)
        enddo
      enddo
    enddo

    if (buoyancy==1) then
      !$hmppcg grid blocksize 512x1
      !$hmppcg gridify (k,j) private(p,xb,DF)
      do k=1,Prnz
        do j=1,Prny
          p = 0
          do i=2*Prnx/3,Prnx-4
            p = p+temperature(i,j,k)
          enddo
          p = p/(Prnx-4-2*Prnx/3+1)
          do i=Prnx-bufn,Prnx+1
            xb=(xPr(i)-xs)/(xe-xs)
            DF = DampF(xb)
            temperature(i,j,k) = p+DF*(temperature(i,j,k)-p)
          enddo
        enddo
      enddo
    endif

  endsubroutine AttenuateOut



  pure function DampF(x)
  implicit none
#ifdef __HMPP
   integer, parameter:: KND = 4,TIM = 4
#endif
  real(KND) DampF
  real(KND),intent(in)::x
  intrinsic exp

  if (x<=0) then
    DampF = 1
  elseif (x>=1) then
    DampF = 0
  else
   DampF=(1-0.1_KND*x**2)*(1-(1-exp(10._KND*x**2))/(1-exp(10._KND)))
  endif
  endfunction Dampf


  !$hmpp <tsteps> NullInterior codelet
  subroutine NullInterior(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                          nUnull,nVnull,nWnull,Unull,Vnull,Wnull,U,V,W)
  implicit none
#ifdef __HMPP
   integer, parameter:: KND = 4,TIM = 4
#endif
  integer,intent(in)      :: Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,nUnull,nVnull,nWnull
  integer,dimension(3,nUnull),intent(in)      :: Unull
  integer,dimension(3,nVnull),intent(in)      :: Vnull
  integer,dimension(3,nWnull),intent(in)      :: Wnull
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(inout) :: U
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(inout) :: V
  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(inout) :: W
  integer i


    do i=1,nUnull
      U(Unull(1,i),Unull(2,i),Unull(3,i)) = 0
    enddo

    do i=1,nVnull
      V(Vnull(1,i),Vnull(2,i),Vnull(3,i)) = 0
    enddo

    do i=1,nWnull
      W(Wnull(1,i),Wnull(2,i),Wnull(3,i)) = 0
    enddo

  endsubroutine NullInterior


  subroutine SubgridStresses(U,V,W,Pr)
  use Geometric, only: ScalFlIBPoints, TIBPoint_Viscosity

  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  integer i

     if (sgstype==SmagorinskyModel) then
                       call Smag(U,V,W)
     elseif (sgstype==VremanModel) then
                       if (GPU>0.and.gridtype==uniformgrid.and. Prnx*Prny*Prnz > 50) then
                           !$hmpp <tsteps> Vreman callsite, args[*].noupdate = true
                           call Vreman_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxmin,dymin,dzmin,dt,Re,U,V,W,Visc)
                           !$hmpp <tsteps> delegatedstore, args[Vreman::Visc]
                       else
                           call Vreman(U,V,W)
                       endif
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
























 !GPU codelets

!$hmpp <tsteps> group, target = CUDA

  !$hmpp <tsteps> mapbyname, Prnx,Prny,Prnz
  !$hmpp <tsteps> mapbyname, Btype,sideU,Re
  !$hmpp <tsteps> mapbyname, dt
  !$hmpp <tsteps> mapbyname, dxPr,dyPr,dzPr,dxU,dyV,dzW
  !$hmpp <tsteps> mapbyname, xPr,zPr,xU,zW
  !$hmpp <tsteps> mapbyname, Pr,Visc

  !$hmpp <tsteps> mapbyname, buoyancy

  !$hmpp <tsteps> map, args[*::Unx,BoundU::nx]
  !$hmpp <tsteps> map, args[*::Vnx,BoundV::nx]
  !$hmpp <tsteps> map, args[*::Wnx,BoundW::nx]

  !$hmpp <tsteps> map, args[*::Uny,BoundU::ny]
  !$hmpp <tsteps> map, args[*::Vny,BoundV::ny]
  !$hmpp <tsteps> map, args[*::Wny,BoundW::ny]

  !$hmpp <tsteps> map, args[*::Unz,BoundU::nz]
  !$hmpp <tsteps> map, args[*::Vnz,BoundV::nz]
  !$hmpp <tsteps> map, args[*::Wnz,BoundW::nz]

  !$hmpp <tsteps> map, args[UnifRedBlack::Uin,BoundU::Uin]
  !$hmpp <tsteps> map, args[UnifRedBlack::Vin,BoundV::Uin]
  !$hmpp <tsteps> map, args[UnifRedBlack::Win,BoundW::Uin]

  !$hmpp <tsteps> map, args[Vreman::dx,*::dxmin]
  !$hmpp <tsteps> map, args[Vreman::dy,*::dymin]
  !$hmpp <tsteps> map, args[Vreman::dz,*::dzmin]

  !$hmpp <tsteps> map, args[PressureGrad::coef,ForwEul::coef,UnifRedBlack::coef]

 !U,V,W
 !$hmpp <tsteps> map, args[Vreman::U,Convection::U,ForwEul::U,UnifRedBlack::U,TimeStepEul::U,&
  !$hmpp <tsteps>   AttenuateOut2::U,NullInterior::U]
 !$hmpp <tsteps> map, args[Vreman::V,Convection::V,ForwEul::V,UnifRedBlack::V,TimeStepEul::V,&
  !$hmpp  <tsteps>  AttenuateOut2::V,NullInterior::V]
 !$hmpp <tsteps> map, args[Vreman::W,Convection::W,ForwEul::W,UnifRedBlack::W,TimeStepEul::W,&
  !$hmpp  <tsteps>  AttenuateOut2::W,NullInterior::W]

 !U2,V2,W2
 !$hmpp <tsteps> map, args[Convection::U2,PressureGrad::U,ForwEul::U2,UnifRedBlack::U2,AttenuateTop::U,AttenuateOut::U]
 !$hmpp <tsteps> map, args[Convection::V2,PressureGrad::V,ForwEul::V2,UnifRedBlack::V2,AttenuateTop::V,AttenuateOut::V]
 !$hmpp <tsteps> map, args[Convection::W2,PressureGrad::W,ForwEul::W2,UnifRedBlack::W2,AttenuateTop::W,AttenuateOut::W]

 !U3,V3,W3 on GPU device mapped also to Ustar,Vstar,Wstar
 !$hmpp <tsteps> map, args[ForwEul::U3,UnifRedBlack::U3]
 !$hmpp <tsteps> map, args[ForwEul::V3,UnifRedBlack::V3]
 !$hmpp <tsteps> map, args[ForwEul::W3,UnifRedBlack::W3]

 !$hmpp <tsteps> map, args[Convection::temperature,AttenuateTop::temperature,AttenuateOut::temperature,AttenuateOut2::temperature]

#include "boundaries_GPU.f90"
#include "cds_GPU.f90"
#include "smagorinsky_GPU.f90"
#include "tsteps_GPU.f90"




end module TSTEPS

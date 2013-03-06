module SCALARS
 use PARAMETERS
 use ArrayUtilities
 use WALLMODELS
 use ImmersedBoundary, only: TIBPoint, ScalFlIBPoints,TIBPoint_ScalFlSource
 use LIMITERS, only: Limiter, limparam
 use BOUNDARIES
 use TILING, only: tilenx, tileny, tilenz
#ifdef __HMPP
 use HMPP_CODELETS
 use HMPP_SCALARS
#endif

implicit none
  private
  public ScalarRK3, constPrt, Rig, partdiam, partrho, percdistrib,&
     Bound_Visc, BoundTemperature, Bound_PassScalar, ComputeTDiff,&
     TTemperatureProfileSection, TTemperatureProfile, TemperatureProfile,&
     InitTemperatureProfile, InitTemperature,&
     SubsidenceProfile, SubsidenceGradient, InitSubsidenceProfile


  real(knd),dimension(:),allocatable :: partdiam,partrho,percdistrib !diameter of particles <=0 for gas

  real(knd),dimension(:),allocatable :: SubsidenceProfile
  real(knd) :: SubsidenceGradient = 0

  type TTemperatureProfileSection
    real(knd) :: top, height
    real(knd) :: jump
    real(knd) :: gradient
  end type TTemperatureProfileSection


  type TTemperatureProfile
    type(TTemperatureProfileSection), allocatable :: sections(:)
    integer   :: randomize
    real(knd) :: randomizeTop
    real(knd) :: randomizeAmplitude
  end type TTemperatureProfile

  type(TTemperatureProfile) :: TemperatureProfile

  real(knd),parameter :: constPrt = 0.6 !constant value of Prt, which may be refined further in this module

contains


  subroutine ScalarRK3(U,V,W,Temperature,Scalar,RKstage,fluxprofile)
    use RK3
    real(knd),intent(in)    :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(knd),intent(inout) :: Temperature(-1:,-1:,-1:),Scalar(-1:,-1:,-1:,-1:)
    real(knd),intent(out)   :: fluxprofile(:)
    integer,intent(in)      :: RKstage

    real(knd),dimension(:,:,:,:),allocatable,save ::Scalar_adv,Scalar_2
    real(knd),dimension(:,:,:),allocatable,save   ::Temperature_adv,Temperature2

    integer :: i
    integer,save :: called = 0
    logical,save :: released=.false.


    if (called==0) then

      called = 1

      allocate(Scalar_adv(lbound(Scalar,1):ubound(Scalar,1),lbound(Scalar,2):ubound(Scalar,2),&
           lbound(Scalar,3):ubound(Scalar,3),lbound(Scalar,4):ubound(Scalar,4)))
      allocate(Scalar_2(lbound(Scalar,1):ubound(Scalar,1),lbound(Scalar,2):ubound(Scalar,2),&
           lbound(Scalar,3):ubound(Scalar,3),lbound(Scalar,4):ubound(Scalar,4)))

      allocate(Temperature_adv(lbound(Temperature,1):ubound(Temperature,1),lbound(Temperature,2):ubound(Temperature,2),&
           lbound(Temperature,3):ubound(Temperature,3)))
      allocate(Temperature2(lbound(Temperature,1):ubound(Temperature,1),lbound(Temperature,2):ubound(Temperature,2),&
           lbound(Temperature,3):ubound(Temperature,3)))
    endif


!     if (num_of_scalars>0.and..not.released) call Release(Scalar,released)


    if (num_of_scalars>0) then

      !$omp parallel workshare
      Scalar_2=0
      !$omp end parallel workshare

      if (RKstage>1) then
        !$omp parallel workshare
        Scalar_2 = Scalar_2+Scalar_adv*RK_rho(RKstage)
        !$omp end parallel workshare
      endif

      do i=1,num_of_scalars

        call Bound_PassScalar(Scalar(:,:,:,i))

        call AdvScalar(Scalar_adv(:,:,:,i),Scalar(:,:,:,i),U,V,W,2,1._knd,SubsidenceProfile)

      enddo

      !$omp parallel workshare
      Scalar_2 = Scalar_2+Scalar_adv*RK_beta(RKstage)
      !$omp end parallel workshare

      if (scalsourcetype==pointsource) then

        do i=1,num_of_scalars
          Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i) = &
                     Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i)+&
                     percdistrib(i)*(RK_rho(RKstage) + &
                     RK_beta(RKstage))*dt*totalscalsource / &
                         (dxPr(scalsrci(i))*dyPr(scalsrcj(i))*dzPr(scalsrck(i)))
        enddo

      elseif (scalsourcetype==volumesource) then

        do i=1,num_of_scalars
          call VolScalSource(Scalar_2,RK_rho(RKstage)+RK_beta(RKstage))
        enddo

      endif

      !$omp parallel
      !$omp workshare
      Scalar = Scalar+Scalar_2
      !$omp end workshare
      !$omp end parallel

      do i=1,num_of_scalars
         call DiffScalar(Scalar_2(:,:,:,i),Scalar(:,:,:,i),2,2._knd*RK_alpha(RKstage))
      enddo

      if (computedeposition>0) call Deposition(Scalar_2,2._knd*RK_alpha(RKstage))

      if (computegravsettling>0) call GravSettling(Scalar_2,2._knd*RK_alpha(RKstage))

      !$omp parallel
      !$omp workshare
      Scalar = Scalar_2
      !$omp end workshare
      !$omp end parallel

    endif


    if (enable_buoyancy>0) then

#ifdef __HMPP

      call HMPP_RKstage_Temperature(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                            dxmin,dymin,dzmin,maxCNiter,epsCN,Re,limparam,&
                            TBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr,&
                            TDiff,Temperature2,Temperature,Temperature_adv,U,V,W,&
                            SubsidenceProfile,fluxProfile,dt,RKstage,RK_alpha,RK_beta,RK_rho)

#else
      call BoundTemperature(temperature)

      if (RKstage>1) then

        call assign(temperature2, temperature_adv)
        call multiply(temperature2, RK_rho(RKstage))

      else

        call set(temperature2,0._knd)

      endif

! #ifdef __HMPP
!       call HMPP_KappaTemperature(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
!                             dxmin,dymin,dzmin,limparam,&
!                             temperature_adv,temperature,U,V,W,&
!                             1._knd,dt,SubsidenceProfile,fluxProfile)
! #else
      call AdvScalar(temperature_adv,temperature,U,V,W,1,1._knd,SubsidenceProfile,fluxProfile)
! #endif

      call add_multiplied(temperature2, temperature_adv, RK_beta(RKstage))

      call add(temperature, temperature2)

      call BoundTemperature(temperature)

! #ifdef __HMPP
!       call HMPP_DiffTemperature(Prnx,Prny,Prnz,dxmin,dymin,dzmin,maxCNiter,epsCN,Re,&
!                              TBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr,TDiff,&
!                              temperature2,temperature,2._knd*RK_alpha(RKstage),dt)
! #else
      call DiffScalar(temperature2,temperature,1,2._knd*RK_alpha(RKstage))
! #endif

      call assign(temperature,temperature2)
      
#endif
    endif

  end subroutine ScalarRK3




  subroutine ADVSCALAR(SCAL2,SCAL,U,V,W,sctype,coef,SubsidenceProfile,fluxProfile)
  real(knd),intent(out)      :: Scal2(-1:,-1:,-1:)
  real(knd),intent(in)       :: Scal(-1:,-1:,-1:)
  real(knd),intent(in)       :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in)         :: sctype
  real(knd),intent(in)       :: SubsidenceProfile(0:)
  real(knd),intent(out),optional :: fluxProfile(0:)
  real(knd),allocatable,save :: SubsidenceProfileLoc(:),fluxProfileLoc(:)
  integer,save :: called=0

    if (called==0) then
      allocate(SubsidenceProfileLoc(0:Prnz))
      allocate(fluxProfileLoc(0:Prnz))
      called = 1
    end if

    if (size(SubsidenceProfile)==size(SubsidenceProfileLoc)) then
      SubsidenceProfileLoc = SubsidenceProfile
    else
      SubsidenceProfileLoc = 0
    endif

    if (gridtype==uniformgrid) then
      call KAPPASCALARUG(SCAL2,SCAL,U,V,W,coef,SubsidenceProfileLoc,fluxProfileLoc)
    else
      call KAPPASCALARGG(SCAL2,SCAL,U,V,W,coef,SubsidenceProfileLoc,fluxProfileLoc)
    endif

    if (present(fluxProfile).and.(size(fluxProfile)==size(fluxProfileLoc)) )&
      fluxProfile(0:Prnz) = fluxProfileLoc

  endsubroutine ADVSCALAR

  subroutine CDSSCALAR(SCAL2,SCAL,U,V,W,coef)
    real(knd),intent(out) :: Scal2(-1:,-1:,-1:)
    real(knd),intent(in)  :: Scal(-1:,-1:,-1:)
    real(knd),intent(in)  :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
    integer nx,ny,nz,i,j,k
    real(knd) Ax,Ay,Az

        nx = Prnx
        ny = Prny
        nz = Prnz

        SCAL2 = 0

        Ax = 0.5_knd*coef*dt/dxmin
        Ay = 0.5_knd*coef*dt/dymin
        Az = 0.5_knd*coef*dt/dzmin


        do k = 1,nz
            do j = 1,ny
                do i = 1,nx
                    SCAL2(i,j,k) = SCAL2(i,j,k)- ((Az*(SCAL(i,j,k+1)+SCAL(i,j,k))*(W(i,j,k))&
                    -Az*(SCAL(i,j,k)+SCAL(i,j,k-1))*(W(i,j,k-1)))&
                    +(Ay*(SCAL(i,j+1,k)+SCAL(i,j,k))*(V(i,j,k))&
                    -Ay*(SCAL(i,j,k)+SCAL(i,j-1,k))*(V(i,j-1,k)))&
                    +(Ax*(SCAL(i+1,j,k)+SCAL(i,j,k))*(U(i,j,k))&
                    -Ax*(SCAL(i,j,k)+SCAL(i-1,j,k))*(U(i-1,j,k))))
                enddo
            enddo
        enddo
  end subroutine CDSSCALAR



  subroutine KAPPASCALARUG(SCAL2,SCAL,U,V,W,coef,SubsidenceProfile,fluxProfile) !Kappa scheme with flux limiter
  real(knd),intent(out) :: Scal2(-1:,-1:,-1:) !Hunsdorfer et al. 1995, JCP
  real(knd),intent(in)  :: Scal(-1:,-1:,-1:) !Hunsdorfer et al. 1995, JCP
  real(knd),intent(in)  :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  real(knd),intent(in)  :: SubsidenceProfile(0:)
  real(knd),intent(out) :: fluxProfile(0:)
  integer i,j,k,l
  real(knd) A,Ax,Ay,Az              !Auxiliary variables to store muliplication constants for efficiency
  real(knd) vel,SL,SR,FLUX
  real(knd),allocatable,save :: SLOPE(:,:,:)
  real(knd),parameter ::eps = 1e-8
  real(knd) :: FluxLimiter,r
  integer, save :: called = 0

  FluxLimiter(r)=max(0._knd,min(2._knd*r,min(limparam,(1+2._knd*r)/3._knd)))

  if (called==0) then
    allocate(SLOPE(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    called = 1
  end if

  A = coef*dt
  Ax = coef*dt/dxmin
  Ay = coef*dt/dymin
  Az = coef*dt/dzmin


  call set(SCAL2,0._knd)
  call set(SLOPE,0._knd)

  !$omp parallel private(i,j,k,l,vel,SL,SR,FLUX) shared(SLOPE,SCAL,SCAL2,fluxProfile)
  !$omp do
  do k = 1,Prnz
   do j = 1,Prny
    do i = 0,Prnx
     if (U(i,j,k)>0) then
      SR = (SCAL(i+1,j,k)-SCAL(i,j,k))
      SL = (SCAL(i,j,k)-SCAL(i-1,j,k))
     else
      SR = (SCAL(i,j,k)-SCAL(i+1,j,k))
      SL = (SCAL(i+1,j,k)-SCAL(i+2,j,k))
     endif
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    enddo
   enddo
  enddo
  !$omp end do

  !$omp do
  do k = 1,Prnz
   do j = 1,Prny
    do i = 0,Prnx
     if (U(i,j,k)>0) then
      FLUX = U(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i-1,j,k))*SLOPE(i,j,k)/2._knd)
     else
      FLUX = U(i,j,k)*(SCAL(i+1,j,k)+(SCAL(i+1,j,k)-SCAL(i+2,j,k))*SLOPE(i,j,k)/2._knd)
     endif
     SCAL2(i,j,k) = SCAL2(i,j,k)-Ax*FLUX
     SCAL2(i+1,j,k) = SCAL2(i+1,j,k)+Ax*FLUX
    enddo
   enddo
  enddo
  !$omp end do nowait
  !$omp end parallel


  call set(SLOPE,0._knd)

  !$omp parallel private(i,j,k,l,vel,SL,SR,FLUX) shared(SLOPE,SCAL,SCAL2,fluxProfile)
  !$omp do
  do k = 1,Prnz
   do j = 0,Prny
    do i = 1,Prnx
     if (V(i,j,k)>0) then
      SR = (SCAL(i,j+1,k)-SCAL(i,j,k))
      SL = (SCAL(i,j,k)-SCAL(i,j-1,k))
     else
      SR = (SCAL(i,j,k)-SCAL(i,j+1,k))
      SL = (SCAL(i,j+1,k)-SCAL(i,j+2,k))
     endif
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    enddo
   enddo
  enddo
  !$omp end do


  !$omp do
  do k = 1,Prnz
   do j = 0,Prny
    do i = 1,Prnx
     if (V(i,j,k)>0) then
      FLUX = V(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j-1,k))*SLOPE(i,j,k)/2._knd)
     else
      FLUX = V(i,j,k)*(SCAL(i,j+1,k)+(SCAL(i,j+1,k)-SCAL(i,j+2,k))*SLOPE(i,j,k)/2._knd)
     endif

     SCAL2(i,j,k) = SCAL2(i,j,k)-Ay*FLUX
     SCAL2(i,j+1,k) = SCAL2(i,j+1,k)+Ay*FLUX
    enddo
   enddo
  enddo
  !$omp end do nowait
  !$omp end parallel


  call set(SLOPE ,0._knd)

  !$omp parallel private(i,j,k,l,vel,SL,SR,FLUX) shared(SLOPE,SCAL,SCAL2,fluxProfile)
  !$omp do
  do k = 0,Prnz
   do j = 1,Prny
    do i = 1,Prnx
     if (W(i,j,k)>0) then
      SR = (SCAL(i,j,k+1)-SCAL(i,j,k))
      SL = (SCAL(i,j,k)-SCAL(i,j,k-1))
     else
      SR = (SCAL(i,j,k)-SCAL(i,j,k+1))
      SL = (SCAL(i,j,k+1)-SCAL(i,j,k+2))
     endif
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    enddo
   enddo
  enddo
  !$omp end do nowait
  !$omp end parallel

  call set(fluxprofile,0._knd)

  !$omp parallel private(i,j,k,l,vel,SL,SR,FLUX) shared(SLOPE,SCAL,SCAL2,fluxProfile)
  do l=0,1  !odd-even separation to avoid a race condition
    !$omp do reduction(+:fluxprofile)
    do k = 0+l,Prnz,2
     do j = 1,Prny
      do i = 1,Prnx

       vel = W(i,j,k) - SubsidenceProfile(k)

       if (vel>0) then
        FLUX = vel*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j,k-1))*SLOPE(i,j,k)/2._knd)
       else
        FLUX = vel*(SCAL(i,j,k+1)+(SCAL(i,j,k+1)-SCAL(i,j,k+2))*SLOPE(i,j,k)/2._knd)
       endif

       if (abs(vel)>=1e-6) fluxprofile(k) = fluxprofile(k) + FLUX/vel*W(i,j,k)

       SCAL2(i,j,k) = SCAL2(i,j,k) - Az*FLUX
       SCAL2(i,j,k+1) = SCAL2(i,j,k+1) + Az*FLUX
      enddo
     enddo
    enddo
    !$omp end do
  enddo

  !$omp workshare
  fluxprofile = fluxprofile / (Prnx*Prny)
  !$omp end workshare

  !$omp end parallel

  endsubroutine KAPPASCALARUG



  subroutine KAPPASCALARGG(SCAL2,SCAL,U,V,W,coef,SubsidenceProfile,fluxProfile) !Kappa scheme with flux limiter
  real(knd),intent(out) :: Scal2(-1:,-1:,-1:)                                   !Hunsdorfer et al. 1995, JCP
  real(knd),intent(in)  :: Scal(-1:,-1:,-1:)
  real(knd),intent(in)  :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  real(knd),intent(in)  :: SubsidenceProfile(0:)
  real(knd),intent(out) :: fluxProfile(0:)
  integer i,j,k,l
  real(knd) A                       !Auxiliary variables to store muliplication constants for efficiency
  real(knd) vel,SL,SR,FLUX
  real(knd),allocatable,save :: SLOPE(:,:,:)
  real(knd),parameter::eps = 1e-8
  real(knd) :: FluxLimiter,r
  integer,save :: called = 0

  FluxLimiter(r)=max(0._knd,min(2._knd*r,min(limparam,(1+2._knd*r)/3._knd)))

  if (called==0) then
    allocate(SLOPE(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    called = 1
  end if

  A = coef*dt

  !$omp parallel private(i,j,k,vel,SL,SR,FLUX)

  !$omp workshare
  SCAL2 = 0
  SLOPE = 0
  !$omp end workshare
  !$omp do
  do k = 1,Prnz
   do j = 1,Prny
    do i = 0,Prnx
     if (U(i,j,k)>0) then
      SR = (SCAL(i+1,j,k)-SCAL(i,j,k))/dxU(i)
      SL = (SCAL(i,j,k)-SCAL(i-1,j,k))/dxU(i-1)
     else
      SR = (SCAL(i,j,k)-SCAL(i+1,j,k))/dxU(i)
      SL = (SCAL(i+1,j,k)-SCAL(i+2,j,k))/dxU(i-1)
     endif
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    enddo
   enddo
  enddo
  !$omp end do


  !$omp do
  do k = 1,Prnz
   do j = 1,Prny
    do i = 0,Prnx
     if (U(i,j,k)>0) then
      FLUX = U(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i-1,j,k))*SLOPE(i,j,k)/2._knd)
     else
      FLUX = U(i,j,k)*(SCAL(i+1,j,k)+(SCAL(i+1,j,k)-SCAL(i+2,j,k))*SLOPE(i,j,k)/2._knd)
     endif
     SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dxPr(i)
     SCAL2(i+1,j,k) = SCAL2(i+1,j,k)+A*FLUX/dxPr(i+1)
    enddo
   enddo
  enddo
  !$omp end do nowait


  !$omp workshare
  SLOPE = 0
  !$omp end workshare
  !$omp do
  do k = 1,Prnz
   do j = 0,Prny
    do i = 1,Prnx
     if (V(i,j,k)>0) then
      SR = (SCAL(i,j+1,k)-SCAL(i,j,k))/dyV(j)
      SL = (SCAL(i,j,k)-SCAL(i,j-1,k))/dyV(j-1)
     else
      SR = (SCAL(i,j,k)-SCAL(i,j+1,k))/dyV(j)
      SL = (SCAL(i,j+1,k)-SCAL(i,j+2,k))/dyV(j-1)
     endif
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    enddo
   enddo
  enddo
  !$omp end do


  !$omp do
  do k = 1,Prnz
   do j = 0,Prny
    do i = 1,Prnx
     if (V(i,j,k)>0) then
      FLUX = V(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j-1,k))*SLOPE(i,j,k)/2._knd)
     else
      FLUX = V(i,j,k)*(SCAL(i,j+1,k)+(SCAL(i,j+1,k)-SCAL(i,j+2,k))*SLOPE(i,j,k)/2._knd)
     endif

     SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dyPr(j)
     SCAL2(i,j+1,k) = SCAL2(i,j+1,k)+A*FLUX/dyPr(j+1)
    enddo
   enddo
  enddo
  !$omp end do nowait


  !$omp workshare
  SLOPE = 0
  !$omp end workshare
  !$omp do
  do k = 0,Prnz
   do j = 1,Prny
    do i = 1,Prnx
     if (W(i,j,k)>0) then
      SR = (SCAL(i,j,k+1)-SCAL(i,j,k))/dzW(k)
      SL = (SCAL(i,j,k)-SCAL(i,j,k-1))/dzW(k-1)
     else
      SR = (SCAL(i,j,k)-SCAL(i,j,k+1))/dzW(k)
      SL = (SCAL(i,j,k+1)-SCAL(i,j,k+2))/dzW(k-1)
     endif
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    enddo
   enddo
  enddo
  !$omp end do nowait

  !$omp workshare
  fluxprofile = 0
  !$omp end workshare

  do l=0,1  !odd-even separation to avoid a race condition
    !$omp do reduction(+:fluxprofile)
    do j = 1,Prny  !loop order due to avoid race condition
     do k = 0,Prnz
      do i = 1,Prnx

       vel = W(i,j,k) - SubsidenceProfile(k)

       if (vel>0) then
        FLUX = vel*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j,k-1))*SLOPE(i,j,k)/2._knd)
       else
        FLUX = vel*(SCAL(i,j,k+1)+(SCAL(i,j,k+1)-SCAL(i,j,k+2))*SLOPE(i,j,k)/2._knd)
       endif

       if (abs(vel)>=1e-6) fluxprofile(k) = fluxprofile(k+1) + FLUX/vel*W(i,j,k)

       SCAL2(i,j,k) = SCAL2(i,j,k) - A*FLUX/dzPr(k)
       SCAL2(i,j,k+1) = SCAL2(i,j,k+1) + A*FLUX/dzPr(k+1)
      enddo
     enddo
    enddo
    !$omp end do
  enddo

  !$omp workshare
  fluxProfile = fluxProfile/(Prnx * Prny)
  !$omp end workshare

  !$omp end parallel

  endsubroutine KAPPASCALARGG




  subroutine DIFFSCALAR(SCAL2,SCAL,sctype,coef)
  real(knd),intent(inout) ::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(knd),intent(in) :: coef
  integer,intent(in) :: sctype
  real(knd),allocatable,save :: Scal3(:,:,:)
  real(knd),allocatable,save :: Ap(:,:,:)
  integer nx,ny,nz,i,j,k,bi,bj,bk,l,xi,yj,zk
  real(knd) p,S
  real(knd) A,Ax,Ay,Az
  integer,parameter :: narr = 3, narr2 = 5
  integer,save :: called = 0

  if (called==0) then
    allocate(Scal3(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    allocate(Ap(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    called = 1
  end if

  nx = Prnx
  ny = Prny
  nz = Prnz


  if (Re>0) then
   if (sctype==1) then
    call BoundTemperature(SCAL)
   else
    call BOUND_PASSSCALAR(SCAL)
   endif

   A = dt*coef
   Ax = 1._knd/(dxmin**2)
   Ay = 1._knd/(dymin**2)
   Az = 1._knd/(dzmin**2)

   !$omp parallel private (i,j,k,bi,bj,bk)

   !initital value using forward Euler
   if (gridtype==uniformgrid) then
     !$omp do
     do bk = 1,Prnz,tilenz(narr)
      do bj = 1,Prny,tileny(narr)
       do bi = 1,Prnx,tilenx(narr)
        do k = bk,min(bk+tilenz(narr)-1,Prnz)
         do j = bj,min(bj+tileny(narr)-1,Prny)
          do i = bi,min(bi+tilenx(narr)-1,Prnx)
            SCAL3(i,j,k) = ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL(i+1,j,k)-SCAL(i,j,k))-&
              (TDiff(i,j,k)+TDiff(i-1,j,k))*(SCAL(i,j,k)-SCAL(i-1,j,k)))*Ax

            SCAL3(i,j,k) = SCAL3(i,j,k) +&
              ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL(i,j+1,k)-SCAL(i,j,k))-&
              (TDiff(i,j,k)+TDiff(i,j-1,k))*(SCAL(i,j,k)-SCAL(i,j-1,k)))*Ay

            SCAL3(i,j,k) = SCAL3(i,j,k) +&
              ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL(i,j,k+1)-SCAL(i,j,k))-&
              (TDiff(i,j,k)+TDiff(i,j,k-1))*(SCAL(i,j,k)-SCAL(i,j,k-1)))*Az
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
     !$omp end do
   else
     !$omp do
     do bk = 1,Prnz,tilenz(narr)
      do bj = 1,Prny,tileny(narr)
       do bi = 1,Prnx,tilenx(narr)
        do k = bk,min(bk+tilenz(narr)-1,Prnz)
         do j = bj,min(bj+tileny(narr)-1,Prny)
          do i = bi,min(bi+tilenx(narr)-1,Prnx)
            SCAL3(i,j,k) = (((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL(i+1,j,k)-SCAL(i,j,k))/dxU(i)-&
              (TDiff(i,j,k)+TDiff(i-1,j,k))*(SCAL(i,j,k)-SCAL(i-1,j,k))/dxU(i-1))/(dxPr(i))+&
             ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL(i,j+1,k)-SCAL(i,j,k))/dyV(j)-&
              (TDiff(i,j,k)+TDiff(i,j-1,k))*(SCAL(i,j,k)-SCAL(i,j-1,k))/dyV(j-1))/(dyPr(j))+&
             ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL(i,j,k+1)-SCAL(i,j,k))/dzW(k)-&
              (TDiff(i,j,k)+TDiff(i,j,k-1))*(SCAL(i,j,k)-SCAL(i,j,k-1))/dzW(k-1))/(dzPr(k)))
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
     !$omp end do
   endif
   !$omp end parallel

   Ax = 1._knd/(4._knd*dxmin**2)
   Ay = 1._knd/(4._knd*dymin**2)
   Az = 1._knd/(4._knd*dzmin**2)

   !SCAL2 = SCAL + SCAL3 * A
   call assign(SCAL2,SCAL)
   call add_multiplied(SCAL2,SCAL3,A)

   if (sctype==1) then
     call BoundTemperature(SCAL2)
   else
     call BOUND_PASSSCALAR(SCAL2)
   endif

   !$omp parallel private(i,j,k,bi,bj,bk,p,xi,yj,zk)
   !$omp do
   do i = 1,ubound(ScalFlIBPoints,1)
     xi = ScalFlIBPoints(i)%xi
     yj = ScalFlIBPoints(i)%yj
     zk = ScalFlIBPoints(i)%zk
     p = TIBPoint_ScalFlSource(ScalFlIBPoints(i),Scal2,sctype) * dt
     SCAL2(xi,yj,zk) = SCAL2(xi,yj,zk) + p
     SCAL3(xi,yj,zk) = SCAL3(xi,yj,zk) + p
   enddo
   !$omp end do nowait

   if (gridtype==uniformgrid) then
    !$omp do
    do bk = 1,Prnz,tilenz(narr)
     do bj = 1,Prny,tileny(narr)
      do bi = 1,Prnx,tilenx(narr)
       do k = bk,min(bk+tilenz(narr)-1,Prnz)
        do j = bj,min(bj+tileny(narr)-1,Prny)
         do i = bi,min(bi+tilenx(narr)-1,Prnx)
           Ap(i,j,k) = 1._knd/(1._knd/A+(((TDiff(i+1,j,k)+TDiff(i,j,k))+&
                                (TDiff(i,j,k)+TDiff(i-1,j,k)))*Ax+&
                                ((TDiff(i,j+1,k)+TDiff(i,j,k))+&
                                (TDiff(i,j,k)+TDiff(i,j-1,k)))*Ay+&
                                ((TDiff(i,j,k+1)+TDiff(i,j,k))+&
                                (TDiff(i,j,k)+TDiff(i,j,k-1)))*Az))
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
    !$omp end do
   else
    !$omp do
    do bk = 1,Prnz,tilenz(narr)
     do bj = 1,Prny,tileny(narr)
      do bi = 1,Prnx,tilenx(narr)
       do k = bk,min(bk+tilenz(narr)-1,Prnz)
        do j = bj,min(bj+tileny(narr)-1,Prny)
         do i = bi,min(bi+tilenx(narr)-1,Prnx)
           Ap(i,j,k) = 1._knd/(1._knd/A+(((TDiff(i+1,j,k)+TDiff(i,j,k))/dxU(i)+&
                                (TDiff(i,j,k)+TDiff(i-1,j,k))/dxU(i-1))/(4._knd*dxPr(i))+&
                                ((TDiff(i,j+1,k)+TDiff(i,j,k))/dyV(j)+&
                                (TDiff(i,j,k)+TDiff(i,j-1,k))/dyV(j-1))/(4._knd*dyPr(j))+&
                                ((TDiff(i,j,k+1)+TDiff(i,j,k))/dzW(k)+&
                                (TDiff(i,j,k)+TDiff(i,j,k-1))/dzW(k-1))/(4._knd*dzPr(k))))
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
    !$omp end do
   endif
   !$omp end parallel

   do l = 1,maxCNiter
    S = 0
    if (sctype==1) then
     call BoundTemperature(SCAL2)
    else
     call BOUND_PASSSCALAR(SCAL2)
    endif

    if (gridtype==uniformgrid) then
     !$omp parallel private(i,j,k,p) reduction(max:S)
     !$omp do
     do bk = 1,Prnz,tilenz(narr2)
      do bj = 1,Prny,tileny(narr2)
       do bi = 1,Prnx,tilenx(narr2)
        do k = bk,min(bk+tilenz(narr2)-1,Prnz)
         do j = bj,min(bj+tileny(narr2)-1,Prny)
          do i = bi+mod(bi+j+k-1,2),min(bi+tilenx(narr2)-1,Prnx),2
            p = (SCAL(i,j,k)/A)+(SCAL3(i,j,k)/4._knd+&
             ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL2(i+1,j,k))-&
              (TDiff(i,j,k)+TDiff(i-1,j,k))*(-SCAL2(i-1,j,k)))*Ax+&
             ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL2(i,j+1,k))-&
              (TDiff(i,j,k)+TDiff(i,j-1,k))*(-SCAL2(i,j-1,k)))*Ay+&
             ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL2(i,j,k+1))-&
              (TDiff(i,j,k)+TDiff(i,j,k-1))*(-SCAL2(i,j,k-1)))*Az&
             )
             p = p*Ap(i,j,k)
             S = max(S,abs(p-SCAL2(i,j,k)))
             SCAL2(i,j,k) = p
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
     !$omp enddo
     !$omp do
     do bk = 1,Prnz,tilenz(narr2)
      do bj = 1,Prny,tileny(narr2)
       do bi = 1,Prnx,tilenx(narr2)
        do k = bk,min(bk+tilenz(narr2)-1,Prnz)
         do j = bj,min(bj+tileny(narr2)-1,Prny)
          do i = bi+mod(bi+j+k,2),min(bi+tilenx(narr2)-1,Prnx),2
            p = (SCAL(i,j,k)/A)+(SCAL3(i,j,k)/4._knd+&
             ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL2(i+1,j,k))-&
              (TDiff(i,j,k)+TDiff(i-1,j,k))*(-SCAL2(i-1,j,k)))*Ax+&
             ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL2(i,j+1,k))-&
              (TDiff(i,j,k)+TDiff(i,j-1,k))*(-SCAL2(i,j-1,k)))*Ay+&
             ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL2(i,j,k+1))-&
              (TDiff(i,j,k)+TDiff(i,j,k-1))*(-SCAL2(i,j,k-1)))*Az&
             )
             p = p*Ap(i,j,k)
             S = max(S,abs(p-SCAL2(i,j,k)))
             SCAL2(i,j,k) = p
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
     !$omp enddo
     !$omp endparallel
    else
     !$omp parallel private(i,j,k,p) reduction(max:S)
     !$omp do
     do bk = 1,Prnz,tilenz(narr2)
      do bj = 1,Prny,tileny(narr2)
       do bi = 1,Prnx,tilenx(narr2)
        do k = bk,min(bk+tilenz(narr2)-1,Prnz)
         do j = bj,min(bj+tileny(narr2)-1,Prny)
          do i = bi+mod(bi+j+k-1,2),min(bi+tilenx(narr2)-1,Prnx),2
            p = (SCAL(i,j,k)/A)+(SCAL3(i,j,k)+&
             ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL2(i+1,j,k))/dxU(i)-&
              (TDiff(i,j,k)+TDiff(i-1,j,k))*(-SCAL2(i-1,j,k))/dxU(i-1))/(dxPr(i))+&
             ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL2(i,j+1,k))/dyV(j)-&
              (TDiff(i,j,k)+TDiff(i,j-1,k))*(-SCAL2(i,j-1,k))/dyV(j-1))/(dyPr(j))+&
             ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL2(i,j,k+1))/dzW(k)-&
              (TDiff(i,j,k)+TDiff(i,j,k-1))*(-SCAL2(i,j,k-1))/dzW(k-1))/(dzPr(k))&
             )/4._knd
             p = p*Ap(i,j,k)
             S = max(S,abs(p-SCAL2(i,j,k)))
             SCAL2(i,j,k) = p
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
     !$omp enddo
     !$omp do
     do bk = 1,Prnz,tilenz(narr2)
      do bj = 1,Prny,tileny(narr2)
       do bi = 1,Prnx,tilenx(narr2)
        do k = bk,min(bk+tilenz(narr2)-1,Prnz)
         do j = bj,min(bj+tileny(narr2)-1,Prny)
          do i = bi+mod(bi+j+k,2),min(bi+tilenx(narr2)-1,Prnx),2
            p = (SCAL(i,j,k)/A)+(SCAL3(i,j,k)+&
             ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL2(i+1,j,k))/dxU(i)-&
              (TDiff(i,j,k)+TDiff(i-1,j,k))*(-SCAL2(i-1,j,k))/dxU(i-1))/(dxPr(i))+&
             ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL2(i,j+1,k))/dyV(j)-&
              (TDiff(i,j,k)+TDiff(i,j-1,k))*(-SCAL2(i,j-1,k))/dyV(j-1))/(dyPr(j))+&
             ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL2(i,j,k+1))/dzW(k)-&
              (TDiff(i,j,k)+TDiff(i,j,k-1))*(-SCAL2(i,j,k-1))/dzW(k-1))/(dzPr(k))&
             )/4._knd
             p = p*Ap(i,j,k)
             S = max(S,abs(p-SCAL2(i,j,k)))
             SCAL2(i,j,k) = p
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
     !$omp enddo
     !$omp endparallel
    endif
    write (*,*) "CNscalar ",l,S
    if (S<=epsCN) exit
   enddo


  else


   call assign(SCAL2,SCAL)

   if (sctype==1) then
     call BoundTemperature(SCAL2)
   else
     call BOUND_PASSSCALAR(SCAL2)
   endif

   !$omp parallel do private(i,j,k,bi,bj,bk,p,xi,yj,zk)
   do i = 1,ubound(ScalFlIBPoints,1)
     xi = ScalFlIBPoints(i)%xi
     yj = ScalFlIBPoints(i)%yj
     zk = ScalFlIBPoints(i)%zk
     SCAL2(xi,yj,zk) = SCAL2(xi,yj,zk) + TIBPoint_ScalFlSource(ScalFlIBPoints(i),Scal2,sctype) * dt
   enddo
   !$omp end parallel do


  endif
  endsubroutine DIFFSCALAR


  subroutine BOUND_PASSSCALAR(SCAL)
  real(knd),intent(inout) :: SCAL(-1:,-1:,-1:)
  integer i,j,k,nx,ny,nz

   nx = Prnx
   ny = Prny
   nz = Prnz
    if (ScalBtype(We)==DIRICHLET) then
      do k = 1,nz
       do j = 1,ny
        SCAL(0,j,k) = sideScal(We)-(SCAL(1,j,k)-sideScal(We))
        SCAL(-1,j,k) = sideScal(We)-(SCAL(2,j,k)-sideScal(We))
       enddo
      enddo
!   if (ScalBtype(We)==DIRICHLET) then
!     do k = 1,nz
!      do j = 1,ny
!       SCAL(0,j,k) = Scalin(j,k)!-(SCAL(1,j,k)-sideScal(We))
!       SCAL(-1,j,k) = Scalin(j,k)!-(SCAL(2,j,k)-sideScal(We))
!      enddo
!     enddo
   else if (ScalBtype(We)==PERIODIC) then
    do k = 1,nz
     do j = 1,ny
      SCAL(0,j,k) = SCAL(nx,j,k)
      SCAL(-1,j,k) = SCAL(nx-1,j,k)
     enddo
    enddo
   else if (ScalBtype(We)==NEUMANN) then
    do k = 1,nz
     do j = 1,ny
      SCAL(0,j,k) = SCAL(1,j,k)-sideScal(We)*dxU(0)
      SCAL(-1,j,k) = SCAL(1,j,k)-sideScal(We)*(dxU(0)+dxU(-1))
     enddo
    enddo
   else if (ScalBtype(We)==CONSTFLUX) then
    do k = 1,nz
     do j = 1,ny
!       SCAL(0,j,k) = (SCAL(1,j,k)*(U(0,j,k)/2._knd-(TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))+sideScal(We))/&
!                                   (U(0,j,k)/2._knd+(TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))
      SCAL(0,j,k) = (SCAL(1,j,k)*((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))+sideScal(We))/&
                                  ((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))
      SCAL(-1,j,k) = SCAL(0,j,k)-(SCAL(1,j,k)-SCAL(0,j,k))
     enddo
    enddo
   else
    do k = 1,nz
     do j = 1,ny
      SCAL(0,j,k) = SCAL(1,j,k)-sideScal(We)*dxU(0)
      SCAL(-1,j,k) = SCAL(1,j,k)-sideScal(We)*(dxU(0)+dxU(-1))
     enddo
    enddo
   endif

   if (ScalBtype(Ea)==DIRICHLET) then
    do k = 1,nz
     do j = 1,ny
      SCAL(nx+1,j,k) = sideScal(Ea)-(SCAL(nx,j,k)-sideScal(Ea))
      SCAL(nx+2,j,k) = sideScal(Ea)-(SCAL(nx-1,j,k)-sideScal(Ea))
     enddo
    enddo
   else if (ScalBtype(Ea)==PERIODIC) then
    do k = 1,nz
     do j = 1,ny
      SCAL(nx+1,j,k) = SCAL(1,j,k)
      SCAL(nx+2,j,k) = SCAL(2,j,k)
     enddo
    enddo
   else if (ScalBtype(Ea)==NEUMANN) then
    do k = 1,nz
     do j = 1,ny
      SCAL(nx+1,j,k) = SCAL(nx,j,k)+sideScal(Ea)*dxU(nx+1)
      SCAL(nx+2,j,k) = SCAL(nx,j,k)+sideScal(Ea)*(dxU(nx+1)+dxU(nx+2))
     enddo
    enddo
   else if (ScalBtype(Ea)==CONSTFLUX) then
    do k = 1,nz
     do j = 1,ny
!       SCAL(nx+1,j,k) = (SCAL(nx,j,k)*(U(nx,j,k)/2._knd+(TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))+sideScal(Ea))/&
!        (-U(nx,j,k)/2._knd+(TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))
      SCAL(nx+1,j,k) = (SCAL(nx,j,k)*((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))+sideScal(Ea))/&
       ((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))
      SCAL(nx+2,j,k) = SCAL(nx+1,j,k)-(SCAL(nx,j,k)-SCAL(nx+1,j,k))
     enddo
    enddo
   else
    do k = 1,nz
     do j = 1,ny
      SCAL(nx+1,j,k) = SCAL(nx,j,k)+sideScal(Ea)*dxU(nx+1)
      SCAL(nx+2,j,k) = SCAL(nx,j,k)+sideScal(Ea)*(dxU(nx+1)+dxU(nx+2))
     enddo
    enddo
   endif


   if (ScalBtype(So)==DIRICHLET) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,0,k) = sideScal(So)-(SCAL(i,1,k)-sideScal(So))
      SCAL(i,-1,k) = sideScal(So)-(SCAL(i,2,k)-sideScal(So))
     enddo
    enddo
   else if (ScalBtype(So)==PERIODIC) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,0,k) = SCAL(i,ny,k)
      SCAL(i,-1,k) = SCAL(i,ny-1,k)
     enddo
    enddo
   else if (ScalBtype(So)==NEUMANN) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,0,k) = SCAL(i,1,k)-sideScal(So)*dyV(0)
      SCAL(i,-1,k) = SCAL(i,1,k)-sideScal(So)*(dyV(0)+dyV(-1))
     enddo
    enddo
   else if (ScalBtype(So)==CONSTFLUX) then
    do k = 1,nz
     do i=-1,nx+2
!       SCAL(i,0,k) = (SCAL(i,1,k)*(V(i,0,k)/2._knd-(TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))+sideScal(So))/&
!                                   (V(i,0,k)/2._knd+(TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))
      SCAL(i,0,k) = (SCAL(i,1,k)*((TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))+sideScal(So))/&
                                  ((TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))
      SCAL(i,-1,k) = SCAL(i,0,k)-(SCAL(i,1,k)-SCAL(i,0,k))
     enddo
    enddo
   else
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,0,k) = SCAL(i,1,k)-sideScal(So)*dyV(0)
      SCAL(i,-1,k) = SCAL(i,1,k)-sideScal(So)*(dyV(0)+dyV(-1))
     enddo
    enddo
   endif

   if (ScalBtype(No)==DIRICHLET) then
    do k = 1,nz
     do i=-1,nx+2
       SCAL(i,ny+1,k) = sideScal(No)-(SCAL(i,ny,k)-sideScal(No))
       SCAL(i,ny+2,k) = sideScal(No)-(SCAL(i,ny-1,k)-sideScal(No))
      enddo
     enddo
   else if (ScalBtype(No)==PERIODIC) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,ny+1,k) = SCAL(i,1,k)
      SCAL(i,ny+2,k) = SCAL(i,2,k)
     enddo
    enddo
   else if (ScalBtype(No)==NEUMANN) then
    do k = 1,nz
     do i=-1,nx+2
      SCAL(i,ny+1,k) = SCAL(i,ny,k)+sideScal(No)*dyV(ny+1)
      SCAL(i,ny+2,k) = SCAL(i,ny,k)+sideScal(No)*(dyV(ny+1)+dyV(ny+2))
     enddo
    enddo
   else if (ScalBtype(No)==CONSTFLUX) then
    do k = 1,nz
     do i=-1,nx+2
!       SCAL(i,ny+1,k) = (SCAL(i,ny,k)*(V(i,ny,k)/2._knd+(TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))+sideScal(No))/&
!        (-V(i,ny,k)/2._knd+(TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))
      SCAL(i,ny+1,k) = (SCAL(i,ny,k)*((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))+sideScal(No))/&
       ((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))
      SCAL(i,ny+2,k) = SCAL(i,ny+1,k)-(SCAL(i,ny,k)-SCAL(i,ny+1,k))
     enddo
    enddo
   else
    do k = 1,nz
     do i=-1,ny+2
      SCAL(i,ny+1,k) = SCAL(i,ny,k)+sideScal(No)*dyV(ny+1)
      SCAL(i,ny+2,k) = SCAL(i,ny,k)+sideScal(No)*(dyV(ny+1)+dyV(ny+2))
     enddo
    enddo
   endif



    if (ScalBtype(Bo)==DIRICHLET)  then
     do j=-1,ny+2
      do i=-1,nx+2
       SCAL(i,j,0) = sideScal(Bo)-(SCAL(i,j,1)-sideScal(Bo))
       SCAL(i,j,-1) = sideScal(Bo)-(SCAL(i,j,2)-sideScal(Bo))
      enddo
     enddo
   else if (ScalBtype(Bo)==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,0) = SCAL(i,j,nz)
      SCAL(i,j,-1) = SCAL(i,j,nz-1)
     enddo
    enddo
   else if (ScalBtype(Bo)==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,0) = SCAL(i,j,1)-sideScal(Bo)*dzW(0)
      SCAL(i,j,-1) = SCAL(i,j,1)-sideScal(Bo)*(dzW(0)+dzW(-1))
     enddo
    enddo
   else if (ScalBtype(Bo)==CONSTFLUX.or.ScalBtype(Bo)==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,0) = SCAL(i,j,1)+sideScal(Bo)*dzW(0)/((TDiff(i,j,1)+TDiff(i,j,0))/(2._knd))
      SCAL(i,j,-1) = SCAL(i,j,0)-(SCAL(i,j,1)-SCAL(i,j,0))
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,0) = SCAL(i,j,1)-sideScal(Bo)*dzW(0)
      SCAL(i,j,-1) = SCAL(i,j,1)-sideScal(Bo)*(dzW(0)+dzW(-1))
     enddo
    enddo
   endif

   if (ScalBtype(To)==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
       SCAL(i,j,nz+1) = sideScal(To)-(SCAL(i,j,nz)-sideScal(To))
       SCAL(i,j,nz+2) = sideScal(To)-(SCAL(i,j,nz-1)-sideScal(To))
      enddo
     enddo
   else if (ScalBtype(To)==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,nz+1) = SCAL(i,j,1)
      SCAL(i,j,nz+2) = SCAL(i,j,2)
     enddo
    enddo
   else if (ScalBtype(To)==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,nz+1) = SCAL(i,j,nz)+sideScal(To)*dzW(nz+1)
      SCAL(i,j,nz+2) = SCAL(i,j,nz)+sideScal(To)*(dzW(nz+1)+dzW(nz+2))
     enddo
    enddo
   else if (ScalBtype(To)==CONSTFLUX) then
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,nz+1) = SCAL(i,j,nz)-sideScal(To)*dzW(nz)/((TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2._knd))
      SCAL(i,j,nz+2) = SCAL(i,j,nz+1)-(SCAL(i,j,nz)-SCAL(i,j,nz+1))
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2
      SCAL(i,j,nz+1) = SCAL(i,j,nz)+sideScal(To)*dzW(nz+1)
      SCAL(i,j,nz+2) = SCAL(i,j,nz)+sideScal(To)*(dzW(nz+1)+dzW(nz+2))
     enddo
    enddo
   endif
  end subroutine BOUND_PASSSCALAR






  subroutine BoundTemperature(Theta)
  real(knd),intent(inout) :: theta(-1:,-1:,-1:)
  integer i,j,k,nx,ny,nz
   nx = Prnx
   ny = Prny
   nz = Prnz
!    if (TBtype(We)==DIRICHLET) then
!      do k = 1,nz
!       do j = 1,ny
!        theta(0,j,k) = sideTemp(We)!-(theta(1,j,k)-sideTemp(We))
!        theta(-1,j,k) = sideTemp(We)!-(theta(2,j,k)-sideTemp(We))
!       enddo
!      enddo
   if (TBtype(We)==DIRICHLET) then
     do k = 1,nz
      do j = 1,ny
       theta(0,j,k) = Tempin(j,k)!-(theta(1,j,k)-sideTemp(We))
       theta(-1,j,k) = Tempin(j,k)!-(theta(2,j,k)-sideTemp(We))
      enddo
     enddo
   else if (TBtype(We)==PERIODIC) then
    do k = 1,nz
     do j = 1,ny
      theta(0,j,k) = theta(nx,j,k)
      theta(-1,j,k) = theta(nx-1,j,k)
     enddo
    enddo
   else if (TBtype(We)==NEUMANN) then
    do k = 1,nz
     do j = 1,ny
      theta(0,j,k) = theta(1,j,k)-sideTemp(We)*dxU(0)
      theta(-1,j,k) = theta(1,j,k)-sideTemp(We)*(dxU(0)+dxU(-1))
     enddo
    enddo
   else if (TBtype(We)==CONSTFLUX) then
    do k = 1,nz
     do j = 1,ny
!       theta(0,j,k) = (theta(1,j,k)*(U(0,j,k)/2._knd-(TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))+sideTemp(We))/&
!                                   (U(0,j,k)/2._knd+(TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))
      theta(0,j,k) = (theta(1,j,k)*((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))+sideTemp(We))/&
                                  ((TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))
      theta(-1,j,k) = theta(0,j,k)-(theta(1,j,k)-theta(0,j,k))
     enddo
    enddo
   else
    do k = 1,nz
     do j = 1,ny
      theta(0,j,k) = theta(1,j,k)-sideTemp(We)*dxU(0)
      theta(-1,j,k) = theta(1,j,k)-sideTemp(We)*(dxU(0)+dxU(-1))
     enddo
    enddo
   endif

   if (TBtype(Ea)==DIRICHLET) then
    do k = 1,nz
     do j = 1,ny
      theta(nx+1,j,k) = sideTemp(Ea)!-(theta(nx,j,k)-sideTemp(Ea))
      theta(nx+2,j,k) = sideTemp(Ea)!-(theta(nx-1,j,k)-sideTemp(Ea))
     enddo
    enddo
   else if (TBtype(Ea)==PERIODIC) then
    do k = 1,nz
     do j = 1,ny
      theta(nx+1,j,k) = theta(1,j,k)
      theta(nx+2,j,k) = theta(2,j,k)
     enddo
    enddo
   else if (TBtype(Ea)==NEUMANN) then
    do k = 1,nz
     do j = 1,ny
      theta(nx+1,j,k) = theta(nx,j,k)+sideTemp(Ea)*dxU(nx+1)
      theta(nx+2,j,k) = theta(nx,j,k)+sideTemp(Ea)*(dxU(nx+1)+dxU(nx+2))
     enddo
    enddo
   else if (TBtype(Ea)==CONSTFLUX) then
    do k = 1,nz
     do j = 1,ny
!       theta(nx+1,j,k) = (theta(nx,j,k)*(U(nx,j,k)/2._knd+(TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))+sideTemp(Ea))/&
!        (-U(nx,j,k)/2._knd+(TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))
      theta(nx+1,j,k) = (theta(nx,j,k)*((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))+sideTemp(Ea))/&
       ((TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))
      theta(nx+2,j,k) = theta(nx+1,j,k)-(theta(nx,j,k)-theta(nx+1,j,k))
     enddo
    enddo
   else
    do k = 1,nz
     do j = 1,ny
      theta(nx+1,j,k) = theta(nx,j,k)+sideTemp(Ea)*dxU(nx+1)
      theta(nx+2,j,k) = theta(nx,j,k)+sideTemp(Ea)*(dxU(nx+1)+dxU(nx+2))
     enddo
    enddo
   endif


   if (TBtype(So)==DIRICHLET) then
    do k = 1,nz
     do i=-1,nx+2
      theta(i,0,k) = sideTemp(So)!-(theta(i,1,k)-sideTemp(So))
      theta(i,-1,k) = sideTemp(So)!-(theta(i,2,k)-sideTemp(So))
     enddo
    enddo
   else if (TBtype(So)==PERIODIC) then
    do k = 1,nz
     do i=-1,nx+2
      theta(i,0,k) = theta(i,ny,k)
      theta(i,-1,k) = theta(i,ny-1,k)
     enddo
    enddo
   else if (TBtype(So)==NEUMANN) then
    do k = 1,nz
     do i=-1,nx+2
      theta(i,0,k) = theta(i,1,k)-sideTemp(So)*dyV(0)
      theta(i,-1,k) = theta(i,1,k)-sideTemp(So)*(dyV(0)+dyV(-1))
     enddo
    enddo
   else if (TBtype(So)==CONSTFLUX) then
    do k = 1,nz
     do i=-1,nx+2
!       theta(i,0,k) = (theta(i,1,k)*(V(i,0,k)/2._knd-(TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))+sideTemp(So))/&
!                                   (V(i,0,k)/2._knd+(TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))
      theta(i,0,k) = (theta(i,1,k)*((TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))+sideTemp(So))/&
                                  ((TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))
      theta(i,-1,k) = theta(i,0,k)-(theta(i,1,k)-theta(i,0,k))
     enddo
    enddo
   else
    do k = 1,nz
     do i=-1,nx+2
      theta(i,0,k) = theta(i,1,k)-sideTemp(So)*dyV(0)
      theta(i,-1,k) = theta(i,1,k)-sideTemp(So)*(dyV(0)+dyV(-1))
     enddo
    enddo
   endif

   if (TBtype(No)==DIRICHLET) then
    do k = 1,nz
     do i=-1,nx+2
       theta(i,ny+1,k) = sideTemp(No)!-(theta(i,ny,k)-sideTemp(No))
       theta(i,ny+2,k) = sideTemp(No)!-(theta(i,ny-1,k)-sideTemp(No))
      enddo
     enddo
   else if (TBtype(No)==PERIODIC) then
    do k = 1,nz
     do i=-1,nx+2
      theta(i,ny+1,k) = theta(i,1,k)
      theta(i,ny+2,k) = theta(i,2,k)
     enddo
    enddo
   else if (TBtype(No)==NEUMANN) then
    do k = 1,nz
     do i=-1,nx+2
      theta(i,ny+1,k) = theta(i,ny,k)+sideTemp(No)*dyV(ny+1)
      theta(i,ny+2,k) = theta(i,ny,k)+sideTemp(No)*(dyV(ny+1)+dyV(ny+2))
     enddo
    enddo
   else if (TBtype(No)==CONSTFLUX) then
    do k = 1,nz
     do i=-1,nx+2
!       theta(i,ny+1,k) = (theta(i,ny,k)*(V(i,ny,k)/2._knd+(TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))+sideTemp(No))/&
!        (-V(i,ny,k)/2._knd+(TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))
      theta(i,ny+1,k) = (theta(i,ny,k)*((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))+sideTemp(No))/&
       ((TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))
      theta(i,ny+2,k) = theta(i,ny+1,k)-(theta(i,ny,k)-theta(i,ny+1,k))
     enddo
    enddo
   else
    do k = 1,nz
     do i=-1,ny+2
      theta(i,ny+1,k) = theta(i,ny,k)+sideTemp(No)*dyV(ny+1)
      theta(i,ny+2,k) = theta(i,ny,k)+sideTemp(No)*(dyV(ny+1)+dyV(ny+2))
     enddo
    enddo
   endif



!    if (TBtype(Bo)==DIRICHLET)  then
!     do j=-1,ny+2
!      do i=-1,nx+2
!       theta(i,j,0) = sideTemp(Bo)!-(theta(i,j,1)-sideTemp(Bo))
!       theta(i,j,-1) = sideTemp(Bo)!-(theta(i,j,2)-sideTemp(Bo))
!      enddo
!     enddo
!    else
   if (TBtype(Bo)==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,0) = theta(i,j,nz)
      theta(i,j,-1) = theta(i,j,nz-1)
     enddo
    enddo
   else if (TBtype(Bo)==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,0) = theta(i,j,1)-sideTemp(Bo)*dzW(0)
      theta(i,j,-1) = theta(i,j,1)-sideTemp(Bo)*(dzW(0)+dzW(-1))
     enddo
    enddo
   else if (TBtype(Bo)==CONSTFLUX.or.TBtype(Bo)==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
      if (abs(BsideTFLArr(i,j))<tiny(1._knd).and.TDiff(i,j,1)<1.1_knd/(Prandtl*Re).and.TBtype(Bo)==DIRICHLET) then
       theta(i,j,0) = BsideTArr(i,j)
      else
       theta(i,j,0) = theta(i,j,1)+BsideTFLArr(i,j)*dzW(0)/((TDiff(i,j,1)+TDiff(i,j,0))/(2._knd))
      endif
      theta(i,j,-1) = theta(i,j,0)-(theta(i,j,1)-theta(i,j,0))
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,0) = theta(i,j,1)-sideTemp(Bo)*dzW(0)
      theta(i,j,-1) = theta(i,j,1)-sideTemp(Bo)*(dzW(0)+dzW(-1))
     enddo
    enddo
   endif

   if (TBtype(To)==DIRICHLET) then
    do j=-1,ny+2
     do i=-1,nx+2
       theta(i,j,nz+1) = sideTemp(To)!-(theta(i,j,nz)-sideTemp(To))
       theta(i,j,nz+2) = sideTemp(To)!-(theta(i,j,nz-1)-sideTemp(To))
      enddo
     enddo
   else if (TBtype(To)==PERIODIC) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,nz+1) = theta(i,j,1)
      theta(i,j,nz+2) = theta(i,j,2)
     enddo
    enddo
   else if (TBtype(To)==NEUMANN) then
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,nz+1) = theta(i,j,nz)+sideTemp(To)*dzW(nz+1)
      theta(i,j,nz+2) = theta(i,j,nz)+sideTemp(To)*(dzW(nz+1)+dzW(nz+2))
     enddo
    enddo
   else if (TBtype(To)==CONSTFLUX) then
    do j=-1,ny+2
     do i=-1,nx+2
!       theta(i,j,nz+1) = (theta(i,j,nz)*(W(i,j,nz)/2._knd+(TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2*dzW(nz)))+sideTemp(To))/&
!        (-W(i,j,nz)/2._knd+(TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2*dzW(nz)))
      theta(i,j,nz+1) = theta(i,j,nz)-sideTemp(To)*dzW(nz)/((TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2._knd))
      theta(i,j,nz+2) = theta(i,j,nz+1)-(theta(i,j,nz)-theta(i,j,nz+1))
     enddo
    enddo
   else
    do j=-1,ny+2
     do i=-1,nx+2
      theta(i,j,nz+1) = theta(i,j,nz)+sideTemp(To)*dzW(nz+1)
      theta(i,j,nz+2) = theta(i,j,nz)+sideTemp(To)*(dzW(nz+1)+dzW(nz+2))
     enddo
    enddo
   endif
  endsubroutine BoundTemperature





  pure subroutine BOUND_Visc(Nu)
  real(knd),intent(inout) :: Nu(-1:,-1:,-1:)
  integer i,j,k,nx,ny,nz

  nx = Prnx
  ny = Prny
  nz = Prnz

  if (Btype(To)==PERIODIC) then
   do j = 1,ny
    do i = 1,nx
      Nu(i,j,0) = Nu(i,j,nz)
      Nu(i,j,nz+1) = Nu(i,j,1)
    enddo
   enddo
  else
   do j = 1,ny
    do i = 1,nx
      Nu(i,j,0) = Nu(i,j,1)
      Nu(i,j,nz+1) = Nu(i,j,nz)
    enddo
   enddo
  endif

  if (Btype(Ea)==PERIODIC) then
   do k = 0,nz+1
    do j = 1,ny
      Nu(0,j,k) = Nu(nx,j,k)
      Nu(nx+1,j,k) = Nu(1,j,k)
    enddo
   enddo
  else
   do k = 0,nz+1
    do j = 1,ny
      Nu(0,j,k) = Nu(1,j,k)
      Nu(nx+1,j,k) = Nu(nx,j,k)
    enddo
   enddo
  endif

  if (Btype(No)==PERIODIC) then
   do k = 0,nz+1
    do i = 0,nx+1
      Nu(i,0,k) = Nu(i,ny,k)
      Nu(i,ny+1,k) = Nu(i,1,k)
    enddo
   enddo
  else
   do k = 0,nz+1
    do i = 0,nx+1
      Nu(i,0,k) = Nu(i,1,k)
      Nu(i,ny+1,k) = Nu(i,ny,k)
    enddo
   enddo
  endif

  endsubroutine BOUND_Visc















  pure real(knd) function AirDensity(press,temp)
  real(knd),intent(in) :: press,temp
  AirDensity = press/(287.05_knd*temp)
  endfunction AirDensity


  pure real(knd) function AirDynVisc(temp)
  real(knd),intent(in) :: temp
  AirDynVisc = 1.85e-5_knd
  endfunction AirDynVisc



  pure real(knd) function CorrFactor(dp,press,temp)
  real(knd),intent(in) :: dp,press,temp
  real(DBL) l
  l = MeanFreePath(press,temp)
  CorrFactor = real(1+(2*l/dp)*(1.257_knd+0.4_knd*exp(-0.55_knd*dp/l)), knd)
  endfunction CorrFactor


  pure real(DBL) function MeanFreePath(press,temp)
  real(knd),intent(in) :: press,temp
  MeanFreePath = 2.24e-5_DBL*temp/press
  endfunction MeanFreePath


  pure real(DBL) function BrownDiffusivity(dp,press,temp)
  real(knd),intent(in) :: dp,press,temp
  real(DBL) C
  C = Corrfactor(dp,press,temp)
  BrownDiffusivity = BoltzC*temp*C/(3._knd*pi*AirDynVisc(temp)*dp)
  endfunction BrownDiffusivity


  pure real(knd) function BrownEff(dp,press,temp)
  real(knd),intent(in) :: dp,press,temp
  real(knd) Sc
  Sc = real((AirDynVisc(temp)/AirDensity(press,temp))/BrownDiffusivity(dp,press,temp), knd)
  BrownEff = Sc**(-0.54_knd)
  endfunction BrownEff


  pure real(knd) function ImpactEff(dp,press,temp,ustar,visc)
  real(knd),intent(in) :: dp,press,temp,ustar,visc
  real(knd) St
  St = SedimVelocity2(dp,press,temp)*ustar**2/visc
  ImpactEff = St**2/(400+St**2)
  endfunction ImpactEff


  pure real(knd) function SedimVelocity2(dp,press,temp)
  real(knd),intent(in) :: dp,temp,press
  real(knd) C
  C = CorrFactor(dp,press,temp)
  SedimVelocity2 = AirDensity(press,temp)*dp**2*9.81_knd*C/(18._knd*AirDynVisc(temp))
  endfunction SedimVelocity2

  pure real(DBL) function SedimVelocity(dp,rhop,press,temp)
  real(knd),intent(in) :: dp,rhop,temp,press
  real(DBL) C,us,rho,mu
  rho = AirDensity(press,temp)
  mu = AirDynVisc(temp)
  C = CorrFactor(dp,press,temp)
  us = 1+(0.42_DBL*C**2*rho*rhop/(108*mu**2))*dp**3*(1-rho/rhop)*9.81_DBL
  us = sqrt(us)
  us = (12._DBL*mu/(0.42_DBL*C*rho*dp))*(us-1._DBL)
  SedimVelocity = us
  endfunction SedimVelocity


  pure real(knd) function DyerH(zL)
  real(knd),intent(in) :: zL
  if (zL>=0) then
   DyerH = 1._knd+5._knd*zl
  else
   DyerH = 1._knd/sqrt(1._knd-16._knd*zl)
  endif
  endfunction DyerH



  pure real(knd) function AerResist(z,z0,zL,ustar,visc)
  real(knd),intent(in) :: z,z0,zL,ustar,visc
  real(knd),parameter :: yplcrit = 11.225

  if (z>z0.and.z0>0) then
   AerResist = (log(z/z0)-DyerH(zl))/(0.4_knd*ustar)
  else
   if ((z*ustar/visc)<yplcrit) then
     AerResist = (z*ustar/visc)
   else
    AerResist = (log(abs(ustar*z/visc))/0.4_knd+5.2_knd)/ustar
   endif
  endif
  endfunction AerResist

  pure real(knd) function DepositionVelocity3(dp,rhop,press,temp,z,z0,zL,ustar) !EMRAS recommended values
  real(knd),intent(in) :: dp,rhop,press,temp,z,z0,zL,ustar

   if (dp<1e-6) then
    DepositionVelocity3 = 0.5e-4
   elseif (dp<2e-6) then
    DepositionVelocity3 = 1.5e-4
   elseif (dp<10e-6) then
    DepositionVelocity3 = 10e-4
   else
    DepositionVelocity3 = 80e-4
   endif
  endfunction DepositionVelocity3



  pure real(knd) function DepositionVelocity2(dp,press,temp,z,z0,zL,ustar)
  real(knd),intent(in) :: dp,press,temp,z,z0,zL,ustar
  real(knd) SurfResist,St,visc

   visc = real(AirDynVisc(temp)/AirDensity(press,temp), knd)
   St = SedimVelocity2(dp,press,temp)*ustar**2/(visc)
   SurfResist = 1._knd/(3._knd*ustar*(BrownEff(dp,press,temp)+ImpactEff(dp,press,temp,ustar,visc)))
   DepositionVelocity2 = SedimVelocity2(dp,press,temp)+1._knd/(AerResist(z,z0,zL,ustar,visc)+SurfResist)
  endfunction DepositionVelocity2


  pure real(knd) function DepositionVelocity(dp,rhop,press,temp,z,z0,zL,ustar) !Kharchenko
  real(knd),intent(in) :: dp,rhop,press,temp,z,z0,zL,ustar
  real(DBL) visc,us,Intz,Intexp,BD,tp
  real(DBL),parameter :: zexp = 0.01

   us = SedimVelocity(dp,rhop,press,temp)
   visc = AirDynVisc(temp)/AirDensity(press,temp)
   tp = (us/9.81_knd)*ustar**2/visc
   BD = BrownDiffusivity(dp,press,temp)

   if (zl>=0) then
    Intz=-log(z/zexp+6*(zl-zexp*zl/z))/Karman
   else
    Intz = ((sqrt(1-9*zl)-1)*(sqrt(1-9*zexp*zl/z)+1))
    Intz = Intz/((sqrt(1-9*zl)+1)*(sqrt(1-9*zexp*zl/z)-1))
    Intz=-Intz/Karman
   endif

   Intexp=-367.8_knd
   Intexp = Intexp+16.4*log(visc/BD)
   Intexp = Intexp-0.73*log(100*dp)*log(1e4*BD)-0.5*(log(100*dp))**2
   Intexp = Intexp+0.13*log(0.03/z0)
   Intexp = Intexp+0.25*log(0.2/ustar)*(1-0.2*log(0.03/z0))
   Intexp = Intexp-0.03*log(tp)*log(0.03/z0)
   Intexp = Intexp-32.7*log(100*dp)
   Intexp=-exp(Intexp)

   DepositionVelocity = real(us/(1._knd-exp((us/ustar)*(Intexp+Intz))), knd)
  endfunction DepositionVelocity


  pure real(knd) function DepositionFlux(WMP,conc,partdiam,rhop)
  type(WMPoint),intent(in) :: WMP
  real(knd),intent(in) :: conc,partdiam,rhop
  real(knd) :: press,temp,depvel

   press = 101300
   temp = temperature_ref
   if ((.2_knd*WMP%distz)**2>(WMP%distx)**2+(WMP%disty)**2.and.WMP%distz>0) then
    depvel = DepositionVelocity(partdiam,rhop,press,temp,WMP%distz,WMP%z0,0._knd,WMP%ustar)
   else
    depvel = 0
   endif
   DepositionFlux = abs(depvel)*abs(conc)

  endfunction DepositionFlux

  subroutine Deposition(SCAL,coef)
  real(knd),dimension(-1:,-1:,-1:,1:),intent(inout) :: SCAL
  real(knd),intent(in) :: coef
  type(WMPoint),pointer :: WMP
  integer i
  real(knd) deptmp
   if (associated(FirstWMPoint)) then
   WMP => FirstWMPoint
   do
    if (allocated(WMP%depscalar)) then
    if (partdistrib>0) then
     do i = 1,partdistrib
      deptmp = abs(DepositionFlux(WMP,SCAL(WMP%xi,WMP%yj,WMP%zk,1)*percdistrib(i),partdiam(i),partrho(i)))&
          *coef*dt*dxPr(WMP%xi)*dyPr(WMP%yj)
      WMP%depscalar(1) = WMP%depscalar(1)+deptmp/(dxPr(WMP%xi)*dyPr(WMP%yj)*dzPr(WMP%zk))
      SCAL(WMP%xi,WMP%yj,WMP%zk,1) = SCAL(WMP%xi,WMP%yj,WMP%zk,1)-deptmp/(dxPr(WMP%xi)*dyPr(WMP%yj)*dzPr(WMP%zk))
     enddo
    else
     do i = 1,num_of_scalars
      deptmp = abs(DepositionFlux(WMP,SCAL(WMP%xi,WMP%yj,WMP%zk,i),partdiam(i),partrho(i)))&
          *coef*dt*dxPr(WMP%xi)*dyPr(WMP%yj)
      WMP%depscalar(i) = WMP%depscalar(i)+deptmp/(dxPr(WMP%xi)*dyPr(WMP%yj)*dzPr(WMP%zk))
      SCAL(WMP%xi,WMP%yj,WMP%zk,i) = SCAL(WMP%xi,WMP%yj,WMP%zk,i)-deptmp/(dxPr(WMP%xi)*dyPr(WMP%yj)*dzPr(WMP%zk))
     enddo
    endif
    endif
    if (associated(WMP%next)) then
     WMP => WMP%next
    else
     exit
    endif
   enddo
   endif
  endsubroutine Deposition


  subroutine Gravsettling(SCAL,coef)
  real(knd),dimension(-1:,-1:,-1:,1:) :: SCAL
  integer i,j,k,l
  real(knd),dimension(Prnx,Prny,Prnz) :: flux
  real(knd) :: coef,press,temp,us

  press = 101300
  temp = temperature_ref
  if (partdistrib==0) then
   do l = 1,num_of_scalars
    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       us = real(SedimVelocity(partdiam(l),partrho(l),press,temp), knd)
       flux(i,j,k) = us*SCAL(i,j,k+1,l)*coef*dt*dxPr(i)*dyPr(j)
      enddo
     enddo
    enddo
    do k = 1,Prnz-1
     do j = 1,Prny
      do i = 1,Prnx
       SCAL(i,j,k+1,l) = SCAL(i,j,k+1,l)-flux(i,j,k)/(dxPr(i)*dyPr(j)*dzPr(k+1))
       SCAL(i,j,k,l) = SCAL(i,j,k,l)+flux(i,j,k)/(dxPr(i)*dyPr(j)*dzPr(k))
      enddo
     enddo
    enddo
   enddo
  endif
  endsubroutine Gravsettling


  pure real(knd) function Rig(i,j,k,U,V,temperature)
  integer,intent(in) :: i,j,k
  real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V
  real(knd),dimension(-1:,-1:,-1:),intent(in) :: temperature
  real(knd) num,denom

  num = (grav_acc/temperature_ref)*(temperature(i,j,k+1)-temperature(i,j,k-1))/(zPr(k+1)-zPr(k-1))
  denom = ((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._knd*(zPr(k+1)-zPr(k-1))))**2
  denom = denom+((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._knd*(zPr(k+1)-zPr(k-1))))**2

  if (abs(denom)>1E-5_knd*abs(num)) then
   Rig = num/denom
  else
   Rig = 100000._knd*sign(1.0_knd,num)*sign(1.0_knd,denom)
  endif
  endfunction Rig



  subroutine ComputeTDiff(U,V,W)
  real(knd),intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer:: i,j,k
  real(knd),parameter :: Prt = constPrt !if variable, then implement as an internal or statement function due to problems with inlining

   if (Re>0) then
    !$omp parallel do private(i,j,k)
    do k=1,Prnz
      do j=1,Prny
        do i=1,Prnx
          TDiff(i,j,k) = (Visc(i,j,k)-1._knd/Re)/Prt + (1._knd/(Re*Prandtl))
        end do
      end do
    end do
    !$omp end parallel do
   else
    !$omp parallel do private(i,j,k)
    do k=1,Prnz
      do j=1,Prny
        do i=1,Prnx
          TDiff(i,j,k) = Visc(i,j,k)/Prt
        end do
      end do
    end do
    !$omp end parallel do
   endif
  end subroutine ComputeTDiff



  subroutine InitTemperatureProfile(TempIn)
    real(knd),intent(out) :: TempIn(-1:,-1:)
    integer   :: SectionToUse(-1:ubound(TempIn,2))
    integer   :: section,nSections,s
    integer   :: i,j,k
    real(knd) :: temp

    nSections = size(TemperatureProfile%Sections)

    if (nSections > 0) then

       if (size(TemperatureProfile%Sections)>0) then
         if (TemperatureProfile%Sections(1)%jump<=0) TemperatureProfile%Sections(1)%jump = temperature_ref
       endif

      TemperatureProfile%Sections(1)%height = TemperatureProfile%Sections(1)%top

      do i = 2,nSections
        if (TemperatureProfile%Sections(i)%top < TemperatureProfile%Sections(i-1)%top)&
          TemperatureProfile%Sections(i)%top = TemperatureProfile%Sections(i-1)%top

        TemperatureProfile%Sections(i)%height = TemperatureProfile%Sections(i)%top - TemperatureProfile%Sections(i-1)%top
      enddo

      if (TemperatureProfile%Sections(nSections)%top < zW(Prnz+2) )&
          TemperatureProfile%Sections(nSections)%top = zW(Prnz+2)

      section = 1

      do k = -1, Prnz+2
        do while (zPr(k) > TemperatureProfile%Sections(section)%top) !should be safe, because last top is adjusted above
          section = section + 1
        enddo
        SectionToUse(k) = section
      enddo

    else

      SectionToUse = 0

    endif

    do k=-1,Prnz+2

      s = SectionToUse(k)

      if (s==0) then

        temp = temperature_ref

      else

        temp = 0

        do i = 1, s-1
          temp = temp + TemperatureProfile%Sections(i)%jump
          temp = temp + TemperatureProfile%Sections(i)%height * TemperatureProfile%Sections(i)%gradient
        enddo

        temp = temp + TemperatureProfile%Sections(s)%jump

        if (s>1) then
          temp = temp + (zPr(k) - TemperatureProfile%Sections(s-1)%top) * TemperatureProfile%Sections(s)%gradient
        else
          temp = temp + zPr(k) * TemperatureProfile%Sections(s)%gradient
        endif

      endif

      do j=-1,Prny+2
          Tempin(j,k) = temp
      enddo

    enddo
  end subroutine InitTemperatureProfile


  subroutine InitTemperature(TempIn,Temperature)
    real(knd),intent(in) :: TempIn(-1:,-1:)
    real(knd),intent(out) :: Temperature(-1:,-1:,-1:)
    real(knd) :: p_x,p_y,p_z,p
    integer   :: i,j,k

    if (TemperatureProfile%randomize==1) then

!       p_x=0.5_knd
!       p_y=0.5_knd
!       p_z=0.5_knd

      do k=0,Prnz+1

!         if (zPr(k) <= TemperatureProfile%randomizeTop) then
!           call random_number(p)
!           p = p - 0.5
!           p_z=(p_z+p)/2
!         endif

        do j=0,Prny+1

!           if (zPr(k) <= TemperatureProfile%randomizeTop) then
!             call random_number(p)
!             p = p - 0.5
!             p_y=(p_y+p)/2
!           endif

           do i=0,Prnx+1

             if (zPr(k) <= TemperatureProfile%randomizeTop) then
               call random_number(p)
               p = p - 0.5
!               p_x=(p_x+p)/2
!               p=p_x+p_y+p_z
             else
               p=0
             endif

             Temperature(i,j,k)=Tempin(j,k) + TemperatureProfile%randomizeAmplitude * 2 * p

           enddo
        enddo
      enddo

    else

      forall(i=0:Prnx+1) Temperature(i,:,:)=Tempin(:,:)

    endif

  end subroutine InitTemperature




  subroutine InitSubsidenceProfile
    integer k

    if (SubsidenceGradient/=0) then
      allocate(SubsidenceProfile(0:Prnz))
      SubsidenceProfile = (/ (zW(k)*SubsidenceGradient, k=0,Prnz) /)
    else
      allocate(SubsidenceProfile(0))
    endif
  end subroutine InitSubsidenceProfile



  subroutine VolScalSource(Scal,coef)
  real(knd),intent(inout) :: Scal(-1:,-1:,-1:,:)
  real(knd),intent(in) :: coef

#include "customvolscalsource.f90"
  end subroutine VolScalSource




  subroutine Release(Scalar,released)
    real(knd),intent(inout) :: Scalar(-1:,-1:,-1:,-1:)
    logical,intent(inout)   :: released
    real(knd) xc,yc,xs,xf,ys,yf,zs,zf,dxp,dyp,dzp,ct,cr,xp,yp,zp,p
    integer i,j,k,xi,yj,zk,nprobx,nproby,nprobz
    ct = 7
    cr = 1.5
    xc = 3*cos((xheading-70)*pi/180.)
    yc = 3*sin((xheading-70)*pi/180.)
    xs = xc-cr
    ys = yc-cr
    zs = 0
    xf = xc+cr
    yf = yc+cr
    zf = ct
    nprobx = 100
    nproby = 100
    nprobz = 100
    dxp=(xf-xs)/nprobx
    dyp=(yf-ys)/nproby
    dzp=(zf-zs)/nprobz
    if (num_of_scalars>=4) then
     if (time>(endtime-starttime)/3._knd) then
      Scalar = 0
      p = 0
      do k=0,nprobz
       zp = zs+k*dzp
       do j=0,nproby
        yp = ys+j*dyp
        do i=0,nprobx
         xp = xs+i*dxp
          call GridCoords(xi,yj,zk,xp,yp,zp)
          if ((xp-xc)**2+(yp-yc)**2<cr**2) then
          if   (zp<ct*0.2) then
           Scalar(xi,yj,zk,:) = Scalar(xi,yj,zk,:) + percdistrib(:)*0.2
           p = p+1
          elseif (zp<ct*0.4) then
           Scalar(xi,yj,zk,:) = Scalar(xi,yj,zk,:) + percdistrib(:)*0.8
           p = p+1
          elseif (zp<ct*0.6) then
           Scalar(xi,yj,zk,:) = Scalar(xi,yj,zk,:) + percdistrib(:)*1.25
           p = p+1
          elseif (zp<ct*0.8) then
           Scalar(xi,yj,zk,:) = Scalar(xi,yj,zk,:) + percdistrib(:)*1.75
           p = p+1
          elseif (zp<=ct) then
           Scalar(xi,yj,zk,:) = Scalar(xi,yj,zk,:)  +percdistrib(:)*1.1
           p = p+1
          endif
         endif
        enddo
       enddo
      enddo
      Scalar = totalscalsource*Scalar/p
      Scalar = Scalar/(dxmin*dymin*dzmin)
      released=.true.
     endif
    endif
  endsubroutine Release

end module SCALARS

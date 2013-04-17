module Scalars
 use Parameters
 use ArrayUtilities
 use Wallmodels
 use ImmersedBoundary, only: TIBPoint, ScalFlIBPoints,TIBPoint_ScalFlSource
 use Limiters, only: Limiter, limparam
 use Boundaries
 use ScalarBoundaries
 use TILING, only: tilenx, tileny, tilenz
#ifdef __HMPP
 use HMPP_Codelets
 use HMPP_Scalars
#endif

implicit none
  private
  public ScalarRK3, constPrt, Rig, partdiam, partrho, percdistrib, &
     ComputeTDiff, &
     TemperatureProfile, MoistureProfile, &
     InitScalarProfile, InitScalar, InitHydrostaticPressure, &
     SubsidenceProfile, SubsidenceGradient, InitSubsidenceProfile


  real(knd),dimension(:),allocatable :: partdiam,partrho,percdistrib !diameter of particles <=0 for gas

  real(knd),dimension(:),allocatable :: SubsidenceProfile
  real(knd) :: SubsidenceGradient = 0

  type TScalarProfileSection
    real(knd) :: top, height
    real(knd) :: jump
    real(knd) :: gradient
  end type


  type TScalarProfile
    type(TScalarProfileSection), allocatable :: sections(:)
    integer   :: randomize
    real(knd) :: randomizeTop
    real(knd) :: randomizeAmplitude
  end type

  type(TScalarProfile) :: TemperatureProfile
  type(TScalarProfile) :: MoistureProfile

  real(knd),parameter :: constPrt = 0.6 !constant value of Prt, which may be refined further in this module

  abstract interface
    subroutine boundary_interface(array)
      use Parameters
      real(knd),intent(inout) :: array(-1:,-1:,-1:)
    end subroutine
    subroutine extra_interface
    end subroutine
  end interface

contains


  subroutine ScalarRK3(U,V,W,Temperature,Moisture,Scalar,RK_stage,fluxprofile)
    use RK3
    use VolumeSources, only: ScalarVolumeSources
    real(knd),intent(in)    :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(knd),intent(inout) :: Temperature(-1:,-1:,-1:),Moisture(-1:,-1:,-1:)
    real(knd),intent(inout) :: Scalar(-1:,-1:,-1:,1:)
    real(knd),intent(out)   :: fluxprofile(:)
    integer,  intent(in)    :: RK_stage

    real(knd),dimension(:,:,:),  allocatable,save :: Temperature_adv, Temperature2
    real(knd),dimension(:,:,:),  allocatable,save :: Moisture_adv,    Moisture2
    real(knd),dimension(:,:,:,:),allocatable,save :: Scalar_adv,      Scalar_2

    integer :: i
    integer,save :: called = 0
    logical,save :: released=.false.


    if (called==0) then

      called = 1

      allocate(Scalar_adv(lbound(Scalar,1):ubound(Scalar,1), &
                          lbound(Scalar,2):ubound(Scalar,2), &
                          lbound(Scalar,3):ubound(Scalar,3), &
                          lbound(Scalar,4):ubound(Scalar,4)))
      allocate(Scalar_2(lbound(Scalar,1):ubound(Scalar,1), &
                        lbound(Scalar,2):ubound(Scalar,2), &
                        lbound(Scalar,3):ubound(Scalar,3), &
                        lbound(Scalar,4):ubound(Scalar,4)))

      allocate(Temperature_adv(lbound(Temperature,1):ubound(Temperature,1), &
                               lbound(Temperature,2):ubound(Temperature,2), &
                               lbound(Temperature,3):ubound(Temperature,3)))
      allocate(Temperature2(lbound(Temperature,1):ubound(Temperature,1), &
                            lbound(Temperature,2):ubound(Temperature,2), &
                            lbound(Temperature,3):ubound(Temperature,3)))

      allocate(Moisture_adv(lbound(Moisture,1):ubound(Moisture,1), &
                               lbound(Moisture,2):ubound(Moisture,2), &
                               lbound(Moisture,3):ubound(Moisture,3)))
      allocate(Moisture2(lbound(Moisture,1):ubound(Moisture,1), &
                            lbound(Moisture,2):ubound(Moisture,2), &
                            lbound(Moisture,3):ubound(Moisture,3)))
    end if


!     if (num_of_scalars>0.and..not.released) call Release(Scalar,released)


    if (num_of_scalars>0) then

      !$omp parallel workshare
      Scalar_2=0
      !$omp end parallel workshare

      if (RK_stage>1) then
        !$omp parallel workshare
        Scalar_2 = Scalar_2+Scalar_adv*RK_rho(RK_stage)
        !$omp end parallel workshare
      end if

      do i=1,num_of_scalars

        call BoundScalar(Scalar(:,:,:,i))

        call AdvScalar(Scalar_adv(:,:,:,i),Scalar(:,:,:,i),U,V,W,1._knd,SubsidenceProfile)

      end do

      call ScalarVolumeSources(Scalar_adv)

      !$omp parallel workshare
      Scalar_2 = Scalar_2+Scalar_adv*RK_beta(RK_stage)
      !$omp end parallel workshare

      if (scalsourcetype==pointsource) then

        do i=1,num_of_scalars
          Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i) = &
                     Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i)+&
                     percdistrib(i)*(RK_rho(RK_stage) + &
                     RK_beta(RK_stage))*dt*totalscalsource / &
                         (dxPr(scalsrci(i))*dyPr(scalsrcj(i))*dzPr(scalsrck(i)))
        end do

      else if (scalsourcetype==volumesource) then

        do i=1,num_of_scalars
          call VolScalSource(Scalar_2,RK_rho(RK_stage)+RK_beta(RK_stage))
        end do

      end if

      !$omp parallel
      !$omp workshare
      Scalar = Scalar+Scalar_2
      !$omp end workshare
      !$omp end parallel

      do i=1,num_of_scalars
         call DiffScalar(Scalar_2(:,:,:,i),Scalar(:,:,:,i), &
                         ScalarTypePassive,BoundScalar,2._knd*RK_alpha(RK_stage))
      end do

      if (computedeposition>0) call Deposition(Scalar_2,2._knd*RK_alpha(RK_stage))

      if (computegravsettling>0) call GravSettling(Scalar_2,2._knd*RK_alpha(RK_stage))

      !$omp parallel
      !$omp workshare
      Scalar = Scalar_2
      !$omp end workshare
      !$omp end parallel

    end if


#ifdef __HMPP
    if (enable_buoyancy>0) then
      call HMPP_RK_stage_Temperature(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                            dxmin,dymin,dzmin,maxCNiter,epsCN,Re,limparam,&
                            TempBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr,&
                            TDiff,Temperature2,Temperature,Temperature_adv,U,V,W,&
                            SubsidenceProfile,fluxProfile,dt,RK_stage,RK_alpha,RK_beta,RK_rho)

    end if
#else
    if (enable_buoyancy>0) then
print *,"temp"
      call stage(Temperature,Temperature2,Temperature_adv, &
                 ScalarTypeTemperature,BoundTemperature,TemperatureExtra,fluxProfile)
        where (Prtype>0) Temperature(0:Prnx+1,0:Prny+1,0:Prnz+1) = temperature_ref
    end if

    if (enable_moisture>0) then
print *,"moist"
      call stage(Moisture,Moisture2,Moisture_adv, &
                 ScalarTypeMoisture,BoundMoisture,MoistureExtra)
    end if
#endif

    contains

      subroutine stage(Array,Array2,Array_adv, &
                       scalar_type,boundary_procedure,extra_procedure, &
                       flux_profile)
        real(knd),dimension(-1:,-1:,-1:),intent(inout) :: Array,Array2,Array_adv
        integer,intent(in) :: scalar_type
        procedure(boundary_interface) :: boundary_procedure
        procedure(extra_interface) :: extra_procedure
        real(knd),intent(out),optional :: flux_profile(0:)
        integer i,j,k


        call BoundTemperature(Array)

        if (RK_stage>1) then

          call assign(Array2, Array_adv)
          call multiply(Array2, RK_rho(RK_stage))

        else

          call set(Array2,0._knd)

        end if

        if (present(flux_profile)) then
          call AdvScalar(Array_adv,Array,U,V,W,1._knd,SubsidenceProfile,flux_profile)
        else
          call AdvScalar(Array_adv,Array,U,V,W,1._knd,SubsidenceProfile)
        end if

        !$omp parallel do private(i,j,k)
        do k=0,Prnz+1
          do j=0,Prny+1
            do i=0,Prnx+1
              if (Prtype(i,j,k)>0) Array_adv(i,j,k) = 0
            end do
          end do
        end do
        !$omp end parallel do

        call extra_procedure

        call add_multiplied(Array2, Array_adv, RK_beta(RK_stage))

        call add(Array, Array2)

        call boundary_procedure(Array)

        call DiffScalar(Array2,Array, &
                        scalar_type,boundary_procedure,2._knd*RK_alpha(RK_stage))

        call assign(Array,Array2)
      end subroutine

      subroutine TemperatureExtra
        use VolumeSources, only: TemperatureVolumeSources
        call TemperatureVolumeSources(Temperature_adv)
      end subroutine

      subroutine MoistureExtra
        use VolumeSources, only: MoistureVolumeSources
        call MoistureVolumeSources(Moisture_adv)
      end subroutine

      subroutine ScalarExtra
        use VolumeSources, only: ScalarVolumeSources
        call ScalarVolumeSources(Scalar_adv)
      end subroutine


  end subroutine ScalarRK3




  subroutine ADVSCALAR(SCAL2,SCAL,U,V,W,coef,SubsidenceProfile,fluxProfile)
  real(knd),intent(out)      :: Scal2(-1:,-1:,-1:)
  real(knd),intent(in)       :: Scal(-1:,-1:,-1:)
  real(knd),intent(in)       :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
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
    end if

    if (gridtype==uniformgrid) then
      call KAPPASCALARUG(SCAL2,SCAL,U,V,W,coef,SubsidenceProfileLoc,fluxProfileLoc)
    else
      call KAPPASCALARGG(SCAL2,SCAL,U,V,W,coef,SubsidenceProfileLoc,fluxProfileLoc)
    end if

    if (present(fluxProfile)) then
      if (size(fluxProfile)==size(fluxProfileLoc)) &
        fluxProfile(0:Prnz) = fluxProfileLoc
    end if

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
                end do
            end do
        end do
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
     end if
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    end do
   end do
  end do
  !$omp end do

  !$omp do
  do k = 1,Prnz
   do j = 1,Prny
    do i = 0,Prnx
     if (U(i,j,k)>0) then
      FLUX = U(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i-1,j,k))*SLOPE(i,j,k)/2._knd)
     else
      FLUX = U(i,j,k)*(SCAL(i+1,j,k)+(SCAL(i+1,j,k)-SCAL(i+2,j,k))*SLOPE(i,j,k)/2._knd)
     end if
     SCAL2(i,j,k) = SCAL2(i,j,k)-Ax*FLUX
     SCAL2(i+1,j,k) = SCAL2(i+1,j,k)+Ax*FLUX
    end do
   end do
  end do
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
     end if
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    end do
   end do
  end do
  !$omp end do


  !$omp do
  do k = 1,Prnz
   do j = 0,Prny
    do i = 1,Prnx
     if (V(i,j,k)>0) then
      FLUX = V(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j-1,k))*SLOPE(i,j,k)/2._knd)
     else
      FLUX = V(i,j,k)*(SCAL(i,j+1,k)+(SCAL(i,j+1,k)-SCAL(i,j+2,k))*SLOPE(i,j,k)/2._knd)
     end if

     SCAL2(i,j,k) = SCAL2(i,j,k)-Ay*FLUX
     SCAL2(i,j+1,k) = SCAL2(i,j+1,k)+Ay*FLUX
    end do
   end do
  end do
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
     end if
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    end do
   end do
  end do
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
       end if

       if (abs(vel)>=1e-6) fluxprofile(k) = fluxprofile(k) + FLUX/vel*W(i,j,k)

       SCAL2(i,j,k) = SCAL2(i,j,k) - Az*FLUX
       SCAL2(i,j,k+1) = SCAL2(i,j,k+1) + Az*FLUX
      end do
     end do
    end do
    !$omp end do
  end do

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
     end if
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    end do
   end do
  end do
  !$omp end do


  !$omp do
  do k = 1,Prnz
   do j = 1,Prny
    do i = 0,Prnx
     if (U(i,j,k)>0) then
      FLUX = U(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i-1,j,k))*SLOPE(i,j,k)/2._knd)
     else
      FLUX = U(i,j,k)*(SCAL(i+1,j,k)+(SCAL(i+1,j,k)-SCAL(i+2,j,k))*SLOPE(i,j,k)/2._knd)
     end if
     SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dxPr(i)
     SCAL2(i+1,j,k) = SCAL2(i+1,j,k)+A*FLUX/dxPr(i+1)
    end do
   end do
  end do
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
     end if
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    end do
   end do
  end do
  !$omp end do


  !$omp do
  do k = 1,Prnz
   do j = 0,Prny
    do i = 1,Prnx
     if (V(i,j,k)>0) then
      FLUX = V(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j-1,k))*SLOPE(i,j,k)/2._knd)
     else
      FLUX = V(i,j,k)*(SCAL(i,j+1,k)+(SCAL(i,j+1,k)-SCAL(i,j+2,k))*SLOPE(i,j,k)/2._knd)
     end if

     SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dyPr(j)
     SCAL2(i,j+1,k) = SCAL2(i,j+1,k)+A*FLUX/dyPr(j+1)
    end do
   end do
  end do
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
     end if
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._knd,SL))/(SL+eps*sign(1._knd,SL)))
    end do
   end do
  end do
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
       end if

       if (abs(vel)>=1e-6) fluxprofile(k) = fluxprofile(k+1) + FLUX/vel*W(i,j,k)

       SCAL2(i,j,k) = SCAL2(i,j,k) - A*FLUX/dzPr(k)
       SCAL2(i,j,k+1) = SCAL2(i,j,k+1) + A*FLUX/dzPr(k+1)
      end do
     end do
    end do
    !$omp end do
  end do

  !$omp workshare
  fluxProfile = fluxProfile/(Prnx * Prny)
  !$omp end workshare

  !$omp end parallel

  endsubroutine KAPPASCALARGG




  subroutine DIFFSCALAR(SCAL2,SCAL,sctype,boundary_procedure,coef)
  real(knd),intent(in)    :: Scal(-1:,-1:,-1:)
  real(knd),intent(inout) :: Scal2(-1:,-1:,-1:)
  real(knd),intent(in) :: coef
  integer,intent(in) :: sctype
  procedure(boundary_interface) :: boundary_procedure
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
            if (Prtype(i,j,k)<=0) then
              SCAL3(i,j,k) = ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL(i+1,j,k)-SCAL(i,j,k))-&
                (TDiff(i,j,k)+TDiff(i-1,j,k))*(SCAL(i,j,k)-SCAL(i-1,j,k)))*Ax

              SCAL3(i,j,k) = SCAL3(i,j,k) +&
                ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL(i,j+1,k)-SCAL(i,j,k))-&
                (TDiff(i,j,k)+TDiff(i,j-1,k))*(SCAL(i,j,k)-SCAL(i,j-1,k)))*Ay

              SCAL3(i,j,k) = SCAL3(i,j,k) +&
                ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL(i,j,k+1)-SCAL(i,j,k))-&
                (TDiff(i,j,k)+TDiff(i,j,k-1))*(SCAL(i,j,k)-SCAL(i,j,k-1)))*Az
            else
              SCAL3(i,j,k) = 0
            end if
          end do
         end do
        end do
       end do
      end do
     end do
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
          end do
         end do
        end do
       end do
      end do
     end do
     !$omp end do
   end if
   !$omp end parallel

   Ax = 1._knd/(4._knd*dxmin**2)
   Ay = 1._knd/(4._knd*dymin**2)
   Az = 1._knd/(4._knd*dzmin**2)

   !SCAL2 = SCAL + SCAL3 * A
   call assign(SCAL2,SCAL)
   call add_multiplied(SCAL2,SCAL3,A)

   call boundary_procedure(SCAL2)

   !$omp parallel private(i,j,k,bi,bj,bk,p,xi,yj,zk)
   !$omp do
   do i = 1,ubound(ScalFlIBPoints,1)
     xi = ScalFlIBPoints(i)%xi
     yj = ScalFlIBPoints(i)%yj
     zk = ScalFlIBPoints(i)%zk
     p = TIBPoint_ScalFlSource(ScalFlIBPoints(i),Scal2,sctype) * dt
     SCAL2(xi,yj,zk) = SCAL2(xi,yj,zk) + p
     SCAL3(xi,yj,zk) = SCAL3(xi,yj,zk) + p
   end do
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
         end do
        end do
       end do
      end do
     end do
    end do
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
         end do
        end do
       end do
      end do
     end do
    end do
    !$omp end do
   end if
   !$omp end parallel

   do l = 1,maxCNiter
    S = 0
    call boundary_procedure(SCAL2)

    if (gridtype==uniformgrid) then
     !$omp parallel private(i,j,k,p) reduction(max:S)
     !$omp do
     do bk = 1,Prnz,tilenz(narr2)
      do bj = 1,Prny,tileny(narr2)
       do bi = 1,Prnx,tilenx(narr2)
        do k = bk,min(bk+tilenz(narr2)-1,Prnz)
         do j = bj,min(bj+tileny(narr2)-1,Prny)
          do i = bi+mod(bi+j+k-1,2),min(bi+tilenx(narr2)-1,Prnx),2
            if (Prtype(i,j,k)<=0) then
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
            end if
          end do
         end do
        end do
       end do
      end do
     end do
     !$omp end do
     !$omp do
     do bk = 1,Prnz,tilenz(narr2)
      do bj = 1,Prny,tileny(narr2)
       do bi = 1,Prnx,tilenx(narr2)
        do k = bk,min(bk+tilenz(narr2)-1,Prnz)
         do j = bj,min(bj+tileny(narr2)-1,Prny)
          do i = bi+mod(bi+j+k,2),min(bi+tilenx(narr2)-1,Prnx),2
            if (Prtype(i,j,k)<=0) then
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
            end if
          end do
         end do
        end do
       end do
      end do
     end do
     !$omp end do
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
          end do
         end do
        end do
       end do
      end do
     end do
     !$omp end do
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
          end do
         end do
        end do
       end do
      end do
     end do
     !$omp end do
     !$omp endparallel
    end if
    write (*,*) "CNscalar ",l,S
    if (S<=epsCN) exit
   end do


  else


   call assign(SCAL2,SCAL)

   call boundary_procedure(SCAL2)

   !$omp parallel do private(i,j,k,bi,bj,bk,p,xi,yj,zk)
   do i = 1,ubound(ScalFlIBPoints,1)
     xi = ScalFlIBPoints(i)%xi
     yj = ScalFlIBPoints(i)%yj
     zk = ScalFlIBPoints(i)%zk
     SCAL2(xi,yj,zk) = SCAL2(xi,yj,zk) + TIBPoint_ScalFlSource(ScalFlIBPoints(i),Scal2,sctype) * dt
   end do
   !$omp end parallel do


  end if
  endsubroutine DIFFSCALAR







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
  end if
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
   end if
  end if
  endfunction AerResist

  pure real(knd) function DepositionVelocity3(dp,rhop,press,temp,z,z0,zL,ustar) !EMRAS recommended values
  real(knd),intent(in) :: dp,rhop,press,temp,z,z0,zL,ustar

   if (dp<1e-6) then
    DepositionVelocity3 = 0.5e-4
   else if (dp<2e-6) then
    DepositionVelocity3 = 1.5e-4
   else if (dp<10e-6) then
    DepositionVelocity3 = 10e-4
   else
    DepositionVelocity3 = 80e-4
   end if
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
   end if

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
   end if
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
     end do
    else
     do i = 1,num_of_scalars
      deptmp = abs(DepositionFlux(WMP,SCAL(WMP%xi,WMP%yj,WMP%zk,i),partdiam(i),partrho(i)))&
          *coef*dt*dxPr(WMP%xi)*dyPr(WMP%yj)
      WMP%depscalar(i) = WMP%depscalar(i)+deptmp/(dxPr(WMP%xi)*dyPr(WMP%yj)*dzPr(WMP%zk))
      SCAL(WMP%xi,WMP%yj,WMP%zk,i) = SCAL(WMP%xi,WMP%yj,WMP%zk,i)-deptmp/(dxPr(WMP%xi)*dyPr(WMP%yj)*dzPr(WMP%zk))
     end do
    end if
    end if
    if (associated(WMP%next)) then
     WMP => WMP%next
    else
     exit
    end if
   end do
   end if
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
      end do
     end do
    end do
    do k = 1,Prnz-1
     do j = 1,Prny
      do i = 1,Prnx
       SCAL(i,j,k+1,l) = SCAL(i,j,k+1,l)-flux(i,j,k)/(dxPr(i)*dyPr(j)*dzPr(k+1))
       SCAL(i,j,k,l) = SCAL(i,j,k,l)+flux(i,j,k)/(dxPr(i)*dyPr(j)*dzPr(k))
      end do
     end do
    end do
   end do
  end if
  endsubroutine Gravsettling


  pure real(knd) function Rig(i,j,k,U,V,temperature)
  integer,intent(in) :: i,j,k
  real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V
  real(knd),dimension(-1:,-1:,-1:),intent(in) :: temperature
  real(knd) num,denom

  num = (grav_acc/temperature_ref)*(temperature(i,j,k+1)-temperature(i,j,k-1))/(zPr(k+1)-zPr(k-1))
  denom = ((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._knd*(zPr(k+1)-zPr(k-1))))**2
  denom = denom+((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._knd*(zPr(k+1)-zPr(k-1))))**2

  if (abs(denom)>1E-5_knd*abs(num).and.abs(denom)>epsilon(1._knd)) then
   Rig = num/denom
  else
   Rig = 100000._knd*sign(1.0_knd,num)*sign(1.0_knd,denom)
  end if
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
   end if
  end subroutine ComputeTDiff



  subroutine InitScalarProfile(ScalarIn,ScalarProfile,default_value)
    real(knd),intent(inout) :: ScalarIn(-1:,-1:)
    type(TScalarProfile),intent(inout) :: ScalarProfile
    real(knd),intent(in) :: default_value
    integer   :: SectionToUse(-1:ubound(ScalarIn,2))
    integer   :: section,nSections,s
    integer   :: i,j,k
    real(knd) :: temp

    nSections = size(ScalarProfile%Sections)

    if (nSections > 0) then

       if (size(ScalarProfile%Sections)>0) then
         if (ScalarProfile%Sections(1)%jump<=0) ScalarProfile%Sections(1)%jump = temperature_ref
       end if

      ScalarProfile%Sections(1)%height = ScalarProfile%Sections(1)%top

      do i = 2,nSections
        if (ScalarProfile%Sections(i)%top < ScalarProfile%Sections(i-1)%top)&
          ScalarProfile%Sections(i)%top = ScalarProfile%Sections(i-1)%top

        ScalarProfile%Sections(i)%height = ScalarProfile%Sections(i)%top - ScalarProfile%Sections(i-1)%top
      end do

      if (ScalarProfile%Sections(nSections)%top < zW(Prnz+2) )&
          ScalarProfile%Sections(nSections)%top = zW(Prnz+2)

      section = 1

      do k = -1, Prnz+2
        do while (zPr(k) > ScalarProfile%Sections(section)%top) !should be safe, because last top is adjusted above
          section = section + 1
        end do
        SectionToUse(k) = section
      end do

    else

      SectionToUse = 0

    end if

    do k=-1,Prnz+2

      s = SectionToUse(k)

      if (s==0) then

        temp = default_value

      else

        temp = 0

        do i = 1, s-1
          temp = temp + ScalarProfile%Sections(i)%jump
          temp = temp + ScalarProfile%Sections(i)%height * ScalarProfile%Sections(i)%gradient
        end do

        temp = temp + ScalarProfile%Sections(s)%jump

        if (s>1) then
          temp = temp + (zPr(k) - ScalarProfile%Sections(s-1)%top) * ScalarProfile%Sections(s)%gradient
        else
          temp = temp + zPr(k) * ScalarProfile%Sections(s)%gradient
        end if

      end if

      do j=-1,Prny+2
          ScalarIn(j,k) = temp
      end do

    end do
  end subroutine InitScalarProfile


  subroutine InitScalar(ScalarIn,ScalarProfile,Sc)
    real(knd),intent(in) :: ScalarIn(-1:,-1:)
    type(TScalarProfile),intent(in) :: ScalarProfile
    real(knd),intent(out) :: Sc(-1:,-1:,-1:)
    real(knd) :: p_x,p_y,p_z,p
    integer   :: i,j,k

    if (ScalarProfile%randomize==1) then

!       p_x=0.5_knd
!       p_y=0.5_knd
!       p_z=0.5_knd

      do k=0,Prnz+1

!         if (zPr(k) <= ScalarProfile%randomizeTop) then
!           call random_number(p)
!           p = p - 0.5
!           p_z=(p_z+p)/2
!         end if

        do j=0,Prny+1

!           if (zPr(k) <= ScalarProfile%randomizeTop) then
!             call random_number(p)
!             p = p - 0.5
!             p_y=(p_y+p)/2
!           end if

           do i=0,Prnx+1

             if (zPr(k) <= ScalarProfile%randomizeTop) then
               call random_number(p)
               p = p - 0.5
!               p_x=(p_x+p)/2
!               p=p_x+p_y+p_z
             else
               p=0
             end if

             Sc(i,j,k)=ScalarIn(j,k) + ScalarProfile%randomizeAmplitude * 2 * p

           end do
        end do
      end do

    else

      forall(i=0:Prnx+1) Sc(i,:,:)=ScalarIn(:,:)

    end if

  end subroutine InitScalar


  subroutine InitHydrostaticPressure(Pr,Temperature,Moisture)
    real(knd),intent(out) :: Pr(1:,1:,1:)
    real(knd),intent(in) :: Temperature(-1:,-1:,-1:),Moisture(-1:,-1:,-1:)
    real(knd) t_virt,t_virt_prev,p
    integer i,j,k

    interface theta_v
      procedure :: theta_v_ijk
      procedure :: theta_v_Tq
    end interface


    p = bottom_pressure - (zW(Prnz+1)-zW(0))

    do k=1,Prnz
     t_virt = sum(TempIn(1:Prny,k))/Prny

     if (enable_moisture==1) t_virt = theta_v(t_virt,sum(MoistIn(1:Prny,k))/Prny)

     p = p + grav_acc*dzPr(k) * &
                  ( t_virt - temperature_ref )&
                  / temperature_ref

    end do
    top_pressure = p


    do j=1,Vny+1
      do i=1,Unx+1
        t_virt_prev = theta_v(i,j,Prnz)
        Pr(i,j,Prnz) = top_pressure - grav_acc*(zW(Prnz+1)-zPr(Prnz)) * &
                ( t_virt_prev - temperature_ref )&
                / temperature_ref
        do k=Prnz-1,1,-1
          t_virt = theta_v(i,j,k)
          Pr(i,j,k) = Pr(i,j,k+1) - &
                 grav_acc*dzW(k) * &
                ( (t_virt+t_virt_prev)/2._knd - temperature_ref )&
                / temperature_ref
          t_virt_prev = t_virt
        end do
      end do
    end do

    contains

      pure function theta_v_ijk(i,j,k) result(res)
        real(knd) :: res
        integer,intent(in) :: i,j,k

        if (enable_moisture==1) then
          res = Temperature(i,j,k) * (1._knd + 0.61_knd * Moisture(i,j,k))
        else
          res = Temperature(i,j,k)
        end if
      end function

      pure function theta_v_tq(T,q) result(res)
        real(knd) :: res
        real(knd),intent(in) :: T,q

        if (enable_moisture==1) then
          res = T * (1._knd + 0.61_knd * q)
        else
          res = T
        end if
      end function

  end subroutine


  subroutine InitSubsidenceProfile
    integer k

    if (SubsidenceGradient/=0) then
      allocate(SubsidenceProfile(0:Prnz))
      SubsidenceProfile = (/ (zW(k)*SubsidenceGradient, k=0,Prnz) /)
    else
      allocate(SubsidenceProfile(0))
    end if
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
    xc = 3*cos((x_axis_azimuth-70)*pi/180.)
    yc = 3*sin((x_axis_azimuth-70)*pi/180.)
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
          else if (zp<ct*0.4) then
           Scalar(xi,yj,zk,:) = Scalar(xi,yj,zk,:) + percdistrib(:)*0.8
           p = p+1
          else if (zp<ct*0.6) then
           Scalar(xi,yj,zk,:) = Scalar(xi,yj,zk,:) + percdistrib(:)*1.25
           p = p+1
          else if (zp<ct*0.8) then
           Scalar(xi,yj,zk,:) = Scalar(xi,yj,zk,:) + percdistrib(:)*1.75
           p = p+1
          else if (zp<=ct) then
           Scalar(xi,yj,zk,:) = Scalar(xi,yj,zk,:)  +percdistrib(:)*1.1
           p = p+1
          end if
         end if
        end do
       end do
      end do
      Scalar = totalscalsource*Scalar/p
      Scalar = Scalar/(dxmin*dymin*dzmin)
      released=.true.
     end if
    end if
  endsubroutine Release

end module Scalars

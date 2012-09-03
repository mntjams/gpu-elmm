module SCALARS
 use PARAMETERS
 use WALLMODELS
 use GEOMETRIC, only: TIBPoint, ScalFlIBPoints,TIBPoint_ScalFlSource
 use LIMITERS, only: Limiter, limparam
 use BOUNDARIES
 use TILING, only: tilenx, tileny, tilenz
#ifdef __HMPP
 use HMPP_CODELETS
#endif

implicit none
  private
  public Scalar, ScalarRK3, constPrt, Rig, partdiam, partrho, percdistrib,&
     Bound_Visc, Bound_Temp, Bound_PassScalar,&
     TTemperatureProfileSection, TTemperatureProfile, TemperatureProfile,&
     InitTemperatureProfile, InitTemperature,&
     SubsidenceProfile, SubsidenceGradient, InitSubsidenceProfile


  real(KND),allocatable,dimension(:,:,:,:) :: SCALAR  !last index is a number of scalar (because of paging)
  real(KND),dimension(:),allocatable :: partdiam,partrho,percdistrib !diameter of particles <=0 for gas

  real(KND),dimension(:),allocatable :: SubsidenceProfile
  real(KND) :: SubsidenceGradient = 0

  type TTemperatureProfileSection
    real(KND) :: top, height
    real(KND) :: jump
    real(KND) :: gradient
  end type TTemperatureProfileSection


  type TTemperatureProfile
    type(TTemperatureProfileSection), allocatable :: sections(:)
    integer   :: randomize
    real(KND) :: randomizeTop
    real(KND) :: randomizeAmplitude
  end type TTemperatureProfile

  type(TTemperatureProfile) :: TemperatureProfile

  real(KND),parameter :: constPrt = 0.6 !constant value of Prt, which may be refined further in this module

contains


  subroutine ScalarRK3(U,V,W,Temperature,Scalar,RKstage,fluxprofile)
    real(KND),intent(in)    :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(KND),intent(inout) :: Temperature(-1:,-1:,-1:),Scalar(-1:,-1:,-1:,-1:)
    real(KND),intent(out),optional   :: fluxprofile(0:)
    integer,intent(in)      :: RKstage

    real(KND),dimension(:,:,:,:),allocatable,save ::Scalar_adv,Scalar_2
    real(KND),dimension(:,:,:),allocatable,save   ::Temperature_adv,Temperature2

    real(KND),dimension(1:3),parameter :: alpha = (/ 4._KND/15._KND, 1._KND/15._KND, 1._KND/6._KND /)
    real(KND),dimension(1:3),parameter :: beta  = (/ 8._KND/15._KND, 5._KND/12._KND, 3._KND/4._KND /)
    real(KND),dimension(1:3),parameter :: rho   = (/       0._KND, -17._KND/60._KND,-5._KND/12._KND/)

    integer :: i
    integer,save :: called = 0
    logical,save :: released=.false.

#ifdef __HMPP
    integer iters
    real(KND) res
#endif

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
      !$hmpp <tsteps> advancedload, args[DiffTemperature::TBtype,DiffTemperature::sideTemp,DiffTemperature::TempIn]
      !$hmpp <tsteps> advancedload, args[DiffTemperature::epsCN,DiffTemperature::sideTemp,DiffTemperature::maxCNiter]
    endif


!     if (computescalars>0.and..not.released) call Release(Scalar,released)


    if (RKstage == 1) then
      !$omp parallel
      !$omp workshare
      temperature_adv = 0
      Scalar_adv = 0
      !$omp end workshare
      !$omp end parallel
    endif

    if (computescalars>0) then

      !$omp parallel
      !$omp workshare
      Scalar_2=0
      !$omp end workshare

      if (RKstage>1) then
        !$omp workshare
        Scalar_2 = Scalar_2+Scalar_adv*rho(RKstage)
        !$omp end workshare
      endif

      !$omp workshare
      Scalar_adv = 0
      !$omp end workshare
      !$omp end parallel

      do i=1,computescalars
        call AdvScalar(Scalar_adv(:,:,:,i),Scalar(:,:,:,i),U,V,W,2,1._KND,SubsidenceProfile)
      enddo

      !$omp parallel
      !$omp workshare
      Scalar_2 = Scalar_2+Scalar_adv*beta(RKstage)
      !$omp end workshare
      !$omp end parallel

      if (scalsourcetype==pointsource) then

        do i=1,computescalars
          Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i) = Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i)+&
           percdistrib(i)*(rho(RKstage)+beta(RKstage))*dt*totalscalsource/(dxPr(scalsrci(i))*dyPr(scalsrcj(i))*dzPr(scalsrck(i)))
        enddo

      elseif (scalsourcetype==volumesource) then

        do i=1,computescalars
          call VolScalSource(Scalar_2,rho(RKstage)+beta(RKstage))
        enddo

      endif

      !$omp parallel
      !$omp workshare
      Scalar = Scalar+Scalar_2
      !$omp end workshare
      !$omp end parallel

      if (sgstype/=StabSmagorinskyModel)  call ComputeTDiff(U,V,W)

      call Bound_Visc(TDiff)

      do i=1,computescalars
         call DiffScalar(Scalar_2(:,:,:,i),Scalar(:,:,:,i),2,2._KND*alpha(RKstage))
      enddo

      if (computedeposition>0) call Deposition(Scalar_2,2._KND*alpha(RKstage))

      if (computegravsettling>0) call GravSettling(Scalar_2,2._KND*alpha(RKstage))

      !$omp parallel
      !$omp workshare
      Scalar = Scalar_2
      !$omp end workshare
      !$omp end parallel

    endif


    if (buoyancy>0) then

      if (sgstype/=StabSmagorinskyModel)  call ComputeTDiff(U,V,W)

      call Bound_Visc(TDiff)

      call Bound_Temp(temperature)

      if (RKstage>1) then
        !$omp parallel
        !$omp workshare
        temperature2 = temperature_adv*rho(RKstage)
        !$omp end workshare
        !$omp end parallel
      else
        !$omp parallel
        !$omp workshare
        temperature2=0
        !$omp end workshare
        !$omp end parallel
      endif

      !$omp parallel
      !$omp workshare
      temperature_adv = 0
      !$omp end workshare
      !$omp end parallel

      if (RKstage==3.and.present(fluxprofile)) then
        call AdvScalar(temperature_adv,temperature,U,V,W,1,1._KND,SubsidenceProfile,FluxProfile)
      else
        call AdvScalar(temperature_adv,temperature,U,V,W,1,1._KND,SubsidenceProfile)
      end if

      !$omp parallel
      !$omp workshare
      temperature2 = temperature2+temperature_adv*beta(RKstage)

      temperature = temperature+temperature2
      !$omp end workshare
      !$omp end parallel

      call Bound_Temp(temperature)

#ifdef __HMPP
        !$hmpp <tsteps> advancedload, args[DiffTemperature::coef]
        if (size(BsideTArr)>0) then
          !$hmpp <tsteps> advancedload, args[DiffTemperature::BsideTArr]
        endif
        if (size(BsideTFLArr)>0) then
          !$hmpp <tsteps> advancedload, args[DiffTemperature::BsideTFLArr]
        endif
        !$hmpp <tsteps> advancedload, args[DiffTemperature::TDiff,DiffTemperature::Temperature]
        !$hmpp <tsteps> DiffTemperature callsite, args[*].noupdate=true
        call DiffTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,maxCNiter,epsCN,Re,&
                            TBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr,TDiff,&
                            temperature2,temperature,2._KND*alpha(RKstage),dt,iters,res)
        !$hmpp <tsteps> delegatedstore, args[DiffTemperature::l,DiffTemperature::res]
        !$hmpp <tsteps> delegatedstore, args[DiffTemperature::Temperature2]
        write (*,*) "CNscalar ",iters,res
#else
        call DiffScalar(temperature2,temperature,1,2._KND*alpha(RKstage))
#endif

      !$omp parallel
      !$omp workshare
      temperature = temperature2
      !$omp end workshare
      !$omp end parallel
    endif

  end subroutine ScalarRK3




  subroutine ADVSCALAR(SCAL2,SCAL,U,V,W,sctype,coef,SubsidenceProfile,fluxProfile)
  real(KND),intent(inout) :: Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in)    :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in)      :: sctype
  real(KND),intent(in),allocatable :: SubsidenceProfile(:)
  real(KND),intent(out),optional :: fluxProfile(0:)
  real(KND) :: SubsidenceProfileLoc(0:Prnz)
  real(KND) :: fluxProfileLoc(0:Prnz)

   if (allocated(SubsidenceProfile)) then
     SubsidenceProfileLoc = SubsidenceProfile(0:Prnz)
   else
     SubsidenceProfileLoc = 0
   endif

   if (gridtype==uniformgrid) then
       call KAPPASCALARUG(SCAL2,SCAL,U,V,W,sctype,coef,SubsidenceProfileLoc,fluxProfileLoc)
   else
       call KAPPASCALARGG(SCAL2,SCAL,U,V,W,sctype,coef,SubsidenceProfileLoc,fluxProfileLoc)
   endif

   if (present(fluxProfile)) fluxProfile(0:Prnz) = fluxProfileLoc

  endsubroutine ADVSCALAR

  subroutine CDSSCALAR(SCAL2,SCAL,U,V,W,sctype,coef)
  real(KND),intent(inout) :: Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in)    :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in)      :: sctype
  integer nx,ny,nz,i,j,k
  real(KND) Ax,Ay,Az

        nx = Prnx
        ny = Prny
        nz = Prnz

        if (sctype==1) then
         call BOUND_Temp(SCAL)
        else
         call BOUND_PASSSCALAR(SCAL)
        endif

        Ax = 0.5_KND*coef*dt/dxmin
        Ay = 0.5_KND*coef*dt/dymin
        Az = 0.5_KND*coef*dt/dzmin


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

  subroutine UDSSCALAR(SCAL2,SCAL,U,V,W,sctype,coef)
  real(KND),intent(inout) :: Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in) :: sctype
  integer nx,ny,nz,i,j,k
  real(KND) Ax,Ay,Az

        nx = Prnx
        ny = Prny
        nz = Prnz

        if (sctype==1) then
         call BOUND_Temp(SCAL)
        else
         call BOUND_PASSSCALAR(SCAL)
        endif

        Ax = coef*dt/dxmin
        Ay = coef*dt/dymin
        Az = coef*dt/dzmin

        do k = 1,nz
            do j = 1,ny
                do i = 1,nx
                    if (U(i,j,k)>0) then
                     SCAL2(i,j,k) = SCAL2(i,j,k)-Ax*SCAL(i,j,k)*(U(i,j,k))
                    else
                     SCAL2(i,j,k) = SCAL2(i,j,k)-Ax*SCAL(i+1,j,k)*(U(i,j,k))
                    endif
                    if (U(i-1,j,k)>0) then
                     SCAL2(i,j,k) = SCAL2(i,j,k)+Ax*SCAL(i-1,j,k)*(U(i-1,j,k))
                    else
                     SCAL2(i,j,k) = SCAL2(i,j,k)+Ax*SCAL(i,j,k)*(U(i-1,j,k))
                    endif

                     if (V(i,j,k)>0) then
                      SCAL2(i,j,k) = SCAL2(i,j,k)-Ay*SCAL(i,j,k)*(V(i,j,k))
                     else
                      SCAL2(i,j,k) = SCAL2(i,j,k)-Ay*SCAL(i,j+1,k)*(V(i,j,k))
                     endif
                     if (V(i,j-1,k)>0) then
                     SCAL2(i,j,k) = SCAL2(i,j,k)+Ay*SCAL(i,j-1,k)*(V(i,j-1,k))
                    else
                     SCAL2(i,j,k) = SCAL2(i,j,k)+Ay*SCAL(i,j,k)*(V(i,j-1,k))
                    endif


                    if (W(i,j,k)>0) then
                     SCAL2(i,j,k) = SCAL2(i,j,k)-Az*SCAL(i,j,k)*(W(i,j,k))
                    else
                     SCAL2(i,j,k) = SCAL2(i,j,k)-Az*SCAL(i,j,k+1)*(W(i,j,k))
                    endif
                    if (W(i,j,k-1)>0) then
                     SCAL2(i,j,k) = SCAL2(i,j,k)+Az*SCAL(i,j,k-1)*(W(i,j,k-1))
                    else
                     SCAL2(i,j,k) = SCAL2(i,j,k)+Az*SCAL(i,j,k)*(W(i,j,k-1))
                    endif

                enddo
            enddo
        enddo
  end subroutine UDSSCALAR



  subroutine QUICKSCALAR(SCAL2,SCAL,U,V,W,sctype,coef)
  real(KND),intent(inout) :: Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in) :: sctype
  real(KND) :: Scal3(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  integer nx,ny,nz,i,j,k
  real(KND) Ax,Ay,Az,F

    if (sctype==1) then
     call BOUND_Temp(SCAL)
    else
     call BOUND_PASSSCALAR(SCAL)
    endif

    Ax = coef*dt/dxmin
    Ay = coef*dt/dymin
    Az = coef*dt/dzmin
    SCAL2 = SCAL
    SCAL3 = SCAL
    do k = 1,Prnz
     do j = 1,Prny
      do i = 0,Prnx
        if (U(i,j,k)>0) then
         F = Ax*U(i,j,k)*((6._KND/8._KND)*SCAL3(i,j,k)+(3._KND/8._KND)*SCAL3(i-1,j,k)-(1._KND/8._KND)*SCAL3(i-1,j,k))
        else
         F = Ax*U(i,j,k)*((6._KND/8._KND)*SCAL3(i+1,j,k)+(3._KND/8._KND)*SCAL3(i,j,k)-(1._KND/8._KND)*SCAL3(i+2,j,k))
        endif
        SCAL2(i,j,k) = SCAL2(i,j,k)-F
        SCAL2(i+1,j,k) = SCAL2(i+1,j,k)+F
      enddo
     enddo
    enddo

    if (sctype==1) then
     call BOUND_Temp(SCAL2)
    else
     call BOUND_PASSSCALAR(SCAL2)
    endif
    SCAL3 = SCAL2
    do k = 1,Prnz
     do j = 0,Prny
      do i = 1,Prnx
        if (V(i,j,k)>0) then
         F = Ay*V(i,j,k)*((6._KND/8._KND)*SCAL3(i,j,k)+(3._KND/8._KND)*SCAL3(i,j-1,k)-(1._KND/8._KND)*SCAL3(i,j-1,k))
        else
         F = Ay*V(i,j,k)*((6._KND/8._KND)*SCAL3(i,j+1,k)+(3._KND/8._KND)*SCAL3(i,j,k)-(1._KND/8._KND)*SCAL3(i,j+1,k))
        endif
        SCAL2(i,j,k) = SCAL2(i,j,k)-F
        SCAL2(i,j+1,k) = SCAL2(i,j+1,k)+F
      enddo
     enddo
    enddo

    if (sctype==1) then
     call BOUND_Temp(SCAL2)
    else
     call BOUND_PASSSCALAR(SCAL2)
    endif
    SCAL3 = SCAL2
    do k = 0,Prnz
     do j = 1,Prny
      do i = 1,Prnx
        if (W(i,j,k)>0) then
         F = Az*W(i,j,k)*((6._KND/8._KND)*SCAL3(i,j,k)+(3._KND/8._KND)*SCAL3(i,j,k-1)-(1._KND/8._KND)*SCAL3(i,j,k-1))
        else
         F = Az*W(i,j,k)*((6._KND/8._KND)*SCAL3(i,j,k+1)+(3._KND/8._KND)*SCAL3(i,j,k)-(1._KND/8._KND)*SCAL3(i,j,k+1))
        endif
        SCAL2(i,j,k) = SCAL2(i,j,k)-F
        SCAL2(i,j,k+1) = SCAL2(i,j,k+1)+F
      enddo
     enddo
    enddo
    SCAL2 = SCAL2-SCAL
  end subroutine QUICKSCALAR





  subroutine PLMSCALAR(SCAL2,SCAL,U,V,W,sctype,coef)
  real(KND),intent(inout) :: Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in) :: sctype
  integer i,j,k
  real(KND) SL,SR,FLUX,A
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2) :: SLOPE
  logical,save :: direction=.true.

   if (sctype==1) then
    call BOUND_Temp(SCAL)
   else
    call BOUND_PASSSCALAR(SCAL)
   endif

  A = coef*dt
  SCAL2 = SCAL2+SCAL

  if (direction) then
   SLOPE = 0
   do i = 0,Prnx+1
    do j = 0,Prny+1
     do k = 0,Prnz+1
      SR = (SCAL2(i+1,j,k)-SCAL2(i,j,k))/dxU(i)
      SL = (SCAL2(i,j,k)-SCAL2(i-1,j,k))/dxU(i-1)
      SLOPE(i,j,k) = LIMITER(SL,SR)*dxPr(i)
     enddo
    enddo
   enddo


   do i = 0,Prnx
    do j = 1,Prny
     do k = 1,Prnz
      if (U(i,j,k)>0) then
       FLUX = U(i,j,k)*SCAL2(i,j,k)+U(i,j,k)*(1-U(i,j,k)*A/dxPr(i))*SLOPE(i,j,k)/2._KND
      else
       FLUX = U(i,j,k)*SCAL2(i+1,j,k)-U(i,j,k)*(1+U(i,j,k)*A/dxPr(i+1))*SLOPE(i+1,j,k)/2._KND
      endif
      SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dxPr(i)
      SCAL2(i+1,j,k) = SCAL2(i+1,j,k)+A*FLUX/dxPr(i+1)
     enddo
    enddo
   enddo

   SLOPE = 0
   do i = 0,Prnx+1
    do j = 0,Prny+1
     do k = 0,Prnz+1
      SR = (SCAL2(i,j+1,k)-SCAL2(i,j,k))/dyV(j)
      SL = (SCAL2(i,j,k)-SCAL2(i,j-1,k))/dyV(j-1)
      SLOPE(i,j,k) = LIMITER(SL,SR)*dyPr(j)
     enddo
    enddo
   enddo


   do i = 1,Prnx
    do j = 0,Prny
     do k = 1,Prnz
      if (V(i,j,k)>0) then
       FLUX = V(i,j,k)*SCAL2(i,j,k)+V(i,j,k)*(1-V(i,j,k)*A/dyPr(j))*SLOPE(i,j,k)/2._KND
      else
       FLUX = V(i,j,k)*SCAL2(i,j+1,k)-V(i,j,k)*(1+V(i,j,k)*A/dyPr(j+1))*SLOPE(i,j+1,k)/2._KND
      endif
      SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dyPr(j)
      SCAL2(i,j+1,k) = SCAL2(i,j+1,k)+A*FLUX/dyPr(j+1)
     enddo
    enddo
   enddo

   SLOPE = 0
   do i = 0,Prnx+1
    do j = 0,Prny+1
     do k = 0,Prnz+1
      SR = (SCAL2(i,j,k+1)-SCAL2(i,j,k))/dzW(k)
      SL = (SCAL2(i,j,k)-SCAL2(i,j,k-1))/dzW(k-1)
      SLOPE(i,j,k) = LIMITER(SL,SR)*dzPr(k)
     enddo
    enddo
   enddo

   do i = 1,Prnx
    do j = 1,Prny
     do k = 0,Prnz
      if (W(i,j,k)>0) then
       FLUX = W(i,j,k)*SCAL2(i,j,k)+W(i,j,k)*(1-W(i,j,k)*A/dzPr(k))*SLOPE(i,j,k)/2._KND
      else
       FLUX = W(i,j,k)*SCAL2(i,j,k+1)-W(i,j,k)*(1+W(i,j,k)*A/dzPr(k+1))*SLOPE(i,j,k+1)/2._KND
      endif
      SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dzPr(k)
      SCAL2(i,j,k+1) = SCAL2(i,j,k+1)+A*FLUX/dzPr(k+1)
     enddo
    enddo
   enddo
  else
   SLOPE = 0
   do i = 0,Prnx+1
    do j = 0,Prny+1
     do k = 0,Prnz+1
      SR = (SCAL2(i,j,k+1)-SCAL2(i,j,k))/dzW(k)
      SL = (SCAL2(i,j,k)-SCAL2(i,j,k-1))/dzW(k-1)
      SLOPE(i,j,k) = LIMITER(SL,SR)*dzPr(k)
     enddo
    enddo
   enddo

   do i = 1,Prnx
    do j = 1,Prny
     do k = 0,Prnz
      if (W(i,j,k)>0) then
       FLUX = W(i,j,k)*SCAL2(i,j,k)+W(i,j,k)*(1-W(i,j,k)*A/dzPr(k))*SLOPE(i,j,k)/2._KND
      else
       FLUX = W(i,j,k)*SCAL2(i,j,k+1)-W(i,j,k)*(1+W(i,j,k)*A/dzPr(k+1))*SLOPE(i,j,k+1)/2._KND
      endif
      SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dzPr(k)
      SCAL2(i,j,k+1) = SCAL2(i,j,k+1)+A*FLUX/dzPr(k+1)
     enddo
    enddo
   enddo

   SLOPE = 0
   do i = 0,Prnx+1
    do j = 0,Prny+1
     do k = 0,Prnz+1
      SR = (SCAL2(i,j+1,k)-SCAL2(i,j,k))/dyV(j)
      SL = (SCAL2(i,j,k)-SCAL2(i,j-1,k))/dyV(j-1)
      SLOPE(i,j,k) = LIMITER(SL,SR)*dyPr(j)
     enddo
    enddo
   enddo


   do i = 1,Prnx
    do j = 0,Prny
     do k = 1,Prnz
      if (V(i,j,k)>0) then
       FLUX = V(i,j,k)*SCAL2(i,j,k)+V(i,j,k)*(1-V(i,j,k)*A/dyPr(j))*SLOPE(i,j,k)/2._KND
      else
       FLUX = V(i,j,k)*SCAL2(i,j+1,k)-V(i,j,k)*(1+V(i,j,k)*A/dyPr(j+1))*SLOPE(i,j+1,k)/2._KND
      endif
      SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dyPr(j)
      SCAL2(i,j+1,k) = SCAL2(i,j+1,k)+A*FLUX/dyPr(j+1)
     enddo
    enddo
   enddo

   SLOPE = 0
   do i = 0,Prnx+1
    do j = 0,Prny+1
     do k = 0,Prnz+1
      SR = (SCAL2(i+1,j,k)-SCAL2(i,j,k))/dxU(i)
      SL = (SCAL2(i,j,k)-SCAL2(i-1,j,k))/dxU(i-1)
      SLOPE(i,j,k) = LIMITER(SL,SR)*dxPr(i)
     enddo
    enddo
   enddo


   do i = 0,Prnx
    do j = 1,Prny
     do k = 1,Prnz
      if (U(i,j,k)>0) then
       FLUX = U(i,j,k)*SCAL2(i,j,k)+U(i,j,k)*(1-U(i,j,k)*A/dxPr(i))*SLOPE(i,j,k)/2._KND
      else
       FLUX = U(i,j,k)*SCAL2(i+1,j,k)-U(i,j,k)*(1+U(i,j,k)*A/dxPr(i+1))*SLOPE(i+1,j,k)/2._KND
      endif
      SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dxPr(i)
      SCAL2(i+1,j,k) = SCAL2(i+1,j,k)+A*FLUX/dxPr(i+1)
     enddo
    enddo
   enddo
  endif

  SCAL2 = SCAL2-SCAL
  endsubroutine PLMSCALAR



  subroutine PLMSCALARNOSPLIT(SCAL2,SCAL,U,V,W,sctype,coef)
  real(KND),intent(inout) :: Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in) :: sctype
  integer i,j,k
  real(KND) SL,SR,FLUX,A
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2) :: SLOPE

   if (sctype==1) then
    call BOUND_Temp(SCAL)
   else
    call BOUND_PASSSCALAR(SCAL)
   endif

  A = coef*dt

  SLOPE = 0
  do i = 0,Prnx+1
   do j = 0,Prny+1
    do k = 0,Prnz+1
     SR = (SCAL(i+1,j,k)-SCAL(i,j,k))/dxU(i)
     SL = (SCAL(i,j,k)-SCAL(i-1,j,k))/dxU(i-1)
     SLOPE(i,j,k) = LIMITER(SL,SR)*dxPr(i)
    enddo
   enddo
  enddo


  do i = 0,Prnx
   do j = 1,Prny
    do k = 1,Prnz
     if (U(i,j,k)>0) then
      FLUX = U(i,j,k)*SCAL(i,j,k)+U(i,j,k)*(1-U(i,j,k)*A/dxPr(i))*SLOPE(i,j,k)/2._KND
     else
      FLUX = U(i,j,k)*SCAL(i+1,j,k)-U(i,j,k)*(1+U(i,j,k)*A/dxPr(i+1))*SLOPE(i+1,j,k)/2._KND
     endif
     SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dxPr(i)
     SCAL2(i+1,j,k) = SCAL2(i+1,j,k)+A*FLUX/dxPr(i+1)
    enddo
   enddo
  enddo

  SLOPE = 0
  do i = 0,Prnx+1
   do j = 0,Prny+1
    do k = 0,Prnz+1
     SR = (SCAL(i,j+1,k)-SCAL(i,j,k))/dyV(j)
     SL = (SCAL(i,j,k)-SCAL(i,j-1,k))/dyV(j-1)
     SLOPE(i,j,k) = LIMITER(SL,SR)*dyPr(j)
    enddo
   enddo
  enddo


  do i = 1,Prnx
   do j = 0,Prny
    do k = 1,Prnz
     if (V(i,j,k)>0) then
      FLUX = V(i,j,k)*SCAL(i,j,k)+V(i,j,k)*(1-V(i,j,k)*A/dyPr(j))*SLOPE(i,j,k)/2._KND
     else
      FLUX = V(i,j,k)*SCAL(i,j+1,k)-V(i,j,k)*(1+V(i,j,k)*A/dyPr(j+1))*SLOPE(i,j+1,k)/2._KND
     endif
     SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dyPr(j)
     SCAL2(i,j+1,k) = SCAL2(i,j+1,k)+A*FLUX/dyPr(j+1)
    enddo
   enddo
  enddo

  SLOPE = 0
  do i = 0,Prnx+1
   do j = 0,Prny+1
    do k = 0,Prnz+1
     SR = (SCAL(i,j,k+1)-SCAL(i,j,k))/dzW(k)
     SL = (SCAL(i,j,k)-SCAL(i,j,k-1))/dzW(k-1)
     SLOPE(i,j,k) = LIMITER(SL,SR)*dzPr(k)
    enddo
   enddo
  enddo

  do i = 1,Prnx
   do j = 1,Prny
    do k = 0,Prnz
     if (W(i,j,k)>0) then
      FLUX = W(i,j,k)*SCAL(i,j,k)+W(i,j,k)*(1-W(i,j,k)*A/dzPr(k))*SLOPE(i,j,k)/2._KND
     else
      FLUX = W(i,j,k)*SCAL(i,j,k+1)-W(i,j,k)*(1+W(i,j,k)*A/dzPr(k+1))*SLOPE(i,j,k+1)/2._KND
     endif
     SCAL2(i,j,k) = SCAL2(i,j,k)-A*FLUX/dzPr(k)
     SCAL2(i,j,k+1) = SCAL2(i,j,k+1)+A*FLUX/dzPr(k+1)
    enddo
   enddo
  enddo
  endsubroutine PLMSCALARNOSPLIT



!   subroutine KAPPASCALAR(SCAL2,SCAL,U,V,W,sctype,coef,fluxProfile,SubsidenceProfile) !Kappa scheme with flux limiter
!   real(KND),intent(inout) :: Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:) !Hunsdorfer et al. 1995, JCP
!   real(KND),intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
!   integer,intent(in) :: sctype
!
!
!   endsubroutine KAPPASCALAR


  subroutine KAPPASCALARUG(SCAL2,SCAL,U,V,W,sctype,coef,SubsidenceProfile,fluxProfile) !Kappa scheme with flux limiter
  real(KND),intent(inout) ::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:) !Hunsdorfer et al. 1995, JCP
  real(KND),intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in) :: sctype
  real(KND),intent(in) :: SubsidenceProfile(0:)
  real(KND),intent(out) :: fluxProfile(0:)
  integer i,j,k
  real(KND) A,Ax,Ay,Az              !Auxiliary variables to store muliplication constants for efficiency
  real(KND) vel,SL,SR,FLUX
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2) :: SLOPE
  real(KND),parameter ::eps = 1e-8
  real(KND) :: FluxLimiter,r

  FluxLimiter(r)=max(0._KND,min(2._KND*r,min(limparam,(1+2._KND*r)/3._KND)))

  if (sctype==1) then
    call BOUND_Temp(SCAL)
  else
    call BOUND_PASSSCALAR(SCAL)
  endif

  A = coef*dt
  Ax = coef*dt/dxmin
  Ay = coef*dt/dymin
  Az = coef*dt/dzmin


  !$omp parallel private(i,j,k,vel,SL,SR,FLUX)

  !$omp workshare
  SLOPE = 0
  !$omp end workshare
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
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
    enddo
   enddo
  enddo
  !$omp end do

  !$omp do
  do k = 1,Prnz
   do j = 1,Prny
    do i = 0,Prnx
     if (U(i,j,k)>0) then
      FLUX = U(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i-1,j,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX = U(i,j,k)*(SCAL(i+1,j,k)+(SCAL(i+1,j,k)-SCAL(i+2,j,k))*SLOPE(i,j,k)/2._KND)
     endif
     SCAL2(i,j,k) = SCAL2(i,j,k)-Ax*FLUX
     SCAL2(i+1,j,k) = SCAL2(i+1,j,k)+Ax*FLUX
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
      SR = (SCAL(i,j+1,k)-SCAL(i,j,k))
      SL = (SCAL(i,j,k)-SCAL(i,j-1,k))
     else
      SR = (SCAL(i,j,k)-SCAL(i,j+1,k))
      SL = (SCAL(i,j+1,k)-SCAL(i,j+2,k))
     endif
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
    enddo
   enddo
  enddo
  !$omp end do


  !$omp do
  do k = 1,Prnz
   do j = 0,Prny
    do i = 1,Prnx
     if (V(i,j,k)>0) then
      FLUX = V(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j-1,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX = V(i,j,k)*(SCAL(i,j+1,k)+(SCAL(i,j+1,k)-SCAL(i,j+2,k))*SLOPE(i,j,k)/2._KND)
     endif

     SCAL2(i,j,k) = SCAL2(i,j,k)-Ay*FLUX
     SCAL2(i,j+1,k) = SCAL2(i,j+1,k)+Ay*FLUX
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
      SR = (SCAL(i,j,k+1)-SCAL(i,j,k))
      SL = (SCAL(i,j,k)-SCAL(i,j,k-1))
     else
      SR = (SCAL(i,j,k)-SCAL(i,j,k+1))
      SL = (SCAL(i,j,k+1)-SCAL(i,j,k+2))
     endif
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
    enddo
   enddo
  enddo
  !$omp end do nowait

  !$omp workshare
  fluxprofile = 0
  !$omp end workshare

  !$omp do reduction(+:fluxprofile)
  do j = 1,Prny  !loop order due to avoid race condition
   do k = 0,Prnz
    do i = 1,Prnx

     vel = W(i,j,k) - SubsidenceProfile(k)

     if (vel>0) then
      FLUX = vel*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j,k-1))*SLOPE(i,j,k)/2._KND)
     else
      FLUX = vel*(SCAL(i,j,k+1)+(SCAL(i,j,k+1)-SCAL(i,j,k+2))*SLOPE(i,j,k)/2._KND)
     endif

     if (vel>=1e-6) fluxprofile(k) = fluxprofile(k) + FLUX/vel*W(i,j,k)

     SCAL2(i,j,k) = SCAL2(i,j,k) - Az*FLUX
     SCAL2(i,j,k+1) = SCAL2(i,j,k+1) + Az*FLUX
    enddo
   enddo
  enddo
  !$omp end do

  !$omp workshare
  fluxprofile = fluxprofile / (Prnx*Prny)
  !$omp end workshare

  !$omp end parallel

  endsubroutine KAPPASCALARUG



  subroutine KAPPASCALARGG(SCAL2,SCAL,U,V,W,sctype,coef,SubsidenceProfile,fluxProfile) !Kappa scheme with flux limiter
  real(KND),intent(inout) ::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:) !Hunsdorfer et al. 1995, JCP
  real(KND),intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),coef
  integer,intent(in) :: sctype
  real(KND),intent(in) :: SubsidenceProfile(0:)
  real(KND),intent(out) :: fluxProfile(0:)
  integer i,j,k
  real(KND) A                       !Auxiliary variables to store muliplication constants for efficiency
  real(KND) vel,SL,SR,FLUX
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2) :: SLOPE
  real(KND),parameter::eps = 1e-8
  real(KND) :: FluxLimiter,r

  FluxLimiter(r)=max(0._KND,min(2._KND*r,min(limparam,(1+2._KND*r)/3._KND)))

  if (sctype==1) then
   call BOUND_Temp(SCAL)
  else
   call BOUND_PASSSCALAR(SCAL)
  endif

  A = coef*dt

  !$omp parallel private(i,j,k,vel,SL,SR,FLUX)

  !$omp workshare
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
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
    enddo
   enddo
  enddo
  !$omp end do


  !$omp do
  do k = 1,Prnz
   do j = 1,Prny
    do i = 0,Prnx
     if (U(i,j,k)>0) then
      FLUX = U(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i-1,j,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX = U(i,j,k)*(SCAL(i+1,j,k)+(SCAL(i+1,j,k)-SCAL(i+2,j,k))*SLOPE(i,j,k)/2._KND)
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
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
    enddo
   enddo
  enddo
  !$omp end do


  !$omp do
  do k = 1,Prnz
   do j = 0,Prny
    do i = 1,Prnx
     if (V(i,j,k)>0) then
      FLUX = V(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j-1,k))*SLOPE(i,j,k)/2._KND)
     else
      FLUX = V(i,j,k)*(SCAL(i,j+1,k)+(SCAL(i,j+1,k)-SCAL(i,j+2,k))*SLOPE(i,j,k)/2._KND)
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
     SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
    enddo
   enddo
  enddo
  !$omp end do nowait

  !$omp workshare
  fluxprofile = 0
  !$omp end workshare

  !$omp do reduction(+:fluxprofile)
  do j = 1,Prny  !loop order due to avoid race condition
   do k = 0,Prnz
    do i = 1,Prnx

     vel = W(i,j,k) - SubsidenceProfile(k)

     if (vel>0) then
      FLUX = vel*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j,k-1))*SLOPE(i,j,k)/2._KND)
     else
      FLUX = vel*(SCAL(i,j,k+1)+(SCAL(i,j,k+1)-SCAL(i,j,k+2))*SLOPE(i,j,k)/2._KND)
     endif

     if (vel>=1e-6) fluxprofile(k) = fluxprofile(k+1) + FLUX/vel*W(i,j,k)

     SCAL2(i,j,k) = SCAL2(i,j,k) - A*FLUX/dzPr(k)
     SCAL2(i,j,k+1) = SCAL2(i,j,k+1) + A*FLUX/dzPr(k+1)
    enddo
   enddo
  enddo
  !$omp end do

  !$omp workshare
  fluxProfile = fluxProfile/(Prnx * Prny)
  !$omp end workshare

  !$omp end parallel

  endsubroutine KAPPASCALARGG




  subroutine DIFFSCALAR(SCAL2,SCAL,sctype,coef)
  real(KND),intent(inout) ::Scal2(-1:,-1:,-1:),Scal(-1:,-1:,-1:)
  real(KND),intent(in) :: coef
  integer,intent(in) :: sctype
  real(KND) Scal3(-1:Prnx+1,-1:Prny+1,-1:Prnz+1)
  integer nx,ny,nz,i,j,k,bi,bj,bk,l,xi,yj,zk
  real(KND) p,S
  real(KND) A,Ax,Ay,Az,Ap(-1:Prnx+1,-1:Prny+1,-1:Prnz+1)
  integer,parameter :: narr = 3, narr2 = 5

   nx = Prnx
   ny = Prny
   nz = Prnz


  if (Re>0) then
   if (sctype==1) then
    call BOUND_Temp(SCAL)
   else
    call BOUND_PASSSCALAR(SCAL)
   endif

   A = dt*coef
   Ax = 1._KND/(4._KND*dxmin**2)
   Ay = 1._KND/(4._KND*dymin**2)
   Az = 1._KND/(4._KND*dzmin**2)

   !$omp parallel private (i,j,k,bi,bj,bk)

   !initital value using forward Euler
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

   !$omp workshare
   SCAL2(1:Prnx,1:Prny,1:Prnz) = SCAL(1:Prnx,1:Prny,1:Prnz) + A * SCAL3(1:Prnx,1:Prny,1:Prnz)
   !$omp end workshare
   !$omp end parallel

   if (sctype==1) then
     call BOUND_Temp(SCAL2)
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
           Ap(i,j,k) = 1._KND/(1._KND/A+(((TDiff(i+1,j,k)+TDiff(i,j,k))+&
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
           Ap(i,j,k) = 1._KND/(1._KND/A+(((TDiff(i+1,j,k)+TDiff(i,j,k))/dxU(i)+&
                                (TDiff(i,j,k)+TDiff(i-1,j,k))/dxU(i-1))/(4._KND*dxPr(i))+&
                                ((TDiff(i,j+1,k)+TDiff(i,j,k))/dyV(j)+&
                                (TDiff(i,j,k)+TDiff(i,j-1,k))/dyV(j-1))/(4._KND*dyPr(j))+&
                                ((TDiff(i,j,k+1)+TDiff(i,j,k))/dzW(k)+&
                                (TDiff(i,j,k)+TDiff(i,j,k-1))/dzW(k-1))/(4._KND*dzPr(k))))
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
     call BOUND_Temp(SCAL2)
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
            p = (SCAL(i,j,k)/A)+(SCAL3(i,j,k)/4._KND+&
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
            p = (SCAL(i,j,k)/A)+(SCAL3(i,j,k)/4._KND+&
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
             )/4._KND
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
             )/4._KND
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


   SCAL2 = SCAL

   if (sctype==1) then
     call BOUND_Temp(SCAL2)
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
  real(KND),intent(inout) :: SCAL(-1:,-1:,-1:)
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
!       SCAL(0,j,k) = (SCAL(1,j,k)*(U(0,j,k)/2._KND-(TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))+sideScal(We))/&
!                                   (U(0,j,k)/2._KND+(TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))
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
!       SCAL(nx+1,j,k) = (SCAL(nx,j,k)*(U(nx,j,k)/2._KND+(TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))+sideScal(Ea))/&
!        (-U(nx,j,k)/2._KND+(TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))
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
!       SCAL(i,0,k) = (SCAL(i,1,k)*(V(i,0,k)/2._KND-(TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))+sideScal(So))/&
!                                   (V(i,0,k)/2._KND+(TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))
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
!       SCAL(i,ny+1,k) = (SCAL(i,ny,k)*(V(i,ny,k)/2._KND+(TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))+sideScal(No))/&
!        (-V(i,ny,k)/2._KND+(TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))
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
      SCAL(i,j,0) = SCAL(i,j,1)+sideScal(Bo)*dzW(0)/((TDiff(i,j,1)+TDiff(i,j,0))/(2._KND))
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
      SCAL(i,j,nz+1) = SCAL(i,j,nz)-sideScal(To)*dzW(nz)/((TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2._KND))
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






  subroutine BOUND_TEMP(Theta)
  real(KND),intent(inout) :: theta(-1:,-1:,-1:)
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
!       theta(0,j,k) = (theta(1,j,k)*(U(0,j,k)/2._KND-(TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))+sideTemp(We))/&
!                                   (U(0,j,k)/2._KND+(TDiff(1,j,k)+TDiff(0,j,k))/(2*dxU(0)))
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
!       theta(nx+1,j,k) = (theta(nx,j,k)*(U(nx,j,k)/2._KND+(TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))+sideTemp(Ea))/&
!        (-U(nx,j,k)/2._KND+(TDiff(nx,j,k)+TDiff(nx+1,j,k))/(2*dxU(nx)))
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
!       theta(i,0,k) = (theta(i,1,k)*(V(i,0,k)/2._KND-(TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))+sideTemp(So))/&
!                                   (V(i,0,k)/2._KND+(TDiff(i,1,k)+TDiff(i,0,k))/(2*dyV(0)))
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
!       theta(i,ny+1,k) = (theta(i,ny,k)*(V(i,ny,k)/2._KND+(TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))+sideTemp(No))/&
!        (-V(i,ny,k)/2._KND+(TDiff(i,ny,k)+TDiff(i,ny+1,k))/(2*dyV(ny)))
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
      if (abs(BsideTFLArr(i,j))<tiny(1._KND).and.TDiff(i,j,1)<1.1_KND/(Prandtl*Re).and.TBtype(Bo)==DIRICHLET) then
       theta(i,j,0) = BsideTArr(i,j)
      else
       theta(i,j,0) = theta(i,j,1)+BsideTFLArr(i,j)*dzW(0)/((TDiff(i,j,1)+TDiff(i,j,0))/(2._KND))
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
!       theta(i,j,nz+1) = (theta(i,j,nz)*(W(i,j,nz)/2._KND+(TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2*dzW(nz)))+sideTemp(To))/&
!        (-W(i,j,nz)/2._KND+(TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2*dzW(nz)))
      theta(i,j,nz+1) = theta(i,j,nz)-sideTemp(To)*dzW(nz)/((TDiff(i,j,nz)+TDiff(i,j,nz+1))/(2._KND))
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
  endsubroutine BOUND_TEMP





  pure subroutine BOUND_Visc(Nu)
  real(KND),intent(inout) :: Nu(-1:,-1:,-1:)
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















  pure real(DBL) function AirDensity(press,temp)
  real(KND),intent(in) :: press,temp
  AirDensity = press/(287.05_DBL*temp)
  endfunction AirDensity


  pure real(DBL) function AirDynVisc(temp)
  real(KND),intent(in) :: temp
  AirDynVisc = 1.85e-5_DBL
  endfunction AirDynVisc



  pure real(DBL) function CorrFactor(dp,press,temp)
  real(KND),intent(in) :: dp,press,temp
  real(DBL) l
  l = MeanFreePath(press,temp)
  CorrFactor = 1+(2*l/dp)*(1.257_KND+0.4_KND*exp(-0.55_KND*dp/l))
  endfunction CorrFactor


  pure real(DBL) function MeanFreePath(press,temp)
  real(KND),intent(in) :: press,temp
  MeanFreePath = 2.24e-5_DBL*temp/press
  endfunction MeanFreePath


  pure real(DBL) function BrownDiffusivity(dp,press,temp)
  real(KND),intent(in) :: dp,press,temp
  real(DBL) C
  C = Corrfactor(dp,press,temp)
  BrownDiffusivity = BoltzC*temp*C/(3._KND*pi*AirDynVisc(temp)*dp)
  endfunction BrownDiffusivity


  pure real(KND) function BrownEff(dp,press,temp)
  real(KND),intent(in) :: dp,press,temp
  real(KND) Sc
  Sc = (AirDynVisc(temp)/AirDensity(press,temp))/BrownDiffusivity(dp,press,temp)
  BrownEff = Sc**(-0.54_KND)
  endfunction BrownEff


  pure real(KND) function ImpactEff(dp,press,temp,ustar,visc)
  real(KND),intent(in) :: dp,press,temp,ustar,visc
  real(KND) St
  St = SedimVelocity2(dp,press,temp)*ustar**2/visc
  ImpactEff = St**2/(400+St**2)
  endfunction ImpactEff


  pure real(KND) function SedimVelocity2(dp,press,temp)
  real(KND),intent(in) :: dp,temp,press
  real(KND) C
  C = CorrFactor(dp,press,temp)
  SedimVelocity2 = AirDensity(press,temp)*dp**2*9.81_KND*C/(18._KND*AirDynVisc(temp))
  endfunction SedimVelocity2

  pure real(DBL) function SedimVelocity(dp,rhop,press,temp)
  real(KND),intent(in) :: dp,rhop,temp,press
  real(DBL) C,us,rho,mu
  rho = AirDensity(press,temp)
  mu = AirDynVisc(temp)
  C = CorrFactor(dp,press,temp)
  us = 1+(0.42_DBL*C**2*rho*rhop/(108*mu**2))*dp**3*(1-rho/rhop)*9.81_DBL
  us = sqrt(us)
  us = (12._DBL*mu/(0.42_DBL*C*rho*dp))*(us-1._DBL)
  SedimVelocity = us
  endfunction SedimVelocity


  pure real(KND) function DyerH(zL)
  real(KND),intent(in) :: zL
  if (zL>=0) then
   DyerH = 1._KND+5._KND*zl
  else
   DyerH = 1._KND/sqrt(1._KND-16._KND*zl)
  endif
  endfunction DyerH



  pure real(KND) function AerResist(z,z0,zL,ustar,visc)
  real(KND),intent(in) :: z,z0,zL,ustar,visc
  real(KND),parameter :: yplcrit = 11.225

  if (z>z0.and.z0>0) then
   AerResist = (log(z/z0)-DyerH(zl))/(0.4_KND*ustar)
  else
   if ((z*ustar/visc)<yplcrit) then
     AerResist = (z*ustar/visc)
   else
    AerResist = (log(abs(ustar*z/visc))/0.4_KND+5.2_KND)/ustar
   endif
  endif
  endfunction AerResist

  pure real(KND) function DepositionVelocity3(dp,rhop,press,temp,z,z0,zL,ustar) !EMRAS recommended values
  real(KND),intent(in) :: dp,rhop,press,temp,z,z0,zL,ustar

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



  pure real(KND) function DepositionVelocity2(dp,press,temp,z,z0,zL,ustar)
  real(KND),intent(in) :: dp,press,temp,z,z0,zL,ustar
  real(KND) SurfResist,St,visc

   visc = AirDynVisc(temp)/AirDensity(press,temp)
   St = SedimVelocity2(dp,press,temp)*ustar**2/(visc)
   SurfResist = 1._KND/(3._KND*ustar*(BrownEff(dp,press,temp)+ImpactEff(dp,press,temp,ustar,visc)))
   DepositionVelocity2 = SedimVelocity2(dp,press,temp)+1._KND/(AerResist(z,z0,zL,ustar,visc)+SurfResist)
  endfunction DepositionVelocity2


  pure real(KND) function DepositionVelocity(dp,rhop,press,temp,z,z0,zL,ustar) !Kharchenko
  real(KND),intent(in) :: dp,rhop,press,temp,z,z0,zL,ustar
  real(DBL) visc,us,Intz,Intexp,BD,tp
  real(DBL),parameter :: zexp = 0.01

   us = SedimVelocity(dp,rhop,press,temp)
   visc = AirDynVisc(temp)/AirDensity(press,temp)
   tp = (us/9.81_KND)*ustar**2/visc
   BD = BrownDiffusivity(dp,press,temp)

   if (zl>=0) then
    Intz=-log(z/zexp+6*(zl-zexp*zl/z))/Karman
   else
    Intz = ((sqrt(1-9*zl)-1)*(sqrt(1-9*zexp*zl/z)+1))
    Intz = Intz/((sqrt(1-9*zl)+1)*(sqrt(1-9*zexp*zl/z)-1))
    Intz=-Intz/Karman
   endif

   Intexp=-367.8_KND
   Intexp = Intexp+16.4*log(visc/BD)
   Intexp = Intexp-0.73*log(100*dp)*log(1e4*BD)-0.5*(log(100*dp))**2
   Intexp = Intexp+0.13*log(0.03/z0)
   Intexp = Intexp+0.25*log(0.2/ustar)*(1-0.2*log(0.03/z0))
   Intexp = Intexp-0.03*log(tp)*log(0.03/z0)
   Intexp = Intexp-32.7*log(100*dp)
   Intexp=-exp(Intexp)

   DepositionVelocity = us/(1._KND-exp((us/ustar)*(Intexp+Intz)))
  endfunction DepositionVelocity


  pure real(KND) function DepositionFlux(WMP,conc,partdiam,rhop)
  type(WMPoint),intent(in) :: WMP
  real(KND),intent(in) :: conc,partdiam,rhop
  real(KND) :: press,temp,depvel

   press = 101300
   temp = temperature_ref
   if ((.2_KND*WMP%distz)**2>(WMP%distx)**2+(WMP%disty)**2.and.WMP%distz>0) then
    depvel = DepositionVelocity3(partdiam,rhop,press,temp,WMP%distz,WMP%z0,0._KND,WMP%ustar)
   else
    depvel = 0
   endif
   DepositionFlux = abs(depvel)*abs(conc)

  endfunction DepositionFlux

  subroutine Deposition(SCAL,coef)
  real(KND),dimension(-1:,-1:,-1:,1:),intent(inout) :: SCAL
  real(KND),intent(in) :: coef
  type(WMPoint),pointer :: WMP
  integer i
  real(KND) deptmp
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
     do i = 1,computescalars
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
  real(KND),dimension(-1:,-1:,-1:,1:) :: SCAL
  integer i,j,k,l
  real(KND),dimension(Prnx,Prny,Prnz) :: flux
  real(KND) :: coef,press,temp,us

  press = 101300
  temp = temperature_ref
  if (partdistrib==0) then
   do l = 1,computescalars
    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       us = SedimVelocity(partdiam(l),partrho(l),press,temp)
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


  pure real(KND) function Rig(i,j,k,U,V,temperature)
  integer,intent(in) :: i,j,k
  real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V
  real(KND),dimension(-1:,-1:,-1:),intent(in) :: temperature
  real(KND) num,denom

  num = (grav_acc/temperature_ref)*(temperature(i,j,k+1)-temperature(i,j,k-1))/(zPr(k+1)-zPr(k-1))
  denom = ((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(2._KND*(zPr(k+1)-zPr(k-1))))**2
  denom = denom+((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(zPr(k+1)-zPr(k-1))))**2

  if (abs(denom)>1E-5_KND*abs(num)) then
   Rig = num/denom
  else
   Rig = 100000._KND*sign(1.0_KND,num)*sign(1.0_KND,denom)
  endif
  endfunction Rig



  subroutine ComputeTDiff(U,V,W)
  real(KND),intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer:: i,j,k
  real(KND),parameter :: Prt = constPrt !if variable, then implement as a statement function due to problems with inlining

   if (Re>0) then
    !$omp workshare
    forall(k=1:Prnz,j=1:Prny,i=1:Prnx)&
          TDiff(i,j,k) = (Visc(i,j,k)-1._KND/Re)/Prt + (1._KND/(Re*Prandtl))
    !$omp end workshare
   else
    !$omp workshare
    forall(k=1:Prnz,j=1:Prny,i=1:Prnx)&
          TDiff(i,j,k) = Visc(i,j,k)/Prt
    !$omp end workshare
   endif
  end subroutine ComputeTDiff



  subroutine InitTemperatureProfile(TempIn)
    real(KND),intent(out) :: TempIn(-1:,-1:)
    integer   :: SectionToUse(-1:ubound(TempIn,2))
    integer   :: section,nSections,s
    integer   :: i,j,k
    real(KND) :: temp

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
    real(KND),intent(out) :: TempIn(-1:,-1:)
    real(KND),intent(out) :: Temperature(-1:,-1:,-1:)
    real(KND) :: p_x,p_y,p_z,p
    integer   :: i,j,k

    if (TemperatureProfile%randomize==1) then

!       p_x=0.5_KND
!       p_y=0.5_KND
!       p_z=0.5_KND

      do k=0,Prnz+1

!         if (zPr(k) <= TemperatureProfile%randomizeTop) then
!           call RANdoM_NUMBER(p)
!           p = p - 0.5
!           p_z=(p_z+p)/2
!         endif

        do j=0,Prny+1

!           if (zPr(k) <= TemperatureProfile%randomizeTop) then
!             call RANdoM_NUMBER(p)
!             p = p - 0.5
!             p_y=(p_y+p)/2
!           endif

           do i=0,Prnx+1

             if (zPr(k) <= TemperatureProfile%randomizeTop) then
               call RANdoM_NUMBER(p)
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

    allocate(SubsidenceProfile(0:Prnz))
    SubsidenceProfile = (/ (zW(k)*SubsidenceGradient, k=0,Prnz) /)
  end subroutine InitSubsidenceProfile



  subroutine VolScalSource(Scal,coef)
  real(KND),intent(inout) :: Scal(-1:,-1:,-1:,:)
  real(KND),intent(in) :: coef

#include "customvolscalsource.f90"
  end subroutine VolScalSource



!  !$hmpp <tsteps> ScalarConvection codelet
 subroutine ScalarConvection(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,buoyancy,convmet,&
                       dxmin,dymin,dzmin,coriolisparam,grav_acc,temperature_ref,&
                       U,V,W,U2,V2,W2,Ustar,Vstar,Wstar,temperature,beta,rho,lev,dt)
  implicit none
! #ifdef __HMPP
!   integer,parameter :: KND = 4
  intrinsic abs
! #endif

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




  end subroutine ScalarConvection


  subroutine Release(Scalar,released)
    real(KND),intent(inout) :: Scalar(-1:,-1:,-1:,-1:)
    logical,intent(inout)   :: released
    real(KND) xc,yc,xs,xf,ys,yf,zs,zf,dxp,dyp,dzp,ct,cr,xp,yp,zp,p
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
    if (computescalars>=4) then
     if (time>(endtime-starttime)/3._KND) then
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

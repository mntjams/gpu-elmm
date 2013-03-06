module HMPP_SCALARS
  use HMPP_CODELETS
  use PARAMETERS

  implicit none

  private
  public HMPP_RKstage_Temperature,GetTemperatureFromGPU

  contains


!     subroutine HMPP_KappaTemperature(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
!                             dxmin,dymin,dzmin,limparam,&
!                             Temperature2,Temperature,U,V,W,coef,dt,SubsidenceProfile,fluxProfile)
!
!       integer,intent(in)         :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz
!       real(knd),intent(out)      :: Temperature2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
!       real(knd),intent(in)       :: Temperature(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
!       real(knd),intent(in)       :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3),V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
!       real(knd),intent(in)       :: dxmin,dymin,dzmin,coef,dt,limparam
!       real(knd),intent(in)       :: SubsidenceProfile(0:)
!       real(knd),intent(out)      :: fluxProfile(0:)
!       real(knd),allocatable,save :: SubsidenceProfileLoc(:),fluxProfileLoc(:)
!       integer,save :: called = 0
!
!       if (called==0) then
!         allocate(SubsidenceProfileLoc(0:Prnz))
!         allocate(fluxProfileLoc(0:Prnz))
!         called = 1
!       !$hmpp <tsteps> advancedload, args[KappaTemperature::limparam]
!       end if
!
!
!       !$hmpp <tsteps> advancedload, args[KappaTemperature::coef]
!
!
!       if (size(SubsidenceProfile)==size(SubsidenceProfileLoc)) then
!         SubsidenceProfileLoc = SubsidenceProfile
!       else
!         SubsidenceProfileLoc = 0
!       end if
!
!       !$hmpp <tsteps> advancedload, args[KappaTemperature::Temperature]
!       !$hmpp <tsteps> advancedload, args[KappaTemperature::SubsidenceProfile]
!
!       !$hmpp <tsteps> KappaTemperature callsite, args[*].noupdate=true
!       call KappaTemperature_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
!                             dxmin,dymin,dzmin,limparam,&
!                             Temperature2,Temperature,U,V,W,&
!                             coef,dt,SubsidenceProfileLoc,fluxProfileLoc)
!
!       !$hmpp <tsteps> delegatedstore, args[KappaTemperature::Temperature2]
!       !$hmpp <tsteps> delegatedstore, args[KappaTemperature::fluxProfile]
!
!       if (size(fluxProfile)==size(fluxProfileLoc)) fluxProfile(0:) = fluxProfileLoc
!
!     end subroutine HMPP_KappaTemperature
!
!
!     subroutine HMPP_DiffTemperature(Prnx,Prny,Prnz,dxmin,dymin,dzmin,maxCNiter,epsCN,Re,&
!                               TBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr,&
!                               TDiff,Temperature2,Temperature,&
!                               coef,dt)
!
!       integer,intent(in)    :: Prnx,Prny,Prnz,maxCNiter,TBtype(6)
!       real(knd),intent(out),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2) :: Temperature2
!       real(knd),intent(in),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)  :: Temperature
!       real(knd),intent(in)  :: dxmin,dymin,dzmin,epsCN,Re,coef,dt
!       real(knd),intent(in)  :: sideTemp(6),Tempin(-1:Prny+2,-1:Prnz+2),TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
!       real(knd),intent(in),dimension(-1:,-1:) :: BsideTArr,BsideTFLArr
!       real(knd),dimension(:,:),allocatable,save :: BsideTArrLoc,BsideTFLArrLoc
!       integer   :: l
!       real(knd) :: res
!       integer,save :: called = 0
!
!
!
!       if (called==0) then
!         allocate(BsideTArrLoc(-1:Prnx+2,-1:Prny+2),BsideTFLArrLoc(-1:Prnx+2,-1:Prny+2))
!         !$hmpp <tsteps> advancedload, args[DiffTemperature::TBtype,DiffTemperature::sideTemp,DiffTemperature::TempIn]
!         !$hmpp <tsteps> advancedload, args[DiffTemperature::epsCN,DiffTemperature::maxCNiter]
!         called = 1
!       end if
!
!
!       !$hmpp <tsteps> advancedload, args[DiffTemperature::coef]
!       if (size(BsideTArr)==size(BsideTArrLoc)) then
!         BsideTArrLoc = BsideTArr
!         !$hmpp <tsteps> advancedload, args[DiffTemperature::BsideTArr]
!       endif
!       if (size(BsideTFLArr)==size(BsideTFlArrLoc)) then
!         BsideTFlArrloc = BsideTFlArr
!         !$hmpp <tsteps> advancedload, args[DiffTemperature::BsideTFLArr]
!       endif
!
!       !$hmpp <tsteps> advancedload, args[DiffTemperature::TDiff,DiffTemperature::Temperature]
!
!       !$hmpp <tsteps> DiffTemperature callsite, args[*].noupdate=true
!       call DiffTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,maxCNiter,epsCN,Re,&
!                             TBtype,sideTemp,TempIn,BsideTArrLoc,BsideTFLArrLoc,&
!                             TDiff,Temperature2,Temperature,&
!                             coef,dt,l,res)
!
!       !$hmpp <tsteps> delegatedstore, args[DiffTemperature::l,DiffTemperature::res]
!       !$hmpp <tsteps> delegatedstore, args[DiffTemperature::Temperature2]
!
!       write (*,*) "CNscalar ",l,res
!     end subroutine HMPP_DiffTemperature


    subroutine HMPP_RKstage_Temperature(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                            dxmin,dymin,dzmin,maxCNiter,epsCN,Re,limparam,&
                            TBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr,&
                            TDiff,Temperature2,Temperature,Temperature_adv,U,V,W,&
                            SubsidenceProfile,fluxProfile,dt,RKstage,alpha,beta,rho)

      integer,intent(in)      :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,maxCNiter,TBtype(6),RKstage
      real(knd),intent(inout) :: Temperature2(-1:,-1:,-1:)  !Only to avoid allocation of an automatic array, memory transfers to codelet should be avoided. Assuming we do not have to save main memory.
      real(knd),intent(inout) :: Temperature(-1:,-1:,-1:)
      real(knd),intent(inout) :: Temperature_adv(-1:,-1:,-1:)
      real(knd),intent(in)    :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
      real(knd),intent(in)    :: dxmin,dymin,dzmin,epsCN,Re,limparam,dt
      real(knd),intent(in)    :: sideTemp(6),Tempin(-1:,-1:),TDiff(-1:,-1:,-1:)
      real(knd),intent(in)    :: BsideTArr(-1:,-1:),BsideTFLArr(-1:,-1:)
      real(knd),intent(in)    :: SubsidenceProfile(0:)
      real(knd),intent(out)   :: fluxProfile(0:)
      real(knd),intent(in)    :: alpha(3),beta(3),rho(3)

      integer   :: CNiters
      real(knd) :: CNres

      real(knd),allocatable,save :: SubsidenceProfileLoc(:),fluxProfileLoc(:)
      real(knd),allocatable,save :: BsideTArrLoc(:,:),BsideTFLArrLoc(:,:)
      integer,save               :: called = 0

      if (called==0) then
        allocate(SubsidenceProfileLoc(0:Prnz))
        allocate(fluxProfileLoc(0:Prnz))
        fluxProfileLoc = 0
        allocate(BsideTArrLoc(-1:Prnx+2,-1:Prny+2),BsideTFLArrLoc(-1:Prnx+2,-1:Prny+2))
       !$hmpp <tsteps> advancedload, args[RKstage_Temperature::limparam]
       !$hmpp <tsteps> advancedload, args[RKstage_Temperature::alpha,RKstage_Temperature::beta,RKstage_Temperature::rho]
       !$hmpp <tsteps> advancedload, args[RKstage_Temperature::TBtype,RKstage_Temperature::sideTemp,RKstage_Temperature::TempIn]
       !$hmpp <tsteps> advancedload, args[RKstage_Temperature::epsCN,RKstage_Temperature::maxCNiter]

       if (size(SubsidenceProfile)==size(SubsidenceProfileLoc)) then
          SubsidenceProfileLoc = SubsidenceProfile
        else
          SubsidenceProfileLoc = 0
        end if

       !$hmpp <tsteps> advancedload, args[RKstage_Temperature::SubsidenceProfile]

        called = 1
      end if

      if (size(BsideTArr)==size(BsideTArrLoc)) then
        BsideTArrLoc = BsideTArr
        !$hmpp <tsteps> advancedload, args[RKstage_Temperature::BsideTArr]
      end if

      if (size(BsideTFLArr)==size(BsideTFlArrLoc)) then
        BsideTFlArrloc = BsideTFlArr
        !$hmpp <tsteps> advancedload, args[RKstage_Temperature::BsideTFLArr]
      end if

      !$hmpp <tsteps> advancedload, args[RKstage_Temperature::RKstage]

      !No need to upload Temperature, it was uploaded for momentum solution
      !No need to upload Temperature2, it is needed only inside
      !No need to upload Temperature_avg, it is needed only inside
!
!       !$hmpp <tsteps> advancedload, args[RKstage_Temperature::TDiff]
!       !$hmpp <tsteps> advancedload, args[RKstage_Temperature::Temperature]

      !$hmpp <tsteps> RKstage_Temperature callsite, args[*].noupdate=true
      call RKstage_Temperature_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                            dxmin,dymin,dzmin,maxCNiter,epsCN,Re,limparam,&
                            TBtype,sideTemp,TempIn,BsideTArrLoc,BsideTFLArrLoc,&
                            TDiff,Temperature2,Temperature,Temperature_adv,U,V,W,&
                            SubsidenceProfileLoc,fluxProfileLoc,dt,RKstage,alpha,beta,rho,CNiters,CNres)

      !No need to download Temperature output procedures can download eplicitly.
      !No need to download Temperature2, it is needed only inside
      !No need to download Temperature_avg, it is needed only inside

      if (RKstage==3) then
!         !$hmpp <tsteps> delegatedstore, args[RKstage_Temperature::Temperature]
        !$hmpp <tsteps> delegatedstore, args[RKstage_Temperature::fluxProfile]
        if (size(fluxProfile)==size(fluxProfileLoc)) fluxProfile(0:) = fluxProfileLoc
      end if

      !$hmpp <tsteps> delegatedstore, args[RKstage_Temperature::CNiters,RKstage_Temperature::CNres]


      write (*,*) "CNscalar ",CNiters,CNres
    end subroutine HMPP_RKstage_Temperature


  subroutine GetTemperatureFromGPU(Temperature)
    real(knd),intent(inout) :: Temperature(-1:,-1:,-1:)
    !$hmpp <tsteps> delegatedstore, args[RKstage_Temperature::Temperature]
  end subroutine GetTemperatureFromGPU


end module HMPP_SCALARS
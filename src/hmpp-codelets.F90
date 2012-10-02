module HMPP_CODELETS
  use parameters
  implicit none


contains
!GPU codelets

!$hmpp <tsteps> group, target = CUDA

  !$hmpp <tsteps> mapbyname, Prnx,Prny,Prnz
  !$hmpp <tsteps> mapbyname, Btype,sideU,Re,sideTemp,TBtype,TempIn,BsideTArr,BsideTFLArr
  !$hmpp <tsteps> mapbyname, dt,dt2
  !$hmpp <tsteps> mapbyname, dxPr,dyPr,dzPr,dxU,dyV,dzW
  !$hmpp <tsteps> mapbyname, xPr,zPr,xU,zW
  !$hmpp <tsteps> mapbyname, Pr,Visc,TDiff

  !$hmpp <tsteps> mapbyname, buoyancy

  !$hmpp <tsteps> mapbyname, temperature_ref,grav_acc

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
  !$hmpp &   AttenuateOut2::U,NullInterior::U,ComputeViscWM::U,RKstage_Temperature::U,UpdateU::U]
 !$hmpp <tsteps> map, args[Vreman::V,Convection::V,ForwEul::V,UnifRedBlack::V,TimeStepEul::V,&
  !$hmpp  &  AttenuateOut2::V,NullInterior::V,ComputeViscWM::V,RKstage_Temperature::V,UpdateU::V]
 !$hmpp <tsteps> map, args[Vreman::W,Convection::W,ForwEul::W,UnifRedBlack::W,TimeStepEul::W,&
  !$hmpp  &  AttenuateOut2::W,NullInterior::W,ComputeViscWM::W,RKstage_Temperature::W,UpdateU::W]


 !U2,V2,W2
 !$hmpp <tsteps> map, args[Convection::U2,PressureGrad::U,ForwEul::U2,UnifRedBlack::U2,AttenuateTop::U,&
   !$hmpp &   AttenuateOut::U,PrePoisson::U,PostPoisson::U,UpdateU::U2]
 !$hmpp <tsteps> map, args[Convection::V2,PressureGrad::V,ForwEul::V2,UnifRedBlack::V2,AttenuateTop::V,&
   !$hmpp &   AttenuateOut::V,PrePoisson::V,PostPoisson::V,UpdateU::V2]
 !$hmpp <tsteps> map, args[Convection::W2,PressureGrad::W,ForwEul::W2,UnifRedBlack::W2,AttenuateTop::W,&
   !$hmpp &   AttenuateOut::W,PrePoisson::W,PostPoisson::W,UpdateU::W2]

 !U3,V3,W3 on GPU device mapped also to Ustar,Vstar,Wstar
 !$hmpp <tsteps> map, args[ForwEul::U3,UnifRedBlack::U3]
 !$hmpp <tsteps> map, args[ForwEul::V3,UnifRedBlack::V3]
 !$hmpp <tsteps> map, args[ForwEul::W3,UnifRedBlack::W3]

 !$hmpp <tsteps> map, args[Convection::temperature,AttenuateTop::temperature,AttenuateOut::temperature,&
 !$hmpp &  AttenuateOut2::temperature,RKstage_Temperature::Temperature]
!,BoundTemperature::theta

!  !$hmpp <tsteps> map, args[PrePoisson::RHS,PostPoisson::Phi]

  !$hmpp <tsteps> map, args[UnifRedBlack::maxCNiter]
  !,DiffTemperature::maxCNiter,DiffScalar::maxCNiter]
  !$hmpp <tsteps> map, args[UnifRedBlack::epsCN]
  !,DiffTemperature::epsCN,DiffScalar::epsCN]

#define myblocksize 512x1

#define myblocksize2 32x16

#include "tsteps-shared.f90"

#include "boundaries_GPU.f90"

#include "tsteps_GPU.f90"

#include "cds_GPU.f90"

#include "smagorinsky_GPU.f90"

#include "scalars_GPU.f90"

#include "scalarsbound_GPU.f90"

#include "wallmodels_GPU.F90"

#include "pressure_GPU.F90"

end module HMPP_CODELETS
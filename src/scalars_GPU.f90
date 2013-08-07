
  pure function FluxLimiter(r,limparam) result(res)
    implicit none
#include "hmpp-include.f90"
    real(knd) res
    real(knd),intent(in) :: r,limparam
    intrinsic max,min
    res = max(0._knd,min(2._knd*r,min(limparam,(1+2._knd*r)/3._knd)))
  end function FluxLimiter

  pure function hsign(r) result(res)  !sign is not supported by HMPP yet
    implicit none
#include "hmpp-include.f90"
    real(knd) res
    real(knd),intent(in) :: r

    if (r>0) then
      res = 1._knd
    else if (r<0) then
      res = -1._knd
    else
      res = 0._knd
    end if
  end function hsign


!   !$hmpp <tsteps> KappaTemperature codelet
  subroutine KappaTemperature_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                            dxmin,dymin,dzmin,limparam,&
                            Temperature2,Temperature,U,V,W,&
                            coef,dt,SubsidenceProfile,fluxProfile) !Kappa scheme with flux limiter
  implicit none                                                                                   !Hunsdorfer et al. 1995, JCP
#include "hmpp-include.f90"

  integer,intent(in)    :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz
  real(knd),intent(in)  :: Temperature(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  real(knd),intent(out) :: Temperature2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  real(knd),intent(in)  :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3),V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(knd),intent(in)  :: dxmin,dymin,dzmin,coef,dt,limparam
  real(knd),intent(in)  :: SubsidenceProfile(0:Prnz)
  real(knd),intent(out) :: fluxProfile(0:Prnz)
  integer i,j,k,l
  real(knd) A,Ax,Ay,Az              !Auxiliary variables to store muliplication constants for efficiency
  real(knd) vel,SL,SR,SP,FLUX
  real(knd),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2) :: SLOPE
  real(knd),parameter ::eps = 1e-8
  intrinsic sign


  A = coef*dt
  Ax = coef*dt/dxmin
  Ay = coef*dt/dymin
  Az = coef*dt/dzmin

  !$hmppcg grid blocksize myblocksize
  !$hmppcg gridify(k,i), private(j)
  do k = -1,Prnz+2
   do i = -1,Prnx+2
    do j = -1,Prny+2
      Temperature2(i,j,k) = 0
    enddo
   enddo
  enddo

  !$hmppcg grid blocksize myblocksize
  !$hmppcg gridify(k,i), private(j,SL,SR,SP)
  do k = 1,Prnz
   do i = 0,Prnx
    do j = 1,Prny
     if (U(i,j,k)>0) then
      SR = (Temperature(i+1,j,k)-Temperature(i,j,k))
      SL = (Temperature(i,j,k)-Temperature(i-1,j,k))
     else
      SR = (Temperature(i,j,k)-Temperature(i+1,j,k))
      SL = (Temperature(i+1,j,k)-Temperature(i+2,j,k))
     endif
     SP = (SR+eps*hsign(SL))/(SL+eps*hsign(SL))
     SLOPE(i,j,k) = FLUXLIMITER(SP,limparam)
    enddo
   enddo
  enddo


  do l=0,1  !odd-even separation to avoid a race condition
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify(k,i), private(j,FLUX)
    do k = 1,Prnz
     do i = 0+l,Prnx,2
      do j = 1,Prny
       if (U(i,j,k)>0) then
        FLUX = U(i,j,k)*(Temperature(i,j,k)+(Temperature(i,j,k)-Temperature(i-1,j,k))*SLOPE(i,j,k)/2._knd)
       else
        FLUX = U(i,j,k)*(Temperature(i+1,j,k)+(Temperature(i+1,j,k)-Temperature(i+2,j,k))*SLOPE(i,j,k)/2._knd)
       endif
       Temperature2(i,j,k) = Temperature2(i,j,k)-Ax*FLUX
       Temperature2(i+1,j,k) = Temperature2(i+1,j,k)+Ax*FLUX
      enddo
     enddo
    enddo
  enddo


  !$hmppcg grid blocksize myblocksize
  !$hmppcg gridify(k,i), private(j,SL,SR,SP)
  do k = 1,Prnz
   do i = 1,Prnx
    do j = 0,Prny
     if (V(i,j,k)>0) then
      SR = (Temperature(i,j+1,k)-Temperature(i,j,k))
      SL = (Temperature(i,j,k)-Temperature(i,j-1,k))
     else
      SR = (Temperature(i,j,k)-Temperature(i,j+1,k))
      SL = (Temperature(i,j+1,k)-Temperature(i,j+2,k))
     endif
     SP = (SR+eps*hsign(SL))/(SL+eps*hsign(SL))
     SLOPE(i,j,k) = FLUXLIMITER(SP,limparam)
    enddo
   enddo
  enddo


  do l=0,1  !odd-even separation to avoid a race condition
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify(k,i), private(j,FLUX)
    do k = 1,Prnz
     do i = 1,Prnx
      do j = 0+l,Prny,2
       if (V(i,j,k)>0) then
        FLUX = V(i,j,k)*(Temperature(i,j,k)+(Temperature(i,j,k)-Temperature(i,j-1,k))*SLOPE(i,j,k)/2._knd)
       else
        FLUX = V(i,j,k)*(Temperature(i,j+1,k)+(Temperature(i,j+1,k)-Temperature(i,j+2,k))*SLOPE(i,j,k)/2._knd)
       endif

       Temperature2(i,j,k) = Temperature2(i,j,k)-Ay*FLUX
       Temperature2(i,j+1,k) = Temperature2(i,j+1,k)+Ay*FLUX
      enddo
     enddo
    enddo
  enddo


  !$hmppcg grid blocksize myblocksize
  !$hmppcg gridify(k,i), private(j,SL,SR,SP)
  do k = 0,Prnz
   do i = 1,Prnx
    do j = 1,Prny
     if (W(i,j,k)>0) then
      SR = (Temperature(i,j,k+1)-Temperature(i,j,k))
      SL = (Temperature(i,j,k)-Temperature(i,j,k-1))
     else
      SR = (Temperature(i,j,k)-Temperature(i,j,k+1))
      SL = (Temperature(i,j,k+1)-Temperature(i,j,k+2))
     endif
     SP = (SR+eps*hsign(SL))/(SL+eps*hsign(SL))
     SLOPE(i,j,k) = FLUXLIMITER(SP,limparam)
    enddo
   enddo
  enddo


  do l=0,1  !odd-even separation to avoid a race condition
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify(k,i), private(j,FLUX,vel)
    do k = 0+l,Prnz,2
     do i = 1,Prnx
      do j = 1,Prny

       vel = W(i,j,k) - SubsidenceProfile(k)

       if (vel>0) then
        FLUX = vel*(Temperature(i,j,k)+(Temperature(i,j,k)-Temperature(i,j,k-1))*SLOPE(i,j,k)/2._knd)
       else
        FLUX = vel*(Temperature(i,j,k+1)+(Temperature(i,j,k+1)-Temperature(i,j,k+2))*SLOPE(i,j,k)/2._knd)
       endif

       if (abs(vel)>=1e-6) then
         fluxprofile(k) = fluxprofile(k) + FLUX/vel*W(i,j,k)
       else
         fluxprofile(k) = 0
       endif

       Temperature2(i,j,k) = Temperature2(i,j,k) - Az*FLUX
       Temperature2(i,j,k+1) = Temperature2(i,j,k+1) + Az*FLUX
      enddo
     enddo
    enddo
  end do

  !$hmppcg gridify(k)
  do k = 0,Prnz
      fluxprofile(k) = fluxprofile(k) / (Prnx*Prny)
  enddo

  endsubroutine KappaTemperature_GPU





!   !$hmpp <tsteps> DiffTemperature codelet
  subroutine DiffTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,maxCNiter,epsCN,Re,&
                            TempBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr,&
                            TDiff,Temperature2,Temperature,&
                            coef,dt,l,res)
  implicit none
#include "hmpp-include.f90"

  integer,intent(in)    :: Prnx,Prny,Prnz,maxCNiter,TempBtype(6)
  real(knd),intent(out),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2) :: Temperature2
  real(knd),intent(in),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)  :: Temperature

  real(knd),intent(in)  :: dxmin,dymin,dzmin,epsCN,Re,coef,dt
  real(knd),intent(in)  :: sideTemp(6),TempIn(-1:Prny+2,-1:Prnz+2),TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  real(knd),intent(in),dimension(-1:Prnx+2,-1:Prny+2) :: BsideTArr,BsideTFLArr
  integer,intent(out)   :: l
  real(knd),intent(out) :: res
  real(knd) Temperature3(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  integer i,j,k,xi,yj,zk
  real(knd) p
  real(knd) A,Ax,Ay,Az,Ap(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)

  intrinsic max, abs, mod


  Ax = 1._knd/(dxmin**2)
  Ay = 1._knd/(dymin**2)
  Az = 1._knd/(dzmin**2)

  if (Re>0) then

  !$hmppcg grid blocksize myblocksize
   !$hmppcg gridify(k,i)
   do k = 1,Prnz  !initital value using forward Euler
    do i = 1,Prnx
     do j = 1,Prny
      Temperature3(i,j,k) = ((TDiff(i+1,j,k)+TDiff(i,j,k))*(Temperature(i+1,j,k)-Temperature(i,j,k))-&
        (TDiff(i,j,k)+TDiff(i-1,j,k))*(Temperature(i,j,k)-Temperature(i-1,j,k)))*Ax

      Temperature3(i,j,k) = Temperature3(i,j,k)+ &
        ((TDiff(i,j+1,k)+TDiff(i,j,k))*(Temperature(i,j+1,k)-Temperature(i,j,k))-&
        (TDiff(i,j,k)+TDiff(i,j-1,k))*(Temperature(i,j,k)-Temperature(i,j-1,k)))*Ay

      Temperature3(i,j,k) = Temperature3(i,j,k)+ &
        ((TDiff(i,j,k+1)+TDiff(i,j,k))*(Temperature(i,j,k+1)-Temperature(i,j,k))-&
        (TDiff(i,j,k)+TDiff(i,j,k-1))*(Temperature(i,j,k)-Temperature(i,j,k-1)))*Az
     enddo
    enddo
   enddo

   A = dt*coef
  !$hmppcg grid blocksize myblocksize
   !$hmppcg gridify(k,i)
   do k = 1,Prnz
    do i = 1,Prnx
     do j = 1,Prny
       Temperature2(i,j,k) = Temperature(i,j,k)+A*Temperature3(i,j,k)
     enddo
    enddo
   enddo

   call BoundTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,TempBtype,sideTemp,BsideTArr,BsideTFLArr,TDiff,TempIn,Temperature2)



   Ax = 1._knd/(4._knd*dxmin**2)
   Ay = 1._knd/(4._knd*dymin**2)
   Az = 1._knd/(4._knd*dzmin**2)

  !$hmppcg grid blocksize myblocksize
   !$hmppcg gridify(k,i)
   do k = 1,Prnz
    do i = 1,Prnx
     do j = 1,Prny
      Ap(i,j,k) = 1._knd/(1._knd/A+(((TDiff(i+1,j,k)+TDiff(i,j,k))+&
                           (TDiff(i,j,k)+TDiff(i-1,j,k)))*Ax+&
                           ((TDiff(i,j+1,k)+TDiff(i,j,k))+&
                           (TDiff(i,j,k)+TDiff(i,j-1,k)))*Ay+&
                           ((TDiff(i,j,k+1)+TDiff(i,j,k))+&
                           (TDiff(i,j,k)+TDiff(i,j,k-1)))*Az))
     enddo
    enddo
   enddo

   l=0
   res = epsCN + 1._knd
   do while (l<maxCNiter)!.and.res>epsCN
    l=l+1
!     res = 0

    call BoundTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,TempBtype,sideTemp,BsideTArr,BsideTFLArr,TDiff,TempIn,Temperature2)


   !$hmppcg grid blocksize myblocksize
   !$hmppcg gridify(k,i), private(j,p)
   !, reduce(max:res)
    do k = 1,Prnz
     do i = 1,Prnx
      do j = 1+mod(i+k,2),Prny,2
        p = (Temperature(i,j,k)/A)+(Temperature3(i,j,k)/4._knd+&
         ((TDiff(i+1,j,k)+TDiff(i,j,k))*(Temperature2(i+1,j,k))-&
          (TDiff(i,j,k)+TDiff(i-1,j,k))*(-Temperature2(i-1,j,k)))*Ax+&
         ((TDiff(i,j+1,k)+TDiff(i,j,k))*(Temperature2(i,j+1,k))-&
          (TDiff(i,j,k)+TDiff(i,j-1,k))*(-Temperature2(i,j-1,k)))*Ay+&
         ((TDiff(i,j,k+1)+TDiff(i,j,k))*(Temperature2(i,j,k+1))-&
          (TDiff(i,j,k)+TDiff(i,j,k-1))*(-Temperature2(i,j,k-1)))*Az&
         )
         p = p*Ap(i,j,k)
!          res = max(res,abs(p-Temperature2(i,j,k)))
         Temperature2(i,j,k) = p
      enddo
     enddo
    enddo

   !$hmppcg grid blocksize myblocksize
   !$hmppcg gridify(k,i), private(j,p)
   !, reduce(max:res)
    do k = 1,Prnz
     do i = 1,Prnx
      do j = 1+mod(i+k+1,2),Prny,2
        p = (Temperature(i,j,k)/A)+(Temperature3(i,j,k)/4._knd+&
         ((TDiff(i+1,j,k)+TDiff(i,j,k))*(Temperature2(i+1,j,k))-&
          (TDiff(i,j,k)+TDiff(i-1,j,k))*(-Temperature2(i-1,j,k)))*Ax+&
         ((TDiff(i,j+1,k)+TDiff(i,j,k))*(Temperature2(i,j+1,k))-&
          (TDiff(i,j,k)+TDiff(i,j-1,k))*(-Temperature2(i,j-1,k)))*Ay+&
         ((TDiff(i,j,k+1)+TDiff(i,j,k))*(Temperature2(i,j,k+1))-&
          (TDiff(i,j,k)+TDiff(i,j,k-1))*(-Temperature2(i,j,k-1)))*Az&
         )
         p = p*Ap(i,j,k)
!          res = max(res,abs(p-Temperature2(i,j,k)))
         Temperature2(i,j,k) = p
      enddo
     enddo
    enddo

   enddo


  else

   !$hmppcg grid blocksize myblocksize
   !$hmppcg gridify(k,i)
   do k = 1,Prnz
    do i = 1,Prnx
     do j = 1,Prny
       Temperature2(i,j,k) = Temperature(i,j,k)
     enddo
    enddo
   enddo

   call BoundTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,TempBtype,sideTemp,BsideTArr,BsideTFLArr,TDiff,TempIn,Temperature2)

  endif
  endsubroutine DiffTemperature_GPU




  !$hmpp <tsteps> RKstage_Temperature codelet
  subroutine RKstage_Temperature_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                            dxmin,dymin,dzmin,maxCNiter,epsCN,Re,limparam,&
                            TempBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr,&
                            TDiff,Temperature2,Temperature,Temperature_adv,U,V,W,&
                            SubsidenceProfile,fluxProfile,dt,RKstage,alpha,beta,rho,CNiters,CNres)
    implicit none                                                                                   !Hunsdorfer et al. 1995, JCP
#include "hmpp-include.f90"

    integer,intent(in)      :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,maxCNiter,TempBtype(6),RKstage
    real(knd),intent(inout) :: Temperature(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
    real(knd),intent(inout) :: Temperature_adv(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
    real(knd),intent(in)    :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3),V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
    real(knd),intent(in)    :: dxmin,dymin,dzmin,epsCN,Re,limparam,dt
    real(knd),intent(in)    :: sideTemp(6),TempIn(-1:Prny+2,-1:Prnz+2),TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
    real(knd),intent(in),dimension(-1:Prnx+2,-1:Prny+2) :: BsideTArr,BsideTFLArr
    real(knd),intent(in)    :: SubsidenceProfile(0:Prnz)
    real(knd),intent(out)   :: fluxProfile(0:Prnz)
    real(knd),intent(in)    :: alpha(3),beta(3),rho(3)
    integer,intent(out)     :: CNiters
    real(knd),intent(out)   :: CNres
    real(knd),intent(out)   :: Temperature2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)

    integer i,j,k

      call BoundTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,TempBtype,sideTemp,BsideTArr,BsideTFLArr,TDiff,TempIn,temperature)

      if (RKstage>1) then

        !$hmppcg grid blocksize myblocksize
        !$hmppcg gridify(k,i), private(j)
        do k=-1,Prnz+2
         do i=-1,Prnx+2
          do j=-1,Prny+2
            temperature2(i,j,k) = temperature_adv(i,j,k)*rho(RKstage)
          enddo
         enddo
        enddo

      else

        !$hmppcg grid blocksize myblocksize
        !$hmppcg gridify(k,i), private(j)
        do k=-1,Prnz+2
         do i=-1,Prnx+2
          do j=-1,Prny+2
            temperature2(i,j,k)=0
          enddo
         enddo
        enddo

      endif


      call KappaTemperature_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                            dxmin,dymin,dzmin,limparam,&
                            temperature_adv,temperature,U,V,W,&
                            1._knd,dt,SubsidenceProfile,fluxProfile)



      !$hmppcg grid blocksize myblocksize
      !$hmppcg gridify(k,i), private(j)
      do k=-1,Prnz+2
       do i=-1,Prnx+2
        do j=-1,Prny+2
          temperature2(i,j,k) = temperature2(i,j,k)+temperature_adv(i,j,k)*beta(RKstage)
        enddo
       enddo
      enddo

      !$hmppcg grid blocksize myblocksize
      !$hmppcg gridify(k,i), private(j)
      do k=-1,Prnz+2
       do i=-1,Prnx+2
        do j=-1,Prny+2
          temperature(i,j,k) = temperature(i,j,k)+temperature2(i,j,k)
        enddo
       enddo
      enddo


      call BoundTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,TempBtype,sideTemp,BsideTArr,BsideTFLArr,TDiff,TempIn,temperature)

      call DiffTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,maxCNiter,epsCN,Re,&
                             TempBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr,TDiff,&
                             temperature2,temperature,2._knd*alpha(RKstage),dt,CNiters,CNres)


      !$hmppcg grid blocksize myblocksize
      !$hmppcg gridify(k,i), private(j)
      do k=-1,Prnz+2
       do i=-1,Prnx+2
        do j=-1,Prny+2
          temperature(i,j,k) = temperature2(i,j,k)
        enddo
       enddo
      enddo


  end subroutine RKstage_Temperature_GPU





  subroutine BoundViscosity_GPU(Prnx,Prny,Prnz,Btype,Nu)
    implicit none
#include "hmpp-include.f90"
    integer, intent(in) :: Prnx,Prny,Prnz,Btype(6)
    real(knd),intent(inout) :: Nu(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
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

  endsubroutine BoundViscosity_GPU


  !$hmpp <tsteps> BoundViscosity_TDiff codelet
  subroutine BoundViscosity_TDiff(Prnx,Prny,Prnz,enable_buoyancy,Btype,Re,Prandtl,Visc,TDiff)
    implicit none
#include "hmpp-include.f90"
    integer,intent(in)      :: Prnx,Prny,Prnz,enable_buoyancy,Btype(6)
    real(knd),intent(in)    :: Re,Prandtl
    real(knd),intent(inout) :: Viscosity(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
    integer                 :: i,j,k
    real(knd),parameter     :: Prt = 0.6

    call BoundViscosity_GPU(Prnx,Prny,Prnz,Btype,Visc)

    if (enable_buoyancy==1) then
      if (Re>0) then

       !$hmppcg grid blocksize myblocksize
       !$hmppcg permute (k,i,j)
       !$hmppcg gridify(k,i), private(j)
       do k=1,Prnz
         do j=1,Prny
           do i=1,Prnx
             TDiff(i,j,k) = (Viscosity(i,j,k)-1._knd/Re)/Prt + (1._knd/(Re*Prandtl))
           end do
         end do
       end do
      else
       !$hmppcg grid blocksize myblocksize
       !$hmppcg permute (k,i,j)
       !$hmppcg gridify(k,i), private(j)
       do k=1,Prnz
         do j=1,Prny
           do i=1,Prnx
             TDiff(i,j,k) = Viscosity(i,j,k)/Prt
           end do
         end do
       end do
      endif

      call BoundViscosity_GPU(Prnx,Prny,Prnz,Btype,TDiff)
    end if

  end subroutine BoundViscosity_TDiff








! !   !$hmpp <tsteps> DiffScalar codelet
!   subroutine DIFFSCALAR_GPU(nscalars,Prnx,Prny,Prnz,dxmin,dymin,dzmin,maxCNiter,epsCN,Re,&
!                             ScalBtype,sideScal,&
!                             TDiff,SCAL2,SCAL,&
!                             coef,dt,l,res)
!   implicit none
! #include "hmpp-include.f90"
!
!   integer,intent(in)    :: nscalars,Prnx,Prny,Prnz,maxCNiter,ScalBtype(6)
!   real(knd),intent(out),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,nscalars) :: Scal2
!   real(knd),intent(in),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,nscalars)  :: Scal
!
!   real(knd),intent(in)  :: dxmin,dymin,dzmin,epsCN,Re,coef,dt
!   real(knd),intent(in)  :: sideScal(6),TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
!   integer,intent(out)   :: l
!   real(knd),intent(out) :: res
!   real(knd) Scal3(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,nscalars)
!   integer nx,ny,nz,i,j,k,m,xi,yj,zk
!   real(knd) p
!   real(knd) A,Ax,Ay,Az,Ap(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
!
!   intrinsic max, abs, mod
!
!
!   Ax = 1._knd/(4._knd*dxmin**2)
!   Ay = 1._knd/(4._knd*dymin**2)
!   Az = 1._knd/(4._knd*dzmin**2)
!
!   !$hmppcg grid blocksize myblocksize
!   !$hmppcg gridify(k,i)
!   do k = 1,Prnz
!    do i = 1,Prnx
!     do j = 1,Prny
!       Ap(i,j,k) = 1._knd/(1._knd/A+(((TDiff(i+1,j,k)+TDiff(i,j,k))+&
!                            (TDiff(i,j,k)+TDiff(i-1,j,k)))*Ax+&
!                            ((TDiff(i,j+1,k)+TDiff(i,j,k))+&
!                            (TDiff(i,j,k)+TDiff(i,j-1,k)))*Ay+&
!                            ((TDiff(i,j,k+1)+TDiff(i,j,k))+&
!                            (TDiff(i,j,k)+TDiff(i,j,k-1)))*Az))
!     enddo
!    enddo
!   enddo
!
!   do m = 1,nscalars
!
!     if (Re>0) then
!
!     !$hmppcg grid blocksize myblocksize
!      !$hmppcg gridify(k,i)
!      do k = 1,Prnz  !initital value using forward Euler
!       do i = 1,Prnx
!        do j = 1,Prny
!         SCAL3(i,j,k,l) = (((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL(i+1,j,k,l)-SCAL(i,j,k,l))/dxmin-&
!           (TDiff(i,j,k)+TDiff(i-1,j,k))*(SCAL(i,j,k,l)-SCAL(i-1,j,k,l))/dxmin)/(dxmin)+&
!          ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL(i,j+1,k,l)-SCAL(i,j,k,l))/dymin-&
!           (TDiff(i,j,k)+TDiff(i,j-1,k))*(SCAL(i,j,k,l)-SCAL(i,j-1,k,l))/dymin)/(dymin)+&
!          ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL(i,j,k+1,l)-SCAL(i,j,k,l))/dzmin-&
!           (TDiff(i,j,k)+TDiff(i,j,k-1))*(SCAL(i,j,k,l)-SCAL(i,j,k-1,l))/dzmin)/(dzmin))
!        enddo
!       enddo
!      enddo
!
!      A = dt*coef
!     !$hmppcg grid blocksize myblocksize
!      !$hmppcg gridify(k,i)
!      do k = 1,Prnz
!       do i = 1,Prnx
!        do j = 1,Prny
!          SCAL2(i,j,k,l) = SCAL(i,j,k,l)+A*SCAL3(i,j,k,l)
!        enddo
!       enddo
!      enddo
!
!
!      call BoundScalar_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,ScalBtype,sideScal,TDiff,SCAL2(:,:,:,l))
!
!
!
!
!
!      do l = 1,maxCNiter
!       res = 0
!
!       call BoundScalar_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,ScalBtype,sideScal,TDiff,SCAL2(:,:,:,l))
!
!
!      !$hmppcg grid blocksize myblocksize
!      !$hmppcg gridify(k,i), reduce(max:res)
!       do k = 1,Prnz
!        do i = 1,Prnx
!         do j = 1+mod(i+k,2),Prny,2
!           p = (SCAL(i,j,k,l)/A)+(SCAL3(i,j,k,l)/4._knd+&
!            ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL2(i+1,j,k,l))-&
!             (TDiff(i,j,k)+TDiff(i-1,j,k))*(-SCAL2(i-1,j,k,l)))*Ax+&
!            ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL2(i,j+1,k,l))-&
!             (TDiff(i,j,k)+TDiff(i,j-1,k))*(-SCAL2(i,j-1,k,l)))*Ay+&
!            ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL2(i,j,k+1,l))-&
!             (TDiff(i,j,k)+TDiff(i,j,k-1))*(-SCAL2(i,j,k-1,l)))*Az&
!            )
!            p = p*Ap(i,j,k)
!            res = max(res,abs(p-SCAL2(i,j,k,l)))
!            SCAL2(i,j,k,l) = p
!         enddo
!        enddo
!       enddo
!
!      !$hmppcg grid blocksize myblocksize
!      !$hmppcg gridify(k,i), reduce(max:res)
!       do k = 1,Prnz
!        do i = 1,Prnx
!         do j = 1+mod(i+k+1,2),Prny,2
!           p = (SCAL(i,j,k,l)/A)+(SCAL3(i,j,k,l)/4._knd+&
!            ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL2(i+1,j,k,l))-&
!             (TDiff(i,j,k)+TDiff(i-1,j,k))*(-SCAL2(i-1,j,k,l)))*Ax+&
!            ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL2(i,j+1,k,l))-&
!             (TDiff(i,j,k)+TDiff(i,j-1,k))*(-SCAL2(i,j-1,k,l)))*Ay+&
!            ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL2(i,j,k+1,l))-&
!             (TDiff(i,j,k)+TDiff(i,j,k-1))*(-SCAL2(i,j,k-1,l)))*Az&
!            )
!            p = p*Ap(i,j,k)
!            res = max(res,abs(p-SCAL2(i,j,k,l)))
!            SCAL2(i,j,k,l) = p
!         enddo
!        enddo
!       enddo
!
!
!       if (res<=epsCN) exit
!      enddo
!
!
!     else
!
!      !$hmppcg grid blocksize myblocksize
!      !$hmppcg gridify(k,i)
!      do k = 1,Prnz
!       do i = 1,Prnx
!        do j = 1,Prny
!          SCAL2(i,j,k,l) = SCAL(i,j,k,l)
!        enddo
!       enddo
!      enddo
!
!
!       call BoundScalar_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,ScalBtype,sideScal,TDiff,SCAL2(:,:,:,l))
!
!     endif
!
!   enddo
!   endsubroutine DIFFSCALAR_GPU


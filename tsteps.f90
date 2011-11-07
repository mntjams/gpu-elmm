module TSTEPS

  use UPWIND
  use CDS, only: CDU, CDV, CDW, KappaU, KappaV, KappaW
  use LAXFRIED
  use LAXWEND
  use PARAMETERS
  use BOUNDARIES, only: Bound_CondU,Bound_CondV,Bound_CondW,Bound_Q
  use POISSON, only: Pr_Correct
  use SMAGORINSKY, only: Smag, Dynsmag, StabSmag, Vreman
  use SCALARS, only:  Bound_Temp, Bound_Visc, Scalar, percdistrib, AdvScalar,&
     DiffScalar, ComputeTDiff, Deposition, GravSettling, ScalFlSourc
  use TURBINLET, only: GetTurbInlet, GetInletFromFile
  use GEOMETRIC, only: IBLinInt,IBBiLinInt,IBTriLinInt
  use Wallmodels, only: ComputeViscsWM

  implicit none


  private
  public TMarchEul,TMarchAB2,TMarchRK2,TMarchRK3,TMarchShiftInlet

  logical:: released=.false.

contains
 subroutine TMarchEul(U,V,W,Pr,delta)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),intent(out):: delta
  real(KND),allocatable,dimension(:,:,:),save:: Q
  real(KND),dimension(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),lbound(U,3):ubound(U,3)):: U2
  real(KND),dimension(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),lbound(V,3):ubound(V,3)):: V2
  real(KND),dimension(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),lbound(W,3):ubound(W,3)):: W2
  real(KND),allocatable,dimension(:,:,:,:),save:: Scalar_2
  real(KND),allocatable,dimension(:,:,:),save:: temperature2
  integer,save:: called=0
  integer i

  
 if (called==0) then
   called=1
   if (computescalars>0) then
    allocate(Scalar_2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
   endif
   if (buoyancy>0) then
    allocate(temperature2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
   endif 
   if (masssourc==1) allocate(Q(0:Prnx+1,0:Prny+1,0:Prnz+1))
  endif 

  if ((BtypeW==TurbulentInlet).or.(BtypeE==TurbulentInlet)) call GetTurbInlet
 
  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)
  if (buoyancy==1)  call Bound_Temp(temperature)

 

  U2=1e19
  V2=1e19
  W2=1e19

  call timestepEUL(U,V,W)

  if (steady==0.and.dt+time>endtime)  dt=endtime-time
  write (*,*) "time:",time,"dt: ",dt



  if (convmet==1) then
   call LF(U2,V2,W2,U,V,W)
  elseif (convmet==2) then
   call MC1(U2,V2,W2,U,V,W)
  elseif (convmet==3) then
   call Tau(U2,V2,W2,U,V,W,Pr)
  else
   U2(1:Unx,1:Uny,1:Unz)=0
   V2(1:Vnx,1:Vny,1:Vnz)=0
   W2(1:Wnx,1:Wny,1:Wnz)=0
  endif
  call CoriolisForce(U2,V2,U,V,1._KND)
  if (buoyancy==1) call BuoyancyForce(W2,temperature,1._KND)


  call OtherTerms(U,V,W,U2,V2,W2,Pr,1._KND)


  if (BtypeT==FreeSlipBuff)  call AttenuateTop(U2,V2,W2,Pr)


  if (masssourc==1) then
      call MASS_SOURC(Q,U2,V2,W2)
  endif
  
 
  
  call Bound_CondU(U2)
  call Bound_CondV(V2)
  call Bound_CondW(W2)

  if (poissmet>0) then
  if (masssourc==1) then
    call Pr_Correct(U2,V2,W2,Pr,1._KND,Q)
  else
    call Pr_Correct(U2,V2,W2,Pr,1._KND)
  endif
  endif


  
  if (computescalars>0) then
   Scalar_2=0
   do i=1,computescalars
    call AdvScalar(Scalar_2(:,:,:,i),Scalar(:,:,:,i),U2,V2,W2,2,1._KND)
   enddo
   Scalar=Scalar+Scalar_2
   if (sgstype/=StabSmagorinskyModel)  call ComputeTDiff(U2,V2,W2)
   call Bound_Visc(TDiff)
   do i=1,computescalars
    call DiffScalar(Scalar_2(:,:,:,i),Scalar(:,:,:,i),2,1._KND)
   enddo
   if (computedeposition>0) call Deposition(Scalar_2,1._KND)
   if (computegravsettling>0) call GravSettling(Scalar_2,1._KND)
   Scalar=Scalar_2
  endif
    
  if (buoyancy>0) then
   if (sgstype/=StabSmagorinskyModel)  call ComputeTDiff(U2,V2,W2)
   call Bound_Visc(TDiff)
   temperature2=0
   call Bound_Temp(temperature)
   call AdvScalar(temperature2,temperature,U2,V2,W2,1,1._KND)
   temperature=temperature+temperature2
   call Bound_Temp(temperature2)
   call DiffScalar(temperature2,temperature,1,1._KND)
!    temperature2=temperature+temperature2
   temperature=temperature2
  endif


  delta=sum(abs(U(1:Unx,1:Uny,1:Unz)-U2(1:Unx,1:Uny,1:Unz)))/(Unx*Uny*Unz)
  delta=delta+sum(abs(V(1:Vnx,1:Vny,1:Vnz)-V2(1:Vnx,1:Vny,1:Vnz)))/(Vnx*Vny*Vnz)
  delta=delta+sum(abs(W(1:Wnx,1:Wny,1:Wnz)-W2(1:Wnx,1:Wny,1:Wnz)))/(Wnx*Wny*Wnz)
  delta=delta/dt


  U=U2
  V=V2
  W=W2
 end subroutine TMarchEul
  



 subroutine TMarchAB2(U,V,W,Pr,delta)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),intent(out):: delta

  real(KND),save:: dtlast
  real(KND),save:: rho,beta

  real(KND),allocatable,dimension(:,:,:),save:: Q,U2,V2,W2,Ustar,Vstar,Wstar

  real(KND),allocatable,dimension(:,:,:),save:: Temperature_adv,Temperature2
  real(KND),allocatable,dimension(:,:,:,:),save:: Scalar_adv,Scalar_2

  integer i,j,k
  real(KND) p
  integer,save:: called=0
  real time1,time2


  if (called==0) then

   called=1
   if (masssourc==1) allocate(Q(0:Prnx+1,0:Prny+1,0:Prnz+1))

   allocate(Ustar(-2:Unx+3,-2:Uny+3,-2:Unz+3))
   allocate(Vstar(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
   allocate(Wstar(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))

   allocate(U2(-2:Unx+3,-2:Uny+3,-2:Unz+3))
   allocate(V2(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
   allocate(W2(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))

   if (buoyancy==1) then
    allocate(temperature_adv(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    allocate(temperature2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
   endif

   if (computescalars>0) then
    allocate(Scalar_adv(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
    allocate(Scalar_2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
   endif

   rho=0
   beta=1._KND

   Ustar=0
   Vstar=0
   Wstar=0
   Scalar_adv=0
   Temperature_adv=0

  else
   rho=-(dt/dtlast)/2._KND
   beta=1._KND+(dt/dtlast)/2._KND

  endif

  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)

  if (buoyancy==1)  call Bound_Temp(temperature)

  if ((BtypeW==TurbulentInlet).or.(BtypeE==TurbulentInlet)) call GetTurbInlet
  if (BtypeW==InletFromFile) call GetInletFromFile(time)

  call timestepEUL(U,V,W)

  if (steady==0.and.dt+time>endtime)  dt=endtime-time

  write (*,*) "time:",time,"dt: ",dt

  temperature_adv=0


    if (debugparam>1) call cpu_time(time1)

    U2=0
    V2=0
    W2=0


    if (convmet>0) then

      U2=U2+Ustar*rho
      V2=V2+Vstar*rho
      W2=W2+Wstar*rho


      Ustar=0
      Vstar=0
      Wstar=0

      if (convmet==1) then
       call LF(Ustar,Vstar,Wstar,U,V,W)
      elseif (convmet==3) then
       call KAPPAU(Ustar,U,V,W,1._KND)
       call KAPPAV(Vstar,U,V,W,1._KND)
       call KAPPAW(Wstar,U,V,W,1._KND)
      else
       call CDU(Ustar,U,V,W,1._KND)
       call CDV(Vstar,U,V,W,1._KND)
       call CDW(Wstar,U,V,W,1._KND)
      endif

      call CoriolisForce(Ustar,Vstar,U,V,1._KND)
      if (buoyancy==1) call BuoyancyForce(Wstar,temperature,1._KND)

      U2=U2+Ustar*beta
      V2=V2+Vstar*beta
      W2=W2+Wstar*beta
    endif

    if (debugparam>1) then
     call cpu_time(time2)
     write (*,*) "ET of part 1", (time2-time1)
     time1=time2
    endif

    call OtherTerms(U,V,W,U2,V2,W2,Pr,1._KND)

    if (debugparam>1) then
     call cpu_time(time2)
     write (*,*) "ET of part 2", (time2-time1)
     time1=time2
    endif

    if (BtypeT==FreeSlipBuff)  call AttenuateTop(U2,V2,W2,Pr)
    if (BtypeE==OutletBuff) then
      if (buoyancy==1) then
        call AttenuateOut(U2,V2,W2,Pr,temperature)
      else
        call AttenuateOut(U2,V2,W2,Pr)
      endif
    endif

    if (masssourc==1) then
        call MASS_SOURC(Q,U2,V2,W2)
    endif


    call Bound_CondU(U2)
    call Bound_CondV(V2)
    call Bound_CondW(W2)


    if (poissmet>0) then
     if (masssourc==1) then
       call Pr_Correct(U2,V2,W2,Pr,1._KND,Q)
     else
       call Pr_Correct(U2,V2,W2,Pr,1._KND)
     endif
    endif


    if (computescalars>0.and..not.released) call Explosion

    if (computescalars>0) then
      Scalar_2=0
 
      Scalar_2=Scalar_2+Scalar_adv*rho

      Scalar_adv=0
      do i=1,computescalars
       call AdvScalar(Scalar_adv(:,:,:,i),Scalar(:,:,:,i),U2,V2,W2,2,1._KND)
      enddo
      Scalar_2=Scalar_2+Scalar_adv*beta

      if (pointscalsource==1) then
       do i=1,computescalars
        Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i)=Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i)+&
         percdistrib(i)*dt*totalscalsource/(dxPr(scalsrci(i))*dyPr(scalsrcj(i))*dzPr(scalsrck(i)))
       enddo
      endif

      Scalar=Scalar+Scalar_2
      if (sgstype/=StabSmagorinskyModel)  call ComputeTDiff(U2,V2,W2)
      call Bound_Visc(TDiff)
      do i=1,computescalars
         call DiffScalar(Scalar_2(:,:,:,i),Scalar(:,:,:,i),2,1._KND)
      enddo
      if (computedeposition>0) call Deposition(Scalar_2,1._KND)
      if (computegravsettling>0) call GravSettling(Scalar_2,1._KND)
      Scalar=Scalar_2
    endif

    if (buoyancy>0) then
      if (sgstype/=StabSmagorinskyModel)  call ComputeTDiff(U2,V2,W2)
      call Bound_Visc(TDiff)
      temperature2=0
      call Bound_Temp(temperature)


      temperature2=temperature2+temperature_adv*rho

      temperature_adv=0
      call AdvScalar(temperature_adv,temperature,U2,V2,W2,1,1._KND)
      temperature2=temperature2+temperature_adv*beta

      temperature=temperature+temperature2
      call Bound_Temp(temperature)
       call DiffScalar(temperature2,temperature,1,1._KND)
      temperature=temperature2
    endif

    if (BtypeE==OutletBuff) then
      if (buoyancy==1) then
        call AttenuateOut(U2,V2,W2,Pr,temperature)
      else
        call AttenuateOut(U2,V2,W2,Pr)
      endif
    endif


    delta=sum(abs(U(1:Unx,1:Uny,1:Unz)-U2(1:Unx,1:Uny,1:Unz)))/(Unx*Uny*Unz)
    delta=delta+sum(abs(V(1:Vnx,1:Vny,1:Vnz)-V2(1:Vnx,1:Vny,1:Vnz)))/(Vnx*Vny*Vnz)
    delta=delta+sum(abs(W(1:Wnx,1:Wny,1:Wnz)-W2(1:Wnx,1:Wny,1:Wnz)))/(Wnx*Wny*Wnz)


    call NullInterior(U2,V2,W2)

    U=U2
    V=V2
    W=W2

    dtlast=dt

 end subroutine TMarchAB2






  
 subroutine TMarchRK2(U,V,W,Pr,delta)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),intent(out):: delta
  real(KND),allocatable,dimension(:,:,:),save:: Q
  real(KND),dimension(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),&
                                              lbound(U,3):ubound(U,3)):: U2,Ustar
  real(KND),dimension(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),&
                                              lbound(V,3):ubound(V,3)):: V2,Vstar
  real(KND),dimension(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),&
                                              lbound(W,3):ubound(W,3)):: W2,Wstar
  real(KND) p
  integer i,j,k
  real(KND),allocatable,dimension(:,:,:,:),save:: Scalar_2
  real(KND),allocatable,dimension(:,:,:),save:: temperature2
  integer,save:: called=0


  if (called==0) then
   called=1
   if (computescalars>0) then
    allocate(Scalar_2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
   endif
   if (buoyancy>0) then
    allocate(temperature2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
   endif 
   if (masssourc==1) allocate(Q(0:Prnx+1,0:Prny+1,0:Prnz+1))
  endif 
  
  if ((BtypeW==TurbulentInlet).or.(BtypeE==TurbulentInlet)) call GetTurbInlet

  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)
  if (buoyancy==1)  call Bound_Temp(temperature)

  


      call timestepEUL(U,V,W)

  if (steady==0.and.dt+time>endtime)  dt=endtime-time
  write (*,*) "time:",time,"dt: ",dt

  U2=1e9
  V2=1e9
  W2=1e9



  if (convmet>0) then
   U2(1:Unx,1:Uny,1:Unz)=0
   V2(1:Vnx,1:Vny,1:Vnz)=0
   W2(1:Wnx,1:Wny,1:Wnz)=0
   call CDU(U2,U,V,W,1._KND)
   call CDV(V2,U,V,W,1._KND)
   call CDW(W2,U,V,W,1._KND)
   call CoriolisForce(Ustar,Vstar,U,V,1._KND)
   if (buoyancy==1) call BuoyancyForce(Wstar,temperature,1._KND)

   Ustar=U+U2
   Vstar=V+V2
   Wstar=W+W2

   U2=0
   V2=0
   W2=0

   call CDU(U2,Ustar,Vstar,Wstar,0.5_KND)
   call CDV(V2,Ustar,Vstar,Wstar,0.5_KND)
   call CDW(W2,Ustar,Vstar,Wstar,0.5_KND)
   call CoriolisForce(U2,V2,Ustar,Vstar,1._KND)
     if (buoyancy==1) call BuoyancyForce(W2,temperature,1._KND)

   
   U2=Ustar/2._KND-U/2._KND+U2
   V2=Vstar/2._KND-V/2._KND+V2
   W2=Wstar/2._KND-W/2._KND+W2
  else
   U2(1:Unx,1:Uny,1:Unz)=0
   V2(1:Vnx,1:Vny,1:Vnz)=0
   W2(1:Wnx,1:Wny,1:Wnz)=0
  endif 
 

   call OtherTerms(U,V,W,U2,V2,W2,Pr,1._KND)


  if (BtypeT==FreeSlipBuff)  call AttenuateTop(U2,V2,W2,Pr)
  if (BtypeE==OutletBuff) call AttenuateOut(U2,V2,W2,Pr)

  
  if (masssourc==1) then
      call MASS_SOURC(Q,U2,V2,W2)
  endif
   

  call Bound_CondU(U2)
  call Bound_CondV(V2)
  call Bound_CondW(W2)
  

  if (poissmet>0) then
  if (masssourc==1) then
    call Pr_Correct(U2,V2,W2,Pr,1._KND,Q)
  else
    call Pr_Correct(U2,V2,W2,Pr,1._KND)
  endif  
  endif

  
  call Bound_CondU(U2)
  call Bound_CondV(V2)
  call Bound_CondW(W2)
  
  

  if (computescalars>=4)then
   if (time>(endtime-starttime)/3._KND.and.maxval(scalar(:,:,:,1))==0) then
    p=0
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (((xPr(i)>-3.5.and.xPr(i)<3.5).or.(xPr(i)>3.5.and.xPr(i-1)<-3.5)).and.&
           ((yPr(j)>-3.5.and.yPr(j)<3.5).or.(yPr(j)>3.5.and.yPr(j-1)<-3.5)).and.&
           (zPr(k)<12.or.(zPr(k)>12.and.zPr(k-1)<0))) then
        Scalar(i,j,k,1)=0.10
        p=p+dxPr(i)*dyPr(j)*dzPr(k)
       endif
      enddo
     enddo
    enddo
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (((xPr(i)>-3.5.and.xPr(i)<3.5).or.(xPr(i)>3.5.and.xPr(i-1)<-3.5)).and.&
           ((yPr(j)>-3.5.and.yPr(j)<3.5).or.(yPr(j)>3.5.and.yPr(j-1)<-3.5)).and.&
           (zPr(k)<12.or.(zPr(k)>12.and.zPr(k-1)<0))) then
        Scalar(i,j,k,2)=0.13
       endif
      enddo
     enddo
    enddo
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (((xPr(i)>-3.5.and.xPr(i)<3.5).or.(xPr(i)>3.5.and.xPr(i-1)<-3.5)).and.&
           ((yPr(j)>-3.5.and.yPr(j)<3.5).or.(yPr(j)>3.5.and.yPr(j-1)<-3.5)).and.&
           (zPr(k)<12.or.(zPr(k)>12.and.zPr(k-1)<0))) then
        Scalar(i,j,k,3)=0.64
       endif
      enddo
     enddo
    enddo
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (((xPr(i)>-3.5.and.xPr(i)<3.5).or.(xPr(i)>3.5.and.xPr(i-1)<-3.5)).and.&
           ((yPr(j)>-3.5.and.yPr(j)<3.5).or.(yPr(j)>3.5.and.yPr(j-1)<-3.5)).and.&
           (zPr(k)<12.or.(zPr(k)>12.and.zPr(k-1)<0))) then
        Scalar(i,j,k,4)=0.13
       endif
      enddo
     enddo
    enddo
    Scalar=Scalar/p
   endif

  elseif (computescalars>0.and.partdistrib>0) then
   if (time>(endtime-starttime)/2._KND.and.maxval(scalar(:,:,:,1))==0) then
    p=0
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (((xPr(i)>-3.5.and.xPr(i)<3.5).or.(xPr(i)>3.5.and.xPr(i-1)<-3.5)).and.&
           ((yPr(j)>-3.5.and.yPr(j)<3.5).or.(yPr(j)>3.5.and.yPr(j-1)<-3.5)).and.&
           (zPr(k)<12.or.(zPr(k)>12.and.zPr(k-1)<0))) then
        Scalar(i,j,k,1)=1._KND
        p=p+dxPr(i)*dyPr(j)*dzPr(k)
       endif
      enddo
     enddo
    enddo
    Scalar=Scalar/p
   endif
  endif
  
  if (computescalars>0) then
   Scalar_2=0
   do i=1,computescalars
    call AdvScalar(Scalar_2(:,:,:,i),Scalar(:,:,:,i),U2,V2,W2,2,1._KND)
   enddo
   Scalar=Scalar+Scalar_2
   if (sgstype/=StabSmagorinskyModel)  call ComputeTDiff(U2,V2,W2)
   call Bound_Visc(TDiff)
    do i=1,computescalars
     call DiffScalar(Scalar_2(:,:,:,i),Scalar(:,:,:,i),2,1._KND)
    enddo
    if (computedeposition>0) call Deposition(Scalar_2,1._KND)
    if (computegravsettling>0) call GravSettling(Scalar_2,1._KND)
   Scalar=Scalar_2
  endif

  if (buoyancy>0) then
   if (sgstype/=StabSmagorinskyModel)  call ComputeTDiff(U2,V2,W2)
   call Bound_Visc(TDiff)
   temperature2=0
   call Bound_Temp(temperature)
   call AdvScalar(temperature2,temperature,U2,V2,W2,1,1._KND)
   temperature=temperature+temperature2
   call Bound_Temp(temperature)
   call DiffScalar(temperature2,temperature,1,1._KND)
   temperature=temperature2
  endif
  
  call NullInterior(U2,V2,W2)
  
  delta=sum(abs(U(1:Unx,1:Uny,1:Unz)-U2(1:Unx,1:Uny,1:Unz)))/(Unx*Uny*Unz)
  delta=delta+sum(abs(V(1:Vnx,1:Vny,1:Vnz)-V2(1:Vnx,1:Vny,1:Vnz)))/(Vnx*Vny*Vnz)
  delta=delta+sum(abs(W(1:Wnx,1:Wny,1:Wnz)-W2(1:Wnx,1:Wny,1:Wnz)))/(Wnx*Wny*Wnz)
  U=U2
  V=V2
  W=W2
    
 endsubroutine TMarchRK2
  
  
  
  












  
 subroutine TMarchRK3(U,V,W,Pr,delta)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),intent(out):: delta

  real(KND),allocatable,dimension(:,:,:),save:: Q
  real(KND),dimension(lbound(U,1):ubound(U,1),lbound(U,2):ubound(U,2),&
                                              lbound(U,3):ubound(U,3))::U2,Ustar
  real(KND),dimension(lbound(V,1):ubound(V,1),lbound(V,2):ubound(V,2),&
                                              lbound(V,3):ubound(V,3))::V2,Vstar
  real(KND),dimension(lbound(W,1):ubound(W,1),lbound(W,2):ubound(W,2),&
                                              lbound(W,3):ubound(W,3))::W2,Wstar
  real(KND),dimension(lbound(Scalar,1):ubound(Scalar,1),lbound(Scalar,2):ubound(Scalar,2),&
                lbound(Scalar,3):ubound(Scalar,3),lbound(Scalar,4):ubound(Scalar,4))::Scalar_adv,Scalar_2
  real(KND),dimension(lbound(Temperature,1):ubound(Temperature,1),lbound(Temperature,2):ubound(Temperature,2),&
   lbound(Temperature,3):ubound(Temperature,3))::Temperature_adv,Temperature2
  real(KND),dimension(1:3),save:: alpha,beta,rho
  integer i,j,k,l
  real(KND) p
  integer,save:: called=0
  real time1,time2

i=0

write(*,*) "Test codelet:", i
  if (called==0) then
   alpha(1)=4._KND/15._KND
   alpha(2)=1._KND/15._KND
   alpha(3)=1._KND/6._KND
   beta(1)=8._KND/15._KND
   beta(2)=5._KND/12._KND
   beta(3)=3._KND/4._KND
   rho(1)=0
   rho(2)=-17._KND/60._KND
   rho(3)=-5._KND/12._KND

   called=1
   if (masssourc==1) allocate(Q(0:Prnx+1,0:Prny+1,0:Prnz+1))
  endif

  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)

  if (buoyancy==1)  call Bound_Temp(temperature)

  if ((BtypeW==TurbulentInlet).or.(BtypeE==TurbulentInlet)) call GetTurbInlet
  if (BtypeW==InletFromFile) call GetInletFromFile(time)

  call timestepEUL(U,V,W)

  if (steady==0.and.dt+time>endtime)  dt=endtime-time
  write (*,*) "time:",time,"dt: ",dt

  temperature_adv=0
  Ustar=0
  Vstar=0
  Wstar=0

  do l=1,3
    if (debugparam>1) call cpu_time(time1)

    U2=0
    V2=0
    W2=0

    if (convmet>0) then
      if (l>1) then
        U2=U2+Ustar*rho(l)
        V2=V2+Vstar*rho(l)
        W2=W2+Wstar*rho(l)
      endif

      Ustar=0
      Vstar=0
      Wstar=0

      if (convmet==1) then
       call LF(Ustar,Vstar,Wstar,U,V,W)
      elseif (convmet==3) then
       call KAPPAU(Ustar,U,V,W,1._KND)
       call KAPPAV(Vstar,U,V,W,1._KND)
       call KAPPAW(Wstar,U,V,W,1._KND)
      else if (GPU>0) then
       write(*,*) "Call GPU CDS"
       call BOUND_CONDU(U)
       call BOUND_CONDV(V)
       call BOUND_CONDW(W)
       !$hmpp CDS_GPU callsite
       call CDS_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxmin,dymin,dzmin,dt,Ustar,Vstar,Wstar,U,V,W)
       write(*,*) "Back from  GPU CDS"
      else
       call CDU(Ustar,U,V,W,1._KND)
       call CDV(Vstar,U,V,W,1._KND)
       call CDW(Wstar,U,V,W,1._KND)
      endif

      call CoriolisForce(Ustar,Vstar,U,V,1._KND)
      if (buoyancy==1) call BuoyancyForce(Wstar,temperature,1._KND)

      U2=U2+Ustar*beta(l)
      V2=V2+Vstar*beta(l)
      W2=W2+Wstar*beta(l)
    endif

    if (debugparam>1) then
     call cpu_time(time2)
     write (*,*) "ET of part 1", (time2-time1)
     time1=time2
    endif

    call OtherTerms(U,V,W,U2,V2,W2,Pr,2._KND*alpha(l))

    if (debugparam>1) then
     call cpu_time(time2)
     write (*,*) "ET of part 2", (time2-time1)
     time1=time2
    endif

    if (BtypeT==FreeSlipBuff)  call AttenuateTop(U2,V2,W2,Pr)
    if (BtypeE==OutletBuff) then
      if (buoyancy==1) then
        call AttenuateOut(U2,V2,W2,Pr,temperature)
      else
        call AttenuateOut(U2,V2,W2,Pr)
      endif
    endif

    if (masssourc==1) then
        call MASS_SOURC(Q,U2,V2,W2)
    endif


    call Bound_CondU(U2)
    call Bound_CondV(V2)
    call Bound_CondW(W2)


    if (poissmet>0) then
     if (masssourc==1) then
       call Pr_Correct(U2,V2,W2,Pr,2._KND*alpha(l),Q)
     else
       call Pr_Correct(U2,V2,W2,Pr,2._KND*alpha(l))
     endif
    endif


    if (computescalars>0.and..not.released) call Explosion

    if (computescalars>0) then
      Scalar_2=0
      if (l>1) then
       Scalar_2=Scalar_2+Scalar_adv*rho(l)
      endif
      Scalar_adv=0
      do i=1,computescalars
       call AdvScalar(Scalar_adv(:,:,:,i),Scalar(:,:,:,i),U2,V2,W2,2,1._KND)
      enddo
      Scalar_2=Scalar_2+Scalar_adv*beta(l)

      if (pointscalsource==1) then
       do i=1,computescalars
        Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i)=Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i)+&
         percdistrib(i)*(rho(l)+beta(l))*dt*totalscalsource/(dxPr(scalsrci(i))*dyPr(scalsrcj(i))*dzPr(scalsrck(i)))
       enddo
      endif

      Scalar=Scalar+Scalar_2
      if (sgstype/=StabSmagorinskyModel)  call ComputeTDiff(U2,V2,W2)
      call Bound_Visc(TDiff)
      do i=1,computescalars
         call DiffScalar(Scalar_2(:,:,:,i),Scalar(:,:,:,i),2,2._KND*alpha(l))
      enddo
      if (computedeposition>0) call Deposition(Scalar_2,2._KND*alpha(l))
      if (computegravsettling>0) call GravSettling(Scalar_2,2._KND*alpha(l))
      Scalar=Scalar_2
    endif

    if (buoyancy>0) then
      if (sgstype/=StabSmagorinskyModel)  call ComputeTDiff(U2,V2,W2)
      call Bound_Visc(TDiff)
      temperature2=0
      call Bound_Temp(temperature)

      if (l>1) then
       temperature2=temperature2+temperature_adv*rho(l)
      endif
      temperature_adv=0
      call AdvScalar(temperature_adv,temperature,U2,V2,W2,1,1._KND)
      temperature2=temperature2+temperature_adv*beta(l)

      temperature=temperature+temperature2
      call Bound_Temp(temperature)
       call DiffScalar(temperature2,temperature,1,2._KND*alpha(l))
      temperature=temperature2
    endif

    if (BtypeE==OutletBuff) then
      if (buoyancy==1) then
        call AttenuateOut(U2,V2,W2,Pr,temperature)
      else
        call AttenuateOut(U2,V2,W2,Pr)
      endif
    endif

    if (l==1) delta=0

    delta=delta+sum(abs(U(1:Unx,1:Uny,1:Unz)-U2(1:Unx,1:Uny,1:Unz)))/(Unx*Uny*Unz)
    delta=delta+sum(abs(V(1:Vnx,1:Vny,1:Vnz)-V2(1:Vnx,1:Vny,1:Vnz)))/(Vnx*Vny*Vnz)
    delta=delta+sum(abs(W(1:Wnx,1:Wny,1:Wnz)-W2(1:Wnx,1:Wny,1:Wnz)))/(Wnx*Wny*Wnz)


    call NullInterior(U2,V2,W2)

    U=U2
    V=V2
    W=W2
   enddo
 end subroutine TMarchRK3
    
  



  
   
 subroutine TMarchShiftInlet(U,V,W,Pr,delta) !Only shifts inlet in the x direction, for debugging purposes
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),intent(out):: delta
  integer i,j,k
  real(KND),allocatable,dimension(:,:,:,:),save:: Scalar_2
  real(KND),allocatable,dimension(:,:,:),save:: temperature2
  integer,save:: called=0


  if (called==0) then
    called=1
    if (computescalars>0) then
     allocate(Scalar_2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
    endif
    if (buoyancy>0) then
     allocate(temperature2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    endif 
  endif 
 
  if ((BtypeW==TurbulentInlet).or.(BtypeE==TurbulentInlet)) call GetTurbInlet
  if (BtypeW==InletFromFile) call GetInletFromFile(time)

  write (*,*) "time:",time,"dt: ",dt
 
  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)

  do k=1,Unz
   do j=1,Uny
    do i=Unx,1,-1
     U(i,j,k)=U(i-1,j,k)
    enddo
   enddo
  enddo
  do k=1,Vnz
   do j=1,Vny
    do i=Vnx,1,-1
     V(i,j,k)=V(i-1,j,k)
    enddo
   enddo
  enddo
  do k=1,Wnz
   do j=1,Wny
    do i=Wnx,1,-1
     W(i,j,k)=W(i-1,j,k)
    enddo
   enddo
  enddo

  if (buoyancy==1) then
    call Bound_Temp(temperature)
    do k=1,Prnz
     do j=1,Prny
      do i=Prnx,1,-1
       temperature(i,j,k)=temperature(i-1,j,k)
      enddo
     enddo
    enddo
  endif

  if (poissmet>0) then
    call Pr_Correct(U,V,W,Pr,1._KND)
  endif

  if (computescalars>0)then
   if (time>5.and.maxval(scalar(:,:,:,1))==0) then
    call GridCoords(i,j,k,0._KND,0._KND,0.5_KND)
    Scalar(i,j,k,1)=1._KND/(dxPr(i)*dyPr(j)*dzPr(k))
   endif
  endif
  


  if (wallmodeltype>0) then
                   call ComputeViscsWM(U,V,W,Pr)
  endif

  delta=1
 endsubroutine TMarchShiftInlet
  
  
  
 subroutine Explosion
  real(KND) xc,yc,xs,xf,ys,yf,zs,zf,dxp,dyp,dzp,ct,cr,xp,yp,zp,p
  integer i,j,k,xi,yj,zk,nprobx,nproby,nprobz
  ct=7
  cr=1.5
  xc=3*cos((xheading-70)*pi/180.)
  yc=3*sin((xheading-70)*pi/180.)
  xs=xc-cr
  ys=yc-cr
  zs=0
  xf=xc+cr
  yf=yc+cr
  zf=ct
  nprobx=100
  nproby=100
  nprobz=100
  dxp=(xf-xs)/nprobx
  dyp=(yf-ys)/nproby
  dzp=(zf-zs)/nprobz
  if (computescalars>=4) then
   if (time>(endtime-starttime)/3._KND) then
    Scalar=0
    p=0
    do k=0,nprobz
     zp=zs+k*dzp
     do j=0,nproby
      yp=ys+j*dyp
      do i=0,nprobx
       xp=xs+i*dxp
        call GridCoords(xi,yj,zk,xp,yp,zp)
        if ((xp-xc)**2+(yp-yc)**2<cr**2) then
        if   (zp<ct*0.2) then
         Scalar(xi,yj,zk,:)=Scalar(xi,yj,zk,:)+percdistrib(:)*0.2
         p=p+1
        elseif (zp<ct*0.4) then
         Scalar(xi,yj,zk,:)=Scalar(xi,yj,zk,:)+percdistrib(:)*0.8
         p=p+1
        elseif (zp<ct*0.6) then
         Scalar(xi,yj,zk,:)=Scalar(xi,yj,zk,:)+percdistrib(:)*1.25
         p=p+1
        elseif (zp<ct*0.8) then
         Scalar(xi,yj,zk,:)=Scalar(xi,yj,zk,:)+percdistrib(:)*1.75
         p=p+1
        elseif (zp<=ct) then
         Scalar(xi,yj,zk,:)=Scalar(xi,yj,zk,:)+percdistrib(:)*1.1
         p=p+1
        endif
       endif
      enddo
     enddo
    enddo
    Scalar=totalscalsource*Scalar/p
    Scalar=Scalar/(dxmin*dymin*dzmin)
    released=.true.
   endif
  endif

  
 endsubroutine Explosion
  
  
  
  


  subroutine MomSourc(U,V,W)
   real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
   real(KND) intvel,p0,p1,p2,p3,p4,p5,vel1,vel2,vel3,vel4,vel5,vel6,vel7
   type(TIBPoint),pointer:: IBP
   integer x,y,z,dirx,diry,dirz,n

   call Bound_CondU(U)
   call Bound_CondV(V)
   call Bound_CondW(W)

  if (associated(FirstIBPoint)) then
   IBP => FirstIBPoint
   do
    x=IBP%x
    y=IBP%y
    z=IBP%z
    dirx=IBP%dirx
    diry=IBP%diry
    dirz=IBP%dirz
    if (IBP%interp<=0) then
     intvel=0
    elseif (IBP%interp==1) then
      n=0
      if (dirx/=0)  n=n+1
      if (diry/=0)  n=n+1
      if (dirz/=0)  n=n+1
      if (n/=1) then
                 write (*,*) "Assert error, more directions.", dirx,diry,dirz
                 stop
       endif

     if (IBP%component==1) then
      vel1=U(x+dirx,y+diry,z+dirz)
      vel2=U(x+2*dirx,y+2*diry,z+2*dirz)
      if (IBP%interpdir==1) then
       p0=abs(IBP%distx)
       p1=abs(xU(x+dirx)-xU(x))
       p2=abs(xU(x+2*dirx)-xU(x))
      elseif (IBP%interpdir==2) then
       p0=abs(IBP%disty)
       p1=abs(yPr(y+diry)-yPr(y))
       p2=abs(yPr(y+2*diry)-yPr(y))
      else
       p0=abs(IBP%distz)
       p1=abs(zPr(z+dirz)-zPr(z))
       p2=abs(zPr(z+2*dirz)-zPr(z))
      endif
    elseif (IBP%component==2) then
      vel1=V(x+dirx,y+diry,z+dirz)
      vel2=V(x+2*dirx,y+2*diry,z+2*dirz)
      if (IBP%interpdir==1) then
       p0=abs(IBP%distx)
       p1=abs(xPr(x+dirx)-xPr(x))
       p2=abs(xPr(x+2*dirx)-xPr(x))
      elseif (IBP%interpdir==2) then
       p0=abs(IBP%disty)
       p1=abs(yV(y+diry)-yV(y))
       p2=abs(yV(y+2*diry)-yV(y))
      else
       p0=abs(IBP%distz)
       p1=abs(zPr(z+dirz)-zPr(z))
       p2=abs(zPr(z+2*dirz)-zPr(z))
     endif
     else
      vel1=W(x+dirx,y+diry,z+dirz)
      vel2=W(x+2*dirx,y+2*diry,z+2*dirz)
      if (IBP%interpdir==1) then
       p0=abs(IBP%distx)
       p1=abs(xPr(x+dirx)-xPr(x))
       p2=abs(xPr(x+2*dirx)-xPr(x))
      elseif (IBP%interpdir==2) then
       p0=abs(IBP%disty)
       p1=abs(yPr(y+diry)-yPr(y))
       p2=abs(yPr(y+2*diry)-yPr(y))
      else
       p0=abs(IBP%distz)
       p1=abs(zW(z+dirz)-zW(z))
       p2=abs(zW(z+2*dirz)-zW(z))
      endif
     endif

     intvel=IBLinInt(p0,p1,p2,vel1,vel2)

    elseif (IBP%interp==2) then
     if (IBP%component==1) then
      if (IBP%interpdir==1) then
       vel1=U(x,y+diry,z)
       vel2=U(x,y,z+dirz)
       vel3=U(x,y+diry,z+dirz)
       p0=abs(IBP%disty)
       p1=abs(IBP%distz)
       p2=abs(yPr(y+diry)-yPr(y))
       p3=abs(zPr(z+dirz)-zPr(z))
      elseif (IBP%interpdir==2) then
       vel1=U(x,y,z+dirz)
       vel2=U(x+dirx,y,z)
       vel3=U(x+dirx,y,z+dirz)
       p0=abs(IBP%distz)
       p1=abs(IBP%distx)
       p2=abs(zPr(z+dirz)-zPr(z))
       p3=abs(xU(x+dirx)-xU(x))
      else
       vel1=U(x+dirx,y,z)
       vel2=U(x,y+diry,z)
       vel3=U(x+dirx,y+diry,z)
       p0=abs(IBP%distx)
       p1=abs(IBP%disty)
       p2=abs(xU(x+dirx)-xU(x))
       p3=abs(yPr(y+diry)-yPr(y))
      endif
     elseif (IBP%component==2) then
      if (IBP%interpdir==1) then
       vel1=V(x,y+diry,z)
       vel2=V(x,y,z+dirz)
       vel3=V(x,y+diry,z+dirz)
       p0=abs(IBP%disty)
       p1=abs(IBP%distz)
       p2=abs(yV(y+diry)-yV(y))
       p3=abs(zPr(z+dirz)-zPr(z))
      elseif (IBP%interpdir==2) then
       vel1=V(x,y,z+dirz)
       vel2=V(x+dirx,y,z)
       vel3=V(x+dirx,y,z+dirz)
       p0=abs(IBP%distz)
       p1=abs(IBP%distx)
       p2=abs(zPr(z+dirz)-zPr(z))
       p3=abs(xPr(x+dirx)-xPr(x))
      else
       vel1=V(x+dirx,y,z)
       vel2=V(x,y+diry,z)
       vel3=V(x+dirx,y+diry,z)
       p0=abs(IBP%distx)
       p1=abs(IBP%disty)
       p2=abs(xPr(x+dirx)-xPr(x))
       p3=abs(yV(y+diry)-yV(y))
      endif
     else
      if (IBP%interpdir==1) then
       vel1=W(x,y+diry,z)
       vel2=W(x,y,z+dirz)
       vel3=W(x,y+diry,z+dirz)
       p0=abs(IBP%disty)
       p1=abs(IBP%distz)
       p2=abs(yPr(y+diry)-yPr(y))
       p3=abs(zW(z+dirz)-zW(z))
      elseif (IBP%interpdir==2) then
       vel1=W(x,y,z+dirz)
       vel2=W(x+dirx,y,z)
       vel3=W(x+dirx,y,z+dirz)
       p0=abs(IBP%distz)
       p1=abs(IBP%distx)
       p2=abs(zW(z+dirz)-zW(z))
       p3=abs(xPr(x+dirx)-xPr(x))
      else
       vel1=W(x+dirx,y,z)
       vel2=W(x,y+diry,z)
       vel3=W(x+dirx,y+diry,z)
       p0=abs(IBP%distx)
       p1=abs(IBP%disty)
       p2=abs(xPr(x+dirx)-xPr(x))
       p3=abs(yPr(y+diry)-yPr(y))
      endif
     endif

     intvel=IBBiLinInt(p0,p1,p2,p3,vel1,vel2,vel3)


    elseif (IBP%interp==3) then
     if (IBP%component==1) then
       vel1=U(x+dirx,y,z)
       vel2=U(x,y+diry,z)
       vel3=U(x+dirx,y+diry,z)
       vel4=U(x,y,z+dirz)
       vel5=U(x+dirx,y,z+dirz)
       vel6=U(x,y+diry,z+dirz)
       vel7=U(x+dirx,y+diry,z+dirz)
       p0=abs(IBP%distx)
       p1=abs(IBP%disty)
       p2=abs(IBP%distz)
       p3=abs(xU(x+dirx)-xU(x))
       p4=abs(yPr(y+diry)-yPr(y))
       p5=abs(zPr(z+dirz)-zPr(z))
     elseif (IBP%component==2) then
       vel1=V(x+dirx,y,z)
       vel2=V(x,y+diry,z)
       vel3=V(x+dirx,y+diry,z)
       vel4=V(x,y,z+dirz)
       vel5=V(x+dirx,y,z+dirz)
       vel6=V(x,y+diry,z+dirz)
       vel7=V(x+dirx,y+diry,z+dirz)
       p0=abs(IBP%distx)
       p1=abs(IBP%disty)
       p2=abs(IBP%distz)
       p3=abs(xPr(x+dirx)-xPr(x))
       p4=abs(yV(y+diry)-yV(y))
       p5=abs(zPr(z+dirz)-zPr(z))
     else
       vel1=W(x+dirx,y,z)
       vel2=W(x,y+diry,z)
       vel3=W(x+dirx,y+diry,z)
       vel4=W(x,y,z+dirz)
       vel5=W(x+dirx,y,z+dirz)
       vel6=W(x,y+diry,z+dirz)
       vel7=W(x+dirx,y+diry,z+dirz)
       p0=abs(IBP%distx)
       p1=abs(IBP%disty)
       p2=abs(IBP%distz)
       p3=abs(xPr(x+dirx)-xPr(x))
       p4=abs(yPr(y+diry)-yPr(y))
       p5=abs(zW(z+dirz)-zW(z))
     endif

     intvel=IBTriLinInt(p0,p1,p2,p3,p4,p5,vel1,vel2,vel3,vel4,vel5,vel6,vel7)

    else
     intvel=0
    endif

    if (IBP%component==1) then
     IBP%MSourc=(intvel-U(x,y,z))/dt
    elseif (IBP%component==2) then
     IBP%MSourc=(intvel-V(x,y,z))/dt
    elseif (IBP%component==3) then
     IBP%MSourc=(intvel-W(x,y,z))/dt
    endif


    if (associated(IBP%next)) then
     IBP=>IBP%next
    else
     exit
    endif
   enddo
  endif

  end subroutine MomSourc
   

  subroutine MASS_SOURC(Q,U,V,W)
  real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
  real(KND),intent(out):: Q(0:,0:,0:)
  type(TIBPoint),pointer:: IBP
  integer x,y,z

   Q=0
   if (associated(FirstIBPoint)) then
    IBP => FirstIBPoint
    do
     x=IBP%x
     y=IBP%y
     z=IBP%z
      if (IBP%component==1) then
       Q(x,y,z)=Q(x,y,z)+U(x,y,z)/(dxPr(x))
       Q(x+1,y,z)=Q(x+1,y,z)-U(x,y,z)/(dxPr(x+1))
      elseif (IBP%component==2) then
       Q(x,y,z)=Q(x,y,z)+V(x,y,z)/(dyPr(y))
       Q(x,y+1,z)=Q(x,y+1,z)-V(x,y,z)/(dyPr(y+1))
      elseif (IBP%component==3) then
       Q(x,y,z)=Q(x,y,z)+W(x,y,z)/(dzPr(z))
       Q(x,y,z+1)=Q(x,y,z+1)-W(x,y,z)/(dzPr(z+1))
      endif
     if (associated(IBP%next)) then
      IBP=>IBP%next
     else
      exit
     endif
    enddo
   endif
   call Bound_Q(Q)
  end subroutine MASS_SOURC







  subroutine TIMESTEPCW(U,V,W)
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer i,j,k
  real(KND) m,p
  m=0

  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
     p=MAX(abs(U(i,j,k)/dxU(i)),abs(U(i-1,j,k)/dxU(i-1)))+&
        MAX(abs(V(i,j,k)/dyV(j)),abs(V(i,j-1,k)/dyV(j-1)))+&
        MAX(abs(W(i,j,k)/dzW(k)),abs(W(i,j,k-1)/dzW(k-1)))
     if (p>m) m=p
    enddo
   enddo
  enddo
  
  if (m>0.1*Uinlet/dxmin) then 
   dt=0.25_KND/(m)
  else
   dt=0.25_KND/(Uinlet/dxmin)
  endif
  endsubroutine TIMESTEPCW
  

  
  subroutine TIMESTEPEUL(U,V,W)
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer i,j,k,maxI,maxJ,maxK
  real(KND) m,p
  m=0
 maxI=0
 maxJ=0
 maxK=0

  do k=1,Prnz
   do j=1,Prny
    do i=1,Prnx
     p=MAX(MAX(abs(U(i,j,k)/dxU(i)),abs(U(i-1,j,k)/dxU(i-1))),MAX(abs(V(i,j,k)/dyV(j)),abs(V(i,j-1,k)/dyV(j-1))))&
         +MAX(abs(W(i,j,k)/dzW(k)),abs(W(i,j,k-1)/dzW(k-1)))
     if (p>m) then
               m=p
               maxI=i
               maxj=j
               maxk=k
     endif
    enddo
   enddo
  enddo


  if (m>0) then 
   dt=MIN(CFL/(m),dxmin/Uref)
  else
   dt=dxmin/Uref
  endif
  endsubroutine TIMESTEPEUL







  subroutine BuoyancyForce(W2,theta,coef)
  real(KND),dimension(-2:,-2:,-2:),intent(inout):: W2
  real(KND),dimension(-1:,-1:,-1:),intent(in):: theta
  real(KND),intent(in):: coef
  real(KND) A
  integer i,j,k

    A=grav_acc*coef*dt/temperature_ref
    do k=1,Wnz
     do j=1,Wny
      do i=1,Wnx
            W2(i,j,k)=W2(i,j,k)+A*((theta(i,j,k+1)+theta(i,j,k))/2._KND-temperature_ref)
      enddo
     enddo
    enddo
  endsubroutine BuoyancyForce





  subroutine CoriolisForce(U2,V2,U,V,coef)
  real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V
  real(KND),dimension(-2:,-2:,-2:),intent(inout):: U2,V2
  real(KND),intent(in):: coef
  real(KND) A
  integer i,j,k

   A=-coef*dt
   if (coriolisparam>0) then
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx   
           U2(i,j,k)=U2(i,j,k)-A*coriolisparam*(V(i,j-1,k)+V(i+1,j-1,k)+V(i,j,k)+V(i+1,j,k))/4._KND
      enddo
     enddo
    enddo   
    do k=1,Vnz
     do j=1,Vny
      do i=1,Vnx
           V2(i,j,k)=V2(i,j,k)+A*coriolisparam*(U(i-1,j,k)+U(i-1,j+1,k)+U(i,j,k)+U(i,j+1,k))/4._KND
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
   real(KND) Af,Ap,Apre,Aprn,Aprt,S

   integer i,j,k,x,y,z

   type(TIBPoint),pointer:: IBP
write(*,*) "Otherterms:"   
   Apre=-coef*dt
   Aprn=-coef*dt
   Aprt=-coef*dt
   Af=dt

   call Bound_Pr(Pr)
   !Pressure gradient terms
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx   
          U2(i,j,k)=U2(i,j,k)+Apre*(Pr(i+1,j,k)-Pr(i,j,k))/dxU(i)+Apre*prgradientx
     enddo
    enddo
   enddo
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
          V2(i,j,k)=V2(i,j,k)+Aprn*(Pr(i,j+1,k)-Pr(i,j,k))/dyV(j)+Aprn*prgradienty
     enddo
    enddo
   enddo
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
          W2(i,j,k)=W2(i,j,k)+Aprt*(Pr(i,j,k+1)-Pr(i,j,k))/dzW(k)
     enddo
    enddo
   enddo


   call Bound_CondU2(U2)
   call Bound_CondV2(V2)
   call Bound_CondW2(W2)

   Re_gt_0: if (Re>0) then

     !Diffusion using Crank Nicolson
     !first approximation using forward Euler
     !iteration SOR or Gauss-Seidel


     call Bound_CondU(U)
     call Bound_CondV(V)
     call Bound_CONDW(W)
     U3=U+U2
     V3=V+V2
     W3=W+W2
     if (sgstype==SmagorinskyModel) then
                       call Smag(U,V,W)
     elseif (sgstype==DynSmagorinskyModel) then
                       call DynSmag(U,V,W)
     elseif (sgstype==VremanModel) then
                       call Vreman(U,V,W)
     elseif (sgstype==StabSmagorinskyModel) then
                       call StabSmag(U,V,W)
     else
                       Visc=1._KND/Re
     endif

     if (debuglevel>0) then
      write(*,*) "NUt", sum(Visc(1:Prnx,1:Prny,1:Prnz))/(Prnx*Prny*Prnz)
      write(*,*) "maxNUt", MAXVAL(Visc(1:Prnx,1:Prny,1:Prnz))
      write(*,*) "minNUt", MINVAL(Visc(1:Prnx,1:Prny,1:Prnz))
     endif

     if (wallmodeltype>0) then
                     call ComputeViscsWM(U,V,W,Pr)
     endif

     call ScalFlSourc(Visc,3)

     Ap=coef*dt

     do k=1,Unz    !Forward Euler for the first approximation
      do j=1,Uny
       do i=1,Unx
        U3(i,j,k)=U(i,j,k)+U2(i,j,k)+Ap*(&
        ((Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))/dxPr(i+1)-&
        Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k))/dxPr(i))/dxU(i)+&
         0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))/dyV(j)-&
         (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k))/dyV(j-1))/dyPr(j)+&
         ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)-&
         (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1))/dzW(k-1))/dzPr(k))))
       enddo
      enddo
     enddo
     do k=1,Vnz
      do j=1,Vny
       do i=1,Vnx
        V3(i,j,k)=V(i,j,k)+V2(i,j,k)+Ap*(&
        ((Visc(i,j+1,k)*(V(i,j+1,k)-V(i,j,k))/dyPr(j+1)-&
         Visc(i,j,k)*(V(i,j,k)-V(i,j-1,k))/dyPr(j))/dyV(j)+&
         0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V(i+1,j,k)-V(i,j,k))/dxU(i)-&
        (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(V(i,j,k)-V(i-1,j,k))/dxU(i-1))/dxPr(i)+&
         ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))/dzW(k)-&
         (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(V(i,j,k)-V(i,j,k-1))/dzW(k-1))/dzPr(k))))
       enddo
      enddo
     enddo
     do k=1,Wnz
      do j=1,Wny
       do i=1,Wnx
        W3(i,j,k)=W(i,j,k)+W2(i,j,k)+Ap*(&
        ((0.25_KND*((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W(i+1,j,k)-W(i,j,k))/dxU(i)-&
        (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(W(i,j,k)-W(i-1,j,k))/dxU(i-1))/dxPr(i)+&
         ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W(i,j+1,k)-W(i,j,k))/dyV(j)-&
         (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(W(i,j,k)-W(i,j-1,k))/dyV(j-1))/dyPr(j))+&
         (Visc(i,j,k+1)*(W(i,j,k+1)-W(i,j,k))/dzPr(k+1)-&
         Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1))/dzPr(k))/dzW(k)))
       enddo
      enddo
     enddo


     call Bound_CondU(U3)
     call Bound_CondV(V3)
     call Bound_CondW(W3)
     call MomSourc(U3,V3,W3)

     if (associated(FirstIBPoint)) then   !Immersed boundary terms, in the future should be in an array
      IBP => FirstIBPoint
      do
       x=IBP%x
       y=IBP%y
       z=IBP%z
        if (IBP%component==1) then
         U3(x,y,z)=U3(x,y,z)+IBP%MSourc*dt
         U2(x,y,z)=U2(x,y,z)+IBP%MSourc*dt
        elseif (IBP%component==2) then
         V3(x,y,z)=V3(x,y,z)+IBP%MSourc*dt
         V2(x,y,z)=V2(x,y,z)+IBP%MSourc*dt
        elseif (IBP%component==3) then
         W3(x,y,z)=W3(x,y,z)+IBP%MSourc*dt
         W2(x,y,z)=W2(x,y,z)+IBP%MSourc*dt
        endif
       if (associated(IBP%next)) then
        IBP=>IBP%next
       else
        exit
       endif
      enddo
     endif

     if (gridtype==UNIFORMGRID.and.GPU>0) then                  !Performs the diffusion terms
      write (*,*) "GPU CN call"
      !$hmpp UNIFREDBLACK_GPU callsite
      call UNIFREDBLACK_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,dt,dxmin,dymin,dzmin,&
                            U,V,W,U2,V2,W2,U3,V3,W3,Visc,coef,maxCNiter,epsCN)
      write(*,*) "back from GPU CN"
     else if (gridtype==UNIFORMGRID) then
      call UNIFREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
     else
      call GENREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
     endif

     U2=U3
     V2=V3
     W2=W3
  
   else  Re_gt_0  !Re<=0

    U2=U+U2
    V2=V+V2
    W2=W+W2

   endif   Re_gt_0

   if (debuglevel>=2) then  !Compute and output the mean friction in the domain.
    S=0
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
       S=S-((Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))/dxPr(i+1)-&
       Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k))/dxPr(i))/dxU(i)+&
         (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))/dyV(j)-&
         0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k))/dyV(j-1))/dyPr(j)+&
          (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)-&
         0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1))/dzW(k-1))/dzPr(k))
      enddo
     enddo
    enddo

    S=S/(Unx*Uny*Unz)
    write(*,*) "Mean friction:", S
   endif
  endsubroutine OtherTerms




  subroutine UNIFREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
   real(KND),dimension(-2:,-2:,-2:),intent(inout):: U,V,W,U2,V2,W2,U3,V3,W3
   real(KND),intent(in):: coef
   real(KND),dimension(1:Unx,1:Uny,1:Unz):: Apu
   real(KND),dimension(1:Vnx,1:Vny,1:Vnz):: ApV
   real(KND),dimension(1:Wnx,1:Wny,1:Wnz):: ApW
   real(KND) recdxmin2,recdymin2,recdzmin2                                                               !reciprocal values of dx**2
   real(KND) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
   integer i,j,k,l

       Ap=coef*dt/(2._KND)
       S=0
       l=0

       recdxmin2=1./dxmin**2
       recdymin2=1./dymin**2
       recdzmin2=1./dzmin**2

       do k=1,Unz    !The explicit part, which doesn't have to be changed inside the loop
        do j=1,Uny
         do i=1,Unx
          U2(i,j,k)=U2(i,j,k)+Ap*(&
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
          V2(i,j,k)=V2(i,j,k)+Ap*(&
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
         W2(i,j,k)=W2(i,j,k)+Ap*(&
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
          ApU(i,j,k)=((Visc(i+1,j,k)+&
                      Visc(i,j,k))*recdxmin2+&
                      0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))+&
                      (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k)))*recdymin2+&
                      ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                      (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1)))*recdzmin2))
         enddo
        enddo
       enddo

       ApU=1._KND/(1._KND+Ap*ApU)

       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          ApV(i,j,k)=((Visc(i,j+1,k)+&
                     Visc(i,j,k))*recdymin2+&
                     0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))+&
                      (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k)))*recdxmin2+&
                     ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                     (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1)))*recdzmin2))
         enddo
        enddo
       enddo

       ApV=1._KND/(1._KND+Ap*ApV)

       do k=1,Wnz
        do j=1,Wny
         do i=1,Wnx
          ApW(i,j,k)=(0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                      (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k)))*recdxmin2+&
                     ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))+&
                     (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k)))*recdymin2)+&
                     (Visc(i,j,k+1)+&
                     Visc(i,j,k))*recdzmin2)
         enddo
        enddo
       enddo

       ApW=1._KND/(1._KND+Ap*ApW)


       Suavg=abs(MAXVAL(U3(1:Unx,1:Uny,1:Unz)))  !maximum values of velocities to norm the residues.
       Svavg=abs(MAXVAL(V3(1:Vnx,1:Vny,1:Vnz)))
       Swavg=abs(MAXVAL(W3(1:Wnx,1:Wny,1:Wnz)))
       if (Suavg<=1e-3_KND) Suavg=1
       if (Svavg<=1e-3_KND) Svavg=1
       if (Swavg<=1e-3_KND) Swavg=1



       do l=1,maxCNiter               !Gauss-Seidel iteration for Crank-Nicolson result
        call Bound_CondU(U3)
        call Bound_CondV(V3)
        call Bound_CondW(W3)
        S=0
        Su=0
        Sv=0
        Sw=0
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
            p=Ap*p+U2(i,j,k)+U(i,j,k)
            p=p*ApU(i,j,k)


            Su=max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
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
            p=Ap*p+V2(i,j,k)+V(i,j,k)
            p=p*ApV(i,j,k)
            Sv=max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
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
            p=Ap*p+W2(i,j,k)+W(i,j,k)
            p=p*ApW(i,j,k)
            Sw=max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
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
            p=Ap*p+U2(i,j,k)+U(i,j,k)
            p=p*ApU(i,j,k)


            Su=max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
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
            p=Ap*p+V2(i,j,k)+V(i,j,k)
            p=p*ApV(i,j,k)
            Sv=max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
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
            p=Ap*p+W2(i,j,k)+W(i,j,k)
            p=p*ApW(i,j,k)
            Sw=max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        !$OMP ENDDO
        !$OMP END PARALLEL
        S=max(Su/Suavg,Sv/Svavg,Sw/Swavg)
        write (*,*) "CN ",l,S
        if (S<=epsCN) exit
       enddo
  endsubroutine UNIFREDBLACK



  subroutine GENREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
   real(KND),dimension(-2:,-2:,-2:),intent(inout):: U,V,W,U2,V2,W2,U3,V3,W3
   real(KND),intent(in):: coef
   real(KND),dimension(1:Unx,1:Uny,1:Unz):: Apu
   real(KND),dimension(1:Vnx,1:Vny,1:Vnz):: ApV
   real(KND),dimension(1:Wnx,1:Wny,1:Wnz):: ApW
   real(KND) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
   integer i,j,k,l
   integer ind(3)


       Ap=coef*dt/(2._KND)
       S=0
       l=0


       do k=1,Unz    !The explicit part, which doesn't have to be changed inside the loop
        do j=1,Uny
         do i=1,Unx
          U2(i,j,k)=U2(i,j,k)+Ap*(&
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
          V2(i,j,k)=V2(i,j,k)+Ap*(&
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
         W2(i,j,k)=W2(i,j,k)+Ap*(&
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
          ApU(i,j,k)=((Visc(i+1,j,k)/dxPr(i+1)+&
                      Visc(i,j,k)/dxPr(i))/dxU(i)+&
                      (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))/dyV(j)+&
                      0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))/dyV(j-1))/dyPr(j)+&
                      (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))/dzW(k)+&
                      0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))/dzW(k-1))/dzPr(k))
         enddo
        enddo
       enddo

      ApU=1._KND/(1._KND+Ap*ApU)


       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          ApV(i,j,k)= ((0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))/dxU(i)+&
                      0.25_KND*(Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))/dxU(i-1))/dxPr(i)+&
                     (Visc(i,j+1,k)/dyPr(j+1)+&
                     Visc(i,j,k)/dyPr(j))/dyV(j)+&
                     (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))/dzW(k)+&
                     0.25_KND*(Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))/dzW(k-1))/dzPr(k))
         enddo
        enddo
       enddo

       ApV=1._KND/(1._KND+Ap*ApV)


       do k=1,Wnz
        do j=1,Wny
         do i=1,Wnx
          ApW(i,j,k)= ((0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))/dxU(i)+&
                      0.25_KND*(Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))/dxU(i-1))/dxPr(i)+&
                     (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))/dyV(j)+&
                     0.25_KND*(Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))/dyV(j-1))/dyPr(j)+&
                     (Visc(i,j,k+1)/dzPr(k+1)+&
                     Visc(i,j,k)/dzPr(k))/dzW(k))
         enddo
        enddo
       enddo

       ApW=1._KND/(1._KND+Ap*ApW)

       Suavg=abs(MAXVAL(U3(1:Unx,1:Uny,1:Unz)))  !maximum values of velocities to norm the residues.
       Svavg=abs(MAXVAL(V3(1:Vnx,1:Vny,1:Vnz)))
       Swavg=abs(MAXVAL(W3(1:Wnx,1:Wny,1:Wnz)))
       if (Suavg<=1e-3_KND) Suavg=1
       if (Svavg<=1e-3_KND) Svavg=1
       if (Swavg<=1e-3_KND) Swavg=1



       do l=1,maxCNiter               !Gauss-Seidel iteration for Crank-Nicolson result
        call Bound_CondU(U3)
        call Bound_CondV(V3)
        call Bound_CondW(W3)
        S=0
        Su=0
        Sv=0
        Sw=0
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
            p=Ap*p+U2(i,j,k)+U(i,j,k)
            p=p*ApU(i,j,k)


            Su=max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
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
            p=Ap*p+V2(i,j,k)+V(i,j,k)
            p=p*ApV(i,j,k)
            Sv=max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
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
            p=Ap*p+W2(i,j,k)+W(i,j,k)
            p=p*ApW(i,j,k)
            Sw=max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
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
            p=Ap*p+U2(i,j,k)+U(i,j,k)
            p=p*ApU(i,j,k)


            Su=max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
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
            p=Ap*p+V2(i,j,k)+V(i,j,k)
            p=p*ApV(i,j,k)
            Sv=max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
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
            p=Ap*p+W2(i,j,k)+W(i,j,k)
            p=p*ApW(i,j,k)
            Sw=max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        !$OMP ENDDO
        !$OMP END PARALLEL
        S=max(Su/Suavg,Sv/Svavg,Sw/Swavg)
        write (*,*) "CN ",l,S
        if (S<=epsCN) exit
       enddo
   endsubroutine GENREDBLACK








  subroutine AttenuateTop(U,V,W,Pr)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  integer i,j,k,bufn
  real(KND) p,ze,zs,zb,DF

  bufn=max(5,Prnz/4)
  zs=zW(Prnz-bufn)
  ze=zW(Prnz)

    do k=Vnz-bufn,Vnz
     p=0
      do i=1,Vnx
        do j=1,Vny
         p=p+V(i,j,k)
        enddo
      enddo
      p=p/(Vnx*Vny)
      zb=(zPr(k)-zs)/(ze-zs)
      DF=DampF(zb)
      do i=-1,Vnx+1
       do j=-1,Vny+1
        V(i,j,k)=p+DF*(V(i,j,k)-p)
       enddo
      enddo
    enddo
    do k=Unz-bufn,Unz
     p=0
      do i=1,Unx
        do j=1,Uny
         p=p+U(i,j,k)
        enddo
      enddo
      p=p/(Unx*Uny)
      zb=(zPr(k)-zs)/(ze-zs)
      DF=DampF(zb)
      do i=-1,Unx+1
       do j=-1,Uny+1
        U(i,j,k)=p+DF*(U(i,j,k)-p)
       enddo
      enddo
    enddo
    if (buoyancy>0) then
     do k=Prnz-bufn,Prnz
       p=0
       do i=1,Prnx
        do j=1,Prny
         p=p+Temperature(i,j,k)
        enddo
       enddo
      p=p/(Prnx*Prny)
      zb=(zPr(k)-zs)/(ze-zs)
      DF=DampF(zb)
      do i=-1,Prnx+1
       do j=-1,Prny+1
        Temperature(i,j,k)=p+DF*(temperature(i,j,k)-p)
       enddo
      enddo
     enddo
    endif
    do k=Wnz-bufn,Wnz
     p=0
     do i=1,Wnx
      do j=1,Wny
        p=p+W(i,j,k)
      enddo
     enddo
     p=p/(Wnx*Wny)
     zb=(zW(k)-zs)/(ze-zs)
     DF=DampF(zb)
     do i=-1,Wnx+1
      do j=-1,Wny+1
        W(i,j,k)=p+DF*(W(i,j,k)-p)
      enddo
     enddo
    enddo
  endsubroutine AttenuateTop


  subroutine AttenuateOut(U,V,W,Pr,temperature)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),optional,intent(inout):: temperature(-1:,-1:,-1:)
  integer i,j,k,bufn
  real(KND) p,xe,xs,xb,DF

  bufn=max(10,Prnx/8)
  xs=xU(Prnx-bufn)
  xe=xU(Prnx)
    do k=1,Unz
     p=0
     do i=2*Unx/3,Unx-4
      do j=1,Uny
       p=p+U(i,j,k)
      enddo
     enddo
     p=p/((Unx-4-2*Unx/3+1)*Uny)
     do i=Unx-bufn,Unx+1
      xb=(xU(i)-xs)/(xe-xs)
      DF=DampF(xb)
      do j=-1,Uny+1
        U(i,j,k)=p+DF*(U(i,j,k)-p)
      enddo
     enddo
    enddo
    do k=1,Vnz
     p=0
     do i=2*Vnx/3,Vnx-4
      do j=1,Vny
       p=p+V(i,j,k)
      enddo
     enddo
     p=p/((Vnx-4-2*Vnx/3+1)*Vny)
     do i=Vnx-bufn,Vnx+1
      xb=(xPr(i)-xs)/(xe-xs)
      DF=DampF(xb)
      do j=-1,Vny+1
        V(i,j,k)=p+DF*(V(i,j,k)-p)
      enddo
     enddo
    enddo
    do k=1,Wnz
     p=0
     do i=2*Wnx/3,Wnx-4
      do j=1,Wny
       p=p+W(i,j,k)
      enddo
     enddo
     p=p/((Wnx-4-2*Wnx/3+1)*Wny)
     do i=Wnx-bufn,Wnx+1
      xb=(xPr(i)-xs)/(xe-xs)
      DF=DampF(xb)
      do j=-1,Wny+1
        W(i,j,k)=p+DF*(W(i,j,k)-p)
      enddo
     enddo
    enddo
   if (present(temperature).and.buoyancy==1) then
    do k=1,Prnz
     p=0
     do i=2*Prnx/3,Prnx-4
      do j=1,Prny
       p=p+temperature(i,j,k)
      enddo
     enddo
     p=p/((Prnx-4-2*Prnx/3+1)*Prny)
     do i=Prnx-bufn,Prnx+1
      xb=(xU(i)-xs)/(xe-xs)
      DF=DampF(xb)
      do j=-1,Prny+1
        temperature(i,j,k)=p+DF*(temperature(i,j,k)-p)
      enddo
     enddo
    enddo
   endif
  endsubroutine AttenuateOut



  pure function DampF(x)
  real(KND) DampF
  real(KND),intent(in)::x
  if (x<=0) then
    DampF=1
  elseif (x>=1) then
    DampF=0
  else
   DampF=(1-0.1_KND*x**2)*(1-(1-exp(10._KND*x**2))/(1-exp(10._KND)))
  endif
  endfunction Dampf


  subroutine NullInterior(U,V,W)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer i,j,k

  
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
      if (Utype(i,j,k)>0.and.Utype(i,j,k+1)>0.and.Utype(i,j,k-1)>0&
          .and.Utype(i,j-1,k)>0.and.Utype(i,j+1,k)>0&
          .and.Utype(i-1,j,k)>0.and.Utype(i+1,j,k)>0)  U(i,j,k)=0
     enddo
    enddo
   enddo

   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
      if (Vtype(i,j,k)>0.and.Vtype(i,j,k+1)>0.and.Vtype(i,j,k-1)>0&
          .and.Vtype(i,j-1,k)>0.and.Vtype(i,j+1,k)>0&
          .and.Vtype(i-1,j,k)>0.and.Vtype(i+1,j,k)>0)  V(i,j,k)=0
     enddo
    enddo
   enddo

   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
      if (Wtype(i,j,k)>0.and.Wtype(i,j,k+1)>0.and.Wtype(i,j,k-1)>0&
          .and.Wtype(i,j-1,k)>0.and.Wtype(i,j+1,k)>0&
          .and.Wtype(i-1,j,k)>0.and.Wtype(i+1,j,k)>0)  W(i,j,k)=0
     enddo
    enddo
   enddo
  endsubroutine NullInterior






























 !GPU codelets

  !$hmpp UNIFREDBLACK_GPU codelet, target=CUDA
  subroutine UNIFREDBLACK_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,dt,dxmin,dymin,dzmin,&
          U,V,W,U2,V2,W2,U3,V3,W3,Visc,coef,maxCNiter,epsCN)
  implicit none

   integer, parameter:: KND=4,TIM=4

   integer,intent(in):: Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz
   real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(inout):: U,U2,U3
   real(KND),dimension(-2:Vny+3,-2:Vny+3,-2:Vnz+3),intent(inout):: V,V2,V3
   real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(inout):: W,W2,W3
   real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(in):: Visc
   real(TIM),intent(in):: dt
   real(KND),intent(in):: dxmin,dymin,dzmin,coef,epsCN
   integer(KND),intent(in):: maxCNiter
   real(KND),dimension(1:Unx,1:Uny,1:Unz):: Apu
   real(KND),dimension(1:Vnx,1:Vny,1:Vnz):: ApV
   real(KND),dimension(1:Wnx,1:Wny,1:Wnz):: ApW
   real(KND) recdxmin2,recdymin2,recdzmin2                                                               !reciprocal values of dx**2
   real(KND) Ap,p,S,Suavg,Svavg,Swavg,Su,Sv,Sw
   integer i,j,k,l
   intrinsic mod,abs,max

       Ap=coef*dt/(2._KND)
       S=0
       l=0

       recdxmin2=1./dxmin**2
       recdymin2=1./dymin**2
       recdzmin2=1./dzmin**2

       !$hmppcg permute (k,i,j)
       do k=1,Unz    !The explicit part, which doesn't have to be changed inside the loop
        do j=1,Uny
         do i=1,Unx
          U2(i,j,k)=U2(i,j,k)+Ap*(&
          ((Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))-&
          Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k)))*recdxmin2+0.25_KND*(&
           ((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))-&
           (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k)))*recdymin2+&
           ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))-&
           (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1)))*recdzmin2)))
         enddo
        enddo
       enddo
       !$hmppcg permute (k,i,j)
       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          V2(i,j,k)=V2(i,j,k)+Ap*(&
          (0.25_KND*((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V(i+1,j,k)-V(i,j,k))-&
          (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(V(i,j,k)-V(i-1,j,k)))*recdxmin2+&
           (Visc(i,j+1,k)*(V(i,j+1,k)-V(i,j,k))-&
           Visc(i,j,k)*(V(i,j,k)-V(i,j-1,k)))*recdymin2+&
           0.25_KND*((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))-&
           (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(V(i,j,k)-V(i,j,k-1)))*recdzmin2))
         enddo
        enddo
       enddo
       !$hmppcg permute (k,i,j)
       do k=1,Wnz
       do j=1,Wny
        do i=1,Wnx
         W2(i,j,k)=W2(i,j,k)+Ap*(&
         (0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W(i+1,j,k)-W(i,j,k))-&
         (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(W(i,j,k)-W(i-1,j,k)))*recdxmin2+&
          ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W(i,j+1,k)-W(i,j,k))-&
          (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(W(i,j,k)-W(i,j-1,k)))*recdymin2)+&
          (Visc(i,j,k+1)*(W(i,j,k+1)-W(i,j,k))-&
          Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1)))*recdzmin2))
        enddo
       enddo
      enddo

       !$hmppcg permute (k,i,j)
       do k=1,Unz         !Auxiliary coefficients to better efficiency in loops
        do j=1,Uny
         do i=1,Unx
          ApU(i,j,k)=((Visc(i+1,j,k)+&
                      Visc(i,j,k))*recdxmin2+&
                      0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))+&
                      (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k)))*recdymin2+&
                      ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                      (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1)))*recdzmin2))
          ApU(i,j,k)=1._KND/(1._KND+Ap*ApU(i,j,k))
         enddo
        enddo
       enddo


       !$hmppcg permute (k,i,j)
       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          ApV(i,j,k)=((Visc(i,j+1,k)+&
                     Visc(i,j,k))*recdymin2+&
                     0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))+&
                      (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k)))*recdxmin2+&
                     ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                     (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1)))*recdzmin2))
          ApV(i,j,k)=1._KND/(1._KND+Ap*ApV(i,j,k))
         enddo
        enddo
       enddo


       !$hmppcg permute (k,i,j)
       do k=1,Wnz
        do j=1,Wny
         do i=1,Wnx
          ApW(i,j,k)=(0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))+&
                      (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k)))*recdxmin2+&
                     ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))+&
                     (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k)))*recdymin2)+&
                     (Visc(i,j,k+1)+&
                     Visc(i,j,k))*recdzmin2)
          ApW(i,j,k)=1._KND/(1._KND+Ap*ApW(i,j,k))
         enddo
        enddo
       enddo


       Suavg=0    !maximum values of velocities to norm the residues.
       !$hmppcg gridify (k,i), reduce (max:Suavg)
       do k=1,Unz
        do i=1,Unx
         do j=1,Uny
          Suavg=max(Suavg,abs(U3(i,j,k)))
         enddo
        enddo
       enddo
       Svavg=0
       !$hmppcg gridify (k,i), reduce (max:Svavg)
       do k=1,Vnz
        do i=1,Vnx
         do j=1,Vny
          Svavg=max(Svavg,abs(V3(i,j,k)))
         enddo
        enddo
       enddo
       Swavg=0
       !$hmppcg gridify (k,i), reduce (max:Swavg)
       do k=1,Wnz
        do i=1,Wnx
         do j=1,Wny
          Swavg=max(Swavg,abs(W3(i,j,k)))
         enddo
        enddo
       enddo
       if (Suavg<=1e-3_KND) Suavg=1
       if (Svavg<=1e-3_KND) Svavg=1
       if (Swavg<=1e-3_KND) Swavg=1


       l=1
       S=epsCN+1.
       do while (S>epsCN.and.l<=maxCNiter)          !Gauss-Seidel iteration for Crank-Nicolson result
        !call Bound_CondU(U3)
        !call Bound_CondV(V3)
        !call Bound_CondW(W3)
        S=0
        Su=0
        Sv=0
        Sw=0
        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Su)
        do k=1,Unz
         do i=1,Unx
          do j=1+mod(i+k,2),Uny,2
            p=((Visc(i+1,j,k)*(U3(i+1,j,k))-&
             Visc(i,j,k)*(-U3(i-1,j,k)))*recdxmin2+&
             0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U3(i,j+1,k))-&
             (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(-U3(i,j-1,k)))*recdymin2+&
             ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U3(i,j,k+1))-&
             (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-U3(i,j,k-1)))*recdzmin2))
            p=Ap*p+U2(i,j,k)+U(i,j,k)
            p=p*ApU(i,j,k)


            Su=max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Sv)
        do k=1,Vnz
         do i=1,Vnx
          do j=1+mod(i+k,2),Vny,2
            p=((Visc(i,j+1,k)*(V3(i,j+1,k))-&
             Visc(i,j,k)*(-V3(i,j-1,k)))*recdymin2+&
             0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V3(i+1,j,k))-&
             (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(-V3(i-1,j,k)))*recdxmin2+&
             ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V3(i,j,k+1))-&
             (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-V3(i,j,k-1)))*recdzmin2))
            p=Ap*p+V2(i,j,k)+V(i,j,k)
            p=p*ApV(i,j,k)
            Sv=max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Sw)
        do k=1,Wnz
         do i=1,Wnx
          do j=1+mod(i+k,2),Wny,2
            p=(0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W3(i+1,j,k))-&
             (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(-W3(i-1,j,k)))*recdxmin2+&
             ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W3(i,j+1,k))-&
             (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(-W3(i,j-1,k)))*recdymin2)+&
             (Visc(i,j,k+1)*(W3(i,j,k+1))-&
             Visc(i,j,k)*(-W3(i,j,k-1)))*recdzmin2)
            p=Ap*p+W2(i,j,k)+W(i,j,k)
            p=p*ApW(i,j,k)
            Sw=max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo


        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Su)
        do k=1,Unz
         do i=1,Unx
          do j=1+mod(i+k+1,2),Uny,2
            p=((Visc(i+1,j,k)*(U3(i+1,j,k))-&
             Visc(i,j,k)*(-U3(i-1,j,k)))*recdxmin2+&
             0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U3(i,j+1,k))-&
             (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(-U3(i,j-1,k)))*recdymin2+&
             ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U3(i,j,k+1))-&
             (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-U3(i,j,k-1)))*recdzmin2))
            p=Ap*p+U2(i,j,k)+U(i,j,k)
            p=p*ApU(i,j,k)


            Su=max(Su,abs(p-U3(i,j,k)))
            U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Sv)
        do k=1,Vnz
         do i=1,Vnx
          do j=1+mod(i+k+1,2),Vny,2
            p=((Visc(i,j+1,k)*(V3(i,j+1,k))-&
             Visc(i,j,k)*(-V3(i,j-1,k)))*recdymin2+&
             0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V3(i+1,j,k))-&
             (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(-V3(i-1,j,k)))*recdxmin2+&
             ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V3(i,j,k+1))-&
             (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-V3(i,j,k-1)))*recdzmin2))
            p=Ap*p+V2(i,j,k)+V(i,j,k)
            p=p*ApV(i,j,k)
            Sv=max(Sv,abs(p-V3(i,j,k)))
            V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        !$hmppcg grid blocksize 512x1
        !$hmppcg gridify(k,i), reduce(max:Sw)
        do k=1,Wnz
         do i=1,Wnx
          do j=1+mod(i+k+1,2),Wny,2
            p=(0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W3(i+1,j,k))-&
             (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(-W3(i-1,j,k)))*recdxmin2+&
             ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W3(i,j+1,k))-&
             (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(-W3(i,j-1,k)))*recdymin2)+&
             (Visc(i,j,k+1)*(W3(i,j,k+1))-&
             Visc(i,j,k)*(-W3(i,j,k-1)))*recdzmin2)
            p=Ap*p+W2(i,j,k)+W(i,j,k)
            p=p*ApW(i,j,k)
            Sw=max(Sw,abs(p-W3(i,j,k)))
            W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
          enddo
         enddo
        enddo
        S=max(Su/Suavg,Sv/Svavg,Sw/Swavg)
        l=l+1
       enddo
  endsubroutine UNIFREDBLACK_GPU


    !$hmpp CDS_GPU codelet, target=CUDA
    subroutine CDS_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dx,dy,dz,dt,U2,V2,W2,U,V,W)
    implicit none
    integer, parameter:: KND=4

    integer,intent(in)    :: Unx, Uny, Unz, Vnx, Vny, Vnz, Wnx, Wny, Wnz
    real(KND),intent(in)  :: dx, dy, dz, dt
    real(KND),intent(out) :: U2(-2:Unx+3,-2:Uny+3,-2:Unz+3)
    real(KND),intent(in)  :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
    real(KND),intent(out) :: V2(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
    real(KND),intent(in)  :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
    real(KND),intent(out) :: W2(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
    real(KND),intent(in)  :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
    integer i,j,k
    real(KND) Ax,Ay,Az


        Ax=0.25_KND*dt/dx
        Ay=0.25_KND*dt/dy
        Az=0.25_KND*dt/dz
       !$hmppcg grid blocksize 512x1        
       !$hmppcg permute (k,i,j)
        do k=1,Unz
            do j=1,Uny
                do i=1,Unx
                    U2(i,j,k)= - ((Ax*(U(i+1,j,k)+U(i,j,k))*(U(i+1,j,k)+U(i,j,k))&
                    -Ax*(U(i,j,k)+U(i-1,j,k))*(U(i,j,k)+U(i-1,j,k)))&
                    +(Ay*(U(i,j+1,k)+U(i,j,k))*(V(i+1,j,k)+V(i,j,k))&
                    -Ay*(U(i,j,k)+U(i,j-1,k))*(V(i+1,j-1,k)+V(i,j-1,k)))&
                    +(Az*(U(i,j,k+1)+U(i,j,k))*(W(i+1,j,k)+W(i,j,k))&
                    -Az*(U(i,j,k)+U(i,j,k-1))*(W(i+1,j,k-1)+W(i,j,k-1))))
                enddo
            enddo
        enddo

       !$hmppcg grid blocksize 512x1        
       !$hmppcg permute (k,i,j)
        do k=1,Vnz
            do j=1,Vny
                do i=1,Vnx
                    V2(i,j,k)= - ((Ay*(V(i,j+1,k)+V(i,j,k))*(V(i,j+1,k)+V(i,j,k))&
                    -Ay*(V(i,j,k)+V(i,j-1,k))*(V(i,j,k)+V(i,j-1,k)))&
                    +(Ax*(V(i+1,j,k)+V(i,j,k))*(U(i,j+1,k)+U(i,j,k))&
                    -Ax*(V(i,j,k)+V(i-1,j,k))*(U(i-1,j+1,k)+U(i-1,j,k)))&
                    +(Az*(V(i,j,k+1)+V(i,j,k))*(W(i,j+1,k)+W(i,j,k))&
                    -Az*(V(i,j,k)+V(i,j,k-1))*(W(i,j+1,k-1)+W(i,j,k-1))))
                enddo
            enddo
        enddo

       !$hmppcg grid blocksize 512x1        
       !$hmppcg permute (k,i,j)
        do k=1,Wnz
            do j=1,Wny
                do i=1,Wnx
                    W2(i,j,k)= - ((Az*(W(i,j,k+1)+W(i,j,k))*(W(i,j,k+1)+W(i,j,k))&
                    -Az*(W(i,j,k)+W(i,j,k-1))*(W(i,j,k)+W(i,j,k-1)))&
                    +(Ay*(W(i,j+1,k)+W(i,j,k))*(V(i,j,k+1)+V(i,j,k))&
                    -Ay*(W(i,j,k)+W(i,j-1,k))*(V(i,j-1,k)+V(i,j-1,k+1)))&
                    +(Ax*(W(i+1,j,k)+W(i,j,k))*(U(i,j,k+1)+U(i,j,k))&
                    -Ax*(W(i,j,k)+W(i-1,j,k))*(U(i-1,j,k+1)+U(i-1,j,k))))
                enddo
            enddo
        enddo

    end subroutine CDS_GPU



end module TSTEPS

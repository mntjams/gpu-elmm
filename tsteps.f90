module TSTEPS

use UPWIND
use CDS
use LAXFRIED
use LAXWEND
use PARAMETERS
use BOUNDARIES
use POISSON
use SMAGORINSKY
use SCALARS

implicit none

real(KND):: odey=1._KND !for debugging
logical:: released=.false.

contains
  subroutine TMARCHEUL(U,V,W,Pr,delta)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),intent(out):: delta
  real(KND),allocatable,dimension(:,:,:),save:: Q
  real(KND),dimension(LBOUND(U,1):UBOUND(U,1),LBOUND(U,2):UBOUND(U,2),LBOUND(U,3):UBOUND(U,3)):: U2
  real(KND),dimension(LBOUND(V,1):UBOUND(V,1),LBOUND(V,2):UBOUND(V,2),LBOUND(V,3):UBOUND(V,3)):: V2
  real(KND),dimension(LBOUND(W,1):UBOUND(W,1),LBOUND(W,2):UBOUND(W,2),LBOUND(W,3):UBOUND(W,3)):: W2
  integer i,j,k
  real(KND),allocatable,dimension(:,:,:,:),save:: SCALAR_2
  real(KND),allocatable,dimension(:,:,:),save:: temperature2
  integer,save:: called=0

  
 if (called==0) then
   called=1
   if (computescalars>0) then
    allocate(SCALAR_2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
   endif
   if (buoyancy>0) then
    allocate(temperature2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
   endif 
   if (masssourc==1) allocate(Q(0:Prnx+1,0:Prny+1,0:Prnz+1))
  endif 

  if ((BtypeW==TURBULENTINLET).or.(BtypeE==TURBULENTINLET)) call GETTURBINLET
 
  call BOUND_CONDU(U)
  call BOUND_CONDV(V)
  call BOUND_CONDW(W)
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


  call OTHERTERMS(U,V,W,U2,V2,W2,Pr,1._KND)


  if (BtypeT==FREESLIPBUFF)  call ATTENUATETOP(U2,V2,W2,Pr)


  if (masssourc==1) then
      call MASS_SOURC(Q,U2,V2,W2)
  endif
  
 
  
  call BOUND_CONDU(U2)
  call BOUND_CONDV(V2)
  call BOUND_CONDW(W2)

  if (poissmet>0) then
  if (masssourc==1) then
    call PR_CORRECT(U2,V2,W2,Pr,1._KND,Q)
  else
    call PR_CORRECT(U2,V2,W2,Pr,1._KND)
  endif
  endif


  
  if (computescalars>0) then
   SCALAR_2=0
   do i=1,computescalars
    call ADVSCALAR(SCALAR_2(:,:,:,i),SCALAR(:,:,:,i),U2,V2,W2,2,1._KND)
   enddo
   SCALAR=SCALAR+SCALAR_2
   call ComputeTDiff(U2,V2,W2)
   call Bound_Visc(TDiff)
   do i=1,computescalars
    call DIFFSCALAR(SCALAR_2(:,:,:,i),SCALAR(:,:,:,i),2,1._KND)
   enddo
   if (computedeposition>0) call Deposition(SCALAR_2,1._KND)
   if (computegravsettling>0) call Gravsettling(SCALAR_2,1._KND)
   SCALAR=SCALAR_2
  endif
    
  if (buoyancy>0) then
   call ComputeTDiff(U2,V2,W2)
   call Bound_Visc(TDiff)
   temperature2=0
   call Bound_Temp(temperature)
   call CDSSCALAR(temperature2,temperature,U2,V2,W2,1,1._KND)
   temperature=temperature+temperature2
   call Bound_Temp(temperature2)
   call DIFFSCALAR(temperature2,temperature,1,1._KND)
!    temperature2=temperature+temperature2
   temperature=temperature2
  endif


  delta=SUM(ABS(U(1:Unx,1:Uny,1:Unz)-U2(1:Unx,1:Uny,1:Unz)))/(Unx*Uny*Unz)
  delta=delta+SUM(ABS(V(1:Vnx,1:Vny,1:Vnz)-V2(1:Vnx,1:Vny,1:Vnz)))/(Vnx*Vny*Vnz)
  delta=delta+SUM(ABS(W(1:Wnx,1:Wny,1:Wnz)-W2(1:Wnx,1:Wny,1:Wnz)))/(Wnx*Wny*Wnz)
  delta=delta/dt


  U=U2
  V=V2
  W=W2
  end subroutine TMARCHEUL
  

  
  subroutine TMARCHRK2(U,V,W,Pr,delta)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),intent(out):: delta
  real(KND),allocatable,dimension(:,:,:),save:: Q
  real(KND),dimension(LBOUND(U,1):UBOUND(U,1),LBOUND(U,2):UBOUND(U,2),&
                                              LBOUND(U,3):UBOUND(U,3)):: U2,Ustar
  real(KND),dimension(LBOUND(V,1):UBOUND(V,1),LBOUND(V,2):UBOUND(V,2),&
                                              LBOUND(V,3):UBOUND(V,3)):: V2,Vstar
  real(KND),dimension(LBOUND(W,1):UBOUND(W,1),LBOUND(W,2):UBOUND(W,2),&
                                              LBOUND(W,3):UBOUND(W,3)):: W2,Wstar
  real(KND) p,odey2,odeystar
  integer i,j,k
  real(KND),allocatable,dimension(:,:,:,:),save:: SCALAR_2
  real(KND),allocatable,dimension(:,:,:),save:: temperature2
  integer,save:: called=0


 if (called==0) then
   called=1
   if (computescalars>0) then
    allocate(SCALAR_2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
   endif
   if (buoyancy>0) then
    allocate(temperature2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
   endif 
   if (masssourc==1) allocate(Q(0:Prnx+1,0:Prny+1,0:Prnz+1))
  endif 
  
  if ((BtypeW==TURBULENTINLET).or.(BtypeE==TURBULENTINLET)) call GETTURBINLET

  call BOUND_CONDU(U)
  call BOUND_CONDV(V)
  call BOUND_CONDW(W)
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
   call CDU2(U2,U,V,W,1._KND)
   call CDV2(V2,U,V,W,1._KND)
   call CDW2(W2,U,V,W,1._KND)
   call CoriolisForce(Ustar,Vstar,U,V,1._KND)
   if (buoyancy==1) call BuoyancyForce(Wstar,temperature,1._KND)

   Ustar=U+U2
   Vstar=V+V2
   Wstar=W+W2

   U2=0
   V2=0
   W2=0

   call CDU2(U2,Ustar,Vstar,Wstar,0.5_KND)
   call CDV2(V2,Ustar,Vstar,Wstar,0.5_KND)
   call CDW2(W2,Ustar,Vstar,Wstar,0.5_KND)
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
 

   call OTHERTERMS(U,V,W,U2,V2,W2,Pr,1._KND)


  if (BtypeT==FREESLIPBUFF)  call ATTENUATETOP(U2,V2,W2,Pr)
  if (BtypeE==NEUMANN) call ATTENUATEOUT(U2,V2,W2,Pr)

  
  if (masssourc==1) then
      call MASS_SOURC(Q,U2,V2,W2)
  endif
   

  call BOUND_CONDU(U2)
  call BOUND_CONDV(V2)
  call BOUND_CONDW(W2)
  

  if (poissmet>0) then
  if (masssourc==1) then
    call PR_CORRECT(U2,V2,W2,Pr,1._KND,Q)
  else
    call PR_CORRECT(U2,V2,W2,Pr,1._KND)
  endif  
  endif

  
  call BOUND_CONDU(U2)
  call BOUND_CONDV(V2)
  call BOUND_CONDW(W2)
  
  

  if (computescalars>=4)then
   if (time>(endtime-starttime)/3._KND.and.maxval(scalar(:,:,:,1))==0) then
    p=0
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (((xPr(i)>-3.5.and.xPr(i)<3.5).or.(xPr(i)>3.5.and.xPr(i-1)<-3.5)).and.&
           ((yPr(j)>-3.5.and.yPr(j)<3.5).or.(yPr(j)>3.5.and.yPr(j-1)<-3.5)).and.&
           (zPr(k)<12.or.(zPr(k)>12.and.zPr(k-1)<0))) then
        SCALAR(i,j,k,1)=0.10
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
        SCALAR(i,j,k,2)=0.13
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
        SCALAR(i,j,k,3)=0.64
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
        SCALAR(i,j,k,4)=0.13
       endif
      enddo
     enddo
    enddo
    SCALAR=SCALAR/p
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
        SCALAR(i,j,k,1)=1._KND
        p=p+dxPr(i)*dyPr(j)*dzPr(k)
       endif
      enddo
     enddo
    enddo
    SCALAR=SCALAR/p
   endif
  endif
  
  if (computescalars>0) then
   SCALAR_2=0
   do i=1,computescalars
    call ADVSCALAR(SCALAR_2(:,:,:,i),SCALAR(:,:,:,i),U2,V2,W2,2,1._KND)
   enddo
   SCALAR=SCALAR+SCALAR_2
   call ComputeTDiff(U2,V2,W2)
   call Bound_Visc(TDiff)
    do i=1,computescalars
     call DIFFSCALAR(SCALAR_2(:,:,:,i),SCALAR(:,:,:,i),2,1._KND)
    enddo
    if (computedeposition>0) call Deposition(SCALAR_2,1._KND)
    if (computegravsettling>0) call Gravsettling(SCALAR_2,1._KND)
   SCALAR=SCALAR_2
  endif

  if (buoyancy>0) then
   call ComputeTDiff(U2,V2,W2)
   call Bound_Visc(TDiff)
   temperature2=0
   call Bound_Temp(temperature)
   call CDSSCALAR(temperature2,temperature,U2,V2,W2,1,1._KND)
   temperature=temperature+temperature2
   call Bound_Temp(temperature)
   call DIFFSCALAR(temperature2,temperature,1,1._KND)
   temperature=temperature2
  endif
  
  
  delta=SUM(ABS(U(1:Unx,1:Uny,1:Unz)-U2(1:Unx,1:Uny,1:Unz)))/(Unx*Uny*Unz)
  delta=delta+SUM(ABS(V(1:Vnx,1:Vny,1:Vnz)-V2(1:Vnx,1:Vny,1:Vnz)))/(Vnx*Vny*Vnz)
  delta=delta+SUM(ABS(W(1:Wnx,1:Wny,1:Wnz)-W2(1:Wnx,1:Wny,1:Wnz)))/(Wnx*Wny*Wnz)
  U=U2
  V=V2
  W=W2
    
  endsubroutine TMARCHRK2
  
  
  
  












  
  subroutine TMARCHRK3(U,V,W,Pr,delta)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),intent(out):: delta

  real(KND),allocatable,dimension(:,:,:),save:: Q
  real(KND),dimension(LBOUND(U,1):UBOUND(U,1),LBOUND(U,2):UBOUND(U,2),&
                                              LBOUND(U,3):UBOUND(U,3))::U2,Ustar
  real(KND),dimension(LBOUND(V,1):UBOUND(V,1),LBOUND(V,2):UBOUND(V,2),&
                                              LBOUND(V,3):UBOUND(V,3))::V2,Vstar
  real(KND),dimension(LBOUND(W,1):UBOUND(W,1),LBOUND(W,2):UBOUND(W,2),&
                                              LBOUND(W,3):UBOUND(W,3))::W2,Wstar
  real(KND),dimension(LBOUND(SCALAR,1):UBOUND(SCALAR,1),LBOUND(SCALAR,2):UBOUND(SCALAR,2),&
                LBOUND(SCALAR,3):UBOUND(SCALAR,3),LBOUND(SCALAR,4):UBOUND(SCALAR,4))::SCALAR_adv,SCALAR_2
  real(KND),dimension(LBOUND(TEMPERATURE,1):UBOUND(TEMPERATURE,1),LBOUND(TEMPERATURE,2):UBOUND(TEMPERATURE,2),&
   LBOUND(TEMPERATURE,3):UBOUND(TEMPERATURE,3))::TEMPERATURE_adv,TEMPERATURE2
  real(KND),dimension(1:3),save:: alpha,beta,rho
  integer i,j,k,l
  real(KND) p
  integer,save:: called=0


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

  call BOUND_CONDU(U)
  call BOUND_CONDV(V)
  call BOUND_CONDW(W)
  if (buoyancy==1)  call Bound_Temp(temperature)

  if ((BtypeW==TURBULENTINLET).or.(BtypeE==TURBULENTINLET)) call GETTURBINLET
  if (BtypeW==INLETFROMFILE) call GETINLETFROMFILE(time)

      call timestepEUL(U,V,W)
  if (steady==0.and.dt+time>endtime)  dt=endtime-time
 
  write (*,*) "time:",time,"dt: ",dt

  temperature_adv=0
  Ustar=0
  Vstar=0
  Wstar=0

  do l=1,3
  
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
   else
    call CDU2(Ustar,U,V,W,1._KND)
    call CDV2(Vstar,U,V,W,1._KND)
    call CDW2(Wstar,U,V,W,1._KND)
   endif
   call CoriolisForce(Ustar,Vstar,U,V,1._KND)
   if (buoyancy==1) call BuoyancyForce(Wstar,temperature,1._KND)

   U2=U2+Ustar*beta(l)
   V2=V2+Vstar*beta(l)
   W2=W2+Wstar*beta(l)
  endif



  call OTHERTERMS(U,V,W,U2,V2,W2,Pr,2.*alpha(l))

  if (BtypeT==FREESLIP)  call ATTENUATETOP(U2,V2,W2,Pr)
  if (BtypeE==NEUMANN) then
    if (buoyancy==1) then
      call ATTENUATEOUT(U2,V2,W2,Pr,temperature)
    else
      call ATTENUATEOUT(U2,V2,W2,Pr)
    endif
  endif

  if (masssourc==1) then
      call MASS_SOURC(Q,U2,V2,W2)
  endif
   
  
  call BOUND_CONDU(U2)
  call BOUND_CONDV(V2)
  call BOUND_CONDW(W2)


  if (poissmet>0) then
   if (masssourc==1) then
     call PR_CORRECT(U2,V2,W2,Pr,2.*alpha(l),Q)
   else
     call PR_CORRECT(U2,V2,W2,Pr,2.*alpha(l))
   endif  
  endif


  if (computescalars>0.and..not.released) call EXPLOSION

  if (computescalars>0) then
   SCALAR_2=0
   if (l>1) then
    SCALAR_2=SCALAR_2+SCALAR_adv*rho(l)
   endif
   SCALAR_adv=0
   do i=1,computescalars
    call KAPPASCALAR(SCALAR_adv(:,:,:,i),SCALAR(:,:,:,i),U2,V2,W2,2,1._KND)
   enddo
   SCALAR_2=SCALAR_2+SCALAR_adv*beta(l)

   if (pointscalsource==1) then
    do i=1,computescalars
     SCALAR_2(scalsrci(i),scalsrcj(i),scalsrck(i),i)=SCALAR_2(scalsrci(i),scalsrcj(i),scalsrck(i),i)+&
      percdistrib(i)*(rho(l)+beta(l))*dt*totalscalsource/(dxPr(scalsrci(i))*dyPr(scalsrcj(i))*dzPr(scalsrck(i)))
    enddo
   endif

   SCALAR=SCALAR+SCALAR_2
   call ComputeTDiff(U2,V2,W2)
   call Bound_Visc(TDiff)
   do i=1,computescalars
      call DIFFSCALAR(SCALAR_2(:,:,:,i),SCALAR(:,:,:,i),2,2._KND*alpha(l))
   enddo
   if (computedeposition>0) call Deposition(SCALAR_2,2._KND*alpha(l))
   if (computegravsettling>0) call Gravsettling(SCALAR_2,2._KND*alpha(l))
   SCALAR=SCALAR_2
  endif

  if (buoyancy>0) then
   call ComputeTDiff(U2,V2,W2)
   call Bound_Visc(TDiff)
   temperature2=0
   call Bound_Temp(temperature)

   if (l>1) then
    temperature2=temperature2+temperature_adv*rho(l)
   endif
   temperature_adv=0
   call KAPPASCALAR(temperature_adv,temperature,U2,V2,W2,1,1._KND)
   temperature2=temperature2+temperature_adv*beta(l)

   temperature=temperature+temperature2
   call Bound_Temp(temperature)
    call DIFFSCALAR(temperature2,temperature,1,2._KND*alpha(l))
   temperature=temperature2
   endif



   if (l==1) delta=0
   delta=delta+SUM(ABS(U(1:Unx,1:Uny,1:Unz)-U2(1:Unx,1:Uny,1:Unz)))/(Unx*Uny*Unz)
   delta=delta+SUM(ABS(V(1:Vnx,1:Vny,1:Vnz)-V2(1:Vnx,1:Vny,1:Vnz)))/(Vnx*Vny*Vnz)
   delta=delta+SUM(ABS(W(1:Wnx,1:Wny,1:Wnz)-W2(1:Wnx,1:Wny,1:Wnz)))/(Wnx*Wny*Wnz)



   where(Utype>0) U2=0
   where(Vtype>0) V2=0
   where(Wtype>0) W2=0

   call NULLINTERIOR(U2,V2,W2)
   U=U2
   V=V2
   W=W2
  enddo
  end subroutine TMARCHRK3   
    
  


  subroutine ComputeTDiff(U,V,W)
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer:: i,j,k
   if (Re>0) then
    forall(k=1:Prnz,j=1:Prny,i=1:Prnx)
     TDiff(i,j,k)=(Visc(i,j,k)-1._KND/Re)/Prt(i,j,k,U,V,temperature)+(1._KND/(Re*Prandtl))
    end forall
   else
    forall(k=1:Prnz,j=1:Prny,i=1:Prnx)
     TDiff(i,j,k)=Visc(i,j,k)/Prt(i,j,k,U,V,temperature)
    end forall
   endif
  end subroutine ComputeTDiff



  
   
  subroutine TMARCHSHIFTINLET(U,V,W,Pr,delta) !Only shifts inlet in the x direction, for debugging purposes
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),intent(out):: delta
  integer i,j,k
  real(KND),allocatable,dimension(:,:,:,:),save:: SCALAR_2
  real(KND),allocatable,dimension(:,:,:),save:: temperature2
  integer,save:: called=0


 if (called==0) then
   called=1
   if (computescalars>0) then
    allocate(SCALAR_2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
   endif
   if (buoyancy>0) then
    allocate(temperature2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
   endif 
  endif 
 
  if ((BtypeW==TURBULENTINLET).or.(BtypeE==TURBULENTINLET)) call GETTURBINLET
  if (BtypeW==INLETFROMFILE) call GETINLETFROMFILE(time)
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
    call PR_CORRECT(U,V,W,Pr,1._KND)
  endif

  if (computescalars>0)then
   if (time>5.and.maxval(scalar(:,:,:,1))==0) then
    call GridCoords(i,j,k,0._KND,0._KND,0.5_KND)
    SCALAR(i,j,k,1)=1._KND/(dxPr(i)*dyPr(j)*dzPr(k))
   endif
  endif
  


   if (wallmodeltype>0) then
                   call ComputeViscsWM(U,V,W,Pr)
   endif

  delta=1
  endsubroutine TMARCHSHIFTINLET
  
  
  
  subroutine EXPLOSION
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
    SCALAR=0
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
         SCALAR(xi,yj,zk,:)=SCALAR(xi,yj,zk,:)+percdistrib(:)*0.2
         p=p+1
        elseif (zp<ct*0.4) then
         SCALAR(xi,yj,zk,:)=SCALAR(xi,yj,zk,:)+percdistrib(:)*0.8
         p=p+1
        elseif (zp<ct*0.6) then
         SCALAR(xi,yj,zk,:)=SCALAR(xi,yj,zk,:)+percdistrib(:)*1.25
         p=p+1
        elseif (zp<ct*0.8) then
         SCALAR(xi,yj,zk,:)=SCALAR(xi,yj,zk,:)+percdistrib(:)*1.75
         p=p+1
        elseif (zp<=ct) then
         SCALAR(xi,yj,zk,:)=SCALAR(xi,yj,zk,:)+percdistrib(:)*1.1
         p=p+1
        endif
       endif
      enddo
     enddo
    enddo
    SCALAR=totalscalsource*SCALAR/p
    SCALAR=SCALAR/(dxmin*dymin*dzmin)
    released=.true.
   endif
  endif

  
  endsubroutine EXPLOSION
  
  
  
  


  subroutine MOMSOURC(U,V,W)
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
!     if (IBP%interp<=0) then
!      intvel=0
!     elseif (IBP%interp==1) then
!       n=0
!       if (dirx/=0)  n=n+1
!       if (diry/=0)  n=n+1
!       if (dirz/=0)  n=n+1
!       if (n/=1) then
!                  write (*,*) dirx,diry,dirz
!                  stop
!        endif
! 
!      if (IBP%component==1) then
!       vel1=U(x+dirx,y+diry,z+dirz)
!       vel2=U(x+2*dirx,y+2*diry,z+2*dirz)
!       if (IBP%interpdir==1) then
!        p0=abs(IBP%distx)
!        p1=abs(xU(x+dirx)-xU(x))
!        p2=abs(xU(x+2*dirx)-xU(x))
!       elseif (IBP%interpdir==2) then
!        p0=abs(IBP%disty)
!        p1=abs(yPr(y+diry)-yPr(y))
!        p2=abs(yPr(y+2*diry)-yPr(y))
!       else
!        p0=abs(IBP%distz)
!        p1=abs(zPr(z+dirz)-zPr(z))
!        p2=abs(zPr(z+2*dirz)-zPr(z))
!       endif 
!     elseif (IBP%component==2) then
!       vel1=V(x+dirx,y+diry,z+dirz)
!       vel2=V(x+2*dirx,y+2*diry,z+2*dirz)
!       if (IBP%interpdir==1) then
!        p0=abs(IBP%distx)
!        p1=abs(xPr(x+dirx)-xPr(x))
!        p2=abs(xPr(x+2*dirx)-xPr(x))
!       elseif (IBP%interpdir==2) then
!        p0=abs(IBP%disty)
!        p1=abs(yV(y+diry)-yV(y))
!        p2=abs(yV(y+2*diry)-yV(y))
!       else
!        p0=abs(IBP%distz)
!        p1=abs(zPr(z+dirz)-zPr(z))
!        p2=abs(zPr(z+2*dirz)-zPr(z))
!      endif 
!      else
!       vel1=W(x+dirx,y+diry,z+dirz)
!       vel2=W(x+2*dirx,y+2*diry,z+2*dirz)
!       if (IBP%interpdir==1) then
!        p0=abs(IBP%distx)
!        p1=abs(xPr(x+dirx)-xPr(x))
!        p2=abs(xPr(x+2*dirx)-xPr(x))
!       elseif (IBP%interpdir==2) then
!        p0=abs(IBP%disty)
!        p1=abs(yPr(y+diry)-yPr(y))
!        p2=abs(yPr(y+2*diry)-yPr(y))
!       else
!        p0=abs(IBP%distz)
!        p1=abs(zW(z+dirz)-zW(z))
!        p2=abs(zW(z+2*dirz)-zW(z))
!       endif
!      endif
! 
!      intvel=IBLinInt(p0,p1,p2,vel1,vel2)*.5
! 
! 
!     elseif (IBP%interp==2) then
!      if (IBP%component==1) then
!       if (IBP%interpdir==1) then
!        vel1=U(x,y+diry,z)
!        vel2=U(x,y,z+dirz)
!        vel3=U(x,y+diry,z+dirz)
!        p0=abs(IBP%disty)
!        p1=abs(IBP%distz)
!        p2=abs(yPr(y+diry)-yPr(y))
!        p3=abs(zPr(z+dirz)-zPr(z))
!       elseif (IBP%interpdir==2) then
!        vel1=U(x,y,z+dirz)
!        vel2=U(x+dirx,y,z)
!        vel3=U(x+dirx,y,z+dirz)
!        p0=abs(IBP%distz)
!        p1=abs(IBP%distx)
!        p2=abs(zPr(z+dirz)-zPr(z))
!        p3=abs(xU(x+dirx)-xU(x))
!       else
!        vel1=U(x+dirx,y,z)
!        vel2=U(x,y+diry,z)
!        vel3=U(x+dirx,y+diry,z)
!        p0=abs(IBP%distx)
!        p1=abs(IBP%disty)
!        p2=abs(xU(x+dirx)-xU(x))
!        p3=abs(yPr(y+diry)-yPr(y))
!       endif
!      elseif (IBP%component==2) then
!       if (IBP%interpdir==1) then
!        vel1=V(x,y+diry,z)
!        vel2=V(x,y,z+dirz)
!        vel3=V(x,y+diry,z+dirz)
!        p0=abs(IBP%disty)
!        p1=abs(IBP%distz)
!        p2=abs(yV(y+diry)-yV(y))
!        p3=abs(zPr(z+dirz)-zPr(z))
!       elseif (IBP%interpdir==2) then
!        vel1=V(x,y,z+dirz)
!        vel2=V(x+dirx,y,z)
!        vel3=V(x+dirx,y,z+dirz)
!        p0=abs(IBP%distz)
!        p1=abs(IBP%distx)
!        p2=abs(zPr(z+dirz)-zPr(z))
!        p3=abs(xPr(x+dirx)-xPr(x))
!       else
!        vel1=V(x+dirx,y,z)
!        vel2=V(x,y+diry,z)
!        vel3=V(x+dirx,y+diry,z)
!        p0=abs(IBP%distx)
!        p1=abs(IBP%disty)
!        p2=abs(xPr(x+dirx)-xPr(x))
!        p3=abs(yV(y+diry)-yV(y))
!       endif
!      else
!       if (IBP%interpdir==1) then
!        vel1=W(x,y+diry,z)
!        vel2=W(x,y,z+dirz)
!        vel3=W(x,y+diry,z+dirz)
!        p0=abs(IBP%disty)
!        p1=abs(IBP%distz)
!        p2=abs(yPr(y+diry)-yPr(y))
!        p3=abs(zW(z+dirz)-zW(z))
!       elseif (IBP%interpdir==2) then
!        vel1=W(x,y,z+dirz)
!        vel2=W(x+dirx,y,z)
!        vel3=W(x+dirx,y,z+dirz)
!        p0=abs(IBP%distz)
!        p1=abs(IBP%distx)
!        p2=abs(zW(z+dirz)-zW(z))
!        p3=abs(xPr(x+dirx)-xPr(x))
!       else
!        vel1=W(x+dirx,y,z)
!        vel2=W(x,y+diry,z)
!        vel3=W(x+dirx,y+diry,z)
!        p0=abs(IBP%distx)
!        p1=abs(IBP%disty)
!        p2=abs(xPr(x+dirx)-xPr(x))
!        p3=abs(yPr(y+diry)-yPr(y))
!       endif
!      endif
! 
!      intvel=IBBiLinInt(p0,p1,p2,p3,vel1,vel2,vel3)*.5
! 
! 
!     elseif (IBP%interp==3) then
!      if (IBP%component==1) then
!        vel1=U(x+dirx,y,z)
!        vel2=U(x,y+diry,z)
!        vel3=U(x+dirx,y+diry,z)
!        vel4=U(x,y,z+dirz)
!        vel5=U(x+dirx,y,z+dirz)
!        vel6=U(x,y+diry,z+dirz)
!        vel7=U(x+dirx,y+diry,z+dirz)
!        p0=abs(IBP%distx)
!        p1=abs(IBP%disty)
!        p2=abs(IBP%distz)
!        p3=abs(xU(x+dirx)-xU(x))
!        p4=abs(yPr(y+diry)-yPr(y))
!        p5=abs(zPr(z+dirz)-zPr(z))
!      elseif (IBP%component==2) then
!        vel1=V(x+dirx,y,z)
!        vel2=V(x,y+diry,z)
!        vel3=V(x+dirx,y+diry,z)
!        vel4=V(x,y,z+dirz)
!        vel5=V(x+dirx,y,z+dirz)
!        vel6=V(x,y+diry,z+dirz)
!        vel7=V(x+dirx,y+diry,z+dirz)
!        p0=abs(IBP%distx)
!        p1=abs(IBP%disty)
!        p2=abs(IBP%distz)
!        p3=abs(xPr(x+dirx)-xPr(x))
!        p4=abs(yV(y+diry)-yV(y))
!        p5=abs(zPr(z+dirz)-zPr(z))
!      else      
!        vel1=W(x+dirx,y,z)
!        vel2=W(x,y+diry,z)
!        vel3=W(x+dirx,y+diry,z)
!        vel4=W(x,y,z+dirz)
!        vel5=W(x+dirx,y,z+dirz)
!        vel6=W(x,y+diry,z+dirz)
!        vel7=W(x+dirx,y+diry,z+dirz)
!        p0=abs(IBP%distx)
!        p1=abs(IBP%disty)
!        p2=abs(IBP%distz)
!        p3=abs(xPr(x+dirx)-xPr(x))
!        p4=abs(yPr(y+diry)-yPr(y))
!        p5=abs(zW(z+dirz)-zW(z))
!      endif
! 
!      intvel=IBTriLinInt(p0,p1,p2,p3,p4,p5,vel1,vel2,vel3,vel4,vel5,vel6,vel7)
! 
!     else
!      intvel=0
!     endif
 intvel=0
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

  end subroutine MOMSOURC
   

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
     Q(x,y,z)=0
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
     p=MAX(ABS(U(i,j,k)/dxU(i)),ABS(U(i-1,j,k)/dxU(i-1)))+&
        MAX(ABS(V(i,j,k)/dyV(j)),ABS(V(i,j-1,k)/dyV(j-1)))+&
        MAX(ABS(W(i,j,k)/dzW(k)),ABS(W(i,j,k-1)/dzW(k-1)))
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
     p=MAX(MAX(ABS(U(i,j,k)/dxU(i)),ABS(U(i-1,j,k)/dxU(i-1))),MAX(ABS(V(i,j,k)/dyV(j)),ABS(V(i,j-1,k)/dyV(j-1))))&
         +MAX(ABS(W(i,j,k)/dzW(k)),ABS(W(i,j,k-1)/dzW(k-1)))
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







  subroutine BUOYANCYFORCE(W2,theta,coef)
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
  endsubroutine BUOYANCYFORCE

  subroutine CORIOLISFORCE(U2,V2,U,V,coef)
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
  endsubroutine CORIOLISFORCE







  subroutine OTHERTERMS(U,V,W,U2,V2,W2,Pr,coef)
   real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
   real(KND),intent(inout):: Pr(1:,1:,1:)
   real(KND),intent(inout):: U2(-2:,-2:,-2:),V2(-2:,-2:,-2:),W2(-2:,-2:,-2:)
   real(KND),intent(in):: coef
       
   real(KND),dimension(LBOUND(U,1):UBOUND(U,1),LBOUND(U,2):UBOUND(U,2),LBOUND(U,3):UBOUND(U,3)):: U3, Apu
   real(KND),dimension(LBOUND(V,1):UBOUND(V,1),LBOUND(V,2):UBOUND(V,2),LBOUND(V,3):UBOUND(V,3)):: V3, ApV
   real(KND),dimension(LBOUND(W,1):UBOUND(W,1),LBOUND(W,2):UBOUND(W,2),LBOUND(W,3):UBOUND(W,3)):: W3, ApW
   real(KND) p,S,Su,Sv,Sw,Suavg,Svavg,Swavg,Af,Ap,Apre,Aprn,Aprt,Ap2,Str(3,3)

   integer i,j,k,l,x,y,z

   type(TIBPoint),pointer:: IBP
   
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


   call BOUND_CONDU2(U2)
   call BOUND_CONDV2(V2)
   call BOUND_CONDW2(W2)

  if (Re>0) then
   !Diffusion using Crank Nicolson 
   !first approximation using forward Euler
   !iteration SOR or Gauss-Seidel
   call BOUND_CONDU(U)
   call BOUND_CONDV(V)
   call Bound_CONDW(W)
   U3=U+U2
   V3=V+V2
   W3=W+W2
   if (sgstype==1) then
                     call Smag(U,V,W)
   elseif (sgstype==2) then
                     call DynSmag(U,V,W)
   elseif (sgstype==3) then
                     call VREMAN(U,V,W)
   else 
                     Visc=1._KND/Re
   endif

   if (debuglevel>0) then               
    write(*,*) "NUt", SUM(Visc(1:Prnx,1:Prny,1:Prnz))/(Prnx*Prny*Prnz)
    write(*,*) "maxNUt", MAXVAL(Visc(1:Prnx,1:Prny,1:Prnz))
    write(*,*) "minNUt", MINVAL(Visc(1:Prnx,1:Prny,1:Prnz))
   endif

   if (wallmodeltype>0) then
                   call ComputeViscsWM(U,V,W,Pr)
   endif
   call ScalFlSourc(Visc,3)

   Ap=coef*dt
   if (fullstress==1) then
    do k=0,Prnz+1
     do j=0,Prny+1
      do i=0,Prnx+1
       call STRAINIJ(i,j,k,U,V,W,Str)
       tstress(:,:,i,j,k)=-2*Str(:,:)*(Visc(i,j,k)-1./Re)*0.5
      enddo
     enddo
    enddo
    Visc=1./Re+0.5*(Visc-1./Re)
    call Bound_Visc(Visc)
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
         U2(i,j,k)=U2(i,j,k)+Ap*(tstress(1,1,i+1,j,k)-tstress(1,1,i,j,k))/dxU(i)
         U2(i,j,k)=U2(i,j,k)+Ap*(tstress(1,2,i+1,j+1,k)+tstress(1,2,i,j+1,k)&
                                -tstress(1,2,i+1,j-1,k)-tstress(1,2,i,j-1,k))/(2._KND*(dyV(j)+dyV(j-1)))
         U2(i,j,k)=U2(i,j,k)+Ap*(tstress(1,3,i+1,j,k+1)+tstress(1,3,i,j,k+1)&
                                -tstress(1,3,i+1,j,k-1)-tstress(1,3,i,j,k-1))/(2._KND*(dzW(k)+dzW(k-1)))
      enddo
     enddo
    enddo
    do k=1,Vnz
     do j=1,Vny
      do i=1,Vnx
         V2(i,j,k)=V2(i,j,k)+Ap*(tstress(2,2,i,j+1,k)-tstress(2,2,i,j,k))/dyV(j)
         V2(i,j,k)=V2(i,j,k)+Ap*(tstress(2,1,i+1,j+1,k)+tstress(2,2,i+1,j,k)&
                                -tstress(2,1,i-1,j+1,k)-tstress(2,2,i-1,j,k))/(2._KND*(dxU(i)+dxU(i-1)))
         V2(i,j,k)=V2(i,j,k)+Ap*(tstress(2,3,i,j+1,k+1)+tstress(2,3,i,j,k+1)&
                                -tstress(2,3,i,j+1,k-1)-tstress(2,3,i,j,k-1))/(2._KND*(dzW(k)+dzW(k-1)))
      enddo
     enddo
    enddo
    do k=1,Wnz
     do j=1,Wny
      do i=1,Wnx
         W2(i,j,k)=W2(i,j,k)+Ap*(tstress(3,3,i,j,k+1)-tstress(3,3,i,j,k))/dzW(k)
         W2(i,j,k)=W2(i,j,k)+Ap*(tstress(3,2,i,j+1,k+1)+tstress(3,2,i,j+1,k)&
                                -tstress(3,2,i,j-1,k+1)-tstress(3,2,i,j-1,k))/(2._KND*(dyV(j)+dyV(j-1)))
         W2(i,j,k)=W2(i,j,k)+Ap*(tstress(3,1,i+1,j,k+1)+tstress(3,1,i+1,j,k)&
                                -tstress(3,1,i-1,j,k+1)-tstress(3,1,i-1,j,k))/(2._KND*(dxU(i)+dxU(i-1)))
      enddo
     enddo
    enddo
    call Bound_CondU(U2)
    call Bound_CondV(V2)
    call Bound_CondW(W2)
   endif

   do k=1,Unz
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
  call MOMSOURC(U3,V3,W3)

   if (associated(FirstIBPoint)) then
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


  !now in W3 is first approximation
  impldiff=1
  if (impldiff==1) then

   if (gridtype==UNIFORMGRID) then
    call UNIFREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
   else

    !Gauss-Seidel iteration for Crank-Nicolson result
    Ap2=coef*dt/(2._KND)
   S=0
   l=0



   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
      U2(i,j,k)=U2(i,j,k)+Ap2*(&
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
      V2(i,j,k)=V2(i,j,k)+Ap2*(&
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
      W2(i,j,k)=W2(i,j,k)+Ap2*(&
      ((0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W(i+1,j,k)-W(i,j,k))/dxU(i)-&
      0.25_KND*(Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(W(i,j,k)-W(i-1,j,k))/dxU(i-1))/dxPr(i)+&
       (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W(i,j+1,k)-W(i,j,k))/dyV(j)-&
       0.25_KND*(Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(W(i,j,k)-W(i,j-1,k))/dyV(j-1))/dyPr(j)+&
       (Visc(i,j,k+1)*(W(i,j,k+1)-W(i,j,k))/dzPr(k+1)-&
       Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1))/dzPr(k))/dzW(k)))
     enddo
    enddo
   enddo

     do k=1,Unz
      do j=1,Uny
       do i=1,Unx
        ApU(i,j,k)=1.+Ap2* ((Visc(i+1,j,k)/dxPr(i+1)+&
                    Visc(i,j,k)/dxPr(i))/dxU(i)+&
                    (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))/dyV(j)+&
                    0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))/dyV(j-1))/dyPr(j)+&
                    (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))/dzW(k)+&
                    0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))/dzW(k-1))/dzPr(k))
       enddo
      enddo
     enddo
     do k=1,Vnz
      do j=1,Vny
       do i=1,Vnx
        ApV(i,j,k)=1.+Ap2* ((0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))/dxU(i)+&
                    0.25_KND*(Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))/dxU(i-1))/dxPr(i)+&
                   (Visc(i,j+1,k)/dyPr(j+1)+&
                   Visc(i,j,k)/dyPr(j))/dyV(j)+&
                   (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))/dzW(k)+&
                   0.25_KND*(Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))/dzW(k-1))/dzPr(k))
       enddo
      enddo
     enddo
     do k=1,Wnz
      do j=1,Wny
       do i=1,Wnx
        ApW(i,j,k)=1.+Ap2* ((0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))/dxU(i)+&
                    0.25_KND*(Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))/dxU(i-1))/dxPr(i)+&
                   (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))/dyV(j)+&
                   0.25_KND*(Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))/dyV(j-1))/dyPr(j)+&
                   (Visc(i,j,k+1)/dzPr(k+1)+&
                   Visc(i,j,k)/dzPr(k))/dzW(k))
       enddo
      enddo
     enddo

    Suavg=abs(MAXVAL(U3(1:Unx,1:Uny,1:Unz)))
    Svavg=abs(MAXVAL(V3(1:Vnx,1:Vny,1:Vnz)))
    Swavg=abs(MAXVAL(W3(1:Wnx,1:Wny,1:Wnz)))
    if (Suavg<=1e-3_KND) Suavg=1
    if (Svavg<=1e-3_KND) Svavg=1
    if (Swavg<=1e-3_KND) Swavg=1
    do l=1,maxCNiter
      call BOUND_CONDU(U3)
      call BOUND_CONDV(V3)
      call BOUND_CONDW(W3)
     S=0
     Su=0
     Sv=0
     Sw=0
     !$OMP PARALLEL PRIVATE(i,j,k,p) REDUCTION(max:Su,Sv,Sw)
     !$OMP DO
     do k=1,Unz
      do j=1,Uny
       do i=1,Unx
        if (REDBLACKU(i,j,k)) then
         p=((Visc(i+1,j,k)*(U3(i+1,j,k))/dxPr(i+1)-&
          Visc(i,j,k)*(-U3(i-1,j,k))/dxPr(i))/dxU(i)+&
          (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U3(i,j+1,k))/dyV(j)-&
          0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(-U3(i,j-1,k))/dyV(j-1))/dyPr(j)+&
          (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U3(i,j,k+1))/dzW(k)-&
          0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-U3(i,j,k-1))/dzW(k-1))/dzPr(k))
         p=Ap2*p+U2(i,j,k)+U(i,j,k)
         p=p/ApU(i,j,k)

        
         Su=max(Su,abs(p-U3(i,j,k)))
         U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
        endif
       enddo
      enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP DO
     do k=1,Vnz
      do j=1,Vny
       do i=1,Vnx
        if (REDBLACKV(i,j,k)) then
         p=((0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V3(i+1,j,k))/dxU(i)-&
          0.25_KND*(Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(-V3(i-1,j,k))/dxU(i-1))/dxPr(i)+&     
          (Visc(i,j+1,k)*(V3(i,j+1,k))/dyPr(j+1)-&
          Visc(i,j,k)*(-V3(i,j-1,k))/dyPr(j))/dyV(j)+&
          (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V3(i,j,k+1))/dzW(k)-&
          0.25_KND*(Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-V3(i,j,k-1))/dzW(k-1))/dzPr(k))
         p=Ap2*p+V2(i,j,k)+V(i,j,k)
         p=p/ApV(i,j,k)
         Sv=max(Sv,abs(p-V3(i,j,k)))
         V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
        endif
       enddo
      enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP DO
     do k=1,Wnz
      do j=1,Wny
       do i=1,Wnx
        if (REDBLACKW(i,j,k)) then
         p=((0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W3(i+1,j,k))/dxU(i)-&
          0.25_KND*(Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(-W3(i-1,j,k))/dxU(i-1))/dxPr(i)+&     
          (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W3(i,j+1,k))/dyV(j)-&
          0.25_KND*(Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(-W3(i,j-1,k))/dyV(j-1))/dyPr(j)+&
          (Visc(i,j,k+1)*(W3(i,j,k+1))/dzPr(k+1)-&
          Visc(i,j,k)*(-W3(i,j,k-1))/dzPr(k))/dzW(k))
         p=Ap2*p+W2(i,j,k)+W(i,j,k)
         p=p/ApW(i,j,k)
         Sw=max(Sw,abs(p-W3(i,j,k)))
         W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
        endif
       enddo
      enddo
     enddo
     !$OMP ENDDO

     !$OMP DO
     do k=1,Unz
      do j=1,Uny
       do i=1,Unx
        if (.not.REDBLACKU(i,j,k)) then
         p=((Visc(i+1,j,k)*(U3(i+1,j,k))/dxPr(i+1)-&
          Visc(i,j,k)*(-U3(i-1,j,k))/dxPr(i))/dxU(i)+&
          (0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U3(i,j+1,k))/dyV(j)-&
          0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(-U3(i,j-1,k))/dyV(j-1))/dyPr(j)+&
          (0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U3(i,j,k+1))/dzW(k)-&
          0.25_KND*(Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-U3(i,j,k-1))/dzW(k-1))/dzPr(k))
         p=Ap2*p+U2(i,j,k)+U(i,j,k)
         p=p/ApU(i,j,k)

        
         Su=max(Su,abs(p-U3(i,j,k)))
         U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))!*1.72_KND
        endif
       enddo
      enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP DO
     do k=1,Vnz
      do j=1,Vny
       do i=1,Vnx
        if (.not.REDBLACKV(i,j,k)) then
         p=((0.25_KND*(Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V3(i+1,j,k))/dxU(i)-&
          0.25_KND*(Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(-V3(i-1,j,k))/dxU(i-1))/dxPr(i)+&     
          (Visc(i,j+1,k)*(V3(i,j+1,k))/dyPr(j+1)-&
          Visc(i,j,k)*(-V3(i,j-1,k))/dyPr(j))/dyV(j)+&
          (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V3(i,j,k+1))/dzW(k)-&
          0.25_KND*(Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(-V3(i,j,k-1))/dzW(k-1))/dzPr(k))
         p=Ap2*p+V2(i,j,k)+V(i,j,k)
         p=p/ApV(i,j,k)
         Sv=max(Sv,abs(p-V3(i,j,k)))
         V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))!*1.72_KND
        endif
       enddo
      enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP DO
     do k=1,Wnz
      do j=1,Wny
       do i=1,Wnx
        if (.not.REDBLACKW(i,j,k)) then
         p=((0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W3(i+1,j,k))/dxU(i)-&
          0.25_KND*(Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(-W3(i-1,j,k))/dxU(i-1))/dxPr(i)+&     
          (0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W3(i,j+1,k))/dyV(j)-&
          0.25_KND*(Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(-W3(i,j-1,k))/dyV(j-1))/dyPr(j)+&
          (Visc(i,j,k+1)*(W3(i,j,k+1))/dzPr(k+1)-&
          Visc(i,j,k)*(-W3(i,j,k-1))/dzPr(k))/dzW(k))
         p=Ap2*p+W2(i,j,k)+W(i,j,k)
         p=p/ApW(i,j,k)
         Sw=max(Sw,abs(p-W3(i,j,k)))
         W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))!*1.72_KND
        endif
       enddo
      enddo
     enddo
     !$OMP ENDDO
     !$OMP END PARALLEL
     S=max(Su/Suavg,Sv/Svavg,Sw/Swavg)
     write (*,*) "CN ",l,S
     if (S<=epsCN) exit
    enddo
   endif

   U2=U3
   V2=V3
   W2=W3
  
   else

    U2=U3
    V2=V3
    W2=W3
   endif

   else

  U2=U+U2
  V2=V+V2
  W2=W+W2
  endif

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
  if (debuglevel>=2) then
   write(*,*) "Mean friction:", S 
  endif
  endsubroutine OTHERTERMS


  subroutine UNIFREDBLACK(U,V,W,U2,V2,W2,U3,V3,W3,coef)
  real(KND),dimension(-2:,-2:,-2:),intent(inout):: U,V,W,U2,V2,W2,U3,V3,W3
  real(KND),intent(in):: coef
  real(KND) dxmin2,dymin2,dzmin2,Ap,Ap2,Ap3,Ap4,p,S
  integer i,j,k,l,Unyz,Vnyz,Wnyz

   dxmin2=dxmin*dxmin
   dymin2=dymin*dymin
   dzmin2=dzmin*dzmin
   Unyz=Uny*Unz
   Vnyz=Vny*Vnz
   Wnyz=Wny*Wnz
   Ap4=coef*dt/(2._KND)!*Re)
   Ap3= (2._KND/(dxmin2)+2._KND/(dymin2)+2._KND/(dzmin2))
   S=0
   l=0
    do l=1,maxCNiter
     call BOUND_CONDU(U3)
     call BOUND_CONDV(V3)
     call BOUND_CONDW(W3)
     S=0
     !$OMP PARALLEL PRIVATE(i,j,k,p) REDUCTION(+:S)
     !$OMP DO
     do k=1,Unz
      do j=1,Uny
       do i=1,Unx
        if (REDBLACKU(i,j,k)) then
         Ap2=Ap4/Re
         Ap=1._KND+ Ap2*Ap3
         p=((U(i+1,j,k)-U(i,j,k))-&
            (U(i,j,k)-U(i-1,j,k)))/(dxmin2)+&
          ((U(i,j+1,k)-U(i,j,k))-&
          (U(i,j,k)-U(i,j-1,k)))/(dymin2)+&
          ((U(i,j,k+1)-U(i,j,k))-&
          (U(i,j,k)-U(i,j,k-1)))/(dzmin2)+&
          (U3(i+1,j,k)+&
          U3(i-1,j,k))/(dxmin2)+&
          (U3(i,j+1,k)+&
          U3(i,j-1,k))/(dymin2)+&
          (U3(i,j,k+1)+&
          U3(i,j,k-1))/(dzmin2)
         p=Ap2*p+U2(i,j,k)+U(i,j,k)
         p=p/Ap
         S=S+(p-U3(i,j,k))*(p-U3(i,j,k))
         U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))
        endif 
       enddo
      enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP DO
     do k=1,Vnz
      do j=1,Vny
       do i=1,Vnx
        if (REDBLACKV(i,j,k)) then
         Ap2=Ap4/Re
         Ap=1._KND+ Ap2*Ap3
         p=(V(i+1,j,k)-2*V(i,j,k)+V(i-1,j,k))/(dxmin2)+&
          (V(i,j+1,k)-2*V(i,j,k)+V(i,j-1,k))/(dymin2)+&
          (V(i,j,k+1)-2*V(i,j,k)+V(i,j,k-1))/(dzmin2)+&
          (V3(i+1,j,k)+V3(i-1,j,k))/(dxmin2)+&
          (V3(i,j+1,k)+V3(i,j-1,k))/(dymin2)+&
          (V3(i,j,k+1)+V3(i,j,k-1))/(dzmin2)
         p=Ap2*p+V2(i,j,k)+V(i,j,k)
         p=p/Ap
         S=S+(p-V3(i,j,k))*(p-V3(i,j,k))
         V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))
        endif
       enddo
      enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP DO
     do k=1,Wnz
      do j=1,Wny
       do i=1,Wnx
        if (REDBLACKW(i,j,k)) then
         Ap2=Ap4/Re
         Ap=1._KND+ Ap2*Ap3
         p=(W(i+1,j,k)-2*W(i,j,k)+W(i-1,j,k))/(dxmin2)+&
          (W(i,j+1,k)-2*W(i,j,k)+W(i,j-1,k))/(dymin2)+&
          (W(i,j,k+1)-2*W(i,j,k)+W(i,j,k-1))/(dzmin2)+&
          (W3(i+1,j,k)+W3(i-1,j,k))/(dxmin2)+&
          (W3(i,j+1,k)+W3(i,j-1,k))/(dymin2)+&
          (W3(i,j,k+1)+W3(i,j,k-1))/(dzmin2)
         p=Ap2*p+W2(i,j,k)+W(i,j,k)
         p=p/Ap
         S=S+(p-W3(i,j,k))*(p-W3(i,j,k))
         W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))
        endif
       enddo
      enddo
     enddo
     !$OMP ENDDO

     !$OMP DO
     do k=1,Unz
      do j=1,Uny
       do i=1,Unx 
        if (.not.REDBLACKU(i,j,k)) then
         Ap2=Ap4/Re
         Ap=1._KND+ Ap2*Ap3
         p=(U(i+1,j,k)-2*U(i,j,k)+U(i-1,j,k))/(dxmin2)+&
          (U(i,j+1,k)-2*U(i,j,k)+U(i,j-1,k))/(dymin2)+&
          (U(i,j,k+1)-2*U(i,j,k)+U(i,j,k-1))/(dzmin2)+&
          (U3(i+1,j,k)+U3(i-1,j,k))/(dxmin2)+&
          (U3(i,j+1,k)+U3(i,j-1,k))/(dymin2)+&
          (U3(i,j,k+1)+U3(i,j,k-1))/(dzmin2)
         p=Ap2*p+U2(i,j,k)+U(i,j,k)
         p=p/Ap
         S=S+(p-U3(i,j,k))*(p-U3(i,j,k))
         U3(i,j,k)=U3(i,j,k)+(p-U3(i,j,k))
        endif 
       enddo
      enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP DO
     do k=1,Vnz
      do j=1,Vny
       do i=1,Vnx
        if (.not.REDBLACKV(i,j,k)) then
         Ap2=Ap4/Re
         Ap=1._KND+ Ap2*Ap3
         p=(V(i+1,j,k)-2*V(i,j,k)+V(i-1,j,k))/(dxmin2)+&
          (V(i,j+1,k)-2*V(i,j,k)+V(i,j-1,k))/(dymin2)+&
          (V(i,j,k+1)-2*V(i,j,k)+V(i,j,k-1))/(dzmin2)+&
          (V3(i+1,j,k)+V3(i-1,j,k))/(dxmin2)+&
          (V3(i,j+1,k)+V3(i,j-1,k))/(dymin2)+&
          (V3(i,j,k+1)+V3(i,j,k-1))/(dzmin2)
         p=Ap2*p+V2(i,j,k)+V(i,j,k)
         p=p/Ap
         S=S+(p-V3(i,j,k))*(p-V3(i,j,k))
         V3(i,j,k)=V3(i,j,k)+(p-V3(i,j,k))
        endif
       enddo
      enddo
     enddo
     !$OMP ENDDO NOWAIT
     !$OMP DO
     do k=1,Wnz
      do j=1,Wny
       do i=1,Wnx
        if (.not.REDBLACKW(i,j,k)) then
         Ap2=Ap4/Re
         Ap=1._KND+ Ap2*Ap3
         p=(W(i+1,j,k)-2*W(i,j,k)+W(i-1,j,k))/(dxmin2)+&
          (W(i,j+1,k)-2*W(i,j,k)+W(i,j-1,k))/(dymin2)+&
          (W(i,j,k+1)-2*W(i,j,k)+W(i,j,k-1))/(dzmin2)+&
          (W3(i+1,j,k)+W3(i-1,j,k))/(dxmin2)+&
          (W3(i,j+1,k)+W3(i,j-1,k))/(dymin2)+&
          (W3(i,j,k+1)+W3(i,j,k-1))/(dzmin2)
         p=Ap2*p+W2(i,j,k)+W(i,j,k)
         p=p/Ap
         S=S+(p-W3(i,j,k))*(p-W3(i,j,k))
         W3(i,j,k)=W3(i,j,k)+(p-W3(i,j,k))
        endif
       enddo
      enddo
     enddo
     !$OMP ENDDO
     !$OMP ENDPARALLEL

     S=S/(Prnx*Prny*Prnz)
     if (Sqrt(S)<=epsCN) exit
     write (*,*) "CN ",l,Sqrt(S)
    enddo
  endsubroutine UNIFREDBLACK

  subroutine ATTENUATETOP(U,V,W,Pr)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  integer i,j,k,bufn
  real(KND) p,ze,zs,zb

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
      do i=-1,Vnx+1
       do j=-1,Vny+1
        V(i,j,k)=p+DampF(zb)*(V(i,j,k)-p)
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
      do i=-1,Unx+1
       do j=-1,Uny+1
        U(i,j,k)=p+DampF(zb)*(U(i,j,k)-p)
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
      do i=-1,Prnx+1
       do j=-1,Prny+1
        Temperature(i,j,k)=p+DampF(zb)*(temperature(i,j,k)-p)
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
     do i=-1,Wnx+1
      do j=-1,Wny+1
        W(i,j,k)=p+DampF(zb)*(W(i,j,k)-p)
      enddo
     enddo
    enddo
  endsubroutine ATTENUATETOP


  subroutine ATTENUATEOUT(U,V,W,Pr,temperature)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(KND),optional,intent(inout):: temperature(-1:,-1:,-1:)
  integer i,j,k,bufn
  real(KND) p,xe,xs,xb

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
      do j=-1,Uny+1
        U(i,j,k)=p+DampF(xb)*(U(i,j,k)-p)
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
      do j=-1,Vny+1
        V(i,j,k)=p+DampF(xb)*(V(i,j,k)-p)
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
      do j=-1,Wny+1
        W(i,j,k)=p+DampF(xb)*(W(i,j,k)-p)
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
      do j=-1,Prny+1
        temperature(i,j,k)=p+DampF(xb)*(temperature(i,j,k)-p)
      enddo
     enddo
    enddo
   endif
  endsubroutine ATTENUATEOUT


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

  subroutine NULLINTERIOR(U,V,W)
  real(KND),intent(inout):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer i,j,k

  
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
      if (Utype(i,j,k)>0.and.Utype(i,j,k+1)==Utype(i,j,k).and.Utype(i,j,k-1)==Utype(i,j,k)&
          .and.Utype(i,j-1,k)==Utype(i,j,k).and.Utype(i,j+1,k)==Utype(i,j,k)&
          .and.Utype(i-1,j,k)==Utype(i,j,k).and.Utype(i+1,j,k)==Utype(i,j,k).and.Utype(i,j,k+1)==Utype(i,j,k))  U(i,j,k)=0
     enddo
    enddo
   enddo

   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
      if (Vtype(i,j,k)>0.and.Vtype(i,j,k+1)==Vtype(i,j,k).and.Vtype(i,j,k-1)==Vtype(i,j,k)&
          .and.Vtype(i,j-1,k)==Vtype(i,j,k).and.Vtype(i,j+1,k)==Vtype(i,j,k)&
          .and.Vtype(i-1,j,k)==Vtype(i,j,k).and.Vtype(i+1,j,k)==Vtype(i,j,k).and.Vtype(i,j,k+1)==Vtype(i,j,k))  V(i,j,k)=0
     enddo
    enddo
   enddo

   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
      if (Wtype(i,j,k)>0.and.Wtype(i,j,k+1)==Wtype(i,j,k).and.Wtype(i,j,k-1)==Wtype(i,j,k)&
          .and.Wtype(i,j-1,k)==Wtype(i,j,k).and.Wtype(i,j+1,k)==Wtype(i,j,k)&
          .and.Wtype(i-1,j,k)==Wtype(i,j,k).and.Wtype(i+1,j,k)==Wtype(i,j,k).and.Wtype(i,j,k+1)==Wtype(i,j,k))  W(i,j,k)=0
     enddo
    enddo
   enddo
  endsubroutine NULLINTERIOR


end module TSTEPS

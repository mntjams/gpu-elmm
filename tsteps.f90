module TSTEPS

  use UPWIND
  use CDS, only: CDU, CDV, CDW, KappaU, KappaV, KappaW
  use LAXFRIED
  use LAXWEND
  use PARAMETERS
  use BOUNDARIES, only: BoundU,Bound_Q
  use POISSON, only: Pr_Correct
  use SMAGORINSKY, only: Smag, StabSmag, Vreman
  use SCALARS, only:  Bound_Temp, Bound_Visc, Scalar, percdistrib, AdvScalar,&
     DiffScalar, ComputeTDiff, Deposition, GravSettling, ScalFlSourc, VolScalSource
  use TURBINLET, only: GetTurbInlet, GetInletFromFile
  use GEOMETRIC, only: IBLinInt,IBBiLinInt,IBTriLinInt
  use Wallmodels, only: ComputeViscsWM

  implicit none


  private
  public TMarchRK3

  logical:: released=.false.

contains






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
  real(KND),dimension(1:3),parameter:: alpha = (/ 4._KND/15._KND, 1._KND/15._KND, 1._KND/6._KND /)
  real(KND),dimension(1:3),parameter:: beta  = (/ 8._KND/15._KND, 5._KND/12._KND, 3._KND/4._KND /)
  real(KND),dimension(1:3),parameter:: rho   = (/       0._KND, -17._KND/60._KND,-5._KND/12._KND/)
  integer i,j,k,l
  real(KND) p
  integer,save:: called=0
  real time1,time2


  if (called==0) then
   called=1
   if (masssourc==1) allocate(Q(0:Prnx+1,0:Prny+1,0:Prnz+1))
  endif

  call BoundU(1,U)
  call BoundU(2,V)
  call BoundU(3,W)

  if (buoyancy==1)  call Bound_Temp(temperature)

  if ((Btype(We)==TurbulentInlet).or.(Btype(Ea)==TurbulentInlet)) call GetTurbInlet
  if (Btype(We)==InletFromFile) call GetInletFromFile(time)

  call timestepEUL(U,V,W)

  if (steady==0.and.dt+time>endtime)  dt=endtime-time
  write (*,*) "time:",time,"dt: ",dt

  temperature_adv=0
  Ustar=0
  Vstar=0
  Wstar=0

  do l=1,3
    !$hmpp <tsteps> allocate
    !$hmpp <tsteps> advancedload, args[Vreman::Prnx,Vreman::Prny,Vreman::Prnz]
    !$hmpp <tsteps> advancedload, args[Vreman::Unx,Vreman::Uny,Vreman::Unz,&
    !$hmpp <tsteps>     Vreman::Vnx,Vreman::Vny,Vreman::Vnz,Vreman::Wnx,Vreman::Wny,Vreman::Wnz]
    !$hmpp <tsteps> advancedload, args[Vreman::dx,Vreman::dy,Vreman::dz,Vreman::dt,Vreman::Re]


    if (debugparam>1) call cpu_time(time1)

    call BoundU(1,U)
    call BoundU(2,V)
    call BoundU(3,W)

    !$hmpp <tsteps> advancedload, args[Vreman::U,Vreman::V,Vreman::W]
    call SubgridStresses(U,V,W,Pr)

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


       !$hmpp <tsteps> CDS callsite, args[*].noupdate=true
       call CDS_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxmin,dymin,dzmin,dt,Ustar,Vstar,Wstar,U,V,W)

       !$hmpp <tsteps> delegatedstore, args[CDS::U2,CDS::V2,CDS::W2]

       write(*,*) "Back from  GPU CDS"
      else
       call CDU(Ustar,U,V,W,1._KND)
       call CDV(Vstar,U,V,W,1._KND)
       call CDW(Wstar,U,V,W,1._KND)
      endif

write (*,*) "sum Ustar",sum(Ustar(1:Unx,1:Uny,1:Unz))
write (*,*) "sum Vstar",sum(Vstar(1:Vnx,1:Vny,1:Vnz))
write (*,*) "sum Wstar",sum(Wstar(1:Wnx,1:Wny,1:Wnz))

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
   !$hmpp <tsteps> release



    if (debugparam>1) then
     call cpu_time(time2)
     write (*,*) "ET of part 2", (time2-time1)
     time1=time2
    endif

    if (Btype(To)==FreeSlipBuff)  call AttenuateTop(U2,V2,W2,Pr)
    if (Btype(Ea)==OutletBuff) then
      if (buoyancy==1) then
        call AttenuateOut(U2,V2,W2,Pr,temperature)
      else
        call AttenuateOut(U2,V2,W2,Pr)
      endif
    endif

    if (masssourc==1) then
        call MASS_SOURC(Q,U2,V2,W2)
    endif


    call BoundU(1,U2)
    call BoundU(2,V2)
    call BoundU(3,W2)



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

      if (scalsourcetype==pointsource) then
       do i=1,computescalars
        Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i)=Scalar_2(scalsrci(i),scalsrcj(i),scalsrck(i),i)+&
         percdistrib(i)*(rho(l)+beta(l))*dt*totalscalsource/(dxPr(scalsrci(i))*dyPr(scalsrcj(i))*dzPr(scalsrck(i)))
       enddo
      endif

      if (scalsourcetype==volumesource) then
       do i=1,computescalars
        call VolScalSource(Scalar_2,rho(l)+beta(l))
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

    if (Btype(Ea)==OutletBuff) then
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

        call BoundU_GPU(1,Unx,Uny,Unz,Prny,Prnz,&
                             Btype,sideU,&
                             Uin,U,0)
        call BoundU_GPU(1,Vnx,Vny,Vnz,Prny,Prnz,&
                             Btype,sideU,&
                             Vin,V,0)
        call BoundU_GPU(1,Wnx,Wny,Wnz,Prny,Prnz,&
                             Btype,sideU,&
                             Win,W,0)

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
   real(KND) Ap,Apre,Aprn,Aprt,S

   integer i,j,k,it,x,y,z

   type(TIBPoint),pointer:: IBP


   write(*,*) "Otherterms:"

   !$hmpp <tsteps> advancedload, args[PressureGrad::dxU,PressureGrad::dyV,PressureGrad::dzW]
   !$hmpp <tsteps> advancedload, args[PressureGrad::prgradientx,PressureGrad::prgradienty]
   !$hmpp <tsteps> advancedload, args[PressureGrad::Pr,PressureGrad::U,PressureGrad::V,PressureGrad::W,PressureGrad::dt]
   !$hmpp <tsteps> advancedload, args[PressureGrad::Btype]

   !$hmpp <tsteps> advancedload, args[PressureGrad::coef]

   !Pressure gradient terms
   !$hmpp <tsteps> PressureGrad callsite, args[*].noupdate=true
   call PressureGrad(Prnx,Prny,Prnz,&
                     Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                     dxU,dyV,dzW,&
                     Btype,prgradientx,prgradienty,&
                     Pr,U2,V2,W2,&
                     dt,coef)

   !$hmpp <tsteps> delegatedstore, args[PressureGrad::U,PressureGrad::V,PressureGrad::W]
write (*,*) "sum U2",sum(U2(1:Unx,1:Uny,1:Unz))
write (*,*) "sum V2",sum(V2(1:Vnx,1:Vny,1:Vnz))
write (*,*) "sum W2",sum(W2(1:Wnx,1:Wny,1:Wnz))

   Re_gt_0: if (Re>0) then

     !Diffusion using Crank Nicolson
     !first approximation using forward Euler
     !iteration SOR or Gauss-Seidel


     !$hmpp <tsteps> advancedload, args[ForwEul::Visc]
     !$hmpp <tsteps> advancedload, args[ForwEul::dxPr,ForwEul::dyPr,ForwEul::dzPr]

     !$hmpp <tsteps> ForwEul callsite, args[*].noupdate=true
     call ForwEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                  dxPr,dyPr,dzPr,dxU,dyV,dzW,&
                  U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                  dt,coef)

   !$hmpp <tsteps> delegatedstore, args[ForwEul::U3,ForwEul::V3,ForwEul::W3]
write (*,*) "sum U3",sum(U3(1:Unx,1:Uny,1:Unz))
write (*,*) "sum V3",sum(V3(1:Vnx,1:Vny,1:Vnz))
write (*,*) "sum W3",sum(W3(1:Wnx,1:Wny,1:Wnz))

!      call MomSourc(U3,V3,W3)


!      do i=1,NIBPointsU
!       xi=IBPijkU(1,i)
!       yj=IBPijkU(2,i)
!       zk=IBPijkU(3,i)
!       U3(xi,yj,zk)=U3(xi,yj,zk)+IBPsourcU*dt
!       U2(xi,yj,zk)=U2(xi,yj,zk)+IBPsourcU*dt
!      enddo
!
!      do i=1,NIBPointsV
!       xi=IBPijkV(1,i)
!       yj=IBPijkV(2,i)
!       zk=IBPijkV(3,i)
!       V3(xi,yj,zk)=V3(xi,yj,zk)+IBPsourcV*dt
!       V2(xi,yj,zk)=V2(xi,yj,zk)+IBPsourcV*dt
!      enddo
!
!      do i=1,NIBPointsW
!       xi=IBPijkW(1,i)
!       yj=IBPijkW(2,i)
!       zk=IBPijkW(3,i)
!       W3(xi,yj,zk)=W3(xi,yj,zk)+IBPsourcW*dt
!       W2(xi,yj,zk)=W2(xi,yj,zk)+IBPsourcW*dt
!      enddo

!      if (associated(FirstIBPoint)) then   !Immersed boundary terms, in the future should be in an array
!       IBP => FirstIBPoint
!       do
!        x=IBP%x
!        y=IBP%y
!        z=IBP%z
!         if (IBP%component==1) then
!          U3(x,y,z)=U3(x,y,z)+IBP%MSourc*dt
!          U2(x,y,z)=U2(x,y,z)+IBP%MSourc*dt
!         elseif (IBP%component==2) then
!          V3(x,y,z)=V3(x,y,z)+IBP%MSourc*dt
!          V2(x,y,z)=V2(x,y,z)+IBP%MSourc*dt
!         elseif (IBP%component==3) then
!          W3(x,y,z)=W3(x,y,z)+IBP%MSourc*dt
!          W2(x,y,z)=W2(x,y,z)+IBP%MSourc*dt
!         endif
!        if (associated(IBP%next)) then
!         IBP=>IBP%next
!        else
!         exit
!        endif
!       enddo
!      endif

     if (gridtype==UNIFORMGRID.and.GPU>0) then                  !Performs the diffusion terms
      write (*,*) "GPU CN call"

      !$hmpp <tsteps> advancedload, args[UnifRedBlack::sideU,UnifRedBlack::Uin,UnifRedBlack::Vin,UnifRedBlack::Win,&
      !$hmpp <tsteps>  UnifRedBlack::maxCNiter,UnifRedBlack::epsCN]


      !$hmpp <tsteps> UNIFREDBLACK callsite, args[*].noupdate=true
      call UNIFREDBLACK_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,&
                            Btype,sideU,&
                            dt,dxmin,dymin,dzmin,&
                            Uin,Vin,Win,&
                            U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                            coef,maxCNiter,epsCN,it,S)
      write(*,*) "back from GPU CN", it,S

     !$hmpp <tsteps> delegatedstore, args[UnifRedBlack::U3,UnifRedBlack::V3,UnifRedBlack::W3,&
     !$hmpp <tsteps>   UnifRedBlack::iters,UnifRedBlack::residuum]

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


  !$hmpp <tsteps> PressureGrad codelet
  subroutine PressureGrad(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxU,dyV,dzW,&
                          Btype,prgradientx,prgradienty,&
                          Pr,U,V,W,&
                          dt,coef)
  implicit none
#ifdef __HMPP
  integer,parameter :: KND=4
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
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
          U(i,j,k)=U(i,j,k)+A*(Pr(i+1,j,k)-Pr(i,j,k))/dxU(i)+A*prgradientx
     enddo
    enddo
   enddo
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
          V(i,j,k)=V(i,j,k)+A*(Pr(i,j+1,k)-Pr(i,j,k))/dyV(j)+A*prgradienty
     enddo
    enddo
   enddo
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
          W(i,j,k)=W(i,j,k)+A*(Pr(i,j,k+1)-Pr(i,j,k))/dzW(k)
     enddo
    enddo
   enddo
  end subroutine PressureGrad



  !$hmpp <tsteps> ForwEul codelet
  subroutine ForwEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                    dxPr,dyPr,dzPr,dxU,dyV,dzW,&
                    U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                    dt,coef)
  implicit none
#ifdef __HMPP
  integer,parameter :: KND=4
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
  end subroutine ForwEul


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
        call BoundU(1,U3)
        call BoundU(2,V3)
        call BoundU(3,W3)
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
        call BoundU(1,U3)
        call BoundU(2,V3)
        call BoundU(3,W3)
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


  subroutine SubgridStresses(U,V,W,Pr)
  real(KND),intent(in):: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)

     if (sgstype==SmagorinskyModel) then
                       call Smag(U,V,W)
     elseif (sgstype==VremanModel) then
                       if (GPU>0.and.gridtype==uniformgrid.and. Prnx*Prny*Prnz > 50) then
                           !$hmpp <tsteps> Vreman callsite, args[*].noupdate=true
     !$hmpp <tsteps> delegatedstore, args[Vreman::Visc]
                           call Vreman_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxmin,dymin,dzmin,dt,Re,U,V,W,Visc)
                       else
                           call Vreman(U,V,W)
                       endif
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
  end subroutine SubgridStresses


























 !GPU codelets

!$hmpp <tsteps> group, target=CUDA

  !$hmpp <tsteps> mapbyname, Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz
  !$hmpp <tsteps> mapbyname, Btype,sideU,Re
  !$hmpp <tsteps> mapbyname, dt
  !$hmpp <tsteps> mapbyname, dxPr,dyPr,dzPr,dxU,dyV,dzW
  !$hmpp <tsteps> mapbyname, Uin,Vin,Win,Pr,Visc
  !$hmpp <tsteps> map, args[Vreman::dx,CDS::dx,UnifRedBlack::dxmin]
  !$hmpp <tsteps> map, args[Vreman::dy,CDS::dy,UnifRedBlack::dymin]
  !$hmpp <tsteps> map, args[Vreman::dz,CDS::dz,UnifRedBlack::dzmin]

  !$hmpp <tsteps> map, args[PressureGrad::coef,ForwEul::coef,UnifRedBlack::coef]

 !U,V,W
 !$hmpp <tsteps> map, args[Vreman::U,CDS::U,ForwEul::U,UnifRedBlack::U]
 !$hmpp <tsteps> map, args[Vreman::V,CDS::V,ForwEul::V,UnifRedBlack::V]
 !$hmpp <tsteps> map, args[Vreman::W,CDS::W,ForwEul::W,UnifRedBlack::W]

 !U2,V2,W2
 !$hmpp <tsteps> map, args[PressureGrad::U,ForwEul::U2,UnifRedBlack::U2]
 !$hmpp <tsteps> map, args[PressureGrad::V,ForwEul::V2,UnifRedBlack::V2]
 !$hmpp <tsteps> map, args[PressureGrad::W,ForwEul::W2,UnifRedBlack::W2]

 !U3,V3,W3 on GPU device mapped also to Ustar,Vstar,Wstar
 !$hmpp <tsteps> map, args[CDS::U2,ForwEul::U3,UnifRedBlack::U3]
 !$hmpp <tsteps> map, args[CDS::V2,ForwEul::V3,UnifRedBlack::V3]
 !$hmpp <tsteps> map, args[CDS::W2,ForwEul::W3,UnifRedBlack::W3]


#include "boundaries_GPU.f90"
#include "cds_GPU.f90"
#include "smagorinsky_GPU.f90"
#include "tsteps_GPU.f90"




end module TSTEPS

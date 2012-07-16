module Wallmodels
use PARAMETERS
! use SCALARS


implicit none

 private
 public WMPoint, AddWMPoint, FirstWMPoint, ComputeViscsWM, MoveWMPointsToArray,&
        InitTempFl, GroundDeposition, GroundUstar

 type WMpoint   !points in which we apply wall model

  integer   :: xi
  integer   :: yj
  integer   :: zk

  real(KND) :: distx
  real(KND) :: disty
  real(KND) :: distz

  real(KND) :: z0=0
  real(KND) :: ustar=1
  real(KND) :: temp=0
  real(KND) :: tempfl=0
  real(KND) :: wallu=0
  real(KND) :: wallv=0
  real(KND) :: wallw=0

  real(KND),allocatable:: depscalar(:)

  type(WMpoint),pointer:: next=>null()

 endtype WMpoint

 type(WMpoint),pointer::FirstWMPoint=>null(), LastWMPoint=>null()

 type(WMPoint),dimension(:),allocatable :: WMPoints

 interface assignment(=)
   module procedure WMPtoWMP
 end interface

 contains


  subroutine WMPtoWMP(ToWMP,FromWMP)
    type(WMPoint),intent(out) :: ToWMP
    type(WMPoint),intent(in)  :: FromWMP

    ToWMP%xi = FromWMP%xi
    ToWMP%yj = FromWMP%yj
    ToWMP%zk = FromWMP%zk

    ToWMP%distx = FromWMP%distx
    ToWMP%disty = FromWMP%disty
    ToWMP%distz = FromWMP%distz

    ToWMP%z0     = FromWMP%z0
    ToWMP%ustar  = FromWMP%ustar
    ToWMP%temp   = FromWMP%temp
    ToWMP%tempfl = FromWMP%tempfl
    ToWMP%wallu  = FromWMP%wallu
    ToWMP%wallv  = FromWMP%wallv
    ToWMP%wallw  = FromWMP%wallw

    allocate(ToWMP%depscalar(size(FromWMP%depscalar)))

    ToWMP%depscalar = FromWMP%depscalar

  end subroutine WMPtoWMP


  subroutine AddWMpoint(WMP)
  type(WMpoint),intent(in):: WMP

   if (.not.associated(LastWMPoint)) then

    allocate(FirstWMPoint)

    if (computedeposition>0) then
      allocate( FirstWMPoint%depscalar(size(WMP%depscalar)) )
      FirstWMPoint%depscalar = 0
    endif

    FirstWMPoint=WMP
    LastWMPoint=>FirstWMPoint

   else

    allocate(LastWMPoint%next)

    if (computedeposition>0) then
      allocate( LastWMPoint%next%depscalar(size(WMP%depscalar)) )
      LastWMPoint%next%depscalar = 0
    endif

    LastWMPoint%next=WMP
    LastWMPoint=>LastWMPoint%next

   endif
  endsubroutine AddWMPoint


  subroutine WMPoint_DeallocateList(WMP)
    type(WMPoint),pointer,intent(inout) :: WMP
    type(WMPoint),pointer :: Aux,Aux2

    if (.not.associated(WMP)) return

    Aux => WMP

    do
      if (allocated(Aux%depscalar)) deallocate(Aux%depscalar)
      if (associated(Aux%next)) then
        Aux => Aux%next
      else
        exit
      endif
    enddo


    Aux =>WMP

    do

      if (associated(Aux%next)) then
        Aux2 => Aux%next
      else
        Aux2 => null()
      endif

      deallocate(Aux)

      if (associated(Aux2)) then
        Aux => Aux2
      else
        exit
      endif

    enddo

    WMP => null()

  end subroutine WMPoint_DeallocateList



  subroutine MoveWMPointsToArray
    type(WMPoint),pointer :: CurrentWMPoint
    integer :: i,nWMPoints

    NWMPoints = 0


    CurrentWMPoint => FirstWMPoint
    do
     if (associated(CurrentWMPoint)) then
       NWMPoints = NWMPoints + 1
     else
       exit
     endif
     CurrentWMPoint => CurrentWMPoint%next
    enddo


    allocate(WMPoints(NWMPoints))


    CurrentWMPoint => FirstWMPoint
    i = 0
    do
     if (associated(CurrentWMPoint)) then
       i = i + 1
       WMPoints(i) = CurrentWMPoint
     else
       exit
     endif
     CurrentWMPoint => CurrentWMPoint%next
    enddo


    call WMPoint_DeallocateList(FirstWMPoint)

  end subroutine MoveWMPointsToArray















  real(KND) function WM1ustar(vel,dist,ustar0,dp,dptrans)
   real(KND),parameter :: eps=1e-4_KND
   real(KND) :: yplcrit=11.225_KND
   real(KND),intent(in) :: vel,dist,ustar0,dp,dptrans
   real(KND) :: ustar,ustar2,kprime
   integer i


   if (wallmodeltype==1) then

    if ((dist*ustar0*Re)<yplcrit) then
      ustar=sqrt(vel/(dist*Re))
    else
      ustar=vel/(log(abs(ustar0*dist*Re))/0.41_KND+5.2_KND)
    endif

    i=1

    if ((dist*ustar*Re)>yplcrit) then
     ustar=ustar0
     do
     i=i+1
      ustar2=ustar
      ustar=vel/(log(abs(ustar2*dist*Re))/0.41_KND+5.2_KND)
      if  (abs(ustar-ustar2)/abs(ustar)<eps) exit
      if (i>=50) then
                  ustar=0
                  exit
      endif
     enddo
    endif

   else
    ustar=ustar0
    kprime=0.41_KND*(1-dptrans*(1./(0.41_KND*ustar0**2)-1./(2*ustar0**2)))
    yplcrit=-LambertW(-kprime*0.119_KND)/kprime

    if ((dist*ustar0*Re)<yplcrit) then
      ustar=sqrt(vel/(dist*Re)-dp/2)
    else
      ustar=vel*(1-Re*dp/0.41_KND)/(log(abs(ustar0*dist*Re))/0.41_KND+5.2_KND)
    endif

    i=1

    if ((dist*ustar*Re)>yplcrit) then
     do
     i=i+1
      ustar2=ustar
      ustar=vel*(1._KND-Re*dp/0.41_KND)/(log(abs(ustar2*dist*Re))/0.41_KND+5.2_KND)
      if  (abs(ustar-ustar2)/abs(ustar)<eps) exit
      if (i>=50) then
                  ustar=0
                  exit
      endif
     enddo

    endif

   endif

   WM1ustar=ustar
  endfunction WM1ustar

  real(KND) function WM1Visc(WMP,U,V,W,Pr)
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND),dimension(1:,1:,1:),intent(in):: Pr
   integer i,j,k
   real(KND) ustar,vel,dist,dp,dptrans,dpx,dpy,dpz
   type(WMPoint):: WMP

    i=WMP%xi
    j=WMP%yj
    k=WMP%zk

    dist=sqrt(WMP%distx**2+WMP%disty**2+WMP%distz**2)

    vel=0
    ustar=0
    dp=0
    dptrans=0

    if (abs(WMP%disty)/dymin<1.e-2.and.abs(WMP%distz)/dzmin<1.e-2) vel=vel+((U(i,j,k)+U(i-1,j,k))/2._KND-WMP%wallu)**2
    if (abs(WMP%distx)/dxmin<1.e-2.and.abs(WMP%distz)/dzmin<1.e-2) vel=vel+((V(i,j,k)+V(i,j-1,k))/2._KND-WMP%wallv)**2
    if (abs(WMP%disty)/dymin<1.e-2.and.abs(WMP%distx)/dxmin<1.e-2) vel=vel+((W(i,j,k)+W(i,j,k-1))/2._KND-WMP%wallw)**2

    vel=sqrt(vel)

    if (vel/=0) then
      if (wallmodeltype>1) then
       dpx=0
       if (i>1) dpx=dpx+(Pr(i,j,k)-Pr(i-1,j,k))/(xPr(i)-xPr(i-1))
       if (i<Unx) dpx=dpx+(Pr(i+1,j,k)-Pr(i,j,k))/(xPr(i+1)-xPr(i))
       if (i==1.and.Btype(We)==PERIODIC) dpx=dpx+(Pr(1,j,k)-Pr(Prnx,j,k))/(xPr(1)-xPr(0))
       dpx=dpx/2

       dpy=0
       if (j>1) dpy=dpy+(Pr(i,j,k)-Pr(i,j-1,k))/(yPr(j)-yPr(j-1))
       if (j<Vny) dpy=dpy+(Pr(i,j+1,k)-Pr(i,j,k))/(yPr(j+1)-yPr(j))
       if (j==1.and.Btype(So)==PERIODIC) dpy=dpy+(Pr(i,1,k)-Pr(i,Prny,k))/(yPr(1)-yPr(0))
       dpy=dpy/2

       dpz=0
       if (k>1) dpz=dpz+(Pr(i,j,k)-Pr(i,j,k-1))/(zPr(k)-zPr(k-1))
       if (k<Wnz) dpz=dpz+(Pr(i,j,k+1)-Pr(i,j,k))/(zPr(k+1)-zPr(k))
       if (k==1.and.Btype(Bo)==PERIODIC) dpz=dpz+(Pr(i,j,1)-Pr(i,j,Prnz))/(zPr(1)-zPr(0))
       dpz=dpz/2

       if (abs(WMP%distx)*1.1>dist)&
          dp=dp+((U(i,j,k)+U(i-1,j,k))/(2._KND*vel))*dpx
       if (abs(WMP%disty)*1.1<dist)&
          dp=dp+((V(i,j,k)+V(i,j-1,k))/(2._KND*vel))*dpy
       if (abs(WMP%distz)*1.1<dist)&
          dp=dp+((W(i,j,k)+W(i,j,k-1))/(2._KND*vel))*dpz
       dptrans=dptrans+WMP%distx/dist*dpx
       dptrans=dptrans+WMP%disty/dist*dpy
       dptrans=dptrans+WMP%distz/dist*dpz
      else
       dp=0
       dptrans=0
      endif
      ustar=WMP%ustar
      ustar=WM1ustar(vel,dist,ustar,dp,dptrans)
      WMP%ustar=ustar
    endif

    if (ustar<0) ustar=0

    if (vel>0) then
     if (dist*ustar*Re>1) then
      WM1Visc=ustar*ustar*dist/vel
     elseif (Re>0) then
      WM1Visc=1._KND/Re
     else
      WM1Visc=0
     endif
    elseif (Re>0) then
      WM1Visc=1._KND/Re
    else
      WM1Visc=0
    endif
  endfunction WM1Visc




  real(KND) function WM2ustar(vel,dist,ustar0,dp,dptrans,z0)
   real(KND),parameter:: eps=1e-4_KND
   real(KND):: yplcrit=11.225_KND
   real(KND),intent(in):: vel,dist,ustar0,z0,dp,dptrans
   real(KND) :: kprime

   if (wallmodeltype==1) then

    if (dist<=z0) then
     if (Re>0) then
      if ((dist*ustar0*Re)<yplcrit) then
        WM2ustar=sqrt(vel/(dist*Re))
      else
        WM2ustar=vel/(log(abs(ustar0*dist*Re))/0.41_KND+5.2_KND)
      endif
     else
      stop "The wall model need positive viscosity under roughness length."
     endif
    else
      WM2ustar=vel*0.41_KND/log(dist/z0)
    endif

   else

    kprime=0.41_KND*(1-dptrans*(1./(0.41*ustar0**2)-1./(2*ustar0**2)))
    yplcrit=-LambertW(-kprime*0.119_KND)/kprime

    if (dist<=z0) then
     if (Re>0) then
      if ((dist*ustar0*Re)<yplcrit) then
        WM2ustar=sqrt(vel/(dist*Re)-dp/2)
      else
        WM2ustar=vel*(1-Re*dp/0.41_KND)/(log(abs(ustar0*dist*Re))/0.41_KND+5.2_KND)
      endif
     else
      stop "The wall model need positive viscosity under roughness length."
     endif
    else
      WM2ustar=vel*(1-Re*dp/0.41_KND)*0.41_KND/log(dist/z0)
    endif

   endif
  endfunction WM2ustar


  pure real(KND) function PsiM_MO(zeta)
   real(KND),intent(in):: zeta
   real(KND) x

   if (zeta<0) then
    x=(1-15._KND*zeta)**(1/4._KND)
    PsiM_MO=log(((1+x**2)/2._KND)*((1+x)/2._KND)**2)-2._KND*atan(x)+pi/2
   else
    PsiM_MO=-4.8_KND*zeta !GABLS recommendation
   endif
  endfunction PsiM_MO


  pure real(KND) function PsiH_MO(zeta)
   real(KND),intent(in):: zeta
   real(KND) x

   if (zeta<0) then
    x=(1-15._KND*zeta)**(1/4._KND)
    PsiH_MO=2._KND*log((1+x**2)/2._KND)
   else
    PsiH_MO=-7.8_KND*zeta !GABLS recommendation
   endif
  endfunction PsiH_MO


 pure real(KND) function Obukhov_zL(ustar,tempfl,tempref,g,z)
   real(KND),intent(in):: ustar,tempfl,tempref,g,z

   Obukhov_zL=z*(0.4_KND*(g/tempref)*tempfl)/(-ustar**3)
  endfunction Obukhov_zL

  real(KND) function WM_MO_FLUX_ustar(vel,dist,ustar0,z0,tempflux)
   real(KND),parameter:: eps=1e-3
   real(KND),parameter:: yplcrit=11.225_KND
   real(KND),intent(in):: vel,dist,ustar0,z0,tempflux
   real(KND) ustar,ustar2,zL,Psi
   integer i

   if (dist<=z0) then

    if (Re>0) then

     if ((dist*ustar0*Re)<yplcrit) then
       ustar=sqrt(vel/(dist*Re))
     else
       ustar=vel/(log(abs(ustar0*dist*Re))/0.4_KND+5.2_KND)
     endif

    else

     stop "ERROR: The wall model needs positive viscosity under roughness length."

    endif

   else

    ustar=ustar0
    i=0
    zL=0
    Psi=0

    do
     i=i+1
     ustar2=ustar

     if (((dist/z0)-Psi)<1E-5) then
       write(*,*) "i", i
       write(*,*) vel,dist,ustar0,z0,tempflux
       write(*,*) Psi, "=", dist, "/", z0
       STOP
     endif

     if (log((dist/z0)-Psi)<1E-5) then
       write(*,*) "i", i
       write(*,*) vel,dist,ustar0,z0,tempflux
       write(*,*) Psi, "=", dist, "/", z0
       STOP
     endif

     ustar=ustar+(max(vel*0.4_KND/(log(max((dist/z0)-Psi,1E-5))),0._KND)-ustar)/2

     if (ustar<1E-4) then
      zL=-10000
     else
      zL=zL+(Obukhov_zL(ustar,tempflux,temperature_ref,grav_acc,dist)-zL)/2
     endif

     Psi=PsiM_MO(zL)

     if  (abs(ustar-ustar2)/max(abs(ustar),1.e-3_KND)<eps) exit

     if (i>=50) then
                 ustar=0
                 exit
     endif

    enddo

   endif

   WM_MO_FLUX_ustar=ustar
  endfunction WM_MO_FLUX_ustar


  subroutine WM_MO_DIRICHLET_ustar_tfl(vel,dist,z0,ustar,tempflux,tempdif)
   real(KND),parameter:: eps=1e-3
   real(KND):: yplcrit=11.225
   real(KND) vel,dist,ustar,z0,tempflux,tempdif,ustar0,tempflux0,zL,zL0,Rib
   integer i

   ustar0=WM2ustar(vel,dist,ustar0,0._KND,0._KND,z0)
   tempflux0=tempflux

   if (dist<=z0) then

     if (Re>0) then
      if ((dist*ustar0*Re)<yplcrit) then
        ustar=sqrt(vel/(dist*Re))
      else
        ustar=vel/(log(abs(ustar0*dist*Re))/0.4_KND+5.2_KND)
      endif
     else
      stop "The wall model needs positive viscosity under roughness length."
     endif

    else

     Rib=-grav_acc*dist*tempdif/(temperature_ref*max(vel**2,1E-6_KND))
     zL=0
     i=0
     if (Rib>0.34_KND) then
                         ustar=0
                         tempflux=0
                         return
     endif
     do
       i=i+1
       zL0=zL
       zL=Rib*(log(dist/z0)-PsiM_MO(zl))**2/(log(dist/z0)-PsiH_MO(zl))
       if  (abs(zL-zL0)/max(abs(zL),1.e-3_KND)<eps) exit
       if (i>=50.or.zL>100) then
                   ustar=0
                   tempflux=0
                   return
       endif
     enddo
     ustar=vel*0.4_KND/(log(dist/z0)-PsiM_MO(zL))
     tempflux=0.4_KND*ustar*tempdif/(log(dist/z0)-PsiH_MO(zL))
   endif
  endsubroutine WM_MO_DIRICHLET_ustar_tfl









  real(KND) function WM2Visc(WMP,U,V,W,Pr)
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND),dimension(1:,1:,1:),intent(in):: Pr
   integer i,j,k
   real(KND) ustar,vel,dist,z0,dp,dptrans,dpx,dpy,dpz
   type(WMPoint):: WMP

   i=WMP%xi
   j=WMP%yj
   k=WMP%zk

   dist=sqrt(WMP%distx**2+WMP%disty**2+WMP%distz**2)

   vel=0
   dp=0
   dptrans=0
   if (abs(WMP%distx)/dymin<0.9_KND*dist) vel=vel+((U(i,j,k)+U(i-1,j,k))/2._KND-WMP%wallu)**2
   if (abs(WMP%disty)/dymin<0.9_KND*dist) vel=vel+((V(i,j,k)+V(i,j-1,k))/2._KND-WMP%wallv)**2
   if (abs(WMP%distz)/dymin<0.9_KND*dist) vel=vel+((W(i,j,k)+W(i,j,k-1))/2._KND-WMP%wallw)**2

   vel=sqrt(vel)

   if (vel/=0) then
     if (wallmodeltype>1) then
      dpx=0
      if (i>1) dpx=dpx+(Pr(i,j,k)-Pr(i-1,j,k))/(xPr(i)-xPr(i-1))
      if (i<Unx) dpx=dpx+(Pr(i+1,j,k)-Pr(i,j,k))/(xPr(i+1)-xPr(i))
      if (i==1.and.Btype(We)==PERIODIC) dpx=dpx+(Pr(1,j,k)-Pr(Prnx,j,k))/(xPr(1)-xPr(0))
      dpx=dpx/2

      dpy=0
      if (j>1) dpy=dpy+(Pr(i,j,k)-Pr(i,j-1,k))/(yPr(j)-yPr(j-1))
      if (j<Vny) dpy=dpy+(Pr(i,j+1,k)-Pr(i,j,k))/(yPr(j+1)-yPr(j))
      if (j==1.and.Btype(So)==PERIODIC) dpy=dpy+(Pr(i,1,k)-Pr(i,Prny,k))/(yPr(1)-yPr(0))
      dpy=dpy/2

      dpz=0
      if (k>1) dpz=dpz+(Pr(i,j,k)-Pr(i,j,k-1))/(zPr(k)-zPr(k-1))
      if (k<Wnz) dpz=dpz+(Pr(i,j,k+1)-Pr(i,j,k))/(zPr(k+1)-zPr(k))
      if (k==1.and.Btype(Bo)==PERIODIC) dpz=dpz+(Pr(i,j,1)-Pr(i,j,Prnz))/(zPr(1)-zPr(0))
      dpz=dpz/2

      if (abs(WMP%distx)*1.1>dist)&
         dp=dp+((U(i,j,k)+U(i-1,j,k))/(2._KND*vel))*dpx
      if (abs(WMP%disty)*1.1<dist)&
         dp=dp+((V(i,j,k)+V(i,j-1,k))/(2._KND*vel))*dpy
      if (abs(WMP%distz)*1.1<dist)&
         dp=dp+((W(i,j,k)+W(i,j,k-1))/(2._KND*vel))*dpz
      dptrans=dptrans+WMP%distx/dist*dpx
      dptrans=dptrans+WMP%disty/dist*dpy
      dptrans=dptrans+WMP%distz/dist*dpz
     else
      dp=0
      dptrans=0
     endif

     ustar=WMP%ustar
     z0=WMP%z0
     ustar=WM2ustar(vel,dist,ustar,dp,dptrans,z0)
     if (ustar<0) ustar=0
     WMP%ustar=ustar
   endif

   if (vel>0.and.ustar*ustar*dist/vel>1._KND/Re) then
     WM2Visc=ustar*ustar*dist/vel
   elseif (Re>0) then
     WM2Visc=1._KND/Re
   else
     WM2Visc=0
   endif

  endfunction WM2Visc




  real(KND) function WM_MO_FLUX(WMP,U,V,W,Pr)
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND),dimension(1:,1:,1:),intent(in):: Pr
   integer i,j,k
   real(KND) ustar,vel,dist,z0
   type(WMPoint):: WMP

   i=WMP%xi
   j=WMP%yj
   k=WMP%zk

   dist=sqrt(WMP%distx**2+WMP%disty**2+WMP%distz**2)

   vel=0
   if (abs(WMP%distx)/dymin<0.9_KND*dist) vel=vel+((U(i,j,k)+U(i-1,j,k))/2._KND-WMP%wallu)**2
   if (abs(WMP%disty)/dymin<0.9_KND*dist) vel=vel+((V(i,j,k)+V(i,j-1,k))/2._KND-WMP%wallv)**2
   if (abs(WMP%distz)/dymin<0.9_KND*dist) vel=vel+((W(i,j,k)+W(i,j,k-1))/2._KND-WMP%wallw)**2

   vel=sqrt(vel)

   if (vel/=0) then
     ustar=WMP%ustar
     z0=WMP%z0
     ustar=WM_MO_FLUX_ustar(vel,dist,ustar,z0,WMP%tempfl)
     if (ustar<0) ustar=0
     WMP%ustar=ustar
   endif

   if (vel>0.and.ustar*ustar*dist/vel>1._KND/Re) then
     WM_MO_FLUX=ustar*ustar*dist/vel
   elseif (Re>0) then
     WM_MO_FLUX=1._KND/Re
   else
     WM_MO_FLUX=0
   endif

  endfunction WM_MO_FLUX



  subroutine WM_MO_DIRICHLET(visc,WMP,U,V,W,Pr)
   real(KND) visc
   real(KND),dimension(-2:,-2:,-2:),intent(in):: U,V,W
   real(KND),dimension(1:,1:,1:),intent(in):: Pr
   integer i,j,k
   real(KND) ustar,vel,dist,z0,tempflux
   type(WMPoint):: WMP

   i=WMP%xi
   j=WMP%yj
   k=WMP%zk
   !provizorni zpusob pro steny orientovane se stenami site
   dist=sqrt(WMP%distx**2+WMP%disty**2+WMP%distz**2)

   vel=0
   if (abs(WMP%distx)/dymin<0.9_KND*dist) vel=vel+((U(i,j,k)+U(i-1,j,k))/2._KND-WMP%wallu)**2
   if (abs(WMP%disty)/dymin<0.9_KND*dist) vel=vel+((V(i,j,k)+V(i,j-1,k))/2._KND-WMP%wallv)**2
   if (abs(WMP%distz)/dymin<0.9_KND*dist) vel=vel+((W(i,j,k)+W(i,j,k-1))/2._KND-WMP%wallw)**2
   vel=sqrt(vel)

   ustar=WMP%ustar
   tempflux=WMP%tempfl

   z0=WMP%z0

   call WM_MO_DIRICHLET_ustar_tfl(vel,dist,z0,ustar,tempflux,WMP%temp-temperature(i,j,k))

   if (ustar<0) ustar=0


   WMP%ustar=ustar
   WMP%tempfl=tempflux

   if ((vel>0.or.tempflux/=0)) then
     visc=MAX(ustar*ustar*dist/vel,1._KND/Re)
   elseif (Re>0) then
     visc=1._KND/Re
   else
     visc=0
   endif

   if (k==1) BsideTFLArr(i,j)=tempflux
  endsubroutine WM_MO_DIRICHLET



  pure subroutine BOUND_tempfl(Nu)
   real(KND),intent(inout):: Nu(-1:,-1:)
   integer i,j,nx,ny

   nx=Prnx
   ny=Prny

   if (Btype(Ea)==PERIODIC) then
     do j=1,ny
       Nu(0,j)=Nu(nx,j)
       Nu(nx+1,j)=Nu(1,j)
     enddo
   else
     do j=1,ny
       Nu(0,j)=Nu(1,j)
       Nu(nx+1,j)=Nu(nx,j)
     enddo
   endif

   if (Btype(No)==PERIODIC) then
     do i=1,nx
       Nu(i,0)=Nu(i,ny)
       Nu(i,ny+1)=Nu(i,1)
     enddo
   else
     do i=1,nx
       Nu(i,0)=Nu(i,1)
       Nu(i,ny+1)=Nu(i,ny)
     enddo
   endif
  endsubroutine BOUND_tempfl





  subroutine InitTempFL
    integer i

    if (buoyancy==1.and.TBtype(Bo)==DIRICHLET) then
      do i=1,size(WMPoints)
        if (WMPoints(i)%zk==1) WMPoints(i)%tempfl = -TDiff(WMPoints(i)%xi,WMPoints(i)%yj,1)*&
                       (temperature(WMPoints(i)%xi,WMPoints(i)%yj,1) - temperature(WMPoints(i)%xi,WMPoints(i)%yj,0))
      enddo
    endif

  endsubroutine InitTempFL


  subroutine ComputeViscsWM(U,V,W,Pr)
   real(KND),dimension(-2:,-2:,-2:):: U,V,W
   real(KND),dimension(1:,1:,1:):: Pr
   type(WMPoint),pointer:: WMP
   integer i,j

   !$omp parallel private(i,j)

   if (buoyancy==1.and.TBtype(Bo)==DIRICHLET) then
     !$omp do
     do j=1,Prny
       do i=1,Prnx
         BsideTArr(i,j)=SurfTemperature(xPr(i),yPr(j), time+dt/2._TIM)
       enddo
      enddo
     !$omp end do
   endif

   !$omp do
   do i = 1,size(WMPoints)

      if (WMPoints(i)%z0>0) then

        if (buoyancy==1 .and. TBtype(Bo)==CONSTFLUX) then

           Visc(WMPoints(i)%xi,WMPoints(i)%yj,WMPoints(i)%zk) = WM_MO_FLUX(WMPoints(i),U,V,W,Pr)

        else if (buoyancy==1 .and. TBtype(Bo)==DIRICHLET) then

          if (WMPoints(i)%zk==1) WMPoints(i)%temp = BsideTArr(WMPoints(i)%xi,WMPoints(i)%yj)

          call WM_MO_DIRICHLET(Visc(WMPoints(i)%xi,WMPoints(i)%yj,WMPoints(i)%zk),WMPoints(i),U,V,W,Pr)

        else

           Visc(WMPoints(i)%xi,WMPoints(i)%yj,WMPoints(i)%zk) = WM2Visc(WMPoints(i),U,V,W,Pr)

        endif

      else

        if (Re<=0) then
         stop "The wall model requires positive viscosity or roughness length."
        endif

        Visc(WMPoints(i)%xi,WMPoints(i)%yj,WMPoints(i)%zk) = WM1Visc(WMPoints(i),U,V,W,Pr)

      endif

   enddo
   !$omp end do
   !$omp end parallel

   if (buoyancy==1.and. TBtype(Bo)==DIRICHLET) call Bound_tempfl(BsideTFLArr)
  endsubroutine ComputeViscsWM

  pure real(KND) function GroundUstar()
    if (any(WMPoints%zk == 1)) then
      GroundUstar = sum(WMPoints%ustar, mask = (WMPoints%zk == 1)) / count(WMPoints%zk == 1)
    else
      GroundUstar = 0
    endif
  end function GroundUstar

  pure real(KND) function TotalUstar()
    if (size(WMPoints) > 0) then
      TotalUstar = sum(WMPoints%ustar) / size(WMPoints%zk)
    else
      TotalUstar = 0
    endif
  end function TotalUstar

  pure function GroundDeposition() result(depos)
    real(KND), dimension(:) :: depos(1:Prnx,1:Prny,computescalars)

    integer :: i, j

    depos = 0

    do j = 1, size(WMPoints)

      if (allocated(WMPoints(j)%depscalar)) then

        do i = 1, computescalars
          depos(WMPoints(j)%xi,WMPoints(j)%yj,i) = depos(WMPoints(j)%xi,WMPoints(j)%yj,i) + WMPoints(j)%depscalar(i)
        enddo

      endif

    enddo
  endfunction GroundDeposition



  pure real(KND) function SurfTemperature(x,y,t)
   real(KND),intent(in):: x,y
   real(TIM),intent(in):: t

   SurfTemperature=sideTemp(Bo)  ! Needs bet to allow time evolution somehow
  endfunction


  real(KND) function LambertW(x)
   real(KND) x

   LambertW=11
   write(*,*) "Warning, Lambert function not defined!"
   STOP
  endfunction LambertW

 endmodule Wallmodels

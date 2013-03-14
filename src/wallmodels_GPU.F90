  subroutine PsiM_MO(Psi,zeta)
   implicit none
#include "hmpp-include.f90"

   real(knd),parameter :: pi = 3.1415926535_knd
   real(knd),intent(out) :: Psi
   real(knd),intent(in):: zeta
   real(knd) x
   intrinsic log,atan

   if (zeta<0) then
    x=(1-15._knd*zeta)**(1/4._knd)
    Psi=log(((1+x**2)/2._knd)*((1+x)/2._knd)**2)-2._knd*atan(x)+pi/2
   else
    Psi=-4.8_knd*zeta !GABLS recommendation
   endif
  endsubroutine PsiM_MO


  subroutine PsiH_MO(Psi,zeta)
   implicit none
#include "hmpp-include.f90"

   real(knd),intent(out) :: Psi
   real(knd),intent(in):: zeta
   real(knd) x
   intrinsic log

   if (zeta<0) then
    x=(1-15._knd*zeta)**(1/4._knd)
    Psi=2._knd*log((1+x**2)/2._knd)
   else
    Psi=-7.8_knd*zeta !GABLS recommendation
   endif
  endsubroutine PsiH_MO


  subroutine Obukhov_zL(zL,ustar,tempfl,tempref,g,z)
   implicit none
#include "hmpp-include.f90"

   real(knd),intent(out) :: zL
   real(knd),intent(in):: ustar,tempfl,tempref,g,z

   zL=z*(0.4_knd*(g/tempref)*tempfl)/(-ustar**3)
  endsubroutine Obukhov_zL

  subroutine WM_MO_FLUX_ustar(vel,dist,ustar,z0,tempflux,Re,temperature_ref,grav_acc)
   implicit none
#include "hmpp-include.f90"

   real(knd),intent(inout) :: ustar
   real(knd),parameter  :: eps=1e-3
   real(knd),parameter  :: yplcrit=11.225_knd
   real(knd),intent(in) :: vel,dist,z0,tempflux
   real(knd),intent(in) :: Re,temperature_ref,grav_acc
   real(knd) ustar2,zL,zL2,Psi
   integer i

   if (dist<=z0) then

    if (Re>0) then

     if ((dist*ustar*Re)<yplcrit) then
       ustar=sqrt(vel/(dist*Re))
     else
       ustar=vel/(log(abs(ustar*dist*Re))/0.4_knd+5.2_knd)
     endif
    else
      ustar = 0
    endif

   else

    i=0
    zL=0
    Psi=0

    do
     i=i+1
     ustar2=ustar

     ustar=ustar+(max(vel*0.4_knd/(log(max((dist/z0)-Psi,1E-5))),0._knd)-ustar)/2

     if (ustar<1E-4) then
      zL=-10000
     else
      call Obukhov_zL(zL2,ustar,tempflux,temperature_ref,grav_acc,dist)
      zL=zL+(zL2-zL)/2
     endif

     call PsiM_MO(Psi,zL)

     if  (abs(ustar-ustar2)/max(abs(ustar),1.e-3_knd)<eps) exit

     if (i>=50) then
                 ustar=0
                 exit
     endif

    enddo

   endif

  endsubroutine WM_MO_FLUX_ustar



  subroutine WM_MO_FLUX_GPU(visc,i,j,k,distx,disty,distz,z0,tempfl,ustar,u,v,w,Re,temperature_ref,grav_acc)
   implicit none
#include "hmpp-include.f90"
   real(knd),intent(out) :: visc

   integer,intent(in)      :: i,j,k
   real(knd),intent(inout) :: ustar
   real(knd),intent(in)    :: distx,disty,distz,z0,tempfl
   real(knd),intent(in)    :: u,v,w
   real(knd),intent(in)    :: Re,temperature_ref,grav_acc
   real(knd) vel,dist


   dist=sqrt(distx**2+disty**2+distz**2)

   vel=0
   vel=vel+(u)**2 !(U(i,j,k)+U(i-1,j,k))/2._knd
   vel=vel+(v)**2 !(V(i,j,k)+V(i,j-1,k))/2._knd
!    if (abs(distz)/dymin<0.9_knd*dist) vel=vel+(w)**2 !(W(i,j,k)+W(i,j,k-1))/2._knd

   vel=sqrt(vel)

   if (vel/=0) then
     call WM_MO_FLUX_ustar(vel,dist,ustar,z0,tempfl,Re,temperature_ref,grav_acc)
     if (ustar<0) ustar=0
   endif

   if (vel>0.and.ustar*ustar*dist/vel>1._knd/Re) then
     visc=ustar*ustar*dist/vel
   elseif (Re>0) then
     visc=1._knd/Re
   else
     visc=0
   endif

  endsubroutine WM_MO_FLUX_GPU


  !$hmpp <tsteps> ComputeViscWM codelet
  subroutine ComputeViscsWM_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                                nWMPoints,WMPoints,&
                                TempBtype,Re,temperature_ref,grav_acc,&
                                U,V,W,Visc)
   implicit none
#include "hmpp-include.f90"

!    type hmppWMpoint   !points in which we apply wall model
!
!     integer   :: xi
!     integer   :: yj
!     integer   :: zk
!
!     real(knd) :: distx
!     real(knd) :: disty
!     real(knd) :: distz
!
!     real(knd) :: z0=0
!     real(knd) :: ustar=1
!     real(knd) :: temp=0
!     real(knd) :: tempfl=0
!
!    endtype hmppWMpoint

   integer,intent(in)    :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz
   integer,intent(in)              :: nWMPoints
   type(hmppWMPoint),intent(inout) :: WMPoints(nWMPoints)
   integer,intent(in)              :: TempBtype(6)
   real(knd),intent(in)            :: Re,temperature_ref,grav_acc
   real(knd),intent(in)            :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
   real(knd),intent(in)            :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
   real(knd),intent(in)            :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
   real(knd),intent(inout)         :: Visc(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
   integer i,xi,yj,zk
   real(knd) up,vp,wp

   !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(i),private(xi,yj,zk,distx,disty,distz,z0,tempfl,ustar,up,vp,wp)
   do i = 1,nWMPoints

        xi = WMPoints(i)%xi
        yj = WMPoints(i)%yj
        zk = WMPoints(i)%zk
        up = (U(xi,yj,zk)+U(xi-1,yj,zk))/2._knd
        vp = (V(xi,yj,zk)+V(xi,yj-1,zk))/2._knd
        wp = (W(xi,yj,zk)+W(xi,yj,zk-1))/2._knd



        if (TempBtype(Bo)==CONSTFLUX) then

           call WM_MO_FLUX_GPU(Visc(xi,yj,zk) ,&
                             xi, yj, zk,&
                             WMPoints(i)%distx, WMPoints(i)%disty, WMPoints(i)%distz,&
                             WMPoints(i)%z0, WMPoints(i)%tempfl, WMPoints(i)%ustar,&
                             up, vp, wp,&
                             Re, temperature_ref, grav_acc)

        endif


   enddo

  endsubroutine ComputeViscsWM_GPU

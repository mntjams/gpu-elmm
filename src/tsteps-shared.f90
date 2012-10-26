
 !$hmpp <tsteps> PressureGrad codelet
  subroutine PressureGrad(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxU,dyV,dzW,&
                          Btype,prgradientx,prgradienty,&
                          Pr,U,V,W,&
                          dt,coef)
  implicit none
#ifdef __HMPP
#include "hmpp-include.f90"
#endif
  integer,intent(in)      :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz
  real(KND),intent(in)    :: prgradientx,prgradienty,dt,coef
  integer,intent(in)      :: Btype(6)
#ifdef __HMPP
  real(KND),intent(in)    :: dxU(-2:Prnx+2),dyV(-2:Prny+2),dzW(-2:Prnz+2)
  real(KND),intent(in)    :: Pr(1:Unx+1,1:Vny+1,1:Wnz+1)
  real(KND),intent(inout) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
  real(KND),intent(inout) :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
  real(KND),intent(inout) :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
#else
  real(KND),intent(in),dimension(-2:) :: dxU,dyV,dzW
  real(KND),intent(inout) :: Pr(1:,1:,1:)
  real(KND),intent(inout),dimension(-2:,-2:,-2:) :: U,V,W
#endif
  real(KND) :: A
  integer i,j,k

   A=-coef*dt
   A=-coef*dt
   A=-coef*dt

   !$omp parallel
   !$omp do
   !$hmppcg grid blocksize myblocksize
   !$hmppcg permute (k,i,j)
   !$hmppcg gridify(k,i)
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
          U(i,j,k) = U(i,j,k)+A*(Pr(i+1,j,k)-Pr(i,j,k))/dxU(i)+A*prgradientx
     enddo
    enddo
   enddo
   !$omp enddo nowait
   !$omp do
   !$hmppcg grid blocksize myblocksize
   !$hmppcg permute (k,i,j)
   !$hmppcg gridify(k,i)
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
          V(i,j,k) = V(i,j,k)+A*(Pr(i,j+1,k)-Pr(i,j,k))/dyV(j)+A*prgradienty
     enddo
    enddo
   enddo
   !$omp enddo nowait
   !$omp do
   !$hmppcg grid blocksize myblocksize
   !$hmppcg permute (k,i,j)
   !$hmppcg gridify(k,i)
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
          W(i,j,k) = W(i,j,k)+A*(Pr(i,j,k+1)-Pr(i,j,k))/dzW(k)
     enddo
    enddo
   enddo
   !$omp enddo nowait
   !$omp end parallel
  end subroutine PressureGrad







 !$hmpp <tsteps> Convection codelet
 subroutine Convection(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,buoyancy,convmet,&
                       dxmin,dymin,dzmin,coriolisparam,grav_acc,temperature_ref,&
                       U,V,W,U2,V2,W2,Ustar,Vstar,Wstar,temperature,beta,rho,lev,dt)
  implicit none
#ifdef __HMPP
#include "hmpp-include.f90"
#endif

  integer,intent(in)   :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,buoyancy,convmet,lev
  real(KND),intent(in) :: dxmin,dymin,dzmin,coriolisparam,grav_acc,temperature_ref
  real(KND),intent(in) :: dt
#ifdef __HMPP
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
#else
  real(KND),dimension(-2:,-2:,-2:),intent(in)    :: U,V,W
  real(KND),dimension(-2:,-2:,-2:),intent(out)   :: U2,V2,W2
  real(KND),dimension(-2:,-2:,-2:),intent(inout) :: Ustar,Vstar,Wstar
  real(KND),dimension(-1:,-1:,-1:),intent(in) :: temperature
#endif
  real(KND),dimension(1:3),intent(in) :: beta,rho
  integer i,j,k
  intrinsic abs

      if (lev>1) then
        !$omp parallel private(i,j,k)
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Unz
         do j=1,Uny
          do i=1,Unx
           U2(i,j,k) = Ustar(i,j,k)*rho(lev)
          enddo
         enddo
        enddo
        !$omp end do nowait
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Vnz
         do j=1,Vny
          do i=1,Vnx
           V2(i,j,k) = Vstar(i,j,k)*rho(lev)
          enddo
         enddo
        enddo
        !$omp end do nowait
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Wnz
         do j=1,Wny
          do i=1,Wnx
           W2(i,j,k) = Wstar(i,j,k)*rho(lev)
          enddo
         enddo
        enddo
        !$omp end do
        !$omp end parallel
      else
        !$omp parallel private(i,j,k)
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Unz
         do j=1,Uny
          do i=1,Unx
           Ustar(i,j,k) = 0
          enddo
         enddo
        enddo
        !$omp end do nowait
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Vnz
         do j=1,Vny
          do i=1,Vnx
           Vstar(i,j,k) = 0
          enddo
         enddo
        enddo
        !$omp end do nowait
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Wnz
         do j=1,Wny
          do i=1,Wnx
           Wstar(i,j,k) = 0
          enddo
         enddo
        enddo
        !$omp end do nowait

        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Unz
         do j=1,Uny
          do i=1,Unx
           U2(i,j,k) = 0
          enddo
         enddo
        enddo
        !$omp end do nowait
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Vnz
         do j=1,Vny
          do i=1,Vnx
           V2(i,j,k) = 0
          enddo
         enddo
        enddo
        !$omp end do nowait
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Wnz
         do j=1,Wny
          do i=1,Wnx
           W2(i,j,k) = 0
          enddo
         enddo
        enddo
        !$omp end do
        !$omp end parallel
      endif

      if (convmet>0) then

#ifdef __HMPP
          call CDS_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxmin,dymin,dzmin,dt,Ustar,Vstar,Wstar,U,V,W)
#else
        if (convmet==2) then
          call CDU(Ustar,U,V,W)
          call CDV(Vstar,U,V,W)
          call CDW(Wstar,U,V,W)
        else if (convmet==4) then
          call CDS4U(Ustar,U,V,W)
          call CDS4V(Vstar,U,V,W)
          call CDS4W(Wstar,U,V,W)
        end if
        write(*,*) "U,V,Wstar",Ustar(Unx/2,Uny,Unz),Wstar(Wnx/2,Wny,Wnz)
        write(*,*) "U,V,Wstar",maxval(Ustar(1:Unx,Uny,Unz)),maxval(Wstar(1:Wnx,Wny,Wnz))
        write(*,*) "U,V,Wstar",maxloc(Ustar(1:Unx,Uny,Unz)),maxloc(Wstar(1:Wnx,Wny,Wnz))
#endif

      else

        !$omp parallel private(i,j,k)
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Unz
         do j=1,Uny
          do i=1,Unx
           Ustar(i,j,k) = 0
          enddo
         enddo
        enddo
        !$omp end do nowait
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Vnz
         do j=1,Vny
          do i=1,Vnx
           Vstar(i,j,k) = 0
          enddo
         enddo
        enddo
        !$omp end do nowait
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Wnz
         do j=1,Wny
          do i=1,Wnx
           Wstar(i,j,k) = 0
          enddo
         enddo
        enddo
        !$omp end do
        !$omp end parallel

      endif


      if (abs(coriolisparam)>1E-10_KND) call CoriolisForce(Unx,Uny,Unz,Vnx,Vny,Vnz,&
                                                    coriolisparam,&
                                                    Ustar,Vstar,U,V,dt)

      if (buoyancy==1) call BuoyancyForce(Prnx,Prny,Prnz,Wnx,Wny,Wnz,&
                                grav_acc,temperature_ref,Wstar,temperature,dt)


      !$omp parallel private(i,j,k)
      !$omp do
      !$hmppcg grid blocksize myblocksize
      !$hmppcg permute (k,i,j)
      !$hmppcg gridify(k,i)
      do k=1,Unz
       do j=1,Uny
        do i=1,Unx
         U2(i,j,k) = U2(i,j,k)+Ustar(i,j,k)*beta(lev)
        enddo
       enddo
      enddo
      !$omp end do nowait
      !$omp do
      !$hmppcg grid blocksize myblocksize
      !$hmppcg permute (k,i,j)
      !$hmppcg gridify(k,i)
      do k=1,Vnz
       do j=1,Vny
        do i=1,Vnx
         V2(i,j,k) = V2(i,j,k)+Vstar(i,j,k)*beta(lev)
        enddo
       enddo
      enddo
      !$omp end do nowait
      !$omp do
      !$hmppcg grid blocksize myblocksize
      !$hmppcg permute (k,i,j)
      !$hmppcg gridify(k,i)
      do k=1,Wnz
       do j=1,Wny
        do i=1,Wnx
         W2(i,j,k) = W2(i,j,k)+Wstar(i,j,k)*beta(lev)
        enddo
       enddo
      enddo
      !$omp end do
      !$omp end parallel

  end subroutine Convection











  !$hmpp <tsteps> TimeStepEul codelet
  subroutine TimeStepEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
               dxmin,dxU,dyV,dzW,CFL,Uref,steady,time,endtime,&
               U,V,W,dt)
  implicit none
#ifdef __HMPP
#include "hmpp-include.f90"
  intrinsic min,max,abs
#endif
  integer,intent(in)   :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,steady
  real(KND),intent(in) :: dxmin
  real(KND),intent(in) :: CFL,Uref
  real(TIM),intent(in) :: time,endtime
#ifdef __HMPP
  real(KND),intent(in) :: dxU(-2:Prnx+2),dyV(-2:Prny+2),dzW(-2:Prnz+2)
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(in)  :: U
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(in)  :: V
  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(in)  :: W
#else
  real(KND),dimension(-2:),intent(in)  :: dxU,dyV,dzW
  real(KND),dimension(-2:,-2:,-2:),intent(in)  :: U,V,W
#endif
  real(TIM),intent(out) :: dt
  integer i,j,k
  real(KND) m,p

    m = 0
   !$omp parallel do private(i,j,k,p) reduction(max:m)
   !$hmppcg grid blocksize myblocksize2
   !$hmppcg permute (k,i,j)
   !$hmppcg gridify (k,i), private(p), reduce (max:m)
   do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       p = MAX(abs(U(i,j,k)/dxU(i)),abs(U(i-1,j,k)/dxU(i-1)))
       p = MAX(p,abs(V(i,j,k)/dyV(j)),abs(V(i,j-1,k)/dyV(j-1)))
       p = MAX(p,abs(W(i,j,k)/dzW(k)),abs(W(i,j,k-1)/dzW(k-1)))
       m = max(m,p)
      enddo
     enddo
    enddo
    !$omp end parallel do

    if (m>0) then
     dt = MIN(CFL/m,dxmin/Uref)
    else
     dt = dxmin/Uref
    endif

    if (steady/=1.and.dt+time>endtime)  dt = endtime-time

  endsubroutine TimeStepEul







  subroutine BuoyancyForce(Prnx,Prny,Prnz,Wnx,Wny,Wnz,grav_acc,temperature_ref,W2,theta,dt)
  implicit none
#ifdef __HMPP
#include "hmpp-include.f90"
#endif
  integer,intent(in) :: Prnx,Prny,Prnz,Wnx,Wny,Wnz
#ifdef __HMPP
  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(inout) :: W2
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(in) :: theta
#else
  real(KND),dimension(-2:,-2:,-2:),intent(inout) :: W2
  real(KND),dimension(-1:,-1:,-1:),intent(in) :: theta
#endif
  real(KND),intent(in) :: grav_acc,temperature_ref,dt
  real(KND) A
  integer i,j,k

    A = grav_acc*dt/temperature_ref
    !$omp parallel do private(i,j,k)
    !$hmppcg grid blocksize myblocksize
    !$hmppcg permute (k,i,j)
    !$hmppcg gridify(k,i)
    do k=1,Wnz
     do j=1,Wny
      do i=1,Wnx
            W2(i,j,k) = W2(i,j,k)+A*((theta(i,j,k+1)+theta(i,j,k))/2._KND-temperature_ref)
      enddo
     enddo
    enddo
    !$omp end parallel do
  endsubroutine BuoyancyForce





  subroutine CoriolisForce(Unx,Uny,Unz,Vnx,Vny,Vnz,&
                                 coriolisparam,&
                                 U2,V2,U,V,dt)
  implicit none
#ifdef __HMPP
#include "hmpp-include.f90"
#endif
  integer,intent(in) :: Unx,Uny,Unz,Vnx,Vny,Vnz
#ifdef __HMPP
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(in)    :: U
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(in)    :: V
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(inout) :: U2
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(inout) :: V2
#else
  real(KND),dimension(-2:,-2:,-2:),intent(in)    :: U,V
  real(KND),dimension(-2:,-2:,-2:),intent(inout) :: U2,V2
#endif
  real(KND),intent(in):: coriolisparam,dt
  real(KND) A
  integer i,j,k

   A=-dt
   if (coriolisparam>0) then
   !$omp parallel private(i,j,k)
   !$omp do
   !$hmppcg grid blocksize myblocksize
   !$hmppcg permute (k,i,j)
   !$hmppcg gridify(k,i)
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
           U2(i,j,k) = U2(i,j,k)-A*coriolisparam*(V(i,j-1,k)+V(i+1,j-1,k)+V(i,j,k)+V(i+1,j,k))/4._KND
      enddo
     enddo
    enddo
    !$omp end do nowait

    !$omp do
    !$hmppcg grid blocksize myblocksize
    !$hmppcg permute (k,i,j)
    !$hmppcg gridify(k,i)
    do k=1,Vnz
     do j=1,Vny
      do i=1,Vnx
           V2(i,j,k) = V2(i,j,k)+A*coriolisparam*(U(i-1,j,k)+U(i-1,j+1,k)+U(i,j,k)+U(i,j+1,k))/4._KND
      enddo
     enddo
    enddo
    !$omp end do
    !$omp end parallel
   endif
  endsubroutine CoriolisForce








#ifdef __HMPP
  !$hmpp <tsteps> ForwEul codelet
  subroutine ForwEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                    dxPr,dyPr,dzPr,dxU,dyV,dzW,&
                    U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                    dt,coef)

  implicit none

#include "hmpp-include.f90"

  integer,intent(in) :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz
  real(KND),intent(in):: U(-2:Unx+3,-2:Uny+3,-2:Unz+3),V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(KND),intent(in):: U2(-2:Unx+3,-2:Uny+3,-2:Unz+3),V2(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),W2(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(KND),intent(out):: U3(-2:Unx+3,-2:Uny+3,-2:Unz+3),V3(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),W3(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(KND),intent(in):: Visc(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  real(KND),intent(in) :: dxU(-2:Prnx+2),dyV(-2:Prny+2),dzW(-2:Prnz+2)
  real(KND),intent(in) :: dxPr(-2:Prnx+3),dyPr(-2:Prny+3),dzPr(-2:Prnz+3),dt,coef
  real(KND) :: Ap,Ax,Ay,Az,dxmin,dymin,dzmin
  integer i,j,k


     dxmin = dxPr(1)
     dymin = dyPr(1)
     dzmin = dzPr(1)

     Ap = coef*dt

     Ax = 1._KND/(dxmin**2)
     Ay = 1._KND/(dymin**2)
     Az = 1._KND/(dzmin**2)


     !$hmppcg grid blocksize myblocksize
     !$hmppcg permute (k,i,j)
     !$hmppcg gridify(k,i)
     do k=1,Unz
      do j=1,Uny
       do i=1,Unx
         U3(i,j,k) =&
             (Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))-&
             Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k)))*Ax
         U3(i,j,k) = U3(i,j,k) +&
             Ay*0.25_KND*((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))-&
             (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k)))
         U3(i,j,k) = U3(i,j,k) +&
             Az*0.25_KND*((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))-&
             (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1)))
         U3(i,j,k) = U3(i,j,k) * Ap
       enddo
      enddo
     enddo
     !$hmppcg grid blocksize myblocksize
     !$hmppcg permute (k,i,j)
     !$hmppcg gridify(k,i)
     do k=1,Unz    !Forward Euler for the first approximation
      do j=1,Uny
       do i=1,Unx
         U3(i,j,k) = U3(i,j,k) + U(i,j,k) + U2(i,j,k)
       enddo
      enddo
     enddo

     !$hmppcg grid blocksize myblocksize
     !$hmppcg permute (k,i,j)
     !$hmppcg gridify(k,i)
     do k=1,Vnz
      do j=1,Vny
       do i=1,Vnx
         V3(i,j,k) =&
             (Visc(i,j+1,k)*(V(i,j+1,k)-V(i,j,k))-&
             Visc(i,j,k)*(V(i,j,k)-V(i,j-1,k)))*Ay
         V3(i,j,k) = V3(i,j,k) +&
             Ax*0.25_KND*((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V(i+1,j,k)-V(i,j,k))-&
             (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(V(i,j,k)-V(i-1,j,k)))
         V3(i,j,k) = V3(i,j,k) +&
             Az*0.25_KND*((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))-&
             (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(V(i,j,k)-V(i,j,k-1)))
         V3(i,j,k) = V3(i,j,k) * Ap
       enddo
      enddo
     enddo
     !$hmppcg grid blocksize myblocksize
     !$hmppcg permute (k,i,j)
     !$hmppcg gridify(k,i)
     do k=1,Vnz
      do j=1,Vny
       do i=1,Vnx
         V3(i,j,k) = V3(i,j,k) + V(i,j,k) + V2(i,j,k)
       enddo
      enddo
     enddo

     !$hmppcg grid blocksize myblocksize
     !$hmppcg permute (k,i,j)
     !$hmppcg gridify(k,i)
     do k=1,Wnz
      do j=1,Wny
       do i=1,Wnx
         W3(i,j,k) =&
             Ax*0.25_KND*((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W(i+1,j,k)-W(i,j,k))-&
             (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(W(i,j,k)-W(i-1,j,k)))
         W3(i,j,k) = W3(i,j,k) +&
             Ay*0.25_KND*((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W(i,j+1,k)-W(i,j,k))-&
             (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(W(i,j,k)-W(i,j-1,k)))
         W3(i,j,k) = W3(i,j,k) +&
             (Visc(i,j,k+1)*(W(i,j,k+1)-W(i,j,k))-&
             Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1)))*Az
         W3(i,j,k) = W3(i,j,k) * Ap
       enddo
      enddo
     enddo
     !$hmppcg grid blocksize myblocksize
     !$hmppcg permute (k,i,j)
     !$hmppcg gridify(k,i)
     do k=1,Wnz
      do j=1,Wny
       do i=1,Wnx
         W3(i,j,k) = W3(i,j,k) + W(i,j,k) + W2(i,j,k)
       enddo
      enddo
     enddo

  end subroutine ForwEul



#else



  subroutine ForwEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                    dxPr,dyPr,dzPr,dxU,dyV,dzW,&
                    U,V,W,U2,V2,W2,U3,V3,W3,Visc,&
                    dt,coef)

  implicit none

  integer,intent(in) :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz


  real(KND),intent(in),dimension(-2:,-2:,-2:) :: U,V,W
  real(KND),intent(in),dimension(-2:,-2:,-2:) :: U2,V2,W2
  real(KND),intent(out),dimension(-2:,-2:,-2:):: U3,V3,W3
  real(KND),intent(in),dimension(-1:,-1:,-1:) :: Visc
  real(KND),intent(in),dimension(-2:) :: dxU,dyV,dzW
  real(KND),intent(in),dimension(-2:) :: dxPr,dyPr,dzPr
  real(KND),intent(in) :: dt,coef

  real(KND) :: Ap,Ax,Ay,Az
  integer i,j,k


     Ap = coef*dt

     Ax = 1._KND/(dxmin**2)
     Ay = 1._KND/(dymin**2)
     Az = 1._KND/(dzmin**2)

     if (gridtype==uniformgrid) then
       !$omp parallel private(i,j,k)

       !$omp do
       do k=1,Unz
        do j=1,Uny
         do i=1,Unx
           U3(i,j,k) =&
               (Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))-&
               Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k)))*Ax
           U3(i,j,k) = U3(i,j,k) +&
               Ay*0.25_KND*((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))-&
               (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k)))
           U3(i,j,k) = U3(i,j,k) +&
               Az*0.25_KND*((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))-&
               (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1)))
           U3(i,j,k) = U3(i,j,k) * Ap
         enddo
        enddo
       enddo
       !$omp end do
       !$omp do
       do k=1,Unz
        do j=1,Uny
         do i=1,Unx
           U3(i,j,k) = U3(i,j,k) + U(i,j,k) + U2(i,j,k)
         enddo
        enddo
       enddo
       !$omp end do nowait


       !$omp do
       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
           V3(i,j,k) =&
               (Visc(i,j+1,k)*(V(i,j+1,k)-V(i,j,k))-&
               Visc(i,j,k)*(V(i,j,k)-V(i,j-1,k)))*Ay
           V3(i,j,k) = V3(i,j,k) +&
               Ax*0.25_KND*((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V(i+1,j,k)-V(i,j,k))-&
               (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(V(i,j,k)-V(i-1,j,k)))
           V3(i,j,k) = V3(i,j,k) +&
               Az*0.25_KND*((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))-&
               (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(V(i,j,k)-V(i,j,k-1)))
           V3(i,j,k) = V3(i,j,k) * Ap
         enddo
        enddo
       enddo
       !$omp end do
       !$omp do
       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
           V3(i,j,k) = V3(i,j,k) + V(i,j,k) + V2(i,j,k)
         enddo
        enddo
       enddo
       !$omp end do nowait


       !$omp do
       do k=1,Wnz
        do j=1,Wny
         do i=1,Wnx
           W3(i,j,k) =+&
               Ax*0.25_KND*((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W(i+1,j,k)-W(i,j,k))-&
               (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(W(i,j,k)-W(i-1,j,k)))
           W3(i,j,k) = W3(i,j,k) +&
               Ay*0.25_KND*((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W(i,j+1,k)-W(i,j,k))-&
               (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(W(i,j,k)-W(i,j-1,k)))
           W3(i,j,k) = W3(i,j,k) +&
               (Visc(i,j,k+1)*(W(i,j,k+1)-W(i,j,k))-&
               Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1)))*Az
           W3(i,j,k) = W3(i,j,k) * Ap
         enddo
        enddo
       enddo
       !$omp end do
       !$omp do
       do k=1,Wnz
        do j=1,Wny
         do i=1,Wnx
           W3(i,j,k) = W3(i,j,k) + W(i,j,k) + W2(i,j,k)
         enddo
        enddo
       enddo
       !$omp end do

       !$omp end parallel
     else
       !$omp parallel private(i,j,k)

       !$omp do
       do k=1,Unz    !Forward Euler for the first approximation
        do j=1,Uny
         do i=1,Unx
          U3(i,j,k) = U(i,j,k)+U2(i,j,k)+Ap*(&
          ((Visc(i+1,j,k)*(U(i+1,j,k)-U(i,j,k))/dxPr(i+1)-&
          Visc(i,j,k)*(U(i,j,k)-U(i-1,j,k))/dxPr(i))/dxU(i)+&
           0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(U(i,j+1,k)-U(i,j,k))/dyV(j)-&
           (Visc(i+1,j,k)+Visc(i+1,j-1,k)+Visc(i,j,k)+Visc(i,j-1,k))*(U(i,j,k)-U(i,j-1,k))/dyV(j-1))/dyPr(j)+&
           ((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)-&
           (Visc(i+1,j,k)+Visc(i+1,j,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(U(i,j,k)-U(i,j,k-1))/dzW(k-1))/dzPr(k))))
         enddo
        enddo
       enddo
       !$omp end do nowait
       !$omp do
       do k=1,Vnz
        do j=1,Vny
         do i=1,Vnx
          V3(i,j,k) = V(i,j,k)+V2(i,j,k)+Ap*(&
          ((Visc(i,j+1,k)*(V(i,j+1,k)-V(i,j,k))/dyPr(j+1)-&
           Visc(i,j,k)*(V(i,j,k)-V(i,j-1,k))/dyPr(j))/dyV(j)+&
           0.25_KND*(((Visc(i+1,j+1,k)+Visc(i+1,j,k)+Visc(i,j+1,k)+Visc(i,j,k))*(V(i+1,j,k)-V(i,j,k))/dxU(i)-&
          (Visc(i,j+1,k)+Visc(i,j,k)+Visc(i-1,j+1,k)+Visc(i-1,j,k))*(V(i,j,k)-V(i-1,j,k))/dxU(i-1))/dxPr(i)+&
           ((Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))/dzW(k)-&
           (Visc(i,j+1,k)+Visc(i,j+1,k-1)+Visc(i,j,k)+Visc(i,j,k-1))*(V(i,j,k)-V(i,j,k-1))/dzW(k-1))/dzPr(k))))
         enddo
        enddo
       enddo
       !$omp end do nowait
       !$omp do
       do k=1,Wnz
        do j=1,Wny
         do i=1,Wnx
          W3(i,j,k) = W(i,j,k)+W2(i,j,k)+Ap*(&
          0.25_KND*(((Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(W(i+1,j,k)-W(i,j,k))/dxU(i)-&
          (Visc(i,j,k+1)+Visc(i,j,k)+Visc(i-1,j,k+1)+Visc(i-1,j,k))*(W(i,j,k)-W(i-1,j,k))/dxU(i-1))/dxPr(i)+&
           ((Visc(i,j+1,k+1)+Visc(i,j,k+1)+Visc(i,j+1,k)+Visc(i,j,k))*(W(i,j+1,k)-W(i,j,k))/dyV(j)-&
           (Visc(i,j,k+1)+Visc(i,j-1,k+1)+Visc(i,j,k)+Visc(i,j-1,k))*(W(i,j,k)-W(i,j-1,k))/dyV(j-1))/dyPr(j))+&
           (Visc(i,j,k+1)*(W(i,j,k+1)-W(i,j,k))/dzPr(k+1)-&
           Visc(i,j,k)*(W(i,j,k)-W(i,j,k-1))/dzPr(k))/dzW(k))
         enddo
        enddo
       enddo
       !$omp end do

       !$omp end parallel
     end if
  end subroutine ForwEul




#endif


















  !Slower using HMPP than CPU, but if we avoid memory transfer, it is still profitable.

  !$hmpp <tsteps> AttenuateTop codelet
  subroutine AttenuateTop(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Btype,zPr,zW,U,V,W,temperature,buoyancy)

  implicit none

#include "hmpp-include.f90"

  integer,intent(in)      :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,buoyancy
  integer,intent(in)      :: Btype(6)
#ifdef __HMPP
  real(KND),intent(in)    :: zPr(-2:Prnz+3)
  real(KND),intent(in)    :: zW(-3:Prnz+3)
  real(KND),intent(inout) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
  real(KND),intent(inout) :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
  real(KND),intent(inout) :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(inout) :: temperature
#else
  real(KND),intent(in)    :: zPr(-2:)
  real(KND),intent(in)    :: zW(-3:)
  real(KND),intent(inout),dimension(-2:,-2:,-2:) :: U,V,W
  real(KND),dimension(-1:,-1:,-1:),intent(inout) :: temperature
#endif
  integer i,j,k,bufn,mini,maxi,maxUi
  real(KND) ze,zs,zb,p
  real(KND),dimension(:),allocatable :: DF,avg
  intrinsic max,min

    if (Btype(We)==DIRICHLET.or.Btype(We)==TURBULENTINLET.or.Btype(We)==INLETFROMFILE) then
      mini = min(5,Unx)
    else
      mini = 1
    endif

    if (Btype(Ea)==DIRICHLET.or.Btype(We)==TURBULENTINLET.or.Btype(We)==OUTLETBUFF) then
      maxi = max(1,Prnx-5)
      maxUi = max(1,Unx-5)
    else
      maxi = Prnx
      maxUi = Unx
    endif

    bufn = max(5,Prnz/4)
    zs = zW(Prnz-bufn)
    ze = zW(Prnz)

    allocate(DF(min(Unz,Vnz,Wnz)-bufn:max(Unz,Vnz,Wnz)))
    allocate(avg(min(Unz,Vnz,Wnz)-bufn:max(Unz,Vnz,Wnz)))



    do k=Unz-bufn,Unz
      avg(k) = 0
    enddo

    !$omp parallel private(i,j,k,p,zb)
    !$omp do
    do k=Unz-bufn,Unz
      p = 0
      !$hmppcg grid blocksize myblocksize
      !$hmppcg gridify (j,i) global(p), reduce(+:p)
      do j=1,Uny
        do i=mini,maxUi
          p = p+U(i,j,k)
        enddo
      enddo
      avg(k) = p
    enddo
    !$omp end do

    !$omp do
    do k=Unz-bufn,Unz
      avg(k) = avg(k)/((maxUi-mini+1)*Uny)
    enddo
    !$omp end do

    !$omp do
    do k=Unz-bufn,Unz
      zb=(zPr(k)-zs)/(ze-zs)
      DF(k) = DampF(zb)
    enddo
    !$omp end do

    !$omp do
    !$hmppcg grid blocksize myblocksize
    !$hmppcg permute(k,i,j)
    !$hmppcg gridify (k,i)
    do k=Unz-bufn,Unz
      do j=-1,Uny+1
        do i=-1,Unx+1
          U(i,j,k) = avg(k)+DF(k)*(U(i,j,k)-avg(k))
        enddo
      enddo
    enddo
    !$omp end do


    !$omp do
    do k=Vnz-bufn,Vnz
      avg(k) = 0
    enddo
    !$omp end do

    !$omp do
    do k=Vnz-bufn,Vnz
      p = 0
      !$hmppcg grid blocksize myblocksize
      !$hmppcg gridify (j,i) global(p), reduce(+:p)
      do j=1,Vny
        do i=mini,maxi
          p = p+V(i,j,k)
        enddo
      enddo
      avg(k) = p
    enddo
    !$omp end do

    !$omp do
    do k=Vnz-bufn,Vnz
      avg(k) = avg(k)/((maxi-mini+1)*Vny)
    enddo
    !$omp end do

    !$omp do
    do k=Vnz-bufn,Vnz
      zb=(zPr(k)-zs)/(ze-zs)
      DF(k) = DampF(zb)
    enddo
    !$omp end do

    !$omp do
    !$hmppcg grid blocksize myblocksize
    !$hmppcg permute(k,i,j)
    !$hmppcg gridify (k,i)
    do k=Vnz-bufn,Vnz
      do j=-1,Vny+1
        do i=-1,Vnx+1
          V(i,j,k) = avg(k)+DF(k)*(V(i,j,k)-avg(k))
        enddo
      enddo
    enddo
    !$omp end do


    !$omp do
    do k=Wnz-bufn,Wnz
      avg(k) = 0
    enddo
    !$omp end do

    !$omp do
    do k=Wnz-bufn,Wnz
      p = 0
      !$hmppcg grid blocksize myblocksize
      !$hmppcg gridify (j,i) global(p), reduce(+:p)
      do j=1,Wny
        do i=mini,maxi
          p = p+W(i,j,k)
        enddo
      enddo
      avg(k) = p
    enddo
    !$omp end do

    !$omp do
    do k=Wnz-bufn,Wnz
      avg(k) = avg(k)/((maxi-mini+1)*Wny)
    enddo
    !$omp end do

    !$omp do
    do k=Wnz-bufn,Wnz
      zb=(zW(k)-zs)/(ze-zs)
      DF(k) = DampF(zb)
    enddo
    !$omp end do

    !$omp do
    !$hmppcg grid blocksize myblocksize
    !$hmppcg permute(k,i,j)
    !$hmppcg gridify (k,i)
    do k=Wnz-bufn,Wnz
      do j=-1,Wny+1
        do i=-1,Wnx+1
          W(i,j,k) = avg(k)+DF(k)*(W(i,j,k)-avg(k))
        enddo
      enddo
    enddo
    !$omp end do



    if (buoyancy==1) then

      !$omp do
      do k=Prnz-bufn,Prnz
        avg(k) = 0
      enddo
      !$omp end do

      !$omp do
      do k=Prnz-bufn,Prnz
        p = 0
        !$hmppcg grid blocksize myblocksize
        !$hmppcg gridify (j,i) global(p), reduce(+:p)
        do j=1,Prny
          do i=mini,maxi
            p = p+temperature(i,j,k)
          enddo
        enddo
        avg(k) = p
      enddo
      !$omp end do

      !$omp do
      do k=Prnz-bufn,Prnz
        avg(k) = avg(k)/((maxi-mini+1)*Prny)
      enddo
      !$omp end do

      !$omp do
      do k=Prnz-bufn,Prnz
        zb=(zPr(k)-zs)/(ze-zs)
        DF(k) = DampF(zb)
      enddo
      !$omp end do

      !$omp do
      !$hmppcg grid blocksize myblocksize
      !$hmppcg permute(k,i,j)
      !$hmppcg gridify (k,i)
      do k=Prnz-bufn,Prnz
        do j=-1,Prny+1
          do i=-1,Prnx+1
            temperature(i,j,k) = avg(k)+DF(k)*(temperature(i,j,k)-avg(k))
          enddo
        enddo
      enddo
      !$omp end do

    endif
    !$omp end parallel

  endsubroutine AttenuateTop



  !$hmpp <tsteps> AttenuateOut codelet
  !$hmpp <tsteps> AttenuateOut2 codelet
  subroutine AttenuateOut(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,xPr,xU,U,V,W,temperature,buoyancy)

  implicit none

#include "hmpp-include.f90"

  integer,intent(in)      :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,buoyancy
#ifdef __HMPP
  real(KND),intent(in)    :: xPr(-2:Prnx+3)
  real(KND),intent(in)    :: xU(-3:Prnx+3)
  real(KND),intent(inout) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
  real(KND),intent(inout) :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
  real(KND),intent(inout) :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(inout) :: temperature
#else
  real(KND),intent(in)    :: xPr(-2:)
  real(KND),intent(in)    :: xU(-3:)
  real(KND),intent(inout),dimension(-2:,-2:,-2:) :: U,V,W
  real(KND),dimension(-1:,-1:,-1:),intent(inout) :: temperature
#endif
  integer i,j,k,bufn
  real(KND) p,xe,xs,xb,DF
  intrinsic max

    bufn = min(max(10,Prnx/8),Prnx/2)
    xs = xU(Prnx-bufn)
    xe = xU(Prnx)

    !$omp parallel private(i,j,k,p,xb,DF)

    !$omp do
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify (k,j) private(p,xb,DF)
    do k=1,Unz
      do j=1,Uny
        p = 0
        do i=2*Unx/3,Unx-4
          p = p+U(i,j,k)
        enddo
        p = p/(Unx-4-2*Unx/3+1)
        do i=Unx-bufn,Unx+1
          xb=(xU(i)-xs)/(xe-xs)
          DF = DampF(xb)
          U(i,j,k) = p+DF*(U(i,j,k)-p)
        enddo
      enddo
    enddo
    !$omp end do

    !$omp do
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify (k,j) private(p,xb,DF)
    do k=1,Vnz
      do j=1,Vny
        p = 0
        do i=2*Vnx/3,Vnx-4
          p = p+V(i,j,k)
        enddo
        p = p/(Vnx-4-2*Vnx/3+1)
        do i=Vnx-bufn,Vnx+1
          xb=(xPr(i)-xs)/(xe-xs)
          DF = DampF(xb)
          V(i,j,k) = p+DF*(V(i,j,k)-p)
        enddo
      enddo
    enddo
    !$omp end do

    !$omp do
    !$hmppcg grid blocksize myblocksize
    !$hmppcg gridify (k,j) private(p,xb,DF)
    do k=1,Wnz
      do j=1,Wny
        p = 0
        do i=2*Wnx/3,Wnx-4
          p = p+W(i,j,k)
        enddo
        p = p/((Wnx-4-2*Wnx/3+1))
        do i=Wnx-bufn,Wnx+1
          xb=(xPr(i)-xs)/(xe-xs)
          DF = DampF(xb)
          W(i,j,k) = p+DF*(W(i,j,k)-p)
        enddo
      enddo
    enddo
    !$omp end do

    if (buoyancy==1) then
      !$omp do
      !$hmppcg grid blocksize myblocksize
      !$hmppcg gridify (k,j) private(p,xb,DF)
      do k=1,Prnz
        do j=1,Prny
          p = 0
          do i=2*Prnx/3,Prnx-4
            p = p+temperature(i,j,k)
          enddo
          p = p/(Prnx-4-2*Prnx/3+1)
          do i=Prnx-bufn,Prnx+1
            xb=(xPr(i)-xs)/(xe-xs)
            DF = DampF(xb)
            temperature(i,j,k) = p+DF*(temperature(i,j,k)-p)
          enddo
        enddo
      enddo
      !$omp end do
    endif
    !$omp end parallel

  endsubroutine AttenuateOut



  pure function DampF(x)

  implicit none

#include "hmpp-include.f90"

  real(KND) DampF
  real(KND),intent(in)::x
  intrinsic exp

  if (x<=0) then
    DampF = 1
  elseif (x>=1) then
    DampF = 0
  else
   DampF=(1-0.04_KND*x**2)*(1-(1-exp(10._KND*x**2))/(1-exp(10._KND)))
  endif
  endfunction Dampf


  !$hmpp <tsteps> NullInterior codelet
  subroutine NullInterior(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                          nUnull,nVnull,nWnull,Unull,Vnull,Wnull,U,V,W)

  implicit none

#include "hmpp-include.f90"

  integer,intent(in)      :: Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,nUnull,nVnull,nWnull
#ifdef __HMPP
  integer,dimension(3,nUnull),intent(in)      :: Unull
  integer,dimension(3,nVnull),intent(in)      :: Vnull
  integer,dimension(3,nWnull),intent(in)      :: Wnull
  real(KND),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(inout) :: U
  real(KND),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(inout) :: V
  real(KND),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(inout) :: W
#else
  integer,dimension(:,:),intent(in)      :: Unull,Vnull,Wnull
  real(KND),dimension(-2:,-2:,-2:),intent(inout) :: U,V,W
#endif
  integer i

    !$omp parallel private(i)
    !$omp do
    !$hmppcg gridify(i)
    do i=1,nUnull
      U(Unull(1,i),Unull(2,i),Unull(3,i)) = 0
    enddo
    !$omp end do nowait

    !$omp do
    !$hmppcg gridify(i)
    do i=1,nVnull
      V(Vnull(1,i),Vnull(2,i),Vnull(3,i)) = 0
    enddo
    !$omp end do nowait

    !$omp do
    !$hmppcg gridify(i)
    do i=1,nWnull
      W(Wnull(1,i),Wnull(2,i),Wnull(3,i)) = 0
    enddo
    !$omp end do
    !$omp end parallel

  endsubroutine NullInterior




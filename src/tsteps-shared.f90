
 !$hmpp <tsteps> PressureGrad codelet
  subroutine PressureGrad(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxU,dyV,dzW,&
                          Btype,prgradientx,prgradienty,&
                          Pr,U,V,W,&
                          coef)
  implicit none
#ifdef __HMPP
#include "hmpp-include.f90"
#endif
  integer,intent(in)      :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz
  real(knd),intent(in)    :: prgradientx,prgradienty,coef
  integer,intent(in)      :: Btype(6)
#ifdef __HMPP
  real(knd),intent(in)    :: dxU(-2:Prnx+2),dyV(-2:Prny+2),dzW(-2:Prnz+2)
  real(knd),intent(in)    :: Pr(1:Unx+1,1:Vny+1,1:Wnz+1)
  real(knd),intent(inout) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
  real(knd),intent(inout) :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
  real(knd),intent(inout) :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
#else
  real(knd),intent(in),contiguous,dimension(-2:) :: dxU,dyV,dzW
  real(knd),intent(inout),contiguous :: Pr(1:,1:,1:)
  real(knd),intent(inout),contiguous,dimension(-2:,-2:,-2:) :: U,V,W
#endif
  real(knd) :: A, Ax, Ay, Az
  integer i,j,k

   A = -coef
   Ax = - coef / dxmin
   Ay = - coef / dymin
   Az = - coef / dzmin

   !$omp parallel
   !$omp do
   !$hmppcg grid blocksize myblocksize
   !$hmppcg permute (k,i,j)
   !$hmppcg gridify(k,i)
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
          U(i,j,k) = U(i,j,k)+Ax * (Pr(i+1,j,k)-Pr(i,j,k)) + A * prgradientx
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
          V(i,j,k) = V(i,j,k)+Ay * (Pr(i,j+1,k)-Pr(i,j,k)) + A * prgradienty
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
          W(i,j,k) = W(i,j,k) + Az * (Pr(i,j,k+1)-Pr(i,j,k))
     enddo
    enddo
   enddo
   !$omp enddo nowait
   !$omp end parallel
  end subroutine PressureGrad


  
  
  subroutine StressBoundaryFlux(U2,V2,dt)
    use ArrayUtilities,only: add
    use Outputs,only: profuw,profvw, profuwsgs, profvwsgs
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout) :: U2,V2
    real(knd),intent(in) :: dt
    real(knd) :: flux
    integer :: first,last
    
    if (Btype(To)==AUTOMATICFLUX) then
      first = min(Prnz*5/6,Prnz-5)
      last = Prnz-5
      
      flux = sum(profuw(first:last)) + sum(profuwsgs(first:last))
      flux = flux / (last-first+1)
      call add(U2(:,:,Unz), -dt*flux/dzPr(Unz))
      
      flux = sum(profvw(first:last)) + sum(profvwsgs(first:last))
      flux = flux / (last-first+1)
      call add(V2(:,:,Vnz), -dt*flux/dzPr(Vnz))
    end if
  end subroutine




#ifdef __HMPP
 !$hmpp <tsteps> Convection codelet
 subroutine Convection(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,enable_buoyancy,convmet,&
                       dxmin,dymin,dzmin,coriolisparam,grav_acc,temperature_ref,&
                       U,V,W,U2,V2,W2,Ustar,Vstar,Wstar,temperature,beta,rho,RK_stage,dt)
  implicit none
#include "hmpp-include.f90"


  integer,intent(in)   :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,enable_buoyancy,convmet
  real(knd),intent(in) :: dxmin,dymin,dzmin,coriolisparam,grav_acc,temperature_ref
  real(knd),intent(in) :: dt

  real(knd),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(in)    :: U
  real(knd),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(out)   :: U2
  real(knd),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(inout) :: Ustar
  real(knd),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(in)    :: V
  real(knd),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(out)   :: V2
  real(knd),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(inout) :: Vstar
  real(knd),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(in)    :: W
  real(knd),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(out)   :: W2
  real(knd),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(inout) :: Wstar
  real(knd),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(in) :: temperature
  intrinsic abs
#else
 subroutine Convection(U,V,W,U2,V2,W2,Ustar,Vstar,Wstar,Temperature,Moisture,beta,rho,RK_stage,dt)
  use VolumeSources, only: ResistanceForce
  use VTKArray
  real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in)    :: U,V,W
  real(knd),dimension(-2:,-2:,-2:),contiguous,intent(out)   :: U2,V2,W2
  real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout) :: Ustar,Vstar,Wstar
  real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in)    :: Temperature
  real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in)    :: Moisture
#endif
  real(knd),dimension(1:3),intent(in) :: beta,rho
  integer,intent(in) :: RK_stage
  real(knd),intent(in) :: dt
  integer i,j,k

      if (RK_stage>1) then
        !$omp parallel private(i,j,k)
        !$omp do
        !$hmppcg grid blocksize myblocksize
        !$hmppcg permute (k,i,j)
        !$hmppcg gridify(k,i)
        do k=1,Unz
         do j=1,Uny
          do i=1,Unx
           U2(i,j,k) = Ustar(i,j,k)*rho(RK_stage)*dt
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
           V2(i,j,k) = Vstar(i,j,k)*rho(RK_stage)*dt
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
           W2(i,j,k) = Wstar(i,j,k)*rho(RK_stage)*dt
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
          call CDS_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,dxmin,dymin,dzmin,Ustar,Vstar,Wstar,U,V,W)
#else
        if (convmet==2) then
          call CDU(Ustar,U,V,W)
          call CDV(Vstar,U,V,W)
          call CDW(Wstar,U,V,W)
        else if (convmet==4) then
          call CD4divU(Ustar,U,V,W)
          call CD4divV(Vstar,U,V,W)
          call CD4divW(Wstar,U,V,W)
        end if
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

      call FilterUstar    

      if (abs(coriolisparam)>tiny(1._knd)) call CoriolisForce(Unx,Uny,Unz,Vnx,Vny,Vnz,&
                                                    coriolisparam,&
                                                    Ustar,Vstar,U,V)

      if (enable_buoyancy) call BuoyancyForce(Prnx,Prny,Prnz,Wnx,Wny,Wnz,&
                                grav_acc,temperature_ref,Wstar,Temperature,Moisture)

      call ResistanceForce(Ustar,Vstar,Wstar,U,V,W)

      call StressBoundaryFlux(Ustar,Vstar,dt)

      if (explicit_diffusion==1) call MomentumDiffusion(Ustar,Vstar,Wstar,U,V,W)


      !$omp parallel private(i,j,k)
      !$omp do
      !$hmppcg grid blocksize myblocksize
      !$hmppcg permute (k,i,j)
      !$hmppcg gridify(k,i)
      do k=1,Unz
       do j=1,Uny
        do i=1,Unx
         U2(i,j,k) = U2(i,j,k)+Ustar(i,j,k)*beta(RK_stage)*dt
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
         V2(i,j,k) = V2(i,j,k)+Vstar(i,j,k)*beta(RK_stage)*dt
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
         W2(i,j,k) = W2(i,j,k)+Wstar(i,j,k)*beta(RK_stage)*dt
        enddo
       enddo
      enddo
      !$omp end do
      !$omp end parallel

      contains

        subroutine FilterUstar

          use Filters, only: filtertype, Filter

          if (filtertype/=0) then
            call BoundU(1,Ustar,Uin,2)
            call BoundU(2,Vstar,Vin,2)
            call BoundU(3,Wstar,Win,2)


            call Filter(Ustar,Utype)

            call Filter(Vstar,Vtype)

            call Filter(Wstar,Wtype)

            call BoundU(1,Ustar,Uin,2)
            call BoundU(2,Vstar,Vin,2)
            call BoundU(3,Wstar,Win,2)
          end if

        end subroutine
  end subroutine Convection











  !$hmpp <tsteps> TimeStepEul codelet
  subroutine TimeStepEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
               dxmin,dxU,dyV,dzW,CFL,Uref,steady,time,endtime,&
               U,V,W,dt)
    use ieee_arithmetic

    integer,intent(in)   :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,steady
    real(knd),intent(in) :: dxmin
    real(knd),intent(in) :: CFL,Uref
    real(TIM),intent(in) :: time,endtime
    real(knd),dimension(-2:),contiguous,intent(in)  :: dxU,dyV,dzW
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in)  :: U,V,W
    real(TIM),intent(out) :: dt
    integer i,j,k
    real(knd) m, p
    logical nan

    nan = .false.
    m = 0
    !$omp parallel do private(i,j,k,p) reduction(max:m) reduction(.or.:nan)
    do k=1,Prnz
      do j=1,Prny
        do i=1,Prnx
          !For scalar advection the sum proved to be necessary when the flow is not aligned to grid.
          p = max(abs(U(i,j,k)/dxU(i)),abs(U(i-1,j,k)/dxU(i-1)))
          p = p + max(abs(V(i,j,k)/dyV(j)),abs(V(i,j-1,k)/dyV(j-1)))
          p = p + max(abs(W(i,j,k)/dzW(k)),abs(W(i,j,k-1)/dzW(k-1)))
          
          m = max(m,p)
          if (ieee_is_nan(p)) nan = .true.
        enddo
      enddo
    enddo
    !$omp end parallel do

    if (nan) then
      dt = tiny(dt)
    else if (m>0) then
      dt = min(CFL/m,max(dxmin,dymin,dzmin)/Uref)
    else
      dt = dxmin/Uref
    endif

    if (steady/=1.and.dt+time>endtime)  dt = endtime-time

  endsubroutine TimeStepEul







  subroutine BuoyancyForce(Prnx,Prny,Prnz,Wnx,Wny,Wnz,grav_acc,temperature_ref,W2,Temperature,Moisture)
    implicit none
#ifdef __HMPP
#include "hmpp-include.f90"
#endif
    integer,intent(in) :: Prnx,Prny,Prnz,Wnx,Wny,Wnz
#ifdef __HMPP
    real(knd),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),contiguous,intent(inout) :: W2
    real(knd),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),contiguous,intent(in) :: Temperature
#else
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout) :: W2
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in) :: Temperature,Moisture
#endif
    real(knd),intent(in) :: grav_acc,temperature_ref
    real(knd) A,A2,temperature_virt
    integer i,j,k


    if (enable_moisture) then
      if (enable_liquid) then
        stop "Liquid water not implemented."
      else
        A = grav_acc / temperature_ref
        A2 = A / 2._KND

        call apply_moist(1)
        call apply_moist(2)
      end if
    else
      A = grav_acc / temperature_ref
      A2 = A / 2._KND
      !$omp parallel do private(i,j,k)
      !$hmppcg grid blocksize myblocksize
      !$hmppcg permute (k,i,j)
      !$hmppcg gridify(k,i)
      do k=1,Wnz
       do j=1,Wny
        do i=1,Wnx
              W2(i,j,k) = W2(i,j,k) + A2 * (Temperature(i,j,k+1)+Temperature(i,j,k)) - A * temperature_ref
        enddo
       enddo
      enddo
      !$omp end parallel do
    end if

    contains

      subroutine apply_moist(start)
        integer,intent(in) :: start
        !$omp parallel do private(i,j,k,temperature_virt)
        do k=start,Wnz+1,2
         do j=1,Wny
          do i=1,Wnx
                temperature_virt = theta_v(i,j,k)
                W2(i,j,k)   = W2(i,j,k)   + A2 * temperature_virt - A * temperature_ref
                W2(i,j,k-1) = W2(i,j,k-1) + A2 * temperature_virt
          enddo
         enddo
        enddo
        !$omp end parallel do
      end subroutine

      real(knd) function theta_v(i,j,k)
        integer :: i,j,k

        theta_v = Temperature(i,j,k) * (1._knd + 0.61_knd * Moisture(i,j,k))
      end function
  endsubroutine BuoyancyForce





  subroutine CoriolisForce(Unx,Uny,Unz,Vnx,Vny,Vnz,&
                                 coriolisparam,&
                                 U2,V2,U,V)
  implicit none
#ifdef __HMPP
#include "hmpp-include.f90"
#endif
  integer,intent(in) :: Unx,Uny,Unz,Vnx,Vny,Vnz
#ifdef __HMPP
  real(knd),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(in)    :: U
  real(knd),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(in)    :: V
  real(knd),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(inout) :: U2
  real(knd),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(inout) :: V2
#else
  real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in)    :: U,V
  real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout) :: U2,V2
#endif
  real(knd),intent(in):: coriolisparam
  integer i,j,k

   if (coriolisparam>0) then
   !$omp parallel private(i,j,k)
   !$omp do
   !$hmppcg grid blocksize myblocksize
   !$hmppcg permute (k,i,j)
   !$hmppcg gridify(k,i)
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
           U2(i,j,k) = U2(i,j,k) + &
              coriolisparam*(V(i,j-1,k)+V(i+1,j-1,k)+V(i,j,k)+V(i+1,j,k))/4._knd
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
           V2(i,j,k) = V2(i,j,k) - &
              coriolisparam*(U(i-1,j,k)+U(i-1,j+1,k)+U(i,j,k)+U(i,j+1,k))/4._knd
      enddo
     enddo
    enddo
    !$omp end do
    !$omp end parallel
   endif
  endsubroutine CoriolisForce





  subroutine ForwEul(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                    dxPr,dyPr,dzPr,dxU,dyV,dzW,&
                    U,V,W,U2,V2,W2,U3,V3,W3,nu,&
                    coef)

  implicit none

  integer,intent(in) :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz


  real(knd),intent(in),dimension(-2:,-2:,-2:),contiguous :: U,V,W
  real(knd),intent(in),dimension(-2:,-2:,-2:),contiguous :: U2,V2,W2
  real(knd),intent(out),dimension(-2:,-2:,-2:),contiguous:: U3,V3,W3
  real(knd),intent(in),dimension(-1:,-1:,-1:),contiguous :: nu
  real(knd),intent(in),dimension(-2:),contiguous :: dxU,dyV,dzW
  real(knd),intent(in),dimension(-2:),contiguous :: dxPr,dyPr,dzPr
  real(knd),intent(in) :: coef

  real(knd) :: Ap, recdxmin2, recdymin2, recdzmin2
  integer i,j,k


     Ap = coef

     recdxmin2 = 1._knd / dxmin**2
     recdymin2 = 1._knd / dymin**2
     recdzmin2 = 1._knd / dzmin**2

     !$omp parallel private(i,j,k)

     !$omp do
     do k=1,Unz
      do j=1,Uny
       do i=1,Unx
            U3(i,j,k) = 0
            if (Uflx_mask(i+1,j,k)) &
              U3(i,j,k) = U3(i,j,k) + &
               nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) *recdxmin2
            if (Uflx_mask(i,j,k)) &
              U3(i,j,k) = U3(i,j,k) - &
                nu(i,j,k) * (U(i,j,k)-U(i-1,j,k)) * recdxmin2 
            if (Ufly_mask(i,j+1,k)) &
              U3(i,j,k) = U3(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (U(i,j+1,k)-U(i,j,k)) * recdymin2
            if (Ufly_mask(i,j,k)) &
              U3(i,j,k) = U3(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * (U(i,j,k)-U(i,j-1,k)) * recdymin2
            if (Uflz_mask(i,j,k+1)) &
              U3(i,j,k) = U3(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (U(i,j,k+1)-U(i,j,k)) * recdzmin2
            if (Uflz_mask(i,j,k)) &
              U3(i,j,k) = U3(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (U(i,j,k)-U(i,j,k-1)) * recdzmin2
       enddo
      enddo
     enddo
     !$omp end do
#define comp 1
#define wrk U3
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
     !$omp do
     do k=1,Unz
      do j=1,Uny
       do i=1,Unx
         U3(i,j,k) = U3(i,j,k) * Ap
         U3(i,j,k) = U3(i,j,k) + U(i,j,k) + U2(i,j,k)
       enddo
      enddo
     enddo
     !$omp end do nowait


     !$omp do
     do k=1,Vnz
      do j=1,Vny
       do i=1,Vnx
            V3(i,j,k) = 0
            if (Vflx_mask(i+1,j,k)) &
              V3(i,j,k) = V3(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (V(i+1,j,k)-V(i,j,k)) * recdxmin2
            if (Vflx_mask(i,j,k)) &
              V3(i,j,k) = V3(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * (V(i,j,k)-V(i-1,j,k)) * recdxmin2
            if (Vfly_mask(i,j+1,k)) &
              V3(i,j,k) = V3(i,j,k) + &
                nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k)) * recdymin2
            if (Vfly_mask(i,j,k)) &
              V3(i,j,k) = V3(i,j,k) - &
                nu(i,j,k) * (V(i,j,k)-V(i,j-1,k)) * recdymin2
            if (Vflz_mask(i,j,k+1)) &
              V3(i,j,k) = V3(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * (V(i,j,k+1)-V(i,j,k)) * recdzmin2
            if (Vflz_mask(i,j,k)) &
              V3(i,j,k) = V3(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (V(i,j,k)-V(i,j,k-1)) * recdzmin2
       enddo
      enddo
     enddo
     !$omp end do
#define comp 2
#define wrk V3
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
     !$omp do
     do k=1,Vnz
      do j=1,Vny
       do i=1,Vnx
         V3(i,j,k) = V3(i,j,k) * Ap
         V3(i,j,k) = V3(i,j,k) + V(i,j,k) + V2(i,j,k)
       enddo
      enddo
     enddo
     !$omp end do nowait


     !$omp do
     do k=1,Wnz
      do j=1,Wny
       do i=1,Wnx
            W3(i,j,k) = 0
            if (Wflx_mask(i+1,j,k)) &
              W3(i,j,k) = W3(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (W(i+1,j,k)-W(i,j,k)) * recdxmin2
            if (Wflx_mask(i,j,k)) &
              W3(i,j,k) = W3(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * (W(i,j,k)-W(i-1,j,k)) * recdxmin2
            if (Wfly_mask(i,j+1,k)) &
              W3(i,j,k) = W3(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * (W(i,j+1,k)-W(i,j,k)) * recdymin2
            if (Wfly_mask(i,j,k)) &
              W3(i,j,k) = W3(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * (W(i,j,k)-W(i,j-1,k)) * recdymin2
            if (Wflz_mask(i,j,k+1)) &
              W3(i,j,k) = W3(i,j,k) + &
                nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k)) * recdzmin2
            if (Wflz_mask(i,j,k)) &
              W3(i,j,k) = W3(i,j,k) - &
                nu(i,j,k) * (W(i,j,k)-W(i,j,k-1)) * recdzmin2
       enddo
      enddo
     enddo
     !$omp end do
#define comp 3
#define wrk W3
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
     !$omp do
     do k=1,Wnz
      do j=1,Wny
       do i=1,Wnx
         W3(i,j,k) = W3(i,j,k) * Ap
         W3(i,j,k) = W3(i,j,k) + W(i,j,k) + W2(i,j,k)
       enddo
      enddo
     enddo
     !$omp end do

     !$omp end parallel

  end subroutine ForwEul










  subroutine MomentumDiffusion(U2,V2,W2,U,V,W)
    use Parameters, nu => Viscosity
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in)    :: U,V,W
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout) :: U2,V2,W2
    real(knd) :: Ap, recdxmin2, recdymin2, recdzmin2
    integer :: i,j,k,bi,bj,bk
    integer, parameter :: narr = 6
       
     recdxmin2 = 1._knd / dxmin**2
     recdymin2 = 1._knd / dymin**2
     recdzmin2 = 1._knd / dzmin**2

     !$omp parallel private(i,j,k,bi,bj,bk)
     !$omp do schedule(runtime) !collapse(3)
     do bk = 1, Unz, tilenz(narr)
      do bj = 1, Uny, tileny(narr)
       do bi = 1, Unx, tilenx(narr)
        do k = bk, min(bk+tilenz(narr)-1,Unz)
         do j = bj, min(bj+tileny(narr)-1,Uny)
          do i = bi, min(bi+tilenx(narr)-1,Unx)
            if (Uflx_mask(i+1,j,k)) &
              U2(i,j,k) = U2(i,j,k) + &
               nu(i+1,j,k) * (U(i+1,j,k)-U(i,j,k)) *recdxmin2
            if (Uflx_mask(i,j,k)) &
              U2(i,j,k) = U2(i,j,k) - &
                nu(i,j,k) * (U(i,j,k)-U(i-1,j,k)) * recdxmin2 
            if (Ufly_mask(i,j+1,k)) &
              U2(i,j,k) = U2(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (U(i,j+1,k)-U(i,j,k)) * recdymin2
            if (Ufly_mask(i,j,k)) &
              U2(i,j,k) = U2(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j-1,k)+nu(i,j,k)+nu(i,j-1,k)) * (U(i,j,k)-U(i,j-1,k)) * recdymin2
            if (Uflz_mask(i,j,k+1)) &
              U2(i,j,k) = U2(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (U(i,j,k+1)-U(i,j,k)) * recdzmin2
            if (Uflz_mask(i,j,k)) &
              U2(i,j,k) = U2(i,j,k) - &
                0.25_knd * (nu(i+1,j,k)+nu(i+1,j,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (U(i,j,k)-U(i,j,k-1)) * recdzmin2
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
     !$omp end do
#define comp 1
#define wrk U2
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp


     !$omp do schedule(runtime) !collapse(3)
     do bk = 1, Vnz, tilenz(narr)
      do bj = 1, Vny, tileny(narr)
       do bi = 1, Vnx, tilenx(narr)
        do k = bk, min(bk+tilenz(narr)-1,Vnz)
         do j = bj, min(bj+tileny(narr)-1,Vny)
          do i = bi, min(bi+tilenx(narr)-1,Vnx)
            if (Vflx_mask(i+1,j,k)) &
              V2(i,j,k) = V2(i,j,k) + &
                0.25_knd * (nu(i+1,j+1,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i,j,k)) * (V(i+1,j,k)-V(i,j,k)) * recdxmin2
            if (Vflx_mask(i,j,k)) &
              V2(i,j,k) = V2(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j,k)+nu(i-1,j+1,k)+nu(i-1,j,k)) * (V(i,j,k)-V(i-1,j,k)) * recdxmin2
            if (Vfly_mask(i,j+1,k)) &
              V2(i,j,k) = V2(i,j,k) + &
                nu(i,j+1,k) * (V(i,j+1,k)-V(i,j,k)) * recdymin2
            if (Vfly_mask(i,j,k)) &
              V2(i,j,k) = V2(i,j,k) - &
                nu(i,j,k) * (V(i,j,k)-V(i,j-1,k)) * recdymin2
            if (Vflz_mask(i,j,k+1)) &
              V2(i,j,k) = V2(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j,k)) * (V(i,j,k+1)-V(i,j,k)) * recdzmin2
            if (Vflz_mask(i,j,k)) &
              V2(i,j,k) = V2(i,j,k) - &
                0.25_knd * (nu(i,j+1,k)+nu(i,j+1,k-1)+nu(i,j,k)+nu(i,j,k-1)) * (V(i,j,k)-V(i,j,k-1)) * recdzmin2
          enddo
         enddo
        enddo
       enddo
      enddo
     enddo
     !$omp end do
#define comp 2
#define wrk V2
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp


     !$omp do schedule(runtime) !collapse(3)
     do bk = 1, Wnz, tilenz(narr)
      do bj = 1, Wny, tileny(narr)
       do bi = 1, Wnx, tilenx(narr)
        do k = bk, min(bk+tilenz(narr)-1,Wnz)
         do j = bj, min(bj+tileny(narr)-1,Wny)
          do i = bi, min(bi+tilenx(narr)-1,Wnx)
            if (Wflx_mask(i+1,j,k)) &
              W2(i,j,k) = W2(i,j,k) + &
                0.25_knd * (nu(i+1,j,k+1)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i,j,k)) * (W(i+1,j,k)-W(i,j,k)) * recdxmin2
            if (Wflx_mask(i,j,k)) &
              W2(i,j,k) = W2(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j,k)+nu(i-1,j,k+1)+nu(i-1,j,k)) * (W(i,j,k)-W(i-1,j,k)) * recdxmin2
            if (Wfly_mask(i,j+1,k)) &
              W2(i,j,k) = W2(i,j,k) + &
                0.25_knd * (nu(i,j+1,k+1)+nu(i,j,k+1)+nu(i,j+1,k)+nu(i,j,k)) * (W(i,j+1,k)-W(i,j,k)) * recdymin2
            if (Wfly_mask(i,j,k)) &
              W2(i,j,k) = W2(i,j,k) - &
                0.25_knd * (nu(i,j,k+1)+nu(i,j-1,k+1)+nu(i,j,k)+nu(i,j-1,k)) * (W(i,j,k)-W(i,j-1,k)) * recdymin2
            if (Wflz_mask(i,j,k+1)) &
              W2(i,j,k) = W2(i,j,k) + &
                nu(i,j,k+1) * (W(i,j,k+1)-W(i,j,k)) * recdzmin2
            if (Wflz_mask(i,j,k)) &
              W2(i,j,k) = W2(i,j,k) - &
                nu(i,j,k) * (W(i,j,k)-W(i,j,k-1)) * recdzmin2
          end do
         end do
        end do
       end do
      end do
     end do
     !$omp end do
#define comp 3
#define wrk V2
#include "wmfluxes-inc.f90"
#undef wrk
#undef comp
     !$omp end parallel

  end subroutine MomentumDiffusion





















  !Slower using HMPP than CPU, but if we avoid memory transfer, it is still profitable.

  !$hmpp <tsteps> AttenuateTop codelet
  subroutine AttenuateTop(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Btype,zPr,zW,U,V,W)

  implicit none

  integer,intent(in)      :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz
  integer,intent(in)      :: Btype(6)
#ifdef __HMPP
#include "hmpp-include.f90"
  real(knd),intent(in)    :: zPr(-2:Prnz+3)
  real(knd),intent(in)    :: zW(-3:Prnz+3)
  real(knd),intent(inout) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
  real(knd),intent(inout) :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
  real(knd),intent(inout) :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
#else
  real(knd),contiguous,intent(in)    :: zPr(-2:)
  real(knd),contiguous,intent(in)    :: zW(-3:)
  real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout) :: U,V,W
#endif
  integer i,j,k,bufn,mini,maxi,maxUi
  real(knd) ze,zs,zb,p
  real(knd),dimension(:),allocatable :: DF,avg
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
    
    !$omp end parallel

  endsubroutine AttenuateTop



  !$hmpp <tsteps> AttenuateOut codelet
  !$hmpp <tsteps> AttenuateOut2 codelet
  subroutine AttenuateOut(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,xPr,xU,U,V,W,temperature,enable_buoyancy)

  implicit none

  integer,intent(in)      :: Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz
  logical, intent(in)     :: enable_buoyancy
#ifdef __HMPP
#include "hmpp-include.f90"
  real(knd),intent(in)    :: xPr(-2:Prnx+3)
  real(knd),intent(in)    :: xU(-3:Prnx+3)
  real(knd),intent(inout) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3)
  real(knd),intent(inout) :: V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3)
  real(knd),intent(inout) :: W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3)
  real(knd),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),intent(inout) :: temperature
#else
  real(knd),contiguous,intent(in)    :: xPr(-2:)
  real(knd),contiguous,intent(in)    :: xU(-3:)
  real(knd),contiguous,intent(inout),dimension(-2:,-2:,-2:) :: U,V,W
  real(knd),dimension(-1:,-1:,-1:),contiguous,intent(inout) :: temperature
#endif
  integer i,j,k,bufn
  real(knd) p,xe,xs,xb,DF
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

    if (enable_buoyancy) then
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
#ifdef __HMPP
#include "hmpp-include.f90"
#endif
    real(knd) DampF
    real(knd),intent(in)::x
    intrinsic exp

    if (x<=0) then
      DampF = 1
    elseif (x>=1) then
      DampF = 0
    else
     DampF = (1-0.04_knd*x**2) * ( 1 - (1-exp(10._knd*x**2)) / (1-exp(10._knd)) )
    endif
  endfunction Dampf


  !$hmpp <tsteps> NullInterior codelet
  subroutine NullInterior(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,&
                          nUnull,nVnull,nWnull,Unull,Vnull,Wnull,U,V,W)

  implicit none

  integer,intent(in)      :: Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,nUnull,nVnull,nWnull
#ifdef __HMPP
#include "hmpp-include.f90"
  integer,dimension(3,nUnull),intent(in)      :: Unull
  integer,dimension(3,nVnull),intent(in)      :: Vnull
  integer,dimension(3,nWnull),intent(in)      :: Wnull
  real(knd),dimension(-2:Unx+3,-2:Uny+3,-2:Unz+3),intent(inout) :: U
  real(knd),dimension(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),intent(inout) :: V
  real(knd),dimension(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),intent(inout) :: W
#else
  integer,dimension(:,:),contiguous,intent(in)              :: Unull,Vnull,Wnull
  real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout) :: U,V,W
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




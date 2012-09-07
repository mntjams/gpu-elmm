
!   !$hmpp <tsteps> KappaTemperature codelet
!   subroutine KappaTemperature_GPU(Prnx,Prny,Prnz,Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz&
!                             dxmin,dymin,dzmin,maxCNiter,epsCN,Re,&
!                             TBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr,&
!                             SCAL2,SCAL,U,V,W,sctype,coef,SubsidenceProfile,fluxProfile) !Kappa scheme with flux limiter
!   real(KND),intent(inout) ::Scal2(-1:Prnx+2,-1:Prny+2,-1:Prnz+2),Scal(-1:Prnx+2,-1:Prny+2,-1:Prnz+2) !Hunsdorfer et al. 1995, JCP
!   real(KND),intent(in) :: U(-2:Unx+3,-2:Uny+3,-2:Unz+3),V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3),W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3),coef
!   integer,intent(in) :: sctype
!   real(KND),intent(in) :: SubsidenceProfile(0:Prnz)
!   real(KND),intent(out) :: fluxProfile(0:Prnz)
!   integer i,j,k
!   real(KND) A,Ax,Ay,Az              !Auxiliary variables to store muliplication constants for efficiency
!   real(KND) vel,SL,SR,FLUX
!   real(KND),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2) :: SLOPE
!   real(KND),parameter ::eps = 1e-8
!   real(KND) :: FluxLimiter,r
!
!   FluxLimiter(r)=max(0._KND,min(2._KND*r,min(limparam,(1+2._KND*r)/3._KND)))
!
!   if (sctype==1) then
!     call BOUND_Temp(SCAL)
!   else
!     call BOUND_PASSSCALAR(SCAL)
!   endif
!
!   A = coef*dt
!   Ax = coef*dt/dxmin
!   Ay = coef*dt/dymin
!   Az = coef*dt/dzmin
!
!
!   !$omp parallel private(i,j,k,vel,SL,SR,FLUX)
!
!   !$omp workshare
!   SLOPE = 0
!   !$omp end workshare
!   !$omp do
!   do k = 1,Prnz
!    do j = 1,Prny
!     do i = 0,Prnx
!      if (U(i,j,k)>0) then
!       SR = (SCAL(i+1,j,k)-SCAL(i,j,k))
!       SL = (SCAL(i,j,k)-SCAL(i-1,j,k))
!      else
!       SR = (SCAL(i,j,k)-SCAL(i+1,j,k))
!       SL = (SCAL(i+1,j,k)-SCAL(i+2,j,k))
!      endif
!      SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
!     enddo
!    enddo
!   enddo
!   !$omp end do
!
!   !$omp do
!   do k = 1,Prnz
!    do j = 1,Prny
!     do i = 0,Prnx
!      if (U(i,j,k)>0) then
!       FLUX = U(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i-1,j,k))*SLOPE(i,j,k)/2._KND)
!      else
!       FLUX = U(i,j,k)*(SCAL(i+1,j,k)+(SCAL(i+1,j,k)-SCAL(i+2,j,k))*SLOPE(i,j,k)/2._KND)
!      endif
!      SCAL2(i,j,k) = SCAL2(i,j,k)-Ax*FLUX
!      SCAL2(i+1,j,k) = SCAL2(i+1,j,k)+Ax*FLUX
!     enddo
!    enddo
!   enddo
!   !$omp end do nowait
!
!
!   !$omp workshare
!   SLOPE = 0
!   !$omp end workshare
!   !$omp do
!   do k = 1,Prnz
!    do j = 0,Prny
!     do i = 1,Prnx
!      if (V(i,j,k)>0) then
!       SR = (SCAL(i,j+1,k)-SCAL(i,j,k))
!       SL = (SCAL(i,j,k)-SCAL(i,j-1,k))
!      else
!       SR = (SCAL(i,j,k)-SCAL(i,j+1,k))
!       SL = (SCAL(i,j+1,k)-SCAL(i,j+2,k))
!      endif
!      SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
!     enddo
!    enddo
!   enddo
!   !$omp end do
!
!
!   !$omp do
!   do k = 1,Prnz
!    do j = 0,Prny
!     do i = 1,Prnx
!      if (V(i,j,k)>0) then
!       FLUX = V(i,j,k)*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j-1,k))*SLOPE(i,j,k)/2._KND)
!      else
!       FLUX = V(i,j,k)*(SCAL(i,j+1,k)+(SCAL(i,j+1,k)-SCAL(i,j+2,k))*SLOPE(i,j,k)/2._KND)
!      endif
!
!      SCAL2(i,j,k) = SCAL2(i,j,k)-Ay*FLUX
!      SCAL2(i,j+1,k) = SCAL2(i,j+1,k)+Ay*FLUX
!     enddo
!    enddo
!   enddo
!   !$omp end do nowait
!
!
!   !$omp workshare
!   SLOPE = 0
!   !$omp end workshare
!   !$omp do
!   do k = 0,Prnz
!    do j = 1,Prny
!     do i = 1,Prnx
!      if (W(i,j,k)>0) then
!       SR = (SCAL(i,j,k+1)-SCAL(i,j,k))
!       SL = (SCAL(i,j,k)-SCAL(i,j,k-1))
!      else
!       SR = (SCAL(i,j,k)-SCAL(i,j,k+1))
!       SL = (SCAL(i,j,k+1)-SCAL(i,j,k+2))
!      endif
!      SLOPE(i,j,k) = FLUXLIMITER((SR+eps*sign(1._KND,SL))/(SL+eps*sign(1._KND,SL)))
!     enddo
!    enddo
!   enddo
!   !$omp end do nowait
!
!   !$omp workshare
!   fluxprofile = 0
!   !$omp end workshare
!
!   !$omp do reduction(+:fluxprofile)
!   do j = 1,Prny  !loop order due to avoid race condition
!    do k = 0,Prnz
!     do i = 1,Prnx
!
!      vel = W(i,j,k) - SubsidenceProfile(k)
!
!      if (vel>0) then
!       FLUX = vel*(SCAL(i,j,k)+(SCAL(i,j,k)-SCAL(i,j,k-1))*SLOPE(i,j,k)/2._KND)
!      else
!       FLUX = vel*(SCAL(i,j,k+1)+(SCAL(i,j,k+1)-SCAL(i,j,k+2))*SLOPE(i,j,k)/2._KND)
!      endif
!
!      if (vel>=1e-6) fluxprofile(k) = fluxprofile(k) + FLUX/vel*W(i,j,k)
!
!      SCAL2(i,j,k) = SCAL2(i,j,k) - Az*FLUX
!      SCAL2(i,j,k+1) = SCAL2(i,j,k+1) + Az*FLUX
!     enddo
!    enddo
!   enddo
!   !$omp end do
!
!   !$omp workshare
!   fluxprofile = fluxprofile / (Prnx*Prny)
!   !$omp end workshare
!
!   !$omp end parallel
!
!   endsubroutine KappaTemperature_GPU








  !$hmpp <tsteps> DiffTemperature codelet
  subroutine DiffTemperature_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,maxCNiter,epsCN,Re,&
                            TBtype,sideTemp,TempIn,BsideTArr,BsideTFLArr,&
                            TDiff,Temperature2,Temperature,&
                            coef,dt,l,res)
  implicit none
#include "hmpp-include.f90"

  integer,intent(in)    :: Prnx,Prny,Prnz,maxCNiter,TBtype(6)
  real(KND),intent(out),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2) :: Temperature2
  real(KND),intent(in),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)  :: Temperature

  real(KND),intent(in)  :: dxmin,dymin,dzmin,epsCN,Re,coef,dt
  real(KND),intent(in)  :: sideTemp(6),Tempin(-1:Prny+2,-1:Prnz+2),TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  real(KND),intent(in),dimension(-1:Prnx+2,-1:Prny+2) :: BsideTArr,BsideTFLArr
  integer,intent(out)   :: l
  real(KND),intent(out) :: res
  real(KND) Temperature3(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  integer nx,ny,nz,i,j,k,xi,yj,zk
  real(KND) p
  real(KND) A,Ax,Ay,Az,Ap(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)

  intrinsic max, abs, mod

   nx = Prnx
   ny = Prny
   nz = Prnz


  if (Re>0) then

  !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
   do k = 1,Prnz  !initital value using forward Euler
    do i = 1,Prnx
     do j = 1,Prny
      Temperature3(i,j,k) = (((TDiff(i+1,j,k)+TDiff(i,j,k))*(Temperature(i+1,j,k)-Temperature(i,j,k))/dxmin-&
        (TDiff(i,j,k)+TDiff(i-1,j,k))*(Temperature(i,j,k)-Temperature(i-1,j,k))/dxmin)/(dxmin)+&
       ((TDiff(i,j+1,k)+TDiff(i,j,k))*(Temperature(i,j+1,k)-Temperature(i,j,k))/dymin-&
        (TDiff(i,j,k)+TDiff(i,j-1,k))*(Temperature(i,j,k)-Temperature(i,j-1,k))/dymin)/(dymin)+&
       ((TDiff(i,j,k+1)+TDiff(i,j,k))*(Temperature(i,j,k+1)-Temperature(i,j,k))/dzmin-&
        (TDiff(i,j,k)+TDiff(i,j,k-1))*(Temperature(i,j,k)-Temperature(i,j,k-1))/dzmin)/(dzmin))
     enddo
    enddo
   enddo

   A = dt*coef
  !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
   do k = 1,Prnz
    do i = 1,Prnx
     do j = 1,Prny
       Temperature2(i,j,k) = Temperature(i,j,k)+A*Temperature3(i,j,k)
     enddo
    enddo
   enddo

   call BOUND_Temp_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,TBtype,sideTemp,BsideTArr,BsideTFLArr,TDiff,TempIn,Temperature2)



   Ax = 1._KND/(4._KND*dxmin**2)
   Ay = 1._KND/(4._KND*dymin**2)
   Az = 1._KND/(4._KND*dzmin**2)

  !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
   do k = 1,Prnz
    do i = 1,Prnx
     do j = 1,Prny
      Ap(i,j,k) = 1._KND/(1._KND/A+(((TDiff(i+1,j,k)+TDiff(i,j,k))+&
                           (TDiff(i,j,k)+TDiff(i-1,j,k)))*Ax+&
                           ((TDiff(i,j+1,k)+TDiff(i,j,k))+&
                           (TDiff(i,j,k)+TDiff(i,j-1,k)))*Ay+&
                           ((TDiff(i,j,k+1)+TDiff(i,j,k))+&
                           (TDiff(i,j,k)+TDiff(i,j,k-1)))*Az))
     enddo
    enddo
   enddo

   l=0
   res = epsCN + 1._KND
   do while (l<maxCNiter)!.and.res>epsCN
    l=l+1
!     res = 0

    call BOUND_Temp_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,TBtype,sideTemp,BsideTArr,BsideTFLArr,TDiff,TempIn,Temperature2)


   !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i), private(j,p)
   !, reduce(max:res)
    do k = 1,Prnz
     do i = 1,Prnx
      do j = 1+mod(i+k,2),Prny,2
        p = (Temperature(i,j,k)/A)+(Temperature3(i,j,k)/4._KND+&
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

   !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i), private(j,p)
   !, reduce(max:res)
    do k = 1,Prnz
     do i = 1,Prnx
      do j = 1+mod(i+k+1,2),Prny,2
        p = (Temperature(i,j,k)/A)+(Temperature3(i,j,k)/4._KND+&
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

   !$hmppcg grid blocksize 512x1
   !$hmppcg gridify(k,i)
   do k = 1,Prnz
    do i = 1,Prnx
     do j = 1,Prny
       Temperature2(i,j,k) = Temperature(i,j,k)
     enddo
    enddo
   enddo

   call BOUND_Temp_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,TBtype,sideTemp,BsideTArr,BsideTFLArr,TDiff,TempIn,Temperature2)

  endif
  endsubroutine DiffTemperature_GPU

!   !$hmpp <tsteps> DiffScalar codelet
  subroutine DIFFSCALAR_GPU(nscalars,Prnx,Prny,Prnz,dxmin,dymin,dzmin,maxCNiter,epsCN,Re,&
                            ScalBtype,sideScal,&
                            TDiff,SCAL2,SCAL,&
                            coef,dt,l,res)
  implicit none
#include "hmpp-include.f90"

  integer,intent(in)    :: nscalars,Prnx,Prny,Prnz,maxCNiter,ScalBtype(6)
  real(KND),intent(out),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,nscalars) :: Scal2
  real(KND),intent(in),dimension(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,nscalars)  :: Scal

  real(KND),intent(in)  :: dxmin,dymin,dzmin,epsCN,Re,coef,dt
  real(KND),intent(in)  :: sideScal(6),TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)
  integer,intent(out)   :: l
  real(KND),intent(out) :: res
  real(KND) Scal3(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,nscalars)
  integer nx,ny,nz,i,j,k,m,xi,yj,zk
  real(KND) p
  real(KND) A,Ax,Ay,Az,Ap(-1:Prnx+2,-1:Prny+2,-1:Prnz+2)

  intrinsic max, abs, mod

  nx = Prnx
  ny = Prny
  nz = Prnz

  Ax = 1._KND/(4._KND*dxmin**2)
  Ay = 1._KND/(4._KND*dymin**2)
  Az = 1._KND/(4._KND*dzmin**2)

  !$hmppcg grid blocksize 512x1
  !$hmppcg gridify(k,i)
  do k = 1,Prnz
   do i = 1,Prnx
    do j = 1,Prny
      Ap(i,j,k) = 1._KND/(1._KND/A+(((TDiff(i+1,j,k)+TDiff(i,j,k))+&
                           (TDiff(i,j,k)+TDiff(i-1,j,k)))*Ax+&
                           ((TDiff(i,j+1,k)+TDiff(i,j,k))+&
                           (TDiff(i,j,k)+TDiff(i,j-1,k)))*Ay+&
                           ((TDiff(i,j,k+1)+TDiff(i,j,k))+&
                           (TDiff(i,j,k)+TDiff(i,j,k-1)))*Az))
    enddo
   enddo
  enddo

  do m = 1,nscalars

    if (Re>0) then

    !$hmppcg grid blocksize 512x1
     !$hmppcg gridify(k,i)
     do k = 1,Prnz  !initital value using forward Euler
      do i = 1,Prnx
       do j = 1,Prny
        SCAL3(i,j,k,l) = (((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL(i+1,j,k,l)-SCAL(i,j,k,l))/dxmin-&
          (TDiff(i,j,k)+TDiff(i-1,j,k))*(SCAL(i,j,k,l)-SCAL(i-1,j,k,l))/dxmin)/(dxmin)+&
         ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL(i,j+1,k,l)-SCAL(i,j,k,l))/dymin-&
          (TDiff(i,j,k)+TDiff(i,j-1,k))*(SCAL(i,j,k,l)-SCAL(i,j-1,k,l))/dymin)/(dymin)+&
         ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL(i,j,k+1,l)-SCAL(i,j,k,l))/dzmin-&
          (TDiff(i,j,k)+TDiff(i,j,k-1))*(SCAL(i,j,k,l)-SCAL(i,j,k-1,l))/dzmin)/(dzmin))
       enddo
      enddo
     enddo

     A = dt*coef
    !$hmppcg grid blocksize 512x1
     !$hmppcg gridify(k,i)
     do k = 1,Prnz
      do i = 1,Prnx
       do j = 1,Prny
         SCAL2(i,j,k,l) = SCAL(i,j,k,l)+A*SCAL3(i,j,k,l)
       enddo
      enddo
     enddo


     call BOUND_PASSSCALAR_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,ScalBtype,sideScal,TDiff,SCAL2(:,:,:,l))





     do l = 1,maxCNiter
      res = 0

      call BOUND_PASSSCALAR_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,ScalBtype,sideScal,TDiff,SCAL2(:,:,:,l))


     !$hmppcg grid blocksize 512x1
     !$hmppcg gridify(k,i), reduce(max:res)
      do k = 1,Prnz
       do i = 1,Prnx
        do j = 1+mod(i+k,2),Prny,2
          p = (SCAL(i,j,k,l)/A)+(SCAL3(i,j,k,l)/4._KND+&
           ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL2(i+1,j,k,l))-&
            (TDiff(i,j,k)+TDiff(i-1,j,k))*(-SCAL2(i-1,j,k,l)))*Ax+&
           ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL2(i,j+1,k,l))-&
            (TDiff(i,j,k)+TDiff(i,j-1,k))*(-SCAL2(i,j-1,k,l)))*Ay+&
           ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL2(i,j,k+1,l))-&
            (TDiff(i,j,k)+TDiff(i,j,k-1))*(-SCAL2(i,j,k-1,l)))*Az&
           )
           p = p*Ap(i,j,k)
           res = max(res,abs(p-SCAL2(i,j,k,l)))
           SCAL2(i,j,k,l) = p
        enddo
       enddo
      enddo

     !$hmppcg grid blocksize 512x1
     !$hmppcg gridify(k,i), reduce(max:res)
      do k = 1,Prnz
       do i = 1,Prnx
        do j = 1+mod(i+k+1,2),Prny,2
          p = (SCAL(i,j,k,l)/A)+(SCAL3(i,j,k,l)/4._KND+&
           ((TDiff(i+1,j,k)+TDiff(i,j,k))*(SCAL2(i+1,j,k,l))-&
            (TDiff(i,j,k)+TDiff(i-1,j,k))*(-SCAL2(i-1,j,k,l)))*Ax+&
           ((TDiff(i,j+1,k)+TDiff(i,j,k))*(SCAL2(i,j+1,k,l))-&
            (TDiff(i,j,k)+TDiff(i,j-1,k))*(-SCAL2(i,j-1,k,l)))*Ay+&
           ((TDiff(i,j,k+1)+TDiff(i,j,k))*(SCAL2(i,j,k+1,l))-&
            (TDiff(i,j,k)+TDiff(i,j,k-1))*(-SCAL2(i,j,k-1,l)))*Az&
           )
           p = p*Ap(i,j,k)
           res = max(res,abs(p-SCAL2(i,j,k,l)))
           SCAL2(i,j,k,l) = p
        enddo
       enddo
      enddo


      if (res<=epsCN) exit
     enddo


    else

     !$hmppcg grid blocksize 512x1
     !$hmppcg gridify(k,i)
     do k = 1,Prnz
      do i = 1,Prnx
       do j = 1,Prny
         SCAL2(i,j,k,l) = SCAL(i,j,k,l)
       enddo
      enddo
     enddo


      call BOUND_PASSSCALAR_GPU(Prnx,Prny,Prnz,dxmin,dymin,dzmin,Re,ScalBtype,sideScal,TDiff,SCAL2(:,:,:,l))

    endif

  enddo
  endsubroutine DIFFSCALAR_GPU


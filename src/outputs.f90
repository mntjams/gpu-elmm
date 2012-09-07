module BigEndian
!   use iso_fortran_env
  implicit none

  private
  public :: GetEndianess, BigEnd

  integer,parameter :: int8=selected_int_kind(1),int32=selected_int_kind(9),int64=selected_int_kind(10)
  integer,parameter :: real32=selected_real_kind(p=6,r=37),real64=selected_real_kind(p=15,r=200)
  logical,save :: littleendian

  interface BigEnd
    module procedure BigEnd32
    module procedure BigEnd64
  end interface

  contains

    subroutine GetEndianess
      integer(int8),dimension(4):: bytes !may not work on some processors

      bytes=transfer(1_int32,bytes,4)
      if (bytes(4)==1) then
        littleendian=.false.
      else
        littleendian=.true.
      endif
    end subroutine GetEndianess

    elemental function BigEnd32(x) result(res)
      real(real32) :: res
      real(real32),intent(in)::x
      integer(int8),dimension(4):: bytes !may not work on some processors

      if (.not.littleendian) then
        res = x
      else
        bytes = transfer(x,bytes,4)
        res = transfer(bytes(4:1:-1),res)
      endif
    end function BigEnd32

    elemental function BigEnd64(x) result(res)
      real(real64) :: res
      real(real64),intent(in)::x
      integer(int8),dimension(8):: bytes !may not work on some processors

      if (.not.littleendian) then
        res = x
      else
        bytes = transfer(x,bytes,8)
        res = transfer(bytes(8:1:-1),res)
      endif
    end function BigEnd64

end module BigEndian


module OUTPUTS
  use PARAMETERS
  use BOUNDARIES
  use SCALARS
  use WALLMODELS, only: GroundDeposition, GroundUstar
  use GEOMETRIC
  use TURBINLET, only: ustarinlet
  use BigEndian

  implicit none


  private
  public OutTStep, Output, store, display, probes, NumProbes, AllocateOutputs,&
         Tprobe, TOutputSwitches, TDisplaySwitches, SetFrameDomain, proftempfl

  real(KND),dimension(:),allocatable :: profuavg,profuavg2,profvavg,profvavg2,profuuavg,profvvavg,profwwavg,&
                                         profU,profV,profuu,profvv,profww,proftauavg,proftau,proftausgs,proftausgsavg,&
                                         proftemp,proftempfl,proftempavg,proftempavg2,proftempflavg,&
                                         proftempflsgs,proftempflsgsavg,proftt,profttavg,&
                                         profuw,profuwavg,profuwsgs,profuwsgsavg,&
                                         profvw,profvwavg,profvwsgs,profvwsgsavg

  real(KND),dimension(:,:),allocatable ::profscal,profscalfl,profscalavg,profscalavg2,profscalflavg,&  !which scalar, height
                                         profscalflsgs,profscalflsgsavg,profss,profssavg

  real(KND),allocatable :: Uavg(:,:,:),Vavg(:,:,:),Wavg(:,:,:),Pravg(:,:,:),ScalarAvg(:,:,:,:)

  real(TIM),allocatable,dimension(:) :: times                                !times of the timesteps

  real(KND),allocatable,dimension(:) :: CDtime,CLtime,deltime,tke,dissip

  real(KND),allocatable,dimension(:,:) :: ustar,tstar                        !first index differentiates flux from friction number
                                                                             !second index is time

  real(KND),allocatable,dimension(:,:) :: Utime,Vtime,Wtime,Prtime,temptime  !position, time

  real(KND),allocatable,dimension(:,:,:) :: scalptime                        !which scalar, position, time
  real(KND),allocatable,dimension(:,:) :: scalsumtime                        !which scalar, time

  integer :: NumProbes                        !number of probes in space to collect timed data

  type TProbe
    integer :: Ui,Uj,Uk,Vi,Vj,Vk,Wi,Wj,Wk    !grid coordinates of probes in the U,V,W grids
    integer :: i,j,k                         !grid coordinates of probes in the scalar grid
    real(KND) :: x,y,z                       !physical coordinates of probes
  end type TProbe

  type(TProbe),allocatable,dimension(:),save :: probes

  type TFrameDomain
    integer   :: dimension,direction
    real(KND) :: position
    integer   :: mini, maxi, minj, maxj, mink, maxk
  endtype TFrameDomain



  type TOutputSwitches
    integer :: U = 0
    integer :: U_interp = 0
    integer :: V = 0
    integer :: V_interp = 0
    integer :: W = 0
    integer :: W_interp = 0

    integer :: out = 1
    integer :: avg = 1

    integer :: scalars = 1
    integer :: scalarsavg = 1

    integer :: deposition = 1


    integer :: out_U = 1
    integer :: out_vort = 1
    integer :: out_Pr = 1
    integer :: out_Prtype = 1
    integer :: out_lambda2 = 1
    integer :: out_T = 1
    integer :: out_div = 0
    integer :: out_visc = 0

    integer :: avg_U = 1
    integer :: avg_vort = 0
    integer :: avg_Pr = 1
    integer :: avg_Prtype = 0
    integer :: avg_T = 1

    integer :: frame_U = 1
    integer :: frame_vort = 0
    integer :: frame_Pr = 0
    integer :: frame_lambda2 = 0
    integer :: frame_scalars = 0
    integer :: frame_sumscalars = 1
    integer :: frame_T = 1
    integer :: frame_tempfl = 0
    integer :: frame_scalfl = 0
    type(TFrameDomain),dimension(:),allocatable :: frame_domains

    integer :: deltime = 0
    integer :: tke = 0
    integer :: dissip = 0
    integer :: scalsumtime = 0
    integer :: scaltotsumtime = 0
    integer :: ustar = 1
    integer :: tstar = 1

    integer :: blprofiles = 0
  end type TOutputSwitches

  type(TOutputSwitches),save :: store

  type TDisplaySwitches
    integer :: delta = 1
    integer :: ustar = 0
    integer :: tstar = 0
  endtype

  type(TDisplaySwitches),save :: display



contains


 subroutine AllocateOutputs
 integer :: k

   if (computescalars>0) then
    if (averaging==1.and.store%scalarsavg>0) then
     allocate(Scalaravg(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
     Scalaravg = 0
    endif
   else
    if (averaging==1) then
     allocate(Scalaravg(0,0,0,0))
    endif
   endif


  if (averaging==1) then
   if (store%avg_U>0) then
     allocate(Uavg(-2:Unx+3,-2:Uny+3,-2:Unz+3))
     allocate(Vavg(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
     allocate(Wavg(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
     Uavg = 0
     Vavg = 0
     Wavg = 0
   endif

   if (store%avg_Pr>0) then
     allocate(Pravg(1:Prnx,1:Prny,1:Prnz))
     Pravg = 0
   endif

   if (buoyancy==1.and.store%avg_T>0) then
    allocate(temperatureavg(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    temperatureavg = 0
   endif
  endif

  allocate(times(0:nt))

  if (NumProbes>0) then
    allocate(Utime(NumProbes,0:nt),Vtime(NumProbes,0:nt),Wtime(NumProbes,0:nt),Prtime(NumProbes,0:nt))
    times = huge(1.0_KND)
    Utime = huge(1.0_KND)
    Vtime = huge(1.0_KND)
    Wtime = huge(1.0_KND)
    Prtime = huge(1.0_KND)
  endif

  if (NumProbes>0.and.buoyancy==1) then
   allocate(temptime(NumProbes,0:nt))
   temptime = huge(1.0_KND)
  endif

  if (store%deltime>0) then
    allocate(deltime(0:nt))
    deltime = huge(1.0_KND)
  endif

  if (store%tke>0) then
    allocate(tke(0:nt))
    tke = huge(1.0_KND)
  endif

  if (store%deltime>0) then
    allocate(dissip(0:nt))
    dissip = huge(1.0_KND)
    dissip(0)=0
  endif

  if (wallmodeltype>0.and.(display%ustar>0.or.store%ustar>0)) then
    allocate(ustar(2,0:nt))
    ustar = huge(1.0)
  endif

  if (wallmodeltype>0.and.buoyancy==1.and.TBtype(Bo)==DIRICHLET.and.(display%tstar>0.or.store%tstar>0)) then
    allocate(tstar(2,0:nt))
    tstar = huge(1.0)
  endif

  if (computescalars>0) then
    if (store%scalsumtime>0) then
      allocate(scalsumtime(1:computescalars,0:nt))
      scalsumtime = huge(1.0_KND)
    endif

    if (NumProbes>0) then
      allocate(scalptime(1:computescalars,1:NumProbes,0:nt))
      scalptime = huge(1.0_KND)
   endif
  endif

  do k = 1,NumProbes
    call GridCoords(probes(k)%i,probes(k)%j,probes(k)%k,probes(k)%x,probes(k)%y,probes(k)%z)

    probes(k)%i = max(probes(k)%i,1)
    probes(k)%j = max(probes(k)%j,1)
    probes(k)%k = max(probes(k)%k,1)
    probes(k)%i = min(probes(k)%i,Prnx)
    probes(k)%j = min(probes(k)%j,Prny)
    probes(k)%k = min(probes(k)%k,Prnz)

    call GridCoordsU(probes(k)%Ui,probes(k)%Uj,probes(k)%Uk,probes(k)%x,probes(k)%y,probes(k)%z)
    call GridCoordsV(probes(k)%Vi,probes(k)%Vj,probes(k)%Vk,probes(k)%x,probes(k)%y,probes(k)%z)
    call GridCoordsW(probes(k)%Wi,probes(k)%Wj,probes(k)%Wk,probes(k)%x,probes(k)%y,probes(k)%z)
  enddo

  if (store%BLprofiles>0.and.averaging==1) then
    allocate(profU(0:Unz+1),profV(0:Vnz+1),profUavg(0:Unz+1),profVavg(0:Vnz+1),profUavg2(0:Unz+1),profVavg2(0:Vnz+1))
    allocate(profuuavg(1:Unz),profvvavg(1:Vnz),profwwavg(0:Wnz),profuu(1:Unz),profvv(1:Vnz),profww(0:Wnz))
    allocate(profuw(0:Prnz),profuwavg(0:Prnz),profuwsgs(0:Prnz),profuwsgsavg(0:Prnz))
    allocate(profvw(0:Prnz),profvwavg(0:Prnz),profvwsgs(0:Prnz),profvwsgsavg(0:Prnz))

    allocate(proftemp(1:Prnz),proftempfl(0:Prnz),proftempavg(1:Prnz),proftempavg2(1:Prnz),proftempflavg(0:Prnz))
    allocate(proftempflsgs(0:Prnz),proftempflsgsavg(0:Prnz),proftt(1:Prnz),profttavg(1:Prnz))

    allocate(profscal(computescalars,1:Prnz),profscalfl(computescalars,0:Prnz),profscalavg(computescalars,1:Prnz))
    allocate(profscalavg2(computescalars,1:Prnz),profscalflavg(computescalars,0:Prnz))
    allocate(profscalflsgs(computescalars,0:Prnz),profscalflsgsavg(computescalars,0:Prnz))
    allocate(profss(computescalars,1:Prnz),profssavg(computescalars,1:Prnz))

    profU = 0
    profV = 0
    profUavg = 0
    profVavg = 0
    profUavg2 = 0
    profVavg2 = 0
    profuu = 0
    profvv = 0
    profww = 0
    profuuavg = 0
    profvvavg = 0
    profwwavg = 0
    profuw = 0
    profvw = 0
    profuwavg = 0
    profvwavg = 0

    proftemp = 0
    proftempavg = 0
    proftempavg2 = 0
    proftempfl = 0
    proftempflavg = 0
    proftempflsgs = 0
    proftempflsgsavg = 0
    proftt = 0
    profttavg = 0

    if (computescalars>0) then
      profscal = 0
      profscalavg = 0
      profscalavg2 = 0
      profscalfl = 0
      profscalflavg = 0
      profscalflsgs = 0
      profscalflsgsavg = 0
      profss = 0
      profssavg = 0
    endif
  endif

  call GetEndianess

 end subroutine AllocateOutputs













 subroutine OutTStep(U,V,W,Pr,delta)
 real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
 real(KND),dimension(1:,1:,1:),intent(in) :: Pr
 real(KND),intent(in) :: delta

 integer :: l,i,j,k,nWM
 real(KND) :: S,S2
 integer, save :: fnum = 0

   times(step)=time

   do l = 1,computescalars
      S = 0
      do k = 1,Prnz
       do j = 1,Prny
        do i = 1,Prnx
         if (Prtype(i,j,k)==0) S = S+scalar(i,j,k,l)*dxPr(i)*dyPr(j)*dzPr(k)
        enddo
       enddo
      enddo
      scalsumtime(l,step)=S
   enddo



   do k = 1,NumProbes
     Utime(k,step)=Trilinint((probes(k)%x-xU(probes(k)%Ui))/(xU(probes(k)%Ui+1)-xU(probes(k)%Ui)),&
                          (probes(k)%y-yPr(probes(k)%Uj))/(yPr(probes(k)%Uj+1)-yPr(probes(k)%Uj)),&
                          (probes(k)%z-zPr(probes(k)%Uk))/(zPr(probes(k)%Uk+1)-zPr(probes(k)%Uk)),&
                          U(probes(k)%Ui,probes(k)%Uj,probes(k)%Uk),U(probes(k)%Ui+1,probes(k)%Uj,probes(k)%Uk),&
                          U(probes(k)%Ui,probes(k)%Uj+1,probes(k)%Uk),U(probes(k)%Ui,probes(k)%Uj,probes(k)%Uk+1),&
                          U(probes(k)%Ui+1,probes(k)%Uj+1,probes(k)%Uk),U(probes(k)%Ui+1,probes(k)%Uj,probes(k)%Uk+1),&
                          U(probes(k)%Ui,probes(k)%Uj+1,probes(k)%Uk+1),U(probes(k)%Ui+1,probes(k)%Uj+1,probes(k)%Uk+1))
     Vtime(k,step)=Trilinint((probes(k)%x-xPr(probes(k)%Vi))/(xPr(probes(k)%Vi+1)-xPr(probes(k)%Vi)),&
                          (probes(k)%y-yV(probes(k)%Vj))/(yV(probes(k)%Vj+1)-yV(probes(k)%Vj)),&
                          (probes(k)%z-zPr(probes(k)%Vk))/(zPr(probes(k)%Vk+1)-zPr(probes(k)%Vk)),&
                          V(probes(k)%Vi,probes(k)%Vj,probes(k)%Vk),V(probes(k)%Vi+1,probes(k)%Vj,probes(k)%Vk),&
                          V(probes(k)%Vi,probes(k)%Vj+1,probes(k)%Vk),V(probes(k)%Vi,probes(k)%Vj,probes(k)%Vk+1),&
                          V(probes(k)%Vi+1,probes(k)%Vj+1,probes(k)%Vk),V(probes(k)%Vi+1,probes(k)%Vj,probes(k)%Vk+1),&
                          V(probes(k)%Vi,probes(k)%Vj+1,probes(k)%Vk+1),V(probes(k)%Vi+1,probes(k)%Vj+1,probes(k)%Vk+1))
     Wtime(k,step)=Trilinint((probes(k)%x-xPr(probes(k)%Wi))/(xPr(probes(k)%Wi+1)-xPr(probes(k)%Wi)),&
                          (probes(k)%y-yPr(probes(k)%Wj))/(yPr(probes(k)%Wj+1)-yPr(probes(k)%Wj)),&
                          (probes(k)%z-zW(probes(k)%Wk))/(zW(probes(k)%Wk+1)-zW(probes(k)%Wk)),&
                          W(probes(k)%Wi,probes(k)%Wj,probes(k)%Wk),W(probes(k)%Wi+1,probes(k)%Wj,probes(k)%Wk),&
                          W(probes(k)%Wi,probes(k)%Wj+1,probes(k)%Wk),W(probes(k)%Wi,probes(k)%Wj,probes(k)%Wk+1),&
                          W(probes(k)%Wi+1,probes(k)%Wj+1,probes(k)%Wk),W(probes(k)%Wi+1,probes(k)%Wj,probes(k)%Wk+1),&
                          W(probes(k)%Wi,probes(k)%Wj+1,probes(k)%Wk+1),W(probes(k)%Wi+1,probes(k)%Wj+1,probes(k)%Wk+1))
     Prtime(k,step)=Trilinint((probes(k)%x-xPr(probes(k)%i))/(xPr(probes(k)%i+1)-xPr(probes(k)%i)),&
                          (probes(k)%y-yPr(probes(k)%j))/(yPr(probes(k)%j+1)-yPr(probes(k)%j)),&
                          (probes(k)%z-zPr(probes(k)%k))/(zPr(probes(k)%k+1)-zPr(probes(k)%k)),&
                           Pr(probes(k)%i,probes(k)%j,probes(k)%k),&
                           Pr(min(probes(k)%i+1,Unx+1),probes(k)%j,probes(k)%k),&
                           Pr(probes(k)%i,min(probes(k)%j+1,Vny+1),probes(k)%k),&
                           Pr(probes(k)%i,probes(k)%j,min(probes(k)%k+1,Wnz+1)),&
                           Pr(min(probes(k)%i+1,Unx+1),min(probes(k)%j+1,Vny+1),probes(k)%k),&
                           Pr(min(probes(k)%i+1,Unx+1),probes(k)%j,min(probes(k)%k+1,Wnz+1)),&
                           Pr(probes(k)%i,min(probes(k)%j+1,Vny+1),min(probes(k)%k+1,Wnz+1)),&
                           Pr(min(probes(k)%i+1,Unx+1),min(probes(k)%j+1,Vny+1),min(probes(k)%k+1,Wnz+1)))

     if (buoyancy>0) then
       temptime(k,step)=Trilinint((probes(k)%x-xPr(probes(k)%i))/(xPr(probes(k)%i+1)-xPr(probes(k)%i)),&
                          (probes(k)%y-yPr(probes(k)%j))/(yPr(probes(k)%j+1)-yPr(probes(k)%j)),&
                          (probes(k)%z-zPr(probes(k)%k))/(zPr(probes(k)%k+1)-zPr(probes(k)%k)),&
                          temperature(probes(k)%i,probes(k)%j,probes(k)%k),&
                          temperature(probes(k)%i+1,probes(k)%j,probes(k)%k),&
                          temperature(probes(k)%i,probes(k)%j+1,probes(k)%k),&
                          temperature(probes(k)%i,probes(k)%j,probes(k)%k+1),&
                          temperature(probes(k)%i+1,probes(k)%j+1,probes(k)%k),&
                          temperature(probes(k)%i+1,probes(k)%j,probes(k)%k+1),&
                          temperature(probes(k)%i,probes(k)%j+1,probes(k)%k+1),&
                          temperature(probes(k)%i+1,probes(k)%j+1,probes(k)%k+1))
     endif

     do l = 1,computescalars
       scalptime(l,k,step)=Trilinint((probes(k)%x-xPr(probes(k)%i))/(xPr(probes(k)%i+1)-xPr(probes(k)%i)),&
                          (probes(k)%y-yPr(probes(k)%j))/(yPr(probes(k)%j+1)-yPr(probes(k)%j)),&
                          (probes(k)%z-zPr(probes(k)%k))/(zPr(probes(k)%k+1)-zPr(probes(k)%k)),&
                          Scalar(probes(k)%i,probes(k)%j,probes(k)%k,l),&
                          Scalar(probes(k)%i+1,probes(k)%j,probes(k)%k,l),&
                          Scalar(probes(k)%i,probes(k)%j+1,probes(k)%k,l),&
                          Scalar(probes(k)%i,probes(k)%j,probes(k)%k+1,l),&
                          Scalar(probes(k)%i+1,probes(k)%j+1,probes(k)%k,l),&
                          Scalar(probes(k)%i+1,probes(k)%j,probes(k)%k+1,l),&
                          Scalar(probes(k)%i,probes(k)%j+1,probes(k)%k+1,l),&
                          Scalar(probes(k)%i+1,probes(k)%j+1,probes(k)%k+1,l))
     enddo
   enddo

   if (store%tke>0) then
     tke(step)=totke(U,V,W)
   endif

   if (store%tke>0.and.store%dissip>0.and.step>0) then
     dissip(step)=(tke(step-1)-tke(step))/(times(step)-times(step-1))
   endif

   if (store%deltime>0) then
     deltime(step)=delta/dt
   endif

   endstep = step

   if (display%delta>0) then
     write (*,*) "delta: ",delta
   endif

   if (frames>0)then
      if ((time>=timefram1).and.(time<=timefram2+(timefram2-timefram1)/(frames-1))&
        .and.(time>=timefram1+fnum*(timefram2-timefram1)/(frames-1))) then
       fnum = fnum+1
       call FRAME(U,V,W,Pr,fnum)
      endif
!    call OUTINLET(U,V,W)
   endif


   if ((averaging==1).and.((time>=timeavg1).and.(time<=timeavg2))) then

     if (store%avg_U>0) then
       Uavg = Uavg+U*dt/(timeavg2-timeavg1)
       Vavg = Vavg+V*dt/(timeavg2-timeavg1)
       Wavg = Wavg+W*dt/(timeavg2-timeavg1)
     endif

     if (store%avg_Pr>0) then
       Pravg = Pravg+Pr(1:Prnx,1:Prny,1:Prnz)*dt/(timeavg2-timeavg1)
     endif

     if (buoyancy==1.and.store%avg_T>0) then
       temperatureavg = temperatureavg+temperature*dt/(timeavg2-timeavg1)
     endif

     if (computescalars>0.and.store%scalarsavg>0) then
       SCALARavg = SCALARavg+SCALAR*dt/(timeavg2-timeavg1)
     endif

  endif


  if (wallmodeltype>0.and.(display%ustar>0.or.store%ustar>0)) then

    S = GroundUstar()

    S2 = S*Re

    if (display%ustar>0) then
      if (allocated(ustarinlet)) then
       write(*,*) "ustar:",S,"Re_tau:",S2,"u*inlet",ustarinlet(1)
      else
       write(*,*) "ustar:",S,"Re_tau:",S2
      endif
    endif
    if (store%ustar>0) then
      ustar(:,step)=(/ S2 , S /)
    endif
  endif


   if (wallmodeltype>0.and.buoyancy==1.and.TBtype(Bo)==DIRICHLET.and.(display%tstar>0.or.store%tstar>0)) then
    S2 = SUM(BsideTFLArr(1:Prnx,1:Prny))/(Prnx*Prny)
    S=-S*S2

    if (display%tstar>0) then
      write(*,*) "Tstar",S,"tflux", S2
    endif

    if (store%tstar>0) then
      tstar(:,step) = (/ S2,S /)
    endif
   endif


   if (store%BLprofiles>0) then
     if ((averaging==1).and.((time>=timeavg1).and.(time<=timeavg2))) then

      if (buoyancy>0.and.computescalars>0) then
        call BLProfiles(U,V,W,temperature,scalar)
      elseif (buoyancy>0) then
        call BLProfiles(U,V,W,temperature = temperature)
      elseif (computescalars>0) then
        call BLProfiles(U,V,W,scalar = scalar)
      else
        call BLProfiles(U,V,W)
      endif

      profuavg = profuavg+profu*dt/(timeavg2-timeavg1)
      profvavg = profvavg+profv*dt/(timeavg2-timeavg1)
      profuwavg = profuwavg+profuw*dt/(timeavg2-timeavg1)
      profuwsgsavg = profuwsgsavg+profuwsgs*dt/(timeavg2-timeavg1)
      profvwavg = profvwavg+profvw*dt/(timeavg2-timeavg1)
      profvwsgsavg = profvwsgsavg+profvwsgs*dt/(timeavg2-timeavg1)
      profuuavg = profuuavg+profuu*dt/(timeavg2-timeavg1)
      profvvavg = profvvavg+profvv*dt/(timeavg2-timeavg1)
      profwwavg = profwwavg+profww*dt/(timeavg2-timeavg1)

      if (buoyancy>0) then
        proftempavg = proftempavg+proftemp*dt/(timeavg2-timeavg1)
        proftempflavg = proftempflavg+proftempfl*dt/(timeavg2-timeavg1)
        proftempflsgsavg = proftempflsgsavg+proftempflsgs*dt/(timeavg2-timeavg1)
        profttavg = profttavg+proftt*dt/(timeavg2-timeavg1)
      endif

      if (computescalars>0) then
        profscalavg = profscalavg+profscal*dt/(timeavg2-timeavg1)
        profscalflavg = profscalflavg+profscalfl*dt/(timeavg2-timeavg1)
        profscalflsgsavg = profscalflsgsavg+profscalflsgs*dt/(timeavg2-timeavg1)
        profssavg = profssavg+profss*dt/(timeavg2-timeavg1)
      endif
    endif
 endif

 end subroutine OutTstep





 subroutine OutputProfiles(U,V,W,Pr)
 real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
 real(KND),dimension(1:,1:,1:),intent(in) :: Pr
 character(70) :: str
 character(2) :: prob
 real(KND) :: S,S2,nom,denom
 integer :: i,j,k

  do k = 1,NumProbes

    write(prob,"(i2.2)") k

    open(11,file="Prtimep"//prob//".txt")
    do j = 0,endstep
     write (11,*) times(j),Prtime(k,j)
    enddo
    close(11)

    open(11,file="Utimep"//prob//".txt")
    do j = 0,endstep
     write (11,*) times(j),Utime(k,j)
    enddo
    close(11)

    open(11,file="Vtimep"//prob//".txt")
    do j = 0,endstep
     write (11,*) times(j),Vtime(k,j)
    enddo
    close(11)

    open(11,file="Wtimep"//prob//".txt")
    do j = 0,endstep
     write (11,*) times(j),Wtime(k,j)
    enddo
    close(11)

    if (buoyancy==1) then
      open(11,file="temptimep"//prob//".txt")
      do j = 0,endstep
       write (11,*) times(j),temptime(k,j)
      enddo
      close(11)
    endif

    if (computescalars>0) then
      open(11,file="scaltimep"//prob//".txt")
      do j = 1,endstep
       write (11,*) times(j),scalptime(:,k,j)
      enddo
      close(11)
    endif

  enddo

  if (store%deltime>0) then
    open(11,file="deltime.txt")
    do j = 1,endstep
     write (11,*) times(j),deltime(j)
    enddo
    close(11)
  endif

  if (store%tke>0) then
    open(11,file="tke.txt")
    do j = 0,endstep
     write (11,*) times(j),tke(j)
    enddo
    close(11)
  endif

  if (store%tke>0.and.store%dissip>0) then
   open(11,file="dissip.txt")
    do j = 1,endstep
     write (11,*) times(j),dissip(j)
    enddo
    close(11)
  endif

  if (wallmodeltype>0.and.display%ustar>0) then
    open(11,file="Retau.txt")
    do j = 0,endstep
     write (11,*) times(j),ustar(:,j)
    enddo
    close(11)
  endif

  if (wallmodeltype>0.and.buoyancy==1.and.TBtype(Bo)==DIRICHLET.and.store%tstar>0) then
    open(11,file="tflux.txt")
    do j = 0,endstep
     write (11,*) times(j),tstar(:,j)
    enddo
    close(11)
  endif


  if (computescalars>0.and.store%scalsumtime>0) then
    open(11,file="scalsumtime.txt")
    do j = 1,endstep
     write (11,*) times(j),scalsumtime(:,j)
    enddo
    close(11)
  endif

  if (computescalars>0.and.store%scaltotsumtime>0) then
    open(11,file="scaltotsumtime.txt")
    do j = 1,endstep
     write (11,*) times(j),sum(scalsumtime(:,j))
    enddo
    close(11)
  endif


  if (store%BLprofiles>0.and.averaging==1) then

     open(11,file="profu.txt")
     do k = 1,Unz
      write (11,*) zPr(k),profuavg(k)
     enddo
     close(11)

     open(11,file="profv.txt")
     do k = 1,Vnz
      write (11,*) zPr(k),profvavg(k)
     enddo
     close(11)

     open(11,file="profuu.txt")
     do k = 1,Unz
      write (11,*) zPr(k),profuuavg(k)
     enddo
     close(11)

     open(11,file="profvv.txt")
     do k = 1,Vnz
      write (11,*) zPr(k),profvvavg(k)
     enddo
     close(11)

     open(11,file="profww.txt")
     do k = 1,Wnz
      write (11,*) zW(k),profwwavg(k)
     enddo
     close(11)

     open(11,file="profuw.txt")
     do k = 0,Prnz
      write (11,*) zW(k),profuwavg(k),profuwsgsavg(k)
     enddo
     close(11)

     open(11,file="profvw.txt")
     do k = 0,Prnz
      write (11,*) zW(k),profvwavg(k),profvwsgsavg(k)
     enddo
     close(11)

     if (buoyancy>0) then

        open(11,file="proftemp.txt")
        do k = 1,Prnz
         write (11,*) zPr(k),proftempavg(k)
        enddo
        close(11)

        open(11,file="proftempfl.txt")
        do k = 0,Prnz
         write (11,*) zW(k),proftempflavg(k),proftempflsgsavg(k)
        enddo
        close(11)

        open(11,file="proftt.txt")
        do k = 1,Prnz
         write (11,*) zPr(k),profttavg(k)
        enddo
        close(11)

        if (allocated(Uavg).and.allocated(Vavg).and.allocated(temperatureavg)) then
          open(11,file="profRig.txt")
          do k = 1,Prnz
           S = 0
           do j = 1,Prny
            do i = 1,Prnx
             S = S+Rig(i,j,k,Uavg,Vavg,temperatureavg)
            enddo
           enddo
           S = S/(Prnx*Prny)
           write (11,*) zPr(k),S
          enddo
          close(11)

          open(11,file="profRf.txt")
          do k = 1,Prnz
           S = 0
           S2 = 0
           do j = 1,Prny
            do i = 1,Prnx
             S = S+(Uavg(i,j,k+1)+Uavg(i-1,j,k+1)-Uavg(i,j,k-1)-Uavg(i-1,j,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
             S2 = S2+(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
            enddo
           enddo

           S = S/(Prnx*Prny)
           S2 = S2/(Prnx*Prny)
           nom=(grav_acc/temperature_ref)*(proftempflavg(k)+proftempflsgsavg(k))
           denom=((profuwavg(k)+profuwsgsavg(k))*S+(profvwavg(k)+profvwsgsavg(k))*S2)

           if (abs(denom)>1E-5_KND*abs(nom)) then
             S = nom/denom
           else
             S = 100000._KND*sign(1.0_KND,nom)*sign(1.0_KND,denom)
           endif

           write (11,*) zPr(k),S
          enddo
          close(11)
        endif !allocated avgs

     endif !buoyancy>0


     if (computescalars>0) then

        open(11,file="profscal.txt")
        do k = 1,Prnz
         write (11,*) zPr(k),(profscalavg(i,k), i = 1,computescalars)
        enddo
        close(11)

        open(11,file="profscalfl.txt")
        do k = 0,Prnz
         write (11,*) zW(k),(profscalflavg(i,k),profscalflsgsavg(i,k), i = 1,computescalars)
        enddo
        close(11)

        open(11,file="profss.txt")
        do k = 1,Prnz
         write (11,*) zPr(k),(profssavg(i,k), i = 1,computescalars)
        enddo
        close(11)

     endif !computescalars>0

  endif !store%BLprofiles

 end subroutine OutputProfiles







 subroutine OutputOut(U,V,W,Pr)
 real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
 real(KND),dimension(1:,1:,1:),intent(in) :: Pr
 character(70) :: str
 integer i,j,k

   if (store%out>0) then
      open(11,file="out.vtk")
      write (11,"(A)") "# vtk DataFile Version 2.0"
      write (11,"(A)") "CLMM output file"
      write (11,"(A)") "ASCII"
      write (11,"(A)") "DATASET RECTILINEAR_GRID"
      str="DIMENSIONS"
      write (str(12:),*) Prnx,Prny,Prnz
      write (11,"(A)") str
      str="X_COORDINATES"
      write (str(15:),'(i5,2x,a)') Prnx,"float"
      write (11,"(A)") str
      write (11,*) xPr(1:Prnx)
      str="Y_COORDINATES"
      write (str(15:),'(i5,2x,a)') Prny,"float"
      write (11,"(A)") str
      write (11,*) yPr(1:Prny)
      str="Z_COORDINATES"
      write (str(15:),'(i5,2x,a)') Prnz,"float"
      write (11,"(A)") str
      write (11,*) zPr(1:Prnz)
      str="POINT_DATA"
      write (str(12:),*) Prnx*Prny*Prnz
      write (11,"(A)") str

      if (store%out_Pr>0) then
        write (11,"(A)") "SCALARS p float"
        write (11,"(A)") "LOOKUP_TABLE default"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            write (11,*) Pr(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (buoyancy>0.and.store%out_T>0) then
        write (11,*) "SCALARS temperature float"
        write (11,*) "LOOKUP_TABLE default"
        do k = 1,Prnz
          do j = 1,Prny
          do i = 1,Prnx
            write (11,*) Temperature(i,j,k)
          enddo
          enddo
        enddo
        write (11,*)
      endif

      if (store%out_Prtype>0) then
        write (11,*) "SCALARS ptype float"
        write (11,*) "LOOKUP_TABLE default"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            Write (11,*) Prtype(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%out_div>0) then
        write (11,*) "SCALARS div float"
        write (11,*) "LOOKUP_TABLE default"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            Write (11,*) (U(i,j,k)-U(i-1,j,k))/(dxPr(i))+(V(i,j,k)-V(i,j-1,k))/(dyPr(j))+(W(i,j,k)-W(i,j,k-1))/(dzPr(k))
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%out_lambda2>0) then
        write (11,*) "SCALARS lambda2 float"
        write (11,*) "LOOKUP_TABLE default"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            write (11,*) Lambda2(i,j,k,U,V,W)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%out_visc>0) then
        write (11,*) "SCALARS visc float"
        write (11,*) "LOOKUP_TABLE default"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            Write (11,*) Visc(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%out_U>0) then
        write (11,"(A)") "VECTORS u float"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            write (11,*) (U(i,j,k)+U(i-1,j,k))/2._KND,(V(i,j,k)+V(i,j-1,k))/2._KND,(W(i,j,k)+W(i,j,k-1))/2._KND
          enddo
         enddo
        enddo
      endif

      if (store%out_vort>0) then
        write (11,"(A)") "VECTORS vort float"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            write (11,*) (W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                            -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),&
                        (U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                            -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),&
                        (V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                            -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin)
          enddo
         enddo
        enddo
        write (11,*)
      endif
      close(11)
   endif  !store%out
 end subroutine OutputOut








 subroutine OutputScalars
 character(70) :: str
 real(KND),dimension(:,:,:),allocatable :: depos
 character(8) ::  scalname="scalar00"
 integer :: i,j,k,l

  if (computescalars>0) then
    if (store%scalars>0) then
        open(11,file="scalars.vtk")
        write (11,"(A)") "# vtk DataFile Version 2.0"
        write (11,"(A)") "CLMM output file"
        write (11,"(A)") "ASCII"
        write (11,"(A)") "DATASET RECTILINEAR_GRID"
        str="DIMENSIONS"
        write (str(12:),*) Prnx,Prny,Prnz
        write (11,"(A)") str
        str="X_COORDINATES"
        write (str(15:),'(i5,2x,a)') Prnx,"float"
        write (11,"(A)") str
        write (11,*) xPr(1:Prnx)
        str="Y_COORDINATES"
        write (str(15:),'(i5,2x,a)') Prny,"float"
        write (11,"(A)") str
        write (11,*) yPr(1:Prny)
        str="Z_COORDINATES"
        write (str(15:),'(i5,2x,a)') Prnz,"float"
        write (11,"(A)") str
        write (11,*) zPr(1:Prnz)
        str="POINT_DATA"
        write (str(12:),*) Prnx*Prny*Prnz
        write (11,"(A)") str

        do l = 1,computescalars
          write(scalname(7:8),"(I2.2)") l
          write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
          write (11,"(A)") "LOOKUP_TABLE default"
          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              write (11,*) SCALAR(i,j,k,l)
            enddo
           enddo
          enddo
          write (11,*)
        enddo
        close(11)
    endif !store%scalars

    if (computedeposition>0.and.store%deposition>0) then

        allocate(depos(1:Prnx,1:Prny,computescalars))
        depos = GroundDeposition()

        open(11,file="deposition.vtk")
        write (11,"(A)") "# vtk DataFile Version 2.0"
        write (11,"(A)") "CLMM output file"
        write (11,"(A)") "ASCII"
        write (11,"(A)") "DATASET RECTILINEAR_GRID"
        str="DIMENSIONS"
        write (str(12:),*) Prnx,Prny,1
        write (11,"(A)") str
        str="X_COORDINATES"
        write (str(15:),'(i5,2x,a)') Prnx,"float"
        write (11,"(A)") str
        write (11,*) xPr(1:Prnx)
        str="Y_COORDINATES"
        write (str(15:),'(i5,2x,a)') Prny,"float"
        write (11,"(A)") str
        write (11,*) yPr(1:Prny)
        str="Z_COORDINATES"
        write (str(15:),'(i5,2x,a)') 1,"float"
        write (11,"(A)") str
        write (11,*) zW(0)
        str="POINT_DATA"
        write (str(12:),*) Prnx*Prny
        write (11,"(A)") str

        do l = 1,computescalars
          write(scalname(7:8),"(I2.2)") l
          write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
          write (11,"(A)") "LOOKUP_TABLE default"

          do j = 1,Prny
            do i = 1,Prnx
              write (11,*) depos(i,j,l)
            enddo
          enddo

          write (11,*)
        enddo

        close(11)

        deallocate(depos)

    endif  !store%deposition
  endif  !compute_scalars
 end subroutine OutputScalars










 subroutine OutputAvg(U,V,W,Pr)
 real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
 real(KND),dimension(1:,1:,1:),intent(in) :: Pr
 character(70) :: str
 character(2) :: sc
 real(KND),dimension(:,:,:),allocatable :: Upat,Vpat,Wpat
 character(8) ::  scalname="scalar00"
 integer i,j,k,l,m

  if (averaging==1.and.store%avg>0) then
      open(11,file="avg.vtk")
      write (11,"(A)") "# vtk DataFile Version 2.0"
      write (11,"(A)") "CLMM output file"
      write (11,"(A)") "ASCII"
      write (11,"(A)") "DATASET RECTILINEAR_GRID"
      str="DIMENSIONS"
      write (str(12:),*) Prnx,Prny,Prnz
      write (11,"(A)") str
      str="X_COORDINATES"
      write (str(15:),'(i5,2x,a)') Prnx,"float"
      write (11,"(A)") str
      write (11,*) xPr(1:Prnx)
      str="Y_COORDINATES"
      write (str(15:),'(i5,2x,a)') Prny,"float"
      write (11,"(A)") str
      write (11,*) yPr(1:Prny)
      str="Z_COORDINATES"
      write (str(15:),'(i5,2x,a)') Prnz,"float"
      write (11,"(A)") str
      write (11,*) zPr(1:Prnz)
      str="POINT_DATA"
      write (str(12:),*) Prnx*Prny*Prnz
      write (11,"(A)") str

      if (store%avg_Pr>0) then
        write (11,"(A)") "SCALARS p float"
        write (11,"(A)") "LOOKUP_TABLE default"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            write (11,*) Pravg(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%avg_Prtype>0) then
        write (11,*) "SCALARS ptype float"
        write (11,*) "LOOKUP_TABLE default"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            write (11,*) Prtype(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (buoyancy>0.and.store%avg_T>0) then
        write (11,*) "SCALARS temperature float"
        write (11,*) "LOOKUP_TABLE default"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            write (11,*) Temperatureavg(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%avg_U>0) then
        write (11,"(A)") "VECTORS u float"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            write (11,*) (Uavg(i,j,k)+Uavg(i-1,j,k))/2._KND,(Vavg(i,j,k)+Vavg(i,j-1,k))/2._KND,(Wavg(i,j,k)+Wavg(i,j,k-1))/2._KND
          enddo
         enddo
        enddo
      endif

      if (store%avg_vort>0) then
        write (11,"(A)") "VECTORS vort float"
        do k = 1,Prnz
         do j = 1,Prny
          do i = 1,Prnx
            write (11,*) (Wavg(i,j+1,k)-Wavg(i,j-1,k)+Wavg(i,j+1,k-1)-Wavg(i,j-1,k-1))/(4*dxmin)&
                            -(Vavg(i,j,k+1)-Vavg(i,j,k-1)+Vavg(i,j-1,k+1)-Vavg(i,j-1,k-1))/(4*dymin),&
                        (Uavg(i,j,k+1)-Uavg(i,j,k-1)+Uavg(i-1,j,k+1)-Uavg(i-1,j,k-1))/(4*dxmin)&
                            -(Wavg(i+1,j,k)-Wavg(i-1,j,k)+Wavg(i+1,j,k-1)-Wavg(i-1,j,k-1))/(4*dymin),&
                        (Vavg(i+1,j,k)-Vavg(i-1,j,k)+Vavg(i+1,j-1,k)-Vavg(i-1,j-1,k))/(4*dxmin)&
                            -(Uavg(i,j+1,k)-Uavg(i,j-1,k)+Uavg(i-1,j+1,k)-Uavg(i-1,j-1,k))/(4*dymin)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      close(11)
  endif !averaging

  if (averaging==1.and.store%scalarsavg>0) then
      if (computescalars>0) then
          open(11,file="scalarsavg.vtk")
          write (11,"(A)") "# vtk DataFile Version 2.0"
          write (11,"(A)") "CLMM output file"
          write (11,"(A)") "ASCII"
          write (11,"(A)") "DATASET RECTILINEAR_GRID"
          str="DIMENSIONS"
          write (str(12:),*) Prnx,Prny,Prnz
          write (11,"(A)") str
          str="X_COORDINATES"
          write (str(15:),'(i5,2x,a)') Prnx,"float"
          write (11,"(A)") str
          write (11,*) xPr(1:Prnx)
          str="Y_COORDINATES"
          write (str(15:),'(i5,2x,a)') Prny,"float"
          write (11,"(A)") str
          write (11,*) yPr(1:Prny)
          str="Z_COORDINATES"
          write (str(15:),'(i5,2x,a)') Prnz,"float"
          write (11,"(A)") str
          write (11,*) zPr(1:Prnz)
          str="POINT_DATA"
          write (str(12:),*) Prnx*Prny*Prnz
          write (11,"(A)") str

          do l = 1,computescalars
            write(scalname(7:8),"(I2.2)") l
            write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
            write (11,"(A)") "LOOKUP_TABLE default"
            do k = 1,Prnz
             do j = 1,Prny
              do i = 1,Prnx
                write (11,*) SCALARavg(i,j,k,l)
              enddo
             enddo
            enddo
            write (11,*)
          enddo
          close(11)
      endif  !computescalars

  endif !averaging
 endsubroutine OutputAvg




 subroutine OutputUVW(U,V,W)
 real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
 character(70) :: str
 integer i,j,k,l


  if (store%U>0) then
      open(11,file="U.vtk")
      write (11,"(A)") "# vtk DataFile Version 2.0"
      write (11,"(A)") "CLMM output file"
      write (11,"(A)") "ASCII"
      write (11,"(A)") "DATASET RECTILINEAR_GRID"
      str="DIMENSIONS"
      write (str(12:),*) Unx,Uny,Unz
      write (11,"(A)") str
      str="X_COORDINATES"
      write (str(15:),'(i5,2x,a)') Unx,"float"
      write (11,"(A)") str
      write (11,*) xU(1:Unx)
      str="Y_COORDINATES"
      write (str(15:),'(i5,2x,a)') Uny,"float"
      write (11,"(A)") str
      write (11,*) yPr(1:Uny)
      str="Z_COORDINATES"
      write (str(15:),'(i5,2x,a)') Unz,"float"
      write (11,"(A)") str
      write (11,*) zPr(1:Unz)
      str="POINT_DATA"
      write (str(12:),*) Unx*Uny*Unz
      write (11,"(A)") str


      write (11,"(A)") "SCALARS U float"
      write (11,"(A)") "LOOKUP_TABLE default"
      do k = 1,Unz
       do j = 1,Uny
        do i = 1,Unx
          write (11,*) U(i,j,k)
        enddo
       enddo
      enddo
      write (11,*)

      close(11)
  endif

  if (store%V>0) then
      open(11,file="V.vtk")
      write (11,"(A)") "# vtk DataFile Version 2.0"
      write (11,"(A)") "CLMM output file"
      write (11,"(A)") "ASCII"
      write (11,"(A)") "DATASET RECTILINEAR_GRID"
      str="DIMENSIONS"
      write (str(12:),*) Vnx,Vny,Vnz
      write (11,"(A)") str
      str="X_COORDINATES"
      write (str(15:),'(i5,2x,a)') Vnx,"float"
      write (11,"(A)") str
      write (11,*) xPr(1:Vnx)
      str="Y_COORDINATES"
      write (str(15:),'(i5,2x,a)') Vny,"float"
      write (11,"(A)") str
      write (11,*) yV(1:Vny)
      str="Z_COORDINATES"
      write (str(15:),'(i5,2x,a)') Vnz,"float"
      write (11,"(A)") str
      write (11,*) zPr(1:Vnz)
      str="POINT_DATA"
      write (str(12:),*) Vnx*Vny*Vnz
      write (11,"(A)") str


      write (11,"(A)") "SCALARS V float"
      write (11,"(A)") "LOOKUP_TABLE default"
      do k = 1,Vnz
       do j = 1,Vny
        do i = 1,Vnx
          write (11,*) V(i,j,k)
        enddo
       enddo
      enddo
      write (11,*)

      close(11)
  endif !store%V

  if (store%W>0) then
      open(11,file="W.vtk")
      write (11,"(A)") "# vtk DataFile Version 2.0"
      write (11,"(A)") "CLMM output file"
      write (11,"(A)") "ASCII"
      write (11,"(A)") "DATASET RECTILINEAR_GRID"
      str="DIMENSIONS"
      write (str(12:),*) Wnx,Wny,Wnz
      write (11,"(A)") str
      str="X_COORDINATES"
      write (str(15:),'(i5,2x,a)') Wnx,"float"
      write (11,"(A)") str
      write (11,*) xPr(1:Wnx)
      str="Y_COORDINATES"
      write (str(15:),'(i5,2x,a)') Wny,"float"
      write (11,"(A)") str
      write (11,*) yPr(1:Wny)
      str="Z_COORDINATES"
      write (str(15:),'(i5,2x,a)') Wnz,"float"
      write (11,"(A)") str
      write (11,*) zW(1:Wnz)
      str="POINT_DATA"
      write (str(12:),*) Wnx*Wny*Wnz
      write (11,"(A)") str


      write (11,"(A)") "SCALARS W float"
      write (11,"(A)") "LOOKUP_TABLE default"
      do k = 1,Wnz
       do j = 1,Wny
        do i = 1,Wnx
          write (11,*) W(i,j,k)
        enddo
       enddo
      enddo
      write (11,*)

      close(11)
  endif !store%W

  if (store%U_interp/=0) then
    open(11,file="Uinterp.txt")
    do i = 1,size(UIBPoints)
      write(11,*) "xi,yj,zk",UIBPoints(i)%xi,UIBPoints(i)%yj,UIBPoints(i)%zk
      write(11,*) "interp",UIBPoints(i)%interp
      write(11,*) "xi",UIBPoints(i)%IntPoints%xi
      write(11,*) "yj",UIBPoints(i)%IntPoints%yj
      write(11,*) "zk",UIBPoints(i)%IntPoints%zk
      write(11,*) "coefs",UIBPoints(i)%IntPoints%coef
    enddo
    close(11)
  endif

  if (store%V_interp/=0) then
    open(11,file="Vinterp.txt")
    do i = 1,size(VIBPoints)
      write(11,*) "xi,yj,zk",VIBPoints(i)%xi,VIBPoints(i)%yj,VIBPoints(i)%zk
      write(11,*) "interp",VIBPoints(i)%interp
      write(11,*) "xi",VIBPoints(i)%IntPoints%xi
      write(11,*) "yj",VIBPoints(i)%IntPoints%yj
      write(11,*) "zk",VIBPoints(i)%IntPoints%zk
      write(11,*) "coefs",VIBPoints(i)%IntPoints%coef
    enddo
    close(11)
  endif

  if (store%W_interp/=0) then
    open(11,file="Winterp.txt")
    do i = 1,size(WIBPoints)
      write(11,*) "xi,yj,zk",WIBPoints(i)%xi,WIBPoints(i)%yj,WIBPoints(i)%zk
      write(11,*) "interp",WIBPoints(i)%interp
      write(11,*) "xi",WIBPoints(i)%IntPoints%xi
      write(11,*) "yj",WIBPoints(i)%IntPoints%yj
      write(11,*) "zk",WIBPoints(i)%IntPoints%zk
      write(11,*) "coefs",WIBPoints(i)%IntPoints%coef
    enddo
    close(11)
  endif

  if (store%U_interp/=0.and.store%V_interp/=0.and.store%W_interp/=0) then
    open(11,file="Scinterp.txt")
    do i = 1,size(ScalFlIBPoints)
      write(11,*) "xi,yj,zk",ScalFlIBPoints(i)%xi,ScalFlIBPoints(i)%yj,ScalFlIBPoints(i)%zk
      write(11,*) "interp",ScalFlIBPoints(i)%interp
      write(11,*) "xi",ScalFlIBPoints(i)%IntPoints%xi
      write(11,*) "yj",ScalFlIBPoints(i)%IntPoints%yj
      write(11,*) "zk",ScalFlIBPoints(i)%IntPoints%zk
      write(11,*) "coefs",ScalFlIBPoints(i)%IntPoints%coef
    enddo
    close(11)
  endif

 end subroutine OutputUVW






 subroutine OUTPUT(U,V,W,Pr)
 real(KND),dimension(-2:,-2:,-2:),intent(inout) :: U,V,W
 real(KND),dimension(1:,1:,1:),intent(inout) :: Pr

  call BoundU(1,U)
  call BoundU(2,V)
  call BoundU(3,W)

  call OutputOut(U,V,W,Pr)

  call OutputScalars

  if (time>=timeavg1) call OutputAvg(U,V,W,Pr)

  call OutputUVW(U,V,W)

  if (time>=timeavg1) call OutputProfiles(U,V,W,Pr)

  write (*,*) "saved"
 end subroutine OUTPUT


  pure real(KND) function ScalarVerticalFlux(i,j,k,Scal,W)
    integer, intent(in)   :: i,j,k
    real(KND), intent(in) :: Scal(-1:,-1:,-1:), W(-2:,-2:,-2:)

    ScalarVerticalFlux = Scal(i,j,k) * (W(i,j,k)+W(i,j,k-1))/2 + &
                       TDiff(i,j,k) * (Scal(i,j,k+1)-Scal(i,j,k-1)) / (zW(k+1)-zW(k-1))
  end function ScalarVerticalFlux

  pure function Vorticity(i,j,k,U,V,W)
    real(KND), dimension(3) :: Vorticity
    integer, intent(in)     :: i,j,k
    real(KND),dimension(-2:,-2:,-2:), intent(in) :: U,V,W

    Vorticity = (/         (W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                              -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),&
                           (U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                             -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),&
                           (V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                             -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin) /)
  end function Vorticity


  elemental subroutine SetFrameDomain(domain)
    type(TFrameDomain),intent(inout) :: domain


        if (domain%dimension==3) then

          domain%mini = 1
          domain%maxi = Prnx
          domain%minj = 1
          domain%maxj = Prny
          domain%mink = 1
          domain%maxk = Prnz

        else

          if (domain%direction==1) then

            call Gridcoords(domain%mini, domain%minj, domain%mink, domain%position, &
                            (yV(Prny+1)+yV(0))/2._KND, (zW(Prnz+1)+zW(0))/2._KND )

            domain%maxi = domain%mini
            domain%minj = 1
            domain%maxj = Prny
            domain%mink = 1
            domain%maxk = Prnz

          elseif (domain%direction==2) then

            call Gridcoords(domain%mini, domain%minj, domain%mink, (xU(Prnx+1)+xU(0))/2._KND, &
                            domain%position, (zW(Prnz+1)+zW(0))/2._KND )

            domain%maxj = domain%minj
            domain%mini = 1
            domain%maxi = Prnx
            domain%mink = 1
            domain%maxk = Prnz

          else

            call Gridcoords(domain%mini, domain%minj, domain%mink, (xU(Prnx+1)+xU(0))/2._KND, &
                            (yV(Prny+1)+yV(0))/2._KND, domain%position )

            domain%maxk = domain%mink
            domain%mini = 1
            domain%maxi = Prnx
            domain%minj = 1
            domain%maxj = Prny

          endif

        endif

  end subroutine SetFrameDomain



  subroutine FRAME(U,V,W,Pr,n)
    real(KND) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
    integer n,i,j,k,l,m
    character(20) :: fname,fnumber
    character(20) :: fsuffix
    character(70) :: str
    character(8) ::  scalname="scalar00"
    character,parameter :: lf = char(10)
    integer mini,maxi,minj,maxj,mink,maxk
    integer unit

    fsuffix=".vtk"

    vtkformat = binaryvtk

    !$omp parallel do default(private) shared(n,store,time,fsuffix,xPr,yPr,zPr,Prtype,Pr,U,V,W,&
    !$omp                                          scalar,temperature,buoyancy,computescalars,vtkformat)
    do m = 1,size(store%frame_domains)

      fname="frame-"//achar(iachar('a')+m-1)//"-"
      write(fnumber,"(I4.4)") n
      write(*,*) "Saving frame:",fnumber,"   time:",time

      unit = 100+m

      if (vtkformat==textvtk) then

        mini = store%frame_domains(m)%mini
        minj = store%frame_domains(m)%minj
        mink = store%frame_domains(m)%mink
        maxi = store%frame_domains(m)%maxi
        maxj = store%frame_domains(m)%maxj
        maxk = store%frame_domains(m)%maxk

        open(unit,file = trim(fname)//trim(fnumber)//trim(fsuffix))

        write (unit,"(A)") "# vtk DataFile Version 2.0"
        write (unit,"(A)") "CLMM output file"
        write (unit,"(A)") "ASCII"
        write (unit,"(A)") "DATASET RECTILINEAR_GRID"
        str="DIMENSIONS"
        write (str(12:),*) maxi-mini+1,maxj-minj+1,maxk-mink+1
        write (unit,"(A)") str
        str="X_COORDINATES"
        write (str(15:),'(i5,2x,a)') maxi-mini+1,"float"
        write (unit,"(A)") str
        write (unit,*) xPr(mini:maxi)
        str="Y_COORDINATES"
        write (str(15:),'(i5,2x,a)') maxj-minj+1,"float"
        write (unit,"(A)") str
        write (unit,*) yPr(minj:maxj)
        str="Z_COORDINATES"
        write (str(15:),'(i5,2x,a)') maxk-mink+1,"float"
        write (unit,"(A)") str
        write (unit,*) zPr(mink:maxk)
        str="POINT_DATA"
        write (str(12:),*) (maxi-mini+1)*(maxj-minj+1)*(maxk-mink+1)
        write (unit,"(A)") str

        if (store%frame_Pr>0) then
          write (unit,"(A)") "SCALARS p float"
          write (unit,"(A)") "LOOKUP_TABLE default"
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
             if (Prtype(i,j,k)==0) then
              write (unit,*) Pr(i,j,k)
             else
              write (unit,*) 0.
             endif
            enddo
           enddo
          enddo
          write (unit,*)
        endif

        if (store%frame_lambda2>0) then
          write (unit,*) "SCALARS lambda2 float"
          write (unit,*) "LOOKUP_TABLE default"
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
             if (Prtype(i,j,k)==0) then
              write (unit,*) Lambda2(i,j,k,U,V,W)
             else
              write (unit,*) 0.
             endif
            enddo
           enddo
          enddo
          write (unit,*)
        endif

        if (store%frame_scalars>0) then
          scalname(1:6)="scalar"
          do l = 1,computescalars
            write(scalname(7:8),"(I2.2)") l
            write (unit,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
            write (unit,"(A)") "LOOKUP_TABLE default"
            do k = mink,maxk
             do j = minj,maxj
              do i = mini,maxi
               if (Prtype(i,j,k)==0) then
                write (unit,*) SCALAR(i,j,k,l)
               else
                write (unit,*) 0.
               endif
              enddo
             enddo
            enddo
            write (unit,*)
          enddo
        elseif (store%frame_sumscalars>0.and.computescalars>0) then
          write (unit,"(A,1X,A,1X,A)") "SCALARS", "scalar" , "float"
          write (unit,"(A)") "LOOKUP_TABLE default"
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
             if (Prtype(i,j,k)==0) then
              write (unit,*) SUM(SCALAR(i,j,k,:))
             else
              write (unit,*) 0.
             endif
            enddo
           enddo
          enddo
          write (unit,*)
        endif

        if (store%frame_T>0) then
          if (buoyancy>0) then
            write (unit,*) "SCALARS temperature float"
            write (unit,*) "LOOKUP_TABLE default"
            do k = mink,maxk
             do j = minj,maxj
              do i = mini,maxi
               if (Prtype(i,j,k)==0) then
                write (unit,*) Temperature(i,j,k)
               else
                write (unit,*) temperature_ref
               endif
              enddo
             enddo
            enddo
            write (unit,*)
          endif
        endif

        if (buoyancy==1.and.store%frame_tempfl>0) then
          write (unit,*) "SCALARS tempfl float"
          write (unit,*) "LOOKUP_TABLE default"
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
             if (Prtype(i,j,k)==0) then
              write (unit,*) ScalarVerticalFlux(i,j,k,Temperature,W)
             else
              write (unit,*) 0.
             endif
            enddo
           enddo
          enddo
          write (unit,*)
        endif

        if (store%frame_scalfl>0) then
          if (store%frame_scalars>0) then
            scalname(1:6)="scalfl"
            do l = 1,computescalars
              write(scalname(7:8),"(I2.2)") l
              write (unit,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
              write (unit,"(A)") "LOOKUP_TABLE default"
              do k = mink,maxk
               do j = minj,maxj
                do i = mini,maxi
                 if (Prtype(i,j,k)==0) then
                  write (unit,*) ScalarVerticalFlux(i,j,k,Scalar(:,:,:,l),W)
                 else
                  write (unit,*) 0.
                 endif
                enddo
               enddo
              enddo
              write (unit,*)
            enddo
          elseif (store%frame_sumscalars>0.and.computescalars>0) then
            write (unit,"(A)") "SCALARS scalfl float"
            write (unit,"(A)") "LOOKUP_TABLE default"
            do k = mink,maxk
             do j = minj,maxj
              do i = mini,maxi
               if (Prtype(i,j,k)==0) then
                write (unit,*) ScalarVerticalFlux(i,j,k,SUM(Scalar,4),W)
               else
                write (unit,*) 0.
               endif
              enddo
             enddo
            enddo
            write (unit,*)
          endif
        endif

        if (store%frame_U>0) then
          write (unit,"(A)") "VECTORS u float"
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
             if (Prtype(i,j,k)==0) then
              write (unit,*) (U(i,j,k)+U(i-1,j,k))/2._KND,(V(i,j,k)+V(i,j-1,k))/2._KND,&
                             (W(i,j,k)+W(i,j,k-1))/2._KND
             else
              write (unit,*) 0.,0.,0.
             endif
            enddo
           enddo
          enddo
          write (unit,*)
        endif

        if (store%frame_vort>0) then
          write (unit,"(A)") "VECTORS vort float"
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
             if (Prtype(i,j,k)==0) then
              write (unit,*) Vorticity(i,j,k,U,V,W)
             else
              write (unit,*) 0.,0.,0.
             endif
            enddo
           enddo
          enddo
        endif
        close(unit)


      else  !Binary VTK format

        mini = store%frame_domains(m)%mini
        minj = store%frame_domains(m)%minj
        mink = store%frame_domains(m)%mink
        maxi = store%frame_domains(m)%maxi
        maxj = store%frame_domains(m)%maxj
        maxk = store%frame_domains(m)%maxk


        open(unit,file = trim(fname)//trim(fnumber)//trim(fsuffix),&
          access='stream',status='replace',form="unformatted",action="write")

        write (unit) "# vtk DataFile Version 2.0",lf
        write (unit) "CLMM output file",lf
        write (unit) "BINARY",lf
        write (unit) "DATASET RECTILINEAR_GRID",lf
        str="DIMENSIONS"
        write (str(12:),*) maxi-mini+1,maxj-minj+1,maxk-mink+1
        write (unit) str,lf
        str="X_COORDINATES"
        write (str(15:),'(i5,2x,a)') maxi-mini+1,"float"
        write (unit) str,lf
        write (unit) (BigEnd(real(xPr(i),SNG)),i = mini,maxi),lf
        str="Y_COORDINATES"
        write (str(15:),'(i5,2x,a)') maxj-minj+1,"float"
        write (unit) str,lf
        write (unit) (BigEnd(real(yPr(j),SNG)),j = minj,maxj),lf
        str="Z_COORDINATES"
        write (str(15:),'(i5,2x,a)') maxk-mink+1,"float"
        write (unit) str,lf
        write (unit) (BigEnd(real(zPr(k),SNG)),k = mink,maxk),lf
        str="POINT_DATA"
        write (str(12:),*) (maxi-mini+1)*(maxj-minj+1)*(maxk-mink+1)
        write (unit) str,lf

        if (store%frame_Pr>0) then
         write (unit) "SCALARS p float",lf
         write (unit) "LOOKUP_TABLE default",lf
         do k = mink,maxk
          do j = minj,maxj
           do i = mini,maxi
            if (Prtype(i,j,k)==0) then
             write (unit) BigEnd(real(Pr(i,j,k),SNG))
            else
             write (unit) BigEnd(0._SNG)
            endif
           enddo
          enddo
         enddo
         write (unit) lf
        endif

        if (store%frame_lambda2>0) then
         write (unit) "SCALARS lambda2 float",lf
         write (unit) "LOOKUP_TABLE default",lf
         do k = mink,maxk
          do j = minj,maxj
           do i = mini,maxi
            if (Prtype(i,j,k)==0) then
             write (unit) BigEnd(real(Lambda2(i,j,k,U,V,W),SNG))
            else
             write (unit) BigEnd(0._SNG)
            endif
           enddo
          enddo
         enddo
         write (unit) lf
        endif

        if (store%frame_scalars>0) then
         scalname(1:6)="scalar"
         do l = 1,computescalars
          write(scalname(7:8),"(I2.2)") l
          write (unit) "SCALARS ", scalname , " float",lf
          write (unit) "LOOKUP_TABLE default",lf
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
             if (Prtype(i,j,k)==0) then
              write (unit) BigEnd(real(SCALAR(i,j,k,l),SNG))
             else
              write (unit) BigEnd(0._SNG)
             endif
            enddo
           enddo
          enddo
          write (unit) lf
         enddo
        elseif (store%frame_sumscalars>0.and.computescalars>0) then
          write (unit) "SCALARS ", "scalar" , " float",lf
          write (unit) "LOOKUP_TABLE default",lf
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
             if (Prtype(i,j,k)==0) then
              write (unit) BigEnd(real(SUM(SCALAR(i,j,k,:)),SNG))
             else
              write (unit) BigEnd(0._SNG)
             endif
            enddo
           enddo
          enddo
          write (unit) lf
        endif

        if (store%frame_T>0) then
         if (buoyancy>0) then
          write (unit) "SCALARS temperature float",lf
          write (unit) "LOOKUP_TABLE default",lf
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
             if (Prtype(i,j,k)==0) then
              write (unit) BigEnd(real(Temperature(i,j,k),SNG))
             else
              write (unit) BigEnd(real(temperature_ref,SNG))
             endif
            enddo
           enddo
          enddo
          write (unit) lf
         endif
        endif

        if (buoyancy==1.and.store%frame_tempfl>0) then
          write (unit) "SCALARS tempfl float", lf
          write (unit) "LOOKUP_TABLE default", lf
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
             if (Prtype(i,j,k)==0) then
              write (unit) BigEnd(real(ScalarVerticalFlux(i,j,k,Temperature,W),SNG))
             else
              write (unit) BigEnd(0._SNG)
             endif
            enddo
           enddo
          enddo
          write (unit) lf
        endif

        if (store%frame_scalfl>0) then
          if (store%frame_scalars>0) then
            scalname(1:6)="scalfl"
            do l = 1,computescalars
              write(scalname(7:8),"(I2.2)") l
              write (unit) "SCALARS ", scalname , " float", lf
              write (unit) "LOOKUP_TABLE default", lf
              do k = mink,maxk
               do j = minj,maxj
                do i = mini,maxi
                 if (Prtype(i,j,k)==0) then
                  write (unit) BigEnd(real(ScalarVerticalFlux(i,j,k,Scalar(:,:,:,l),W),SNG))
                 else
                  write (unit) BigEnd(0._SNG)
                 endif
                enddo
               enddo
              enddo
              write (unit) lf
            enddo
          elseif (store%frame_sumscalars>0.and.computescalars>0) then
            write (unit) "SCALARS scalfl float", lf
            write (unit) "LOOKUP_TABLE default", lf
            do k = mink,maxk
             do j = minj,maxj
              do i = mini,maxi
               if (Prtype(i,j,k)==0) then
                write (unit) BigEnd(real(ScalarVerticalFlux(i,j,k,SUM(Scalar(:,:,:,:),4),W),SNG))
               else
                write (unit) BigEnd(0._SNG)
               endif
              enddo
             enddo
            enddo
            write (unit) lf
          endif
        endif

        if (store%frame_U>0) then
         write (unit) "VECTORS u float",lf
         do k = mink,maxk
          do j = minj,maxj
           do i = mini,maxi
            if (Prtype(i,j,k)==0) then
             write (unit) BigEnd(real((U(i,j,k)+U(i-1,j,k))/2._KND,SNG)),BigEnd(real((V(i,j,k)+V(i,j-1,k))/2._KND,SNG))&
              ,BigEnd(real((W(i,j,k)+W(i,j,k-1))/2._KND,SNG))
            else
             write (unit) BigEnd(0._SNG),BigEnd(0._SNG),BigEnd(0._SNG)
            endif
           enddo
          enddo
         enddo
         write (unit) lf
        endif

        if (store%frame_vort>0) then
         write (unit) "VECTORS vort float",lf
         do k = mink,maxk
          do j = minj,maxj
           do i = mini,maxi
            if (Prtype(i,j,k)==0) then
             write (unit) BigEnd(real((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                             -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),SNG)),&
                          BigEnd(real((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                            -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),SNG)),&
                          BigEnd(real((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                            -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin),SNG))
            else
             write (unit) BigEnd(0._SNG),BigEnd(0._SNG),BigEnd(0._SNG)
            endif
           enddo
          enddo
         enddo
        endif
        close(unit)

      endif !textvtk/binaryvtk
    enddo   !frame_domains
    !$omp end parallel do
  end subroutine FRAME



  subroutine BLProfiles(U,V,W,temperature,scalar)
  real(KND),dimension(-2:,-2:,-2:) :: U,V,W
  real(KND),dimension(-1:,-1:,-1:),optional:: temperature
  real(KND),dimension(-1:,-1:,-1:,1:),optional:: scalar
  real(KND) :: S
  real(KND) :: Str(1:3,1:3)
  real(KND),allocatable,save ::fp(:),ht(:),gp(:)
  integer   :: i,j,k,l,n
  integer,save :: called = 0

    !$omp parallel private(i,j,k,n,S)
    !$omp do
    do k = 0,Unz+1
     S = 0
     n = 0
     do j = 1,Uny
      do i = 1,Unx
       if (Utype(i,j,k)==0) then
        S = S+U(i,j,k)
        n = n+1
       endif
      enddo
     enddo
     profU(k)=S/max(n,1)
    enddo
    !$omp end do nowait

    !$omp do
    do k = 1,Vnz+1
     S = 0
     n = 0
     do j = 1,Vny
      do i = 1,Vnx
       if (Vtype(i,j,k)==0) then
        S = S+V(i,j,k)
        n = n+1
       endif
      enddo
     enddo
     profV(k)=S/max(n,1)
    enddo
    !$omp end do nowait
    !$omp end parallel

    if (present(temperature)) then
     !$omp parallel do private(i,j,k,n,S)
     do k = 1,Prnz
      S = 0
      n = 0
      do j = 1,Prny
       do i = 1,Prnx
        if (Prtype(i,j,k)==0) then
         S = S+temperature(i,j,k)
         n = n+1
        endif
       enddo
      enddo
      proftemp(k)=S/max(n,1)
     enddo
     !$omp end parallel do
    endif

    if (present(scalar)) then
     do l = 1,computescalars
      !$omp parallel do private(i,j,k,n,S)
      do k = 1,Prnz
       S = 0
       n = 0
       do j = 1,Prny
        do i = 1,Prnx
         if (Prtype(i,j,k)==0) then
          S = S+scalar(i,j,k,l)
          n = n+1
         endif
        enddo
       enddo
       profscal(l,k)=S/max(n,1)
      enddo
     !$omp end parallel do
     enddo
    endif

    if (called==0) then
     allocate(fp(0:Prnx+1),ht(0:Prnz+1),gp(0:Prny+1))
     !$omp parallel workshare
     forall (i = 0:Prnx+1)      fp(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
     forall (k = 0:Prnz+1)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
     forall (j = 0:Prny+1)      gp(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
     !$omp end parallel workshare
    endif

    !$omp parallel private(i,j,k,n,S)
    !$omp do
    do k = 0,Prnz
     S = 0
     n = 0
     do j = 1,Uny
      do i = 1,Unx
       if ((Utype(i,j,k+1)==0.or.Utype(i,j,k)==0).and.(Wtype(i+1,j,k)==0.or.Wtype(i,j,k)==0)) then
        S = S+((ht(k)*U(i,j,k+1)+(1-ht(k))*U(i,j,k))-((1-ht(k))*profU(k)+ht(k)*profU(k+1)))*(fp(i)*W(i+1,j,k)+(1-fp(i))*W(i,j,k))
        n = n+1
       endif
      enddo
     enddo
     profuw(k)=S/max(n,1)
    enddo
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
     S = 0
     n = 0
     do j = 1,Vny
      do i = 1,Vnx
       if ((Vtype(i,j,k+1)==0.or.Vtype(i,j,k)==0).and.(Wtype(i,j+1,k)==0.or.Wtype(i,j,k)==0)) then
        S = S+(ht(k)*V(i,j,k+1)+(1-ht(k))*V(i,j,k)-(ht(k)*profV(k+1)+(1-ht(k))*profV(k)))*(gp(j)*W(i,j+1,k)+(1-gp(j))*W(i,j,k))
        n = n+1
       endif
      enddo
     enddo
     profvw(k)=S/max(n,1)
    enddo
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
     S = 0
     n = 0
     do j = 1,Uny
      do i = 1,Unx
       if (Utype(i,j,k+1)==0.or.Utype(i,j,k)==0) then
        S = S-0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)
        n = n+1
       endif
      enddo
     enddo
     profuwsgs(k)=S/max(n,1)
    enddo
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
     S = 0
     n = 0
     do j = 1,Vny
      do i = 1,Vnx
       if (Vtype(i,j,k+1)==0.or.Vtype(i,j,k)==0) then
        S = S-0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))/dzW(k)
        n = n+1
       endif
      enddo
     enddo
     profvwsgs(k)=S/max(n,1)
    enddo
    !$omp end do nowait

    !$omp do
    do k = 1,Unz
     S = 0
     n = 0
     do j = 1,Uny
      do i = 1,Unx
       if (Utype(i,j,k)==0) then
        S = S+(U(i,j,k)-profU(k))**2
        n = n+1
       endif
      enddo
     enddo
     profuu(k)=S/max(n,1)
    enddo
    !$omp end do nowait

    !$omp do
    do k = 1,Vnz
     S = 0
     n = 0
     do j = 1,Vny
      do i = 1,Vnx
       if (Vtype(i,j,k)==0) then
        S = S+(V(i,j,k)-profV(k))**2
        n = n+1
       endif
      enddo
     enddo
     profvv(k)=S/max(n,1)
    enddo
    !$omp end do nowait

    !$omp do
    do k = 0,Wnz
     S = 0
     n = 0
     do j = 1,Wny
      do i = 1,Wnx
       if (Wtype(i,j,k)==0) then
        S = S+(W(i,j,k))**2
        n = n+1
       endif
      enddo
     enddo
     profww(k)=S/max(n,1)
    enddo
    !$omp end do
    !$omp end parallel

    if (present(temperature)) then
     !proftempfl is computed directly during advection step

     !$omp parallel private(i,j,k,n,S)
     !$omp do
     do k = 1,Prnz
      S = 0
      n = 0
      do j = 1,Prny
       do i = 1,Prnx
        if (Prtype(i,j,k)==0) then
         S = S+(temperature(i,j,k)-profTemp(k))**2
         n = n+1
        endif
       enddo
      enddo
      proftt(k)=S/max(n,1)
     enddo
     !$omp end do nowait

     !$omp do
     do k = 0,Prnz
      S = 0
      n = 0
      do j = 1,Prny
       do i = 1,Prnx
        if (Prtype(i,j,k+1)==0.or.Prtype(i,j,k)==0) then
          S = S-(0.5_KND*(TDiff(i,j,k+1)+TDiff(i,j,k))*(temperature(i,j,k+1)-temperature(i,j,k)))/dzW(k)
          n = n+1
        endif
       enddo
      enddo
      proftempflsgs(k)=S/max(n,1)
     enddo
     !$omp end do
     !$omp end parallel
    endif ! present(temperature)


    if (present(scalar)) then
     do l = 1,computescalars
      !$omp parallel private(i,j,k,n,S)
      !$omp do
      do k = 0,Prnz
       S = 0
       n = 0
       do j = 1,Prny
        do i = 1,Prnx
         if (Prtype(i,j,k+1)==0.or.Prtype(i,j,k)==0) then
          S = S+0.5_KND*(scalar(i,j,k+1,l)+scalar(i,j,k,l))*(W(i,j,k))
          n = n+1
         endif
        enddo
       enddo
       profscalfl(l,k)=S/max(n,1)
      enddo
      !$omp end do nowait

      !$omp do
      do k = 1,Prnz
       S = 0
       n = 0
       do j = 1,Prny
        do i = 1,Prnx
         if (Prtype(i,j,k)==0) then
          S = S+(scalar(i,j,k,l)-profscal(l,k))**2
          n = n+1
         endif
        enddo
       enddo
       profss(l,k)=S/max(n,1)
      enddo
      !$omp end do nowait

      !$omp do
      do k = 0,Prnz
       S = 0
       n = 0
       do j = 1,Prny
        do i = 1,Prnx
         if (Prtype(i,j,k+1)==0.or.Prtype(i,j,k)==0) then
           S = S-(0.5_KND*(TDiff(i,j,k+1)+TDiff(i,j,k))*(scalar(i,j,k+1,l)-scalar(i,j,k,l)))/dzW(k)
           n = n+1
         endif
        enddo
       enddo
       profscalflsgs(l,k)=S/max(n,1)
      enddo
      !$omp end do
      !$omp end parallel
     enddo
    endif ! present(scalar)

    called = 1
  end subroutine BLProfiles


  function TotKE(U,V,W) result(res)
  real(KND) :: res
  real(KND),dimension(-2:,-2:,-2:) :: U,V,W
  real(KND) Um,Vm,Wm
  integer i,j,k
   res = 0
   Um = 0
   Vm = 0
   Wm = 0
   !$omp parallel do private(i,j,k) reduction(+:res)
   do k = 1,Prnz
    do j = 1,Prny
     do i = 1,Prnx
      res = res + (((U(i-1,j,k)+U(i,j,k))/2._KND-Um)**2+&
                      ((V(i,j-1,k)+V(i,j,k))/2._KND-Vm)**2+&
                      ((W(i,j,k-1)+W(i,j,k))/2._KND-Wm)**2)
     enddo
    enddo
   enddo
   !$omp end parallel do
   res = res*lx*lz*lz/2
  end function TotKE

  pure real(KND) function VorticityMag(i,j,k,U,V,W) result(res)
    integer, intent(in) :: i,j,k
    real(KND),dimension(-2:,-2:,-2:), intent(in) :: U,V,W

    res = sum(Vorticity(i,j,k,U,V,W)**2)
    res = Sqrt(res)
  end function VorticityMag


  pure real(KND) function Lambda2(i,j,k,U,V,W) result(res)
    integer, intent(in) :: i,j,k
    real(KND),dimension(-2:,-2:,-2:), intent(in) :: U,V,W

    res = ((U(i,j,k)-U(i-1,j,k))/dxmin)**2
    res = res + ((V(i,j,k)-V(i,j-1,k))/dymin)**2
    res = res + ((W(i,j,k)-W(i,j,k-1))/dzmin)**2
    res = res + ((V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(4*dxmin))**2
    res = res + ((W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(4*dxmin))**2
    res = res + ((U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(4*dymin))**2
    res = res + ((W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(4*dymin))**2
    res = res + ((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(4*dzmin))**2
    res = res + ((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(4*dzmin))**2
    res = -Sqrt(res)
    res = VorticityMag(i,j,k,U,V,W) + res
  end function Lambda2





 subroutine OUTPUTU2(U,V,W)
 real(KND),dimension(-2:,-2:,-2:) :: U,V,W
 integer i,j,k
 character(70) :: str

  open(11,file="U2.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Unx,Uny,Unz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Unx,"float"
  write (11,"(A)") str
  write (11,*) xU(1:Unx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Uny,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:Uny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Unz,"float"
  write (11,"(A)") str
  write (11,*) zPr(1:Unz)
  str="POINT_DATA"
  write (str(12:),*) Unx*Uny*Unz
  write (11,"(A)") str


  write (11,"(A)") "SCALARS U float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k = 1,Unz
   do j = 1,Uny
    do i = 1,Unx
      write (11,*) U(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  close(11)


  open(11,file="V2.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Vnx,Vny,Vnz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Vnx,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:Vnx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Vny,"float"
  write (11,"(A)") str
  write (11,*) yV(1:Vny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Vnz,"float"
  write (11,"(A)") str
  write (11,*) zPr(1:Vnz)
  str="POINT_DATA"
  write (str(12:),*) Vnx*Vny*Vnz
  write (11,"(A)") str


  write (11,"(A)") "SCALARS V float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k = 1,Vnz
   do j = 1,Vny
    do i = 1,Vnx
      write (11,*) V(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  close(11)



  open(11,file="W2.vtk")
  write (11,"(A)") "# vtk DataFile Version 2.0"
  write (11,"(A)") "CLMM output file"
  write (11,"(A)") "ASCII"
  write (11,"(A)") "DATASET RECTILINEAR_GRID"
  str="DIMENSIONS"
  write (str(12:),*) Wnx,Wny,Wnz
  write (11,"(A)") str
  str="X_COORDINATES"
  write (str(15:),'(i5,2x,a)') Wnx,"float"
  write (11,"(A)") str
  write (11,*) xPr(1:Wnx)
  str="Y_COORDINATES"
  write (str(15:),'(i5,2x,a)') Wny,"float"
  write (11,"(A)") str
  write (11,*) yPr(1:Wny)
  str="Z_COORDINATES"
  write (str(15:),'(i5,2x,a)') Wnz,"float"
  write (11,"(A)") str
  write (11,*) zW(1:Wnz)
  str="POINT_DATA"
  write (str(12:),*) Wnx*Wny*Wnz
  write (11,"(A)") str


  write (11,"(A)") "SCALARS W float"
  write (11,"(A)") "LOOKUP_TABLE default"
  do k = 1,Wnz
   do j = 1,Wny
    do i = 1,Wnx
      write (11,*) W(i,j,k)
    enddo
   enddo
  enddo
  write (11,*)
  close(11)
  end subroutine OUTPUTU2


  subroutine OUTINLET(U,V,W)
  !for output of 2d data for use as an inilet condition later
  real(KND),dimension(-2:,-2:,-2:) :: U,V,W
  integer i,j,k
  integer,save ::fnum
  integer,save :: called = 0

   if ((time>=timefram1).and.(time<=timefram2+(timefram2-timefram1)/(frames-1))&
       .and.(time>=timefram1+fnum*(timefram2-timefram1)/(frames-1))) then
    if (called==0) then
     open(101,file="inletframeinfo.unf",form='unformatted',status='replace',action='write')
     write(101) Prny,Prnz  !for check of consistency of grids before use
     write(101) Vny
     write(101) Wnz
     write(101) dxPr(0)
     called = 1
     fnum = 0
    endif
    fnum = fnum+1
    write(101) time-timefram1
    call OUTINLETFRAME(U,V,W,fnum)
   elseif (time>timefram2+(timefram2-timefram1)/(frames-1).and.called==1) then
     close(101)
     called = 2
   endif

  end subroutine OUTINLET


  subroutine OUTINLETFRAME(U,V,W,n)
    real(KND) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    integer n,i,j,k,l,m
    character(12) :: fname
    integer mini,maxi,minj,maxj,mink,maxk

    call Gridcoords(mini,minj,mink,store%frame_domains(1)%position,(yV(Prny+1)+yV(0))/2._KND,(zW(Prnz+1)+zW(0))/2._KND)
    maxi = mini
    minj = 1
    maxj = Prny
    mink = 1
    maxk = Prnz


    fname(1:5)="frame"
    write(fname(6:8),"(I3.3)") n
    fname(9:12)=".unf"
    write(*,*) "Saving frame:",fname(1:6),"   time:",time

    open(11,file = fname,form='unformatted',access='sequential',status='replace',action='write')


    write(11) U(mini,1:Uny,1:Unz)
    write(11) V(mini,1:Vny,1:Vnz)
    write(11) W(mini,1:Wny,1:Wnz)
    if (buoyancy>0) then
         write(11) Temperature(mini,1:Prny,1:Prnz)
    endif
    close(11)

  end subroutine OUTINLETFRAME



  pure real(KND) function TriLinInt(a,b,c,vel000,vel100,vel010,vel001,vel110,vel101,vel011,vel111)
  real(KND),intent(in) :: a,b,c,vel000,vel100,vel010,vel001,vel110,vel101,vel011,vel111

    TriLinInt=   (1-a)*(1-b)*(1-c)*vel000+&
                 a*(1-b)*(1-c)*vel100+&
                 (1-a)*b*(1-c)*vel010+&
                 (1-a)*(1-b)*c*vel001+&
                 a*b*(1-c)*vel110+&
                 a*(1-b)*c*vel101+&
                 (1-a)*b*c*vel011+&
                 a*b*c*vel111

  end function TriLinInt


end module OUTPUTS


module vtkarray
 !Simple module to output arrays for visualization. No physical coordinates are used, only the position in the array.
 !Mostly only for debugging.

!   use iso_fortran_env, only: real32, real64
  implicit none

  integer,parameter :: real32=selected_real_kind(p=6,r=37)
  integer,parameter :: real64=selected_real_kind(p=15,r=200)

  interface VtkArraySimple
    module procedure SVtkArraySimple
    module procedure DVtkArraySimple
  end interface

  contains

    subroutine SVtkArraySimple(fname,A)
      character(len=*)                 :: fname
      real(real32),dimension(1:,1:,1:) :: A
      integer                          :: nx,ny,nz
      integer                          :: i,j,k
      character(len=40)                :: str

      nx=Ubound(A,1)
      ny=Ubound(A,2)
      nz=Ubound(A,3)

      open(11,file=fname)
      write (11,"(A)") "# vtk DataFile Version 2.0"
      write (11,"(A)") "CLMM output file"
      write (11,"(A)") "ASCII"
      write (11,"(A)") "DATASET RECTILINEAR_GRID"
      str="DIMENSIONS"
      write (str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
      write (11,"(A)") trim(str)
      str="X_COORDINATES"
      write (str(15:),'(i5,2x,a)') nx,"float"
      write (11,"(A)") str
      write (11,*) (i, i=1,nx)
      str="Y_COORDINATES"
      write (str(15:),'(i5,2x,a)') ny,"float"
      write (11,"(A)") trim(str)
      write (11,*) (i, i=1,ny)
      str="Z_COORDINATES"
      write (str(15:),'(i5,2x,a)') nz,"float"
      write (11,"(A)") trim(str)
      write (11,*) (i, i=1,nz)
      str="POINT_DATA"
      write (str(12:),*) nx*ny*nz
      write (11,"(A)") trim(str)


      write (11,"(A)") "SCALARS array float"
      write (11,"(A)") "LOOKUP_TABLE default"
      do k=1,nz
       do j=1,ny
        do i=1,nx
          write (11,*) A(i,j,k)
        enddo
        enddo
      enddo
      write (11,*)

    end subroutine SVtkArraySimple

    subroutine DVtkArraySimple(fname,A)
      character(len=*)                 :: fname
      real(real64),dimension(1:,1:,1:) :: A
      integer                          :: nx,ny,nz
      integer                          :: i,j,k
      character(len=40)                :: str

      nx=Ubound(A,1)
      ny=Ubound(A,2)
      nz=Ubound(A,3)

      open(11,file=fname)
      write (11,"(A)") "# vtk DataFile Version 2.0"
      write (11,"(A)") "CLMM output file"
      write (11,"(A)") "ASCII"
      write (11,"(A)") "DATASET RECTILINEAR_GRID"
      str="DIMENSIONS"
      write (str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
      write (11,"(A)") trim(str)
      str="X_COORDINATES"
      write (str(15:),'(i5,2x,a)') nx,"float"
      write (11,"(A)") str
      write (11,*) (i, i=1,nx)
      str="Y_COORDINATES"
      write (str(15:),'(i5,2x,a)') ny,"float"
      write (11,"(A)") trim(str)
      write (11,*) (i, i=1,ny)
      str="Z_COORDINATES"
      write (str(15:),'(i5,2x,a)') nz,"float"
      write (11,"(A)") trim(str)
      write (11,*) (i, i=1,nz)
      str="POINT_DATA"
      write (str(12:),*) nx*ny*nz
      write (11,"(A)") trim(str)


      write (11,"(A)") "SCALARS array float"
      write (11,"(A)") "LOOKUP_TABLE default"
      do k=1,nz
       do j=1,ny
        do i=1,nx
          write (11,*) A(i,j,k)
        enddo
        enddo
      enddo
      write (11,*)

    end subroutine DVtkArraySimple
end module vtkarray



module OUTPUTS
 use PARAMETERS
 use BOUNDARIES
 use SCALARS
 use WALLMODELS, only: WMPoint,FirstWMPoint
 use GEOMETRIC
 use TURBINLET, only: ustarinlet

 implicit none


 private
 public OutTStep, Output, store, display, probes, NumProbes, AllocateOutputs,&
        Tprobe, TOutputSwitches, TDisplaySwitches

 real(KND),dimension(:),allocatable :: profuavg,profuavg2,profvavg,profvavg2,profuuavg,profvvavg,profwwavg,&
                                        profU,profV,profuu,profvv,profww,proftauavg,proftau,proftausgs,proftausgsavg,&
                                        proftemp,proftempfl,proftempavg,proftempavg2,proftempflavg,&
                                        proftempflsgs,proftempflsgsavg,proftt,profttavg,&
                                        profuw,profuwavg,profuwsgs,profuwsgsavg,&
                                        profvw,profvwavg,profvwsgs,profvwsgsavg

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
     Scalaravg=0
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
     Uavg=0
     Vavg=0
     Wavg=0
   endif

   if (store%avg_Pr>0) then
     allocate(Pravg(1:Prnx,1:Prny,1:Prnz))
     Pravg=0
   endif

   if (buoyancy==1.and.store%avg_T>0) then
    allocate(temperatureavg(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    temperatureavg=0
   endif
  endif

  allocate(times(0:nt))

  if (NumProbes>0) then
    allocate(Utime(NumProbes,0:nt),Vtime(NumProbes,0:nt),Wtime(NumProbes,0:nt),Prtime(NumProbes,0:nt))
    times=huge(1.0_KND)
    Utime=huge(1.0_KND)
    Vtime=huge(1.0_KND)
    Wtime=huge(1.0_KND)
    Prtime=huge(1.0_KND)
  endif

  if (NumProbes>0.and.buoyancy==1) then
   allocate(temptime(NumProbes,0:nt))
   temptime=huge(1.0_KND)
  endif

  if (store%deltime>0) then
    allocate(deltime(0:nt))
    deltime=huge(1.0_KND)
  endif

  if (store%tke>0) then
    allocate(tke(0:nt))
    tke=huge(1.0_KND)
  endif

  if (store%deltime>0) then
    allocate(dissip(0:nt))
    dissip=huge(1.0_KND)
    dissip(0)=0
  endif

  if (wallmodeltype>0.and.display%ustar>0) then
    allocate(ustar(2,0:nt))
    ustar=huge(1.0)
  endif

  if (wallmodeltype>0.and.buoyancy==1.and.TBtype(Bo)==DIRICHLET.and.store%tstar>0) then
    allocate(tstar(2,0:nt))
    tstar=huge(1.0)
  endif

  if (computescalars>0) then
    if (store%scalsumtime>0) then
      allocate(scalsumtime(1:computescalars,0:nt))
      scalsumtime=huge(1.0_KND)
    endif

    if (NumProbes>0) then
      allocate(scalptime(1:computescalars,1:3,0:nt))
      scalptime=huge(1.0_KND)
   endif
  endif

  do k=1,NumProbes
    call GridCoords(probes(k)%i,probes(k)%j,probes(k)%k,probes(k)%x,probes(k)%y,probes(k)%z)

    probes(k)%i=max(probes(k)%i,1)
    probes(k)%j=max(probes(k)%j,1)
    probes(k)%k=max(probes(k)%k,1)
    probes(k)%i=min(probes(k)%i,Prnx)
    probes(k)%j=min(probes(k)%j,Prny)
    probes(k)%k=min(probes(k)%k,Prnz)

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
    profU(0:Unz+1)=0
    profV(0:Vnz+1)=0
    profUavg2(1:Unz)=0
    profVavg2(1:Vnz)=0
    profuu(1:Unz)=0
    profvv(1:Vnz)=0
    profww(0:Wnz)=0
    profuuavg(1:Unz)=0
    profvvavg(1:Vnz)=0
    profwwavg(1:Wnz)=0
    profuw(1:Prnz)=0
    profvw(1:Prnz)=0
    proftemp(1:Prnz)=0
    proftempfl(1:Prnz)=0
    proftempflavg(1:Prnz)=0
    proftempavg(1:Prnz)=0
    proftempavg2(1:Prnz)=0
    proftempflsgs(1:Prnz)=0
    proftempflsgsavg(1:Prnz)=0
    proftt(1:Prnz)=0
    profttavg(1:Prnz)=0
  endif

 end subroutine AllocateOutputs













 subroutine OutTStep(U,V,W,Pr,delta)
 real(KND),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
 real(KND),dimension(1:,1:,1:),intent(in) :: Pr
 real(KND),intent(in) :: delta

 integer :: l,i,j,k,nWM
 real(KND) :: S,S2
 type(WMPoint),pointer :: WMP
 integer, save :: fnum=0

   times(step)=time

   do l=1,computescalars
      S=0
      do k=1,Prnz
       do j=1,Prny
        do i=1,Prnx
         if (Prtype(i,j,k)==0) S=S+scalar(i,j,k,l)*dxPr(i)*dyPr(j)*dzPr(k)
        enddo
       enddo
      enddo
      scalsumtime(l,step)=S
   enddo



   do k=1,NumProbes
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

     do l=1,computescalars
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

   endstep=step

   if (display%delta>0) then
     write (*,*) "delta: ",delta
   endif

   if (frames>0)then
      if ((time>=timefram1).and.(time<=timefram2+(timefram2-timefram1)/(frames-1))&
        .and.(time>=timefram1+fnum*(timefram2-timefram1)/(frames-1))) then
       fnum=fnum+1
       call FRAME(U,V,W,Pr,fnum)
      endif
!    call OUTINLET(U,V,W)
   endif


   if ((averaging==1).and.((time>=timeavg1).and.(time<=timeavg2))) then

     if (store%avg_U>0) then
       Uavg=Uavg+U*dt/(timeavg2-timeavg1)
       Vavg=Vavg+V*dt/(timeavg2-timeavg1)
       Wavg=Wavg+W*dt/(timeavg2-timeavg1)
     endif

     if (store%avg_Pr>0) then
       Pravg=Pravg+Pr(1:Prnx,1:Prny,1:Prnz)*dt/(timeavg2-timeavg1)
     endif

     if (buoyancy==1.and.store%avg_T>0) then
       temperatureavg=temperatureavg+temperature*dt/(timeavg2-timeavg1)
     endif

     if (computescalars>0.and.store%scalarsavg>0) then
       SCALARavg=SCALARavg+SCALAR*dt/(timeavg2-timeavg1)
     endif

  endif


  if (wallmodeltype>0.and.(display%ustar>0.or.store%ustar>0)) then
    S=0
    nWM=0
    if (associated(FirstWMPoint)) then
    WMP => FirstWMPoint
    do
     S=S+WMP%ustar
     nWM=nWM+1
     if (associated(WMP%next)) then
      WMP=>WMP%next
     else
      exit
     endif
    enddo
    endif
    if (nWM>0) S=S/nWM
    S2=S*Re

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
    S2=SUM(BsideTFLArr(1:Prnx,1:Prny))/(Prnx*Prny)
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

      if (buoyancy>0) then
       call BLProfiles(U,V,W,temperature)
      else
       call BLProfiles(U,V,W)
      endif

      profuavg=profuavg+profu*dt/(timeavg2-timeavg1)
      profvavg=profvavg+profv*dt/(timeavg2-timeavg1)
      proftempavg=proftempavg+proftemp*dt/(timeavg2-timeavg1)
      profuwavg=profuwavg+profuw*dt/(timeavg2-timeavg1)
      profuwsgsavg=profuwsgsavg+profuwsgs*dt/(timeavg2-timeavg1)
      profvwavg=profvwavg+profvw*dt/(timeavg2-timeavg1)
      profvwsgsavg=profvwsgsavg+profvwsgs*dt/(timeavg2-timeavg1)
      profuuavg=profuuavg+profuu*dt/(timeavg2-timeavg1)
      profvvavg=profvvavg+profvv*dt/(timeavg2-timeavg1)
      profwwavg=profwwavg+profww*dt/(timeavg2-timeavg1)

      if (buoyancy>0) then
        proftempflavg=proftempflavg+proftempfl*dt/(timeavg2-timeavg1)
        proftempflsgsavg=proftempflsgsavg+proftempflsgs*dt/(timeavg2-timeavg1)
        profttavg=profttavg+proftt*dt/(timeavg2-timeavg1)
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

  do k=1,NumProbes

    write(prob,"(i2.2)") k

    open(11,file="Prtimep"//prob//".txt")
    do j=0,endstep
     write (11,*) times(j),Prtime(k,j)
    enddo
    close(11)

    open(11,file="Utimep"//prob//".txt")
    do j=0,endstep
     write (11,*) times(j),Utime(k,j)
    enddo
    close(11)

    open(11,file="Vtimep"//prob//".txt")
    do j=0,endstep
     write (11,*) times(j),Vtime(k,j)
    enddo
    close(11)

    open(11,file="Wtimep"//prob//".txt")
    do j=0,endstep
     write (11,*) times(j),Wtime(k,j)
    enddo
    close(11)

    if (buoyancy==1) then
      open(11,file="temptimep"//prob//".txt")
      do j=0,endstep
       write (11,*) times(j),temptime(k,j)
      enddo
      close(11)
    endif

    if (computescalars>0) then
      open(11,file="scaltimep"//prob//".txt")
      do j=1,endstep
       write (11,*) times(j),scalptime(:,k,j)
      enddo
      close(11)
    endif

  enddo

  if (store%deltime>0) then
    open(11,file="deltime.txt")
    do j=1,endstep
     write (11,*) times(j),deltime(j)
    enddo
    close(11)
  endif

  if (store%tke>0) then
    open(11,file="tke.txt")
    do j=0,endstep
     write (11,*) times(j),tke(j)
    enddo
    close(11)
  endif

  if (store%tke>0.and.store%dissip>0) then
   open(11,file="dissip.txt")
    do j=1,endstep
     write (11,*) times(j),dissip(j)
    enddo
    close(11)
  endif

  if (wallmodeltype>0.and.display%ustar>0) then
    open(11,file="Retau.txt")
    do j=0,endstep
     write (11,*) times(j),ustar(:,j)
    enddo
    close(11)
  endif

  if (wallmodeltype>0.and.buoyancy==1.and.TBtype(Bo)==DIRICHLET.and.store%tstar>0) then
    open(11,file="tflux.txt")
    do j=0,endstep
     write (11,*) times(j),tstar(:,j)
    enddo
    close(11)
  endif


  if (computescalars>0.and.store%scalsumtime>0) then
    open(11,file="scalsumtime.txt")
    do j=1,endstep
     write (11,*) times(j),scalsumtime(:,j)
    enddo
    close(11)
  endif

  if (computescalars>0.and.store%scaltotsumtime>0) then
    open(11,file="scaltotsumtime.txt")
    do j=1,endstep
     write (11,*) times(j),sum(scalsumtime(:,j))
    enddo
    close(11)
  endif


  if (store%BLprofiles>0.and.averaging==1) then

     open(11,file="profu.txt")
     do k=1,Unz
      write (11,*) zPr(k),profuavg(k)
     enddo
     close(11)

     open(11,file="profv.txt")
     do k=1,Vnz
      write (11,*) zPr(k),profvavg(k)
     enddo
     close(11)

     open(11,file="profuu.txt")
     do k=1,Unz
      write (11,*) zPr(k),profuuavg(k)
     enddo
     close(11)

     open(11,file="profvv.txt")
     do k=1,Vnz
      write (11,*) zPr(k),profvvavg(k)
     enddo
     close(11)

     open(11,file="profww.txt")
     do k=1,Wnz
      write (11,*) zW(k),profwwavg(k)
     enddo
     close(11)

     open(11,file="profuw.txt")
     do k=0,Prnz
      write (11,*) zW(k),profuwavg(k),profuwsgsavg(k)
     enddo
     close(11)

     open(11,file="profvw.txt")
     do k=0,Prnz
      write (11,*) zW(k),profvwavg(k),profvwsgsavg(k)
     enddo
     close(11)

     if (buoyancy>0) then
        open(11,file="proftemp.txt")
        do k=1,Prnz
         write (11,*) zPr(k),proftempavg(k)
        enddo
        close(11)

        open(11,file="proftempfl.txt")
        do k=0,Prnz
         write (11,*) zW(k),proftempflavg(k),proftempflsgsavg(k)
        enddo
        close(11)

        open(11,file="proftt.txt")
        do k=1,Prnz
         write (11,*) zPr(k),profttavg(k)
        enddo
        close(11)

        open(11,file="profRig.txt")
        do k=1,Prnz
         S=0
         do j=1,Prny
          do i=1,Prnx
           S=S+Rig(i,j,k,Uavg,Vavg,temperatureavg)
          enddo
         enddo
         S=S/(Prnx*Prny)
         write (11,*) zPr(k),S
        enddo
        close(11)

        open(11,file="profRf.txt")
        do k=1,Prnz
         S=0
         S2=0
         do j=1,Prny
          do i=1,Prnx
           S=S+(Uavg(i,j,k+1)+Uavg(i-1,j,k+1)-Uavg(i,j,k-1)-Uavg(i-1,j,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
           S2=S2+(V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(2._KND*(zPr(k+1)-zPr(k-1)))
          enddo
         enddo

         S=S/(Prnx*Prny)
         S2=S2/(Prnx*Prny)
         nom=(grav_acc/temperature_ref)*(proftempflavg(k)+proftempflsgsavg(k))
         denom=((profuwavg(k)+profuwsgsavg(k))*S+(profvwavg(k)+profvwsgsavg(k))*S2)

         if (abs(denom)>1E-5_KND*abs(nom)) then
           S=nom/denom
         else
           S=100000._KND*sign(1.0_KND,nom)*sign(1.0_KND,denom)
         endif

         write (11,*) zPr(k),S
        enddo
        close(11)

     endif !buoyancy>0

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
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
            write (11,*) Pr(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (buoyancy>0.and.store%out_T>0) then
        write (11,*) "SCALARS temperature float"
        write (11,*) "LOOKUP_TABLE default"
        do k=1,Prnz
          do j=1,Prny
          do i=1,Prnx
            write (11,*) Temperature(i,j,k)
          enddo
          enddo
        enddo
        write (11,*)
      endif

      if (store%out_Prtype>0) then
        write (11,*) "SCALARS ptype float"
        write (11,*) "LOOKUP_TABLE default"
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
            Write (11,*) Prtype(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%out_div>0) then
        write (11,*) "SCALARS div float"
        write (11,*) "LOOKUP_TABLE default"
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
            Write (11,*) (U(i,j,k)-U(i-1,j,k))/(dxPr(i))+(V(i,j,k)-V(i,j-1,k))/(dyPr(j))+(W(i,j,k)-W(i,j,k-1))/(dzPr(k))
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%out_lambda2>0) then
        write (11,*) "SCALARS lambda2 float"
        write (11,*) "LOOKUP_TABLE default"
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
            write (11,*) Lambda2(i,j,k,U,V,W)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%out_visc>0) then
        write (11,*) "SCALARS visc float"
        write (11,*) "LOOKUP_TABLE default"
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
            Write (11,*) Visc(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%out_U>0) then
        write (11,"(A)") "VECTORS u float"
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
            write (11,*) (U(i,j,k)+U(i-1,j,k))/2._KND,(V(i,j,k)+V(i,j-1,k))/2._KND,(W(i,j,k)+W(i,j,k-1))/2._KND
          enddo
         enddo
        enddo
      endif

      if (store%out_vort>0) then
        write (11,"(A)") "VECTORS vort float"
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
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
 type(WMPoint),pointer :: WMP
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

        do l=1,computescalars
          write(scalname(7:8),"(I2.2)") l
          write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
          write (11,"(A)") "LOOKUP_TABLE default"
          do k=1,Prnz
           do j=1,Prny
            do i=1,Prnx
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
        depos=0

        if (associated(FirstWMPoint)) then
          WMP => FirstWMPoint

          do
            if (allocated(WMP%depscalar)) then
              do i=1,computescalars
              depos(WMP%x,WMP%y,i)=depos(WMP%x,WMP%y,i)+WMP%depscalar(i)
              enddo
            endif

            if (associated(WMP%next)) then
              WMP=>WMP%next
            else
              exit
            endif
          enddo
        endif

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

        do l=1,computescalars
          write(scalname(7:8),"(I2.2)") l
          write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
          write (11,"(A)") "LOOKUP_TABLE default"

          do j=1,Prny
            do i=1,Prnx
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
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
            write (11,*) Pravg(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%avg_Prtype>0) then
        write (11,*) "SCALARS ptype float"
        write (11,*) "LOOKUP_TABLE default"
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
            write (11,*) Prtype(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (buoyancy>0.and.store%avg_T>0) then
        write (11,*) "SCALARS temperature float"
        write (11,*) "LOOKUP_TABLE default"
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
            write (11,*) Temperatureavg(i,j,k)
          enddo
         enddo
        enddo
        write (11,*)
      endif

      if (store%avg_vort>0) then
        write (11,"(A)") "VECTORS vort float"
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
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

      if (store%avg_U>0) then
        write (11,"(A)") "VECTORS u float"
        do k=1,Prnz
         do j=1,Prny
          do i=1,Prnx
            write (11,*) (Uavg(i,j,k)+Uavg(i-1,j,k))/2._KND,(Vavg(i,j,k)+Vavg(i,j-1,k))/2._KND,(Wavg(i,j,k)+Wavg(i,j,k-1))/2._KND
          enddo
         enddo
        enddo
        close(11)
      endif
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

          do l=1,computescalars
            write(sc,"(I2.2)") l
            write (11,"(A,1X,A,1X,A)") "SCALARS", sc , "float"
            write (11,"(A)") "LOOKUP_TABLE default"
            do k=1,Prnz
             do j=1,Prny
              do i=1,Prnx
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
 real(KND),dimension(1:Unx,1:Uny,1:Unz) :: Uinterp,Uinterpdir
 real(KND),dimension(1:Vnx,1:Vny,1:Vnz) :: Vinterp,Vinterpdir
 real(KND),dimension(1:Wnx,1:Wny,1:Wnz) :: Winterp,Winterpdir
 type(TIBPoint),pointer :: IBP
 character(70) :: str
 integer i,j,k

  if ((store%U>0.or.store%V>0.or.store%W>0) .and. (store%U_interp>0.or.store%V_interp>0.or.store%W_interp>0)) then
      Uinterp=-1
      Uinterpdir=-1
      Vinterp=-1
      Vinterpdir=-1
      Winterp=-1
      Winterpdir=-1

      if (associated(FirstIBPoint)) then
        IBP => FirstIBPoint
        do
          i=IBP%x
          j=IBP%y
          k=IBP%z

          if (IBP%component==1) then
            Uinterp(i,j,k)=IBP%interp
          elseif (IBP%component==2) then
            Vinterp(i,j,k)=IBP%interp
          elseif (IBP%component==3) then
            Winterp(i,j,k)=IBP%interp
          endif

          if (IBP%component==1) then
            Uinterpdir(i,j,k)=IBP%interpdir
          elseif (IBP%component==2) then
            Vinterpdir(i,j,k)=IBP%interpdir
          elseif (IBP%component==3) then
            Winterpdir(i,j,k)=IBP%interpdir
          endif

          if (associated(IBP%next)) then
            IBP=>IBP%next
          else
            exit
          endif

        enddo
      endif
  endif ! U & interp

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
      do k=1,Unz
       do j=1,Uny
        do i=1,Unx
          write (11,*) U(i,j,k)
        enddo
       enddo
      enddo
      write (11,*)

      if (store%U_interp>0) then
          write (11,"(A)") "SCALARS Uinterp float"
          write (11,"(A)") "LOOKUP_TABLE default"
          do k=1,Unz
           do j=1,Uny
            do i=1,Unx
              write (11,*) Uinterp(i,j,k)
            enddo
           enddo
          enddo
          write (11,*)
          write (11,"(A)") "SCALARS Uinterpdir float"
          write (11,"(A)") "LOOKUP_TABLE default"
          do k=1,Unz
           do j=1,Uny
            do i=1,Unx
              write (11,*) Uinterpdir(i,j,k)
            enddo
           enddo
          enddo
          write (11,*)
      endif

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
      do k=1,Vnz
       do j=1,Vny
        do i=1,Vnx
          write (11,*) V(i,j,k)
        enddo
       enddo
      enddo
      write (11,*)

      if (store%V_interp>0) then
          write (11,"(A)") "SCALARS Vinterp float"
          write (11,"(A)") "LOOKUP_TABLE default"
          do k=1,Vnz
           do j=1,Vny
            do i=1,Vnx
              write (11,*) Vinterp(i,j,k)
            enddo
           enddo
          enddo
          write (11,*)
          write (11,"(A)") "SCALARS Vinterpdir float"
          write (11,"(A)") "LOOKUP_TABLE default"
          do k=1,Vnz
           do j=1,Vny
            do i=1,Vnx
              write (11,*) Vinterpdir(i,j,k)
            enddo
           enddo
          enddo
          write (11,*)
      endif

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
      do k=1,Wnz
       do j=1,Wny
        do i=1,Wnx
          write (11,*) W(i,j,k)
        enddo
       enddo
      enddo
      write (11,*)

      if (store%W_interp>0) then
          write (11,"(A)") "SCALARS Winterp float"
          write (11,"(A)") "LOOKUP_TABLE default"
          do k=1,Wnz
           do j=1,Wny
            do i=1,Wnx
              write (11,*) Winterp(i,j,k)
            enddo
           enddo
          enddo
          write (11,*)
          write (11,"(A)") "SCALARS Winterpdir float"
          write (11,"(A)") "LOOKUP_TABLE default"
          do k=1,Wnz
           do j=1,Wny
            do i=1,Wnx
              write (11,*) Winterpdir(i,j,k)
            enddo
           enddo
          enddo
          write (11,*)
      close(11)

      endif
  endif !store%W
 end subroutine OutputUVW






 subroutine OUTPUT(U,V,W,Pr)
 real(KND),dimension(-2:,-2:,-2:),intent(inout) :: U,V,W
 real(KND),dimension(1:,1:,1:),intent(inout) :: Pr

  call BoundU(1,U)
  call BoundU(2,V)
  call BoundU(3,W)

  call OutputProfiles(U,V,W,Pr)

  call OutputOut(U,V,W,Pr)

  call OutputScalars

  call OutputAvg(U,V,W,Pr)

  call OutputUVW(U,V,W)

  write (*,*) "saved"
 end subroutine OUTPUT



  subroutine FRAME(U,V,W,Pr,n)
  real(KND) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  integer n,i,j,k,l
  character(13) :: fname
  character(70) :: str
  character(8) ::  scalname="scalar00"
  character,parameter :: lf=char(10)
  integer mini,maxi,minj,maxj,mink,maxk


  vtkformat=binaryvtk
  if (vtkformat==textvtk) then

   if (framedimension==3) then
   fname(1:5)="frame"
   write(fname(6:9),"(I4.4)") n
   fname(10:13)=".vtk"
   write(*,*) "Saving frame:",fname(6:9),"   time:",time

   open(11,file=fname)
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

   if (store%frame_Pr>0) then
    write (11,"(A)") "SCALARS p float"
    write (11,"(A)") "LOOKUP_TABLE default"
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (11,*) Pr(i,j,k)
       else
        write (11,*) 0.
       endif
      enddo
     enddo
    enddo
    write (11,*)
   endif

   if (store%frame_lambda2>0) then
    write (11,*) "SCALARS lambda2 float"
    write (11,*) "LOOKUP_TABLE default"
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (11,*) Lambda2(i,j,k,U,V,W)
       else
        write (11,*) 0.
       endif
      enddo
     enddo
    enddo
    write (11,*)
   endif

   if (store%frame_scalars>0) then
    do l=1,computescalars
     write(scalname(7:8),"(I2.2)") l
     write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
     write (11,"(A)") "LOOKUP_TABLE default"
     do k=1,Prnz
      do j=1,Prny
       do i=1,Prnx
        if (Prtype(i,j,k)==0) then
         write (11,*) SCALAR(i,j,k,l)
        else
         write (11,*) 0.
        endif
       enddo
      enddo
     enddo
     write (11,*)
    enddo
   elseif (store%frame_sumscalars>0.and.computescalars>0) then
     write (11,"(A,1X,A,1X,A)") "SCALARS", "scalar" , "float"
     write (11,"(A)") "LOOKUP_TABLE default"
     do k=1,Prnz
      do j=1,Prny
       do i=1,Prnx
        if (Prtype(i,j,k)==0) then
         write (11,*) SUM(SCALAR(i,j,k,:))
        else
         write (11,*) 0.
        endif
       enddo
      enddo
     enddo
     write (11,*)
   endif

   if (store%frame_T>0) then
    if (buoyancy>0) then
      write (11,*) "SCALARS temperature float"
      write (11,*) "LOOKUP_TABLE default"
      do k=1,Prnz
       do j=1,Prny
        do i=1,Prnx
         if (Prtype(i,j,k)==0) then
          write (11,*) Temperature(i,j,k)
         else
          write (11,*) 0.
         endif
        enddo
       enddo
      enddo
      write (11,*)
    endif
   endif

   if (store%frame_U>0) then
    write (11,"(A)") "VECTORS u float"
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (11,*) (U(i,j,k)+U(i-1,j,k))/2._KND,(V(i,j,k)+V(i,j-1,k))/2._KND,(W(i,j,k)+W(i,j,k-1))/2._KND
       else
        write (11,*) 0.,0.,0.
       endif
      enddo
     enddo
    enddo
    write (11,*)
   endif

   if (store%frame_vort>0) then
    write (11,"(A)") "VECTORS vort float"
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (11,*) (W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                       -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),&
                     (U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                       -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),&
                     (V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                        -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin)
       else
        write (11,*) 0.,0.,0.
       endif
      enddo
     enddo
    enddo
    close(11)
   endif


  else
   if (slicedir==1) then
     call Gridcoords(mini,minj,mink,slicex,(yV(Prny+1)+yV(0))/2._KND,(zW(Prnz+1)+zW(0))/2._KND)
     maxi=mini
     minj=1
     maxj=Prny
     mink=1
     maxk=Prnz
    elseif (slicedir==2) then
     call Gridcoords(mini,minj,mink,(xU(Prnx+1)+xU(0))/2._KND,slicex,(zW(Prnz+1)+zW(0))/2._KND)
     maxj=minj
     mini=1
     maxi=Prnx
     mink=1
     maxk=Prnz
    else
     call Gridcoords(mini,minj,mink,(xU(Prnx+1)+xU(0))/2._KND,(yV(Prny+1)+yV(0))/2._KND,slicex)
     maxk=mink
     mini=1
     maxi=Prnx
     minj=1
     maxj=Prny
    endif

   fname(1:5)="frame"
   write(fname(6:9),"(I4.4)") n
   fname(10:13)=".vtk"
   write(*,*) "Saving frame:",fname(6:9),"   time:",time

   open(11,file=fname)
   write (11,"(A)") "# vtk DataFile Version 2.0"
   write (11,"(A)") "CLMM output file"
   write (11,"(A)") "ASCII"
   write (11,"(A)") "DATASET RECTILINEAR_GRID"
   str="DIMENSIONS"
   write (str(12:),*) maxi-mini+1,maxj-minj+1,maxk-mink+1
   write (11,"(A)") str
   str="X_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxi-mini+1,"float"
   write (11,"(A)") str
   write (11,*) xPr(mini:maxi)
   str="Y_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxj-minj+1,"float"
   write (11,"(A)") str
   write (11,*) yPr(minj:maxj)
   str="Z_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxk-mink+1,"float"
   write (11,"(A)") str
   write (11,*) zPr(mink:maxk)
   str="POINT_DATA"
   write (str(12:),*) (maxi-mini+1)*(maxj-minj+1)*(maxk-mink+1)
   write (11,"(A)") str

   if (store%frame_Pr>0) then
    write (11,"(A)") "SCALARS p float"
    write (11,"(A)") "LOOKUP_TABLE default"
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (11,*) Pr(i,j,k)
       else
        write (11,*) 0.
       endif
      enddo
     enddo
    enddo
    write (11,*)
   endif

   if (store%frame_lambda2>0) then
    write (11,*) "SCALARS lambda2 float"
    write (11,*) "LOOKUP_TABLE default"
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (11,*) Lambda2(i,j,k,U,V,W)
       else
        write (11,*) 0.
       endif
      enddo
     enddo
    enddo
    write (11,*)
   endif

   if (store%frame_scalars>0) then
    do l=1,computescalars
     write(scalname(7:8),"(I2.2)") l
     write (11,"(A,1X,A,1X,A)") "SCALARS", scalname , "float"
     write (11,"(A)") "LOOKUP_TABLE default"
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (11,*) SCALAR(i,j,k,l)
        else
         write (11,*) 0.
        endif
       enddo
      enddo
     enddo
     write (11,*)
    enddo
   elseif (store%frame_sumscalars>0.and.computescalars>0) then
     write (11,"(A,1X,A,1X,A)") "SCALARS", "scalar" , "float"
     write (11,"(A)") "LOOKUP_TABLE default"
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (11,*) SUM(SCALAR(i,j,k,:))
        else
         write (11,*) 0.
        endif
       enddo
      enddo
     enddo
     write (11,*)
   endif

   if (store%frame_T>0) then
    if (buoyancy>0) then
     write (11,*) "SCALARS temperature float"
     write (11,*) "LOOKUP_TABLE default"
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (11,*) Temperature(i,j,k)
        else
         write (11,*) 0.
        endif
       enddo
      enddo
     enddo
     write (11,*)
    endif
   endif

   if (store%frame_U>0) then
    write (11,"(A)") "VECTORS u float"
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (11,*) (U(i,j,k)+U(i-1,j,k))/2._KND,(V(i,j,k)+V(i,j-1,k))/2._KND,(W(i,j,k)+W(i,j,k-1))/2._KND
       else
        write (11,*) 0.,0.,0.
       endif
      enddo
     enddo
    enddo
    write (11,*)
   endif

   if (store%frame_vort>0) then
    write (11,"(A)") "VECTORS vort float"
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (11,*) (W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                        -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),&
                     (U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                       -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),&
                     (V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                       -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin)
       else
        write (11,*) 0.,0.,0.
       endif
      enddo
     enddo
    enddo
   endif
   close(11)
   endif


  else  !Binary VTK format



   if (framedimension==3) then
   fname(1:5)="frame"
   write(fname(6:9),"(I4.4)") n
   fname(10:13)=".vtk"
   write(*,*) "Saving frame:",fname(6:9),"   time:",time

   open(20,file=fname,access='stream',status='replace',form="unformatted",action="write")
   write (20) "# vtk DataFile Version 2.0",lf
   write (20) "CLMM output file",lf
   write (20) "BINARY",lf
   write (20) "DATASET RECTILINEAR_GRID",lf
   str="DIMENSIONS"
   write (str(12:),*) Prnx,Prny,Prnz
   write (20) str,lf
   str="X_COORDINATES"
   write (str(15:),'(i5,2x,a)') Prnx,"float"
   write (20) str,lf
   write (20) (BigEnd(real(xPr(i),SNG)),i=1,Prnx),lf
   str="Y_COORDINATES"
   write (str(15:),'(i5,2x,a)') Prny,"float"
   write (20) str,lf
   write (20) (BigEnd(real(yPr(j),SNG)),j=1,Prny),lf
   str="Z_COORDINATES"
   write (str(15:),'(i5,2x,a)') Prnz,"float"
   write (20) str,lf
   write (20) (BigEnd(real(zPr(k),SNG)),k=1,Prnz),lf
   str="POINT_DATA"
   write (str(12:),*) Prnx*Prny*Prnz
   write (20) str,lf

   if (store%frame_Pr>0) then
    write (20) "SCALARS p float",lf
    write (20) "LOOKUP_TABLE default",lf
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (20) real(Pr(1:Prnx,1:Prny,1:Prnz),SNG)
       else
        write (20) BigEnd(0._SNG)
       endif
      enddo
     enddo
    enddo
    write (20) lf
   endif

   if (store%frame_lambda2>0) then
    write (20) "SCALARS lambda2 float",lf
    write (20) "LOOKUP_TABLE default",lf
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (20) real(Lambda2(i,j,k,U,V,W),SNG)
       else
        write (20) BigEnd(0._SNG)
       endif
      enddo
     enddo
    enddo
    write (20) lf
   endif

   if (store%frame_scalars>0) then
    do l=1,computescalars
     write(scalname(7:8),"(I2.2)") l
     write (20) "SCALARS ", scalname , " float",lf
     write (20) "LOOKUP_TABLE default",lf
     do k=1,Prnz
      do j=1,Prny
       do i=1,Prnx
        if (Prtype(i,j,k)==0) then
         write (20) real(SCALAR(i,j,k,l),SNG)
        else
         write (20) BigEnd(0._SNG)
        endif
       enddo
      enddo
     enddo
     write (20) lf
    enddo
   elseif (store%frame_sumscalars>0.and.computescalars>0) then
     write (20) "SCALARS ", "scalar" , " float",lf
     write (20) "LOOKUP_TABLE default",lf
     do k=1,Prnz
      do j=1,Prny
       do i=1,Prnx
        if (Prtype(i,j,k)==0) then
         write (20) BigEnd(real(SUM(SCALAR(i,j,k,:)),SNG))
        else
         write (20) BigEnd(0._SNG)
        endif
       enddo
      enddo
     enddo
     write (20) lf
   endif

   if (store%frame_T>0) then
    if (buoyancy>0) then
      write (20) "SCALARS temperature float",lf
      write (20) "LOOKUP_TABLE default",lf
      do k=1,Prnz
       do j=1,Prny
        do i=1,Prnx
         if (Prtype(i,j,k)==0) then
          write (20) BigEnd(real(Temperature(i,j,k),SNG))
         else
          write (20) BigEnd(0._SNG)
         endif
        enddo
       enddo
      enddo
      write (20) lf
    endif
   endif

   if (store%frame_U>0) then
    write (20) "VECTORS u float",lf
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (20) BigEnd(real((U(i,j,k)+U(i-1,j,k))/2._KND,SNG)),BigEnd(real((V(i,j,k)+V(i,j-1,k))/2._KND,SNG))&
         ,BigEnd(real((W(i,j,k)+W(i,j,k-1))/2._KND,SNG))
       else
        write (20) BigEnd(0._SNG),BigEnd(0._SNG),BigEnd(0._SNG)
       endif
      enddo
     enddo
    enddo
    write (20) lf
   endif

   if (store%frame_vort>0) then
    write (20) "VECTORS vort float",lf
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,k)==0) then
        write (20) BigEnd(real((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                       -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),SNG)),&
                     BigEnd(real((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                       -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),SNG)),&
                     BigEnd(real((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                        -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin),SNG))
       else
        write (20) BigEnd(0._SNG),BigEnd(0._SNG),BigEnd(0._SNG)
       endif
      enddo
     enddo
    enddo
    close(20)
   endif


  else
   if (slicedir==1) then
     call Gridcoords(mini,minj,mink,slicex,(yV(Prny+1)+yV(0))/2._KND,(zW(Prnz+1)+zW(0))/2._KND)
     maxi=mini
     minj=1
     maxj=Prny
     mink=1
     maxk=Prnz
    elseif (slicedir==2) then
     call Gridcoords(mini,minj,mink,(xU(Prnx+1)+xU(0))/2._KND,slicex,(zW(Prnz+1)+zW(0))/2._KND)
     maxj=minj
     mini=1
     maxi=Prnx
     mink=1
     maxk=Prnz
    else
     call Gridcoords(mini,minj,mink,(xU(Prnx+1)+xU(0))/2._KND,(yV(Prny+1)+yV(0))/2._KND,slicex)
     maxk=mink
     mini=1
     maxi=Prnx
     minj=1
     maxj=Prny
    endif

   fname(1:5)="frame"
   write(fname(6:9),"(I4.4)") n
   fname(10:13)=".vtk"
   write(*,*) "Saving frame:",fname(6:9),"   time:",time

   open(20,file=fname,access='stream',status='replace',form="unformatted",action="write")
   write (20) "# vtk DataFile Version 2.0",lf
   write (20) "CLMM output file",lf
   write (20) "BINARY",lf
   write (20) "DATASET RECTILINEAR_GRID",lf
   str="DIMENSIONS"
   write (str(12:),*) maxi-mini+1,maxj-minj+1,maxk-mink+1
   write (20) str,lf
   str="X_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxi-mini+1,"float"
   write (20) str,lf
   write (20) (BigEnd(real(xPr(i),SNG)),i=mini,maxi),lf
   str="Y_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxj-minj+1,"float"
   write (20) str,lf
   write (20) (BigEnd(real(yPr(j),SNG)),j=minj,maxj),lf
   str="Z_COORDINATES"
   write (str(15:),'(i5,2x,a)') maxk-mink+1,"float"
   write (20) str,lf
   write (20) (BigEnd(real(zPr(k),SNG)),k=mink,maxk),lf
   str="POINT_DATA"
   write (str(12:),*) (maxi-mini+1)*(maxj-minj+1)*(maxk-mink+1)
   write (20) str,lf

   if (store%frame_Pr>0) then
    write (20) "SCALARS p float",lf
    write (20) "LOOKUP_TABLE default",lf
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (20) BigEnd(real(Pr(i,j,k),SNG))
       else
        write (20) BigEnd(0._SNG)
       endif
      enddo
     enddo
    enddo
    write (20) lf
   endif

   if (store%frame_lambda2>0) then
    write (20) "SCALARS lambda2 float",lf
    write (20) "LOOKUP_TABLE default",lf
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (20) BigEnd(real(Lambda2(i,j,k,U,V,W),SNG))
       else
        write (20) BigEnd(0._SNG)
       endif
      enddo
     enddo
    enddo
    write (20) lf
   endif

   if (store%frame_scalars>0) then
    do l=1,computescalars
     write(scalname(7:8),"(I2.2)") l
     write (20) "SCALARS ", scalname , " float",lf
     write (20) "LOOKUP_TABLE default",lf
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (20) BigEnd(real(SCALAR(i,j,k,l),SNG))
        else
         write (20) BigEnd(0._SNG)
        endif
       enddo
      enddo
     enddo
     write (20) lf
    enddo
   elseif (store%frame_sumscalars>0.and.computescalars>0) then
     write (20) "SCALARS ", "scalar" , " float",lf
     write (20) "LOOKUP_TABLE default",lf
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (20) BigEnd(real(SUM(SCALAR(i,j,k,:)),SNG))
        else
         write (20) BigEnd(0._SNG)
        endif
       enddo
      enddo
     enddo
     write (20) lf
   endif

   if (store%frame_T>0) then
    if (buoyancy>0) then
     write (20) "SCALARS temperature float",lf
     write (20) "LOOKUP_TABLE default",lf
     do k=mink,maxk
      do j=minj,maxj
       do i=mini,maxi
        if (Prtype(i,j,k)==0) then
         write (20) BigEnd(real(Temperature(i,j,k),SNG))
        else
         write (20) BigEnd(0._SNG)
        endif
       enddo
      enddo
     enddo
     write (20) lf
    endif
   endif

   if (store%frame_U>0) then
    write (20) "VECTORS u float",lf
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (20) BigEnd(real((U(i,j,k)+U(i-1,j,k))/2._KND,SNG)),BigEnd(real((V(i,j,k)+V(i,j-1,k))/2._KND,SNG))&
         ,BigEnd(real((W(i,j,k)+W(i,j,k-1))/2._KND,SNG))
       else
        write (20) BigEnd(0._SNG),BigEnd(0._SNG),BigEnd(0._SNG)
       endif
      enddo
     enddo
    enddo
    write (20) lf
   endif

   if (store%frame_vort>0) then
    write (20) "VECTORS vort float",lf
    do k=mink,maxk
     do j=minj,maxj
      do i=mini,maxi
       if (Prtype(i,j,k)==0) then
        write (20) BigEnd(real((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin)&
                        -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin),SNG)),&
                     BigEnd(real((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin)&
                       -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin),SNG)),&
                     BigEnd(real((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                       -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin),SNG))
       else
        write (20) BigEnd(0._SNG),BigEnd(0._SNG),BigEnd(0._SNG)
       endif
      enddo
     enddo
    enddo
   endif
   close(20)
   endif

  endif
  end subroutine FRAME



  subroutine BLProfiles(U,V,W,temperature)
  real(KND),dimension(-2:,-2:,-2:) :: U,V,W
  real(KND),dimension(-1:,-1:,-1:),optional:: temperature
  real(KND) :: S, S2
  real(KND) Str(1:3,1:3)
  real(KND),allocatable,save ::fp(:),ht(:),gp(:)
  integer i,j,k,n
  integer,save :: called=0

   do k=0,Unz+1
    S=0
    n=0
    do j=1,Uny
     do i=1,Unx
      if (Utype(i,j,k)==0) then
       S=S+U(i,j,k)
       n=n+1
      endif
     enddo
    enddo
    profU(k)=S/n
   enddo

   do k=1,Vnz+1
    S=0
    n=0
    do j=1,Vny
     do i=1,Vnx
      if (Vtype(i,j,k)==0) then
       S=S+V(i,j,k)
       n=n+1
      endif
     enddo
    enddo
    profV(k)=S/n
   enddo

  if (present(temperature)) then
   do k=1,Prnz
    S=0
    n=0
    do j=1,Prny
     do i=1,Prnx
      if (Prtype(i,j,k)==0) then
       S=S+temperature(i,j,k)
       n=n+1
      endif
     enddo
    enddo
    profTemp(k)=S/n
   enddo
  endif

   if (called==0) then
    allocate(fp(0:Prnx+1),ht(0:Prnz+1),gp(0:Prny+1))
    forall (i=0:Prnx+1)      fp(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
    forall (k=0:Prnz+1)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
    forall (j=0:Prny+1)      gp(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
   endif

   do k=0,Prnz
    S=0
    S2=0
    n=0
    do j=1,Uny
     do i=1,Unx
      if ((Utype(i,j,k+1)==0.or.Utype(i,j,k)==0).and.(Wtype(i+1,j,k)==0.or.Wtype(i,j,k)==0)) then
       S=S+((ht(k)*U(i,j,k+1)+(1-ht(k))*U(i,j,k))-((1-ht(k))*profU(k)+ht(k)*profU(k+1)))*(fp(i)*W(i+1,j,k)+(1-fp(i))*W(i,j,k))
       n=n+1
      endif
     enddo
    enddo
    profuw(k)=S/n
   enddo

   do k=0,Prnz
    S=0
    n=0
    do j=1,Vny
     do i=1,Vnx
      if ((Vtype(i,j,k+1)==0.or.Vtype(i,j,k)==0).and.(Wtype(i,j+1,k)==0.or.Wtype(i,j,k)==0)) then
       S=S+(ht(k)*V(i,j,k+1)+(1-ht(k))*V(i,j,k)-(ht(k)*profV(k+1)+(1-ht(k))*profV(k)))*(gp(j)*W(i,j+1,k)+(1-gp(j))*W(i,j,k))
       n=n+1
      endif
     enddo
    enddo
    profvw(k)=S/n
   enddo

   do k=0,Prnz
    S=0
    n=0
    do j=1,Uny
     do i=1,Unx
      if (Utype(i,j,k+1)==0.or.Utype(i,j,k)==0) then
       S=S-0.25_KND*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)
       n=n+1
      endif
     enddo
    enddo
    profuwsgs(k)=S/n
   enddo

   do k=0,Prnz
    S=0
    n=0
    do j=1,Vny
     do i=1,Vnx
      if (Vtype(i,j,k+1)==0.or.Vtype(i,j,k)==0) then
       S=S-0.25_KND*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))/dzW(k)
       n=n+1
      endif
     enddo
    enddo
    profvwsgs(k)=S/n
   enddo

   do k=1,Unz
    S=0
    S2=0
    n=0
    do j=1,Uny
     do i=1,Unx
      if (Utype(i,j,k)==0) then
       S=S+(U(i,j,k)-profU(k))**2
       n=n+1
      endif
     enddo
    enddo
    profuu(k)=S/n
   enddo

   do k=1,Vnz
    S=0
    S2=0
    n=0
    do j=1,Vny
     do i=1,Vnx
      if (Vtype(i,j,k)==0) then
       S=S+(V(i,j,k)-profV(k))**2
       n=n+1
      endif
     enddo
    enddo
    profvv(k)=S/n
   enddo

   do k=0,Wnz
    S=0
    S2=0
    n=0
    do j=1,Wny
     do i=1,Wnx
      if (Wtype(i,j,k)==0) then
       S=S+(W(i,j,k))**2
       n=n+1
      endif
     enddo
    enddo
    profww(k)=S/n
   enddo

  if (present(temperature)) then
   do k=0,Prnz
    S=0
    S2=0
    n=0
    do j=1,Prny
     do i=1,Prnx
      if (Prtype(i,j,k+1)==0.or.Prtype(i,j,k)==0) then
       S=S+0.5_KND*(temperature(i,j,k+1)+temperature(i,j,k))*(W(i,j,k))
       n=n+1
      endif
     enddo
    enddo
    proftempfl(k)=S/n
   enddo

   do k=1,Prnz
    S=0
    S2=0
    n=0
    do j=1,Prny
     do i=1,Prnx
      if (Prtype(i,j,k)==0) then
       S=S+(temperature(i,j,k)-profTemp(k))**2
       n=n+1
      endif
     enddo
    enddo
    proftt(k)=S/n
   enddo

   do k=0,Prnz
    S=0
    S2=0
    n=0
    do j=1,Prny
     do i=1,Prnx
      if (Prtype(i,j,k+1)==0.or.Prtype(i,j,k)==0) then
        S=S-(0.5_KND*(TDiff(i,j,k+1)+TDiff(i,j,k))*(temperature(i,j,k+1)-temperature(i,j,k)))/dzW(k)
        n=n+1
      endif
     enddo
    enddo
    proftempflsgs(k)=S/n
   enddo
  endif
  called=1
  end subroutine BLProfiles


  real(KND) function TotKE(U,V,W)
  real(KND),dimension(-2:,-2:,-2:) :: U,V,W
  real(KND) Um,Vm,Wm
  integer i,j,k
   TotKE=0
   Um=0
   Vm=0
   Wm=0
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
      TotKE=TotKE+(((U(i-1,j,k)+U(i,j,k))/2._KND-Um)**2+&
                    ((V(i,j-1,k)+V(i,j,k))/2._KND-Vm)**2+&
                    ((W(i,j,k-1)+W(i,j,k))/2._KND-Wm)**2)
     enddo
    enddo
   enddo
   TotKE=TotKE*lx*lz*lz/2
  endfunction TotKE

  real(KND) function Vorticity(i,j,k,U,V,W)
  integer i,j,k
  real(KND),dimension(-2:,-2:,-2:) :: U,V,W

    Vorticity=((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dymin)&
                      -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dzmin))**2+&
              ((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dzmin)&
                      -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dxmin))**2+&
              ((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin)&
                      -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin))**2
    Vorticity=Sqrt(Vorticity)
  endfunction Vorticity


  real(KND) function Lambda2(i,j,k,U,V,W)
  integer i,j,k
  real(KND),dimension(-2:,-2:,-2:) :: U,V,W

   Lambda2=((U(i,j,k)-U(i-1,j,k))/dxmin)**2
   Lambda2=Lambda2+((V(i,j,k)-V(i,j-1,k))/dymin)**2
   Lambda2=Lambda2+((W(i,j,k)-W(i,j,k-1))/dzmin)**2
   Lambda2=Lambda2+((V(i+1,j,k)+V(i+1,j-1,k)-V(i-1,j,k)-V(i-1,j-1,k))/(4*dxmin))**2
   Lambda2=Lambda2+((W(i+1,j,k)+W(i+1,j,k-1)-W(i-1,j,k)-W(i-1,j,k-1))/(4*dxmin))**2
   Lambda2=Lambda2+((U(i,j+1,k)+U(i-1,j+1,k)-U(i,j-1,k)-U(i-1,j-1,k))/(4*dymin))**2
   Lambda2=Lambda2+((W(i,j+1,k)+W(i,j+1,k-1)-W(i,j-1,k)-W(i,j-1,k-1))/(4*dymin))**2
   Lambda2=Lambda2+((U(i,j,k+1)+U(i-1,j,k+1)-U(i,j,k-1)-U(i-1,j,k-1))/(4*dzmin))**2
   Lambda2=Lambda2+((V(i,j,k+1)+V(i,j-1,k+1)-V(i,j,k-1)-V(i,j-1,k-1))/(4*dzmin))**2
   Lambda2=-Sqrt(Lambda2)
   Lambda2=Vorticity(i,j,k,U,V,W)+Lambda2


  endfunction Lambda2





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
  do k=1,Unz
   do j=1,Uny
    do i=1,Unx
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
  do k=1,Vnz
   do j=1,Vny
    do i=1,Vnx
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
  do k=1,Wnz
   do j=1,Wny
    do i=1,Wnx
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
  integer,save :: called=0

   if ((time>=timefram1).and.(time<=timefram2+(timefram2-timefram1)/(frames-1))&
       .and.(time>=timefram1+fnum*(timefram2-timefram1)/(frames-1))) then
    if (called==0) then
     open(101,file="inletframeinfo.unf",form='unformatted',status='replace',action='write')
     write(101) Prny,Prnz  !for check of consistency of grids before use
     write(101) Vny
     write(101) Wnz
     write(101) dxPr(0)
     called=1
     fnum=0
    endif
    fnum=fnum+1
    write(101) time-timefram1
    call OUTINLETFRAME(U,V,W,fnum)
   elseif (time>timefram2+(timefram2-timefram1)/(frames-1).and.called==1) then
     close(101)
     called=2
   endif

  end subroutine OUTINLET


  subroutine OUTINLETFRAME(U,V,W,n)
  real(KND) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
  integer n,i,j,k,l
  character(12) :: fname
  integer mini,maxi,minj,maxj,mink,maxk



  call Gridcoords(mini,minj,mink,slicex,(yV(Prny+1)+yV(0))/2._KND,(zW(Prnz+1)+zW(0))/2._KND)
  maxi=mini
  minj=1
  maxj=Prny
  mink=1
  maxk=Prnz


  fname(1:5)="frame"
  write(fname(6:8),"(I3.3)") n
  fname(9:12)=".unf"
  write(*,*) "Saving frame:",fname(1:6),"   time:",time

  open(11,file=fname,form='unformatted',access='sequential',status='replace',action='write')


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

  endfunction TriLinInt


  real(SNG) function BigEnd(x)
  real(SNG),intent(in)::x
  integer,save:: called=0
  logical,save:: littleendian !endianess of the machine
  integer(selected_int_kind(1)),dimension(4):: bytes !may not work on some processors

    if (called==0) then
      bytes=transfer(1,bytes,4)
      if (bytes(4)==1) then
       littleendian=.false.
      else
       littleendian=.true.
      endif
      called=1
    endif

    if (.not.littleendian) then
     BigEnd=x
    else
     bytes=transfer(x,bytes,4)
     BigEnd=transfer(bytes(4:1:-1),x)
    endif
  end function BigEnd

end module OUTPUTS

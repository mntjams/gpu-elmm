program DIPLOMKA
use TSTEPS
use PARAMETERS
use BOUNDARIES
use INITIAL
use OUTPUTS
use SCALARS
use GEOMETRIC
use WALLMODELS
!use fish
implicit none

!definice promennych, souboru...

  
  real(KND),allocatable:: U(:,:,:),V(:,:,:),W(:,:,:),Pr(:,:,:)
  real(KND),allocatable:: Ulat(:,:,:),Vlat(:,:,:),Wlat(:,:,:)
  integer l,fnum,i,j,k,nWM
  integer probei,probej,probek,pUi,pUj,pUk,pVi,pVj,pVk,pWi,pWj,pWk
  integer procs,omp_get_num_procs!,nthreads,omp_get_num_threads
  real(KND) delta,S,S2
  TYPE(WMPoint),pointer:: WMP

  pi=2.0_KND*acos(0.0_KND)
!inicializace parametru
  !S=0
  !!$OMP PARALLEL PRIVATE(i) SHARED(S)
  !nthreads=OMP_GET_NUM_THREADS()
  !write (*,*) "Threads:",nthreads
  !!$OMP DO
  !do i=1,10000
  ! nthreads=OMP_GET_NUM_THREADS()
  ! S=S+1
  !enddo
  !!$OMP ENDDO
  !write(*,*) "Threads in do:",nthreads,S
  !!$OMP ENDPARALLEL
 
  write (*,*) "Reading parameters..."
  call readparams

  
  write (*,*) "Setting up boundary conditions..."
  call readBounds
  
  write (*,*) "Allocating arrays..."
  allocate(U(-2:Unx+3,-2:Uny+3,-2:Unz+3))
  allocate(V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
  allocate(W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
  allocate(Pr(1:Unx+1,1:Vny+1,1:Wnz+1))
  U=100000
  V=100000
  W=100000
  Pr=0000
  if (averaging==1) then
   allocate(Uavg(-2:Unx+3,-2:Uny+3,-2:Unz+3))
   allocate(Vavg(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
   allocate(Wavg(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
   allocate(Pravg(1:Prnx,1:Prny,1:Prnz))
   Uavg=0
   Vavg=0
   Wavg=0
   Pravg=0
   if (buoyancy==1) then
    allocate(temperatureavg(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    temperatureavg=0
   endif 
  endif
  allocate(times(0:nt),Utime(0:nt),Vtime(0:nt),Wtime(0:nt),Prtime(0:nt),temptime(0:nt))
  allocate(deltime(0:nt),tke(0:nt),dissip(0:nt),dissip2(0:nt))
  times=huge(1.0_KND)
  Utime=huge(1.0_KND)
  Vtime=huge(1.0_KND)
  Wtime=huge(1.0_KND)
  Prtime=huge(1.0_KND)
  temptime=huge(1.0_KND)
  tke=huge(1.0_KND)
  dissip=huge(1.0_KND)
  dissip2=huge(1.0_KND)

  if (computescalars>0) then
    allocate(scaltime(1:computescalars,0:nt))
    allocate(scalptime(1:computescalars,1:3,0:nt))
    scaltime=huge(1.0_KND)
    scalptime=huge(1.0_KND)
  endif

  if (tempmet==2) then
      allocate(Ulat(-2:Unx+3,-2:Uny+3,-2:Unz+3))
      allocate(Vlat(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
      allocate(Wlat(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
      Ulat=huge(1.0_KND)
      Vlat=huge(1.0_KND)
      Wlat=huge(1.0_KND)
  endif
  allocate(Visc(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
  if (buoyancy>0.or.computescalars>0) then
    allocate(TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
    TDiff=huge(1.0_KND)
  endif
  if (fullstress==1) then
   allocate(tstress(3,3,-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
  endif

  if ((tasktype==4.or.tasktype==7).and.averaging==1) then
   write(*,*) tasktype
    allocate(profUavg(1:Uny),profUavg2(1:Uny),profuuavg(1:Uny),profvvavg(1:Vny),profwwavg(1:Wny),proftauavg(1:Prny))
    allocate(profU(1:Uny),profuu(1:Uny),profvv(1:Vny),profww(1:Wny),proftau(1:Prny))
   profUavg=0
   profUavg2=0
   profuuavg=0
   profvvavg=0
   profwwavg=0
   profU=0
   profuu=0
   profvv=0
   profww=0
   proftauavg=0
   proftau=0
  elseif (tasktype==8) then
  write(*,*) tasktype
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

  call GridCoordsPr(probei,probej,probek,probex,probey,probez)
  probei=max(probei,1)
  probej=max(probej,1)
  probek=max(probek,1)
  probei=min(probei,Prnx)
  probej=min(probej,Prny)
  probek=min(probek,Prnz)
  call GridCoordsU(pUi,pUj,pUk,probex,probey,probez)
  call GridCoordsV(pVi,pVj,pVk,probex,probey,probez)
  call GridCoordsW(pWi,pWj,pWk,probex,probey,probez)

!pocatecni podminky
 Pr=0
 write (*,*) "Setting up initial conditions..."
 call INITCONDS(U,V,W,Pr)
! write (*,*) "Saving initial condiions..."
! call OUTPUT(U,V,W,Pr,"in.vtk")
 !casove kroky
 !v kazdem podkroku vsechny cleny -adv. cleny PPM - diff. cleny C-N
 fnum=0
 time=starttime
 step=0
 times(0)=0


 if (wallmodeltype>0) OPEN(21,file="Retau.txt")
 if (tasktype==5) OPEN(22,file="momthick.txt")
 if (wallmodeltype>0.and.TBtypeB==DIRICHLET) OPEN(23,file="tfluxtime.txt")
 
   Utime(step)=Trilinint((probex-xU(pUi))/(xU(pUi+1)-xU(pUi)),&
                          (probey-yPr(pUj))/(yPr(pUj+1)-yPr(pUj)),&
                          (probez-zPr(pUk))/(zPr(pUk+1)-zPr(pUk)),&
                          U(pUi,pUj,pUk),U(pUi+1,pUj,pUk),U(pUi,pUj+1,pUk),U(pUi,pUj,pUk+1),&
                          U(pUi+1,pUj+1,pUk),U(pUi+1,pUj,pUk+1),U(pUi,pUj+1,pUk+1),U(pUi+1,pUj+1,pUk+1))
   Vtime(step)=Trilinint((probex-xPr(pVi))/(xPr(pVi+1)-xPr(pVi)),&
                          (probey-yV(pVj))/(yV(pVj+1)-yV(pVj)),&
                          (probez-zPr(pVk))/(zPr(pVk+1)-zPr(pVk)),&
                          V(pVi,pVj,pVk),V(pVi+1,pVj,pVk),V(pVi,pVj+1,pVk),V(pVi,pVj,pVk+1),&
                          V(pVi+1,pVj+1,pVk),V(pVi+1,pVj,pVk+1),V(pVi,pVj+1,pVk+1),V(pVi+1,pVj+1,pVk+1))
   Wtime(step)=Trilinint((probex-xPr(pWi))/(xPr(pWi+1)-xPr(pWi)),&
                          (probey-yPr(pWj))/(yPr(pWj+1)-yPr(pWj)),&
                          (probez-zW(pWk))/(zW(pWk+1)-zW(pWk)),&
                          W(pWi,pWj,pWk),W(pWi+1,pWj,pWk),W(pWi,pWj+1,pWk),W(pWi,pWj,pWk+1),&
                          W(pWi+1,pWj+1,pWk),W(pWi+1,pWj,pWk+1),W(pWi,pWj+1,pWk+1),W(pWi+1,pWj+1,pWk+1))
   write (*,*) probek+1,Wnz+1,size(PR,3)
   Prtime(step)=Trilinint((probex-xPr(probei))/(xPr(probei+1)-xPr(probei)),&
                          (probey-yPr(probej))/(yPr(probej+1)-yPr(probej)),&
                          (probez-zPr(probek))/(zPr(probek+1)-zPr(probek)),&
                           Pr(probei,probej,probek),&
                           Pr(min(probei+1,Unx+1),probej,probek),&
                           Pr(probei,min(probej+1,Vny+1),probek),&
                           Pr(probei,probej,min(probek+1,Wnz+1)),&
                           Pr(min(probei+1,Unx+1),min(probej+1,Vny+1),probek),&
                           Pr(min(probei+1,Unx+1),probej,min(probek+1,Wnz+1)),&
                           Pr(probei,min(probej+1,Vny+1),min(probek+1,Wnz+1)),&
                           Pr(min(probei+1,Unx+1),min(probej+1,Vny+1),min(probek+1,Wnz+1)))


   if (buoyancy==1)  temptime(step)=Trilinint((probex-xPr(probei))/(xPr(probei+1)-xPr(probei)),&
                          (probey-yPr(probej))/(yPr(probej+1)-yPr(probej)),&
                          (probez-zPr(probek))/(zPr(probek+1)-zPr(probek)),&
                           temperature(probei,probej,probek),&
                           temperature(min(probei+1,Unx+1),probej,probek),&
                           temperature(probei,min(probej+1,Vny+1),probek),&
                           temperature(probei,probej,min(probek+1,Wnz+1)),&
                           temperature(min(probei+1,Unx+1),min(probej+1,Vny+1),probek),&
                           temperature(min(probei+1,Unx+1),probej,min(probek+1,Wnz+1)),&
                           temperature(probei,min(probej+1,Vny+1),min(probek+1,Wnz+1)),&
                           temperature(min(probei+1,Unx+1),min(probej+1,Vny+1),min(probek+1,Wnz+1)))
   tke(step)=totke(U,V,W)










 write (*,*) "Computing..."
  do step=1,nt
   write (*,*) "-----------"
   write (*,*) " "   
   write (*,*) " "
   write (*,*) "tstep:",step



   if (tempmet==1) then
       call TMARCHEUL(U,V,W,Pr,delta)
   elseif (tempmet==2) then
       if (step==1) then 
         Ulat=U
         Vlat=V
         Wlat=W
         call TMARCHEUL(U,V,W,Pr,delta)
       else
!          call TMARCHAB2(U,V,W,Pr,Ulat,Vlat,Wlat,delta)
       endif
    elseif (tempmet==3) then
       call TMARCHRK2(U,V,W,Pr,delta)
    elseif (tempmet==4) then
       call TMARCHRK3(U,V,W,Pr,delta)
    elseif (tempmet==5) then
       call TMARCHSHIFTINLET(U,V,W,Pr,delta)
   endif




   time=time+dt
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
      scaltime(l,step)=S
   enddo
   if (tasktype==1.or.tasktype==8) then 
    do l=1,computescalars
      k=INT((Prnz-1)*(2._KND-zPr(1))/(zPr(Prnz)-zPr(1)))
      j=INT((Prny-1)*(0._KND-yPr(1))/(yPr(Prny)-yPr(1)))

      i=INT((Prnx-1)*(11._KND-xPr(1))/(xPr(Prnx)-xPr(1)))
      scalptime(l,1,step)=SCALAR(i,j,k,l)

      i=INT((Prnx-1)*(18._KND-xPr(1))/(xPr(Prnx)-xPr(1)))
      scalptime(l,2,step)=SCALAR(i,j,k,l)

      i=INT((Prnx-1)*(25._KND-xPr(1))/(xPr(Prnx)-xPr(1)))
      scalptime(l,3,step)=SCALAR(i,j,k,l)
    enddo
   endif

   Utime(step)=Trilinint((probex-xU(pUi))/(xU(pUi+1)-xU(pUi)),&
                          (probey-yPr(pUj))/(yPr(pUj+1)-yPr(pUj)),&
                          (probez-zPr(pUk))/(zPr(pUk+1)-zPr(pUk)),&
                          U(pUi,pUj,pUk),U(pUi+1,pUj,pUk),U(pUi,pUj+1,pUk),U(pUi,pUj,pUk+1),&
                          U(pUi+1,pUj+1,pUk),U(pUi+1,pUj,pUk+1),U(pUi,pUj+1,pUk+1),U(pUi+1,pUj+1,pUk+1))
   Vtime(step)=Trilinint((probex-xPr(pVi))/(xPr(pVi+1)-xPr(pVi)),&
                          (probey-yV(pVj))/(yV(pVj+1)-yV(pVj)),&
                          (probez-zPr(pVk))/(zPr(pVk+1)-zPr(pVk)),&
                          V(pVi,pVj,pVk),V(pVi+1,pVj,pVk),V(pVi,pVj+1,pVk),V(pVi,pVj,pVk+1),&
                          V(pVi+1,pVj+1,pVk),V(pVi+1,pVj,pVk+1),V(pVi,pVj+1,pVk+1),V(pVi+1,pVj+1,pVk+1))
   Wtime(step)=Trilinint((probex-xPr(pWi))/(xPr(pWi+1)-xPr(pWi)),&
                          (probey-yPr(pWj))/(yPr(pWj+1)-yPr(pWj)),&
                          (probez-zW(pWk))/(zW(pWk+1)-zW(pWk)),&
                          W(pWi,pWj,pWk),W(pWi+1,pWj,pWk),W(pWi,pWj+1,pWk),W(pWi,pWj,pWk+1),&
                          W(pWi+1,pWj+1,pWk),W(pWi+1,pWj,pWk+1),W(pWi,pWj+1,pWk+1),W(pWi+1,pWj+1,pWk+1))
   Prtime(step)=Trilinint((probex-xPr(probei))/(xPr(probei+1)-xPr(probei)),&
                          (probey-yPr(probej))/(yPr(probej+1)-yPr(probej)),&
                          (probez-zPr(probek))/(zPr(probek+1)-zPr(probek)),&
                           Pr(probei,probej,probek),&
                           Pr(min(probei+1,Unx+1),probej,probek),&
                           Pr(probei,min(probej+1,Vny+1),probek),&
                           Pr(probei,probej,min(probek+1,Wnz+1)),&
                           Pr(min(probei+1,Unx+1),min(probej+1,Vny+1),probek),&
                           Pr(min(probei+1,Unx+1),probej,min(probek+1,Wnz+1)),&
                           Pr(probei,min(probej+1,Vny+1),min(probek+1,Wnz+1)),&
                           Pr(min(probei+1,Unx+1),min(probej+1,Vny+1),min(probek+1,Wnz+1)))


   if (buoyancy==1)  temptime(step)=Trilinint((probex-xPr(probei))/(xPr(probei+1)-xPr(probei)),&
                          (probey-yPr(probej))/(yPr(probej+1)-yPr(probej)),&
                          (probez-zPr(probek))/(zPr(probek+1)-zPr(probek)),&
        temperature(probei,probej,probek),temperature(probei+1,probej,probek),temperature(probei,probej+1,probek),&
        temperature(probei,probej,probek+1),&
        temperature(probei+1,probej+1,probek),temperature(probei+1,probej,probek+1),temperature(probei,probej+1,probek+1),&
        temperature(probei+1,probej+1,probek+1))
   tke(step)=totke(U,V,W)
   dissip(step)=(tke(step-1)-tke(step))/(times(step)-times(step-1))
   !dissip2(step)=0
   !do k=1,Prnz
   ! do j=1,Prny
   !  do i=1,Prnx
   !   dissip2(step)=dissip2(step)+Vorticity(i,j,k,U,V,W)
   !  enddo
   ! enddo
   !enddo
   !dissip2(step)=(1._KND/Re)*dissip2(step)/(Prnx*Prny*Prnz)
   deltime(step)=delta/dt
   endstep=step
   write (*,*) "delta: ",delta
   if (frames>0)then
      if ((time>=timefram1).and.(time<=timefram2+(timefram2-timefram1)/(frames-1))&
        .and.(time>=timefram1+fnum*(timefram2-timefram1)/(frames-1))) then
       fnum=fnum+1
       call FRAME(U,V,W,Pr,fnum)
      endif
!    call OUTINLET(U,V,W)
   endif
   if ((averaging==1).and.((time>=timeavg1).and.(time<=timeavg2))) then
     Uavg=Uavg+U*dt/(timeavg2-timeavg1)
     Vavg=Vavg+V*dt/(timeavg2-timeavg1)
     Wavg=Wavg+W*dt/(timeavg2-timeavg1)
     Pravg=Pravg+Pr(1:Prnx,1:Prny,1:Prnz)*dt/(timeavg2-timeavg1)
     if (buoyancy==1) then
      temperatureavg=temperatureavg+temperature*dt/(timeavg2-timeavg1)
     endif
     if (computescalars>0) then
      SCALARavg=SCALARavg+SCALAR*dt/(timeavg2-timeavg1)
     endif
   endif
  if (wallmodeltype>0) then
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
      EXIT
     endif
    enddo
    endif
    if (nWM>0) S=S/nWM
    S2=S*Re
    meanustar=S
    if (allocated(ustarinlet)) then
     write(*,*) "ustar:",S,"Re_tau:",S2,"u*inlet",ustarinlet(1)
    else
     write(*,*) "ustar:",S,"Re_tau:",S2
    endif
    write(21,*) time,S2,S
   endif
   if (wallmodeltype>0.and.TBtypeB==DIRICHLET) then
    S2=SUM(BsideTFLArr(1:Prnx,1:Prny))/(Prnx*Prny)
    S=-S*S2
    write(*,*) "Tstar",S,"tflux", S2
    write (23,*) time,S2,S
   endif
   if (tasktype==5) then
    write(22,*) time,MOMTHICK(U)
   endif
   if ((tasktype==4.or.tasktype==7)) then
    if ((averaging==1).and.((time>=timeavg1).and.(time<=timeavg2))) then
     if (time>=(timeavg1+timeavg2)/2._KND) then
      profuavg2=profuavg*(timeavg2-timeavg1)/(time-timeavg1)
     endif

     call PROFILES(U,V,W,profu,proftau,profuu, profvv, profww)

     profuavg=profuavg+profu*dt/(timeavg2-timeavg1)
     profvvavg=profvvavg+profvv*dt/(timeavg2-timeavg1)
     profwwavg=profwwavg+profww*dt/(timeavg2-timeavg1)
     if (time>=(timeavg1+timeavg2)/2._KND) then
      proftauavg=proftauavg+2._KND*proftau*dt/(timeavg2-timeavg1)
      profuuavg=profuuavg+2._KND*profuu*dt/(timeavg2-timeavg1)
     endif
    endif
   elseif ((tasktype==8)) then
    if ((averaging==1).and.((time>=timeavg1).and.(time<=timeavg2))) then
!      if (time>=(timeavg1+timeavg2)/2._KND) then
!       profuavg2=profuavg*(timeavg2-timeavg1)/(time-timeavg1)
!       profvavg2=profvavg*(timeavg2-timeavg1)/(time-timeavg1)
!       if (buoyancy>0) proftempavg2=proftempavg*(timeavg2-timeavg1)/(time-timeavg1)
!      endif

     if (buoyancy>0) then
      call CBLPROFILES(U,V,W,temperature)
     else
      call CBLPROFILES(U,V,W)
     endif

     profuavg=profuavg+profu*dt/(timeavg2-timeavg1)
     profvavg=profvavg+profv*dt/(timeavg2-timeavg1)
     proftempavg=proftempavg+proftemp*dt/(timeavg2-timeavg1)
!      if (time>=(timeavg1+timeavg2)/2._KND) then
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
!      endif
    endif
   endif

   if (time+dt>=endtime.and.time<endtime-0.01*dt.and.tasktype==3) then
    call OUTPUTU2(U,V,W)
    write (*,*) "U2 saved."
   endif
   
    if ((steady==1) .and. (delta<eps)) then
      write (*,*) "Steady state reached."
      EXIT
    endif  
    if ((steady==0) .and. (time>=endtime)) then
      write (*,*) "Time limit reached."
      EXIT
    endif

    if (dt<abs(CFL*dxmin/Uinlet/100._KND)) then
      write (*,*) "Solution diverged."
      EXIT
    endif  
  enddo 
  if ((wallmodeltype>0)) CLOSE(21)
  if (tasktype==5) CLOSE(22)
 if (wallmodeltype>0.and.TBtypeB==DIRICHLET) CLOSE(23)



 write (*,*) "Saving results..."

 call OUTPUT(U,V,W,Pr)



 deallocate(U,V,W,Pr,xU,yV,zW,xPr,yPr,zPr,dxU,dyV,dzW,dxPr,dyPr,dzPr)
 if (allocated(Utype))  deallocate(Utype)
 if (allocated(Vtype))  deallocate(Vtype)
 if (allocated(Wtype))  deallocate(Wtype)
 if (allocated(Prtype))  deallocate(Prtype)
 if (allocated(U))  deallocate(U)
 if (allocated(V))  deallocate(V)
 if (allocated(W))  deallocate(W)
 if (allocated(Pr))  deallocate(Pr)
 if (allocated(Uavg))  deallocate(Uavg)
 if (allocated(Vavg))  deallocate(Vavg)
 if (allocated(Wavg))  deallocate(Wavg)
 if (allocated(Pravg))  deallocate(Pravg)
 if (allocated(temperature))  deallocate(temperature)
 if (allocated(temperatureavg))  deallocate(temperatureavg)
 if (allocated(Scalar))  deallocate(Scalar)
 if (allocated(Scalaravg))  deallocate(Scalaravg)
 if (allocated(Ulat))  deallocate(Ulat)
 if (allocated(Vlat))  deallocate(Vlat)
 if (allocated(times)) deallocate(times)
 if (allocated(Utime)) deallocate(Utime)
 if (allocated(Vtime)) deallocate(Vtime)
 if (allocated(Wtime)) deallocate(Wtime)
 if (allocated(Prtime)) deallocate(Prtime)
 if (allocated(temptime)) deallocate(temptime)
 !if (allocated(Wtime)) deallocate(Wtime)
 if (allocated(CDtime)) deallocate(CDtime)
 if (allocated(CLtime)) deallocate(CLtime)
 if (allocated(deltime)) deallocate(deltime)
 if (allocated(tke)) deallocate(tke)
 if (allocated(dissip)) deallocate(dissip)
 if (allocated(dissip2)) deallocate(dissip2)
 if (allocated(Uin)) deallocate(Uin)
 if (allocated(Vin)) deallocate(Vin)
 if (allocated(Win)) deallocate(Win)
 if (allocated(ustarinlet))  deallocate(ustarinlet)
 if (allocated(transformtensor))  deallocate(transformtensor)
 if (allocated(bfilt)) deallocate(bfilt)
 if (allocated(Apx)) deallocate(Apx)
 if (allocated(Apy)) deallocate(Apy)
 if (allocated(Apz)) deallocate(Apz)
 if (allocated(Aw)) deallocate(Aw)
 if (allocated(Ae)) deallocate(Ae)
 if (allocated(As)) deallocate(As)
 if (allocated(An)) deallocate(An)
 if (allocated(Ab)) deallocate(Ab)
 if (allocated(At)) deallocate(At)
 if (allocated(Ru)) deallocate(Ru)
 if (allocated(Rv)) deallocate(Rv)
 if (allocated(Rw)) deallocate(Rw)
 if (allocated(Uinavg))  deallocate(Uinavg)
 if (allocated(Vinavg))  deallocate(Vinavg)
 if (allocated(Winavg))  deallocate(Winavg)
 if (allocated(RedBlackU))  deallocate(RedBlackU)
 if (allocated(RedBlackV))  deallocate(RedBlackV)
 if (allocated(RedBlackW))  deallocate(RedBlackW)
 if (allocated(RedBlackPr))  deallocate(RedBlackPr)
 if (allocated(Visc))  deallocate(Visc)
 if (allocated(TDiff))  deallocate(TDiff)
 if (allocated(partrho))  deallocate(partrho)
 if (allocated(partdiam))  deallocate(partdiam)
 if (allocated(percdistrib))  deallocate(percdistrib)
 if (allocated(profU))  deallocate(profU)
 if (allocated(profV))  deallocate(profV)
 if (allocated(profUavg))  deallocate(profUavg)
 if (allocated(profVavg))  deallocate(profVavg)
 if (allocated(profUavg2))  deallocate(profUavg2)
 if (allocated(profVavg2))  deallocate(profVavg2)
 if (allocated(profuuavg))  deallocate(profuuavg)
 if (allocated(profvvavg))  deallocate(profvvavg)
 if (allocated(profwwavg))  deallocate(profwwavg)
 if (allocated(proftauavg))  deallocate(proftauavg)
 if (allocated(profuu))  deallocate(profuu)
 if (allocated(profvv))  deallocate(profvv)
 if (allocated(profww))  deallocate(profww)
 if (allocated(proftau))  deallocate(proftau)
 if (allocated(proftausgs))  deallocate(proftausgs)
 if (allocated(proftausgsavg))  deallocate(proftausgsavg)

 if (allocated(proftemp))  deallocate(proftemp)
 if (allocated(proftempavg))  deallocate(proftempavg)
 if (allocated(proftempavg2))  deallocate(proftempavg2)
 if (allocated(proftemp))  deallocate(proftempfl)
 if (allocated(proftempflavg))  deallocate(proftempflavg)
 if (allocated(proftempflsgs))  deallocate(proftempflsgs)
 if (allocated(proftempflsgsavg))  deallocate(proftempflsgsavg)
 if (allocated(proftt))  deallocate(proftt)
 if (allocated(profttavg))  deallocate(profttavg)

 if (allocated(scaltime))  deallocate(scaltime)
 if (allocated(scalptime))  deallocate(scalptime)


 if (associated(FirstSB)) call DeallSB(FirstSB)
 if (associated(FirstIBPoint)) call DeallIBP(FirstIBPoint)
 if (associated(FirstWMPoint)) call DeallWMP(FirstWMPoint)

 !if (poissmet==2) call deallfish
end

module Outputs
  use Parameters
  use Boundaries
  use Scalars
  use Wallmodels, only: GroundDeposition, GroundUstar, wallmodeltype
  use ImmersedBoundary
  use Turbinlet, only: Ustar_inlet
  use Endianness
  use FreeUnit
  use StaggeredFrames

  implicit none


  private
  public store, display, probes, scalar_probes, number_of_probes, number_of_scalar_probes, &
         OutTStep, Output, AllocateOutputs, SetFrameDomain, ReadProbes, proftempfl,&
         StaggeredFrameDomains

  real(knd),dimension(:),allocatable :: profuavg,profuavg2,profvavg,profvavg2,profuuavg,profvvavg,profwwavg,&
                                         profU,profV,profuu,profvv,profww,proftauavg,proftau,proftausgs,proftausgsavg,&
                                         proftemp,proftempfl,proftempavg,proftempavg2,proftempflavg,&
                                         proftempflsgs,proftempflsgsavg,proftt,profttavg,&
                                         profmoist,profmoistfl,profmoistavg,profmoistavg2,profmoistflavg,&
                                         profmoistflsgs,profmoistflsgsavg,profmm,profmmavg,&
                                         profuw,profuwavg,profuwsgs,profuwsgsavg,&
                                         profvw,profvwavg,profvwsgs,profvwsgsavg

  real(knd),dimension(:,:),allocatable ::profscal,profscalfl,profscalavg,profscalavg2,profscalflavg,&  !which scalar, height
                                         profscalflsgs,profscalflsgsavg,profss,profssavg

  real(knd),allocatable :: U_avg(:,:,:),V_avg(:,:,:),W_avg(:,:,:)
  real(knd),allocatable :: Pr_avg(:,:,:)
  real(knd),allocatable :: Temperature_avg(:,:,:)
  real(knd),allocatable :: Moisture_avg(:,:,:)
  real(knd),allocatable :: Scalar_avg(:,:,:,:)

  real(TIM),allocatable,dimension(:) :: times                                !times of the timesteps

  real(knd),allocatable,dimension(:) :: CDtime,CLtime,deltime,tke,dissip

  real(knd),allocatable,dimension(:,:) :: ustar,tstar                        !first index differentiates flux from friction number
                                                                             !second index is time

  real(knd),allocatable,dimension(:,:) :: Utime,Vtime,Wtime,Prtime,temptime,moisttime  !position, time

  real(knd),allocatable,dimension(:,:,:) :: scalptime                        !which scalar, position, time
  real(knd),allocatable,dimension(:,:) :: scalsumtime                        !which scalar, time

  !number of probes in space to collect timed data
  integer :: number_of_probes = 0, number_of_scalar_probes = 0

  type TProbe
    integer :: Ui,Uj,Uk,Vi,Vj,Vk,Wi,Wj,Wk    !grid coordinates of probes in the U,V,W grids
    integer :: i,j,k                         !grid coordinates of probes in the scalar grid
    real(knd) :: x,y,z                       !physical coordinates of probes
  end type TProbe

  !for flow variables including scalar ones (temperature, moisture)
  type(TProbe),allocatable,dimension(:),save :: probes

  !for passive scalars
  type(TProbe),allocatable,dimension(:),save :: scalar_probes

  type TFrameDomain
    integer   :: dimension,direction
    real(knd) :: position
    integer   :: mini, maxi, minj, maxj, mink, maxk
  end type TFrameDomain


  type TOutputSwitches
    integer :: U = 0
    integer :: U_interp = 0
    integer :: V = 0
    integer :: V_interp = 0
    integer :: W = 0
    integer :: W_interp = 0

    integer :: out = 1
    integer :: avg = 1 !1..one avg.vtk, 2.. velocity only in separate Xavg.vtk, 3.. both

    integer :: scalars = 1
    integer :: scalarsavg = 1

    integer :: deposition = 1


    integer :: out_U = 1
    integer :: out_vort = 1
    integer :: out_Pr = 1
    integer :: out_Prtype = 1
    integer :: out_lambda2 = 1
    integer :: out_temperature = 1
    integer :: out_moisture = 1
    integer :: out_div = 0
    integer :: out_visc = 0

    integer :: avg_U = 1
    integer :: avg_vort = 0
    integer :: avg_Pr = 1
    integer :: avg_Prtype = 0
    integer :: avg_temperature = 1
    integer :: avg_moisture = 1

    integer :: frame_U = 1
    integer :: frame_vort = 0
    integer :: frame_Pr = 0
    integer :: frame_lambda2 = 0
    integer :: frame_scalars = 0
    integer :: frame_sumscalars = 1
    integer :: frame_temperature = 1
    integer :: frame_moisture = 1
    integer :: frame_temperature_flux = 0
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
  end type

  type(TDisplaySwitches),save :: display

  type(TStaggeredFrameDomain),allocatable :: StaggeredFrameDomains(:)

  !line feed
  character,parameter :: lf = achar(10)


contains


 subroutine ReadProbes(ps,nps,pfile)
    type(TProbe),allocatable,intent(out) :: ps(:)
    integer,intent(out)     :: nps
    character(*),intent(in) ::pfile
    integer i,io,unit
    real(knd) :: tmp3(3)

    open(newunit=unit, file=pfile, status="old", action="read",iostat=io)
    if (io/=0) then
      write(*,*) "Error: File ",pfile," could not be opened!"
      stop
    else
      nps = 0
      do
        read(unit,*,iostat=io) tmp3
        if (io/=0) exit
        nps = nps + 1
      end do

      allocate(ps(nps))

      rewind(unit)
      do i=1,nps
        read(unit,*) ps(i)%x,ps(i)%y,ps(i)%z
      end do
      close(unit)
    end if
  end subroutine

  subroutine AllocateOutputs
    integer :: k

    if (num_of_scalars>0) then
     if (averaging==1.and.store%scalarsavg==1) then
      allocate(Scalar_avg(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,num_of_scalars))
      Scalar_avg = 0
     end if
    else
     if (averaging==1) then
      allocate(Scalar_avg(0,0,0,0))
     end if
    end if


   if (averaging==1) then
    if (store%avg_U==1) then
      allocate(U_avg(-2:Unx+3,-2:Uny+3,-2:Unz+3))
      allocate(V_avg(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
      allocate(W_avg(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
      U_avg = 0
      V_avg = 0
      W_avg = 0
    end if

    if (store%avg_Pr==1) then
      allocate(Pr_avg(1:Prnx,1:Prny,1:Prnz))
      Pr_avg = 0
    end if

    if (enable_buoyancy==1.and.store%avg_temperature==1) then
     allocate(Temperature_avg(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
     Temperature_avg = 0
    end if

    if (enable_moisture==1.and.store%avg_moisture==1) then
     allocate(Moisture_avg(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
     Moisture_avg = 0
    end if
   end if

   allocate(times(0:nt))

   if (number_of_probes>0) then
     allocate(Utime(number_of_probes,0:nt),Vtime(number_of_probes,0:nt),Wtime(number_of_probes,0:nt),Prtime(number_of_probes,0:nt))
     times = huge(1.0_knd)
     Utime = huge(1.0_knd)
     Vtime = huge(1.0_knd)
     Wtime = huge(1.0_knd)
     Prtime = huge(1.0_knd)

     if (enable_buoyancy==1) then
       allocate(temptime(number_of_probes,0:nt))
       temptime = huge(1.0_knd)
     end if

      if (enable_moisture==1) then
       allocate(moisttime(number_of_probes,0:nt))
       moisttime = huge(1.0_knd)
     end if

   end if


   if (store%deltime==1) then
     allocate(deltime(0:nt))
     deltime = huge(1.0_knd)
   end if

   if (store%tke==1) then
     allocate(tke(0:nt))
     tke = huge(1.0_knd)
   end if

   if (store%deltime==1) then
     allocate(dissip(0:nt))
     dissip = huge(1.0_knd)
     dissip(0)=0
   end if

   if (wallmodeltype>0.and.(display%ustar==1.or.store%ustar==1)) then
     allocate(ustar(2,0:nt))
     ustar = huge(1.0)
   end if

   if (wallmodeltype>0.and.enable_buoyancy==1.and.TempBtype(Bo)==DIRICHLET.and.(display%tstar==1.or.store%tstar==1)) then
     allocate(tstar(2,0:nt))
     tstar = huge(1.0)
   end if

   if (num_of_scalars>0) then
     if (store%scalsumtime==1) then
       allocate(scalsumtime(1:num_of_scalars,0:nt))
       scalsumtime = huge(1.0_knd)
     end if

     if (number_of_scalar_probes>0) then
       allocate(scalptime(1:num_of_scalars,1:number_of_scalar_probes,0:nt))
       scalptime = huge(1.0_knd)
    end if
   end if

   do k = 1,number_of_probes
     associate(p => probes(k))
       call GridCoords(p%i,p%j,p%k,p%x,p%y,p%z)

       p%i = max(p%i,1)
       p%j = max(p%j,1)
       p%k = max(p%k,1)
       p%i = min(p%i,Prnx)
       p%j = min(p%j,Prny)
       p%k = min(p%k,Prnz)

       call GridCoordsU(p%Ui,p%Uj,p%Uk,p%x,p%y,p%z)
       call GridCoordsV(p%Vi,p%Vj,p%Vk,p%x,p%y,p%z)
       call GridCoordsW(p%Wi,p%Wj,p%Wk,p%x,p%y,p%z)
     end associate
   end do

   do k = 1,number_of_scalar_probes
     associate(p => scalar_probes(k))
       call GridCoords(p%i,p%j,p%k,p%x,p%y,p%z)

       p%i = max(p%i,1)
       p%j = max(p%j,1)
       p%k = max(p%k,1)
       p%i = min(p%i,Prnx)
       p%j = min(p%j,Prny)
       p%k = min(p%k,Prnz)
     end associate
   end do

   if (store%BLprofiles==1.and.averaging==1) then
     allocate(profU(0:Unz+1),profV(0:Vnz+1),profUavg(0:Unz+1),profVavg(0:Vnz+1),profUavg2(0:Unz+1),profVavg2(0:Vnz+1))
     allocate(profuuavg(1:Unz),profvvavg(1:Vnz),profwwavg(0:Wnz),profuu(1:Unz),profvv(1:Vnz),profww(0:Wnz))
     allocate(profuw(0:Prnz),profuwavg(0:Prnz),profuwsgs(0:Prnz),profuwsgsavg(0:Prnz))
     allocate(profvw(0:Prnz),profvwavg(0:Prnz),profvwsgs(0:Prnz),profvwsgsavg(0:Prnz))

     allocate(proftemp(1:Prnz),proftempfl(0:Prnz),proftempavg(1:Prnz),proftempavg2(1:Prnz),proftempflavg(0:Prnz))
     allocate(proftempflsgs(0:Prnz),proftempflsgsavg(0:Prnz),proftt(1:Prnz),profttavg(1:Prnz))

     allocate(profmoist(1:Prnz),profmoistfl(0:Prnz),profmoistavg(1:Prnz),profmoistavg2(1:Prnz),profmoistflavg(0:Prnz))
     allocate(profmoistflsgs(0:Prnz),profmoistflsgsavg(0:Prnz),profmm(1:Prnz),profmmavg(1:Prnz))

     allocate(profscal(num_of_scalars,1:Prnz),profscalfl(num_of_scalars,0:Prnz),profscalavg(num_of_scalars,1:Prnz))
     allocate(profscalavg2(num_of_scalars,1:Prnz),profscalflavg(num_of_scalars,0:Prnz))
     allocate(profscalflsgs(num_of_scalars,0:Prnz),profscalflsgsavg(num_of_scalars,0:Prnz))
     allocate(profss(num_of_scalars,1:Prnz),profssavg(num_of_scalars,1:Prnz))

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
     profuwsgs = 0
     profvwsgs = 0
     profuwsgsavg = 0
     profvwsgsavg = 0

     proftemp = 0
     proftempavg = 0
     proftempavg2 = 0
     proftempfl = 0
     proftempflavg = 0
     proftempflsgs = 0
     proftempflsgsavg = 0
     proftt = 0
     profttavg = 0

     profmoist = 0
     profmoistavg = 0
     profmoistavg2 = 0
     profmoistfl = 0
     profmoistflavg = 0
     profmoistflsgs = 0
     profmoistflsgsavg = 0
     profmm = 0
     profmmavg = 0

     if (num_of_scalars>0) then
       profscal = 0
       profscalavg = 0
       profscalavg2 = 0
       profscalfl = 0
       profscalflavg = 0
       profscalflsgs = 0
       profscalflsgsavg = 0
       profss = 0
       profssavg = 0
     end if
   else
     !to avoid the necessity of an allocatable dummy argument
     allocate(proftempfl(0))

     allocate(profmoistfl(0))

   end if


   call GetEndianness

#if _WIN32 || _WIN64
   call execute_command_line("mkdir output")
#else
   call execute_command_line("mkdir -p output")
#endif

   contains 
     !FIXME: delete this as soon as supported by ifort
     subroutine execute_command_line(cmd)
       character(*) :: cmd
       call system(cmd)
     end subroutine
  end subroutine AllocateOutputs













  subroutine OutTStep(U,V,W,Pr,Temperature,Moisture,Scalar,delta)
    real(knd),dimension(-2:,-2:,-2:),intent(in)   :: U,V,W
    real(knd),dimension(1:,1:,1:),intent(in)      :: Pr
    real(knd),dimension(-1:,-1:,-1:),intent(in)   :: Temperature
    real(knd),dimension(-1:,-1:,-1:),intent(in)   :: Moisture
    real(knd),dimension(-1:,-1:,-1:,:),intent(in) :: Scalar
    real(knd),intent(in) :: delta

    integer :: l,i,j,k
    real(knd) :: S,S2
    real(knd) :: time_weight
    integer, save :: fnum = 0

    times(step)=time

    if (store%scalsumtime==1) then
      do l = 1,num_of_scalars
         S = 0
         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
            if (Prtype(i,j,k)<=0) S = S + Scalar(i,j,k,l)*dxPr(i)*dyPr(j)*dzPr(k)
           end do
          end do
         end do
         scalsumtime(l,step) = S
      end do
    end if



    do k = 1,number_of_probes
      associate (p=> probes(k))
        Utime(k,step)=Trilinint((p%x-xU(p%Ui))/(xU(p%Ui+1)-xU(p%Ui)),&
                             (p%y-yPr(p%Uj))/(yPr(p%Uj+1)-yPr(p%Uj)),&
                             (p%z-zPr(p%Uk))/(zPr(p%Uk+1)-zPr(p%Uk)),&
                             U(p%Ui,p%Uj,p%Uk),U(p%Ui+1,p%Uj,p%Uk),&
                             U(p%Ui,p%Uj+1,p%Uk),U(p%Ui,p%Uj,p%Uk+1),&
                             U(p%Ui+1,p%Uj+1,p%Uk),U(p%Ui+1,p%Uj,p%Uk+1),&
                             U(p%Ui,p%Uj+1,p%Uk+1),U(p%Ui+1,p%Uj+1,p%Uk+1))

        Vtime(k,step)=Trilinint((p%x-xPr(p%Vi))/(xPr(p%Vi+1)-xPr(p%Vi)),&
                             (p%y-yV(p%Vj))/(yV(p%Vj+1)-yV(p%Vj)),&
                             (p%z-zPr(p%Vk))/(zPr(p%Vk+1)-zPr(p%Vk)),&
                             V(p%Vi,p%Vj,p%Vk),V(p%Vi+1,p%Vj,p%Vk),&
                             V(p%Vi,p%Vj+1,p%Vk),V(p%Vi,p%Vj,p%Vk+1),&
                             V(p%Vi+1,p%Vj+1,p%Vk),V(p%Vi+1,p%Vj,p%Vk+1),&
                             V(p%Vi,p%Vj+1,p%Vk+1),V(p%Vi+1,p%Vj+1,p%Vk+1))

        Wtime(k,step)=Trilinint((p%x-xPr(p%Wi))/(xPr(p%Wi+1)-xPr(p%Wi)),&
                             (p%y-yPr(p%Wj))/(yPr(p%Wj+1)-yPr(p%Wj)),&
                             (p%z-zW(p%Wk))/(zW(p%Wk+1)-zW(p%Wk)),&
                             W(p%Wi,p%Wj,p%Wk),W(p%Wi+1,p%Wj,p%Wk),&
                             W(p%Wi,p%Wj+1,p%Wk),W(p%Wi,p%Wj,p%Wk+1),&
                             W(p%Wi+1,p%Wj+1,p%Wk),W(p%Wi+1,p%Wj,p%Wk+1),&
                             W(p%Wi,p%Wj+1,p%Wk+1),W(p%Wi+1,p%Wj+1,p%Wk+1))

        Prtime(k,step)=Trilinint((p%x-xPr(p%i))/(xPr(p%i+1)-xPr(p%i)),&
                             (p%y-yPr(p%j))/(yPr(p%j+1)-yPr(p%j)),&
                             (p%z-zPr(p%k))/(zPr(p%k+1)-zPr(p%k)),&
                              Pr(p%i,p%j,p%k),&
                              Pr(min(p%i+1,Unx+1),p%j,p%k),&
                              Pr(p%i,min(p%j+1,Vny+1),p%k),&
                              Pr(p%i,p%j,min(p%k+1,Wnz+1)),&
                              Pr(min(p%i+1,Unx+1),min(p%j+1,Vny+1),p%k),&
                              Pr(min(p%i+1,Unx+1),p%j,min(p%k+1,Wnz+1)),&
                              Pr(p%i,min(p%j+1,Vny+1),min(p%k+1,Wnz+1)),&
                              Pr(min(p%i+1,Unx+1),min(p%j+1,Vny+1),min(p%k+1,Wnz+1)))

        if (enable_buoyancy==1) then
          temptime(k,step)=Trilinint((p%x-xPr(p%i))/(xPr(p%i+1)-xPr(p%i)),&
                             (p%y-yPr(p%j))/(yPr(p%j+1)-yPr(p%j)),&
                             (p%z-zPr(p%k))/(zPr(p%k+1)-zPr(p%k)),&
                             Temperature(p%i,p%j,p%k),&
                             Temperature(p%i+1,p%j,p%k),&
                             Temperature(p%i,p%j+1,p%k),&
                             Temperature(p%i,p%j,p%k+1),&
                             Temperature(p%i+1,p%j+1,p%k),&
                             Temperature(p%i+1,p%j,p%k+1),&
                             Temperature(p%i,p%j+1,p%k+1),&
                             Temperature(p%i+1,p%j+1,p%k+1))
        end if

        if (enable_moisture==1) then
          moisttime(k,step)=Trilinint((p%x-xPr(p%i))/(xPr(p%i+1)-xPr(p%i)),&
                             (p%y-yPr(p%j))/(yPr(p%j+1)-yPr(p%j)),&
                             (p%z-zPr(p%k))/(zPr(p%k+1)-zPr(p%k)),&
                             Moisture(p%i,p%j,p%k),&
                             Moisture(p%i+1,p%j,p%k),&
                             Moisture(p%i,p%j+1,p%k),&
                             Moisture(p%i,p%j,p%k+1),&
                             Moisture(p%i+1,p%j+1,p%k),&
                             Moisture(p%i+1,p%j,p%k+1),&
                             Moisture(p%i,p%j+1,p%k+1),&
                             Moisture(p%i+1,p%j+1,p%k+1))
        end if
      end associate
    end do

    do k = 1,number_of_scalar_probes    
      associate (p=> scalar_probes(k))
        do l = 1,num_of_scalars
          scalptime(l,k,step)=Trilinint((p%x-xPr(p%i))/(xPr(p%i+1)-xPr(p%i)),&
                             (p%y-yPr(p%j))/(yPr(p%j+1)-yPr(p%j)),&
                             (p%z-zPr(p%k))/(zPr(p%k+1)-zPr(p%k)),&
                             Scalar(p%i,p%j,p%k,l),&
                             Scalar(p%i+1,p%j,p%k,l),&
                             Scalar(p%i,p%j+1,p%k,l),&
                             Scalar(p%i,p%j,p%k+1,l),&
                             Scalar(p%i+1,p%j+1,p%k,l),&
                             Scalar(p%i+1,p%j,p%k+1,l),&
                             Scalar(p%i,p%j+1,p%k+1,l),&
                             Scalar(p%i+1,p%j+1,p%k+1,l))
        end do
      end associate
    end do

    if (store%tke==1) then
      tke(step)=totke(U,V,W)
    end if

    if (store%tke==1.and.store%dissip==1.and.step>0) then
      dissip(step)=(tke(step-1)-tke(step))/(times(step)-times(step-1))
    end if

    if (store%deltime==1) then
      deltime(step)=delta/dt
    end if

    endstep = step






    if (display%delta==1) then
      write (*,*) "delta: ",delta
    end if





    time_weight = dt/(timeavg2-timeavg1)

    if ((averaging==1).and.((time>=timeavg1).and.(time<=timeavg2))) then

      if (store%avg_U==1) then
        U_avg = U_avg + U * time_weight
        V_avg = V_avg + V * time_weight
        W_avg = W_avg + W * time_weight
      end if

      if (store%avg_Pr==1) then
        Pr_avg = Pr_avg + Pr(1:Prnx,1:Prny,1:Prnz) * time_weight
      end if

      if (enable_buoyancy==1.and.store%avg_temperature==1) then
        Temperature_avg = Temperature_avg + temperature * time_weight
      end if

      if (enable_moisture==1.and.store%avg_moisture==1) then
        Moisture_avg = Moisture_avg + moisture * time_weight
      end if

      if (num_of_scalars>0.and.store%scalarsavg==1) then
        Scalar_avg = Scalar_avg + Scalar * time_weight
      end if

   end if


   if (wallmodeltype>0.and.(display%ustar==1.or.store%ustar==1)) then

     S = GroundUstar()

     S2 = S*Re

     if (display%ustar==1) then
       if (allocated(Ustar_inlet)) then
        write(*,*) "ustar:",S,"Re_tau:",S2,"u*inlet",Ustar_inlet(1)
       else
        write(*,*) "ustar:",S,"Re_tau:",S2
       end if
     end if
     if (store%ustar==1) then
       ustar(:,step)=(/ S2 , S /)
     end if
   end if


    if (wallmodeltype>0.and.enable_buoyancy==1.and.TempBtype(Bo)==DIRICHLET.and.(display%tstar==1.or.store%tstar==1)) then
     S2 = SUM(BsideTFLArr(1:Prnx,1:Prny))/(Prnx*Prny)
     S=-S*S2

     if (display%tstar==1) then
       write(*,*) "Tstar",S,"tflux", S2
     end if

     if (store%tstar==1) then
       tstar(:,step) = (/ S2,S /)
     end if
    end if


    if (store%BLprofiles==1) then
      if ((averaging==1).and.((time>=timeavg1).and.(time<=timeavg2))) then

       call BLProfiles(U,V,W,Temperature,Moisture,Scalar)

       profuavg = profuavg + profu * time_weight
       profvavg = profvavg + profv * time_weight
       profuwavg = profuwavg + profuw * time_weight
       profuwsgsavg = profuwsgsavg + profuwsgs * time_weight
       profvwavg = profvwavg + profvw * time_weight
       profvwsgsavg = profvwsgsavg + profvwsgs * time_weight
       profuuavg = profuuavg + profuu * time_weight
       profvvavg = profvvavg + profvv * time_weight
       profwwavg = profwwavg + profww * time_weight

       if (enable_buoyancy==1) then
         proftempavg = proftempavg + proftemp * time_weight
         proftempflavg = proftempflavg + proftempfl * time_weight
         proftempflsgsavg = proftempflsgsavg + proftempflsgs * time_weight
         profttavg = profttavg + proftt * time_weight
       end if

       if (enable_buoyancy==1) then
         profmoistavg = profmoistavg + profmoist * time_weight
         profmoistflavg = profmoistflavg + profmoistfl * time_weight
         profmoistflsgsavg = profmoistflsgsavg + profmoistflsgs * time_weight
         profmmavg = profmmavg + profmm * time_weight
       end if

       if (num_of_scalars>0) then
         profscalavg = profscalavg + profscal * time_weight
         profscalflavg = profscalflavg + profscalfl * time_weight
         profscalflsgsavg = profscalflsgsavg + profscalflsgs * time_weight
         profssavg = profssavg + profss * time_weight
       end if
     end if
  end if



    if (frames>0)then
       if ((time>=timefram1).and.(time<=timefram2+(timefram2-timefram1)/(frames-1))&
         .and.(time>=timefram1+fnum*(timefram2-timefram1)/(frames-1))) then
        fnum = fnum+1
        call Frame(U,V,W,Pr,Temperature,Moisture,Scalar,fnum)
       end if
    end if

    if (allocated(StaggeredFrameDomains)) then
      do i=1,size(StaggeredFrameDomains)
        call Save(StaggeredFrameDomains(i), time, U, V, W, Pr, Temperature, Scalar)
      end do
    end if




  end subroutine OutTstep





  subroutine OutputProfiles
    character(2) :: prob
    real(knd) :: S,S2,nom,denom
    integer :: i,j,k,unit

    call newunit(unit)

    do k = 1,number_of_probes

      write(prob,"(i2.2)") k

      open(unit,file="output/Prtimep"//prob//".txt")
      do j = 0,endstep
       write (unit,*) times(j),Prtime(k,j)
      end do
      close(unit)

      open(unit,file="output/Utimep"//prob//".txt")
      do j = 0,endstep
       write (unit,*) times(j),Utime(k,j)
      end do
      close(unit)

      open(unit,file="output/Vtimep"//prob//".txt")
      do j = 0,endstep
       write (unit,*) times(j),Vtime(k,j)
      end do
      close(unit)

      open(unit,file="output/Wtimep"//prob//".txt")
      do j = 0,endstep
       write (unit,*) times(j),Wtime(k,j)
      end do
      close(unit)

      if (enable_buoyancy==1) then
        open(unit,file="output/temptimep"//prob//".txt")
        do j = 0,endstep
         write (unit,*) times(j),temptime(k,j)
        end do
        close(unit)
      end if

      if (enable_moisture==1) then
        open(unit,file="output/moisttimep"//prob//".txt")
        do j = 0,endstep
         write (unit,*) times(j),moisttime(k,j)
        end do
        close(unit)
      end if

    end do

    do k = 1,number_of_scalar_probes

      write(prob,"(i2.2)") k

      if (num_of_scalars>0) then
        open(unit,file="output/scaltimep"//prob//".txt")
        do j = 1,endstep
         write (unit,'(99(g0,tr3))') times(j),scalptime(:,k,j)
        end do
        close(unit)
      end if

    end do

    if (store%deltime==1) then
      open(unit,file="output/deltime.txt")
      do j = 1,endstep
       write (unit,*) times(j),deltime(j)
      end do
      close(unit)
    end if

    if (store%tke==1) then
      open(unit,file="output/tke.txt")
      do j = 0,endstep
       write (unit,*) times(j),tke(j)
      end do
      close(unit)
    end if

    if (store%tke==1.and.store%dissip==1) then
     open(unit,file="output/dissip.txt")
      do j = 1,endstep
       write (unit,*) times(j),dissip(j)
      end do
      close(unit)
    end if

    if (wallmodeltype>0.and.display%ustar==1) then
      open(unit,file="output/Retau.txt")
      do j = 0,endstep
       write (unit,*) times(j),ustar(:,j)
      end do
      close(unit)
    end if

    if (wallmodeltype>0.and.enable_buoyancy==1.and.TempBtype(Bo)==DIRICHLET.and.store%tstar==1) then
      open(unit,file="output/tflux.txt")
      do j = 0,endstep
       write (unit,*) times(j),tstar(:,j)
      end do
      close(unit)
    end if


    if (num_of_scalars>0.and.store%scalsumtime==1) then
      open(unit,file="output/scalsumtime.txt")
      do j = 1,endstep
       write (unit,*) times(j),scalsumtime(:,j)
      end do
      close(unit)
    end if

    if (num_of_scalars>0.and.store%scaltotsumtime==1) then
      open(unit,file="output/scaltotsumtime.txt")
      do j = 1,endstep
       write (unit,*) times(j),sum(scalsumtime(:,j))
      end do
      close(unit)
    end if


    if (store%BLprofiles==1.and.averaging==1) then

       open(unit,file="output/profu.txt")
       do k = 1,Unz
        write (unit,*) zPr(k),profuavg(k)
       end do
       close(unit)

       open(unit,file="output/profv.txt")
       do k = 1,Vnz
        write (unit,*) zPr(k),profvavg(k)
       end do
       close(unit)

       open(unit,file="output/profuu.txt")
       do k = 1,Unz
        write (unit,*) zPr(k),profuuavg(k)
       end do
       close(unit)

       open(unit,file="output/profvv.txt")
       do k = 1,Vnz
        write (unit,*) zPr(k),profvvavg(k)
       end do
       close(unit)

       open(unit,file="output/profww.txt")
       do k = 1,Wnz
        write (unit,*) zW(k),profwwavg(k)
       end do
       close(unit)

       open(unit,file="output/profuw.txt")
       do k = 0,Prnz
        write (unit,*) zW(k),profuwavg(k),profuwsgsavg(k)
       end do
       close(unit)

       open(unit,file="output/profvw.txt")
       do k = 0,Prnz
        write (unit,*) zW(k),profvwavg(k),profvwsgsavg(k)
       end do
       close(unit)

       if (enable_buoyancy==1) then

          open(unit,file="output/proftemp.txt")
          do k = 1,Prnz
           write (unit,*) zPr(k),proftempavg(k)
          end do
          close(unit)

          open(unit,file="output/proftempfl.txt")
          do k = 0,Prnz
           write (unit,*) zW(k),proftempflavg(k),proftempflsgsavg(k)
          end do
          close(unit)

          open(unit,file="output/proftt.txt")
          do k = 1,Prnz
           write (unit,*) zPr(k),profttavg(k)
          end do
          close(unit)

          if (allocated(U_avg).and.allocated(V_avg).and.allocated(Temperature_avg)) then
            open(unit,file="output/profRig.txt")
            do k = 1,Prnz
             S = 0
             do j = 1,Prny
              do i = 1,Prnx
               S = S + Rig(i,j,k,U_avg,V_avg,Temperature_avg)
              end do
             end do
             S = S / (Prnx*Prny)
             write (unit,*) zPr(k),S
            end do
            close(unit)

            open(unit,file="output/profRf.txt")
            do k = 1,Prnz
             S = 0
             S2 = 0
             do j = 1,Prny
              do i = 1,Prnx
               S = S + (U_avg(i,j,k+1)+U_avg(i-1,j,k+1)-U_avg(i,j,k-1)-U_avg(i-1,j,k-1))/(2._knd*(zPr(k+1)-zPr(k-1)))
               S2 = S2+(V_avg(i,j,k+1)+V_avg(i,j-1,k+1)-V_avg(i,j,k-1)-V_avg(i,j-1,k-1))/(2._knd*(zPr(k+1)-zPr(k-1)))
              end do
             end do

             S = S / (Prnx*Prny)
             S2 = S2/(Prnx*Prny)
             nom=(grav_acc/temperature_ref)*(proftempflavg(k)+proftempflsgsavg(k))
             denom=((profuwavg(k)+profuwsgsavg(k))*S + (profvwavg(k)+profvwsgsavg(k))*S2)

             if (abs(denom)>1E-5_knd*abs(nom).and.abs(nom)>epsilon(1._knd)) then
               S = nom/denom
             else
               S = 100000._knd*sign(1.0_knd,nom)*sign(1.0_knd,denom)
             end if

             write (unit,*) zPr(k),S
            end do
            close(unit)
          end if !allocated avgs

       end if !enable_buoyancy == 1


       if (enable_moisture==1) then

          open(unit,file="output/profmoist.txt")
          do k = 1,Prnz
           write (unit,*) zPr(k),profmoistavg(k)
          end do
          close(unit)

          open(unit,file="output/profmoistfl.txt")
          do k = 0,Prnz
           write (unit,*) zW(k),profmoistflavg(k),profmoistflsgsavg(k)
          end do
          close(unit)

          open(unit,file="output/profmm.txt")
          do k = 1,Prnz
           write (unit,*) zPr(k),profmmavg(k)
          end do
          close(unit)

       end if !enable_moisture == 1


       if (num_of_scalars>0) then

          open(unit,file="output/profscal.txt")
          do k = 1,Prnz
           write (unit,*) zPr(k),(profscalavg(i,k), i = 1,num_of_scalars)
          end do
          close(unit)

          open(unit,file="output/profscalfl.txt")
          do k = 0,Prnz
           write (unit,*) zW(k),(profscalflavg(i,k),profscalflsgsavg(i,k), i = 1,num_of_scalars)
          end do
          close(unit)

          open(unit,file="output/profss.txt")
          do k = 1,Prnz
           write (unit,*) zPr(k),(profssavg(i,k), i = 1,num_of_scalars)
          end do
          close(unit)

       end if !num_of_scalars>0

    end if !store%BLprofiles

  end subroutine OutputProfiles







  subroutine OutputOut(U,V,W,Pr,Temperature,Moisture)
    real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    real(knd),dimension(1:,1:,1:),intent(in) :: Pr
    real(knd),intent(in) :: Temperature(-1:,-1:,-1:)
    real(knd),intent(in) :: Moisture(-1:,-1:,-1:)
    character(70) :: str
    integer i,j,k,unit
    real(real32),allocatable :: tmp(:,:,:,:)

    if (store%out==1) then

       call newunit(unit);

       open(unit,file="output/out.vtk", &
         access='stream',status='replace',form="unformatted",action="write")

       write (unit) "# vtk DataFile Version 2.0",lf
       write (unit) "CLMM output file",lf
       write (unit) "BINARY",lf
       write (unit) "DATASET RECTILINEAR_GRID",lf
       str="DIMENSIONS"
       write (str(12:),*) Prnx,Prny,Prnz
       write (unit) str,lf
       str="X_COORDINATES"
       write (str(15:),'(i5,2x,a)') Prnx,"float"
       write (unit) str,lf
       write (unit) BigEnd(real(xPr(1:Prnx), real32)),lf
       str="Y_COORDINATES"
       write (str(15:),'(i5,2x,a)') Prny,"float"
       write (unit) str,lf
       write (unit) BigEnd(real(yPr(1:Prny), real32)),lf
       str="Z_COORDINATES"
       write (str(15:),'(i5,2x,a)') Prnz,"float"
       write (unit) str,lf
       write (unit) BigEnd(real(zPr(1:Prnz), real32)),lf
       str="POINT_DATA"
       write (str(12:),*) Prnx*Prny*Prnz
       write (unit) str,lf

       if (store%out_Pr==1) then
         write (unit) "SCALARS p float",lf
         write (unit) "LOOKUP_TABLE default",lf

         write (unit) BigEnd(real(Pr(1:Prnx,1:Prny,1:Prnz), real32))

         write (unit) lf
       end if

       if (enable_buoyancy>0.and.store%out_temperature==1) then
         write (unit) "SCALARS temperature float",lf
         write (unit) "LOOKUP_TABLE default",lf

         write (unit) BigEnd(real(Temperature(1:Prnx,1:Prny,1:Prnz), real32))

         write (unit) lf
       end if

       if (enable_moisture>0.and.store%out_moisture==1) then
         write (unit) "SCALARS moisture float",lf
         write (unit) "LOOKUP_TABLE default",lf

         write (unit) BigEnd(real(Moisture(1:Prnx,1:Prny,1:Prnz), real32))

         write (unit) lf
       end if

       if (store%out_Prtype==1) then
         write (unit) "SCALARS ptype float",lf
         write (unit) "LOOKUP_TABLE default",lf

         write (unit) BigEnd(real(Prtype(1:Prnx,1:Prny,1:Prnz), real32))

         write (unit) lf
       end if

       if (store%out_div==1) then
         write (unit) "SCALARS div float",lf
         write (unit) "LOOKUP_TABLE default",lf

         if (gridtype==uniformgrid) then
           write (unit) BigEnd(real((U(1:Prnx,1:Prny,1:Prnz) - &
                                       U(0:Prnx-1,1:Prny,1:Prnz))/dxmin + &
                                    (V(1:Prnx,1:Prny,1:Prnz) - &
                                       V(1:Vnx,0:Vny-1,1:Vnz))/dymin + &
                                    (W(1:Prnx,1:Prny,1:Prnz) - &
                                       W(1:Prnx,1:Prny,0:Prnz-1))/dzmin &
                               , real32))
         else
           do k = 1,Prnz
            do j = 1,Prny
             do i = 1,Prnx
               write (unit) BigEnd(real((U(i,j,k)-U(i-1,j,k))/(dxPr(i)) + &
                                        (V(i,j,k)-V(i,j-1,k))/(dyPr(j)) + &
                                        (W(i,j,k)-W(i,j,k-1))/(dzPr(k)), real32))
             end do
            end do
           end do
         end if

         write (unit) lf
       end if

       if (store%out_lambda2==1) then
         allocate(tmp(1:Prnx,1:Prny,1:Prnz,1))

         write (unit) "SCALARS lambda2 float",lf
         write (unit) "LOOKUP_TABLE default",lf

         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
             tmp(i,j,k,1) = BigEnd(real(Lambda2(i,j,k,U,V,W), real32))
           end do
          end do
         end do

         write (unit) tmp

         write (unit) lf
         deallocate(tmp)
       end if

       if (store%out_visc==1) then
         write (unit) "SCALARS visc float",lf
         write (unit) "LOOKUP_TABLE default",lf

         write (unit) BigEnd(real(Visc(1:Prnx,1:Prny,1:Prnz), real32))

         write (unit) lf
       end if

       if (store%out_U==1.or.store%out_vort==1) allocate(tmp(1:3,1:Prnx,1:Prny,1:Prnz))

       if (store%out_U==1) then
         write (unit) "VECTORS u float",lf

         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
             tmp(1:3,i,j,k) = [ BigEnd(real((U(i,j,k)+U(i-1,j,k))/2._knd, real32)), &
                                BigEnd(real((V(i,j,k)+V(i,j-1,k))/2._knd, real32)), &
                                BigEnd(real((W(i,j,k)+W(i,j,k-1))/2._knd, real32)) ]
           end do
          end do
         end do

         write (unit) tmp

         write (unit) lf
       end if

       if (store%out_vort==1) then
         write (unit) "VECTORS vort float",lf

         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
             tmp(1:3,i,j,k) = [ BigEnd(real((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin) &
                                          -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin), real32)), &
                                BigEnd(real((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin) &
                                          -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin), real32)), &
                                BigEnd(real((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin) &
                                          -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin), real32)) ]
           end do
          end do
         end do

         write (unit) tmp

         write (unit) lf
       end if
       close(unit)
    end if  !store%out
  end subroutine OutputOut








  subroutine OutputScalars(Scalar)
    real(knd),dimension(-1:,-1:,-1:,:),intent(in) :: Scalar
    character(70) :: str
    real(knd),dimension(:,:,:),allocatable :: depos
    character(8) ::  scalname="scalar00"
    integer :: l,unit

    if (num_of_scalars>0) then
      if (store%scalars==1) then

          call newunit(unit)

          open(unit,file="output/scalars.vtk",&
            access='stream',status='replace',form="unformatted",action="write")

          write (unit) "# vtk DataFile Version 2.0",lf
          write (unit) "CLMM output file",lf
          write (unit) "BINARY",lf
          write (unit) "DATASET RECTILINEAR_GRID",lf
          str="DIMENSIONS"
          write (str(12:),*) Prnx,Prny,Prnz
          write (unit) str,lf
          str="X_COORDINATES"
          write (str(15:),'(i5,2x,a)') Prnx,"float"
          write (unit) str,lf
          write (unit) BigEnd(real(xPr(1:Prnx), real32)),lf
          str="Y_COORDINATES"
          write (str(15:),'(i5,2x,a)') Prny,"float"
          write (unit) str,lf
          write (unit) BigEnd(real(yPr(1:Prny), real32)),lf
          str="Z_COORDINATES"
          write (str(15:),'(i5,2x,a)') Prnz,"float"
          write (unit) str,lf
          write (unit) BigEnd(real(zPr(1:Prnz), real32)),lf
          str="POINT_DATA"
          write (str(12:),*) Prnx*Prny*Prnz
          write (unit) str,lf

          do l = 1,num_of_scalars
            write(scalname(7:8),"(I2.2)") l
            write (unit) "SCALARS ", scalname , " float",lf
            write (unit) "LOOKUP_TABLE default",lf

            write (unit) BigEnd(real(Scalar(1:Prnx,1:Prny,1:Prnz,l), real32))

            write (unit) lf
          end do
          close(unit)
      end if !store%scalars

      if (computedeposition>0.and.store%deposition==1) then

          allocate(depos(1:Prnx,1:Prny,num_of_scalars))
          depos = GroundDeposition()

          call newunit(unit)

          open(unit,file="output/deposition.vtk",&
            access='stream',status='replace',form="unformatted",action="write")

          write (unit) "CLMM output file",lf
          write (unit) "BINARY",lf
          write (unit) "DATASET RECTILINEAR_GRID",lf
          str="DIMENSIONS"
          write (str(12:),*) Prnx,Prny,1
          write (unit) str,lf
          str="X_COORDINATES"
          write (str(15:),'(i5,2x,a)') Prnx,"float"
          write (unit) str,lf
          write (unit) BigEnd(real(xPr(1:Prnx), real32)),lf
          str="Y_COORDINATES"
          write (str(15:),'(i5,2x,a)') Prny,"float"
          write (unit) str,lf
          write (unit) BigEnd(real(yPr(1:Prny), real32))
          str="Z_COORDINATES"
          write (str(15:),'(i5,2x,a)') 1,"float"
          write (unit) str,lf
          write (unit) BigEnd(real(zW(0), real32)),lf
          str="POINT_DATA"
          write (str(12:),*) Prnx*Prny
          write (unit) str,lf

          do l = 1,num_of_scalars
            write(scalname(7:8),"(I2.2)") l
            write (unit) "SCALARS ", scalname , " float",lf
            write (unit) "LOOKUP_TABLE default",lf

            write (unit) BigEnd(real(depos(1:Prnx,1:Prny,l), real32))

            write (unit) lf
          end do

          close(unit)

          deallocate(depos)

      end if  !store%deposition
    end if  !num_of_scalars
  end subroutine OutputScalars










  subroutine OutputAvg(U,V,W,Pr,Temperature,Moisture,Scalar)
    real(knd),dimension(-2:,-2:,-2:),intent(in)   :: U,V,W
    real(knd),dimension(1:,1:,1:),intent(in)      :: Pr
    real(knd),dimension(-1:,-1:,-1:),intent(in)   :: Temperature
    real(knd),dimension(-1:,-1:,-1:),intent(in)   :: Moisture
    real(knd),dimension(-1:,-1:,-1:,:),intent(in) :: Scalar
    character(70) :: str
    character(8) ::  scalname="scalar00"
    integer i,j,k,l,unit
    real(real32),allocatable :: tmp(:,:,:,:)

    if (averaging==1.and.btest(store%avg,0)) then
        allocate(tmp(1:3,1:Prnx,1:Prny,1:Prnz))

        call newunit(unit)

        open(unit,file="output/avg.vtk",&
          access='stream',status='replace',form="unformatted",action="write")

        write (unit) "# vtk DataFile Version 2.0",lf
        write (unit) "CLMM output file",lf
        write (unit) "BINARY",lf
        write (unit) "DATASET RECTILINEAR_GRID",lf
        str="DIMENSIONS"
        write (str(12:),*) Prnx,Prny,Prnz
        write (unit) str,lf
        str="X_COORDINATES"
        write (str(15:),'(i5,2x,a)') Prnx,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(xPr(1:Prnx), real32)),lf
        str="Y_COORDINATES"
        write (str(15:),'(i5,2x,a)') Prny,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(yPr(1:Prny), real32)),lf
        str="Z_COORDINATES"
        write (str(15:),'(i5,2x,a)') Prnz,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(zPr(1:Prnz), real32)),lf
        str="POINT_DATA"
        write (str(12:),*) Prnx*Prny*Prnz
        write (unit) str,lf

        if (store%avg_Pr==1) then
          write (unit) "SCALARS p float",lf
          write (unit) "LOOKUP_TABLE default",lf

          write (unit) BigEnd(real(Pr(1:Prnx,1:Prny,1:Prnz), real32))

          write (unit) lf
        end if

        if (store%avg_Prtype==1) then
          write (unit) "SCALARS ptype float",lf
          write (unit) "LOOKUP_TABLE default",lf

          write (unit) BigEnd(real(Prtype(1:Prnx,1:Prny,1:Prnz), real32))

          write (unit) lf
        end if

        if (enable_buoyancy==1.and.store%avg_temperature==1) then
          write (unit) "SCALARS temperature float",lf
          write (unit) "LOOKUP_TABLE default",lf

          write (unit) BigEnd(real(Temperature(1:Prnx,1:Prny,1:Prnz), real32))

          write (unit) lf
        end if

        if (enable_moisture==1.and.store%avg_moisture==1) then
          write (unit) "SCALARS moisture float",lf
          write (unit) "LOOKUP_TABLE default",lf

          write (unit) BigEnd(real(Moisture(1:Prnx,1:Prny,1:Prnz), real32))

          write (unit) lf
        end if

        if (store%avg_U==1) then
          write (unit) "VECTORS u float",lf

          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              tmp(:,i,j,k) = [ BigEnd(real((U(i,j,k)+U(i-1,j,k))/2._knd, real32)), &
                               BigEnd(real((V(i,j,k)+V(i,j-1,k))/2._knd, real32)), &
                               BigEnd(real((W(i,j,k)+W(i,j,k-1))/2._knd, real32)) ]
            end do
           end do
          end do

          write (unit) tmp

          write (unit) lf
        end if

        if (store%avg_vort==1) then
          write (unit) "VECTORS vort float",lf

          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              tmp(:,i,j,k) = [ BigEnd(real((W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin) &
                                         -(V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin), real32)), &
                               BigEnd(real((U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin) &
                                         -(W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin), real32)), &
                               BigEnd(real((V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin) &
                                         -(U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin), real32)) ]
            end do
           end do
          end do

          write (unit) tmp

          write (unit) lf
        end if

        close(unit)
    end if !averaging

    if (averaging==1.and.btest(store%avg,1)) call OutputUVW(U,V,W,"output/Uavg.vtk","output/Vavg.vtk","output/Wavg.vtk",.true.)

    if (averaging==1.and.store%scalarsavg==1) then
        if (num_of_scalars>0) then
            open(unit,file="output/scalarsavg.vtk",&
              access='stream',status='replace',form="unformatted",action="write")

            write (unit) "# vtk DataFile Version 2.0",lf
            write (unit) "CLMM output file",lf
            write (unit) "BINARY",lf
            write (unit) "DATASET RECTILINEAR_GRID",lf
            str="DIMENSIONS"
            write (str(12:),*) Prnx,Prny,Prnz
            write (unit) str,lf
            str="X_COORDINATES"
            write (str(15:),'(i5,2x,a)') Prnx,"float"
            write (unit) str,lf
            write (unit) BigEnd(real(xPr(1:Prnx), real32)),lf
            str="Y_COORDINATES"
            write (str(15:),'(i5,2x,a)') Prny,"float"
            write (unit) str,lf
            write (unit) BigEnd(real(yPr(1:Prny), real32)),lf
            str="Z_COORDINATES"
            write (str(15:),'(i5,2x,a)') Prnz,"float"
            write (unit) str,lf
            write (unit) BigEnd(real(zPr(1:Prnz), real32)),lf
            str="POINT_DATA"
            write (str(12:),*) Prnx*Prny*Prnz
            write (unit) str,lf

            do l = 1,num_of_scalars
              write(scalname(7:8),"(I2.2)") l
              write (unit) "SCALARS ", scalname , " float",lf
              write (unit) "LOOKUP_TABLE default",lf

              write (unit) BigEnd(real(Scalar(1:Prnx,1:Prny,1:Prnz,l), real32))

              write (unit) lf
            end do
            close(unit)
        end if  !num_of_scalars

    end if !averaging
  end subroutine OutputAvg




  subroutine OutputUVW(U,V,W,fnameU,fnameV,fnameW,avg_mode_arg)
    real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    character(*),intent(in) :: fnameU,fnameV,fnameW
    logical,optional,intent(in) :: avg_mode_arg
    character(70) :: str
    integer :: i,unit
    logical :: avg_mode

    if (present(avg_mode_arg)) then
      avg_mode = avg_mode_arg
    else
      avg_mode = .false.
    end if

    if (store%U==1.or.avg_mode) then

        call newunit(unit)

        open(unit,file=fnameU,&
          access='stream',status='replace',form="unformatted",action="write")

        write (unit) "# vtk DataFile Version 2.0",lf
        write (unit) "CLMM output file",lf
        write (unit) "BINARY",lf
        write (unit) "DATASET RECTILINEAR_GRID",lf
        str="DIMENSIONS"
        write (str(12:),*) Unx,Uny,Unz
        write (unit) str,lf
        str="X_COORDINATES"
        write (str(15:),'(i5,2x,a)') Unx,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(xU(1:Unx), real32)),lf
        str="Y_COORDINATES"
        write (str(15:),'(i5,2x,a)') Uny,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(yPr(1:Uny), real32)),lf
        str="Z_COORDINATES"
        write (str(15:),'(i5,2x,a)') Unz,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(zPr(1:Unz), real32)),lf
        str="POINT_DATA"
        write (str(12:),*) Unx*Uny*Unz
        write (unit) str,lf


        write (unit) "SCALARS U float",lf
        write (unit) "LOOKUP_TABLE default",lf

        write (unit) BigEnd(real(U(1:Unx,1:Uny,1:Unz), real32))

        write (unit) lf

        if (store%U_interp==1) then
          write (unit) "SCALARS Utype float",lf
          write (unit) "LOOKUP_TABLE default",lf

          write (unit) BigEnd(real(Utype(1:Unx,1:Uny,1:Unz), real32))

          write (unit) lf
        end if

        close(unit)
    end if

    if (store%V==1.or.avg_mode) then

        call newunit(unit)

        open(unit,file=fnameV,&
          access='stream',status='replace',form="unformatted",action="write")

        write (unit) "# vtk DataFile Version 2.0",lf
        write (unit) "CLMM output file",lf
        write (unit) "BINARY",lf
        write (unit) "DATASET RECTILINEAR_GRID",lf
        str="DIMENSIONS"
        write (str(12:),*) Vnx,Vny,Vnz
        write (unit) str,lf
        str="X_COORDINATES"
        write (str(15:),'(i5,2x,a)') Vnx,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(xPr(1:Vnx), real32)),lf
        str="Y_COORDINATES"
        write (str(15:),'(i5,2x,a)') Vny,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(yV(1:Vny), real32)),lf
        str="Z_COORDINATES"
        write (str(15:),'(i5,2x,a)') Vnz,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(zPr(1:Vnz), real32)),lf
        str="POINT_DATA"
        write (str(12:),*) Vnx*Vny*Vnz
        write (unit) str,lf


        write (unit) "SCALARS V float",lf
        write (unit) "LOOKUP_TABLE default",lf

        write (unit) BigEnd(real(V(1:Vnx,1:Vny,1:Vnz), real32))

        write (unit) lf

        if (store%V_interp==1) then
          write (unit) "SCALARS Vtype float",lf
          write (unit) "LOOKUP_TABLE default",lf

          write (unit) BigEnd(real(Vtype(1:Vnx,1:Vny,1:Vnz), real32))

          write (unit) lf
        end if

        close(unit)
    end if !store%V

    if (store%W==1.or.avg_mode) then

        call newunit(unit)

        open(unit,file=fnameW,&
          access='stream',status='replace',form="unformatted",action="write")

        write (unit) "# vtk DataFile Version 2.0",lf
        write (unit) "CLMM output file",lf
        write (unit) "BINARY",lf
        write (unit) "DATASET RECTILINEAR_GRID",lf
        str="DIMENSIONS"
        write (str(12:),*) Wnx,Wny,Wnz
        write (unit) str,lf
        str="X_COORDINATES"
        write (str(15:),'(i5,2x,a)') Wnx,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(xPr(1:Wnx), real32)),lf
        str="Y_COORDINATES"
        write (str(15:),'(i5,2x,a)') Wny,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(yPr(1:Wny), real32)),lf
        str="Z_COORDINATES"
        write (str(15:),'(i5,2x,a)') Wnz,"float"
        write (unit) str,lf
        write (unit) BigEnd(real(zW(1:Wnz), real32)),lf
        str="POINT_DATA"
        write (str(12:),*) Wnx*Wny*Wnz
        write (unit) str,lf


        write (unit) "SCALARS W float",lf
        write (unit) "LOOKUP_TABLE default",lf

        write (unit) BigEnd(real(W(1:Wnx,1:Wny,1:Wnz), real32))

        write (unit) lf

        if (store%W_interp==1) then
          write (unit) "SCALARS Wtype float",lf
          write (unit) "LOOKUP_TABLE default",lf

          write (unit) BigEnd(real(Wtype(1:Wnx,1:Wny,1:Wnz), real32))

          write (unit) lf
        end if

        close(unit)
    end if !store%W

    if (store%U_interp/=0 .and. .not.avg_mode) then

      call newunit(unit)

      open(unit,file="output/Uinterp.txt")
      do i = 1,size(UIBPoints)
        write(unit,*) "xi,yj,zk",UIBPoints(i)%xi,UIBPoints(i)%yj,UIBPoints(i)%zk
        write(unit,*) "interp",UIBPoints(i)%interp
        write(unit,*) "xi",UIBPoints(i)%IntPoints%xi
        write(unit,*) "yj",UIBPoints(i)%IntPoints%yj
        write(unit,*) "zk",UIBPoints(i)%IntPoints%zk
        write(unit,*) "coefs",UIBPoints(i)%IntPoints%coef
      end do
      close(unit)
    end if

    if (store%V_interp/=0 .and. .not.avg_mode) then

      call newunit(unit)

      open(unit,file="output/Vinterp.txt")
      do i = 1,size(VIBPoints)
        write(unit,*) "xi,yj,zk",VIBPoints(i)%xi,VIBPoints(i)%yj,VIBPoints(i)%zk
        write(unit,*) "interp",VIBPoints(i)%interp
        write(unit,*) "xi",VIBPoints(i)%IntPoints%xi
        write(unit,*) "yj",VIBPoints(i)%IntPoints%yj
        write(unit,*) "zk",VIBPoints(i)%IntPoints%zk
        write(unit,*) "coefs",VIBPoints(i)%IntPoints%coef
      end do
      close(unit)
    end if

    if (store%W_interp/=0) then

      call newunit(unit)

      open(unit,file="output/Winterp.txt")
      do i = 1,size(WIBPoints)
        write(unit,*) "xi,yj,zk",WIBPoints(i)%xi,WIBPoints(i)%yj,WIBPoints(i)%zk
        write(unit,*) "interp",WIBPoints(i)%interp
        write(unit,*) "xi",WIBPoints(i)%IntPoints%xi
        write(unit,*) "yj",WIBPoints(i)%IntPoints%yj
        write(unit,*) "zk",WIBPoints(i)%IntPoints%zk
        write(unit,*) "coefs",WIBPoints(i)%IntPoints%coef
      end do
      close(unit)
    end if

    if (store%U_interp/=0.and.store%V_interp/=0.and.store%W_interp/=0) then

      call newunit(unit)

      open(unit,file="output/Scinterp.txt")
      do i = 1,size(ScalFlIBPoints)
        write(unit,*) "xi,yj,zk",ScalFlIBPoints(i)%xi,ScalFlIBPoints(i)%yj,ScalFlIBPoints(i)%zk
        write(unit,*) "interp",ScalFlIBPoints(i)%interp
        write(unit,*) "xi",ScalFlIBPoints(i)%IntPoints%xi
        write(unit,*) "yj",ScalFlIBPoints(i)%IntPoints%yj
        write(unit,*) "zk",ScalFlIBPoints(i)%IntPoints%zk
        write(unit,*) "coefs",ScalFlIBPoints(i)%IntPoints%coef
      end do
      close(unit)
    end if

  end subroutine OutputUVW






  subroutine Output(U,V,W,Pr,Temperature,Moisture,Scalar)
    real(knd),dimension(-2:,-2:,-2:),intent(inout) :: U,V,W
    real(knd),intent(inout) :: Pr(1:,1:,1:)
    real(knd),intent(in) :: Temperature(-1:,-1:,-1:)
    real(knd),intent(in) :: Moisture(-1:,-1:,-1:)
    real(knd),intent(in) :: Scalar(-1:,-1:,-1:,:)
    integer i

    call BoundU(1,U)
    call BoundU(2,V)
    call BoundU(3,W)

    call OutputOut(U,V,W,Pr,Temperature,Moisture)

    call OutputScalars(Scalar)

    if (time>=timeavg1) call OutputAvg(U_avg,V_avg,W_avg,Pr_avg,Temperature_avg,Moisture_avg,Scalar_avg)

    call OutputUVW(U,V,W,"U.vtk","V.vtk","W.vtk")

    if (time>=timeavg1) call OutputProfiles

    if (allocated(StaggeredFrameDomains)) then
      do i=1,size(StaggeredFrameDomains)
        call Finalize(StaggeredFrameDomains(i))
      end do
    end if

    write (*,*) "saved"
  end subroutine Output


  pure real(knd) function ScalarVerticalFlux(i,j,k,Scal,W)
    integer, intent(in)   :: i,j,k
    real(knd), intent(in) :: Scal(-1:,-1:,-1:), W(-2:,-2:,-2:)

    ScalarVerticalFlux = Scal(i,j,k) * (W(i,j,k)+W(i,j,k-1))/2 + &
                       TDiff(i,j,k) * (Scal(i,j,k+1)-Scal(i,j,k-1)) / (zW(k+1)-zW(k-1))
  end function ScalarVerticalFlux

  pure function Vorticity(i,j,k,U,V,W)
    real(knd), dimension(3) :: Vorticity
    integer, intent(in)     :: i,j,k
    real(knd),dimension(-2:,-2:,-2:), intent(in) :: U,V,W

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
                            (yV(Prny+1)+yV(0))/2._knd, (zW(Prnz+1)+zW(0))/2._knd )

            domain%maxi = domain%mini
            domain%minj = 1
            domain%maxj = Prny
            domain%mink = 1
            domain%maxk = Prnz

          elseif (domain%direction==2) then

            call Gridcoords(domain%mini, domain%minj, domain%mink, (xU(Prnx+1)+xU(0))/2._knd, &
                            domain%position, (zW(Prnz+1)+zW(0))/2._knd )

            domain%maxj = domain%minj
            domain%mini = 1
            domain%maxi = Prnx
            domain%mink = 1
            domain%maxk = Prnz

          else

            call Gridcoords(domain%mini, domain%minj, domain%mink, (xU(Prnx+1)+xU(0))/2._knd, &
                            (yV(Prny+1)+yV(0))/2._knd, domain%position )

            domain%maxk = domain%mink
            domain%mini = 1
            domain%maxi = Prnx
            domain%minj = 1
            domain%maxj = Prny

          end if

        end if

  end subroutine SetFrameDomain



  subroutine Frame(U,V,W,Pr,Temperature,Moisture,Scalar,n)
    real(knd),intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
    real(knd),intent(in) :: Temperature(-1:,-1:,-1:)
    real(knd),intent(in) :: Moisture(-1:,-1:,-1:)
    real(knd),intent(in) :: Scalar(-1:,-1:,-1:,:)
    integer,intent(in)   :: n
    integer i,j,k,l,m
    character(40) :: fname
    character(4) :: fsuffix,fnumber
    character(70) :: str
    character(8) ::  scalname="scalar00"
    integer mini,maxi,minj,maxj,mink,maxk
    integer unit
    real(real32),allocatable :: buffer(:,:,:),vbuffer(:,:,:,:)

    fsuffix=".vtk"


    !$omp parallel do default(private) shared(n,store,time,fsuffix,xPr,yPr,zPr,Prtype,Pr,U,V,W,&
    !$omp    scalar,Temperature,Moisture,enable_moisture,enable_buoyancy,num_of_scalars,&
    !$omp    dxmin,dymin,dzmin,temperature_ref)
    do m = 1,size(store%frame_domains)

      fname="output/frame-"//achar(iachar('a')+m-1)//"-"
      write(fnumber,"(I4.4)") n
      write(*,*) "Saving frame:",fnumber,"   time:",time


      mini = store%frame_domains(m)%mini
      minj = store%frame_domains(m)%minj
      mink = store%frame_domains(m)%mink
      maxi = store%frame_domains(m)%maxi
      maxj = store%frame_domains(m)%maxj
      maxk = store%frame_domains(m)%maxk

      allocate(buffer(mini:maxi,minj:maxj,mink:maxk))

      if (store%frame_U==1 .or. store%frame_vort==1) then
        allocate(vbuffer(3,mini:maxi,minj:maxj,mink:maxk))
      end if


      !$omp critical
      call newunit(unit)
      open(unit,file = trim(fname)//trim(fnumber)//trim(fsuffix),&
        access='stream',status='replace',form="unformatted",action="write")
      !$omp end critical

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
      write (unit) (BigEnd(real(xPr(i), real32)), i = mini,maxi),lf
      str="Y_COORDINATES"
      write (str(15:),'(i5,2x,a)') maxj-minj+1,"float"
      write (unit) str,lf
      write (unit) (BigEnd(real(yPr(j), real32)), j = minj,maxj),lf
      str="Z_COORDINATES"
      write (str(15:),'(i5,2x,a)') maxk-mink+1,"float"
      write (unit) str,lf
      write (unit) (BigEnd(real(zPr(k), real32)), k = mink,maxk),lf
      str="POINT_DATA"
      write (str(12:),*) (maxi-mini+1)*(maxj-minj+1)*(maxk-mink+1)
      write (unit) str,lf

      if (store%frame_Pr==1) then
        write (unit) "SCALARS p float",lf
        write (unit) "LOOKUP_TABLE default",lf
        
        write (unit) BigEnd(real(Pr(mini:maxi,minj:maxj,mink:maxk), real32))

        write (unit) lf
      end if

      if (store%frame_lambda2==1) then
        write (unit) "SCALARS lambda2 float",lf
        write (unit) "LOOKUP_TABLE default",lf
        do k = mink,maxk
         do j = minj,maxj
          do i = mini,maxi
            if (Prtype(i,j,k)<=0) then
              buffer(i,j,k) = real(Lambda2(i,j,k,U,V,W), real32)
            else
              buffer(i,j,k) = 0
            end if
          end do
         end do
        end do

        write (unit) BigEnd(buffer)

        write (unit) lf
      end if

      if (store%frame_scalars==1) then
        scalname(1:6)="scalar"
        do l = 1,num_of_scalars
         write(scalname(7:8),"(I2.2)") l
         write (unit) "SCALARS ", scalname , " float",lf
         write (unit) "LOOKUP_TABLE default",lf

         do k = mink,maxk
          do j = minj,maxj
           do i = mini,maxi
             if (Prtype(i,j,k)<=0) then
               buffer(i,j,k) = real(SCALAR(i,j,k,l), real32)
             else
               buffer(i,j,k) = 0
             end if
           end do
          end do
         end do

         write (unit) BigEnd(buffer)

         write (unit) lf
        end do
      elseif (store%frame_sumscalars==1.and.num_of_scalars==1) then
        write (unit) "SCALARS ", "scalar" , " float",lf
        write (unit) "LOOKUP_TABLE default",lf
        do k = mink,maxk
         do j = minj,maxj
          do i = mini,maxi
            if (Prtype(i,j,k)<=0) then
              buffer(i,j,k) =  real(SUM(SCALAR(i,j,k,:)), real32)
            else
              buffer(i,j,k) = 0._real32
            end if
          end do
         end do
        end do

        write (unit) BigEnd(buffer)

        write (unit) lf
      end if

      if (enable_buoyancy==1.and.store%frame_temperature==1) then
        write (unit) "SCALARS temperature float",lf
        write (unit) "LOOKUP_TABLE default",lf
        do k = mink,maxk
         do j = minj,maxj
          do i = mini,maxi
            if (Prtype(i,j,k)<=0) then
              buffer(i,j,k) = real(Temperature(i,j,k), real32)
            else
              buffer(i,j,k) = real(temperature_ref, real32)
            end if
          end do
         end do
        end do

        write (unit) BigEnd(buffer)

        write (unit) lf
      end if

      if (enable_moisture==1.and.store%frame_moisture==1) then
        write (unit) "SCALARS moisture float",lf
        write (unit) "LOOKUP_TABLE default",lf
        do k = mink,maxk
         do j = minj,maxj
          do i = mini,maxi
            if (Prtype(i,j,k)<=0) then
              buffer(i,j,k) = real(Moisture(i,j,k), real32)
            else
              buffer(i,j,k) = 0._real32
            end if
          end do
         end do
        end do

        write (unit) BigEnd(buffer)

        write (unit) lf
      end if

      if (enable_buoyancy==1.and.store%frame_temperature_flux==1) then
        write (unit) "SCALARS temperature_flux float", lf
        write (unit) "LOOKUP_TABLE default", lf

        do k = mink,maxk
         do j = minj,maxj
          do i = mini,maxi
            if (Prtype(i,j,k)<=0) then
              buffer(i,j,k) = real(ScalarVerticalFlux(i,j,k,Temperature,W), real32)
            else
              buffer(i,j,k) = 0
            end if
          end do
         end do
        end do

        write (unit) BigEnd(buffer)

        write (unit) lf
      end if

      if (store%frame_scalfl==1) then
        if (store%frame_scalars==1) then
          scalname(1:6)="scalfl"
          do l = 1,num_of_scalars
            write(scalname(7:8),"(I2.2)") l
            write (unit) "SCALARS ", scalname , " float", lf
            write (unit) "LOOKUP_TABLE default", lf
            do k = mink,maxk
             do j = minj,maxj
              do i = mini,maxi
                if (Prtype(i,j,k)<=0) then
                  buffer(i,j,k) = real( ScalarVerticalFlux(i,j,k,Scalar(:,:,:,l),W) , real32)
                else
                  buffer(i,j,k) = 0
                end if
              end do
             end do
            end do

            write (unit) BigEnd(buffer)

            write (unit) lf
          end do
        elseif (store%frame_sumscalars==1.and.num_of_scalars>0) then
          write (unit) "SCALARS scalfl float", lf
          write (unit) "LOOKUP_TABLE default", lf
          do k = mink,maxk
           do j = minj,maxj
            do i = mini,maxi
              if (Prtype(i,j,k)<=0) then
                buffer(i,j,k) = real(ScalarVerticalFlux(i,j,k,SUM(Scalar(:,:,:,:),4),W) , real32)
              else
                buffer(i,j,k) = 0
              end if
            end do
           end do
          end do

          write (unit) BigEnd(buffer)

          write (unit) lf
        end if
      end if

      if (store%frame_U==1) then
        write (unit) "VECTORS u float",lf
        do k = mink,maxk
         do j = minj,maxj
          do i = mini,maxi
            if (Prtype(i,j,k)<=0) then
              vbuffer(:,i,j,k) = real([ (U(i,j,k)+U(i-1,j,k))/2, &
                                        (V(i,j,k)+V(i,j-1,k))/2, &
                                        (W(i,j,k)+W(i,j,k-1))/2 ], real32)
            else
              vbuffer(:,i,j,k) = 0
            end if
          end do
         end do
        end do

        write (unit) BigEnd(vbuffer)

        write (unit) lf
      end if

      if (store%frame_vort==1) then
        write (unit) "VECTORS vort float",lf
        do k = mink,maxk
         do j = minj,maxj
          do i = mini,maxi
            if (Prtype(i,j,k)<=0) then
              vbuffer(:,i,j,k) = real([ (W(i,j+1,k)-W(i,j-1,k)+W(i,j+1,k-1)-W(i,j-1,k-1))/(4*dxmin) &
                                       - (V(i,j,k+1)-V(i,j,k-1)+V(i,j-1,k+1)-V(i,j-1,k-1))/(4*dymin), &
                                        (U(i,j,k+1)-U(i,j,k-1)+U(i-1,j,k+1)-U(i-1,j,k-1))/(4*dxmin) &
                                       - (W(i+1,j,k)-W(i-1,j,k)+W(i+1,j,k-1)-W(i-1,j,k-1))/(4*dymin), &
                                        (V(i+1,j,k)-V(i-1,j,k)+V(i+1,j-1,k)-V(i-1,j-1,k))/(4*dxmin) &
                                       - (U(i,j+1,k)-U(i,j-1,k)+U(i-1,j+1,k)-U(i-1,j-1,k))/(4*dymin) ], real32)
            else
              vbuffer(:,i,j,k) = 0
            end if
          end do
         end do
        end do

        write (unit) BigEnd(vbuffer)

      end if
      close(unit)
      if (allocated(buffer)) deallocate(buffer)
      if (allocated(vbuffer)) deallocate(vbuffer)
    end do   !frame_domains
    !$omp end parallel do
  end subroutine Frame



  subroutine BLProfiles(U,V,W,Temperature,Moisture,Scalar)
    real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: Temperature
    real(knd),dimension(-1:,-1:,-1:),intent(in) :: Moisture
    real(knd),dimension(-1:,-1:,-1:,1:),intent(in) :: Scalar
    real(knd) :: S
    real(knd),allocatable,save ::fp(:),ht(:),gp(:)
    integer   :: i,j,k,l,n
    integer,save :: called = 0

    !$omp parallel private(i,j,k,n,S)
    !$omp do
    do k = 0,Unz+1
      S = 0
      n = 0
      do j = 1,Uny
       do i = 1,Unx
         if (Utype(i,j,k)<=0) then
           S = S + U(i,j,k)
           n = n + 1
         end if
       end do
      end do
      profU(k) = S / max(n,1)
    end do
    !$omp end do nowait

    !$omp do
    do k = 1,Vnz+1
      S = 0
      n = 0
      do j = 1,Vny
       do i = 1,Vnx
         if (Vtype(i,j,k)<=0) then
           S = S + V(i,j,k)
           n = n + 1
         end if
       end do
      end do
      profV(k) = S / max(n,1)
    end do
    !$omp end do nowait
    !$omp end parallel

    if (size(Temperature)>0) then
      !$omp parallel do private(i,j,k,n,S)
      do k = 1,Prnz
        S = 0
        n = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k)<=0) then
             S = S + Temperature(i,j,k)
             n = n + 1
           end if
         end do
        end do
        proftemp(k) = S / max(n,1)
      end do
      !$omp end parallel do
    end if

    if (size(Moisture)>0) then
      !$omp parallel do private(i,j,k,n,S)
      do k = 1,Prnz
        S = 0
        n = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k)<=0) then
             S = S + Moisture(i,j,k)
             n = n + 1
           end if
         end do
        end do
        profmoist(k) = S / max(n,1)
      end do
      !$omp end parallel do
    end if

    if (size(Scalar)>0) then
      do l = 1,num_of_scalars
        !$omp parallel do private(i,j,k,n,S)
        do k = 1,Prnz
          S = 0
          n = 0
          do j = 1,Prny
           do i = 1,Prnx
             if (Prtype(i,j,k)<=0) then
               S = S + Scalar(i,j,k,l)
               n = n + 1
             end if
           end do
          end do
          profscal(l,k) = S / max(n,1)
        end do
        !$omp end parallel do
      end do
    end if

    if (called==0) then
      allocate(fp(0:Prnx+1),ht(0:Prnz+1),gp(0:Prny+1))
      !$omp parallel workshare
      forall (i = 0:Prnx+1)      fp(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
      forall (k = 0:Prnz+1)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
      forall (j = 0:Prny+1)      gp(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
      !$omp end parallel workshare
    end if

    !$omp parallel private(i,j,k,n,S)
    !$omp do
    do k = 0,Prnz
      S = 0
      n = 0
      do j = 1,Uny
       do i = 1,Unx
         if ((Utype(i,j,k+1)<=0.or.Utype(i,j,k)<=0).and.(Wtype(i+1,j,k)<=0.or.Wtype(i,j,k)<=0)) then
           S = S + ((ht(k)*U(i,j,k+1)+(1-ht(k))*U(i,j,k))-((1-ht(k))*profU(k)+ht(k)*profU(k+1))) * &
                   (fp(i)*W(i+1,j,k)+(1-fp(i))*W(i,j,k))
           n = n + 1
         end if
       end do
      end do
      profuw(k) = S / max(n,1)
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
      S = 0
      n = 0
      do j = 1,Vny
       do i = 1,Vnx
         if ((Vtype(i,j,k+1)<=0.or.Vtype(i,j,k)<=0).and.(Wtype(i,j+1,k)<=0.or.Wtype(i,j,k)<=0)) then
           S = S + (ht(k)*V(i,j,k+1)+(1-ht(k))*V(i,j,k)-(ht(k)*profV(k+1)+(1-ht(k))*profV(k))) * &
                   (gp(j)*W(i,j+1,k)+(1-gp(j))*W(i,j,k))
           n = n + 1
         end if
       end do
      end do
      profvw(k) = S / max(n,1)
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
      S = 0
      n = 0
      do j = 1,Uny
       do i = 1,Unx
         if (Utype(i,j,k+1)<=0.or.Utype(i,j,k)<=0) then
           S = S-0.25_knd*(Visc(i+1,j,k+1)+Visc(i+1,j,k)+Visc(i,j,k+1)+Visc(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)
           n = n + 1
         end if
       end do
      end do
      profuwsgs(k) = S / max(n,1)
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Prnz
      S = 0
      n = 0
      do j = 1,Vny
       do i = 1,Vnx
         if (Vtype(i,j,k+1)<=0.or.Vtype(i,j,k)<=0) then
           S = S-0.25_knd*(Visc(i,j+1,k+1)+Visc(i,j+1,k)+Visc(i,j,k+1)+Visc(i,j,k))*(V(i,j,k+1)-V(i,j,k))/dzW(k)
           n = n + 1
         end if
       end do
      end do
      profvwsgs(k) = S / max(n,1)
    end do
    !$omp end do nowait

    !$omp do
    do k = 1,Unz
      S = 0
      n = 0
      do j = 1,Uny
       do i = 1,Unx
         if (Utype(i,j,k)<=0) then
           S = S + (U(i,j,k)-profU(k))**2
           n = n + 1
         end if
       end do
      end do
      profuu(k) = S / max(n,1)
    end do
    !$omp end do nowait

    !$omp do
    do k = 1,Vnz
      S = 0
      n = 0
      do j = 1,Vny
       do i = 1,Vnx
         if (Vtype(i,j,k)<=0) then
          S = S + (V(i,j,k)-profV(k))**2
          n = n + 1
         end if
       end do
      end do
      profvv(k) = S / max(n,1)
    end do
    !$omp end do nowait

    !$omp do
    do k = 0,Wnz
      S = 0
      n = 0
      do j = 1,Wny
       do i = 1,Wnx
         if (Wtype(i,j,k)<=0) then
           S = S + (W(i,j,k))**2
           n = n + 1
         end if
       end do
      end do
      profww(k) = S / max(n,1)
    end do
    !$omp end do
    !$omp end parallel

    if (size(Temperature)>0) then
      !proftempfl is computed directly during advection step

      !$omp parallel private(i,j,k,n,S)
      !$omp do
      do k = 1,Prnz
        S = 0
        n = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k)<=0) then
             S = S + (Temperature(i,j,k)-profTemp(k))**2
             n = n + 1
           end if
         end do
        end do
        proftt(k) = S / max(n,1)
      end do
      !$omp end do nowait

      !$omp do
      do k = 0,Prnz
        S = 0
        n = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k+1)<=0.or.Prtype(i,j,k)<=0) then
             S = S-(0.5_knd*(TDiff(i,j,k+1)+TDiff(i,j,k))*(Temperature(i,j,k+1)-Temperature(i,j,k)))/dzW(k)
             n = n + 1
           end if
         end do
        end do
        proftempflsgs(k) = S / max(n,1)
      end do
      !$omp end do
      !$omp end parallel
    end if ! size(Temperature)


    if (size(Moisture)>0) then
      !proftempfl is computed directly during advection step

      !$omp parallel private(i,j,k,n,S)
      !$omp do
      do k = 1,Prnz
        S = 0
        n = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k)<=0) then
             S = S + (Moisture(i,j,k)-profMoist(k))**2
             n = n + 1
           end if
         end do
        end do
        profmm(k) = S / max(n,1)
      end do
      !$omp end do nowait

      !$omp do
      do k = 0,Prnz
        S = 0
        n = 0
        do j = 1,Prny
         do i = 1,Prnx
           if (Prtype(i,j,k+1)<=0.or.Prtype(i,j,k)<=0) then
             S = S-(0.5_knd*(TDiff(i,j,k+1)+TDiff(i,j,k))*(Moisture(i,j,k+1)-Moisture(i,j,k)))/dzW(k)
             n = n + 1
           end if
         end do
        end do
        profmoistflsgs(k) = S / max(n,1)
      end do
      !$omp end do
      !$omp end parallel
    end if ! size(Moisture)


    if (size(Scalar)>0) then
      do l = 1,num_of_scalars
        !$omp parallel private(i,j,k,n,S)
        !$omp do
        do k = 0,Prnz
          S = 0
          n = 0
          do j = 1,Prny
           do i = 1,Prnx
             if (Prtype(i,j,k+1)<=0.or.Prtype(i,j,k)<=0) then
              S = S + 0.5_knd*(Scalar(i,j,k+1,l)+Scalar(i,j,k,l))*(W(i,j,k))
              n = n + 1
             end if
           end do
          end do
          profscalfl(l,k) = S / max(n,1)
        end do
        !$omp end do nowait

        !$omp do
        do k = 1,Prnz
          S = 0
          n = 0
          do j = 1,Prny
           do i = 1,Prnx
             if (Prtype(i,j,k)<=0) then
              S = S + (Scalar(i,j,k,l)-profscal(l,k))**2
              n = n + 1
             end if
           end do
          end do
          profss(l,k) = S / max(n,1)
        end do
        !$omp end do nowait

        !$omp do
        do k = 0,Prnz
          S = 0
          n = 0
          do j = 1,Prny
           do i = 1,Prnx
             if (Prtype(i,j,k+1)<=0.or.Prtype(i,j,k)<=0) then
               S = S-(0.5_knd*(TDiff(i,j,k+1)+TDiff(i,j,k))*(Scalar(i,j,k+1,l)-Scalar(i,j,k,l)))/dzW(k)
               n = n + 1
             end if
           end do
          end do
          profscalflsgs(l,k) = S / max(n,1)
        end do
        !$omp end do
        !$omp end parallel
      end do
    end if ! size(Scalar)

    called = 1
  end subroutine BLProfiles


  function TotKE(U,V,W) result(res)
    real(knd) :: res
    real(knd),dimension(-2:,-2:,-2:) :: U,V,W
    real(knd) Um,Vm,Wm
    integer i,j,k
    real(knd) lx,ly,lz

    lx = xU(Prnx) - xU(0)
    ly = yV(Prny) - yV(0)
    lz = zW(Prnz) - zW(0)

    res = 0
    Um = 0
    Vm = 0
    Wm = 0
    !$omp parallel do private(i,j,k) reduction(+:res)
    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       res = res + (((U(i-1,j,k)+U(i,j,k))/2._knd-Um)**2+&
                       ((V(i,j-1,k)+V(i,j,k))/2._knd-Vm)**2+&
                       ((W(i,j,k-1)+W(i,j,k))/2._knd-Wm)**2)
      end do
     end do
    end do
    !$omp end parallel do
    res = res*lx*ly*lz/2
  end function TotKE

  pure real(knd) function VorticityMag(i,j,k,U,V,W) result(res)
    integer,intent(in) :: i,j,k
    real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W

    res = sum(Vorticity(i,j,k,U,V,W)**2)
    res = Sqrt(res)
  end function VorticityMag


  pure real(knd) function Lambda2(i,j,k,U,V,W) result(res)
    integer,intent(in) :: i,j,k
    real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W

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



  subroutine OutputU2(U,V,W)
    real(knd),dimension(-2:,-2:,-2:) :: U,V,W
    integer unit
    character(70) :: str

    call newunit(unit)

    open(unit,file="output/U2.vtk")
    write (unit) "# vtk DataFile Version 2.0",lf
    write (unit) "CLMM output file",lf
    write (unit) "BINARY",lf
    write (unit) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write (str(12:),*) Unx,Uny,Unz
    write (unit) str,lf
    str="X_COORDINATES"
    write (str(15:),'(i5,2x,a)') Unx,"float"
    write (unit) str,lf
    write (unit) BigEnd(real(xU(1:Unx), real32)),lf
    str="Y_COORDINATES"
    write (str(15:),'(i5,2x,a)') Uny,"float"
    write (unit) str,lf
    write (unit) BigEnd(real(yPr(1:Uny), real32)),lf
    str="Z_COORDINATES"
    write (str(15:),'(i5,2x,a)') Unz,"float"
    write (unit) str,lf
    write (unit) BigEnd(real(zPr(1:Unz), real32)),lf
    str="POINT_DATA"
    write (str(12:),*) Unx*Uny*Unz
    write (unit) str,lf


    write (unit) "SCALARS U float",lf
    write (unit) "LOOKUP_TABLE default",lf

    write (unit) BigEnd(real(U(1:Unx,1:Uny,1:Unz), real32)),lf

    write (unit) lf
    close(unit)


    open(unit,file="output/V2.vtk")
    write (unit) "# vtk DataFile Version 2.0",lf
    write (unit) "CLMM output file",lf
    write (unit) "BINARY",lf
    write (unit) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write (str(12:),*) Vnx,Vny,Vnz
    write (unit) str,lf
    str="X_COORDINATES"
    write (str(15:),'(i5,2x,a)') Vnx,"float"
    write (unit) str,lf
    write (unit) BigEnd(real(xPr(1:Vnx), real32)),lf
    str="Y_COORDINATES"
    write (str(15:),'(i5,2x,a)') Vny,"float"
    write (unit) str,lf
    write (unit) BigEnd(real(yV(1:Vny), real32)),lf
    str="Z_COORDINATES"
    write (str(15:),'(i5,2x,a)') Vnz,"float"
    write (unit) str,lf
    write (unit) BigEnd(real(zPr(1:Vnz), real32)),lf
    str="POINT_DATA"
    write (str(12:),*) Vnx*Vny*Vnz
    write (unit) str,lf


    write (unit) "SCALARS V float",lf
    write (unit) "LOOKUP_TABLE default",lf

    write (unit) BigEnd(real(V(1:Vnx,1:Vny,1:Vnz), real32)),lf

    write (unit) lf
    close(unit)


    open(unit,file="output/W2.vtk")
    write (unit) "# vtk DataFile Version 2.0",lf
    write (unit) "CLMM output file",lf
    write (unit) "BINARY",lf
    write (unit) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write (str(12:),*) Wnx,Wny,Wnz
    write (unit) str,lf
    str="X_COORDINATES"
    write (str(15:),'(i5,2x,a)') Wnx,"float"
    write (unit) str,lf
    write (unit) BigEnd(real(xPr(1:Wnx), real32)),lf
    str="Y_COORDINATES"
    write (str(15:),'(i5,2x,a)') Wny,"float"
    write (unit) str,lf
    write (unit) BigEnd(real(yPr(1:Wny), real32)),lf
    str="Z_COORDINATES"
    write (str(15:),'(i5,2x,a)') Wnz,"float"
    write (unit) str,lf
    write (unit) BigEnd(real(zW(1:Wnz), real32)),lf
    str="POINT_DATA"
    write (str(12:),*) Wnx*Wny*Wnz
    write (unit) str,lf


    write (unit) "SCALARS W float",lf
    write (unit) "LOOKUP_TABLE default",lf

    write (unit) BigEnd(real(W(1:Wnx,1:Wny,1:Wnz), real32)),lf

    write (unit) lf
    close(unit)
  end subroutine OutputU2


  subroutine OUTINLET(U,V,W,Temperature)
    !for output of 2d data for use as an inilet condition later
    real(knd),dimension(-2:,-2:,-2:) :: U,V,W
    real(knd),dimension(-1:,-1:,-1:) :: Temperature
    integer,save ::fnum
    integer,save :: called = 0

    if ((time>=timefram1).and.(time<=timefram2+(timefram2-timefram1)/(frames-1))&
        .and.(time>=timefram1+fnum*(timefram2-timefram1)/(frames-1))) then
     if (called==0) then
      open(101,file="output/inletframeinfo.unf",form='unformatted',status='replace',action='write')
      write(101) Prny,Prnz  !for check of consistency of grids before use
      write(101) Vny
      write(101) Wnz
      write(101) dxPr(0)
      called = 1
      fnum = 0
     end if
     fnum = fnum+1
     write(101) time-timefram1
     call OUTINLETFrame(U,V,W,Temperature,fnum)
    elseif (time>timefram2+(timefram2-timefram1)/(frames-1).and.called==1) then
      close(101)
      called = 2
    end if

  end subroutine OUTINLET


  subroutine OUTINLETFrame(U,V,W,Temperature,n)
    real(knd),intent(in) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(knd),dimension(-1:,-1:,-1:),intent(in)   :: Temperature
    integer n
    character(12) :: fname
    integer mini,maxi,minj,maxj,mink,maxk,unit

    call Gridcoords(mini,minj,mink,store%frame_domains(1)%position,(yV(Prny+1)+yV(0))/2._knd,(zW(Prnz+1)+zW(0))/2._knd)
    maxi = mini
    minj = 1
    maxj = Prny
    mink = 1
    maxk = Prnz


    fname(1:5)="frame"
    write(fname(6:8),"(I3.3)") n
    fname(9:12)=".unf"
    write(*,*) "Saving frame:",fname(1:6),"   time:",time

    call newunit(unit)

    open(unit,file = fname,form='unformatted',access='sequential',status='replace',action='write')


    write(unit) U(mini,1:Uny,1:Unz)
    write(unit) V(mini,1:Vny,1:Vnz)
    write(unit) W(mini,1:Wny,1:Wnz)
    if (enable_buoyancy==1) then
         write(unit) Temperature(mini,1:Prny,1:Prnz)
    end if
    close(unit)

  end subroutine OUTINLETFrame



  pure real(knd) function TriLinInt(a,b,c,vel000,vel100,vel010,vel001,vel110,vel101,vel011,vel111)
  real(knd),intent(in) :: a,b,c,vel000,vel100,vel010,vel001,vel110,vel101,vel011,vel111

    TriLinInt=   (1-a)*(1-b)*(1-c)*vel000+&
                 a*(1-b)*(1-c)*vel100+&
                 (1-a)*b*(1-c)*vel010+&
                 (1-a)*(1-b)*c*vel001+&
                 a*b*(1-c)*vel110+&
                 a*(1-b)*c*vel101+&
                 (1-a)*b*c*vel011+&
                 a*b*c*vel111

  end function TriLinInt


end module Outputs

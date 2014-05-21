module Outputs
  use Parameters
  use Boundaries
  use Scalars
  use Wallmodels, only: GroundDeposition, GroundUstar, wallmodeltype
  use ImmersedBoundary
  use Turbinlet, only: Ustar_inlet
  use Endianness
  use FreeUnit
  use VTKFrames
  use StaggeredFrames
  use Output_helpers

  implicit none
 

  private
  public store, display, probes, scalar_probes,  &
         OutTStep, Output, AllocateOutputs, ReadProbes,  &
          proftempfl, profmoistfl, profuw, profvw, profuwsgs, profvwsgs

  real(knd),dimension(:),allocatable :: profuavg,profuavg2,profvavg,profvavg2,profuuavg,profvvavg,profwwavg, &
                                         profU,profV,profuu,profvv,profww,proftauavg,proftau,proftausgs,proftausgsavg, &
                                         proftemp,proftempfl,proftempavg,proftempavg2,proftempflavg, &
                                         proftempflsgs,proftempflsgsavg,proftt,profttavg, &
                                         profmoist,profmoistfl,profmoistavg,profmoistavg2,profmoistflavg, &
                                         profmoistflsgs,profmoistflsgsavg,profmm,profmmavg, &
                                         profuw,profuwavg,profuwsgs,profuwsgsavg, &
                                         profvw,profvwavg,profvwsgs,profvwsgsavg

  real(knd),dimension(:,:),allocatable ::profscal,profscalfl,profscalavg,profscalavg2,profscalflavg, &  !which scalar, height
                                         profscalflsgs,profscalflsgsavg,profss,profssavg

  real(knd),allocatable :: U_avg(:,:,:),V_avg(:,:,:),W_avg(:,:,:) !<u>
  real(knd),allocatable :: U_rms(:,:,:),V_rms(:,:,:),W_rms(:,:,:) !<uu>, <u'u'> must be computed before saving
  real(knd),allocatable :: Pr_avg(:,:,:) !<p>
  real(knd),allocatable :: Temperature_avg(:,:,:) !<theta>
  real(knd),allocatable :: Moisture_avg(:,:,:) !<q>
  real(knd),allocatable :: Scalar_avg(:,:,:,:) !<c>

  real(knd),allocatable :: Scalar_max(:,:,:,:)

  real(knd),allocatable :: Scalar_intermitency(:,:,:,:)

  real(knd),allocatable :: Scalar_fl_U_avg(:,:,:,:) !<cu>, <c'u'> must be computed before saving
  real(knd),allocatable :: Scalar_fl_V_avg(:,:,:,:)
  real(knd),allocatable :: Scalar_fl_W_avg(:,:,:,:)

  real(TIM),allocatable,dimension(:) :: times                                !times of the timesteps

  real(knd),allocatable,dimension(:) :: delta_time,tke,dissip

  real(knd),allocatable,dimension(:,:) :: ustar,tstar                        !first index differentiates flux from friction number
                                                                             !second index is time

  real(knd),allocatable,dimension(:,:) :: U_time,V_time,W_time,Pr_time,temp_time,moist_time  !position, time

  real(knd),allocatable,dimension(:,:,:) :: scalp_time                        !which scalar, position, time
  real(knd),allocatable,dimension(:,:) :: scalsum_time                        !which scalar, time

  real(knd),allocatable,dimension(:,:,:) :: momentum_fluxes_time, momentum_fluxes_sgs_time   !component, position, time
                                                       !components: 1,1; 1,2; 1,3; 2,2; 2,3; 3,3 for 1..6

  type TProbe
    integer :: Ui,Uj,Uk,Vi,Vj,Vk,Wi,Wj,Wk    !grid coordinates of probes in the U,V,W grids
    integer :: i,j,k                         !grid coordinates of probes in the scalar grid
    real(knd) :: x,y,z                       !physical coordinates of probes
    logical :: inside
    integer :: number
  end type TProbe

  !for flow variables including scalar ones (temperature, moisture)
  type(TProbe),allocatable,dimension(:),save :: probes

  !for fluxes
  type(TProbe),allocatable,dimension(:),save :: flux_probes

  !for passive scalars
  type(TProbe),allocatable,dimension(:),save :: scalar_probes

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
    integer :: scalars_avg = 1
    integer :: scalars_max = 0
    integer :: scalars_intermitency = 0
    
    real(knd) :: scalars_intermitency_threshold = epsilon(1._knd)

    integer :: deposition = 0


    integer :: out_U = 1
    integer :: out_vorticity = 0
    integer :: out_Pr = 0
    integer :: out_Prtype = 0
    integer :: out_lambda2 = 0
    integer :: out_temperature = 1
    integer :: out_moisture = 1
    integer :: out_divergence = 0
    integer :: out_viscosity = 0

    integer :: avg_U = 1  !1..only in avg.vtk, 2..only in separate Xavg.vtk, 3..both
    integer :: avg_vorticity = 0
    integer :: avg_Pr = 1
    integer :: avg_Prtype = 1
    integer :: avg_temperature = 1
    integer :: avg_moisture = 1

    integer :: avg_U_rms = 0 !1..only in  avg.vtk, 2..only in separate Xavg.vtk, 3..both
    integer :: avg_flux_scalar = 0

    integer :: delta_time = 0
    integer :: tke = 0
    integer :: dissip = 0
    integer :: scalsum_time = 1
    integer :: scaltotsum_time = 0
    integer :: ustar = 1
    integer :: tstar = 1

    integer :: blprofiles = 0
    
    integer :: probes_fluxes = 0
  end type TOutputSwitches

  type(TOutputSwitches),save :: store

  type TDisplaySwitches
    integer :: delta = 0
    integer :: ustar = 0
    integer :: tstar = 0
  end type

  type(TDisplaySwitches),save :: display

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
      call error_stop("Error: File "//trim(pfile)//" could not be opened!")
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


   if (store%avg_U==0.and.store%avg_U_rms>1) store%avg_U = 1

   if (averaging==1) then
     if (store%avg_U>0) then
       allocate(U_avg(-2:Unx+3,-2:Uny+3,-2:Unz+3))
       allocate(V_avg(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
       allocate(W_avg(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
       U_avg = 0
       V_avg = 0
       W_avg = 0
     end if

     if (store%avg_U_rms>0) then
       allocate(U_rms(-2:Unx+3,-2:Uny+3,-2:Unz+3))
       allocate(V_rms(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
       allocate(W_rms(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
       U_rms = 0
       V_rms = 0
       W_rms = 0
     end if

     if (store%avg_Pr==1) then
       allocate(Pr_avg(1:Prnx,1:Prny,1:Prnz))
       Pr_avg = 0
     end if

     if (enable_buoyancy.and.store%avg_temperature==1) then
       allocate(Temperature_avg(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
       Temperature_avg = 0
     end if

     if (enable_moisture.and.store%avg_moisture==1) then
       allocate(Moisture_avg(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
       Moisture_avg = 0
     end if
   end if

   if (num_of_scalars>0.and.store%scalars_avg==1) then
     allocate(Scalar_avg(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,num_of_scalars))
     Scalar_avg = 0
   else
     allocate(Scalar_avg(0,0,0,0))
   end if
   
   if (num_of_scalars>0.and.store%scalars_max==1) then
     allocate(Scalar_max(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,num_of_scalars))
     Scalar_max = 0
   else
     allocate(Scalar_max(0,0,0,0))
   end if

   if (num_of_scalars>0.and.store%scalars_intermitency==1) then
     allocate(Scalar_intermitency(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,num_of_scalars))
     Scalar_intermitency = 0
   else
     allocate(Scalar_intermitency(0,0,0,0))
   end if

   if (num_of_scalars>0.and.store%avg_flux_scalar==1) then
     allocate(Scalar_fl_U_avg(Unx,Uny,Unz,num_of_scalars))
     allocate(Scalar_fl_V_avg(Vnx,Vny,Vnz,num_of_scalars))
     allocate(Scalar_fl_W_avg(Wnx,Wny,Wnz,num_of_scalars))
     Scalar_fl_U_avg = 0
     Scalar_fl_V_avg = 0
     Scalar_fl_W_avg = 0
   end if


   allocate(times(0:max_number_of_time_steps))

   
   do k = 1,size(probes)
     associate(p => probes(k))
       p%number = k
       
       p%inside = InDomain(p%x,p%y,p%z)
       
       call GridCoords(p%i,p%j,p%k,p%x,p%y,p%z)

       p%i = max(p%i,1)
       p%j = max(p%j,1)
       p%k = max(p%k,1)
       p%i = min(p%i,Prnx)
       p%j = min(p%j,Prny)
       p%k = min(p%k,Prnz)

       call GridCoords_U(p%Ui,p%Uj,p%Uk,p%x,p%y,p%z)
       call GridCoords_V(p%Vi,p%Vj,p%Vk,p%x,p%y,p%z)
       call GridCoords_W(p%Wi,p%Wj,p%Wk,p%x,p%y,p%z)
     end associate
   end do
   
   probes = pack(probes, probes%inside)

   do k = 1,size(scalar_probes)
     associate(p => scalar_probes(k))
       p%number = k
       
       p%inside = InDomain(p%x,p%y,p%z)
       
       call GridCoords(p%i,p%j,p%k,p%x,p%y,p%z)

       p%i = max(p%i,1)
       p%j = max(p%j,1)
       p%k = max(p%k,1)
       p%i = min(p%i,Prnx)
       p%j = min(p%j,Prny)
       p%k = min(p%k,Prnz)
     end associate
   end do

   scalar_probes = pack(scalar_probes, scalar_probes%inside)

   if (size(probes)>0) then

     allocate(U_time(size(probes),0:max_number_of_time_steps), &
              V_time(size(probes),0:max_number_of_time_steps), &
              W_time(size(probes),0:max_number_of_time_steps), &
              Pr_time(size(probes),0:max_number_of_time_steps))
     times = huge(1.0_knd)
     U_time = huge(1.0_knd)
     V_time = huge(1.0_knd)
     W_time = huge(1.0_knd)
     Pr_time = huge(1.0_knd)

     if (store%probes_fluxes==1) then
       allocate(momentum_fluxes_time(6,size(probes),0:max_number_of_time_steps))
       momentum_fluxes_time = huge(1.0)
       allocate(momentum_fluxes_sgs_time(6,size(probes),0:max_number_of_time_steps))
       momentum_fluxes_sgs_time = huge(1.0)
     end if

     if (enable_buoyancy) then
       allocate(temp_time(size(probes),0:max_number_of_time_steps))
       temp_time = huge(1.0_knd)
     end if

     if (enable_moisture) then
       allocate(moist_time(size(probes),0:max_number_of_time_steps))
       moist_time = huge(1.0_knd)
     end if

   end if


   if (store%delta_time==1) then
     allocate(delta_time(0:max_number_of_time_steps))
     delta_time = huge(1.0_knd)
   end if

   if (store%tke==1) then
     allocate(tke(0:max_number_of_time_steps))
     tke = huge(1.0_knd)
   end if

   if (store%delta_time==1) then
     allocate(dissip(0:max_number_of_time_steps))
     dissip = huge(1.0_knd)
     dissip(0)=0
   end if

   if (wallmodeltype>0.and.(display%ustar==1.or.store%ustar==1)) then
     allocate(ustar(2,0:max_number_of_time_steps))
     ustar = huge(1.0)
   end if

   if (wallmodeltype>0.and.enable_buoyancy.and.TempBtype(Bo)==DIRICHLET.and.(display%tstar==1.or.store%tstar==1)) then
     allocate(tstar(2,0:max_number_of_time_steps))
     tstar = huge(1.0)
   end if

   if (num_of_scalars>0) then
     if (store%scalsum_time==1.or.store%scaltotsum_time==1) then
       allocate(scalsum_time(1:num_of_scalars,0:max_number_of_time_steps))
       scalsum_time = huge(1.0_knd)
     end if

     if (size(scalar_probes)>0) then
       allocate(scalp_time(1:num_of_scalars,1:size(scalar_probes),0:max_number_of_time_steps))
       scalp_time = huge(1.0_knd)
    end if
   end if

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
   
     if (Btype(To)==AUTOMATICFLUX) then
       allocate(profuw(0:Prnz),profuwsgs(0:Prnz))
       allocate(profvw(0:Prnz),profvwsgs(0:Prnz))
     end if
     
     if (TempBtype(To)==AUTOMATICFLUX) then
       allocate(proftempfl(0:Prnz),proftempflsgs(0:Prnz))
     else
       !to avoid the necessity of an allocatable dummy argument
       allocate(proftempfl(0))
     end if

     if (MoistBtype(To)==AUTOMATICFLUX) then
       allocate(profmoistfl(0:Prnz),profmoistflsgs(0:Prnz))
     else
       !to avoid the necessity of an allocatable dummy argument
       allocate(profmoistfl(0))
     end if

   end if


   call GetEndianness

#if _WIN32 || _WIN64
   call execute_command_line("mkdir "//trim(output_dir))
#else
   call execute_command_line("mkdir -p "//trim(output_dir))
#endif

   contains 
     !FIXME: delete this as soon as supported by ifort
     subroutine execute_command_line(cmd)
       character(*) :: cmd
       call system(cmd)
     end subroutine
  end subroutine AllocateOutputs













  subroutine OutTStep(U,V,W,Pr,Temperature,Moisture,Scalar,dt,delta)
    use Wallmodels, only: ComputeViscsWM
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in)   :: U,V,W
    real(knd),dimension(1:,1:,1:),contiguous,intent(in)      :: Pr
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in)   :: Temperature
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in)   :: Moisture
    real(knd),dimension(-1:,-1:,-1:,:),contiguous,intent(in) :: Scalar
    real(knd),intent(in) :: dt
    real(knd),intent(in) :: delta

    integer :: l,i,j,k
    real(knd) :: S,S2
    real(knd) :: time_weight
    real(knd) :: fl_L, fl_R
    integer, save :: fnum = 0

    ! We compute the subgrid stresses from the eddy viscosity even at the walls
    ! We should not set also TDiff
    if (wallmodeltype>0.and.store%probes_fluxes==1.and.store%BLprofiles==1) then
      call ComputeViscsWM(U,V,W,Pr,Temperature)
    end if

    times(step)=time

    if (store%scalsum_time==1.or.store%scaltotsum_time==1) then
      do l = 1,num_of_scalars
         S = 0
         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
            if (Prtype(i,j,k)<=0) S = S + Scalar(i,j,k,l)*dxPr(i)*dyPr(j)*dzPr(k)
           end do
          end do
         end do
         scalsum_time(l,step) = S
      end do
    end if



    do k = 1,size(probes)
      associate (p=> probes(k))
        U_time(k,step)=Trilinint((p%x-xU(p%Ui))/(xU(p%Ui+1)-xU(p%Ui)), &
                             (p%y-yPr(p%Uj))/(yPr(p%Uj+1)-yPr(p%Uj)), &
                             (p%z-zPr(p%Uk))/(zPr(p%Uk+1)-zPr(p%Uk)), &
                             U(p%Ui,p%Uj,p%Uk),U(p%Ui+1,p%Uj,p%Uk), &
                             U(p%Ui,p%Uj+1,p%Uk),U(p%Ui,p%Uj,p%Uk+1), &
                             U(p%Ui+1,p%Uj+1,p%Uk),U(p%Ui+1,p%Uj,p%Uk+1), &
                             U(p%Ui,p%Uj+1,p%Uk+1),U(p%Ui+1,p%Uj+1,p%Uk+1))

        V_time(k,step)=Trilinint((p%x-xPr(p%Vi))/(xPr(p%Vi+1)-xPr(p%Vi)), &
                             (p%y-yV(p%Vj))/(yV(p%Vj+1)-yV(p%Vj)), &
                             (p%z-zPr(p%Vk))/(zPr(p%Vk+1)-zPr(p%Vk)), &
                             V(p%Vi,p%Vj,p%Vk),V(p%Vi+1,p%Vj,p%Vk), &
                             V(p%Vi,p%Vj+1,p%Vk),V(p%Vi,p%Vj,p%Vk+1), &
                             V(p%Vi+1,p%Vj+1,p%Vk),V(p%Vi+1,p%Vj,p%Vk+1), &
                             V(p%Vi,p%Vj+1,p%Vk+1),V(p%Vi+1,p%Vj+1,p%Vk+1))

        W_time(k,step)=Trilinint((p%x-xPr(p%Wi))/(xPr(p%Wi+1)-xPr(p%Wi)), &
                             (p%y-yPr(p%Wj))/(yPr(p%Wj+1)-yPr(p%Wj)), &
                             (p%z-zW(p%Wk))/(zW(p%Wk+1)-zW(p%Wk)), &
                             W(p%Wi,p%Wj,p%Wk),W(p%Wi+1,p%Wj,p%Wk), &
                             W(p%Wi,p%Wj+1,p%Wk),W(p%Wi,p%Wj,p%Wk+1), &
                             W(p%Wi+1,p%Wj+1,p%Wk),W(p%Wi+1,p%Wj,p%Wk+1), &
                             W(p%Wi,p%Wj+1,p%Wk+1),W(p%Wi+1,p%Wj+1,p%Wk+1))

        Pr_time(k,step)=Trilinint((p%x-xPr(p%i))/(xPr(p%i+1)-xPr(p%i)), &
                             (p%y-yPr(p%j))/(yPr(p%j+1)-yPr(p%j)), &
                             (p%z-zPr(p%k))/(zPr(p%k+1)-zPr(p%k)), &
                              Pr(p%i,p%j,p%k), &
                              Pr(min(p%i+1,Unx+1),p%j,p%k), &
                              Pr(p%i,min(p%j+1,Vny+1),p%k), &
                              Pr(p%i,p%j,min(p%k+1,Wnz+1)), &
                              Pr(min(p%i+1,Unx+1),min(p%j+1,Vny+1),p%k), &
                              Pr(min(p%i+1,Unx+1),p%j,min(p%k+1,Wnz+1)), &
                              Pr(p%i,min(p%j+1,Vny+1),min(p%k+1,Wnz+1)), &
                              Pr(min(p%i+1,Unx+1),min(p%j+1,Vny+1),min(p%k+1,Wnz+1)))

        if (enable_buoyancy) then
          temp_time(k,step)=Trilinint((p%x-xPr(p%i))/(xPr(p%i+1)-xPr(p%i)), &
                             (p%y-yPr(p%j))/(yPr(p%j+1)-yPr(p%j)), &
                             (p%z-zPr(p%k))/(zPr(p%k+1)-zPr(p%k)), &
                             Temperature(p%i,p%j,p%k), &
                             Temperature(p%i+1,p%j,p%k), &
                             Temperature(p%i,p%j+1,p%k), &
                             Temperature(p%i,p%j,p%k+1), &
                             Temperature(p%i+1,p%j+1,p%k), &
                             Temperature(p%i+1,p%j,p%k+1), &
                             Temperature(p%i,p%j+1,p%k+1), &
                             Temperature(p%i+1,p%j+1,p%k+1))
        end if

        if (enable_moisture) then
          moist_time(k,step)=Trilinint((p%x-xPr(p%i))/(xPr(p%i+1)-xPr(p%i)), &
                             (p%y-yPr(p%j))/(yPr(p%j+1)-yPr(p%j)), &
                             (p%z-zPr(p%k))/(zPr(p%k+1)-zPr(p%k)), &
                             Moisture(p%i,p%j,p%k), &
                             Moisture(p%i+1,p%j,p%k), &
                             Moisture(p%i,p%j+1,p%k), &
                             Moisture(p%i,p%j,p%k+1), &
                             Moisture(p%i+1,p%j+1,p%k), &
                             Moisture(p%i+1,p%j,p%k+1), &
                             Moisture(p%i,p%j+1,p%k+1), &
                             Moisture(p%i+1,p%j+1,p%k+1))
        end if
        
        if (store%probes_fluxes==1) then 
          !NOTE: valid only for central schemes
          momentum_fluxes_time(1,k,step) = (U_time(k,step))**2
          momentum_fluxes_time(2,k,step) = U_time(k,step) * V_time(k,step)
          momentum_fluxes_time(3,k,step) = U_time(k,step) * W_time(k,step)
          momentum_fluxes_time(4,k,step) = (V_time(k,step))**2
          momentum_fluxes_time(5,k,step) = V_time(k,step) * W_time(k,step)
          momentum_fluxes_time(6,k,step) = (W_time(k,step))**2
          
          !NOTE: subgrid part evaluated in the center of the cell
          momentum_fluxes_sgs_time(1,k,step) = (U(p%i, p%j, p%k) - U(p%i-1, p%j, p%k)) / &
                                             dxPr(p%i)
          momentum_fluxes_sgs_time(4,k,step) = (V(p%i, p%j, p%k) - V(p%i, p%j-1, p%k)) / &
                                             dyPr(p%j)
          momentum_fluxes_sgs_time(6,k,step) = (W(p%i, p%j, p%k) - W(p%i, p%j, p%k-1)) / &
                                             dzPr(p%k)

          fl_L = ( V(p%i+1, p%j, p%k)+V(p%i+1, p%j-1, p%k) - &
                   V(p%i-1, p%j, p%k)+V(p%i-1, p%j-1, p%k) ) / &
                 (2*(xPr(i+1)-xPr(i-1)))
          fl_R = ( U(p%i, p%j+1, p%k)+U(p%i-1, p%j+1, p%k) - &
                   U(p%i, p%j-1, p%k)+U(p%i-1, p%j-1, p%k) ) / &
                 (2*(yPr(j+1)-yPr(j-1)))
          momentum_fluxes_sgs_time(2,k,step) = (fl_L + fl_R) / 2
                 
          fl_L = ( W(p%i+1, p%j, p%k)+W(p%i+1, p%j, p%k-1) - &
                   W(p%i-1, p%j, p%k)+W(p%i-1, p%j, p%k-1) ) / &
                 (2*(xPr(i+1)-xPr(i-1)))
          fl_R = ( U(p%i, p%j, p%k+1)+U(p%i-1, p%j, p%k+1) - &
                   U(p%i, p%j, p%k-1)+U(p%i-1, p%j, p%k-1) ) / &
                 (2*(zPr(k+1)-zPr(k-1)))
          momentum_fluxes_sgs_time(3,k,step) = (fl_L + fl_R) / 2
                                  
          fl_L = ( W(p%i, p%j+1, p%k)+W(p%i, p%j+1, p%k-1) - &
                   W(p%i, p%j-1, p%k)+W(p%i, p%j-1, p%k-1) ) / &
                 (2*(yPr(j+1)-yPr(j-1)))
          fl_R = ( V(p%i, p%j, p%k+1)+V(p%i, p%j-1, p%k+1) - &
                   V(p%i, p%j, p%k-1)+V(p%i, p%j-1, p%k-1) ) / &
                 (2*(zPr(k+1)-zPr(k-1)))
          momentum_fluxes_sgs_time(5,k,step) = (fl_L + fl_R) / 2
                 
          momentum_fluxes_sgs_time(:,k,step) = Viscosity(p%i, p%j, p%k) * &
                                             momentum_fluxes_sgs_time(:,k,step)
        end if
      end associate
    end do

    do k = 1,size(scalar_probes)    
      associate (p=> scalar_probes(k))
        do l = 1,num_of_scalars
          scalp_time(l,k,step)=Trilinint((p%x-xPr(p%i))/(xPr(p%i+1)-xPr(p%i)), &
                             (p%y-yPr(p%j))/(yPr(p%j+1)-yPr(p%j)), &
                             (p%z-zPr(p%k))/(zPr(p%k+1)-zPr(p%k)), &
                             Scalar(p%i,p%j,p%k,l), &
                             Scalar(p%i+1,p%j,p%k,l), &
                             Scalar(p%i,p%j+1,p%k,l), &
                             Scalar(p%i,p%j,p%k+1,l), &
                             Scalar(p%i+1,p%j+1,p%k,l), &
                             Scalar(p%i+1,p%j,p%k+1,l), &
                             Scalar(p%i,p%j+1,p%k+1,l), &
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

    if (store%delta_time==1.and.dt>0) then
      delta_time(step)=delta/dt
    end if

    endstep = step






    if (display%delta==1) then
      if (master) write(*,*) "delta: ",delta
    end if





    if ((averaging==1).and.((time>=timeavg1).and.(time-dt<timeavg2))) then

      time_weight = min(dt, time-timeavg1, timeavg2-(time-dt), timeavg2-timeavg1) / &
                    (timeavg2-timeavg1)

      if (store%avg_U>0) then
         U_avg = U_avg + U * time_weight
         V_avg = V_avg + V * time_weight
         W_avg = W_avg + W * time_weight
      end if

      if (store%avg_U_rms>0) then
        U_rms = U_rms + U**2 * time_weight
        V_rms = V_rms + V**2 * time_weight
        W_rms = W_rms + W**2 * time_weight
        call AddSubgridRMS(U_rms,V_rms,W_rms,U,V,W,time_weight)
      end if

      if (store%avg_Pr==1) then
        Pr_avg = Pr_avg + Pr(1:Prnx,1:Prny,1:Prnz) * time_weight
      end if

      if (enable_buoyancy.and.store%avg_temperature==1) then
        Temperature_avg = Temperature_avg + temperature * time_weight
      end if

      if (enable_moisture.and.store%avg_moisture==1) then
        Moisture_avg = Moisture_avg + moisture * time_weight
      end if

      if (num_of_scalars>0.and.store%scalars_avg==1) then
        Scalar_avg = Scalar_avg + Scalar * time_weight
      end if

      if (num_of_scalars>0.and.store%scalars_max==1) then
        Scalar_max = max(Scalar_max, Scalar)
      end if

      if (num_of_scalars>0.and.store%scalars_intermitency==1) then
        where (Scalar>=store%scalars_intermitency_threshold)  &
          Scalar_intermitency = Scalar_intermitency + time_weight
      end if

      if (num_of_scalars>0.and.store%avg_flux_scalar==1) then
        do i=1,num_of_scalars
          call AddScalarAdvVector(Scalar_fl_U_avg(:,:,:,i), &
                               Scalar_fl_V_avg(:,:,:,i), &
                               Scalar_fl_W_avg(:,:,:,i), &
                               Scalar(:,:,:,i), &
                               U,V,W,time_weight)
          call AddScalarDiffVector(Scalar_fl_U_avg(:,:,:,i), &
                                Scalar_fl_V_avg(:,:,:,i), &
                                Scalar_fl_W_avg(:,:,:,i), &
                                Scalar(:,:,:,i),time_weight)
        end do
      end if
   end if


   if (wallmodeltype>0.and.(display%ustar==1.or.store%ustar==1)) then

     S = GroundUstar()

     S2 = S*Re

     if (display%ustar==1) then
       if (allocated(Ustar_inlet)) then
        if (master) write(*,*) "ustar:",S,"Re_tau:",S2,"u*inlet",Ustar_inlet(1)
       else
        if (master) write(*,*) "ustar:",S,"Re_tau:",S2
       end if
     end if
     if (store%ustar==1) then
       ustar(:,step)=[ S2 , S ]
     end if
   end if


    if (wallmodeltype>0.and.enable_buoyancy.and.TempBtype(Bo)==DIRICHLET.and.(display%tstar==1.or.store%tstar==1)) then
     S2 = SUM(BsideTFLArr(1:Prnx,1:Prny))/(Prnx*Prny)
     S=-S*S2

     if (display%tstar==1) then
       if (master) write(*,*) "Tstar",S,"tflux", S2
     end if

     if (store%tstar==1) then
       tstar(:,step) = [ S2,S ]
     end if
    end if

    if ((store%BLprofiles==1 .and. &
         (averaging==1 .and. &
         ( time>=timeavg1 .and. &
         time<=timeavg2))) .or. &
       Btype(To)==AUTOMATICFLUX) then
      call StressProfiles(U,V,W)
    end if

    if (store%BLprofiles==1 .and. &
         (averaging==1 .and. &
         ( time>=timeavg1 .and. &
         time<=timeavg2))) then

      call FluxSGSProfiles(W,Temperature,Moisture,Scalar)

    else if (TempBtype(To)==AUTOMATICFLUX .or. &
             MoistBtype(To)==AUTOMATICFLUX .or. &
             ScalBType(To)==AUTOMATICFLUX) then

      if (TempBtype(To)==AUTOMATICFLUX.and.enable_buoyancy) &
        call TemperatureFluxSGSProfile(W,Temperature)

      if (MoistBtype(To)==AUTOMATICFLUX.and.enable_moisture) &
        call MoistureFluxSGSProfile(W,Moisture)

      if (ScalBtype(To)==AUTOMATICFLUX.and.num_of_scalars>0) &
        call ScalarFluxSGSProfile(W,Scalar)

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

       if (enable_buoyancy) then
         proftempavg = proftempavg + proftemp * time_weight
         proftempflavg = proftempflavg + proftempfl * time_weight
         proftempflsgsavg = proftempflsgsavg + proftempflsgs * time_weight
         profttavg = profttavg + proftt * time_weight
       end if

       if (enable_buoyancy) then
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


    call SaveVTKFrames(time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)

    call SaveStaggeredFrames(time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
    
    if (mod(endstep,1000) == 0) then
      call OutputTimeSeries
    end if

  end subroutine OutTstep


  subroutine AddSubgridRMS(U_r,V_r,W_r,U,V,W,weight)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout)   :: U_r,V_r,W_r
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in)   :: U,V,W
    real(knd),intent(in) :: weight
    real(knd) :: Ax,Ay,Az
    integer :: i,j,k
    
    Ax = weight/(2*dxmin)
    Ay = weight/(2*dymin)
    Az = weight/(2*dzmin)

    !NOTE:neglecting part caused by molecular viscosity
    !$omp parallel private(i,j,k)
    !$omp do
    do k=1,Unz
      do j=1,Uny
        do i=1,Unx
          U_r(i,j,k) =  U_r(i,j,k) + &
                        Ax * (Viscosity(i+1,j,k)*(U(i+1,j,k)-U(i,j,k)) + &
                              Viscosity(i,j,k)*(U(i,j,k)-U(i-1,j,k)))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k=1,Vnz
      do j=1,Vny
        do i=1,Vnx
          V_r(i,j,k) = V_r(i,j,k) + &
                       Ay * (Viscosity(i,j+1,k)*(V(i,j+1,k)-V(i,j,k)) + &
                             Viscosity(i,j,k)*(V(i,j,k)-V(i,j-1,k)))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp do
    do k=1,Wnz
      do j=1,Wny
        do i=1,Wnx
          W_r(i,j,k) = W_r(i,j,k) + &
                       Az * (Viscosity(i,j,k+1)*(W(i,j,k+1)-W(i,j,k)) + &
                             Viscosity(i,j,k)*(W(i,j,k)-W(i,j,k-1)))
        end do
      end do
    end do
    !$omp end do nowait
    !$omp end parallel
  end subroutine AddSubgridRMS


  subroutine OutputTimeSeries
    character(5) :: prob
    real(knd) :: S,S2,nom,denom
    integer :: i,j,k,unit

    call newunit(unit)
    
    do k = 1,size(probes)

      write(prob,"(i0)") probes(k)%number

      open(unit,file=trim(output_dir)//"Prtimep"//trim(prob)//".txt")
      do j = 0,endstep
       write(unit,*) times(j),Pr_time(k,j)
      end do
      close(unit)

      open(unit,file=trim(output_dir)//"Utimep"//trim(prob)//".txt")
      do j = 0,endstep
       write(unit,*) times(j),U_time(k,j)
      end do
      close(unit)

      open(unit,file=trim(output_dir)//"Vtimep"//trim(prob)//".txt")
      do j = 0,endstep
       write(unit,*) times(j),V_time(k,j)
      end do
      close(unit)

      open(unit,file=trim(output_dir)//"Wtimep"//trim(prob)//".txt")
      do j = 0,endstep
       write(unit,*) times(j),W_time(k,j)
      end do
      close(unit)
      
      if (store%probes_fluxes==1) then
        open(unit,file=trim(output_dir)//"stresstimep"//trim(prob)//".txt")
        do j = 0,endstep
         write(unit,'(7(2x,g13.6))') times(j),momentum_fluxes_time(:,k,j)
        end do
        close(unit)
        open(unit,file=trim(output_dir)//"sgstresstimep"//trim(prob)//".txt")
        do j = 0,endstep
         write(unit,'(7(2x,g13.6))') times(j),momentum_fluxes_time(:,k,j)
        end do
        close(unit)
      end if

      if (enable_buoyancy) then
        open(unit,file=trim(output_dir)//"temptimep"//trim(prob)//".txt")
        do j = 0,endstep
         write(unit,*) times(j),temp_time(k,j)
        end do
        close(unit)
      end if

      if (enable_moisture) then
        open(unit,file=trim(output_dir)//"moisttimep"//trim(prob)//".txt")
        do j = 0,endstep
         write(unit,*) times(j),moist_time(k,j)
        end do
        close(unit)
      end if

    end do

    do k = 1,size(scalar_probes)

      write(prob,"(i0)") scalar_probes(k)%number

      if (num_of_scalars>0) then
        open(unit,file=trim(output_dir)//"scaltimep"//trim(prob)//".txt")
        do j = 1,endstep
         write(unit,'(99(g0,tr3))') times(j),scalp_time(:,k,j)
        end do
        close(unit)
      end if

    end do

    if (store%delta_time==1) then
      open(unit,file=trim(output_dir)//"delta_time.txt")
      do j = 1,endstep
       write(unit,*) times(j),delta_time(j)
      end do
      close(unit)
    end if

    if (store%tke==1) then
      open(unit,file=trim(output_dir)//"tke.txt")
      do j = 0,endstep
       write(unit,*) times(j),tke(j)
      end do
      close(unit)
    end if

    if (store%tke==1.and.store%dissip==1) then
     open(unit,file=trim(output_dir)//"dissip.txt")
      do j = 1,endstep
       write(unit,*) times(j),dissip(j)
      end do
      close(unit)
    end if

    if (wallmodeltype>0.and.display%ustar==1) then
      open(unit,file=trim(output_dir)//"Retau.txt")
      do j = 0,endstep
       write(unit,*) times(j),ustar(:,j)
      end do
      close(unit)
    end if

    if (wallmodeltype>0.and.enable_buoyancy.and.TempBtype(Bo)==DIRICHLET.and.store%tstar==1) then
      open(unit,file=trim(output_dir)//"tflux.txt")
      do j = 0,endstep
       write(unit,*) times(j),tstar(:,j)
      end do
      close(unit)
    end if


    if (num_of_scalars>0.and.store%scalsum_time==1) then
      open(unit,file=trim(output_dir)//"scalsumtime.txt")
      do j = 1,endstep
       write(unit,*) times(j),scalsum_time(:,j)
      end do
      close(unit)
    end if

    if (num_of_scalars>0.and.store%scaltotsum_time==1) then
      open(unit,file=trim(output_dir)//"scaltotsumtime.txt")
      do j = 1,endstep
       write(unit,*) times(j),sum(scalsum_time(:,j))
      end do
      close(unit)
    end if
    
  end subroutine OutputTimeSeries
  
  
  

  subroutine OutputProfiles
    character(5) :: prob
    real(knd) :: S,S2,nom,denom
    integer :: i,j,k,unit

    call newunit(unit)

    if (store%BLprofiles==1.and.averaging==1) then

       open(unit,file=trim(output_dir)//"profu.txt")
       do k = 1,Unz
        write(unit,*) zPr(k),profuavg(k)
       end do
       close(unit)

       open(unit,file=trim(output_dir)//"profv.txt")
       do k = 1,Vnz
        write(unit,*) zPr(k),profvavg(k)
       end do
       close(unit)

       open(unit,file=trim(output_dir)//"profuu.txt")
       do k = 1,Unz
        write(unit,*) zPr(k),profuuavg(k)
       end do
       close(unit)

       open(unit,file=trim(output_dir)//"profvv.txt")
       do k = 1,Vnz
        write(unit,*) zPr(k),profvvavg(k)
       end do
       close(unit)

       open(unit,file=trim(output_dir)//"profww.txt")
       do k = 1,Wnz
        write(unit,*) zW(k),profwwavg(k)
       end do
       close(unit)

       open(unit,file=trim(output_dir)//"profuw.txt")
       do k = 0,Prnz
        write(unit,*) zW(k),profuwavg(k),profuwsgsavg(k)
       end do
       close(unit)

       open(unit,file=trim(output_dir)//"profvw.txt")
       do k = 0,Prnz
        write(unit,*) zW(k),profvwavg(k),profvwsgsavg(k)
       end do
       close(unit)

       if (enable_buoyancy) then

          open(unit,file=trim(output_dir)//"proftemp.txt")
          do k = 1,Prnz
           write(unit,*) zPr(k),proftempavg(k)
          end do
          close(unit)

          open(unit,file=trim(output_dir)//"proftempfl.txt")
          do k = 0,Prnz
           write(unit,*) zW(k),proftempflavg(k),proftempflsgsavg(k)
          end do
          close(unit)

          open(unit,file=trim(output_dir)//"proftt.txt")
          do k = 1,Prnz
           write(unit,*) zPr(k),profttavg(k)
          end do
          close(unit)

          if (allocated(U_avg).and.allocated(V_avg).and.allocated(Temperature_avg)) then
            open(unit,file=trim(output_dir)//"profRig.txt")
            do k = 1,Prnz
             S = 0
             do j = 1,Prny
              do i = 1,Prnx
               S = S + Rig(i,j,k,U_avg,V_avg,Temperature_avg)
              end do
             end do
             S = S / (Prnx*Prny)
             write(unit,*) zPr(k),S
            end do
            close(unit)

            open(unit,file=trim(output_dir)//"profRf.txt")
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

             write(unit,*) zPr(k),S
            end do
            close(unit)
          end if !allocated avgs

       end if !enable_buoyancy == 1


       if (enable_moisture) then

          open(unit,file=trim(output_dir)//"profmoist.txt")
          do k = 1,Prnz
           write(unit,*) zPr(k),profmoistavg(k)
          end do
          close(unit)

          open(unit,file=trim(output_dir)//"profmoistfl.txt")
          do k = 0,Prnz
           write(unit,*) zW(k),profmoistflavg(k),profmoistflsgsavg(k)
          end do
          close(unit)

          open(unit,file=trim(output_dir)//"profmm.txt")
          do k = 1,Prnz
           write(unit,*) zPr(k),profmmavg(k)
          end do
          close(unit)

       end if !enable_moisture == 1


       if (num_of_scalars>0) then

          open(unit,file=trim(output_dir)//"profscal.txt")
          do k = 1,Prnz
           write(unit,*) zPr(k),(profscalavg(i,k), i = 1,num_of_scalars)
          end do
          close(unit)

          open(unit,file=trim(output_dir)//"profscalfl.txt")
          do k = 0,Prnz
           write(unit,*) zW(k),(profscalflavg(i,k),profscalflsgsavg(i,k), i = 1,num_of_scalars)
          end do
          close(unit)

          open(unit,file=trim(output_dir)//"profss.txt")
          do k = 1,Prnz
           write(unit,*) zPr(k),(profssavg(i,k), i = 1,num_of_scalars)
          end do
          close(unit)

       end if !num_of_scalars>0

    end if !store%BLprofiles

  end subroutine OutputProfiles







  subroutine OutputOut(U,V,W,Pr,Temperature,Moisture)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd),dimension(1:,1:,1:),contiguous,intent(in) :: Pr
    real(knd),contiguous,intent(in) :: Temperature(-1:,-1:,-1:)
    real(knd),contiguous,intent(in) :: Moisture(-1:,-1:,-1:)
    character(70) :: str
    integer i,j,k,unit
    real(real32),allocatable :: tmp(:,:,:,:)

    if (store%out==1) then

       call newunit(unit);

       open(unit,file=trim(output_dir)//"out.vtk", &
         access='stream',status='replace',form="unformatted",action="write")

       write(unit) "# vtk DataFile Version 2.0",lf
       write(unit) "CLMM output file",lf
       write(unit) "BINARY",lf
       write(unit) "DATASET RECTILINEAR_GRID",lf
       str="DIMENSIONS"
       write(str(12:),*) Prnx,Prny,Prnz
       write(unit) str,lf
       str="X_COORDINATES"
       write(str(15:),'(i5,2x,a)') Prnx,"float"
       write(unit) str,lf
       write(unit) BigEnd(real(xPr(1:Prnx), real32)),lf
       str="Y_COORDINATES"
       write(str(15:),'(i5,2x,a)') Prny,"float"
       write(unit) str,lf
       write(unit) BigEnd(real(yPr(1:Prny), real32)),lf
       str="Z_COORDINATES"
       write(str(15:),'(i5,2x,a)') Prnz,"float"
       write(unit) str,lf
       write(unit) BigEnd(real(zPr(1:Prnz), real32)),lf
       str="POINT_DATA"
       write(str(12:),*) Prnx*Prny*Prnz
       write(unit) str,lf

       if (store%out_Pr==1) then
         write(unit) "SCALARS p float",lf
         write(unit) "LOOKUP_TABLE default",lf

         write(unit) BigEnd(real(Pr(1:Prnx,1:Prny,1:Prnz), real32))

         write(unit) lf
       end if

       if (enable_buoyancy.and.store%out_temperature==1) then
         write(unit) "SCALARS temperature float",lf
         write(unit) "LOOKUP_TABLE default",lf

         write(unit) BigEnd(real(Temperature(1:Prnx,1:Prny,1:Prnz), real32))

         write(unit) lf
       end if

       if (enable_moisture.and.store%out_moisture==1) then
         write(unit) "SCALARS moisture float",lf
         write(unit) "LOOKUP_TABLE default",lf

         write(unit) BigEnd(real(Moisture(1:Prnx,1:Prny,1:Prnz), real32))

         write(unit) lf
       end if

       if (store%out_Prtype==1) then
         write(unit) "SCALARS ptype float",lf
         write(unit) "LOOKUP_TABLE default",lf

         write(unit) BigEnd(real(Prtype(1:Prnx,1:Prny,1:Prnz), real32))

         write(unit) lf
       end if

       if (store%out_divergence==1) then
         write(unit) "SCALARS div float",lf
         write(unit) "LOOKUP_TABLE default",lf

         if (gridtype==uniformgrid) then
           write(unit) BigEnd(real((U(1:Prnx,1:Prny,1:Prnz) - &
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
               write(unit) BigEnd(real((U(i,j,k)-U(i-1,j,k))/(dxPr(i)) + &
                                        (V(i,j,k)-V(i,j-1,k))/(dyPr(j)) + &
                                        (W(i,j,k)-W(i,j,k-1))/(dzPr(k)), real32))
             end do
            end do
           end do
         end if

         write(unit) lf
       end if

       if (store%out_lambda2==1) then
         allocate(tmp(1:Prnx,1:Prny,1:Prnz,1))

         write(unit) "SCALARS lambda2 float",lf
         write(unit) "LOOKUP_TABLE default",lf

         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
             tmp(i,j,k,1) = BigEnd(real(Lambda2(i,j,k,U,V,W), real32))
           end do
          end do
         end do

         write(unit) tmp

         write(unit) lf
         deallocate(tmp)
       end if

       if (store%out_viscosity==1) then
         write(unit) "SCALARS visc float",lf
         write(unit) "LOOKUP_TABLE default",lf

         write(unit) BigEnd(real(Viscosity(1:Prnx,1:Prny,1:Prnz), real32))

         write(unit) lf
       end if

       if (store%out_U==1.or.store%out_vorticity==1) allocate(tmp(1:3,1:Prnx,1:Prny,1:Prnz))

       if (store%out_U==1) then
         write(unit) "VECTORS u float",lf

         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
             tmp(1:3,i,j,k) = [ BigEnd(real((U(i,j,k)+U(i-1,j,k))/2._knd, real32)), &
                                BigEnd(real((V(i,j,k)+V(i,j-1,k))/2._knd, real32)), &
                                BigEnd(real((W(i,j,k)+W(i,j,k-1))/2._knd, real32)) ]
           end do
          end do
         end do

         write(unit) tmp

         write(unit) lf
       end if

       if (store%out_vorticity==1) then
         write(unit) "VECTORS vort float",lf

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

         write(unit) tmp

         write(unit) lf
       end if
       close(unit)
    end if  !store%out
  end subroutine OutputOut








  subroutine OutputScalars(Scalar)
    real(knd),dimension(-1:,-1:,-1:,:),contiguous,intent(in) :: Scalar
    character(70) :: str
    real(knd),dimension(:,:,:),allocatable :: depos
    character(8) ::  scalname
    integer :: l,unit

    scalname = "scalar"

    if (num_of_scalars>0) then
      if (store%scalars==1) then

          call newunit(unit)

          open(unit,file=trim(output_dir)//"scalars.vtk", &
            access='stream',status='replace',form="unformatted",action="write")

          write(unit) "# vtk DataFile Version 2.0",lf
          write(unit) "CLMM output file",lf
          write(unit) "BINARY",lf
          write(unit) "DATASET RECTILINEAR_GRID",lf
          str="DIMENSIONS"
          write(str(12:),*) Prnx,Prny,Prnz
          write(unit) str,lf
          str="X_COORDINATES"
          write(str(15:),'(i5,2x,a)') Prnx,"float"
          write(unit) str,lf
          write(unit) BigEnd(real(xPr(1:Prnx), real32)),lf
          str="Y_COORDINATES"
          write(str(15:),'(i5,2x,a)') Prny,"float"
          write(unit) str,lf
          write(unit) BigEnd(real(yPr(1:Prny), real32)),lf
          str="Z_COORDINATES"
          write(str(15:),'(i5,2x,a)') Prnz,"float"
          write(unit) str,lf
          write(unit) BigEnd(real(zPr(1:Prnz), real32)),lf
          str="POINT_DATA"
          write(str(12:),*) Prnx*Prny*Prnz
          write(unit) str,lf

          do l = 1,num_of_scalars
            write(scalname(7:8),"(I2.2)") l
            write(unit) "SCALARS ", scalname , " float",lf
            write(unit) "LOOKUP_TABLE default",lf

            write(unit) BigEnd(real(Scalar(1:Prnx,1:Prny,1:Prnz,l), real32))

            write(unit) lf
          end do
          close(unit)
      end if !store%scalars

      if (computedeposition>0.and.store%deposition==1) then

          allocate(depos(1:Prnx,1:Prny,num_of_scalars))
          depos = GroundDeposition()

          call newunit(unit)

          open(unit,file=trim(output_dir)//"deposition.vtk", &
            access='stream',status='replace',form="unformatted",action="write")

          write(unit) "CLMM output file",lf
          write(unit) "BINARY",lf
          write(unit) "DATASET RECTILINEAR_GRID",lf
          str="DIMENSIONS"
          write(str(12:),*) Prnx,Prny,1
          write(unit) str,lf
          str="X_COORDINATES"
          write(str(15:),'(i5,2x,a)') Prnx,"float"
          write(unit) str,lf
          write(unit) BigEnd(real(xPr(1:Prnx), real32)),lf
          str="Y_COORDINATES"
          write(str(15:),'(i5,2x,a)') Prny,"float"
          write(unit) str,lf
          write(unit) BigEnd(real(yPr(1:Prny), real32))
          str="Z_COORDINATES"
          write(str(15:),'(i5,2x,a)') 1,"float"
          write(unit) str,lf
          write(unit) BigEnd(real(zW(0), real32)),lf
          str="POINT_DATA"
          write(str(12:),*) Prnx*Prny
          write(unit) str,lf

          do l = 1,num_of_scalars
            write(scalname(7:8),"(I2.2)") l
            write(unit) "SCALARS ", scalname , " float",lf
            write(unit) "LOOKUP_TABLE default",lf

            write(unit) BigEnd(real(depos(1:Prnx,1:Prny,l), real32))

            write(unit) lf
          end do

          close(unit)

          deallocate(depos)

      end if  !store%deposition
    end if  !num_of_scalars
  end subroutine OutputScalars










  subroutine OutputAvg(U,V,W,U_r,V_r,W_r,Pr,Temperature,Moisture)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in)    :: U,V,W
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout) :: U_r,V_r,W_r
    real(knd),dimension(1:,1:,1:),contiguous,intent(in)    :: Pr
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in)    :: Temperature
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in)    :: Moisture
    character(70) :: str
    character(8) ::  scalname="scalar00"
    integer i,j,k,l,unit
    real(real32),allocatable :: tmp(:,:,:,:)

    if (averaging==1) then
      if (store%avg_U_rms>0) then
        U_r = U_r - U**2
        V_r = V_r - V**2
        W_r = W_r - W**2
      end if

      if (store%avg==1) then
        allocate(tmp(1:3,1:Prnx,1:Prny,1:Prnz))

        call newunit(unit)

        open(unit,file=trim(output_dir)//"avg.vtk", &
          access='stream',status='replace',form="unformatted",action="write")

        write(unit) "# vtk DataFile Version 2.0",lf
        write(unit) "CLMM output file",lf
        write(unit) "BINARY",lf
        write(unit) "DATASET RECTILINEAR_GRID",lf
        str="DIMENSIONS"
        write(str(12:),*) Prnx,Prny,Prnz
        write(unit) str,lf
        str="X_COORDINATES"
        write(str(15:),'(i5,2x,a)') Prnx,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(xPr(1:Prnx), real32)),lf
        str="Y_COORDINATES"
        write(str(15:),'(i5,2x,a)') Prny,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(yPr(1:Prny), real32)),lf
        str="Z_COORDINATES"
        write(str(15:),'(i5,2x,a)') Prnz,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(zPr(1:Prnz), real32)),lf
        str="POINT_DATA"
        write(str(12:),*) Prnx*Prny*Prnz
        write(unit) str,lf

        if (store%avg_Pr==1) then
          write(unit) "SCALARS p float",lf
          write(unit) "LOOKUP_TABLE default",lf

          write(unit) BigEnd(real(Pr(1:Prnx,1:Prny,1:Prnz), real32))

          write(unit) lf
        end if

        if (store%avg_Prtype==1) then
          write(unit) "SCALARS ptype float",lf
          write(unit) "LOOKUP_TABLE default",lf

          write(unit) BigEnd(real(Prtype(1:Prnx,1:Prny,1:Prnz), real32))

          write(unit) lf
        end if

        if (enable_buoyancy.and.store%avg_temperature==1) then
          write(unit) "SCALARS temperature float",lf
          write(unit) "LOOKUP_TABLE default",lf

          write(unit) BigEnd(real(Temperature(1:Prnx,1:Prny,1:Prnz), real32))

          write(unit) lf
        end if

        if (enable_moisture.and.store%avg_moisture==1) then
          write(unit) "SCALARS moisture float",lf
          write(unit) "LOOKUP_TABLE default",lf

          write(unit) BigEnd(real(Moisture(1:Prnx,1:Prny,1:Prnz), real32))

          write(unit) lf
        end if

        if (btest(store%avg_U,0)) then
          write(unit) "VECTORS u float",lf

          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              tmp(:,i,j,k) = [ BigEnd(real((U(i,j,k)+U(i-1,j,k))/2._knd, real32)), &
                               BigEnd(real((V(i,j,k)+V(i,j-1,k))/2._knd, real32)), &
                               BigEnd(real((W(i,j,k)+W(i,j,k-1))/2._knd, real32)) ]
            end do
           end do
          end do

          write(unit) tmp

          write(unit) lf
        end if

        if (btest(store%avg_U_rms,0)) then
          write(unit) "VECTORS u_rms float",lf

          do k = 1,Prnz
           do j = 1,Prny
            do i = 1,Prnx
              tmp(:,i,j,k) = [ BigEnd(real((U_r(i,j,k)+U_r(i-1,j,k))/2._knd, real32)), &
                               BigEnd(real((V_r(i,j,k)+V_r(i,j-1,k))/2._knd, real32)), &
                               BigEnd(real((W_r(i,j,k)+W_r(i,j,k-1))/2._knd, real32)) ]
            end do
           end do
          end do

          write(unit) tmp

          write(unit) lf
        end if

        if (store%avg_vorticity==1) then
          write(unit) "VECTORS vort float",lf

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

          write(unit) tmp

          write(unit) lf
        end if

        close(unit)

      end if !store%avg

      if (btest(store%avg_U,1)) &
        call OutputUVW(U,V,W,trim(output_dir)//"Uavg.vtk",trim(output_dir)//"Vavg.vtk",trim(output_dir)//"Wavg.vtk",.true.)

      if (btest(store%avg_U_rms,1)) &
        call OutputUVW(U_r,V_r,W_r,trim(output_dir)//"Urms.vtk",trim(output_dir)//"Vrms.vtk",trim(output_dir)//"Wrms.vtk",.true.)

    end if !averaging

  end subroutine OutputAvg

  subroutine OutputScalarStats(S_avg,S_max,S_int)
    real(knd),dimension(-1:,-1:,-1:,:),contiguous,intent(in) :: S_avg,S_max,S_int
    character(70) :: str
    integer unit
    
    if (averaging==1.and. &
         (store%scalars_avg==1.or. &
          store%scalars_max==1.or. &
          store%scalars_intermitency==1) &
        .and. num_of_scalars>0) then
       
      open(newunit=unit,file=trim(output_dir)//"scalars_avg.vtk", &
        access='stream',status='replace',form="unformatted",action="write")

      write(unit) "# vtk DataFile Version 2.0",lf
      write(unit) "CLMM output file",lf
      write(unit) "BINARY",lf
      write(unit) "DATASET RECTILINEAR_GRID",lf
      str="DIMENSIONS"
      write(str(12:),*) Prnx,Prny,Prnz
      write(unit) str,lf
      str="X_COORDINATES"
      write(str(15:),'(i5,2x,a)') Prnx,"float"
      write(unit) str,lf
      write(unit) BigEnd(real(xPr(1:Prnx), real32)),lf
      str="Y_COORDINATES"
      write(str(15:),'(i5,2x,a)') Prny,"float"
      write(unit) str,lf
      write(unit) BigEnd(real(yPr(1:Prny), real32)),lf
      str="Z_COORDINATES"
      write(str(15:),'(i5,2x,a)') Prnz,"float"
      write(unit) str,lf
      write(unit) BigEnd(real(zPr(1:Prnz), real32)),lf
      str="POINT_DATA"
      write(str(12:),*) Prnx*Prny*Prnz
      write(unit) str,lf

      if (store%scalars_avg==1) then
        call aux(S_avg,'_avg')
      end if
      
      if (store%scalars_max==1) then
        call aux(S_max,'_max')
      end if
      
      if (store%scalars_intermitency==1) then
        call aux(S_int,'_intermitency')
      end if
      
      close(unit)

    end if
    
  contains
  
    subroutine aux(S,suff)
      real(knd),dimension(-1:,-1:,-1:,:),contiguous,intent(in) :: S
      character(*),intent(in) :: suff
      integer :: l
      character(8) ::  scalname="scalar00"

      do l = 1,num_of_scalars
        write(scalname(7:8),"(I2.2)") l
        write(unit) "SCALARS ", scalname//suff, " float",lf
        write(unit) "LOOKUP_TABLE default",lf

        write(unit) BigEnd(real(S(1:Prnx,1:Prny,1:Prnz,l), real32))

        write(unit) lf
      end do
    end subroutine
  end subroutine OutputScalarStats




  subroutine OutputUVW(U,V,W,fnameU,fnameV,fnameW,avg_mode_arg)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
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

        open(unit,file=fnameU, &
          access='stream',status='replace',form="unformatted",action="write")

        write(unit) "# vtk DataFile Version 2.0",lf
        write(unit) "CLMM output file",lf
        write(unit) "BINARY",lf
        write(unit) "DATASET RECTILINEAR_GRID",lf
        str="DIMENSIONS"
        write(str(12:),*) Unx,Uny,Unz
        write(unit) str,lf
        str="X_COORDINATES"
        write(str(15:),'(i5,2x,a)') Unx,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(xU(1:Unx), real32)),lf
        str="Y_COORDINATES"
        write(str(15:),'(i5,2x,a)') Uny,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(yPr(1:Uny), real32)),lf
        str="Z_COORDINATES"
        write(str(15:),'(i5,2x,a)') Unz,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(zPr(1:Unz), real32)),lf
        str="POINT_DATA"
        write(str(12:),*) Unx*Uny*Unz
        write(unit) str,lf


        write(unit) "SCALARS U float",lf
        write(unit) "LOOKUP_TABLE default",lf

        write(unit) BigEnd(real(U(1:Unx,1:Uny,1:Unz), real32))

        write(unit) lf

        if (store%U_interp==1) then
          write(unit) "SCALARS Utype float",lf
          write(unit) "LOOKUP_TABLE default",lf

          write(unit) BigEnd(real(Utype(1:Unx,1:Uny,1:Unz), real32))

          write(unit) lf
        end if

        close(unit)
    end if

    if (store%V==1.or.avg_mode) then

        call newunit(unit)

        open(unit,file=fnameV, &
          access='stream',status='replace',form="unformatted",action="write")

        write(unit) "# vtk DataFile Version 2.0",lf
        write(unit) "CLMM output file",lf
        write(unit) "BINARY",lf
        write(unit) "DATASET RECTILINEAR_GRID",lf
        str="DIMENSIONS"
        write(str(12:),*) Vnx,Vny,Vnz
        write(unit) str,lf
        str="X_COORDINATES"
        write(str(15:),'(i5,2x,a)') Vnx,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(xPr(1:Vnx), real32)),lf
        str="Y_COORDINATES"
        write(str(15:),'(i5,2x,a)') Vny,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(yV(1:Vny), real32)),lf
        str="Z_COORDINATES"
        write(str(15:),'(i5,2x,a)') Vnz,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(zPr(1:Vnz), real32)),lf
        str="POINT_DATA"
        write(str(12:),*) Vnx*Vny*Vnz
        write(unit) str,lf


        write(unit) "SCALARS V float",lf
        write(unit) "LOOKUP_TABLE default",lf

        write(unit) BigEnd(real(V(1:Vnx,1:Vny,1:Vnz), real32))

        write(unit) lf

        if (store%V_interp==1) then
          write(unit) "SCALARS Vtype float",lf
          write(unit) "LOOKUP_TABLE default",lf

          write(unit) BigEnd(real(Vtype(1:Vnx,1:Vny,1:Vnz), real32))

          write(unit) lf
        end if

        close(unit)
    end if !store%V

    if (store%W==1.or.avg_mode) then

        call newunit(unit)

        open(unit,file=fnameW, &
          access='stream',status='replace',form="unformatted",action="write")

        write(unit) "# vtk DataFile Version 2.0",lf
        write(unit) "CLMM output file",lf
        write(unit) "BINARY",lf
        write(unit) "DATASET RECTILINEAR_GRID",lf
        str="DIMENSIONS"
        write(str(12:),*) Wnx,Wny,Wnz
        write(unit) str,lf
        str="X_COORDINATES"
        write(str(15:),'(i5,2x,a)') Wnx,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(xPr(1:Wnx), real32)),lf
        str="Y_COORDINATES"
        write(str(15:),'(i5,2x,a)') Wny,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(yPr(1:Wny), real32)),lf
        str="Z_COORDINATES"
        write(str(15:),'(i5,2x,a)') Wnz,"float"
        write(unit) str,lf
        write(unit) BigEnd(real(zW(1:Wnz), real32)),lf
        str="POINT_DATA"
        write(str(12:),*) Wnx*Wny*Wnz
        write(unit) str,lf


        write(unit) "SCALARS W float",lf
        write(unit) "LOOKUP_TABLE default",lf

        write(unit) BigEnd(real(W(1:Wnx,1:Wny,1:Wnz), real32))

        write(unit) lf

        if (store%W_interp==1) then
          write(unit) "SCALARS Wtype float",lf
          write(unit) "LOOKUP_TABLE default",lf

          write(unit) BigEnd(real(Wtype(1:Wnx,1:Wny,1:Wnz), real32))

          write(unit) lf
        end if

        close(unit)
    end if !store%W

    if (store%U_interp/=0 .and. .not.avg_mode) then

      call newunit(unit)

      open(unit,file=trim(output_dir)//"Uinterp.txt")
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

      open(unit,file=trim(output_dir)//"Vinterp.txt")
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

      open(unit,file=trim(output_dir)//"Winterp.txt")
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

      open(unit,file=trim(output_dir)//"Scinterp.txt")
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


  subroutine OutputAvgFluxes
    real(knd),allocatable :: Scalar_fl_U_adv(:,:,:,:)
    real(knd),allocatable :: Scalar_fl_V_adv(:,:,:,:)
    real(knd),allocatable :: Scalar_fl_W_adv(:,:,:,:)
    real(knd),allocatable :: Scalar_fl_U_turb(:,:,:,:)
    real(knd),allocatable :: Scalar_fl_V_turb(:,:,:,:)
    real(knd),allocatable :: Scalar_fl_W_turb(:,:,:,:)
    integer :: i
   
    if (num_of_scalars>0.and.store%avg_flux_scalar==1) then

      !FIXME: just to allocate, mold= does not work correctly as of gcc 4.8
!       allocate(Scalar_fl_U_adv,mold=Scalar_fl_U_avg)
      Scalar_fl_U_adv = Scalar_fl_U_avg
      Scalar_fl_V_adv = Scalar_fl_V_avg
      Scalar_fl_W_adv = Scalar_fl_W_avg
      Scalar_fl_U_adv = 0
      Scalar_fl_V_adv = 0
      Scalar_fl_W_adv = 0

      do i=1,num_of_scalars
        call AddScalarAdvVector(Scalar_fl_U_adv(:,:,:,i), &
                             Scalar_fl_V_adv(:,:,:,i), &
                             Scalar_fl_W_adv(:,:,:,i), &
                             Scalar_avg(:,:,:,i), &
                             U_avg,V_avg,W_avg, &
                             1._knd)
      end do

      Scalar_fl_U_turb = Scalar_fl_U_avg - Scalar_fl_U_adv
      Scalar_fl_V_turb = Scalar_fl_V_avg - Scalar_fl_V_adv
      Scalar_fl_W_turb = Scalar_fl_W_avg - Scalar_fl_W_adv
      

      call SaveScalarVTKFluxes(Scalar_fl_U_avg, &
                               Scalar_fl_U_adv, &
                               Scalar_fl_U_turb, &
                               trim(output_dir)//"scalflu.vtk",xU,yPr,zPr)
      call SaveScalarVTKFluxes(Scalar_fl_V_avg, &
                               Scalar_fl_V_adv, &
                               Scalar_fl_V_turb, &
                               trim(output_dir)//"scalflv.vtk",xPr,yV,zPr)
      call SaveScalarVTKFluxes(Scalar_fl_W_avg, &
                               Scalar_fl_W_adv, &
                               Scalar_fl_W_turb, &
                               trim(output_dir)//"scalflw.vtk",xPr,yPr,zW)

    end if

  end subroutine

  subroutine SaveScalarVTKFluxes(Avg,Adv,Turb,file_name,x,y,z)
    real(knd),dimension(:,:,:,:),contiguous,intent(in) :: Avg,Adv,Turb
    character(*),intent(in) :: file_name
    real(knd),allocatable,intent(in) :: x(:),y(:),z(:)
    integer :: nx,ny,nz
    character(70) :: str
    character(13) ::  scalname
    integer :: l,unit

    nx = size(Avg,1)
    ny = size(Avg,2)
    nz = size(Avg,3)

    call newunit(unit)

    open(unit,file=trim(file_name), &
      access='stream',status='replace',form="unformatted",action="write")

    write(unit) "# vtk DataFile Version 2.0",lf
    write(unit) "CLMM output file",lf
    write(unit) "BINARY",lf
    write(unit) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write(str(12:),*) nx,ny,nz
    write(unit) str,lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(x(1:nx), real32)),lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(y(1:ny), real32)),lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(z(1:nz), real32)),lf
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(unit) str,lf

    scalname = "scalfl-avg-"

    do l = 1,num_of_scalars
      write(scalname(12:13),"(I2.2)") l
      write(unit) "SCALARS ", scalname , " float",lf
      write(unit) "LOOKUP_TABLE default",lf

      write(unit) BigEnd(real(Avg(:,:,:,l), real32))

      write(unit) lf
    end do

    scalname = "scalfl-adv-"

    do l = 1,num_of_scalars
      write(scalname(12:13),"(I2.2)") l
      write(unit) "SCALARS ", scalname , " float",lf
      write(unit) "LOOKUP_TABLE default",lf

      write(unit) BigEnd(real(Adv(:,:,:,l), real32))

      write(unit) lf
    end do

    scalname = "scalfl-trb-"

    do l = 1,num_of_scalars
      write(scalname(12:13),"(I2.2)") l
      write(unit) "SCALARS ", scalname , " float",lf
      write(unit) "LOOKUP_TABLE default",lf

      write(unit) BigEnd(real(Turb(:,:,:,l), real32))

      write(unit) lf
    end do
    close(unit)
  end subroutine SaveScalarVTKFluxes



  subroutine Output(U,V,W,Pr,Temperature,Moisture,Scalar)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(inout) :: U,V,W
    real(knd),contiguous,intent(inout) :: Pr(1:,1:,1:)
    real(knd),contiguous,intent(in) :: Temperature(-1:,-1:,-1:)
    real(knd),contiguous,intent(in) :: Moisture(-1:,-1:,-1:)
    real(knd),contiguous,intent(in) :: Scalar(-1:,-1:,-1:,:)
    integer i

    call BoundU(1,U,Uin)
    call BoundU(2,V,Vin)
    call BoundU(3,W,Win)

    call OutputOut(U,V,W,Pr,Temperature,Moisture)

    call OutputScalars(Scalar)

    call OutputUVW(U,V,W,trim(output_dir)//"U.vtk",trim(output_dir)//"V.vtk",trim(output_dir)//"W.vtk")
    
    call OutputTimeSeries

    if (averaging==1.and.time>=timeavg1) then

      call OutputAvg(U_avg,V_avg,W_avg,U_rms,V_rms,W_rms,Pr_avg,Temperature_avg,Moisture_avg)
      
      call OutputScalarStats(Scalar_avg,Scalar_max,Scalar_intermitency)

      call OutputProfiles

      call OutputAvgFluxes

    end if

    call FinalizeVTKFrames

    call FinalizeStaggeredFrames

    if (master) write(*,*) "saved"
  end subroutine Output





  subroutine StressProfiles(U,V,W)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd) :: S
    real(knd),allocatable,save ::fp(:),ht(:),gp(:)
    integer   :: i,j,k,l,n
    integer,save :: called = 0

    if (called==0) then
      allocate(fp(0:Prnx+1),ht(0:Prnz+1),gp(0:Prny+1))
      !$omp parallel workshare
      forall (i = 0:Prnx+1)      fp(i)=(xU(i)-xPr(i))/(xPr(i+1)-xPr(i))
      forall (k = 0:Prnz+1)      ht(k)=(zW(k)-zPr(k))/(zPr(k+1)-zPr(k))
      forall (j = 0:Prny+1)      gp(j)=(yV(j)-yPr(j))/(yPr(j+1)-yPr(j))
      !$omp end parallel workshare
      called = 1
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
           S = S-0.25_knd*(Viscosity(i+1,j,k+1)+Viscosity(i+1,j,k)+Viscosity(i,j,k+1)+Viscosity(i,j,k))*(U(i,j,k+1)-U(i,j,k))/dzW(k)
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
           S = S-0.25_knd*(Viscosity(i,j+1,k+1)+Viscosity(i,j+1,k)+Viscosity(i,j,k+1)+Viscosity(i,j,k))*(V(i,j,k+1)-V(i,j,k))/dzW(k)
           n = n + 1
         end if
       end do
      end do
      profvwsgs(k) = S / max(n,1)
    end do
    !$omp end do nowait
    !$omp end parallel
    
  end subroutine StressProfiles
  
  
  subroutine FluxSGSProfiles(W,Temperature,Moisture,Scalar)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: W
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in) :: Temperature
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in) :: Moisture
    real(knd),dimension(-1:,-1:,-1:,1:),contiguous,intent(in) :: Scalar
    real(knd) :: S
    integer   :: i,j,k,l,n

    
    if (enable_buoyancy) call TemperatureFluxSGSProfile(W,Temperature)
    if (enable_moisture) call MoistureFluxSGSProfile(W,Moisture)
    if (num_of_scalars>0) call ScalarFluxSGSProfile(W,Scalar)
  end subroutine


  subroutine TemperatureFluxSGSProfile(W,Temperature)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: W
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in) :: Temperature
    real(knd) :: S
    integer   :: i,j,k,l,n
    !proftempfl is computed directly during advection step

    !$omp parallel do private(i,j,k,n,S)
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
    !$omp end parallel do
  end subroutine

  subroutine MoistureFluxSGSProfile(W,Moisture)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: W
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in) :: Moisture
    real(knd) :: S
    integer   :: i,j,k,l,n
    !profmoistfl is computed directly during advection step

    !$omp parallel do private(i,j,k,n,S)
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
    !$omp end parallel do
  end subroutine
      
  subroutine ScalarFluxSGSProfile(W,Scalar)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: W
    real(knd),dimension(-1:,-1:,-1:,1:),contiguous,intent(in) :: Scalar
    real(knd) :: S
    integer   :: i,j,k,l,n

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
  end subroutine

  


  subroutine BLProfiles(U,V,W,Temperature,Moisture,Scalar)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in) :: Temperature
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in) :: Moisture
    real(knd),dimension(-1:,-1:,-1:,1:),contiguous,intent(in) :: Scalar
    real(knd) :: S
    integer   :: i,j,k,l,n

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

    !$omp parallel private(i,j,k,n,S)
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

    if (enable_buoyancy) then
      !$omp parallel do private(i,j,k,n,S)
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
      !$omp end parallel do
    end if ! size(Temperature)


    if (enable_moisture) then
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



    do l = 1,num_of_scalars
      !$omp parallel do private(i,j,k,n,S)
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
      !$omp end parallel do
    end do

  end subroutine BLProfiles




  subroutine OutputU2(U,V,W)
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    integer unit
    character(70) :: str

    call newunit(unit)

    open(unit,file=trim(output_dir)//"U2.vtk")
    write(unit) "# vtk DataFile Version 2.0",lf
    write(unit) "CLMM output file",lf
    write(unit) "BINARY",lf
    write(unit) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write(str(12:),*) Unx,Uny,Unz
    write(unit) str,lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') Unx,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(xU(1:Unx), real32)),lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') Uny,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(yPr(1:Uny), real32)),lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') Unz,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(zPr(1:Unz), real32)),lf
    str="POINT_DATA"
    write(str(12:),*) Unx*Uny*Unz
    write(unit) str,lf


    write(unit) "SCALARS U float",lf
    write(unit) "LOOKUP_TABLE default",lf

    write(unit) BigEnd(real(U(1:Unx,1:Uny,1:Unz), real32)),lf

    write(unit) lf
    close(unit)


    open(unit,file=trim(output_dir)//"V2.vtk")
    write(unit) "# vtk DataFile Version 2.0",lf
    write(unit) "CLMM output file",lf
    write(unit) "BINARY",lf
    write(unit) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write(str(12:),*) Vnx,Vny,Vnz
    write(unit) str,lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') Vnx,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(xPr(1:Vnx), real32)),lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') Vny,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(yV(1:Vny), real32)),lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') Vnz,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(zPr(1:Vnz), real32)),lf
    str="POINT_DATA"
    write(str(12:),*) Vnx*Vny*Vnz
    write(unit) str,lf


    write(unit) "SCALARS V float",lf
    write(unit) "LOOKUP_TABLE default",lf

    write(unit) BigEnd(real(V(1:Vnx,1:Vny,1:Vnz), real32)),lf

    write(unit) lf
    close(unit)


    open(unit,file=trim(output_dir)//"W2.vtk")
    write(unit) "# vtk DataFile Version 2.0",lf
    write(unit) "CLMM output file",lf
    write(unit) "BINARY",lf
    write(unit) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write(str(12:),*) Wnx,Wny,Wnz
    write(unit) str,lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') Wnx,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(xPr(1:Wnx), real32)),lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') Wny,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(yPr(1:Wny), real32)),lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') Wnz,"float"
    write(unit) str,lf
    write(unit) BigEnd(real(zW(1:Wnz), real32)),lf
    str="POINT_DATA"
    write(str(12:),*) Wnx*Wny*Wnz
    write(unit) str,lf


    write(unit) "SCALARS W float",lf
    write(unit) "LOOKUP_TABLE default",lf

    write(unit) BigEnd(real(W(1:Wnx,1:Wny,1:Wnz), real32)),lf

    write(unit) lf
    close(unit)
  end subroutine OutputU2


  subroutine OUTINLET(U,V,W,Temperature)
    !for output of 2d data for use as an inilet condition later
    real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
    real(knd),dimension(-1:,-1:,-1:),contiguous,intent(in) :: Temperature
    integer,save ::fnum
    integer,save :: called = 0

    if ((time>=timefram1).and.(time<=timefram2+(timefram2-timefram1)/(frames-1))&
        .and.(time>=timefram1+fnum*(timefram2-timefram1)/(frames-1))) then
     if (called==0) then
      open(101,file=trim(output_dir)//"inletframeinfo.unf",form='unformatted',status='replace',action='write')
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

    call GridCoords(mini,minj,mink,(xU(Prnx+1)+xU(0))/2._knd,(yV(Prny+1)+yV(0))/2._knd,(zW(Prnz+1)+zW(0))/2._knd)
    maxi = mini
    minj = 1
    maxj = Prny
    mink = 1
    maxk = Prnz


    fname(1:5)="frame"
    write(fname(6:8),"(I3.3)") n
    fname(9:12)=".unf"
    if (master) write(*,*) "Saving frame:",fname(1:6),"   time:",time

    call newunit(unit)

    open(unit,file = fname,form='unformatted',access='sequential',status='replace',action='write')


    write(unit) U(mini,1:Uny,1:Unz)
    write(unit) V(mini,1:Vny,1:Vnz)
    write(unit) W(mini,1:Wny,1:Wnz)
    if (enable_buoyancy) then
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

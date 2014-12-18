program CLMM

  use Parameters
  use Initial
  use TSteps, only: TMarchRK3, TSteps_Deallocate
  use Scalars, only: Scalars_Deallocate
  use Outputs, only: AllocateOutputs, OutTStep, Output
  use Endianness, only: GetEndianness
#ifdef MPI
  use custom_mpi
#endif

  implicit none




  real(knd),allocatable:: U(:,:,:),V(:,:,:),W(:,:,:),Pr(:,:,:)
  real(knd),allocatable,dimension(:,:,:) :: Temperature !If buoyancy enabled then liquid potential temperature, othrewise potential temperature.
  real(knd),allocatable,dimension(:,:,:) :: Moisture !Total specific humidity q_t.
  real(knd),allocatable,dimension(:,:,:,:) :: Scalar  !last index is a number of scalar (because of paging)

  real(knd) :: delta = 0

  real(tim) :: dt = 0 !time_step

  integer   :: time_step

  real(knd) :: time_step_time

  integer(dbl) :: time_steps_timer_count_1, time_steps_timer_count_2

  call GetEndianness

#ifdef MPI
  call init_custom_mpi
#endif

  if (master) write (*,*) "Reading parameters..."
  call ReadConfiguration


  if (master) write (*,*) "Setting up boundary conditions..."
  call InitBoundaryConditions


  if (master) write (*,*) "Allocating arrays..."
  call AllocateGlobals

  call AllocateOutputs

  if (master) write (*,*) "Setting up initial conditions..."
  call InitialConditions(U,V,W,Pr,Temperature,Moisture,Scalar)

  time = start_time
  time_step = 0

  call OutTStep(U,V,W,Pr,Temperature,Moisture,Scalar,dt,delta)
  
  init_phase = .false.
  run_phase = .true.

  call system_clock(count_rate = timer_rate)

  if (end_time > start_time) then

    if (master) write (*,*) "Computing..."

    do time_step = 1, max_number_of_time_steps

      call system_clock(count = time_steps_timer_count_1)

      if (master) then
!         write (*,*) "-----------"
!         write (*,*)
!         write (*,*)
        write (*,*) "tstep:",time_step
      end if




      call TMarchRK3(U,V,W,Pr,Temperature,Moisture,Scalar,dt,delta)




      time = time + dt


      call OutTStep(U,V,W,Pr,Temperature,Moisture,Scalar,dt,delta)


      call system_clock(count = time_steps_timer_count_2)
      if (master) then
        time_step_time = real(time_steps_timer_count_2 - time_steps_timer_count_1, knd) / &
          real(timer_rate, knd)
        time_steps_time = time_steps_time + time_step_time
!         print *,"time step wall clock time", time_step_time
      end if

      if ((steady==1) .and. (delta<eps)) then
        if (master) write (*,*) "Steady state reached."
        exit
      endif

      if ((steady==0) .and. (time>=end_time)) then
        if (master) write (*,*) "Time limit reached."
        exit
      endif

      if (time_step>=3 .and. dt < abs(CFL*min(dxmin,dymin,dzmin)/Uinlet/10._knd)) then
        if (master) write (*,*) "Solution diverged."
        exit
      endif

    enddo

  endif

  if (master) write(*,*) "Total wall clock time for time steps", time_steps_time
  if (master) write(*,*) "Wall clock time for poisson solver", poisson_solver_time

#ifdef __HMPP
  call GetDataFromGPU(.true.,.true.,.true.,.true.,enable_buoyancy,&
                      U,     V,     W,     Pr,    Temperature)
#endif

  call TSteps_Deallocate
  call Scalars_Deallocate

  if (master) write (*,*) "Saving results..."


  call Output(U,V,W,Pr,Temperature,Moisture,Scalar)


  call DeallocateGlobals
  
#ifdef MPI
  call finalize_custom_mpi
#endif
  contains


   subroutine AllocateGlobals

     allocate(U(-2:Unx+3,-2:Uny+3,-2:Unz+3))
     allocate(V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
     allocate(W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
     allocate(Pr(1:Unx+1,1:Vny+1,1:Wnz+1))

     U = 0
     V = 0
     W = 0
     Pr = 0


     if (enable_buoyancy) then
       allocate(Temperature(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
       Temperature = 0
     else
       allocate(Temperature(0,0,0))
     endif

     if (enable_moisture) then
       allocate(Moisture(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
       Moisture = 0
     else
       allocate(Moisture(0,0,0))
     endif

     if (num_of_scalars>0) then
       allocate(Scalar(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,num_of_scalars))
       Scalar = 0
     else
       allocate(Scalar(0,0,0,0))
     endif


     allocate(Viscosity(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))

     if (enable_buoyancy .or. &
         enable_moisture .or. &
         num_of_scalars>0)      then

       allocate(TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
       TDiff = 0

     else
       allocate(TDiff(0,0,0))
     endif

   end subroutine AllocateGlobals


   subroutine DeallocateGlobals
    !To assist memory leak checking because some compilers do not
    !deallocate main program and module variables automatically
    !(in F2008 they have the SAVE attribute).

    deallocate(U, V, W, Pr)

    deallocate(Temperature, Moisture, Scalar)

    deallocate(Viscosity, TDiff)

    deallocate(xU, yV, zW, dxU, dyV, dzW)
    deallocate(xPr, yPr, zPr, dxPr, dyPr, dzPr)

    deallocate(Utype, Vtype, Wtype, Prtype)

    deallocate(Uin, Vin, Win)

    if (allocated(TempIn)) deallocate(TempIn)
    if (allocated(MoistIn)) deallocate(MoistIn)

    if (allocated(BsideTArr)) deallocate(BsideTArr)
    if (allocated(BsideTFlArr)) deallocate(BsideTFlArr)

    if (allocated(BsideMArr)) deallocate(BsideMArr)
    if (allocated(BsideMFlArr)) deallocate(BsideMFlArr)


   end subroutine


end program CLMM


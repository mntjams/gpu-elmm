program CLMM

  use Parameters
  use Initial
  use TSteps, only: TMarchRK3, TSteps_Deallocate
  use Scalars, only: Scalars_Deallocate
  use Outputs, only: AllocateOutputs, OutTStep, Output

  implicit none




  real(knd),allocatable:: U(:,:,:),V(:,:,:),W(:,:,:),Pr(:,:,:)
  real(knd),allocatable,dimension(:,:,:) :: Temperature !If buoyancy enabled then liquid potential temperature, othrewise potential temperature.
  real(knd),allocatable,dimension(:,:,:) :: Moisture !Total specific humidity q_t.
  real(knd),allocatable,dimension(:,:,:,:) :: Scalar  !last index is a number of scalar (because of paging)

  real(knd) :: delta = 0

  real(knd) :: dt = 0 !time_step

  pi=2.0_knd*acos(0.0_knd)


  write (*,*) "Reading parameters..."
  call ReadConfiguration


  write (*,*) "Setting up boundary conditions..."
  call InitBoundaryConditions


  write (*,*) "Allocating arrays..."
  call AllocateGlobals

  call AllocateOutputs

  write (*,*) "Setting up initial conditions..."
  call InitialConditions(U,V,W,Pr,Temperature,Moisture,Scalar)

  time = start_time
  step = 0

  call OutTStep(U,V,W,Pr,Temperature,Moisture,Scalar,dt,delta)



  if (end_time > start_time) then

    write (*,*) "Computing..."

    do step = 1, max_number_of_time_steps

      write (*,*) "-----------"
      write (*,*) " "
      write (*,*) " "
      write (*,*) "tstep:",step




      call TMarchRK3(U,V,W,Pr,Temperature,Moisture,Scalar,dt,delta)




      time = time + dt


      call OutTStep(U,V,W,Pr,Temperature,Moisture,Scalar,dt,delta)


      if ((steady==1) .and. (delta<eps)) then
        write (*,*) "Steady state reached."
        exit
      endif

      if ((steady==0) .and. (time>=end_time)) then
        write (*,*) "Time limit reached."
        exit
      endif

      if (step>=3 .and. dt < abs(CFL*min(dxmin,dymin,dzmin)/Uinlet/10._knd)) then
        write (*,*) "Solution diverged."
        exit
      endif

    enddo

  endif

#ifdef __HMPP
  call GetDataFromGPU(.true.,.true.,.true.,.true.,enable_buoyancy,&
                      U,     V,     W,     Pr,    Temperature)
#endif

  call TSteps_Deallocate
  call Scalars_Deallocate

  write (*,*) "Saving results..."


  call OUTPUT(U,V,W,Pr,Temperature,Moisture,Scalar)


  call DeallocateGlobals


  contains


   subroutine AllocateGlobals

     allocate(U(-2:Unx+3,-2:Uny+3,-2:Unz+3))
     allocate(V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
     allocate(W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
     allocate(Pr(1:Unx+1,1:Vny+1,1:Wnz+1))

     U = 0
     V = 0
     W = 0
     Pr = 0000


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


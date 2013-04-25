program CLMM

  use TSTEPS, only: TMarchRK3
  use PARAMETERS
  use Initial
  use OUTPUTS, only: AllocateOutputs, OutTStep, Output

  implicit none




  real(knd),allocatable:: U(:,:,:),V(:,:,:),W(:,:,:),Pr(:,:,:)
  real(knd),allocatable,dimension(:,:,:) :: Temperature !If buoyancy enabled then liquid potential temperature, othrewise potential temperature.
  real(knd),allocatable,dimension(:,:,:) :: Moisture !Total specific humidity q_t.
  real(knd),allocatable,dimension(:,:,:,:) :: Scalar  !last index is a number of scalar (because of paging)

  real(knd) :: delta = 0

  pi=2.0_knd*acos(0.0_knd)


  write (*,*) "Reading parameters..."
  call ReadConfiguration


  write (*,*) "Setting up boundary conditions..."
  call InitBoundaryConditions


  write (*,*) "Allocating arrays..."
  call AllocateGlobals


  write (*,*) "Setting up initial conditions..."
  call InitialConditions(U,V,W,Pr,Temperature,Moisture,Scalar)

  time=starttime
  step=0

  call AllocateOutputs

  call OutTStep(U,V,W,Pr,Temperature,Moisture,Scalar,delta)



  if (endtime>starttime) then

    write (*,*) "Computing..."

    do step = 1,nt

      write (*,*) "-----------"
      write (*,*) " "
      write (*,*) " "
      write (*,*) "tstep:",step




      call TMarchRK3(U,V,W,Pr,Temperature,Moisture,Scalar,delta)




      time = time + dt


      call OutTStep(U,V,W,Pr,Temperature,Moisture,Scalar,delta)


      if ((steady==1) .and. (delta<eps)) then
        write (*,*) "Steady state reached."
        exit
      endif

      if ((steady==0) .and. (time>=endtime)) then
        write (*,*) "Time limit reached."
        exit
      endif

      if (step>=3 .and. dt < abs(CFL*min(dxmin,dymin,dzmin)/Uinlet/50._knd)) then
        write (*,*) "Solution diverged."
        exit
      endif

    enddo

  endif

#ifdef __HMPP
  call GetDataFromGPU(.true.,.true.,.true.,.true.,enable_buoyancy==1,&
                      U,     V,     W,     Pr,    Temperature)
#endif


  write (*,*) "Saving results..."

  call OUTPUT(U,V,W,Pr,Temperature,Moisture,Scalar)


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


     if (enable_buoyancy==1) then
       allocate(Temperature(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
       Temperature = 0
     else
       allocate(Temperature(0,0,0))
     endif

     if (enable_moisture==1) then
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


     allocate(Visc(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))

     if (enable_buoyancy==1 .or. &
         enable_moisture==1 .or. &
         num_of_scalars>0)      then

       allocate(TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
       TDiff = 0

     else
       allocate(TDiff(0,0,0))
     endif

   end subroutine AllocateGlobals


end program CLMM


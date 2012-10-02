program CLMM

  use TSTEPS
  use PARAMETERS
  use BOUNDARIES
  use INITIAL
  use OUTPUTS, only: AllocateOutputs, OutTStep, Output
  use GEOMETRIC
  use SCALARS
  use WALLMODELS

  implicit none




  real(KND),allocatable:: U(:,:,:),V(:,:,:),W(:,:,:),Pr(:,:,:)
  real(KND) :: delta = 0

  pi=2.0_KND*acos(0.0_KND)


  write (*,*) "Reading parameters..."
  call ReadParams


  write (*,*) "Setting up boundary conditions..."
  call ReadBounds


  write (*,*) "Allocating arrays..."
  call AllocateGlobals


  write (*,*) "Setting up initial conditions..."
  call InitConds(U,V,W,Pr)

  time=starttime
  step=0

  call AllocateOutputs

  call OutTStep(U,V,W,Pr,delta)



  if (endtime>starttime) then

    write (*,*) "Computing..."

    do step=1,nt

      write (*,*) "-----------"
      write (*,*) " "
      write (*,*) " "
      write (*,*) "tstep:",step




      call TMarchRK3(U,V,W,Pr,Temperature,Scalar,delta)




      time=time+dt


      call OutTStep(U,V,W,Pr,delta)


      if ((steady==1) .and. (delta<eps)) then
        write (*,*) "Steady state reached."
        exit
      endif

      if ((steady==0) .and. (time>=endtime)) then
        write (*,*) "Time limit reached."
        exit
      endif

      if (dt<abs(CFL*dxmin/Uinlet/100._KND)) then
        write (*,*) "Solution diverged."
        exit
      endif

    enddo

  endif

#ifdef __HMPP
  call GetDataFromGPU(.true.,.true.,.true.,.true.,buoyancy==1,&
                      U,     V,     W,     Pr,    Temperature)
#endif

write(*,*) maxval(U(1:Unx,1:Uny,1:Unz)),maxval(V(1:Vnx,1:Vny,1:Vnz)),maxval(W(1:Wnx,1:Wny,1:Wnz)),&
maxval(temperature(1:Prnx,1:Prny,1:Prnz))
write(*,*) sum(U(1:Unx,1:Uny,1:Unz)),sum(V(1:Vnx,1:Vny,1:Vnz)),sum(W(1:Wnx,1:Wny,1:Wnz)),&
sum(temperature(1:Prnx,1:Prny,1:Prnz))
write(*,*) U(Unx/2,Uny/2,1),V(Vnx/2,Vny/2,1),W(Wnx/2,Wny/2,1),&
temperature(Prnx/2,Prny/2,1)



  write (*,*) "Saving results..."

  call OUTPUT(U,V,W,Pr)


  contains


   subroutine AllocateGlobals

     allocate(U(-2:Unx+3,-2:Uny+3,-2:Unz+3))
     allocate(V(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
     allocate(W(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
     allocate(Pr(1:Unx+1,1:Vny+1,1:Wnz+1))
     U=huge(1._KND)
     V=huge(1._KND)
     W=huge(1._KND)
     Pr=0000


     if (buoyancy==1) then
       allocate(temperature(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
       temperature=0
     else
       allocate(temperature(0,0,0))
     endif

     if (computescalars>0) then
       allocate(Scalar(-1:Prnx+2,-1:Prny+2,-1:Prnz+2,computescalars))
       Scalar=huge(1.0_KND)
     else
       allocate(Scalar(0,0,0,0))
     endif


     allocate(Visc(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))

     if (buoyancy>0.or.computescalars>0) then
       allocate(TDiff(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
       TDiff=huge(1.0_KND)
     endif

     if (fullstress==1) then
       allocate(tstress(3,3,-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
     endif
   end subroutine AllocateGlobals


end program CLMM


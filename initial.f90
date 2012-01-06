module INITIAL

  use PARAMETERS
  use LIMITERS, only: limparam, limitertype
  use MULTIGRID, only: SetMGParams
  use MULTIGRID2d, only: SetMGParams2d
  use POISSON
  use BOUNDARIES
  use OUTPUTS, only: store, display, probes, NumProbes
  use SCALARS
  use SMAGORINSKY
  use TURBINLET, only: GetTurbInlet, GetInletFromFile, TLag, Lturby, Lturbz, ustarinlet, transformtensor
  use GEOMETRIC
  use WALLMODELS

  implicit none

  private
  public  ReadParams, Initconds, ReadBounds



contains
 subroutine ReadParams
  integer   lmg,minmglevel,bnx,bny,bnz,mgncgc,mgnpre,mgnpost,mgmaxinnerGSiter,minGPUlevel
  real(KND) mgepsinnerGS
  integer i,io

  open(11,file="main.conf",status="old",action="read")
  read(11,fmt='(/)')
  read(11,*) tempmet
  read(11,fmt='(/)')
  read(11,*) CFL
  read(11,fmt='(/)')
  read(11,*) Uref
  read(11,fmt='(/)')
  read(11,*) poissmet
  read(11,fmt='(/)')
  read(11,*) convmet
  read(11,fmt='(/)')
  read(11,*) limitertype
  read(11,fmt='(/)')
  read(11,*) limparam
  read(11,fmt='(/)')
  read(11,*) wallmodeltype
  read(11,fmt='(/)')
  read(11,*) sgstype
  read(11,fmt='(/)')
  read(11,*) masssourc
  read(11,fmt='(/)')
  read(11,*) steady
  read(11,fmt='(/)')
  read(11,*) tasktype
  write(*,*) "tasktype=",tasktype
  read(11,fmt='(/)')
  read(11,*) initcondsfromfile
  read(11,fmt='(/)')
  read(11,*) timeavg1
  read(11,fmt='(/)')
  read(11,*) timeavg2
  read(11,fmt='(/)')
  read(11,*) Re
  write(*,*) "Re=",Re
  read(11,fmt='(/)')
  read(11,*) coriolisparam
  write(*,*) "coriolisparam=",coriolisparam
  read(11,fmt='(/)')
  read(11,*) prgradientx
  write(*,*) "prgradientx=",prgradientx
  read(11,fmt='(/)')
  read(11,*) prgradienty
  write(*,*) "prgradienty=",prgradienty
  read(11,fmt='(/)')
  read(11,*) starttime
  write(*,*) "starttime=",starttime
  read(11,fmt='(/)')
  read(11,*) endtime
  write(*,*) "endtime=",endtime
  read(11,fmt='(/)')
  read(11,*) maxiter
  write(*,*) "maxiter=",maxiter
  read(11,fmt='(/)')
  read(11,*) eps
  write(*,*) "eps=",eps
  read(11,fmt='(/)')
  read(11,*) maxCNiter
  write(*,*) "maxCNiter=",maxCNiter
  read(11,fmt='(/)')
  read(11,*) epsCN
  write(*,*) "epsCN=",epsCN
  read(11,fmt='(/)')
  read(11,*) maxPOISSONiter
  write(*,*) "maxPOISSONiter=",maxPOISSONiter
  read(11,fmt='(/)')
  read(11,*) epsPOISSON
  write(*,*) "epsPOISSON=",epsPOISSON
  read(11,fmt='(/)')
  read(11,*) debugparam
  write(*,*) "debug parameter=",debugparam
  close(11)


  open(11,file="grid.conf",status="old",action="read")
  read(11,fmt='(/)')
  read(11,*) xgridfromfile
  read(11,fmt='(/)')
  read(11,*) ygridfromfile
  read(11,fmt='(/)')
  read(11,*) zgridfromfile
  read(11,fmt='(/)')
  read(11,*) x0
  write(*,*) "x0=",x0
  read(11,fmt='(/)')
  read(11,*) y0
  write(*,*) "y0=",y0
  read(11,fmt='(/)')
  read(11,*) z0
  write(*,*) "z0=",z0
  read(11,fmt='(/)')
  read(11,*) lx
  write(*,*) "lx=",lx
  read(11,fmt='(/)')
  read(11,*) ly
  write(*,*) "ly=",ly
  read(11,fmt='(/)')
  read(11,*) lz
  write(*,*) "lz=",lz
  read(11,fmt='(/)')
  read(11,*) Prnx
  write(*,*) "nx=",Prnx
  read(11,fmt='(/)')
  read(11,*) Prny
  write(*,*) "ny=",Prny
  read(11,fmt='(/)')
  read(11,*) Prnz
  write(*,*) "ny=",Prnz
  close(11)

  open(11,file="boundconds.conf",status="old",action="read")
  read(11,fmt='(/)')
  read(11,*) Btype(We)
  read(11,fmt='(/)')
  read(11,*) Btype(Ea)
  read(11,fmt='(/)')
  read(11,*) Btype(So)
  read(11,fmt='(/)')
  read(11,*) Btype(No)
  read(11,fmt='(/)')
  read(11,*) Btype(Bo)
  read(11,fmt='(/)')
  read(11,*) Btype(To)
  read(11,fmt='(/)')
  read(11,*) sideU(1,So)
  read(11,fmt='(/)')
  read(11,*) sideU(2,So)
  read(11,fmt='(/)')
  read(11,*) sideU(3,So)
  read(11,fmt='(/)')
  read(11,*) sideU(1,No)
  read(11,fmt='(/)')
  read(11,*) sideU(2,No)
  read(11,fmt='(/)')
  read(11,*) sideU(3,No)
  read(11,fmt='(/)')
  read(11,*) sideU(1,Bo)
  read(11,fmt='(/)')
  read(11,*) sideU(2,Bo)
  read(11,fmt='(/)')
  read(11,*) sideU(3,Bo)
  read(11,fmt='(/)')
  read(11,*) sideU(1,To)
  read(11,fmt='(/)')
  read(11,*) sideU(2,To)
  read(11,fmt='(/)')
  read(11,*) sideU(3,To)
  read(11,fmt='(/)')
  read(11,*) z0W
  read(11,fmt='(/)')
  read(11,*) z0E
  read(11,fmt='(/)')
  read(11,*) z0S
  read(11,fmt='(/)')
  read(11,*) z0N
  read(11,fmt='(/)')
  read(11,*) z0B
  read(11,fmt='(/)')
  read(11,*) z0T
  close(11)


  open(11,file="thermal.conf",status="old",action="read")
  read(11,fmt='(/)')
  read(11,*) buoyancy
  read(11,fmt='(/)')
  read(11,*) Prandtl
  read(11,fmt='(/)')
  read(11,*) grav_acc
  read(11,fmt='(/)')
  read(11,*) temperature_ref
  read(11,fmt='(/)')
  read(11,*) TBtype(We)
  read(11,fmt='(/)')
  read(11,*) TBtype(Ea)
  read(11,fmt='(/)')
  read(11,*) TBtype(So)
  read(11,fmt='(/)')
  read(11,*) TBtype(No)
  read(11,fmt='(/)')
  read(11,*) TBtype(Bo)
  read(11,fmt='(/)')
  read(11,*) TBtype(To)
  read(11,fmt='(/)')
  read(11,*) sideTemp(We)
  read(11,fmt='(/)')
  read(11,*) sideTemp(Ea)
  read(11,fmt='(/)')
  read(11,*) sideTemp(So)
  read(11,fmt='(/)')
  read(11,*) sideTemp(No)
  read(11,fmt='(/)')
  read(11,*) sideTemp(Bo)
  read(11,fmt='(/)')
  read(11,*) sideTemp(To)
  close(11)

  open(11,file="inlet.conf",status="old",action="read")
  read(11,fmt='(/)')
  read(11,*) inlettype
  read(11,fmt='(/)')
  read(11,*) profiletype
  read(11,fmt='(/)')
  read(11,*) SHEARG
  write(*,*) "G=",SHEARG
  read(11,fmt='(/)')
  read(11,*) Uinlet
  write(*,*) "Uinlet=",Uinlet
  read(11,fmt='(/)')
  read(11,*) ustarsurfin  !-<u'w'>
  read(11,fmt='(/)')
  read(11,*) stressgradin !in relative part per 1m
  read(11,fmt='(/)')
  read(11,*) z0inlet
  read(11,fmt='(/)')
  read(11,*) powerexpin
  read(11,fmt='(/)')
  read(11,*) zrefin
  read(11,fmt='(/)')
  read(11,*) Urefin
  read(11,fmt='(/)')
  read(11,*)  relativestress(1,1)
  read(11,fmt='(/)')
  read(11,*)  relativestress(2,2)
  read(11,fmt='(/)')
  read(11,*) relativestress(3,3)
  read(11,fmt='(/)')
  read(11,*) relativestress(1,2)
  read(11,fmt='(/)')
  read(11,*) relativestress(1,3)
  read(11,fmt='(/)')
  read(11,*) relativestress(2,3)
  read(11,fmt='(/)')
  read(11,*) TLag
  read(11,fmt='(/)')
  read(11,*) Lturby
  read(11,fmt='(/)')
  read(11,*) Lturbz
  close(11)

  relativestress(2,1)=relativestress(1,2)
  relativestress(3,1)=relativestress(1,3)
  relativestress(3,2)=relativestress(2,3)

  open(11,file="scalars.conf",status="old",action="read")
  read(11,fmt='(/)')
  read(11,*) computescalars
  read(11,fmt='(/)')
  read(11,*) computedeposition
  read(11,fmt='(/)')
  read(11,*) computegravsettling
  read(11,fmt='(/)')
  read(11,*) partdistrib
  read(11,fmt='(/)')
  read(11,*) totalscalsource
  read(11,fmt='(/)')
  read(11,*) scalsourcetype

  read(11,fmt='(/)')
  read(11,*) ScalBtype(We)
  read(11,fmt='(/)')
  read(11,*) ScalBtype(Ea)
  read(11,fmt='(/)')
  read(11,*) ScalBtype(So)
  read(11,fmt='(/)')
  read(11,*) ScalBtype(No)
  read(11,fmt='(/)')
  read(11,*) ScalBtype(Bo)
  read(11,fmt='(/)')
  read(11,*) ScalBtype(To)
  read(11,fmt='(/)')
  read(11,*) sideScal(We)
  read(11,fmt='(/)')
  read(11,*) sideScal(Ea)
  read(11,fmt='(/)')
  read(11,*) sideScal(So)
  read(11,fmt='(/)')
  read(11,*) sideScal(No)
  read(11,fmt='(/)')
  read(11,*) sideScal(Bo)
  read(11,fmt='(/)')
  read(11,*) sideScal(To)

  if (partdistrib>0) then
     allocate(partdiam(partdistrib),partrho(partdistrib),percdistrib(partdistrib))
     allocate(scalsrcx(partdistrib),scalsrcy(partdistrib),scalsrcz(partdistrib))
     allocate(scalsrci(partdistrib),scalsrcj(partdistrib),scalsrck(partdistrib))
     do i=1,partdistrib
      read(11,fmt='(/)')
      read(11,*) partdiam(i)
      read(11,fmt='(/)')
      read(11,*) partrho(i)
      read(11,fmt='(/)')
      read(11,*) percdistrib(i)
      read(11,fmt='(/)')
      read(11,*) scalsrcx(i)
      read(11,fmt='(/)')
      read(11,*) scalsrcy(i)
      read(11,fmt='(/)')
      read(11,*) scalsrcz(i)
     enddo

  else
     allocate(partdiam(computescalars),partrho(computescalars),percdistrib(computescalars))
     allocate(scalsrcx(computescalars),scalsrcy(computescalars),scalsrcz(computescalars))
     allocate(scalsrci(computescalars),scalsrcj(computescalars),scalsrck(computescalars))
     do i=1,computescalars
      read(11,fmt='(/)')
      read(11,*) partdiam(i)
      read(11,fmt='(/)')
      read(11,*) partrho(i)
      read(11,fmt='(/)')
      read(11,*) percdistrib(i)
      read(11,fmt='(/)')
      read(11,*) scalsrcx(i)
      read(11,fmt='(/)')
      read(11,*) scalsrcy(i)
      read(11,fmt='(/)')
      read(11,*) scalsrcz(i)
     enddo
  endif
  close(11)

  if (poissmet==3.or.poissmet==4.or.poissmet==5) then
    open(11,file="mgopts.conf",status="old",action="read")
    read(11,fmt='(/)')
    read(11,*) lmg
    read(11,fmt='(/)')
    read(11,*) minmglevel
    read(11,fmt='(/)')
    read(11,*) bnx
    read(11,fmt='(/)')
    read(11,*) bny
    read(11,fmt='(/)')
    read(11,*) bnz
    read(11,fmt='(/)')
    read(11,*) mgncgc
    read(11,fmt='(/)')
    read(11,*) mgnpre
    read(11,fmt='(/)')
    read(11,*) mgnpost
    read(11,fmt='(/)')
    read(11,*) mgmaxinnerGSiter
    read(11,fmt='(/)')
    read(11,*) mgepsinnerGS
    read(11,fmt='(/)',iostat=io)
    read(11,*,iostat=io) minGPUlevel
    close(11)

    if (poissmet==3.or.poissmet==4) then
      if (Prny==1) then
       call SetMGParams2d(llmg=lmg,lminmglevel=minmglevel,lbnx=bnx,lbnz=bnz,&
                          lmgncgc=mgncgc,lmgnpre=mgnpre,lmgnpost=mgnpost,&
                          lmgmaxinnerGSiter=mgmaxinnerGSiter,lmgepsinnerGS=mgepsinnerGS)
      else
       call SetMGParams(llmg=lmg,lminmglevel=minmglevel,lminGPUlevel=minGPUlevel,&
                          lbnx=bnx,lbny=bny,lbnz=bnz,&
                          lmgncgc=mgncgc,lmgnpre=mgnpre,lmgnpost=mgnpost,&
                          lmgmaxinnerGSiter=mgmaxinnerGSiter,lmgepsinnerGS=mgepsinnerGS)
      endif
    elseif (poissmet==5) then
     MUDbnx=bnx
     MUDbny=bny
     MUDbnz=bnz
     MUDlmg=MUDlmg
    endif
  endif

  open(11,file="frames.conf",status="old",action="read")
  read(11,fmt='(/)')
  read(11,*) frames
  read(11,fmt='(/)')
  read(11,*) timefram1
  read(11,fmt='(/)')
  read(11,*) timefram2
  read(11,fmt='(/)')

  read(11,*) store%frame_U
  read(11,fmt='(/)')
  read(11,*) store%frame_vort
  read(11,fmt='(/)')
  read(11,*) store%frame_Pr
  read(11,fmt='(/)')
  read(11,*) store%frame_lambda2
  read(11,fmt='(/)')
  read(11,*) store%frame_scalars
  read(11,fmt='(/)')
  read(11,*) store%frame_sumscalars
  read(11,fmt='(/)')
  read(11,*) store%frame_T


  read(11,fmt='(/)')
  read(11,*) framedimension
  read(11,fmt='(/)')
  read(11,*) slicedir
  read(11,fmt='(/)')
  read(11,*) slicex
  close(11)

  open(11,file="output.conf",status="old",action="read",iostat=io)
  if (io==0) then
    read(11,fmt='(/)')
    read(11,*) display%delta
    read(11,fmt='(/)')
    read(11,*) display%ustar
    read(11,fmt='(/)')
    read(11,*) display%tstar

    read(11,fmt='(/)')
    read(11,*) store%U
    read(11,fmt='(/)')
    read(11,*) store%U_interp
    read(11,fmt='(/)')
    read(11,*) store%V
    read(11,fmt='(/)')
    read(11,*) store%V_interp
    read(11,fmt='(/)')
    read(11,*) store%W
    read(11,fmt='(/)')
    read(11,*) store%W_interp

    read(11,fmt='(/)')
    read(11,*) store%out

    read(11,fmt='(/)')
    read(11,*) store%out_U
    read(11,fmt='(/)')
    read(11,*) store%out_vort
    read(11,fmt='(/)')
    read(11,*) store%out_Pr
    read(11,fmt='(/)')
    read(11,*) store%out_Prtype
    read(11,fmt='(/)')
    read(11,*) store%out_lambda2
    read(11,fmt='(/)')
    read(11,*) store%out_T
    read(11,fmt='(/)')
    read(11,*) store%out_div
    read(11,fmt='(/)')
    read(11,*) store%out_visc

    read(11,fmt='(/)')
    read(11,*) store%avg

    read(11,fmt='(/)')
    read(11,*) store%avg_U
    read(11,fmt='(/)')
    read(11,*) store%avg_vort
    read(11,fmt='(/)')
    read(11,*) store%avg_Pr
    read(11,fmt='(/)')
    read(11,*) store%avg_Prtype
    read(11,fmt='(/)')
    read(11,*) store%avg_T

    read(11,fmt='(/)')
    read(11,*) store%scalars
    read(11,fmt='(/)')
    read(11,*) store%scalarsavg

    read(11,fmt='(/)')
    read(11,*) store%deposition

    read(11,fmt='(/)')
    read(11,*) store%deltime
    read(11,fmt='(/)')
    read(11,*) store%tke
    read(11,fmt='(/)')
    read(11,*) store%dissip
    read(11,fmt='(/)')
    read(11,*) store%scalsumtime
    read(11,fmt='(/)')
    read(11,*) store%scaltotsumtime
    read(11,fmt='(/)')
    read(11,*) store%ustar
    read(11,fmt='(/)')
    read(11,*) store%tstar
    read(11,fmt='(/)')
    read(11,*) store%blprofiles

    read(11,fmt='(/)')
    read(11,*) NumProbes
    allocate(probes(Numprobes))

    do i=1,NumProbes
      read(11,fmt='(/)')
      read(11,*) probes(i)%x
      read(11,fmt='(/)')
      read(11,*) probes(i)%y
      read(11,fmt='(/)')
      read(11,*) probes(i)%z
    enddo
    close(11)
  endif


  write(*,*) "computescalars",computescalars
  write(*,*) "partdiam",partdiam




  windangle=0._KND

  projectiontype=1

  fullstress=0

  if (CFL<=0)  CFL=0.5


  write(*,*) "Boundaries:"
  write(*,'(a2)',advance='no') "W "
  select case (Btype(We))
   case (NOSLIP)
    write(*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect
  write(*,'(a2)',advance='no') "E "
  select case (Btype(Ea))
   case (NOSLIP)
    write(*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect
  write(*,'(a2)',advance='no') "S "
  select case (Btype(So))
   case (NOSLIP)
    write(*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect
  write(*,'(a2)',advance='no') "N "
  select case (Btype(No))
   case (NOSLIP)
    write(*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect
  write(*,'(a2)',advance='no') "B "
  select case (Btype(Bo))
   case (NOSLIP)
    write(*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect
  write(*,'(a2)',advance='no') "T "
  select case (Btype(To))
   case (NOSLIP)
    write(*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect





  dxmin=lx/(Prnx)
  dymin=ly/(Prny)
  dzmin=lz/(Prnz)

  if (Prnx==1) then
   dxmin=sqrt(dymin*dzmin)
  elseif (Prny==1) then
   dymin=sqrt(dxmin*dzmin)
  elseif (Prnz==1) then
   dzmin=sqrt(dxmin*dymin)
  endif
  write(*,*) "dxmin ",dxmin
  write(*,*) "dymin ",dymin
  write(*,*) "dzmin ",dzmin

  write(*,*) "lx:",lx
  write(*,*) "ly:",ly
  write(*,*) "lz:",lz


  nt=maxiter

  if (Btype(Ea)==PERIODIC) then
                         Unx=Prnx
                        else
                         Unx=Prnx-1
  endif
  Uny=Prny
  Unz=Prnz

  Vnx=Prnx
  if (Btype(No)==PERIODIC) then
                         Vny=Prny
                        else
                         Vny=Prny-1
  endif
  Vnz=Prnz

  Wnx=Prnx
  Wny=Prny
  if (Btype(To)==PERIODIC) then
                         Wnz=Prnz
                        else
                         Wnz=Prnz-1
  endif

  if (Btype(We)==TURBULENTINLET) inlettype=TURBULENTINLET
  if (Btype(We)==INLETFROMFILE) inlettype=INLETFROMFILE


  if (Abs(Uinlet)>0) then
   dt=Abs(dxmin/Uinlet)
  else
   dt=dxmin
  endif
  deb=1

  if ((timeavg1>=0).and.(timeavg2>=timeavg1)) then
    averaging=1
  else
    averaging=0
  endif

  if (.not.xgridfromfile.and..not.ygridfromfile.and..not.zgridfromfile) then
   gridtype=UNIFORMGRID
   write(*,*) "Uniform grid"
  else
   gridtype=GENERALGRID
   write(*,*) "General grid"
  endif


  write(*,*) "set"
 end subroutine ReadParams

 subroutine ReadIC(U,V,W,Pr)
 real(KND) U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
 integer i,j,k
   open(11,file="in.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read(11,*)
   enddo
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
      read(11,*) Pr(i,j,k)
     enddo
    enddo
   enddo
   if (buoyancy==1) then
    do i=1,3
     read(11,*)
    enddo
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       read(11,*) temperature(i,j,k)
      enddo
     enddo
    enddo
   endif
   close(11)

   open(11,file="Uin.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read(11,*)
   enddo
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
      read(11,*) U(i,j,k)
     enddo
    enddo
   enddo
   close(11)

   open(11,file="Vin.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read(11,*)
   enddo
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
      read(11,*) V(i,j,k)
     enddo
    enddo
   enddo
   close(11)

   open(11,file="Win.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read(11,*)
   enddo
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
      read(11,*) W(i,j,k)
     enddo
    enddo
   enddo
 endsubroutine ReadIC



 subroutine INITCONDS(U,V,W,Pr)
 real(KND) U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
 integer i,j,k
 real(KND) p,x,y,z,x1,x2,y1,y2,z1,z2

 call init_random_seed

 Pr(1:Prnx,1:Prny,1:Prnz)=0

 U=100000
 V=100000
 W=100000

 V(1:Vnx,1:Vny,1:Vnz)=0
 W(1:Wnx,1:Wny,1:Wnz)=0

 if (initcondsfromfile==1) then
  call ReadIC(U,V,W,Pr)

  if (Re>0) then
   Visc=1._KND/Re
  else
   Visc=0
  endif
  if (buoyancy==1.or.computescalars>0) TDiff=1._KND/(Re*Prandtl)

  call BoundU(1,U)
  call BoundU(2,V)
  call BoundU(3,W)
  call Bound_Pr(Pr)

 else

  if (tasktype==2) then
   U(1:Unx,1:Uny,1:Unz)=0
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
      !if (Utype(i,j,k)==0) then
            !call RANDOM_NUMBER(p)
            x=xU(i)
            y=yPr(j)
            z=zPr(k)
            U(i,j,k)=-2*pi*y!*(1+0.1*(p-0.5))!-Uinlet*cos(z)*sin(x)*cos(y)!
        !0.5_KND*(p-0.5_KND)z
       !else
       !  V(i,j,k)=0
      !endif
     enddo
    enddo
   enddo
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
            !call RANDOM_NUMBER(p)
            x=xPr(i)
            y=yV(j)
            z=zPr(k)
            V(i,j,k)=2*pi*x!*(1+0.1*(p-0.5))!0
     enddo
    enddo
   enddo
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
            x=xPr(i)
            y=yPr(j)
            z=zW(k)
            W(i,j,k)=0
     enddo
    enddo
   enddo
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
            x=xPr(i)
            y=yPr(j)
            z=zPr(k)
            Pr(i,j,k)=0!(Uinlet/16._KND)*((2+cos(2*z))*(cos(2*(x))+cos(2*(y)))-2)
     enddo
    enddo
   enddo

  elseif (tasktype==3) then
   U(1:Unx,1:Uny,1:Unz)=0
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
      !if (Utype(i,j,k)==0) then
            !call RANDOM_NUMBER(p)
            x=xU(i)
            y=yPr(j)
            z=zPr(k)
            U(i,j,k)=Uinlet*sin(x)*cos(z)*cos(-y)!*(1+0.1*(p-0.5))!-Uinlet*cos(z)*sin(x)*cos(y)!
        !0.5_KND*(p-0.5_KND)z
       !else
       !  V(i,j,k)=0
      !endif
     enddo
    enddo
   enddo
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
            !call RANDOM_NUMBER(p)
            x=xPr(i)
            y=yV(j)
            z=zPr(k)
            V(i,j,k)=0!*(1+0.1*(p-0.5))!0
     enddo
    enddo
   enddo
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
            x=xPr(i)
            y=yPr(j)
            z=zW(k)
            W(i,j,k)=-Uinlet*cos(x)*sin(z)*cos(-y)!Uinlet*sin(z)*cos(x)*cos(y)!
     enddo
    enddo
   enddo
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
            x=xPr(i)
            y=yPr(j)
            z=zPr(k)
            Pr(i,j,k)=(Uinlet/16._KND)*((2+cos(2*z))*(cos(2*(-y))+cos(2*(x)))-2)!(Uinlet/16._KND)*((2+cos(2*z))*(cos(2*(x))+cos(2*(y)))-2)
     enddo
    enddo
   enddo

  elseif (tasktype==5) then!temporal mixing layer
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
            y=yPr(j)
            if ((abs(y-yPr(Uny/2)))<1.) then
             call RANDOM_NUMBER(p)
            else
             p=0
            endif
             x1=xPr(i)
             x2=xPr(i+1)
             y1=yV(j-1)
             y2=yV(j)
             z1=zW(k-1)
             z2=zW(k)
             p=0
               !(20/(pi*pi))*(cos((pi/10)*x1)-cos((pi/10)*x2))*((pi/2)*(z2-z1)-0.8*(cos((pi/2)*z2)-cos((pi/2)*z1)))
              !sin(pi*xU(i)/10)*(1+0.8*sin(pi*zPr(k)/2))
             x2=cosh(y2+0.4_KND*p-yPr(Uny/2))
             x1=cosh(y1+0.4_KND*p-yPr(Uny/2))
             if (abs(x2-x1)>0.000001_KND) then
              U(i,j,k)=(log(x2)-log(x1))/(y2-y1)
             else
              U(i,j,k)=0
             endif
                !Uinlet*tanh(y+0.1*sin(pi*(xU(i)/10)*sin(pi*zPr(k)/10)-yPr(Uny/2)) by mean of integrals
    enddo
   enddo
  enddo
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
            y=yV(j)
            !if ((abs(y-yPr(Uny/2)))<1.) then
            ! call RANDOM_NUMBER(p)
            !else
             p=sin(2*pi*xPr(i)/10)*(1+0.3*sin(3*pi*zPr(k)/10))*exp(-(yV(j)-yPr(Uny/2))*(yV(j)-yPr(Uny/2)))
            !endif

             V(i,j,k)=0.1*p!Uinlet*tanh(y-yPr(Uny/2))+Uinlet*0.1_KND*(p-0.5_KND)
     enddo
    enddo
   enddo
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
            y=yPr(j)
            !if ((abs(y-yPr(Uny/2)))<1.) then
            ! call RANDOM_NUMBER(p)
            !else
             p=cos(2*pi*xPr(i)/10)*(1+0.3*cos(3*pi*zW(k)/10))*exp(-(yPr(j)-yPr(Uny/2))*(yPr(j)-yPr(Uny/2)))
            !endif

             W(i,j,k)=0.1*p!Uinlet*tanh(y-yPr(Uny/2))+Uinlet*0.1_KND*(p-0.5_KND)
     enddo
    enddo
   enddo
!   elseif (tasktype==8) then
!    do k=1,Unz
!     do j=1,Uny
!      do i=1,Unx
!             call RANDOM_NUMBER(p)
!       U(i,j,k)=-prgradienty/(coriolisparam)*(1+0.1_KND*(p-0.5_KND))
!      enddo
!     enddo
!    enddo
!    do k=1,Vnz
!     do j=1,Vny
!      do i=1,Vnx
!               call RANDOM_NUMBER(p)
!      V(i,j,k)=prgradientx/(coriolisparam)*(1+0.1_KND*(p-0.5_KND))
!      enddo
!     enddo
!    enddo
!    do k=1,Wnz
!     do j=1,Wny
!      do i=1,Wnx
!       W(i,j,k)=0
!      enddo
!     enddo
!    enddo

  else
   U(1:Unx,1:Uny,1:Unz)=Uin(1,1)
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
      if (Utype(i,j,k)==0) then
            call RANDOM_NUMBER(p)
            U(i,j,k)=Uin(j,k)!+(Sqrt((Uin(j,k))**2+(Vin(j,k))**2))*0.1_KND*(p-0.5_KND)!sin(2.*pi*xU(i)+1)*cos(2.*pi*yPr(j)-2)!Uin(j,k)!*(1+0.03_KND*(p-0.5_KND))
       else
         U(i,j,k)=0
      endif
     enddo
    enddo
   enddo
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
      if (Vtype(i,j,k)==0) then
            call RANDOM_NUMBER(p)
            V(i,j,k)=Vin(j,k)!+(Sqrt((Uin(j,k))**2+(Vin(j,k))**2))*0.1_KND*(p-0.5_KND)!-cos(2.*pi*xPr(i)+1)*sin(2.*pi*yV(j)-2)!Uinlet*(0.3_KND*(p-0.5_KND))
       else
         V(i,j,k)=0
      endif
     enddo
    enddo
   enddo
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
      if (Wtype(i,j,k)==0) then
            call RANDOM_NUMBER(p)
            W(i,j,k)=Win(j,k)!Uinlet*(0.00001_KND*(p-0.5_KND))
       else
         W(i,j,k)=0
      endif
     enddo
    enddo
   enddo
  endif



  do i=1,computescalars
   SCALAR(1:Prnx,1:Prny,1:Prnz,i)=0
  enddo

  if (buoyancy==1.and.tasktype==2) then
   do k=0,Prnz+1
    do j=0,Prny+1
     do i=0,Prnx+1
      x=xPr(i)
      y=yPr(j)
      z=zPr(k)
      if ((x)**2+(y-0.5)**2<0.2_KND**2) then
       temperature(i,j,k)=cos(sqrt(x**2+(y-0.5)**2)*pi/2/0.2_KND)**2
      else
       temperature(i,j,k)=0
      endif
     enddo
    enddo
   enddo
  elseif (buoyancy==1.and.tasktype==3) then
   do k=0,Prnz+1
    do j=0,Prny+1
     do i=0,Prnx+1
      x=xPr(i)
      y=yPr(j)
      z=zPr(k)
      temperature(i,j,k)=temperature_ref+(temperature_ref/100._KND)*((2+cos(2*z))*(cos(2*(-y))+cos(2*(x)))-2)
     enddo
    enddo
   enddo
  elseif (buoyancy==1) then
   freetempgradient=0.02!0.0035 !K/m
   !inversionTjump=2 !K
   do k=-2,Prnz+3
    do j=-2,Prny+3
!         if (zPr(k)>(zW(Wnz+1)+zW(0))/4) then
        Tempin(j,k)=(zPr(k)-zW(0))*freetempgradient+temperature_ref!(zPr(k)-(zW(Wnz+1)+zW(0))/4)*freetempgradient+265!
!         else
!          Tempin(j,k)=265!temperature_ref
!         endif
!        Tempin(j,k)=temperature_ref
!
    enddo
   enddo
!    Tempin(j,1)=temperature_ref*1.01
   do k=0,Prnz+1
    do j=0,Prny+1
     do i=0,Prnx+1
      if (zPr(k)<=lz/8._KND) then
       call RANDOM_NUMBER(p)
      else
       p=0.5_KND
      endif
       temperature(i,j,k)=Tempin(j,k)!+0.2_KND*(p-0.5_KND)
!      if (yPr(j)<=yPr(Uny/2)) then
!        temperature(i,j,k)=0
!      else
!        temperature(i,j,k)=1
!      endif
     enddo
    enddo
   enddo
   Pr(1:Prnx,1:Prny,1)=0
   do k=2,Prnz
    do j=1,Prny
     do i=1,Prnx
      Pr(i,j,k)=Pr(i,j,k-1)+grav_acc*dzW(k-1)*((temperature(i,j,k-1)+temperature(i,j,k))/2._KND-temperature_ref)/temperature_ref
     enddo
    enddo
   enddo


  endif


  if (Re>0) then
   Visc=1._KND/Re
  else
   Visc=0
  endif

 if (Re>0) then
  if (buoyancy==1.or.computescalars>0) TDiff=1._KND/(Re*Prandtl)

  call BoundU(1,U)
  call BoundU(2,V)
  call BoundU(3,W)
  call Bound_Pr(Pr)
  call Pr_Correct(U,V,W,Pr,1._KND)


   if (sgstype==1) then
                     call Smag(U,V,W)
   elseif (sgstype==3) then
                     call VREMAN(U,V,W)
   else
    if (Re>0) then
     Visc=1._KND/Re
    else
     Visc=0
    endif
   endif

   call Bound_Visc(Visc)

  if (buoyancy>0) then
   forall(k=1:Prnz,j=1:Prny,i=1:Prnx)
     TDiff(i,j,k)=1.35*(Visc(i,j,k)-1._KND/Re)+(1._KND/(Re*Prt(i,j,k,U,V,temperature)))
    endforall
     call Bound_Visc(TDiff)
  endif
 endif


  call InitTempFL
   if (wallmodeltype>0) then
                   call ComputeViscsWM(U,V,W,Pr)
   endif
   call Bound_Visc(Visc)
  write(*,*) "set"
 endif
 endsubroutine INITCONDS






  subroutine READBOUNDS
  real(KND),allocatable:: xU2(:),yV2(:),zW2(:)
  integer i,j,k,nx,ny,nz,nxup,nxdown,nyup,nydown,nzup,nzdown,io
  real(KND) P
  type(WMPoint):: WMP


  nx=Prnx-1
  ny=Prny-1
  nz=Prnz-1

  if (xgridfromfile) then
   open(11,file="xgrid.txt")
   j=-1
   do
    read (11,*,iostat=io) P
    if (io==0) then
      j=j+1
    else
      exit
    endif
   enddo
   nx=j
   Prnx=nx
   Vnx=Prnx
   Wnx=Prnx
   if (Btype(Ea)==PERIODIC) then
                         Unx=Prnx
                        else
                         Unx=Prnx-1
   endif
   close(11)
  endif

  if (ygridfromfile) then
   open(11,file="ygrid.txt")
   j=-1
   do
    read (11,*,iostat=io) P
    if (io==0) then
      j=j+1
    else
      exit
    endif
   enddo
   ny=j
   Prny=ny
   Uny=Prny
   Wny=Prny
   if (Btype(No)==PERIODIC) then
                         Vny=Prny
                        else
                         Vny=Prny-1
   endif
   close(11)
  endif

  if (zgridfromfile) then
   open(11,file="zgrid.txt")
   j=-1
   do
    read (11,*,iostat=io) P
    if (io==0) then
      j=j+1
      write(*,*) j
    else
      exit
    endif
   enddo
   nz=j
   Prnz=nz
   Unz=Prnz
   Vnz=Prnz
   if (Btype(To)==PERIODIC) then
                         Wnz=Prnz
                        else
                         Wnz=Prnz-1
   endif
   close(11)
  endif





  allocate(xU2(-3:nx+4))
  allocate(yV2(-3:ny+4))
  allocate(zW2(-3:nz+4))



  if (xgridfromfile) then
   open(11,file="xgrid.txt")
   do j=0,nx
    write(*,*) j
    read(11,*) xU2(j)
   enddo
   close(11)

   if (Btype(We)==PERIODIC) then
    do j=-1,-3,-1
     xU2(j)=xU2(0)-(xU2(nx)-xU2(nx+j))
    enddo
   else
    do j=-1,-3,-1
     xU2(j)=xU2(0)-(xU2(0-j)-xU2(0))
    enddo
   endif

   if (Btype(Ea)==PERIODIC) then
    do j=nx+1,nx+4
     xU2(j)=xU2(nx)+(xU2(j-nx)-xU2(0))
    enddo
   else
    do j=nx+1,nx+4
     xU2(j)=xU2(nx)+(xU2(nx)-xU2(nx-(j-nx)))
    enddo
   endif

   x0=xU2(0)
  else
   forall (i=-3:nx+4)
    xU2(i)=(i)*dxmin+x0
   endforall
  endif


  if (ygridfromfile) then
   open(11,file="ygrid.txt")
   do j=0,ny
    read(11,*) yV2(j)
   enddo
   close(11)

   if (Btype(So)==PERIODIC) then
    do j=-1,-3,-1
     yV2(j)=yV2(0)-(yV2(ny)-yV2(ny+j))
    enddo
   else
    do j=-1,-3,-1
     yV2(j)=yV2(0)-(yV2(0-j)-yV2(0))
    enddo
   endif

   if (Btype(No)==PERIODIC) then
    do j=ny+1,ny+4
     yV2(j)=yV2(ny)+(yV2(j-ny)-yV2(0))
    enddo
   else
    do j=ny+1,ny+4
     yV2(j)=yV2(ny)+(yV2(ny)-yV2(ny-(j-ny)))
    enddo
   endif

   y0=yV2(0)
  else
   forall (j=-3:ny+4)
     yV2(j)=(j)*dymin+y0
   endforall
  endif


  if (zgridfromfile) then
   open(11,file="zgrid.txt")
   do j=0,nz
    read(11,*) zW2(j)
   enddo
   close(11)

   if (Btype(Bo)==PERIODIC) then
    do j=-1,-3,-1
     zW2(j)=zW2(0)-(zW2(nz)-zW2(nz+j))
    enddo
   else
    do j=-1,-3,-1
     zW2(j)=zW2(0)-(zW2(0-j)-zW2(0))
    enddo
   endif

   if (Btype(To)==PERIODIC) then
    do j=nz+1,nz+4
     zW2(j)=zW2(nz)+(zW2(j-nz)-zW2(0))
    enddo
   else
    do j=nz+1,nz+4
     zW2(j)=zW2(nz)+(zW2(nz)-zW2(nz-(j-nz)))
    enddo
   endif

   z0=zW2(0)
  else
   forall (k=-3:nz+4)
    zW2(k)=(k)*dzmin+z0
   endforall
  endif


  nxup=nx+1
  nxdown=0
  nyup=ny+1
  nydown=0
  nzup=nz+1
  nzdown=0


  allocate(xU(-3:nx+4))
  allocate(yV(-3:ny+4))
  allocate(zW(-3:nz+4))
  allocate(dxU(-2:nx+3))
  allocate(dyV(-2:ny+3))
  allocate(dzW(-2:nz+3))
  allocate(xPr(-2:nx+4),dxPr(-2:nx+4))
  allocate(yPr(-2:ny+4),dyPr(-2:ny+4))
  allocate(zPr(-2:nz+4),dzPr(-2:nz+4))

  xU=xU2(nxdown-3:nxup+3)
  yV=yV2(nydown-3:nyup+3)
  zW=zW2(nzdown-3:nzup+3)

  forall (i=-2:nx+4)
   xPr(i)=(xU(i-1)+xU(i))/2._KND
   dxPr(i)=xU(i)-xU(i-1)
  endforall
  forall (j=-2:ny+4)
   yPr(j)=(yV(j-1)+yV(j))/2._KND
   dyPr(j)=yV(j)-yV(j-1)
  endforall
  forall (k=-2:nz+4)
   zPr(k)=(zW(k-1)+zW(k))/2._KND
   dzPr(k)=zW(k)-zW(k-1)
  endforall
  forall (i=-2:nx+3)
   dxU(i)=xPr(i+1)-xPr(i)
  endforall
  forall (j=-2:ny+3)
   dyV(j)=yPr(j+1)-yPr(j)
  endforall
  forall (k=-2:nz+3)
   dzW(k)=zPr(k+1)-zPr(k)
  endforall

  deallocate(xU2)
  deallocate(yV2)
  deallocate(zW2)


  allocate(Utype(-2:Unx+3,-2:Uny+3,-2:Unz+3))
  allocate(Vtype(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
  allocate(Wtype(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
  allocate(Prtype(0:Prnx+1,0:Prny+1,0:Prnz+1))
  Utype=0
  Vtype=0
  Wtype=0
  Prtype=0

!   write (*,*) Prnx
!   write (*,*) Prny
!   write (*,*) Prnz
!   write (*,*) Unx
!   write (*,*) Uny
!   write (*,*) Unz


  allocate(Uin(-2:Uny+3,-2:Unz+3),Vin(-2:Vny+3,-2:Vnz+3),Win(-2:Wny+3,-2:Wnz+3))
  if (buoyancy>0) allocate(Tempin(-2:Prny+3,-2:Prnz+3))

  select case (inlettype)
   case (NOINLET)
    Uin=0
    Vin=0
    Win=0
   case (SHEAR)
    call SHEARINLET(SHEARG)
   case (PARABOLIC)
    call PARINLET
   case (TURBULENTINLET)
    call GETTURBINLET
   case (INLETFROMFILE)
    call GETINLETFROMFILE(starttime)
   case default
    call CONSTINLET
  endselect


  if (buoyancy==1) then
     if (TBtype(Bo)==CONSTFLUX.or.TBtype(Bo)==DIRICHLET) then
       allocate(BsideTFLArr(-1:Prnx+2,-1:Prny+2))
       if (TBtype(Bo)==CONSTFLUX) then
        BsideTFLArr=sideTemp(Bo)
       else
        BsideTFLArr=0
       endif
       if (TBtype(Bo)==DIRICHLET) then
        allocate(BsideTArr(-1:Prnx+2,-1:Prny+2))
        BsideTArr=sideTemp(Bo)
       endif
      endif
  endif





    allocate(WMP%depscalar(computescalars))
    WMP%depscalar=0

    if (Btype(We)==NOSLIP) then
     do k=1,Prnz
      do j=1,Prny
       WMP%x=1
       WMP%y=j
       WMP%z=k
       WMP%distx=(xPr(1)-xU(0))
       WMP%disty=0
       WMP%distz=0
       WMP%ustar=1
       WMP%z0=z0W
       call AddWMPoint(WMP)
      enddo
     enddo
    endif

    if (Btype(Ea)==NOSLIP) then
     do k=1,Prnz
      do j=1,Prny
       WMP%x=Prnx
       WMP%y=j
       WMP%z=k
       WMP%distx=(xPr(Prnx)-xU(Unx+1))
       WMP%disty=0
       WMP%distz=0
       WMP%ustar=1
       WMP%z0=z0E
       call AddWMPoint(WMP)
      enddo
     enddo
    endif

    if (Btype(So)==NOSLIP) then
     do k=1,Prnz
      do i=1,Prnx
       WMP%x=i
       WMP%y=1
       WMP%z=k
       WMP%distx=0
       WMP%disty=(yPr(1)-yV(0))
       WMP%distz=0
       WMP%ustar=1
       WMP%z0=z0S
       call AddWMPoint(WMP)
      enddo
     enddo
    elseif (Btype(So)==DIRICHLET) then
     do k=1,Prnz
      do i=1,Prnx
       WMP%x=i
       WMP%y=1
       WMP%z=k
       WMP%distx=0
       WMP%disty=(yPr(1)-yV(0))
       WMP%distz=0
       WMP%ustar=1
       WMP%wallu=sideU(1,So)
       WMP%wallv=sideU(2,So)
       WMP%wallw=sideU(3,So)
       WMP%z0=z0S
       call AddWMPoint(WMP)
      enddo
     enddo
    endif

    if (Btype(No)==NOSLIP) then
     do k=1,Prnz
      do i=1,Prnx
       WMP%x=i
       WMP%y=Prny
       WMP%z=k
       WMP%distx=0
       WMP%disty=(yPr(Prny)-yV(Vny+1))
       WMP%distz=0
       WMP%ustar=1
       WMP%z0=z0N
       call AddWMPoint(WMP)
      enddo
     enddo
    elseif (Btype(No)==DIRICHLET) then
     do k=1,Prnz
      do i=1,Prnx
       WMP%x=i
       WMP%y=Prny
       WMP%z=k
       WMP%distx=0
       WMP%disty=(yPr(Prny)-yV(Vny+1))
       WMP%distz=0
       WMP%ustar=1
       WMP%wallu=sideU(1,No)
       WMP%wallv=sideU(2,No)
       WMP%wallw=sideU(3,No)
       WMP%z0=z0N
       call AddWMPoint(WMP)
      enddo
     enddo
    endif

    if (Btype(Bo)==NOSLIP) then
     do j=1,Prny
      do i=1,Prnx
       if (Prtype(i,j,0)==0) then
        WMP%x=i
        WMP%y=j
        WMP%z=1
        WMP%distx=0
        WMP%disty=0
        WMP%distz=(zPr(1)-zW(0))
        WMP%ustar=1
        WMP%z0=z0B
       if (TBtype(Bo)==CONSTFLUX) then
        WMP%tempfl=sideTemp(Bo)
       else
        WMP%temp=0
       endif
       if (TBtype(Bo)==DIRICHLET) then
        WMP%temp=sideTemp(Bo)
       endif
        call AddWMPoint(WMP)
       endif
      enddo
     enddo
    elseif (Btype(Bo)==DIRICHLET) then
     do j=1,Prny
      do i=1,Prnx
       WMP%x=i
       WMP%y=j
       WMP%z=1
       WMP%distx=0
       WMP%disty=0
       WMP%distz=(zPr(1)-zW(0))
       WMP%ustar=1
       WMP%wallu=sideU(1,Bo)
       WMP%wallv=sideU(2,Bo)
       WMP%wallw=sideU(3,Bo)
       WMP%z0=z0B
       if (TBtype(Bo)==CONSTFLUX) then
        WMP%tempfl=sideTemp(Bo)
       else
        WMP%temp=0
       endif
       if (TBtype(Bo)==DIRICHLET) then
        WMP%temp=sideTemp(Bo)
       endif

       call AddWMPoint(WMP)
      enddo
     enddo
    endif

    if (Btype(To)==NOSLIP) then
     do j=1,Prny
      do i=1,Prnx
       WMP%x=i
       WMP%y=j
       WMP%z=Prnz
       WMP%distx=0
       WMP%disty=0
       WMP%distz=(zPr(Prnz)-zW(Wnz+1))
       WMP%ustar=1
       WMP%z0=z0T
       call AddWMPoint(WMP)
      enddo
     enddo
    elseif (Btype(To)==DIRICHLET) then
     do j=1,Prny
      do i=1,Prnx
       WMP%x=i
       WMP%y=j
       WMP%z=Prnz
       WMP%distx=0
       WMP%disty=0
       WMP%distz=(zPr(Prnz)-zW(Wnz+1))
       WMP%ustar=1
       WMP%wallu=sideU(1,To)
       WMP%wallv=sideU(2,To)
       WMP%wallw=sideU(3,To)
       WMP%z0=z0T
       call AddWMPoint(WMP)
      enddo
     enddo
    endif


   if (computescalars>0.and.scalsourcetype==pointsource) then
        call Gridcoords(scalsrci(:),scalsrcj(:),scalsrck(:),scalsrcx(:),scalsrcy(:),scalsrcz(:))
   endif

   call InitSolidBodies
   call GetSolidBodiesBC

  write (*,*) "set"
 end subroutine READBOUNDS












  subroutine INIT_RANDOM_SEED()
  integer:: i, n, clock
  integer,dimension(:),allocatable:: seed

   call RANDOM_SEED(size = n)
   allocate(seed(n))

   call SYSTEM_CLOCK(COUNT=clock)

   seed=clock+37*(/(i-1,i=1,n)/)
   call RANDOM_SEED(PUT=seed)

   deallocate(seed)
  endsubroutine


endmodule INITIAL

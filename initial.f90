module INITIAL
use PARAMETERS
use MULTIGRID
use MULTIGRID2d
use POISSON
use BOUNDARIES
use SCALARS
use SMAGORINSKY
use WALLMODELS

implicit none

integer   lmg,minmglevel,bnx,bny,bnz,mgncgc,mgnpre,mgnpost,mgmaxinnerGSiter
real(KND) mgepsinnerGS


contains
 subroutine READPARAMS
 integer i
 
  open(11,file="main.conf",status="OLD",action="READ")
  read (11,fmt='(/)')
  read (11,*) tempmet
  read (11,fmt='(/)')
  read (11,*) CFL
  read (11,fmt='(/)')
  read (11,*) Uref
  read (11,fmt='(/)')
  read (11,*) poissmet
  read (11,fmt='(/)')
  read (11,*) convmet
  read (11,fmt='(/)')
  read (11,*) limitertype
  read (11,fmt='(/)')
  read (11,*) limparam
  read (11,fmt='(/)')
  read (11,*) wallmodeltype
  read (11,fmt='(/)')
  read (11,*) sgstype
  read (11,fmt='(/)')
  read (11,*) masssourc
  read (11,fmt='(/)')
  read (11,*) steady
  read (11,fmt='(/)')
  read (11,*) tasktype
  write(*,*) "tasktype=",tasktype
  read (11,fmt='(/)')
  read (11,*) initcondsfromfile
  read (11,fmt='(/)')
  read (11,*) timeavg1
  read (11,fmt='(/)')
  read (11,*) timeavg2
  read (11,fmt='(/)')
  read (11,*) Re
  write (*,*) "Re=",Re
  read (11,fmt='(/)')
  read (11,*) coriolisparam
  write (*,*) "coriolisparam=",coriolisparam
  read (11,fmt='(/)')
  read (11,*) prgradientx
  write (*,*) "prgradientx=",prgradientx
  read (11,fmt='(/)')
  read (11,*) prgradienty
  write (*,*) "prgradienty=",prgradienty
  read (11,fmt='(/)')
  read (11,*) starttime
  write (*,*) "starttime=",starttime
  read (11,fmt='(/)')
  read (11,*) endtime
  write (*,*) "endtime=",endtime
  read (11,fmt='(/)')
  read (11,*) maxiter
  write (*,*) "maxiter=",maxiter
  read (11,fmt='(/)')
  read (11,*) eps
  write (*,*) "eps=",eps
  read (11,fmt='(/)')
  read (11,*) maxCNiter
  write (*,*) "maxCNiter=",maxCNiter
  read (11,fmt='(/)')
  read (11,*) epsCN
  write (*,*) "epsCN=",epsCN
  read (11,fmt='(/)')
  read (11,*) maxPOISSONiter
  write (*,*) "maxPOISSONiter=",maxPOISSONiter
  read (11,fmt='(/)')
  read (11,*) epsPOISSON
  write (*,*) "epsPOISSON=",epsPOISSON
  read (11,fmt='(/)')
  read (11,*) debugparam
  write (*,*) "debug parameter=",debugparam
  read (11,fmt='(/)')
  read (11,*) probex
  read (11,fmt='(/)')
  read (11,*) probey
  read (11,fmt='(/)')
  read (11,*) probez
  close(11)


  open(11,file="grid.conf",status="OLD",action="READ")
  read (11,fmt='(/)')
  read (11,*) xgridfromfile
  read (11,fmt='(/)')
  read (11,*) ygridfromfile
  read (11,fmt='(/)')
  read (11,*) zgridfromfile
  read (11,fmt='(/)')
  read (11,*) x0
  write (*,*) "x0=",x0
  read (11,fmt='(/)')
  read (11,*) y0
  write (*,*) "y0=",y0
  read (11,fmt='(/)')
  read (11,*) z0
  write (*,*) "z0=",z0
  read (11,fmt='(/)')
  read (11,*) lx
  write (*,*) "lx=",lx
  read (11,fmt='(/)')
  read (11,*) ly
  write (*,*) "ly=",ly
  read (11,fmt='(/)')
  read (11,*) lz
  write (*,*) "lz=",lz
  read (11,fmt='(/)')
  read (11,*) Prnx
  write (*,*) "nx=",Prnx
  read (11,fmt='(/)')
  read (11,*) Prny
  write (*,*) "ny=",Prny
  read (11,fmt='(/)')
  read (11,*) Prnz
  write (*,*) "ny=",Prnz
  close(11)

  open(11,file="boundconds.conf",status="OLD",action="READ")
  read (11,fmt='(/)')
  read (11,*) BtypeW
  read (11,fmt='(/)')
  read (11,*) BtypeE
  read (11,fmt='(/)')
  read (11,*) BtypeS
  read (11,fmt='(/)')
  read (11,*) BtypeN
  read (11,fmt='(/)')
  read (11,*) BtypeB
  read (11,fmt='(/)')
  read (11,*) BtypeT
  read (11,fmt='(/)')
  read (11,*) SsideU
  read (11,fmt='(/)')
  read (11,*) SsideV
  read (11,fmt='(/)')
  read (11,*) SsideW
  read (11,fmt='(/)')
  read (11,*) NsideU
  read (11,fmt='(/)')
  read (11,*) NsideV
  read (11,fmt='(/)')
  read (11,*) NsideW
  read (11,fmt='(/)')
  read (11,*) BsideU
  read (11,fmt='(/)')
  read (11,*) BsideV
  read (11,fmt='(/)')
  read (11,*) BsideW
  read (11,fmt='(/)')
  read (11,*) TsideU
  read (11,fmt='(/)')
  read (11,*) TsideV
  read (11,fmt='(/)')
  read (11,*) TsideW
  read (11,fmt='(/)')
  read (11,*) z0W
  read (11,fmt='(/)')
  read (11,*) z0E
  read (11,fmt='(/)')
  read (11,*) z0S
  read (11,fmt='(/)')
  read (11,*) z0N
  read (11,fmt='(/)')
  read (11,*) z0B
  read (11,fmt='(/)')
  read (11,*) z0T
  close(11)

  
  open(11,file="thermal.conf",status="OLD",action="READ")
  read (11,fmt='(/)')
  read (11,*) buoyancy
  read (11,fmt='(/)')
  read (11,*) Prandtl
  read (11,fmt='(/)')
  read (11,*) grav_acc
  read (11,fmt='(/)')
  read (11,*) temperature_ref
  read (11,fmt='(/)')
  read (11,*) TBtypeW
  read (11,fmt='(/)')
  read (11,*) TBtypeE
  read (11,fmt='(/)')
  read (11,*) TBtypeS
  read (11,fmt='(/)')
  read (11,*) TBtypeN
  read (11,fmt='(/)')
  read (11,*) TBtypeB
  read (11,fmt='(/)')
  read (11,*) TBtypeT
  read (11,fmt='(/)')
  read (11,*) WsideTemp
  read (11,fmt='(/)')
  read (11,*) EsideTemp
  read (11,fmt='(/)')
  read (11,*) SsideTemp
  read (11,fmt='(/)')
  read (11,*) NsideTemp
  read (11,fmt='(/)')
  read (11,*) BsideTemp
  read (11,fmt='(/)')
  read (11,*) TsideTemp
  close(11)

  open(11,file="inlet.conf",status="OLD",action="READ")
  read (11,fmt='(/)')
  read (11,*) inlettype
  read (11,fmt='(/)')
  read (11,*) profiletype
  read (11,fmt='(/)')
  read (11,*) SHEARG
  write (*,*) "G=",SHEARG
  read (11,fmt='(/)')
  read (11,*) Uinlet
  write (*,*) "Uinlet=",Uinlet
  read (11,fmt='(/)')
  read (11,*) ustarsurfin  !-<u'w'>
  read (11,fmt='(/)')
  read (11,*) stressgradin !in relative part per 1m
  read (11,fmt='(/)')
  read (11,*) z0inlet
  read (11,fmt='(/)')
  read (11,*) powerexpin
  read (11,fmt='(/)')
  read (11,*) zrefin
  read (11,fmt='(/)')
  read (11,*) Urefin
  read (11,fmt='(/)')
  read (11,*)  relativestress(1,1)
  read (11,fmt='(/)')
  read (11,*)  relativestress(2,2)
  read (11,fmt='(/)')
  read (11,*) relativestress(3,3)
  read (11,fmt='(/)')
  read (11,*) relativestress(1,2)
  read (11,fmt='(/)')
  read (11,*) relativestress(1,3)
  read (11,fmt='(/)')
  read (11,*) relativestress(2,3)
  close(11)
 
  relativestress(2,1)=relativestress(1,2)
  relativestress(3,1)=relativestress(1,3)
  relativestress(3,2)=relativestress(2,3)

  open(11,file="scalars.conf",status="OLD",action="READ")
  read (11,fmt='(/)')
  read (11,*) computescalars
  read (11,fmt='(/)')
  read (11,*) computedeposition
  read (11,fmt='(/)')
  read (11,*) computegravsettling
  read (11,fmt='(/)')
  read (11,*) partdistrib
  read (11,fmt='(/)')
  read (11,*) totalscalsource
  read (11,fmt='(/)')
  read (11,*) pointscalsource

  read (11,fmt='(/)')
  read (11,*) ScalBtypeW
  read (11,fmt='(/)')
  read (11,*) ScalBtypeE
  read (11,fmt='(/)')
  read (11,*) ScalBtypeS
  read (11,fmt='(/)')
  read (11,*) ScalBtypeN
  read (11,fmt='(/)')
  read (11,*) ScalBtypeB
  read (11,fmt='(/)')
  read (11,*) ScalBtypeT
  read (11,fmt='(/)')
  read (11,*) WsideScal
  read (11,fmt='(/)')
  read (11,*) EsideScal
  read (11,fmt='(/)')
  read (11,*) SsideScal
  read (11,fmt='(/)')
  read (11,*) NsideScal
  read (11,fmt='(/)')
  read (11,*) BsideScal
  read (11,fmt='(/)')
  read (11,*) TsideScal
  
  if (partdistrib>0) then
     allocate(partdiam(partdistrib),partrho(partdistrib),percdistrib(partdistrib))
     allocate(scalsrcx(partdistrib),scalsrcy(partdistrib),scalsrcz(partdistrib))
     allocate(scalsrci(partdistrib),scalsrcj(partdistrib),scalsrck(partdistrib))
     do i=1,partdistrib
      read (11,fmt='(/)')
      read (11,*) partdiam(i)
      read (11,fmt='(/)')
      read (11,*) partrho(i)
      read (11,fmt='(/)')
      read (11,*) percdistrib(i)
      read (11,fmt='(/)')
      read (11,*) scalsrcx(i)
      read (11,fmt='(/)')
      read (11,*) scalsrcy(i)
      read (11,fmt='(/)')
      read (11,*) scalsrcz(i)
     enddo
   
  else
     allocate(partdiam(computescalars),partrho(computescalars),percdistrib(computescalars))
     allocate(scalsrcx(computescalars),scalsrcy(computescalars),scalsrcz(computescalars))
     allocate(scalsrci(computescalars),scalsrcj(computescalars),scalsrck(computescalars))
     do i=1,computescalars
      read (11,fmt='(/)')
      read (11,*) partdiam(i)
      read (11,fmt='(/)')
      read (11,*) partrho(i)
      read (11,fmt='(/)')
      read (11,*) percdistrib(i)
      read (11,fmt='(/)')
      read (11,*) scalsrcx(i)
      read (11,fmt='(/)')
      read (11,*) scalsrcy(i)
      read (11,fmt='(/)')
      read (11,*) scalsrcz(i)
     enddo
  endif
  close(11)

  if (poissmet==3.or.poissmet==4.or.poissmet==5) then
    open(11,file="mgopts.conf",status="OLD",action="READ")
    read (11,fmt='(/)')
    read (11,*) lmg
    read (11,fmt='(/)')
    read (11,*) minmglevel
    read (11,fmt='(/)')
    read (11,*) bnx
    read (11,fmt='(/)')
    read (11,*) bny
    read (11,fmt='(/)')
    read (11,*) bnz
    read (11,fmt='(/)')
    read (11,*) mgncgc
    read (11,fmt='(/)')
    read (11,*) mgnpre
    read (11,fmt='(/)')
    read (11,*) mgnpost
    read (11,fmt='(/)')
    read (11,*) mgmaxinnerGSiter
    read (11,fmt='(/)')
    read (11,*) mgepsinnerGS
    close(11)

    if (poissmet==3.or.poissmet==4) then
      if (Prny==1) then
       call SetMGParams2d(llmg=lmg,lminmglevel=minmglevel,lbnx=bnx,lbnz=bnz,&
                          lmgncgc=mgncgc,lmgnpre=mgnpre,lmgnpost=mgnpost,&
                          lmgmaxinnerGSiter=mgmaxinnerGSiter,lmgepsinnerGS=mgepsinnerGS)
      else
       call SetMGParams(llmg=lmg,lminmglevel=minmglevel,lbnx=bnx,lbny=bny,lbnz=bnz,&
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

  open(11,file="frames.conf",status="OLD",action="READ")
  read (11,fmt='(/)')
  read (11,*) frames
  read (11,fmt='(/)')
  read (11,*) timefram1
  read (11,fmt='(/)')
  read (11,*) timefram2
  read (11,fmt='(/)')
  read (11,*) framedimension
  read (11,fmt='(/)')
  read (11,*) slicedir
  read (11,fmt='(/)')
  read (11,*) slicex
  close(11)



  write(*,*) "computescalars",computescalars
  write(*,*) "partdiam",partdiam




  windangle=0._KND
  projectiontype=1
  fullstress=0

  if (CFL<=0)  CFL=0.5


  write (*,*) "Boundaries:"
  write (*,*) "W"
  select case (BtypeW)
   case (NOSLIP)
    write (*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect
  write (*,*) "E"
  select case (BtypeE)
   case (NOSLIP)
    write (*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect
  write (*,*) "S"
  select case (BtypeS)
   case (NOSLIP)
    write (*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect
  write (*,*) "N"
  select case (BtypeN)
   case (NOSLIP)
    write (*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect
  write (*,*) "B"
  select case (BtypeB)
   case (NOSLIP)
    write (*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect
  write (*,*) "T"
  select case (BtypeT)
   case (NOSLIP)
    write (*,*) "noslip"
   case (FREESLIP)
    write(*,*) "freeslip"
   case (PERIODIC)
    write(*,*) "periodic"
   case (DIRICHLET)
    write(*,*) "dirichlet"
   case (NEUMANN)
    write(*,*) "neumann"
  endselect



  if (tasktype==3) then
   x0=-pi
   y0=-pi
   z0=-pi
  endif
  patternnx=1
  patternny=1
  patternny=1
  if (tasktype==6) then
   patternnx=4
   patternny=4
   patternnz=6
   boxnx=Prnx
   boxny=Prny
   boxnz=Prnz
   tilenx=boxnx*2
   tileny=(boxny*3)/2
   tilenz=boxnz
   Prnx=tilenx*patternnx
   Prny=tileny*patternny
   Prnz=tilenz*patternnz
   lx=2*patternnx
   ly=(3./2.)*patternny
   lz=patternnz
  endif
  
  dxmin=lx/(Prnx)
  write (*,*) "dxmin ",dxmin
  dymin=ly/(Prny)
  write (*,*) "dymin ",dymin
  dzmin=lz/(Prnz)
  write (*,*) "dzmin ",dzmin

  write(*,*) "lx:",lx
  write(*,*) "ly:",ly
  write(*,*) "lz:",lz
  
  
  nt=maxiter

  if (BtypeE==PERIODIC) then
                         Unx=Prnx
                        else
                         Unx=Prnx-1
  endif
  Uny=Prny
  Unz=Prnz

  Vnx=Prnx
  if (BtypeN==PERIODIC) then
                         Vny=Prny
                        else
                         Vny=Prny-1
  endif
  Vnz=Prnz

  Wnx=Prnx
  Wny=Prny
  if (BtypeT==PERIODIC) then
                         Wnz=Prnz
                        else
                         Wnz=Prnz-1
  endif

  if (BtypeW==TURBULENTINLET) inlettype=TURBULENTINLET
  if (BtypeW==INLETFROMFILE) inlettype=INLETFROMFILE
  

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
   write (*,*) "Uniform grid"
  else
   gridtype=GENERALGRID
   write (*,*) "General grid"
  endif

  if (probex<x0) then
   probex=x0
  elseif (probex>x0+lx) then
   probex=x0+lx
  endif

  if (probey<y0) then
   probey=y0
  elseif (probey>y0+ly) then
   probey=y0+ly
  endif

  if (probez<z0) then
   probez=z0
  elseif (probez>z0+lz) then
   probez=z0+lz
  endif


  meanustar=ustarsurfin

  write (*,*) "set"
 end subroutine READPARAMS

 subroutine ReadIC(U,V,W,Pr)
 real(KND) U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
 integer i,j,k
   open(11,file="in.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read (11,*)
   enddo
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
      read (11,*) Pr(i,j,k)
     enddo
    enddo
   enddo
   if (buoyancy==1) then
    do i=1,3
     read (11,*)
    enddo
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       read (11,*) temperature(i,j,k)
      enddo
     enddo
    enddo
   endif
   close(11)

   open(11,file="Uin.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read (11,*)
   enddo
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
      read (11,*) U(i,j,k)
     enddo
    enddo
   enddo
   close(11)

   open(11,file="Vin.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read (11,*)
   enddo
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
      read (11,*) V(i,j,k)
     enddo
    enddo
   enddo
   close(11)

   open(11,file="Win.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read (11,*)
   enddo
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
      read (11,*) W(i,j,k)
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
  
  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)
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
   freetempgradient=0.01!0.0035 !K/m
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
  
  call Bound_CondU(U)
  call Bound_CondV(V)
  call Bound_CondW(W)
  call Bound_Pr(Pr)
  call Pr_Correct(U,V,W,Pr,1._KND)


   if (sgstype==1) then
                     call Smag(U,V,W)
   elseif (sgstype==2) then
                     call DynSmag(U,V,W)
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
  write (*,*) "set"
 endif
 endsubroutine INITCONDS

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

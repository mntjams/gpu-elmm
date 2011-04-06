module PARAMETERS
implicit none
  integer,parameter :: DBL=KIND(1.0D0),SNG=KIND(1.0E0),KND=SNG,SINT=selected_int_kind(2),SLOG=SINT !KND is default real kind for the whole program
  real(KND) :: pi !computed at the first line of main
  real(KND),parameter :: Karman=0.41
  real(KND),parameter :: BoltzC=1.3806503E-23
  !integer nx,ny,nz
  real(DBL) dt,starttime,endtime        !active time step
  real(KND) dxmin,dymin,dzmin,lx,ly,lz,CFL,Uref  !minimum grid spacing, dimensions of the domain
  real(KND) Re,Prandtl !1/molecular viscosity, viscosity/thermal diffusivity
  real(KND) prgradientx,prgradienty,temperature_ref,grav_acc,freetempgradient,coriolisparam
  real(KND) SHEARG,Uinlet,ustarsurfin,z0inlet,z0W,z0E,z0S,z0N,z0B,z0T,stressgradin,urefin,zrefin,powerexpin,windangle
  real(KND),dimension(:),allocatable:: scalsrcx,scalsrcy,scalsrcz
  real(KND) totalscalsource

  integer,dimension(:),allocatable:: scalsrci,scalsrcj,scalsrck
  integer pointscalsource

  real(KND),dimension(1:3,1:3):: relativestress


  real(KND) epsCN,epsPoisson,eps,limparam,debugparam
  real(DBL) time
  real(KND) timefram1,timefram2,framedimension,slicedir,slicex,timeavg1,timeavg2

  integer tempmet,poissmet,convmet,masssourc,frames,steady
  integer tasktype,averaging,projectiontype,limitertype,impldiff,wallmodeltype,sgstype,fullstress
  integer buoyancy,computescalars,partdistrib,computedeposition,computegravsettling
  integer maxCNiter,maxPOISSONiter,maxiter,endstep
  integer Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,nt,Xup,a,Xd,Yh,Yd,Yu,Zu,Zd
  integer inlettype,gridtype,profiletype


  real(KND),allocatable:: Uin(:,:),Vin(:,:),Win(:,:),Uoutb(:,:),Tempin(:,:) !1:Unz
  real(KND),allocatable:: xU(:),xPr(:),yV(:),yPr(:),zW(:),zPr(:)          !coordinates of grid points
  real(KND),allocatable:: dxU(:),dxPr(:),dyV(:),dyPr(:),dzW(:),dzPr(:)    !dxPr(i)=xU(i)-xU(i-1), dxU(i)=xPr(i+1)-xPr(i)

  real(KND) xheading !true geographic heading of the x axis

  real(KND),allocatable,dimension(:,:,:):: temperature,temperatureavg


  real(KND),allocatable:: Visc(:,:,:),TDiff(:,:,:),tstress(:,:,:,:,:)  !(turbulent) vicosity, (turbulent) thermal diffusivity
  integer(sint),allocatable,dimension(:,:,:):: Utype,Vtype,Wtype,Prtype !number of solid body inside the point is or 0

  logical(slog),allocatable,dimension(:,:,:):: REDBLACKU,REDBLACKV,REDBLACKW,REDBLACKPR !Fs and Ts for red-black iteration

  real(KND) x0,y0,z0
  integer step,deb !for debugging purposes
  real(KND) meanustar
  integer patternnx,patternny,patternnz,boxnx,boxny,boxnz,tilenx,tileny,tilenz
  logical xgridfromfile,ygridfromfile,zgridfromfile
  integer initcondsfromfile 
 
  integer BtypeW,BtypeE,BtypeS,BtypeN,BtypeB,BtypeT       !boundary condition types see below for values
  integer TBtypeW,TBtypeE,TBtypeS,TBtypeN,TBtypeB,TBtypeT !boundary condition types for temperature
  integer ScalBtypeW,ScalBtypeE,ScalBtypeS,ScalBtypeN,ScalBtypeB,ScalBtypeT !boundary condition types for scalars


  real(KND) SsideU,NsideU,BsideU,TsideU                     !velocities on boundaries in case of Dirichlet BC
  real(KND) SsideV,NsideV,BsideV,TsideV
  real(KND) SsideW,NsideW,BsideW,TsideW
  real(KND) SsideTemp,NsideTemp,BsideTemp,TsideTemp,WsideTemp,EsideTemp  !temperatures or temperature fluxes on boundaries
  real(KND),allocatable:: BsideTArr(:,:),BsideTFLArr(:,:)

  real(KND) SsideScal,NsideScal,BsideScal,TsideScal,WsideScal,EsideScal  !scalarss or scalar fluxes on boundaries

  integer vtkformat

  integer,parameter:: NOSLIP=1, FREESLIP=2, PERIODIC=3, DIRICHLET=4, NEUMANN=5, CONSTFLUX=6,&  !boundary condition types
                        TURBULENTINLET=7, FREESLIPBUFF=8, OUTLETBUFF=9, INLETFROMFILE=10
  integer,parameter:: NOINLET=0,CONSTANT=1,SHEAR=2,PARABOLIC=3,CONSTPROF=1,LOGPROF=2,POWERPROF=3  !inlet profile types
  integer,parameter:: GENERALGRID=1,UNIFORMGRID=2
  integer,parameter:: minmodlim=1,extminmodlim=2,gammalim=3,vanalbadalim=4,vanleerlim=5,superbeelim=6
  integer,parameter:: textvtk=1,binaryvtk=2
  real(KND) probex,probey,probez

endmodule PARAMETERS

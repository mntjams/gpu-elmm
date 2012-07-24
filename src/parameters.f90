module PARAMETERS
implicit none

  save

  integer,parameter :: DBL=selected_real_kind(p=15,r=200),SNG=selected_real_kind(p=6,r=37)
  integer,parameter :: KND=DBL                                         !KND is default real kind for the whole program

  integer,parameter :: TIM=KND                                         !Kind for time variables, can be double for very small timesteps.
                                                                       !It may affect performance

  integer,parameter :: SINT=kind(1)!selected_int_kind(2)               !To save memory a smaller type can be used for some integer
  integer,parameter :: SLOG=SINT                                       ! and logical arrays. Note the same KIND value is guaranteed
                                                                       ! the default intrinsic types.
                                                                       !This can have some negative effect on speed however
  real(KND)           :: pi !computed at the first line of main
  real(KND),parameter :: Karman=0.41_KND
  real(KND),parameter :: BoltzC=1.3806503e-23_KND

  integer :: GPU = 0   !HMPP support


  integer :: Unx,Uny,Unz     !dimensions of grid for velocity component U

  integer :: Vnx,Vny,Vnz     !dimensions of grid for velocity component V

  integer :: Wnx,Wny,Wnz     !dimensions of grid for velocity component W

  integer :: Prnx,Prny,Prnz  !dimensions of grid for pressure


  integer   :: nt              !maximum number of time steps

  real(TIM) :: dt,starttime,endtime        !active time step

  real(KND) :: dxmin,dymin,dzmin,lx,ly,lz,CFL,Uref  !minimum grid spacing, dimensions of the domain


  real(KND) :: Re,Prandtl !1/molecular viscosity, viscosity/thermal diffusivity


  real(KND) :: prgradientx,prgradienty,temperature_ref,grav_acc,coriolisparam

  real(KND) :: SHEARG,Uinlet,ustarsurfin

  real(KND) :: z0inlet,z0W,z0E,z0S,z0N,z0B,z0T

  real(KND) :: stressgradin,urefin,zrefin,powerexpin

  real(KND) :: windangle

  real(KND),dimension(:),allocatable :: scalsrcx,scalsrcy,scalsrcz

  real(KND) :: totalscalsource

  integer,dimension(:),allocatable :: scalsrci,scalsrcj,scalsrck

  integer :: scalsourcetype

  real(KND),dimension(1:3,1:3) :: relativestress


  real(KND) :: epsCN,epsPoisson,eps,debugparam
  real(TIM) :: time
  real(TIM) :: timefram1,timefram2,timeavg1,timeavg2

  integer :: tempmet,poissmet,convmet,masssourc,frames,steady
  integer :: tasktype,averaging,projectiontype,impldiff,wallmodeltype,sgstype,fullstress
  integer :: buoyancy,computescalars,partdistrib,computedeposition,computegravsettling
  integer :: maxCNiter,maxPOISSONiter,maxiter,endstep


  integer :: inlettype,gridtype,profiletype


  real(KND),allocatable :: Uin(:,:),Vin(:,:),Win(:,:),Uoutb(:,:),Tempin(:,:) !1:Unz
  real(KND),allocatable :: xU(:),xPr(:),yV(:),yPr(:),zW(:),zPr(:)          !coordinates of grid points
  real(KND),allocatable :: dxU(:),dxPr(:),dyV(:),dyPr(:),dzW(:),dzPr(:)    !dxPr(i)=xU(i)-xU(i-1), dxU(i)=xPr(i+1)-xPr(i)

  real(KND) :: xheading !true geographic heading of the x axis

  real(KND),allocatable,dimension(:,:,:) :: temperature,temperatureavg


  real(KND),allocatable :: Visc(:,:,:),TDiff(:,:,:),tstress(:,:,:,:,:)  !(turbulent) viscosity, (turbulent) thermal diffusivity

  integer(sint),allocatable,dimension(:,:,:) :: Utype,Vtype,Wtype,Prtype !number of solid body inside which the point is or 0
  integer,allocatable,dimension(:,:) :: Unull,Vnull,Wnull                !indexes of points to be nulled every timestep

  integer   :: nUnull,nVnull,nWnull  !second dimension of arrays above (number of points)

  real(KND) :: x0,y0,z0
  integer   :: step
  integer   :: deb !for debugging purposes

  logical   :: xgridfromfile,ygridfromfile,zgridfromfile
  integer   :: initcondsfromfile

  integer,parameter :: Ea=1,We=2,So=3,No=4,Bo=5,To=6

  integer,dimension(6) :: Btype     !boundary condition types see below for values
  integer,dimension(6) :: TBtype    !boundary condition types for temperature
  integer,dimension(6) :: ScalBtype !boundary condition types for scalars


  real(KND),dimension(3,6)   :: sideU    !velocities on boundaries in case of Dirichlet BC
  real(KND),dimension(6)     :: sideTemp !temperatures or temperature fluxes on boundaries
  real(KND),allocatable      :: BsideTArr(:,:),BsideTFLArr(:,:)

  real(KND),dimension(6) :: sideScal  !scalars or scalar fluxes on boundaries

  integer :: vtkformat

  integer,parameter :: NOSLIP=1, FREESLIP=2, PERIODIC=3, DIRICHLET=4, NEUMANN=5, CONSTFLUX=6,&  !boundary condition types
                        TURBULENTINLET=7, FREESLIPBUFF=8, OUTLETBUFF=9, INLETFROMFILE=10
  integer,parameter :: NOINLET=0,CONSTANT=1,SHEAR=2,PARABOLIC=3,CONSTPROF=1,LOGPROF=2,POWERPROF=3  !inlet profile types
  integer,parameter :: GENERALGRID=1,UNIFORMGRID=2
  integer,parameter :: textvtk=1,binaryvtk=2
  integer,parameter :: SmagorinskyModel=1,DynSmagorinskyModel=2,VremanModel=3,StabSmagorinskyModel=4
  integer,parameter :: PointSource=1,LineSource=2,AreaSource=3,VolumeSource=4

  integer :: debuglevel !amount of information to write out

endmodule PARAMETERS

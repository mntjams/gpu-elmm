module KINDS
!   use iso_fortran_env
  implicit none

  integer,parameter :: int8   = selected_int_kind(1)
  integer,parameter :: int32  = selected_int_kind(9)
  integer,parameter :: int64  = selected_int_kind(10)
  integer,parameter :: real32 = selected_real_kind(p=6,r=37)
  integer,parameter :: real64 = selected_real_kind(p=15,r=200)

  integer,parameter :: DBL = real64, SNG = real32

#ifdef DPREC
  integer,parameter :: knd = DBL                                       !knd is default real kind for the whole program, choosing double
#else
  integer,parameter :: knd = SNG                                       !knd is default real kind for the whole program, choosing single
#endif

  integer,parameter :: TIM = knd                                       !Kind for time variables, can be double for very small timesteps.
                                                                       !It may affect performance

  integer,parameter :: SINT = kind(1)                                  !To save memory a smaller type can be used for some integer
  integer,parameter :: SLOG = SINT                                     ! and logical arrays. Note the same KIND value is guaranteed
                                                                       ! the default intrinsic types.
                                                                       !This can have some negative effect on speed however
end module KINDS

module PARAMETERS

  use KINDS

  implicit none

  save

  real(knd)           :: pi !computed at the first line of main
  real(knd),parameter :: Karman=0.41_knd
  real(knd),parameter :: BoltzC=1.3806503e-23_knd

  integer :: Unx,Uny,Unz     !dimensions of grid for velocity component U

  integer :: Vnx,Vny,Vnz     !dimensions of grid for velocity component V

  integer :: Wnx,Wny,Wnz     !dimensions of grid for velocity component W

  integer :: Prnx,Prny,Prnz  !dimensions of grid for pressure


  integer   :: nt              !maximum number of time steps

  real(TIM) :: dt,starttime,endtime        !active time step

  real(knd) :: dxmin,dymin,dzmin,CFL,Uref  !minimum grid spacing, dimensions of the domain


  real(knd) :: Re = 70000, Prandtl = 0.7!1/molecular viscosity, viscosity/thermal diffusivity


  real(knd) :: prgradientx = 0, prgradienty = 0, temperature_ref = 295, grav_acc = 0, coriolisparam = 0

  real(knd) :: SHEARG,Uinlet,ustarsurfin

  real(knd) :: z0inlet,z0W,z0E,z0S,z0N,z0B,z0T

  real(knd) :: stressgradin,urefin,zrefin,powerexpin

  real(knd) :: windangle

  real(knd),dimension(:),allocatable :: scalsrcx,scalsrcy,scalsrcz

  real(knd) :: totalscalsource

  integer,dimension(:),allocatable :: scalsrci,scalsrcj,scalsrck

  integer :: scalsourcetype

  real(knd),dimension(1:3,1:3) :: relativestress


  real(knd) :: epsCN,epsPoisson,eps,debugparam
  real(TIM) :: time
  real(TIM) :: timefram1,timefram2,timeavg1,timeavg2

  integer :: tempmet,poissmet,convmet,masssourc,frames,steady
  integer :: tasktype,averaging,impldiff

  integer :: enable_buoyancy = 0 !1 if enabled, zero otherwise
  integer :: enable_moisture = 0
  integer :: enable_liquid   = 0 !enable condensation of water vapor
  integer :: num_of_scalars  = 0

  integer :: partdistrib,computedeposition,computegravsettling
  integer :: maxCNiter,maxPOISSONiter,maxiter,endstep
  integer :: projectiontype,correctcompatibility = 0

  integer :: inlettype,gridtype,profiletype


  real(knd),allocatable :: Uin(:,:),Vin(:,:),Win(:,:),Uoutb(:,:),Tempin(:,:) !1:Unz
  real(knd),allocatable :: xU(:),xPr(:),yV(:),yPr(:),zW(:),zPr(:)          !coordinates of grid points
  real(knd),allocatable :: dxU(:),dxPr(:),dyV(:),dyPr(:),dzW(:),dzPr(:)    !dxPr(i)=xU(i)-xU(i-1), dxU(i)=xPr(i+1)-xPr(i)

  real(knd) :: xheading !true geographic heading of the x axis

  real(knd),allocatable :: Visc(:,:,:),TDiff(:,:,:)  !(turbulent) viscosity, (turbulent) thermal diffusivity

  integer(sint),allocatable,dimension(:,:,:) :: Utype,Vtype,Wtype,Prtype !number of solid body inside which the point is or 0
  integer,allocatable,dimension(:,:) :: Unull,Vnull,Wnull                !indexes of points to be nulled every timestep

  integer   :: nUnull,nVnull,nWnull  !second dimension of arrays above (number of points)

  integer   :: step

  logical   :: xgridfromfile,ygridfromfile,zgridfromfile
  integer   :: initcondsfromfile

  integer,parameter :: Ea=1,We=2,So=3,No=4,Bo=5,To=6

  integer,dimension(6) :: Btype     !boundary condition types see below for values
  integer,dimension(6) :: TBtype    !boundary condition types for temperature
  integer,dimension(6) :: ScalBtype !boundary condition types for scalars


  real(knd),dimension(3,6)   :: sideU    !velocities on boundaries in case of Dirichlet BC
  real(knd),dimension(6)     :: sideTemp !temperatures or temperature fluxes on boundaries
  real(knd),allocatable      :: BsideTArr(:,:),BsideTFLArr(:,:)

  real(knd),dimension(6) :: sideScal  !scalars or scalar fluxes on boundaries

  integer,parameter :: NOSLIP=1, FREESLIP=2, PERIODIC=3, DIRICHLET=4, NEUMANN=5, CONSTFLUX=6,&  !boundary condition types
                        TURBULENTINLET=7, FREESLIPBUFF=8, OUTLETBUFF=9, INLETFROMFILE=10
  integer,parameter :: NOINLET=0,CONSTANT=1,SHEAR=2,PARABOLIC=3,CONSTPROF=1,LOGPROF=2,POWERPROF=3  !inlet profile types
  integer,parameter :: GENERALGRID=1,UNIFORMGRID=2
  integer,parameter :: SubgridModel=1,SigmaModel=2,VremanModel=3,StabSubgridModel=4
  integer,parameter :: PointSource=1,LineSource=2,AreaSource=3,VolumeSource=4
  integer :: debuglevel = 0 !amount of information to write out

! #ifdef __HMPP
 integer :: GPU = 0   !If the GPU is allocated

 type hmppWMpoint   !points in which we apply wall model

  integer   :: xi
  integer   :: yj
  integer   :: zk

  real(knd) :: distx
  real(knd) :: disty
  real(knd) :: distz

  real(knd) :: z0=0
  real(knd) :: ustar=1
  real(knd) :: temp=0
  real(knd) :: tempfl=0

 endtype hmppWMpoint
! #endif

endmodule PARAMETERS

module RK3
  use Kinds

  real(knd),parameter :: RK_alpha(3) = [ 4._knd/15._knd, 1._knd/15._knd,  1._knd/6._knd  ]
  real(knd),parameter :: RK_beta(3)  = [ 8._knd/15._knd, 5._knd/12._knd,  3._knd/4._knd  ]
  real(knd),parameter :: RK_rho(3)   = [ 0._knd,       -17._knd/60._knd, -5._knd/12._knd ]
  integer,parameter   :: RKstages = 3

end module RK3

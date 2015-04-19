module Kinds

  use iso_fortran_env
  use iso_c_binding, only: c_int, c_float, c_double, c_size_t

  implicit none

!for reference when iso_fortran_env is not available:
!   integer, parameter :: int8   = selected_int_kind(1)
!   integer, parameter :: int32  = selected_int_kind(9)
!   integer, parameter :: int64  = selected_int_kind(10)
!   integer, parameter :: real32 = selected_real_kind(p=6, r=37)
!   integer, parameter :: real64 = selected_real_kind(p=15, r=200)

  integer, parameter :: dbl = real64, sng = real32

#ifdef DPREC
  integer, parameter :: knd = DBL                                       !knd is default real kind for the whole program, choosing double
#else
  integer, parameter :: knd = SNG                                       !knd is default real kind for the whole program, choosing single
#endif

  integer, parameter :: tim = knd                                       !Kind for time variables, can be double for very small timesteps.
                                                                       !It may affect performance

  integer, parameter :: sint = kind(1)                                  !To save memory a smaller type can be used for some integer
  integer, parameter :: slog = sint                                     ! and logical arrays. Note the same KIND value is guaranteed
                                                                       ! the default intrinsic types.
                                                                       !This can have some negative effect on speed however
end module Kinds



module Parameters

  use Kinds
  use Stop_procedures

  implicit none

  save

  real(knd), parameter :: pi = acos(-1.0_knd)
  real(knd), parameter :: Karman = 0.41_knd
  real(knd), parameter :: BoltzC = 1.3806503e-23_knd

  integer :: Unx, Uny, Unz     !dimensions of grid for velocity component U

  integer :: Vnx, Vny, Vnz     !dimensions of grid for velocity component V

  integer :: Wnx, Wny, Wnz     !dimensions of grid for velocity component W

  integer :: Prnx, Prny, Prnz  !dimensions of grid for pressure
  
  !relevant for MPI
  integer :: gPrnx, gPrny, gPrnz !global grid dimensions
  integer :: gUnx, gUny, gUnz
  integer :: gVnx, gVny, gVnz
  integer :: gWnx, gWny, gWnz
  
  integer :: offset_to_global_x = 0, offset_to_global_y = 0, offset_to_global_z = 0
  integer :: gPrns(3), offsets_to_global(3) = 0
  

  integer   :: max_number_of_time_steps    !maximum number of time steps

  real(TIM) :: start_time, end_time

  real(knd) :: dxmin, dymin, dzmin, CFL, Uref  !minimum grid spacing, dimensions of the domain


  real(knd) :: Re = 70000, Prandtl = 0.7!1/molecular viscosity, viscosity/thermal diffusivity


  real(knd) :: pr_gradient_x = 0, pr_gradient_y = 0
  real(knd), allocatable :: pr_gradient_profile_x(:), pr_gradient_profile_y(:)

  logical :: enable_pr_gradient_x_uniform = .false.
  logical :: enable_pr_gradient_y_uniform = .false.
  logical :: enable_pr_gradient_x_profile = .false.
  logical :: enable_pr_gradient_y_profile = .false.


  real(knd) :: grav_acc = 0, Coriolis_parameter = 0

  real(knd) :: top_pressure !mean pressure at the top boundary - calculated
  real(knd) :: bottom_pressure = 101325.0

  real(knd) :: ShearInletTypeParameter, Uinlet

  real(knd) :: z0W, z0E, z0S, z0N, z0B, z0T

  real(knd) :: windangle = 0

  real(knd), dimension(:), allocatable :: scalsrcx, scalsrcy, scalsrcz

  real(knd) :: totalscalsource

  integer, dimension(:), allocatable :: scalsrci, scalsrcj, scalsrck

  integer :: scalsourcetype


  real(knd) :: epsCN, epsPoisson, eps, debugparam
  real(TIM) :: time
  real(TIM) :: timefram1, timefram2, timeavg1, timeavg2

  integer :: tempmet, poissmet, convmet, masssourc, frames, steady
  integer :: tasktype, averaging

  logical :: explicit_diffusion = .true.

  logical :: enable_buoyancy = .false. !1 if enabled, zero otherwise
  logical :: enable_moisture = .false.
  logical :: enable_liquid   = .false. !enable condensation of water vapor
  integer :: num_of_scalars  = 0

  logical :: enable_radiation = .false.

  integer :: partdistrib, computedeposition, computegravsettling
  integer :: maxCNiter, maxPOISSONiter, endstep
  integer :: projectiontype, correctcompatibility = 0

  integer :: inlettype, gridtype, profiletype


  real(knd), allocatable :: Uin(:,:), Vin(:,:), Win(:,:), Uoutb(:,:)

  real(knd), allocatable :: TempIn(:,:), MoistIn(:,:)

  real(knd) :: gxmin, gxmax, gymin, gymax, gzmin, gzmax

  real(knd), allocatable :: xU(:), xPr(:), yV(:), yPr(:), zW(:), zPr(:)          !coordinates of grid points
  real(knd), allocatable :: dxU(:), dxPr(:), dyV(:), dyPr(:), dzW(:), dzPr(:)    !dxPr(i)=xU(i)-xU(i-1), dxU(i)=xPr(i+1)-xPr(i)

  real(knd) :: x_axis_azimuth = 90!true geographic heading of the x axis in degrees

  real(knd), allocatable :: Viscosity(:,:,:), TDiff(:,:,:)  !(turbulent) viscosity, (turbulent) thermal diffusivity

  integer(sint), allocatable, dimension(:,:,:) :: Utype, Vtype, Wtype, Prtype !number of solid body inside which the point is or 0
  integer, allocatable, dimension(:,:) :: Unull, Vnull, Wnull                !indexes of points to be nulled every timestep

  integer   :: nUnull, nVnull, nWnull  !second dimension of arrays above (number of points)

  logical   :: xgridfromfile, ygridfromfile, zgridfromfile
  integer   :: initcondsfromfile

  integer, parameter :: We=1, Ea=2, So=3, No=4, Bo=5, To=6

  real(knd) :: temperature_ref = 295
  real(knd) :: moisture_ref = 0.001 !TODO: compute from relative humidity

  integer, dimension(6) :: Btype      !boundary condition types for velocity, see below for values
  integer, dimension(6) :: TempBtype  !boundary condition types for temperature
  integer, dimension(6) :: MoistBtype !boundary condition types for temperature
  integer, dimension(6) :: ScalBtype  !boundary condition types for scalars
  
  integer, dimension(6) :: PoissonBtype !boundary conditions for the pressure solver


  real(knd), dimension(3,6)   :: sideU     !velocities on boundaries in case of Dirichlet BC
  real(knd), dimension(6)     :: sideTemp  !temperatures or temperature fluxes on boundaries
  real(knd), dimension(6)     :: sideMoist !moistures or moisture fluxes on boundaries
  real(knd), dimension(6)     :: sideScal  !scalars or scalar fluxes on boundaries

  real(knd), allocatable      :: BsideTArr(:,:), BsideTFLArr(:,:)
  real(knd), allocatable      :: BsideMArr(:,:), BsideMFLArr(:,:)


  integer, parameter :: ScalarTypeTemperature = 1, &
                       ScalarTypeMoisture = 2, &
                       ScalarTypePassive = 3

  integer, parameter :: NOSLIP=1, FREESLIP=2, PERIODIC=3, DIRICHLET=4, NEUMANN=5, CONSTFLUX=6, &  !boundary condition types
                        TURBULENTINLET=7, FREESLIPBUFF=8, OUTLETBUFF=9, INLETFROMFILE=10, RADIATION=7, &
                        MO_TEMPERATURE=10, &
                        NEUMANN_BUFF=8, AUTOMATICFLUX=11, WALL_DIRICHLET=12, WALL_FLUX=13, &
                        MPI_BOUNDS=1000, MPI_BOUNDARY=1000, MPI_PERIODIC=1001
  !inlet types
  integer, parameter :: ZeroInletType=0, ConstantInletType=1, ShearInletType=2, &
                        ParabolicInletType=3, TurbulentInletType=4, FromFileInletType=5, &
                        GeostrophicInletType=6
  !inlet profile types
  integer, parameter :: CONSTPROF=1, LOGPROF=2, POWERPROF=3 
  integer, parameter :: GENERALGRID=1, UNIFORMGRID=2
  integer, parameter :: SmagorinskyModel=1, SigmaModel=2, VremanModel=3, StabSubgridModel=4
  integer, parameter :: PointSource=1, LineSource=2, AreaSource=3, VolumeSource=4

  integer(c_int), bind(C, name="debuglevel") :: debuglevel = 0 !amount of information to write out
  
  logical :: init_phase = .true., run_phase = .false.
  
  character(:), allocatable :: output_dir, image_input_dir
  character(2048) :: scratch_dir = ""

  integer(dbl) :: timer_rate
  real(knd) :: time_steps_time = 0
  real(knd) :: poisson_solver_time = 0

end module Parameters


module PhysicalProperties
  use Parameters

  implicit none

  real(knd), parameter :: rho_air_ref = 1.196 !kg.m^-3

  real(knd), parameter :: Cp_air_ref = 1005 !J.kg^-1.K^-1

  real(knd), parameter :: Cv_air_ref = 718 !J.kg^-1.K^-1

  real(knd), parameter :: Lv_water_ref = 2442000 !J.kg^-1

end module


module RK3
  use Kinds

  implicit none

  real(knd), parameter :: RK_alpha(3) = [ 4._knd/15._knd, 1._knd/15._knd,  1._knd/6._knd  ]
  real(knd), parameter :: RK_beta(3)  = [ 8._knd/15._knd, 5._knd/12._knd,  3._knd/4._knd  ]
  real(knd), parameter :: RK_rho(3)   = [ 0._knd,       -17._knd/60._knd, -5._knd/12._knd ]
  integer, parameter   :: RK_stages = 3

end module RK3

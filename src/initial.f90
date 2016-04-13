module Initial

  use Parameters
  use ArrayUtilities, only: avg
  use Limiters, only: limiter_parameter, limiter_type
  use Multigrid, only: SetMGParams
  use Multigrid2d, only: SetMGParams2d
  use Pressure
  use Boundaries
  use ScalarBoundaries
  use Outputs, only: store, display, probes, scalar_probes, ReadProbes
  use Scalars
  use Filters, only: filtertype, filter_ratios
  use Subgrid
  use Turbinlet, only: default_turbulence_generator, GetInletFromFile
  use SolarRadiation, only: InitSolarRadiation
  use SolidBodies, only: obstacles_file, roughness_file, displacement_file, InitSolidBodies
  use ImmersedBoundary, only: GetSolidBodiesBC, InitIBPFluxes!, SetIBPFluxes
  use VolumeSources!, only: InitVolumeSources, InitVolumeSourceBodies, ScalarFlVolume, ScalarFlVolumesContainer
  use AreaSources, only: InitAreaSources
  use LineSources, only: InitLineSources
  use PointSources, only: InitPointSources
  use Wallmodels
  use Tiling, only: tilesize,InitTiles
  use FreeUnit, only: newunit
  use Puffs, only: InitPuffSources
  use custom_par
#ifdef PAR
  use exchange_par
  use domains_bc_par
#endif

  implicit none

  private
  public  ReadConfiguration, InitialConditions, InitBoundaryConditions, probes_file, scalar_probes_file

  real(knd) :: x0,y0,z0 !domain boundaries, will become xU(0), yV(0), zW(0)
  real(knd) :: lx,ly,lz !domain extents

  character(80) :: probes_file = ""
  character(80) :: scalar_probes_file = ""

  type spline_coefs
    real(knd), allocatable :: z(:)
    real(knd), allocatable :: cu(:,:), cv(:,:)
  end type
    

contains


 subroutine ReadConfiguration
   use StaggeredFrames, only: rrange, TFrameTimes, TSaveFlags, &
                              TStaggeredFrameDomain,  AddDomain
   use VTKFrames, only: TFrameFlags, &
                              TFrameDomain,  AddDomain
   use PoisFFT, only: PoisFFT_NeumannStag, PoisFFT_Periodic
   integer ::  lmg,minmglevel,bnx,bny,bnz,mgncgc,mgnpre,mgnpost,mgmaxinnerGSiter
   real(knd) :: mgepsinnerGS
   integer ::  i,io,io2,itmp
   integer :: numframeslices

   character(len = 1024) :: command_line, msg
   integer :: exenamelength
   integer :: unit

   type(rrange) :: range
   type(TFrameTimes) :: frame_times
   type(TSaveFlags) :: frame_save_flags
   type(TFrameFlags) :: frame_flags
   character(10) :: domain_label
   integer :: num_staggered_domains
   integer :: number_of_probes, number_of_scalar_probes

   integer :: dimension,direction
   real(knd) :: position

   integer :: masssourc

   real(knd) :: Re = 70000, Prandtl = 0.7!1/molecular viscosity, viscosity/thermal diffusivity

   namelist /obstacles/ obstacles_file, roughness_file, displacement_file

   interface get
     procedure chget1
     procedure lget1, lget2, lget3
     procedure iget1, iget2, iget3
     procedure rget1, rget2, rget3
     procedure rgetv3
   end interface

   interface

     subroutine CustomConfiguration_First
     end subroutine

     subroutine CustomConfiguration_Last
     end subroutine

   end interface

   image_input_dir = "input/"
   output_dir = "output/"


   call init_random_seed


   call newunit(unit)
  
   !try read scratch_dir from environment, command line has priority
   call get_scratch

!    the command line arguments can also be specified in cmd.conf
   call read_cmd_conf

   call parse_command_line


   !the actual command line has priority over cmd.conf
   call read_command_line
   !NOTE: it is parsed one more time lower in this subroutine
   call parse_command_line

#ifdef CUSTOM_CONFIG
   call CustomConfiguration_First
#endif

   open(unit,file="main.conf",status="old",action="read")
   call get(advection_method)
   call get(limiter_type)
   call get(limiter_parameter)

   call get(masssourc)
   masssourc = 1 !Seems to be necessary for stability, change with caution.
   enable_ibm_mass_sources = (masssourc>0)

   call get(steady)
   call get(task_type)
   call get(initcondsfromfile)
   call get(timeavg1)
   call get(timeavg2)
   call get(Re)
   if (master) write(*,*) "Re=",Re
   call get(eps)
   if (master) write(*,*) "eps=",eps
   call get(maxCNiter)
   if (master) write(*,*) "maxCNiter=",maxCNiter
   call get(epsCN)
   if (master) write(*,*) "epsCN=",epsCN
   call get(maxPOISSONiter)
   if (master) write(*,*) "maxPOISSONiter=",maxPOISSONiter
   call get(epsPOISSON)
   if (master) write(*,*) "epsPOISSON=",epsPOISSON
   call get(debugparam)
   if (master) write(*,*) "debug parameter=",debugparam
   close(unit)


   open(unit,file="les.conf",status="old",action="read")
   call get(sgstype)
   call get(filtertype)

   if (filtertype > size(filter_ratios)) then
     if (master) write(*,*) "Chosen filter type does not exist. Maximum index is:",size(filter_ratios)
     call error_stop
   end if

   call get(wallmodeltype)
   close(unit)


   open(unit,file="grid.conf",status="old",action="read")
   call get(xgridfromfile)
   call get(ygridfromfile)
   call get(zgridfromfile)
   read(unit,fmt='(/)')

   read(unit,*) x0
   if (master) write(*,*) "x0=",x0
   call get(y0)
   if (master)write(*,*) "y0=",y0
   call get(z0)
   if (master) write(*,*) "z0=",z0
   read(unit,fmt='(/)')

   read(unit,*) lx
   if (lx>0) then
     if (master) write(*,*) "lx=",lx
   else
     if (master) write (*,*) "Domain length in x direction must be positive."
     call error_stop
   end if
   read(unit,fmt='(/)')

   read(unit,*) ly
   if (ly>0) then
     if (master) write(*,*) "ly=",ly
   else
     if (master) write (*,*) "Domain length in y direction must be positive."
     call error_stop
   end if
   read(unit,fmt='(/)')

   read(unit,*) lz
   if (lz>0) then
     if (master) write(*,*) "lz=",lz
   else
     if (master) write (*,*) "Domain length in z direction must be positive."
     call error_stop
   end if
   read(unit,fmt='(/)')

   read(unit,*) Prnx
   if (Prnx>0) then
     if (master) write(*,*) "nx=",Prnx
   else
     if (master) write (*,*) "Number of cells in x direction must be positive."
     call error_stop
   end if
   read(unit,fmt='(/)')

   read(unit,*) Prny
   if (Prny>0) then
     if (master) write(*,*) "ny=",Prny
   else
     if (master) write (*,*) "Number of cells in y direction must be positive."
     call error_stop
   end if
   read(unit,fmt='(/)')

   read(unit,*) Prnz
   if (Prnz>0) then
     if (master) write(*,*) "nz=",Prnz
   else
     if (master) write (*,*) "Number of cells in z direction must be positive."
     call error_stop
   end if
   close(unit)

   open(unit,file="boundconds.conf",status="old",action="read")
   call get(Btype(We))
   call get(Btype(Ea))
   call get(Btype(So))
   call get(Btype(No))
   call get(Btype(Bo))
   call get(Btype(To))
   call get(sideU(1,So))
   call get(sideU(2,So))
   call get(sideU(3,So))
   call get(sideU(1,No))
   call get(sideU(2,No))
   call get(sideU(3,No))
   call get(sideU(1,Bo))
   call get(sideU(2,Bo))
   call get(sideU(3,Bo))
   call get(sideU(1,To))
   call get(sideU(2,To))
   call get(sideU(3,To))
   call get(z0W)
   call get(z0E)
   call get(z0S)
   call get(z0N)
   call get(z0B)
   call get(z0T)
   close(unit)

   open(unit,file="large_scale.conf",status="old",action="read",iostat = io)
   if (io==0) then

     call get(Coriolis_parameter)
     if (master) write(*,*) "Coriolis_parameter=", Coriolis_parameter

     call get(pr_gradient_x)
     if (master) write(*,*) "pr_gradient_x=", pr_gradient_x
     enable_pr_gradient_x_uniform = pr_gradient_x /= 0

     call get(pr_gradient_y)
     if (master) write(*,*) "pr_gradient_y=", pr_gradient_y
     enable_pr_gradient_y_uniform = pr_gradient_y /= 0

     call get(SubsidenceGradient)
     if (master) write(*,*) "SubsidenceGradient=",SubsidenceGradient
     close(unit)

   else

     if (master) write(*,*) "Warning! Could not open file large_scale.conf. Using defaults."

   end if


   open(unit,file="thermal.conf",status="old",action="read",iostat = io)
   if (io==0) then
     call get(enable_buoyancy)
     call get(Prandtl)
     call get(grav_acc)
     call get(temperature_ref)
     call get(TempBtype(We))
     call get(TempBtype(Ea))
     call get(TempBtype(So))
     call get(TempBtype(No))
     call get(TempBtype(Bo))
     call get(TempBtype(To))
     call get(sideTemp(We))
     call get(sideTemp(Ea))
     call get(sideTemp(So))
     call get(sideTemp(No))
     call get(sideTemp(Bo))
     call get(sideTemp(To))
     close(unit)
   else
     enable_buoyancy = .false.
   end if

   if (enable_buoyancy) then

     open(unit,file="temp_profile.conf",status="old",action="read",iostat = io)

     if (io==0) then
       call get(TemperatureProfile%randomize)
       call get(TemperatureProfile%randomizeTop)
       call get(TemperatureProfile%randomizeAmplitude)
       call get(itmp)

       allocate(TemperatureProfile%Sections(max(itmp,0)))

       do i = 1, size(TemperatureProfile%Sections)
         call get(TemperatureProfile%Sections(i)%top)
         call get(TemperatureProfile%Sections(i)%jump)
         call get(TemperatureProfile%Sections(i)%gradient)
       end do

       close(unit)

     else

       if (master) write(*,*) "Warning! Could not open file temp_profile.conf. Using defaults."
       TemperatureProfile%randomize = 0

       allocate(TemperatureProfile%Sections(0))

     end if

   else

     TempBtype = 0

   end if



   open(unit,file="moisture.conf",status="old",action="read",iostat = io)
   if (io==0) then
     call get(enable_moisture)
     call get(moisture_ref)
     call get(MoistBtype(We))
     call get(MoistBtype(Ea))
     call get(MoistBtype(So))
     call get(MoistBtype(No))
     call get(MoistBtype(Bo))
     call get(MoistBtype(To))
     call get(sideMoist(We))
     call get(sideMoist(Ea))
     call get(sideMoist(So))
     call get(sideMoist(No))
     call get(sideMoist(Bo))
     call get(sideMoist(To))
     close(unit)
   else
     enable_moisture = .false.
   end if

   if (enable_moisture) then

     open(unit,file="moist_profile.conf",status="old",action="read",iostat = io)

     if (io==0) then
       call get(MoistureProfile%randomize)
       call get(MoistureProfile%randomizeTop)
       call get(MoistureProfile%randomizeAmplitude)
       call get(itmp)

       allocate(MoistureProfile%Sections(max(itmp,0)))

       do i = 1, size(MoistureProfile%Sections)
         call get(MoistureProfile%Sections(i)%top)
         call get(MoistureProfile%Sections(i)%jump)
         call get(MoistureProfile%Sections(i)%gradient)
       end do

       close(unit)

     else

       if (master) write(*,*) "Warning! Could not open file moist_profile.conf. Using defaults."
       MoistureProfile%randomize = 0

       allocate(MoistureProfile%Sections(0))

     end if

   else

     MoistBtype = 0

   end if

   open(unit,file="inlet.conf",status="old",action="read")
   call get(inlettype)
   call get(profiletype)
   call get(ShearInletTypeParameter)
   if (master) write(*,*) "G=",ShearInletTypeParameter
   call get(Uinlet)
   if (master) write(*,*) "Uinlet=",Uinlet
   call get(default_turbulence_generator%Ustar_surf_inlet)  !-<u'w'>
   call get(default_turbulence_generator%stress_gradient_inlet) !in relative part per 1m
   call get(default_turbulence_generator%z0_inlet)
   call get(default_turbulence_generator%power_exponent_inlet)
   call get(default_turbulence_generator%z_ref_inlet)
   call get(default_turbulence_generator%U_ref_inlet)
   call get(default_turbulence_generator%relative_stress(1,1))
   call get(default_turbulence_generator%relative_stress(2,2))
   call get(default_turbulence_generator%relative_stress(3,3))
   call get(default_turbulence_generator%relative_stress(1,2))
   call get(default_turbulence_generator%relative_stress(1,3))
   call get(default_turbulence_generator%relative_stress(2,3))
   call get(default_turbulence_generator%T_Lag)
   call get(default_turbulence_generator%L_y)
   call get(default_turbulence_generator%L_z)
   close(unit)

   default_turbulence_generator%relative_stress(2,1) = default_turbulence_generator%relative_stress(1,2)
   default_turbulence_generator%relative_stress(3,1) = default_turbulence_generator%relative_stress(1,3)
   default_turbulence_generator%relative_stress(3,2) = default_turbulence_generator%relative_stress(2,3)

   open(unit,file="scalars.conf",status="old",action="read",iostat = io)
   if (io==0) then
     call get(num_of_scalars)
     call get(computedeposition)
     call get(computegravsettling)
     call get(partdistrib)
     call get(totalscalsource)
     call get(scalsourcetype)

     call get(ScalBtype(We))
     call get(ScalBtype(Ea))
     call get(ScalBtype(So))
     call get(ScalBtype(No))
     call get(ScalBtype(Bo))
     call get(ScalBtype(To))
     call get(sideScal(We))
     call get(sideScal(Ea))
     call get(sideScal(So))
     call get(sideScal(No))
     call get(sideScal(Bo))
     call get(sideScal(To))

     if (scalsourcetype==1) then
       if (partdistrib>0) then

          allocate(partdiam(partdistrib),partrho(partdistrib),percdistrib(partdistrib))

          do i = 1, partdistrib
            call get(partdiam(i))
            call get(partrho(i))
            call get(percdistrib(i))
          end do

       else

          allocate(partdiam(num_of_scalars),partrho(num_of_scalars),percdistrib(num_of_scalars))

          do i = 1, num_of_scalars
            call get(partdiam(i))
            call get(partrho(i))
            call get(percdistrib(i))
          end do
       end if
     end if
     close(unit)
   else
     num_of_scalars = 0
     if (master) write (*,*) "scalars.conf not found, no passive scalars for computation."
   end if


   if (num_of_scalars>0) then
      open(unit, file="line_sources.conf",status="old",action="read",iostat=io)
      
      if (io==0) then
        call get_line_sources
        close(unit)
      end if
   end if

   call get_area_sources("area_sources.conf")

   if (num_of_scalars>0) then
      open(unit, file="point_sources.conf",status="old",action="read",iostat=io)
      
      if (io==0) then
        call get_point_sources
        close(unit)
      end if
   end if

   if (pressure_solution%poisson_solver==POISSON_SOLVER_MULTIGRID) then
     open(unit,file="mgopts.conf",status="old",action="read")
     call get(lmg)
     call get(minmglevel)
     call get(bnx)
     call get(bny)
     call get(bnz)
     call get(mgncgc)
     call get(mgnpre)
     call get(mgnpost)
     call get(mgmaxinnerGSiter)
     call get(mgepsinnerGS)
     close(unit)

     if (Prny==1) then
      call SetMGParams2d(llmg = lmg,lminmglevel = minmglevel,lbnx = bnx,lbnz = bnz,&
                         lmgncgc = mgncgc,lmgnpre = mgnpre,lmgnpost = mgnpost,&
                         lmgmaxinnerGSiter = mgmaxinnerGSiter,lmgepsinnerGS = mgepsinnerGS)
     else
      call SetMGParams(llmg = lmg,lminmglevel = minmglevel,&
                         lbnx = bnx,lbny = bny,lbnz = bnz,&
                         lmgncgc = mgncgc,lmgnpre = mgnpre,lmgnpost = mgnpost,&
                         lmgmaxinnerGSiter = mgmaxinnerGSiter,lmgepsinnerGS = mgepsinnerGS)
     end if
   end if


  

   open(unit,file="output.conf",status="old",action="read",iostat = io)
   if (io==0) then
     call read_namelist_output

     if (.not.enable_buoyancy) then
       store%out_temperature = 0
       store%out_moisture = 0
     end if

     if (.not.enable_buoyancy) then
       store%avg_temperature = 0
       store%avg_moisture = 0
     end if

     if (num_of_scalars < 1) then
       store%scalars = 0
       store%scalars_avg = 0
       store%scalsum_time = 0
       store%scaltotsum_time = 0
     end if
     close(unit)
   else
     if (master) write(*,*) "No output.conf found, defaults will be used."
   end if
   
   
   call get_profiles("profiles.conf")

   !probes_file and scalar_probes_file read from command line
   if (probes_file == "" .and. scalar_probes_file == "") then

     open(unit,file="probes.conf",status="old",action="read",iostat = io)
     if (io==0) then
       call get(number_of_probes)

       allocate(probes(number_of_probes))

       do i = 1, number_of_probes
         call get(probes(i)%x, probes(i)%y, probes(i)%z)
       end do

       scalar_probes = probes

       close(unit)

     else
       allocate(probes(0))
       allocate(scalar_probes(0))
     end if

   else

     if (len_trim(probes_file)>0) then
       call ReadProbes(probes,number_of_probes,probes_file)
     else
       allocate(probes(0))
     end if

     if (len_trim(scalar_probes_file)>0) then
       call ReadProbes(scalar_probes,number_of_scalar_probes,scalar_probes_file)
     else
       allocate(scalar_probes(0))
     end if
   end if


   open(unit,file="obstacles.conf",status="old",action="read",iostat = io)
   if (io==0) then
     read(unit,nml=obstacles,iostat=io,iomsg=msg)
     if (io/=0) then
       if (master) write(*,*) io,"Error reading from obstacles.conf."
       if (master) write(*,*) msg
       if (master) write(*,*) command_line
     end if
     close(unit)
   end if

   if (master) write(*,*) "num_of_scalars",num_of_scalars

   call parse_command_line

   if (master) then

     write(*,*) "Boundaries:"

     write(*,'(a2)',advance='no') " W "
     select case (Btype(We))
       case (BC_NOSLIP)
         write(*,*) "noslip"
       case (BC_FREESLIP)
         write(*,*) "freeslip"
       case (BC_PERIODIC)
         write(*,*) "periodic"
       case (BC_DIRICHLET)
         write(*,*) "dirichlet"
       case (BC_NEUMANN)
         write(*,*) "neumann"
     endselect

     write(*,'(a2)',advance='no') " E "
     select case (Btype(Ea))
       case (BC_NOSLIP)
         write(*,*) "noslip"
       case (BC_FREESLIP)
         write(*,*) "freeslip"
       case (BC_PERIODIC)
         write(*,*) "periodic"
       case (BC_DIRICHLET)
         write(*,*) "dirichlet"
       case (BC_NEUMANN)
         write(*,*) "neumann"
     endselect

     write(*,'(a2)',advance='no') " S "
     select case (Btype(So))
       case (BC_NOSLIP)
         write(*,*) "noslip"
       case (BC_FREESLIP)
         write(*,*) "freeslip"
       case (BC_PERIODIC)
         write(*,*) "periodic"
       case (BC_DIRICHLET)
         write(*,*) "dirichlet"
       case (BC_NEUMANN)
         write(*,*) "neumann"
     endselect

     write(*,'(a2)',advance='no') " N "
     select case (Btype(No))
       case (BC_NOSLIP)
         write(*,*) "noslip"
       case (BC_FREESLIP)
         write(*,*) "freeslip"
       case (BC_PERIODIC)
         write(*,*) "periodic"
       case (BC_DIRICHLET)
         write(*,*) "dirichlet"
       case (BC_NEUMANN)
         write(*,*) "neumann"
     endselect

     write(*,'(a2)',advance='no') " B "
     select case (Btype(Bo))
       case (BC_NOSLIP)
         write(*,*) "noslip"
       case (BC_FREESLIP)
         write(*,*) "freeslip"
       case (BC_PERIODIC)
         write(*,*) "periodic"
       case (BC_DIRICHLET)
         write(*,*) "dirichlet"
       case (BC_NEUMANN)
         write(*,*) "neumann"
     endselect

     write(*,'(a2)',advance='no') " T "
     select case (Btype(To))
       case (BC_NOSLIP)
         write(*,*) "noslip"
       case (BC_FREESLIP)
         write(*,*) "freeslip"
       case (BC_PERIODIC)
         write(*,*) "periodic"
       case (BC_DIRICHLET)
         write(*,*) "dirichlet"
       case (BC_NEUMANN)
         write(*,*) "neumann"
     endselect
     
     if ((Btype(We)==BC_PERIODIC.or.Btype(Ea)==BC_PERIODIC).and.Btype(We)/=Btype(Ea)) &
       call error_stop("Error: Both X boundary conditions must be periodic or not periodic.")
     if ((Btype(So)==BC_PERIODIC.or.Btype(No)==BC_PERIODIC).and.Btype(So)/=Btype(No)) &
       call error_stop("Error: Both Y boundary conditions must be periodic or not periodic.")
     if ((Btype(Bo)==BC_PERIODIC.or.Btype(To)==BC_PERIODIC).and.Btype(Bo)/=Btype(To)) &
       call error_stop("Error: Both Z boundary conditions must be periodic or not periodic.")

   end if



   if (Re > 0) then
     molecular_viscosity = 1 / Re
   else
     molecular_viscosity = 0
   end if

   molecular_diffusivity = molecular_viscosity / Prandtl


   dxmin = lx/(Prnx)
   dymin = ly/(Prny)
   dzmin = lz/(Prnz)

   if (Prnx==1) then
     dxmin = sqrt(dymin*dzmin)
     lx = dxmin
   elseif (Prny==1) then
     dymin = sqrt(dxmin*dzmin)
     ly = dymin
   elseif (Prnz==1) then
     dzmin = sqrt(dxmin*dymin)
     lz = dzmin
   end if

   !can be overwritten after MPI decomposition
   gPrnx = Prnx
   gPrny = Prny
   gPrnz = Prnz   


   if (master) then
      write(*,*) "dxmin ",dxmin
      write(*,*) "dymin ",dymin
      write(*,*) "dzmin ",dzmin

      write(*,*) "lx:",lx
      write(*,*) "ly:",ly
      write(*,*) "lz:",lz
   end if


   if (Btype(We)==BC_TURBULENT_INLET) inlettype = TurbulentInletType
   if (Btype(We)==BC_INLET_FROM_FILE) inlettype = FromFileInletType


   call get_time_stepping("time_stepping.conf")

   if (timeavg2>=timeavg1) then
     averaging = 1
   else
     averaging = 0
   end if

   if (.not.xgridfromfile.and..not.ygridfromfile.and..not.zgridfromfile) then
     gridtype = UNIFORMGRID
     if (master) write(*,*) "Uniform grid"
   else
     gridtype = GENERALGRID
     if (master) write(*,*) "General grid not supported."; call error_stop
   end if

   !Btype might get overwritten by MPI procedures
   do i = We, To
     if (Btype(i)==BC_PERIODIC) then
        PoissonBtype(i) = PoisFFT_PERIODIC
     else
        PoissonBtype(i) = PoisFFT_NeumannStag
     end if
   end do
   
   gxmin = x0
   gxmax = x0 + lx
   gymin = y0
   gymax = y0 + ly
   gzmin = z0
   gzmax = z0 + lz


   if (len_trim(scratch_dir)>0) then
     if (scratch_dir(len_trim(scratch_dir):len_trim(scratch_dir))/="/") &
       scratch_dir = trim(scratch_dir) // "/"
     output_dir = trim(scratch_dir) // output_dir
   end if


   im_xmin = gxmin
   im_ymin = gymin
   im_zmin = gzmin

   im_xmax = gxmax
   im_ymax = gymax
   im_zmax = gzmax

#ifdef PAR
   call get_domains("domains.conf")

   call par_init_grid
   
   call par_init_boundaries
   
   x0 = im_xmin
   y0 = im_ymin
   z0 = im_zmin
#endif


   call get_pressure_solution("pressure_solution.conf")


   !both procedures below use the name of the output directory (affected by MPI)
   call read_frames
   
   call read_staggered_frames


   if (Btype(Ea)==BC_PERIODIC .or. &
       (Btype(Ea)>=BC_MPI_BOUNDS_MIN.and.Btype(Ea)<=BC_MPI_BOUNDS_MAX)) then
                          Unx = Prnx
   else
                          Unx = Prnx-1
   end if
   Uny = Prny
   Unz = Prnz

   Vnx = Prnx
   if (Btype(No)==BC_PERIODIC .or. &
       (Btype(No)>=BC_MPI_BOUNDS_MIN.and.Btype(Ea)<=BC_MPI_BOUNDS_MAX)) then
                          Vny = Prny
   else
                          Vny = Prny-1
   end if
   Vnz = Prnz

   Wnx = Prnx
   Wny = Prny
   if (Btype(To)==BC_PERIODIC .or. &
       (Btype(To)>=BC_MPI_BOUNDS_MIN.and.Btype(Ea)<=BC_MPI_BOUNDS_MAX)) then
                          Wnz = Prnz
   else
                          Wnz = Prnz-1
   end if
   
   
#ifdef PAR
   call par_init_exchange
      
   gUnx = par_co_sum(Unx, comm_row_x)
   gUny = par_co_sum(Uny, comm_row_y)
   gUnz = par_co_sum(Unz, comm_row_z)
   
   gVnx = par_co_sum(Vnx, comm_row_x)
   gVny = par_co_sum(Vny, comm_row_y)
   gVnz = par_co_sum(Vnz, comm_row_z)
   
   gWnx = par_co_sum(Wnx, comm_row_x)
   gWny = par_co_sum(Wny, comm_row_y)
   gWnz = par_co_sum(Wnz, comm_row_z)
#else
   gUnx = Unx
   gUny = Uny
   gUnz = Unz
   
   gVnx = Vnx
   gVny = Vny
   gVnz = Vnz
   
   gWnx = Wnx
   gWny = Wny
   gWnz = Wnz 
#endif


#ifdef PAR
   call par_init_domain_boundary_conditions
#endif
   

#ifdef CUSTOM_CONFIG
   call CustomConfiguration_Last
#endif

   if (master) write(*,*) "set"
   
   
   

 contains
 
 
 

     subroutine chget1(x)
       character(*),intent(out) :: x
       read(unit,fmt='(/)')
       read(unit,*) x
     end subroutine
     subroutine lget1(x)
       logical,intent(out) :: x
       character(120) :: line, fname
       integer :: ierr, tmp
       read(unit,fmt='(/)')
       read(unit,'(a)') line
       read(line,*,iostat=ierr) x
       if (ierr/=0) then
         read(line,*,iostat=ierr) tmp
         if (ierr/=0) then
           inquire(unit,name=fname)
           if (master) write(*,*) "Stop expected a boolean flag in file "//trim(fname)
           if (master) write(*,*) "Received '"//trim(line)//"' instead."
           call error_stop
         end if
         x = tmp /=0
       end if
     end subroutine
     subroutine lget2(x,y)
       logical,intent(out) :: x,y
       read(unit,fmt='(/)')
       read(unit,*) x,y
     end subroutine
     subroutine lget3(x,y,z)
       logical,intent(out) :: x,y,z
       read(unit,fmt='(/)')
       read(unit,*) x,y,z
     end subroutine
     subroutine iget1(x)
       integer,intent(out) :: x
       read(unit,fmt='(/)')
       read(unit,*) x
     end subroutine
     subroutine iget2(x,y)
       integer,intent(out) :: x,y
       read(unit,fmt='(/)')
       read(unit,*) x,y
     end subroutine
     subroutine iget3(x,y,z)
       integer,intent(out) :: x,y,z
       read(unit,fmt='(/)')
       read(unit,*) x,y,z
     end subroutine
     subroutine rget1(x)
       real(knd),intent(out) :: x
       read(unit,fmt='(/)')
       read(unit,*) x
     end subroutine
     subroutine rget2(x,y)
       real(knd),intent(out) :: x,y
       read(unit,fmt='(/)')
       read(unit,*) x,y
     end subroutine
     subroutine rget3(x,y,z)
       real(knd),intent(out) :: x,y,z
       read(unit,fmt='(/)')
       read(unit,*) x,y,z
     end subroutine
     subroutine rgetv3(v)
       real(knd),intent(out) :: v(3)
       read(unit,fmt='(/)')
       read(unit,*) v
     end subroutine

     subroutine read_cmd_conf
        command_line = ""
        open(unit,file="cmd.conf",status="old",action="read",iostat = io)
        if (io==0) then
          read(unit,'(a)') command_line
        end if
        command_line = "&cmd "//adjustl(trim(command_line))//" /"
     end subroutine

     subroutine read_command_line
       command_line = ""
       call get_command(command = command_line,status = io)
       if (io==0) then
         call get_command_argument(0,length = exenamelength,status = io2)
         if (io2==0) then
           command_line = "&cmd "//adjustl(trim(command_line(exenamelength+1:)))//" /"
         else
           command_line = "&cmd "//adjustl(trim(command_line))//" /"
         end if
       else
         write(*,*) io,"Error getting command line."
       end if
     end subroutine

     subroutine parse_command_line
#ifdef PAR
       use custom_par
#endif
       namelist /cmd/ tilesize, debugparam, debuglevel, windangle, &
                       Prnx, Prny, Prnz,&
#ifdef PAR
                       npxyz, domain_index, number_of_domains, &
#endif
                       obstacles_file, probes_file, scalar_probes_file, scratch_dir, &
                       enable_fixed_flow_rate

       if (len_trim(command_line)>0) then
         msg = ''
         read(command_line,nml = cmd,iostat = io,iomsg = msg)
         if (io/=0) then
           if (master) write(*,*) io,"Error parsing command line."
           if (master) write(*,*) msg
           if (master) write(*,*) command_line
         end if
       else
         if (master) write(*,*) io,"Error getting command line."
       end if
     end subroutine


     subroutine get_line_sources
        use Strings, only: itoa
        use LineSources, only: ScalarLineSource, ScalarLineSources
        type(ScalarLineSource) :: src
        integer :: n

        call get(n)
        allocate(ScalarLineSources(0))
        do i = 1, n
          read(unit,fmt=*)
          call get(src%scalar_number)
          if (src%scalar_number<0) call error_stop("Error: Scalar number of line source "//itoa(i)//" negative.")
          if (src%scalar_number>num_of_scalars) call error_stop("Error: Scalar number of line source "//itoa(i)//" too large.")
          
          call get(src%start)
          call get(src%end)
          call get(src%flux)
          src%number_of_points = 20 * max(Prnx,Prny,Prnz)
          ScalarLineSources = [ScalarLineSources, src]
        end do
     end subroutine

     subroutine get_point_sources
        use Strings, only: itoa
        use PointSources, only: ScalarPointSource, ScalarPointSources
        type(ScalarPointSource) :: src
        integer :: n

        call get(n)
        allocate(ScalarPointSources(0))
        do i = 1, n
          read(unit,fmt=*)
          call get(src%scalar_number)
          if (src%scalar_number<0) call error_stop("Error: Scalar number of point source "//itoa(i)//" negative.")
          if (src%scalar_number>num_of_scalars) call error_stop("Error: Scalar number of point source "//itoa(i)//" too large.")
          
          call get(src%position)
          call get(src%flux)
          ScalarPointSources = [ScalarPointSources, src]
        end do
     end subroutine

     subroutine read_namelist_output
       namelist /output/ store, display

       read(unit,nml = output,iostat = io,iomsg = msg)
     end subroutine
     
     subroutine read_frames  
       open(unit,file="frames.conf",status="old",action="read",iostat = io)
       if (io==0) then
       
         call get(frame_times%nframes)
         
         if (frame_times%nframes>0) then
           call get(frame_times%start)
           call get(frame_times%end)
           read(unit,fmt='(/)')

           read(unit,*) frame_flags%U
           call get(frame_flags%vorticity)
           call get(frame_flags%Pr)
           call get(frame_flags%lambda2)
           call get(frame_flags%scalars)
           if (num_of_scalars < 1) frame_flags%scalars = 0
           call get(frame_flags%sumscalars)
           if (num_of_scalars < 1) frame_flags%sumscalars = 0
           call get(frame_flags%temperature)
           if (.not.enable_buoyancy) frame_flags%temperature = 0
           call get(frame_flags%moisture)
           if (.not.enable_buoyancy) frame_flags%moisture = 0
           call get(frame_flags%temperature_flux)
           if (.not.enable_buoyancy) frame_flags%temperature_flux = 0
           call get(frame_flags%scalar_flux)
           if (num_of_scalars < 1) frame_flags%scalar_flux = 0

           call get(numframeslices)

           do i = 1, numframeslices
             call get(dimension)
             call get(direction)
             call get(position)
             call AddDomain(TFrameDomain(achar(iachar('a')+i-1), &
                            dimension, direction, position, &
                            frame_times, frame_flags))
           end do
         end if

         close(unit)
       else
         frames = 0
         if (master) write (*,*) "frames.conf not found, no vtk frames will be saved."
       end if
     end subroutine read_frames

     subroutine read_staggered_frames
       open(unit,file="stagframes.conf",status="old",action="read",iostat = io)
       if (io==0) then
         call get(num_staggered_domains)

         do i = 1, num_staggered_domains
           read(unit,fmt=*)
           call get(domain_label)
           call get(range%min%x,range%max%x)
           call get(range%min%y,range%max%y)
           call get(range%min%z,range%max%z)
           call get(frame_times%nframes)
           call get(frame_times%start, frame_times%end)
           call get(frame_save_flags%U, frame_save_flags%V, frame_save_flags%W)
           call get(frame_save_flags%Pr)
           call get(frame_save_flags%Viscosity)
           call get(frame_save_flags%Temperature)
           if (.not.enable_buoyancy) frame_save_flags%Temperature = .false.
           call get(frame_save_flags%Moisture)
           if (.not.enable_moisture) frame_save_flags%Moisture = .false.
           call get(frame_save_flags%Scalar)
           if (num_of_scalars < 1) frame_save_flags%Scalar = .false.

           call AddDomain(TStaggeredFrameDomain(trim(domain_label), &
                                                range, &
                                                frame_times, &
                                                frame_save_flags))
         end do
         close(unit)
       else
         if (master) write (*,*) "stagframes.conf not found, no staggered frames will be saved."
       end if
     end subroutine read_staggered_frames

  end subroutine ReadConfiguration






  subroutine get_area_sources(fname)
    use Strings
    use ParseTrees
    use GeometricShapes2D
    use AreaSources
    character(*), intent(in) :: fname
    type(tree_object), allocatable :: tree(:)
    logical :: ex
    integer :: i, stat

    if (num_of_scalars==0) return

    inquire(file=fname, exist=ex)

    if (.not.ex) return

    call parse_file(tree, fname, stat)

    if (stat==0) then

      if (allocated(tree)) then
        do i = 1, size(tree)
          call get_area_source(tree(i))
          call tree(i)%finalize
        end do
      end if

    else

      write(*,*) "Error parsing file " // fname
      call error_stop

    end if

  contains

    subroutine get_area_source(obj)
      type(tree_object), intent(in) :: obj
      class(GeometricShape2D), allocatable :: shp
      type(ScalarAreaSource) :: src
      real(knd) :: flux
      integer :: scnum
      integer :: j

      if (downcase(obj%name)=='area_source') then

        if (allocated(obj%fields%array)) then

          associate(fields => obj%fields%array)

            do j = 1, size(fields)

              if (downcase(fields(j)%name)=='scalar_number') then
                read(fields(j)%value, *) scnum
              else if (downcase(fields(j)%name)=='flux') then
                read(fields(j)%value, *) flux
              else if (downcase(fields(j)%name)=='geometric_shape') then
                if (fields(j)%is_object .and. associated(fields(j)%object_value)) then
                  call get_geometric_shape(shp, fields(j)%object_value)
                else 
                  write(*,*) "Invalid geometric shape for area source in " // fname
                  call error_stop
                end if
              end if

            end do

          end associate

        else

          write(*,*) "No fields in the AreaSource object in " // fname
          call error_stop

        end if

        src = ScalarAreaSource(shp, scnum, flux)
        call add_element(ScalarAreaSources, src)

      else
        write(*,*) "Unknown object type " // downcase(obj%name) // " in " // fname
        call error_stop
      end if
    end subroutine

    subroutine get_geometric_shape(res, obj)
      class(GeometricShape2D), allocatable :: res
      type(tree_object), intent(in) :: obj
      real(knd) :: xc, yc, r
      integer :: j

      if (downcase(obj%name)=='circle') then

        if (allocated(obj%fields%array)) then

          associate(fields => obj%fields%array)

            do j = 1, size(fields)

              if (downcase(fields(j)%name)=='xc') then
                read(fields(j)%value, *) xc
              else if (downcase(fields(j)%name)=='yc') then
                read(fields(j)%value, *) yc
              else if (downcase(fields(j)%name)=='r') then
                read(fields(j)%value, *) r
              end if

            end do

          end associate

        else

          write(*,*) "No fields in the Circle object in " // fname
          call error_stop

        end if

        allocate(res, source = Circle(xc, yc, r))

      else

        write(*,*) "Invalid geometric shape for area source in " // fname
        write(*,*) "Supported variants: Circle"
        write(*,*) "Received:", obj%name

        call error_stop

      end if

    end subroutine

    subroutine add_element(a,e)
      type(ScalarAreaSource), allocatable, intent(inout) :: a(:)
      type(ScalarAreaSource), intent(inout) :: e
      type(ScalarAreaSource), allocatable :: tmp(:)
      integer :: i, n

      if (.not.allocated(a)) then
        a = [e]
      else
        n = size(a)
        call move_alloc(a,tmp)
        allocate(a(n+1))

        do i = 1, n
          call assign(a(i), tmp(i))
        end do
        call assign(a(n+1), e)
      end if
    end subroutine
    
    subroutine assign(l, r)
      type(ScalarAreaSource), intent(out) :: l
      type(ScalarAreaSource), intent(inout) :: r
      l%flux = r%flux
      l%scalar_number = r%scalar_number
      call move_alloc(r%GeometricShape, l%GeometricShape)
    end subroutine

  end subroutine get_area_sources


  subroutine get_geostrophic_wind(fname, g)
    use Interpolation
    character(*), intent(in) :: fname
    type(spline_coefs), intent(out) :: g
    character(256) :: line
    real(knd) :: r3(3)
    real(knd), allocatable :: ug(:), vg(:)
    integer :: unit, io, n, i, j

    open(newunit=unit,file=fname,status="old",action="read",iostat = io)
    if (io/=0) return

    n = 0
    do
      read(unit,'(a)',iostat=io) line
      if (io/=0) exit
      read(line,*,iostat=io) r3
      if (io/=0) exit
      n = n + 1
    end do
    
    rewind(unit)
    
    if (n>0) then
      allocate(g%z(n), ug(n), vg(n))
      allocate(g%cu(0:1,n), g%cv(0:1,n))
      do i = 1, n
        read(unit,'(a)',iostat=io) line
        read(line, *) g%z(i), ug(i), vg(i)
      end do
    else
      stop "Geostrophic profile empty."
    end if

    if (n > 1) then
      call linear_interpolation(g%z, ug, g%cu)
      call linear_interpolation(g%z, vg, g%cv)

      allocate(pr_gradient_profile_x(1:Prnz))
      allocate(pr_gradient_profile_y(1:Prnz))

      j = 1
      do i = 1, Prnz
        pr_gradient_profile_x(i) =   Coriolis_parameter * linear_interpolation_eval(zPr(i), g%z, g%cv, j)
      end do

      j = 1
      do i = 1, Prnz
        pr_gradient_profile_y(i) = - Coriolis_parameter * linear_interpolation_eval(zPr(i), g%z, g%cu, j)
      end do

      enable_pr_gradient_x_profile = any(pr_gradient_profile_x/=0)
      enable_pr_gradient_y_profile = any(pr_gradient_profile_y/=0)

      if (enable_pr_gradient_x_profile) enable_pr_gradient_x_uniform = .false.
      if (enable_pr_gradient_y_profile) enable_pr_gradient_y_uniform = .false.

    else

      pr_gradient_x =   Coriolis_parameter * vg(1)
      pr_gradient_y = - Coriolis_parameter * ug(1)

      enable_pr_gradient_x_uniform = pr_gradient_x /= 0
      enable_pr_gradient_y_uniform = pr_gradient_y /= 0

    end if
    

  end subroutine get_geostrophic_wind


  subroutine get_profiles(fname)
    use Strings
    use ParseTrees
    use Outputs, only: profiles_config, enable_profiles
    
    character(*), intent(in) :: fname
    type(tree_object), allocatable :: tree(:)
    logical :: ex
    integer :: i, stat

    inquire(file=fname, exist=ex)

    if (.not.ex) return

    call parse_file(tree, fname, stat)

    if (stat==0) then

      if (allocated(tree)) then
        do i = 1, size(tree)
          call get_profile(tree(i))
          call tree(i)%finalize
        end do
      end if

    else

      write(*,*) "Error parsing file " // fname
      call error_stop

    end if
    
  contains
  
    subroutine get_profile(obj)
      type(tree_object), intent(in) :: obj
      integer :: j
      real(knd) :: rval

      if (downcase(obj%name)=='average_profiles') then

        if (allocated(obj%fields%array)) then

          associate(fields => obj%fields%array)

            do j = 1, size(fields)
              read(fields(j)%value, *) rval
             
              select case (downcase(fields(j)%name))
                case ('start')
                  profiles_config%average_start = rval
                case ('end')
                  profiles_config%average_end = rval
              end select

            end do

          end associate

        else

          write(*,*) "No fields in the profile object in " // fname
          call error_stop

        end if

        enable_profiles = .true.

      else if (downcase(obj%name)=='instantaneous_profiles') then

        if (allocated(obj%fields%array)) then

          associate(fields => obj%fields%array)

            do j = 1, size(fields)

              read(fields(j)%value, *) rval
              
              select case (downcase(fields(j)%name))
                case ('start')
                  profiles_config%instant_start = rval
                case ('end')
                  profiles_config%instant_end = rval
                case ('interval')
                  profiles_config%instant_interval = rval
              end select

            end do

          end associate

        else

          write(*,*) "No fields in the profile object in " // fname
          call error_stop

        end if

        enable_profiles = .true.

      else if (downcase(obj%name)=='running_average_profiles') then

        if (allocated(obj%fields%array)) then

          associate(fields => obj%fields%array)

            do j = 1, size(fields)

              read(fields(j)%value, *) rval
              
              select case (downcase(fields(j)%name))
                case ('start')
                  profiles_config%running_start = rval
                case ('end')
                  profiles_config%running_end = rval
                case ('interval')
                  profiles_config%running_interval = rval
              end select

            end do

          end associate
          
        else

          write(*,*) "No fields in the profile object in " // fname
          call error_stop

        end if

        enable_profiles = .true.

      else
        write(*,*) "Unknown object type " // downcase(obj%name) // " in " // fname
        call error_stop
      end if
    end subroutine

    
  end subroutine get_profiles


  subroutine get_time_stepping(fname)
    use Strings
    use ParseTrees
    character(*), intent(in) :: fname
    type(tree_object), allocatable :: tree(:)
    logical :: ex
    integer :: obj, stat

    type field_names
      character(char_len) :: name
      class(*), pointer :: var
    end type

    type field_names_a
      character(char_len) :: name
      class(*), pointer :: var(:)
    end type

    interface field_names
      procedure field_names_init
    end interface

    interface field_names_a
      procedure field_names_a_init
    end interface

    inquire(file=fname, exist=ex)

    if (.not.ex) return

    call parse_file(tree, fname, stat)

    if (stat==0) then

      if (allocated(tree)) then
        do obj = 1, size(tree)
          call get_object(tree(obj), time_stepping)
          call tree(obj)%finalize
        end do
      else
        write(*,*) "Error, no content in " // fname
        call error_stop
      end if

    else

      write(*,*) "Error parsing file " // fname
      call error_stop

    end if

  contains

    subroutine get_object(obj, t_s)
      type(tree_object), intent(in) :: obj
      type(time_step_control), intent(out), target :: t_s
      integer :: i, j
      logical, target :: constant_time_steps = .false.

      type(field_names) :: names(12)
      type(field_names_a) :: names_a(3)

      names = [field_names_init("max_number_of_time_steps",   t_s%max_number_of_time_steps), &
               field_names_init("variable_time_steps",        t_s%variable_time_steps), &
               field_names_init("constant_time_steps",        constant_time_steps), &
               field_names_init("enable_U_scaling",           t_s%enable_U_scaling), &
               field_names_init("enable_CFL_check",           t_s%enable_CFL_check), &
               field_names_init("dt",                         t_s%dt_constant), &
               field_names_init("dt_max",                     t_s%dt_max), &
               field_names_init("dt_min",                     t_s%dt_min), &
               field_names_init("CFL",                        t_s%CFL), &
               field_names_init("CFL_max",                    t_s%CFL_max), &
               field_names_init("start_time",                 t_s%start_time), &
               field_names_init("end_time",                   t_s%end_time)]

      names_a = [field_names_a_init("U_scaling", t_s%U_scaling), &
                 field_names_a_init("U_max",     t_s%U_max), &
                 field_names_a_init("U_min",     t_s%U_min)]

      if (downcase(obj%name)=='time_stepping') then

        if (allocated(obj%fields%array)) then

          associate(fields => obj%fields%array)

fields_do:  do j = 1, size(fields)
             
              do i = 1, size(names)
                if (fields(j)%name == names(i)%name) then
                  select type (var => names(i)%var)
                    type is (integer)
                      read(fields(j)%value, *) var
                    type is (real(real32))
                      read(fields(j)%value, *) var
                    type is (real(real64))
                      read(fields(j)%value, *) var
                    type is (logical)
                      read(fields(j)%value, *) var
                    class default
                      call error_stop("Unexpected type in time_stepping.")
                  end select
                  cycle fields_do
                end if
              end do

              do i = 1, size(names_a)
                if (fields(j)%name == names_a(i)%name) then

                  if (size(names_a(i)%var)/=size(fields(j)%array_value)) then
                    write(*,*) "Error, expecting", &
                                size(names_a(i)%var), &
                                "vector components", &
                                "but", &
                                size(fields(j)%array_value), &
                                "components present."
                    call error_stop()
                  end if

                  select type (var => names_a(i)%var)
                    type is (integer)
                      read(fields(j)%array_value, *) var
                    type is (real(real32))
                      read(fields(j)%array_value, *) var
                    type is (real(real64))
                      read(fields(j)%array_value, *) var
                    type is (logical)
                      read(fields(j)%array_value, *) var
                    class default
                      call error_stop("Unexpected type in time_stepping.")
                  end select
                  cycle fields_do
                end if
              end do

            end do fields_do

          end associate

        end if

      end if

      if (constant_time_steps) t_s%variable_time_steps = .false.

      if (t_s%variable_time_steps) then

        t_s%U_max = abs(t_s%U_max)
        t_s%U_min = abs(t_s%U_min)

        if (maxval(t_s%U_max) <=  0) &
          call error_stop("Error, time_stepping%U_max must have at least one non-zero component.")

        if (maxval(t_s%U_min) <=  0) &
          call error_stop("Error, time_stepping%U_min must have at least one non-zero component.")

        if (t_s%CFL <= 0) &
          call error_stop("Error, time_stepping%CFL must be positive.")


        t_s%dt_min = t_s%CFL / (t_s%U_max(1) / dxmin + &
                                t_s%U_max(2) / dymin + &
                                t_s%U_max(3) / dzmin)
        
        t_s%dt_max = t_s%CFL / (t_s%U_min(1) / dxmin + &
                                t_s%U_min(2) / dymin + &
                                t_s%U_min(3) / dzmin)

        t_s%dt = t_s%dt_min

      else

        t_s%U_scaling = abs(t_s%U_scaling)

        
        if (maxval(t_s%U_scaling) > 0) then

          t_s%enable_U_scaling = .true.

          if (t_s%CFL<=0) &
            call error_stop("Error, time_stepping%CFL must be positive when using U_scaling.")

          t_s%dt_constant = t_s%CFL / (t_s%U_scaling(1) / dxmin + &
                                       t_s%U_scaling(2) / dymin + &
                                       t_s%U_scaling(3) / dzmin)
        else
          if (t_s%dt_constant <=  0) &
            call error_stop("Error, time_stepping%dt or t_s%U_scaling must be positive.")
        end if

        t_s%enable_CFL_check = (t_s%CFL_max > 0)

        t_s%dt = t_s%dt_constant

      end if
      
    
    end subroutine

    function field_names_init(name, var) result(res)
      type(field_names) :: res
      character(*) :: name
      class(*), target, intent(in) :: var

      res%name = name
      res%var => var
    end function

    function field_names_a_init(name, var) result(res)
      type(field_names_a) :: res
      character(*) :: name
      class(*), target, intent(in) :: var(:)

      res%name = name
      res%var => var
    end function

  end subroutine get_time_stepping
  
  
  
  

  subroutine get_pressure_solution(fname)
    use Strings
    use ParseTrees
    use Pressure
    character(*), intent(in) :: fname
    type(tree_object), allocatable :: tree(:)
    logical :: ex
    integer :: obj, stat

    type field_names
      character(char_len) :: name
      class(*), pointer :: var
    end type

    type field_names_a
      character(char_len) :: name
      class(*), pointer :: var(:)
    end type

    interface field_names
      procedure field_names_init
    end interface

    interface field_names_a
      procedure field_names_a_init
    end interface

    inquire(file=fname, exist=ex)

    if (.not.ex) return

    call parse_file(tree, fname, stat)

    if (stat==0) then

      if (allocated(tree)) then
        do obj = 1, size(tree)
          call get_object(tree(obj), pressure_solution)
          call tree(obj)%finalize
        end do
      else
        write(*,*) "Error, no content in " // fname
        call error_stop
      end if

    else

      write(*,*) "Error parsing file " // fname
      call error_stop

    end if

  contains

    subroutine get_object(obj, p_s)
      type(tree_object), intent(in) :: obj
      type(pressure_solution_control), intent(out), target :: p_s
      integer :: i, j
      logical, target :: constant_time_steps = .false.

      type(field_names) :: names(10)
      type(field_names_a) :: names_a(1)

      names = [field_names_init("check_mass_flux", p_s%check_mass_flux), &
               field_names_init("correct_mass_flux_west",   p_s%correct_mass_flux(We)), &
               field_names_init("correct_mass_flux_east",   p_s%correct_mass_flux(Ea)), &
               field_names_init("correct_mass_flux_south",  p_s%correct_mass_flux(So)), &
               field_names_init("correct_mass_flux_north",  p_s%correct_mass_flux(No)), &
               field_names_init("correct_mass_flux_bottom", p_s%correct_mass_flux(Bo)), &
               field_names_init("correct_mass_flux_top",    p_s%correct_mass_flux(To)), &
               field_names_init("poisson_solver",      p_s%poisson_solver), &
               field_names_init("check_divergence",    p_s%check_divergence), &
               field_names_init("bottom_pressure",     p_s%bottom_pressure)]

      names_a = [field_names_a_init("correct_mass_flux", p_s%correct_mass_flux)]

      if (downcase(obj%name)=='pressure_solution') then

        if (allocated(obj%fields%array)) then

          associate(fields => obj%fields%array)

fields_do:  do j = 1, size(fields)
             
              do i = 1, size(names)
                if (fields(j)%name == names(i)%name) then
                  select type (var => names(i)%var)
                    type is (integer)
                      read(fields(j)%value, *) var
                    type is (real(real32))
                      read(fields(j)%value, *) var
                    type is (real(real64))
                      read(fields(j)%value, *) var
                    type is (logical)
                      read(fields(j)%value, *) var
                    class default
                      call error_stop("Unexpected type in pressure_solution.")
                  end select
                  cycle fields_do
                end if
              end do

              do i = 1, size(names_a)
                if (fields(j)%name == names_a(i)%name) then

                  if (size(names_a(i)%var)/=size(fields(j)%array_value)) then
                    write(*,*) "Error, expecting", &
                                size(names_a(i)%var), &
                                "array components", &
                                "but", &
                                size(fields(j)%array_value), &
                                "components present."
                    call error_stop()
                  end if

                  select type (var => names_a(i)%var)
                    type is (integer)
                      read(fields(j)%array_value, *) var
                    type is (real(real32))
                      read(fields(j)%array_value, *) var
                    type is (real(real64))
                      read(fields(j)%array_value, *) var
                    type is (logical)
                      read(fields(j)%array_value, *) var
                    class default
                      call error_stop("Unexpected type in pressure_solution.")
                  end select
                  cycle fields_do
                end if
              end do

            end do fields_do

          end associate

        end if

      end if

      if (any(p_s%correct_mass_flux)) p_s%check_mass_flux = .true.

      where(Btype>=BC_MPI_BOUNDS_MIN .and. Btype<=BC_MPI_BOUNDS_MAX) p_s%correct_mass_flux = .false.
    
      where(Btype==BC_NOSLIP) p_s%correct_mass_flux = .false.

      where(Btype==BC_PERIODIC) p_s%correct_mass_flux = .false.
    
    end subroutine

    function field_names_init(name, var) result(res)
      type(field_names) :: res
      character(*) :: name
      class(*), target, intent(in) :: var

      res%name = name
      res%var => var
    end function

    function field_names_a_init(name, var) result(res)
      type(field_names_a) :: res
      character(*) :: name
      class(*), target, intent(in) :: var(:)

      res%name = name
      res%var => var
    end function

  end subroutine get_pressure_solution
  
  
  
  
#ifdef PAR
  subroutine get_domains(fname)
    use Strings
    use ParseTrees
    use custom_par
    character(*), intent(in) :: fname
    type(tree_object), allocatable :: tree(:)
    logical :: ex
    integer :: obj, stat

    type field_names
      character(char_len) :: name
      class(*), pointer :: var
    end type

    interface field_names
      procedure field_names_init
    end interface

    inquire(file=fname, exist=ex)

    if (.not.ex) return

    call parse_file(tree, fname, stat)

    if (stat==0) then

      if (allocated(tree)) then
        do obj = 1, size(tree)
          call get_object(tree(obj))
          call tree(obj)%finalize
        end do
      else
        write(*,*) "Error, no content in " // fname
        call error_stop
      end if

    else

      write(*,*) "Error parsing file " // fname
      call error_stop

    end if

  contains

    subroutine get_object(obj)
      type(tree_object), intent(in) :: obj
      integer :: i, j
      logical, target :: constant_time_steps = .false.

      type(field_names) :: names(5)
      
      logical, target :: enable_multiple_domains_l = .false.
      integer, target :: number_of_domains_l = -99
      integer, target :: domain_index_l = -99
      integer, target :: spatial_ratio_l = -99
      integer, target :: time_step_ratio_l = -99

      names = [field_names_init("enable_multiple_domains", enable_multiple_domains_l), &
               field_names_init("domain_index",   domain_index_l), &
               field_names_init("number_of_domains",   number_of_domains_l), &
               field_names_init("spatial_ratio",   spatial_ratio_l), &
               field_names_init("time_step_ratio",   time_step_ratio_l)]

      if (downcase(obj%name)=='domains') then

        if (allocated(obj%fields%array)) then

          associate(fields => obj%fields%array)

fields_do:  do j = 1, size(fields)
             
              do i = 1, size(names)
                if (fields(j)%name == names(i)%name) then
                  select type (var => names(i)%var)
                    type is (integer)
                      read(fields(j)%value, *) var
                    type is (real(real32))
                      read(fields(j)%value, *) var
                    type is (real(real64))
                      read(fields(j)%value, *) var
                    type is (logical)
                      read(fields(j)%value, *) var
                    class default
                      call error_stop("Unexpected type in domains.")
                  end select
                  cycle fields_do
                end if
              end do

            end do fields_do

          end associate

        end if

      end if

      enable_multiple_domains = enable_multiple_domains_l
      
      if (enable_multiple_domains) then
          
        number_of_domains = number_of_domains_l
        if (number_of_domains<1) &
          call error_stop("Error, positive number_of_domains must be specified in domains.")
    
        domain_index = domain_index_l
        if (domain_index<1) &
          call error_stop("Error, positive domain_index must be specified in domains.")
        if (number_of_domains<domain_index) &
          call error_stop("Error, domain_index must be smaller or equal to number_of_domains.")

        if (spatial_ratio_l>0) domain_spatial_ratio = spatial_ratio_l
        if (time_step_ratio_l>0) domain_time_step_ratio = time_step_ratio_l
      end if
      
    end subroutine

    function field_names_init(name, var) result(res)
      type(field_names) :: res
      character(*) :: name
      class(*), target, intent(in) :: var

      res%name = name
      res%var => var
    end function

  end subroutine get_domains
#endif  
  
  
  

  subroutine ReadInitialConditions(U,V,W,Pr,Temperature,Moisture,Scalar,scalars_optional)
    use Endianness
    real(knd), intent(inout) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
    real(knd), intent(inout) :: Pr(1:,1:,1:)
    real(knd), intent(inout) :: Temperature(-1:,-1:,-1:)
    real(knd), intent(inout) :: Moisture(-1:,-1:,-1:)
    real(knd), intent(inout) :: Scalar(-1:,-1:,-1:,:)
    logical, intent(in) :: scalars_optional
    real(real32), allocatable :: buffer(:,:,:), UVWbuffer(:,:,:,:)
    integer :: i, unit, stat, file_stat
    logical :: exU(3)
    character(2) :: scalnum

    allocate(buffer(Prnx,Prny,Prnz))

    open(newunit=unit,file=image_input_dir//"out.vtk",access="stream",status="old",action="read",iostat=file_stat)

    if (file_stat/=0) call error_stop("Error opening "//image_input_dir//"out.vtk")

    call skip_to("SCALARS p float",stat)
    if (stat/=0) then
      Pr = 0
      rewind(unit)
    else
      call skip_line
      call skip_line
      read(unit) buffer
      Pr(1:Prnx,1:Prny,1:Prnz) =  real(BigEnd(buffer),knd)
    end if

    if (enable_buoyancy) then
      call skip_to("SCALARS temperature float",stat)
      if (stat/=0) then
        call error_stop("No temperature field found in the initial conditions file " // &
                        image_input_dir // "out.vtk")
      else
        call skip_line
        call skip_line
        read(unit) buffer
        Temperature(1:Prnx,1:Prny,1:Prnz) =  real(BigEnd(buffer),knd)
      end if
    end if

    if (enable_moisture) then
      call skip_to("SCALARS moisture float",stat)
      if (stat/=0) then
        call error_stop("No moisture field found in the initial conditions file " // &
                        image_input_dir // "out.vtk")
      else
        call skip_line
        call skip_line
        read(unit) buffer
        Moisture(1:Prnx,1:Prny,1:Prnz) =  real(BigEnd(buffer),knd)
      end if
    end if

    close(unit)


    if (num_of_scalars > 0) then
      open(newunit=unit,file=image_input_dir//"scalars.vtk",access="stream",status="old",action="read",iostat=file_stat)

      if (file_stat/=0 .and. scalars_optional) then
        !TODO: check consistency between images?
        if (master) write(*,*) "Warning, no scalar initial conditions found."
        if (master) write(*,*) "Scalar fields set to zero."
        
        Scalar = 0

      else if (file_stat/=0) then

        call error_stop("Error opening "//image_input_dir//"scalars.vtk")

      else

        do i = 1, num_of_scalars
          write(scalnum,"(I2.2)") i
          call skip_to("SCALARS scalar"//scalnum//" float",stat)
          call skip_line
          if (stat/=0) then
            call error_stop("scalar"//scalnum//" field not found in the initial conditions file " // &
                            image_input_dir // "scalars.vtk")
          else
            call skip_line
            read(unit) buffer
            Scalar(1:Prnx,1:Prny,1:Prnz,i) =  real(BigEnd(buffer),knd)
          end if
        end do

        close(unit)

      end if

    end if


    deallocate(buffer)


    inquire(file=image_input_dir//"U.vtk", exist=exU(1))
    inquire(file=image_input_dir//"V.vtk", exist=exU(2))
    inquire(file=image_input_dir//"W.vtk", exist=exU(3))
    
    if (all(exU)) then
    
      open(newunit=unit,file=image_input_dir//"U.vtk",access="stream",status="old",action="read")
      call skip_to("SCALARS U float",stat)
      call skip_line
      call skip_line
      allocate(buffer(1:Unx, 1:Uny, 1:Unz))
      read(unit) buffer
      U(1:Unx, 1:Uny, 1:Unz) = real(BigEnd(buffer), knd)
      deallocate(buffer)
      close(unit)

      open(newunit=unit,file=image_input_dir//"V.vtk",access="stream",status="old",action="read")
      call skip_to("SCALARS V float",stat)
      call skip_line
      call skip_line
      allocate(buffer(1:Vnx, 1:Vny, 1:Vnz))
      read(unit) buffer
      V(1:Vnx, 1:Vny, 1:Vnz) = real(BigEnd(buffer), knd)
      deallocate(buffer)
      close(unit)

      open(newunit=unit,file=image_input_dir//"W.vtk",access="stream",status="old",action="read")
      call skip_to("SCALARS W float",stat)
      call skip_line
      call skip_line
      allocate(buffer(1:Wnx, 1:Wny, 1:Wnz))
      read(unit) buffer
      W(1:Wnx, 1:Wny, 1:Wnz) = real(BigEnd(buffer), knd)
      deallocate(buffer)
      close(unit)
      
    else
    
      open(newunit=unit,file=image_input_dir//"out.vtk",access="stream",status="old",action="read")
      call skip_to("VECTORS u float",stat)
      call skip_line
      allocate(UVWbuffer(3, 1:Prnx, 1:Prny, 1:Prnz))
      read(unit) UVWbuffer
      UVWbuffer = BigEnd(UVWbuffer)
      close(unit)
      
      U = 0
      U(1:Prnx,1:Prny,1:Prnz) = real(UVWbuffer(1,1:Prnx,1:Prny,1:Prnz),knd)
      U(0:Prnx-1,1:Prny,1:Prnz) = U(0:Prnx-1,1:Prny,1:Prnz) + U(1:Prnx,1:Prny,1:Prnz)
    
      V = 0
      V(1:Prnx,1:Prny,1:Prnz) = real(UVWbuffer(2,1:Prnx,1:Prny,1:Prnz),knd)
      V(1:Prnx,0:Prny-1,1:Prnz) = V(1:Prnx,0:Prny-1,1:Prnz) + V(1:Prnx,1:Prny,1:Prnz)
    
      W = 0
      W(1:Prnx,1:Prny,1:Prnz) = real(UVWbuffer(3,1:Prnx,1:Prny,1:Prnz),knd)
      W(0:Prnx-1,1:Prny,1:Prnz) = W(0:Prnx-1,1:Prny,1:Prnz) + W(1:Prnx,1:Prny,1:Prnz)

#ifdef PAR
      call par_exchange_U_x(U, 1)
      call par_exchange_U_y(V, 2)
      call par_exchange_U_z(W, 3)
#endif
      if ((Btype(Ea)>=BC_MPI_BOUNDS_MIN .and. (Btype(Ea)<=BC_MPI_BOUNDS_MAX)) .or. &
          Btype(Ea)==BC_PERIODIC) &
        U(Prnx,1:Prny,1:Prnz) = U(Prnx,1:Prny,1:Prnz) + U(0,1:Prny,1:Prnz)
      if ((Btype(No)>=BC_MPI_BOUNDS_MIN .and. (Btype(No)<=BC_MPI_BOUNDS_MAX)) .or. &
          Btype(No)==BC_PERIODIC) &
        V(1:Prnx,Prny,1:Prnz) = V(1:Prnx,Prny,1:Prnz) + V(1:Prnx,0,1:Prnz)
      if ((Btype(To)>=BC_MPI_BOUNDS_MIN .and. (Btype(To)<=BC_MPI_BOUNDS_MAX)) .or. &
          Btype(To)==BC_PERIODIC) &
        W(1:Prnx,1:Prny,Prnz) = W(1:Prnx,1:Prny,Prnz) + W(1:Prnx,1:Prny,0)

      U = U / 2
      V = V / 2
      W = W / 2
    
    end if


  contains

    subroutine skip_line
      character :: ch
      do
        read(unit) ch
        if (ch==new_line("a")) return
      end do
    end subroutine

    subroutine skip_to(str, stat)
      character(*), intent(in) :: str
      integer, intent(out) :: stat
      character :: ch
      integer :: io

      do
        read(unit, iostat=io) ch

        if (io/=0) then
          stat = 1
          return
        end if

        if (ch==str(1:1)) then
          call check(str(2:), stat)
          if (stat == 0) return
        end if

      end do
    end subroutine

    subroutine check(str, stat)
      character(*), intent(in) :: str
      integer, intent(out) :: stat
      character :: ch
      integer :: i, io

      stat = 1
      i = 0

      do
        i = i + 1

        read(unit, iostat=io) ch

        if (io/=0) return

        if (ch/=str(i:i)) return

        if (i==len(str)) then
          stat = 0
          return
        end if
      end do
    end subroutine

  end subroutine ReadInitialConditions



  subroutine InitialConditions(U,V,W,Pr,Temperature,Moisture,Scalar,dt)
    use custom_par
    use ArrayUtilities
    real(knd),contiguous,intent(inout) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
    real(knd),contiguous,intent(inout) :: Temperature(-1:,-1:,-1:)
    real(knd),contiguous,intent(inout) :: Moisture(-1:,-1:,-1:)
    real(knd),contiguous,intent(inout) :: Scalar(-1:,-1:,-1:,:)
    real(knd), intent(out) :: dt
    integer :: i,j,k
    real(knd) :: p,x,y,z,x1,x2,y1,y2,z1,z2
    real(knd),allocatable :: Q(:,:,:)

#ifdef CUSTOM_INITIAL_CONDITIONS
    interface
      subroutine CustomInitialConditions(U,V,W,Pr,Temperature,Moisture,Scalar)
        use Parameters
        real(knd),contiguous,intent(inout) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
        real(knd),contiguous,intent(inout) :: Temperature(-1:,-1:,-1:)
        real(knd),contiguous,intent(inout) :: Moisture(-1:,-1:,-1:)
        real(knd),contiguous,intent(inout) :: Scalar(-1:,-1:,-1:,:)
      end subroutine
    end interface
#endif

    call par_sync_out("  ...setting initial dummy values.")

    if (abs(Uinlet)>0) then
      dt = min(abs(dxmin/Uinlet), abs(dymin/Uinlet), abs(dzmin/Uinlet))
    else if (maxval(abs(time_stepping%U_min))>0) then
      dt = min(abs(dxmin/time_stepping%U_min(1)), &
               abs(dymin/time_stepping%U_min(2)), &
               abs(dzmin/time_stepping%U_min(3)))
    else
      dt = dxmin
    end if
        
    Pr(1:Prnx,1:Prny,1:Prnz) = 0

    U = huge(1._knd)/2
    V = huge(1._knd)/2
    W = huge(1._knd)/2

    V(1:Vnx,1:Vny,1:Vnz) = 0
    W(1:Wnx,1:Wny,1:Wnz) = 0

    if (initcondsfromfile>0) then

      call par_sync_out("  ...reading initial conditions from input files.")

      if (initcondsfromfile==2) then
        call ReadInitialConditions(U,V,W,Pr,Temperature,Moisture,Scalar,scalars_optional=.true.)
      else
        call ReadInitialConditions(U,V,W,Pr,Temperature,Moisture,Scalar,scalars_optional=.false.)
      end if

      Viscosity = molecular_viscosity

      if (enable_buoyancy.or. &
          enable_moisture.or. &
          num_of_scalars>0)        TDiff = molecular_diffusivity

       call BoundUVW(U, V, W)
       call Bound_Pr(Pr)

    else   !init conditions not from file

       call par_sync_out("  ...computing initial conditions.")

#ifdef CUSTOM_INITIAL_CONDITIONS
       call CustomInitialConditions(U,V,W,Pr,Temperature,Moisture,Scalar)
#else
       if (task_type==2) then
         U(1:Unx,1:Uny,1:Unz) = 0
         do k = 1, Unz
          do j = 1, Uny
           do i = 1, Unx
                  x = xU(i)
                  y = yPr(j)
                  z = zPr(k)
                  U(i,j,k) = -cos(pi*x)*sin(pi*y)
           end do
          end do
         end do
         do k = 1, Vnz
          do j = 1, Vny
           do i = 1, Vnx
                  x = xPr(i)
                  y = yV(j)
                  z = zPr(k)
                  V(i,j,k) = sin(pi*x)*cos(pi*y)
           end do
          end do
         end do
         do k = 1, Wnz
          do j = 1, Wny
           do i = 1, Wnx
                  x = xPr(i)
                  y = yPr(j)
                  z = zW(k)
                  W(i,j,k) = 0
           end do
          end do
         end do
         do k = 1, Prnz
          do j = 1, Prny
           do i = 1, Prnx
                  x = xPr(i)
                  y = yPr(j)
                  z = zPr(k)
                  Pr(i,j,k) = -(1._knd/4._knd)*((cos(2*pi*y)+cos(2*pi*x)))
           end do
          end do
         end do

       elseif (task_type==3) then
         U(1:Unx,1:Uny,1:Unz) = 0
         do k = 1, Unz
          do j = 1, Uny
           do i = 1, Unx

                  x = xU(i)
                  y = yPr(j)
                  z = zPr(k)
                  U(i,j,k) = Uinlet*sin(x)*cos(z)*cos(-y)
           end do
          end do
         end do
         do k = 1, Vnz
          do j = 1, Vny
           do i = 1, Vnx
                  x = xPr(i)
                  y = yV(j)
                  z = zPr(k)
                  V(i,j,k) = 0
           end do
          end do
         end do
         do k = 1, Wnz
          do j = 1, Wny
           do i = 1, Wnx
                  x = xPr(i)
                  y = yPr(j)
                  z = zW(k)
                  W(i,j,k) = -Uinlet*cos(x)*sin(z)*cos(-y)
           end do
          end do
         end do
         do k = 1, Prnz
          do j = 1, Prny
           do i = 1, Prnx
                  x = xPr(i)
                  y = yPr(j)
                  z = zPr(k)
                  Pr(i,j,k) = (Uinlet/16._knd)*((2+cos(2*z))*(cos(2*(-y))+cos(2*(x)))-2)
           end do
          end do
         end do

       elseif (InletType==TurbulentInletType) then

         call par_sync_out("  ...computing turbulent initial conditions.")

         dt = hypot(dxmin,dymin) / hypot(avg(Uin(1:Uny,1:Unz)),avg(Vin(1:Vny,1:Vnz)))

         do i = 1, Prnx
           
           call default_turbulence_generator%time_step(Uin, Vin, Win, time_stepping%dt)
           
           !$omp parallel private(j,k)
           !$omp do collapse(2)
           do k = 1, Unz
            do j = 1, Uny
              if (Utype(i,j,k)<=0) then
                 U(i,j,k) = Uin(j,k)
              else
                 U(i,j,k) = 0
              end if
            end do
           end do
           !$omp end do nowait
           !$omp do collapse(2)
           do k = 1, Vnz
            do j = 1, Vny
              if (Vtype(i,j,k)<=0) then
                 V(i,j,k) = Vin(j,k)
              else
                 V(i,j,k) = 0
              end if
            end do
           end do
           !$omp end do nowait
           !$omp do collapse(2)
           do k = 1, Wnz
            do j = 1, Wny
              if (Wtype(i,j,k)<=0) then
                 W(i,j,k) = Win(j,k)
              else
                 W(i,j,k) = 0
              end if
            end do
           end do
           !$omp end do
           !$omp end parallel
         end do

       else

         call par_sync_out("  ...setting initial conditions.")

         !$omp parallel private(i,j,k)
         !$omp do collapse(3)
         do k = 1, Unz
          do j = 1, Uny
           do i = 1, Unx
            if (Utype(i,j,k)<=0) then
               U(i,j,k) = Uin(j,k)
             else
               U(i,j,k) = 0
            end if
           end do
          end do
         end do
         !$omp end do nowait
         !$omp do collapse(3)
         do k = 1, Vnz
          do j = 1, Vny
           do i = 1, Vnx
            if (Vtype(i,j,k)<=0) then
               V(i,j,k) = Vin(j,k)
             else
               V(i,j,k) = 0
            end if
           end do
          end do
         end do
         !$omp end do nowait
         !$omp do collapse(3)
         do k = 1, Wnz
          do j = 1, Wny
           do i = 1, Wnx
            if (Wtype(i,j,k)<=0) then
               W(i,j,k) = Win(j,k)
             else
               W(i,j,k) = 0
            end if
           end do
          end do
         end do
         !$omp end do
         !$omp end parallel
       end if  !task_type



       if (num_of_scalars>0) then
         call par_sync_out("  ...setting initial scalar values.")
         !$omp parallel
         !$omp workshare
         SCALAR(1:Prnx,1:Prny,1:Prnz,:) = 0
         !$omp end workshare
         !$omp end parallel
       end if

       if (enable_buoyancy.and.task_type==2) then
         call par_sync_out("  ...setting initial temperature values.")

         do k = 0, Prnz+1
          do j = 0, Prny+1
           do i = 0, Prnx+1
            x = xPr(i)
            y = yPr(j)
            z = zPr(k)
            if ((x)**2+(y-0.5)**2<0.2_knd**2) then
             temperature(i,j,k) = cos(sqrt(x**2+(y-0.5)**2)*pi/2/0.2_knd)**2
            else
             temperature(i,j,k) = 0
            end if
           end do
          end do
         end do

       elseif (enable_buoyancy.and.task_type==3) then
         call par_sync_out("  ...setting initial temperature values.")

         do k = 0, Prnz+1
          do j = 0, Prny+1
           do i = 0, Prnx+1
            x = xPr(i)
            y = yPr(j)
            z = zPr(k)
            temperature(i,j,k) = temperature_ref + &
               (temperature_ref/100._knd) * ((2+cos(2*z))*(cos(2*(-y))+cos(2*(x)))-2)
           end do
          end do
         end do

       elseif (enable_buoyancy) then
         call par_sync_out("  ...setting initial temperature values.")

         call InitScalarProfile(TempIn,TemperatureProfile,temperature_ref)

         call InitScalar(TempIn,TemperatureProfile,Temperature)

       end if !buoyancy and task_type

       if (enable_moisture) then
         call par_sync_out("  ...setting initial moisture values.")

         call InitScalarProfile(MoistIn,MoistureProfile,moisture_ref)

         call InitScalar(MoistIn,MoistureProfile,Moisture)

       end if
#endif


       if (enable_buoyancy) then
         call par_sync_out("  ...setting hydrostatic pressure.")

         call InitHydrostaticPressure(Pr,Temperature,Moisture)

       end if

       call par_sync_out("  ...setting initial viscosity and diffusivity.")

       call set(Viscosity, molecular_viscosity)

       if (molecular_viscosity > 0 .and. &
             (enable_buoyancy .or. &
              enable_moisture .or. &
              num_of_scalars > 0))     then

         call set(TDiff, molecular_diffusivity)

       end if

       call par_sync_out("  ...setting ghost cell values.")

#ifdef PAR
       call BoundUVW(U, V, W)

       call par_exchange_domain_bounds(U, V, W, Temperature, Moisture, Scalar, time_stepping%time, 1._knd)

       call par_update_domain_bounds(U, V, W, Temperature, Moisture, Scalar, time_stepping%start_time)
#endif

       call BoundUVW(U, V, W)

       call Bound_Pr(Pr)

       call par_sync_out("  ...computing pressure correction.")
       call PressureCorrection(U,V,W,Pr,Q,1._knd)

       call par_sync_out("  ...computing initial eddy viscosity.")
       call SubgridModel(U, V, W)

       call par_sync_out("  ...setting viscosity in ghost cells.")
       call BoundViscosity(Viscosity)

       if (enable_buoyancy.or. &
           enable_moisture.or. &
           num_of_scalars>0)     then
         call par_sync_out("  ...setting subgrid diffusivity.")

         !$omp parallel do private(i,j,k)
         do k = 1, Prnz
           do j = 1, Prny
             do i = 1, Prnx
               TDiff(i,j,k) = (Viscosity(i,j,k) - molecular_viscosity) / constPr_sgs + &
                              molecular_diffusivity
             end do
           end do
         end do
         !$omp end parallel do

         call BoundViscosity(TDiff)
       end if

       if (enable_buoyancy) then
         call par_sync_out("  ...setting temperature in ghost cells.")
         call BoundTemperature(Temperature)
       end if

       if (enable_moisture) then
         call par_sync_out("  ...setting moisture in ghost cells.")
         call BoundMoisture(Moisture)
       end if

       if (num_of_scalars>0)  call par_sync_out("  ...setting scalars in ghost cells.")
       do i = 1, num_of_scalars
         call BoundScalar(Scalar(:,:,:,i))
       end do

       call par_sync_out("  ...setting initial temperature flux.")
       call InitTempFl(Temperature)

       if (wallmodeltype>0) then
                      call par_sync_out("  ...computing wall model in scalar points.")
                      call ComputeViscsWM(U,V,W,Pr,Temperature)
                      call par_sync_out("  ...computing wall model in velocity points.")
                      call ComputeUVWFluxesWM(U,V,W,Pr,Temperature)
       end if

       call par_sync_out("  ...setting viscosity in ghost cells.")
       call BoundViscosity(Viscosity)

       call par_sync_out("  ...initializing fixed flow_rates.")
       call InitFlowRates(U, V)

    end if !init conditions not from file

    call par_sync_out("initial conditions set.")


  end subroutine InitialConditions






  subroutine InitBoundaryConditions
    use VTKFrames, only: InitVTKFrames
    use SurfaceFrames, only: InitSurfaceFrames
    use custom_par
    
    real(knd), allocatable:: xU2(:), yV2(:), zW2(:)
    integer   :: i, j, k
    integer   :: nx, ny, nz
    integer   :: nxup, nxdown, nyup, nydown, nzup, nzdown
    real(knd) :: P, dt
    integer   :: unit, io
    
    type(spline_coefs) :: geostrophic_wind
    
#ifdef CUSTOM_BOUNDARY_CONDITIONS
    interface
      subroutine CustomBoundaryConditions
      end subroutine
    end interface
#endif

    !Important to have some defined value before the first call to GetTurbInlet.
    !The value can be quite arbitrary.
    dt = min( min(dxmin,dymin,dzmin) / Uinlet, &
              min(dxmin / time_stepping%U_max(1), &
                  dymin / time_stepping%U_max(2), &
                  dzmin / time_stepping%U_max(3)))


    call par_sync_out("  ...computing grid coordinates.")


    nx = Prnx-1
    ny = Prny-1
    nz = Prnz-1

    if (xgridfromfile) then

      call newunit(unit)

      open(unit,file="xgrid.txt")
      j = -1
      do
        read (unit,*,iostat = io) P
        if (io==0) then
          j = j+1
        else
          exit
        end if
      end do

      nx = j
      Prnx = nx
      Vnx = Prnx
      Wnx = Prnx

      if (Btype(Ea)==BC_PERIODIC) then
                            Unx = Prnx
      else
                            Unx = Prnx-1
      end if

      close(unit)

    end if

    if (ygridfromfile) then

      call newunit(unit)

      open(unit,file="ygrid.txt")
      j = -1
      do
        read (unit,*,iostat = io) P
        if (io==0) then
          j = j+1
        else
          exit
        end if
      end do

      ny = j
      Prny = ny
      Uny = Prny
      Wny = Prny

      if (Btype(No)==BC_PERIODIC) then
                            Vny = Prny
      else
                            Vny = Prny-1
      end if

      close(unit)

    end if

    if (zgridfromfile) then

      call newunit(unit)

      open(unit,file="zgrid.txt")
      j = -1
      do
        read (unit,*,iostat = io) P
        if (io==0) then
          j = j+1
          if (master) write(*,*) j
        else
          exit
        end if
      end do

      nz = j
      Prnz = nz
      Unz = Prnz
      Vnz = Prnz

      if (Btype(To)==BC_PERIODIC) then
                            Wnz = Prnz
      else
                            Wnz = Prnz-1
      end if

      close(unit)

    end if





    allocate(xU2(-3:nx+4))
    allocate(yV2(-3:ny+4))
    allocate(zW2(-3:nz+4))



    if (xgridfromfile) then

      call newunit(unit)

      open(unit,file="xgrid.txt")
      do j = 0, nx
        if (master) write(*,*) j
        read(unit,*) xU2(j)
      end do
      close(unit)

      if (Btype(We)==BC_PERIODIC) then
        do j = -1, -3, -1
          xU2(j) = xU2(0) - (xU2(nx)-xU2(nx+j))
        end do
      else
        do j = -1, -3, -1
          xU2(j) = xU2(0) - (xU2(0-j)-xU2(0))
        end do
      end if

      if (Btype(Ea)==BC_PERIODIC) then
        do j = nx+1, nx+4
          xU2(j) = xU2(nx) + (xU2(j-nx)-xU2(0))
        end do
      else
        do j = nx+1, nx+4
          xU2(j) = xU2(nx) + (xU2(nx)-xU2(nx-(j-nx)))
        end do
      end if

    else

      forall (i=-3:nx+4)
         xU2(i) = i * dxmin + im_xmin
      end forall

    end if


    if (ygridfromfile) then

      call newunit(unit)

      open(unit,file="ygrid.txt")
      do j = 0, ny
        read(unit,*) yV2(j)
      end do
      close(unit)

      if (Btype(So)==BC_PERIODIC) then
        do j = -1, -3, -1
          yV2(j) = yV2(0) - (yV2(ny)-yV2(ny+j))
        end do
      else
        do j = -1, -3, -1
          yV2(j) = yV2(0) - (yV2(0-j)-yV2(0))
        end do
      end if

      if (Btype(No)==BC_PERIODIC) then
        do j = ny+1, ny+4
          yV2(j) = yV2(ny) + (yV2(j-ny)-yV2(0))
        end do
      else
        do j = ny+1, ny+4
          yV2(j) = yV2(ny) + (yV2(ny)-yV2(ny-(j-ny)))
        end do
      end if

      y0 = yV2(0)

    else

      forall (j=-3:ny+4)
        yV2(j) = j * dymin + im_ymin
      end forall

    end if


    if (zgridfromfile) then

      call newunit(unit)

      open(unit,file="zgrid.txt")
      do j = 0, nz
        read(unit,*) zW2(j)
      end do
      close(unit)

      if (Btype(Bo)==BC_PERIODIC) then
        do j = -1, -3, -1
          zW2(j) = zW2(0) - (zW2(nz)-zW2(nz+j))
        end do
      else
        do j = -1, -3, -1
          zW2(j) = zW2(0) - (zW2(0-j)-zW2(0))
        end do
      end if

      if (Btype(To)==BC_PERIODIC) then
        do j = nz+1, nz+4
          zW2(j) = zW2(nz) + (zW2(j-nz)-zW2(0))
        end do
      else
        do j = nz+1, nz+4
          zW2(j) = zW2(nz) + (zW2(nz)-zW2(nz-(j-nz)))
        end do
      end if

      z0 = zW2(0)

    else

       forall (k=-3:nz+4)
         zW2(k) = k * dzmin + im_zmin
       end forall

    end if


    nxup = nx+1
    nxdown = 0
    nyup = ny+1
    nydown = 0
    nzup = nz+1
    nzdown = 0


    allocate(xU(-3:nx+4))
    allocate(yV(-3:ny+4))
    allocate(zW(-3:nz+4))
    allocate(dxU(-2:nx+3))
    allocate(dyV(-2:ny+3))
    allocate(dzW(-2:nz+3))
    allocate(xPr(-2:nx+4),dxPr(-2:nx+4))
    allocate(yPr(-2:ny+4),dyPr(-2:ny+4))
    allocate(zPr(-2:nz+4),dzPr(-2:nz+4))

    xU = xU2(nxdown-3:nxup+3)
    yV = yV2(nydown-3:nyup+3)
    zW = zW2(nzdown-3:nzup+3)

    forall (i=-2:nx+4)
      xPr(i) = (xU(i-1)+xU(i))/2._knd
      dxPr(i) = xU(i)-xU(i-1)
    end forall

    forall (j=-2:ny+4)
      yPr(j) = (yV(j-1)+yV(j))/2._knd
      dyPr(j) = yV(j)-yV(j-1)
    end forall

    forall (k=-2:nz+4)
      zPr(k) = (zW(k-1)+zW(k))/2._knd
      dzPr(k) = zW(k)-zW(k-1)
    end forall

    forall (i=-2:nx+3)
      dxU(i) = xPr(i+1)-xPr(i)
    end forall

    forall (j=-2:ny+3)
      dyV(j) = yPr(j+1)-yPr(j)
    end forall

    forall (k=-2:nz+3)
      dzW(k) = zPr(k+1)-zPr(k)
    end forall

    deallocate(xU2)
    deallocate(yV2)
    deallocate(zW2)

    call par_sync_out("  ...creating grid cell type arrays.")


    allocate(Utype(-2:Unx+3,-2:Uny+3,-2:Unz+3))
    allocate(Vtype(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
    allocate(Wtype(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
    allocate(Prtype(0:Prnx+1,0:Prny+1,0:Prnz+1))

    Utype = 0
    Vtype = 0
    Wtype = 0
    Prtype = 0


    call par_sync_out("  ...reading geostrophic wind.")

    !Requires grid coordinates
    call get_geostrophic_wind("geostrophic_wind_profile.conf", geostrophic_wind)


    call par_sync_out("  ...getting inlet conditions.")


    allocate(Uin(-2:Uny+3,-2:Unz+3),Vin(-2:Vny+3,-2:Vnz+3),Win(-2:Wny+3,-2:Wnz+3))
    Uin = 0
    Vin = 0
    Win = 0

    if (enable_buoyancy) allocate(TempIn(-1:Prny+2,-1:Prnz+2))

    if (enable_moisture) allocate(MoistIn(-1:Prny+2,-1:Prnz+2))

    select case (inlettype)
      case (ZeroInletType)
        Uin = 0
        Vin = 0
        Win = 0
      case (ShearInletType)
        call ShearInlet(ShearInletTypeParameter)
      case (ParabolicInletType)
        call ParabolicInlet
      case (TurbulentInletType)
        call default_turbulence_generator%init()
        call default_turbulence_generator%time_step(Uin, Vin, Win, time_stepping%dt)
      case (FromFileInletType)
        call GetInletFromFile(time_stepping%start_time)
      case (GeostrophicInletType)
        call GeostrophicWindInlet(geostrophic_wind)
      case default
        call ConstantInlet
    endselect


    if (enable_buoyancy) then
       
       call par_sync_out("  ...setting boundary temperature and temperature flux.")

       if (TempBtype(Bo)==BC_CONSTFLUX.or.TempBtype(Bo)==BC_DIRICHLET) then

         allocate(BsideTFlArr(-1:Prnx+2,-1:Prny+2))

         if (enable_radiation) then
           BsideTFlArr = 0
         else if (TempBtype(Bo)==BC_CONSTFLUX) then
           BsideTFlArr = sideTemp(Bo)
         else
           BsideTFlArr = 0
         end if

         if (TempBtype(Bo)==BC_DIRICHLET) then
           allocate(BsideTArr(-1:Prnx+2,-1:Prny+2))
           BsideTArr = sideTemp(Bo)
         end if

        end if
    end if

    if (.not.allocated(BsideTArr))  allocate(BsideTArr(0,0))
    if (.not.allocated(BsideTFlArr))  allocate(BsideTFlArr(0,0))


    if (enable_moisture) then

       call par_sync_out("  ...setting boundary moisture and moisture flux.")

       if (MoistBtype(Bo)==BC_CONSTFLUX.or.MoistBtype(Bo)==BC_DIRICHLET) then

         allocate(BsideMFlArr(-1:Prnx+2,-1:Prny+2))

         if (enable_radiation) then
           BsideMFlArr = 0
         else if (MoistBtype(Bo)==BC_CONSTFLUX) then
           BsideMFlArr = sideMoist(Bo)
         else
           BsideMFlArr = 0
         end if

         if (MoistBtype(Bo)==BC_DIRICHLET) then
           allocate(BsideMArr(-1:Prnx+2,-1:Prny+2))
           BsideMArr = sideMoist(Bo)
         end if

        end if
    end if

    if (.not.allocated(BsideMArr))  allocate(BsideMArr(0,0))
    if (.not.allocated(BsideMFlArr))  allocate(BsideMFlArr(0,0))



    call par_sync_out("  ...pressure correction.")
    call InitPressureCorrection

    call par_sync_out("  ...initializing subsidence profile.")
    call InitSubsidenceProfile

    call par_sync_out("  ...initializing loop tiles.")
    call InitTiles(Prnx,Prny,Prnz)

    if (enable_radiation) then
      call par_sync_out("  ...solar radiation.")
      call InitSolarRadiation
    end if

    call par_sync_out("  ...preparing solid bodies.")
   !prepare the geometry of the solid bodies
    call InitSolidBodies

    call par_sync_out("  ...volume source bodies.")
    !prepare the geometry of the plant bodies
    call InitVolumeSourceBodies

    call par_sync_out("  ...preparing solid bodies boundary conditions.")
    !create actual immersed boundary and wall model points, i.m. points end up in an array
    call GetSolidBodiesBC

    call par_sync_out("  ...getting outside boundaries wall model points.")
    !add also the wall model points from domain boundaries to the list
    call GetOutsideBoundariesWM(num_of_scalars)

    call par_sync_out("  ...creating wall model points arrays.")
    !create arrays of w.m. points from the list
    call MoveWMPointsToArray

    call par_sync_out("  ...creating wall model mask arrays.")
    !creates masks for computation of diffusive fluxes
    call InitWMMasks

    if (enable_radiation) then
       call par_sync_out("  ...computing initial radiation balance.")
      !compute the radiation balance and prepare the wall fluxes 
      call InitIBPFluxes
      !set the immersed boundary values for fluxes
 !      call SetIBPFluxes
    end if

    call par_sync_out("  ...setting nullified points.")
    call SetNullifiedPoints

    call par_sync_out("  ...getting volume scalar sources.")
    !create actual arrays of the source points
    call InitVolumeSources

    call par_sync_out("  ...getting area scalar sources.")
    !add arrays of the line source points
    call InitAreaSources

    call par_sync_out("  ...getting line scalar sources.")
    !add arrays of the line source points
    call InitLineSources

    call par_sync_out("  ...getting point scalar sources.")
    !add arrays of the point source points
    call InitPointSources

    call par_sync_out("  ...getting puff scalar sources.")
    !add puff sources, each containing one or more points
    call InitPuffSources
    
    !filter out frames outside the domain
    call par_sync_out("  ...preparing VTK frames.")
    call InitVTKFrames
    call par_sync_out("  ...preparing constant height frames.")
    call InitSurfaceFrames
   
#ifdef CUSTOM_BOUNDARY_CONDITIONS
    call par_sync_out("  ...executing custom boundary conditions.")
    call CustomBoundaryConditions
#endif

        call par_sync_out("boundary conditions set.")
    
  contains
  
    subroutine GeostrophicWindInlet(g)
      use Interpolation
      type(spline_coefs), intent(in) :: g
      real(knd) :: ug, vg
      integer :: j, k
      
      Win = 0
      
      j = 0
      do k = -2, Unz+3
        ug = linear_interpolation_eval(zPr(k), g%z, g%cu, j)
        Uin(:,k) = ug
      end do
      
      j = 0
      do k = -2, Vnz+3
        vg = linear_interpolation_eval(zPr(k), g%z, g%cv, j)
        Vin(:,k) = vg
      end do
    end subroutine
    
  end subroutine InitBoundaryConditions


  subroutine InitFlowRates(U, V)
#ifdef PAR
    use custom_par
#endif
    real(knd), intent(in), contiguous, dimension(-2:,-2:,-2:) :: U, V

    if (enable_fixed_flow_rate) then

      if (iim==nxims .and. Btype(Ea)==BC_PERIODIC) then
        flow_rate_x_fixed = .true.
        flow_rate_x = sum(U(Unx, 1:Uny, 1:Unz)) * dymin * dzmin
#ifdef PAR
        flow_rate_x = par_co_sum_plane_yz(flow_rate_x)
#endif
      end if

#ifdef PAR        
      flow_rate_x_fixed = par_co_any(flow_rate_x_fixed)     
      if (flow_rate_x_fixed) call par_broadcast_from_last_x(flow_rate_x)
#endif


      if (jim==nyims .and. Btype(No)==BC_PERIODIC) then
        flow_rate_y_fixed = .true.
        flow_rate_y = sum(V(1:Vnx, Vny, 1:Vnz)) * dxmin * dzmin
#ifdef PAR
        flow_rate_y = par_co_sum_plane_xz(flow_rate_y)
#endif
      end if

#ifdef PAR        
      flow_rate_y_fixed = par_co_any(flow_rate_y_fixed)     
      if (flow_rate_y_fixed) call par_broadcast_from_last_y(flow_rate_y)
#endif

    end if

  end subroutine




  subroutine SetNullifiedPoints
    integer :: i,j,k,n

    !$omp parallel do reduction(+:nUnull)
    do k = 1, Unz
     do j = 1, Uny
      do i = 1, Unx
       if (Utype(i,j,k)>0.and.Utype(i,j,k+1)>0.and.Utype(i,j,k-1)>0&
           .and.Utype(i,j-1,k)>0.and.Utype(i,j+1,k)>0&
           .and.Utype(i-1,j,k)>0.and.Utype(i+1,j,k)>0)  nUnull = nUnull+1
      end do
     end do
    end do
    !$omp end parallel do

    allocate(Unull(3,nUnull))

    n = 0

    do k = 1, Unz
     do j = 1, Uny
      do i = 1, Unx
       if (Utype(i,j,k)>0.and.Utype(i,j,k+1)>0.and.Utype(i,j,k-1)>0&
           .and.Utype(i,j-1,k)>0.and.Utype(i,j+1,k)>0&
           .and.Utype(i-1,j,k)>0.and.Utype(i+1,j,k)>0)  then
            !$omp atomic
            n = n+1
            Unull(:,n) = [ i,j,k ]

       end if
      end do
     end do
    end do

    nVnull = 0

    !$omp parallel do reduction(+:nVnull)
    do k = 1, Vnz
     do j = 1, Vny
      do i = 1, Vnx
       if (Vtype(i,j,k)>0.and.Vtype(i,j,k+1)>0.and.Vtype(i,j,k-1)>0&
           .and.Vtype(i,j-1,k)>0.and.Vtype(i,j+1,k)>0&
           .and.Vtype(i-1,j,k)>0.and.Vtype(i+1,j,k)>0)  nVnull = nVnull+1
      end do
     end do
    end do
    !$omp end parallel do

    allocate(Vnull(3,nVnull))

    n = 0

    do k = 1, Vnz
     do j = 1, Vny
      do i = 1, Vnx
       if (Vtype(i,j,k)>0.and.Vtype(i,j,k+1)>0.and.Vtype(i,j,k-1)>0&
           .and.Vtype(i,j-1,k)>0.and.Vtype(i,j+1,k)>0&
           .and.Vtype(i-1,j,k)>0.and.Vtype(i+1,j,k)>0)  then

            n = n+1
            Vnull(:,n) = [ i,j,k ]

       end if
      end do
     end do
    end do

    nWnull = 0

    !$omp parallel do reduction(+:nWnull)
    do k = 1, Wnz
     do j = 1, Wny
      do i = 1, Wnx
       if (Wtype(i,j,k)>0.and.Wtype(i,j,k+1)>0.and.Wtype(i,j,k-1)>0&
           .and.Wtype(i,j-1,k)>0.and.Wtype(i,j+1,k)>0&
           .and.Wtype(i-1,j,k)>0.and.Wtype(i+1,j,k)>0)  nWnull = nWnull+1
      end do
     end do
    end do
    !$omp end parallel do

    allocate(Wnull(3,nWnull))

    n = 0


    do k = 1, Wnz
     do j = 1, Wny
      do i = 1, Wnx
       if (Wtype(i,j,k)>0.and.Wtype(i,j,k+1)>0.and.Wtype(i,j,k-1)>0&
           .and.Wtype(i,j-1,k)>0.and.Wtype(i,j+1,k)>0&
           .and.Wtype(i-1,j,k)>0.and.Wtype(i+1,j,k)>0)  then

            n = n+1
            Wnull(:,n) = [ i,j,k ]

       end if
      end do
     end do
    end do

  end subroutine SetNullifiedPoints



  subroutine get_scratch
    use strings, only: itoa
    integer :: l, stat

    call get_environment_variable(name="SCRATCHDIR", length=l, status=stat)

    if (stat==0 .and. l>0 .and. l<=len(scratch_dir)) then
      call get_environment_variable(name="SCRATCHDIR", value=scratch_dir)
    else if (stat==0 .and. l>len(scratch_dir)) then
      call error_stop("SCRATCHDIR length exceeds variable length "//itoa(len(scratch_dir)))
    end if
  end subroutine


  subroutine init_random_seed()
    use rng_par_zig
    !$ use omp_lib
#ifdef PAR
    use custom_par
#endif
    integer :: i
    integer(int32),dimension(:),allocatable :: seed
    integer :: nt
    
    nt = 1
    !$omp parallel
    !$ nt = omp_get_num_threads()
    !$omp end parallel
    
    allocate(seed(nt))

#ifdef PAR
    seed = [( 1000 * i + myrank, i=1,size(seed) )]
#else
    seed = [( 1000 * i, i=1,size(seed) )]
#endif

    call rng_init(nt, seed)

  end subroutine


endmodule Initial

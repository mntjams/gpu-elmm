module Initial

  use PARAMETERS
  use ArrayUtilities, only: avg
  use LIMITERS, only: limparam, limitertype
  use MULTIGRID, only: SetMGParams
  use MULTIGRID2d, only: SetMGParams2d
  use POISSON
  use BOUNDARIES
  use ScalarBoundaries
  use OUTPUTS, only: store, display, probes, scalar_probes, frame_flags, &
                     AddFrameDomain, SetFrameDomains, StaggeredFrameDomains, ReadProbes
  use SCALARS
  use Filters, only: filtertype, filter_ratios
  use Subgrid
  use TURBINLET, only: GetTurbulentInlet, GetInletFromFile, TLag, Lturby, Lturbz, Ustar_inlet, relative_stress, &
               Ustar_surf_inlet, stress_gradient_inlet, U_ref_inlet, z_ref_inlet, z0_inlet, power_exponent_inlet
  use SolarRadiation, only: InitSolarRadiation
  use SolidBodies, only: obstacles_file, InitSolidBodies
  use ImmersedBoundary, only: GetSolidBodiesBC, InitIBPFluxes!, SetIBPFluxes
  use VolumeSources!, only: InitVolumeSources, InitVolumeSourceBodies, ScalarFlVolume, ScalarFlVolumesContainer
  use LineSources, only: InitLineSources
  use WALLMODELS
  use TILING, only: tilesize,InitTiles
  use FreeUnit, only: newunit
  use Puffs, only: InitPuffSources

  implicit none

  private
  public  ReadConfiguration, InitialConditions, InitBoundaryConditions

  real(knd) x0,y0,z0 !domain boundaries, will become xU(0), yV(0), zW(0)
  real(knd) lx,ly,lz !domain extents

contains


 subroutine ReadConfiguration
   use StaggeredFrames, only: rrange, TFrameTimes, TSaveFlags, Init
   integer   lmg,minmglevel,bnx,bny,bnz,mgncgc,mgnpre,mgnpost,mgmaxinnerGSiter,minGPUlevel
   real(knd) mgepsinnerGS
   integer   i,io,io2,itmp
   integer numframeslices

   character(80), save :: probes_file = ""
   character(80), save :: scalar_probes_file = ""

   character(len = 1024) :: commandline,msg
   integer :: exenamelength
   integer :: unit

   type(rrange) :: range
   type(TFrameTimes) :: frame_times
   type(TSaveFlags) :: frame_save_flags
   character(10) :: domain_label
   integer :: num_staggered_domains
   integer :: number_of_probes, number_of_scalar_probes

   integer :: dimension,direction
   real(knd) :: position

   interface get
     procedure chget1
     procedure lget1, lget2, lget3
     procedure iget1, iget2, iget3
     procedure rget1, rget2, rget3
     procedure rgetv3
   end interface

   call read_command_line

   call parse_command_line

   call newunit(unit)

   open(unit,file="main.conf",status="old",action="read")
   call get(CFL)
   call get(Uref)
   call get(poissmet)
   call get(convmet)
   call get(limitertype)
   call get(limparam)
   call get(masssourc)
   call get(steady)
   call get(tasktype)
   write(*,*) "tasktype=",tasktype
   call get(initcondsfromfile)
   call get(timeavg1)
   call get(timeavg2)
   call get(Re)
   write(*,*) "Re=",Re
   call get(start_time)
   write(*,*) "start_time=",start_time
   call get(end_time)
   write(*,*) "end_time=",end_time
   call get(max_number_of_time_steps)
   write(*,*) "max_number_of_time_steps=",max_number_of_time_steps
   call get(eps)
   write(*,*) "eps=",eps
   call get(maxCNiter)
   write(*,*) "maxCNiter=",maxCNiter
   call get(epsCN)
   write(*,*) "epsCN=",epsCN
   call get(maxPOISSONiter)
   write(*,*) "maxPOISSONiter=",maxPOISSONiter
   call get(epsPOISSON)
   write(*,*) "epsPOISSON=",epsPOISSON
   call get(debugparam)
   write(*,*) "debug parameter=",debugparam
   close(unit)


   open(unit,file="les.conf",status="old",action="read")
   call get(sgstype)
   call get(filtertype)

   if (filtertype > size(filter_ratios)) then
     write(*,*) "Chosen filter type does not exist. Maximum index is:",size(filter_ratios)
     stop
   end if

   call get(wallmodeltype)
   close(unit)


   open(unit,file="grid.conf",status="old",action="read")
   call get(xgridfromfile)
   call get(ygridfromfile)
   call get(zgridfromfile)
   read(unit,fmt='(/)')

   read(unit,*) x0
   write(*,*) "x0=",x0
   call get(y0)
   write(*,*) "y0=",y0
   call get(z0)
   write(*,*) "z0=",z0
   read(unit,fmt='(/)')

   read(unit,*) lx
   if (lx>0) then
     write(*,*) "lx=",lx
   else
     write (*,*) "Domain length in x direction must be positive."
     stop
   end if
   read(unit,fmt='(/)')

   read(unit,*) ly
   if (ly>0) then
     write(*,*) "ly=",ly
   else
     write (*,*) "Domain length in y direction must be positive."
     stop
   end if
   read(unit,fmt='(/)')

   read(unit,*) lz
   if (lz>0) then
     write(*,*) "lz=",lz
   else
     write (*,*) "Domain length in z direction must be positive."
     stop
   end if
   read(unit,fmt='(/)')

   read(unit,*) Prnx
   if (Prnx>0) then
     write(*,*) "nx=",Prnx
   else
     write (*,*) "Number of cells in x direction must be positive."
     stop
   end if
   read(unit,fmt='(/)')

   read(unit,*) Prny
   if (Prny>0) then
     write(*,*) "ny=",Prny
   else
     write (*,*) "Number of cells in y direction must be positive."
     stop
   end if
   read(unit,fmt='(/)')

   read(unit,*) Prnz
   if (Prnz>0) then
     write(*,*) "nz=",Prnz
   else
     write (*,*) "Number of cells in z direction must be positive."
     stop
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
     call get(CoriolisParam)
     write(*,*) "coriolisparam=",CoriolisParam
     call get(PrGradientX)
     write(*,*) "prgradientx=",PrGradientX
     call get(PrGradientY)
     write(*,*) "prgradienty=",PrGradientY
     call get(SubsidenceGradient)
     write(*,*) "SubsidenceGradient=",SubsidenceGradient
     close(unit)

   else

     write(*,*) "Warning! Could not open file large_scale.conf. Using defaults."

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
     enable_buoyancy = 0
   end if

   if (enable_buoyancy==1) then

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

       write(*,*) "Warning! Could not open file temp_profile.conf. Using defaults."
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
     enable_moisture = 0
   end if

   if (enable_moisture==1) then

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

       write(*,*) "Warning! Could not open file moist_profile.conf. Using defaults."
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
   write(*,*) "G=",ShearInletTypeParameter
   call get(Uinlet)
   write(*,*) "Uinlet=",Uinlet
   call get(Ustar_surf_inlet)  !-<u'w'>
   call get(stress_gradient_inlet) !in relative part per 1m
   call get(z0_inlet)
   call get(power_exponent_inlet)
   call get(z_ref_inlet)
   call get(U_ref_inlet)
   call get(relative_stress(1,1))
   call get(relative_stress(2,2))
   call get(relative_stress(3,3))
   call get(relative_stress(1,2))
   call get(relative_stress(1,3))
   call get(relative_stress(2,3))
   call get(TLag)
   call get(Lturby)
   call get(Lturbz)
   close(unit)

   relative_stress(2,1) = relative_stress(1,2)
   relative_stress(3,1) = relative_stress(1,3)
   relative_stress(3,2) = relative_stress(2,3)

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
          allocate(scalsrcx(partdistrib),scalsrcy(partdistrib),scalsrcz(partdistrib))
          allocate(scalsrci(partdistrib),scalsrcj(partdistrib),scalsrck(partdistrib))

          do i = 1,partdistrib
            call get(partdiam(i))
            call get(partrho(i))
            call get(percdistrib(i))
            call get(scalsrcx(i))
            call get(scalsrcy(i))
            call get(scalsrcz(i))
          end do

       else

          allocate(partdiam(num_of_scalars),partrho(num_of_scalars),percdistrib(num_of_scalars))
          allocate(scalsrcx(num_of_scalars),scalsrcy(num_of_scalars),scalsrcz(num_of_scalars))
          allocate(scalsrci(num_of_scalars),scalsrcj(num_of_scalars),scalsrck(num_of_scalars))

          do i = 1,num_of_scalars
            call get(partdiam(i))
            call get(partrho(i))
            call get(percdistrib(i))
            call get(scalsrcx(i))
            call get(scalsrcy(i))
            call get(scalsrcz(i))
          end do
       end if
     end if
     close(unit)
   else
     num_of_scalars = 0
     write (*,*) "scalars.conf not found, no passive scalars for computation."
   end if


   if (num_of_scalars>0) then
      open(unit, file="line_sources.conf",status="old",action="read",iostat=io)
      
      if (io==0) then
        call get_line_sources
        close(unit)
      end if
   end if

   if (poissmet==3.or.poissmet==4.or.poissmet==5) then
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
     read(unit,fmt='(/)',iostat = io)
     read(unit,*,iostat = io) minGPUlevel
     close(unit)

     if (poissmet==3.or.poissmet==4) then
       if (Prny==1) then
        call SetMGParams2d(llmg = lmg,lminmglevel = minmglevel,lbnx = bnx,lbnz = bnz,&
                           lmgncgc = mgncgc,lmgnpre = mgnpre,lmgnpost = mgnpost,&
                           lmgmaxinnerGSiter = mgmaxinnerGSiter,lmgepsinnerGS = mgepsinnerGS)
       else
        call SetMGParams(llmg = lmg,lminmglevel = minmglevel,lminGPUlevel = minGPUlevel,&
                           lbnx = bnx,lbny = bny,lbnz = bnz,&
                           lmgncgc = mgncgc,lmgnpre = mgnpre,lmgnpost = mgnpost,&
                           lmgmaxinnerGSiter = mgmaxinnerGSiter,lmgepsinnerGS = mgepsinnerGS)
       end if
     end if
   end if

   open(unit,file="frames.conf",status="old",action="read",iostat = io)
   if (io==0) then
     call get(frames)
     call get(timefram1)
     call get(timefram2)
     read(unit,fmt='(/)')

     read(unit,*) frame_flags%U
     call get(frame_flags%vort)
     call get(frame_flags%Pr)
     call get(frame_flags%lambda2)
     call get(frame_flags%scalars)
     if (num_of_scalars < 1) frame_flags%scalars = 0
     call get(frame_flags%sumscalars)
     if (num_of_scalars < 1) frame_flags%sumscalars = 0
     call get(frame_flags%temperature)
     if (enable_buoyancy /= 1) frame_flags%temperature = 0
     call get(frame_flags%moisture)
     if (enable_buoyancy /= 1) frame_flags%moisture = 0
     call get(frame_flags%temperature_flux)
     if (enable_buoyancy /= 1) frame_flags%temperature_flux = 0
     call get(frame_flags%scalfl)
     if (num_of_scalars < 1) frame_flags%scalfl = 0

     call get(numframeslices)

     do i = 1,numframeslices
       call get(dimension)
       call get(direction)
       call get(position)
       call AddFrameDomain(dimension,direction,position)
     end do

   else
     frames = 0
     write (*,*) "frames.conf not found, no vtk frames will be saved."
   end if
   close(unit)


   open(unit,file="stagframes.conf",status="old",action="read",iostat = io)
   if (io==0) then
     call get(num_staggered_domains)

     allocate(StaggeredFrameDomains(num_staggered_domains))

     do i = 1,num_staggered_domains
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
       if (enable_buoyancy /= 1) frame_save_flags%Temperature = .false.
       call get(frame_save_flags%Moisture)
       if (enable_moisture /= 1) frame_save_flags%Moisture = .false.
       call get(frame_save_flags%Scalar)
       if (num_of_scalars < 1) frame_save_flags%Scalar = .false.

       call Init(StaggeredFrameDomains(i), trim(domain_label), &
                 range, &
                 frame_times, &
                 frame_save_flags )
     end do
     close(unit)
   else
     allocate(StaggeredFrameDomains(0))
     write (*,*) "stagframes.conf not found, no staggered frames will be saved."
   end if




   open(unit,file="output.conf",status="old",action="read",iostat = io)
   if (io==0) then
     call read_namelist_output

     if (enable_buoyancy /= 1) then
       store%out_temperature = 0
       store%out_moisture = 0
     end if

     if (enable_buoyancy/=1) then
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
     write(*,*) "No output.conf found, defaults will be used."
   end if

   !probes_file and scalar_probes_file read from command line
   if (probes_file == "" .and. scalar_probes_file == "") then

     open(unit,file="probes.conf",status="old",action="read",iostat = io)
     if (io==0) then
       call get(number_of_probes)

       allocate(probes(number_of_probes))

       do i = 1,number_of_probes
         call get(probes(i)%x)
         call get(probes(i)%y)
         call get(probes(i)%z)
       end do

       scalar_probes = probes
     else
       allocate(probes(0))
       allocate(scalar_probes(0))
     end if

   else

     close(unit)

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
     read(unit,fmt='(/)')
     read(unit,'(a)') obstacles_file
     close(unit)
   end if

   write(*,*) "num_of_scalars",num_of_scalars

   windangle = 0._knd

   projectiontype = 1

   call parse_command_line

   if (CFL<=0)  CFL = 0.5


   write(*,*) "Boundaries:"

   write(*,'(a2)',advance='no') " W "
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

   write(*,'(a2)',advance='no') " E "
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

   write(*,'(a2)',advance='no') " S "
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

   write(*,'(a2)',advance='no') " N "
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

   write(*,'(a2)',advance='no') " B "
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

   write(*,'(a2)',advance='no') " T "
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

   write(*,*) "dxmin ",dxmin
   write(*,*) "dymin ",dymin
   write(*,*) "dzmin ",dzmin

   write(*,*) "lx:",lx
   write(*,*) "ly:",ly
   write(*,*) "lz:",lz


   if (Btype(Ea)==PERIODIC) then
                          Unx = Prnx
   else
                          Unx = Prnx-1
   end if
   Uny = Prny
   Unz = Prnz

   Vnx = Prnx
   if (Btype(No)==PERIODIC) then
                          Vny = Prny
   else
                          Vny = Prny-1
   end if
   Vnz = Prnz

   Wnx = Prnx
   Wny = Prny
   if (Btype(To)==PERIODIC) then
                          Wnz = Prnz
   else
                          Wnz = Prnz-1
   end if

   if (Btype(We)==TURBULENTINLET) inlettype = TurbulentInletType
   if (Btype(We)==INLETFROMFILE) inlettype = FromFileInletType


   if (Abs(Uinlet)>0) then
     dt = Abs(dxmin/Uinlet)
   else
     dt = dxmin
   end if

   if ((timeavg1>=0).and.(timeavg2>=timeavg1)) then
     averaging = 1
   else
     averaging = 0
   end if

   if (.not.xgridfromfile.and..not.ygridfromfile.and..not.zgridfromfile) then
     gridtype = UNIFORMGRID
     write(*,*) "Uniform grid"
   else
     gridtype = GENERALGRID
     write(*,*) "General grid"
   end if


   write(*,*) "set"

   contains

     subroutine chget1(x)
       character(*),intent(out) :: x
       read(unit,fmt='(/)')
       read(unit,*) x
     end subroutine
     subroutine lget1(x)
       logical,intent(out) :: x
       read(unit,fmt='(/)')
       read(unit,*) x
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

     subroutine read_command_line
       commandline = ""
       call get_command(command = commandline,status = io)
       if (io==0) then
         call get_command_argument(0,length = exenamelength,status = io2)
         if (io2==0) then
           commandline="&cmd "//adjustl(trim(commandline(exenamelength+1:)))//" /"
         else
           commandline="&cmd "//adjustl(trim(commandline))//" /"
         end if
       else
         write(*,*) io,"Error getting command line."
       end if
     end subroutine

     subroutine parse_command_line
       namelist /cmd/ tilesize, debugparam, debuglevel, windangle, projectiontype, &
                       Prnx, Prny, Prnz,&
                       obstacles_file, probes_file, scalar_probes_file

       if (len_trim(commandline)>0) then
         msg = ''
         read(commandline,nml = cmd,iostat = io,iomsg = msg)
         if (io/=0) then
           write(*,*) io,"Error parsing command line."
           write(*,*) msg
           write(*,*) commandline
         end if
       else
         write(*,*) io,"Error getting command line."
       end if
     end subroutine

     subroutine get_line_sources
        use LineSources, only: ScalarLineSource, ScalarLineSources
        type(ScalarLineSource) :: src
        integer n

        call get(n)
        allocate(ScalarLineSources(0))
        do i=1,n
          read(unit,fmt=*)
          call get(src%scalar_number)
          call get(src%start)
          call get(src%end)
          call get(src%flux)
          src%number_of_points = 20 * max(Prnx,Prny,Prnz)
          ScalarLineSources = [ScalarLineSources, src]
        end do
     end subroutine

     subroutine read_namelist_output
       namelist /output/ store, display

       read(unit,nml = output,iostat = io,iomsg = msg)
     end subroutine

 end subroutine ReadConfiguration


 subroutine ReadIC(U,V,W,Pr,Temperature)
   real(knd),intent(inout) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:)
   real(knd),intent(inout) :: Pr(1:,1:,1:)
   real(knd),intent(inout) :: Temperature(-1:,-1:,-1:)
   integer i,j,k,unit

   call newunit(unit)

   open(unit,file="in.vtk",position="rewind",status="old",action="read")
   do i = 1,14
    read(unit,*)
   end do
   do k = 1,Prnz
    do j = 1,Prny
     do i = 1,Prnx
      read(unit,*) Pr(i,j,k)
     end do
    end do
   end do
   if (enable_buoyancy==1) then
    do i = 1,3
     read(unit,*)
    end do
    do k = 1,Prnz
     do j = 1,Prny
      do i = 1,Prnx
       read(unit,*) temperature(i,j,k)
      end do
     end do
    end do
   end if
   close(unit)

   open(unit,file="Uin.vtk",position="rewind",status="old",action="read")
   do i = 1,14
    read(unit,*)
   end do
   do k = 1,Unz
    do j = 1,Uny
     do i = 1,Unx
      read(unit,*) U(i,j,k)
     end do
    end do
   end do
   close(unit)

   open(unit,file="Vin.vtk",position="rewind",status="old",action="read")
   do i = 1,14
    read(unit,*)
   end do
   do k = 1,Vnz
    do j = 1,Vny
     do i = 1,Vnx
      read(unit,*) V(i,j,k)
     end do
    end do
   end do
   close(unit)

   open(unit,file="Win.vtk",position="rewind",status="old",action="read")
   do i = 1,14
    read(unit,*)
   end do
   do k = 1,Wnz
    do j = 1,Wny
     do i = 1,Wnx
      read(unit,*) W(i,j,k)
     end do
    end do
   end do
 endsubroutine ReadIC



  subroutine InitialConditions(U,V,W,Pr,Temperature,Moisture,Scalar)
  real(knd),intent(inout) :: U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  real(knd),intent(inout) :: Temperature(-1:,-1:,-1:)
  real(knd),intent(inout) :: Moisture(-1:,-1:,-1:)
  real(knd),intent(inout) :: Scalar(-1:,-1:,-1:,:)
  integer i,j,k
  real(knd) p,x,y,z,x1,x2,y1,y2,z1,z2
  real(knd),allocatable :: Q(:,:,:)

    call init_random_seed

    Pr(1:Prnx,1:Prny,1:Prnz) = 0

    U = huge(1._knd)/2
    V = huge(1._knd)/2
    W = huge(1._knd)/2

    V(1:Vnx,1:Vny,1:Vnz) = 0
    W(1:Wnx,1:Wny,1:Wnz) = 0

    if (initcondsfromfile==1) then

       call ReadIC(U,V,W,Pr,Temperature)

       if (Re>0) then
         Viscosity = 1._knd/Re
       else
         Viscosity = 0
       end if
       if (enable_buoyancy==1.or. &
           enable_moisture==1.or. &
           num_of_scalars>0)        TDiff = 1._knd/(Re*Prandtl)

       call BoundU(1,U,Uin)
       call BoundU(2,V,Vin)
       call BoundU(3,W,Win)
       call Bound_Pr(Pr)

    else   !init conditions not from file

       if (tasktype==2) then
         U(1:Unx,1:Uny,1:Unz) = 0
         do k = 1,Unz
          do j = 1,Uny
           do i = 1,Unx
            !if (Utype(i,j,k)<=0) then
                  !call RANDOM_NUMBER(p)
                  x = xU(i)
                  y = yPr(j)
                  z = zPr(k)
                  U(i,j,k)=-2*pi*y!*(1+0.1*(p-0.5))!-Uinlet*cos(z)*sin(x)*cos(y)!
              !0.5_knd*(p-0.5_knd)z
             !else
             !  V(i,j,k) = 0
            !end if
           end do
          end do
         end do
         do k = 1,Vnz
          do j = 1,Vny
           do i = 1,Vnx
                  !call RANDOM_NUMBER(p)
                  x = xPr(i)
                  y = yV(j)
                  z = zPr(k)
                  V(i,j,k) = 2*pi*x!*(1+0.1*(p-0.5))!0
           end do
          end do
         end do
         do k = 1,Wnz
          do j = 1,Wny
           do i = 1,Wnx
                  x = xPr(i)
                  y = yPr(j)
                  z = zW(k)
                  W(i,j,k) = 0
           end do
          end do
         end do
         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
                  x = xPr(i)
                  y = yPr(j)
                  z = zPr(k)
                  Pr(i,j,k) = 0!(Uinlet/16._knd)*((2+cos(2*z))*(cos(2*(x))+cos(2*(y)))-2)
           end do
          end do
         end do

       elseif (tasktype==3) then
         U(1:Unx,1:Uny,1:Unz) = 0
         do k = 1,Unz
          do j = 1,Uny
           do i = 1,Unx
            !if (Utype(i,j,k)<=0) then
                  !call RANDOM_NUMBER(p)
                  x = xU(i)
                  y = yPr(j)
                  z = zPr(k)
                  U(i,j,k) = Uinlet*sin(x)*cos(z)*cos(-y)!*(1+0.1*(p-0.5))!-Uinlet*cos(z)*sin(x)*cos(y)!
              !0.5_knd*(p-0.5_knd)z
             !else
             !  V(i,j,k) = 0
            !end if
           end do
          end do
         end do
         do k = 1,Vnz
          do j = 1,Vny
           do i = 1,Vnx
                  !call RANDOM_NUMBER(p)
                  x = xPr(i)
                  y = yV(j)
                  z = zPr(k)
                  V(i,j,k) = 0!*(1+0.1*(p-0.5))!0
           end do
          end do
         end do
         do k = 1,Wnz
          do j = 1,Wny
           do i = 1,Wnx
                  x = xPr(i)
                  y = yPr(j)
                  z = zW(k)
                  W(i,j,k)=-Uinlet*cos(x)*sin(z)*cos(-y)!Uinlet*sin(z)*cos(x)*cos(y)!
           end do
          end do
         end do
         do k = 1,Prnz
          do j = 1,Prny
           do i = 1,Prnx
                  x = xPr(i)
                  y = yPr(j)
                  z = zPr(k)
                  Pr(i,j,k)=(Uinlet/16._knd)*((2+cos(2*z))*(cos(2*(-y))+cos(2*(x)))-2)!(Uinlet/16._knd)*((2+cos(2*z))*(cos(2*(x))+cos(2*(y)))-2)
           end do
          end do
         end do

       elseif (tasktype==5) then!temporal mixing layer
         do k = 1,Unz
          do j = 1,Uny
           do i = 1,Unx
                  y = yPr(j)
                  if ((abs(y-yPr(Uny/2)))<1.) then
                   call RANDOM_NUMBER(p)
                  else
                   p = 0
                  end if
                   x1 = xPr(i)
                   x2 = xPr(i+1)
                   y1 = yV(j-1)
                   y2 = yV(j)
                   z1 = zW(k-1)
                   z2 = zW(k)
                   p = 0
                     !(20/(pi*pi))*(cos((pi/10)*x1)-cos((pi/10)*x2))*((pi/2)*(z2-z1)-0.8*(cos((pi/2)*z2)-cos((pi/2)*z1)))
                    !sin(pi*xU(i)/10)*(1+0.8*sin(pi*zPr(k)/2))
                   x2 = cosh(y2+0.4_knd*p-yPr(Uny/2))
                   x1 = cosh(y1+0.4_knd*p-yPr(Uny/2))
                   if (abs(x2-x1)>0.000001_knd) then
                    U(i,j,k)=(log(x2)-log(x1))/(y2-y1)
                   else
                    U(i,j,k) = 0
                   end if
                      !Uinlet*tanh(y+0.1*sin(pi*(xU(i)/10)*sin(pi*zPr(k)/10)-yPr(Uny/2)) by mean of integrals
           end do
          end do
         end do
         do k = 1,Vnz
          do j = 1,Vny
           do i = 1,Vnx
                  y = yV(j)
                  !if ((abs(y-yPr(Uny/2)))<1.) then
                  ! call RANDOM_NUMBER(p)
                  !else
                   p = sin(2*pi*xPr(i)/10)*(1+0.3*sin(3*pi*zPr(k)/10))*exp(-(yV(j)-yPr(Uny/2))*(yV(j)-yPr(Uny/2)))
                  !end if

                   V(i,j,k) = 0.1*p!Uinlet*tanh(y-yPr(Uny/2))+Uinlet*0.1_knd*(p-0.5_knd)
           end do
          end do
         end do
         do k = 1,Wnz
          do j = 1,Wny
           do i = 1,Wnx
                  y = yPr(j)
                  !if ((abs(y-yPr(Uny/2)))<1.) then
                  ! call RANDOM_NUMBER(p)
                  !else
                   p = cos(2*pi*xPr(i)/10)*(1+0.3*cos(3*pi*zW(k)/10))*exp(-(yPr(j)-yPr(Uny/2))*(yPr(j)-yPr(Uny/2)))
                  !end if

                   W(i,j,k) = 0.1*p!Uinlet*tanh(y-yPr(Uny/2))+Uinlet*0.1_knd*(p-0.5_knd)
           end do
          end do
         end do
      !   elseif (tasktype==8) then
      !    do k = 1,Unz
      !     do j = 1,Uny
      !      do i = 1,Unx
      !             call RANDOM_NUMBER(p)
      !       U(i,j,k)=-prgradienty/(coriolisparam)*(1+0.1_knd*(p-0.5_knd))
      !      end do
      !     end do
      !    end do
      !    do k = 1,Vnz
      !     do j = 1,Vny
      !      do i = 1,Vnx
      !               call RANDOM_NUMBER(p)
      !      V(i,j,k) = prgradientx/(coriolisparam)*(1+0.1_knd*(p-0.5_knd))
      !      end do
      !     end do
      !    end do
      !    do k = 1,Wnz
      !     do j = 1,Wny
      !      do i = 1,Wnx
      !       W(i,j,k) = 0
      !      end do
      !     end do
      !    end do

       elseif (InletType==TurbulentInletType) then

         dt = hypot(dxmin,dymin) / hypot(avg(Uin(1:Uny,1:Unz)),avg(Vin(1:Vny,1:Vnz)))

         do i = 1,Prnx
           
           call GetTurbulentInlet
           !$omp parallel private(j,k)
           !$omp do
           do k = 1,Unz
            do j = 1,Uny
              if (Utype(i,j,k)<=0) then
                    U(i,j,k) = Uin(j,k)
              else
                 U(i,j,k) = 0
              end if
            end do
           end do
           !$omp end do nowait
           !$omp do
           do k = 1,Vnz
            do j = 1,Vny
              if (Vtype(i,j,k)<=0) then
                    V(i,j,k) = Vin(j,k)
              else
                 V(i,j,k) = 0
              end if
            end do
           end do
           !$omp end do nowait
           !$omp do
           do k = 1,Wnz
            do j = 1,Wny
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

!          !$omp parallel private(i,j,k)
!          !$omp do
         do k = 1,Unz
          do j = 1,Uny
           do i = 1,Unx
            if (Utype(i,j,k)<=0) then
                  call RANDOM_NUMBER(p)
                  U(i,j,k) = Uin(j,k)!+(Sqrt((Uin(j,k))**2+(Vin(j,k))**2))*0.1_knd*(p-0.5_knd)!sin(2.*pi*xU(i)+1)*cos(2.*pi*yPr(j)-2)!Uin(j,k)!*(1+0.03_knd*(p-0.5_knd))
             else
               U(i,j,k) = 0
            end if
           end do
          end do
         end do
!          !$omp end do
!          !$omp do
         do k = 1,Vnz
          do j = 1,Vny
           do i = 1,Vnx
            if (Vtype(i,j,k)<=0) then
                  call RANDOM_NUMBER(p)
                  V(i,j,k) = Vin(j,k)!+(Sqrt((Uin(j,k))**2+(Vin(j,k))**2))*0.1_knd*(p-0.5_knd)!-cos(2.*pi*xPr(i)+1)*sin(2.*pi*yV(j)-2)!Uinlet*(0.3_knd*(p-0.5_knd))
             else
               V(i,j,k) = 0
            end if
           end do
          end do
         end do
!          !$omp end do
!          !$omp do
         do k = 1,Wnz
          do j = 1,Wny
           do i = 1,Wnx
            if (Wtype(i,j,k)<=0) then
                  call RANDOM_NUMBER(p)
                  W(i,j,k) = Win(j,k)!Uinlet*(0.00001_knd*(p-0.5_knd))
             else
               W(i,j,k) = 0
            end if
           end do
          end do
         end do
!          !$omp end do
!          !$omp end parallel
       end if  !tasktype



       if (num_of_scalars>0) then
         !$omp parallel
         !$omp workshare
         SCALAR(1:Prnx,1:Prny,1:Prnz,:) = 0
         !$omp end workshare
         !$omp end parallel
       end if

       if (enable_buoyancy==1.and.tasktype==2) then

         do k = 0,Prnz+1
          do j = 0,Prny+1
           do i = 0,Prnx+1
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

       elseif (enable_buoyancy==1.and.tasktype==3) then

         do k = 0,Prnz+1
          do j = 0,Prny+1
           do i = 0,Prnx+1
            x = xPr(i)
            y = yPr(j)
            z = zPr(k)
            temperature(i,j,k) = temperature_ref + &
               (temperature_ref/100._knd) * ((2+cos(2*z))*(cos(2*(-y))+cos(2*(x)))-2)
           end do
          end do
         end do

       elseif (enable_buoyancy==1) then

         call InitScalarProfile(TempIn,TemperatureProfile,temperature_ref)

         call InitScalar(TempIn,TemperatureProfile,Temperature)

       end if !byoyancy and tasktype

       if (enable_moisture==1) then

         call InitScalarProfile(MoistIn,MoistureProfile,moisture_ref)

         call InitScalar(MoistIn,MoistureProfile,Moisture)

       end if

       if (enable_buoyancy==1) then

         call InitHydrostaticPressure(Pr,Temperature,Moisture)

       end if


       if (Re>0) then
         !$omp parallel
         !$omp workshare
         Viscosity = 1._knd/Re
         !$omp end workshare
         !$omp end parallel
       else
         !$omp parallel
         !$omp workshare
         Viscosity = 0
         !$omp end workshare
         !$omp end parallel
       end if

       if (Re>0 .and.&
             (enable_buoyancy==1.or. &
              enable_moisture==1.or. &
              num_of_scalars>0))     then
         !$omp parallel
         !$omp workshare
         TDiff = 1._knd/(Re*Prandtl)
         !$omp end workshare
         !$omp end parallel
       end if  !Re>0

       !$omp parallel
       !$omp sections
       !$omp section
       call BoundU(1,U,Uin)
       !$omp section
       call BoundU(2,V,Vin)
       !$omp section
       call BoundU(3,W,Win)
       !$omp section
       call Bound_Pr(Pr)
       !$omp end sections
       !$omp end parallel


       call Pr_Correct(U,V,W,Pr,Q,1._knd)


       if (sgstype==SubgridModel) then
                         call SGS_Smag(U,V,W,2._knd)
       elseif (sgstype==SigmaModel) then
                         call SGS_Sigma(U,V,W,2._knd)
       elseif (sgstype==VremanModel) then
                         call SGS_Vreman(U,V,W,2._knd)
       elseif (sgstype==StabSubgridModel) then
                         call SGS_StabSmag(U,V,W,Temperature,2._knd)
       else
         if (Re>0) then
           Viscosity = 1._knd/Re
         else
           Viscosity = 0
         end if
       end if

       call BoundViscosity(Viscosity)

       if (enable_buoyancy==1.or. &
           enable_moisture==1.or. &
           num_of_scalars>0)     then

         !$omp parallel
         !$omp workshare
         forall(k = 1:Prnz,j = 1:Prny,i = 1:Prnx)
           TDiff(i,j,k) = 1.35*(Viscosity(i,j,k)-1._knd/Re)+(1._knd/(Re*constPrt))
         endforall
         !$omp end workshare
         !$omp end parallel

         call BoundViscosity(TDiff)
       end if

       if (enable_buoyancy==1) then
         call BoundTemperature(Temperature)
       end if

       if (enable_moisture==1) then
         call BoundMoisture(Moisture)
       end if

       do i = 1,num_of_scalars
         call BoundScalar(Scalar(:,:,:,i))
       end do

       call InitTempFl(Temperature)

       if (wallmodeltype>0) then
                      call ComputeViscsWM(U,V,W,Pr,Temperature)
       end if

       call BoundViscosity(Viscosity)

    end if !init conditions not from file


    !prepare arrays with indexes of points to be nulled every timestep

    nUnull = 0

    !$omp parallel do reduction(+:nUnull)
    do k = 1,Unz
     do j = 1,Uny
      do i = 1,Unx
       if (Utype(i,j,k)>0.and.Utype(i,j,k+1)>0.and.Utype(i,j,k-1)>0&
           .and.Utype(i,j-1,k)>0.and.Utype(i,j+1,k)>0&
           .and.Utype(i-1,j,k)>0.and.Utype(i+1,j,k)>0)  nUnull = nUnull+1
      end do
     end do
    end do
    !$omp end parallel do


    write(*,*) "set"
  endsubroutine InitialConditions






  subroutine InitBoundaryConditions
  real(knd),allocatable:: xU2(:),yV2(:),zW2(:)
  integer i,j,k,nx,ny,nz,nxup,nxdown,nyup,nydown,nzup,nzdown,io
  real(knd) P
  integer unit

    nx = Prnx-1
    ny = Prny-1
    nz = Prnz-1

    if (xgridfromfile) then

      call newunit(unit)

      open(unit,file="xgrid.txt")
      j=-1
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

      if (Btype(Ea)==PERIODIC) then
                            Unx = Prnx
      else
                            Unx = Prnx-1
      end if

      close(unit)

    end if

    if (ygridfromfile) then

      call newunit(unit)

      open(unit,file="ygrid.txt")
      j=-1
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

      if (Btype(No)==PERIODIC) then
                            Vny = Prny
      else
                            Vny = Prny-1
      end if

      close(unit)

    end if

    if (zgridfromfile) then

      call newunit(unit)

      open(unit,file="zgrid.txt")
      j=-1
      do
        read (unit,*,iostat = io) P
        if (io==0) then
          j = j+1
          write(*,*) j
        else
          exit
        end if
      end do

      nz = j
      Prnz = nz
      Unz = Prnz
      Vnz = Prnz

      if (Btype(To)==PERIODIC) then
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
      do j = 0,nx
        write(*,*) j
        read(unit,*) xU2(j)
      end do
      close(unit)

      if (Btype(We)==PERIODIC) then
        do j=-1,-3,-1
          xU2(j) = xU2(0)-(xU2(nx)-xU2(nx+j))
        end do
      else
        do j=-1,-3,-1
          xU2(j) = xU2(0)-(xU2(0-j)-xU2(0))
        end do
      end if

      if (Btype(Ea)==PERIODIC) then
        do j = nx+1,nx+4
          xU2(j) = xU2(nx)+(xU2(j-nx)-xU2(0))
        end do
      else
        do j = nx+1,nx+4
          xU2(j) = xU2(nx)+(xU2(nx)-xU2(nx-(j-nx)))
        end do
      end if

      x0 = xU2(0)

    else

      forall (i=-3:nx+4)
         xU2(i)=(i)*dxmin+x0
      endforall

    end if


    if (ygridfromfile) then

      call newunit(unit)

      open(unit,file="ygrid.txt")
      do j = 0,ny
        read(unit,*) yV2(j)
      end do
      close(unit)

      if (Btype(So)==PERIODIC) then
        do j=-1,-3,-1
          yV2(j) = yV2(0)-(yV2(ny)-yV2(ny+j))
        end do
      else
        do j=-1,-3,-1
          yV2(j) = yV2(0)-(yV2(0-j)-yV2(0))
        end do
      end if

      if (Btype(No)==PERIODIC) then
        do j = ny+1,ny+4
          yV2(j) = yV2(ny)+(yV2(j-ny)-yV2(0))
        end do
      else
        do j = ny+1,ny+4
          yV2(j) = yV2(ny)+(yV2(ny)-yV2(ny-(j-ny)))
        end do
      end if

      y0 = yV2(0)

    else

      forall (j=-3:ny+4)
        yV2(j)=(j)*dymin+y0
      endforall

    end if


    if (zgridfromfile) then

      call newunit(unit)

      open(unit,file="zgrid.txt")
      do j = 0,nz
        read(unit,*) zW2(j)
      end do
      close(unit)

      if (Btype(Bo)==PERIODIC) then
        do j=-1,-3,-1
          zW2(j) = zW2(0)-(zW2(nz)-zW2(nz+j))
        end do
      else
        do j=-1,-3,-1
          zW2(j) = zW2(0)-(zW2(0-j)-zW2(0))
        end do
      end if

      if (Btype(To)==PERIODIC) then
        do j = nz+1,nz+4
          zW2(j) = zW2(nz)+(zW2(j-nz)-zW2(0))
        end do
      else
        do j = nz+1,nz+4
          zW2(j) = zW2(nz)+(zW2(nz)-zW2(nz-(j-nz)))
        end do
      end if

      z0 = zW2(0)

    else

       forall (k=-3:nz+4)
         zW2(k)=(k)*dzmin+z0
       endforall

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
      xPr(i)=(xU(i-1)+xU(i))/2._knd
      dxPr(i) = xU(i)-xU(i-1)
    endforall

    forall (j=-2:ny+4)
      yPr(j)=(yV(j-1)+yV(j))/2._knd
      dyPr(j) = yV(j)-yV(j-1)
    endforall

    forall (k=-2:nz+4)
      zPr(k)=(zW(k-1)+zW(k))/2._knd
      dzPr(k) = zW(k)-zW(k-1)
    endforall

    forall (i=-2:nx+3)
      dxU(i) = xPr(i+1)-xPr(i)
    endforall

    forall (j=-2:ny+3)
      dyV(j) = yPr(j+1)-yPr(j)
    endforall

    forall (k=-2:nz+3)
      dzW(k) = zPr(k+1)-zPr(k)
    endforall

    deallocate(xU2)
    deallocate(yV2)
    deallocate(zW2)


    allocate(Utype(-2:Unx+3,-2:Uny+3,-2:Unz+3))
    allocate(Vtype(-2:Vnx+3,-2:Vny+3,-2:Vnz+3))
    allocate(Wtype(-2:Wnx+3,-2:Wny+3,-2:Wnz+3))
    allocate(Prtype(0:Prnx+1,0:Prny+1,0:Prnz+1))

    Utype = 0
    Vtype = 0
    Wtype = 0
    Prtype = 0


    allocate(Uin(-2:Uny+3,-2:Unz+3),Vin(-2:Vny+3,-2:Vnz+3),Win(-2:Wny+3,-2:Wnz+3))
    Uin = 0
    Vin = 0
    Win = 0

    if (enable_buoyancy>0) allocate(TempIn(-1:Prny+2,-1:Prnz+2))

    if (enable_moisture>0) allocate(MoistIn(-1:Prny+2,-1:Prnz+2))

    select case (inlettype)
      case (ZeroInletType)
        Uin = 0
        Vin = 0
        Win = 0
      case (ShearInletType)
        call SHEARINLET(ShearInletTypeParameter)
      case (ParabolicInletType)
        call PARINLET
      case (TurbulentInletType)
        call GetTurbulentInlet
      case (FromFileInletType)
        call GetInletFromFile(start_time)
      case default
        call CONSTINLET
    endselect


    if (enable_buoyancy==1) then
       if (TempBtype(Bo)==CONSTFLUX.or.TempBtype(Bo)==DIRICHLET) then

         allocate(BsideTFlArr(-1:Prnx+2,-1:Prny+2))

         if (enable_radiation==1) then
           BsideTFlArr = 0
         else if (TempBtype(Bo)==CONSTFLUX) then
           BsideTFlArr = sideTemp(Bo)
         else
           BsideTFlArr = 0
         end if

         if (TempBtype(Bo)==DIRICHLET) then
           allocate(BsideTArr(-1:Prnx+2,-1:Prny+2))
           BsideTArr = sideTemp(Bo)
         end if

        end if
    end if

    if (.not.allocated(BsideTArr))  allocate(BsideTArr(0,0))
    if (.not.allocated(BsideTFlArr))  allocate(BsideTFlArr(0,0))


    if (enable_moisture==1) then
       if (MoistBtype(Bo)==CONSTFLUX.or.MoistBtype(Bo)==DIRICHLET) then

         allocate(BsideMFlArr(-1:Prnx+2,-1:Prny+2))

         if (enable_radiation==1) then
           BsideMFlArr = 0
         else if (MoistBtype(Bo)==CONSTFLUX) then
           BsideMFlArr = sideMoist(Bo)
         else
           BsideMFlArr = 0
         end if

         if (MoistBtype(Bo)==DIRICHLET) then
           allocate(BsideMArr(-1:Prnx+2,-1:Prny+2))
           BsideMArr = sideMoist(Bo)
         end if

        end if
    end if

    if (.not.allocated(BsideMArr))  allocate(BsideMArr(0,0))
    if (.not.allocated(BsideMFlArr))  allocate(BsideMFlArr(0,0))


    call InitSubsidenceProfile


   if (num_of_scalars>0.and.scalsourcetype==pointsource) then
        call GridCoords(scalsrci(:),scalsrcj(:),scalsrck(:),scalsrcx(:),scalsrcy(:),scalsrcz(:))
   end if

   call InitTiles(Prnx,Prny,Prnz)

   call InitSolarRadiation

   !prepare the geometry of the solid bodies
   call InitSolidBodies

   !prepare the geometry of the plant bodies
   call InitVolumeSourceBodies

   !create actual immersed boundary and wall model points, i.m. points end up in an array
   call GetSolidBodiesBC

   !add also the wall model points from domain boundaries to the list
   call GetOutsideBoundariesWM(num_of_scalars)

   !create arrays of w.m. points from the list
   call MoveWMPointsToArray

   if (enable_radiation==1) then
     !compute the radiation balance and prepare the wall fluxes 
     call InitIBPFluxes
     !set the immersed boundary values for fluxes
!      call SetIBPFluxes
   end if

   call SetNullifiedPoints

   !create actual arrays of the source points
   call InitVolumeSources

   !add arrays of the line source points
   call InitLineSources

   !add puff sources, each containing one or more points
   call InitPuffSources

   call SetFrameDomains

    write (*,*) "set"
  end subroutine InitBoundaryConditions




  subroutine SetNullifiedPoints
    integer i,j,k,n

    !$omp parallel do reduction(+:nUnull)
    do k = 1,Unz
     do j = 1,Uny
      do i = 1,Unx
       if (Utype(i,j,k)>0.and.Utype(i,j,k+1)>0.and.Utype(i,j,k-1)>0&
           .and.Utype(i,j-1,k)>0.and.Utype(i,j+1,k)>0&
           .and.Utype(i-1,j,k)>0.and.Utype(i+1,j,k)>0)  nUnull = nUnull+1
      end do
     end do
    end do
    !$omp end parallel do

    allocate(Unull(3,nUnull))

    n = 0

    do k = 1,Unz
     do j = 1,Uny
      do i = 1,Unx
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
    do k = 1,Vnz
     do j = 1,Vny
      do i = 1,Vnx
       if (Vtype(i,j,k)>0.and.Vtype(i,j,k+1)>0.and.Vtype(i,j,k-1)>0&
           .and.Vtype(i,j-1,k)>0.and.Vtype(i,j+1,k)>0&
           .and.Vtype(i-1,j,k)>0.and.Vtype(i+1,j,k)>0)  nVnull = nVnull+1
      end do
     end do
    end do
    !$omp end parallel do

    allocate(Vnull(3,nVnull))

    n = 0

    do k = 1,Vnz
     do j = 1,Vny
      do i = 1,Vnx
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
    do k = 1,Wnz
     do j = 1,Wny
      do i = 1,Wnx
       if (Wtype(i,j,k)>0.and.Wtype(i,j,k+1)>0.and.Wtype(i,j,k-1)>0&
           .and.Wtype(i,j-1,k)>0.and.Wtype(i,j+1,k)>0&
           .and.Wtype(i-1,j,k)>0.and.Wtype(i+1,j,k)>0)  nWnull = nWnull+1
      end do
     end do
    end do
    !$omp end parallel do

    allocate(Wnull(3,nWnull))

    n = 0


    do k = 1,Wnz
     do j = 1,Wny
      do i = 1,Wnx
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




  subroutine INIT_RANDOM_SEED()
    integer :: i, n, clock
    integer,dimension(:),allocatable :: seed

    call RANDOM_SEED(size = n)
    allocate(seed(n))

    call SYSTEM_CLOCK(COUNT = clock)

    seed = 0!clock+37*[(i-1,i = 1,n)]
    call RANDOM_SEED(PUT = seed)

    deallocate(seed)
  endsubroutine


endmodule Initial

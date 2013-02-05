module INITIAL

  use PARAMETERS
  use LIMITERS, only: limparam, limitertype
  use MULTIGRID, only: SetMGParams
  use MULTIGRID2d, only: SetMGParams2d
  use POISSON
  use BOUNDARIES
  use OUTPUTS, only: store, display, probes, NumProbes, SetFrameDomain, StaggeredFrameDomains
  use SCALARS
  use Filters, only: filtertype, filter_ratios
  use Subgrid
  use TURBINLET, only: GetTurbInlet, GetInletFromFile, TLag, Lturby, Lturbz, ustarinlet, transformtensor
  use GEOMETRIC
  use WALLMODELS
  use TILING, only: tilesize,InitTiles
  use FreeUnit, only: newunit

  implicit none

  private
  public  ReadParams, Initconds, ReadBounds

  real(KND) x0,y0,z0 !domain boundaries, will become xU(0), yV(0), zW(0)
  real(KND) lx,ly,lz !domain extents

contains


 subroutine ReadParams
   use StaggeredFrames, only: rrange, TFrameTimes, TSaveFlags, Init
   integer   lmg,minmglevel,bnx,bny,bnz,mgncgc,mgnpre,mgnpost,mgmaxinnerGSiter,minGPUlevel
   real(KND) mgepsinnerGS
   integer   i,io,io2,itmp
   integer numframeslices
   namelist /cmd/ tilesize, debugparam, debuglevel, windangle, projectiontype, Prnx, Prny, Prnz,&
                   obstaclefile
   character(len=1024) :: commandline,msg
   integer :: exenamelength
   integer :: unit

   type(rrange) :: range
   type(TFrameTimes) :: frame_times
   type(TSaveFlags) :: frame_save_flags
   character(10) :: domain_label
   integer :: num_staggered_domains

   call newunit(unit)

   open(unit,file="main.conf",status="old",action="read")
   read(unit,fmt='(/)')
   read(unit,*) tempmet
   read(unit,fmt='(/)')
   read(unit,*) CFL
   read(unit,fmt='(/)')
   read(unit,*) Uref
   read(unit,fmt='(/)')
   read(unit,*) poissmet
   read(unit,fmt='(/)')
   read(unit,*) convmet
   read(unit,fmt='(/)')
   read(unit,*) limitertype
   read(unit,fmt='(/)')
   read(unit,*) limparam
   read(unit,fmt='(/)')
   read(unit,*) masssourc
   read(unit,fmt='(/)')
   read(unit,*) steady
   read(unit,fmt='(/)')
   read(unit,*) tasktype
   write(*,*) "tasktype=",tasktype
   read(unit,fmt='(/)')
   read(unit,*) initcondsfromfile
   read(unit,fmt='(/)')
   read(unit,*) timeavg1
   read(unit,fmt='(/)')
   read(unit,*) timeavg2
   read(unit,fmt='(/)')
   read(unit,*) Re
   write(*,*) "Re=",Re
   read(unit,fmt='(/)')
   read(unit,*) starttime
   write(*,*) "starttime=",starttime
   read(unit,fmt='(/)')
   read(unit,*) endtime
   write(*,*) "endtime=",endtime
   read(unit,fmt='(/)')
   read(unit,*) maxiter
   write(*,*) "maxiter=",maxiter
   read(unit,fmt='(/)')
   read(unit,*) eps
   write(*,*) "eps=",eps
   read(unit,fmt='(/)')
   read(unit,*) maxCNiter
   write(*,*) "maxCNiter=",maxCNiter
   read(unit,fmt='(/)')
   read(unit,*) epsCN
   write(*,*) "epsCN=",epsCN
   read(unit,fmt='(/)')
   read(unit,*) maxPOISSONiter
   write(*,*) "maxPOISSONiter=",maxPOISSONiter
   read(unit,fmt='(/)')
   read(unit,*) epsPOISSON
   write(*,*) "epsPOISSON=",epsPOISSON
   read(unit,fmt='(/)')
   read(unit,*) debugparam
   write(*,*) "debug parameter=",debugparam
   close(unit)


   open(unit,file="les.conf",status="old",action="read")
   read(unit,fmt='(/)')
   read(unit,*) sgstype
   read(unit,fmt='(/)')
   read(unit,*) filtertype

   if (filtertype > size(filter_ratios)) then
     write(*,*) "Chosen filter type does not exist. Maximum index is:",size(filter_ratios)
     stop
   end if

   read(unit,fmt='(/)')
   read(unit,*) wallmodeltype
   close(unit)


   open(unit,file="grid.conf",status="old",action="read")
   read(unit,fmt='(/)')
   read(unit,*) xgridfromfile
   read(unit,fmt='(/)')
   read(unit,*) ygridfromfile
   read(unit,fmt='(/)')
   read(unit,*) zgridfromfile
   read(unit,fmt='(/)')

   read(unit,*) x0
   write(*,*) "x0=",x0
   read(unit,fmt='(/)')
   read(unit,*) y0
   write(*,*) "y0=",y0
   read(unit,fmt='(/)')
   read(unit,*) z0
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
   read(unit,fmt='(/)')
   read(unit,*) Btype(We)
   read(unit,fmt='(/)')
   read(unit,*) Btype(Ea)
   read(unit,fmt='(/)')
   read(unit,*) Btype(So)
   read(unit,fmt='(/)')
   read(unit,*) Btype(No)
   read(unit,fmt='(/)')
   read(unit,*) Btype(Bo)
   read(unit,fmt='(/)')
   read(unit,*) Btype(To)
   read(unit,fmt='(/)')
   read(unit,*) sideU(1,So)
   read(unit,fmt='(/)')
   read(unit,*) sideU(2,So)
   read(unit,fmt='(/)')
   read(unit,*) sideU(3,So)
   read(unit,fmt='(/)')
   read(unit,*) sideU(1,No)
   read(unit,fmt='(/)')
   read(unit,*) sideU(2,No)
   read(unit,fmt='(/)')
   read(unit,*) sideU(3,No)
   read(unit,fmt='(/)')
   read(unit,*) sideU(1,Bo)
   read(unit,fmt='(/)')
   read(unit,*) sideU(2,Bo)
   read(unit,fmt='(/)')
   read(unit,*) sideU(3,Bo)
   read(unit,fmt='(/)')
   read(unit,*) sideU(1,To)
   read(unit,fmt='(/)')
   read(unit,*) sideU(2,To)
   read(unit,fmt='(/)')
   read(unit,*) sideU(3,To)
   read(unit,fmt='(/)')
   read(unit,*) z0W
   read(unit,fmt='(/)')
   read(unit,*) z0E
   read(unit,fmt='(/)')
   read(unit,*) z0S
   read(unit,fmt='(/)')
   read(unit,*) z0N
   read(unit,fmt='(/)')
   read(unit,*) z0B
   read(unit,fmt='(/)')
   read(unit,*) z0T
   close(unit)

   open(unit,file="large_scale.conf",status="old",action="read",iostat=io)

   if (io==0) then
     read(unit,fmt='(/)')
     read(unit,*) CoriolisParam
     write(*,*) "coriolisparam=",CoriolisParam
     read(unit,fmt='(/)')
     read(unit,*) PrGradientX
     write(*,*) "prgradientx=",PrGradientX
     read(unit,fmt='(/)')
     read(unit,*) PrGradientY
     write(*,*) "prgradienty=",PrGradientY
     read(unit,fmt='(/)')
     read(unit,*) SubsidenceGradient
     write(*,*) "SubsidenceGradient=",SubsidenceGradient
     close(unit)

   else

     write(*,*) "Warning! Could not open file large_scale.conf. Using defaults."

   endif


   open(unit,file="thermal.conf",status="old",action="read")
   read(unit,fmt='(/)')
   read(unit,*) buoyancy
   read(unit,fmt='(/)')
   read(unit,*) Prandtl
   read(unit,fmt='(/)')
   read(unit,*) grav_acc
   read(unit,fmt='(/)')
   read(unit,*) temperature_ref
   read(unit,fmt='(/)')
   read(unit,*) TBtype(We)
   read(unit,fmt='(/)')
   read(unit,*) TBtype(Ea)
   read(unit,fmt='(/)')
   read(unit,*) TBtype(So)
   read(unit,fmt='(/)')
   read(unit,*) TBtype(No)
   read(unit,fmt='(/)')
   read(unit,*) TBtype(Bo)
   read(unit,fmt='(/)')
   read(unit,*) TBtype(To)
   read(unit,fmt='(/)')
   read(unit,*) sideTemp(We)
   read(unit,fmt='(/)')
   read(unit,*) sideTemp(Ea)
   read(unit,fmt='(/)')
   read(unit,*) sideTemp(So)
   read(unit,fmt='(/)')
   read(unit,*) sideTemp(No)
   read(unit,fmt='(/)')
   read(unit,*) sideTemp(Bo)
   read(unit,fmt='(/)')
   read(unit,*) sideTemp(To)
   close(unit)

   if (buoyancy==1) then

     open(unit,file="temp_profile.conf",status="old",action="read",iostat=io)

     if (io==0) then
       read(unit,fmt='(/)')
       read(unit,*) TemperatureProfile%randomize
       read(unit,fmt='(/)')
       read(unit,*) TemperatureProfile%randomizeTop
       read(unit,fmt='(/)')
       read(unit,*) TemperatureProfile%randomizeAmplitude
       read(unit,fmt='(/)')
       read(unit,*) itmp

       allocate(TemperatureProfile%Sections(max(itmp,0)))

       do i = 1, size(TemperatureProfile%Sections)
         read(unit,fmt='(/)')
         read(unit,*) TemperatureProfile%Sections(i)%top
         read(unit,fmt='(/)')
         read(unit,*) TemperatureProfile%Sections(i)%jump
         read(unit,fmt='(/)')
         read(unit,*) TemperatureProfile%Sections(i)%gradient
       enddo

       close(unit)

     else

       write(*,*) "Warning! Could not open file temp_profile.conf. Using defaults."
       TemperatureProfile%randomize = 0

       allocate(TemperatureProfile%Sections(0))

     endif

   else

     TBtype = 0

   endif

   open(unit,file="inlet.conf",status="old",action="read")
   read(unit,fmt='(/)')
   read(unit,*) inlettype
   read(unit,fmt='(/)')
   read(unit,*) profiletype
   read(unit,fmt='(/)')
   read(unit,*) SHEARG
   write(*,*) "G=",SHEARG
   read(unit,fmt='(/)')
   read(unit,*) Uinlet
   write(*,*) "Uinlet=",Uinlet
   read(unit,fmt='(/)')
   read(unit,*) ustarsurfin  !-<u'w'>
   read(unit,fmt='(/)')
   read(unit,*) stressgradin !in relative part per 1m
   read(unit,fmt='(/)')
   read(unit,*) z0inlet
   read(unit,fmt='(/)')
   read(unit,*) powerexpin
   read(unit,fmt='(/)')
   read(unit,*) zrefin
   read(unit,fmt='(/)')
   read(unit,*) Urefin
   read(unit,fmt='(/)')
   read(unit,*)  relativestress(1,1)
   read(unit,fmt='(/)')
   read(unit,*)  relativestress(2,2)
   read(unit,fmt='(/)')
   read(unit,*) relativestress(3,3)
   read(unit,fmt='(/)')
   read(unit,*) relativestress(1,2)
   read(unit,fmt='(/)')
   read(unit,*) relativestress(1,3)
   read(unit,fmt='(/)')
   read(unit,*) relativestress(2,3)
   read(unit,fmt='(/)')
   read(unit,*) TLag
   read(unit,fmt='(/)')
   read(unit,*) Lturby
   read(unit,fmt='(/)')
   read(unit,*) Lturbz
   close(unit)

   relativestress(2,1)=relativestress(1,2)
   relativestress(3,1)=relativestress(1,3)
   relativestress(3,2)=relativestress(2,3)

   open(unit,file="scalars.conf",status="old",action="read")
   read(unit,fmt='(/)')
   read(unit,*) computescalars
   read(unit,fmt='(/)')
   read(unit,*) computedeposition
   read(unit,fmt='(/)')
   read(unit,*) computegravsettling
   read(unit,fmt='(/)')
   read(unit,*) partdistrib
   read(unit,fmt='(/)')
   read(unit,*) totalscalsource
   read(unit,fmt='(/)')
   read(unit,*) scalsourcetype

   read(unit,fmt='(/)')
   read(unit,*) ScalBtype(We)
   read(unit,fmt='(/)')
   read(unit,*) ScalBtype(Ea)
   read(unit,fmt='(/)')
   read(unit,*) ScalBtype(So)
   read(unit,fmt='(/)')
   read(unit,*) ScalBtype(No)
   read(unit,fmt='(/)')
   read(unit,*) ScalBtype(Bo)
   read(unit,fmt='(/)')
   read(unit,*) ScalBtype(To)
   read(unit,fmt='(/)')
   read(unit,*) sideScal(We)
   read(unit,fmt='(/)')
   read(unit,*) sideScal(Ea)
   read(unit,fmt='(/)')
   read(unit,*) sideScal(So)
   read(unit,fmt='(/)')
   read(unit,*) sideScal(No)
   read(unit,fmt='(/)')
   read(unit,*) sideScal(Bo)
   read(unit,fmt='(/)')
   read(unit,*) sideScal(To)

   if (partdistrib>0) then

      allocate(partdiam(partdistrib),partrho(partdistrib),percdistrib(partdistrib))
      allocate(scalsrcx(partdistrib),scalsrcy(partdistrib),scalsrcz(partdistrib))
      allocate(scalsrci(partdistrib),scalsrcj(partdistrib),scalsrck(partdistrib))

      do i=1,partdistrib
        read(unit,fmt='(/)')
        read(unit,*) partdiam(i)
        read(unit,fmt='(/)')
        read(unit,*) partrho(i)
        read(unit,fmt='(/)')
        read(unit,*) percdistrib(i)
        read(unit,fmt='(/)')
        read(unit,*) scalsrcx(i)
        read(unit,fmt='(/)')
        read(unit,*) scalsrcy(i)
        read(unit,fmt='(/)')
        read(unit,*) scalsrcz(i)
      enddo

   else

      allocate(partdiam(computescalars),partrho(computescalars),percdistrib(computescalars))
      allocate(scalsrcx(computescalars),scalsrcy(computescalars),scalsrcz(computescalars))
      allocate(scalsrci(computescalars),scalsrcj(computescalars),scalsrck(computescalars))

      do i=1,computescalars
        read(unit,fmt='(/)')
        read(unit,*) partdiam(i)
        read(unit,fmt='(/)')
        read(unit,*) partrho(i)
        read(unit,fmt='(/)')
        read(unit,*) percdistrib(i)
        read(unit,fmt='(/)')
        read(unit,*) scalsrcx(i)
        read(unit,fmt='(/)')
        read(unit,*) scalsrcy(i)
        read(unit,fmt='(/)')
        read(unit,*) scalsrcz(i)
      enddo

   endif
   close(unit)

   if (poissmet==3.or.poissmet==4.or.poissmet==5) then
     open(unit,file="mgopts.conf",status="old",action="read")
     read(unit,fmt='(/)')
     read(unit,*) lmg
     read(unit,fmt='(/)')
     read(unit,*) minmglevel
     read(unit,fmt='(/)')
     read(unit,*) bnx
     read(unit,fmt='(/)')
     read(unit,*) bny
     read(unit,fmt='(/)')
     read(unit,*) bnz
     read(unit,fmt='(/)')
     read(unit,*) mgncgc
     read(unit,fmt='(/)')
     read(unit,*) mgnpre
     read(unit,fmt='(/)')
     read(unit,*) mgnpost
     read(unit,fmt='(/)')
     read(unit,*) mgmaxinnerGSiter
     read(unit,fmt='(/)')
     read(unit,*) mgepsinnerGS
     read(unit,fmt='(/)',iostat=io)
     read(unit,*,iostat=io) minGPUlevel
     close(unit)

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
     endif
   endif

   open(unit,file="frames.conf",status="old",action="read")
   read(unit,fmt='(/)')
   read(unit,*) frames
   read(unit,fmt='(/)')
   read(unit,*) timefram1
   read(unit,fmt='(/)')
   read(unit,*) timefram2
   read(unit,fmt='(/)')

   read(unit,*) store%frame_U
   read(unit,fmt='(/)')
   read(unit,*) store%frame_vort
   read(unit,fmt='(/)')
   read(unit,*) store%frame_Pr
   read(unit,fmt='(/)')
   read(unit,*) store%frame_lambda2
   read(unit,fmt='(/)')
   read(unit,*) store%frame_scalars
   read(unit,fmt='(/)')
   read(unit,*) store%frame_sumscalars
   read(unit,fmt='(/)')
   read(unit,*) store%frame_T
   read(unit,fmt='(/)')
   read(unit,*) store%frame_tempfl
   read(unit,fmt='(/)')
   read(unit,*) store%frame_scalfl


   read(unit,fmt='(/)')
   read(unit,*) numframeslices

   allocate(store%frame_domains(numframeslices))

   do i=1,numframeslices
     read(unit,fmt='(/)')
     read(unit,*) store%frame_domains(i)%dimension
     read(unit,fmt='(/)')
     read(unit,*) store%frame_domains(i)%direction
     read(unit,fmt='(/)')
     read(unit,*) store%frame_domains(i)%position
   enddo

   close(unit)


   open(unit,file="stagframes.conf",status="old",action="read",iostat=io)
   if (io==0) then
     read(unit,fmt='(/)')
     read(unit,*) num_staggered_domains

     allocate(StaggeredFrameDomains(num_staggered_domains))

     do i=1,num_staggered_domains
       read(unit,fmt='(//)')
       read(unit,*) domain_label
       read(unit,fmt='(/)')
       read(unit,*) range%min%x,range%max%x
       read(unit,fmt='(/)')
       read(unit,*) range%min%y,range%max%y
       read(unit,fmt='(/)')
       read(unit,*) range%min%z,range%max%z
       read(unit,fmt='(/)')
       read(unit,*) frame_times%nframes
       read(unit,fmt='(/)')
       read(unit,*) frame_times%start, frame_times%end
       read(unit,fmt='(/)')
       read(unit,*) frame_save_flags%U, frame_save_flags%V, frame_save_flags%W
       read(unit,fmt='(/)')
       read(unit,*) frame_save_flags%Pr
       read(unit,fmt='(/)')
       read(unit,*) frame_save_flags%Temperature
       if (buoyancy == 0) frame_save_flags%Temperature = .false.
       read(unit,fmt='(/)')
       read(unit,*) frame_save_flags%Scalar
       if (computescalars == 0) frame_save_flags%Scalar = .false.

       call Init(StaggeredFrameDomains(i), trim(domain_label), &
                 range, &
                 frame_times, &
                 frame_save_flags )
     end do
     close(unit)
   end if




   open(unit,file="output.conf",status="old",action="read",iostat=io)
   if (io==0) then
     read(unit,fmt='(/)')
     read(unit,*) display%delta
     read(unit,fmt='(/)')
     read(unit,*) display%ustar
     read(unit,fmt='(/)')
     read(unit,*) display%tstar

     read(unit,fmt='(/)')
     read(unit,*) store%U
     read(unit,fmt='(/)')
     read(unit,*) store%U_interp
     read(unit,fmt='(/)')
     read(unit,*) store%V
     read(unit,fmt='(/)')
     read(unit,*) store%V_interp
     read(unit,fmt='(/)')
     read(unit,*) store%W
     read(unit,fmt='(/)')
     read(unit,*) store%W_interp

     read(unit,fmt='(/)')
     read(unit,*) store%out

     read(unit,fmt='(/)')
     read(unit,*) store%out_U
     read(unit,fmt='(/)')
     read(unit,*) store%out_vort
     read(unit,fmt='(/)')
     read(unit,*) store%out_Pr
     read(unit,fmt='(/)')
     read(unit,*) store%out_Prtype
     read(unit,fmt='(/)')
     read(unit,*) store%out_lambda2
     read(unit,fmt='(/)')
     read(unit,*) store%out_T
     read(unit,fmt='(/)')
     read(unit,*) store%out_div
     read(unit,fmt='(/)')
     read(unit,*) store%out_visc

     read(unit,fmt='(/)')
     read(unit,*) store%avg

     read(unit,fmt='(/)')
     read(unit,*) store%avg_U
     read(unit,fmt='(/)')
     read(unit,*) store%avg_vort
     read(unit,fmt='(/)')
     read(unit,*) store%avg_Pr
     read(unit,fmt='(/)')
     read(unit,*) store%avg_Prtype
     read(unit,fmt='(/)')
     read(unit,*) store%avg_T

     read(unit,fmt='(/)')
     read(unit,*) store%scalars
     read(unit,fmt='(/)')
     read(unit,*) store%scalarsavg

     read(unit,fmt='(/)')
     read(unit,*) store%deposition

     read(unit,fmt='(/)')
     read(unit,*) store%deltime
     read(unit,fmt='(/)')
     read(unit,*) store%tke
     read(unit,fmt='(/)')
     read(unit,*) store%dissip
     read(unit,fmt='(/)')
     read(unit,*) store%scalsumtime
     read(unit,fmt='(/)')
     read(unit,*) store%scaltotsumtime
     read(unit,fmt='(/)')
     read(unit,*) store%ustar
     read(unit,fmt='(/)')
     read(unit,*) store%tstar
     read(unit,fmt='(/)')
     read(unit,*) store%blprofiles

     read(unit,fmt='(/)')
     read(unit,*) NumProbes
     allocate(probes(Numprobes))

     do i=1,NumProbes
       read(unit,fmt='(/)')
       read(unit,*) probes(i)%x
       read(unit,fmt='(/)')
       read(unit,*) probes(i)%y
       read(unit,fmt='(/)')
       read(unit,*) probes(i)%z
     enddo
     close(unit)
   endif

   open(unit,file="obstacles.conf",status="old",action="read",iostat=io)
   if (io==0) then
     read(unit,fmt='(/)')
     read(unit,'(a)') obstaclefile
     close(unit)
   end if

   write(*,*) "computescalars",computescalars
   write(*,*) "partdiam",partdiam




   windangle=0._KND

   projectiontype=1

#ifndef NO_INTERNAL_NML
   !Parsing of command line uses namelist input on an internal file (Fortran 2003).
   !If not supported by the processor yet, define macro NO_INTERNAL_NML.
   call get_command(command=commandline,status=io)
   call get_command_argument(0,length=exenamelength,status=io2)
   if (io2==0) then
     commandline="&cmd "//adjustl(trim(commandline(exenamelength+1:)))//" /"
   else
     commandline="&cmd "//adjustl(trim(commandline))//" /"
   end if

   if (io==0) then
     msg = ''
     read(commandline,nml=cmd,iostat=io,iomsg=msg)
     if (io/=0) then
       write(*,*) io,"Error parsing command line."
       write(*,*) msg
       write(*,*) commandline
     end if
   else
     write(*,*) io,"Error getting command line."
   end if
#endif

   if (CFL<=0)  CFL=0.5


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





   dxmin=lx/(Prnx)
   dymin=ly/(Prny)
   dzmin=lz/(Prnz)

   if (Prnx==1) then
     dxmin=sqrt(dymin*dzmin)
     lx=dxmin
   elseif (Prny==1) then
     dymin=sqrt(dxmin*dzmin)
     ly=dymin
   elseif (Prnz==1) then
     dzmin=sqrt(dxmin*dymin)
     lz=dzmin
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
 integer i,j,k,unit

   call newunit(unit)

   open(unit,file="in.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read(unit,*)
   enddo
   do k=1,Prnz
    do j=1,Prny
     do i=1,Prnx
      read(unit,*) Pr(i,j,k)
     enddo
    enddo
   enddo
   if (buoyancy==1) then
    do i=1,3
     read(unit,*)
    enddo
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
       read(unit,*) temperature(i,j,k)
      enddo
     enddo
    enddo
   endif
   close(unit)

   open(unit,file="Uin.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read(unit,*)
   enddo
   do k=1,Unz
    do j=1,Uny
     do i=1,Unx
      read(unit,*) U(i,j,k)
     enddo
    enddo
   enddo
   close(unit)

   open(unit,file="Vin.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read(unit,*)
   enddo
   do k=1,Vnz
    do j=1,Vny
     do i=1,Vnx
      read(unit,*) V(i,j,k)
     enddo
    enddo
   enddo
   close(unit)

   open(unit,file="Win.vtk",position="rewind",status="old",action="read")
   do i=1,14
    read(unit,*)
   enddo
   do k=1,Wnz
    do j=1,Wny
     do i=1,Wnx
      read(unit,*) W(i,j,k)
     enddo
    enddo
   enddo
 endsubroutine ReadIC



  subroutine INITCONDS(U,V,W,Pr)
  real(KND) U(-2:,-2:,-2:),V(-2:,-2:,-2:),W(-2:,-2:,-2:),Pr(1:,1:,1:)
  integer i,j,k
  real(KND) p,x,y,z,x1,x2,y1,y2,z1,z2
  real(KND),allocatable :: Q(:,:,:)

    call init_random_seed

    Pr(1:Prnx,1:Prny,1:Prnz)=0

    U=huge(1._KND)
    V=huge(1._KND)
    W=huge(1._KND)

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

    else   !init conditions not from file

       if (tasktype==2) then
         U(1:Unx,1:Uny,1:Unz)=0
         do k=1,Unz
          do j=1,Uny
           do i=1,Unx
            !if (Utype(i,j,k)<=0) then
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
            !if (Utype(i,j,k)<=0) then
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

       elseif (Btype(We)==TURBULENTINLET.or.Btype(Ea)==TURBULENTINLET) then

         do i=1,Prnx

           call GetTurbInlet
           !$omp parallel private(j,k)
           !$omp do
           do k=1,Unz
            do j=1,Uny
              if (Utype(i,j,k)<=0) then
                    U(i,j,k)=Uin(j,k)
              else
                 U(i,j,k)=0
              endif
            enddo
           enddo
           !$omp end do nowait
           !$omp do
           do k=1,Vnz
            do j=1,Vny
              if (Vtype(i,j,k)<=0) then
                    V(i,j,k)=Vin(j,k)
              else
                 V(i,j,k)=0
              endif
            enddo
           enddo
           !$omp end do nowait
           !$omp do
           do k=1,Wnz
            do j=1,Wny
              if (Wtype(i,j,k)<=0) then
                    W(i,j,k)=Win(j,k)
              else
                 W(i,j,k)=0
              endif
            enddo
           enddo
           !$omp end do
           !$omp end parallel
         enddo

       else

!          !$omp parallel private(i,j,k)
!          !$omp do
         do k=1,Unz
          do j=1,Uny
           do i=1,Unx
            if (Utype(i,j,k)<=0) then
                  call RANDOM_NUMBER(p)
                  U(i,j,k)=Uin(j,k)!+(Sqrt((Uin(j,k))**2+(Vin(j,k))**2))*0.1_KND*(p-0.5_KND)!sin(2.*pi*xU(i)+1)*cos(2.*pi*yPr(j)-2)!Uin(j,k)!*(1+0.03_KND*(p-0.5_KND))
             else
               U(i,j,k)=0
            endif
           enddo
          enddo
         enddo
!          !$omp end do
!          !$omp do
         do k=1,Vnz
          do j=1,Vny
           do i=1,Vnx
            if (Vtype(i,j,k)<=0) then
                  call RANDOM_NUMBER(p)
                  V(i,j,k)=Vin(j,k)!+(Sqrt((Uin(j,k))**2+(Vin(j,k))**2))*0.1_KND*(p-0.5_KND)!-cos(2.*pi*xPr(i)+1)*sin(2.*pi*yV(j)-2)!Uinlet*(0.3_KND*(p-0.5_KND))
             else
               V(i,j,k)=0
            endif
           enddo
          enddo
         enddo
!          !$omp end do
!          !$omp do
         do k=1,Wnz
          do j=1,Wny
           do i=1,Wnx
            if (Wtype(i,j,k)<=0) then
                  call RANDOM_NUMBER(p)
                  W(i,j,k)=Win(j,k)!Uinlet*(0.00001_KND*(p-0.5_KND))
             else
               W(i,j,k)=0
            endif
           enddo
          enddo
         enddo
!          !$omp end do
!          !$omp end parallel
       endif  !tasktype



       if (computescalars>0) then
         !$omp parallel
         !$omp workshare
         SCALAR(1:Prnx,1:Prny,1:Prnz,:)=0
         !$omp end workshare
         !$omp end parallel
       end if

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

         call InitTemperatureProfile(TempIn)

         call InitTemperature(TempIn,Temperature)

         !$omp parallel private(i,j,k)
         !$omp workshare
         Pr(:,:,1)=0
         !$omp end workshare
         !$omp end parallel
         do k=2,Prnz
          !$omp parallel
          !$omp do
          do j=1,Vny+1
           do i=1,Unx+1
            Pr(i,j,k)=Pr(i,j,k-1) + &
                   grav_acc*dzW(k-1) * &
                  ( (temperature(i,j,k-1)+temperature(i,j,k))/2._KND - temperature_ref )&
                  / temperature_ref
           enddo
          enddo
          !$omp end do
          !$omp end parallel
         enddo

       endif !byoyancy and tasktype

       if (Re>0) then
         !$omp parallel
         !$omp workshare
         Visc=1._KND/Re
         !$omp end workshare
         !$omp end parallel
       else
         !$omp parallel
         !$omp workshare
         Visc=0
         !$omp end workshare
         !$omp end parallel
       endif

       if (Re>0.and.(buoyancy==1.or.computescalars>0)) then
         !$omp parallel
         !$omp workshare
         TDiff=1._KND/(Re*Prandtl)
         !$omp end workshare
         !$omp end parallel
       endif  !Re>0

       !$omp parallel
       !$omp sections
       !$omp section
       call BoundU(1,U)
       !$omp section
       call BoundU(2,V)
       !$omp section
       call BoundU(3,W)
       !$omp section
       call Bound_Pr(Pr)
       !$omp end sections
       !$omp end parallel


       call Pr_Correct(U,V,W,Pr,Q,1._KND)


       if (sgstype==SubgridModel) then
                         call SGS_Smag(U,V,W,2._KND)
       elseif (sgstype==SigmaModel) then
                         call SGS_Sigma(U,V,W,2._KND)
       elseif (sgstype==VremanModel) then
                         call SGS_Vreman(U,V,W,2._KND)
       elseif (sgstype==StabSubgridModel) then
                         call SGS_StabSmag(U,V,W,2._KND)
       else
         if (Re>0) then
           Visc=1._KND/Re
         else
           Visc=0
         endif
       endif

       call Bound_Visc(Visc)

       if (buoyancy>0) then
         !$omp parallel
         !$omp workshare
         forall(k=1:Prnz,j=1:Prny,i=1:Prnx)
           TDiff(i,j,k)=1.35*(Visc(i,j,k)-1._KND/Re)+(1._KND/(Re*constPrt))
         endforall
         !$omp end workshare
         !$omp end parallel

         call Bound_Visc(TDiff)
         call BoundTemperature(temperature)
       endif

       do i=1,computescalars
         call Bound_PassScalar(Scalar(:,:,:,i))
       enddo

       call InitTempFL

       if (wallmodeltype>0) then
                      call ComputeViscsWM(U,V,W,Pr)
       endif

       call Bound_Visc(Visc)

    endif !init conditions not from file


    !prepare arrays with indexes of points to be nulled every timestep

    nUnull=0

    !$omp parallel do reduction(+:nUnull)
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
       if (Utype(i,j,k)>0.and.Utype(i,j,k+1)>0.and.Utype(i,j,k-1)>0&
           .and.Utype(i,j-1,k)>0.and.Utype(i,j+1,k)>0&
           .and.Utype(i-1,j,k)>0.and.Utype(i+1,j,k)>0)  nUnull=nUnull+1
      enddo
     enddo
    enddo
    !$omp end parallel do


    write(*,*) "set"
  endsubroutine INITCONDS






  subroutine READBOUNDS
  real(KND),allocatable:: xU2(:),yV2(:),zW2(:)
  integer i,j,k,nx,ny,nz,nxup,nxdown,nyup,nydown,nzup,nzdown,io
  real(KND) P
  integer unit

    nx=Prnx-1
    ny=Prny-1
    nz=Prnz-1

    if (xgridfromfile) then

      call newunit(unit)

      open(unit,file="xgrid.txt")
      j=-1
      do
        read (unit,*,iostat=io) P
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

      close(unit)

    endif

    if (ygridfromfile) then

      call newunit(unit)

      open(unit,file="ygrid.txt")
      j=-1
      do
        read (unit,*,iostat=io) P
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

      close(unit)

    endif

    if (zgridfromfile) then

      call newunit(unit)

      open(unit,file="zgrid.txt")
      j=-1
      do
        read (unit,*,iostat=io) P
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

      close(unit)

    endif





    allocate(xU2(-3:nx+4))
    allocate(yV2(-3:ny+4))
    allocate(zW2(-3:nz+4))



    if (xgridfromfile) then

      call newunit(unit)

      open(unit,file="xgrid.txt")
      do j=0,nx
        write(*,*) j
        read(unit,*) xU2(j)
      enddo
      close(unit)

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

      call newunit(unit)

      open(unit,file="ygrid.txt")
      do j=0,ny
        read(unit,*) yV2(j)
      enddo
      close(unit)

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

      call newunit(unit)

      open(unit,file="zgrid.txt")
      do j=0,nz
        read(unit,*) zW2(j)
      enddo
      close(unit)

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


    allocate(Uin(-2:Uny+3,-2:Unz+3),Vin(-2:Vny+3,-2:Vnz+3),Win(-2:Wny+3,-2:Wnz+3))
    Uin=huge(1._KND)
    Vin=huge(1._KND)
    Win=huge(1._KND)

    if (buoyancy>0) allocate(Tempin(-1:Prny+2,-1:Prnz+2))

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

    if (.not.allocated(BsideTArr))  allocate(BsideTArr(0,0))
    if (.not.allocated(BsideTFLArr))  allocate(BsideTFLArr(0,0))


    call InitSubsidenceProfile


   if (computescalars>0.and.scalsourcetype==pointsource) then
        call Gridcoords(scalsrci(:),scalsrcj(:),scalsrck(:),scalsrcx(:),scalsrcy(:),scalsrcz(:))
   endif

   call InitTiles(Prnx,Prny,Prnz)

   call InitSolidBodies

   call GetSolidBodiesBC

   call GetOutsideBoundariesWM(computescalars)

   call MoveWMPointsToArray

   call SetNullifiedPoints

   call SetFrameDomain(store%frame_domains)

    write (*,*) "set"
  end subroutine READBOUNDS




  subroutine SetNullifiedPoints
    integer i,j,k,n

    !$omp parallel do reduction(+:nUnull)
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
       if (Utype(i,j,k)>0.and.Utype(i,j,k+1)>0.and.Utype(i,j,k-1)>0&
           .and.Utype(i,j-1,k)>0.and.Utype(i,j+1,k)>0&
           .and.Utype(i-1,j,k)>0.and.Utype(i+1,j,k)>0)  nUnull=nUnull+1
      enddo
     enddo
    enddo
    !$omp end parallel do

    allocate(Unull(3,nUnull))

    n=0

    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
       if (Utype(i,j,k)>0.and.Utype(i,j,k+1)>0.and.Utype(i,j,k-1)>0&
           .and.Utype(i,j-1,k)>0.and.Utype(i,j+1,k)>0&
           .and.Utype(i-1,j,k)>0.and.Utype(i+1,j,k)>0)  then
            !$omp atomic
            n=n+1
            Unull(:,n)=(/ i,j,k /)

       endif
      enddo
     enddo
    enddo

    nVnull=0

    !$omp parallel do reduction(+:nVnull)
    do k=1,Vnz
     do j=1,Vny
      do i=1,Vnx
       if (Vtype(i,j,k)>0.and.Vtype(i,j,k+1)>0.and.Vtype(i,j,k-1)>0&
           .and.Vtype(i,j-1,k)>0.and.Vtype(i,j+1,k)>0&
           .and.Vtype(i-1,j,k)>0.and.Vtype(i+1,j,k)>0)  nVnull=nVnull+1
      enddo
     enddo
    enddo
    !$omp end parallel do

    allocate(Vnull(3,nVnull))

    n=0

    do k=1,Vnz
     do j=1,Vny
      do i=1,Vnx
       if (Vtype(i,j,k)>0.and.Vtype(i,j,k+1)>0.and.Vtype(i,j,k-1)>0&
           .and.Vtype(i,j-1,k)>0.and.Vtype(i,j+1,k)>0&
           .and.Vtype(i-1,j,k)>0.and.Vtype(i+1,j,k)>0)  then

            n=n+1
            Vnull(:,n)=(/ i,j,k /)

       endif
      enddo
     enddo
    enddo

    nWnull=0

    !$omp parallel do reduction(+:nWnull)
    do k=1,Wnz
     do j=1,Wny
      do i=1,Wnx
       if (Wtype(i,j,k)>0.and.Wtype(i,j,k+1)>0.and.Wtype(i,j,k-1)>0&
           .and.Wtype(i,j-1,k)>0.and.Wtype(i,j+1,k)>0&
           .and.Wtype(i-1,j,k)>0.and.Wtype(i+1,j,k)>0)  nWnull=nWnull+1
      enddo
     enddo
    enddo
    !$omp end parallel do

    allocate(Wnull(3,nWnull))

    n=0


    do k=1,Wnz
     do j=1,Wny
      do i=1,Wnx
       if (Wtype(i,j,k)>0.and.Wtype(i,j,k+1)>0.and.Wtype(i,j,k-1)>0&
           .and.Wtype(i,j-1,k)>0.and.Wtype(i,j+1,k)>0&
           .and.Wtype(i-1,j,k)>0.and.Wtype(i+1,j,k)>0)  then

            n=n+1
            Wnull(:,n)=(/ i,j,k /)

       endif
      enddo
     enddo
    enddo

  end subroutine SetNullifiedPoints




  subroutine INIT_RANDOM_SEED()
  integer:: i, n, clock
  integer,dimension(:),allocatable:: seed

    call RANDOM_SEED(size = n)
    allocate(seed(n))

    call SYSTEM_CLOCK(COUNT=clock)

    seed=0!clock+37*(/(i-1,i=1,n)/)
    call RANDOM_SEED(PUT=seed)

    deallocate(seed)
  endsubroutine


endmodule INITIAL

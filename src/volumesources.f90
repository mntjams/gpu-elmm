module VolumeSources
  use Parameters
  use Lists
  use Body_class

  implicit none
! 
!   private
! 
!   public UResistanceVolumes, VResistanceVolumes, WResistanceVolumes, &
!          TemperatureFlVolumes, MoistureFlVolumes, ScalarFlVolumes, VolumeSourceBody

  type :: GridVolume
    integer   :: xi, yj, zk
  end type GridVolume

  type,extends(GridVolume) :: FluxVolume
    real(knd) :: flux
  end type FluxVolume
  
  interface FluxVolume
    module procedure FluxVolume_v3
  end interface

  type,extends(FluxVolume) :: TResistanceVolume
    ! f_drag_i = - Cd * a * V * u_i , V = sqrt(sum(u_i**2))
    !flux = Cd * a
  end type TResistanceVolume

  type,extends(FluxVolume) :: TemperatureFlVolume
    ! <T'w'> = temperature_flux, Qh = rho*Cp*temperature_flux
    !flux = temperature_flux
  end type TemperatureFlVolume

  type,extends(FluxVolume) :: MoistureFlVolume
    ! <q'w'> = moisture_flux, Qe = rho*Lv*moisture_flux
    !flux = moisture_flux
  end type MoistureFlVolume

  type,extends(FluxVolume) :: ScalarFlVolume
    ! <c'w'> = scalar_flux
    !flux = flux
  end type ScalarFlVolume

  interface ScalarFlVolume
    module procedure ScalarFlVolume_3i
    module procedure ScalarFlVolume_v3
  end interface

  type ScalarFlVolumesContainer
    integer :: scalar_number
    type(ScalarFlVolume), allocatable :: Volumes(:)
  end type

  interface Add
    module procedure AddScalarFlVolumesContainer
  end interface

  type ScalarFlVolumesListContainer
    type(List) :: list
  end type
  
  type, extends(Body) :: VolumeSourceBody
    procedure(resistance_interface),pointer  :: get_resistance       => null()
    procedure(resistance_interface),pointer  :: get_temperature_flux => null()
    procedure(resistance_interface),pointer  :: get_moisture_flux    => null()
    procedure(scalar_flux_interface),pointer :: get_scalar_flux      => null()
  end type VolumeSourceBody

  abstract interface
    function resistance_interface(self,x,y,z) result(res)
      import
      real(knd) :: res
      class(VolumeSourceBody),intent(in) :: self
      real(knd),intent(in) :: x,y,z
    end function
    function scalar_flux_interface(self,x,y,z,num_of_scalar) result(res)
      import
      real(knd) :: res
      class(VolumeSourceBody),intent(in) :: self
      real(knd),intent(in) :: x,y,z
      integer,intent(in) :: num_of_scalar
    end function
  end interface

  type(List) :: SourceBodiesList

  type(List) :: UResistanceVolumesList, &
                 VResistanceVolumesList, &
                 WResistanceVolumesList, &
                 TemperatureFlVolumesList, &
                 MoistureFlVolumesList
                 
  type(ScalarFlVolumesListContainer),allocatable :: ScalarFlVolumesLists(:)

  !final volume momentum sinks - resistances
  type(TResistanceVolume),allocatable :: UResistanceVolumes(:), &
                                         VResistanceVolumes(:), &
                                         WResistanceVolumes(:)

  !final scalar quantities sources/sinks
  type(TemperatureFlVolume),allocatable :: TemperatureFlVolumes(:)
  type(MoistureFlVolume)   ,allocatable :: MoistureFlVolumes(:)
  type(ScalarFlVolumesContainer),allocatable :: ScalarFlVolumes(:)
  
  type, extends(VolumeSourceBody) :: PlantBody
    integer :: plant_type !problem specific, used by custom routines
    real(knd) :: albedo = 0.3_knd
    real(knd) :: emmissivity = 0.7_knd
    real(knd) :: evaporative_fraction = 0.6_knd
  end type PlantBody
  
  contains

      subroutine InitVolumeSourceBodies
#ifdef CUSTOMPB
        interface
          subroutine CustomVolumeSourceBodies
          end subroutine
        end interface

      call CustomVolumeSourceBodies
#endif
    end subroutine InitVolumeSourceBodies
    
    
    function FluxVolume_v3(pos,flux) result(res)
      type(FluxVolume) :: res
      integer,intent(in) :: pos(3)
      real(knd),intent(in) :: flux
      res%xi = pos(1)
      res%yj = pos(2)
      res%zk = pos(3)
      res%flux = flux
    end function

    
    function ScalarFlVolume_v3(pos,flux) result(res)
      type(ScalarFlVolume) :: res
      integer,intent(in) :: pos(3)
      real(knd),intent(in) :: flux
      res%FluxVolume = FluxVolume(pos,flux)
    end function

    
    function ScalarFlVolume_3i(x,y,z,flux) result(res)
      type(ScalarFlVolume) :: res
      integer,intent(in) :: x,y,z
      real(knd),intent(in) :: flux
      res%FluxVolume = FluxVolume([x,y,z],flux)
    end function

    
    subroutine InsideCellsToLists
      !find cells inside the canopy and store them in a list
      logical,allocatable :: InsidePr(:,:,:)
      
!       !$omp parallel sections
!       !$omp section
      call SourceBodiesList%for_each(GetUCells)

!       !$omp section
      call SourceBodiesList%for_each(GetVCells)

!       !$omp section
      call SourceBodiesList%for_each(GetWCells)

!       !$omp section
      allocate(ScalarFlVolumesLists(num_of_scalars))

      allocate(InsidePr(-1:Prnx+2,-1:Prny+2,-1:Prnz+2))
      
      call SourceBodiesList%for_each(GetOtherCells)
!       !$omp end parallel sections

      contains

        subroutine GetUCells(PB)
          class(*) :: PB
          type(TResistanceVolume) :: elem
          integer i,j,k

          select type (PB)
            class is (VolumeSourceBody)
              if (associated(PB%get_resistance)) then
                do k = 0,Unz+1
                 do j = 0,Uny+1
                  do i = 0,Unx+1
                     if (Inside(PB,xU(i),yPr(j),zPr(k))) then
                       elem%xi = i
                       elem%yj = j
                       elem%zk = k
                       elem%flux = PB%get_resistance(xU(i),yPr(j),zPr(k))
                       call UResistanceVolumesList%add(elem)
                     end if
                  enddo
                 enddo
                enddo
              end if
            class default
              stop "Not a VolumeSourceBody."
          end select
        end subroutine

        subroutine GetVCells(PB)
          class(*) :: PB
          type(TResistanceVolume) :: elem
          integer i,j,k

          select type (PB)
            class is (VolumeSourceBody)
              if (associated(PB%get_resistance)) then
                do k = 0,Vnz+1
                 do j = 0,Vny+1
                  do i = 0,Vnx+1
                     if (Inside(PB,xPr(i),yV(j),zPr(k))) then
                       elem%xi = i
                       elem%yj = j
                       elem%zk = k
                       elem%flux = PB%get_resistance(xPr(i),yV(j),zPr(k))
                       call VResistanceVolumesList%add(elem)
                     end if
                  enddo
                 enddo
                enddo
              end if
            class default
              stop "Not a VolumeSourceBody."
          end select
        end subroutine

        subroutine GetWCells(PB)
          class(*) :: PB
          type(TResistanceVolume) :: elem
          integer i,j,k

          select type (PB)
            class is (VolumeSourceBody)
              if (associated(PB%get_resistance)) then
                do k = 0,Wnz+1
                 do j = 0,Wny+1
                  do i = 0,Wnx+1
                     if (Inside(PB,xPr(i),yPr(j),zW(k))) then
                       elem%xi = i
                       elem%yj = j
                       elem%zk = k
                       elem%flux = PB%get_resistance(xPr(i),yPr(j),zW(k))
                       call WResistanceVolumesList%add(elem)
                     end if
                  enddo
                 enddo
                enddo
              end if
            class default
              stop "Not a VolumeSourceBody."
          end select
        end subroutine

        subroutine GetOtherCells(PB)
          use SolarRadiation, only: enable_radiation
          class(*) :: PB
          type(TemperatureFlVolume) :: Telem
          type(MoistureFlVolume) :: Melem
          type(ScalarFlVolume) :: Selem
          integer i,j,k,sc
          
          select type (PB)
            class is (VolumeSourceBody)
              InsidePr = .false.
              
              do k = 0,Prnz+1
               do j = 0,Prny+1
                do i = 0,Prnx+1
                   InsidePr(i,j,k) = PB%Inside(xPr(i),yPr(j),zPr(k))
                end do
               end do
              end do
              
            
              do k = 0,Prnz+1
               do j = 0,Prny+1
                do i = 0,Prnx+1
                   if (InsidePr(i,j,k)) then
                   
                     Telem%xi = i
                     Telem%yj = j
                     Telem%zk = k
                     Melem%xi = i
                     Melem%yj = j
                     Melem%zk = k
                       
                     if (enable_buoyancy==1 .and. enable_radiation==1) then

                       associate (p=>PB) !workaround for gfortran 4.8 bug
                       select type (p)
                         class is (PlantBody)
                           if (on_border(i,j,k)) then
                             call GetRadiationFluxes(p,Telem,Melem,xPr(i),yPr(j),zPr(k))
                           else
                             Telem%flux = 0
                             Melem%flux = 0
                           end if
                         class default
                           Telem%flux = 0
                           Melem%flux = 0
                       end select
                       end associate
                     else
                       Telem%flux = 0
                       Melem%flux = 0
                     end if
                     
                     if (enable_buoyancy==1 .and. associated(PB%get_temperature_flux)) then
                     
                         Telem%flux = Telem%flux + &
                                      PB%get_temperature_flux(xPr(i),yPr(j),zPr(k))
                     end if
                     
                     if (enable_moisture==1 .and. associated(PB%get_moisture_flux)) then

                         Melem%flux = Melem%flux + &
                                      PB%get_moisture_flux(xPr(i),yPr(j),zPr(k))
                       
                     end if

                     if (Telem%flux/=0) call TemperatureFlVolumesList%add(Telem)

                     if (Melem%flux/=0) call MoistureFlVolumesList%add(Melem)
                      
                     
                     if (associated(PB%get_scalar_flux)) then
                       do sc = 1,num_of_scalars
                         Selem%xi = i
                         Selem%yj = j
                         Selem%zk = k
                         Selem%flux = PB%get_scalar_flux(xPr(i),yPr(j),zPr(k),sc)
                         if (Selem%flux/=0) call ScalarFlVolumesLists(sc)%list%add(Selem)
                       end do
                     end if
                   end if
                enddo
               enddo
              enddo
            class default
              stop "Not a VolumeSourceBody."
          end select
        end subroutine

        
        
        logical function on_border(xi,yj,zk)
          integer xi,yj,zk
          on_border = InsidePr(xi,yj,zk) .and. &
                      .not. all([InsidePr(xi-1,yj,zk), &
                                 InsidePr(xi+1,yj,zk), &
                                 InsidePr(xi,yj-1,zk), &
                                 InsidePr(xi,yj+1,zk), &
                                 InsidePr(xi,yj,zk-1), &
                                 InsidePr(xi,yj,zk+1)])
        end function
    end subroutine InsideCellsToLists


    
    
    
    subroutine GetRadiationFluxes(PB,Telem,Melem,x,y,z)
      use SolarRadiation
      use PhysicalProperties
      use GeometricShapes, only: Ray
      use SolidBodies, only: SolidBodiesList, SolidBody
      class(PlantBody),intent(inout) :: PB
      type(TemperatureFlVolume),intent(inout) :: Telem
      type(MoistureFlVolume),intent(inout)    :: Melem
      real(knd),intent(in) :: x, y, z
      
      type(Ray) :: r
      real(knd) :: out_norm(3)
      real(knd) :: inc_radiation_flux, radiation_balance, angle_to_sun, &
                   total_heat_flux, sensible_heat_flux, latent_heat_flux,&
                   xnear, ynear, znear, distx, disty, distz
      real(knd) :: xr, yr, zr, svf
                   
      call PB%ClosestOut(xnear,ynear,znear,x,y,z)

      distx = xnear - x
      disty = ynear - y
      distz = znear - z
      
      out_norm = [distx, disty, distz] / &
                 hypot(hypot(distx,disty),distz)
  
      angle_to_sun = max(0._knd, dot_product(out_norm, vector_to_sun))

      xr = x + distx * 1.1_knd
      yr = y + disty * 1.1_knd
      zr = z + distz * 1.1_knd
      if (angle_to_sun>0._knd) then
        r = Ray([xr, yr, zr], &
                vector_to_sun)

        if ((SolidBodiesList%any(DoIntersect)) .or. &
            (SourceBodiesList%any(DoIntersectPlants))) then
          angle_to_sun = 0
        end if
      end if
      
      !temporary
      svf = sky_view_factor(xr,yr,zr)
      
      inc_radiation_flux = angle_to_sun * solar_direct_flux() + &
                           solar_diffuse_flux()*svf + &
                           in_lw_radiation()

      radiation_balance = inc_radiation_flux * (1-PB%albedo) - &
                          out_lw_radiation(PB%emmissivity, temperature_ref) * svf

      total_heat_flux = radiation_balance - &
                        (0.1 * radiation_balance) !crude guess of the storage flux

      !FIXME: surface flux to volume flux APPROXIMATE!                  
      total_heat_flux = total_heat_flux/dot_product([dxmin,dymin,dzmin],abs(out_norm))
                        
      if (enable_moisture==1) then
        latent_heat_flux = total_heat_flux * PB%evaporative_fraction
        sensible_heat_flux = total_heat_flux - latent_heat_flux
        Melem%flux  = latent_heat_flux / (rho_air_ref * Lv_water_ref)
      else
        sensible_heat_flux = total_heat_flux
        Melem%flux  = 0
      end if
      
      Telem%flux = sensible_heat_flux / (rho_air_ref * Cp_air_ref)

      contains
          
        real(knd) function sky_view_factor(x,y,z)
          real(knd),intent(in) :: x,y,z
          integer :: i,nfree

          
          nfree = 0
          
          do i=1,svf_nrays
            r = Ray([x,y,z],svf_vecs(:,i))
          
            if (.not.((SolidBodiesList%any(DoIntersect)) .or. &
                      (SourceBodiesList%any(DoIntersectPlants)))) then
              nfree = nfree + 1
            end if
          end do
          
          sky_view_factor = real(nfree,knd) / svf_nrays
          
        end function
  
        logical function DoIntersect(item) result(res)
          class(*) :: item
          
          select type (item)
            class is (SolidBody)
              res = item%IntersectsRay(r)
            class default
              stop "Non-SolidBody i the list."
          end select
        end function
        
        logical function DoIntersectPlants(item) result(res)
          class(*) :: item
          
          select type (item)
            class is (PlantBody)
              res = item%IntersectsRay(r)
            class is (VolumeSourceBody)
              res = .false.
            class default
              stop "Non-SolidBody i the list."
          end select
        end function
      

    end subroutine

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    subroutine MovePointsToArray
    !It would be posible call a generic procedure with the list and array
    ! as an argument, but we want the final arrays not polymorphic.
      integer :: i, j, component

      allocate(UResistanceVolumes(UResistanceVolumesList%Len()))
      allocate(VResistanceVolumes(VResistanceVolumesList%Len()))
      allocate(WResistanceVolumes(WResistanceVolumesList%Len()))
      allocate(TemperatureFlVolumes(TemperatureFlVolumesList%Len()))
      allocate(MoistureFlVolumes(MoistureFlVolumesList%Len()))

      allocate(ScalarFlVolumes(num_of_scalars))
      
      do j=1,num_of_scalars
        allocate(ScalarFlVolumes(j)%volumes(ScalarFlVolumesLists(j)%list%Len()))
      end do

      i = 0
      component = 1
      call UResistanceVolumesList%for_each(CopyPoint)
      i = 0
      component = 2
      call VResistanceVolumesList%for_each(CopyPoint)
      i = 0
      component = 3
      call WResistanceVolumesList%for_each(CopyPoint)
      i = 0
      call TemperatureFlVolumesList%for_each(CopyPoint)
      i = 0
      call MoistureFlVolumesList%for_each(CopyPoint)

      do j=1,num_of_scalars
        i = 0
        call ScalarFlVolumesLists(j)%list%for_each(CopyPoint)
      end do
      
      if (enable_radiation==1) call SaveFluxes

      contains

        subroutine CopyPoint(elem)
          class(*) :: elem

          i = i + 1

          select type (elem)
            type is (TResistanceVolume)          
              if (component==1) then
                UResistanceVolumes(i) = elem
              elseif (component==2) then
                VResistanceVolumes(i) = elem
              else
                WResistanceVolumes(i) = elem
              endif
            type is (TemperatureFlVolume)
              TemperatureFlVolumes(i) = elem
            type is (MoistureFlVolume)
              MoistureFlVolumes(i) = elem
            type is (ScalarFlVolume)
              ScalarFlVolumes(j)%volumes(i) = elem
            class default
              stop "Type error in volume source list."
          end select
        end subroutine

        subroutine SaveFluxes
          use VTKArray
          real(knd),allocatable :: temperature_flux(:,:,:)
          allocate(temperature_flux(0:Prnx+1,0:Prny+1,0:Prnz+1))
          
          temperature_flux = 0
          
          do i=1,size(TemperatureFlVolumes)
            associate(p => TemperatureFlVolumes(i))
              temperature_flux(p%xi,p%yj,p%zk) = p%flux
            end associate
          end do
          
          call VtkArraySimple("output/tempflplants.vtk",temperature_flux)
          
        end subroutine
        
    end subroutine MovePointsToArray



    subroutine InitVolumeSources
 
      call InsideCellsToLists
 
      call MovePointsToArray

    end subroutine InitVolumeSources

    

    subroutine ResistanceForce(U2,V2,W2,U,V,W)
      real(knd),dimension(-2:,-2:,-2:),intent(inout) :: U2,V2,W2
      real(knd),dimension(-2:,-2:,-2:),intent(in)    :: U,V,W

      call apply(U2,U,UResistanceVolumes,totU_u)

      call apply(V2,V,VResistanceVolumes,totU_v)

      call apply(W2,W,WResistanceVolumes,totU_w)

      contains

         pure real(knd) function totU_u(i,j,k)
           integer,intent(in) :: i,j,k

           totU_u = hypot( hypot( U(i,j,k), &
                                  (V(i,j,k)+V(i,j-1,k)+V(i-1,j,k)+V(i-1,j-1,k))/4._knd ) , &
                                  (W(i,j,k)+W(i,j,k-1)+W(i-1,j,k)+W(i-1,j,k-1))/4._knd )
         end function

         pure real(knd) function totU_v(i,j,k)
           integer,intent(in) :: i,j,k

           totU_v = hypot( hypot( (U(i,j,k)+U(i-1,j,k)+U(i,j-1,k)+U(i-1,j-1,k))/4._knd, &
                                  V(i,j,k) ) , &
                                  (W(i,j,k)+W(i,j,k-1)+W(i,j-1,k)+W(i,j-1,k-1))/4._knd )
         end function

         pure real(knd) function totU_w(i,j,k)
           integer,intent(in) :: i,j,k

           totU_w = hypot( hypot( (U(i,j,k)+U(i-1,j,k)+U(i,j,k-1)+U(i-1,j,k-1))/4._knd, &
                                  (V(i,j,k)+V(i,j-1,k)+V(i,j,k-1)+V(i,j-1,k-1))/4._knd ) , &
                                  W(i,j,k) )
         end function

         subroutine apply(X2,X,src,fun)
           real(knd),dimension(-2:,-2:,-2:),intent(inout)  :: X2
           real(knd),dimension(-2:,-2:,-2:),intent(in)     :: X
           type(TResistanceVolume),allocatable,intent(in)  :: src(:)
           procedure(totU_u) :: fun
           integer i

           if (allocated(src)) then
             do i=1,size(src)
               associate (xi => src(i)%xi, yj => src(i)%yj, zk => src(i)%zk)
                 X2(xi,yj,zk) = X2(xi,yj,zk) - dt * src(i)%flux * fun(xi,yj,zk) * X(xi,yj,zk)              
               end associate
             end do
           end if
         end subroutine
         
    end subroutine ResistanceForce

    subroutine TemperatureVolumeSources(Temperature)
      real(knd),dimension(-1:,-1:,-1:),intent(inout) :: Temperature
      integer i
!       call FluxKernel(Temperature,TemperatureFlVolumes)
      associate (X=>Temperature, src=> TemperatureFlVolumes)
       do i=1,size(src)
         associate (xi => src(i)%xi, yj => src(i)%yj, zk => src(i)%zk)
           X(xi,yj,zk) = X(xi,yj,zk) + dt * src(i)%flux
         end associate
       end do
      end associate      
    end subroutine TemperatureVolumeSources

    subroutine MoistureVolumeSources(Moisture)
      real(knd),dimension(-1:,-1:,-1:),intent(inout) :: Moisture
      integer i
!       call FluxKernel(Moisture,MoistureFlVolumes)
      associate (X=>Moisture, src=> MoistureFlVolumes)
       do i=1,size(src)
         associate (xi => src(i)%xi, yj => src(i)%yj, zk => src(i)%zk)
           X(xi,yj,zk) = X(xi,yj,zk) + dt * src(i)%flux
         end associate
       end do
      end associate
    end subroutine MoistureVolumeSources

    subroutine ScalarVolumeSources(Scalar, sc)
      real(knd),dimension(-1:,-1:,-1:),intent(inout) :: Scalar
      integer, intent(in) :: sc
      integer j

      if (allocated(ScalarFlVolumes)) then
        call FluxKernel(Scalar(:,:,:), ScalarFlVolumes(sc)%volumes)
      end if

    end subroutine ScalarVolumeSources

    subroutine FluxKernel(X,src)
      real(knd),intent(inout) :: X(-1:,-1:,-1:)
      class(FluxVolume),intent(in) :: src(:)
      integer i
      !Assume src is allocated. It must hold if we called MovePointsToArray properly.
       do i=1,size(src)
         associate (xi => src(i)%xi, yj => src(i)%yj, zk => src(i)%zk)
           X(xi,yj,zk) = X(xi,yj,zk) + dt * src(i)%flux
         end associate
       end do

    end subroutine FluxKernel
    
    
    
    
    subroutine AddScalarFlVolumesContainer(l,r)
      type(ScalarFlVolumesContainer),intent(inout) :: l(:)
      type(ScalarFlVolumesContainer),intent(in)    :: r
      type(ScalarFlVolume),allocatable :: tmp(:)
      integer i
      associate (sn => r%scalar_number)
        if (size(r%volumes)>0) then
          !NOTE: the shorter version problematic in gfortran4.8 and ifort 14
          !l(sn)%volumes = [l(sn)%volumes, r%volumes]
          allocate(tmp( size(l(sn)%volumes) + size(r%volumes) ))
          
          tmp(1:size(l(sn)%volumes)) = l(sn)%volumes
          tmp(size(l(sn)%volumes)+1:) = r%volumes
          
          call move_alloc(tmp,l(sn)%volumes)
        end if
      end associate
    end subroutine

end module VolumeSources




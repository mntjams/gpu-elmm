module VolumeSources
  use Parameters
  use Lists
  use TBody_class
  use ImmersedBoundary

  implicit none
! 
!   private
! 
!   public UResistanceVolumes, VResistanceVolumes, WResistanceVolumes, &
!          TemperatureFlVolumes, MoistureFlVolumes, ScalarFlVolumes, TVolumeSourceBody

  type,extends(TListable) :: TGridVolume
    integer   :: xi, yj, zk
  end type TGridVolume

  type,extends(TGridVolume) :: TResistanceVolume
    ! f_drag_i = - Cd * a * V * u_i , V = sqrt(sum(u_i**2))
    real(KND) :: resistance  !Cd * a
  end type TResistanceVolume

  type,extends(TGridVolume) :: TTemperatureFlVolume
    ! <T'w'> = temperature_flux, Qh = rho*Cp*temperature_flux
    real(KND) :: temperature_flux
  end type TTemperatureFlVolume

  type,extends(TGridVolume) :: TMoistureFlVolume
    ! <q'w'> = moisture_flux, Qe = rho*Lv*moisture_flux
    real(KND) :: moisture_flux
  end type TMoistureFlVolume

  type,extends(TGridVolume) :: TScalarFlVolume
    ! <c'w'> = scalar_flux
    real(KND),allocatable :: scalar_flux(:)
  end type TScalarFlVolume

  type, extends(TBody) :: TVolumeSourceBody
    procedure(resistance_interface),pointer  :: get_resistance       => null()
    procedure(resistance_interface),pointer  :: get_temperature_flux => null()
    procedure(resistance_interface),pointer  :: get_moisture_flux    => null()
!     procedure(scalar_flux_interface),pointer :: get_scalar_flux      => null()
  end type TVolumeSourceBody

  abstract interface
    function resistance_interface(self,x,y,z) result(res)
      import
      real(KND) :: res
      class(TVolumeSourceBody),intent(in) :: self
      real(KND),intent(in) :: x,y,z
    end function
!     function scalar_flux_interface(self,x,y,z) result(res)
!       import
!       real(KND) :: res(computescalars) !HACK, allocatable causes ICE in gfortran 4.7.1
!       class(TVolumeSourceBody),intent(in) :: self
!       real(KND),intent(in) :: x,y,z
!     end function
  end interface

  type(TList) :: SourceBodiesList

  type(TList) :: UResistanceVolumesList, &
                 VResistanceVolumesList, &
                 WResistanceVolumesList, &
                 TemperatureFlVolumesList, &
                 MoistureFlVolumesList, &
                 ScalarFlVolumesList

  !final volume momentum sinks - resistances
  type(TResistanceVolume),allocatable :: UResistanceVolumes(:), &
                                         VResistanceVolumes(:), &
                                         WResistanceVolumes(:)

  !final scalar quantities sources/sinks
  type(TTemperatureFlVolume),allocatable :: TemperatureFlVolumes(:)
  type(TMoistureFlVolume)   ,allocatable :: MoistureFlVolumes(:)
  type(TScalarFlVolume)       ,allocatable :: ScalarFlVolumes(:)

  contains

    subroutine InsideCellsToLists
      !find cells inside the canopy and store them in a list

      call SourceBodiesList%ForEach(GetUCells)

      call SourceBodiesList%ForEach(GetVCells)

      call SourceBodiesList%ForEach(GetWCells)

      call SourceBodiesList%ForEach(GetOtherCells)
      contains

        subroutine GetUCells(PB)
          class(TListable) :: PB
          type(TResistanceVolume) :: elem
          integer i,j,k

          select type (PB)
            class is (TVolumeSourceBody)
              if (associated(PB%get_resistance)) then
                !$omp parallel do private(i,j,k,elem)
                do k = 0,Unz+1
                 do j = 0,Uny+1
                  do i = 0,Unx+1
                     if (PB%Inside(xU(i),yPr(j),zPr(k))) then
                       elem%resistance = PB%get_resistance(xU(i),yPr(j),zPr(k))
                       call UResistanceVolumesList%Add(elem)
                     end if
                  enddo
                 enddo
                enddo
                !$omp end parallel do
              end if
            class default
              stop "Not a TVolumeSourceBody."
          end select
        end subroutine

        subroutine GetVCells(PB)
          class(TListable) :: PB
          type(TResistanceVolume) :: elem
          integer i,j,k

          select type (PB)
            class is (TVolumeSourceBody)
              if (associated(PB%get_resistance)) then
                !$omp parallel do private(i,j,k,elem)
                do k = 0,Vnz+1
                 do j = 0,Vny+1
                  do i = 0,Vnx+1
                     if (PB%Inside(xPr(i),yV(j),zPr(k))) then
                       elem%resistance = PB%get_resistance(xPr(i),yV(j),zPr(k))
                       call VResistanceVolumesList%Add(elem)
                     end if
                  enddo
                 enddo
                enddo
                !$omp end parallel do
              end if
            class default
              stop "Not a TVolumeSourceBody."
          end select
        end subroutine

        subroutine GetWCells(PB)
          class(TListable) :: PB
          type(TResistanceVolume) :: elem
          integer i,j,k

          select type (PB)
            class is (TVolumeSourceBody)
              if (associated(PB%get_resistance)) then
                !$omp parallel do private(i,j,k,elem)
                do k = 0,Wnz+1
                 do j = 0,Wny+1
                  do i = 0,Wnx+1
                     if (PB%Inside(xPr(i),yPr(j),zW(k))) then
                       elem%resistance = PB%get_resistance(xPr(i),yPr(j),zW(k))
                       call WResistanceVolumesList%Add(elem)
                     end if
                  enddo
                 enddo
                enddo
                !$omp end parallel do
              end if
            class default
              stop "Not a TVolumeSourceBody."
          end select
        end subroutine

        subroutine GetOtherCells(PB)
          class(TListable) :: PB
          type(TTemperatureFlVolume) :: Telem
          type(TMoistureFlVolume) :: Melem
          type(TScalarFlVolume) :: Selem
          integer i,j,k

          select type (PB)
            class is (TVolumeSourceBody)
              allocate(Selem%scalar_flux(computescalars))  !to avoid allocation inside the loop
              !$omp parallel do private(i,j,k,Telem,Melem,Selem)
              do k = 0,Prnz+1
               do j = 0,Prny+1
                do i = 0,Prnx+1
                   if (PB%Inside(xPr(i),yPr(j),zPr(k))) then
                     if (associated(PB%get_temperature_flux)) then
                       Telem%temperature_flux = PB%get_temperature_flux(xPr(i),yPr(j),zPr(k))
                       call TemperatureFlVolumesList%Add(Telem)
                     end if
                     if (associated(PB%get_moisture_flux)) then
                       Melem%moisture_flux = PB%get_moisture_flux(xPr(i),yPr(j),zPr(k))
                       call MoistureFlVolumesList%Add(Melem)
                     end if
!                      if (associated(PB%get_scalar_flux)) then
!                        Selem%scalar_flux(:) = PB%get_scalar_flux(xPr(i),yPr(j),zPr(k)) !(:) avoiding reallocation
!                        call ScalarFlVolumesList%Add(Selem)
!                      end if
                   end if
                enddo
               enddo
              enddo
              !$omp end parallel do
            class default
              stop "Not a TVolumeSourceBody."
          end select
        end subroutine

    end subroutine InsideCellsToLists


    subroutine MovePointsToArray
    !It would be posible call a generic procedure with the list and array
    ! as an argument, but we want the final arrays not polymorphic.
      integer :: i, component

      component = 1
      allocate(UResistanceVolumes(UResistanceVolumesList%Len()))
      component = 2
      allocate(VResistanceVolumes(VResistanceVolumesList%Len()))
      component = 3
      allocate(WResistanceVolumes(WResistanceVolumesList%Len()))
      allocate(TemperatureFlVolumes(TemperatureFlVolumesList%Len()))
      allocate(MoistureFlVolumes(MoistureFlVolumesList%Len()))
      allocate(ScalarFlVolumes(ScalarFlVolumesList%Len()))

      i = 0
      call UResistanceVolumesList%ForEach(CopyIBPoint)
      i = 0
      call VResistanceVolumesList%ForEach(CopyIBPoint)
      i = 0
      call WResistanceVolumesList%ForEach(CopyIBPoint)
      i = 0
      call TemperatureFlVolumesList%ForEach(CopyIBPoint)
      i = 0
      call MoistureFlVolumesList%ForEach(CopyIBPoint)
      i = 0
      call ScalarFlVolumesList%ForEach(CopyIBPoint)

      contains

        subroutine CopyIBPoint(elem)
          class(TListable) :: elem

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
            type is (TTemperatureFlVolume)
              TemperatureFlVolumes(i) = elem
            type is (TMoistureFlVolume)
              MoistureFlVolumes(i) = elem
            type is (TScalarFlVolume)
              ScalarFlVolumes(i) = elem
            class default
              stop "Type error in volume source list."
          end select
        end subroutine

    end subroutine MovePointsToArray


end module VolumeSources




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
    procedure(scalar_flux_interface),pointer :: get_scalar_flux      => null()
  end type TVolumeSourceBody

  abstract interface
    function resistance_interface(self,x,y,z) result(res)
      import
      real(KND) :: res
      class(TVolumeSourceBody),intent(in) :: self
      real(KND),intent(in) :: x,y,z
    end function
    function scalar_flux_interface(self,x,y,z) result(res)
      import
      real(KND),allocatable :: res(:) !allocatable causes ICE in gfortran 4.7.1, ok in repository since March 2013
      class(TVolumeSourceBody),intent(in) :: self
      real(KND),intent(in) :: x,y,z
    end function
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
                do k = 0,Unz+1
                 do j = 0,Uny+1
                  do i = 0,Unx+1
                     if (Inside(PB,xU(i),yPr(j),zPr(k))) then
                       elem%xi = i
                       elem%yj = j
                       elem%zk = k
                       elem%resistance = PB%get_resistance(xU(i),yPr(j),zPr(k))
                       call UResistanceVolumesList%Add(elem)
                     end if
                  enddo
                 enddo
                enddo
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
                do k = 0,Vnz+1
                 do j = 0,Vny+1
                  do i = 0,Vnx+1
                     if (Inside(PB,xPr(i),yV(j),zPr(k))) then
                       elem%xi = i
                       elem%yj = j
                       elem%zk = k
                       elem%resistance = PB%get_resistance(xPr(i),yV(j),zPr(k))
                       call VResistanceVolumesList%Add(elem)
                     end if
                  enddo
                 enddo
                enddo
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
                do k = 0,Wnz+1
                 do j = 0,Wny+1
                  do i = 0,Wnx+1
                     if (Inside(PB,xPr(i),yPr(j),zW(k))) then
                       elem%xi = i
                       elem%yj = j
                       elem%zk = k
                       call WResistanceVolumesList%Add(elem)
                     end if
                  enddo
                 enddo
                enddo
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
              do k = 0,Prnz+1
               do j = 0,Prny+1
                do i = 0,Prnx+1
                   if (Inside(PB,xPr(i),yPr(j),zPr(k))) then
                     if (associated(PB%get_temperature_flux)) then
                       Telem%xi = i
                       Telem%yj = j
                       Telem%zk = k
                       Telem%temperature_flux = PB%get_temperature_flux(xPr(i),yPr(j),zPr(k))
                       call TemperatureFlVolumesList%Add(Telem)
                     end if
                     if (associated(PB%get_moisture_flux)) then
                       Melem%xi = i
                       Melem%yj = j
                       Melem%zk = k
                       Melem%moisture_flux = PB%get_moisture_flux(xPr(i),yPr(j),zPr(k))
                       call MoistureFlVolumesList%Add(Melem)
                     end if
                     if (associated(PB%get_scalar_flux)) then
                       Selem%xi = i
                       Selem%yj = j
                       Selem%zk = k
                       Selem%scalar_flux(:) = PB%get_scalar_flux(xPr(i),yPr(j),zPr(k)) !(:) avoiding reallocation
                       call ScalarFlVolumesList%Add(Selem)
                     end if
                   end if
                enddo
               enddo
              enddo
            class default
              stop "Not a TVolumeSourceBody."
          end select
        end subroutine

    end subroutine InsideCellsToLists


    subroutine MovePointsToArray
    !It would be posible call a generic procedure with the list and array
    ! as an argument, but we want the final arrays not polymorphic.
      integer :: i, component

      allocate(UResistanceVolumes(UResistanceVolumesList%Len()))
      allocate(VResistanceVolumes(VResistanceVolumesList%Len()))
      allocate(WResistanceVolumes(WResistanceVolumesList%Len()))
      allocate(TemperatureFlVolumes(TemperatureFlVolumesList%Len()))
      allocate(MoistureFlVolumes(MoistureFlVolumesList%Len()))
      allocate(ScalarFlVolumes(ScalarFlVolumesList%Len()))

      i = 0
      component = 1
      call UResistanceVolumesList%ForEach(CopyPoint)
      i = 0
      component = 2
      call VResistanceVolumesList%ForEach(CopyPoint)
      i = 0
      component = 3
      call WResistanceVolumesList%ForEach(CopyPoint)
      i = 0
      call TemperatureFlVolumesList%ForEach(CopyPoint)
      i = 0
      call MoistureFlVolumesList%ForEach(CopyPoint)
      i = 0
      call ScalarFlVolumesList%ForEach(CopyPoint)

      contains

        subroutine CopyPoint(elem)
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



    subroutine InitVolumeSources
    
      call InsideCellsToLists
      call MovePointsToArray
    end subroutine InitVolumeSources

    

    subroutine ResistanceForce(U2,V2,W2,U,V,W)
      real(KND),dimension(-2:,-2:,-2:),intent(inout) :: U2,V2,W2
      real(KND),dimension(-2:,-2:,-2:),intent(in)    :: U,V,W

      call apply(U2,U,UResistanceVolumes,totU_u)
      call apply(V2,V,VResistanceVolumes,totU_v)
      call apply(W2,W,WResistanceVolumes,totU_w)

      contains

         pure real(KND) function totU_u(i,j,k)
           integer,intent(in) :: i,j,k

           totU_u = hypot( hypot( U(i,j,k), &
                                  (V(i,j,k)+V(i,j-1,k)+V(i-1,j,k)+V(i-1,j-1,k))/4._KND ) , &
                                  (W(i,j,k)+W(i,j,k-1)+W(i-1,j,k)+W(i-1,j,k-1))/4._KND )
         end function

         pure real(KND) function totU_v(i,j,k)
           integer,intent(in) :: i,j,k

           totU_v = hypot( hypot( (U(i,j,k)+U(i-1,j,k)+U(i,j-1,k)+U(i-1,j-1,k))/4._KND, &
                                  V(i,j,k) ) , &
                                  (W(i,j,k)+W(i,j,k-1)+W(i,j-1,k)+W(i,j-1,k-1))/4._KND )
         end function

         pure real(KND) function totU_w(i,j,k)
           integer,intent(in) :: i,j,k

           totU_w = hypot( hypot( (U(i,j,k)+U(i-1,j,k)+U(i,j,k-1)+U(i-1,j,k-1))/4._KND, &
                                  (V(i,j,k)+V(i,j-1,k)+V(i,j,k-1)+V(i,j-1,k-1))/4._KND ) , &
                                  W(i,j,k) )
         end function

         subroutine apply(X2,X,src,fun)
           real(KND),dimension(-2:,-2:,-2:),intent(inout)  :: X2
           real(KND),dimension(-2:,-2:,-2:),intent(in)     :: X
           type(TResistanceVolume),allocatable,intent(in)  :: src(:)
           procedure(totU_u) :: fun
           integer i

           if (allocated(src)) then
             do i=1,size(src)
               associate (xi => src(i)%xi, yj => src(i)%yj, zk => src(i)%zk)
                 X2(xi,yj,zk) = X2(xi,yj,zk) - dt * src(i)%resistance * fun(xi,yj,zk) * X(xi,yj,zk)
               end associate
             end do
           end if
         end subroutine
         
    end subroutine ResistanceForce



end module VolumeSources




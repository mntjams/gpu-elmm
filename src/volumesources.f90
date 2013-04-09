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

  type,extends(TGridVolume) :: TFluxVolume
    real(knd) :: flux
  end type TFluxVolume

  type,extends(TFluxVolume) :: TResistanceVolume
    ! f_drag_i = - Cd * a * V * u_i , V = sqrt(sum(u_i**2))
    !flux = Cd * a
  end type TResistanceVolume

  type,extends(TFluxVolume) :: TTemperatureFlVolume
    ! <T'w'> = temperature_flux, Qh = rho*Cp*temperature_flux
    !flux = temperature_flux
  end type TTemperatureFlVolume

  type,extends(TFluxVolume) :: TMoistureFlVolume
    ! <q'w'> = moisture_flux, Qe = rho*Lv*moisture_flux
    !flux = moisture_flux
  end type TMoistureFlVolume

  type,extends(TFluxVolume) :: TScalarFlVolume
    ! <c'w'> = scalar_flux
    !flux = flux
  end type TScalarFlVolume

  type TScalarFlVolumesContainer
    type(TScalarFlVolume), allocatable :: Volumes(:)
  end type

  type TScalarFlVolumesListContainer
    type(TList) :: list
  end type

  type, extends(TBody) :: TVolumeSourceBody
    procedure(resistance_interface),pointer  :: get_resistance       => null()
    procedure(resistance_interface),pointer  :: get_temperature_flux => null()
    procedure(resistance_interface),pointer  :: get_moisture_flux    => null()
    procedure(scalar_flux_interface),pointer :: get_scalar_flux      => null()
  end type TVolumeSourceBody

  abstract interface
    function resistance_interface(self,x,y,z) result(res)
      import
      real(knd) :: res
      class(TVolumeSourceBody),intent(in) :: self
      real(knd),intent(in) :: x,y,z
    end function
    function scalar_flux_interface(self,x,y,z,num_of_scalar) result(res)
      import
      real(knd) :: res !allocatable causes ICE in gfortran 4.7.1, ok in repository since March 2013
      class(TVolumeSourceBody),intent(in) :: self
      real(knd),intent(in) :: x,y,z
      integer,intent(in) :: num_of_scalar
    end function
  end interface

  type(TList) :: SourceBodiesList

  type(TList) :: UResistanceVolumesList, &
                 VResistanceVolumesList, &
                 WResistanceVolumesList, &
                 TemperatureFlVolumesList, &
                 MoistureFlVolumesList
                 
  type(TScalarFlVolumesListContainer),allocatable :: ScalarFlVolumesLists(:)

  !final volume momentum sinks - resistances
  type(TResistanceVolume),allocatable :: UResistanceVolumes(:), &
                                         VResistanceVolumes(:), &
                                         WResistanceVolumes(:)

  !final scalar quantities sources/sinks
  type(TTemperatureFlVolume),allocatable :: TemperatureFlVolumes(:)
  type(TMoistureFlVolume)   ,allocatable :: MoistureFlVolumes(:)
  type(TScalarFlVolumesContainer) ,allocatable :: ScalarFlVolumes(:)
  

  contains

    subroutine InsideCellsToLists
      !find cells inside the canopy and store them in a list

      !$omp parallel sections
      !$omp section
      call SourceBodiesList%ForEach(GetUCells)

      !$omp section
      call SourceBodiesList%ForEach(GetVCells)

      !$omp section
      call SourceBodiesList%ForEach(GetWCells)

      !$omp section
      allocate(ScalarFlVolumesLists(num_of_scalars))

      call SourceBodiesList%ForEach(GetOtherCells)
      !$omp end parallel sections

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
                       elem%flux = PB%get_resistance(xU(i),yPr(j),zPr(k))
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
                       elem%flux = PB%get_resistance(xPr(i),yV(j),zPr(k))
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
                       elem%flux = PB%get_resistance(xPr(i),yPr(j),zW(k))
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
          integer i,j,k,sc

          select type (PB)
            class is (TVolumeSourceBody)
              do k = 0,Prnz+1
               do j = 0,Prny+1
                do i = 0,Prnx+1
                   if (Inside(PB,xPr(i),yPr(j),zPr(k))) then
                     if (associated(PB%get_temperature_flux)) then
                       Telem%xi = i
                       Telem%yj = j
                       Telem%zk = k
                       Telem%flux = PB%get_temperature_flux(xPr(i),yPr(j),zPr(k))
                       if (Telem%flux/=0) call TemperatureFlVolumesList%Add(Telem)
                     end if
                     if (associated(PB%get_moisture_flux)) then
                       Melem%xi = i
                       Melem%yj = j
                       Melem%zk = k
                       Melem%flux = PB%get_moisture_flux(xPr(i),yPr(j),zPr(k))
                       if (Melem%flux/=0) call MoistureFlVolumesList%Add(Melem)
                     end if
                     if (associated(PB%get_scalar_flux)) then
                       do sc = 1,num_of_scalars
                         Selem%xi = i
                         Selem%yj = j
                         Selem%zk = k
                         Selem%flux = PB%get_scalar_flux(xPr(i),yPr(j),zPr(k),sc) !(:) avoiding reallocation
                         if (Selem%flux/=0) call ScalarFlVolumesLists(sc)%list%Add(Selem)
                       end do
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
      integer :: i, j, component

      allocate(UResistanceVolumes(UResistanceVolumesList%Len()))
      allocate(VResistanceVolumes(VResistanceVolumesList%Len()))
      allocate(WResistanceVolumes(WResistanceVolumesList%Len()))
      allocate(TemperatureFlVolumes(TemperatureFlVolumesList%Len()))
      allocate(MoistureFlVolumes(MoistureFlVolumesList%Len()))

      allocate(ScalarFlVolumes(num_of_scalars))
      
      do j=1,num_of_scalars
        allocate(ScalarFlVolumes(j)%Volumes(ScalarFlVolumesLists(j)%list%Len()))
      end do

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

      do j=1,num_of_scalars
        i = 0
        call ScalarFlVolumesLists(j)%list%ForEach(CopyPoint)
      end do

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
              ScalarFlVolumes(j)%Volumes(i) = elem
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

    subroutine ScalarVolumeSources(Scalar)
      real(knd),dimension(-1:,-1:,-1:,-1:),intent(inout) :: Scalar
      integer j

      if (allocated(ScalarFlVolumes)) then
        do j=1,num_of_scalars
          call FluxKernel(Scalar(:,:,:,j),ScalarFlVolumes(j)%Volumes)
        end do
      end if

    end subroutine ScalarVolumeSources

    subroutine FluxKernel(X,src)
      real(knd),intent(inout) :: X(-1:,-1:,-1:)
      class(TFluxVolume),intent(in) :: src(:)
      integer i
      !Assume src is allocated. It must hold if we called MovePointsToArray properly.
       do i=1,size(src)
         associate (xi => src(i)%xi, yj => src(i)%yj, zk => src(i)%zk)
           X(xi,yj,zk) = X(xi,yj,zk) + dt * src(i)%flux
         end associate
       end do

    end subroutine FluxKernel

end module VolumeSources




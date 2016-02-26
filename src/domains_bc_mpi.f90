module domains_bc_par

  use mpi, only: MPI_PROC_NULL, MPI_COMM_NULL, MPI_STATUS_SIZE
  use Parameters
  use custom_mpi
  use custom_par

  implicit none

  private

  public par_init_domain_boundary_conditions, &
         par_exchange_domain_bounds, &
         par_update_domain_bounds, &
         par_update_domain_bounds_UVW, &
         par_update_domain_bounds_temperature, &
         par_update_domain_bounds_moisture, &
         par_domain_bound_relaxation

  type dom_bc_buffer_copy
    logical :: enabled = .false.

    logical :: rescale_compatibility = .false.
    logical :: rescale_dU_dt = .true.

    logical :: relaxation = .true.

    integer :: remote_domain
    !This is the index of the cell boundary corresponding to the domain boundary 
    ! or to the nested domain boundary.
    integer :: position
    !The direction of the buffer from position, 1 to 6 (We to To).
    integer :: direction = 0

    !grid coordinates of the buffer (from 1 to 2)
    integer :: Ui1, Ui2, Vi1, Vi2, Wi1, Wi2, Pri1, Pri2
    integer :: Uj1, Uj2, Vj1, Vj2, Wj1, Wj2, Prj1, Prj2
    integer :: Uk1, Uk2, Vk1, Vk2, Wk1, Wk2, Prk1, Prk2

    !grid coordinates of the boundary region (from 1 to 2)
    integer :: bUi1, bUi2, bVi1, bVi2, bWi1, bWi2, bPri1, bPri2
    integer :: bUj1, bUj2, bVj1, bVj2, bWj1, bWj2, bPrj1, bPrj2
    integer :: bUk1, bUk2, bVk1, bVk2, bWk1, bWk2, bPrk1, bPrk2

    !the most simple type, just a copy, no interpolation
    !the grid resolution must be identical
    !the boundary position must always exactly coincide with the grid cell boundaries
    real(knd), allocatable, dimension(:,:,:) :: U, V, W, Pr, Temperature, Moisture
    real(knd), allocatable, dimension(:,:,:,:) :: Scalar
    real(knd), allocatable, dimension(:,:,:) :: dU_dt, dV_dt, dW_dt, dPr_dt, &
                                                dTemperature_dt, dMoisture_dt
    real(knd), allocatable, dimension(:,:,:,:) :: dScalar_dt

    integer :: U_mpi_type, V_mpi_type, W_mpi_type, Pr_mpi_type

    !time at which the data is valid
    real(tim) :: time = -tiny(1.0_tim)

    !MPI communicator for the remote communication
    integer :: comm = MPI_COMM_NULL
    integer :: remote_rank = MPI_PROC_NULL
    !The buffers are transferred every `time_step_ratio` time steps
    integer :: time_step_ratio = 1 
  end type


  !send buffers should be indexed from 1, there may be more of them from several nested domains
  !they are not accessed by the boundary conditions routines
  type(dom_bc_buffer_copy), allocatable :: domain_bc_send_buffers(:)
  !receive buffers should be indexed with the index of the domain side (West to Top)
  ! to ease finding the right buffer to a given nested boundary 
  type(dom_bc_buffer_copy), allocatable :: domain_bc_recv_buffers_copy(:)

contains



  subroutine par_init_domain_boundary_conditions
    integer :: Ui1, Ui2, Vi1, Vi2, Wi1, Wi2, Pri1, Pri2
    integer :: Uj1, Uj2, Vj1, Vj2, Wj1, Wj2, Prj1, Prj2
    integer :: Uk1, Uk2, Vk1, Vk2, Wk1, Wk2, Prk1, Prk2
    !temporary testing hack

#include "custom_domains_setup.f90"

  end subroutine


  function U_mpi_subarray_type(i1, i2, j1, j2, k1, k2) result(res)
    integer :: res
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    integer :: ie

    call MPI_Type_create_subarray(3, &
                                  [Unx+6, Uny+6, Unz+6], &          !whole array shape
                                  [i2-i1+1, j2-j1+1, k2-k1+1], & !subarray shape
                                  [i1+2, j1+2, k1+2], &          !offsets, index + lbound
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  res, &
                                  ie)                                  
    if (ie/=0) call error_stop("Error creating MPI subarray derived type in U_mpi_subarray_type().")

    call MPI_Type_commit(res, ie)
    if (ie/=0) call error_stop("Error commiting MPI subarray derived type in U_mpi_subarray_type().")

  end function

  function V_mpi_subarray_type(i1, i2, j1, j2, k1, k2) result(res)
    integer :: res
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    integer :: ie

    call MPI_Type_create_subarray(3, &
                                  [Vnx+6, Vny+6, Vnz+6], &          !whole array shape
                                  [i2-i1+1, j2-j1+1, k2-k1+1], & !subarray shape
                                  [i1+2, j1+2, k1+2], &          !offsets, index + lbound
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  res, &
                                  ie)                                  
    if (ie/=0) call error_stop("Error creating MPI subarray derived type in V_mpi_subarray_type().")

    call MPI_Type_commit(res, ie)
    if (ie/=0) call error_stop("Error commiting MPI subarray derived type in V_mpi_subarray_type().")

  end function

  function W_mpi_subarray_type(i1, i2, j1, j2, k1, k2) result(res)
    integer :: res
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    integer :: ie

    call MPI_Type_create_subarray(3, &
                                  [Wnx+6, Wny+6, Wnz+6], &          !whole array shape
                                  [i2-i1+1, j2-j1+1, k2-k1+1], & !subarray shape
                                  [i1+2, j1+2, k1+2], &          !offsets, index + lbound
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  res, &
                                  ie)                                  
    if (ie/=0) call error_stop("Error creating MPI subarray derived type in W_mpi_subarray_type().")

    call MPI_Type_commit(res, ie)
    if (ie/=0) call error_stop("Error commiting MPI subarray derived type in W_mpi_subarray_type().")

  end function

  function Pr_mpi_subarray_type(i1, i2, j1, j2, k1, k2) result(res)
    integer :: res
    integer, intent(in) :: i1, i2, j1, j2, k1, k2
    integer :: ie

    call MPI_Type_create_subarray(3, &
                                  [Prnx+4, Prny+4, Prnz+4], &          !whole array shape
                                  [i2-i1+1, j2-j1+1, k2-k1+1], & !subarray shape
                                  [i1+1, j1+1, k1+1], &          !offsets, index + lbound
                                  MPI_ORDER_FORTRAN, &
                                  MPI_KND, &
                                  res, &
                                  ie)                                  
    if (ie/=0) call error_stop("Error creating MPI subarray derived type in Pr_mpi_subarray_type().")

    call MPI_Type_commit(res, ie)
    if (ie/=0) call error_stop("Error commiting MPI subarray derived type in Pr_mpi_subarray_type().")

  end function




  subroutine par_exchange_domain_bounds(U, V, W, Temperature, Moisture, Scalar, time, dt)
    !if necessary send or receive the new boundary conditions
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(in) :: U, V ,W 
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(in) :: Temperature, Moisture
    real(knd), dimension(-1:,-1:,-1:,1:), contiguous, intent(in) :: Scalar
    real(TIM), intent(in) :: time, dt

    integer, allocatable :: requests(:), statuses(:,:)
    integer :: i
    integer :: ie

    allocate(requests(0))

    if (enable_multiple_domains) then

      if (allocated(domain_bc_send_buffers)) then
        do i = lbound(domain_bc_send_buffers,1), &
               ubound(domain_bc_send_buffers,1)
          associate(b => domain_bc_send_buffers(i))

            if (b%enabled) then
              !derivatives approximated by the backward difference

              if (b%time<=-tiny(-1.0_tim)) then
                b%dU_dt = 0
                b%dV_dt = 0
                b%dW_dt = 0
                if (enable_buoyancy) then
                  b%dTemperature_dt = 0
                end if
                if (enable_moisture) then
                  b%dMoisture_dt = 0
                end if
              else
                b%dU_dt = (U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2) - b%U) / dt
                b%dV_dt = (V(b%Vi1:b%Vi2,b%Vj1:b%Vj2,b%Vk1:b%Vk2) - b%V) / dt
                b%dW_dt = (W(b%Wi1:b%Wi2,b%Wj1:b%Wj2,b%Wk1:b%Wk2) - b%W) / dt
                if (enable_buoyancy) then
                  b%dTemperature_dt = (Temperature(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2) - b%Temperature) / dt
                end if
                if (enable_moisture) then
                  b%dMoisture_dt = (Moisture(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2) - b%Moisture) / dt
                end if
              end if

              b%U = U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2)
              b%V = V(b%Vi1:b%Vi2,b%Vj1:b%Vj2,b%Vk1:b%Vk2)
              b%W = W(b%Wi1:b%Wi2,b%Wj1:b%Wj2,b%Wk1:b%Wk2)
              if (enable_buoyancy) then
                b%Temperature = Temperature(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2)
              end if
              if (enable_moisture) then
                b%Moisture = Moisture(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2)
              end if

              b%time = time

              call send_arrays(b)
            end if

          end associate
        end do
      end if


      if (allocated(domain_bc_recv_buffers_copy)) then
        do i = lbound(domain_bc_recv_buffers_copy,1), &
               ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(i))
  
            if (b%enabled) then
              call recv_arrays(b)
              b%time = time
            end if

          end associate
        end do
      end if


    end if

    call MPI_Waitall(size(requests), requests, MPI_STATUSES_IGNORE, ie)

  contains

    subroutine send_arrays(b)
      type(dom_bc_buffer_copy), intent(inout) :: b
      
      call send(b%U, b%U_mpi_type, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 1)
      call send(b%V, b%V_mpi_type, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 2)
      call send(b%W, b%W_mpi_type, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 3)

      call send(b%dU_dt, b%U_mpi_type, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 1)
      call send(b%dV_dt, b%V_mpi_type, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 2)
      call send(b%dW_dt, b%W_mpi_type, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 3)

      if (enable_buoyancy) then
        call send(b%Temperature, b%Pr_mpi_type, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 4)
        call send(b%dTemperature_dt, b%Pr_mpi_type, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 4)
      end if

      if (enable_moisture) then
        call send(b%Moisture, b%Pr_mpi_type, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 5)
        call send(b%dMoisture_dt, b%Pr_mpi_type, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 5)
      end if

    end subroutine

    subroutine recv_arrays(b)
      type(dom_bc_buffer_copy), intent(inout) :: b
      real(knd) :: avg
    
      call recv(b%U, b%U_mpi_type, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 1)   
      call recv(b%V, b%V_mpi_type, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 2)     
      call recv(b%W, b%W_mpi_type, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 3)    

      call recv(b%dU_dt, b%U_mpi_type, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 1)   
      call recv(b%dV_dt, b%V_mpi_type, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 2)     
      call recv(b%dW_dt, b%W_mpi_type, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 3)

      if (b%rescale_dU_dt) then
        !rescales the time derivative to be zero on average
        !makes sense only in certain situations
        select case (b%direction)
          case(We)
            avg = sum(b%dU_dt(0, 1:Uny, 1:Unz)) / (Uny*Unz)
            b%dU_dt = b%dU_dt - avg
          case(Ea)
            avg = sum(b%dU_dt(Unx+1, 1:Uny, 1:Unz)) / (Uny*Unz)
            b%dU_dt = b%dU_dt - avg
          case(So)
            avg = sum(b%dV_dt(1:Vnx, 0, 1:Vnz)) / (Vnx*Vnz)
            b%dV_dt = b%dV_dt - avg
          case(No)
            avg = sum(b%dV_dt(1:Vnx, Vny+1, 1:Vnz)) / (Vnx*Vnz)
            b%dV_dt = b%dV_dt - avg
          case(Bo)
            avg = sum(b%dW_dt(1:Wnx, 1:Wny, 0)) / (Wnx*Wny)
            b%dW_dt = b%dW_dt - avg
          case(To)
            avg = sum(b%dW_dt(1:Wnx, 1:Wny, Wnz+1)) / (Wnx*Wny)
            b%dW_dt = b%dW_dt - avg
        end select
      end if

      if (enable_buoyancy) then
        call recv(b%Temperature, b%Pr_mpi_type, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 4)
        call recv(b%dTemperature_dt, b%Pr_mpi_type, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 4)
      end if

      if (enable_moisture) then
        call recv(b%Moisture, b%Pr_mpi_type, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 5)
        call recv(b%dMoisture_dt, b%Pr_mpi_type, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 5)
      end if
    end subroutine

    subroutine send(a, dtype, rank, comm, tag)
      real(knd), intent(in), contiguous :: a(:,:,:)
      integer, intent(in) :: dtype, rank, comm, tag
      integer :: request, ie

      call MPI_ISsend(a, size(a), MPI_KND, rank, tag, comm, request, ie)
      if (ie/=0) stop "error sending MPI message."

      requests = [requests, request]
    end subroutine

    subroutine recv(a, dtype, rank, comm, tag)
      real(knd), intent(inout), contiguous :: a(:,:,:)
      integer, intent(in) :: dtype, rank, comm, tag
      integer :: request, ie

      call MPI_IRecv(a, size(a), MPI_KND, rank, tag, comm, request, ie)
      if (ie/=0) stop "error receiving MPI message."

      requests = [requests, request]
    end subroutine

  end subroutine par_exchange_domain_bounds



  subroutine par_update_domain_bounds_UVW(U, V, W, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V, W
    real(knd), intent(in) :: eff_time
    integer :: bi
    real(knd) :: S, t_diff

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do bi = lbound(domain_bc_recv_buffers_copy,1), &
                ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(bi))


             U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = b%U (b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2)                 

             V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) = b%V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2)

             W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) = b%W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2)


            if (eff_time > b%time) then
              t_diff = eff_time - b%time

              U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) + &
                                                       b%dU_dt(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) * t_diff
                

              V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) = V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) + &
                                                       b%dV_dt(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) * t_diff

              W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) = W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) + &
                                                       b%dW_dt(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) * t_diff
            end if
          end associate
        end do
      end if

    end if

  end subroutine



  subroutine par_update_domain_bounds_temperature(Temperature, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(inout) :: Temperature
    real(knd), intent(in) :: eff_time
    real(knd) :: t_diff
    integer :: i

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do i = lbound(domain_bc_recv_buffers_copy,1), &
               ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(i))

              Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                b%Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)

              t_diff = eff_time - b%time
              if (t_diff > 0) then
                  Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                  Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) + &
                    b%dTemperature_dt(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) * t_diff
              end if

          end associate
        end do
      end if

    end if

  end subroutine

  subroutine par_update_domain_bounds_moisture(Moisture, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(inout) :: Moisture
    real(knd), intent(in) :: eff_time
    real(knd) :: t_diff
    integer :: i

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do i = lbound(domain_bc_recv_buffers_copy,1), &
               ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(i))

              Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                b%Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)

              t_diff = eff_time - b%time
              if (t_diff > 0) then
                  Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                  Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) + &
                    b%dMoisture_dt(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) * t_diff
              end if

          end associate
        end do
      end if

    end if

  end subroutine

  subroutine par_update_domain_bounds(U, V, W, Temperature, Moisture, Scalar, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V ,W 
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(inout) :: Temperature, Moisture
    real(knd), dimension(-1:,-1:,-1:,1:), contiguous, intent(inout) :: Scalar
    real(knd), intent(in) :: eff_time
    real(knd) :: t_diff
    integer :: bi

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do bi = lbound(domain_bc_recv_buffers_copy,1), &
                ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(bi))

              U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = b%U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2)

              V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) = b%V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2)

              W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) = b%W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2)

              if (enable_buoyancy) then
                Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                  b%Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)
              end if
              if (enable_moisture) then
                Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                  b%Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2)
              end if

              if (eff_time > b%time) then
                t_diff = eff_time - b%time

                U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) + &
                                                         b%dU_dt(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) * t_diff
                V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) = V(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) + &
                                                         b%dV_dt(b%bVi1:b%bVi2,b%bVj1:b%bVj2,b%bVk1:b%bVk2) * t_diff
                W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) = W(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) + &
                                                         b%dW_dt(b%bWi1:b%bWi2,b%bWj1:b%bWj2,b%bWk1:b%bWk2) * t_diff

                if (enable_buoyancy) then
                  Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                  Temperature(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) + &
                    b%dTemperature_dt(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) * t_diff
                end if
                if (enable_moisture) then
                  Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) = &
                  Moisture(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) + &
                    b%dMoisture_dt(b%bPri1:b%bPri2,b%bPrj1:b%bPrj2,b%bPrk1:b%bPrk2) * t_diff
                end if
              end if

          end associate
        end do
      end if

    end if

  end subroutine


  subroutine par_domain_bound_relaxation(U, V, W, Temperature, Moisture, Scalar, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V ,W 
    real(knd), dimension(-1:,-1:,-1:), contiguous, intent(inout) :: Temperature, Moisture
    real(knd), dimension(-1:,-1:,-1:,1:), contiguous, intent(inout) :: Scalar
    real(knd), intent(in) :: eff_time
    real(knd) :: t_diff
    integer :: bi

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do bi = lbound(domain_bc_recv_buffers_copy,1), &
                ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(bi))

            if (b%relaxation) then

              t_diff = eff_time - b%time

              select case(b%direction)
                case (We)
                  call relax_domain_We(b)
                case (Ea)
                  call relax_domain_Ea(b)
                case (So)
                case (No)
                case (Bo)
                case (To)
              end select
            end if
          end associate
        end do
      end if

    end if
  contains

    subroutine relax_domain_We(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real, parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

      U(1,1:Uny,1:Unz) = ca(1) * U(1,1:Uny,1:Unz) + &
              cb(1) * (b%U(1,1:Uny,1:Unz) + t_diff * b%dU_dt(1,1:Uny,1:Unz))
      U(2,1:Uny,1:Unz) = ca(2) * U(2,1:Uny,1:Unz) + &
              cb(2) * (b%U(2,1:Uny,1:Unz) + t_diff * b%dU_dt(2,1:Uny,1:Unz))

      V(1,1:Vny,1:Vnz) = ca(1) * V(1,1:Vny,1:Vnz) + &
              cb(1) * (b%V(1,1:Vny,1:Vnz) + t_diff * b%dV_dt(1,1:Vny,1:Vnz))
      V(2,1:Vny,1:Vnz) = ca(2) * V(2,1:Vny,1:Vnz) + &
              cb(2) * (b%V(2,1:Vny,1:Vnz) + t_diff * b%dV_dt(2,1:Vny,1:Vnz))

      W(1,1:Wny,1:Wnz) = ca(1) * W(1,1:Wny,1:Wnz) + &
              cb(1) * (b%W(1,1:Wny,1:Wnz) + t_diff * b%dW_dt(1,1:Wny,1:Wnz))
      W(2,1:Wny,1:Wnz) = ca(2) * W(2,1:Wny,1:Wnz) + &
              cb(2) * (b%W(2,1:Wny,1:Wnz) + t_diff * b%dW_dt(2,1:Wny,1:Wnz))
    end subroutine

    subroutine relax_domain_Ea(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real, parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

      U(Unx,1:Uny,1:Unz) = ca(1) * U(Unx,1:Uny,1:Unz) + &
              cb(1) * (b%U(Unx,1:Uny,1:Unz) + t_diff * b%dU_dt(Unx,1:Uny,1:Unz))
      U(Unx-1,1:Uny,1:Unz) = ca(2) * U(Unx-1,1:Uny,1:Unz) + &
              cb(2) * (b%U(Unx-1,1:Uny,1:Unz) + t_diff * b%dU_dt(Unx-1,1:Uny,1:Unz))

      V(Vnx,1:Vny,1:Vnz) = ca(1) * V(Vnx,1:Vny,1:Vnz) + &
              cb(1) * (b%V(Vnx,1:Vny,1:Vnz) + t_diff * b%dV_dt(Vnx,1:Vny,1:Vnz))
      V(Vnx-1,1:Vny,1:Vnz) = ca(2) * V(Vnx-1,1:Vny,1:Vnz) + &
              cb(2) * (b%V(Vnx-1,1:Vny,1:Vnz) + t_diff * b%dV_dt(Vnx-1,1:Vny,1:Vnz))

      W(Wnx,1:Wny,1:Wnz) = ca(1) * W(Wnx,1:Wny,1:Wnz) + &
              cb(1) * (b%W(Wnx,1:Wny,1:Wnz) + t_diff * b%dW_dt(Wnx,1:Wny,1:Wnz))
      W(Wnx-1,1:Wny,1:Wnz) = ca(2) * W(Wnx-1,1:Wny,1:Wnz) + &
              cb(2) * (b%W(Wnx-1,1:Wny,1:Wnz) + t_diff * b%dW_dt(Wnx-1,1:Wny,1:Wnz))
    end subroutine

  end subroutine par_domain_bound_relaxation

end module domains_bc_par
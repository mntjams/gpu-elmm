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

    !time at which the data is valid
    real(tim) :: time = -tiny(1.0_tim)

    !MPI communicator for the remote communication
    integer :: comm = MPI_COMM_NULL
    integer :: remote_rank = MPI_PROC_NULL
    !The buffers are transferred every `time_step_ratio` time steps
    integer :: time_step_ratio = 1 
  end type

  type, extends(dom_bc_buffer_copy) :: dom_bc_buffer_refined
    integer :: spatial_ratio = 3

    !receive buffer grid, fields are interpolated from this grid to the finer grid
    integer :: r_Ui1, r_Ui2, r_Vi1, r_Vi2, r_Wi1, r_Wi2, r_Pri1, r_Pri2
    integer :: r_Uj1, r_Uj2, r_Vj1, r_Vj2, r_Wj1, r_Wj2, r_Prj1, r_Prj2
    integer :: r_Uk1, r_Uk2, r_Vk1, r_Vk2, r_Wk1, r_Wk2, r_Prk1, r_Prk2

    real(knd) :: r_dx, r_dy, r_dz
    real(knd) :: r_x0, r_y0, r_z0

    real(knd), allocatable :: r_xU(:), r_yV(:), r_zW(:)
    real(knd), allocatable :: r_x(:), r_y(:), r_z(:)

    !receive buffers
    real(knd), allocatable, dimension(:,:,:) :: r_U, r_V, r_W, r_Pr, r_Temperature, r_Moisture
    real(knd), allocatable, dimension(:,:,:,:) :: r_Scalar
    real(knd), allocatable, dimension(:,:,:) :: r_dU_dt, r_dV_dt, r_dW_dt, r_dPr_dt, &
                                                r_dTemperature_dt, r_dMoisture_dt
    real(knd), allocatable, dimension(:,:,:,:) :: r_dScalar_dt    
  end type


  !send buffers should be indexed from 1, there may be more of them from several nested domains
  !they are not accessed by the boundary conditions routines
  type(dom_bc_buffer_copy), allocatable :: domain_bc_send_buffers(:)
  !receive buffers should be indexed with the index of the domain side (West to Top)
  ! to ease finding the right buffer to a given nested boundary 
  type(dom_bc_buffer_copy), allocatable :: domain_bc_recv_buffers_copy(:)

  !receive buffers should be indexed with the index of the domain side (West to Top)
  ! to ease finding the right buffer to a given nested boundary 
  type(dom_bc_buffer_refined), allocatable :: domain_bc_recv_buffers(:)

contains



  subroutine par_init_domain_boundary_conditions
    integer :: Ui1, Ui2, Vi1, Vi2, Wi1, Wi2, Pri1, Pri2
    integer :: Uj1, Uj2, Vj1, Vj2, Wj1, Wj2, Prj1, Prj2
    integer :: Uk1, Uk2, Vk1, Vk2, Wk1, Wk2, Prk1, Prk2
    !temporary testing hack

#include "custom_domains_setup_refined.f90"

  end subroutine



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
              call recv_arrays_copy(b)
              b%time = time
            end if

          end associate
        end do
      end if


      if (allocated(domain_bc_recv_buffers)) then
        do i = lbound(domain_bc_recv_buffers,1), &
               ubound(domain_bc_recv_buffers,1)
          associate(b => domain_bc_recv_buffers(i))
  
            if (b%enabled) then
              call recv_arrays(b)
              b%time = time
            end if

          end associate
        end do
      end if


      call MPI_Waitall(size(requests), requests, MPI_STATUSES_IGNORE, ie)
      if (ie/=0) call error_stop("Error, MPI_Waitall in par_exchange_domain_bounds returns", ie)


      if (allocated(domain_bc_recv_buffers)) then
        do i = lbound(domain_bc_recv_buffers,1), &
               ubound(domain_bc_recv_buffers,1)
          associate(b => domain_bc_recv_buffers(i))
  
            if (b%enabled) then
              call par_interpolate_buffers(b)
            end if

          end associate
        end do
      end if
      
    end if

  contains

    subroutine send_arrays(b)
      type(dom_bc_buffer_copy), intent(inout) :: b
      
      call send(b%U, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 1)
      call send(b%V, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 2)
      call send(b%W, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 3)

      call send(b%dU_dt, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 1)
      call send(b%dV_dt, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 2)
      call send(b%dW_dt, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 3)

      if (enable_buoyancy) then
        call send(b%Temperature, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 4)
        call send(b%dTemperature_dt, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 4)
      end if

      if (enable_moisture) then
        call send(b%Moisture, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 5)
        call send(b%dMoisture_dt, b%remote_rank, b%comm, b%remote_domain*100 + b%direction*10 + 5)
      end if

    end subroutine

    subroutine recv_arrays_copy(b)
      type(dom_bc_buffer_copy), intent(inout) :: b
      real(knd) :: avg
    
      call recv(b%U, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 1)   
      call recv(b%V, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 2)     
      call recv(b%W, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 3)    

      call recv(b%dU_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 1)   
      call recv(b%dV_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 2)     
      call recv(b%dW_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 3)

      if (enable_buoyancy) then
        call recv(b%Temperature, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 4)
        call recv(b%dTemperature_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 4)
      end if

      if (enable_moisture) then
        call recv(b%Moisture, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 5)
        call recv(b%dMoisture_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 5)
      end if
    end subroutine

    subroutine recv_arrays(b)
      type(dom_bc_buffer_refined), intent(inout) :: b
      real(knd) :: avg
    
      call recv(b%r_U, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 1)   
      call recv(b%r_V, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 2)     
      call recv(b%r_W, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 3)    

      call recv(b%r_dU_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 1)   
      call recv(b%r_dV_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 2)     
      call recv(b%r_dW_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 3)

      if (enable_buoyancy) then
        call recv(b%r_Temperature, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 4)
        call recv(b%r_dTemperature_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 4)
      end if

      if (enable_moisture) then
        call recv(b%r_Moisture, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 5)
        call recv(b%r_dMoisture_dt, b%remote_rank, b%comm, domain_index*100 + b%direction*10 + 5)
      end if

    end subroutine

    subroutine send(a, rank, comm, tag)
      real(knd), intent(in), contiguous :: a(:,:,:)
      integer, intent(in) :: rank, comm, tag
      integer :: request, ie

      call MPI_ISend(a, size(a), MPI_KND, rank, tag, comm, request, ie)
      if (ie/=0) call error_stop("error sending MPI message.")

      requests = [requests, request]
    end subroutine

    subroutine recv(a, rank, comm, tag)
      real(knd), intent(inout), contiguous :: a(:,:,:)
      integer, intent(in) :: rank, comm, tag
      integer :: request, ie

      call MPI_IRecv(a, size(a), MPI_KND, rank, tag, comm, request, ie)
      if (ie/=0) call error_stop("error receiving MPI message.")

      requests = [requests, request]
    end subroutine

  end subroutine par_exchange_domain_bounds




  subroutine par_interpolate_buffers(b)
    type(dom_bc_buffer_refined), intent(inout) :: b

    call interpolate_U_spline(b%r_U, b%U)
    call interpolate_U_spline(b%r_dU_dt, b%dU_dt)

    call interpolate_V_spline(b%r_V, b%V)
    call interpolate_V_spline(b%r_dV_dt, b%dV_dt)

    call interpolate_W_spline(b%r_W, b%W)
    call interpolate_W_spline(b%r_dW_dt, b%dW_dt)


  contains

    subroutine interpolate_U_spline(in, out)
      use bspline_oo_module
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout), allocatable :: out(:,:,:)
      integer :: i, j, k
      type(bspline_3d) :: spl
      integer :: iflag, idx, idy, idz
      real(real64) :: val

      idx = 0; idy = 0; idz = 0
      
      call spl%initialize(b%r_xU(b%r_Ui1:b%r_Ui2), &
                          b%r_y(b%r_Uj1:b%r_Uj2), &
                          b%r_z(b%r_Uk1:b%r_Uk2), &
                          in, &
                          3, 3, 3, iflag)
      if (iflag/=1) then
        write(*,*) "Error in spline initialization, iflag:",iflag
        stop
      end if
      
      
      do k = b%Uk1, b%Uk2
        do j = b%Uj1, b%Uj2
          do i = b%Ui1, b%Ui2
            call spl%evaluate(xU(i), yPr(j), zPr(k), &
                              idx, idy, idz, val, iflag)
            if (iflag/=0) then
              write(*,*) "Error in spline interpolation at", &
                          xU(i), yPr(j), zPr(k)
              call error_stop()
            end if

            out(i,j,k) = val
          end do
        end do
      end do

      call spl%destroy()
    end subroutine

    subroutine interpolate_V_spline(in, out)
      use bspline_oo_module
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout), allocatable :: out(:,:,:)
      integer :: i, j, k
      type(bspline_3d) :: spl
      integer :: iflag, idx, idy, idz
      real(real64) :: val

      idx = 0; idy = 0; idz = 0
      
      call spl%initialize(b%r_x(b%r_Vi1:b%r_Vi2), &
                          b%r_yV(b%r_Vj1:b%r_Vj2), &
                          b%r_z(b%r_Vk1:b%r_Vk2), &
                          in, &
                          3, 3, 3, iflag)
      if (iflag/=1) then
        write(*,*) "Error in spline initialization, iflag:",iflag
        stop
      end if
      
      
      do k = b%Vk1, b%Vk2
        do j = b%Vj1, b%Vj2
          do i = b%Vi1, b%Vi2
            call spl%evaluate(xPr(i), yV(j), zPr(k), &
                              idx, idy, idz, val, iflag)
            if (iflag/=0) then
              write(*,*) "Error in spline interpolation at", &
                          xPr(i), yV(j), zPr(k)
              call error_stop()
            end if

            out(i,j,k) = val
          end do
        end do
      end do

      call spl%destroy()
    end subroutine

    subroutine interpolate_W_spline(in, out)
      use bspline_oo_module
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout), allocatable :: out(:,:,:)
      integer :: i, j, k
      type(bspline_3d) :: spl
      integer :: iflag, idx, idy, idz
      real(real64) :: val

      idx = 0; idy = 0; idz = 0
      
      call spl%initialize(b%r_x(b%r_Wi1:b%r_Wi2), &
                          b%r_y(b%r_Wj1:b%r_Wj2), &
                          b%r_zW(b%r_Wk1:b%r_Wk2), &
                          b%r_W, &
                          3, 3, 3, iflag)
      if (iflag/=1) then
        write(*,*) "Error in spline initialization, iflag:",iflag
        stop
      end if
      
      
      do k = b%Wk1, b%Wk2
        do j = b%Wj1, b%Wj2
          do i = b%Wi1, b%Wi2
            call spl%evaluate(xPr(i), yPr(j), zW(k), &
                              idx, idy, idz, val, iflag)
            if (iflag/=0) then
              write(*,*) "Error in spline interpolation at", &
                          xPr(i), yPr(j), zW(k)
              call error_stop()
            end if

            out(i,j,k) = val
          end do
        end do
      end do

      call spl%destroy()
    end subroutine

  end subroutine











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
            if (b%enabled) then

              U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = b%U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2)                 

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

            end if
          end associate
        end do
      end if

      if (allocated(domain_bc_recv_buffers)) then
        do bi = lbound(domain_bc_recv_buffers,1), &
                ubound(domain_bc_recv_buffers,1)
          associate(b => domain_bc_recv_buffers(bi))
            if (b%enabled) then

              U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2) = b%U(b%bUi1:b%bUi2,b%bUj1:b%bUj2,b%bUk1:b%bUk2)                 

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

    call par_update_domain_bounds_UVW(U, V, W, eff_time)

    call par_update_domain_bounds_temperature(Temperature, eff_time)

    call par_update_domain_bounds_moisture(Moisture, eff_time)

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
                  call relax_domain_copy_We(b)
                case (Ea)
                  call relax_domain_copy_Ea(b)
                case (So)
                  call relax_domain_copy_So(b)
                case (No)
                  call relax_domain_copy_No(b)
                case (Bo)
                case (To)
                  call relax_domain_copy_To(b)
              end select
            end if
          end associate
        end do
      end if

     if (allocated(domain_bc_recv_buffers)) then
        do bi = lbound(domain_bc_recv_buffers,1), &
                ubound(domain_bc_recv_buffers,1)
          associate(b => domain_bc_recv_buffers(bi))

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

    subroutine relax_domain_copy_We(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real(knd), parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

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

    subroutine relax_domain_copy_Ea(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real(knd), parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

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

    subroutine relax_domain_copy_So(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real(knd), parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

      U(1:Unx,1,1:Unz) = ca(1) * U(1:Unx,1,1:Unz) + &
              cb(1) * (b%U(1:Unx,1,1:Unz) + t_diff * b%dU_dt(1:Unx,1,1:Unz))
      U(1:Unx,2,1:Unz) = ca(2) * U(1:Unx,2,1:Unz) + &
              cb(2) * (b%U(1:Unx,2,1:Unz) + t_diff * b%dU_dt(1:Unx,2,1:Unz))

      V(1:Vnx,1,1:Vnz) = ca(1) * V(1:Vnx,1,1:Vnz) + &
              cb(1) * (b%V(1:Vnx,1,1:Vnz) + t_diff * b%dV_dt(1:Vnx,1,1:Vnz))
      V(1:Vnx,2,1:Vnz) = ca(2) * V(1:Vnx,2,1:Vnz) + &
              cb(2) * (b%V(1:Vnx,2,1:Vnz) + t_diff * b%dV_dt(1:Vnx,2,1:Vnz))

      W(1:Wnx,1,1:Wnz) = ca(1) * W(1:Wnx,1,1:Wnz) + &
              cb(1) * (b%W(1:Wnx,1,1:Wnz) + t_diff * b%dW_dt(1:Wnx,1,1:Wnz))
      W(1:Wnx,2,1:Wnz) = ca(2) * W(1:Wnx,2,1:Wnz) + &
              cb(2) * (b%W(1:Wnx,2,1:Wnz) + t_diff * b%dW_dt(1:Wnx,2,1:Wnz))
    end subroutine

    subroutine relax_domain_copy_No(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real(knd), parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

      U(1:Unx,Uny,1:Unz) = ca(1) * U(1:Unx,Uny,1:Unz) + &
              cb(1) * (b%U(1:Unx,Uny,1:Unz) + t_diff * b%dU_dt(1:Unx,Uny,1:Unz))
      U(1:Unx,Uny-1,1:Unz) = ca(2) * U(1:Unx,Uny-1,1:Unz) + &
              cb(2) * (b%U(1:Unx,Uny-1,1:Unz) + t_diff * b%dU_dt(1:Unx,Uny-1,1:Unz))

      V(1:Vnx,Vny,1:Vnz) = ca(1) * V(1:Vnx,Vny,1:Vnz) + &
              cb(1) * (b%V(1:Vnx,Vny,1:Vnz) + t_diff * b%dV_dt(1:Vnx,Vny,1:Vnz))
      V(1:Vnx,Vny-1,1:Vnz) = ca(2) * V(1:Vnx,Vny-1,1:Vnz) + &
              cb(2) * (b%V(1:Vnx,Vny-1,1:Vnz) + t_diff * b%dV_dt(1:Vnx,Vny-1,1:Vnz))

      W(1:Wnx,Wny,1:Wnz) = ca(1) * W(1:Wnx,Wny,1:Wnz) + &
              cb(1) * (b%W(1:Wnx,Wny,1:Wnz) + t_diff * b%dW_dt(1:Wnx,Wny,1:Wnz))
      W(1:Wnx,Wny-1,1:Wnz) = ca(2) * W(1:Wnx,Wny-1,1:Wnz) + &
              cb(2) * (b%W(1:Wnx,Wny-1,1:Wnz) + t_diff * b%dW_dt(1:Wnx,Wny-1,1:Wnz))
    end subroutine

    subroutine relax_domain_copy_To(b)
      type(dom_bc_buffer_copy), intent(in) :: b
      real(knd), parameter, dimension(*) :: ca = [.3_knd,.6_knd], cb = 1._knd - ca

      U(1:Unx,1:Uny,Unz) = ca(1) * U(1:Unx,1:Uny,Unz) + &
              cb(1) * (b%U(1:Unx,1:Uny,Unz) + t_diff * b%dU_dt(1:Unx,1:Uny,Unz))
      U(1:Unx,1:Uny,Unz-1) = ca(2) * U(1:Unx,1:Uny,Unz-1) + &
              cb(2) * (b%U(1:Unx,1:Uny,Unz-1) + t_diff * b%dU_dt(1:Unx,1:Uny,Unz-1))

      V(1:Vnx,1:Vny,Vnz) = ca(1) * V(1:Vnx,1:Vny,Vnz) + &
              cb(1) * (b%V(1:Vnx,1:Vny,Vnz) + t_diff * b%dV_dt(1:Vnx,1:Vny,Vnz))
      V(1:Vnx,1:Vny,Vnz-1) = ca(2) * V(1:Vnx,1:Vny,Vnz-1) + &
              cb(2) * (b%V(1:Vnx,1:Vny,Vnz-1) + t_diff * b%dV_dt(1:Vnx,1:Vny,Vnz-1))

      W(1:Wnx,1:Wny,Wnz) = ca(1) * W(1:Wnx,1:Wny,Wnz) + &
              cb(1) * (b%W(1:Wnx,1:Wny,Wnz) + t_diff * b%dW_dt(1:Wnx,1:Wny,Wnz))
      W(1:Wnx,1:Wny,Wnz-1) = ca(2) * W(1:Wnx,1:Wny,Wnz-1) + &
              cb(2) * (b%W(1:Wnx,1:Wny,Wnz-1) + t_diff * b%dW_dt(1:Wnx,1:Wny,Wnz-1))
    end subroutine



    subroutine relax_domain_We(b)
      type(dom_bc_buffer_refined), intent(in) :: b
      real(knd), parameter, dimension(*) :: cb = [.12_knd,.08_knd,.05_knd,.03_knd,.02_knd,.01_knd], &
                                            ca = 1._knd - cb
      integer, parameter :: width = size(ca)
      integer :: i

      do i = 1, width
        U(i,1:Uny,1:Unz) = ca(i) * U(i,1:Uny,1:Unz) + &
                cb(i) * (b%U(i,1:Uny,1:Unz) + t_diff * b%dU_dt(i,1:Uny,1:Unz))


        V(i,1:Vny,1:Vnz) = ca(i) * V(i,1:Vny,1:Vnz) + &
                cb(i) * (b%V(i,1:Vny,1:Vnz) + t_diff * b%dV_dt(i,1:Vny,1:Vnz))


        W(i,1:Wny,1:Wnz) = ca(i) * W(i,1:Wny,1:Wnz) + &
                cb(i) * (b%W(i,1:Wny,1:Wnz) + t_diff * b%dW_dt(i,1:Wny,1:Wnz))
      end do
    end subroutine

    subroutine relax_domain_Ea(b)
      type(dom_bc_buffer_refined), intent(in) :: b
      real(knd), parameter, dimension(*) :: cb = [.12_knd,.08_knd,.05_knd,.03_knd,.02_knd,.01_knd], &
                                            ca = 1._knd - cb
      integer, parameter :: width = size(ca)
      integer :: i, bi

      do bi = 1, width
        i = Unx - bi + 1
        U(i,1:Uny,1:Unz) = ca(bi) * U(i,1:Uny,1:Unz) + &
                cb(bi) * (b%U(i,1:Uny,1:Unz) + t_diff * b%dU_dt(i,1:Uny,1:Unz))
      end do

      do i = Vnx-width+1, Vnx
        V(i,1:Vny,1:Vnz) = ca(bi) * V(i,1:Vny,1:Vnz) + &
                cb(i) * (b%V(i,1:Vny,1:Vnz) + t_diff * b%dV_dt(i,1:Vny,1:Vnz))
      end do

      do i = Wnx-width+1, Wnx
        W(i,1:Wny,1:Wnz) = ca(bi) * W(i,1:Wny,1:Wnz) + &
                cb(bi) * (b%W(i,1:Wny,1:Wnz) + t_diff * b%dW_dt(i,1:Wny,1:Wnz))
      end do
    end subroutine


  end subroutine par_domain_bound_relaxation

end module domains_bc_par
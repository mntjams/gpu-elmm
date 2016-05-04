module domains_bc_par

  use mpi, only: MPI_PROC_NULL, MPI_COMM_NULL, MPI_STATUS_SIZE
  use Parameters
  use custom_mpi
  use custom_par
  use TurbInlet

  implicit none

  private

  public par_init_domain_boundary_conditions, &
         par_exchange_domain_bounds, &
         par_update_domain_bounds, &
         par_update_pr_gradient, &
         par_update_domain_bounds_UVW, &
         par_update_domain_bounds_temperature, &
         par_update_domain_bounds_moisture, &
         par_domain_bound_relaxation, &
         domain_spatial_ratio, &
         domain_time_step_ratio

  type dom_bc_buffer_copy
    logical :: enabled = .false.

    logical :: rescale_compatibility = .false.
 
    logical :: relaxation = .true.
    real(knd) :: relax_factor = 1

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

    !whether child gets the pressure gradient from parent
    !makes sence when parent changes pressure-gradient dynamically
    logical :: exchange_pr_gradient_x = .false.
    logical :: exchange_pr_gradient_y = .false.

    real(knd) :: pr_gradient_x = 0
    real(knd) :: pr_gradient_y = 0
    real(knd) :: pr_gradient_x_dt = 0
    real(knd) :: pr_gradient_y_dt = 0

    !time at which the data is valid
    real(tim) :: time = -tiny(1.0_tim)

    !MPI communicator for the remote communication
    integer :: comm = MPI_COMM_NULL
    integer :: remote_rank = MPI_PROC_NULL
    !The buffers are transferred every `time_step_ratio` time steps
    integer :: time_step_ratio = 1
    !whether the interpolation from the receive buffers is necessary in this tim-step
    logical :: interpolate = .true.

    integer :: time_step = 0
  end type

  type, extends(dom_bc_buffer_copy) :: dom_bc_buffer_refined
    !the ratio of grid resolutions between parent and child domains
    integer :: spatial_ratio = 3
    !order of accuracy of the spatial interpolation in the respective directions (2..linear, 3..parabolic)
    integer :: interp_order(3) = 2

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


  type, extends(dom_bc_buffer_refined) :: dom_bc_buffer_turbulence_generator
    logical :: turb_generator_enabled = .false.
    type(turbulence_generator_nesting), allocatable :: turb_generator
    real(knd), allocatable, dimension(:,:) :: U_turb, V_turb, W_turb
  contains
    procedure :: compute_sgs_tke => dom_bc_buffer_turbulence_generator_compute_sgs_tke
  end type


  integer :: domain_spatial_ratio = 3

  integer :: domain_time_step_ratio = 3


  !send buffers should be indexed from 1, there may be more of them from several nested domains
  !they are not accessed by the boundary conditions routines
  type(dom_bc_buffer_copy), allocatable :: domain_bc_send_buffers(:)
  !receive buffers should be indexed with the index of the domain side (West to Top)
  ! to ease finding the right buffer to a given nested boundary 
  type(dom_bc_buffer_copy), allocatable :: domain_bc_recv_buffers_copy(:)

  !receive buffers should be indexed with the index of the domain side (West to Top)
  ! to ease finding the right buffer to a given nested boundary 
  type(dom_bc_buffer_turbulence_generator), allocatable :: domain_bc_recv_buffers(:)

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

    integer, allocatable :: requests(:)
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

              if (b%exchange_pr_gradient_x) then
                b%pr_gradient_x_dt = (pr_gradient_x - b%pr_gradient_x) / dt
                b%pr_gradient_x = pr_gradient_x
                call send_scalar(b, 11, b%pr_gradient_x)
                call send_scalar(b, 12, b%pr_gradient_x_dt)
              end if

              if (b%exchange_pr_gradient_x) then
                b%pr_gradient_y_dt = (pr_gradient_y - b%pr_gradient_y) / dt
                b%pr_gradient_y = pr_gradient_y
                call send_scalar(b, 13, b%pr_gradient_y)
                call send_scalar(b, 14, b%pr_gradient_y_dt)
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

              b%interpolate = .false.

              !for 0 and 1, initialization and the first time-step
              if (b%time_step<=1) then
                call recv_arrays(b)

                if (b%exchange_pr_gradient_x) then
                  call recv_scalar(b, 11, b%pr_gradient_x)
                  call recv_scalar(b, 12, b%pr_gradient_x_dt)
                end if

                if (b%exchange_pr_gradient_x) then
                  call recv_scalar(b, 13, b%pr_gradient_y)
                  call recv_scalar(b, 14, b%pr_gradient_y_dt)
                end if

                b%time = time
                b%interpolate = .true.
              end if
              if (b%time_step==b%time_step_ratio) then
                b%time_step=1
              else
                b%time_step = b%time_step + 1
              end if
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
  
            if (b%enabled.and.b%interpolate) then
              call par_interpolate_buffers(b)
            end if

            if (b%turb_generator_enabled) then
              if (b%interpolate) call b%compute_sgs_tke
              call b%turb_generator%time_step(b%U_turb, b%V_turb, b%V_turb, dt)

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
      class(dom_bc_buffer_refined), intent(inout) :: b
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

      !ISsend because we do want synchronisation here. Otherwise a fast domain
      !could run too quickly, buffer too many requests and fill up the memory
      call MPI_ISsend(a, size(a), MPI_KND, rank, tag, comm, request, ie)
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

    subroutine send_scalar(b, tag_base, a)
      type(dom_bc_buffer_copy), intent(inout) :: b
      integer, intent(in) :: tag_base
      real(knd), intent(in) :: a
      integer :: request, ie
      
      call MPI_ISend(a, 1, MPI_KND, &
                     b%remote_rank, &
                     b%remote_domain*100 + b%direction*10 + tag_base, &
                     b%comm, &
                     request, ie)
      if (ie/=0) call error_stop("error sending MPI message.")

      requests = [requests, request]
    end subroutine

    subroutine recv_scalar(b, tag_base, a)
      class(dom_bc_buffer_refined), intent(inout) :: b
      integer, intent(in) :: tag_base
      real(knd), intent(inout) :: a
      integer :: request, ie
      
      call MPI_IRecv(a, 1, MPI_KND, &
                      b%remote_rank, &
                      domain_index*100 + b%direction*10 + tag_base, &
                      b%comm, &
                      request, ie)
      if (ie/=0) call error_stop("error sending MPI message.")

      requests = [requests, request]
    end subroutine

  end subroutine par_exchange_domain_bounds




  subroutine par_interpolate_buffers(b)
    class(dom_bc_buffer_refined), intent(inout) :: b

    if (all(b%interp_order==2)) then
      call interpolate_U_trilinear(b%r_U, b%U)
      call interpolate_U_trilinear(b%r_dU_dt, b%dU_dt)

      call interpolate_V_trilinear(b%r_V, b%V)
      call interpolate_V_trilinear(b%r_dV_dt, b%dV_dt)

      call interpolate_W_trilinear(b%r_W, b%W)
      call interpolate_W_trilinear(b%r_dW_dt, b%dW_dt)
    else
      call interpolate_U_spline(b%r_U, b%U)
      call interpolate_U_spline(b%r_dU_dt, b%dU_dt)

      call interpolate_V_spline(b%r_V, b%V)
      call interpolate_V_spline(b%r_dV_dt, b%dV_dt)

      call interpolate_W_spline(b%r_W, b%W)
      call interpolate_W_spline(b%r_dW_dt, b%dW_dt)
    end if

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
                          b%interp_order(1), b%interp_order(2), b%interp_order(3), &
                          iflag)
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
              write(*,*) "Error in spline interpolation of U at", &
                          xU(i), yPr(j), zPr(k)
              write(*,*) "flag:",iflag
              write(*,*) "xs:",b%r_xU(b%r_Ui1:b%r_Ui2)
              write(*,*) "ys:",b%r_y(b%r_Uj1:b%r_Uj2)
              write(*,*) "zs:",b%r_z(b%r_Uk1:b%r_Uk2)
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
                          b%interp_order(1), b%interp_order(2), b%interp_order(3), &
                          iflag)
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
              write(*,*) "Error in spline interpolation of V at", &
                          xPr(i), yV(j), zPr(k)
              write(*,*) "flag:",iflag
              write(*,*) "xs:",b%r_x(b%r_Vi1:b%r_Vi2)
              write(*,*) "ys:",b%r_yV(b%r_Vj1:b%r_Vj2)
              write(*,*) "zs:",b%r_z(b%r_Vk1:b%r_Vk2)
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
                          b%interp_order(1), b%interp_order(2), b%interp_order(3), &
                          iflag)
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
              write(*,*) "Error in spline interpolation of W at", &
                          xPr(i), yPr(j), zW(k)
              write(*,*) "flag:",iflag
              write(*,*) "xs:",b%r_x(b%r_Wi1:b%r_Wi2)
              write(*,*) "ys:",b%r_y(b%r_Wj1:b%r_Wj2)
              write(*,*) "zs:",b%r_zW(b%r_Wk1:b%r_Wk2)
              call error_stop()
            end if

            out(i,j,k) = val
          end do
        end do
      end do

      call spl%destroy()
    end subroutine

    subroutine interpolate_U_trilinear(in, out)
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout), allocatable :: out(:,:,:)
      integer :: xi, yj, zk
      integer :: i, j, k
      
      !$omp parallel do private(i, j, k, xi, yj, zk)
      do k = b%Uk1, b%Uk2
        do j = b%Uj1, b%Uj2
          do i = b%Ui1, b%Ui2
            call U_r_index(xU(i), yPr(j), zPr(k), xi, yj, zk)

            out(i,j,k) = &
              TriLinInt((xU(i)   - b%r_xU(xi)) / b%r_dx, &
                        (yPr(j)  - b%r_y(yj) ) / b%r_dy, &
                        (zPr(k)  - b%r_z(zk) ) / b%r_dz, &
                        in(xi  , yj  , zk  ), &
                        in(xi+1, yj  , zk  ), &
                        in(xi  , yj+1, zk  ), &
                        in(xi  , yj  , zk+1), &
                        in(xi+1, yj+1, zk  ), &
                        in(xi+1, yj  , zk+1), &
                        in(xi  , yj+1, zk+1), &
                        in(xi+1, yj+1, zk+1))
          end do
        end do
      end do
    end subroutine

    pure subroutine U_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_xU,1)
      ly = lbound(b%r_y,1)
      lz = lbound(b%r_z,1)

      xi = min(max(floor( (x - b%r_xU(lx))/b%r_dx ), 0) + lx, ubound(b%r_xU,1)-1)
      yj = min(max(floor( (y - b%r_y(ly))/b%r_dy ), 0) + ly, ubound(b%r_y,1)-1)
      zk = min(max(floor( (z - b%r_z(lz))/b%r_dz ), 0) + lz, ubound(b%r_z,1)-1)
    end subroutine

    subroutine interpolate_V_trilinear(in, out)
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout), allocatable :: out(:,:,:)
      integer :: xi, yj, zk
      integer :: i, j, k
      
      !$omp parallel do private(i, j, k, xi, yj, zk)
      do k = b%Vk1, b%Vk2
        do j = b%Vj1, b%Vj2
          do i = b%Vi1, b%Vi2
            call V_r_index(xPr(i), yV(j), zPr(k), xi, yj, zk)

            out(i,j,k) = &
              TriLinInt((xPr(i) - b%r_x(xi) ) / b%r_dx, &
                        (yV(j)  - b%r_yV(yj)) / b%r_dy, &
                        (zPr(k) - b%r_z(zk) ) / b%r_dz, &
                        in(xi  , yj  , zk  ), &
                        in(xi+1, yj  , zk  ), &
                        in(xi  , yj+1, zk  ), &
                        in(xi  , yj  , zk+1), &
                        in(xi+1, yj+1, zk  ), &
                        in(xi+1, yj  , zk+1), &
                        in(xi  , yj+1, zk+1), &
                        in(xi+1, yj+1, zk+1))
          end do
        end do
      end do
    end subroutine

    pure subroutine V_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_x,1)
      ly = lbound(b%r_yV,1)
      lz = lbound(b%r_z,1)

      xi = min(max(floor( (x - b%r_x(lx))/b%r_dx ), 0) + lx, ubound(b%r_x,1)-1)
      yj = min(max(floor( (y - b%r_yV(ly))/b%r_dy ), 0) + ly, ubound(b%r_yV,1)-1)
      zk = min(max(floor( (z - b%r_z(lz))/b%r_dz ), 0) + lz, ubound(b%r_z,1)-1)
    end subroutine

    subroutine interpolate_W_trilinear(in, out)
      real(knd), intent(in), allocatable :: in(:,:,:)
      real(knd), intent(inout), allocatable :: out(:,:,:)
      integer :: xi, yj, zk
      integer :: i, j, k
      
      !$omp parallel do private(i, j, k, xi, yj, zk)
      do k = b%Wk1, b%Wk2
        do j = b%Wj1, b%Wj2
          do i = b%Wi1, b%Wi2
            call W_r_index(xPr(i), yPr(j), zW(k), xi, yj, zk)

            out(i,j,k) = &
              TriLinInt((xPr(i) - b%r_x(xi) ) / b%r_dx, &
                        (yV(j)  - b%r_y(yj) ) / b%r_dy, &
                        (zPr(k) - b%r_zW(zk)) / b%r_dz, &
                        in(xi  , yj  , zk  ), &
                        in(xi+1, yj  , zk  ), &
                        in(xi  , yj+1, zk  ), &
                        in(xi  , yj  , zk+1), &
                        in(xi+1, yj+1, zk  ), &
                        in(xi+1, yj  , zk+1), &
                        in(xi  , yj+1, zk+1), &
                        in(xi+1, yj+1, zk+1))
          end do
        end do
      end do
    end subroutine

    pure subroutine W_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_x,1)
      ly = lbound(b%r_y,1)
      lz = lbound(b%r_zW,1)

      xi = min(max(floor( (x - b%r_x(lx))/b%r_dx ), 0) + lx, ubound(b%r_x,1)-1)
      yj = min(max(floor( (y - b%r_y(ly))/b%r_dy ), 0) + ly, ubound(b%r_y,1)-1)
      zk = min(max(floor( (z - b%r_zW(lz))/b%r_dz ), 0) + lz, ubound(b%r_zW,1)-1)
    end subroutine

    pure real(knd) function TriLinInt(a, b, c, &
                                     val000, val100, val010, val001, val110, val101, val011, val111)
      real(knd), intent(in) :: a, b, c
      real(knd), intent(in) :: val000, val100, val010, val001, val110, val101, val011, val111

      TriLinInt =  (1-a) * (1-b) * (1-c) * val000 + &
                   a     * (1-b) * (1-c) * val100 + &
                   (1-a) * b     * (1-c) * val010 + &
                   (1-a) * (1-b) * c     * val001 + &
                   a     * b     * (1-c) * val110 + &
                   a     * (1-b) * c     * val101 + &
                   (1-a) * b     * c     * val011 + &
                   a     * b     * c     * val111
    end function TriLinInt


  end subroutine











  subroutine par_update_domain_bounds_UVW(U, V, W, eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U, V, W
    real(knd), intent(in) :: eff_time
    integer :: bi
    real(knd) :: S, t_diff
    integer :: i, j, k

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

      if (allocated(domain_bc_recv_buffers)) then
        do bi = lbound(domain_bc_recv_buffers,1), &
                ubound(domain_bc_recv_buffers,1)
          associate(b => domain_bc_recv_buffers(bi))
            if (b%enabled.and.b%turb_generator_enabled) then
              select case (b%direction)
                case (We, Ea)
                  do k = b%bUk1, b%bUk2
                    do j = b%bUj1, b%bUj2
                      do i = b%bUi1, b%bUi2
                        U(i,j,k) = U(i,j,k) + b%U_turb(j,k)
                      end do
                    end do
                  end do
                  do k = b%bVk1, b%bVk2
                    do j = b%bVj1, b%bVj2
                      do i = b%bVi1, b%bVi2
                        V(i,j,k) = V(i,j,k) + b%V_turb(j,k)
                      end do
                    end do
                  end do
                  do k = b%bWk1, b%bWk2
                    do j = b%bWj1, b%bWj2
                      do i = b%bWi1, b%bWi2
                        W(i,j,k) = W(i,j,k) + b%W_turb(j,k)
                      end do
                    end do
                  end do
                case (So, No)
                  do k = b%bUk1, b%bUk2
                    do j = b%bUj1, b%bUj2
                      do i = b%bUi1, b%bUi2
                        U(i,j,k) = U(i,j,k) + b%U_turb(i,k)
                      end do
                    end do
                  end do
                  do k = b%bVk1, b%bVk2
                    do j = b%bVj1, b%bVj2
                      do i = b%bVi1, b%bVi2
                        V(i,j,k) = V(i,j,k) + b%V_turb(i,k)
                      end do
                    end do
                  end do
                  do k = b%bWk1, b%bWk2
                    do j = b%bWj1, b%bWj2
                      do i = b%bWi1, b%bWi2
                        W(i,j,k) = W(i,j,k) + b%W_turb(i,k)
                      end do
                    end do
                  end do
                case default
                  continue
              end select
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

  subroutine par_update_pr_gradient(eff_time)
    !effective time, because it can also reflect individual RK stages
    real(knd), intent(in) :: eff_time
    real(knd) :: t_diff
    integer :: bi

    if (enable_multiple_domains) then

     if (allocated(domain_bc_recv_buffers)) then
        do bi = lbound(domain_bc_recv_buffers,1), &
                ubound(domain_bc_recv_buffers,1)
          associate(b => domain_bc_recv_buffers(bi))
            if (b%enabled) then

              t_diff = eff_time - b%time

              if (b%exchange_pr_gradient_x) then
                pr_gradient_x = b%pr_gradient_x + b%pr_gradient_x_dt * t_diff
              end if
              if (b%exchange_pr_gradient_y) then
                pr_gradient_y = b%pr_gradient_y + b%pr_gradient_y_dt * t_diff
              end if

              if (b%exchange_pr_gradient_x.or.b%exchange_pr_gradient_y) exit
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

    call par_update_pr_gradient(eff_time)

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
    integer :: bi, width
    real(knd), allocatable :: ca(:), cb(:)

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

              width = b%spatial_ratio * 2

              cb = cb_table(width) * b%relax_factor
              where(cb>1) cb = 1
              ca = 1 - cb

              t_diff = eff_time - b%time

              select case(b%direction)
                case (We)
                  call relax_domain_We(b)
                case (Ea)
                  call relax_domain_Ea(b)
                case (So)
                  call relax_domain_So(b)
                case (No)
                  call relax_domain_No(b)
                case (Bo)
                  call relax_domain_Bo(b)
                case (To)
                  call relax_domain_To(b)
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
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i

      width = b%spatial_ratio * 2

      !$omp parallel do private(i)
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
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i, wi

      width = b%spatial_ratio * 2

      !$omp parallel do private(i, wi)
      do wi = 1, width
        i = Unx - wi + 1
        U(i,1:Uny,1:Unz) = ca(wi) * U(i,1:Uny,1:Unz) + &
                           cb(wi) * (b%U(i,1:Uny,1:Unz) + &
                                     t_diff * b%dU_dt(i,1:Uny,1:Unz))

        i = Vnx - wi + 1
        V(i,1:Vny,1:Vnz) = ca(wi) * V(i,1:Vny,1:Vnz) + &
                           cb(wi) * (b%V(i,1:Vny,1:Vnz) + &
                                     t_diff * b%dV_dt(i,1:Vny,1:Vnz))

        i = Wnx - wi + 1
        W(i,1:Wny,1:Wnz) = ca(wi) * W(i,1:Wny,1:Wnz) + &
                           cb(wi) * (b%W(i,1:Wny,1:Wnz) + &
                                     t_diff * b%dW_dt(i,1:Wny,1:Wnz))
      end do
    end subroutine

    subroutine relax_domain_So(b)
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i

      width = b%spatial_ratio * 2

      !$omp parallel do private(i)
      do i = 1, width
        U(1:Unx,i,1:Unz) = ca(i) * U(1:Unx,i,1:Unz) + &
                cb(i) * (b%U(1:Unx,i,1:Unz) + t_diff * b%dU_dt(1:Unx,i,1:Unz))


        V(1:Vnx,i,1:Vnz) = ca(i) * V(1:Vnx,i,1:Vnz) + &
                cb(i) * (b%V(1:Vnx,i,1:Vnz) + t_diff * b%dV_dt(1:Vnx,i,1:Vnz))


        W(1:Wnx,i,1:Wnz) = ca(i) * W(1:Wnx,i,1:Wnz) + &
                cb(i) * (b%W(1:Wnx,i,1:Wnz) + t_diff * b%dW_dt(1:Wnx,i,1:Wnz))
      end do
    end subroutine

    subroutine relax_domain_No(b)
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i, wi

      width = b%spatial_ratio * 2

      !$omp parallel do private(i, wi)
      do wi = 1, width
        i = Uny - wi + 1
        U(1:Unx,i,1:Unz) = ca(wi) * U(1:Unx,i,1:Unz) + &
                           cb(wi) * (b%U(1:Unx,i,1:Unz) + &
                                     t_diff * b%dU_dt(1:Unx,i,1:Unz))

        i = Vny - wi + 1
        V(1:Vnx,i,1:Vnz) = ca(wi) * V(1:Vnx,i,1:Vnz) + &
                           cb(wi) * (b%V(1:Vnx,i,1:Vnz) + &
                                     t_diff * b%dV_dt(1:Vnx,i,1:Vnz))

        i = Wny - wi + 1
        W(1:Wnx,i,1:Wnz) = ca(wi) * W(1:Wnx,i,1:Wnz) + &
                           cb(wi) * (b%W(1:Wnx,i,1:Wnz) + &
                                     t_diff * b%dW_dt(1:Wnx,i,1:Wnz))
      end do
    end subroutine

    subroutine relax_domain_Bo(b)
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i

      width = b%spatial_ratio * 2

      !$omp parallel do private(i)
      do i = 1, width
        U(1:Unx,1:Uny,i) = ca(i) * U(1:Unx,1:Uny,i) + &
                cb(i) * (b%U(1:Unx,1:Uny,i) + t_diff * b%dU_dt(1:Unx,1:Uny,i))


        V(1:Vnx,1:Vny,i) = ca(i) * V(1:Vnx,1:Vny,i) + &
                cb(i) * (b%V(1:Vnx,1:Vny,i) + t_diff * b%dV_dt(1:Vnx,1:Vny,i))


        W(1:Wnx,1:Wny,i) = ca(i) * W(1:Wnx,1:Wny,i) + &
                cb(i) * (b%W(1:Wnx,1:Wny,i) + t_diff * b%dW_dt(1:Wnx,1:Wny,i))
      end do
    end subroutine

    subroutine relax_domain_To(b)
      class(dom_bc_buffer_refined), intent(in) :: b
      integer :: width
      integer :: i, wi

      width = b%spatial_ratio * 2

      !$omp parallel do private(i, wi)
      do wi = 1, width
        i = Unz - wi + 1
        U(1:Unx,1:Uny,i) = ca(wi) * U(1:Unx,1:Uny,i) + &
                           cb(wi) * (b%U(1:Unx,1:Uny,i) + &
                                     t_diff * b%dU_dt(1:Unx,1:Uny,i))

        i = Vnz - wi + 1
        V(1:Vnx,1:Vny,i) = ca(wi) * V(1:Vnx,1:Vny,i) + &
                           cb(wi) * (b%V(1:Vnx,1:Vny,i) + &
                                     t_diff * b%dV_dt(1:Vnx,1:Vny,i))

        i = Wnz - wi + 1
        W(1:Wnx,1:Wny,i) = ca(wi) * W(1:Wnx,1:Wny,i) + &
                           cb(wi) * (b%W(1:Wnx,1:Wny,i) + &
                                     t_diff * b%dW_dt(1:Wnx,1:Wny,i))
      end do
    end subroutine

    pure function ca_table(width) result(res)
      integer, intent(in) :: width
      real(knd) :: res(width)
      integer :: i
     
      res = 1 - cb_table(width)
    end function

    pure function cb_table(width) result(res)
      integer, intent(in) :: width
      real(knd) :: res(width)
      integer :: i
     
      do i = 1, width
        res(i) = (1 - DampF((width - i + 0.5_knd)/width)) / 10
      end do
   end function

    pure function DampF(x) result(res)
      real(knd) :: res
      real(knd), intent(in) :: x

      if (x<=0) then
        res = 1
      else if (x>=1) then
        res = 0
      else
        res = (1 - 0.04_knd*x**2) * &
                ( 1 - (1 - exp(10._knd*x**2)) / (1 - exp(10._knd)) )
      end if
    end function

  end subroutine par_domain_bound_relaxation


  subroutine dom_bc_buffer_turbulence_generator_compute_sgs_tke(b)
    use Filters
    class(dom_bc_buffer_turbulence_generator), intent(inout) :: b
    real(knd), allocatable, dimension(:,:,:) :: Uf, Vf, Wf
    real(knd) :: uu
    integer :: i, j, k, xi, yj, zk

    select case (b%direction)
      case(We)
        allocate(Uf(b%r_Ui1+2-1:b%r_Ui1+2+1,b%r_Uj1:b%r_Uj2,b%r_Uk1:b%r_Uk2))
        allocate(Vf(b%r_Vi1+2-1:b%r_Vi1+2+1,b%r_Vj1:b%r_Vj2,b%r_Vk1:b%r_Vk2))
        allocate(Wf(b%r_Wi1+2-1:b%r_Wi1+2+1,b%r_Wj1:b%r_Wj2,b%r_Wk1:b%r_Wk2))

        call FilterTopHatSimple(Uf, b%r_U(b%r_Ui1+2-1:b%r_Ui1+2+1,b%r_Uj1:b%r_Uj2,b%r_Uk1:b%r_Uk2))
        call FilterTopHatSimple(Vf, b%r_V(b%r_Vi1+2-1:b%r_Vi1+2+1,b%r_Vj1:b%r_Vj2,b%r_Vk1:b%r_Vk2))
        call FilterTopHatSimple(Wf, b%r_W(b%r_Wi1+2-1:b%r_Wi1+2+1,b%r_Wj1:b%r_Wj2,b%r_Wk1:b%r_Wk2))

        Uf(b%r_Ui1+2,:,:) = (b%r_U(b%r_Ui1+2,:,:) - Uf(b%r_Ui1+2,:,:))**2
        Vf(b%r_Ui1+2,:,:) = (b%r_V(b%r_Vi1+2,:,:) - Vf(b%r_Vi1+2,:,:))**2
        Wf(b%r_Ui1+2,:,:) = (b%r_W(b%r_Wi1+2,:,:) - Wf(b%r_Wi1+2,:,:))**2

        do k = 1, Prnz
          do j = 1, Prny
            b%turb_generator%sgs_tke = 0

            call U_r_index(im_xmin, yPr(j), zPr(k), xi, yj, zk)
            uu = BiLinInt((yPr(j)  - b%r_y(yj) ) / b%r_dy, &
                          (zPr(k)  - b%r_z(zk) ) / b%r_dz, &
                          Uf(b%r_Ui1+2  , yj  , zk  ), &
                          Uf(b%r_Ui1+2  , yj+1, zk  ), &
                          Uf(b%r_Ui1+2  , yj  , zk+1), &
                          Uf(b%r_Ui1+2  , yj+1, zk+1))
            uu = max(uu,0._knd)
            b%turb_generator%sgs_tke = b%turb_generator%sgs_tke + uu

            

            call V_r_index(im_xmin, yV(j), zPr(k), xi, yj, zk)
            uu = BiLinInt((yV(j)   - b%r_y(yj) ) / b%r_dy, &
                          (zPr(k)  - b%r_z(zk) ) / b%r_dz, &
                          Vf(b%r_Vi1+2  , yj  , zk  ), &
                          Vf(b%r_Vi1+2  , yj+1, zk  ), &
                          Vf(b%r_Vi1+2  , yj  , zk+1), &
                          Vf(b%r_Vi1+2  , yj+1, zk+1))
            uu = max(uu,0._knd)
            b%turb_generator%sgs_tke = b%turb_generator%sgs_tke + uu

            call W_r_index(im_xmin, yPr(j), zW(k), xi, yj, zk)
            uu = BiLinInt((yPr(j)  - b%r_y(yj) ) / b%r_dy, &
                          (zW(k)   - b%r_z(zk) ) / b%r_dz, &
                          Wf(b%r_Wi1+2  , yj  , zk  ), &
                          Wf(b%r_Wi1+2  , yj+1, zk  ), &
                          Wf(b%r_Wi1+2  , yj  , zk+1), &
                          Wf(b%r_Wi1+2  , yj+1, zk+1))
            uu = max(uu,0._knd)
            b%turb_generator%sgs_tke = b%turb_generator%sgs_tke + uu

          end do
        end do

        call b%turb_generator%update_turbulence_profiles
            
    end select

  contains

    pure subroutine U_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_xU,1)
      ly = lbound(b%r_y,1)
      lz = lbound(b%r_z,1)

      xi = min(max(floor( (x - b%r_xU(lx))/b%r_dx ), 0) + lx, ubound(b%r_xU,1)-1)
      yj = min(max(floor( (y - b%r_y(ly))/b%r_dy ), 0) + ly, ubound(b%r_y,1)-1)
      zk = min(max(floor( (z - b%r_z(lz))/b%r_dz ), 0) + lz, ubound(b%r_z,1)-1)
    end subroutine

    pure subroutine V_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_x,1)
      ly = lbound(b%r_yV,1)
      lz = lbound(b%r_z,1)

      xi = min(max(floor( (x - b%r_x(lx))/b%r_dx ), 0) + lx, ubound(b%r_x,1)-1)
      yj = min(max(floor( (y - b%r_yV(ly))/b%r_dy ), 0) + ly, ubound(b%r_yV,1)-1)
      zk = min(max(floor( (z - b%r_z(lz))/b%r_dz ), 0) + lz, ubound(b%r_z,1)-1)
    end subroutine

    pure subroutine W_r_index(x, y, z, xi, yj, zk)
      real(knd), intent(in) :: x, y, z
      integer, intent(out) :: xi, yj, zk
      integer :: lx, ly, lz

      lx = lbound(b%r_x,1)
      ly = lbound(b%r_y,1)
      lz = lbound(b%r_zW,1)

      xi = min(max(floor( (x - b%r_x(lx))/b%r_dx ), 0) + lx, ubound(b%r_x,1)-1)
      yj = min(max(floor( (y - b%r_y(ly))/b%r_dy ), 0) + ly, ubound(b%r_y,1)-1)
      zk = min(max(floor( (z - b%r_zW(lz))/b%r_dz ), 0) + lz, ubound(b%r_zW,1)-1)
    end subroutine

    pure real(knd) function BiLinInt(a, b, &
                                     val00, val10, val01, val11)
      real(knd), intent(in) :: a, b
      real(knd), intent(in) :: val00, val10, val01, val11

      BiLinInt =  (1-a) * (1-b) * val00 + &
                   a     * (1-b) * val10 + &
                   (1-a) * b     * val01 + &
                   a     * b     * val11

    end function BiLinInt

  end subroutine


end module domains_bc_par
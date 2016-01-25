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
         par_update_domain_bounds_U

  type dom_bc_buffer_copy
    logical :: enabled = .false.

    logical :: rescale_compatibility = .false.
    logical :: rescale_dU_dt = .true.

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
  type(dom_bc_buffer_copy), allocatable :: domain_bc_send_buffers_copy(:)
  !receive buffers should be indexed with the index of the domain side (West to Top)
  ! to ease finding the right buffer to a given nested boundary 
  type(dom_bc_buffer_copy), allocatable :: domain_bc_recv_buffers_copy(:)

contains



  subroutine par_init_domain_boundary_conditions
    integer :: Ui1, Ui2, Vi1, Vi2, Wi1, Wi2, Pri1, Pri2
    integer :: Uj1, Uj2, Vj1, Vj2, Wj1, Wj2, Prj1, Prj2
    integer :: Uk1, Uk2, Vk1, Vk2, Wk1, Wk2, Prk1, Prk2
    !temporary testing hack

    if (enable_multiple_domains) then

      if (domain_index==1) then

        if (iim==nxims) then
          !send buffers should be counted from 1, there may be more of them from several nested domains
          allocate(domain_bc_send_buffers_copy(1))
          domain_bc_send_buffers_copy(1)%comm = world_comm
          domain_bc_send_buffers_copy(1)%remote_rank = domain_ranks_grid(2)%arr(1,jim,kim)
          domain_bc_send_buffers_copy(1)%remote_domain = 2
          domain_bc_send_buffers_copy(1)%direction = Ea
          domain_bc_send_buffers_copy(1)%enabled = .true.

          Ui1 = Unx-2
          Ui2 = Unx
          Uj1 = -2
          Uj2 = Uny+3
          Uk1 = -2
          Uk2 = Unz+3

          Vi1 = Vnx-2
          Vi2 = Vnx
          Vj1 = -2
          Vj2 = Vny+3
          Vk1 = -2
          Vk2 = Vnz+3

          Wi1 = Wnx-2
          Wi2 = Wnx
          Wj1 = -2
          Wj2 = Wny+3
          Wk1 = -2
          Wk2 = Wnz+3

          Pri1 = Prnx-2
          Pri2 = Prnx
          Prj1 = -2
          Prj2 = Prny+3
          Prk1 = -2
          Prk2 = Prnz+3

          domain_bc_send_buffers_copy(1)%Ui1 = Ui1 
          domain_bc_send_buffers_copy(1)%Ui2 = Ui2 
          domain_bc_send_buffers_copy(1)%Uj1 = Uj1 
          domain_bc_send_buffers_copy(1)%Uj2 = Uj2 
          domain_bc_send_buffers_copy(1)%Uk1 = Uk1 
          domain_bc_send_buffers_copy(1)%Uk2 = Uk2 

          domain_bc_send_buffers_copy(1)%Vi1 = Vi1 
          domain_bc_send_buffers_copy(1)%Vi2 = Vi2 
          domain_bc_send_buffers_copy(1)%Vj1 = Vj1 
          domain_bc_send_buffers_copy(1)%Vj2 = Vj2 
          domain_bc_send_buffers_copy(1)%Vk1 = Vk1 
          domain_bc_send_buffers_copy(1)%Vk2 = Vk2 

          domain_bc_send_buffers_copy(1)%Wi1 = Wi1 
          domain_bc_send_buffers_copy(1)%Wi2 = Wi2 
          domain_bc_send_buffers_copy(1)%Wj1 = Wj1 
          domain_bc_send_buffers_copy(1)%Wj2 = Wj2 
          domain_bc_send_buffers_copy(1)%Wk1 = Wk1 
          domain_bc_send_buffers_copy(1)%Wk2 = Wk2 

          domain_bc_send_buffers_copy(1)%Pri1 = Pri1
          domain_bc_send_buffers_copy(1)%Pri2 = Pri2
          domain_bc_send_buffers_copy(1)%Prj1 = Prj1
          domain_bc_send_buffers_copy(1)%Prj2 = Prj2
          domain_bc_send_buffers_copy(1)%Prk1 = Prk1
          domain_bc_send_buffers_copy(1)%Prk2 = Prk2


          allocate(domain_bc_send_buffers_copy(1)%U(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
          allocate(domain_bc_send_buffers_copy(1)%V(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
          allocate(domain_bc_send_buffers_copy(1)%W(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

          allocate(domain_bc_send_buffers_copy(1)%dU_dt(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
          allocate(domain_bc_send_buffers_copy(1)%dV_dt(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
          allocate(domain_bc_send_buffers_copy(1)%dW_dt(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

          if (enable_buoyancy) then
            allocate(domain_bc_send_buffers_copy(1)%Temperature(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
            allocate(domain_bc_send_buffers_copy(1)%dTemperature_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          end if
          if (enable_moisture) then
            allocate(domain_bc_send_buffers_copy(1)%Moisture(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
            allocate(domain_bc_send_buffers_copy(1)%dMoisture_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          end if

          domain_bc_send_buffers_copy(1)%U_mpi_type = U_mpi_subarray_type(Ui1,Ui2,Uj1,Uj2,Uk1,Uk2)
          domain_bc_send_buffers_copy(1)%V_mpi_type = V_mpi_subarray_type(Vi1,Vi2,Vj1,Vj2,Vk1,Vk2)
          domain_bc_send_buffers_copy(1)%W_mpi_type = W_mpi_subarray_type(Wi1,Wi2,Wj1,Wj2,Wk1,Wk2)
          !for arrays defined from -1 to Prn+2
          if (enable_buoyancy.or.enable_moisture) then
            domain_bc_send_buffers_copy(1)%Pr_mpi_type = Pr_mpi_subarray_type(Pri1,Pri2,Prj1,Prj2,Prk1,Prk2)
          end if
        endif

      else if (domain_index==2) then

        if (iim==1) then
          Btype(We) = BC_DOMAIN_COPY
          allocate(domain_bc_recv_buffers_copy(We:We))
          domain_bc_recv_buffers_copy(We)%comm = world_comm
          domain_bc_recv_buffers_copy(We)%remote_rank = domain_ranks_grid(1)%arr(domain_nxims(1),jim,kim)        
          domain_bc_recv_buffers_copy(We)%remote_domain = 1
          domain_bc_recv_buffers_copy(We)%direction = We
          domain_bc_recv_buffers_copy(We)%enabled = .true.

          Ui1 = -2
          Ui2 = 0
          Uj1 = -2
          Uj2 = Uny+3
          Uk1 = -2
          Uk2 = Unz+3

          Vi1 = -2
          Vi2 = 0
          Vj1 = -2
          Vj2 = Vny+3
          Vk1 = -2
          Vk2 = Vnz+3

          Wi1 = -2
          Wi2 = 0
          Wj1 = -2
          Wj2 = Wny+3
          Wk1 = -2
          Wk2 = Wnz+3

          Pri1 = -2
          Pri2 = 0
          Prj1 = -2
          Prj2 = Prny+3
          Prk1 = -2
          Prk2 = Prnz+3

          domain_bc_recv_buffers_copy(We)%Ui1 = Ui1 
          domain_bc_recv_buffers_copy(We)%Ui2 = Ui2 
          domain_bc_recv_buffers_copy(We)%Uj1 = Uj1 
          domain_bc_recv_buffers_copy(We)%Uj2 = Uj2 
          domain_bc_recv_buffers_copy(We)%Uk1 = Uk1 
          domain_bc_recv_buffers_copy(We)%Uk2 = Uk2 

          domain_bc_recv_buffers_copy(We)%Vi1 = Vi1 
          domain_bc_recv_buffers_copy(We)%Vi2 = Vi2 
          domain_bc_recv_buffers_copy(We)%Vj1 = Vj1 
          domain_bc_recv_buffers_copy(We)%Vj2 = Vj2 
          domain_bc_recv_buffers_copy(We)%Vk1 = Vk1 
          domain_bc_recv_buffers_copy(We)%Vk2 = Vk2 

          domain_bc_recv_buffers_copy(We)%Wi1 = Wi1 
          domain_bc_recv_buffers_copy(We)%Wi2 = Wi2 
          domain_bc_recv_buffers_copy(We)%Wj1 = Wj1 
          domain_bc_recv_buffers_copy(We)%Wj2 = Wj2 
          domain_bc_recv_buffers_copy(We)%Wk1 = Wk1 
          domain_bc_recv_buffers_copy(We)%Wk2 = Wk2 

          domain_bc_recv_buffers_copy(We)%Pri1 = Pri1
          domain_bc_recv_buffers_copy(We)%Pri2 = Pri2
          domain_bc_recv_buffers_copy(We)%Prj1 = Prj1
          domain_bc_recv_buffers_copy(We)%Prj2 = Prj2
          domain_bc_recv_buffers_copy(We)%Prk1 = Prk1
          domain_bc_recv_buffers_copy(We)%Prk2 = Prk2


          allocate(domain_bc_recv_buffers_copy(We)%U(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
          allocate(domain_bc_recv_buffers_copy(We)%V(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
          allocate(domain_bc_recv_buffers_copy(We)%W(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

          allocate(domain_bc_recv_buffers_copy(We)%dU_dt(Ui1:Ui2,Uj1:Uj2,Uk1:Uk2))
          allocate(domain_bc_recv_buffers_copy(We)%dV_dt(Vi1:Vi2,Vj1:Vj2,Vk1:Vk2))
          allocate(domain_bc_recv_buffers_copy(We)%dW_dt(Wi1:Wi2,Wj1:Wj2,Wk1:Wk2))

          if (enable_buoyancy) then
            allocate(domain_bc_recv_buffers_copy(We)%Temperature(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
            allocate(domain_bc_recv_buffers_copy(We)%dTemperature_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          end if
          if (enable_moisture) then
            allocate(domain_bc_recv_buffers_copy(We)%Moisture(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
            allocate(domain_bc_recv_buffers_copy(We)%dMoisture_dt(Pri1:Pri2,Prj1:Prj2,Prk1:Prk2))
          end if

          domain_bc_recv_buffers_copy(We)%U_mpi_type = U_mpi_subarray_type(Ui1,Ui2,Uj1,Uj2,Uk1,Uk2)
          domain_bc_recv_buffers_copy(We)%V_mpi_type = V_mpi_subarray_type(Vi1,Vi2,Vj1,Vj2,Vk1,Vk2)
          domain_bc_recv_buffers_copy(We)%W_mpi_type = W_mpi_subarray_type(Wi1,Wi2,Wj1,Wj2,Wk1,Wk2)
          !for arrays defined from -1 to Prn+2
          if (enable_buoyancy.or.enable_moisture) then
            domain_bc_recv_buffers_copy(We)%Pr_mpi_type = Pr_mpi_subarray_type(Pri1,Pri2,Prj1,Prj2,Prk1,Prk2)
          end if
        endif

      end if

    end if !enable_multiple_domains

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

      if (allocated(domain_bc_send_buffers_copy)) then
        do i = lbound(domain_bc_send_buffers_copy,1), &
               ubound(domain_bc_send_buffers_copy,1)
          associate(b => domain_bc_send_buffers_copy(i))

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
      
      call send(b%U, b%U_mpi_type, b%remote_rank, b%comm, b%remote_domain*10 + 1)
      call send(b%V, b%V_mpi_type, b%remote_rank, b%comm, b%remote_domain*10 + 2)
      call send(b%W, b%W_mpi_type, b%remote_rank, b%comm, b%remote_domain*10 + 3)

      call send(b%dU_dt, b%U_mpi_type, b%remote_rank, b%comm, b%remote_domain*10 + 1)
      call send(b%dV_dt, b%V_mpi_type, b%remote_rank, b%comm, b%remote_domain*10 + 2)
      call send(b%dW_dt, b%W_mpi_type, b%remote_rank, b%comm, b%remote_domain*10 + 3)

      if (enable_buoyancy) then
        call send(b%Temperature, b%Pr_mpi_type, b%remote_rank, b%comm, b%remote_domain*10 + 4)
        call send(b%dTemperature_dt, b%Pr_mpi_type, b%remote_rank, b%comm, b%remote_domain*10 + 4)
      end if

      if (enable_moisture) then
        call send(b%Moisture, b%Pr_mpi_type, b%remote_rank, b%comm, b%remote_domain*10 + 5)
        call send(b%dMoisture_dt, b%Pr_mpi_type, b%remote_rank, b%comm, b%remote_domain*10 + 5)
      end if

    end subroutine

    subroutine recv_arrays(b)
      type(dom_bc_buffer_copy), intent(inout) :: b
      real(knd) :: avg
    
      call recv(b%U, b%U_mpi_type, b%remote_rank, b%comm, domain_index*10 + 1)   
      call recv(b%V, b%V_mpi_type, b%remote_rank, b%comm, domain_index*10 + 2)     
      call recv(b%W, b%W_mpi_type, b%remote_rank, b%comm, domain_index*10 + 3)    

      call recv(b%dU_dt, b%U_mpi_type, b%remote_rank, b%comm, domain_index*10 + 1)   
      call recv(b%dV_dt, b%V_mpi_type, b%remote_rank, b%comm, domain_index*10 + 2)     
      call recv(b%dW_dt, b%W_mpi_type, b%remote_rank, b%comm, domain_index*10 + 3)

      if (b%rescale_dU_dt) then
        !rescales the time derivative to be zero on average
        !makes sense only in certain situations
        select case (b%direction)
          case(We)
            avg = sum(b%dU_dt(0, 1:Uny, 1:Unz)) / (Uny*Unz)
            b%dU_dt = b%dU_dt - avg
          case(Ea)
            avg = sum(b%dU_dt(Unz+1, 1:Uny, 1:Unz)) / (Uny*Unz)
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
        call recv(b%Temperature, b%Pr_mpi_type, b%remote_rank, b%comm, domain_index*10 + 4)
        call recv(b%dTemperature_dt, b%Pr_mpi_type, b%remote_rank, b%comm, domain_index*10 + 4)
      end if

      if (enable_moisture) then
        call recv(b%Moisture, b%Pr_mpi_type, b%remote_rank, b%comm, domain_index*10 + 5)
        call recv(b%dMoisture_dt, b%Pr_mpi_type, b%remote_rank, b%comm, domain_index*10 + 5)
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



  subroutine par_update_domain_bounds_U(U, eff_time, component)
    !effective time, because it can also reflect individual RK stages
    real(knd), dimension(-2:,-2:,-2:), contiguous, intent(inout) :: U
    real(knd), intent(in) :: eff_time
    integer, intent(in) :: component
    integer :: i
    real(knd) :: S

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do i = lbound(domain_bc_recv_buffers_copy,1), &
               ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(i))

            select case (component)
              case (1)
                U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2) = b%U
                 U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2) = U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2) + &
                                                         b%dU_dt * (eff_time - b%time)

                if (b%direction==We .and. b%rescale_compatibility) then
                    S = sum(U(1,1:Prny,1:Prnz))
                    S = S / sum(b%U(0,1:Prny,1:Prnz))
                    U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2) = S * U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2)  
                end if
                  
              case (2)
                U(b%Vi1:b%Vi2,b%Vj1:b%Vj2,b%Vk1:b%Vk2) = b%V
                U(b%Vi1:b%Vi2,b%Vj1:b%Vj2,b%Vk1:b%Vk2) = U(b%Vi1:b%Vi2,b%Vj1:b%Vj2,b%Vk1:b%Vk2) + &
                                                         b%dV_dt * (eff_time - b%time)
              case (3)
                U(b%Wi1:b%Wi2,b%Wj1:b%Wj2,b%Wk1:b%Wk2) = b%W
                U(b%Wi1:b%Wi2,b%Wj1:b%Wj2,b%Wk1:b%Wk2) = U(b%Wi1:b%Wi2,b%Wj1:b%Wj2,b%Wk1:b%Wk2) + &
                                                         b%dW_dt * (eff_time - b%time)
            end select

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
    integer :: i

    if (enable_multiple_domains) then

      if (allocated(domain_bc_recv_buffers_copy)) then
        do i = lbound(domain_bc_recv_buffers_copy,1), &
               ubound(domain_bc_recv_buffers_copy,1)
          associate(b => domain_bc_recv_buffers_copy(i))

              U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2) = b%U
              U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2) = U(b%Ui1:b%Ui2,b%Uj1:b%Uj2,b%Uk1:b%Uk2) + &
                                                       b%dU_dt * (eff_time - b%time)

              V(b%Vi1:b%Vi2,b%Vj1:b%Vj2,b%Vk1:b%Vk2) = b%V
              V(b%Vi1:b%Vi2,b%Vj1:b%Vj2,b%Vk1:b%Vk2) = V(b%Vi1:b%Vi2,b%Vj1:b%Vj2,b%Vk1:b%Vk2) + &
                                                       b%dV_dt * (eff_time - b%time)

              W(b%Wi1:b%Wi2,b%Wj1:b%Wj2,b%Wk1:b%Wk2) = b%W
              W(b%Wi1:b%Wi2,b%Wj1:b%Wj2,b%Wk1:b%Wk2) = W(b%Wi1:b%Wi2,b%Wj1:b%Wj2,b%Wk1:b%Wk2) + &
                                                       b%dW_dt * (eff_time - b%time)


              if (enable_buoyancy) then
                Temperature(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2) = b%Temperature
              end if
              if (enable_moisture) then
                Moisture(b%Pri1:b%Pri2,b%Prj1:b%Prj2,b%Prk1:b%Prk2) = b%Moisture
              end if

              !TODO: the derivatives
          end associate
        end do
      end if

    end if

  end subroutine

end module domains_bc_par
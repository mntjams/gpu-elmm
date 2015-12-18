module domains_bc_par

  use custom_par


  type dom_bc_buffer_copy
    !the most simple type, just a copy, no interpolation
    !the grid resolution must be identical
    !the boundary position must always exactly coincide with the grid cell boundaries
    real(knd), allocatable, dimension(:,:,:) :: U, V, W, Pr, Temperature, Moisture
    real(knd), allocatable, dimension(:,:,:,:) :: Scalar
    real(knd), allocatable, dimension(:,:,:) :: dU_dt, dV_dt, dW_dt, dPr_dt, &
                                                dTemperature_dt, dMoisture_dt
    real(knd), allocatable, dimension(:,:,:,:) :: dScalar_dt
    integer :: comm !communicator connecting
    integer :: remote_rank
  end type

  type(dom_bc_buffer_copy), allocatable :: domain_bc_buffers_copy(:)
contains



  subroutine par_init_domain_boundary_conditions

    !temporary testing hack

    if (domain_index==1) then
      if (iim==nxims) then
        Btype(Ea) = BC_DOMAIN_COPY
        allocate(domain_bc_buffers_copy(Ea:Ea)
        domain_bc_buffers_copy(Ea)%comm = world_comm
        domain_bc_buffers_copy(Ea)%remote = domain_ranks_grid(2)%arr(1,jim,kim)
        
      endif
    else if (domain_ranks_grid==2) then
      if (iim==1) then
        Btype(We) = BC_DOMAIN_COPY
        allocate(domain_bc_buffers_copy(We:We)
        domain_bc_buffers_copy(We)%comm = world_comm
        domain_bc_buffers_copy(We)%remote = domain_ranks_grid(1)%arr(domain_nxims(1),jim,kim)        
      endif

    end if


  end subroutine



end module domains_bc_par
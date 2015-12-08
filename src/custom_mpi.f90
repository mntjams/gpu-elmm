module custom_mpi
  use mpi
  use Kinds

  integer, parameter :: MPI_real32 = MPI_REAL4, MPI_real64 = MPI_REAL8

#ifdef DPREC
  integer, parameter :: MPI_knd = MPI_real64
#else
  integer, parameter :: MPI_knd = MPI_real32
#endif

  real(real32), pointer, contiguous :: MPI_IN_PLACE_real32(:)
  real(real64), pointer, contiguous :: MPI_IN_PLACE_real64(:)

end module

module custom_par
  use Kinds
  use custom_mpi
  use stop_procedures
  
  implicit none

  private :: MPI_knd, MPI_real32, MPI_real64, MPI_IN_PLACE_real32, MPI_IN_PLACE_real64  

  logical :: enable_multiple_domains = .true.
  
  !number of the domain
  integer :: domain_index = 1
  
  !number of the domain
  integer :: number_of_domains = 1
  
  !numbers of images in individual domains
  integer, allocatable :: domain_nims(:), &
                          domain_nxims(:), &
                          domain_nyims(:), &
                          domain_nzims(:), &
                          domain_comms(:), &
                          domain_groups(:), &
                          domain_comms_with_nth(:), &
                          domain_groups_with_nth(:)
  
  !at which number starts numbering of this domain?
  integer :: first_domain_rank_in_world = 0
  integer :: first_domain_im_in_world = 1
  
  integer :: world_comm_size
  
  
  integer :: nims, npx = -1, npy = -1, npz = -1
  integer :: npxyz(3) = -1, pxyz(3)
  integer :: nxims, nyims, nzims
  integer :: myim, myrank, iim, jim, kim
  integer :: w_im, e_im, s_im, n_im, b_im, t_im
  integer :: w_rank, e_rank, s_rank, n_rank, b_rank, t_rank
  integer :: neigh_ims(6), neigh_ranks(6)
  
  
  integer, parameter :: world_comm = MPI_COMM_WORLD
  integer :: domain_comm = MPI_COMM_WORLD
  integer :: poisfft_comm = MPI_COMM_NULL
  integer :: cart_comm = MPI_COMM_NULL

  integer :: world_group = MPI_GROUP_NULL
  
  !MPI communicators which include the inner or the outer domain
  integer :: inner_comm = MPI_COMM_NULL, outer_comm = MPI_COMM_NULL
  
  integer :: comm_plane_yz = MPI_COMM_NULL, comm_plane_xz = MPI_COMM_NULL, comm_plane_xy = MPI_COMM_NULL
  integer :: comm_row_x = MPI_COMM_NULL, comm_row_y = MPI_COMM_NULL, comm_row_z = MPI_COMM_NULL
  integer :: cart_comm_dim = -1
  integer, allocatable :: images_grid(:,:,:), ranks_grid(:,:,:)
  
  interface par_co_reduce
    module procedure par_co_reduce_32
    module procedure par_co_reduce_64
    module procedure par_co_reduce_32_1d
    module procedure par_co_reduce_64_1d
    module procedure par_co_reduce_logical
    module procedure par_co_reduce_int
  end interface
  
  interface par_co_min
    module procedure par_co_min_32
    module procedure par_co_min_64
    module procedure par_co_min_int
  end interface
  
  interface par_co_max
    module procedure par_co_max_32
    module procedure par_co_max_64
  end interface
  
  interface par_co_sum
    module procedure par_co_sum_32
    module procedure par_co_sum_32_comm
    module procedure par_co_sum_32_comm_1d
    module procedure par_co_sum_64
    module procedure par_co_sum_64_comm
    module procedure par_co_sum_64_comm_1d
    module procedure par_co_sum_int
    module procedure par_co_sum_int_comm
  end interface

  interface par_sum_to_master_horizontal
    module procedure par_sum_to_master_horizontal_32_1d
    module procedure par_sum_to_master_horizontal_64_1d
  end interface

  interface par_broadcast_from_top
    module procedure par_broadcast_from_top_real32
    module procedure par_broadcast_from_top_real64
  end interface

contains

  subroutine par_sync_all()
    integer ie
    call MPI_Barrier(domain_comm, ie)
    if (ie/=0) call error_stop("Error in MPI Barrier.")
  end subroutine

  subroutine par_sync_out(str)
    character(*), intent(in) :: str
    call par_sync_all()
    if (master) write(*,*) str
    call par_sync_all()
  end subroutine

 
  integer function par_this_image() result(res)
    integer ie
    call MPI_Comm_rank(domain_comm, res, ie)
    res = res + 1
    if (ie/=0) call error_stop("MPI_Comm_rank ERROR")
  end function
  

  integer function par_num_images() result(res)
    integer ie
    call MPI_Comm_size(domain_comm, res, ie)  
    if (ie/=0) call error_stop("MPI_Comm_size ERROR")
  end function

  
  integer function par_image_index(sub) result(res)
    integer, intent(in) :: sub(3)
    integer ie
    
    if (cart_comm_dim==-1) then
      call MPI_Cartdim_get(cart_comm, cart_comm_dim, ie)
      if (ie/=0) call error_stop("MPI_Cartdim_get")
    end if
    call MPI_Cart_rank(cart_comm, sub(3:1:-1)-1, res, ie)  
    if (ie/=0) call error_stop("MPI_Cart_rank")
    res = res + 1
  end function
  

  subroutine par_get_image_coords()
    integer :: ie
    integer :: i, j, k
    
    call MPI_Cart_coords(cart_comm, myrank, 3, pxyz, ie)
    if (ie/=0) call error_stop("MPI_Cart_coords")
           
    pxyz = pxyz(3:1:-1)
           
    iim = pxyz(1) + 1
    jim = pxyz(2) + 1
    kim = pxyz(3) + 1
    
    write(*,*) myim, "coords:",iim,jim,kim
    
    if (iim>1) then
      w_im = par_image_index([iim-1, jim, kim])
    else
      w_im = par_image_index([nxims, jim, kim])
    end if
    if (iim<nxims) then
      e_im = par_image_index([iim+1, jim, kim])
    else
      e_im = par_image_index([1, jim, kim])
    end if
    
    if (jim>1) then
      s_im = par_image_index([iim, jim-1, kim])
    else
      s_im = par_image_index([iim, nyims, kim])
    end if
    if (jim<nyims) then
      n_im = par_image_index([iim, jim+1, kim])
    else
      n_im = par_image_index([iim, 1, kim])
    end if

    if (kim>1) then
      b_im = par_image_index([iim, jim, kim-1])
    else
      b_im = par_image_index([iim, jim, nzims])
    end if
    if (kim<nzims) then
      t_im = par_image_index([iim, jim, kim+1])
    else
      t_im = par_image_index([iim, jim, 1])
    end if
    
    w_rank = w_im - 1
    e_rank = e_im - 1
    s_rank = s_im - 1
    n_rank = n_im - 1
    b_rank = b_im - 1
    t_rank = t_im - 1
    
    neigh_ranks = [w_rank, e_rank, s_rank, n_rank, b_rank, t_rank]
    
    neigh_ims = [w_im, e_im, s_im, n_im, b_im, t_im]
    
    allocate(images_grid(1:nxims, 1:nyims, 1:nzims))
    do k = 1, nzims
      do j = 1, nyims
        do i = 1, nxims
          images_grid(i,j,k) = par_image_index([i, j, k])
        end do
      end do
    end do
    
    ranks_grid = images_grid - 1
     
  end subroutine
  
  
  
  subroutine par_init
    use Kinds
    use iso_c_binding

    integer :: ie
    integer :: required, provided

    required = MPI_THREAD_SERIALIZED
  
    call MPI_Init_thread(required, provided, ie)
    if (ie/=0) call error_stop("Error in MPI_Init")

    if (provided<required) then
      write(*,*) "------------------------------"
      write(*,*) "Error, the provided MPI threading support smaller than required!"
      write(*,*) "required:", required
      write(*,*) "provided:", provided
      write(*,*) "Trying to continue anyway, but a crash is likely and the results will be questionable."
      write(*,*) "------------------------------"

      call par_sync_all()
    end if

    call c_f_pointer(my_loc(MPI_IN_PLACE), MPI_IN_PLACE_real32, [1])
    call c_f_pointer(my_loc(MPI_IN_PLACE), MPI_IN_PLACE_real64, [1])
    
    call MPI_Errhandler_set(world_comm, MPI_ERRORS_RETURN, ie)
    if (ie/=0) call error_stop("Error in MPI_Errhandler_set")
    
    
    call MPI_Comm_size(world_comm, world_comm_size, ie)
    if (ie/=0) call error_stop("Error calling MPI_Comm_size")
    
    
  contains

    type(c_ptr) function my_loc(t)
      integer, intent(in), target :: t
      my_loc = c_loc(t)
    end function
      
  end subroutine par_init
  
  
  
  subroutine par_init_domains
    integer :: ie
    integer, allocatable :: check_n(:)
    integer :: dom, i
    
    interface
      subroutine MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        integer :: SENDBUF, RECVBUF(*)
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
      end subroutine
      subroutine MPI_ALLGATHER(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, &
                 RECVTYPE, COMM, IERROR)
        integer ::  SENDBUF, RECVBUF(*)
        integer :: SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE, COMM
        integer :: IERROR
      end subroutine
    end interface
    
    if (number_of_domains < domain_index) then
      write(*,*) "Error, domain index", domain_index, " larger than the number of domains", number_of_domains
      call error_stop()
    end if
    
    
    allocate(check_n(world_comm_size))
    
    call MPI_AllGather(number_of_domains, 1, MPI_INTEGER, &
                       check_n, 1, MPI_INTEGER, &
                       world_comm, ie)
                       
    if (any(check_n /= number_of_domains)) then
      write(*,*) "Error, number of domains is not defined equally for all images."
      call error_stop()
    end if

  
    allocate(domain_nims(number_of_domains))
    allocate(domain_nxims(number_of_domains))
    allocate(domain_nyims(number_of_domains))
    allocate(domain_nzims(number_of_domains))
    
    allocate(domain_comms(number_of_domains))
    allocate(domain_groups(number_of_domains))
    allocate(domain_comms_with_nth(number_of_domains))
    allocate(domain_groups_with_nth(number_of_domains))
    
    domain_nims = 0
    
    domain_nims(domain_index) = 1
    
    call MPI_AllReduce(MPI_IN_PLACE, domain_nims, number_of_domains, &
                       MPI_INTEGER, MPI_SUM, world_comm, ie)
                       
    if (domain_nims(domain_index) /= product(npxyz)) then
      write(*,*) "Error, npxyz must be specified and equal to the number of MPI processes for each domain."
      write(*,*) domain_nims(domain_index), " /= ", npxyz(1) * npxyz(2) * npxyz(3), " for domain ", domain_index
      call error_stop()
    end if
    
    if (any(domain_nims <= 0)) then
      dom = minloc(domain_nims, 1)
      write(*,*) "Error, ", domain_nims(dom), " images in domain", dom
      call error_stop()
    end if

    call MPI_Comm_group(world_comm, world_group, ie)
    if (ie/=0) call error_stop("Error calling MPI_Comm_group.")


    first_domain_rank_in_world = sum(domain_nims(1:domain_index-1))
    first_domain_im_in_world = first_domain_rank_in_world + 1

    do dom = 1, number_of_domains
      call MPI_Group_incl(world_group, domain_nims(dom), &
             [( sum(domain_nims(1:dom-1)) + i - 1, i = 1,  domain_nims(dom) )], &
             domain_groups(dom), &
             ie)
      if (ie/=0) call error_stop("Error calling MPI_Group_incl.")

      call MPI_Comm_create(world_comm, domain_groups(dom), domain_comms(dom), ie)
      if (ie/=0) call error_stop("Error calling MPI_Comm_create.")
    end do

    domain_comm = domain_comms(domain_index)  
    
    do dom = 1, number_of_domains
      call MPI_Group_union(domain_groups(domain_index), domain_groups(dom), domain_groups_with_nth(dom), ie)
      if (ie/=0) call error_stop("Error calling MPI_Group_union.")

      call MPI_Comm_create(world_comm, domain_groups_with_nth(dom), domain_comms_with_nth(dom), ie)
      if (ie/=0) call error_stop("Error calling MPI_Comm_create.")
    end do

  end subroutine par_init_domains
  
  
  
  subroutine par_init_grid
    use PoisFFT
    use iso_c_binding
    use Parameters, only: gPrnx, gPrny, gPrnz, gPrns, Prnx, Prny, Prnz, &
                          offset_to_global_x, offset_to_global_y, offset_to_global_z, &
                          offsets_to_global, output_dir, image_input_dir
    use strings, only: itoa
    
    integer(c_intptr_t) :: ng(3), nxyz(3), off(3), nxyz2(3), nsxyz2(3)
    integer :: ie
    integer :: pos
    character(80) :: str_dir

    if (enable_multiple_domains) then
      call par_init_domains
    else
      domain_comm = world_comm
    end if

    nims = par_num_images()

    myim = par_this_image()
    
    myrank = myim - 1

    master = (myrank==0)

!     npxyz = [npx, npy, npz]
    
    if (any(npxyz<1)) then
      npxyz(1) = 1
      npxyz(2) = nint(sqrt(real(nims)))
      npxyz(3) = nims / npxyz(2)
      if (master)    write (*,*) "Trying to decompose in",npxyz,"process grid."
    end if

    nxims = npxyz(1)
    nyims = npxyz(2)
    nzims = npxyz(3)
    
    if (product(npxyz)/=nims) then
      if (master) then
        write(*,*) "Could not decompose the processes to N x N grid."
        write(*,*) "Try a perfect square for number of processes."
      end if
      call error_stop(25)
    end if

    !The number of images does not have to be 1 necesarilly, but it is 
    ! for the optimal preformance of the FFT in Poisson solver PoisFFT.
    if (npxyz(1)/=1) call error_stop("The number of images in x direction must be one.")
    
    !Create the 3D Cartesian communicator "cart_comm" for use in CLMM
    call MPI_Cart_create(domain_comm, 3, int(npxyz(3:1:-1)), &
                         [.false.,.false.,.false.], &
                         .true., cart_comm, ie)
    if (ie/=0) call error_stop("Error calling MPI_Cart_create.")
    
    domain_comm = cart_comm
    
    myim = par_this_image()
    
    myrank = myim - 1
    !very unlikely to be changed
    master = (myrank==0)
    
    !creates a 2D! Cartesian communicator "poisfft_comm"
    !2D because of the PFFT library
    !the dimensions in the poisfft_comm are z,y
    call PoisFFT_InitMPIGrid(domain_comm, npxyz(3:2:-1), poisfft_comm, ie)
    
    call par_init_sub_comms
    
    call par_get_image_coords
    
    gPrnx = Prnx
    gPrny = Prny
    gPrnz = Prnz
    
    gPrns = [gPrnx, gPrny, gPrnz]
    
    ng = gPrns

    call PoisFFT_LocalGridSize(3,ng,cart_comm,nxyz,off,nxyz2,nsxyz2)
    if (any(nxyz/=nxyz2).or.any(off/=nsxyz2)) call error_stop(40)
    
    call MPI_Barrier(domain_comm, ie)
    if (ie/=0) call error_stop("Error in MPI Barrier.")
    
    write(*,*) iim,jim,kim, "nxyz:", nxyz
    
    call MPI_Barrier(domain_comm, ie)
    if (ie/=0) call error_stop("Error in MPI Barrier.")
    
    if (par_co_any(any(nxyz<=0))) then
      call error_stop("Zero or negative grid size on one or more images. Try different process grid.")
    end if
    
    Prnx = int(nxyz(1))
    Prny = int(nxyz(2))
    Prnz = int(nxyz(3))
    
    offset_to_global_x = int(off(1))
    offset_to_global_y = int(off(2))
    offset_to_global_z = int(off(3))
    
    offsets_to_global = int(off)
    
    output_dir = output_dir // "im-" // itoa(iim) // &
                                 "-" // itoa(jim) // &
                                 "-" // itoa(kim) // "/"

    image_input_dir = image_input_dir // "im-" // itoa(iim) // &
                                           "-" // itoa(jim) // &
                                           "-" // itoa(kim) // "/"

  end subroutine par_init_grid

  
  subroutine par_init_boundaries
    use Parameters
    
    call helper(Btype)
    call helper(TempBtype)
    call helper(MoistBtype)
    call helper(ScalBtype)
    
  contains
  
    subroutine helper(Bt)
      integer, intent(inout) :: Bt(:)
      
      if (nxims>1) then
        if (iim>1) then
          Bt(We) = MPI_BOUNDARY
        end if
        if (iim<nxims) then
          Bt(Ea) = MPI_BOUNDARY
        end if
        if (iim==1.and.Bt(We)==PERIODIC) Bt(We) = MPI_PERIODIC
        if (iim==nxims.and.Bt(Ea)==PERIODIC) Bt(Ea) = MPI_PERIODIC
      end if
    
      if (nyims>1) then
        if (jim>1) then
          Bt(So) = MPI_BOUNDARY
        end if
        if (jim<nyims) then
          Bt(No) = MPI_BOUNDARY
        end if
        if (jim==1.and.Bt(So)==PERIODIC) Bt(So) = MPI_PERIODIC
        if (jim==nyims.and.Bt(No)==PERIODIC) Bt(No) = MPI_PERIODIC
      end if
    
      if (nzims>1) then
        if (kim>1) then
          Bt(Bo) = MPI_BOUNDARY
        end if
        if (kim<nzims) then
          Bt(To) = MPI_BOUNDARY
        end if
        if (kim==1.and.Bt(Bo)==PERIODIC) Bt(Bo) = MPI_PERIODIC
        if (kim==nzims.and.Bt(To)==PERIODIC) Bt(To) = MPI_PERIODIC
      end if
    end subroutine
    
  end subroutine par_init_boundaries
  
  subroutine par_finalize
    integer :: ie
    call MPI_Barrier(domain_comm, ie)
    if (ie/=0) call error_stop("Error when waiting before finalizing MPI.")
    call MPI_Finalize(ie)
    if (ie/=0) call error_stop("Error finalizing MPI.")
  end subroutine
  
  
  
  subroutine par_init_sub_comms
    integer :: ie
  
    !note the order of the dimensions is reversed
    call MPI_Cart_sub(cart_comm, [.true., .true., .false.], comm_plane_yz, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")
  
    call MPI_Cart_sub(cart_comm, [.true., .false., .true.], comm_plane_xz, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")
  
    call MPI_Cart_sub(cart_comm, [.false., .true., .true.], comm_plane_xy, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")
  
  
    call MPI_Cart_sub(cart_comm, [.false., .false., .true.], comm_row_x, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")
  
    call MPI_Cart_sub(cart_comm, [.false., .true., .false.], comm_row_y, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")
  
    call MPI_Cart_sub(cart_comm, [.true., .false., .false.], comm_row_z, ie)
    if (ie/=0) call error_stop("Error creating comm_plane_x.")

  end subroutine
  
  
  function par_co_reduce_32(x,op,comm) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer, intent(in) :: op, comm
    integer ie
    
    interface
      subroutine MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        real(real32) :: SENDBUF, RECVBUF
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
      end subroutine
    end interface
    
    call MPI_AllReduce(x, res, &
                       count=1, datatype=MPI_KND, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_32.")
  end function

  function par_co_reduce_64(x,op,comm) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer, intent(in) :: op, comm
    integer ie
    
    interface
      subroutine MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        real(real64) :: SENDBUF, RECVBUF
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
      end subroutine
    end interface
    
    call MPI_AllReduce(x, res, &
                       count=1, datatype=MPI_KND, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_64.")
  end function

  function par_co_reduce_32_1d(x,op,comm) result(res)
    real(real32),intent(in) :: x(:)
    real(real32) :: res(size(x))
    integer, intent(in) :: op, comm
    integer ie
    
    interface
      subroutine MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
        real(real32) :: SENDBUF(count), RECVBUF(count)
      end subroutine
    end interface
    
    call MPI_AllReduce(x, res, &
                       count=size(x), datatype=MPI_KND, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_32.")
  end function

  function par_co_reduce_64_1d(x,op,comm) result(res)
    real(real64),intent(in) :: x(:)
    real(real64) :: res(size(x))
    integer, intent(in) :: op, comm
    integer ie
    
    interface
      subroutine MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
        real(real64) :: SENDBUF(count), RECVBUF(count)
      end subroutine
    end interface
    
    call MPI_AllReduce(x, res, &
                       count=size(x), datatype=MPI_KND, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_64.")
  end function

  function par_co_reduce_logical(x,op,comm) result(res)
    logical :: res
    logical,intent(in) :: x
    integer, intent(in) :: op, comm
    integer ie
    
    interface
      subroutine MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        logical :: SENDBUF, RECVBUF
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
      end subroutine
    end interface
    
    call MPI_AllReduce(x, res, &
                       count=1, datatype=MPI_LOGICAL, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_logical.")
  end function

  function par_co_reduce_int(x,op,comm) result(res)
    integer :: res
    integer,intent(in) :: x
    integer, intent(in) :: op, comm
    integer ie
    
    interface
      subroutine MPI_ALLREDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, COMM, IERROR)
        import
        integer :: SENDBUF, RECVBUF
        integer :: COUNT, DATATYPE, OP, COMM, IERROR
      end subroutine
    end interface
    
    call MPI_AllReduce(x, res, &
                       count=1, datatype=MPI_INTEGER, op=op, &
                       comm=comm, ierror=ie)
    if (ie/=0) call error_stop("Error in par_co_reduce_logical.")
  end function

  function par_co_min_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer ie
    
    res = par_co_reduce(x, MPI_MIN, domain_comm)
  end function

  function par_co_min_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer ie
    
    res = par_co_reduce(x, MPI_MIN, domain_comm)
  end function

  function par_co_min_int(x) result(res)
    integer :: res
    integer,intent(in) :: x
    integer ie
    
    res = par_co_reduce(x, MPI_MIN, domain_comm)
  end function

  function par_co_max_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer ie
    
    res = par_co_reduce(x, MPI_MAX, domain_comm)
  end function

  function par_co_max_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer ie
    
    res = par_co_reduce(x, MPI_MAX, domain_comm)
  end function

  function par_co_sum_int(x) result(res)
    integer :: res
    integer,intent(in) :: x
    integer ie
    
    res = par_co_reduce(x, MPI_SUM, domain_comm)
  end function

  function par_co_sum_int_comm(x, comm) result(res)
    integer :: res
    integer,intent(in) :: x
    integer, intent(in) :: comm
    integer ie
    
    res = par_co_reduce(x, MPI_SUM, comm)
  end function

  function par_co_sum_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer ie
    
    res = par_co_reduce(x, MPI_SUM, domain_comm)
  end function

  function par_co_sum_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer ie
    
    res = par_co_reduce(x, MPI_SUM, domain_comm)
  end function

  function par_co_sum_32_comm(x, comm) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer, intent(in) :: comm
    integer ie
    
    res = par_co_reduce(x, MPI_SUM, comm)
  end function

  function par_co_sum_64_comm(x, comm) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer, intent(in) :: comm
    integer ie
    
    res = par_co_reduce(x, MPI_SUM, comm)
  end function

  function par_co_sum_32_comm_1d(x, comm) result(res)
    real(real32),intent(in) :: x(:)
    real(real32) :: res(size(x))
    integer, intent(in) :: comm
    integer ie
    
    res = par_co_reduce(x, MPI_SUM, comm)
  end function

  function par_co_sum_64_comm_1d(x, comm) result(res)
    real(real64),intent(in) :: x(:)
    real(real64) :: res(size(x))
    integer, intent(in) :: comm
    integer ie
    
    res = par_co_reduce(x, MPI_SUM, comm)
  end function

  function par_co_any(x) result(res)
    logical :: res
    logical,intent(in) :: x
    integer ie
    
    res = par_co_reduce(x, MPI_LOR, domain_comm)
  end function


  function par_co_all(x) result(res)
    logical :: res
    logical,intent(in) :: x
    integer ie
    
    res = par_co_reduce(x, MPI_LAND, domain_comm)
  end function

  subroutine par_sum_to_master_horizontal_32_1d(x)
    real(real32), intent(inout) :: x(:)
    integer ie
    interface
      subroutine MPI_REDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR)
        import
        INTEGER    COUNT, DATATYPE, OP, ROOT, COMM, IERROR
        real(real32)  SENDBUF(*), RECVBUF(count)
      end subroutine
    end interface
    !sums the values across horizontal plane of images to the master image (iim==1 and jim==1)
    if (iim==1.and.jim==1) then
      call MPI_Reduce(MPI_IN_PLACE_real32, x, size(x), MPI_real32, MPI_SUM, 0, comm_plane_xy, ie)
    else
      call MPI_Reduce(x, x, size(x), MPI_real32, MPI_SUM, 0, comm_plane_xy, ie)
    end if
  end subroutine

  subroutine par_sum_to_master_horizontal_64_1d(x)
    real(real64), intent(inout) :: x(:)
    integer ie
    interface
      subroutine MPI_REDUCE(SENDBUF, RECVBUF, COUNT, DATATYPE, OP, ROOT, COMM, IERROR)
        import
        INTEGER    COUNT, DATATYPE, OP, ROOT, COMM, IERROR
        real(real64)  SENDBUF(*), RECVBUF(count)
      end subroutine
    end interface
    !sums the values across horizontal plane of images to the master image (iim==1 and jim==1)
    if (iim==1.and.jim==1) then
      call MPI_Reduce(MPI_IN_PLACE_real64, x, size(x), MPI_real64, MPI_SUM, 0, comm_plane_xy, ie)
    else
      call MPI_Reduce(x, x, size(x), MPI_real64, MPI_SUM, 0, comm_plane_xy, ie)
    end if
  end subroutine

  subroutine par_broadcast_from_top_real32(x)
    real(real32), intent(inout) :: x
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        real(real32) ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_real32, nzims-1, comm_row_z, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine

  subroutine par_broadcast_from_top_real64(x)
    real(real64), intent(inout) :: x
    integer :: ie
    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        real(real64) ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    call MPI_Bcast(x, 1, MPI_real64, nzims-1, comm_row_z, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)
  end subroutine

end module custom_par

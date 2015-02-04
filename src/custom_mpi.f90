#ifdef MPI
module custom_mpi
  use Kinds
  use mpi
  use stop_procedures
  
  implicit none

  
  integer :: nims, npx = -1, npy = -1, npz = -1
  integer :: npxyz(3) = -1, pxyz(3)
  integer :: nxims, nyims, nzims
  integer :: myim, myrank, iim, jim, kim
  integer :: w_im, e_im, s_im, n_im, b_im, t_im
  integer :: w_rank, e_rank, s_rank, n_rank, b_rank, t_rank
  integer :: global_comm, poisfft_comm, cart_comm
  integer :: comm_plane_yz = -1, comm_plane_xz = -1, comm_plane_xy = -1
  integer :: comm_row_x = -1, comm_row_y = -1, comm_row_z = -1
  integer :: cart_comm_dim = -1
  integer :: MPI_knd = -huge(1)
  integer, allocatable :: images_grid(:,:,:), ranks_grid(:,:,:)
  
  interface mpi_co_reduce
    module procedure mpi_co_reduce_32
    module procedure mpi_co_reduce_64
    module procedure mpi_co_reduce_logical
    module procedure mpi_co_reduce_int
  end interface
  
  interface mpi_co_min
    module procedure mpi_co_min_32
    module procedure mpi_co_min_64
    module procedure mpi_co_min_int
  end interface
  
  interface mpi_co_max
    module procedure mpi_co_max_32
    module procedure mpi_co_max_64
  end interface
  
  interface mpi_co_sum
    module procedure mpi_co_sum_32
    module procedure mpi_co_sum_32_comm
    module procedure mpi_co_sum_64
    module procedure mpi_co_sum_64_comm
  end interface

contains

 
  integer function mpi_this_image() result(res)
    integer ie
    call MPI_Comm_rank(global_comm, res, ie)
    res = res + 1
    if (ie/=0) call error_stop("MPI_Comm_rank ERROR")
  end function
  

  integer function mpi_num_images() result(res)
    integer ie
    call MPI_Comm_size(global_comm, res, ie)  
    if (ie/=0) call error_stop("MPI_Comm_size ERROR")
  end function

  
  integer function mpi_image_index(sub) result(res)
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
  

  subroutine get_image_coords()
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
      w_im = mpi_image_index([iim-1, jim, kim])
    else
      w_im = mpi_image_index([nxims, jim, kim])
    end if
    if (iim<nxims) then
      e_im = mpi_image_index([iim+1, jim, kim])
    else
      e_im = mpi_image_index([1, jim, kim])
    end if
    
    if (jim>1) then
      s_im = mpi_image_index([iim, jim-1, kim])
    else
      s_im = mpi_image_index([iim, nyims, kim])
    end if
    if (jim<nyims) then
      n_im = mpi_image_index([iim, jim+1, kim])
    else
      n_im = mpi_image_index([iim, 1, kim])
    end if

    if (kim>1) then
      b_im = mpi_image_index([iim, jim, kim-1])
    else
      b_im = mpi_image_index([iim, jim, nzims])
    end if
    if (kim<nzims) then
      t_im = mpi_image_index([iim, jim, kim+1])
    else
      t_im = mpi_image_index([iim, jim, 1])
    end if
    
    w_rank = w_im - 1
    e_rank = e_im - 1
    s_rank = s_im - 1
    n_rank = n_im - 1
    b_rank = b_im - 1
    t_rank = t_im - 1
    
    allocate(images_grid(1:nxims, 1:nyims, 1:nzims))
    do k = 1, nzims
      do j = 1, nyims
        do i = 1, nxims
          images_grid(i,j,k) = mpi_image_index([i, j, k])
        end do
      end do
    end do
    
    ranks_grid = images_grid - 1
     
  end subroutine

  
  
  subroutine init_custom_mpi
    use Kinds, only: knd
    integer :: ie
  
    if (knd == kind(1.)) then
      MPI_knd = MPI_REAL
    else if (knd == kind(1.D0)) then
      MPI_knd = MPI_DOUBLE_PRECISION
    end if
    
    global_comm = MPI_COMM_WORLD

    call MPI_Init(ie)
    if (ie/=0) call error_stop("Error in MPI_Init")
    
    call MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL, ie)
    if (ie/=0) call error_stop("Error in MPI_Errhandler_set")
    
    nims = mpi_num_images()

    myim = mpi_this_image()
    
    myrank = myim - 1

    master = (myrank==0)
    
  end subroutine
  
  subroutine init_mpi_grid
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
    call MPI_Cart_create(global_comm, 3, int(npxyz(3:1:-1)), &
                         [.false.,.false.,.false.], &
                         .false., cart_comm, ie)
    if (ie/=0) call error_stop("Error calling MPI_Cart_create.")

    !creates a 2D! Cartesian communicator "poisfft_comm"
    !2D because of the PFFT library
    !the dimnsions in the poisfft_comm are z,y
    call PoisFFT_InitMPIGrid(global_comm, npxyz(3:2:-1), poisfft_comm, ie)
    !This uses the 3D communicator, slower, but necessary for 3D decomposition
!     poisfft_comm = cart_comm
    
    call init_sub_comms
    
    call get_image_coords
    
    gPrnx = Prnx
    gPrny = Prny
    gPrnz = Prnz
    
    gPrns = [gPrnx, gPrny, gPrnz]
    
    ng = gPrns

    call PoisFFT_LocalGridSize(3,ng,cart_comm,nxyz,off,nxyz2,nsxyz2)
    if (any(nxyz/=nxyz2).or.any(off/=nsxyz2)) call error_stop(40)
    
    call MPI_Barrier(global_comm, ie)
    if (ie/=0) call error_stop("Error in MPI Barrier.")
    
    write(*,*) iim,jim,kim, "nxyz:", nxyz
    
    call MPI_Barrier(global_comm, ie)
    if (ie/=0) call error_stop("Error in MPI Barrier.")
    
    if (mpi_co_any(any(nxyz<=0))) then
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

  end subroutine

  
  subroutine init_mpi_boundaries
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
    
  end subroutine
  
  subroutine finalize_custom_mpi
    integer :: ie
    call MPI_Barrier(global_comm, ie)
    if (ie/=0) call error_stop("Error when waiting before finalizing MPI.")
    call MPI_Finalize(ie)
    if (ie/=0) call error_stop("Error finalizing MPI.")
  end subroutine
  
  
  
  subroutine init_sub_comms
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
  
  
  function mpi_co_reduce_32(x,op,comm) result(res)
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
    if (ie/=0) call error_stop("Error in mpi_co_reduce_32.")
  end function

  function mpi_co_reduce_64(x,op,comm) result(res)
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
    if (ie/=0) call error_stop("Error in mpi_co_reduce_64.")
  end function

  function mpi_co_reduce_logical(x,op,comm) result(res)
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
    if (ie/=0) call error_stop("Error in mpi_co_reduce_logical.")
  end function

  function mpi_co_reduce_int(x,op,comm) result(res)
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
    if (ie/=0) call error_stop("Error in mpi_co_reduce_logical.")
  end function

  function mpi_co_min_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer ie
    
    res = mpi_co_reduce(x, MPI_MIN, global_comm)
  end function

  function mpi_co_min_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer ie
    
    res = mpi_co_reduce(x, MPI_MIN, global_comm)
  end function

  function mpi_co_min_int(x) result(res)
    integer :: res
    integer,intent(in) :: x
    integer ie
    
    res = mpi_co_reduce(x, MPI_MIN, global_comm)
  end function

  function mpi_co_max_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer ie
    
    res = mpi_co_reduce(x, MPI_MAX, global_comm)
  end function

  function mpi_co_max_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer ie
    
    res = mpi_co_reduce(x, MPI_MAX, global_comm)
  end function

  function mpi_co_sum_32(x) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer ie
    
    res = mpi_co_reduce(x, MPI_SUM, global_comm)
  end function

  function mpi_co_sum_64(x) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer ie
    
    res = mpi_co_reduce(x, MPI_SUM, global_comm)
  end function

  function mpi_co_sum_32_comm(x, comm) result(res)
    real(real32) :: res
    real(real32),intent(in) :: x
    integer, intent(in) :: comm
    integer ie
    
    res = mpi_co_reduce(x, MPI_SUM, comm)
  end function

  function mpi_co_sum_64_comm(x, comm) result(res)
    real(real64) :: res
    real(real64),intent(in) :: x
    integer, intent(in) :: comm
    integer ie
    
    res = mpi_co_reduce(x, MPI_SUM, comm)
  end function

  function mpi_co_any(x) result(res)
    logical :: res
    logical,intent(in) :: x
    integer ie
    
    res = mpi_co_reduce(x, MPI_LOR, global_comm)
  end function


  function mpi_co_all(x) result(res)
    logical :: res
    logical,intent(in) :: x
    integer ie
    
    res = mpi_co_reduce(x, MPI_LAND, global_comm)
  end function


end module
#endif
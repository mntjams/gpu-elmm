module exchange_par
  use custom_par, only: global_comm, w_rank, e_rank, s_rank, n_rank, b_rank, t_rank, &
                        nxims, nyims, nzims, iim, jim, kim

#ifdef MPI
  use custom_mpi, only: MPI_knd, MPI_STATUS_SIZE
#endif

  use Kinds
  
  implicit none
  
  private
  
  public par_exchange_boundaries, par_exchange_boundaries_yz, par_exchange_Pr, &
         par_exchange_Q, par_exchange_Sc_x, par_exchange_Sc_y, par_exchange_Sc_z, &
         par_exchange_U_x, par_exchange_U_y, par_exchange_U_z

  interface
    subroutine MPI_Recv(BUF, COUNT, DATATYPE, SOURCE, TAG, COMM, STATUS, IERROR)
      import
      real(knd) :: BUF(*)
      integer   :: COUNT, DATATYPE, SOURCE, TAG, COMM
      integer   :: STATUS(MPI_STATUS_SIZE), IERROR
    end subroutine
    subroutine MPI_Send(BUF, COUNT, DATATYPE, DEST, TAG, COMM, IERROR)
      import
      real(knd) :: BUF(*)
      integer   :: COUNT, DATATYPE, DEST, TAG, COMM, IERROR
    end subroutine
  end interface
      
contains

  subroutine par_exchange_boundaries(Phi, nx, ny, nz, Btype, lb, width, dir)
    use Parameters, only: We, Ea, So, No, Bo, To, MPI_PERIODIC
    real(knd), intent(inout),contiguous :: Phi(lb:,lb:,lb:)
    integer, intent(in) :: nx, ny, nz
    integer, intent(in) :: Btype(6)
    integer, intent(in) :: lb, width
    integer, intent(in), optional :: dir
    logical :: oddx, oddy, oddz, evenx, eveny, evenz
    integer :: ierr,tag,status(MPI_STATUS_SIZE)
    logical :: xdir, ydir, zdir
    
    ierr = 0; status=0; tag = 0
    
    if (present(dir)) then
      xdir = .false.
      ydir = .false.
      zdir = .false.

      select case(dir)
        case(1)
          xdir = .true.
        case(2)
          ydir = .true.
        case(3)
          zdir = .true.
      end select
    else
      xdir = .true.
      ydir = .true.
      zdir = .true.
    end if
    
    oddx = mod(iim,2) == 1
    evenx = .not. oddx
    
    oddy = mod(jim,2) == 1
    eveny = .not. oddy
    
    oddz = mod(kim,2) == 1
    evenz = .not. oddz

    !internal boundaries
    if (xdir) then
      tag = tag + 1
      if (oddx) then
        call send_w
      else
        call recv_e
      end if
      if (evenx) then
        call send_w
      else
        call recv_e
      end if

      tag = tag + 1
      if (oddx) then
        call send_e
      else
        call recv_w
      end if
      if (evenx) then
        call send_e
      else
        call recv_w
      end if
    end if

    if (ydir) then
      tag = tag + 1
      if (oddy) then
        call send_s
      else
        call recv_n
      end if
      if (eveny) then
        call send_s
      else
        call recv_n
      end if

      tag = tag + 1
      if (oddy) then
        call send_n
      else
        call recv_s
      end if
      if (eveny) then
        call send_n
      else
        call recv_s
      end if
    end if

    if (zdir) then
      tag = tag + 1
      if (oddz) then
        call send_b
      else
        call recv_t
      end if
      if (evenz) then
        call send_b
      else
        call recv_t
      end if

      tag = tag + 1
      if (oddz) then
        call send_t
      else
        call recv_b
      end if
      if (evenz) then
        call send_t
      else
        call recv_b
      end if
    end if

    !global domain boundaries
    if (xdir) then
      tag = tag + 1
      if (Btype(We)==MPI_PERIODIC.or.Btype(Ea)==MPI_PERIODIC) then
        if (iim==1) then
          call send(Phi(1:0+width,1:ny,1:nz), w_rank)
        else if (iim==nxims) then
          call recv(Phi(nx+1:nx+width,1:ny,1:nz), e_rank)
        end if     
        if (iim==nxims) then
          call send(Phi(nx+1-width:nx,1:ny,1:nz), e_rank)
        else if (iim==1) then
          call recv(Phi(1-width:0,1:ny,1:nz), w_rank)
        end if
      end if
    end if

    if (ydir) then  
      tag = tag + 1
      if (Btype(So)==MPI_PERIODIC.or.Btype(No)==MPI_PERIODIC) then
        if (jim==1) then
          call send(Phi(1-width:nx+width,1:0+width,1:nz), s_rank)
        else if (jim==nyims) then
          call recv(Phi(1-width:nx+width,ny+1:ny+width,1:nz), n_rank)
        end if
        if (jim==nyims) then
          call send(Phi(1-width:nx+width,ny+1-width:ny,1:nz), n_rank)
        else if (jim==1) then
          call recv(Phi(1-width:nx+width,1-width:0,1:nz), s_rank)
        end if
      end if
    end if

    if (zdir) then
      tag = tag + 1
      if (Btype(Bo)==MPI_PERIODIC.or.Btype(To)==MPI_PERIODIC) then
        if (kim==1) then
          call send(Phi(1-width:nx+width,1-width:ny+width,1:0+width), b_rank)
        else if (kim==nzims) then
          call recv(Phi(1-width:nx+width,1-width:ny+width,nz+1:nz+width), t_rank)
        end if
        if (kim==nzims) then
          call send(Phi(1-width:nx+width,1-width:ny+width,nz+1-width:nz), t_rank)
        else if (kim==1) then
          call recv(Phi(1-width:nx+width,1-width:ny+width,1-width:0), b_rank)
        end if
      end if
    end if
    
  contains
  
  
    subroutine send(a,to)
      real(knd), intent(in) :: a(:,:,:)
      integer, intent(in) :: to

      call MPI_Send(a, size(a) , MPI_KND, to, tag, global_comm, ierr)
      if (ierr/=0) stop "error sending MPI message."
    end subroutine
    
    subroutine recv(a,from)
      real(knd), intent(out) :: a(:,:,:)
      integer, intent(in) :: from

      call MPI_Recv(a, size(a) , MPI_KND, from, tag, global_comm, status, ierr)
      if (ierr/=0) stop "error sending MPI message."
    end subroutine
    

    subroutine send_w
      if (iim>1) then
        call send(Phi(1:0+width,1:ny,1:nz), w_rank)
      end if
    end subroutine
    subroutine recv_w
      if (iim>1) then
        call recv(Phi(1-width:0,1:ny,1:nz), w_rank)
      end if
    end subroutine
    subroutine send_e
      if (iim<nxims) then
        call send(Phi(nx+1-width:nx,1:ny,1:nz), e_rank)
      end if
    end subroutine       
    subroutine recv_e
      if (iim<nxims) then
        call recv(Phi(nx+1:nx+width,1:ny,1:nz), e_rank)
      end if
     end subroutine
    subroutine send_s
      if (jim>1) then
        call send(Phi(1-width:nx+width,1:0+width,1:nz), s_rank)
      end if
    end subroutine
    subroutine recv_s
      if (jim>1) then
        call recv(Phi(1-width:nx+width,1-width:0,1:nz), s_rank)
      end if
    end subroutine
    subroutine send_n
      if (jim<nyims) then
        call send(Phi(1-width:nx+width,ny+1-width:ny,1:nz), n_rank)
      end if
    end subroutine
    subroutine recv_n
      if (jim<nyims) then
        call recv(Phi(1-width:nx+width,ny+1:ny+width,1:nz), n_rank)
      end if
    end subroutine
    subroutine send_b
      if (kim>1) then
        call send(Phi(1-width:nx+width,1-width:ny+width,1:0+width), b_rank)
      end if
    end subroutine
    subroutine recv_b
      if (kim>1) then
        call recv(Phi(1-width:nx+width,1-width:ny+width,1-width:0), b_rank)
      end if
    end subroutine
    subroutine send_t
      if (kim<nzims) then
        call send(Phi(1-width:nx+width,1-width:ny+width,nz+1-width:nz), t_rank)
      end if
    end subroutine
    subroutine recv_t
      if (kim<nzims) then
        call recv(Phi(1-width:nx+width,1-width:ny+width,nz+1:nz+width), t_rank)
      end if
    end subroutine
    
  end subroutine par_exchange_boundaries
  
  
  

  subroutine par_exchange_U_x(U, nx, ny, nz)
    use Parameters, only: Btype
    real(knd), intent(inout) :: U(-2:,-2:,-2:)
    integer, intent(in) :: nx, ny, nz
    call par_exchange_boundaries(U, nx, ny, nz, Btype, -2, 3, dir=1)
  end subroutine
  
  subroutine par_exchange_U_y(U, nx, ny, nz)
    use Parameters, only: Btype
    real(knd), intent(inout) :: U(-2:,-2:,-2:)
    integer, intent(in) :: nx, ny, nz
    call par_exchange_boundaries(U, nx, ny, nz, Btype, -2, 3, dir=2)
  end subroutine
  
  subroutine par_exchange_U_z(U, nx, ny, nz)
    use Parameters, only: Btype
    real(knd), intent(inout) :: U(-2:,-2:,-2:)
    integer, intent(in) :: nx, ny, nz
    call par_exchange_boundaries(U, nx, ny, nz, Btype, -2, 3, dir=3)
  end subroutine
  
  subroutine par_exchange_Sc_x(U, SBtype)
    use Parameters, only: Prnx, Prny, Prnz
    real(knd), intent(inout) :: U(-1:,-1:,-1:)
    integer, intent(in) :: SBtype(6)
    call par_exchange_boundaries(U, Prnx, Prny, Prnz, SBtype, -1, 2, dir=1)
  end subroutine
  
  subroutine par_exchange_Sc_y(U, SBType)
    use Parameters, only: Prnx, Prny, Prnz
    real(knd), intent(inout) :: U(-1:,-1:,-1:)
    integer, intent(in) :: SBtype(6)
    call par_exchange_boundaries(U, Prnx, Prny, Prnz, SBtype, -1, 2, dir=2)
  end subroutine
  
  subroutine par_exchange_Sc_z(U, SBType)
    use Parameters, only: Prnx, Prny, Prnz
    real(knd), intent(inout) :: U(-1:,-1:,-1:)
    integer, intent(in) :: SBtype(6)
    call par_exchange_boundaries(U, Prnx, Prny, Prnz, SBtype, -1, 2, dir=3)
  end subroutine
   
   
   
  subroutine par_exchange_Pr(Phi)
    use Parameters, only: We, Ea, So, No, Bo, To, MPI_PERIODIC, Prnx, Prny, Prnz, Btype
    real(knd), intent(inout), contiguous :: Phi(1:,1:,1:)
    logical :: oddx, oddy, oddz, evenx, eveny, evenz
    integer :: ierr, tag, status(MPI_STATUS_SIZE)
    integer :: nx, ny, nz
    
    nx = Prnx
    ny = Prny
    nz = Prnz
    
    oddx = mod(iim,2) == 1
    evenx = .not. oddx
    
    oddy = mod(jim,2) == 1
    eveny = .not. oddy
    
    oddz = mod(kim,2) == 1
    evenz = .not. oddz
    

    !internal boundaries
    if (oddx) then
      call send_w
    else
      call recv_e
    end if
    if (evenx) then
      call send_w
    else
      call recv_e
    end if

    if (oddy) then
      call send_s
    else
      call recv_n
    end if
    if (eveny) then
      call send_s
    else
      call recv_n
    end if

    if (oddz) then
      call send_b
    else
      call recv_t
    end if
    if (evenz) then
      call send_b
    else
      call recv_t
    end if
   

    !global domain boundaries
    if (Btype(We)==MPI_PERIODIC.or.Btype(Ea)==MPI_PERIODIC) then
      if (iim==1) then
        call send(Phi(1,1:ny,1:nz), w_rank)
      else if (iim==nxims) then
        call recv(Phi(nx+1,1:ny,1:nz), e_rank)
      end if     
    end if

    if (Btype(So)==MPI_PERIODIC.or.Btype(No)==MPI_PERIODIC) then
      if (jim==1) then
        call send(Phi(1:nx,1,1:nz), s_rank)
      else if (jim==nyims) then
        call recv(Phi(1:nx,ny+1,1:nz), n_rank)
      end if
    end if
          
    if (Btype(Bo)==MPI_PERIODIC.or.Btype(To)==MPI_PERIODIC) then
      if (kim==1) then
        call send(Phi(1:nx,1:ny,1), b_rank)
      else if (kim==nzims) then
        call recv(Phi(1:nx,1:ny,nz+1), t_rank)
      end if
    end if
    
    
  contains
  
  
    subroutine send(a,to)
      real(knd), intent(in) :: a(:,:)
      integer, intent(in) :: to

      call MPI_Send(a, size(a) , MPI_KND, to, 1, global_comm, ierr)
      if (ierr/=0) stop "error sending MPI message."
    end subroutine
    
    subroutine recv(a,from)
      real(knd), intent(out) :: a(:,:)
      integer, intent(in) :: from

      call MPI_Recv(a, size(a) , MPI_KND, from, 1, global_comm, status, ierr)
      if (ierr/=0) stop "error sending MPI message."
    end subroutine
    

    subroutine send_w
      if (iim>1) then
        call send(Phi(1,1:ny,1:nz), w_rank)
      end if
    end subroutine
     
    subroutine recv_e
      if (iim<nxims) then
        call recv(Phi(nx+1,1:ny,1:nz), e_rank)
      end if
    end subroutine
    
    subroutine send_s
      if (jim>1) then
        call send(Phi(1:nx,1,1:nz), s_rank)
      end if
    end subroutine

    subroutine recv_n
      if (jim<nyims) then
        call recv(Phi(1:nx,ny+1,1:nz), n_rank)
      end if
    end subroutine
    
    subroutine send_b
      if (kim>1) then
        call send(Phi(1:nx,1:ny,1), b_rank)
      end if
    end subroutine
    subroutine recv_t
      if (kim<nzims) then
        call recv(Phi(1:nx,1:ny,nz+1), t_rank)
      end if
    end subroutine
    
  end subroutine par_exchange_Pr
  
  
  
  
  
  
  
  
  
  
  
  
  
  subroutine par_exchange_Q(Phi)
    use Parameters, only: We, Ea, So, No, Bo, To, MPI_PERIODIC, Prnx, Prny, Prnz, Btype
    real(knd), intent(inout), contiguous :: Phi(0:,0:,0:)
    logical :: oddx, oddy, oddz, evenx, eveny, evenz
    integer :: ierr, tag, status(MPI_STATUS_SIZE)
    logical :: xdir, ydir, zdir
    integer :: nx, ny, nz
    
    ierr = 0; status=0; tag = 0
    
    nx = Prnx
    ny = Prny
    nz = Prnz
    
    oddx = mod(iim,2) == 1
    evenx = .not. oddx
    
    oddy = mod(jim,2) == 1
    eveny = .not. oddy
    
    oddz = mod(kim,2) == 1
    evenz = .not. oddz

    !internal boundaries
    tag = tag + 1
    if (oddx) then
      call recv_w
    else
      call send_e
    end if
    if (evenx) then
      call recv_w
    else
      call send_e
    end if

    tag = tag + 1
    if (oddx) then
      call recv_e
    else
      call send_w
    end if
    if (evenx) then
      call recv_e
    else
      call send_w
    end if


    tag = tag + 1
    if (oddy) then
      call recv_s
    else
      call send_n
    end if
    if (eveny) then
      call recv_s
    else
      call send_n
    end if

    tag = tag + 1
    if (oddy) then
      call recv_n
    else
      call send_s
    end if
    if (eveny) then
      call recv_n
    else
      call send_s
    end if


    tag = tag + 1
    if (oddz) then
      call recv_b
    else
      call send_t
    end if
    if (evenz) then
      call recv_b
    else
      call send_t
    end if

    tag = tag + 1
    if (oddz) then
      call recv_t
    else
      call send_b
    end if
    if (evenz) then
      call recv_t
    else
      call send_b
    end if

    !global domain boundaries
    tag = tag + 1
    if (Btype(We)==MPI_PERIODIC.or.Btype(Ea)==MPI_PERIODIC) then
      if (iim==1) then
        call recv(Phi(1,1:ny,1:nz), w_rank)
      else if (iim==nxims) then
        call send(Phi(nx+1,1:ny,1:nz), e_rank)
      end if     
      if (iim==nxims) then
        call recv(Phi(nx,1:ny,1:nz), e_rank)
      else if (iim==1) then
        call send(Phi(0,1:ny,1:nz), w_rank)
      end if
    end if

    tag = tag + 1
    if (Btype(So)==MPI_PERIODIC.or.Btype(No)==MPI_PERIODIC) then
      if (jim==1) then
        call recv(Phi(0:nx+1,1,1:nz), s_rank)
      else if (jim==nyims) then
        call send(Phi(0:nx+1,ny+1,1:nz), n_rank)
      end if
      if (jim==nyims) then
        call recv(Phi(0:nx+1,ny,1:nz), n_rank)
      else if (jim==1) then
        call send(Phi(0:nx+1,0,1:nz), s_rank)
      end if
    end if

    tag = tag + 1
    if (Btype(Bo)==MPI_PERIODIC.or.Btype(To)==MPI_PERIODIC) then
      if (kim==1) then
        call recv(Phi(0:nx+1,0:ny+1,1), b_rank)
      else if (kim==nzims) then
        call send(Phi(0:nx+1,0:ny+1,nz+1), t_rank)
      end if
      if (kim==nzims) then
        call recv(Phi(0:nx+1,0:ny+1,nz), t_rank)
      else if (kim==1) then
        call send(Phi(0:nx+1,0:ny+1,0), b_rank)
      end if
    end if
    
  contains
  
  
    subroutine send(a,to)
      real(knd), intent(in) :: a(:,:)
      integer, intent(in) :: to

      call MPI_Send(a, size(a) , MPI_KND, to, tag, global_comm, ierr)
      if (ierr/=0) stop "error sending MPI message."
    end subroutine
    
    subroutine recv(a,from)
      real(knd), intent(out) :: a(:,:)
      real(knd) :: tmp(size(a,1),size(a,2))
      integer, intent(in) :: from

      call MPI_Recv(tmp, size(a) , MPI_KND, from, tag, global_comm, status, ierr)
      if (ierr/=0) stop "error sending MPI message."
      a = a + tmp
    end subroutine
    

    subroutine recv_w
      if (iim>1) then
        call recv(Phi(1,1:ny,1:nz), w_rank)
      end if
    end subroutine
    subroutine send_w
      if (iim>1) then
        call send(Phi(0,1:ny,1:nz), w_rank)
      end if
    end subroutine
    subroutine recv_e
      if (iim<nxims) then
        call recv(Phi(nx,1:ny,1:nz), e_rank)
      end if
    end subroutine       
    subroutine send_e
      if (iim<nxims) then
        call send(Phi(nx+1,1:ny,1:nz), e_rank)
      end if
     end subroutine
    subroutine recv_s
      if (jim>1) then
        call recv(Phi(0:nx+1,1,1:nz), s_rank)
      end if
    end subroutine
    subroutine send_s
      if (jim>1) then
        call send(Phi(0:nx+1,0,1:nz), s_rank)
      end if
    end subroutine
    subroutine recv_n
      if (jim<nyims) then
        call recv(Phi(0:nx+1,ny,1:nz), n_rank)
      end if
    end subroutine
    subroutine send_n
      if (jim<nyims) then
        call send(Phi(0:nx+1,ny+1,1:nz), n_rank)
      end if
    end subroutine
    subroutine recv_b
      if (kim>1) then
        call recv(Phi(0:nx+1,0:ny+1,1), b_rank)
      end if
    end subroutine
    subroutine send_b
      if (kim>1) then
        call send(Phi(0:nx+1,0:ny+1,0), b_rank)
      end if
    end subroutine
    subroutine recv_t
      if (kim<nzims) then
        call recv(Phi(0:nx+1,0:ny+1,nz), t_rank)
      end if
    end subroutine
    subroutine send_t
      if (kim<nzims) then
        call send(Phi(0:nx+1,0:ny+1,nz+1), t_rank)
      end if
    end subroutine
    
  end subroutine par_exchange_Q
  
  
  
  
  
  
  

  
  
  subroutine par_exchange_boundaries_yz(Phi, ny, nz, Btype, lby, lbz, widthy, widthz)
    use Parameters, only: So, No, Bo, To, MPI_PERIODIC
    real(knd), intent(inout), contiguous :: Phi(lby:, lbz:)
    integer, intent(in) :: ny, nz
    integer, intent(in) :: Btype(6)
    integer, intent(in) :: lby, lbz, widthy, widthz
    logical :: oddy, oddz, eveny, evenz
    integer :: ierr, tag, status(MPI_STATUS_SIZE)
    logical :: ydir, zdir
    
    ierr = 0; status=0; tag = 0
    
    oddy = mod(jim,2) == 1
    eveny = .not. oddy
    
    oddz = mod(kim,2) == 1
    evenz = .not. oddz

    !internal boundaries

    tag = tag + 1
    if (oddy) then
      call send_s
    else
      call recv_n
    end if
    if (eveny) then
      call send_s
    else
      call recv_n
    end if

    tag = tag + 1
    if (oddy) then
      call send_n
    else
      call recv_s
    end if
    if (eveny) then
      call send_n
    else
      call recv_s
    end if


    tag = tag + 1
    if (oddz) then
      call send_b
    else
      call recv_t
    end if
    if (evenz) then
      call send_b
    else
      call recv_t
    end if

    tag = tag + 1
    if (oddz) then
      call send_t
    else
      call recv_b
    end if
    if (evenz) then
      call send_t
    else
      call recv_b
    end if

    !global domain boundaries

    tag = tag + 1
    if (Btype(So)==MPI_PERIODIC.or.Btype(No)==MPI_PERIODIC) then
      if (jim==1) then
        call send(Phi(1:0+widthy,1:nz), s_rank)
      else if (jim==nyims) then
        call recv(Phi(ny+1:ny+widthy,1:nz), n_rank)
      end if
      if (jim==nyims) then
        call send(Phi(ny+1-widthy:ny,1:nz), n_rank)
      else if (jim==1) then
        call recv(Phi(1-widthy:0,1:nz), s_rank)
      end if
    end if



    tag = tag + 1
    if (Btype(Bo)==MPI_PERIODIC.or.Btype(To)==MPI_PERIODIC) then
      if (kim==1) then
        call send(Phi(1-widthy:ny+widthy,1:0+widthz), b_rank)
      else if (kim==nzims) then
        call recv(Phi(1-widthy:ny+widthy,nz+1:nz+widthz), t_rank)
      end if
      if (kim==nzims) then
        call send(Phi(1-widthy:ny+widthy,nz+1-widthz:nz), t_rank)
      else if (kim==1) then
        call recv(Phi(1-widthy:ny+widthy,1-widthz:0), b_rank)
      end if
    end if

    
  contains
  
  
    subroutine send(a,to)
      real(knd), intent(in) :: a(:,:)
      integer, intent(in) :: to
      !Must be global_comm, for `to` and `from` derived from `x_rank` to be valid!
      call MPI_Send(a, size(a) , MPI_KND, to, tag, global_comm, ierr)
      if (ierr/=0) stop "error sending MPI message."
    end subroutine
    
    subroutine recv(a,from)
      real(knd), intent(out) :: a(:,:)
      integer, intent(in) :: from

      call MPI_Recv(a, size(a) , MPI_KND, from, tag, global_comm, status, ierr)
      if (ierr/=0) stop "error sending MPI message."
    end subroutine
    

    subroutine send_s
      if (jim>1) then
        call send(Phi(1:0+widthy,1:nz), s_rank)
      end if
    end subroutine
    subroutine recv_s
      if (jim>1) then
        call recv(Phi(1-widthy:0,1:nz), s_rank)
      end if
    end subroutine
    subroutine send_n
      if (jim<nyims) then
        call send(Phi(ny+1-widthy:ny,1:nz), n_rank)
      end if
    end subroutine
    subroutine recv_n
      if (jim<nyims) then
        call recv(Phi(ny+1:ny+widthy,1:nz), n_rank)
      end if
    end subroutine
    subroutine send_b
      if (kim>1) then
        call send(Phi(1-widthy:ny+widthy,1:0+widthz), b_rank)
      end if
    end subroutine
    subroutine recv_b
      if (kim>1) then
        call recv(Phi(1-widthy:ny+widthy,1-widthz:0), b_rank)
      end if
    end subroutine
    subroutine send_t
      if (kim<nzims) then
        call send(Phi(1-widthy:ny+widthy,nz+1-widthz:nz), t_rank)
      end if
    end subroutine
    subroutine recv_t
      if (kim<nzims) then
        call recv(Phi(1-widthy:ny+widthy,nz+1:nz+widthz), t_rank)
      end if
    end subroutine
    
  end subroutine par_exchange_boundaries_yz
  
  
end module exchange_par
module Kinds
  use iso_fortran_env
  
  integer,parameter :: rp = real32
end module

module Endianness

  use iso_fortran_env
  
  implicit none

  private
  public :: GetEndianness, BigEnd, SwapB, littleendian

  logical, save :: littleendian = .true.

  interface BigEnd
    module procedure BigEnd16
    module procedure BigEnd32
    module procedure BigEnd32_a1
    module procedure BigEnd32_a2
    module procedure BigEnd32_a3
    module procedure BigEnd32_a4
    module procedure BigEnd64
    module procedure BigEnd64_a1
    module procedure BigEnd64_a2
    module procedure BigEnd64_a3
    module procedure BigEnd64_a4
  end interface

  interface SwapB
    module procedure SwapB32
    module procedure SwapB64
  end interface

  contains

    subroutine GetEndianness
      character(4) :: bytes !may not work on some processors

      bytes = transfer(1_int32,bytes)
      if (ichar(bytes(4:4))==1) then
        littleendian=.false.
      else
        littleendian=.true.
      endif
    end subroutine GetEndianness

   
    elemental function BigEnd16(x) result(res)
      integer(int16) :: res
      integer(int16),intent(in)::x
      character(2) :: bytes
      
      if (.not.littleendian) then
        res = x
      else
        bytes = transfer(x,bytes)
        res = ichar(bytes(2:2),int16)
        res = ior( ishft(ichar(bytes(1:1),int16),8), res )
      endif
    end function
    
    function BigEnd32(x) result(res)
      real(real32),intent(in) :: x
      real(real32) :: res
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a1(x) result(res)
      real(real32),intent(in) :: x(:)
      real(real32) :: res(size(x,1))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a2(x) result(res)
      real(real32),intent(in) :: x(:,:)
      real(real32) :: res(size(x,1),size(x,2))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a3(x) result(res)
      real(real32),intent(in) :: x(:,:,:)
      real(real32) :: res(size(x,1),size(x,2),size(x,3))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd32_a4(x) result(res)
      real(real32),intent(in) :: x(:,:,:,:)
      real(real32) :: res(size(x,1),size(x,2),size(x,3),size(x,4))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64(x) result(res)
      real(real64),intent(in) :: x
      real(real64) :: res
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a1(x) result(res)
      real(real64),intent(in) :: x(:)
      real(real64) :: res(size(x,1))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a2(x) result(res)
      real(real64),intent(in) :: x(:,:)
      real(real64) :: res(size(x,1),size(x,2))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a3(x) result(res)
      real(real64),intent(in) :: x(:,:,:)
      real(real64) :: res(size(x,1),size(x,2),size(x,3))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    function BigEnd64_a4(x) result(res)
      real(real64),intent(in) :: x(:,:,:,:)
      real(real64) :: res(size(x,1),size(x,2),size(x,3),size(x,4))
      
      if (.not.littleendian) then
        res = x
      else
        res = SwapB(x)
      endif
    end function
    
    elemental function SwapB32(x) result(res)
      real(real32) :: res
      real(real32),intent(in) :: x
      character(4) :: bytes
      integer(int32) :: t
      real(real32) :: rbytes, rt
      !nicer looking TRANSFER is problematic to optimize by ifort
      ! gfortran would be OK with that
      equivalence (rbytes, bytes)
      equivalence (t, rt)
      
      rbytes = x

      t = ichar(bytes(4:4),int32)

      t = ior( ishftc(ichar(bytes(3:3),int32),8),  t )

      t = ior( ishftc(ichar(bytes(2:2),int32),16), t )

      t = ior( ishftc(ichar(bytes(1:1),int32),24), t )

      res = rt
        
    end function
    
    elemental function SwapB64(x) result(res)
      real(real64) :: res
      real(real64),intent(in) :: x
      character(8) :: bytes
      integer(int64) :: t
      real(real64) :: rbytes, rt
      equivalence (rbytes, bytes)
      equivalence (t, rt)
      
      rbytes = x

      t = ichar(bytes(8:8),int64)

      t = ior( ishftc(ichar(bytes(7:7),int64),8),  t )

      t = ior( ishftc(ichar(bytes(6:6),int64),16), t )

      t = ior( ishftc(ichar(bytes(5:5),int64),24), t )

      t = ior( ishftc(ichar(bytes(4:4),int64),32), t )

      t = ior( ishftc(ichar(bytes(3:3),int64),40), t )

      t = ior( ishftc(ichar(bytes(2:2),int64),48), t )

      t = ior( ishftc(ichar(bytes(1:1),int64),56), t )

      res = rt

    end function

end module Endianness




module Types
  use Kinds
  
  implicit none

  type grid
    integer :: nx,ny,nz
    integer :: offx,offy,offz
    integer :: unit
    real(rp),allocatable :: x(:),y(:),z(:)
    character(100) :: fname
  contains
    procedure :: read_header
    procedure :: read_title
    procedure :: read_scalar
    procedure :: read_vector
  end type
  
  integer, parameter :: SCALAR = 7, VECTOR = 8
  character,parameter :: lf = achar(10)

contains
  

  subroutine read_header(g)
    use Endianness
    class(grid),intent(inout) :: g
    integer :: io

    open(newunit=g%unit,file=g%fname,access="stream",form="unformatted",status="old",action="read",iostat=io)

    if (io/=0) then
      g%nx = 0
      g%ny = 0
      g%nz = 0
      return
    end if

    call skip_to("X_COORDINATES", io)
    if (io/=0) then
      write(*,*) "Error, cannot find X_COORDINATES in '",g%fname,"'"
      stop
    end if
    call get_number(g%nx)
    call skip_line
    
    allocate(g%x(g%nx))
    read(g%unit,iostat=io) g%x
    g%x=BigEnd(g%x)

    call skip_to("Y_COORDINATES", io)
    if (io/=0) then
      write(*,*) "Error, cannot find Y_COORDINATES in '",g%fname,"'"
      stop
    end if
    call get_number(g%ny)
    call skip_line
    
    allocate(g%y(g%ny))
    read(g%unit,iostat=io) g%y
    g%y=BigEnd(g%y)

    call skip_to("Z_COORDINATES", io)
    if (io/=0) then
      write(*,*) "Error, cannot find Z_COORDINATES in '",g%fname,"'"
      stop
    end if
    call get_number(g%nz)
    call skip_line
    
    allocate(g%z(g%nz))
    read(g%unit,iostat=io) g%z
    g%z=BigEnd(g%z)
    
    call skip_line
    call skip_line
    
  contains
    subroutine skip_line
      character :: ch
      do
        read(g%unit) ch
        if (ch==new_line("a")) return
      end do
    end subroutine

    subroutine get_number(n)
      integer, intent(out) :: n
      character :: ch
      integer :: io
      character(:), allocatable :: num_str
      
      do
        read(g%unit, iostat=io) ch
        if (io/=0) then
          write(*,*) "Error reading from '",g%fname,"'."
          stop 3
        end if
        if (ch/=' ') exit
      end do
      
      num_str=ch
      do
        read(g%unit) ch
        if (ch<'0' .or. ch>'9') exit
        num_str = num_str // ch
      end do

      read(num_str,*) n
    end subroutine

    subroutine skip_to(str, stat)
      character(*), intent(in) :: str
      integer, intent(out) :: stat
      character :: ch
      integer :: io

      do
        read(g%unit, iostat=io) ch

        if (io/=0) then
          stat = 1
          return
        end if

        if (ch==str(1:1)) then
          call check(str(2:), stat)
          if (stat == 0) return
        end if

      end do
    end subroutine
    
    subroutine check(str, stat)
      character(*), intent(in) :: str
      integer, intent(out) :: stat
      character :: ch
      integer :: i, io

      stat = 1
      i = 0

      do
        i = i + 1

        read(g%unit, iostat=io) ch

        if (io/=0) return

        if (ch/=str(i:i)) return

        if (i==len(str)) then
          stat = 0
          return
        end if
      end do
    end subroutine
   
  end subroutine
  
  subroutine read_title(g,status,title)
    use iso_fortran_env
    class(grid),intent(in) :: g
    integer :: status
    character(:),allocatable,intent(out) :: title
    character(7) :: vtype
    character :: ch
    integer :: io, n, lf_pos
    character(256) :: msg
  
    do
      vtype=""
      read(g%unit,iostat=io,iomsg=msg) vtype
      
      if (io==iostat_end) then
        status = 0
        return
      else if (io/=0) then
        write(*,'(*(g0))') io, trim(msg), "  ", g%fname, " ", "'",vtype,"'"
        stop "Error reading title."
      end if
      
      if (vtype=='VECTORS') then
      
        status = VECTOR
        title = vtype
        do
          read(g%unit) ch
          title = title // ch
          if (ch==lf) exit
        end do
        
        exit
        
      else if (vtype=='SCALARS') then
      
        status = SCALAR
        title = vtype
        n = 0
        do
          read(g%unit) ch
          title = title // ch
          if (ch==lf) n = n + 1
          if (n==2) exit
        end do
        
        exit
        
      else
        !if our vtype buffer contains the end-of-line LF character, move back just behind it
        !otherwise skip to the nearest one
        lf_pos = scan(vtype, achar(10),.true.)
        if (lf_pos>0) then
          call stream_back(len(vtype)-lf_pos+1)
        else
          call skip_line
        end if
        
      end if
    end do
    
  contains
  
    subroutine skip_line
      character :: ch
      do
        read(g%unit) ch
        if (ch==new_line("a")) return
      end do
    end subroutine
    
    subroutine stream_back(n)
      integer, intent(in) :: n
      integer :: pos
      character :: dummy
      inquire(unit=g%unit,pos=pos)
      read(g%unit,pos=pos-n) dummy
    end subroutine

  end subroutine
  
  subroutine read_scalar(g,buf)
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,:,:)
    character :: ch
    
    read(g%unit) buf(1:g%nx,1:g%ny,1:g%nz)
    read(g%unit) ch
  end subroutine
  
  subroutine read_vector(g,buf)
  use Endianness
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,:,:,:)
    character :: ch

    read(g%unit) buf(:,1:g%nx,1:g%ny,1:g%nz)
    read(g%unit) ch
  end subroutine
  
end module Types






program spectrum_1D
  use fftw3
  use Kinds
  use Types
  use Endianness
  
  implicit none
  
  integer :: nx, ny, nz
  real(rp), allocatable :: sc(:,:,:),vec(:,:,:,:)
 
  character(:), allocatable :: title
  integer :: status, var_type, arg_len, title_pos
  
  character(:), allocatable :: in_name, out_name, arg
  
  type(grid) :: in_g
  
  integer :: dir
  
  call GetEndianness

  if (command_argument_count()<1) then
    write(*,*) "Usage: 1D_spectrum U.vtk spectrum.txt [direction]"
    stop 1
  end if
 
  call get_command_argument(1, length=arg_len)
  allocate(character(arg_len) :: in_name)
  call get_command_argument(1, value=in_name)

  call get_command_argument(2, length=arg_len)
  allocate(character(arg_len) :: out_name)
  call get_command_argument(2, value=out_name)
  
  if (command_argument_count()>2) then
    call get_command_argument(3, length=arg_len)
    allocate(character(arg_len) :: arg)
    call get_command_argument(3, value=arg)
    read(arg,*) dir
  else
    dir = 1
  end if
  
  in_g%fname = in_name
  call in_g%read_header
  
  
  nx = in_g%nx
  ny = in_g%ny
  nz = in_g%nz
  write(*,*) in_g%nx, in_g%ny, in_g%nz
  
  allocate(sc(1:nx,1:ny,1:nz))
  allocate(vec(1:3,1:nx,1:ny,1:nz))
  
  do
    call in_g%read_title(status,title)

    if (status==SCALAR) then
      call get_buffer(status,sc,vec)
    
      sc = BigEnd(sc)
    
      call do_transform(sc, dir)
    else if (status==VECTOR) then
      call get_buffer(status,sc,vec)
    
      sc = BigEnd(vec(1,:,:,:))
    
      call do_transform(sc, dir)
    end if
    exit
  end do
  
contains

  subroutine do_transform(u3d, dir)
    real(rp), intent(in) :: u3d(:,:,:)
    integer, intent(in) :: dir    
    real(rp), allocatable :: u(:), sp(:), sp_avg(:)
    complex(rp), allocatable :: u_hat(:)
    integer :: n, nj, nk
    type(c_ptr) :: forw
    integer :: i, j, k, iu
    real(rp) :: dx, lx, u_var, u_m
    real(rp), parameter :: pi = acos(-1.0_rp)
    
    n = size(u3d,dir)
    
    select case (dir)
      case (1)
        nj = size(u3d, 2)
        nk = size(u3d, 3)
      case (2)
        nj = size(u3d, 1)
        nk = size(u3d, 3)
      case (3)
        nj = size(u3d, 1)
        nk = size(u3d, 3)
    end select
    
    allocate(u(n))
    allocate(u_hat(0:n/2))
    allocate(sp(n/2))
    allocate(sp_avg(n/2))
    
    sp_avg = 0
    u_var = 0
    
    forw = fftw_plan_gen(size(u), &
                u, u_hat, FFTW_UNALIGNED)
                
    do k = 1, nk
      do j = 1, nj
        select case (dir)
          case (1)
            u = u3d(:,j,k)
          case (2)
            u = u3d(j,:,k)
          case (3)
            u = u3d(j,k,:)
        end select
        
#ifdef DPREC
        call fftw_execute_dft_r2c(forw, u, u_hat)
#else
        call fftwf_execute_dft_r2c(forw, u, u_hat)
#endif
        sp = real(u_hat(1:) * conjg(u_hat(1:)))
        sp = sp / n
        
        sp_avg = sp_avg + sp
        
        u_m = sum(u)/n
        do i = 1, n
          u_var = u_var + (u(i)-u_m)**2
        end do
      end do
    end do
    
    sp_avg = sp_avg / (nj*nk)
    
    u_var = u_var / size(u3d)
    
    write(*,*) "variance:", u_var

    select case (dir)
      case(1)
        lx = (in_g%x(n)-in_g%x(1))
        dx = (in_g%x(n)-in_g%x(1))/n
      case(2)
        lx = (in_g%y(n)-in_g%y(1))
        dx = (in_g%y(n)-in_g%y(1))/n
      case(3)
        lx = (in_g%z(n)-in_g%z(1))
        dx = (in_g%z(n)-in_g%z(1))/n
    end select
    
    open(newunit=iu,file=out_name)
    do j = 1, size(sp_avg)
      write(iu,*) 2*pi*j/lx, dx*sp_avg(j)/pi
    end do 
    close(iu)
    
    write(*,*) "integral:", 2*pi*n/lx   *  dx*sum(sp_avg)/n/pi
  end subroutine
  
  subroutine get_buffer(status,sc,vec)
    integer,intent(in) :: status
    real(rp),intent(out) :: sc(:,:,:),vec(:,:,:,:)

    if (status==SCALAR) then
      call in_g%read_scalar(sc)
    else
      call in_g%read_vector(vec)
    end if
  end subroutine
      
end program

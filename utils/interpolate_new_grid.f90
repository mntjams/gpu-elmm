module WorkKinds
  use iso_fortran_env
  
  integer,parameter :: rp = real32
end module

module Endianness
  use iso_fortran_env
  
  implicit none

  private
  public :: GetEndianness, BigEnd

  logical,save :: littleendian = .true.

  interface BigEnd
    module procedure BigEnd32
    module procedure BigEnd64
  end interface

contains

  subroutine GetEndianness
    integer(int8),dimension(4):: bytes !may not work on some processors

    bytes=transfer(1_int32,bytes,4)
    if (bytes(4)==1) then
      littleendian=.false.
    else
      littleendian=.true.
    endif
  end subroutine GetEndianness

  elemental function BigEnd32(x) result(res)
    real(real32) :: res
    real(real32),intent(in)::x
    integer(int8),dimension(4):: bytes !may not work on some processors

    if (.not.littleendian) then
      res = x
    else
      bytes = transfer(x,bytes,4)
      res = transfer(bytes(4:1:-1),res)
    endif
  end function BigEnd32

  elemental function BigEnd64(x) result(res)
    real(real64) :: res
    real(real64),intent(in)::x
    integer(int8),dimension(8):: bytes !may not work on some processors

    if (.not.littleendian) then
      res = x
    else
      bytes = transfer(x,bytes,8)
      res = transfer(bytes(8:1:-1),res)
    endif
  end function BigEnd64

end module Endianness




module Types
  use WorkKinds
  
  implicit none

  type grid
    integer nx,ny,nz
    integer offx,offy,offz
    integer unit
    real(rp),allocatable :: x(:),y(:),z(:)
    real(rp) :: x1, y1, z1
    real(rp) :: dx, dy, dz
    character(100) :: fname
  contains
    procedure :: read_header
    procedure :: read_title
    procedure :: read_scalar
    procedure :: read_vector
  end type
  
  integer, parameter :: SCALAR = 1, VECTOR = 2
  character,parameter :: lf = achar(10)

contains
  

  subroutine read_header(g)
    use Endianness
    class(grid),intent(inout) :: g
    integer :: io,nx,ny,nz
    character(5) :: ch5
    character(70) :: str
    character :: ch

    open(newunit=g%unit,file=g%fname,access="stream",form="unformatted",status="old",action="read",iostat=io)

    if (io/=0) then
      g%nx = 0
      g%ny = 0
      g%nz = 0
      return
    end if

    read(g%unit,pos=162,iostat=io) ch5
    if (io/=0) call error

    read(ch5,'(i5)') g%nx
    nx=g%nx
    allocate(g%x(g%nx))
    read(g%unit,pos=219,iostat=io) g%x
    g%x=BigEnd(g%x)

    read(g%unit,pos=234+nx*4,iostat=io) ch5
    read(ch5,'(i5)') g%ny
    ny=g%ny

    allocate(g%y(g%ny))
    read(g%unit,pos=291+nx*4,iostat=io) g%y
    g%y=BigEnd(g%y)

    read(g%unit,pos=306+nx*4+ny*4,iostat=io) ch5
    read(ch5,'(i5)') g%nz
    nz=g%nz
    allocate(g%z(g%nz))
    read(g%unit,pos=363+nx*4+ny*4,iostat=io) g%z
    g%z=BigEnd(g%z)
    read(g%unit) ch,str
    read(g%unit) ch
    
    g%x1 = g%x(1)
    g%y1 = g%y(1)
    g%z1 = g%z(1)
    
    g%dx = (g%x(nx) - g%x(1)) / (nx-1)
    g%dy = (g%x(ny) - g%y(1)) / (ny-1)
    g%dz = (g%z(nz) - g%z(1)) / (nz-1)
    
  contains
    subroutine error()
      write(*,*) "Error reading from ",trim(g%fname)
      stop 2
    end subroutine    
  end subroutine
  
  subroutine read_title(g,status,title)
    use iso_fortran_env
    class(grid),intent(in) :: g
    integer :: status
    character(:),allocatable,intent(out) :: title
    character(7) :: vtype
    character :: ch
    integer :: io, n
    character(256) :: msg
    
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
    else
      write(*,'(*(g0))') "'",vtype,"'"
      stop "Error reading title."
    end if
  end subroutine
  
  subroutine read_scalar(g,buf)
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,:,:)
    character :: ch
 
    read(g%unit) buf(g%offx+1:g%offx+g%nx, g%offy+1:g%offy+g%ny, g%offz+1:g%offz+g%nz)
    read(g%unit) ch
  end subroutine
  
  subroutine read_vector(g,buf)
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,:,:,:)
    character :: ch
    
    read(g%unit) buf(:, g%offx+1:g%offx+g%nx, g%offy+1:g%offy+g%ny, g%offz+1:g%offz+g%nz)
    read(g%unit) ch
  end subroutine
  
end module Types



program interpolate%new%grid
  use WorkKinds
  use Types
  use Endianness
  use Parameters
  use Strings
  
  implicit none
  
  type(grid) :: old, new

  real(rp), allocatable :: sc(:,:,:),vec(:,:,:,:)
  
  integer :: it
 
  call GetEndianness

  if (command_argument_count()<2) then
    write(*,*) "Usage: joinvtk file_name npx,npy,npz"
    stop 1
  end if
 
  call get_command_argument(1, value=file_name)
  base_name = trim(file_name(scan(file_name,'/',back=.true.) + 1 :  ))

  call get_command_argument(2, value=arg)
  read(arg,*,iostat=io) nx_new, ny_new, nz_new
  
  call execute_command_line("mkdir -p grid-"//itoa(nx_new)//"-"//itoa(ny_new)//"-"//itoa(nz_new))
  
  call old%read_header
  
  lo = old%x(1) - old%dx / 2
  up = old%x(old%nx) + old%dx / 2
  
  new%dx = (up - lo) / nx_new
  new%x = [ (lo + new%dx * (it - 0.5_rp), it = 1, nx) ]
  new%x1 = new%x(1)
  
  lo = old%y(1) - old%dy / 2
  up = old%y(old%ny) + old%dy / 2
  
  new%dy = (up - lo) / ny_new
  new%y = [ (lo + new%dy * (it - 0.5_rp), it = 1, ny) ]
  new%y1 = new%y(1)
  
  lo = old%z(1) - old%dz / 2
  up = old%z(old%nz) + old%dz / 2
  
  new%dz = (up - lo) / nz_new
  new%z = [ (lo + new%dz * (it - 0.5_rp), it = 1, nz) ]
  new%z1 = new%z(1)
  
  call save_header


  do
    call get_next(status,title)
    if (status==0) exit

    call get_buffer(status,sc,vec)

    call save_buffer(status,title,sc,vec)
  end do
  
  close(old%unit)
  close(new%unit)

  
contains

  subroutine get_next(status,title)
    integer :: status
    character(:),allocatable,intent(out) :: title
    integer :: i,j,k


    call old%read_title(status,title)

  end subroutine
  
  subroutine get_buffer(status,sc,vec)
    integer,intent(in) :: status
    real(rp),intent(out) :: sc(:,:,:),vec(:,:,:,:)
    integer :: i,j,k
    

    if (status==SCALAR) then
      call old%read_scalar(sc)
    else
      call old%read_vector(vec)
    end if
  end subroutine
  
  subroutine save_header
    integer i
    character(70) :: str
   
    open(newunit=new%unit,file=file_name, &
      access='stream',status='replace',form="unformatted",action="write")

    write(new%unit) "# vtk DataFile Version 2.0",lf
    write(new%unit) "CLMM output file",lf
    write(new%unit) "BINARY",lf
    write(new%unit) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write(str(12:),*) new%nx, new%ny, new%nz
    write(new%unit) str,lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') new%nx,"float"
    write(new%unit) str,lf
    write(new%unit) BigEnd(real(x, real32)),lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') new%ny,"float"
    write(new%unit) str,lf
    write(new%unit) BigEnd(real(y, real32)),lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') new%nz,"float"
    write(new%unit) str,lf
    write(new%unit) BigEnd(real(z, real32)),lf
    str="POINT_DATA"
    write(str(12:),*) new%nx * new%ny * new%nz
    write(new%unit) str,lf 
  end subroutine
  
  subroutine save_buffer(status,title,sc,vec)
    integer,intent(in) :: status
    character(*),intent(in) :: title
    real(rp),intent(in) :: sc(:,:,:),vec(:,:,:,:)
    
    write(new%unit) title
    if (status==SCALAR) then
      write(new%unit) sc,lf
    else
      write(new%unit) vec,lf
    end if
  end subroutine

  subroutine old_index(x, y, z, xi, yj, zk)
    xi = min(max(nint( (x - old%x(1))/old%dx ), 1), old%nx-1)
    yj = min(max(nint( (y - old%y(1))/old%dy ), 1), old%nx-1)
    zk = min(max(nint( (z - old%z(1))/old%dz ), 1), old%nx-1)
  end subroutine


end program

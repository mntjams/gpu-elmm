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


module Interpolation
  use WorkKinds

  implicit none

contains

  pure real(rp) function TriLinInt(a, b, c, &
                                   vel000, vel100, vel010, vel001, vel110, vel101, vel011, vel111)
    real(rp), intent(in) :: a, b, c
    real(rp), intent(in) :: vel000, vel100, vel010, vel001, vel110, vel101, vel011, vel111

    TriLinInt =  (1-a) * (1-b) * (1-c) * vel000 + &
                 a     * (1-b) * (1-c) * vel100 + &
                 (1-a) * b     * (1-c) * vel010 + &
                 (1-a) * (1-b) * c     * vel001 + &
                 a     * b     * (1-c) * vel110 + &
                 a     * (1-b) * c     * vel101 + &
                 (1-a) * b     * c     * vel011 + &
                 a     * b     * c     * vel111

  end function TriLinInt

end module


module Types
  use WorkKinds
  use Endianness
  
  implicit none

  type grid
    integer nx,ny,nz
    integer unit
    real(rp),allocatable :: x(:),y(:),z(:)
    real(rp) :: x1, y1, z1
    real(rp) :: dx, dy, dz
    character(1024) :: fname
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
    g%dy = (g%y(ny) - g%y(1)) / (ny-1)
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
    real(rp),intent(out) :: buf(-1:,-1:,-1:)
    character :: ch
 
    read(g%unit) buf(1:g%nx, 1:g%ny, 1:g%nz)
    buf = BigEnd(buf)
    read(g%unit) ch
  end subroutine
  
  subroutine read_vector(g,buf)
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,-2:,-2:,-2:)
    character :: ch
    
    read(g%unit) buf(:, 1:g%nx, 1:g%ny, 1:g%nz)
    buf = BigEnd(buf)
    read(g%unit) ch
  end subroutine
  
end module Types



program interpolate_new_grid
  use WorkKinds
  use Types
  use Endianness
  use Interpolation
  use Strings
  use Boundaries
  
  implicit none
  
  type(grid) :: old, new

  real(rp), allocatable :: old_sc(:,:,:), old_vec(:,:,:,:)
  real(rp), allocatable :: new_sc(:,:,:), new_vec(:,:,:,:)

  character(1024) :: file_name, arg

  character(:), allocatable :: base_name, dir_name, title
  
  integer :: it, io, status

  real(rp) :: lo, up
 
  call GetEndianness

  if (command_argument_count()<2) then
    write(*,*) "Usage: joinvtk file_name npx,npy,npz"
    stop 1
  end if
 
  call get_command_argument(1, value=file_name)
  base_name = trim(file_name(scan(file_name,'/',back=.true.) + 1 :  ))

  call get_command_argument(2, value=arg)
  read(arg,*,iostat=io) new%nx, new%ny, new%nz
  
  dir_name = "grid-"//itoa(new%nx)//"-"//itoa(new%ny)//"-"//itoa(new%nz)

  call execute_command_line("mkdir -p "//dir_name)
  
  old%fname = file_name
  call old%read_header

  new%fname = dir_name // '/' // base_name
  
  lo = old%x(1) - old%dx / 2
  up = old%x(old%nx) + old%dx / 2
  
  new%dx = (up - lo) / new%nx
  new%x = [ (lo + new%dx * (it - 0.5_rp), it = 1, new%nx) ]
  new%x1 = new%x(1)
  
  lo = old%y(1) - old%dy / 2
  up = old%y(old%ny) + old%dy / 2
  
  new%dy = (up - lo) / new%ny
  new%y = [ (lo + new%dy * (it - 0.5_rp), it = 1, new%ny) ]
  new%y1 = new%y(1)
  
  lo = old%z(1) - old%dz / 2
  up = old%z(old%nz) + old%dz / 2
  
  new%dz = (up - lo) / new%nz
  new%z = [ (lo + new%dz * (it - 0.5_rp), it = 1, new%nz) ]
  new%z1 = new%z(1)
  
  call save_header


  allocate(old_sc(-2:old%nx+3,-2:old%ny+3,-2:old%nz+3), &
           old_vec(3,-2:old%nx+3,-2:old%ny+3,-2:old%nz+3))
  allocate(new_sc(-2:new%nx+3,-2:new%ny+3,-2:new%nz+3), &
           new_vec(3,-2:new%nx+3,-2:new%ny+3,-2:new%nz+3))

  do
    call get_next(status,title)
    if (status==0) exit

    call get_buffer(status, old_sc, old_vec)

    call interpolate_buffer(status, new_sc, new_vec, old_sc, old_vec)

    call save_buffer(status,title, new_sc, new_vec)
  end do
  
  close(old%unit)
  close(new%unit)

  
contains

  subroutine get_next(status,title)
    integer :: status
    character(:),allocatable,intent(out) :: title

    call old%read_title(status,title)

  end subroutine
  
  subroutine get_buffer(status,sc,vec)
    integer,intent(in) :: status
    real(rp),intent(out) :: sc(:,:,:),vec(:,:,:,:)
    
    if (status==SCALAR) then
      call old%read_scalar(sc)
    else
      call old%read_vector(vec)
    end if
  end subroutine
  
  subroutine save_header
    character(70) :: str
   
    open(newunit=new%unit,file=new%fname, &
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
    write(new%unit) BigEnd(real(new%x, real32)),lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') new%ny,"float"
    write(new%unit) str,lf
    write(new%unit) BigEnd(real(new%y, real32)),lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') new%nz,"float"
    write(new%unit) str,lf
    write(new%unit) BigEnd(real(new%z, real32)),lf
    str="POINT_DATA"
    write(str(12:),*) new%nx * new%ny * new%nz
    write(new%unit) str,lf 
  end subroutine
  
  subroutine save_buffer(status,title,sc,vec)
    integer, intent(in) :: status
    character(*), intent(in) :: title
    real(rp), intent(inout) :: sc(-1:,-1:,-1:), vec(:,-2:,-2:,-2:)
    
    write(new%unit) title
    if (status==SCALAR) then
      sc = BigEnd(sc)
      write(new%unit) sc(1:new%nx, 1:new%ny, 1:new%nz),lf
    else
      vec = BigEnd(vec)
      write(new%unit) vec(:,1:new%nx, 1:new%ny, 1:new%nz),lf
    end if
  end subroutine

  subroutine old_index(x, y, z, xi, yj, zk)
    real(rp), intent(in) :: x, y, z
    integer, intent(out) :: xi, yj, zk

    xi = min(max(floor( (x - old%x(1))/old%dx )+1, 1), old%nx-1)
    yj = min(max(floor( (y - old%y(1))/old%dy )+1, 1), old%ny-1)
    zk = min(max(floor( (z - old%z(1))/old%dz )+1, 1), old%nz-1)
  end subroutine

  function interpolate_trilinear(arr, lb, x, y, z) result(res)
    real(rp) :: res
    real(rp), intent(in) :: arr(lb:,lb:,lb:)
    integer, intent(in) :: lb
    real(rp), intent(in) :: x, y, z
    integer :: xi, yj, zk

    call old_index(x, y, z, xi, yj, zk)

    res = TriLinInt((x - old%x(xi)) / old%dx, &
                    (y - old%y(yj)) / old%dy, &
                    (z - old%z(zk)) / old%dz, &
                    arr(xi, yj, zk), &
                    arr(xi+1, yj  , zk  ), &
                    arr(xi  , yj+1, zk  ), &
                    arr(xi  , yj  , zk+1), &
                    arr(xi+1, yj+1, zk  ), &
                    arr(xi+1, yj  , zk+1), &
                    arr(xi  , yj+1, zk+1), &
                    arr(xi+1, yj+1, zk+1))
  end function

  subroutine interpolate_scalar(new_sc, old_sc)
    real(rp), intent(out) :: new_sc(-1:,-1:,-1:)
    real(rp), intent(in)  :: old_sc(-2:,-2:,-2:)
    integer :: i, j, k
    do k = 1, new%nz
      do j = 1, new%ny
        do i = 1, new%nx
          new_sc(i,j,k) = interpolate_trilinear(old_sc, -1, new%x(i), new%y(j), new%z(k))
        end do
      end do
    end do
  end subroutine

  subroutine interpolate_vector(new_vec, old_vec)
    real(rp),intent(out) :: new_vec(:,-2:,-2:,-2:)
    real(rp),intent(in)  :: old_vec(:,-2:,-2:,-2:)
    integer :: i, j, k, comp
    do k = 1, new%nz
      do j = 1, new%ny
        do i = 1, new%nx
          do comp = 1, 3
            new_vec(comp,i,j,k) = interpolate_trilinear(old_vec(comp,:,:,:), -2, new%x(i), new%y(j), new%z(k))
          end do
        end do
      end do
    end do
  end subroutine

  subroutine interpolate_buffer(status, new_sc, new_vec, old_sc, old_vec)
    integer,intent(in) :: status
    real(rp),intent(out) :: new_sc(:,:,:), new_vec(:,:,:,:)
    real(rp),intent(in) :: old_sc(:,:,:), old_vec(:,:,:,:)
    
    if (status==SCALAR) then
      call interpolate_scalar(new_sc, old_sc)
    else
      call interpolate_vector(new_vec, old_vec)
    end if
  end subroutine

end program

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
    return
    
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

    read(g%unit) buf(1:g%nx,1:g%ny,1:g%nz)
    read(g%unit) ch
  end subroutine
  
  subroutine read_vector(g,buf)
    class(grid),intent(in) :: g
    real(rp),intent(out) :: buf(:,:,:,:)
    character :: ch
    
    read(g%unit) buf(:,1:g%nx,1:g%ny,1:g%nz)
    read(g%unit) ch
  end subroutine
  
end module Types





program vtk_subset
  use Kinds
  use Types
  use Endianness
  
  implicit none
  
  integer :: out_unit
  integer :: nx_orig, ny_orig, nz_orig
  integer :: nx, ny, nz
  real(rp), allocatable :: sc(:,:,:),vec(:,:,:,:)
 
  character(:), allocatable :: title
  integer :: status, var_type, arg_len, title_pos
  
  character(:), allocatable :: in_name, out_name, extent, fields
  
  character(256), allocatable :: fields_arr(:)
  
  integer :: i1, i2, j1, j2, k1, k2, ifield
  
  type(grid) :: in_g
  
  call GetEndianness

  if (command_argument_count()<2) then
    write(*,*) "Usage: vtk_subset in_name out_name extent [field]"
    stop 1
  end if
 
  call get_command_argument(1, length=arg_len)
  allocate(character(arg_len) :: in_name)
  call get_command_argument(1, value=in_name)

  call get_command_argument(2, length=arg_len)
  allocate(character(arg_len) :: out_name)
  call get_command_argument(2, value=out_name)

  call get_command_argument(3, length=arg_len)
  allocate(character(arg_len) :: extent)
  call get_command_argument(3, value=extent)
  read(extent,*) i1,i2,j1,j2,k1,k2
  
  nx = i2 - i1 + 1
  ny = j2 - j1 + 1
  nz = k2 - k1 + 1

  if (command_argument_count()==4) then
    call get_command_argument(4, length=arg_len)
    allocate(character(arg_len) :: fields)
    call get_command_argument(4, value=fields)
    fields_arr = split(fields,",")
  else
    allocate(fields_arr(0))
  end if
  
  in_g%fname = in_name
  call in_g%read_header
  

  nx_orig = in_g%nx
  ny_orig = in_g%ny
  nz_orig = in_g%nz
  allocate(sc(1:nx_orig,1:ny_orig,1:nz_orig))
  allocate(vec(1:3,1:nx_orig,1:ny_orig,1:nz_orig))
  
  call save_header
  
  if (size(fields_arr)>0) then
    do ifield = 1, size(fields_arr)
      call find(trim(fields_arr(ifield))//" float", status, var_type)   
      call in_g%read_title(status,title)
          
      call get_buffer(var_type,sc,vec)
      call save_buffer(var_type,title,sc,vec)
    end do

  else
  
    do
      call in_g%read_title(status,title)
      if (status==0) exit
      call get_buffer(status,sc,vec)
      call save_buffer(status,title,sc,vec)
    end do
 
  end if
  
  close(in_g%unit)
  close(out_unit)

contains


  subroutine get_buffer(status,sc,vec)
    integer,intent(in) :: status
    real(rp),intent(out) :: sc(:,:,:),vec(:,:,:,:)

    if (status==SCALAR) then
      call in_g%read_scalar(sc)
    else
      call in_g%read_vector(vec)
    end if
  end subroutine
  
  subroutine save_header
    character(70) :: str
    real(rp), allocatable :: x(:),y(:),z(:)
  
    allocate(x(0), y(0), z(0))
    x = in_g%x(i1:i2)
    y = in_g%y(j1:j2)
    z = in_g%z(k1:k2)    

    open(newunit=out_unit,file=out_name, &
      access='stream',status='replace',form="unformatted",action="write")

    write(out_unit) "# vtk DataFile Version 2.0",lf
    write(out_unit) "CLMM output file",lf
    write(out_unit) "BINARY",lf
    write(out_unit) "DATASET RECTILINEAR_GRID",lf
    str="DIMENSIONS"
    write(str(12:),*) nx,ny,nz
    write(out_unit) str,lf
    str="X_COORDINATES"
    write(str(15:),'(i5,2x,a)') nx,"float"
    write(out_unit) str,lf
    write(out_unit) BigEnd(real(x, real32)),lf
    str="Y_COORDINATES"
    write(str(15:),'(i5,2x,a)') ny,"float"
    write(out_unit) str,lf
    write(out_unit) BigEnd(real(y, real32)),lf
    str="Z_COORDINATES"
    write(str(15:),'(i5,2x,a)') nz,"float"
    write(out_unit) str,lf
    write(out_unit) BigEnd(real(z, real32)),lf
    str="POINT_DATA"
    write(str(12:),*) nx*ny*nz
    write(out_unit) str,lf 
  end subroutine
  
  subroutine save_buffer(status,title,sc,vec)
    integer,intent(in) :: status
    character(*),intent(in) :: title
    real(rp),intent(in) :: sc(:,:,:),vec(:,:,:,:)
    
    write(out_unit) title
    if (status==SCALAR) then
      write(out_unit) sc(i1:i2,j1:j2,k1:k2),lf
    else
      write(out_unit) vec(:,i1:i2,j1:j2,k1:k2),lf
    end if
  end subroutine
  
  subroutine find(str, stat, var_type)
    character(*), intent(in) :: str
    character :: ch
    integer, intent(out) :: stat, var_type
    
    rewind(in_g%unit)
    call skip_to("SCALARS "//str, stat)
    if (stat==0) then
      read(in_g%unit, pos=title_pos) ch
      var_type = SCALAR    
      return
    end if
    rewind(in_g%unit)
    call skip_to("VECTORS "//str, stat)
    if (stat==0) then
      read(in_g%unit, pos=title_pos) ch
      var_type = VECTOR
      return
    else
      write(*,*) "Variable declaration '",str,"' not found in the vtk file."
      stop
    end if
  end subroutine
  
  subroutine skip_to(str, stat)
    character(*), intent(in) :: str
    integer, intent(out) :: stat
    character :: ch
    integer :: io

    do
      read(in_g%unit, iostat=io) ch

      if (io/=0) then
        stat = 1
        return
      end if

      if (ch==str(1:1)) then
        inquire(unit=in_g%unit, pos=title_pos)
        title_pos = title_pos - 2
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

      read(in_g%unit, iostat=io) ch

      if (io/=0) return

      if (ch/=str(i:i)) return

      if (i==len(str)) then
        stat = 0
        return
      end if
    end do
  end subroutine
  
  function split(str,sep) result(res)
    character(256), allocatable :: res(:)
    character(*), intent(in) :: str
    character, intent(in) :: sep
    integer :: i, n, pos, istr
    n = 1
    do i = 1, len(str)
      if (str(i:i)==sep) then
        n = n + 1
      end if
    end do
    
    allocate(res(n))
    res = ""
    
    pos = 1
    do istr = 1,n
      i = pos
      do
        if (str(i:i)==sep) then
          res(istr) = str(pos:i-1)
          pos = i + 1
          exit
        else if (i==len(str)) then
          res(istr) = str(pos:i)
          pos = i
          exit
        end if
        i = i + 1
      end do
    end do
  end function
    
end program

module ElevationModels
  use Kinds
  
  implicit none

  type,abstract :: DEM
  contains
    procedure(elevation_interface),deferred :: elevation
  end type
  
  interface
    function elevation_interface(self, x, y) result(res)
      import
      real(knd) :: res
      class(DEM), intent(in) :: self
      real(knd), intent(in) :: x, y
    end function
  end interface
  
  type, extends(DEM) :: uniform_DEM
    private
    real(knd), allocatable :: z(:,:)
    real(knd) :: dx, dy, x1, y1, xn, yn    
  contains
    procedure :: elevation => uniform_DEM_elevation
  end type
  
  interface uniform_DEM
    module procedure uniform_DEM_init
  end interface
   
contains

  function uniform_DEM_XYZ(filename) result(res)
    type(uniform_DEM) :: res
    character(*), intent(in) :: filename
    real(knd),allocatable :: points(:,:) !coordinate, index
    logical :: cmajor, xincreasing, yincreasing
    integer :: startx, endx, stepx, starty, endy, stepy
    integer :: unit, io
    integer :: i,j,k, current, top
    integer :: nx, ny
    real(knd), parameter :: eps = 1E-3
    
    top = 1000
    current = 0
    allocate(points(3,top))
       
    open(newunit=unit, file=filename, status='old', action='read', iostat=io)
    if (io/=0) then
      write(*,'(*(g0))') "Could not open file ",filename
      stop
    else
      do
        read(unit,*,iostat=io) points(:,current+1)
        if (io/=0) exit
        current = current + 1
        if (current==top) call realloc    
      end do
      close(unit)
    end if

    call trim_points
    
    res%x1 = minval(points(1,:))
    res%y1 = minval(points(2,:))

    !get ordering
    cmajor = (abs(points(1,2)-points(1,1))>eps)

    if (cmajor) then
      xincreasing = points(1,2)>points(1,1)
      res%dx = abs(points(1,2)-points(1,1))
      i = 1
      do
        i = i + 1
        if (abs(points(2,i)-points(2,1))>eps) then
          res%dy = abs(points(2,i)-points(2,1))
          nx = i - 1
          
          if (points(2,i) > points(2,1)) then
            yincreasing = .true.
          else
            yincreasing = .false.
          end if
          
          exit
        end if
      end do
      ny = size(points,2) / nx
      
      if (nx*ny/=size(points,2)) stop "Error decomposing uniform DEM terrain grid."
    else
      yincreasing = points(2,2)>points(2,1)
      res%dy = abs(points(2,2)-points(2,1))
      i = 1
      do
        i = i + 1
        if (abs(points(1,i)-points(1,1))>eps) then
          res%dx = abs(points(1,i)-points(1,1))
          ny = i - 1
            
          if (points(1,i) > points(1,1)) then
            xincreasing = .true.
          else
            xincreasing = .false.
          end if

          exit
        end if
      end do
      nx = size(points,2) / ny
      
      if (nx*ny/=size(points,2)) stop "Error decomposing uniform DEM terrain grid."
    end if

    allocate(res%z(nx,ny))
      
    if (xincreasing) then
      startx = 1
      endx = nx
      stepx = 1
    else
      startx = nx
      endx = 1
      stepx = -1
    end if
    if (yincreasing) then
      starty = 1
      endy = ny
      stepy = 1
    else
      starty = ny
      endy = 1
      stepy = -1
    end if
    
    if (cmajor) then
      k = 0
      do j = starty, endy, stepy
        do i = startx, endx, stepx
          k = k + 1
          res%z(i,j) = points(3,k)
        end do
      end do
    else
      k = 0
      do i = startx, endx, stepx
        do j = starty, endy, stepy
          k = k + 1
          res%z(i,j) = points(3,k)
        end do
      end do
    end if
    
    res%xn = res%x1 + (size(res%z,1)-1) * res%dx
    res%yn = res%y1 + (size(res%z,2)-1) * res%dy

  contains
    subroutine realloc
      real(knd),allocatable :: tmp(:,:) !coordinate, index
      call move_alloc(points, tmp)
      top = top * 2
      allocate(points(3,top))
      points(:,:size(tmp,2)) = tmp
    end subroutine
    subroutine trim_points
      real(knd),allocatable :: tmp(:,:) !coordinate, index
      call move_alloc(points, tmp)
      allocate(points(3,current))
      points = tmp(:,:current)
    end subroutine
  end function
  
  
  
  
  function uniform_DEM_ELV(filename) result(res)
    type(uniform_DEM) :: res
    character(*), intent(in) :: filename
    real(real32),allocatable :: z(:,:)
    integer(int32) :: i
    real(real32) :: r
    integer :: unit, io
    integer :: nx, ny

  
    !orthogonal grid of elevations in a custom format .elv
    !unformatted stream file 32-bit kinds, header: 1, nx, x1, dx, ny, y1, dy
    !body: whole array z in Fortran native order
  
  
    open(newunit=unit, file=filename, access="stream", form="unformatted", status="old", action="read", iostat=io)
    if (io/=0) then
      write(*,'(*(g0,1x))') "Error opening file ",filename
      stop
    end if
    read(unit) i
    if (i/=1) then
      !TODO: implement automatic conversion using module Endianness
      write(*,'(*(g0,1x))') "Error reading file ",filename,", wrong endianness?"
    end if
    read(unit) i
    nx = i
    read(unit) r
    res%x1 = r
    read(unit) r
    res%dx = r
    read(unit) i
    ny = i
    read(unit) r
    res%y1 = r
    read(unit) r
    res%dy = r
    allocate(z(nx,ny))
    read(unit) z
    close(unit)

    res%xn = res%x1 + (nx-1) * res%dx
    res%yn = res%y1 + (ny-1) * res%dy
    !NOTE: must be here because of http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58861 ?
    allocate(res%z(nx,ny))
    res%z = real(z, kind=knd)
  end function
  
  
  
  function uniform_DEM_init(filename) result(res)
    type(uniform_DEM) :: res
    character(*), intent(in) :: filename
    
    if (filename(len(filename)-2:) == "xyz") then
      res = uniform_DEM_XYZ(filename)
    else if (filename(len(filename)-2:) == "elv") then
      res = uniform_DEM_ELV(filename)
    else
      write(*,'(*(g0,1x))')  "Error, unknown file extension in file ",filename
      stop
    end if
  end function
  
  
  
  
  
  function uniform_DEM_elevation(self, x, y) result(res)
    real(knd) :: res
    class(uniform_DEM), intent(in) :: self
    real(knd), intent(in) :: x, y
    real(knd) :: x1, x2, y1, y2
    integer :: i1, i2, j1, j2
    integer :: nx, ny
  
    nx = size(self%z,1)
    ny = size(self%z,2)
  
    i1 = int((x-self%x1)/self%dx) + 1
    if (i1>=nx) then
      i2 = nx
      i1 = nx - 1
    else if (i1<1) then
      i1 = 1
      i2 = 2
    else
      i2 = i1 + 1
    end if
    
    j1 = int((y-self%y1)/self%dy) + 1
    if (j1>=ny) then
      j2 = ny
      j1 = ny - 1
    else if (j1<1) then
      j1 = 1
      j2 = 2
    else
      j2 = j1 + 1
    end if
    
    x1 = self%x1 + (i1-1)*self%dx
    x2 = self%x1 + (i2-1)*self%dx
    y1 = self%y1 + (j1-1)*self%dy
    y2 = self%y1 + (j2-1)*self%dy

    res = ( self%z(i1,j1) * (x2 - x) * (y2 - y) + &
            self%z(i2,j1) * (x - x1) * (y2 - y) + &
            self%z(i1,j2) * (x2 - x) * (y - y1) + &
            self%z(i2,j2) * (x - x1) * (y - y1) ) / &
          ((x2 - x1) * (y2 - y1))
  end function
  
end module ElevationModels

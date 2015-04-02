module Custom_gabls4
  use Kinds
  use Interpolation
  
  implicit none
  
  logical :: initialized = .false.
  
  integer :: last = 1
  real(knd) :: surf_temp
  
  real(knd), allocatable :: coefs(:,:), table(:,:)
  
contains

  subroutine read_table(fname)
    character(*), intent(in) :: fname
    integer :: io, u, n, i
    real(knd) :: a
    character(80) :: line
    
    open(newunit=u, file=fname, status="old", action="read")
    read(u,*)
    n = 0
    do
      read(u,'(a)',iostat=io) line
      if (io/=0) exit
      n = n + 1
    end do

    allocate(table(n,2))

    rewind(u)
    read(u,*)
    
    do i = 1, n
      read(u,'(a)',iostat=io) line
      read(line, *) table(i,:)
    end do

    close(u)

  end subroutine
    
  subroutine initialize_surface_temp
  
    call read_table("input/Tg.txt")

    allocate(coefs(0:3,size(table,1)))
 
    call linear_interpolation(table(:,1), table(:,2), coefs)
    
    initialized = .true.
  
  end subroutine

end module


subroutine CustomTimeStepProcedure
  use Custom_gabls4
  use Parameters

  surf_temp = cubic_spline_eval(time, table(:,1), coefs, last)
end subroutine


function CustomSurfaceTemperature(x,y,z,t) result(res)
   use Kinds
   use Custom_gabls4
   
   implicit none
   
   real(knd) :: res
   real(knd), intent(in) :: x, y, z
   real(tim), intent(in) :: t
   
   res = surf_temp

end function

subroutine CustomConfiguration_first
end subroutine

subroutine CustomConfiguration_last
  use Kinds
  use Scalars
  use Custom_gabls4
  
  implicit none
  
  integer :: n, u, i, io
  real(knd), allocatable :: theta(:,:)
  character(80) :: line
  
  call initialize_surface_temp
  
  open(newunit=u, file="input/theta_profile.txt", status="old", action="read")
  n = 0
  do
    read(u,'(a)',iostat=io) line
    if (io/=0) exit
    n = n + 1
  end do

  allocate(theta(n,2))

  rewind(u)
  
  do i = 1, n
    read(u,'(a)',iostat=io) line
    read(line, *) theta(i,:)
  end do

  close(u)
  
  call move_alloc(theta, TemperatureProfile%points)
end subroutine

subroutine CustomBoundaryConditions
  use Wallmodels
  
  implicit none
  
  integer :: i, j
  
  WMPoints%z0 = 0.01
  WMPoints%z0H = 0.001
  
  do j = 1, 3
    do i = 1, 6
      WMPointsUVW(i,j)%points%z0 = 0.01
      WMPointsUVW(i,j)%points%z0H = 0.001
    end do
  end do
end subroutine

subroutine  CustomSolidBodies
end subroutine  CustomSolidBodies

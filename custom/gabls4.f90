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
    
  subroutine initialize
  
    call read_table("input/Tg.txt")

    allocate(coefs(0:3,size(table,1)))
 
    call cubic_spline(table(:,1), table(:,2), coefs)
    
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


subroutine  CustomSolidBodies
end subroutine  CustomSolidBodies

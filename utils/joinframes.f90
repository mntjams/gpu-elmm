program joinvtkframes

  character(:), allocatable :: domain, filename, cmd
  character(80) :: arg
  integer :: i, stat
  
  if (command_argument_count()<2) then
    write(*,*) "Usage: joinvtkframes label npx,npy,npz"
    stop 1
  end if
 
  call get_command_argument(1, value=arg)
  
  domain = trim(arg)

  call get_command_argument(2, value=arg)

  i = 0
  do
    filename = 'frame-'//domain//'-'//itoa(i)//'.vtk'
    cmd = 'joinvtk '//filename//' '//arg

    call execute_command_line(cmd, exitstat=stat)

    if (stat/=0) exit
    i = i + 1
  end do
    
contains
  
  function itoa(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
  end function
  
end program

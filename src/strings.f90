module Strings
  !Auxiliary string processing procedures.

  implicit none

  private

  public upcase, count_multispaces, itoa

  contains

    function upcase(string) result(res)
      character(*),intent(in) :: string
      character(len(string))  :: res
      integer :: j

      do j = 1,len(string)
        if(string(j:j) >= "a" .and. string(j:j) <= "z") then
             res(j:j) = achar(iachar(string(j:j)) - 32)
        else
             res(j:j) = string(j:j)
        end if
      end do
    end function upcase

    function count_multispaces(string) result(res)
      !count groups of spaces between groups of other characters
      integer :: res
      character(*),intent(in) :: string
      character(len(string))  :: loc_string
      integer :: i
      logical :: lastspace

      loc_string = adjustl(string)

      lastspace = .false.
      res = 0
      do i=1,len_trim(loc_string)
        if (loc_string(i:i) == ' ') then
          if (.not.lastspace) res = res + 1
          lastspace = .true.
        else
          lastspace = .false.
        end if
      end do

    end function count_multispaces
    
    function itoa(i) result(res)
      character(:),allocatable :: res
      integer,intent(in) :: i
      character(range(i)+2) :: tmp
      write(tmp,'(i0)') i
      res = trim(tmp)
    end function

end module Strings

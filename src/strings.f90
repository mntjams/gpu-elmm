module Strings
  !Auxiliary string processing procedures.

  implicit none

  private

  public :: upcase, downcase, count_multispaces, itoa, tokenize

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

  function downcase(string) result(res)
    character(*),intent(in) :: string
    character(len(string))  :: res
    integer :: j

    do j = 1,len(string)
      if(string(j:j) >= "A" .and. string(j:j) <= "Z") then
           res(j:j) = achar(iachar(string(j:j)) + 32)
      else
           res(j:j) = string(j:j)
      end if
    end do
  end function downcase

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


  subroutine tokenize(stream,res,ierr)
    use str_lists, only: str_list => list, char_len
    type(str_list),intent(out) :: res
    character(*),intent(in),target :: stream
    character,pointer :: s
    integer, intent(out) :: ierr
    character(char_len) :: token
    integer :: i,pos
    logical :: in_string

    in_string = .false.
    pos = 1
    token = ''
    
    do i = 1, len_trim(stream)

      s => stream(i:i)

      if (in_string) then
       
        if (s=='"'.and.stream(i-1:i-1)/='\') then
          token(pos:pos) = s
          call send
          in_string = .false.
        else
          token(pos:pos) = s
          pos = pos + 1
        end if
        
      else
       
        if (s=='('.or.s==')') then
          if (pos>1) call res%add(token)
          token = s
          call send
        else if (s=='"') then
          if (pos>1) then
            ierr = 1
            return
          else
            token(1:1) = s
            pos = pos + 1
            in_string = .true.
          end if
        else if (s==' ') then
          if (pos>1) call send
        else if (s==';') then
          if (pos>1) call send
          exit
        else
          token(pos:pos) = s
          pos = pos + 1
        end if
        
      end if

    end do

    if (pos>1) call send

    ierr = 0
    
  contains

    subroutine send
      call res%add(token)
      token = ''
      pos = 1
    end subroutine

  end subroutine tokenize

end module Strings



module ParseTrees

  use Strings

  use str_lists, only: str_list => list, char_len

  implicit none

  private

  public tree_object, tree_object_field, tree_object_fields, parse_file

  type tree_object_field
    character(char_len) :: name
    logical :: is_object = .false.
    character(char_len) :: value
    type(tree_object), pointer :: object_value => null()
  end type

  type tree_object_fields
    type(tree_object_field), allocatable :: array(:)
  end type

  type tree_object
    character(char_len) :: name
    type(tree_object_fields) :: fields
  contains
    procedure :: finalize => tree_object_finalize
  end type

contains

  recursive subroutine tree_object_finalize(obj)
    class(tree_object) :: obj
    integer :: i

    if (allocated(obj%fields%array)) then
      do i = 1, size(obj%fields%array)
        associate(field => obj%fields%array(i))
          if (associated(field%object_value)) then
            call field%object_value%finalize
            deallocate(field%object_value)
          end if
        end associate
      end do
      deallocate(obj%fields%array)
    end if
  end subroutine

  subroutine to_array(ch_list, ch_array)
    type(str_list) :: ch_list
    character(char_len), allocatable :: ch_array(:)
    allocate(ch_array(0))
    call ch_list%for_each(aux)
  contains
    subroutine aux(item)
      character(char_len) :: item

      ch_array = [ch_array, item]
    end subroutine
  end subroutine

  subroutine parse_file(tree, fname, stat)
    use Strings, only: itoa, downcase
    type(tree_object), allocatable, intent(out) :: tree(:)
    character(*), intent(in) :: fname
    integer, intent(out) :: stat

    character(:), allocatable :: stream
    type(tree_object), allocatable :: tmp(:)
    type(str_list) :: token_list
    character(char_len), allocatable :: tokens(:)
    type(tree_object) :: object

    integer :: pos, io, unit
    character(256) :: line


    open(newunit=unit, file=fname,status="old",action="read",iostat=stat)
       
    if (stat/=0) return

    stream = ''

    do
      read(unit,'(a)',iostat=io) line
      if (io/=0) exit
      stream = stream // ' ' // trim(adjustl(line))
    end do

    call tokenize(stream, token_list, stat)

    if (stat/=0) then
      call token_list%finalize
      return
    end if

    call to_array(token_list, tokens)
    call token_list%finalize

    allocate(tree(0))
    pos = 1
    do
      call get_object(object, pos, stat)
      if (stat /=0) exit

      call move_alloc(tree, tmp)
      allocate(tree(size(tmp)+1))
      tree(:size(tmp)) = tmp

      tree(size(tmp)+1) = object
      call object%finalize

      if (pos > size(tokens)) exit
    end do

  contains

    recursive subroutine get_object(object, pos, stat)
      type(tree_object) :: object
      integer, intent(inout) :: pos
      integer, intent(out) :: stat

      if (pos+2 > size(tokens)) then
        !need 3 tokens for an empty object
        stat = 2
        return
      end if

      if (tokens(pos+1)=='(') then
        call get_string(object%name, pos, stat)
        if (stat /= 0) return
        pos = pos + 1 !checked above
        call get_object_fields(object%fields, pos, stat)
        if (stat /= 0) return
        call check_string(')', pos, stat)
      end if
    end subroutine

    subroutine check_string(str, pos, stat)
      character(*), intent(in) :: str
      integer, intent(inout) :: pos
      integer, intent(out) :: stat

      if (downcase(str)==downcase(tokens(pos))) then
        stat = 0
        pos = pos + 1
      else
        write(*,*) "Error in '" // &
                   fname // &
                   "', expected '" // &
                   str // &
                   "' read '" // &
                   downcase(tokens(pos)) // &
                   "' instead."
        stat = 1
      end if
    end subroutine

    recursive subroutine get_object_fields(fields, pos, stat)
      type(tree_object_fields) :: fields
      integer, intent(inout) :: pos
      integer, intent(out) :: stat
      type(tree_object_field), allocatable :: tmp(:)

      allocate(fields%array(0))
      do
        if (tokens(pos) == ')') return

        call move_alloc(fields%array, tmp)
        allocate(fields%array(size(tmp)+1))
        fields%array(:size(tmp)) = tmp

        call get_field(fields%array(size(tmp)+1), pos, stat)
        if (stat /= 0) return
      end do
    end subroutine

    recursive subroutine get_field(field, pos, stat)
      type(tree_object_field) :: field
      integer, intent(inout) :: pos
      integer, intent(out) :: stat

      call get_string(field%name, pos, stat)
      call check_string('=', pos, stat)

      if (tokens(pos+1) == '(') then
        field%is_object = .true.
        allocate(field%object_value)
        call get_object(field%object_value, pos, stat)
      else
        call get_string(field%value, pos, stat)
      end if
    end subroutine

    recursive subroutine get_string(str, pos, stat)
      character(char_len) :: str
      integer, intent(inout) :: pos
      integer, intent(out) :: stat
      if (pos>size(tokens)) then
        stat = 1
        return
      else
        str = tokens(pos)
        pos = pos + 1
      end if
    end subroutine

  end subroutine parse_file

end module ParseTrees

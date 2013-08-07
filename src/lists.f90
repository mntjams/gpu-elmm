module Lists

  implicit none

  type,abstract :: Listable  !A placeholder, unlimited polymorphism still problematic in compilers as of January 2013.
  end type

  type ListNode
    class(Listable),allocatable :: item
    type(ListNode),pointer :: next =>null()
  end type ListNode

  type List
    type(ListNode),pointer,private :: first => null()
    type(ListNode),pointer,private :: last => null()
    type(ListNode),pointer,private :: iter => null()
    integer                ,private :: length = 0
    contains
      procedure :: Deallocate => List_Deallocate
      procedure :: Add => List_Add
      procedure :: IterNext => List_IterNext
      procedure :: IterRestart => List__IterRestart
      procedure :: ForEach => List_ForEach
      procedure :: All => List_All
      procedure :: Any => List_Any
      procedure :: Len => List_Len
  end type List

  abstract interface
    subroutine foreach_sub(item)
      import
      class(Listable) :: item
    end subroutine
    logical function elem_logical_fun(item)
      import
      class(Listable) :: item
    end function
  end interface


  contains


    subroutine List_Deallocate(self)
      class(List),intent(inout) :: self
      type(ListNode),pointer :: node,tmp

      node => self%first

      do while (associated(node))
        tmp => node
        node => node%next

        deallocate(tmp)
      end do

      self%last => null()
      self%first => null()

      self%length = 0
      

    end subroutine List_Deallocate



    subroutine List_Add(self,item)
      use iso_c_binding
      class(List),intent(inout) :: self
      class(Listable),intent(in) :: item

      if (.not.associated(self%last)) then
        allocate(self%first)
        self%last => self%first
      else
        allocate(self%last%next)
        self%last => self%last%next
      endif

      allocate(self%last%item, source=item)

      self%length = self%length + 1

    end subroutine List_Add


    subroutine List__IterRestart(self)
      class(List),intent(inout) :: self

      self%iter => self%first

    end subroutine List__IterRestart


    subroutine List_IterNext(self,res)
      class(List),intent(inout) :: self
      class(Listable),pointer,intent(out) :: res

      if (associated(self%iter)) then
        res => self%iter%item
        self%iter => self%iter%next
      else
        res => null()
      end if
    end subroutine List_IterNext

    subroutine List_ForEach(self,proc)
      class(List),intent(inout) :: self
      procedure(foreach_sub) :: proc
      type(ListNode),pointer :: node

      node => self%first

      do while (associated(node))
        if (allocated(node%item)) call proc(node%item)
        node => node%next
      end do

    end subroutine

    logical function List_All(self,proc) result(res)
      class(List),intent(inout) :: self
      procedure(elem_logical_fun) :: proc
      type(ListNode),pointer :: node

      res = .true.
      
      node => self%first

      do while (associated(node))
        if (allocated(node%item)) then
          res =  proc(node%item)
          if (.not.res) return
        end if
        node => node%next
      end do

    end function

    logical function List_Any(self,proc) result(res)
      class(List),intent(inout) :: self
      procedure(elem_logical_fun) :: proc
      type(ListNode),pointer :: node

      res = .false.
      
      node => self%first

      do while (associated(node))
        if (allocated(node%item)) then
          res =  proc(node%item)
          if (res) return
        end if
        node => node%next
      end do

    end function

    pure integer function List_Len(self)
      class(List),intent(in) :: self

      List_Len = self%length
    end function
    
end module Lists




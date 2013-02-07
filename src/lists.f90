module Lists

  implicit none

  type,abstract :: TListable  !A placeholder, unlimited polymorphism still problematic in compilers as of January 2013.
  end type

  type TListNode
    class(TListable),allocatable :: item
    type(TListNode),pointer :: next =>null()
  end type TListNode

  type TList
    type(TListNode),pointer,private :: first => null()
    type(TListNode),pointer,private :: last => null()
    type(TListNode),pointer,private :: iter => null()
    integer                ,private :: length = 0
    contains
      procedure :: Deallocate => TList_Deallocate
      procedure :: Add => TList_Add
      procedure :: IterNext => TList_IterNext
      procedure :: IterRestart => TList__IterRestart
      procedure :: ForEach => TList_ForEach
      procedure :: Len => TList_Len
  end type TList

  abstract interface
    subroutine foreach_sub(item)
      import
      class(TListable) :: item
    end subroutine
  end interface


  contains


    subroutine TList_Deallocate(self)
      class(Tlist),intent(inout) :: self
      type(TListNode),pointer :: node,tmp

      node => self%first

      do while (associated(node))
        tmp => node
        node => node%next

        deallocate(tmp)
      end do

      self%last => null()
      self%first => null()

      self%length = 0
      

    end subroutine TList_Deallocate



    subroutine TList_Add(self,item)
      use iso_c_binding
      class(TList),intent(inout) :: self
      class(TListable),intent(in) :: item

      if (.not.associated(self%last)) then
        allocate(self%first)
        self%last => self%first
      else
        allocate(self%last%next)
        self%last => self%last%next
      endif

      allocate(self%last%item, source=item)

      self%length = self%length + 1

    end subroutine TList_Add


    subroutine TList__IterRestart(self)
      class(TList),intent(inout) :: self

      self%iter => self%first

    end subroutine TList__IterRestart


    subroutine TList_IterNext(self,res)
      class(TList),intent(inout) :: self
      class(TListable),pointer,intent(out) :: res

      if (associated(self%iter)) then
        res => self%iter%item
        self%iter => self%iter%next
      else
        res => null()
      end if
    end subroutine TList_IterNext

    subroutine TList_ForEach(self,proc)
      class(TList),intent(inout) :: self
      procedure(foreach_sub) :: proc
      type(TListNode),pointer :: node

      node => self%first

      do while (associated(node))
        if (allocated(node%item)) call proc(node%item)
        node => node%next
      end do

    end subroutine

    pure integer function TList_Len(self)
      class(TList),intent(in) :: self

      TList_Len = self%length
    end function
    
end module Lists




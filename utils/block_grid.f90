module kinds
  integer, parameter :: rp = kind(1.d0)
end module

module grid_blocks
  use kinds

  implicit none

  type block
    real(rp) :: lb, ub
    integer :: ncells
    real(rp) :: ratio
    real(rp), allocatable :: z_u(:)
  end type

contains

  function block_uniform(lb, ub, ncells) result(res)
    type(block) :: res
    real(rp), intent(in) :: lb, ub
    integer, intent(in) :: ncells
    real(rp) :: d
    integer :: i
    res%ub = ub
    res%lb = lb
    
    res%ratio = 1
    res%ncells = ncells
    
    d = (ub - lb) / ncells
    res%z_u = [(d*i, i = 1, ncells-1)]

  end function

  function block_ratio_ncells(lb, ub, ratio, ncells) result(res)
    type(block) :: res
    real(rp), intent(in) :: lb, ub, ratio
    integer, intent(in) :: ncells
    real(rp) :: start
    
    start = (ub - lb) / gsum(1._rp, ratio, ncells)
    
    res%ub = ub
    res%lb = lb
    res%ratio = ratio
    res%ncells = ncells
    res%z_u = lb + divisions(start, ratio, ncells)
  end function

  function block_start_end(lb, ub, start, end, start_priority) result(res)
    type(block) :: res
    real(rp), intent(in) :: lb, ub, start, end
    logical, optional :: start_priority
    logical :: s_p
    integer :: n
    real(rp) :: l, r
    
    if (present(start_priority)) then
      s_p = start_priority
    else
      s_p = .true.
    end if
    
    res%ub = ub
    res%lb = lb
    
    l = ub - lb

    r = (l - start) / (l - end)
    n = nint(1 + log((end/start)) / log(r))

    if (s_p) then
      r = gratio(start, l, n, r)
      res%z_u = lb + divisions(start, r, n)
    else
      r = 1/gratio(end, l, n, 1/r)
      res%z_u = ub - divisions(end, r, n)
      res%z_u = res%z_u(size(res%z_u):1:-1)
    end if
      
    res%ratio = r
    res%ncells = n
  end function

  function block_start_ratio(lb, ub, start, ratio) result(res)
    type(block) :: res
    real(rp), intent(in) :: lb, ub, start, ratio
    real(rp) :: l, r1, r2, r, gs
    integer :: i, n
    integer, parameter :: maxn = 200
    
    res%ub = ub
    res%lb = lb
    
    l = ub - lb
    
    if (gsum(start, ratio, maxn) < l) then
     r = gratio(start, l, maxn, ratio)
    else
     r = ratio
    end if
  
    i = 1
    do
      gs = gsum(start, r, i)
      if (gs>l) then
        r1 = gratio(start, l, i, ratio)
        r2 = gratio(start, l, i-1, ratio)
         
        if (abs(ratio-r1)<=abs(ratio-r2)) then
          r = r1
          n  = i
        else
          r = r2
          n = i - 1
        end if
        
        exit
      end if
      i = i + 1
    end do
    
    res%ratio = r
    res%ncells = n
    res%z_u = lb + divisions(start, r, n)
  end function

  function block_end_ratio(lb, ub, end, ratio) result(res)
    type(block) :: res
    real(rp), intent(in) :: lb, ub, end, ratio
    
    res = block_start_ratio(lb, ub, end, ratio)
    res%z_u = ub - divisions(end, res%ratio, res%ncells)
    res%z_u = res%z_u(size(res%z_u):1:-1)
  end function
  
  function gsum(a, r, n)
    real(rp) :: gsum
    real(rp), intent(in) :: a, r
    integer, intent(in) :: n
    gsum = a * (1 - r**n) / (1 - r)
  end function
  
  function gratio(a, s, n, r0) result(r)
    real(rp) :: r
    real(rp), intent(in) :: a, s
    real(rp), value :: r0
    integer, intent(in) :: n
    real(rp), parameter :: eps = epsilon(r0)
    integer :: i
    ! fixed point iteration
    if (a*n<=s) then
      if (r0<1) r0 = 1.05
      do i = 1, 100
        r = f(r0)       
        if (abs(r-r0)<=eps) exit
        r0 = r
      end do
    else
      if (r0>1.or.r0<=0) r0 = 0.95
      do i = 1, 100
        r = g(r0)      
        if (abs(r-r0)<=eps) exit
        r0 = r
      end do
    end if
  contains
    function f(r)
      real(rp) :: f, r
      f = (1-(s/a)*(1-r))**(1._rp/n)
    end function
    function g(r)
      real(rp) :: g, r
      g = (r**n-1)/(s/a) + 1    
    end function
  end function
  
  function divisions(a, r, n)
    real(rp) :: divisions(1:n-1)
    real(rp), intent(in) :: a, r
    integer, intent(in) :: n
    integer :: i
    
    divisions(1) = a
    do i = 2, n-1
      divisions(i) = divisions(i-1) + a * r**(i-1)
    end do
  end function
end module grid_blocks

program block_grid
  use grid_blocks
  
  implicit none
  
  type(block), allocatable :: blocks(:)
  integer :: i, iblock
  
  allocate(blocks(1))
  
  blocks(1) = block_ratio_ncells(0._rp, 1._rp, 0.95_rp, 20)
  
  do iblock = 1, size(blocks)
    write(*,*) blocks(iblock)%lb
    do i=1, size(blocks(iblock)%z_u) 
      write(*,*) blocks(iblock)%z_u(i)
    end do
  end do
  write(*,*) blocks(size(blocks))%ub
end program 

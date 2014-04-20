module Frames_common

  use Kinds
  use iso_c_binding, only: c_ptr
  use Pthreads

  type TFrameTimes
    real(knd) :: start, end
    integer :: nframes
  end type

  type, abstract :: TFrameBase

    type(TFrameTimes) :: frame_times

    integer   :: frame_number = -1 !number of the previous frame

    real(knd), allocatable :: times(:)

    integer   :: minPri, maxPri, minPrj, maxPrj, minPrk, maxPrk

    integer   :: sizePr

    real(knd),allocatable,dimension(:) :: xPr,yPr,zPr

    logical :: ranges_set = .false.

    character(40) :: base_name

    integer :: unit

    logical :: in_progress = .false.

    type(c_ptr) :: threadptr
  contains
    procedure(save_interface), deferred :: save
    procedure(fill_interface), deferred :: fill
    procedure :: SaveTimes => TFrameBase_SaveTimes
    procedure :: Wait => TFrameBase_Wait
    procedure :: Finalize => TFrameBase_Finalize
  end type
  
  abstract interface
    subroutine fill_interface(D, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
      import
      class(TFrameBase),intent(inout) :: D
      real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
      real(knd),contiguous,intent(in) :: Pr(1:,1:,1:), Viscosity(-1:,-1:,-1:), &
                              Temperature(-1:,-1:,-1:), Moisture(-1:,-1:,-1:), &
                              Scalar(-1:,-1:,-1:,1:)
    end subroutine
    subroutine save_interface(D, time, U, V, W, Pr, Viscosity, Temperature, Moisture, Scalar)
      import
      class(TFrameBase),target,asynchronous,intent(inout) :: D
      real(knd),intent(in) :: time
      real(knd),dimension(-2:,-2:,-2:),contiguous,intent(in) :: U,V,W
      real(knd),contiguous,intent(in) :: Pr(1:,1:,1:), Viscosity(-1:,-1:,-1:), &
                                         Temperature(-1:,-1:,-1:), Moisture(-1:,-1:,-1:), &
                                         Scalar(-1:,-1:,-1:,1:)
    end subroutine
  end interface

  interface Add
    module procedure TimeSeries_Add
  end interface
  
contains

  subroutine TimeSeries_Add(array,idx,val)
    real(knd),allocatable,intent(inout) :: array(:)
    integer,intent(in)   :: idx
    real(knd),intent(in) :: val
    real(knd),allocatable :: tmp(:)

    !assume allocated, otherwise error by the runtime library

    if (idx < lbound(array,1)) then
      write(*,*) "Error. Tried to write under the start of the series."
      stop
    end if
    if (idx > ubound(array,1)) then
      allocate( tmp(lbound(array,1):lbound(array,1)+(idx-lbound(array,1))*2) )
      tmp(lbound(array,1):ubound(array,1)) = array
      deallocate(array)
      call move_alloc(tmp, array)
    end if
    
    array(idx) = val

  end subroutine TimeSeries_Add
  
  
  subroutine TFrameBase_SaveTimes(D)
    class(TFrameBase),intent(in) :: D
    character(40) :: file_name
    integer i

    file_name = trim(D%base_name)//"-times.txt"

    open(unit=D%unit, file=file_name, access='sequential', form='formatted', action='write', status='replace')

    do i = 0, D%frame_number
      write(D%unit,'(i8,2x,es12.5)') i, D%times(i)
    end do
    close(D%unit)

  end subroutine
  
  
  subroutine TFrameBase_Wait(D)
    class(TFrameBase),asynchronous,intent(inout) :: D
    integer err

    call pthread_join_opaque(D%threadptr,err)

    if (err/=0) then
      write (*,*) "Error in joining staggered frame thread. Will try to continue anyway. Code:",err
    end if

    D%in_progress = .false.
  end subroutine



  subroutine TFrameBase_Finalize(D)
    class(TFrameBase),intent(inout) :: D

    if (D%in_progress) call D%Wait

    call D%SaveTimes

  end subroutine


end module
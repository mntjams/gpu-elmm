module Pthreads
  use iso_c_binding

  implicit none

  interface
    subroutine pthread_create_opaque(threadptr, procptr, dataptr, err) bind(C,name="pthread_create_opaque")
      use iso_c_binding
      type(c_ptr) :: threadptr
      type(c_funptr),value :: procptr
      type(c_ptr),value :: dataptr
      integer(c_int),intent(out) :: err
    end subroutine

    subroutine pthread_join_opaque(thread, err) bind(C,name="pthread_join_opaque")
      use iso_c_binding
      type(c_ptr),value :: thread
      integer(c_int),intent(out) :: err
    end subroutine
  end interface
end module Pthreads


module StaggeredFrames
  use iso_c_binding, only: c_ptr
  use Kinds
  use Pthreads

  implicit none

  private
  public i3, r3, irange, rrange, TFrameTimes, TSaveFlags, TStaggeredFrameDomain, &
         Init, Save, Finalize


  type i3
    integer i,j,k
  end type

  type r3
    real(knd) x,y,z
  end type

  type irange
    type(i3) min,max
  end type

  type rrange
    type(r3) min,max
  end type

  type TSaveFlags
    logical :: U = .false.
    logical :: V = .false.
    logical :: W = .false.
    logical :: Pr = .false.
    logical :: Temperature = .false.
    logical :: Moisture = .false.
    logical :: Scalar = .false.
    integer :: num_scalars = 0
  end type

  type TFrameTimes
    real(knd) :: start, end
    integer :: nframes
  end type

  type :: TStaggeredFrameDomain
    !domain for saving full timed data
    !fluxes should be computed by postprocessing, preferably using the same numerical methods
    !as in the primary solver
    private

    type(TFrameTimes) :: frame_times

    integer   :: frame_number = -1 !number of the previous frame

    real(knd), allocatable :: times(:)

    type(rrange) :: range
    integer   :: minUi, maxUi, minUj, maxUj, minUk, maxUk
    integer   :: minVi, maxVi, minVj, maxVj, minVk, maxVk
    integer   :: minWi, maxWi, minWj, maxWj, minWk, maxWk
    integer   :: minPri, maxPri, minPrj, maxPrj, minPrk, maxPrk
    integer   :: sizeU, sizeV, sizeW, sizePr
    real(knd),allocatable,dimension(:) :: xPr,yPr,zPr,xU,yV,zW
      !temperature and scalars use Prsize and min/maxPr_
    logical :: ranges_set = .false.

    type(TSaveFlags) :: save_flags

    integer :: buffer_size = 0
    real(knd),allocatable :: buffer(:)

    character(20) :: base_name
    character(4)  :: suffix = ".unf"
    integer       :: ndigits

    integer :: unit

    logical :: in_progress = .false.

    type(c_ptr) :: threadptr

!     contains
    !  procedure :: Init => TStaggeredFrameDomain_Init
    !  procedure :: Fill => TStaggeredFrameDomain_Fill
    !  procedure :: Save => TStaggeredFrameDomain_Save
    !  procedure :: Wait => TStaggeredFrameDomain_Wait
    !  procedure :: Finalize => TStaggeredFrameDomain_Finalize
  endtype TStaggeredFrameDomain

  interface Init
    module procedure TStaggeredFrameDomain_Init
  end interface
  interface SetRest
    module procedure TStaggeredFrameDomain_SetRest
  end interface
  interface Fill
    module procedure TStaggeredFrameDomain_Fill
  end interface
  interface Wait
    module procedure TStaggeredFrameDomain_Wait
  end interface
  interface Save
    module procedure TStaggeredFrameDomain_Save
  end interface
  interface SaveTimes
    module procedure TStaggeredFrameDomain_SaveTimes
  end interface
  interface SaveHeader
    module procedure TStaggeredFrameDomain_SaveHeader
  end interface
  interface Finalize
    module procedure TStaggeredFrameDomain_Finalize
  end interface

  interface Add
    module procedure TimeSeries_Add
  end interface

  interface assignment(=)
    module procedure assign3to1
  end interface

  integer :: stagframe_unit = 200

  contains



    subroutine TimeSeries_Add(array,idx,val)
      real(knd),allocatable,intent(inout) :: array(:)
      integer,intent(in)   :: idx
      real(knd),intent(in) :: val
      real(knd),allocatable :: tmp(:)

      !assume allocated, otherwise error by the runtime library

      if (idx < lbound(array,1)) then
        write(*,*) "Error. Tried to write under the start of the series. Ignoring request."
        return
      end if
      if (idx > ubound(array,1)) then
        allocate( tmp(lbound(array,1):lbound(array,1)+(idx-lbound(array,1))*2) )
      end if

      array(idx) = val

    end subroutine TimeSeries_Add



    subroutine assign3to1(A1,A3)
      real(knd),intent(out) :: A1(1:)
      real(knd),intent(in)  :: A3(1:,1:,1:)
      integer j,k
      integer s1,s2,s3 !sizes
      integer off2,off3 !offsets

      !if s1*s2*s3 /= size(A1) let it fail by the Fortran runtime library

      s1 = size(A3,1)
      s2 = size(A3,2)
      s3 = size(A3,3)

      off3 = 0
      do k=1,size(A3,3)
        off2 = 0
        do j=1,size(A3,2)
          A1(off2+off3+1:off2+off3+s1) = A3(:,j,k)
          off2 = off2 + s1
        end do
        off3 = off3 + s1*s2
      end do
    end subroutine assign3to1



    subroutine TStaggeredFrameDomain_Init(D,label,range,time_params,save_flags)
      type(TStaggeredFrameDomain),intent(out) :: D
      character(*) :: label
      type(rrange),intent(in) :: range
      type(TFrameTimes),intent(in) :: time_params
      type(TSaveFlags),intent(in) :: save_flags

      D%base_name = "stagframe-"//label

      D%frame_times = time_params

      D%range = range

      D%save_flags = save_flags

      D%unit = stagframe_unit + 1
      stagframe_unit = D%unit

    end subroutine TStaggeredFrameDomain_Init


    subroutine TStaggeredFrameDomain_SetRest(D, num_scalars)

      use Boundaries, only: GridCoords
      use Parameters, only: xU, yV, zW, xPr, yPr, zPr

      type(TStaggeredFrameDomain),intent(inout) :: D
      integer,intent(in) :: num_scalars

      call GridCoords(D%minPri, D%minPrj, D%minPrk, D%range%min%x, D%range%min%y, D%range%min%z)

      call GridCoords(D%maxPri, D%maxPrj, D%maxPrk, D%range%max%x, D%range%max%y, D%range%max%z)

      D%minUi = D%minPri - 1
      D%minUj = D%minPrj
      D%minUk = D%minPrk

      D%minVi = D%minPri
      D%minVj = D%minPrj - 1
      D%minVk = D%minPrk

      D%minWi = D%minPri
      D%minWj = D%minPrj
      D%minWk = D%minPrk - 1

      D%maxUi = D%maxPri
      D%maxUj = D%maxPrj
      D%maxUk = D%maxPrk

      D%maxVi = D%maxPri
      D%maxVj = D%maxPrj
      D%maxVk = D%maxPrk

      D%maxWi = D%maxPri
      D%maxWj = D%maxPrj
      D%maxWk = D%maxPrk

      D%sizeU = max( (D%maxUi - D%minUi + 1) * (D%maxUj - D%minUj + 1) * (D%maxUk - D%minUk + 1) , 0 )
      D%sizeV = max( (D%maxVi - D%minVi + 1) * (D%maxVj - D%minVj + 1) * (D%maxVk - D%minVk + 1) , 0 )
      D%sizeW = max( (D%maxWi - D%minWi + 1) * (D%maxWj - D%minWj + 1) * (D%maxWk - D%minWk + 1) , 0 )
      D%sizePr = max( (D%maxPri - D%minPri + 1) * (D%maxPrj - D%minPrj + 1) * (D%maxPrk - D%minPrk + 1) , 0 )

      allocate(D%xU(D%maxUi-D%minUi+1))  !not needed if Fortran 2003 reallocation enabled, but problems with zero sized arrays in gfortran 4.7.1
      D%xU(:) = xU(D%minUi:D%maxUi)
      allocate(D%yV(D%maxVj-D%minVj+1))
      D%yV(:) = yV(D%minVj:D%maxVj)
      allocate(D%zW(D%maxWk-D%minWk+1))
      D%zW(:) = zW(D%minWk:D%maxWk)
      allocate(D%xPr(D%maxPri-D%minPri+1))
      D%xPr(:) = xPr(D%minPri:D%maxPri)
      allocate(D%yPr(D%maxPrj-D%minPrj+1))
      D%yPr(:) = yPr(D%minPrj:D%maxPrj)
      allocate(D%zPr(D%maxPrk-D%minPrk+1))
      D%zPr(:) = zPr(D%minPrk:D%maxPrk)

      D%buffer_size = 0

      if (D%save_flags%U) D%buffer_size = D%buffer_size + D%sizeU
      if (D%save_flags%V) D%buffer_size = D%buffer_size + D%sizeV
      if (D%save_flags%W) D%buffer_size = D%buffer_size + D%sizeW
      if (D%save_flags%Pr) D%buffer_size = D%buffer_size + D%sizePr
      if (D%save_flags%Temperature) D%buffer_size = D%buffer_size + D%sizePr
      if (D%save_flags%Moisture) D%buffer_size = D%buffer_size + D%sizePr
      if (D%save_flags%Scalar) D%buffer_size = D%buffer_size + D%sizePr * num_scalars

      if (D%save_flags%Scalar) then
        D%save_flags%num_scalars = num_scalars
      else !probably redundant
        D%save_flags%num_scalars = 0
      end if

      call SaveHeader(D)

      allocate(D%buffer(0:D%buffer_size-1))

      D%ranges_set = .true.

    end subroutine TStaggeredFrameDomain_SetRest



    subroutine TStaggeredFrameDomain_Fill(D, U, V, W, Pr, Temperature, Scalar)
      !Fill the output buffer for asynchronous output
      type(TStaggeredFrameDomain),intent(inout) :: D
      real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
      real(knd),intent(in) :: Pr(1:,1:,1:), Temperature(-1:,-1:,-1:), Scalar(-1:,-1:,-1:,-1:)
      integer :: offset
      integer :: i

      if (.not.D%ranges_set) call SetRest(D, num_scalars = size(Scalar,4))

      !assume allocated by SetRest, otherwise error by the runtime library

      offset = 0

      if (D%save_flags%U) then
        D%buffer(offset:offset+D%sizeU-1) = &
            U(D%minUi:D%maxUi, D%minUj:D%maxUj, D%minUk:D%maxUk)

        offset = offset + D%sizeU
      end if

      if (D%save_flags%V) then
        D%buffer(offset:offset+D%sizeV-1) = &
            V(D%minVi:D%maxVi, D%minVj:D%maxVj, D%minVk:D%maxVk)

        offset = offset + D%sizeV
      end if

      if (D%save_flags%W) then
        D%buffer(offset:offset+D%sizeW-1) = &
            W(D%minWi:D%maxWi, D%minWj:D%maxWj, D%minWk:D%maxWk)

        offset = offset + D%sizeW
      end if

      if (D%save_flags%Pr) then
        D%buffer(offset:offset+D%sizePr-1) = &
            Pr(D%minPri:D%maxPri, D%minPrj:D%maxPrj, D%minPrk:D%maxPrk)

        offset = offset + D%sizePr
      end if

      if (D%save_flags%Temperature) then
        D%buffer(offset:offset+D%sizePr-1) = &
            Temperature(D%minPri:D%maxPri, D%minPrj:D%maxPrj, D%minPrk:D%maxPrk)

        offset = offset + D%sizePr
      end if

      if (D%save_flags%Scalar) then
        do i=1,size(Scalar,4)
          D%buffer(offset:offset+D%sizePr-1) = &
              Scalar(D%minPri:D%maxPri, D%minPrj:D%maxPrj, D%minPrk:D%maxPrk,i)

          offset = offset + D%sizePr
        end do
      end if

    end subroutine TStaggeredFrameDomain_Fill



    subroutine TStaggeredFrameDomain_Wait(D)
      type(TStaggeredFrameDomain),asynchronous,intent(inout) :: D
      integer err

      call pthread_join_opaque(D%threadptr,err)

      if (err/=0) then
        write (*,*) "Error in joining staggered frame thread. Will try to continue anyway. Code:",err
      end if

      D%in_progress = .false.
    end subroutine TStaggeredFrameDomain_Wait


    
    subroutine SaveBuffer(Dptr) bind(C)
      type(c_ptr),value :: Dptr
      type(TStaggeredFrameDomain),pointer :: D
      character(40) :: file_name

      call c_f_pointer(Dptr, D)

      write(file_name,'(a,i0,a)') trim(D%base_name)//"-",D%frame_number,trim(D%suffix)

      write(*,*) "Writing frame: ",file_name, " size: ",D%buffer_size

      open(unit=D%unit, file=file_name, access='stream', form='unformatted', action='write', status='replace')

      write(D%unit) D%times(D%frame_number)          !current time

      write(D%unit) D%buffer

      close(D%unit)

    end subroutine SaveBuffer



    subroutine TStaggeredFrameDomain_Save(D, time, U, V, W, Pr, Temperature, Scalar)
      type(TStaggeredFrameDomain),target,asynchronous,intent(inout) :: D
      real(knd),intent(in) :: time
      real(knd),dimension(-2:,-2:,-2:),intent(in) :: U,V,W
      real(knd),intent(in) :: Pr(1:,1:,1:), Temperature(-1:,-1:,-1:), Scalar(-1:,-1:,-1:,-1:)
      integer err

      associate(start   => D%frame_times%start,&
                end     => D%frame_times%end,&
                nframes => D%frame_times%nframes)

        if ( (time >= start) .and. (time <= end + (end-start)/(nframes-1)) &
          .and. (time >= start + (D%frame_number+1)*(end-start) / (nframes-1)) ) then

          D%frame_number = D%frame_number + 1

          if (.not.allocated(D%times)) allocate(D%times(D%frame_number:D%frame_number+1000))

          call Add(D%times,D%frame_number,time)

          if (D%in_progress) call Wait(D)

          call Fill(D, U, V, W, Pr, Temperature, Scalar)

          call pthread_create_opaque(D%threadptr, c_funloc(SaveBuffer), c_loc(D), err)

          if (err==0) then
            D%in_progress = .true.
          else
            write (*,*) "Error in creating frame thread. Will run again synchronously. Code:",err
            call SaveBuffer(c_loc(D))
          end if

        end if

      end associate
    end subroutine TStaggeredFrameDomain_Save



    subroutine TStaggeredFrameDomain_SaveTimes(D)
      type(TStaggeredFrameDomain),target,intent(in) :: D
      character(40) :: file_name
      integer i

      file_name = trim(D%base_name)//"-times.txt"

      open(unit=D%unit, file=file_name, access='sequential', form='formatted', action='write', status='replace')

      do i = 0, D%frame_number
        write(D%unit,'(i8,2x,es12.5)') i, D%times(i)
      end do
      close(D%unit)

    end subroutine TStaggeredFrameDomain_SaveTimes



    subroutine TStaggeredFrameDomain_SaveHeader(D)
      type(TStaggeredFrameDomain),target,intent(in) :: D
      character(40) :: file_name

      file_name = trim(D%base_name)//"-header"//trim(D%suffix)

      open(unit=D%unit, file=file_name, access='stream', form='unformatted', action='write', status='replace')

      write(D%unit) 1_int32 !endianess can be infered from this
      write(D%unit) int(storage_size(1._knd),int32)  !save number of bits of the used real kind

      write(D%unit) D%save_flags%U, D%save_flags%V, D%save_flags%W, D%save_flags%Pr
      write(D%unit) D%save_flags%Temperature, D%save_flags%Scalar, D%save_flags%num_scalars

      write(D%unit) D%minUi, D%maxUi, D%minUj, D%maxUj, D%minUk, D%maxUk
      write(D%unit) D%minVi, D%maxVi, D%minVj, D%maxVj, D%minVk, D%maxVk
      write(D%unit) D%minWi, D%maxWi, D%minWj, D%maxWj, D%minWk, D%maxWk
      write(D%unit) D%minPri, D%maxPri, D%minPrj, D%maxPrj, D%minPrk, D%maxPrk

      write(D%unit) D%xPr,D%yPr,D%zPr,D%xU,D%yV,D%zW

      close(D%unit)

    end subroutine TStaggeredFrameDomain_SaveHeader



    subroutine TStaggeredFrameDomain_Finalize(D)
      type(TStaggeredFrameDomain),intent(inout) :: D

      if (D%in_progress) call Wait(D)

      call SaveTimes(D)

      if (allocated(D%buffer)) deallocate(D%buffer)
      if (allocated(D%times)) deallocate(D%times)
      if (allocated(D%xPr)) deallocate(D%xPr)
      if (allocated(D%yPr)) deallocate(D%yPr)
      if (allocated(D%zPr)) deallocate(D%zPr)
      if (allocated(D%xU)) deallocate(D%xU)
      if (allocated(D%yV)) deallocate(D%yV)
      if (allocated(D%zW)) deallocate(D%zW)

    end subroutine TStaggeredFrameDomain_Finalize

end module StaggeredFrames







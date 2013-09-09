module Puffs

  use Parameters
  
  use ArrayUtilities
  
  use VolumeSources, only: ScalarFlVolumesContainer, ScalarFlVolume
  
  implicit none
  
  type PuffSource
    real(knd) :: start_time,end_time,period, duration, last_start = -huge(1._knd)
    logical   :: clear !clear the scalar field? (all scalars for which source below exists)
    type(ScalarFlVolumesContainer),allocatable :: sources(:) !must be allocated to 0 if empty
  end type

  type(PuffSource),allocatable :: PuffSources(:) !index means individual puff source objects
  
  integer :: enable_puffs = 1

  
contains

  subroutine DoPuffs(Scalar, Scalar_derivative, stage, RK_stages)
    real(knd),dimension(-1:,-1:,-1:,1:),intent(inout) :: Scalar, Scalar_derivative
    integer,intent(in) :: stage, RK_stages
    integer :: i
      
    if (allocated(PuffSources)) then
      do i = 1,size(PuffSources)
        call DoPuff(PuffSources(i), Scalar, Scalar_derivative, stage, RK_stages)
      end do
    end if
  end subroutine
  
  
  subroutine DoPuff(p, Scalar, Scalar_derivative, stage, RK_stages)
    type(PuffSource),intent(inout) :: p
    real(knd),dimension(-1:,-1:,-1:,1:),intent(inout) :: Scalar, Scalar_derivative
    integer,intent(in) :: stage, RK_stages
    logical :: restart
    real(knd) :: effective_dt
    integer :: i,j
  
    !we assume period >> dt and duration >> dt in all following
    if (time+dt > p%last_start + p%period .and. time+dt >= p%start_time .and. time <= p%end_time) then
      restart = .true.
    else
      restart = .false.
    end if
    

    if ((.not.restart) .and. (time <= p%last_start + p%duration) .and. (time >= p%start_time)) then
      effective_dt = max(min(p%last_start + p%duration, time + dt, p%end_time) - time, 0._knd)
    else if (restart) then
      effective_dt = time + dt - max(p%last_start + p%period, time, p%start_time)
    else
      effective_dt = 0
    end if

    do i=1,size(p%sources)

      associate(src => p%sources(i))
      
        if (p%clear.and.restart.and.stage==1) then
          call set(Scalar(:,:,:,src%scalar_number), 0)
        end if
        
        if (effective_dt>0) then
      
          do j=1,size(src%volumes)
            associate (xi => src%volumes(j)%xi, &
                       yj => src%volumes(j)%yj, &
                       zk => src%volumes(j)%zk)
                Scalar_derivative(xi,yj,zk,src%scalar_number) = &
                  Scalar_derivative(xi,yj,zk,src%scalar_number) + effective_dt * src%volumes(j)%flux
            end associate
          end do
        
        end if
        
      end associate

    end do
    
    if (restart.and.stage==RK_stages) then
      if (p%last_start>=p%start_time) then
        p%last_start = p%last_start + p%period
      else
        p%last_start = p%start_time
      end if
    end if

  end subroutine
  
  
  subroutine InitPuffSources
  end subroutine



end module Puffs

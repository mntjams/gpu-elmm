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
   use Boundaries, only: GridCoords
   !!!FIXME: Should come up with reading from a configuration file.
!    integer :: xi(3),yj(3),zk(3)
!    if (.not.allocated(PuffSources)) then
!       allocate(PuffSources(3))
!       call GridCoords(xi(1),yj(1),zk(1), &
!                  -361.9_knd,125.1_knd,0.1_knd)
!       call GridCoords(xi(2),yj(2),zk(2), &
!                  -64.4_knd,38.3_knd,0.1_knd)
!       call GridCoords(xi(3),yj(3),zk(3), &
!                  -0._knd,0._knd,0.1_knd)
!       PuffSources(1) = PuffSource(start_time=0.,&
!                                   end_time=10000.0,&
!                                   period=900.,&
!                                   duration=29.,&
!                                   clear=.true., &
!                                   sources=[ScalarFlVolumesContainer(1, &
!                                             [ScalarFlVolume(xi(1),yj(1),zk(2),1./(dxmin*dymin*dzmin))]), &
!                                            ScalarFlVolumesContainer(2, &
!                                             [ScalarFlVolume(xi(2),yj(2),zk(3),1./(dxmin*dymin*dzmin))]), &
!                                            ScalarFlVolumesContainer(3, &
!                                             [ScalarFlVolume(xi(3),yj(3),zk(3),1./(dxmin*dymin*dzmin))]) &
!                                           ] &
!                                  )
!       PuffSources(2) = PuffSource(300.,10000.0,900.,29,clear=.true., &
!                                   sources=[ScalarFlVolumesContainer(4, &
!                                             [ScalarFlVolume(xi(1),yj(1),zk(2),1./(dxmin*dymin*dzmin))]), &
!                                            ScalarFlVolumesContainer(5, &
!                                             [ScalarFlVolume(xi(2),yj(2),zk(3),1./(dxmin*dymin*dzmin))]), &
!                                            ScalarFlVolumesContainer(6, &
!                                             [ScalarFlVolume(xi(3),yj(3),zk(3),1./(dxmin*dymin*dzmin))]) &
!                                           ] &
!                                  )
!       PuffSources(3) = PuffSource(600.,10000.0,900.,29,clear=.true., &
!                                   sources=[ScalarFlVolumesContainer(7, &
!                                             [ScalarFlVolume(xi(1),yj(1),zk(2),1./(dxmin*dymin*dzmin))]), &
!                                            ScalarFlVolumesContainer(8, &
!                                             [ScalarFlVolume(xi(2),yj(2),zk(3),1./(dxmin*dymin*dzmin))]), &
!                                            ScalarFlVolumesContainer(9, &
!                                             [ScalarFlVolume(xi(3),yj(3),zk(3),1./(dxmin*dymin*dzmin))]) &
!                                           ] &
!                                  )
!     end if

  end subroutine



end module Puffs

module SolidBodies

  use Parameters
  use Body_class
  use GeometricShapes

  implicit none

  private

  public SolidBody, InitSolidBodies, SetCurrentSB, FindInsideCells, &
         obstacles_file, roughness_file, displacement_file
#ifdef CUSTOMSB
  public  AddSolidBody, SolidBodiesList
#endif

  type, extends(Body) :: SolidBody
    logical   :: rough = .false.                             !T rough surface, F flat surface
    real(knd) :: z0 = 0                                      !roughness parameter
    real(knd) :: temperature_flux = 0
    real(knd) :: moisture_flux = 0
    !other scalar fluxes assumed zero
    
    !building brick - Santamouris - Environmental Design of Urban Buildings
    real(knd) :: emissivity = 0.45 
    real(knd) :: albedo = 0.3
  end type SolidBody
  
#define TYPEPARAM type(SolidBody)
#include "list-inc-def.f90"

  type(List) :: SolidBodiesList

  character(1024) :: obstacles_file = ''
  
  character(1024) :: roughness_file = ''
  
  character(1024) :: displacement_file = ''
  
  interface AddSolidBody
    module procedure AddSolidBody_scalar
    module procedure AddSolidBody_array
  end interface
  
  interface SolidBody
    module procedure SolidBody_Init
  end interface

contains

#include "list-inc-proc.f90"
#undef TYPEPARAM

  subroutine SetCurrentSB(SB,n)
    type(SolidBody),pointer,intent(out) :: SB
    integer,intent(in) :: n
    integer i

    call SolidBodiesList%iter_restart

    do i = 1,n
      SB => SolidBodiesList%iter_next()
    enddo

    if (.not.associated(SB)) then
      write (*,*) "Error, no SolidBody with number",n,"in the list."
      call error_stop
    end if

  end subroutine SetCurrentSB








  subroutine FindInsideCells
    !find if the gridpoints lie inside a solid body and write it's number
    !do not nullify the .type arrays, they could have been made nonzero by other unit

    call SolidBodiesList%for_each(SetPrtype)
    call SolidBodiesList%for_each(SetUtype)
    call SolidBodiesList%for_each(SetVtype)
    call SolidBodiesList%for_each(SetWtype)

    contains

      subroutine SetPrtype(CurrentSB)
        type(SolidBody) :: CurrentSB
        integer i,j,k

        if (CurrentSB%numofbody==0) call error_stop("Error, numofbody==0, did you use AddSolidBody()?")
    
        !$omp parallel do private(i,j,k) schedule(dynamic)
        do k = 0,Prnz+1
         do j = 0,Prny+1
          do i = 0,Prnx+1
             if (CurrentSB%Inside(xPr(i),yPr(j),zPr(k),(dxmin*dymin*dzmin)**(1._knd/3)/20)) &
                      Prtype(i,j,k) = CurrentSB%numofbody
          enddo
         enddo
        enddo
        !$omp end parallel do
      end subroutine

      subroutine SetUtype(CurrentSB)
        type(SolidBody) :: CurrentSB
        integer i,j,k

        if (CurrentSB%numofbody==0) call error_stop("Error, numofbody==0, did you use AddSolidBody()?")
    
        !$omp parallel do private(i,j,k) schedule(dynamic)
        do k = 0,Unz+1
         do j = 0,Uny+1
          do i = 0,Unx+1
             if (CurrentSB%Inside(xU(i),yPr(j),zPr(k),(dxmin*dymin*dzmin)**(1._knd/3)/20))&
                       Utype(i,j,k) = CurrentSB%numofbody
          enddo
         enddo
        enddo
        !$omp end parallel do
      end subroutine

      subroutine SetVtype(CurrentSB)
        type(SolidBody) :: CurrentSB
        integer i,j,k

        if (CurrentSB%numofbody==0) call error_stop("Error, numofbody==0, did you use AddSolidBody()?")
    
        !$omp parallel do private(i,j,k) schedule(dynamic)
        do k = 0,Vnz+1
         do j = 0,Vny+1
          do i = 0,Vnx+1
             if (CurrentSB%Inside(xPr(i),yV(j),zPr(k),(dxmin*dymin*dzmin)**(1._knd/3)/20))&
                       Vtype(i,j,k) = CurrentSB%numofbody
          enddo
         enddo
        enddo
        !$omp end parallel do
      end subroutine

      subroutine SetWtype(CurrentSB)
        type(SolidBody) :: CurrentSB
        integer i,j,k

        !$omp parallel do private(i,j,k) schedule(dynamic)
        do k = 0,Wnz+1
         do j = 0,Wny+1
          do i = 0,Wnx+1
             if (CurrentSB%Inside(xPr(i),yPr(j),zW(k),(dxmin*dymin*dzmin)**(1._knd/3)/20))&
                       Wtype(i,j,k) = CurrentSB%numofbody
          enddo
         enddo
        enddo
        !$omp end parallel do
      end subroutine

  end subroutine FindInsideCells




  subroutine InitSolidBodies

#ifdef CUSTOMSB
    interface
      subroutine CustomSolidBodies
      end subroutine
    end interface
    !An external subroutine, it should use this module and use AddSolidBody to supply
    ! pointers to the new solid bodies.
    call CustomSolidBodies
#endif

    if (len_trim(obstacles_file)>0) then
      call ReadSolidBodiesFromFile(trim(obstacles_file))
    end if

  end subroutine InitSolidBodies


  subroutine ReadSolidBodiesFromFile(filename)
    use Strings
    character(*),intent(in) :: filename
    integer :: l
    character(5) :: suffix
    
    l = len(filename)
    
    suffix = filename(index(filename,'.',back=.true.):)

    if (suffix=='.obst'.or.suffix=='.geom') then
      call ReadUnion(filename)
    else if (suffix=='.off') then
      call ReadOff(filename)
    else if (suffix=='.xyz'.or.suffix=='.elv') then
      call ReadTerrain(filename)
    else if (suffix=='.ltop') then
      call ReadTopPoints(filename,.false.)
    else if (suffix=='.rtop') then
      call ReadTopPoints(filename,.true.)
    else
      write(*,*) "Unknown file format "//suffix
      call error_stop
    end if

  end subroutine ReadSolidBodiesFromFile
  
  

  subroutine ReadUnion(filename)
    character(*),intent(in) :: filename
    
    !assume z0 the same as for the lower boundary
    call AddSolidBody(SolidBody(Union(filename), z0 = z0B))
  end subroutine ReadUnion
  
  
  
  subroutine ReadOff(filename)
    character(*),intent(in) :: filename
    logical :: ex

    inquire(file=filename,exist=ex)

    if (ex) then
      !assume z0 the same as for the lower boundary
      call AddSolidBody(SolidBody(Polyhedron(filename), z0 = z0B))
    else
      call error_stop("Error, file "//filename//" does not exist.")
    end if

  end subroutine ReadOff
  
  
  subroutine ReadTopPoints(filename, right)
    character(*),intent(in) :: filename
    logical, intent(in) :: right
    logical :: ex, body_opened
    real(knd), allocatable :: points(:,:)
    integer :: u, stat

    inquire(file=filename,exist=ex)

    if (ex) then
      
      open(newunit=u, file=filename)
      
      body_opened = .false.
      do
        call next_line(stat)
        if (stat>0) exit
      end do
      
      close(u)
    else
      call error_stop("Error, file "//filename//" does not exist.")
    end if

  contains
    subroutine next_line(stat)
      integer, intent(out) :: stat
      character(1024) :: line
      integer :: io
      
      read(u,'(a)', iostat=io) line
      
      if (io/=0) then
        if (body_opened) call close_body
        stat = 1
        return
      end if
      
      if (len_trim(line)==0) then  
        if (body_opened) call close_body
        stat = 3
        return
      end if
        
      call next_point(line, stat)
    end subroutine
    
    subroutine next_point(line, stat)
      character(*), intent(in) :: line
      integer, intent(out) :: stat
      real(knd) :: point(3)
      real(knd), allocatable :: tmp(:,:)
      integer :: io
      
      read(line,*,iostat=io) point
      if (io/=0) then
        stat = 2
        return
      end if
     
      if (.not.body_opened) then
        body_opened = .true.
        allocate(points(3,1))
      else
        tmp = points
        deallocate(points)
        allocate(points(3,size(tmp,2)+1))
        points(:,1:size(tmp,2)) = tmp
      end if
      points(:,size(points,2)) = point
      
      stat = 0
    end subroutine
    
    subroutine close_body
      !assume z0 the same as for the lower boundary
      if (right) then
        call AddSolidBody(SolidBody(ConvexPolyhedron_FromTopPoints(points), z0 = z0B))
      else
        call AddSolidBody(SolidBody(ConvexPolyhedron_FromTopPoints(points(:,size(points,2):1:-1)), z0 = z0B))
      end if

      deallocate(points)
      body_opened = .false.
    end subroutine

  end subroutine ReadTopPoints
  
  
  
  subroutine ReadTerrain(filename)
    use ElevationModels, only: uniform_map
    character(*), intent(in) :: filename

    !assume z0 the same as for the lower boundary if not specified otherwise
    
    if (len_trim(roughness_file)>0) then
      if (len_trim(displacement_file)>0) then
        call AddSolidBody(SolidBody(Terrain(uniform_map(filename), &
                                            uniform_map(trim(roughness_file), default = z0B), &
                                            uniform_map(trim(displacement_file), default = 0._knd))))
      else
        call AddSolidBody(SolidBody(Terrain(uniform_map(trim(filename)), &
                                            uniform_map(trim(roughness_file), default = z0B))))
      end if
    else
      call AddSolidBody(SolidBody(Terrain(uniform_map(filename), z0 = z0B)))
    end if

  end subroutine ReadTerrain

  
  
  subroutine AddSolidBody_scalar(SB)
    type(SolidBody),intent(in) :: SB
    !NOTE: expecting numofbody value unspecified and not important in the calling code
    type(SolidBody) :: tmp

    tmp = SB
    tmp%numofbody = SolidBodiesList%Len() + 1

    call SolidBodiesList%add(tmp)

  end subroutine

  subroutine AddSolidBody_array(SB)
    type(SolidBody),intent(in) :: SB(:)
    !NOTE: expecting numofbody value unspecified and not important in the calling code
    integer :: i
    type(SolidBody) :: tmp
    
    do i = 1,size(SB)
      tmp = SB(i)
      tmp%numofbody = SolidBodiesList%Len() + 1

      call SolidBodiesList%add(tmp)
    end do

  end subroutine
  
  function SolidBody_Init(gs, z0, temperature_flux, moisture_flux) result(res)
    type(SolidBody) :: res
    class(GeometricShape), intent(in) :: gs
    real(knd), optional, intent(in) :: z0, temperature_flux, moisture_flux
    
    allocate(res%GeometricShape, source=gs)
    if (present(z0)) res%z0 = z0
    res%rough = res%z0 > 0
    if (present(temperature_flux)) res%temperature_flux = temperature_flux
    if (present(moisture_flux)) res%moisture_flux = moisture_flux
  end function

end module SolidBodies

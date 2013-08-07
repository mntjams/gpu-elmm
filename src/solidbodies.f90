module SolidBodies

  use Parameters
  use Lists
  use Body_class
  use GeometricShapes

  implicit none

  private

  public SolidBody, InitSolidBodies, SetCurrentSB, FindInsideCells, obstacles_file
#ifdef CUSTOMSB
  public  AddSolidBody, SolidBodiesList
#endif

  type, extends(Body) :: SolidBody
    logical   :: rough = .false.                             !T rough surface, F flat surface
    real(knd) :: z0 = 0                                      !roughness parameter
    real(knd) :: temperature_flux = 0
    real(knd) :: moisture_flux = 0
    !other scalar fluxes assumed zero
  end type SolidBody

  type(List) :: SolidBodiesList

  character(80) :: obstacles_file = ''
  
  interface AddSolidBody
    module procedure AddSolidBody_scalar
    module procedure AddSolidBody_array
  end interface

contains



  subroutine SetCurrentSB(SB,n)
    type(SolidBody),pointer,intent(out) :: SB
    integer,intent(in) :: n
    class(Listable),pointer :: ptr
    integer i

    call SolidBodiesList%IterRestart

    do i = 1,n
      call SolidBodiesList%IterNext(ptr)
    enddo

    if (associated(ptr)) then
      select type (ptr)
        type is (SolidBody)
          SB => ptr
        class default
          write (*,*) "Error, wrong dynamic type of item in SolidBodiesList."
          stop
      end select
    else
      write (*,*) "Error, no SolidBody with number",n,"in the list."
      stop
    end if

  end subroutine SetCurrentSB








  subroutine FindInsideCells
    !find if the gridpoints lie inside a solid body and write it's number
    !do not nullify the .type arrays, they could have been made nonzero by other unit

    call SolidBodiesList%ForEach(SetPrtype)
    call SolidBodiesList%ForEach(SetUtype)
    call SolidBodiesList%ForEach(SetVtype)
    call SolidBodiesList%ForEach(SetWtype)

    contains

      subroutine SetPrtype(CurrentSB)
        class(Listable) :: CurrentSB
        integer i,j,k

        select type (CurrentSB)
          type is (SolidBody)

            if (CurrentSB%numofbody==0) stop "Error, numofbody==0, did you use AddSolidBody()?"
        
            !$omp parallel do private(i,j,k) schedule(dynamic)
            do k = 0,Prnz+1
             do j = 0,Prny+1
              do i = 0,Prnx+1
                 if (CurrentSB%Inside(xPr(i),yPr(j),zPr(k),(dxmin*dymin*dzmin)**(1._knd/3)/1000)) &
                          Prtype(i,j,k) = CurrentSB%numofbody
              enddo
             enddo
            enddo
            !$omp end parallel do
        end select
      end subroutine

      subroutine SetUtype(CurrentSB)
        class(Listable) :: CurrentSB
        integer i,j,k

        select type (CurrentSB)
          type is (SolidBody)

            if (CurrentSB%numofbody==0) stop "Error, numofbody==0, did you use AddSolidBody()?"
        
            !$omp parallel do private(i,j,k) schedule(dynamic)
            do k = 0,Unz+1
             do j = 0,Uny+1
               do i = 0,Unx+1
                  if (CurrentSB%Inside(xU(i),yPr(j),zPr(k),(dxmin*dymin*dzmin)**(1._knd/3)/1000))&
                           Utype(i,j,k) = CurrentSB%numofbody
               enddo
              enddo
             enddo
            !$omp end parallel do
        end select
      end subroutine

      subroutine SetVtype(CurrentSB)
        class(Listable) :: CurrentSB
        integer i,j,k

        select type (CurrentSB)
          type is (SolidBody)

            if (CurrentSB%numofbody==0) stop "Error, numofbody==0, did you use AddSolidBody()?"
        
            !$omp parallel do private(i,j,k) schedule(dynamic)
            do k = 0,Vnz+1
             do j = 0,Vny+1
              do i = 0,Vnx+1
                 if (CurrentSB%Inside(xPr(i),yV(j),zPr(k),(dxmin*dymin*dzmin)**(1._knd/3)/1000))&
                           Vtype(i,j,k) = CurrentSB%numofbody
              enddo
             enddo
            enddo
            !$omp end parallel do
        end select
      end subroutine

      subroutine SetWtype(CurrentSB)
        class(Listable) :: CurrentSB
        integer i,j,k

         select type (CurrentSB)
          type is (SolidBody)
            !$omp parallel do private(i,j,k) schedule(dynamic)
            do k = 0,Wnz+1
             do j = 0,Wny+1
              do i = 0,Wnx+1
                 if (CurrentSB%Inside(xPr(i),yPr(j),zW(k),(dxmin*dymin*dzmin)**(1._knd/3)/1000))&
                           Wtype(i,j,k) = CurrentSB%numofbody
              enddo
             enddo
            enddo
            !$omp end parallel do
        end select
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

    if (filename(len(filename)-4:len(filename))=='.obst') then
      call ReadObst(filename)
    else if (filename(len(filename)-3:len(filename))=='.off') then
      call ReadOff(filename)
    end if
    
  end subroutine ReadSolidBodiesFromFile
  
  

  subroutine ReadObst(filename)
    use Strings
    character(*),intent(in) :: filename
    integer unit,io
    character(180) :: line
    type(SolidBody) :: SB

    allocate(SB%GeometricShape, source=Union(filename))
    
    call AddSolidBody(SB)
  end subroutine ReadObst
  
  
  
  subroutine ReadOff(filename)
    character(*),intent(in) :: filename
    logical :: ex
    type(SolidBody) :: SB

    inquire(file=filename,exist=ex)

    if (ex) then
      allocate(Polyhedron :: SB%GeometricShape)
      
      select type (geomshape => SB%GeometricShape)
        type is (Polyhedron)
       
          geomshape = Polyhedron(filename)
      
          call AddSolidBody(SB)
      end select
      
    else
      write(*,*) "Error, file ",filename," does not exist."
      stop
    end if

  end subroutine ReadOff



  subroutine AddSolidBody_scalar(SB)
    class(Body),intent(inout) :: SB
    !NOTE: expecting numofbody value unspecified and not important in the calling code

    SB%numofbody = SolidBodiesList%Len() + 1

    call SolidBodiesList%Add(SB)

  end subroutine

  subroutine AddSolidBody_array(SB)
    class(Body),intent(inout) :: SB(:)
    !NOTE: expecting numofbody value unspecified and not important in the calling code
    integer :: i
    
    do i = 1,size(SB)
      SB(i)%numofbody = SolidBodiesList%Len() + 1

      call SolidBodiesList%Add(SB(i))
    end do

  end subroutine

end module SolidBodies

module SolidBodies

  use Parameters
  use Lists
  use TBody_class
  use GeometricShapes

  implicit none

  private

  public TSolidBody, InitSolidBodies, SetCurrentSB, FindInsideCells, obstaclefile
#ifdef CUSTOMSB
  public  AddBody, SolidBodiesList
#endif

  type, extends(TBody) :: TSolidBody
    logical   :: rough = .false.                             !T rough surface, F flat surface
    real(KND) :: z0 = 0                                      !roughness parameter
    real(KND) :: temperatureflux = 0
  end type TSolidbody

  type(TList) :: SolidBodiesList

  character(80) :: obstaclefile = ''

contains



  subroutine SetCurrentSB(SB,n)
    type(TSolidBody),pointer,intent(out) :: SB
    integer,intent(in) :: n
    class(TListable),pointer :: ptr
    integer i

    call SolidBodiesList%IterRestart

    do i = 1,n
      call SolidBodiesList%IterNext(ptr)
    enddo

    if (associated(ptr)) then
      select type (ptr)
        type is (TSolidBody)
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
        class(TListable) :: CurrentSB
        integer i,j,k

        select type (CurrentSB)
          type is (TSolidBody)
            !$omp parallel do private(i,j,k)
            do k = 0,Prnz+1
             do j = 0,Prny+1
              do i = 0,Prnx+1
                 if (CurrentSB%Inside(xPr(i),yPr(j),zPr(k))) Prtype(i,j,k) = CurrentSB%numofbody
              enddo
             enddo
            enddo
            !$omp end parallel do
        end select
      end subroutine

      subroutine SetUtype(CurrentSB)
        class(TListable) :: CurrentSB
        integer i,j,k

        select type (CurrentSB)
          type is (TSolidBody)
            !$omp parallel do private(i,j,k)
            do k = 0,Unz+1
             do j = 0,Uny+1
               do i = 0,Unx+1
                  if (CurrentSB%Inside(xU(i),yPr(j),zPr(k),(dxmin*dymin*dzmin)**(1._KND/3)/1000))&
                           Utype(i,j,k) = CurrentSB%numofbody
               enddo
              enddo
             enddo
            !$omp end parallel do
        end select
      end subroutine

      subroutine SetVtype(CurrentSB)
        class(TListable) :: CurrentSB
        integer i,j,k

        select type (CurrentSB)
          type is (TSolidBody)
            !$omp parallel do private(i,j,k)
            do k = 0,Vnz+1
             do j = 0,Vny+1
              do i = 0,Vnx+1
                 if (CurrentSB%Inside(xPr(i),yV(j),zPr(k),(dxmin*dymin*dzmin)**(1._KND/3)/1000))&
                           Vtype(i,j,k) = CurrentSB%numofbody
              enddo
             enddo
            enddo
            !$omp end parallel do
        end select
      end subroutine

      subroutine SetWtype(CurrentSB)
        class(TListable) :: CurrentSB
        integer i,j,k

         select type (CurrentSB)
          type is (TSolidBody)
            !$omp parallel do private(i,j,k)
            do k = 0,Wnz+1
             do j = 0,Wny+1
              do i = 0,Wnx+1
                 if (CurrentSB%Inside(xPr(i),yPr(j),zW(k),(dxmin*dymin*dzmin)**(1._KND/3)/1000))&
                           Wtype(i,j,k) = CurrentSB%numofbody
              enddo
             enddo
            enddo
            !$omp end parallel do
        end select
      end subroutine

  end subroutine FindInsideCells




  subroutine InitSolidBodies

    if (len_trim(obstaclefile)>0) then
      call ReadSolidBodiesFromFile(obstaclefile)
    end if

#ifdef CUSTOMSB
    !An external subroutine, it should use this module and use AddBody to supply
    ! pointers to the new solid bodies.
    call CustomSolidBodies
#endif
  end subroutine InitSolidBodies




  subroutine ReadSolidBodiesFromFile(filename)
    use Strings
    character(*),intent(in) :: filename
    integer unit,io
    character(180) :: line
    type(TSolidBody),pointer :: SB => null()

    unit=20
    !provisional, Solaris Studio still doesn't support newunit.
    open(unit=unit,file=filename,action='read',status='old',iostat=io)

    if (io==0) then
      do
        read(unit,'(a)',iostat=io) line

        if (io/=0) exit

        line = adjustl(line)

        if (len_trim(line)>0) then

          if (upcase(line(1:10))=='POLYHEDRON') then

            call ReadPolyhedron(unit,line(11:),SB)
          end if

          if (associated(SB)) then

            call AddBody(SolidBodiesList,SB)
            deallocate(SB)

          end if

        end if

      end do

      close(unit)

    else

      write(*,*) "Could not open file",filename
      stop

    end if
  end subroutine ReadSolidBodiesFromFile


  subroutine ReadPolyhedron(unit,restline,SB)
    integer,intent(in)       :: unit
    character(*),intent(in)  :: restline
    type(TSolidBody),pointer :: SB
    integer nPlanes,i,io

    read(restline,*,iostat=io) nPlanes

    if (io/=0) then
      write(*,*) "Expected number of planes in polyhedron, received '",trim(restline),"' instead."
      stop
    end if

    allocate(SB)
    allocate(TPolyhedron :: SB%GeometricShape)

    select type (geomshape => SB%GeometricShape)
      type is (TPolyhedron)
        allocate(geomshape%Planes(nPlanes))

        geomshape%nplanes = nPlanes

        do i=1,nPlanes
          call ReadPlane(unit,geomshape%Planes(i))
        end do
    end select

  end subroutine ReadPolyhedron


  subroutine ReadPlane(unit,Pl)
    use Strings
    integer,intent(in) :: unit
    type(TPlane),intent(inout) :: Pl
    character(180) :: line
    integer io

    read(unit,'(a)',iostat=io) line
    if (io/=0) then
      write(*,*) "Error reading the line with the plane definition."
      stop
    end if

    if (count_multispaces(line) == 4) then
      read(line,*,iostat=io) Pl%a,Pl%b,Pl%c,Pl%d,Pl%gl
    else if (count_multispaces(line) == 6) then
      read(line,*,iostat=io) Pl%a,Pl%b,Pl%c,Pl%d,Pl%gl,Pl%rough, Pl%z0
    else
      io = 999
    end if
    if (io/=0) then
      write(*,*) "Error parsing the line with the plane definition."
      stop
    end if
  end subroutine ReadPlane



  subroutine AddBody(List,SB)
    type(TList),intent(inout) :: List
    class(TBody) :: SB

    SB%numofbody = List%Len() + 1

    call List%Add(SB)

  end subroutine AddBody


end module SolidBodies

module vtkarray
 !Simple module to output arrays for visualization. No physical coordinates are used, only the position in the array.
 !Mostly only for debugging.

!   use iso_fortran_env, only: real32, real64
  use FreeUnit

  implicit none

  integer,parameter :: real32=selected_real_kind(p=6,r=37)
  integer,parameter :: real64=selected_real_kind(p=15,r=200)

  interface VtkArraySimple
    module procedure SVtkArraySimple
    module procedure DVtkArraySimple
  end interface

  contains

    subroutine SVtkArraySimple(fname,A)
      character(len=*)              :: fname
      real(real32),dimension(:,:,:) :: A
      integer                       :: nx,ny,nz
      integer                       :: i
      integer                       :: unit
      character(len=40)             :: str

      nx=Ubound(A,1)
      ny=Ubound(A,2)
      nz=Ubound(A,3)

      call newunit(unit)

      open(unit,file=fname,status="replace",action="write")
      write (unit,"(A)") "# vtk DataFile Version 2.0"
      write (unit,"(A)") "CLMM output file"
      write (unit,"(A)") "ASCII"
      write (unit,"(A)") "DATASET RECTILINEAR_GRID"
      str="DIMENSIONS"
      write (str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
      write (unit,"(A)") trim(str)
      str="X_COORDINATES"
      write (str(15:),'(i5,2x,a)') nx,"float"
      write (unit,"(A)") str
      write (unit,*) (i, i=1,nx)
      str="Y_COORDINATES"
      write (str(15:),'(i5,2x,a)') ny,"float"
      write (unit,"(A)") trim(str)
      write (unit,*) (i, i=1,ny)
      str="Z_COORDINATES"
      write (str(15:),'(i5,2x,a)') nz,"float"
      write (unit,"(A)") trim(str)
      write (unit,*) (i, i=1,nz)
      str="POINT_DATA"
      write (str(12:),*) nx*ny*nz
      write (unit,"(A)") trim(str)


      write (unit,"(A)") "SCALARS array float"
      write (unit,"(A)") "LOOKUP_TABLE default"

      write (unit,'(*(g0,/))') A(1:nx,1:ny,1:nz)

      write (unit,*)
      close(unit)

    end subroutine SVtkArraySimple

    subroutine DVtkArraySimple(fname,A)
      character(len=*)              :: fname
      real(real64),dimension(:,:,:) :: A
      integer                       :: nx,ny,nz
      integer                       :: i
      integer                       :: unit
      character(len=40)             :: str

      nx=Ubound(A,1)
      ny=Ubound(A,2)
      nz=Ubound(A,3)

      call newunit(unit)

      open(unit,file=fname,status="replace",action="write")
      write (unit,"(A)") "# vtk DataFile Version 2.0"
      write (unit,"(A)") "CLMM output file"
      write (unit,"(A)") "ASCII"
      write (unit,"(A)") "DATASET RECTILINEAR_GRID"
      str="DIMENSIONS"
      write (str(12:),'(i0,1x,i0,1x,i0)') nx,ny,nz
      write (unit,"(A)") trim(str)
      str="X_COORDINATES"
      write (str(15:),'(i5,2x,a)') nx,"float"
      write (unit,"(A)") str
      write (unit,*) (i, i=1,nx)
      str="Y_COORDINATES"
      write (str(15:),'(i5,2x,a)') ny,"float"
      write (unit,"(A)") trim(str)
      write (unit,*) (i, i=1,ny)
      str="Z_COORDINATES"
      write (str(15:),'(i5,2x,a)') nz,"float"
      write (unit,"(A)") trim(str)
      write (unit,*) (i, i=1,nz)
      str="POINT_DATA"
      write (str(12:),*) nx*ny*nz
      write (unit,"(A)") trim(str)


      write (unit,"(A)") "SCALARS array float"
      write (unit,"(A)") "LOOKUP_TABLE default"

      write (unit,'(*(g0,/))') A(1:nx,1:ny,1:nz)

      write (unit,*)
      close(unit)

    end subroutine DVtkArraySimple
end module vtkarray



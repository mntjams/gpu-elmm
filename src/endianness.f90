module Endianness
!   use iso_fortran_env
  implicit none

  private
  public :: GetEndianness, BigEnd

  integer,parameter :: int8=selected_int_kind(1),int32=selected_int_kind(9),int64=selected_int_kind(10)
  integer,parameter :: real32=selected_real_kind(p=6,r=37),real64=selected_real_kind(p=15,r=200)
  logical,save :: littleendian

  interface BigEnd
    module procedure BigEnd32
    module procedure BigEnd64
  end interface

  contains

    subroutine GetEndianness
      integer(int8),dimension(4):: bytes !may not work on some processors

      bytes=transfer(1_int32,bytes,4)
      if (bytes(4)==1) then
        littleendian=.false.
      else
        littleendian=.true.
      endif
    end subroutine GetEndianness

    elemental function BigEnd32(x) result(res)
      real(real32) :: res
      real(real32),intent(in)::x
      integer(int8),dimension(4):: bytes !may not work on some processors

      if (.not.littleendian) then
        res = x
      else
        bytes = transfer(x,bytes,4)
        res = transfer(bytes(4:1:-1),res)
      endif
    end function BigEnd32

    elemental function BigEnd64(x) result(res)
      real(real64) :: res
      real(real64),intent(in)::x
      integer(int8),dimension(8):: bytes !may not work on some processors

      if (.not.littleendian) then
        res = x
      else
        bytes = transfer(x,bytes,8)
        res = transfer(bytes(8:1:-1),res)
      endif
    end function BigEnd64

end module Endianness


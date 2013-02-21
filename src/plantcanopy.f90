module PlantCanopy
  use VolumeSources

  implicit none

  type, extends(TVolumeSourceBody) :: TPlantBody
    integer :: plant_type !problem specific, used by custom routines
  end type TPlantBody

  contains

      subroutine InitPlantBodies
#ifdef CUSTOMPB
      !An external subroutine, it should use this module and use AddBody to supply
      ! pointers to the new solid bodies.
      call CustomPlantBodies
#endif
    end subroutine InitPlantBodies


end module PlantCanopy

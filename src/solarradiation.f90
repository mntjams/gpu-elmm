module SolarRadiation
  use Parameters
  
  implicit none
  
  real(knd),parameter :: SB_sigma = 5.6704E-8_knd

  real(knd) :: sun_azimuth, sun_elevation !in degrees
  real(knd) :: vector_to_sun(3) !unit vector pointing to sun in grid coordinates
  
  integer :: svf_nrays = 50
  
  real(knd),allocatable :: svf_vecs(:,:)
  
contains
  
  pure function solar_direct_flux() result(res)
    real(knd) :: res
    res = 800._knd
  end function

  pure function solar_diffuse_flux() result(res)
    real(knd) :: res
    res = 200._knd
  end function

  subroutine InitSolarRadiation
    real(knd) :: horiz_component, horiz_angle
    real(knd) :: azimuth,elevation
    integer i
   
    sun_azimuth = 237
    sun_elevation = 53
    
    vector_to_sun(3) = sin(sun_elevation * pi / 180)
    horiz_component = cos(sun_elevation * pi / 180)
    horiz_angle = (x_axis_azimuth - sun_azimuth) * pi / 180
    
    vector_to_sun(1) = horiz_component * cos(horiz_angle)
    vector_to_sun(2) = horiz_component * sin(horiz_angle)
    
    allocate(svf_vecs(3,svf_nrays))
    
    do i=1,svf_nrays
    
      call random_number(azimuth)
      call random_number(elevation)
      elevation = elevation * pi / 2
      azimuth = azimuth * pi * 2
      horiz_component = cos(elevation)

      svf_vecs(1,i) = horiz_component * cos(azimuth)
      svf_vecs(2,i) = horiz_component * sin(azimuth)
      svf_vecs(3,i) = sin(elevation)
    end do
    
  end subroutine
  
  
  
  pure function longwave_radiation(T) result(res)
    real(knd) :: res
    real(knd),intent(in) :: T
    
    res = SB_sigma * T**4
  end function
  
  
end module SolarRadiation

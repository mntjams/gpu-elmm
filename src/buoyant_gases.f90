module BuoyantGases
  use Parameters
  use PhysicalPropertis
  
  ! At the moment buoyant gases should not be combined with air moisture.

  ! Stores information necessary for computation of density
  ! of mixture of gases based on their concantration stored in Scalar.
  ! The concentrations shall mass per volume kg . m^-3.
  
  ! The conversion to mass fraction w = m_g / (m_g + m_d) shall
  ! assume pressure and temperature that would give the reference value of density rho_air_ref.

  ! The array does not have to cover all scalars simulated. Th simulation shall safely proceed
  ! assuming the  restscalars are neutral.
  
  type buoyant_scalar
    real(knd) :: m_mol = m_mol_air_ref
    real(knd) :: R = Rd_air_ref  ! R_gas_universal / m_mol
    real(knd) :: eps = 0         ! m_mol_air_ref/m_mol - 1
  end type
  
  ! Index of the array correspends to the number of the scalar.
  type(buoyant_scalar), allocatable :: buoyant_scalars(:)
  
contains

   
  function buoyant_scalar_init(m_mol) result(res)
    type(buoyant_scalar) :: res
    real(knd), intent(in) :: m_mol
    
    res%m_mol = m_mol
    res%R = R_gas_universal / m_mol
    res%eps = m_mol_air_ref / m_mol - 1
  end function
   
end module


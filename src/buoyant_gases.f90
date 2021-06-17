module BuoyantGases

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
    real(knd) :: R = Rd_air_ref
    real(knd) :: eps = 1
  end type
  
  type(buoyant_scalar), allocatable :: buoyant_scalars(:)
end module


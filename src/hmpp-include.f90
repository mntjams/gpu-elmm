  integer,parameter :: KND=4, TIM=4
  integer,parameter       :: NOSLIP=1, FREESLIP=2, PERIODIC=3, DIRICHLET=4, NEUMANN=5, CONSTFLUX=6,&  !boundary condition types
                               TURBULENTINLET=7, FREESLIPBUFF=8, OUTLETBUFF=9, INLETFROMFILE=10
  integer, parameter :: Ea=1,We=2,So=3,No=4,Bo=5,To=6
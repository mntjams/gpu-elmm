! Marsaglia & Tsang generator for random normals & random exponentials.
! Translated from C by Alan Miller (amiller@bigpond.net.au)

! Marsaglia, G. & Tsang, W.W. (2000) `The ziggurat method for generating
! random variables', J. Statist. Software, v5(8).

! This is an electronic journal which can be downloaded from:
! http://www.jstatsoft.org/v05/i08

! N.B. It is assumed that all integer(int32)s are 32-bit.
! N.B. The value of M2 has been halved to compensate for the lack of
!      unsigned integer(int32)s in Fortran.

! Latest version - 1 January 2001

! Parallel version - October 2006

! This version has been customised for parallel processing use,
! specifically with OpenMP.  Each thread uses its own pseudo-random 
! sequence. (Gib Bogle)

! Changed to subroutines, cleaned tab characters, 
! unified to lower case, new interface (Vladimir Fuka) -  May 2014
!--------------------------------------------------------------------------

module rng_par_zig

   use iso_fortran_env

   implicit none

   private

   public  :: rng_init, rng_int, rng_uni, rng_norm, rng_exp
   
   integer(int32),  parameter  ::  dp = real64
   real(dp), parameter  ::  m1=2147483648.0_dp,   m2=2147483648.0_dp,      &
                            half=0.5_dp
   real(dp)             ::  dn0=3.442619855899_dp, tn0=3.442619855899_dp,    &
                            vn=0.00991256303526217_dp,                     &
                            q,                    de0=7.697117470131487_dp, &
                            te0=7.697117470131487_dp,                       &
                            ve=0.003949659822581572_dp

   integer, save :: par_n = 0, par_step
   integer(int32), allocatable, save ::  par_jsr(:), par_kn(:,:), par_ke(:,:)
   real(dp), allocatable, save :: par_wn(:,:), par_fn(:,:), par_we(:,:), par_fe(:,:)

   interface rng_uni
     module procedure par_uni
     module procedure par_uni_32
     module procedure par_uni_ser
     module procedure par_uni_32_ser
   end interface

   interface rng_norm
     module procedure par_rnor
     module procedure par_rnor_32
     module procedure par_rnor_ser
     module procedure par_rnor_32_ser
   end interface

   interface rng_exp
     module procedure par_rexp
     module procedure par_rexp_32
     module procedure par_rexp_ser
     module procedure par_rexp_32_ser
   end interface

   interface rng_int
     module procedure par_shr3
   end interface
   
   interface rng_init
     module procedure par_zigset
   end interface
   

contains


  subroutine par_zigset( npar, par_jsrseed, grainsize)

    integer(int32), intent(in)  :: npar
    integer(int32), intent(in)  :: par_jsrseed(0:npar-1)
    integer, intent(in), optional  :: grainsize

    integer :: i, kpar
    real(dp) dn, tn, de, te

    par_n = npar

    if (present(grainsize)) then
      par_step = grainsize
    else
      par_step = 32
    end if

    ! First we need to allocate all the non-volatile arrays with the size npar
    allocate(par_jsr(0:npar*par_step))
    allocate(par_kn(0:127,0:npar-1))
    allocate(par_ke(0:255,0:npar-1))
    allocate(par_wn(0:127,0:npar-1))
    allocate(par_fn(0:127,0:npar-1))
    allocate(par_we(0:255,0:npar-1))
    allocate(par_fe(0:255,0:npar-1))

    ! Now treat each instance separately
    do kpar = 0,npar-1
      !  Set the seed
      par_jsr(kpar*par_step) = par_jsrseed(kpar)

      !  Tables for RNOR
      dn = dn0
      tn = tn0
      q = vn*exp(half*dn*dn)
      par_kn(0,kpar) = int((dn/q)*m1, int32)
      par_kn(1,kpar) = 0
      par_wn(0,kpar) = q/m1
      par_wn(127,kpar) = dn/m1
      par_fn(0,kpar) = 1.0_dp
      par_fn(127,kpar) = exp( -half*dn*dn )
      do  i = 126, 1, -1
          dn = sqrt( -2.0_dp * log( vn/dn + exp( -half*dn*dn ) ) )        ! dn
          par_kn(i+1,kpar) = int((dn/tn)*m1, int32)
          tn = dn                                                        ! tn
          par_fn(i,kpar) = exp(-half*dn*dn)
          par_wn(i,kpar) = dn/m1
      end do

      !  Tables for Rexp
      de = de0
      te = te0
      q = ve*exp( de )
      par_ke(0,kpar) = int((de/q)*m2, int32)
      par_ke(1,kpar) = 0
      par_we(0,kpar) = q/m2
      par_we(255,kpar) = de/m2
      par_fe(0,kpar) = 1.0_dp
      par_fe(255,kpar) = exp( -de )
      do  i = 254, 1, -1
          de = -log( ve/de + exp( -de ) )                                ! de
          par_ke(i+1,kpar) = int(m2 * (de/te), int32)
          te = de                                                        ! te
          par_fe(i,kpar) = exp( -de )
          par_we(i,kpar) = de/m2
      end do
    enddo
  end subroutine par_zigset



  !  Generate random 32-bit integers
  subroutine par_shr3(ival, kpar)
    integer(int32), intent(out)  ::  ival
    integer, intent(in)  ::  kpar
    integer(int32) :: jz, jsr
    interface
      function sum_and_overflow(a, b) result(res) bind(C, name="sum_and_overflow")
        use, intrinsic :: iso_c_binding
        integer(c_int32_t) :: res
        integer(c_int32_t), value :: a, b
      end function
    end interface

    jsr = par_jsr(kpar*par_step)
    jz = jsr
    jsr = ieor( jsr, ishft( jsr,  13 ) )
    jsr = ieor( jsr, ishft( jsr, -17 ) )
    jsr = ieor( jsr, ishft( jsr,   5 ) )
    par_jsr(kpar*par_step) = jsr
    !ival = jz + jsr causes overflow of signed integer (against standard, undefined behaviour)
    ival = sum_and_overflow(jz, jsr)
  end subroutine par_shr3



  !  Generate uniformly distributed random numbers, sequence kpar
  subroutine par_uni(fn_val, kpar)
    integer :: kpar
    real(dp), intent(out) ::  fn_val
    integer(int32) :: x

    if (kpar >= par_n) then
        write(*,*) 'thread number exceeds initialized max: ',kpar,par_n-1
        stop
    endif
    
    call par_shr3(x, kpar)
    
    fn_val = half + 0.2328306e-9_dp * x
  end subroutine par_uni



  !  Generate random normals, sequence kpar
  subroutine par_rnor(fn_val, kpar)
    real(dp), intent(out) ::  fn_val
    integer :: kpar

    real(dp), parameter  ::  r = 3.442620_dp
    real(dp)             ::  x, y, z
    integer(int32) :: iz, hz

    if (kpar >= par_n) then
        write(*,*) 'thread number exceeds initialized max: ',kpar,par_n
        stop
    endif
    
    call par_shr3(hz, kpar)
    iz = iand( hz, 127 )
    
    if( abs( hz ) < par_kn(iz,kpar) ) then
        fn_val = hz * par_wn(iz,kpar)
    else
        do
          if( iz == 0 ) then
              do
                call par_uni(z, kpar)
                x = -0.2904764_dp * log( z )
                call par_uni(z, kpar)
                y = -log( z )
                if( y+y >= x*x ) exit
              end do
              fn_val = r+x
              if( hz <= 0 ) fn_val = -fn_val
              return
          end if
          
          x = hz * par_wn(iz,kpar)
          
          call par_uni(z, kpar)
          if( par_fn(iz,kpar) + z*(par_fn(iz-1,kpar)-par_fn(iz,kpar)) < exp(-half*x*x) ) then
              fn_val = x
              return
          end if
          call par_shr3(hz, kpar)
          iz = iand( hz, 127 )
          if( abs( hz ) < par_kn(iz,kpar) ) then
              fn_val = hz * par_wn(iz,kpar)
              return
          end if
        end do
    end if
  end subroutine par_rnor



  !  Generate random exponentials, sequence kpar
  subroutine par_rexp(fn_val, kpar)
    real(dp), intent(out)  ::  fn_val
    integer :: kpar

    real(dp)  ::  x, y
    integer(int32) :: iz, jz

    if (kpar >= par_n) then
        write(*,*) 'thread number exceeds initialized max: ',kpar,par_n-1
        stop
    endif

    call par_shr3(jz, kpar)

    iz = iand( jz, 255 )
    if( abs( jz ) < par_ke(iz,kpar) ) then
        fn_val = abs(jz) * par_we(iz,kpar)
        return
    end if

    do
        if( iz == 0 ) then
          call par_uni(y, kpar)
          fn_val = 7.69711_dp - log( y )
          return
        end if
        x = abs( jz ) * par_we(iz,kpar)
        
        call par_uni(y, kpar)
        
        if( par_fe(iz,kpar) + y * (par_fe(iz-1,kpar) - par_fe(iz,kpar)) < exp( -x ) ) then
          fn_val = x
          return
        end if
        
        call par_shr3(jz, kpar)
        iz = iand( jz, 255 )
        if( abs( jz ) < par_ke(iz,kpar) ) then
          fn_val = abs( jz ) * par_we(iz,kpar)
          return
        end if
    end do
  end subroutine par_rexp

  subroutine par_uni_32(x, kpar)
    real(real32), intent(out) :: x
    integer :: kpar
    real(dp) :: x64
    call par_uni(x64, kpar)
    x = real(x64, real32)
  end subroutine
  
  subroutine par_rnor_32(x, kpar)
    real(real32), intent(out) :: x
    integer :: kpar
    real(dp) :: x64
    call par_rnor(x64, kpar)
    x = real(x64, real32)
  end subroutine
  
  subroutine par_rexp_32(x, kpar)
    real(real32), intent(out) :: x
    integer :: kpar
    real(dp) :: x64
    call par_rexp(x64, kpar)
    x = real(x64, real32)
  end subroutine
  
  subroutine par_uni_ser(x)
    real(dp), intent(out) :: x

    call par_uni(x, 0_int32)
  end subroutine
  
  subroutine par_rnor_ser(x)
    real(dp), intent(out) :: x

    call par_rnor(x, 0_int32)
  end subroutine
  
  subroutine par_rexp_ser(x)
    real(dp), intent(out) :: x

    call par_rexp(x, 0_int32)
  end subroutine
  
  subroutine par_uni_32_ser(x)
    real(real32), intent(out) :: x
    real(dp) :: x64
    call par_uni(x64, 0_int32)
    x = real(x64, real32)
  end subroutine
  
  subroutine par_rnor_32_ser(x)
    real(real32), intent(out) :: x
    real(dp) :: x64
    call par_rnor(x64, 0_int32)
    x = real(x64, real32)
  end subroutine
  
  subroutine par_rexp_32_ser(x)
    real(real32), intent(out) :: x
    real(dp) :: x64
    call par_rexp(x64, 0_int32)
    x = real(x64, real32)
  end subroutine
  
end module rng_par_zig




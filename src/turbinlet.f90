module TurbInlet

  use Parameters
  use ArrayUtilities
  use rng_par_zig
  !$ use omp_lib
  
  implicit none

  private
  public :: turbulence_generator, default_turbulence_generator, GetInletFromFile

  type turbulence_generator
    integer :: direction

    real(knd) :: T_lag
    real(knd) :: L_y
    real(knd) :: L_z

    real(knd) :: Ustar_surf_inlet, &
                 stress_gradient_inlet, &
                 U_ref_inlet, &
                 z_ref_inlet, &
                 z0_inlet, &
                 power_exponent_inlet

    integer :: filtny, filtnz, bigNy, bigNz

    real(knd), allocatable, dimension(:) :: expsy, expsz
    real(knd) :: compat
    real(knd), allocatable, dimension(:,:,:) :: Ru, Rv, Rw !arrays of randoms
    real(knd), allocatable, dimension(:,:,:) :: Psiu, Psiv, Psiw
    real(knd), allocatable, dimension(:,:,:) :: bfilt !filter coefficients (ii,jj,kk,kz)


    real(knd),allocatable,dimension(:)     :: Ustar_inlet !friction velocity profile at inlet
    real(knd),allocatable,dimension(:,:)   :: Uinavg, Vinavg, Winavg !mean values of U,V,W at inflow
    real(knd),allocatable,dimension(:,:,:) :: transform_tensor
    real(knd),dimension(1:3,1:3) :: relative_stress
  contains
    procedure :: init => turbulence_generator_init
    procedure :: time_step => turbulence_generator_time_step
    procedure, private :: init_turbulence_profiles => turbulence_generator_init_turbulence_profiles
    procedure, private :: init_mean_profiles => turbulence_generator_init_mean_profiles
  end type

  type TInlet
    real(knd),allocatable,dimension(:,:) :: U, V, W, temperature
  end type

  type(turbulence_generator) :: default_turbulence_generator

contains

  subroutine turbulence_generator_init(g)
#ifdef PAR
    use custom_par, only: iim, jim, kim, nxims, nyims, nzims, par_co_sum, par_co_min
    use exchange_par
#endif
    class(turbulence_generator), intent(inout) :: g
    integer :: i, j, k
    integer :: jlo, jup, klo, kup
    real(knd) :: bysum, bzsum
    integer :: tid

#ifdef PAR
    !should not happen for 2D decomposition in Y and Z
    if (.not. (iim==1 .or. iim==nxims)) return
    if (nxims > 1 .and. iim==nxims .and. Btype(Ea) /= TurbulentInletType) return
#endif

    call g%init_turbulence_profiles

    call g%init_mean_profiles

#ifdef PAR
    g%filtny = min(max(nint(g%L_y / dymin),1), ceiling(1._knd * gPrny / 3),Prny)
    g%filtnz = min(max(nint(g%L_z / dzmin),1), ceiling(1._knd * gPrnz / 3),Prnz)

    g%filtny = par_co_min(g%filtny)
    g%filtnz = par_co_min(g%filtnz)
#else
    g%filtny = min(max(nint(g%L_y / dymin),1), ceiling(1._knd * Prny / 3))
    g%filtnz = min(max(nint(g%L_z / dzmin),1), ceiling(1._knd * Prnz / 3))
#endif

    g%bigNy = 2 * g%filtny
    g%bigNz = 2 * g%filtnz

    jlo = -g%bigNy + 1
    jup = Prny + g%bigNy

    klo = -g%bigNz + 1
    kup = Prnz + g%bigNz

    allocate(g%Ru(jlo:jup,klo:kup,2))
    allocate(g%Rv(jlo:jup,klo:kup,2))
    allocate(g%Rw(jlo:jup,klo:kup,2))
    allocate(g%Psiu(1:Prny,1:Prnz,1:2))
    allocate(g%Psiv(1:Prny,1:Prnz,1:2))
    allocate(g%Psiw(1:Prny,1:Prnz,1:2))

    !$omp parallel private(j,k,tid)
    tid = 0
    !$ tid = omp_get_thread_num()
    
    !$omp do collapse(2)
    do k = klo, kup
     do j = jlo, jup
       call rng_norm(g%Ru(j,k,1), tid)
       call rng_norm(g%Rv(j,k,1), tid)
       call rng_norm(g%Rw(j,k,1), tid)
     end do
    end do

    !$omp do collapse(2)
    do k = klo, kup
     do j = jlo, jup
       call rng_norm(g%Ru(j,k,2), tid)
       call rng_norm(g%Rv(j,k,2), tid)
       call rng_norm(g%Rw(j,k,2), tid)
     end do
    end do

    !$omp end parallel


    if ((Btype(So)==BC_PERIODIC).or.(Btype(No)==BC_PERIODIC)) then
      !$omp parallel workshare
      forall(k = klo:kup, j = jlo:0)
         g%Ru(j,k,1:2) = g%Ru(j+Prny,k,1:2)
         g%Rv(j,k,1:2) = g%Rv(j+Prny,k,1:2)
         g%Rw(j,k,1:2) = g%Rw(j+Prny,k,1:2)
      end forall
      forall(k = klo:kup,j = Prny+1:jup)
         g%Ru(j,k,1:2) = g%Ru(j-Prny,k,1:2)
         g%Rv(j,k,1:2) = g%Rv(j-Prny,k,1:2)
         g%Rw(j,k,1:2) = g%Rw(j-Prny,k,1:2)
      end forall
      !$omp end  parallel workshare
    end if

    if  ((Btype(Bo)==BC_PERIODIC).or.(Btype(To)==BC_PERIODIC)) then
      !$omp parallel workshare
      forall(k = klo:0, j = jlo:jup)
         g%Ru(j,k,1:2) = g%Ru(j,k+Prnz,1:2)
         g%Rv(j,k,1:2) = g%Rv(j,k+Prnz,1:2)
         g%Rw(j,k,1:2) = g%Rw(j,k+Prnz,1:2)
      end forall
      forall(k = Prnz+1:kup,j = jlo:jup)
         g%Ru(j,k,1:2) = g%Ru(j,k-Prnz,1:2)
         g%Rv(j,k,1:2) = g%Rv(j,k-Prnz,1:2)
         g%Rw(j,k,1:2) = g%Rw(j,k-Prnz,1:2)
      end forall
      !$omp end  parallel workshare
    end if

#ifdef PAR
    do i=1,2
      call par_exchange_boundaries_yz(g%Ru(:,:,i), Prny, Prnz, Btype, &
                                      jlo, klo, g%bigNy, g%bigNz)
      call par_exchange_boundaries_yz(g%Rv(:,:,i), Prny, Prnz, Btype, &
                                      jlo, klo, g%bigNy, g%bigNz)
      call par_exchange_boundaries_yz(g%Rw(:,:,i), Prny, Prnz, Btype, &
                                      jlo, klo, g%bigNy, g%bigNz)
    end do
#endif


    allocate(g%bfilt(-g%bigNy:g%bigNy, -g%bigNz:g%bigNz,1))
    allocate(g%expsy(-g%bigNy:g%bigNy), g%expsz(-g%bigNz:g%bigNz))

    bysum = 0
    do i = -g%bigNy, g%bigNy
     g%expsy(i) = exp(-pi * abs(i) / (g%bigNy))
     bysum = bysum + g%expsy(i)**2
    end do
    bysum = sqrt(bysum)
    g%expsy = g%expsy / bysum


    bzsum = 0
    do i = -g%bigNz, g%bigNz
     g%expsz(i) = exp(-pi * abs(i) / (g%bigNz))
     bzsum = bzsum + g%expsz(i)**2
    end do
    bzsum = sqrt(bzsum)
    g%expsz = g%expsz / bzsum

    !$omp parallel private(j,k)
    !$omp do
    do k = -g%bigNz, g%bigNz
      do j = -g%bigNy, g%bigNy
         g%bfilt(j,k,1) = g%expsy(j) * g%expsz(k)
      end do
    end do
    !$omp end  do
    !$omp workshare
    forall(k = 1:Prnz, j = 1:Prny)
         g%Psiu(j,k,1) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Ru(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,1))
         g%Psiv(j,k,1) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Rv(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,1))
         g%Psiw(j,k,1) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Rw(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,1))
    end forall

    g%compat = sum(g%Uinavg(1:Prny,1:Prnz))
    !$omp end  workshare
    !$omp end  parallel
#ifdef PAR
    g%compat = par_co_sum(g%compat)
#endif

  end subroutine


  subroutine turbulence_generator_time_step(g, Uin, Vin, Win, dt)
#ifdef PAR
    use custom_par, only: iim, jim, kim, nxims, nyims, nzims, par_co_sum, par_co_min
    use exchange_par
#endif
    class(turbulence_generator), intent(inout) :: g
    real(knd), intent(out) :: Uin(-2:,-2:), Vin(-2:,-2:), Win(-2:,-2:)
    real(knd), intent(in) :: dt
    integer :: i, j, k
    integer :: jlo, jup, klo, kup
    real(knd) :: Ui, Vi, Wi, p
    integer :: tid
    
#ifdef PAR
    !should not happen for 2D decomposition in X and Y
    if (.not. (iim==1 .or. iim==nxims)) return
    if (nxims > 1 .and. iim==nxims .and. Btype(Ea) /= TurbulentInletType) return
#endif

    jlo = -g%bigNy + 1
    jup = Prny + g%bigNy

    klo = -g%bigNz + 1
    kup = Prnz + g%bigNz


    !$omp parallel do private(i,j,k,Ui,Vi,Wi)  
    do k = 1, Prnz
      do j = 1, Prny
         g%Psiu(j,k,2) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Ru(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,2))
         g%Psiv(j,k,2) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Rv(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,2))
         g%Psiw(j,k,2) = sum(g%bfilt(-g%bigNy:g%bigNy,-g%bigNz:g%bigNz,1) * &
                             g%Rw(j-g%bigNy:j+g%bigNy,k-g%bigNz:k+g%bigNz,2))
      end  do
    end  do
    !$omp end  parallel do

    call multiply(g%Psiu(:,:,1), exp(-pi * dt / (2._knd * g%T_lag)))
    call add_multiplied(g%Psiu(:,:,1), g%Psiu(:,:,2), sqrt(1 - exp(-pi * dt / (g%T_lag))))

    call multiply(g%Psiv(:,:,1),exp(-pi * dt / (2._knd * g%T_lag)))
    call add_multiplied(g%Psiv(:,:,1), g%Psiv(:,:,2), sqrt(1 - exp(-pi * dt / (g%T_lag))))

    call multiply(g%Psiw(:,:,1),exp(-pi * dt / (2._knd * g%T_lag)))
    call add_multiplied(g%Psiw(:,:,1), g%Psiw(:,:,2), sqrt(1 - exp(-pi * dt / (g%T_lag))))
    
    call set(Uin, 0._knd)
    call set(Vin, 0._knd)
    call set(Win, 0._knd)
    !$omp parallel do private(i,j,k,Ui,Vi,Wi)
    do k = 1, Prnz
     do j = 1, Prny
      Ui = g%Psiu(j,k,1)
      Vi = g%Psiv(j,k,1)
      Wi = g%Psiw(j,k,1)

      Uin(j,k) = g%Uinavg(j,k) + g%transform_tensor(1,j,k) * Ui   !a12,a13,a23 = 0

      Vin(j,k) = g%Vinavg(j,k) + g%transform_tensor(2,j,k) * Ui&
                            +g%transform_tensor(3,j,k) * Vi

      Win(j,k) = g%Winavg(j,k) + g%transform_tensor(4,j,k) * Ui&
                            +g%transform_tensor(5,j,k) * Vi&
                            +g%transform_tensor(6,j,k) * Wi
     end do
    end do
    !$omp end parallel do

    if (g%direction==We .or. g%direction==Ea) then
      !$omp parallel workshare
      p = sum(Uin(1:Prny,1:Prnz))
      !$omp end parallel workshare
    else if (g%direction==So .or. g%direction==No) then
      !$omp parallel workshare
      p = sum(Vin(1:Prny,1:Prnz))
      !$omp end parallel workshare
    end if

    if (g%direction==We .or. g%direction==Ea .or. &
        g%direction==So .or. g%direction==No) then
#ifdef PAR
      p = par_co_sum(p)
#endif
      p = g%compat / p                  !To ensure the g%compatibility condition.

      call multiply(Uin, p)
      call multiply(Vin, p)
      call multiply(Win, p)
    end if


    call BoundUin(1, Uin)

    call BoundUin(2, Vin)

    call BoundUin(3, Win)

    !$omp parallel private(j,k,tid)
    tid = 0
    !$ tid = omp_get_thread_num()
    
    !$omp do collapse(2) 
    do k = klo, kup
     do j = jlo, jup
       call rng_norm(g%Ru(j,k,2), tid)
       call rng_norm(g%Rv(j,k,2), tid)
       call rng_norm(g%Rw(j,k,2), tid)
     end do
    end do
    !$omp end do

    if ((Btype(So)==BC_PERIODIC) .or. (Btype(No)==BC_PERIODIC)) then
      !workaround to a bug in Cray Fortran
      !$omp single
      i = Prny + 1
      !$omp end single

      !$omp do collapse(2)
      do k = klo, kup
        do j = jlo, 0
          g%Ru(j,k,1:2) = g%Ru(j+Prny,k,1:2)
          g%Rv(j,k,1:2) = g%Rv(j+Prny,k,1:2)
          g%Rw(j,k,1:2) = g%Rw(j+Prny,k,1:2)
        end do
      end do
      !$omp end do nowait
      !$omp do collapse(2)
      do k = klo, kup
        do j = i, jup
          g%Ru(j,k,1:2) = g%Ru(j-Prny,k,1:2)
          g%Rv(j,k,1:2) = g%Rv(j-Prny,k,1:2)
          g%Rw(j,k,1:2) = g%Rw(j-Prny,k,1:2)
        end do
      end do
      !$omp end do
    end if
    if  ((Btype(Bo)==BC_PERIODIC) .or. (Btype(To)==BC_PERIODIC)) then
      !$omp do collapse(2)
      do k = klo, 0
        do j = jlo, jup
          g%Ru(j,k,1:2) = g%Ru(j,k+Prnz,1:2)
          g%Rv(j,k,1:2) = g%Rv(j,k+Prnz,1:2)
          g%Rw(j,k,1:2) = g%Rw(j,k+Prnz,1:2)
        end do
      end do
      !$omp end do nowait
      !$omp do collapse(2)
      do k = Prnz+1, kup
        do j = jlo, jup
          g%Ru(j,k,1:2) = g%Ru(j,k-Prnz,1:2)
          g%Rv(j,k,1:2) = g%Rv(j,k-Prnz,1:2)
          g%Rw(j,k,1:2) = g%Rw(j,k-Prnz,1:2)
        end do
      end do
      !$omp end do
    end if
    !$omp end parallel

#ifdef PAR
    do i=1,2
      call par_exchange_boundaries_yz(g%Ru(:,:,i), Prny, Prnz, Btype, &
                                      jlo, klo, g%bigNy, g%bigNz)
      call par_exchange_boundaries_yz(g%Rv(:,:,i), Prny, Prnz, Btype, &
                                      jlo, klo, g%bigNy, g%bigNz)
      call par_exchange_boundaries_yz(g%Rw(:,:,i), Prny, Prnz, Btype, &
                                      jlo, klo, g%bigNy, g%bigNz)
    end do
#endif

  end subroutine turbulence_generator_time_step












  subroutine turbulence_generator_init_turbulence_profiles(g)
    class(turbulence_generator), intent(inout) :: g
    integer :: j,k

    allocate(g%Ustar_inlet(1:Prnz))
    allocate(g%transform_tensor(1:6,1:Prny,1:Prnz))
    
    !constant stress assumption
    if ((profiletype==LOGPROF.and.g%Ustar_surf_inlet<=0).or.(profiletype==POWERPROF.and.g%Ustar_surf_inlet<=0)) then
      g%Ustar_surf_inlet = abs(Karman * g%U_ref_inlet / log(g%z0_inlet / g%z_ref_inlet))
    end if


    do k = 1, Prnz
      g%Ustar_inlet(k) = g%Ustar_surf_inlet * sqrt(max(1 + g%stress_gradient_inlet * zPr(k),1E-5_knd))
    end do


    do k = 1, Prnz
     do j = 1, Prny  ! tt1 = a11,tt2 = a21,tt3 = a22, tt4 = a31, tt5 = a32, tt6 = a33
        g%transform_tensor(1,j,k) = sqrt((g%Ustar_inlet(k)**2) * g%relative_stress(1,1))
        g%transform_tensor(2,j,k) = (g%Ustar_inlet(k)**2) * g%relative_stress(2,1) / g%transform_tensor(1,j,k)
        g%transform_tensor(3,j,k) = sqrt((g%Ustar_inlet(k)**2) * g%relative_stress(2,2) - g%transform_tensor(2,j,k)**2)
        g%transform_tensor(4,j,k) = (g%Ustar_inlet(k)**2) * g%relative_stress(3,1) / g%transform_tensor(1,j,k)
        g%transform_tensor(5,j,k) = ((g%Ustar_inlet(k)**2) * g%relative_stress(3,2)-&
                                 g%transform_tensor(2,j,k) * g%transform_tensor(4,j,k)) / g%transform_tensor(3,j,k)
        g%transform_tensor(6,j,k) = sqrt((g%Ustar_inlet(k)**2) * g%relative_stress(3,3)-&
                                 g%transform_tensor(4,j,k)**2 - g%transform_tensor(5,j,k)**2)

     end do
    end do
  end subroutine

  subroutine turbulence_generator_init_mean_profiles(g)
    class(turbulence_generator), intent(inout) :: g
    real(knd) :: Ustar_prof, utmp
    integer :: k

    allocate(g%Uinavg(-2:Uny+3,-2:Unz+3), g%Vinavg(-2:Vny+3,-2:Vnz+3), g%Winavg(-2:Wny+3,-2:Wnz+3))
    g%Vinavg = 0
    g%Winavg = 0

    if  (profiletype==CONSTPROF) then

      do k = 1, Prnz

        g%Uinavg(:,k) = Uinlet

      end do

    else if (profiletype==LOGPROF) then

      if (g%U_ref_inlet/=0.and.g%z_ref_inlet>0) then

        Ustar_prof = g%U_ref_inlet * Karman / log(g%z_ref_inlet / g%z0_inlet)

        do k = 1, Prnz
          utmp = (Ustar_prof / Karman) * log(zPr(k) / g%z0_inlet)
          if (sign(1._knd,Ustar_prof) * utmp<abs(Ustar_prof) / 2) &
            utmp = (Ustar_prof / 2) * zPr(k) / (g%z0_inlet * 1.22)
          g%Uinavg(:,k) = utmp
        end do

      else

        utmp = (g%Ustar_inlet(1) / Karman) * log(zPr(1) / g%z0_inlet)
        g%Uinavg(:,1) = utmp
        do k = 2, Prnz
          utmp = (g%Ustar_inlet(k) / Karman) * log(zPr(k) / zPr(k-1)) + utmp
          if (sign(1._knd,g%Ustar_inlet(1)) * utmp < abs(g%Ustar_inlet(1)) / 2) &
            utmp = (g%Ustar_inlet(1) / 2) * zPr(k) / (g%z0_inlet * 1.22)
          g%Uinavg(:,k) = utmp
        end do

      end  if

    else if (profiletype==POWERPROF) then

      do k = 1, Prnz
        g%Uinavg(:,k) = g%U_ref_inlet * (zPr(k) / g%z_ref_inlet)**g%power_exponent_inlet
      end do

    else

      g%Uinavg = 0

    end if

    if (windangle/=0) then
        g%Vinavg = g%Uinavg * sin(windangle / 180._knd * pi)
        g%Uinavg = g%Uinavg * cos(windangle / 180._knd * pi)
    end if

  end subroutine turbulence_generator_init_mean_profiles


  subroutine GetInletFromFile(t)
    real(TIM),intent(in):: t
    integer,save:: called = 0
    integer :: Prny2, Prnz2, Vny2, Wnz2
    real(knd) :: dx2
    character(12):: fname
    integer,save:: inletfnum

    type(TInlet),pointer,save:: In1=>null(),In2=>null(),Inp=>null() !pointer to inlets to avoid transfers, time(In1)<time(In2)
    real(TIM),save:: t1,t2 !time if file1, file 2
    real(TIM) :: tp
    integer,save:: io

    real(knd) :: c1,c2

    if (called==0) then
       open(102,file="inletframeinfo.unf",form='unformatted',status='old',action='read',iostat = io)
       if (io/=0) then
        write(*,*) 'Error while opening file inletframeinfo.unf'
       end if
       read(102) Prny2, Prnz2  !for check of consistency of grids before use
       read(102) Vny2
       read(102) Wnz2
       read(102) dx2

       allocate(In1,In2)
       allocate(In1%U(Uny,Unz),In1%V(Vny,Vnz),In1%W(Wny,Wnz))
       allocate(In2%U(Uny,Unz),In2%V(Vny,Vnz),In2%W(Wny,Wnz))
       if (enable_buoyancy) allocate(In1%temperature(Prny,Prnz),In2%temperature(Prny,Prnz))

       if ((Prny/=Prny2).or.(Prnz/=Prnz2).or.(Vny/=Vny2).or.(Wnz/=Wnz2).or.((dx2-dxPr(0))/dx2>0.1)) then
        call error_stop("Mismatch of computational grid and inlet file.")
       end if
       called = 1
       inletfnum = 1

       fname(1:5)="frame"
       write(fname(6:8),"(I3.3)") inletfnum
       fname(9:12)=".unf"

       open(11,file = fname,form='unformatted',status='old',action='read',iostat = io)
       if (io/=0) then
        write(*,*) "Error while opening file ", fname
        call error_stop
       end if
       read(102) t1
       call ReadBC_INLET_FROM_FILE(11,In1)
       close(11)

       inletfnum = inletfnum+1

       fname(1:5)="frame"
       write(fname(6:8),"(I3.3)") inletfnum
       fname(9:12)=".unf"

        open(11,file = fname,form='unformatted',status='old',action='read',iostat = io)
       if (io/=0) then
        write(*,*) "Error while opening file ", fname
        call error_stop
       end if
       read(102) t2
       call ReadBC_INLET_FROM_FILE(11,In2)
       close(11)
    end if

    io = 0

    do while (t2<t.and.io==0)
      read(102,iostat = io) tp
      if (io==0) then
       t1 = t2
       t2 = tp

       inletfnum = inletfnum+1
       fname(1:5)="frame"
       write(fname(6:8),"(I3.3)") inletfnum
       fname(9:12)=".unf"

       open(11,file = fname,form='unformatted',status='old',action='read',iostat = io)
       if (io/=0) then
        write(*,*) "Error while opening file ", fname
        call error_stop
       end if
       Inp=>In1
       In1=>In2
       In2=>Inp
       call ReadBC_INLET_FROM_FILE(11,In2)
       close(11)
      end if
    end do


    if (t<=t1) then
     c1 = 1
     c2 = 0
    else if (t>=t2) then
     c1 = 0
     c2 = 1
    else
     c1=(t2-t)/(t2-t1)
     c2=(t-t1)/(t2-t1)
    end if
    write(*,*) "t",t1,t2,t
    write(*,*) "c",c1,c2

    Uin(1:Uny,1:Unz) = c1*In1%U+c2*In2%U
    Vin(1:Vny,1:Vnz) = c1*In1%V+c2*In2%V
    Win(1:Wny,1:Wnz) = c1*In1%W+c2*In2%W
    if (enable_buoyancy) TempIn(1:Prny,1:Prnz) = c1*In1%temperature+c2*In2%temperature

  end subroutine GetInletFromFile

  subroutine  ReadBC_INLET_FROM_FILE(unitnum,In)
    integer,intent(in):: unitnum
    type(TInlet),intent(inout)::In

    read(unitnum) In%U
    read(unitnum) In%V
    read(unitnum) In%W
    if (enable_buoyancy) then
         read(unitnum) In%temperature
    end if
  end subroutine ReadBC_INLET_FROM_FILE


  subroutine BoundUin(component,Uin)
#ifdef PAR
    use exchange_par, only: par_exchange_boundaries_yz
#endif
    integer,  intent(in)    :: component
    real(knd),intent(inout) :: Uin(-2:,-2:)
    integer :: nx, ny, nz

    if (component==1) then
      nx = Unx
      ny = Uny
      nz = Unz
    else if (component==2) then
      nx = Vnx
      ny = Vny
      nz = Vnz
    else
      nx = Wnx
      ny = Wny
      nz = Wnz
    end if

#ifdef PAR
     call par_exchange_boundaries_yz(Uin, ny, nz, Btype, &
                                     -2, -2, 2, 2)
#endif

    if (Btype(Bo)==BC_DIRICHLET) then
      if (component==3) then
        Uin(1:ny,0) = sideU(component,Bo)
        Uin(1:ny,-1) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,1))
      else
        Uin(1:ny,0) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,1))
        Uin(1:ny,-1) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,2))
      end  if
    else if (Btype(Bo)==BC_PERIODIC) then
        Uin(1:ny,-1:0) = Uin(1:ny,nz-1:nz)
    else if (Btype(Bo)==BC_NOSLIP.or.(component==3.and.Btype(Bo)==BC_FREESLIP)) then
      if (component==3) then
        Uin(1:ny,0) = 0
        Uin(1:ny,-1)=-Uin(1:ny,1)
      else
        Uin(1:ny,-1:0)=-Uin(1:ny,2:1:-1)
      end if
    else if (Btype(Bo)==BC_NEUMANN.or.(component/=3.and.Btype(Bo)==BC_FREESLIP)) then
        Uin(1:ny,0) = Uin(1:ny,1)
        Uin(1:ny,-1) = Uin(1:ny,1)
    else if (Btype(Bo)==BC_PERIODIC) then
        Uin(1:ny,-1:0) = Uin(1:ny,nz-1:nz)
    end if

    if (Btype(To)==BC_DIRICHLET) then
      if (component==3) then
        Uin(1:ny,nz+1) = sideU(component,To)
        Uin(1:ny,nz+2) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz))
      else
        Uin(1:ny,nz+1) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz))
        Uin(1:ny,nz+2) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz-1))
      end  if
    else if (Btype(To)==BC_PERIODIC) then
        Uin(1:ny,nz+1:nz+2) = Uin(1:ny,nz-1:nz)
    else if (Btype(To)==BC_NOSLIP.or.(component==3.and.(Btype(To)==BC_FREESLIP))) then
      if (component==3) then
        Uin(1:ny,nz+1) = 0
        Uin(1:ny,nz+2)=-Uin(1:ny,nz)
      else
        Uin(1:ny,nz+1:nz+2)=-Uin(1:ny,nz:nz-1:-1)
      end if
    else if (Btype(To)==BC_NEUMANN.or.(component/=3.and.(Btype(To)==BC_FREESLIP))) then
        Uin(1:ny,nz+1) = Uin(1:ny,nz)
        Uin(1:ny,nz+2) = Uin(1:ny,nz)
    else if (Btype(To)==BC_PERIODIC) then
        Uin(1:ny,nz+1:nz+2) = Uin(1:ny,1:2)
    end if

    if (Btype(So)==BC_DIRICHLET) then
      if (component==2) then
        Uin(0,-1:nz+2) = sideU(component,So)
        Uin(-1,-1:nz+2) = sideU(component,So)+(sideU(component,So)-Uin(1,-1:nz+2))
      else
        Uin(0,-1:nz+2) = sideU(component,So)+(sideU(component,So)-Uin(1,-1:nz+2))
        Uin(-1,-1:nz+2) = sideU(component,So)+(sideU(component,So)-Uin(2,-1:nz+2))
      end if
    else if (Btype(So)==BC_NOSLIP.or.(component==2.and.Btype(So)==BC_FREESLIP)) then
      if (component==2) then
        Uin(0,-1:nz+2) = 0
        Uin(-1,-1:nz+2)=-Uin(1,-1:nz+2)
      else
        Uin(0,-1:nz+2)=-Uin(1,-1:nz+2)
        Uin(-1,-1:nz+2)=-Uin(2,-1:nz+2)
      end if
    else if (Btype(So)==BC_NEUMANN.or.(component/=2.and.Btype(So)==BC_FREESLIP)) then
        Uin(0,-1:nz+2) = Uin(1,-1:nz+2)
        Uin(-1,-1:nz+2) = Uin(1,-1:nz+2)
    else if (Btype(So)==BC_PERIODIC) then  !Periodic BC
        Uin(0,-1:nz+2) = Uin(ny,-1:nz+2)
        Uin(-1,-1:nz+2) = Uin(ny-1,-1:nz+2)
    end if

    if (Btype(No)==BC_DIRICHLET) then
      if (component==2) then
        Uin(ny+1,-1:nz+2) = sideU(component,No)
        Uin(ny+2,-1:nz+2) = sideU(component,No)+(sideU(component,No)-Uin(ny,-1:nz+2))
      else
        Uin(ny+1,-1:nz+2) = sideU(component,No)+(sideU(component,No)-Uin(ny,-1:nz+2))
        Uin(ny+2,-1:nz+2) = sideU(component,No)+(sideU(component,No)-Uin(ny-1,-1:nz+2))
      end if
    else if (Btype(No)==BC_NOSLIP.or.(component==2.and.Btype(So)==BC_FREESLIP)) then
      if (component==2) then
        Uin(ny+1,-1:nz+2) = 0
        Uin(ny+2,-1:nz+2)=-Uin(ny,-1:nz+2)
      else
        Uin(ny+1:ny+2,-1:nz+2)=-Uin(ny:ny-1:-1,-1:nz+2)
      end if
    else if (Btype(No)==BC_NEUMANN.or.(component/=2.and.Btype(No)==BC_FREESLIP)) then
        Uin(ny+1,-1:nz+2) = Uin(ny,-1:nz+2)
        Uin(ny+2,-1:nz+2) = Uin(ny,-1:nz+2)
    else if (Btype(No)==BC_PERIODIC) then  !Periodic BC
        Uin(ny+1:ny+2,-1:nz+2) = Uin(1:2,-1:nz+2)
    end if


  end  subroutine BoundUin





end  module TurbInlet

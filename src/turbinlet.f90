module TurbInlet

  use Parameters
  use ArrayUtilities
  use rng_par_zig
  !$ use omp_lib
  
  implicit none

  private
  public GetTurbulentInlet, GetInletFromFile, TLag, Lturby, Lturbz, Ustar_inlet, relative_stress, &
         Ustar_surf_inlet, stress_gradient_inlet, U_ref_inlet, z_ref_inlet, z0_inlet, power_exponent_inlet

  real(knd) :: TLag
  real(knd) :: Lturby
  real(knd) :: Lturbz

  real(knd) :: Ustar_surf_inlet, stress_gradient_inlet, U_ref_inlet, z_ref_inlet, z0_inlet, power_exponent_inlet


  real(knd),allocatable,dimension(:)     :: Ustar_inlet !friction velocity profile at inlet
  real(knd),allocatable,dimension(:,:)   :: Uinavg,Vinavg,Winavg !mean values of U,V,W at inflow
  real(knd),allocatable,dimension(:,:,:) :: transform_tensor
  real(knd),dimension(1:3,1:3) :: relative_stress


   type TInlet
     real(knd),allocatable,dimension(:,:) :: U,V,W,temperature
   end type

  interface GetTurbulentInlet
    module procedure GetTurbInletXie
  end interface

contains



 subroutine GetTurbInletXie(dt)
#ifdef PAR
   use custom_par, only: iim, jim, kim, nxims, nyims, nzims, par_co_sum, par_co_min
   use exchange_par
#endif
   real(knd),intent(in) :: dt
   logical,save:: called=.false.
   integer i,j,k
   integer,save:: filtny,filtnz,bigNy,bigNz
   real(knd) Ui,Vi,Wi,bysum,bzsum,p
   real(knd),allocatable,dimension(:):: expsy,expsz
   real(knd),save:: compat
   real(knd),allocatable,dimension(:,:,:),save :: Ru,Rv,Rw !arrays of randoms
   real(knd),allocatable,dimension(:,:,:),save :: Psiu,Psiv,Psiw
   real(knd),allocatable,dimension(:,:,:),save :: bfilt !filter coefficients (ii,jj,kk,kz)
   integer :: tid
   
#ifdef PAR
   !should not happen for 2D decomposition in X and Y
   if (.not.(iim==1.or.iim==nxims)) return
   if (nxims>1.and.iim==nxims.and.Btype(Ea)/=TurbulentInletType) return
#endif

   if (.not. called) then

      call InitTurbulenceProfiles

      call InitMeanProfiles

#ifdef PAR
      filtny = min(max(nint(lturby/dymin),1),ceiling(1._knd*gPrny/3),Prny)
      filtnz = min(max(nint(lturbz/dzmin),1),ceiling(1._knd*gPrnz/3),Prnz)

      filtny = par_co_min(filtny)
      filtnz = par_co_min(filtnz)
#else
      filtny = min(max(nint(lturby/dymin),1),ceiling(1._knd*Prny/3))
      filtnz = min(max(nint(lturbz/dzmin),1),ceiling(1._knd*Prnz/3))
#endif

      bigNy = 2*filtny
      bigNz = 2*filtnz


      allocate(Ru(-bigNy+1:Prny+bigNy,-bigNz+1:Prnz+bigNz,2))
      allocate(Rv(-bigNy+1:Prny+bigNy,-bigNz+1:Prnz+bigNz,2))
      allocate(Rw(-bigNy+1:Prny+bigNy,-bigNz+1:Prnz+bigNz,2))
      allocate(Psiu(1:Prny,1:Prnz,1:2))
      allocate(Psiv(1:Prny,1:Prnz,1:2))
      allocate(Psiw(1:Prny,1:Prnz,1:2))

      !$omp parallel private(j,k,tid)
      tid = 0
      !$ tid = omp_get_thread_num()
      
      !$omp do collapse(2)
      do k = -bigNz+1, Prnz+bigNz
       do j = -bigNy+1, Prny+bigNy
         call rng_norm(Ru(j,k,1), tid)
         call rng_norm(Rv(j,k,1), tid)
         call rng_norm(Rw(j,k,1), tid)
       end do
      end do

      !$omp do collapse(2)
      do k = -bigNz+1, Prnz+bigNz
       do j = -bigNy+1, Prny+bigNy
         call rng_norm(Ru(j,k,2), tid)
         call rng_norm(Rv(j,k,2), tid)
         call rng_norm(Rw(j,k,2), tid)
       end do
      end do
 
      !$omp end parallel


      if ((Btype(So)==PERIODIC).or.(Btype(No)==PERIODIC)) then
        !$omp parallel workshare
        forall(k = -bigNz+1:Prnz+bigNz,j = -bigNy+1:0)
           Ru(j,k,1:2) = Ru(j+Prny,k,1:2)
           Rv(j,k,1:2) = Rv(j+Prny,k,1:2)
           Rw(j,k,1:2) = Rw(j+Prny,k,1:2)
        end forall
        forall(k = -bigNz+1:Prnz+bigNz,j = Prny+1:Prny+bigNy)
           Ru(j,k,1:2) = Ru(j-Prny,k,1:2)
           Rv(j,k,1:2) = Rv(j-Prny,k,1:2)
           Rw(j,k,1:2) = Rw(j-Prny,k,1:2)
        end forall
        !$omp end  parallel workshare
      end if

      if  ((Btype(Bo)==PERIODIC).or.(Btype(To)==PERIODIC)) then
        !$omp parallel workshare
        forall(k = -bigNz+1:0,j = -bigNy+1:Prny+bigNy)
           Ru(j,k,1:2) = Ru(j,k+Prnz,1:2)
           Rv(j,k,1:2) = Rv(j,k+Prnz,1:2)
           Rw(j,k,1:2) = Rw(j,k+Prnz,1:2)
        end forall
        forall(k = Prnz+1:Prnz+bigNz,j = -bigNy+1:Prny+bigNy)
           Ru(j,k,1:2) = Ru(j,k-Prnz,1:2)
           Rv(j,k,1:2) = Rv(j,k-Prnz,1:2)
           Rw(j,k,1:2) = Rw(j,k-Prnz,1:2)
        end forall
        !$omp end  parallel workshare
      end if

#ifdef PAR
      do i=1,2
        call par_exchange_boundaries_yz(Ru(:,:,i), Prny, Prnz, Btype, &
                                        -bigNy+1, -bigNz+1, bigNy, bigNz)
        call par_exchange_boundaries_yz(Rv(:,:,i), Prny, Prnz, Btype, &
                                        -bigNy+1, -bigNz+1, bigNy, bigNz)
        call par_exchange_boundaries_yz(Rw(:,:,i), Prny, Prnz, Btype, &
                                        -bigNy+1, -bigNz+1, bigNy, bigNz)
      end do
#endif


      allocate(bfilt(-bigNy:bigNy,-bigNz:bigNz,1))
      allocate(expsy(-bigNy:bigNy),expsz(-bigNz:bigNz))

      bysum = 0
      do i = -bigNy,bigNy
       expsy(i) = exp(-pi*abs(i)/(bigNy))
       bysum = bysum+expsy(i)**2
      end do
      bysum = sqrt(bysum)
      expsy = expsy/bysum


      bzsum = 0
      do i = -bigNz,bigNz
       expsz(i) = exp(-pi*abs(i)/(bigNz))
       bzsum = bzsum+expsz(i)**2
      end do
      bzsum = sqrt(bzsum)
      expsz = expsz/bzsum

      !$omp parallel private(j,k)
      !$omp do
      do k = -bigNz,bigNz
        do j = -bigNy,bigNy
           bfilt(j,k,1) = expsy(j)*expsz(k)
        end do
      end do
      !$omp end  do
      !$omp workshare
      forall(k = 1:Prnz,j = 1:Prny)
           Psiu(j,k,1) = sum(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Ru(j-bigNy:j+bigNy,k-bigNz:k+bigNz,1))
           Psiv(j,k,1) = sum(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Rv(j-bigNy:j+bigNy,k-bigNz:k+bigNz,1))
           Psiw(j,k,1) = sum(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Rw(j-bigNy:j+bigNy,k-bigNz:k+bigNz,1))
      end forall

      compat = sum(Uinavg(1:Prny,1:Prnz))
      !$omp end  workshare
      !$omp end  parallel
#ifdef PAR
      compat = par_co_sum(compat)
#endif

      called=.true.

   end if

   !$omp parallel do private(i,j,k,Ui,Vi,Wi)  
   do k = 1, Prnz
     do j = 1, Prny
        Psiu(j,k,2) = sum(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Ru(j-bigNy:j+bigNy,k-bigNz:k+bigNz,2))
        Psiv(j,k,2) = sum(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Rv(j-bigNy:j+bigNy,k-bigNz:k+bigNz,2))
        Psiw(j,k,2) = sum(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Rw(j-bigNy:j+bigNy,k-bigNz:k+bigNz,2))
     end  do
   end  do
   !$omp end  parallel do

   call multiply(Psiu(:,:,1), exp(-pi*dt/(2._knd*TLag)))
   call add_multiplied(Psiu(:,:,1), Psiu(:,:,2), sqrt(1-exp(-pi*dt/(TLag))))

   call multiply(Psiv(:,:,1),exp(-pi*dt/(2._knd*TLag)))
   call add_multiplied(Psiv(:,:,1), Psiv(:,:,2), sqrt(1-exp(-pi*dt/(TLag))))

   call multiply(Psiw(:,:,1),exp(-pi*dt/(2._knd*TLag)))
   call add_multiplied(Psiw(:,:,1), Psiw(:,:,2), sqrt(1-exp(-pi*dt/(TLag))))
   
   call set(Uin,0._knd)
   call set(Vin,0._knd)
   call set(Win,0._knd)
   !$omp parallel do private(i,j,k,Ui,Vi,Wi)
   do k = 1, Prnz
    do j = 1, Prny
     Ui = Psiu(j,k,1)
     Vi = Psiv(j,k,1)
     Wi = Psiw(j,k,1)

     Uin(j,k) = Uinavg(j,k)+transform_tensor(1,j,k)*Ui   !a12,a13,a23 = 0

     Vin(j,k) = Vinavg(j,k)+transform_tensor(2,j,k)*Ui&
                           +transform_tensor(3,j,k)*Vi

     Win(j,k) = Winavg(j,k)+transform_tensor(4,j,k)*Ui&
                           +transform_tensor(5,j,k)*Vi&
                           +transform_tensor(6,j,k)*Wi
    end do
   end do
   !$omp end parallel do

   if (Btype(We)==TURBULENTINLET .or. Btype(Ea)==TURBULENTINLET) then
     !$omp parallel workshare
     p = sum(Uin(1:Prny,1:Prnz))
     !$omp end parallel workshare
   else if (Btype(So)==TURBULENTINLET .or. Btype(No)==TURBULENTINLET) then
     !$omp parallel workshare
     p = sum(Vin(1:Prny,1:Prnz))
     !$omp end parallel workshare
   end if

   if (Btype(We)==TURBULENTINLET .or. Btype(Ea)==TURBULENTINLET .or. &
       Btype(So)==TURBULENTINLET .or. Btype(No)==TURBULENTINLET) then
#ifdef PAR
     p = par_co_sum(p)
#endif
     p = compat/p                  !To ensure the compatibility condition.

     call multiply(Uin,p)
     call multiply(Vin,p)
     call multiply(Win,p)
   end if


   call BoundUin(1,Uin)

   call BoundUin(2,Vin)

   call BoundUin(3,Win)

   !$omp parallel private(j,k,tid)
   tid = 0
   !$ tid = omp_get_thread_num()
   
   !$omp do collapse(2) 
   do k = -bigNz+1, Prnz+bigNz
    do j = -bigNy+1, Prny+bigNy
      call rng_norm(Ru(j,k,2), tid)
      call rng_norm(Rv(j,k,2), tid)
      call rng_norm(Rw(j,k,2), tid)
    end do
   end do
   !$omp end do

   if ((Btype(So)==PERIODIC).or.(Btype(No)==PERIODIC)) then
     !$omp do collapse(2)
     do k = -bigNz+1, Prnz+bigNz
       do j = -bigNy+1, 0
         Ru(j,k,1:2) = Ru(j+Prny,k,1:2)
         Rv(j,k,1:2) = Rv(j+Prny,k,1:2)
         Rw(j,k,1:2) = Rw(j+Prny,k,1:2)
       end do
     end do
     !$omp end do nowait
     !$omp do collapse(2)
     do k = -bigNz+1, Prnz+bigNz
       do j = Prny+1, Prny+bigNy
         Ru(j,k,1:2) = Ru(j-Prny,k,1:2)
         Rv(j,k,1:2) = Rv(j-Prny,k,1:2)
         Rw(j,k,1:2) = Rw(j-Prny,k,1:2)
       end do
     end do
     !$omp end do
   end if
   if  ((Btype(Bo)==PERIODIC).or.(Btype(To)==PERIODIC)) then
     !$omp do collapse(2)
     do k = -bigNz+1, 0
       do j = -bigNy+1, Prny+bigNy
         Ru(j,k,1:2) = Ru(j,k+Prnz,1:2)
         Rv(j,k,1:2) = Rv(j,k+Prnz,1:2)
         Rw(j,k,1:2) = Rw(j,k+Prnz,1:2)
       end do
     end do
     !$omp end do nowait
     !$omp do collapse(2)
     do k = Prnz+1, Prnz+bigNz
       do j = -bigNy+1, Prny+bigNy
         Ru(j,k,1:2) = Ru(j,k-Prnz,1:2)
         Rv(j,k,1:2) = Rv(j,k-Prnz,1:2)
         Rw(j,k,1:2) = Rw(j,k-Prnz,1:2)
       end do
     end do
     !$omp end do
   end if
   !$omp end parallel

#ifdef PAR
   do i=1,2
     call par_exchange_boundaries_yz(Ru(:,:,i), Prny, Prnz, Btype, &
                                     -bigNy+1, -bigNz+1, bigNy, bigNz)
     call par_exchange_boundaries_yz(Rv(:,:,i), Prny, Prnz, Btype, &
                                     -bigNy+1, -bigNz+1, bigNy, bigNz)
     call par_exchange_boundaries_yz(Rw(:,:,i), Prny, Prnz, Btype, &
                                     -bigNy+1, -bigNz+1, bigNy, bigNz)
   end do
#endif

 end subroutine GetTurbInletXie












  subroutine InitTurbulenceProfiles
    integer j,k

    allocate(Ustar_inlet(1:Prnz))
    allocate(transform_tensor(1:6,1:Prny,1:Prnz))
    
    !constant stress assumption
    if ((profiletype==LOGPROF.and.Ustar_surf_inlet<=0).or.(profiletype==POWERPROF.and.Ustar_surf_inlet<=0)) then
      Ustar_surf_inlet = abs(Karman*U_ref_inlet/log(z0_inlet/z_ref_inlet))
    end if


    do k = 1, Prnz
      Ustar_inlet(k) = Ustar_surf_inlet*sqrt(max(1+stress_gradient_inlet*zPr(k),1E-5_knd))
    end do


    do k = 1, Prnz
     do j = 1, Prny  ! tt1 = a11,tt2 = a21,tt3 = a22, tt4 = a31, tt5 = a32, tt6 = a33
        transform_tensor(1,j,k) = sqrt((Ustar_inlet(k)**2)*relative_stress(1,1))
        transform_tensor(2,j,k) = (Ustar_inlet(k)**2)*relative_stress(2,1)/transform_tensor(1,j,k)
        transform_tensor(3,j,k) = sqrt((Ustar_inlet(k)**2)*relative_stress(2,2)-transform_tensor(2,j,k)**2)
        transform_tensor(4,j,k) = (Ustar_inlet(k)**2)*relative_stress(3,1)/transform_tensor(1,j,k)
        transform_tensor(5,j,k) = ((Ustar_inlet(k)**2)*relative_stress(3,2)-&
                                 transform_tensor(2,j,k)*transform_tensor(4,j,k))/transform_tensor(3,j,k)
        transform_tensor(6,j,k) = sqrt((Ustar_inlet(k)**2)*relative_stress(3,3)-&
                                 transform_tensor(4,j,k)**2-transform_tensor(5,j,k)**2)

     end do
    end do
  end subroutine InitTurbulenceProfiles

  subroutine InitMeanProfiles
    real(knd) Ustar_prof
    integer k

    allocate(Uinavg(-2:Uny+3,-2:Unz+3),Vinavg(-2:Vny+3,-2:Vnz+3),Winavg(-2:Wny+3,-2:Wnz+3))
    Vinavg = 0
    Winavg = 0

    if  (profiletype==CONSTPROF) then

      do k = 1, Prnz

        Uinavg(:,k) = Uinlet

      end do

    else if (profiletype==LOGPROF) then

      if (U_ref_inlet/=0.and.z_ref_inlet>0) then

        Ustar_prof = U_ref_inlet * Karman / log(z_ref_inlet/z0_inlet)
        do k = 1, Prnz
          Uinavg(:,k) = (Ustar_prof/Karman)*log(zPr(k)/z0_inlet)
        end do

      else

        Uinavg(:,1) = (Ustar_inlet(1)/Karman)*log(zPr(1)/z0_inlet)
        do k = 2, Prnz
          Uinavg(:,k) = (Ustar_inlet(k)/Karman)*log(zPr(k)/zPr(k-1)) + Uinavg(:,k-1)
        end do

      end  if

    else if (profiletype==POWERPROF) then

      do k = 1, Prnz
        Uinavg(:,k) = U_ref_inlet*(zPr(k)/z_ref_inlet)**power_exponent_inlet
      end do

    else

      Uinavg = 0

    end if

    if (windangle/=0) then
        Vinavg = Uinavg * sin(windangle/180._knd*pi)
        Uinavg = Uinavg * cos(windangle/180._knd*pi)
    end if

  end subroutine InitMeanProfiles


  subroutine GetInletFromFile(t)
    real(TIM),intent(in):: t
    integer,save:: called = 0
    integer Prny2, Prnz2, Vny2, Wnz2
    real(knd) dx2
    character(12):: fname
    integer,save:: inletfnum

    type(TInlet),pointer,save:: In1=>null(),In2=>null(),Inp=>null() !pointer to inlets to avoid transfers, time(In1)<time(In2)
    real(TIM),save:: t1,t2 !time if file1, file 2
    real(TIM) tp
    integer,save:: io

    real(knd) c1,c2

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
       call ReadInletFromFile(11,In1)
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
       call ReadInletFromFile(11,In2)
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
       call ReadInletFromFile(11,In2)
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

  subroutine  ReadInletFromFile(unitnum,In)
    integer,intent(in):: unitnum
    type(TInlet),intent(inout)::In

    read(unitnum) In%U
    read(unitnum) In%V
    read(unitnum) In%W
    if (enable_buoyancy) then
         read(unitnum) In%temperature
    end if
  end subroutine ReadInletFromFile


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

    if (Btype(Bo)==DIRICHLET) then
      if (component==3) then
        Uin(1:ny,0) = sideU(component,Bo)
        Uin(1:ny,-1) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,1))
      else
        Uin(1:ny,0) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,1))
        Uin(1:ny,-1) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,2))
      end  if
    else if (Btype(Bo)==PERIODIC) then
        Uin(1:ny,-1:0) = Uin(1:ny,nz-1:nz)
    else if (Btype(Bo)==NOSLIP.or.(component==3.and.Btype(Bo)==FREESLIP)) then
      if (component==3) then
        Uin(1:ny,0) = 0
        Uin(1:ny,-1)=-Uin(1:ny,1)
      else
        Uin(1:ny,-1:0)=-Uin(1:ny,2:1:-1)
      end if
    else if (Btype(Bo)==NEUMANN.or.(component/=3.and.Btype(Bo)==FREESLIP)) then
        Uin(1:ny,0) = Uin(1:ny,1)
        Uin(1:ny,-1) = Uin(1:ny,1)
    else if (Btype(Bo)==PERIODIC) then
        Uin(1:ny,-1:0) = Uin(1:ny,nz-1:nz)
    end if

    if (Btype(To)==DIRICHLET) then
      if (component==3) then
        Uin(1:ny,nz+1) = sideU(component,To)
        Uin(1:ny,nz+2) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz))
      else
        Uin(1:ny,nz+1) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz))
        Uin(1:ny,nz+2) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz-1))
      end  if
    else if (Btype(To)==PERIODIC) then
        Uin(1:ny,nz+1:nz+2) = Uin(1:ny,nz-1:nz)
    else if (Btype(To)==NOSLIP.or.(component==3.and.(Btype(To)==FREESLIP))) then
      if (component==3) then
        Uin(1:ny,nz+1) = 0
        Uin(1:ny,nz+2)=-Uin(1:ny,nz)
      else
        Uin(1:ny,nz+1:nz+2)=-Uin(1:ny,nz:nz-1:-1)
      end if
    else if (Btype(To)==NEUMANN.or.(component/=3.and.(Btype(To)==FREESLIP))) then
        Uin(1:ny,nz+1) = Uin(1:ny,nz)
        Uin(1:ny,nz+2) = Uin(1:ny,nz)
    else if (Btype(To)==PERIODIC) then
        Uin(1:ny,nz+1:nz+2) = Uin(1:ny,1:2)
    end if

    if (Btype(So)==DIRICHLET) then
      if (component==2) then
        Uin(0,-1:nz+2) = sideU(component,So)
        Uin(-1,-1:nz+2) = sideU(component,So)+(sideU(component,So)-Uin(1,-1:nz+2))
      else
        Uin(0,-1:nz+2) = sideU(component,So)+(sideU(component,So)-Uin(1,-1:nz+2))
        Uin(-1,-1:nz+2) = sideU(component,So)+(sideU(component,So)-Uin(2,-1:nz+2))
      end if
    else if (Btype(So)==NOSLIP.or.(component==2.and.Btype(So)==FREESLIP)) then
      if (component==2) then
        Uin(0,-1:nz+2) = 0
        Uin(-1,-1:nz+2)=-Uin(1,-1:nz+2)
      else
        Uin(0,-1:nz+2)=-Uin(1,-1:nz+2)
        Uin(-1,-1:nz+2)=-Uin(2,-1:nz+2)
      end if
    else if (Btype(So)==NEUMANN.or.(component/=2.and.Btype(So)==FREESLIP)) then
        Uin(0,-1:nz+2) = Uin(1,-1:nz+2)
        Uin(-1,-1:nz+2) = Uin(1,-1:nz+2)
    else if (Btype(So)==PERIODIC) then  !Periodic BC
        Uin(0,-1:nz+2) = Uin(ny,-1:nz+2)
        Uin(-1,-1:nz+2) = Uin(ny-1,-1:nz+2)
    end if

    if (Btype(No)==DIRICHLET) then
      if (component==2) then
        Uin(ny+1,-1:nz+2) = sideU(component,No)
        Uin(ny+2,-1:nz+2) = sideU(component,No)+(sideU(component,No)-Uin(ny,-1:nz+2))
      else
        Uin(ny+1,-1:nz+2) = sideU(component,No)+(sideU(component,No)-Uin(ny,-1:nz+2))
        Uin(ny+2,-1:nz+2) = sideU(component,No)+(sideU(component,No)-Uin(ny-1,-1:nz+2))
      end if
    else if (Btype(No)==NOSLIP.or.(component==2.and.Btype(So)==FREESLIP)) then
      if (component==2) then
        Uin(ny+1,-1:nz+2) = 0
        Uin(ny+2,-1:nz+2)=-Uin(ny,-1:nz+2)
      else
        Uin(ny+1:ny+2,-1:nz+2)=-Uin(ny:ny-1:-1,-1:nz+2)
      end if
    else if (Btype(No)==NEUMANN.or.(component/=2.and.Btype(No)==FREESLIP)) then
        Uin(ny+1,-1:nz+2) = Uin(ny,-1:nz+2)
        Uin(ny+2,-1:nz+2) = Uin(ny,-1:nz+2)
    else if (Btype(No)==PERIODIC) then  !Periodic BC
        Uin(ny+1:ny+2,-1:nz+2) = Uin(1:2,-1:nz+2)
    end if


  end  subroutine BoundUin





end  module TurbInlet

module POISSON

 use PARAMETERS
 use BOUNDARIES
 use MULTIGRID,   only: POISSMG
 use MULTIGRID2d, only: POISSMG2d
#ifdef __HMPP
 use HMPP_CODELETS
#endif

implicit none

 private
 public PR_CORRECT

#ifdef __HMPP
 public GetPrFromGPU
#endif

contains


  subroutine PR_CORRECT(U,V,W,Pr,Q,coef)                       !Pressure correction
    real(knd),dimension(-2:,-2:,-2:),intent(inout)    :: U,V,W !Phi is computed in Poisson eq. with div of U in RHS
    real(knd),dimension(1:,1:,1:),intent(inout)       :: Pr    !Depend ing on active projection method Phi becomes new pressure
    real(knd),dimension(:,:,:),allocatable,intent(in) :: Q     !or is added to last pressure
    real(knd),intent(in) :: coef
                                                           !U,V,W velocity field for correction
    real(knd),save,allocatable :: Phi(:,:,:), RHS(:,:,:)   !Pr pressure
                                                           !coef cofficient from Runge Kutta, Q mass sources from immersed boundary
    real(TIM) dt2,dt3                                      !RHS right hand side of eq. with divergence of U
                                                           !Phi computed pseudopressure, saved as first guess for next time
    real(knd) :: divergence,uncompatibility
    character(70) :: str
    integer,save :: called = 0
    integer(DBL), save :: trate
    integer(DBL), save :: time1, time2, time3, time4
    integer(DBL), save :: timem1, timem2, timem3, timem4
correctcompatibility = 2
    if (called==0) then
      allocate(Phi(0:Prnx+1,0:Prny+1,0:Prnz+1))
      allocate(RHS(0:Prnx+1,0:Prny+1,0:Prnz+1))
      Phi = 0

      if (debugparam>1) call system_clock(count_rate=trate)
    end if
    called = called + 1


    if (debugparam>1.and.called>1) call system_clock(count=time1)

    write (*,*) "Preparing pressure correction"


    dt2 = coef*dt
    if (Re>0)  then
      dt3 = coef*dt / (2._TIM*Re)
    else
      dt3 = 0
    end if

#ifdef __HMPP
    if (GPU>0) then !If GPU already allocated. typically bot the case on first call.

      !$hmpp <tsteps> advancedload, args[PrePoisson::dt2,PrePoisson::correctcompatibility]

      !$hmpp <tsteps> PrePoisson callsite, args[*].noupdate=true
      call PrePoisson_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,correctcompatibility,&
                          Btype,sideU,&
                          dt2,dxmin,dymin,dzmin,&
                          Uin,Vin,Win,U,V,W,RHS,uncompatibility,divergence)

      if (correctcompatibility==1) then
        !$hmpp <tsteps> delegatedstore, args[PrePoisson::uncompatibility]
      end if
      if (debugparam>1.and.called>1) call system_clock(count=timem1)
      !$hmpp <tsteps> delegatedstore, args[PrePoisson::divergence,PrePoisson::RHS]
      if (debugparam>1.and.called>1) call system_clock(count=timem2)
    else
      call PrePoisson(U,V,W,RHS,dt2,uncompatibility,divergence)
    end if
#else
    call PrePoisson(U,V,W,Q,RHS,dt2,uncompatibility,divergence)
#endif

    if (debugparam>1.and.called>1) call system_clock(count=time2)

    if (correctcompatibility>=1) write(*,*) "Uncompatibility:",uncompatibility

    write(*,*) "avgRHS:",divergence


    if (poissmet==1) then

        call POISSSOR(Phi,RHS)

    else if (poissmet==2) then

        call POISS_POISFFT(Phi,RHS)

    else if (poissmet==3) then

        if (Prny==1) then
          call POISSMG2d(Phi,RHS)
        else
          call POISSMG(Phi,RHS)
        end if

    end if


    call Bound_Phi(Phi)


    if (debugparam>1.and.called>1) then
     call system_clock(count=time3)
     write (*,*) "ET of part 2", real(time3-time2)/real(trate)
    endif

#ifdef __HMPP
    if (GPU>0) then
      if (debugparam>1.and.called>1) call system_clock(count=timem3)
      !$hmpp <tsteps> advancedload, args[PostPoisson::dt3,PostPoisson::Phi]
      if (debugparam>1.and.called>1) then
       call system_clock(count=timem4)
       write (*,*) "ET of part 4", real((timem4-timem3)+(timem2-timem1))/real(trate)
      endif

      !$hmpp <tsteps> PostPoisson callsite, args[*].noupdate=true
      call PostPoisson_GPU(Unx,Uny,Unz,Vnx,Vny,Vnz,Wnx,Wny,Wnz,Prnx,Prny,Prnz,&
                           Btype,sideU,&
                           dt2,dt3,dxmin,dymin,dzmin,&
                           Uin,Vin,Win,U,V,W,Pr,Phi)
    else
      call PostPoisson(U,V,W,Pr,Phi,dt2,dt3)
    end if
#else
    call PostPoisson(U,V,W,Pr,Phi,dt2,dt3)
#endif
    if (debugparam>1.and.called>1) then
     call system_clock(count=time4)
     write (*,*) "ET of part 3", real((time4-time1)-(time3-time2))/real(trate)
    endif
  end subroutine PR_CORRECT



#ifdef __HMPP
  subroutine GetPrFromGPU(Pr)
    real(knd),dimension(1:,1:,1:),intent(in) :: Pr
    !$hmpp <tsteps> delegatedstore, args[PostPoisson::Pr]
  end subroutine GetPrFromGPU
#endif





  subroutine PrePoisson(U,V,W,Q,RHS,dt2,uncompatibility,divergence)
    real(knd),intent(inout) :: U(-2:,-2:,-2:)
    real(knd),intent(inout) :: V(-2:,-2:,-2:)
    real(knd),intent(inout) :: W(-2:,-2:,-2:)
    real(knd),intent(out)   :: RHS(0:,0:,0:)
    real(knd),allocatable,intent(in) :: Q(:,:,:)
    real(knd),intent(out)   :: uncompatibility,divergence
    real(knd),intent(in)    :: dt2
    real(knd) :: S,S2
    integer   :: i,j,k

    !$omp parallel
    !$omp sections
    !$omp section
    call BoundU(1,U,Uin)
    !$omp section
    call BoundU(2,V,Vin)
    !$omp section
    call BoundU(3,W,Win)
    !$omp end sections
    !$omp end parallel

    if (correctcompatibility>=1) then
      S=0
      !$omp parallel do private(i,j,k) reduction(+:S)
      do k=1,Prnz
       do j=1,Prny
        do i=1,Prnx
         S = S + (-((U(i,j,k)-U(i-1,j,k))/(dxmin)+(V(i,j,k)-V(i,j-1,k))/(dymin)+(W(i,j,k)-W(i,j,k-1))&
                      /(dzmin)))
        end do
       end do
      end do
      !$omp end parallel do
      
      if (abs(windangle-90)<1.or.abs(windangle+90)<1) then
        S=S*dymin/(Prnx*Prnz)

        uncompatibility = S
        
        if (correctcompatibility==1) then
          !$omp parallel workshare
          V(:,Vny+1,:) = V(:,Vny+1,:)+S
          !$omp end parallel workshare
        end if
      else
        S=S*dxmin/(Prny*Prnz)

        uncompatibility = S
        
        if (correctcompatibility==1) then
          !$omp parallel workshare
          U(Unx+1,:,:) = U(Unx+1,:,:)+S
          !$omp end parallel workshare
        end if
      end if
    end if

    S=0
    S2=0

    if (allocated(Q)) then
      !$omp parallel do private(i,j,k) reduction(+:s2)
      do k=1,Prnz            !divergence of U -> RHS
       do j=1,Prny
        do i=1,Prnx
             RHS(i,j,k) = (U(i,j,k)-U(i-1,j,k))/(dxmin)&
                         +(V(i,j,k)-V(i,j-1,k))/(dymin)&
                         +(W(i,j,k)-W(i,j,k-1))/(dzmin)&
                         -Q(i,j,k)
             S2=S2+abs(RHS(i,j,k))
             RHS(i,j,k) = RHS(i,j,k)/(dt2)
        end do
       end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(i,j,k) reduction(+:s2)
      do k=1,Prnz            !divergence of U -> RHS
       do j=1,Prny
        do i=1,Prnx
             RHS(i,j,k) = (U(i,j,k)-U(i-1,j,k))/(dxmin)&
                         +(V(i,j,k)-V(i,j-1,k))/(dymin)&
                         +(W(i,j,k)-W(i,j,k-1))/(dzmin)
             S2=S2+abs(RHS(i,j,k))
             RHS(i,j,k) = RHS(i,j,k)/(dt2)
        end do
       end do
      end do
      !$omp end parallel do
    end if

    divergence = S2/(Prnx*Prny*Prnz)

   end subroutine PrePoisson




  subroutine PostPoisson(U,V,W,Pr,Phi,dt2,dt3)
    real(knd),intent(inout) :: U(-2:,-2:,-2:)
    real(knd),intent(inout) :: V(-2:,-2:,-2:)
    real(knd),intent(inout) :: W(-2:,-2:,-2:)
    real(knd),intent(inout) :: Pr(1:,1:,1:)
    real(knd),intent(inout) :: Phi(0:,0:,0:)
    real(knd),intent(in)    :: dt2,dt3
    real(knd) :: Phiref,Au,Av,Aw,dxmin2,dymin2,dzmin2
    integer   :: i,j,k


    Au = dt2/dxmin
    Av = dt2/dymin
    Aw = dt2/dzmin

    dxmin2 = dxmin**2
    dymin2 = dymin**2
    dzmin2 = dzmin**2


    !$omp parallel
    !$omp do private(k,j,i)
    do k=1,Unz
     do j=1,Uny
      do i=1,Unx
        U(i,j,k)=U(i,j,k)-Au*(Phi(i+1,j,k)-Phi(i,j,k))
      end do
     end do
    end do
    !$omp end do nowait
    !$omp do private(k,j,i)
    do k=1,Vnz
     do j=1,Vny
      do i=1,Vnx
        V(i,j,k)=V(i,j,k)-Av*(Phi(i,j+1,k)-Phi(i,j,k))
      end do
     end do
    end do
    !$omp end do nowait
    !$omp do private(k,j,i)
    do k=1,Wnz
     do j=1,Wny
      do i=1,Wnx
        W(i,j,k)=W(i,j,k)-Aw*(Phi(i,j,k+1)-Phi(i,j,k))
      end do
     end do
    end do
    !$omp end do nowait

    !$omp do private(k,j,i)
    do k=1,Prnz
     do j=1,Prny
      do i=1,Prnx
        Pr(i,j,k)=Pr(i,j,k)+Phi(i,j,k)-dt3*(((Phi(i+1,j,k)-Phi(i,j,k))-(Phi(i,j,k)-Phi(i-1,j,k)))/dxmin2+&
                                            ((Phi(i,j+1,k)-Phi(i,j,k))-(Phi(i,j,k)-Phi(i,j-1,k)))/dymin2+&
                                            ((Phi(i,j,k+1)-Phi(i,j,k))-(Phi(i,j,k)-Phi(i,j,k-1)))/dzmin2)
      end do
     end do
    end do
    !$omp end do

    !$omp workshare
    Phiref = sum(Pr(1:Prnx,1:Prny,Prnz))/(Prnx*Prny)
    Pr = Pr - (Phiref - top_pressure)
    !$omp end workshare
    

    !$omp sections
    !$omp section
    call BoundU(1,U,Uin)
    !$omp section
    call BoundU(2,V,Vin)
    !$omp section
    call BoundU(3,W,Win)
    !$omp section
    call Bound_Pr(Pr)
    !$omp end sections
    !$omp end parallel



   end subroutine PostPoisson




























  subroutine POISS_POISFFT(Phi,RHS)
    use poisfft
#ifdef DPREC
    type(PoisFFT_Solver3D_DP),save :: Solver
#else
    type(PoisFFT_Solver3D_SP),save :: Solver
#endif
    real(knd),dimension(0:,0:,0:),intent(inout) :: Phi
    real(knd),dimension(0:,0:,0:),intent(in) :: RHS
    integer i
    logical, save :: called = .false.
    integer bounds(6)
    integer(DBL), save :: trate
    integer(DBL)       :: t1, t2

    do i=1,6
      if (Btype(i)==PERIODIC) then
        bounds(i) = PoisFFT_PERIODIC
      else
        bounds(i) = PoisFFT_NeumannStag
      end if
    end do

    if (.not.called) then
      call New(Solver,Prnx,Prny,Prnz,dxmin,dymin,dzmin,bounds)
      called = .true.
      call system_clock(count_rate=trate)
    end if

    call system_clock(count=t1)

    call Execute(Solver,Phi,RHS)

    call system_clock(count=t2)
    write(*,*) "solver cpu time", real(t2-t1)/real(trate)

  end subroutine POISS_POISFFT





  subroutine POISSSOR(Phi,RHS)        !Solves Poisson equation using Successive over-relaxation

  real(knd),dimension(0:,0:,0:),intent(inout) :: Phi
  real(knd),dimension(0:,0:,0:),intent(in) :: RHS
  integer,save :: called=0
  integer nx,ny,nz,i,j,k,l,blockit
  real(knd) S,P,Ap
  real(knd),dimension(:),allocatable,save :: Aw,Ae,As,An,Ab,At


   write (*,*) "Computing Poisson equation"
   S=0
   nx=Prnx
   ny=Prny
   nz=Prnz
   if (called==0) then                             !coefficients based on grid spacing computed only once
       allocate(Aw(1:nx),Ae(1:nx))
       allocate(As(1:ny),An(1:ny))
       allocate(Ab(1:nz),At(1:nz))
       forall(i=1:nx)
        Ae(i)=1._knd/(dxU(i)*dxPr(i))
        Aw(i)=1._knd/(dxU(i-1)*dxPr(i))
       end forall
       forall(j=1:ny)
        An(j)=1._knd/(dyV(j)*dyPr(j))
        As(j)=1._knd/(dyV(j-1)*dyPr(j))
       end forall
       forall(k=1:nz)
        At(k)=1._knd/(dzW(k)*dzPr(k))
        Ab(k)=1._knd/(dzW(k-1)*dzPr(k))
       end forall
       called=1
   end if

    l=0
    S=huge(1.0_knd)
        do while (l<=maxPoissoniter.and.S>epsPoisson)
          l=l+1
          S=0
          !$OMP PARALLEL PRIVATE(i,j,k,p) REDUCTION(max:S)
          !$OMP DO
          do k=1,nz
             do j=1,ny
                do i=1+mod(j+k,2),nx,2
                  p=0
                  Ap=0
                  if (i>1) then
                            p=p+Phi(i-1,j,k)*Aw(i)
                            Ap=Ap+Aw(i)
                  else if (Btype(We)==PERIODIC) then
                            p=p+Phi(nx,j,k)*Aw(i)
                            Ap=Ap+Aw(i)
                  end if
                  if (i<nx) then
                            p=p+Phi(i+1,j,k)*Ae(i)
                            Ap=Ap+Ae(i)
                  else if (Btype(We)==PERIODIC) then
                            p=p+Phi(1,j,k)*Ae(i)
                            Ap=Ap+Ae(i)
                  end if
                  if (j>1) then
                            p=p+Phi(i,j-1,k)*As(j)
                            Ap=Ap+As(j)
                  else if (Btype(No)==PERIODIC) then
                            p=p+Phi(i,ny,k)*As(j)
                            Ap=Ap+As(j)
                  end if
                  if (j<ny) then
                            p=p+Phi(i,j+1,k)*An(j)
                            Ap=Ap+An(j)
                  else if (Btype(No)==PERIODIC) then
                            p=p+Phi(i,1,k)*An(j)
                            Ap=Ap+An(j)
                  end if
                  if (k>1) then
                            p=p+Phi(i,j,k-1)*Ab(k)
                            Ap=Ap+Ab(k)
                  else if (Btype(To)==PERIODIC) then
                            p=p+Phi(i,j,nz)*Ab(k)
                            Ap=Ap+Ab(k)
                  end if
                  if (k<nz) then
                            p=p+Phi(i,j,k+1)*At(k)
                            Ap=Ap+At(k)
                  else if (Btype(To)==PERIODIC) then
                            p=p+Phi(i,j,1)*At(k)
                            Ap=Ap+At(k)
                  end if
                  p=p-RHS(i,j,k)
                  p=p/Ap
                  S=max(abs(p-Phi(i,j,k)),S)
                  Phi(i,j,k)=p
                end do
             end do
          end do
          !$OMP ENDDO
          !$OMP DO
          do k=1,nz
             do j=1,ny
                do i=1+mod(j+k+1,2),nx,2
                  p=0
                  Ap=0
                  if (i>1) then
                            p=p+Phi(i-1,j,k)*Aw(i)
                            Ap=Ap+Aw(i)
                  else if (Btype(We)==PERIODIC) then
                            p=p+Phi(nx,j,k)*Aw(i)
                            Ap=Ap+Aw(i)
                  end if
                  if (i<nx) then
                            p=p+Phi(i+1,j,k)*Ae(i)
                            Ap=Ap+Ae(i)
                  else if (Btype(We)==PERIODIC) then
                            p=p+Phi(1,j,k)*Ae(i)
                            Ap=Ap+Ae(i)
                  end if
                  if (j>1) then
                            p=p+Phi(i,j-1,k)*As(j)
                            Ap=Ap+As(j)
                  else if (Btype(No)==PERIODIC) then
                            p=p+Phi(i,ny,k)*As(j)
                            Ap=Ap+As(j)
                  end if
                  if (j<ny) then
                            p=p+Phi(i,j+1,k)*An(j)
                            Ap=Ap+An(j)
                  else if (Btype(No)==PERIODIC) then
                            p=p+Phi(i,1,k)*An(j)
                            Ap=Ap+An(j)
                  end if
                  if (k>1) then
                            p=p+Phi(i,j,k-1)*Ab(k)
                            Ap=Ap+Ab(k)
                  else if (Btype(To)==PERIODIC) then
                            p=p+Phi(i,j,nz)*Ab(k)
                            Ap=Ap+Ab(k)
                  end if
                  if (k<nz) then
                            p=p+Phi(i,j,k+1)*At(k)
                            Ap=Ap+At(k)
                  else if (Btype(To)==PERIODIC) then
                            p=p+Phi(i,j,1)*At(k)
                            Ap=Ap+At(k)
                  end if
                  p=p-RHS(i,j,k)
                  p=p/Ap
                  S=max(abs(p-Phi(i,j,k)),S)
                  Phi(i,j,k)=p
                end do
             end do
          end do
          !$OMP ENDDO
        !$OMP ENDPARALLEL
          p=abs(maxval(Phi(1:nx,1:ny,1:nz)))
          if (p>0) S=S/p
          if (MOD(l,10)==0)  write (*,*) "   Poisson iter: ",l,S
        end do

 end subroutine POISSSOR



end module POISSON

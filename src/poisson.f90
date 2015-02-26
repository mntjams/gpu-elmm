module PoissonSolvers

  use Parameters
  use Multigrid,   only: PoissMG
  use Multigrid2D, only: PoissMG2d


contains

  subroutine Poiss_PoisFFT(Phi,RHS)
#ifdef DPREC
    use PoisFFT, PoisFFT_Solver => PoisFFT_Solver3D_DP
#else
    use PoisFFT, PoisFFT_Solver => PoisFFT_Solver3D_SP
#endif
#ifdef MPI
    use custom_mpi
#endif

    type(PoisFFT_Solver),save :: Solver
    real(knd),dimension(0:,0:,0:),intent(inout) :: Phi
    real(knd),dimension(0:,0:,0:),intent(in) :: RHS
    logical, save :: called = .false.
#ifdef POIS_SOLVER_TIME
    integer(DBL), save :: trate
    integer(DBL)       :: t1, t2
#endif


    if (.not.called) then
#ifdef MPI
      Solver =  PoisFFT_Solver([Prnx,Prny,Prnz], &
                               [gxmax-gxmin,gymax-gymin,gzmax-gzmin], &
                               PoissonBtype, &
                               PoisFFT_FiniteDifference2, &
                               gPrns,offsets_to_global,poisfft_comm)
#else
      Solver =  PoisFFT_Solver([Prnx,Prny,Prnz], &
                               [gxmax-gxmin,gymax-gymin,gzmax-gzmin], &
                               PoissonBtype, &
                               PoisFFT_FiniteDifference2)
#endif
      called = .true.
#ifdef POIS_SOLVER_TIME
      call system_clock(count_rate=trate)
#endif
    end if

#ifdef POIS_SOLVER_TIME
    call system_clock(count=t1)
#endif

    call Execute(Solver,Phi,RHS)

#ifdef POIS_SOLVER_TIME
    call system_clock(count=t2)
    if (master) then
      poisson_solver_time = poisson_solver_time + real(t2-t1,knd)/real(trate,knd)
      if (master) write(*,*) "solver cpu time", real(t2-t1)/real(trate)
    end if
#endif

  end subroutine





  subroutine PoissSOR(Phi,RHS) 
#ifdef MPI
    use Boundaries
    use custom_mpi
#endif
    !Solves Poisson equation using Successive over-relaxation

    real(knd),dimension(0:,0:,0:),intent(inout) :: Phi
    real(knd),dimension(0:,0:,0:),intent(in) :: RHS
    integer,save :: called=0
    integer :: nx,ny,nz,i,j,k,l
    real(knd) :: S,P,Ap
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
#ifdef MPI
      call Bound_Phi(Phi)
#endif
      !$OMP PARALLEL PRIVATE(i,j,k,p) REDUCTION(max:S)
      !$OMP DO
      do k=1,nz
        do j=1,ny
          do i=1+mod(j+k,2),nx,2
            p=0
            Ap=0
            if (i>1.or.Btype(We)>=MPI_BOUNDS) then
                      p=p+Phi(i-1,j,k)*Aw(i)
                      Ap=Ap+Aw(i)
            else if (Btype(We)==PERIODIC) then
                      p=p+Phi(nx,j,k)*Aw(i)
                      Ap=Ap+Aw(i)
            end if
            if (i<nx.or.Btype(Ea)>=MPI_BOUNDS) then
                      p=p+Phi(i+1,j,k)*Ae(i)
                      Ap=Ap+Ae(i)
            else if (Btype(We)==PERIODIC) then
                      p=p+Phi(1,j,k)*Ae(i)
                      Ap=Ap+Ae(i)
            end if
            if (j>1.or.Btype(So)>=MPI_BOUNDS) then
                      p=p+Phi(i,j-1,k)*As(j)
                      Ap=Ap+As(j)
            else if (Btype(No)==PERIODIC) then
                      p=p+Phi(i,ny,k)*As(j)
                      Ap=Ap+As(j)
            end if
            if (j<ny.or.Btype(No)>=MPI_BOUNDS) then
                      p=p+Phi(i,j+1,k)*An(j)
                      Ap=Ap+An(j)
            else if (Btype(No)==PERIODIC) then
                      p=p+Phi(i,1,k)*An(j)
                      Ap=Ap+An(j)
            end if
            if (k>1.or.Btype(Bo)>=MPI_BOUNDS) then
                      p=p+Phi(i,j,k-1)*Ab(k)
                      Ap=Ap+Ab(k)
            else if (Btype(To)==PERIODIC) then
                      p=p+Phi(i,j,nz)*Ab(k)
                      Ap=Ap+Ab(k)
            end if
            if (k<nz.or.Btype(To)>=MPI_BOUNDS) then
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
      
#ifdef MPI
      S = mpi_co_max(S)
#endif

      if (MOD(l,10)==0)  write (*,*) "   Poisson iter: ",l,S
    end do

  end subroutine PoissSOR


end module PoissonSolvers
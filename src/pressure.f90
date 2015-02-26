module Pressure

  use Parameters
  use Boundaries
  use PoissonSolvers

  implicit none

  private
  public PressureCorrection


contains


  subroutine PressureCorrection(U,V,W,Pr,Q,coef)                    !Pressure correction
    real(knd), dimension(-2:,-2:,-2:), intent(inout)     :: U, V, W !Phi is computed in Poisson eq. with div of U in RHS
    real(knd), dimension(1:,1:,1:), intent(inout)        :: Pr      !Depend ing on active projection method Phi becomes new pressure
    real(knd), dimension(:,:,:), allocatable, intent(in) :: Q       !or is added to last pressure
    real(knd), intent(in) :: coef
                                                           !U,V,W velocity field for correction
    real(knd), save, allocatable :: Phi(:,:,:), RHS(:,:,:) !Pr pressure
                                                           !coef cofficient from Runge Kutta, Q mass sources from immersed boundary
    real(TIM) dt2,dt3                                      !RHS right hand side of eq. with divergence of U
                                                           !Phi computed pseudopressure, saved as first guess for next time
    real(knd) :: uncompatibility
    integer, save :: called = 0
    integer(DBL), save :: trate
    integer(DBL), save :: time1, time2, time3, time4
#ifdef MPI
    correctcompatibility = 0
#else
    correctcompatibility = 0
#endif
    if (called==0) then
      allocate(Phi(0:Prnx+1,0:Prny+1,0:Prnz+1))
      allocate(RHS(0:Prnx+1,0:Prny+1,0:Prnz+1))
      Phi = 0

      if (debugparam>1) call system_clock(count_rate=trate)
    end if
    called = called + 1


    if (debugparam>1.and.called>1) call system_clock(count=time1)

    dt2 = coef
    if (Re>0)  then
      dt3 = coef / (2._TIM*Re)
    else
      dt3 = 0
    end if


    call PrePoisson(U,V,W,Q,RHS,dt2,uncompatibility)


    if (debugparam>1.and.called>1) call system_clock(count=time2)

    if (correctcompatibility>=1) then
      if (master) write(*,*) "Uncompatibility:",uncompatibility
    end if

    if (poissmet==1) then

        call PoissSOR(Phi,RHS)

    else if (poissmet==2) then

        call Poiss_PoisFFT(Phi,RHS)

    else if (poissmet==3) then

        if (Prny==1) then
          call PoissMG2d(Phi,RHS)
        else
          call PoissMG(Phi,RHS)
        end if

    end if


    call Bound_Phi(Phi)


    if (debugparam>1.and.called>1) then
     call system_clock(count=time3)
     if (master) write(*,*) "ET of part 2", real(time3-time2)/real(trate)
    endif


    call PostPoisson(U,V,W,Pr,Q,Phi,dt2,dt3)


    if (debugparam>1.and.called>1) then
     call system_clock(count=time4)
     if (master) write(*,*) "ET of part 3", real((time4-time1)-(time3-time2))/real(trate)
    endif
    
  end subroutine PressureCorrection





  subroutine PrePoisson(U,V,W,Q,RHS,dt2,uncompatibility)
#ifdef MPI
    use custom_mpi
#endif
    real(knd), intent(inout) :: U(-2:,-2:,-2:)
    real(knd), intent(inout) :: V(-2:,-2:,-2:)
    real(knd), intent(inout) :: W(-2:,-2:,-2:)
    real(knd), intent(out)   :: RHS(0:,0:,0:)
    real(knd), allocatable, intent(in) :: Q(:,:,:)
    real(knd), intent(out)   :: uncompatibility
    real(knd), intent(in)    :: dt2
    real(knd) :: S
    integer   :: i,j,k

!     !$omp parallel
!     !$omp sections
!     !$omp section
    call BoundU(1,U,Uin)
!     !$omp section
    call BoundU(2,V,Vin)
!     !$omp section
    call BoundU(3,W,Win)
!     !$omp end sections
!     !$omp end parallel

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

#ifdef MPI
      if (abs(windangle-90)<1.or.abs(windangle+90)<1) then
        S=S*dymin/(gPrnx*gPrnz)

        uncompatibility = S

        uncompatibility = mpi_co_sum(uncompatibility)

        if (correctcompatibility==1.and.jim==nyims) then
          !$omp parallel workshare
          V(:,Vny+1,:) = V(:,Vny+1,:)+S
          !$omp end parallel workshare
        end if
      else
        S=S*dxmin/(gPrny*gPrnz)

        uncompatibility = S

        uncompatibility = mpi_co_sum(uncompatibility)

        if (correctcompatibility==1.and.iim==nxims) then
          !$omp parallel workshare
          U(Unx+1,:,:) = U(Unx+1,:,:)+S
          !$omp end parallel workshare
        end if
      end if
#else
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
#endif
    end if

    S=0

    if (allocated(Q)) then
      !$omp parallel do private(i,j,k)
      do k=1,Prnz            !divergence of U -> RHS
       do j=1,Prny
        do i=1,Prnx
             RHS(i,j,k) = (U(i,j,k)-U(i-1,j,k))/(dxmin)&
                         +(V(i,j,k)-V(i,j-1,k))/(dymin)&
                         +(W(i,j,k)-W(i,j,k-1))/(dzmin)&
                         -Q(i,j,k)
             RHS(i,j,k) = RHS(i,j,k)/(dt2)
        end do
       end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do private(i,j,k)
      do k=1,Prnz            !divergence of U -> RHS
       do j=1,Prny
        do i=1,Prnx
             RHS(i,j,k) = (U(i,j,k)-U(i-1,j,k))/(dxmin)&
                         +(V(i,j,k)-V(i,j-1,k))/(dymin)&
                         +(W(i,j,k)-W(i,j,k-1))/(dzmin)
             RHS(i,j,k) = RHS(i,j,k)/(dt2)
        end do
       end do
      end do
      !$omp end parallel do
    end if

  end subroutine PrePoisson




  subroutine PostPoisson(U,V,W,Pr,Q,Phi,dt2,dt3)
#ifdef MPI
    use custom_mpi, only: kim, nzims, &
                          mpi_co_max, mpi_co_sum, &
                          comm_plane_xy, comm_row_z, &
                          MPI_KND

    interface
      subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
        import
        real(knd) ::  BUFFER
        INTEGER   COUNT, DATATYPE, ROOT, COMM, IERROR
      end subroutine
    end interface

    integer :: ie
#endif
    real(knd), intent(inout) :: U(-2:,-2:,-2:)
    real(knd), intent(inout) :: V(-2:,-2:,-2:)
    real(knd), intent(inout) :: W(-2:,-2:,-2:)
    real(knd), intent(inout) :: Pr(1:,1:,1:)
    real(knd), allocatable, intent(in) :: Q(:,:,:)
    real(knd), intent(inout) :: Phi(0:,0:,0:)
    real(knd), intent(in)    :: dt2,dt3
    real(knd) :: Phi_ref,Au,Av,Aw,dxmin2,dymin2,dzmin2,S,p
    integer   :: i,j,k
#ifdef CHECK_DIVERGENCE
    logical, parameter :: check_divergence = .true.
#else
    logical, parameter :: check_divergence = .false.
#endif


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

    if (explicit_diffusion) then
      !$omp do private(k,j,i)
      do k=1,Prnz
       do j=1,Prny
        do i=1,Prnx
          Pr(i,j,k) = Pr(i,j,k) + Phi(i,j,k)
        end do
       end do
      end do
      !$omp end do
    else
      !$omp do private(k,j,i)
      do k=1,Prnz
       do j=1,Prny
        do i=1,Prnx
          Pr(i,j,k) = Pr(i,j,k) + Phi(i,j,k) - &
                       dt3*(((Phi(i+1,j,k)-Phi(i,j,k)) - &
                             (Phi(i,j,k)-Phi(i-1,j,k)))/dxmin2 + &
                            ((Phi(i,j+1,k)-Phi(i,j,k)) - &
                             (Phi(i,j,k)-Phi(i,j-1,k)))/dymin2 + &
                            ((Phi(i,j,k+1)-Phi(i,j,k)) - &
                             (Phi(i,j,k)-Phi(i,j,k-1)))/dzmin2)
        end do
       end do
      end do
      !$omp end do
    end if

#ifdef MPI
    !images in top plane compute the reference pressure
    if (kim==nzims) then
      !$omp workshare
      Phi_ref = sum(Pr(1:Prnx,1:Prny,Prnz))
      !$omp end workshare
      Phi_ref = mpi_co_sum(Phi_ref, comm = comm_plane_xy)
      Phi_ref = Phi_ref / (gPrnx * gPrny)
    end if

    !all kim==nzims broadcast to images with smaller kim
    call MPI_Bcast(Phi_ref, 1, MPI_KND, nzims-1, comm_row_z, ie)
    if (ie/=0) call error_stop("Error calling MPI_Bcast, in "//__FILE__//" line",__LINE__)

#else
    !$omp workshare
    Phi_ref = sum(Pr(1:Prnx,1:Prny,Prnz)) / (Prnx*Prny)
    !$omp end workshare
#endif
    !$omp workshare
    Pr = Pr - (Phi_ref - top_pressure)
    !$omp end workshare
    

    !$omp end parallel
!     !$omp sections
!     !$omp section
    call BoundU(1,U,Uin)
!     !$omp section
    call BoundU(2,V,Vin)
!     !$omp section
    call BoundU(3,W,Win)
!     !$omp section
    call Bound_Pr(Pr)
!     !$omp end sections


    if (check_divergence) then
      S = 0
      if (allocated(Q)) then
        !$omp parallel do private(i,j,k,p) reduction(max:S)
        do k=1,Prnz            !divergence of U -> RHS
         do j=1,Prny
          do i=1,Prnx
               p = (U(i,j,k)-U(i-1,j,k))/(dxmin)&
                           +(V(i,j,k)-V(i,j-1,k))/(dymin)&
                           +(W(i,j,k)-W(i,j,k-1))/(dzmin)&
                           -Q(i,j,k)
               S = max(S,abs(p))
          end do
         end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do private(i,j,k,p) reduction(max:S)
        do k=1,Prnz            !divergence of U -> RHS
         do j=1,Prny
          do i=1,Prnx
               p = (U(i,j,k)-U(i-1,j,k))/(dxmin)&
                           +(V(i,j,k)-V(i,j-1,k))/(dymin)&
                           +(W(i,j,k)-W(i,j,k-1))/(dzmin)
               S = max(S,abs(p))
          end do
         end do
        end do
        !$omp end parallel do
      end if
#ifdef MPI
      S = mpi_co_max(S)
#endif
      
      if (master) write(*,*) "max divergence:", S
    end if

  end subroutine PostPoisson


end module Pressure

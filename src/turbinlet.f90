module TURBINLET

  use PARAMETERS

  implicit none

  private
  public GetTurbInlet, GetInletFromFile, TLag, Lturby, Lturbz, ustarinlet, transformtensor

  real(knd),save :: TLag
  real(knd),save :: Lturby
  real(knd),save :: Lturbz


  real(knd),allocatable,dimension(:),save     :: ustarinlet !friction velocity profile at inlet
  real(knd),allocatable,dimension(:,:),save   :: Uinavg,Vinavg,Winavg !mean values of U,V,W at inflow
  real(knd),allocatable,dimension(:,:,:),save :: transformtensor

   type TInlet
    real(knd),allocatable,dimension(:,:) :: U,V,W,temperature
   endtype

  interface GETTURBINLET
   module procedure GETTURBINLETXIE
  endinterface

contains


 subroutine GETTURBINLETKLEIN
 logical,save:: called=.false.
 integer i,j,k,ii,jj,kk
 integer,save:: filtnx,filtny,filtnz,bigNx,bigNy,bigNz
 real(knd) Ui,Vi,Wi,bx,by,bz,bxsum,bysum,bzsum
 real(knd),allocatable,dimension(:,:,:),save :: Ru,Rv,Rw !arrays of randoms
 real(knd),allocatable,dimension(:,:,:),save :: bfilt !filter coefficients (ii,jj,kk,kz)

 if (.not. called) then
    write(*,*) "INITTURBULENCEPROFILES"
    call INITTURBULENCEPROFILES
    write(*,*) "INITMEANPROFILES"
    call INITMEANPROFILES
 write(*,*) "GETTURBINLET-cont"

    filtnz=1!max(Prnz/20,2)
    filtny=1!max((Prnz*dzmin/dymin)/20,2)
    filtnx=1!max((Prnz*dzmin/dxmin)/20,2)

    write(*,*) "filtn",filtnx,filtny,filtnz
    bigNx=2*filtnx
    bigNy=2*filtny
    bigNz=2*filtnz
    write(*,*) "bigN",bigNx,bigNy,bigNz

    allocate(Ru(-bigNx:bigNx,-bigNy+1:Prny+bigNy,-bigNz+1:Prnz+bigNz))
    allocate(Rv(-bigNx:bigNx,-bigNy+1:Prny+bigNy,-bigNz+1:Prnz+bigNz))
    allocate(Rw(-bigNx:bigNx,-bigNy+1:Prny+bigNy,-bigNz+1:Prnz+bigNz))

    do k=-bigNz+1,Prnz+bigNz
     do j=-bigNy+1,Prny+bigNy
      do i=-bigNx,bigNx
       call RandomGauss(Ru(i,j,k))
       call RandomGauss(Rv(i,j,k))
       call RandomGauss(Rw(i,j,k))
      enddo
     enddo
    enddo

    if ((Btype(So)==PERIODIC).or.(Btype(No)==PERIODIC)) then
      do k=-bigNz+1,Prnz+bigNz
       do j=-bigNy+1,0
        do i=-bigNx,bigNx
         Ru(i,j,k)=Ru(i,j+Prny,k)
         Rv(i,j,k)=Rv(i,j+Prny,k)
         Rw(i,j,k)=Rw(i,j+Prny,k)
        enddo
       enddo
      enddo
      do k=-bigNz+1,Prnz+bigNz
       do j=Prny+1,Prny+bigNy
        do i=-bigNx,bigNx
         Ru(i,j,k)=Ru(i,j-Prny,k)
         Rv(i,j,k)=Rv(i,j-Prny,k)
         Rw(i,j,k)=Rw(i,j-Prny,k)
        enddo
       enddo
      enddo
    endif
    if  ((Btype(Bo)==PERIODIC).or.(Btype(To)==PERIODIC)) then
      do k=-bigNz+1,0
       do j=-bigNy+1,Prny+bigNy
        do i=-bigNx,bigNx
         Ru(i,j,k)=Ru(i,j,k+Prnz)
         Rv(i,j,k)=Rv(i,j,k+Prnz)
         Rw(i,j,k)=Rw(i,j,k+Prnz)
        enddo
       enddo
      enddo
      do k=Prnz+1,Prnz+bigNz
       do j=-bigNy+1,Prny+bigNy
        do i=-bigNx,bigNx
         Ru(i,j,k)=Ru(i,j,k-Prnz)
         Rv(i,j,k)=Rv(i,j,k-Prnz)
         Rw(i,j,k)=Rw(i,j,k-Prnz)
        enddo
       enddo
      enddo
    endif


    allocate(bfilt(-bigNx:bigNx,-bigNy:bigNy,-bigNz:bigNz))

    bxsum=0
    do i=-bigNx,bigNx
     bxsum=bxsum+exp(-pi*i*i/(bigNx*bigNx))
    enddo
    bxsum=SQRT(bxsum)

    bysum=0
    do i=-bigNy,bigNy
     bysum=bysum+exp(-pi*i*i/(bigNy*bigNy))
    enddo
    bysum=SQRT(bysum)

    bzsum=0
    do i=-bigNz,bigNz
     bzsum=bzsum+exp(-pi*i*i/(bigNz*bigNz))
    enddo
    bzsum=SQRT(bzsum)

    do k=-bigNz,bigNz
      bz=exp(-pi*k*k/(2*bigNz*bigNz))/bzsum

      do j=-bigNy,bigNy
        by=exp(-pi*j*j/(2*bigNy*bigNy))/bysum

        do i=-bigNx,bigNx
         bx=exp(-pi*i*i/(2*bigNx*bigNx))/bxsum
         bfilt(i,j,k)=bx*by*bz
        enddo
      enddo
    enddo

    called=.true.
 endif
    Uin=0
    Vin=0
    Win=0
    do k=1,Prnz
     do j=1,Prny
      do kk=-bigNz,bigNz
       do jj=-bigNy,bigNy
        do ii=-bigNx,bigNx
         Uin(j,k)=Uin(j,k)+bfilt(ii,jj,kk)*Ru(ii,j+jj,k+kk)
         Vin(j,k)=Vin(j,k)+bfilt(ii,jj,kk)*RV(ii,j+jj,k+kk)
         Win(j,k)=Win(j,k)+bfilt(ii,jj,kk)*RW(ii,j+jj,k+kk)
        enddo
       enddo
      enddo
     enddo
    enddo

    do k=1,Prnz
     do j=1,Prny
      Ui=Uin(j,k)
      Vi=Vin(j,k)
      Wi=Win(j,k)

      Uin(j,k)=Uinavg(j,k)+transformtensor(1,j,k)*Ui   !a12,a13,a23=0

      Vin(j,k)=Vinavg(j,k)+transformtensor(2,j,k)*Ui&
                            +transformtensor(3,j,k)*Vi

      Win(j,k)=Winavg(j,k)+transformtensor(4,j,k)*Ui&
                            +transformtensor(5,j,k)*Vi&
                            +transformtensor(6,j,k)*Wi
     enddo
    enddo

    do k=-bigNz+1,Prnz+bigNz
     do j=-bigNy+1,Prny+bigNy
      do i=-bigNx,bigNx-1
       Ru(i,j,k)=Ru(i+1,j,k)
       Rv(i,j,k)=Rv(i+1,j,k)
       Rw(i,j,k)=Rw(i+1,j,k)
      enddo
     enddo
    enddo


    do k=-bigNz+1,Prnz+bigNz
     do j=-bigNy+1,Prny+bigNy
       call RandomGauss(Ru(bigNx,j,k))
       call RandomGauss(Rv(bigNx,j,k))
       call RandomGauss(Rw(bigNx,j,k))
      enddo
     enddo


    if ((Btype(So)==PERIODIC).or.(Btype(No)==PERIODIC)) then
      do k=-bigNz+1,Prnz+bigNz
       do j=-bigNy+1,0
         Ru(bigNx,j,k)=Ru(bigNx,j+Prny,k)
         Rv(bigNx,j,k)=Rv(bigNx,j+Prny,k)
         Rw(bigNx,j,k)=Rw(bigNx,j+Prny,k)
       enddo
      enddo
      do k=-bigNz+1,Prnz+bigNz
       do j=Prny+1,Prny+bigNy
         Ru(bigNx,j,k)=Ru(bigNx,j-Prny,k)
         Rv(bigNx,j,k)=Rv(bigNx,j-Prny,k)
         Rw(bigNx,j,k)=Rw(bigNx,j-Prny,k)
       enddo
      enddo
    endif
    if  ((Btype(Bo)==PERIODIC).or.(Btype(To)==PERIODIC)) then
      do k=-bigNz+1,0
       do j=-bigNy+1,Prny+bigNy
         Ru(bigNx,j,k)=Ru(bigNx,j,k+Prnz)
         Rv(bigNx,j,k)=Rv(bigNx,j,k+Prnz)
         Rw(bigNx,j,k)=Rw(bigNx,j,k+Prnz)
       enddo
      enddo
      do k=Prnz+1,Prnz+bigNz
       do j=-bigNy+1,Prny+bigNy
         Ru(bigNx,j,k)=Ru(bigNx,j,k-Prnz)
         Rv(bigNx,j,k)=Rv(bigNx,j,k-Prnz)
         Rw(bigNx,j,k)=Rw(bigNx,j,k-Prnz)
       enddo
      enddo
    endif

!     write(*,*) "Inlets:"
!     do k=1,Unz
!      write(*,*) "U",k,Uin(1,k),Uinavg(1,k)
!     enddo
!     do k=1,Vnz
!      write(*,*) "V",k,Vin(1,k),Vinavg(1,k)
!     enddo
!     do k=1,Wnz
!      write(*,*) "W",k,Win(1,k),Winavg(1,k)
!     enddo
!
!   OPEN(11,file="RandU.vtk")
!   write (11,"(A)") "# vtk DataFile Version 2.0"
!   write (11,"(A)") "diplomka output file"
!   write (11,"(A)") "ASCII"
!   write (11,"(A)") "DATASET RECTILINEAR_GRID"
!   str="DIMENSIONS"
!   write (str(12:),*) (2*bigNx+1),(Prny+2*bigNy),(Prnz+2*bigNz)
!   write (11,"(A)") str
!   str="X_COORDINATES"
!   write (str(15:),*) (2*bigNx+1),"float"
!   write (11,"(A)") str
!   write (11,*) (/(i,i=-bigNx,bigNx)/)
!   str="Y_COORDINATES"
!   write (str(15:),*) (Prny+2*bigNy),"float"
!   write (11,"(A)") str
!   write (11,*) (/(i,i=-bigNy+1,Prny+bigNy)/)
!   str="Z_COORDINATES"
!   write (str(15:),*) (Prnz+2*bigNz),"float"
!   write (11,"(A)") str
!   write (11,*) (/(i,i=-bigNz+1,Prnz+bigNz)/)
!   str="POINT_DATA"
!   write (str(12:),*) (2*bigNx+1)*(Prny+2*bigNy)*(Prnz+2*bigNz)
!   write (11,"(A)") str
!
!
!   write (11,"(A)") "SCALARS RandU float"
!   write (11,"(A)") "LOOKUP_TABLE default"
!     do k=-bigNz+1,Prnz+bigNz
!      do j=-bigNy+1,Prny+bigNy
!       do i=-bigNx,bigNx
!       Write (11,*) Ru(i,j,k)
!     enddo
!    enddo
!   enddo
!   write (11,*)
!   CLOSE(11)

 endsubroutine GETTURBINLETKLEIN













 subroutine GETTURBINLETXIE
 logical,save:: called=.false.
 integer i,j,k
 integer,save:: filtny,filtnz,bigNy,bigNz
 real(knd) Ui,Vi,Wi,bysum,bzsum,p
 real(knd),allocatable,dimension(:):: expsy,expsz
 real(knd),save:: compat
 real(knd),allocatable,dimension(:,:,:),save :: Ru,Rv,Rw !arrays of randoms
 real(knd),allocatable,dimension(:,:,:),save :: Psiu,Psiv,Psiw
 real(knd),allocatable,dimension(:,:,:),save :: bfilt !filter coefficients (ii,jj,kk,kz)

 if (.not. called) then
    write(*,*) "INITTURBULENCEPROFILES"
    call INITTURBULENCEPROFILES
    write(*,*) "INITMEANPROFILES"
    call INITMEANPROFILES
 write(*,*) "GETTURBINLET-cont"


    filtnz=min(max(NINT(lturbz/dzmin),1),ceiling(1._knd*Prnz/3))
    filtny=min(max(NINT(lturby/dymin),1),ceiling(1._knd*Prny/3))


    bigNy=2*filtny
    bigNz=2*filtnz


    allocate(Ru(-bigNy+1:Prny+bigNy,-bigNz+1:Prnz+bigNz,2))
    allocate(Rv(-bigNy+1:Prny+bigNy,-bigNz+1:Prnz+bigNz,2))
    allocate(Rw(-bigNy+1:Prny+bigNy,-bigNz+1:Prnz+bigNz,2))
    allocate(Psiu(1:Prny,1:Prnz,1:2))
    allocate(Psiv(1:Prny,1:Prnz,1:2))
    allocate(Psiw(1:Prny,1:Prnz,1:2))


    do k=-bigNz+1,Prnz+bigNz
     do j=-bigNy+1,Prny+bigNy
       call RandomGauss(Ru(j,k,1))
       call RandomGauss(Rv(j,k,1))
       call RandomGauss(Rw(j,k,1))
     enddo
    enddo

    do k=-bigNz+1,Prnz+bigNz
     do j=-bigNy+1,Prny+bigNy
       call RandomGauss(Ru(j,k,2))
       call RandomGauss(Rv(j,k,2))
       call RandomGauss(Rw(j,k,2))
     enddo
    enddo



    if ((Btype(So)==PERIODIC).or.(Btype(No)==PERIODIC)) then
      !$omp parallel workshare
      forall(k=-bigNz+1:Prnz+bigNz,j=-bigNy+1:0)
         Ru(j,k,1:2)=Ru(j+Prny,k,1:2)
         Rv(j,k,1:2)=Rv(j+Prny,k,1:2)
         Rw(j,k,1:2)=Rw(j+Prny,k,1:2)
      endforall
      forall(k=-bigNz+1:Prnz+bigNz,j=Prny+1:Prny+bigNy)
         Ru(j,k,1:2)=Ru(j-Prny,k,1:2)
         Rv(j,k,1:2)=Rv(j-Prny,k,1:2)
         Rw(j,k,1:2)=Rw(j-Prny,k,1:2)
      endforall
      !$omp end parallel workshare
    endif
    if  ((Btype(Bo)==PERIODIC).or.(Btype(To)==PERIODIC)) then
      !$omp parallel workshare
      forall(k=-bigNz+1:0,j=-bigNy+1:Prny+bigNy)
         Ru(j,k,1:2)=Ru(j,k+Prnz,1:2)
         Rv(j,k,1:2)=Rv(j,k+Prnz,1:2)
         Rw(j,k,1:2)=Rw(j,k+Prnz,1:2)
      endforall
      forall(k=Prnz+1:Prnz+bigNz,j=-bigNy+1:Prny+bigNy)
         Ru(j,k,1:2)=Ru(j,k-Prnz,1:2)
         Rv(j,k,1:2)=Rv(j,k-Prnz,1:2)
         Rw(j,k,1:2)=Rw(j,k-Prnz,1:2)
      endforall
      !$omp end parallel workshare
    endif


    allocate(bfilt(-bigNy:bigNy,-bigNz:bigNz,1))
    allocate(expsy(-bigNy:bigNy),expsz(-bigNz:bigNz))

    bysum=0
    do i=-bigNy,bigNy
     expsy(i)=exp(-pi*abs(i)/(bigNy))
     bysum=bysum+expsy(i)**2
    enddo
    bysum=SQRT(bysum)
    expsy=expsy/bysum


    bzsum=0
    do i=-bigNz,bigNz
     expsz(i)=exp(-pi*abs(i)/(bigNz))
     bzsum=bzsum+expsz(i)**2
    enddo
    bzsum=SQRT(bzsum)
    expsz=expsz/bzsum

    !$omp parallel private(j,k)
    !$omp do
    do k=-bigNz,bigNz
      do j=-bigNy,bigNy
         bfilt(j,k,1)=expsy(j)*expsz(k)
      enddo
    enddo
    !$omp end do
    !$omp workshare
    forall(k=1:Prnz,j=1:Prny)
         Psiu(j,k,1)=SUM(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Ru(j-bigNy:j+bigNy,k-bigNz:k+bigNz,1))
         Psiv(j,k,1)=SUM(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Rv(j-bigNy:j+bigNy,k-bigNz:k+bigNz,1))
         Psiw(j,k,1)=SUM(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Rw(j-bigNy:j+bigNy,k-bigNz:k+bigNz,1))
    endforall

    compat=SUM(Uinavg(1:Prny,1:Prnz))
    !$omp end workshare
    !$omp end parallel

    called=.true.

 endif

    !$omp parallel private(i,j,k,Ui,Vi,Wi)

    !$omp workshare
    forall(k=1:Prnz,j=1:Prny)
         Psiu(j,k,2)=SUM(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Ru(j-bigNy:j+bigNy,k-bigNz:k+bigNz,2))
         Psiv(j,k,2)=SUM(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Rv(j-bigNy:j+bigNy,k-bigNz:k+bigNz,2))
         Psiw(j,k,2)=SUM(bfilt(-bigNy:bigNy,-bigNz:bigNz,1)*Rw(j-bigNy:j+bigNy,k-bigNz:k+bigNz,2))
    endforall

    Psiu(:,:,1)=Psiu(:,:,1)*exp(-pi*dt/(2._knd*TLag))+Psiu(:,:,2)*sqrt(1-exp(-pi*dt/(TLag)))
    Psiv(:,:,1)=Psiv(:,:,1)*exp(-pi*dt/(2._knd*TLag))+Psiv(:,:,2)*sqrt(1-exp(-pi*dt/(TLag)))
    Psiw(:,:,1)=Psiw(:,:,1)*exp(-pi*dt/(2._knd*TLag))+Psiw(:,:,2)*sqrt(1-exp(-pi*dt/(TLag)))
    !$omp end workshare

    !$omp do
    do k=1,Prnz
     do j=1,Prny
      Ui=Psiu(j,k,1)
      Vi=Psiv(j,k,1)
      Wi=Psiw(j,k,1)

      Uin(j,k)=Uinavg(j,k)+transformtensor(1,j,k)*Ui   !a12,a13,a23=0

      Vin(j,k)=Vinavg(j,k)+transformtensor(2,j,k)*Ui&
                            +transformtensor(3,j,k)*Vi

      Win(j,k)=Winavg(j,k)+transformtensor(4,j,k)*Ui&
                            +transformtensor(5,j,k)*Vi&
                            +transformtensor(6,j,k)*Wi
     enddo
    enddo
    !$omp end do

    !$omp workshare
    p=compat/SUM(Uin(1:Prny,1:Prnz))                  !To ensure the compatibility condition.
    !$omp end workshare

    !gfortran-4.7 signaled SIGFPE here when array assignments were used - compiler bug?
    !$omp do
    do k=1,Prnz
     do j=1,Prny
       Uin(j,k)=Uin(j,k)*p
     end do
    end do
    !$omp end do
    !$omp do
    do k=1,Prnz
     do j=1,Prny
       Vin(j,k)=Vin(j,k)*p
     end do
    end do
    !$omp end do
    !$omp do
    do k=1,Prnz
     do j=1,Prny
       Win(j,k)=Win(j,k)*p
     end do
    end do
    !$omp end do

    !$omp sections
    !$omp section
    call BoundUin(1,Uin)
    !$omp section
    call BoundUin(1,Vin)
    !$omp section
    call BoundUin(1,Win)
    !$omp end sections
    !$omp end parallel

    do k=-bigNz+1,Prnz+bigNz
     do j=-bigNy+1,Prny+bigNy
       call RandomGauss(Ru(j,k,2))
       call RandomGauss(Rv(j,k,2))
       call RandomGauss(Rw(j,k,2))
     enddo
    enddo


    if ((Btype(So)==PERIODIC).or.(Btype(No)==PERIODIC)) then
      !$omp parallel workshare
      forall(k=-bigNz+1:Prnz+bigNz,j=-bigNy+1:0)
         Ru(j,k,1:2)=Ru(j+Prny,k,1:2)
         Rv(j,k,1:2)=Rv(j+Prny,k,1:2)
         Rw(j,k,1:2)=Rw(j+Prny,k,1:2)
      endforall
      forall(k=-bigNz+1:Prnz+bigNz,j=Prny+1:Prny+bigNy)
         Ru(j,k,1:2)=Ru(j-Prny,k,1:2)
         Rv(j,k,1:2)=Rv(j-Prny,k,1:2)
         Rw(j,k,1:2)=Rw(j-Prny,k,1:2)
      endforall
      !$omp end parallel workshare
    endif
    if  ((Btype(Bo)==PERIODIC).or.(Btype(To)==PERIODIC)) then
      !$omp parallel workshare
      forall(k=-bigNz+1:0,j=-bigNy+1:Prny+bigNy)
         Ru(j,k,1:2)=Ru(j,k+Prnz,1:2)
         Rv(j,k,1:2)=Rv(j,k+Prnz,1:2)
         Rw(j,k,1:2)=Rw(j,k+Prnz,1:2)
      endforall
      forall(k=Prnz+1:Prnz+bigNz,j=-bigNy+1:Prny+bigNy)
         Ru(j,k,1:2)=Ru(j,k-Prnz,1:2)
         Rv(j,k,1:2)=Rv(j,k-Prnz,1:2)
         Rw(j,k,1:2)=Rw(j,k-Prnz,1:2)
      endforall
      !$omp end parallel workshare
    endif

 endsubroutine GETTURBINLETXIE












  subroutine INITTURBULENCEPROFILES
  integer j,k
   allocate(ustarinlet(1:Prnz))
   allocate(transformtensor(1:6,1:Prny,1:Prnz))
  !constant stress assumption
    if ((profiletype==LOGPROF.and.ustarsurfin<=0).or.(profiletype==POWERPROF.and.ustarsurfin<=0)) then
     ustarsurfin=abs(Karman*Urefin/log(z0inlet/zrefin))
    endif


     do k=1,Prnz
      ustarinlet(k)=ustarsurfin*SQRT(max(1+stressgradin*zPr(k),1E-5_knd))
     enddo

!    !Stress tensor divided by ustar**2
!    relativestress(1,1)=2.4**2 !Garrat
!    relativestress(2,2)=1.9**2
!    relativestress(3,3)=1.25**2
!    relativestress(1,3)=-1
!    relativestress(3,1)=-1
!    relativestress(2,3)=0 !No corriolis effects
!    relativestress(3,2)=0
!    relativestress(1,2)=0
!    relativestress(2,1)=0
   do k=1,Prnz
    do j=1,Prny  ! tt1=a11,tt2=a21,tt3=a22, tt4=a31, tt5=a32, tt6=a33
       transformtensor(1,j,k)=SQRT((ustarinlet(k)**2)*relativestress(1,1))
       transformtensor(2,j,k)=(ustarinlet(k)**2)*relativestress(2,1)/transformtensor(1,j,k)
       transformtensor(3,j,k)=SQRT((ustarinlet(k)**2)*relativestress(2,2)-transformtensor(2,j,k)**2)
       transformtensor(4,j,k)=(ustarinlet(k)**2)*relativestress(3,1)/transformtensor(1,j,k)
       transformtensor(5,j,k)=((ustarinlet(k)**2)*relativestress(3,2)-&
                                transformtensor(2,j,k)*transformtensor(4,j,k))/transformtensor(3,j,k)
       transformtensor(6,j,k)=SQRT((ustarinlet(k)**2)*relativestress(3,3)-&
                                transformtensor(4,j,k)**2-transformtensor(5,j,k)**2)

    enddo
   enddo
  endsubroutine INITTURBULENCEPROFILES

  subroutine INITMEANPROFILES
  integer j,k

  allocate(Uinavg(-2:Uny+3,-2:Unz+3),Vinavg(-2:Vny+3,-2:Vnz+3),Winavg(-2:Wny+3,-2:Wnz+3))
  Vinavg=0
  Winavg=0

  if  (profiletype==CONSTPROF) then
   do k=1,Prnz
    do j=1,Prny
     Uinavg(j,k)=Uinlet
    enddo
   enddo
  elseif (profiletype==LOGPROF) then
   do j=1,Prny
    Uinavg(j,1)=(ustarinlet(1)/Karman)*log(zPr(1)/z0inlet)
   enddo
   do k=2,Prnz
    do j=1,Prny
     Uinavg(j,k)=(ustarinlet(k)/Karman)*log(zPr(k)/zPr(k-1))+Uinavg(j,k-1)
    enddo
   enddo
  elseif (profiletype==POWERPROF) then
   do k=1,Prnz
    do j=1,Prny
     Uinavg(j,k)=Urefin*(zPr(k)/zrefin)**powerexpin
    enddo
   enddo
  else
   Uinavg=0
  endif
  if (windangle/=0) then
   do k=1,Prnz
    do j=1,Prny
     Vinavg(j,k)=Uinavg(j,k)*sin(windangle/180._knd*pi)
     Uinavg(j,k)=Uinavg(j,k)*cos(windangle/180._knd*pi)
    enddo
   enddo
  endif
  endsubroutine INITMEANPROFILES

  subroutine RANDOMGAUSS(res)  !more effective way: http://www.taygeta.com/random/gaussian.html
  real(knd), intent(out) :: res
  real(knd) p,S
  integer,parameter:: n=6
  integer i

  S=0
  do i=1,n
   call RANDOM_NUMBER(p)
   S=S+p
  enddo
  S=S-0.5*n
  S=S/SQRT(n/12._knd)
  res=(p-0.5)*sqrt(12._knd)!S

  end subroutine RANDOMGAUSS


  subroutine GETINLETFROMFILE(t)
  real(TIM),intent(in):: t
  integer,save:: called=0
  integer Prny2,Prnz2,Vny2,Wnz2
  real(knd) dx2
  character(12):: fname
  integer,save:: inletfnum

  type(TInlet),pointer,save:: In1=>null(),In2=>null(),Inp=>null() !pointer to inlets to avoid transfers, time(In1)<time(In2)
  real(TIM),save:: t1,t2 !time if file1, file 2
  real(TIM) tp
  integer,save:: io

  real(knd) c1,c2

  if (called==0) then
     open(102,file="inletframeinfo.unf",form='unformatted',status='old',action='read',iostat=io)
     if (io/=0) then
      write(*,*) 'Error while opening file inletframeinfo.unf'
     endif
     read(102) Prny2,Prnz2  !for check of consistency of grids before use
     read(102) Vny2
     read(102) Wnz2
     read(102) dx2

     allocate(In1,In2)
     allocate(In1%U(Uny,Unz),In1%V(Vny,Vnz),In1%W(Wny,Wnz))
     allocate(In2%U(Uny,Unz),In2%V(Vny,Vnz),In2%W(Wny,Wnz))
     if (enable_buoyancy==1) allocate(In1%temperature(Prny,Prnz),In2%temperature(Prny,Prnz))

     if ((Prny/=Prny2).or.(Prnz/=Prnz2).or.(Vny/=Vny2).or.(Wnz/=Wnz2).or.((dx2-dxPr(0))/dx2>0.1)) then
      write(*,*) "Mismatch of computational grid and inlet file."
      stop
     endif
     called=1
     inletfnum=1

     fname(1:5)="frame"
     write(fname(6:8),"(I3.3)") inletfnum
     fname(9:12)=".unf"

     open(11,file=fname,form='unformatted',status='old',action='read',iostat=io)
     if (io/=0) then
      write(*,*) "Error while opening file ", fname
      stop
     endif
     read(102) t1
     call READINLETFROMFILE(11,In1)
     close(11)

     inletfnum=inletfnum+1

     fname(1:5)="frame"
     write(fname(6:8),"(I3.3)") inletfnum
     fname(9:12)=".unf"

      open(11,file=fname,form='unformatted',status='old',action='read',iostat=io)
     if (io/=0) then
      write(*,*) "Error while opening file ", fname
      stop
     endif
     read(102) t2
     call READINLETFROMFILE(11,In2)
     close(11)
  endif

  io=0

  do while (t2<t.and.io==0)
    read(102,iostat=io) tp
    if (io==0) then
     t1=t2
     t2=tp

     inletfnum=inletfnum+1
     fname(1:5)="frame"
     write(fname(6:8),"(I3.3)") inletfnum
     fname(9:12)=".unf"

     open(11,file=fname,form='unformatted',status='old',action='read',iostat=io)
     if (io/=0) then
      write(*,*) "Error while opening file ", fname
      stop
     endif
     Inp=>In1
     In1=>In2
     In2=>Inp
     call READINLETFROMFILE(11,In2)
     close(11)
    endif
  enddo


  if (t<=t1) then
   c1=1
   c2=0
  elseif (t>=t2) then
   c1=0
   c2=1
  else
   c1=(t2-t)/(t2-t1)
   c2=(t-t1)/(t2-t1)
  endif
  write(*,*) "t",t1,t2,t
  write(*,*) "c",c1,c2

  Uin(1:Uny,1:Unz)=c1*In1%U+c2*In2%U
  Vin(1:Vny,1:Vnz)=c1*In1%V+c2*In2%V
  Win(1:Wny,1:Wnz)=c1*In1%W+c2*In2%W
  if (enable_buoyancy==1) Tempin(1:Prny,1:Prnz)=c1*In1%temperature+c2*In2%temperature

  endsubroutine GETINLETFROMFILE

  subroutine  READINLETFROMFILE(unitnum,In)
  integer,intent(in):: unitnum
  type(TInlet),intent(inout)::In

  read(unitnum) In%U
  read(unitnum) In%V
  read(unitnum) In%W
  if (enable_buoyancy>0) then
       read(unitnum) In%temperature
  endif
  endsubroutine READINLETFROMFILE


  pure subroutine BoundUin(component,Uin)
    integer,  intent(in)    :: component
    real(knd),intent(inout) :: Uin(-2:,-2:)
    integer :: nx, ny, nz

    if (component==1) then
      nx=Unx
      ny=Uny
      nz=Unz
    elseif (component==2) then
      nx=Vnx
      ny=Vny
      nz=Vnz
    else
      nx=Wnx
      ny=Wny
      nz=Wnz
    endif

    if (Btype(Bo)==DIRICHLET) then
      if (component==3) then
        Uin(1:ny,0) = sideU(component,Bo)
        Uin(1:ny,-1) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,1))
      else
        Uin(1:ny,0) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,1))
        Uin(1:ny,-1) = sideU(component,Bo)+(sideU(component,Bo)-Uin(1:ny,2))
      end if
    else if (Btype(Bo)==PERIODIC) then
        Uin(1:ny,-1:0) = Uin(1:ny,nz-1:nz)
    elseif (Btype(Bo)==NOSLIP.or.(component==3.and.Btype(Bo)==FREESLIP)) then
      if (component==3) then
        Uin(1:ny,0)=0
        Uin(1:ny,-1)=-Uin(1:ny,1)
      else
        Uin(1:ny,-1:0)=-Uin(1:ny,2:1:-1)
      endif
    elseif (Btype(Bo)==NEUMANN.or.(component/=3.and.Btype(Bo)==FREESLIP)) then
        Uin(1:ny,0)=Uin(1:ny,1)
        Uin(1:ny,-1)=Uin(1:ny,1)
    elseif (Btype(Bo)==PERIODIC) then
        Uin(1:ny,-1:0) = Uin(1:ny,nz-1:nz)
    endif

    if (Btype(To)==DIRICHLET) then
      if (component==3) then
        Uin(1:ny,nz+1) = sideU(component,To)
        Uin(1:ny,nz+2) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz))
      else
        Uin(1:ny,nz+1) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz))
        Uin(1:ny,nz+2) = sideU(component,To)+(sideU(component,To)-Uin(1:ny,nz-1))
      end if
    else if (Btype(To)==PERIODIC) then
        Uin(1:ny,nz+1:nz+2) = Uin(1:ny,nz-1:nz)
    elseif (Btype(To)==NOSLIP.or.(component==3.and.(Btype(To)==FREESLIP.or.Btype(To)==FREESLIPBUFF))) then
      if (component==3) then
        Uin(1:ny,nz+1)=0
        Uin(1:ny,nz+2)=-Uin(1:ny,nz)
      else
        Uin(1:ny,nz+1:nz+2)=-Uin(1:ny,nz:nz-1:-1)
      endif
    elseif (Btype(To)==NEUMANN.or.(component/=3.and.(Btype(To)==FREESLIP.or.Btype(To)==FREESLIPBUFF))) then
        Uin(1:ny,nz+1)=Uin(1:ny,nz)
        Uin(1:ny,nz+2)=Uin(1:ny,nz)
    elseif (Btype(To)==PERIODIC) then
        Uin(1:ny,nz+1:nz+2) = Uin(1:ny,1:2)
    endif

    if (Btype(So)==DIRICHLET) then
      if (component==2) then
        Uin(0,-1:nz+2)=sideU(component,So)
        Uin(-1,-1:nz+2)=sideU(component,So)+(sideU(component,So)-Uin(1,-1:nz+2))
      else
        Uin(0,-1:nz+2)=sideU(component,So)+(sideU(component,So)-Uin(1,-1:nz+2))
        Uin(-1,-1:nz+2)=sideU(component,So)+(sideU(component,So)-Uin(2,-1:nz+2))
      endif
    elseif (Btype(So)==NOSLIP.or.(component==2.and.Btype(So)==FREESLIP)) then
      if (component==2) then
        Uin(0,-1:nz+2)=0
        Uin(-1,-1:nz+2)=-Uin(1,-1:nz+2)
      else
        Uin(0,-1:nz+2)=-Uin(1,-1:nz+2)
        Uin(-1,-1:nz+2)=-Uin(2,-1:nz+2)
      endif
    elseif (Btype(So)==NEUMANN.or.(component/=2.and.Btype(So)==FREESLIP)) then
        Uin(0,-1:nz+2)=Uin(1,-1:nz+2)
        Uin(-1,-1:nz+2)=Uin(1,-1:nz+2)
    elseif (Btype(So)==PERIODIC) then  !Periodic BC
        Uin(0,-1:nz+2)=Uin(ny,-1:nz+2)
        Uin(-1,-1:nz+2)=Uin(ny-1,-1:nz+2)
    endif

    if (Btype(No)==DIRICHLET) then
      if (component==2) then
        Uin(ny+1,-1:nz+2)=sideU(component,No)
        Uin(ny+2,-1:nz+2)=sideU(component,No)+(sideU(component,No)-Uin(ny,-1:nz+2))
      else
        Uin(ny+1,-1:nz+2)=sideU(component,No)+(sideU(component,No)-Uin(ny,-1:nz+2))
        Uin(ny+2,-1:nz+2)=sideU(component,No)+(sideU(component,No)-Uin(ny-1,-1:nz+2))
      endif
    elseif (Btype(No)==NOSLIP.or.(component==2.and.Btype(So)==FREESLIP)) then
      if (component==2) then
        Uin(ny+1,-1:nz+2)=0
        Uin(ny+2,-1:nz+2)=-Uin(ny,-1:nz+2)
      else
        Uin(ny+1:ny+2,-1:nz+2)=-Uin(ny:ny-1:-1,-1:nz+2)
      endif
    elseif (Btype(No)==NEUMANN.or.(component/=2.and.Btype(No)==FREESLIP)) then
        Uin(ny+1,-1:nz+2)=Uin(ny,-1:nz+2)
        Uin(ny+2,-1:nz+2)=Uin(ny,-1:nz+2)
    elseif (Btype(No)==PERIODIC) then  !Periodic BC
        Uin(ny+1:ny+2,-1:nz+2)=Uin(1:2,-1:nz+2)
    endif

  end subroutine BoundUin





end module TURBINLET

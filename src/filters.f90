module Filters

  use Parameters

  implicit none

  private

  public filtertype, filter_ratios, Filter

  integer filtertype

  real(KND), parameter :: filter_ratios(0:1) = [ 1.0_KND, 2._KND]

 
  contains


    subroutine Filter3ord(U,dir)
      !3D filter with 3 vanishing moments width 2 delta x
      !Vasilyev, Lund, Moin, 1998, http://dx.doi.org/10.1006/jcph.1998.6060
      real(KND),dimension(-2:,-2:,-2:),intent(inout)  :: U
      integer,intent(in) :: dir
      real(KND),parameter :: w(-3:3) = (/ -1._KND/32, 0._KND, 9._KND/32, 1._KND/2, 9._KND/32, 0._KND, -1._KND/32 /)
      integer i,j,k
      integer mini,minj,mink
      integer maxi,maxj,maxk

      real(KND) :: tmp(-2:ubound(U,dir))
      
      mini = lbound(U,1)+3
      maxi = ubound(U,1)-3
      minj = lbound(U,2)+3
      maxj = ubound(U,2)-3
      mink = lbound(U,3)+3
      maxk = ubound(U,3)-3
      
      if (dir==1) then
        !$omp parallel do private(i,j,k,tmp)
        do k = mink,maxk
         do j = minj,maxj
           tmp = U(:,j,k)
           do i = 1,maxi
             U(i,j,k) = sum(tmp(i-3:i+3) * w)
           end do
         end do
        end do
        !$omp end parallel do
      else if (dir==2) then
        !$omp parallel do private(i,j,k,tmp)
        do k = mink,maxk
         do i = mini,maxi
           tmp = U(i,:,k)
           do j = 1,maxj
             U(i,j,k) = sum(tmp(j-3:j+3) * w)
           end do
         end do
        end do
        !$omp end parallel do
      else
        !$omp parallel do private(i,j,k,tmp)
        do j = mink,maxk
         do i = mini,maxi
           tmp = U(i,j,:)
           do k = 1,maxk
             U(i,j,k) = sum(tmp(k-3:k+3) * w)
           end do
         end do
        end do
        !$omp end parallel do
      end if
    end subroutine Filter3ord



    subroutine Filter(U)    !Calls a selected filter  U2 <- filt(U1)
      real(KND),dimension(:,:,:),intent(inout) :: U

      if (filtertype==1) then
        if (Prnx==1) then
          call Filter3ord(U,2)
          call Filter3ord(U,3)
        else if (Prny==1) then
          call Filter3ord(U,1)
          call Filter3ord(U,3)
        else if (Prnz==1) then
          call Filter3ord(U,3)
          call Filter3ord(U,3)
        else
          call Filter3ord(U,1)
          call Filter3ord(U,2)
          call Filter3ord(U,3)
        end if
      end if

    end subroutine Filter

end module Filters

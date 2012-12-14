module Filters

  use Parameters

  implicit none

  private

  public filtertype, filter_ratios, Filter

  integer filtertype

  real(KND), parameter :: filter_ratios(0:1) = [ 1.0_KND, 2._KND]


  interface Filter3ord
    module procedure Filter3ord_2D
    module procedure Filter3ord_3D
  end interface
  
  contains

    subroutine Filter3ord_3D(U2,U1)
      !3D filter with 3 vanishing moments width 2 delta x
      !Vasilyev, Lund, Moin, 1998, http://dx.doi.org/10.1006/jcph.1998.6060
      real(KND),dimension(:,:,:),intent(in)  :: U1
      real(KND),dimension(:,:,:),intent(out) :: U2
      real(KND),parameter :: w(-3:3) = (/ -1._KND/32, 0._KND, 9._KND/32, 1._KND/2, 9._KND/32, 0._KND, -1._KND/32 /)
      real(KND),save :: w3(-3:3,-3:3,-3:3)
      real(KND) S
      integer ii,jj,i,j,k,kk
      integer,save :: called = 0

      if (called ==0) then
        do kk=-3,3
         do jj=-3,3
          do ii=-3,3
            w3(ii,jj,kk) = w(ii)*w(jj)*w(kk)
          end do
         end do
        end do
        
        called = 1
        
      end if

      !$omp parallel do private(i,j,k,ii,jj,kk,S)
      do k = lbound(U1,3)+3,ubound(U1,3)-3
       do j = lbound(U1,2)+3,ubound(U1,2)-3
        do i = lbound(U1,1)+3,ubound(U1,1)-3
        
          S = 0

          do kk = -3,3
           do jj = -3,3
            do ii = -3,3
              S = S + W3(ii,jj,kk) * U1(i+ii,j+jj,k+kk)
            end do
           end do
          end do

          U2(i,j,k) = S
         
        end do
       end do
      end do
      !$omp end parallel do
    end subroutine Filter3ord_3D


    subroutine Filter3ord_2D(U2,U1,dir)
      !3D filter with 3 vanishing moments width 2 delta x
      !Vasilyev, Lund, Moin, 1998, http://dx.doi.org/10.1006/jcph.1998.6060
      real(KND),dimension(-2:,-2:,-2:),intent(in)  :: U1
      real(KND),dimension(-2:,-2:,-2:),intent(out) :: U2
      integer,intent(in) :: dir
      real(KND),parameter :: w(-3:3) = (/ -1._KND/32, 0._KND, 9._KND/32, 1._KND/2, 9._KND/32, 0._KND, -1._KND/32 /)
      real(KND) :: w2(-3:3,-3:3,-3:3)
      real(KND) S
      integer ii,jj,i,j,k,kk
      integer mini,minii,minj,minjj,mink,minkk
      integer maxi,maxii,maxj,maxjj,maxk,maxkk

      minii = -3
      maxii = 3
      minjj = -3
      maxjj = 3
      minkk = -3
      maxkk = 3

      mini = lbound(U1,1)+3
      maxi = ubound(U1,1)-3
      minj = lbound(U1,2)+3
      maxj = ubound(U1,2)-3
      mink = lbound(U1,3)+3
      maxk = ubound(U1,3)-3

      if (dir==1) then
        minii=0
        maxii=0
        mini=1
        maxi=1
      else if (dir==2) then
        minjj=0
        maxjj=0
        minj=1
        maxj=1
      else
        minkk=0
        maxkk=0
        mink=1
        maxk=1
      end if

      do kk = minkk,maxkk
       do jj = minjj,maxjj
        do ii = minii,maxii  
          w2(ii,jj,kk) = 2 * w(ii)*w(jj)*w(kk)         
        end do
       end do
      end do

      !$omp parallel do  private(i,j,k,ii,jj,kk,S)
      do k = mink,maxk
       do j = minj,maxj
        do i = mini,maxi
        
          S = 0
          
          do kk = minkk,maxkk
           do jj = minjj,maxjj
            do ii = minii,maxii
              S = S + W2(ii,jj,kk) * U1(i+ii, j+jj, k+kk)
            end do
           end do
          end do
          
          U2(i,j,k) = S
          
        end do
       end do
      end do
      !$omp end parallel do
    end subroutine Filter3ord_2D



    subroutine Filter(U2,U1)    !Calls a selected filter  U2 <- filt(U1)
      real(KND),dimension(:,:,:),intent(in) :: U1
      real(KND),dimension(:,:,:),intent(out) :: U2

      if (filtertype==1) then
        if (Prnx==1) then
          call Filter3ord(U2,U1,1)
        else if (Prny==1) then
          call Filter3ord(U2,U1,2)
        else if (Prnz==1) then
          call Filter3ord(U2,U1,3)
        else
          call Filter3ord(U2,U1)
        end if
      else   !should not happen, just to be safe
        U2 = U1
      end if

    end subroutine Filter

end module Filters
